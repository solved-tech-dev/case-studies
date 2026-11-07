// SPDX-License-Identifier: LGPL-2.1-or-later
#include "Xianzhi45Adapter.hpp"
#include "XianzhiKrylovAPI.h"
#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <dlfcn.h>
#include <execinfo.h>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <csignal>
#include <unistd.h>

namespace SU2Xianzhi45 {
namespace {
// Diagnostic watchdog: when SU2_HANG_BT=1, a one-shot alarm dumps the current
// stack (including the library frames) and exits if a solve call is stuck.
void hang_backtrace_handler(int) {
  const char* enabled = std::getenv("SU2_HANG_BT");
  if (!enabled || std::strcmp(enabled, "1") != 0) return;
  void* frames[64];
  const int count = backtrace(frames, 64);
  std::fprintf(stderr, "SU2_HANG_BT: watchdog fired, stack trace follows\n");
  backtrace_symbols_fd(frames, count, 2);
  _exit(86);
}

struct HangWatchdog {
  HangWatchdog() {
    const char* enabled = std::getenv("SU2_HANG_BT");
    if (enabled && std::strcmp(enabled, "1") == 0) {
      const char* secs = std::getenv("SU2_HANG_BT_SECS");
      const unsigned delay = secs ? std::atoi(secs) : 25;
      std::signal(SIGALRM, hang_backtrace_handler);
      alarm(delay ? delay : 25);
    }
  }
};
static HangWatchdog g_hang_watchdog;

// Internal library entry points living on the private
// XIANZHI_KRYLOV_INTERNAL_26_45 version node. They reproduce SU2's
// CSysMatrix::HtDTransfer bulk uploads inside the library so the adapter can
// hand the SDK a device-resident BCSR descriptor instead of copying blocks
// through host memory on every linear call.
extern "C" {
int xk_internal_su2_htd_alloc_f64(double** device, std::uint64_t count);
int xk_internal_su2_htd_upload_f64(double* device, const double* host, std::uint64_t count);
int xk_internal_su2_htd_alloc_u64(std::uint64_t** device, std::uint64_t count);
int xk_internal_su2_htd_upload_u64(std::uint64_t* device, const std::uint64_t* host, std::uint64_t count);
int xk_internal_su2_htd_alloc_i32(std::int32_t** device, std::uint64_t count);
int xk_internal_su2_htd_upload_i32(std::int32_t* device, const std::int32_t* host, std::uint64_t count);
void xk_internal_su2_htd_free(void* device);
int xk_internal_get_last_stats(void* handle, std::uint32_t* iterations, double* residual, double* rhs_norm);
}
using Clock = std::chrono::steady_clock;
void check(int status, const char* operation) {
  if (status != XIANZHIKRYLOV_SUCCESS)
    throw std::runtime_error(std::string(operation) + " status=" + std::to_string(status));
}
void trace(const char* operation) {
  std::cerr << "SU2_XIANZHI45 " << operation << std::endl;
}
void verify_library() {
  Dl_info info{};
  char executable[PATH_MAX]{}, actual[PATH_MAX]{}, expected[PATH_MAX]{};
  const auto length = readlink("/proc/self/exe", executable, sizeof(executable) - 1);
  if (length <= 0 || static_cast<size_t>(length) >= sizeof(executable) - 1)
    throw std::runtime_error("Cannot find the SU2 executable path.");
  executable[length] = '\0';
  const auto package = std::filesystem::path(executable).parent_path().parent_path().parent_path();
  const auto expected_library = package / "sdk/lib/libXianzhiKrylov.so.2";
  if (!dladdr(reinterpret_cast<void*>(&xianzhikrylov_create), &info) ||
      !realpath(info.dli_fname, actual) || !realpath(expected_library.c_str(), expected) ||
      std::strcmp(actual, expected))
    throw std::runtime_error("The loaded XianzhiKrylov library is not the package library.");
  static bool printed = false;
  if (!printed) {
    std::cerr << "SU2_XIANZHI45_LIBRARY " << actual << std::endl;
    printed = true;
  }
}
double norm(const std::vector<double>& vector) {
  long double sum = 0;
  for (double value : vector) {
    if (!std::isfinite(value)) throw std::runtime_error("Non-finite vector data.");
    sum += static_cast<long double>(value) * value;
  }
  return std::sqrt(sum);
}
double residual(const Matrix& matrix, const std::vector<double>& rhs, const std::vector<double>& x) {
  auto ax = product(matrix, x);
  for (std::size_t i = 0; i < ax.size(); ++i) ax[i] -= rhs[i];
  return norm(ax);
}
void validate(const Matrix& matrix, std::size_t length) {
  const auto bs = matrix.block_size;
  if (!bs || matrix.rows.size() < 2 || matrix.rows.front() != 0 ||
      matrix.rows.back() != matrix.columns.size() ||
      (matrix.rows.size() - 1) * bs != length ||
      !matrix.values ||
      matrix.value_count != matrix.columns.size() * std::size_t(bs) * bs)
    throw std::runtime_error("Invalid native block CSR dimensions.");
  const auto n = matrix.rows.size() - 1;
  for (std::size_t row = 0; row < n; ++row) {
    if (matrix.rows[row] >= matrix.rows[row + 1]) throw std::runtime_error("Empty or unordered block row.");
    bool diagonal = false;
    for (auto k = matrix.rows[row]; k < matrix.rows[row + 1]; ++k) {
      if (matrix.columns[k] >= n || (k > matrix.rows[row] && matrix.columns[k - 1] >= matrix.columns[k]))
        throw std::runtime_error("Invalid or unsorted native block column.");
      diagonal |= matrix.columns[k] == row;
    }
    if (!diagonal) throw std::runtime_error("Native block row has no diagonal.");
  }
  const double* finite_values = matrix.host_values ? matrix.host_values : matrix.values;
  for (std::uint64_t k = 0; k < matrix.value_count; ++k)
    if (!std::isfinite(finite_values[k])) throw std::runtime_error("Non-finite native matrix value.");
}
std::string statistics(xianzhikrylov_handle_t handle) {
  std::cout.flush();
  std::fflush(stdout);
  FILE* file = std::tmpfile();
  if (!file) throw std::runtime_error("Cannot allocate statistics capture.");
  const int saved = dup(STDOUT_FILENO);
  if (saved < 0 || dup2(fileno(file), STDOUT_FILENO) < 0) {
    if (saved >= 0) close(saved);
    std::fclose(file);
    throw std::runtime_error("Cannot capture SDK statistics.");
  }
  const int status = xianzhikrylov_print_info(handle);
  std::cout.flush();
  std::fflush(stdout);
  const int restored = dup2(saved, STDOUT_FILENO);
  close(saved);
  std::rewind(file);
  std::string output;
  char data[4096];
  while (const auto count = std::fread(data, 1, sizeof(data), file)) output.append(data, count);
  std::fclose(file);
  if (restored < 0) throw std::runtime_error("Cannot restore stdout.");
  check(status, "print_info");
  return output;
}
double number(const std::string& stats, const char* name) {
  std::smatch match;
  if (!std::regex_search(stats, match, std::regex(std::string("(?:^|[ \\n])") + name + "=([-+0-9.eE]+)")))
    throw std::runtime_error(std::string("Missing SDK statistic: ") + name);
  return std::stod(match[1]);
}
struct Session {
  Matrix matrix;
  std::vector<int32_t> rows32, columns32;
  std::vector<uint64_t> diagonal;
  std::vector<double> rhs, x;
  double* device_values = nullptr;
  uint64_t* device_rows = nullptr;
  uint64_t* device_columns = nullptr;
  uint64_t* device_diagonal = nullptr;
  int32_t* device_rows32 = nullptr;
  int32_t* device_columns32 = nullptr;
  xianzhikrylov_handle_t handle = nullptr;
  xianzhikrylov_matrix_t descriptor = nullptr;
  xianzhikrylov_vector_t rhs_object = nullptr, x_object = nullptr;
  unsigned long calls = 0;
  bool configured = false;
  bool structure_ready = false;
  bool numeric_ready = false;
  int device = -1;
  std::string config_path;
  ~Session() {
    if (handle) xianzhikrylov_destroy(handle);
    if (x_object) xianzhikrylov_vector_destroy(x_object);
    if (rhs_object) xianzhikrylov_vector_destroy(rhs_object);
    if (descriptor) xianzhikrylov_matrix_destroy(descriptor);
    xk_internal_su2_htd_free(device_values);
    xk_internal_su2_htd_free(device_rows);
    xk_internal_su2_htd_free(device_columns);
    xk_internal_su2_htd_free(device_diagonal);
    xk_internal_su2_htd_free(device_rows32);
    xk_internal_su2_htd_free(device_columns32);
  }
};
std::map<const void*, std::unique_ptr<Session>> sessions;
std::mutex session_mutex;
void initialize(Session& session, const Matrix& matrix, const std::vector<double>& rhs,
                const std::vector<double>& x, const Options& options) {
  session.matrix = matrix;
  session.rhs = rhs;
  session.x = x;
  session.device = options.device;
  const auto bs = matrix.block_size;
  const auto n = matrix.rows.size() - 1;
  std::filesystem::create_directories(options.output_directory);
  std::ostringstream filename;
  filename << "linear-" << reinterpret_cast<std::uintptr_t>(&session) << ".json";
  session.config_path = (std::filesystem::path(options.output_directory) / filename.str()).string();
  if (bs == 1) {
    if (n > INT_MAX || matrix.columns.size() > INT_MAX) throw std::runtime_error("CSR exceeds public SDK int32 limits.");
    session.rows32.assign(matrix.rows.begin(), matrix.rows.end());
    session.columns32.assign(matrix.columns.begin(), matrix.columns.end());
  } else {
    for (std::size_t row = 0; row < n; ++row)
      for (auto k = matrix.rows[row]; k < matrix.rows[row + 1]; ++k)
        if (matrix.columns[k] == row) { session.diagonal.push_back(k); break; }
  }
  const auto nnz = matrix.columns.size();
  const bool borrow = bs > 1 && matrix.device_rows != nullptr &&
      matrix.device_columns != nullptr;
  // Upload the SU2 contiguous block array once per topology; later calls only
  // refresh the value buffer at the same device address.
  if (bs == 1) {
    check(xk_internal_su2_htd_alloc_i32(&session.device_rows32, n + 1), "htd_rows32_alloc");
    check(xk_internal_su2_htd_upload_i32(session.device_rows32, session.rows32.data(), n + 1), "htd_rows32_upload");
    check(xk_internal_su2_htd_alloc_i32(&session.device_columns32, nnz), "htd_columns32_alloc");
    check(xk_internal_su2_htd_upload_i32(session.device_columns32, session.columns32.data(), nnz), "htd_columns32_upload");
  } else {
    if (!borrow) {
      check(xk_internal_su2_htd_alloc_u64(&session.device_rows, n + 1), "htd_rows_alloc");
      check(xk_internal_su2_htd_upload_u64(session.device_rows, matrix.rows.data(), n + 1), "htd_rows_upload");
      check(xk_internal_su2_htd_alloc_u64(&session.device_columns, nnz), "htd_columns_alloc");
      check(xk_internal_su2_htd_upload_u64(session.device_columns, matrix.columns.data(), nnz), "htd_columns_upload");
    }
    check(xk_internal_su2_htd_alloc_u64(&session.device_diagonal, n), "htd_diagonal_alloc");
    check(xk_internal_su2_htd_upload_u64(session.device_diagonal, session.diagonal.data(), n), "htd_diagonal_upload");
  }
  if (!borrow) {
    check(xk_internal_su2_htd_alloc_f64(&session.device_values, matrix.value_count), "htd_values_alloc");
    check(xk_internal_su2_htd_upload_f64(session.device_values, matrix.values, matrix.value_count), "htd_values_upload");
  }
  trace("create_begin");
  check(xianzhikrylov_create(options.device, &session.handle), "create");
  trace("create_end");
  xianzhikrylov_matrix_data_t d{};
  d.abi_version = XIANZHIKRYLOV_ABI_VERSION;
  d.struct_size = sizeof(d);
  d.format = bs == 1 ? XIANZHIKRYLOV_CSR : XIANZHIKRYLOV_BSR;
  d.value_type = XIANZHIKRYLOV_FP64;
  d.index_type = bs == 1 ? XIANZHIKRYLOV_INDEX_I32 : XIANZHIKRYLOV_INDEX_I64;
  d.block_size = bs;
  d.block_layout = XIANZHIKRYLOV_BLOCK_ROW_MAJOR;
  d.memory_space = XIANZHIKRYLOV_DEVICE;
  d.device_id = options.device;
  d.rows = d.columns = n;
  d.entries = nnz;
  d.row_offsets = {bs == 1 ? static_cast<void*>(session.device_rows32)
                   : static_cast<void*>(borrow ? const_cast<std::uint64_t*>(matrix.device_rows) : session.device_rows),
                   (bs == 1 ? n + 1 : matrix.rows.size()) * (bs == 1 ? 4 : 8)};
  d.column_indices = {bs == 1 ? static_cast<void*>(session.device_columns32)
                      : static_cast<void*>(borrow ? const_cast<std::uint64_t*>(matrix.device_columns) : session.device_columns),
                      nnz * (bs == 1 ? 4 : 8)};
  if (bs > 1) d.diagonal_positions = {session.device_diagonal, session.diagonal.size() * 8};
  d.values = {borrow ? const_cast<double*>(matrix.values) : session.device_values, matrix.value_count * 8};
  check(xianzhikrylov_matrix_create(&session.descriptor, &d), "matrix_create");
  check(xianzhikrylov_vector_create(&session.rhs_object, rhs.size(), session.rhs.data(), XIANZHIKRYLOV_FP64, XIANZHIKRYLOV_HOST), "rhs_create");
  check(xianzhikrylov_vector_create(&session.x_object, x.size(), session.x.data(), XIANZHIKRYLOV_FP64, XIANZHIKRYLOV_HOST), "solution_create");
}
}

std::vector<double> product(const Matrix& matrix, const std::vector<double>& x) {
  const auto bs = matrix.block_size;
  const double* values = matrix.host_values ? matrix.host_values : matrix.values;
  std::vector<double> output(x.size(), 0);
  for (std::size_t row = 0; row + 1 < matrix.rows.size(); ++row)
    for (std::uint32_t i = 0; i < bs; ++i) {
      long double sum = 0;
      for (auto k = matrix.rows[row]; k < matrix.rows[row + 1]; ++k)
        for (std::uint32_t j = 0; j < bs; ++j)
          sum += static_cast<long double>(values[(k * bs + i) * bs + j]) * x[matrix.columns[k] * bs + j];
      output[row * bs + i] = sum;
    }
  return output;
}

Result solve(const void* owner, const Matrix& matrix, const std::vector<double>& rhs,
             std::vector<double>& solution, const Options& options) {
  std::lock_guard<std::mutex> lock(session_mutex);
  verify_library();
  if (!owner || rhs.empty() || rhs.size() != solution.size() || options.output_directory.empty() ||
      !std::isfinite(options.tolerance) || options.tolerance <= 0 || !options.max_iterations)
    throw std::runtime_error("Invalid adapter input/options.");
  auto& entry = sessions[owner];
  // SU2 keeps one CSysMatrix per solver and only refreshes coefficient values
  // between outer iterations, so the topology stays constant. Check its size
  // in O(1) and run the full host-side validation only when the session is
  // created; scanning the whole value array on every call was the dominant
  // per-call host cost for large blocks.
  const bool same_structure = entry &&
      entry->matrix.block_size == matrix.block_size &&
      entry->matrix.rows.size() == matrix.rows.size() &&
      (matrix.rows.size() < 2 || entry->matrix.rows.back() == matrix.rows.back()) &&
      entry->device == options.device;
  if (entry && !same_structure) entry.reset();
  if (!entry) {
    validate(matrix, rhs.size());
    auto created = std::make_unique<Session>();
    initialize(*created, matrix, rhs, solution, options);
    entry = std::move(created);
  }
  auto& session = *entry;
  Result result;
  const bool borrow = matrix.device_rows != nullptr;
  if (!borrow) {
    const auto upload_started = Clock::now();
    check(xk_internal_su2_htd_upload_f64(session.device_values, matrix.values, matrix.value_count),
          "htd_values_update");
    result.upload_seconds = std::chrono::duration<double>(Clock::now() - upload_started).count();
  }
  std::copy(rhs.begin(), rhs.end(), session.rhs.begin());
  if (options.zero_initial_guess) std::fill(session.x.begin(), session.x.end(), 0.0);
  else std::copy(solution.begin(), solution.end(), session.x.begin());
  result.rhs_norm = norm(rhs);
  // The initial residual only matters when the stopping test is relative to
  // the initial guess; otherwise this host matvec is wasted work per call.
  result.initial_norm = (options.relative_to_initial && !options.zero_initial_guess)
      ? residual(matrix, rhs, session.x) : result.rhs_norm;
  result.normalizer = options.relative_to_initial ? result.initial_norm : result.rhs_norm;
  result.target = options.tolerance * result.normalizer;
  if (result.initial_norm == 0) {
    result.sdk_status = XIANZHIKRYLOV_CONVERGED;
    std::copy(session.x.begin(), session.x.end(), solution.begin());
    return result;
  }
  // Preserve SU2's stopping threshold when it is relative to the initial residual.
  const double sdk_rtol = result.rhs_norm > 0
      ? ((!options.relative_to_initial || options.zero_initial_guess) ? options.tolerance : result.target / result.rhs_norm)
      : 0.0;
  const double sdk_atol = result.rhs_norm > 0 ? 0.0 : result.target;
  if (!session.configured) {
    const int level = options.preconditioner_level;
    std::ofstream config(session.config_path);
    config << std::setprecision(17)
           << "{\"schema_version\":1,\"matrix_format\":\"" << (matrix.block_size == 1 ? "csr" : "bsr")
           << "\",\"matrix_symmetry\":\"general\",\"block_size\":" << matrix.block_size
           << ",\"solver\":{\"method\":\"fgmres\",\"level\":" << level
           << ",\"max_iterations\":" << options.max_iterations << ",\"restart\":" << options.max_iterations
           << ",\"relative_tolerance\":" << sdk_rtol << ",\"absolute_tolerance\":" << sdk_atol
           << ",\"initial_guess\":\"" << (options.zero_initial_guess ? "zero" : "provided") << "\"}";
    config << "}\n";
    config.close();
    if (!config) throw std::runtime_error("Cannot write SDK configuration.");
    check(xianzhikrylov_config_set(session.handle, XIANZHIKRYLOV_CONFIG_FILE, session.config_path.c_str()), "config_set");
    session.configured = true;
  }
  trace("all_begin");
  const auto started = Clock::now();
  // SU2 reuses one CSysMatrix per solver across outer iterations: only the
  // coefficient values change. Build the symbolic structure once, then refresh
  // values + factor and solve. This avoids rebuilding the CSR/BSR structure,
  // re-uploading row/column arrays and rerunning symbolic analysis on every
  // linear system, which dominated the execute(ALL) time.
  if (!session.structure_ready) {
    const auto analysis_started = Clock::now();
    trace("analysis_begin");
    check(xianzhikrylov_execute(session.handle, XIANZHIKRYLOV_PHASE_ANALYSIS,
                                session.descriptor, nullptr, nullptr), "analysis");
    trace("analysis_end");
    result.analysis_seconds = std::chrono::duration<double>(Clock::now() - analysis_started).count();
    session.structure_ready = true;
  }
  // The archived v26.42.13 client built the symbolic structure and factor
  // once per topology and only refreshed the RHS per linear call.  Rebuilding
  // the ILU factors on every call re-pays the level-scheduled factorization
  // and is what makes the 26.45 adapter slower than the archived run.
  // SU2 rebuilds the coefficient values on every outer iteration, so the ILU
  // factors must be refreshed each linear call to stay numerically valid.
  // Default to refresh-every-call; SU2_XK_NUMERIC_REFRESH=0 opts into the
  // faster stale-factor mode for benchmarks that can accept the residual loss.
  const char* refresh_env = std::getenv("SU2_XK_NUMERIC_REFRESH");
  const long refresh_every = refresh_env ? std::atol(refresh_env) : 1;
  if (!session.numeric_ready || (refresh_every > 0 &&
                                 session.calls % refresh_every == 0)) {
    const auto numeric_started = Clock::now();
    trace("numeric_begin");
    check(xianzhikrylov_execute(session.handle, XIANZHIKRYLOV_PHASE_NUMERIC_PREPARE,
                                session.descriptor, nullptr, nullptr), "numeric_prepare");
    trace("numeric_end");
    result.numeric_prepare_seconds = std::chrono::duration<double>(Clock::now() - numeric_started).count();
    session.numeric_ready = true;
  }
  const auto solve_started = Clock::now();
  trace("solve_begin");
  result.sdk_status = xianzhikrylov_execute(session.handle, XIANZHIKRYLOV_PHASE_SOLVE,
                                            session.descriptor, session.x_object, session.rhs_object);
  trace("solve_end");
  result.solve_seconds = std::chrono::duration<double>(Clock::now() - solve_started).count();
  result.execute_seconds = std::chrono::duration<double>(Clock::now() - started).count();
  trace("all_end");
  if (result.sdk_status != XIANZHIKRYLOV_CONVERGED && result.sdk_status != XIANZHIKRYLOV_MAX_ITERATIONS)
    throw std::runtime_error("execute(ALL) status=" + std::to_string(result.sdk_status));
  // Read the last-result counters through the private internal accessor
  // instead of capturing print_info output on every call.
  std::uint32_t iterations = 0;
  double reported_residual = 0, reported_rhs_norm = 0;
  check(xk_internal_get_last_stats(session.handle, &iterations, &reported_residual, &reported_rhs_norm), "get_stats");
  if (iterations > options.max_iterations) throw std::runtime_error("Invalid SDK iteration count.");
  result.iterations = iterations;
  result.residual_norm = reported_residual;
  if (!std::isfinite(result.residual_norm)) throw std::runtime_error("Non-finite SDK residual.");
  norm(session.x);
  const char* verify_every_env = std::getenv("SU2_XK_VERIFY_EVERY");
  const int verify_every = verify_every_env ? std::atoi(verify_every_env) : 1;
  const bool do_verify = options.verify && (verify_every <= 1 || session.calls % verify_every == 0);
  if (do_verify) {
    const auto verify_started = Clock::now();
    // One host matvec serves both the independent residual and the wrapper's
    // native-product cross check; the previous code computed it twice.
    result.matvec_product = product(matrix, session.x);
    long double sum = 0;
    for (std::size_t i = 0; i < result.matvec_product.size(); ++i) {
      const long double difference =
          static_cast<long double>(result.matvec_product[i]) - rhs[i];
      sum += difference * difference;
    }
    const double actual = std::sqrt(sum);
    const double rounding = 100 * std::numeric_limits<double>::epsilon() * std::max(result.rhs_norm, result.initial_norm);
    if (result.sdk_status == XIANZHIKRYLOV_CONVERGED && actual > result.target * (1 + 1e-6) + rounding)
      throw std::runtime_error("The independent native-matrix residual fails the requested tolerance.");
    result.residual_norm = actual;
    result.verify_seconds = std::chrono::duration<double>(Clock::now() - verify_started).count();
  }
  std::copy(session.x.begin(), session.x.end(), solution.begin());
  std::ofstream log(std::filesystem::path(options.output_directory) / "linear-solves.jsonl", std::ios::app);
  log << std::setprecision(17) << "{\"call\":" << session.calls++ << ",\"block_size\":" << matrix.block_size
      << ",\"iterations\":" << result.iterations << ",\"sdk_status\":" << result.sdk_status
      << ",\"execute_seconds\":" << result.execute_seconds << ",\"residual\":" << result.residual_norm
      << ",\"upload_seconds\":" << result.upload_seconds << ",\"stats_seconds\":" << result.stats_seconds
      << ",\"analysis_seconds\":" << result.analysis_seconds
      << ",\"numeric_prepare_seconds\":" << result.numeric_prepare_seconds
      << ",\"solve_seconds\":" << result.solve_seconds
      << ",\"verify_seconds\":" << result.verify_seconds << ",\"target\":" << result.target << ",\"phase\":\"ALL\"}\n";
  std::cout << "SU2_XIANZHI45_SOLVE block_size=" << matrix.block_size << " iterations=" << result.iterations
            << " status=" << result.sdk_status << " residual=" << result.residual_norm << std::endl;
  return result;
}

void clear() { std::lock_guard<std::mutex> lock(session_mutex); sessions.clear(); }
}
