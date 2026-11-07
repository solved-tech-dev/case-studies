// SPDX-License-Identifier: LGPL-2.1-or-later
// This translation unit is added at final link time with --wrap. Original
// SU2 files, class layout, and native CPU implementation remain unchanged.
#include "CSysSolve.hpp"
#define private public
#include "CSysMatrix.hpp"
#undef private
#include "CMatrixVectorProduct.hpp"
#include "Xianzhi45Adapter.hpp"
#include "wrapped_symbol.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <map>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
struct Verification {
  SU2Xianzhi45::Result result;
  std::vector<double> expected_product;
};
struct StructureCache {
  std::vector<unsigned long> rows;
  std::vector<unsigned long> columns;
};
struct CallTimer {
  const CSysSolve<double>* solver;
  std::string backend;
  unsigned long block_size;
  std::chrono::steady_clock::time_point started = std::chrono::steady_clock::now();
  ~CallTimer() {
    if (omp_get_thread_num() != 0) return;
    const char* directory = std::getenv("SU2_XIANZHI_OUTPUT_DIR");
    if (!directory) return;
    const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
    std::ofstream log(std::string(directory) + "/whole-linear-calls.jsonl", std::ios::app);
    log << std::setprecision(17) << "{\"backend\":\"" << backend << "\",\"block_size\":" << block_size
        << ",\"iterations\":" << solver->GetIterations() << ",\"residual\":" << solver->GetResidual()
        << ",\"wall_seconds\":" << wall << "}\n";
  }
};
std::map<const void*, Verification> verification;
std::map<const CSysMatrix<double>*, StructureCache> structure_cache;
bool verify_native() {
  const char* value = std::getenv("SU2_XIANZHI_VERIFY");
  return !value || std::string(value) != "0";
}
}

extern "C" unsigned long original_su2_solve(CSysSolve<double>*, CSysMatrix<double>&,
    const CSysVector<double>&, CSysVector<double>&, CGeometry*, const CConfig*) asm(SU2_REAL_SOLVE_SYMBOL);

template <>
unsigned long CSysSolve<double>::Solve(CSysMatrix<double>&, const CSysVector<double>&,
    CSysVector<double>&, CGeometry*, const CConfig*) asm(SU2_WRAP_SOLVE_SYMBOL);

template <>
unsigned long CSysSolve<double>::Solve(CSysMatrix<double>& matrix, const CSysVector<double>& rhs,
    CSysVector<double>& solution, CGeometry* geometry, const CConfig* config) {
  const char* selected = std::getenv("SU2_XIANZHI_BACKEND");
  const std::string backend = selected ? selected : "cpu";
  CallTimer timer{this, lin_sol_mode == LINEAR_SOLVER_MODE::STANDARD ? backend : "cpu", rhs.GetNVar()};
  if (backend == "cpu" || lin_sol_mode != LINEAR_SOLVER_MODE::STANDARD)
    return original_su2_solve(this, matrix, rhs, solution, geometry, config);
  if (backend != "gpu") SU2_MPI::Error("SU2_XIANZHI_BACKEND must be cpu or gpu.", CURRENT_FUNCTION);

  SU2_OMP_BARRIER
  SU2_OMP_MASTER {
    try {
      if (SU2_MPI::GetSize() != 1 || rhs.GetNBlk() != rhs.GetNBlkDomain())
        throw std::runtime_error("The new adapter requires one MPI process and no halo rows.");
      if (config->GetKind_Linear_Solver() != FGMRES || config->GetKind_Linear_Solver_Prec() != ILU ||
          config->GetLinear_Solver_ILU_n() != 0 || config->GetKind_Linear_Solver_Inner() != LINEAR_SOLVER_INNER::NONE)
        throw std::runtime_error("The adapter requires the original FGMRES + ILU(0) recipe without a nested solver.");
      if (config->GetDiscrete_Adjoint())
        throw std::runtime_error("Only ordinary primal flow is supported by the SDK adapter.");
      const auto n = rhs.GetNBlk();
      const auto bs = rhs.GetNVar();
      const auto& pattern = geometry->GetSparsePattern(ConnectivityType::FiniteVolume, 0);
      const auto* rows = pattern.outerPtr();
      const auto* cols = pattern.innerIdx();
      // The compressed pattern is built once per matrix; cache the two index
      // arrays instead of re-copying them on every linear solve.
      auto& cached = structure_cache[&matrix];
      if (cached.rows.size() != n + 1 || cached.rows.empty() || cached.rows.back() != rows[n]) {
        cached.rows.assign(rows, rows + n + 1);
        cached.columns.assign(cols, cols + rows[n]);
      }
      SU2Xianzhi45::Matrix exported;
      exported.block_size = bs;
      exported.rows = cached.rows;
      exported.columns = cached.columns;
      // CSysMatrix keeps every block in one contiguous array indexed by the
      // compressed-pattern position, i.e. the same layout as BCSR blocks.
      // GetBlock(0, cols[0]) resolves &matrix[row_ptr[0]]; the pointer is
      // borrowed and the library uploads it straight to the GPU.
      const double* host_base = matrix.GetBlock(0, cols[0]);
      if (!host_base) throw std::runtime_error("Cannot resolve the contiguous CSysMatrix block storage.");
      exported.host_values = host_base;
      exported.value_count = rows[n] * bs * bs;
#if defined(HAVE_CUDA)
      if (config->GetCUDA() && bs > 1) {
        // SU2 CUDA build: the matrix already lives on the GPU. Refresh the
        // value array with SU2's own HtDTransfer and hand the device pointers
        // straight to the solver instead of copying through host memory.
        matrix.HtDTransfer();
        exported.values = matrix.d_matrix;
        exported.device_rows = matrix.d_row_ptr;
        exported.device_columns = matrix.d_col_ind;
      } else {
        exported.values = host_base;
      }
#else
      exported.values = host_base;
#endif
      std::vector<double> b(rhs.GetLocSize()), x(solution.GetLocSize());
      for (unsigned long i = 0; i < rhs.GetLocSize(); ++i) { b[i] = rhs[i]; x[i] = solution[i]; }
      SU2Xianzhi45::Options options;
      options.tolerance = SU2_TYPE::GetValue(config->GetLinear_Solver_Error());
      options.max_iterations = config->GetLinear_Solver_Iter();
      options.preconditioner_level = config->GetXianzhiPreconditionerLevel();
      if (options.preconditioner_level > 3)
        throw std::runtime_error("XIANZHI_PRECONDITIONER_LEVEL must be between 0 and 3.");
      options.zero_initial_guess = xIsZero;
      // Non-ILU0 preconditioner levels do not accept a caller-supplied BSR
      // guess; SU2_XK_ZERO=1 forces a zero start so every level can be compared.
      if (const char* zero = std::getenv("SU2_XK_ZERO"); zero && std::string(zero) != "0")
        options.zero_initial_guess = true;
      options.relative_to_initial = tol_type == LinearToleranceType::RELATIVE;
      options.verify = verify_native();
      const char* directory = std::getenv("SU2_XIANZHI_OUTPUT_DIR");
      if (!directory) throw std::runtime_error("Set SU2_XIANZHI_OUTPUT_DIR for generated SDK descriptors and logs.");
      options.output_directory = directory;
      const auto result = SU2Xianzhi45::solve(this, exported, b, x, options);
      LinSysSol_tmp.Initialize(n, rhs.GetNBlkDomain(), bs, static_cast<const double*>(nullptr));
      LinSysRes_tmp.Initialize(n, rhs.GetNBlkDomain(), bs, static_cast<const double*>(nullptr));
      for (unsigned long i = 0; i < rhs.GetLocSize(); ++i) LinSysSol_tmp[i] = x[i];
      Iterations = result.iterations;
      Residual = result.normalizer > 0 ? result.residual_norm / result.normalizer : result.residual_norm;
    } catch (const std::exception& error) {
      SU2_MPI::Error(std::string("Original SDK45 adapter: ") + error.what(), CURRENT_FUNCTION);
    }
  }
  END_SU2_OMP_MASTER
  SU2_OMP_BARRIER
  solution = LinSysSol_tmp;
  SU2_OMP_BARRIER
  return Iterations;
}
