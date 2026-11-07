// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once
#include <cstdint>
#include <string>
#include <vector>

namespace SU2Xianzhi45 {
struct Matrix {
  std::uint32_t block_size = 1;
  std::vector<std::uint64_t> rows, columns;
  const double* values = nullptr;
  std::uint64_t value_count = 0;
  // SU2 CUDA build: device-resident arrays owned by CSysMatrix itself.
  const std::uint64_t* device_rows = nullptr;
  const std::uint64_t* device_columns = nullptr;
  // Host block array of the same matrix, used only for residual verification.
  const double* host_values = nullptr;
};
struct Options {
  double tolerance = 1e-6;
  unsigned long max_iterations = 100;
  unsigned short preconditioner_level = 2;
  bool zero_initial_guess = false;
  bool relative_to_initial = false;
  bool verify = true;
  int device = 0;
  std::string output_directory;
};
struct Result {
  unsigned long iterations = 0;
  int sdk_status = 0;
  double rhs_norm = 0, initial_norm = 0, normalizer = 0;
  double residual_norm = 0, target = 0, execute_seconds = 0;
  double upload_seconds = 0, stats_seconds = 0, verify_seconds = 0;
  double analysis_seconds = 0, numeric_prepare_seconds = 0, solve_seconds = 0;
  std::vector<double> matvec_product;
};

// Called by the OpenMP master. Only public SDK ABI 3 functions are used.
// The library does not alter any SU2_* / XIANZHI_* environment setting.
Result solve(const void* owner, const Matrix& matrix, const std::vector<double>& rhs,
             std::vector<double>& solution, const Options& options);
std::vector<double> product(const Matrix& matrix, const std::vector<double>& vector);
void clear();
}
