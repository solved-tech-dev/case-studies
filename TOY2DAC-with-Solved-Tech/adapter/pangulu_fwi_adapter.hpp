// SPDX-License-Identifier: AGPL-3.0-only
// Copyright 2026 Solved Tech
#pragma once

#include <complex>
#include <cstdint>

using SolverComplex = std::complex<double>;
constexpr int kRhs = 1;
constexpr double kTolerance = 1e-8;

extern "C" {
int case_factor(int n, int64_t nnz, const int *rows, const int *cols,
                const double *real, const double *imag);
int case_solve_many(const SolverComplex *rhs, SolverComplex *solution);
int case_cleanup();
}
