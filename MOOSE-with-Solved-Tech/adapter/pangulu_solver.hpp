// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Solved Tech
#pragma once
#include "sparse.hpp"
#include <chrono>
#include <cstdio>
#include <cuda_runtime_api.h>
#include <dlfcn.h>
#include <memory>
#include <pangulu.h>
#include <string>
static_assert(PANGULU_ABI_VERSION == 2, "PanguLU ABI 2 is required");
namespace solved {
using Clock = std::chrono::steady_clock;
inline void check(int status, const char *op) {
  if (status)
    throw std::runtime_error(std::string(op) + " returned status " + std::to_string(status));
}
struct Solver {
  Sparse matrix;
  pangulu_handle_t handle = nullptr;
  pangulu_matrix_t descriptor = nullptr;
  bool factored = false;
  unsigned analyses = 0, numerics = 0, solves = 0;
  ~Solver() { clear(); }
  void clear() noexcept {
    if (descriptor)
      pangulu_matrix_destroy(descriptor);
    if (handle)
      pangulu_destroy(handle);
    descriptor = nullptr;
    handle = nullptr;
    factored = false;
  }
  template <class F> void phase(const char *name, F fn) {
    check(cudaDeviceSynchronize(), "CUDA synchronize before native phase");
    const auto begin = Clock::now();
    const int status = fn();
    check(cudaDeviceSynchronize(), "CUDA synchronize after native phase");
    const double seconds = std::chrono::duration<double>(Clock::now() - begin).count();
    std::printf("MOOSE_PANGULU_NATIVE phase=%s seconds=%.12g status=%d\n", name, seconds, status);
    std::fflush(stdout);
    check(status, name);
  }
  void factor(int64_t n, const std::vector<int64_t> &p, const std::vector<int32_t> &c,
              const std::vector<double> &a) {
    require(std::string(PANGULU_VERSION_STRING) == "26.46", "Build against PanguLU 26.46");
    const bool graph_changed = !matrix.same_graph(n, p, c);
    if (graph_changed) {
      clear();
      matrix.reset(n, p, c);
    }
    const bool changed = matrix.refresh(a);
    if (!handle) {
      check(pangulu_create(0, &handle), "pangulu_create");
      check(pangulu_config_set(handle, PANGULU_REORDERING_METHOD, L1), "ordering L1");
      check(pangulu_config_set(handle, PANGULU_REUSE_LEVEL, L2), "reuse L2");
      Dl_info info{};
      dladdr(reinterpret_cast<void *>(pangulu_execute), &info);
      cudaDeviceProp props{};
      int device = 0;
      check(cudaGetDevice(&device), "cudaGetDevice");
      check(cudaGetDeviceProperties(&props, device), "cudaGetDeviceProperties");
      std::printf(
          "MOOSE_PANGULU_PROVIDER version=%s abi=%d ordering=L1 reuse=L2 gpu=%s library=%s\n",
          PANGULU_VERSION_STRING, PANGULU_ABI_VERSION, props.name,
          info.dli_fname ? info.dli_fname : "unknown");
      check(pangulu_matrix_create(&descriptor, PANGULU_MATRIX_LAYOUT_CSC, n, n, int64_t(c.size()),
                                  1, matrix.cp.data(), matrix.ri.data(), matrix.csc.data(),
                                  PANGULU_R_64I, PANGULU_R_32I, PANGULU_R_64F, PANGULU_LU,
                                  PANGULU_MEMORY_HOST),
            "matrix create");
      phase("analysis_order_symbolic", [&] {
        return pangulu_execute(handle, PANGULU_PHASE_ANALYSIS, descriptor, nullptr, nullptr);
      });
      ++analyses;
    }
    if (!factored || changed) {
      const auto stage =
          factored ? PANGULU_PHASE_REFACTORISATION : PANGULU_PHASE_NUMERIC_FACTORISATION;
      phase(factored ? "refactor_numeric" : "numeric",
            [&] { return pangulu_execute(handle, stage, descriptor, nullptr, nullptr); });
      factored = true;
      ++numerics;
    } else
      std::printf("MOOSE_PANGULU_REUSE numeric=retained coefficients=exactly_equal\n");
  }
  // A block solve is used only when PETSc actually supplies multiple independent RHS.
  void solve(const double *b, double *x, int64_t nrhs = 1, double tolerance = 1e-8) {
    require(factored && nrhs > 0, "Solve before factorization or empty RHS block");
    require(nrhs <= INT64_MAX / matrix.n && size_t(matrix.n * nrhs) <= SIZE_MAX / sizeof(double),
            "RHS block too large");
    const size_t bytes = size_t(matrix.n * nrhs) * sizeof(double);
    struct Buffers {
      void *b = nullptr, *x = nullptr;
      pangulu_vector_t bv = nullptr, xv = nullptr;
      ~Buffers() {
        if (xv)
          pangulu_vector_destroy(xv);
        if (bv)
          pangulu_vector_destroy(bv);
        if (x)
          cudaFree(x);
        if (b)
          cudaFree(b);
      }
    } d;
    check(cudaMalloc(&d.b, bytes), "allocate RHS");
    check(cudaMalloc(&d.x, bytes), "allocate solution");
    check(cudaMemcpy(d.b, b, bytes, cudaMemcpyHostToDevice), "upload RHS");
    check(cudaMemset(d.x, 0, bytes), "clear solution");
    check(pangulu_multiple_vector_create(&d.bv, matrix.n, nrhs, matrix.n, d.b, PANGULU_R_64F,
                                         PANGULU_MEMORY_DEVICE),
          "RHS descriptor");
    check(pangulu_multiple_vector_create(&d.xv, matrix.n, nrhs, matrix.n, d.x, PANGULU_R_64F,
                                         PANGULU_MEMORY_DEVICE),
          "solution descriptor");
    phase("solve",
          [&] { return pangulu_execute(handle, PANGULU_PHASE_SOLVE, descriptor, d.xv, d.bv); });
    ++solves;
    check(cudaMemcpy(x, d.x, bytes, cudaMemcpyDeviceToHost), "download solution");
    for (int64_t k = 0; k < nrhs; ++k) {
      const auto begin = Clock::now();
      const double r = matrix.residual(b + k * matrix.n, x + k * matrix.n);
      std::printf("MOOSE_PANGULU_CHECK solve=%u rhs=%lld n=%lld stored_nnz=%zu relative_l2=%.17g "
                  "evidence_seconds=%.12g analyses=%u numerics=%u\n",
                  solves, static_cast<long long>(k), static_cast<long long>(matrix.n),
                  matrix.ci.size(), r, std::chrono::duration<double>(Clock::now() - begin).count(),
                  analyses, numerics);
      require(r <= tolerance, "Original full CSR residual exceeds tolerance");
    }
  }
};
} // namespace solved
