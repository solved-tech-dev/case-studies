// SPDX-License-Identifier: AGPL-3.0-only
// Copyright 2026 Solved Tech
#include "pangulu_adapter.h"
#include <pangulu.h>
#include <cuda_runtime_api.h>
#include <algorithm>
#include <chrono>
#include <complex>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

using Complex = std::complex<double>;
static_assert(sizeof(Complex) == 16 && sizeof(int) == 4);

namespace {
struct System {
    pangulu_handle_t handle = nullptr;
    pangulu_matrix_t matrix = nullptr;
    std::vector<int64_t> pointers;
    std::vector<int32_t> rows;
    std::vector<Complex> coefficients;
    int n = 0;
};

System system_state;

void check(int status, const char *operation) {
    if (status) throw std::runtime_error(std::string(operation) + " status=" + std::to_string(status));
}

int report(const char *entry, const std::exception &error) {
    std::fprintf(stderr, "PANGULU_ERROR entry=%s %s\n", entry, error.what());
    return 70;
}

void phase(pangulu_phase_t operation, const char *name, pangulu_vector_t x = nullptr,
           pangulu_vector_t b = nullptr) {
    check(cudaDeviceSynchronize(), "synchronize before phase");
    const auto start = std::chrono::steady_clock::now();
    const int status = pangulu_execute(system_state.handle, operation, system_state.matrix, x, b);
    const int sync = cudaDeviceSynchronize();
    const auto elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start);
    std::printf("PANGULU_PHASE phase=%s seconds=%.9f status=%d\n", name, elapsed.count(), status ? status : sync);
    std::fflush(stdout);
    check(status, name);
    check(sync, "synchronize after phase");
}

void factor(int n, int64_t nnz, const int *rows, const int *cols, const double *real, const double *imag) {
    if (n <= 0 || nnz <= 0) throw std::runtime_error("invalid matrix size");
    check(case_cleanup(), "release previous factors");
    auto &s = system_state;
    s.n = n;
    s.pointers.assign(static_cast<size_t>(n) + 1, 0);
    s.rows.resize(nnz);
    s.coefficients.resize(nnz);
    for (int64_t k = 0; k < nnz; ++k) {
        if (rows[k] < 0 || rows[k] >= n || cols[k] < 0 || cols[k] >= n)
            throw std::runtime_error("COO coordinate outside matrix");
        ++s.pointers[cols[k] + 1];
    }
    for (int col = 0; col < n; ++col) s.pointers[col + 1] += s.pointers[col];
    auto next = s.pointers;
    for (int64_t k = 0; k < nnz; ++k) {
        const auto position = next[cols[k]]++;
        s.rows[position] = rows[k];
        s.coefficients[position] = Complex(real[k], imag[k]);
    }
    std::vector<std::pair<int32_t, Complex>> column;
    for (int col = 0; col < n; ++col) {
        const auto start = s.pointers[col];
        column.clear();
        for (auto k = start; k < s.pointers[col + 1]; ++k)
            column.emplace_back(s.rows[k], s.coefficients[k]);
        std::sort(column.begin(), column.end(), [](const auto &a, const auto &b) { return a.first < b.first; });
        for (size_t k = 0; k < column.size(); ++k) {
            if (k && column[k].first == column[k - 1].first)
                throw std::runtime_error("duplicate COO entry");
            s.rows[start + k] = column[k].first;
            s.coefficients[start + k] = column[k].second;
        }
    }
    std::printf("PANGULU_CONFIG version=%s ordering=L0 reuse=L1 precision=complex128 nrhs=1 n=%d nnz=%lld\n",
                PANGULU_VERSION_STRING, n, static_cast<long long>(nnz));
    check(pangulu_create(0, &s.handle), "create (verify the SDK and GPU license)");
    check(pangulu_config_set(s.handle, PANGULU_REORDERING_METHOD, L0), "configure ordering");
    check(pangulu_config_set(s.handle, PANGULU_REUSE_LEVEL, L1), "configure reuse");
    check(pangulu_matrix_create(&s.matrix, PANGULU_MATRIX_LAYOUT_CSC, n, n, nnz, 1,
                               s.pointers.data(), s.rows.data(), s.coefficients.data(), PANGULU_R_64I,
                               PANGULU_R_32I, PANGULU_C_64F, PANGULU_LU, PANGULU_MEMORY_HOST), "create CSC matrix");
    phase(PANGULU_PHASE_ANALYSIS, "analysis");
    phase(PANGULU_PHASE_NUMERIC_FACTORISATION, "numeric");
}

struct Vectors {
    Complex *db = nullptr, *dx = nullptr;
    pangulu_vector_t b = nullptr, x = nullptr;
    ~Vectors() {
        if (x) pangulu_vector_destroy(x);
        if (b) pangulu_vector_destroy(b);
        if (dx) cudaFree(dx);
        if (db) cudaFree(db);
    }
};

std::vector<Complex> solve(int n, const Complex *rhs) {
    const auto &s = system_state;
    if (n != s.n || !s.handle || !s.matrix) throw std::runtime_error("solve requires matching factors");
    const size_t bytes = static_cast<size_t>(n) * sizeof(Complex);
    Vectors device;
    check(cudaMalloc(reinterpret_cast<void **>(&device.db), bytes), "allocate RHS");
    check(cudaMalloc(reinterpret_cast<void **>(&device.dx), bytes), "allocate solution");
    check(cudaMemset(device.dx, 0, bytes), "initialize solution");
    check(pangulu_multiple_vector_create(&device.b, n, 1, n, device.db, PANGULU_C_64F,
                                         PANGULU_MEMORY_DEVICE), "create RHS vector");
    check(pangulu_multiple_vector_create(&device.x, n, 1, n, device.dx, PANGULU_C_64F,
                                         PANGULU_MEMORY_DEVICE), "create solution vector");
    check(cudaMemcpy(device.db, rhs, bytes, cudaMemcpyHostToDevice), "upload RHS");
    phase(PANGULU_PHASE_SOLVE, "solve", device.x, device.b);
    std::vector<Complex> x(n);
    check(cudaMemcpy(x.data(), device.dx, bytes, cudaMemcpyDeviceToHost), "download solution");
    return x;
}

}  // namespace

extern "C" int case_cleanup() {
    auto &s = system_state;
    const int matrix_status = s.matrix ? pangulu_matrix_destroy(s.matrix) : 0;
    const int handle_status = s.handle ? pangulu_destroy(s.handle) : 0;
    s = System{};
    return matrix_status ? matrix_status : handle_status;
}

extern "C" int case_factor(int n, int64_t nnz, const int *rows, const int *cols,
                           const double *ar, const double *ai) {
    try {
        factor(n, nnz, rows, cols, ar, ai);
        return 0;
    } catch (const std::exception &error) { return report("factor", error); }
}

extern "C" int case_solve(int n, const double *br, const double *bi, double *xr, double *xi) {
    try {
        std::vector<Complex> rhs(n);
        for (int i = 0; i < n; ++i) rhs[i] = Complex(br[i], bi[i]);
        const auto x = solve(n, rhs.data());
        for (int i = 0; i < n; ++i) {
            xr[i] = x[i].real();
            xi[i] = x[i].imag();
        }
        return 0;
    } catch (const std::exception &error) { return report("solve", error); }
}
