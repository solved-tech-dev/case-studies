// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (c) 2026 Solved Tech
#include <cuda_runtime_api.h>
#include <pangulu.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

static_assert(PANGULU_ABI_VERSION == 2);

namespace {
using Clock = std::chrono::steady_clock;

void require(bool ok, const char *message)
{
    if (!ok)
        throw std::runtime_error(message);
}

void check(int code, const char *phase)
{
    if (code)
        throw std::runtime_error(std::string(phase) + " status=" + std::to_string(code));
}

template<class Execute>
void native_phase(int call, const char *phase, Execute execute)
{
    check(cudaDeviceSynchronize(), "pre-phase synchronization");
    const auto start = Clock::now();
    const int status = execute();
    const int sync = cudaDeviceSynchronize();
    const double seconds = std::chrono::duration<double>(Clock::now() - start).count();
    std::cout << "NATIVE call=" << call << " phase=" << phase
              << " seconds=" << seconds << " status=" << (status ? status : sync)
              << std::endl;
    check(status, phase);
    check(sync, "post-phase synchronization");
}

template<class T>
std::vector<T> read_exact(const std::string &path, size_t count)
{
    std::ifstream file(path, std::ios::binary);
    std::vector<T> data(count);
    if (!file.read(reinterpret_cast<char *>(data.data()), count * sizeof(T)) ||
        file.peek() != std::ifstream::traits_type::eof() || file.bad())
        throw std::runtime_error("unexpected size or read error: " + path);
    return data;
}

struct Matrix {
    std::vector<int64_t> rows;
    std::vector<int32_t> columns;
    std::vector<double> values;
    std::vector<int64_t> csc_ptr;
    std::vector<int32_t> csc_rows;
    std::vector<double> csc_values;

    Matrix(const std::string &stem, int32_t n, int32_t nz)
        : rows(read_exact<int64_t>(stem + ".rp.i64", size_t(n) + 1)),
          columns(read_exact<int32_t>(stem + ".ci.i32", nz)),
          values(read_exact<double>(stem + ".a.f64", nz)),
          csc_ptr(size_t(n) + 1), csc_rows(nz), csc_values(nz)
    {
        require(rows.front() == 0 && rows.back() == nz, "invalid CSR endpoints");
        for (int32_t row = 0; row < n; ++row) {
            require(rows[row] >= 0 && rows[row] <= rows[row+1] && rows[row+1] <= nz,
                    "invalid CSR row pointer");
            int32_t previous = -1;
            for (auto k = rows[row]; k < rows[row+1]; ++k) {
                const int32_t col = columns[k];
                require(col > previous && col < n && std::isfinite(values[k]),
                        "invalid, duplicate or unsorted CSR entry");
                ++csc_ptr[col+1];
                previous = col;
            }
        }
        for (int32_t col = 0; col < n; ++col)
            csc_ptr[col+1] += csc_ptr[col];
        auto next = csc_ptr;
        for (int32_t row = 0; row < n; ++row) {
            for (auto k = rows[row]; k < rows[row+1]; ++k) {
                auto slot = next[columns[k]]++;
                csc_rows[slot] = row;
                csc_values[slot] = values[k];
            }
        }
    }
};

struct DeviceBuffer {
    void *data = nullptr;
    explicit DeviceBuffer(size_t bytes) { check(cudaMalloc(&data, bytes), "cudaMalloc"); }
    ~DeviceBuffer() { cudaFree(data); }
    DeviceBuffer(const DeviceBuffer &) = delete;
    DeviceBuffer &operator=(const DeviceBuffer &) = delete;
};

struct Solver {
    const size_t bytes;
    DeviceBuffer rhs, solution;
    pangulu_handle_t handle = nullptr;
    pangulu_matrix_t matrix = nullptr;
    pangulu_vector_t b = nullptr, x = nullptr;

    explicit Solver(int32_t n) : bytes(size_t(n) * sizeof(double)), rhs(bytes), solution(bytes) {}
    ~Solver()
    {
        if (x) pangulu_vector_destroy(x);
        if (b) pangulu_vector_destroy(b);
        if (matrix) pangulu_matrix_destroy(matrix);
        if (handle) pangulu_destroy(handle);
    }

    void factor(Matrix &a, int32_t n, int32_t nz)
    {
        check(pangulu_create(0, &handle), "PanguLU create");
        check(pangulu_config_set(handle, PANGULU_REORDERING_METHOD, L0), "PanguLU ordering");
        check(pangulu_config_set(handle, PANGULU_REUSE_LEVEL, L1), "PanguLU reuse");
        check(pangulu_matrix_create(&matrix, PANGULU_MATRIX_LAYOUT_CSC, n, n, nz, 1,
              a.csc_ptr.data(), a.csc_rows.data(), a.csc_values.data(),
              PANGULU_R_64I, PANGULU_R_32I, PANGULU_R_64F,
              PANGULU_LU, PANGULU_MEMORY_HOST), "PanguLU matrix");
        check(pangulu_multiple_vector_create(&b, n, 1, n, rhs.data,
              PANGULU_R_64F, PANGULU_MEMORY_DEVICE), "PanguLU RHS");
        check(pangulu_multiple_vector_create(&x, n, 1, n, solution.data,
              PANGULU_R_64F, PANGULU_MEMORY_DEVICE), "PanguLU solution");
        native_phase(1, "analysis", [&] {
            return pangulu_execute(handle, PANGULU_PHASE_ANALYSIS, matrix, nullptr, nullptr);
        });
        native_phase(1, "numeric", [&] {
            return pangulu_execute(handle, PANGULU_PHASE_NUMERIC_FACTORISATION, matrix, nullptr, nullptr);
        });
    }

    void solve(int call, const std::vector<double> &host_rhs, std::vector<double> &host_x)
    {
        check(cudaMemcpy(rhs.data, host_rhs.data(), bytes, cudaMemcpyHostToDevice), "RHS H2D");
        check(cudaMemset(solution.data, 0, bytes), "initialize solution");
        native_phase(call, "solve", [&] {
            return pangulu_execute(handle, PANGULU_PHASE_SOLVE, matrix, x, b);
        });
        check(cudaMemcpy(host_x.data(), solution.data, bytes, cudaMemcpyDeviceToHost), "solution D2H");
    }
};

void verify(int call, const Matrix &a, const std::vector<double> &rhs,
            const std::vector<double> &x, const std::vector<double> &reference)
{
    long double residual = 0, rhs_norm = 0, difference = 0, reference_norm = 0;
    for (size_t i = 0; i < x.size(); ++i) {
        require(std::isfinite(x[i]) && std::isfinite(rhs[i]) && std::isfinite(reference[i]),
                "nonfinite vector");
        long double r = -rhs[i];
        for (auto k = a.rows[i]; k < a.rows[i+1]; ++k)
            r += static_cast<long double>(a.values[k]) * x[a.columns[k]];
        const long double d = static_cast<long double>(x[i]) - reference[i];
        residual += r*r;
        rhs_norm += static_cast<long double>(rhs[i]) * rhs[i];
        difference += d*d;
        reference_norm += static_cast<long double>(reference[i]) * reference[i];
    }
    const auto rel_res = std::sqrt(residual) / (rhs_norm > 0 ? std::sqrt(rhs_norm) : 1.0L);
    const auto rel_difference = std::sqrt(difference) / std::max(std::sqrt(reference_norm), 1e-300L);
    std::cout << "VERIFY call=" << call << " residual=" << rel_res
              << " klu_relative_l2=" << rel_difference << std::endl;
    require(std::isfinite(rel_res) && rel_res <= 1e-8, "original CSR residual");
    require(std::isfinite(rel_difference) && rel_difference <= 1e-6, "KLU solution mismatch");
}

int32_t positive_index(const char *text)
{
    size_t consumed = 0;
    auto number = std::stoll(text, &consumed);
    require(consumed == std::strlen(text) && number > 0 &&
            number < std::numeric_limits<int32_t>::max(), "invalid LP64 dimension");
    return int32_t(number);
}
}

int main(int argc, char **argv) try {
    require(argc >= 4, "usage: replay ORDER CSR_ENTRIES CAPTURE_STEM...");
    require(std::string(PANGULU_VERSION_STRING) == "26.46", "expected PanguLU 26.46 headers");
    const auto n = positive_index(argv[1]), nz = positive_index(argv[2]);
    std::cout << std::setprecision(17);
    int device = 0, runtime = 0;
    cudaDeviceProp properties{};
    check(cudaGetDevice(&device), "CUDA device");
    check(cudaGetDeviceProperties(&properties, device), "CUDA device properties");
    check(cudaRuntimeGetVersion(&runtime), "CUDA runtime");
    std::cout << "GPU name=" << properties.name << " runtime=" << runtime << std::endl;
    Matrix matrix(argv[3], n, nz);
    Solver solver(n);
    solver.factor(matrix, n, nz);
    std::vector<double> x(n);
    for (int i = 3; i < argc; ++i) {
        const int call = i - 2;
        const std::string stem = argv[i];
        auto rhs = read_exact<double>(stem + ".b.f64", n);
        auto reference = read_exact<double>(stem + ".x.f64", n);
        solver.solve(call, rhs, x);
        verify(call, matrix, rhs, x, reference);
    }
    std::cout << "REPLAY_COMPLETE calls=" << argc - 3 << std::endl;
    return 0;
} catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 2;
}
