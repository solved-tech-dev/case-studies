// SPDX-License-Identifier: AGPL-3.0-only
#include <cuda_runtime_api.h>
#include <pangulu.h>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
static_assert(PANGULU_ABI_VERSION == 2, "PanguLU ABI 2 is required");
#ifdef FEALPY_COMPLEX
using Scalar = std::complex<double>;
constexpr auto value_type = PANGULU_C_64F;
#else
using Scalar = double;
constexpr auto value_type = PANGULU_R_64F;
#endif

void check(int status, const char* operation) {
    if (status != 0)
        throw std::runtime_error(std::string(operation) + " returned " + std::to_string(status));
}

struct CscMatrix {
    std::vector<int64_t> pointers;
    std::vector<int32_t> rows;
    std::vector<Scalar> values;

    CscMatrix(int n, const int64_t* rp, const int32_t* ci, const Scalar* a)
        : pointers(static_cast<size_t>(n) + 1) {
        if (rp[0] != 0) throw std::invalid_argument("CSR must start at zero");
        for (int i = 0; i < n; ++i) {
            if (rp[i] > rp[i+1]) throw std::invalid_argument("invalid CSR pointers");
            for (auto k = rp[i]; k < rp[i+1]; ++k) {
                if (ci[k] < 0 || ci[k] >= n) throw std::invalid_argument("invalid column");
                ++pointers[static_cast<size_t>(ci[k]) + 1];
            }
        }
        for (size_t i = 1; i < pointers.size(); ++i) pointers[i] += pointers[i-1];
        rows.resize(pointers.back());
        values.resize(pointers.back());
        auto next = pointers;
        for (int i = 0; i < n; ++i) {
            for (auto k = rp[i]; k < rp[i+1]; ++k) {
                const auto dest = next[ci[k]]++;
                rows[dest] = i;
                values[dest] = a[k];
            }
        }
    }
};

struct SolveResources {
    pangulu_handle_t handle = nullptr;
    pangulu_matrix_t matrix = nullptr;
    pangulu_vector_t b = nullptr, x = nullptr;
    void* device_b = nullptr;
    void* device_x = nullptr;
    ~SolveResources() {
        if (x) pangulu_vector_destroy(x);
        if (b) pangulu_vector_destroy(b);
        if (matrix) pangulu_matrix_destroy(matrix);
        if (handle) pangulu_destroy(handle);
        if (device_x) cudaFree(device_x);
        if (device_b) cudaFree(device_b);
    }
};

void solve_system(int n, const int64_t* rp, const int32_t* ci, const Scalar* a,
                  const Scalar* b, Scalar* x) {
    CscMatrix input(n, rp, ci, a);
    SolveResources ctx;
    const size_t bytes = static_cast<size_t>(n) * sizeof(Scalar);
    check(pangulu_create(0, &ctx.handle), "pangulu_create");
    check(pangulu_config_set(ctx.handle, PANGULU_REORDERING_METHOD, L0), "ordering");
    check(pangulu_config_set(ctx.handle, PANGULU_REUSE_LEVEL, L2), "reuse level");
    check(pangulu_matrix_create(&ctx.matrix, PANGULU_MATRIX_LAYOUT_CSC, n, n, rp[n], 1,
                                input.pointers.data(), input.rows.data(), input.values.data(),
                                PANGULU_R_64I, PANGULU_R_32I, value_type, PANGULU_LU,
                                PANGULU_MEMORY_HOST), "matrix creation");
    check(cudaMalloc(&ctx.device_b, bytes), "allocate RHS");
    check(cudaMalloc(&ctx.device_x, bytes), "allocate solution");
    check(pangulu_multiple_vector_create(&ctx.b, n, 1, n, ctx.device_b, value_type,
                                         PANGULU_MEMORY_DEVICE), "RHS creation");
    check(pangulu_multiple_vector_create(&ctx.x, n, 1, n, ctx.device_x, value_type,
                                         PANGULU_MEMORY_DEVICE), "solution creation");
    check(pangulu_execute(ctx.handle, PANGULU_PHASE_ANALYSIS, ctx.matrix, nullptr, nullptr), "analysis");
    check(cudaDeviceSynchronize(), "analysis synchronization");
    check(pangulu_execute(ctx.handle, PANGULU_PHASE_NUMERIC_FACTORISATION, ctx.matrix, nullptr, nullptr), "factorization");
    check(cudaDeviceSynchronize(), "factorization synchronization");
    check(cudaMemcpy(ctx.device_b, b, bytes, cudaMemcpyHostToDevice), "upload RHS");
    check(cudaMemset(ctx.device_x, 0, bytes), "initialize solution");
    check(pangulu_execute(ctx.handle, PANGULU_PHASE_SOLVE, ctx.matrix, ctx.x, ctx.b), "solve");
    check(cudaDeviceSynchronize(), "solve synchronization");
    check(cudaMemcpy(x, ctx.device_x, bytes, cudaMemcpyDeviceToHost), "download solution");
}
}  // namespace
