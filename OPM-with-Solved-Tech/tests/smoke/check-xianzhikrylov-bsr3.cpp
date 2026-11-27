#include <XianzhiKrylovAPI.h>
#include <cuda_runtime_api.h>

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void checkCuda(cudaError_t err, const char* what)
{
    if (err != cudaSuccess) {
        std::fprintf(stderr, "CUDA error (%s): %s\n", what, cudaGetErrorString(err));
        std::exit(1);
    }
}

void checkXk(xianzhikrylov_status_t status, xianzhikrylov_status_t expected, const char* what,
             xianzhikrylov_handle_t handle = nullptr)
{
    if (status == expected) {
        return;
    }
    if (handle) {
        xianzhikrylov_print_info(handle);
    }
    std::fprintf(stderr, "XianzhiKrylov error (%s): expected %d got %d\n",
                 what, static_cast<int>(expected), static_cast<int>(status));
    std::exit(1);
}

void assertNoLinglongLoaded()
{
    std::ifstream maps("/proc/self/maps");
    std::string line;
    while (std::getline(maps, line)) {
        if (line.find("liblinglongamg") != std::string::npos) {
            std::fprintf(stderr, "Isolation failure: unexpected liblinglongamg mapping: %s\n",
                         line.c_str());
            std::exit(1);
        }
    }
}

double computeTrueRelativeL2(int64_t blockRows,
                             const std::vector<uint64_t>& rowOffsets,
                             const std::vector<uint64_t>& cols,
                             const std::vector<double>& vals,
                             const std::vector<double>& x,
                             const std::vector<double>& b)
{
    constexpr int B = 3;
    double errSq = 0.0;
    double rhsSq = 0.0;
    for (int64_t r = 0; r < blockRows; ++r) {
        double ax[B] = {0.0, 0.0, 0.0};
        for (uint64_t p = rowOffsets[r]; p < rowOffsets[r + 1]; ++p) {
            const uint64_t c = cols[p];
            for (int i = 0; i < B; ++i) {
                for (int j = 0; j < B; ++j) {
                    ax[i] += vals[p * B * B + i * B + j] * x[c * B + j];
                }
            }
        }
        for (int i = 0; i < B; ++i) {
            const double diff = b[r * B + i] - ax[i];
            errSq += diff * diff;
            rhsSq += b[r * B + i] * b[r * B + i];
        }
    }
    return std::sqrt(errSq) / std::sqrt(rhsSq);
}

} // namespace

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::fprintf(stderr, "Usage: %s <xk-bsr3-config.json>\n", argv[0]);
        return 2;
    }
    const char* configPath = argv[1];

    constexpr int64_t kRows = 2;
    constexpr int64_t kBlocks = 4;
    constexpr int B = 3;
    constexpr int64_t kScalars = kRows * B;

    const std::vector<uint64_t> rowOffsets = {0, 2, 4};
    const std::vector<uint64_t> cols = {0, 1, 0, 1};
    const std::vector<uint64_t> diag = {0, 3};

    // Pass 1: nonsymmetric diagonally dominant 2-block BSR3 system.
    std::vector<double> vals = {
        // block (0, 0)
        6.0, -0.5, 0.2,
        0.1,  5.5, -0.3,
        -0.2, 0.4, 7.0,
        // block (0, 1)
        -0.8, 0.1, 0.0,
        0.2, -0.6, 0.1,
        0.0,  0.2, -0.5,
        // block (1, 0)
        -0.4, 0.0, 0.2,
        0.1, -0.7, 0.1,
        -0.1, 0.2, -0.6,
        // block (1, 1)
        5.8, 0.3, -0.1,
        -0.2, 6.2, 0.4,
        0.1, -0.3, 6.5
    };
    const std::vector<double> rhs = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    std::vector<double> sol(kScalars, 0.0);

    uint64_t *dRows = nullptr, *dCols = nullptr, *dDiag = nullptr;
    double *dVals = nullptr, *dRhs = nullptr, *dSol = nullptr;

    checkCuda(cudaMalloc(reinterpret_cast<void**>(&dRows), rowOffsets.size() * sizeof(uint64_t)), "dRows");
    checkCuda(cudaMalloc(reinterpret_cast<void**>(&dCols), cols.size() * sizeof(uint64_t)), "dCols");
    checkCuda(cudaMalloc(reinterpret_cast<void**>(&dDiag), diag.size() * sizeof(uint64_t)), "dDiag");
    checkCuda(cudaMalloc(reinterpret_cast<void**>(&dVals), vals.size() * sizeof(double)), "dVals");
    checkCuda(cudaMalloc(reinterpret_cast<void**>(&dRhs), rhs.size() * sizeof(double)), "dRhs");
    checkCuda(cudaMalloc(reinterpret_cast<void**>(&dSol), sol.size() * sizeof(double)), "dSol");

    checkCuda(cudaMemcpy(dRows, rowOffsets.data(), rowOffsets.size() * sizeof(uint64_t), cudaMemcpyHostToDevice), "copy rows");
    checkCuda(cudaMemcpy(dCols, cols.data(), cols.size() * sizeof(uint64_t), cudaMemcpyHostToDevice), "copy cols");
    checkCuda(cudaMemcpy(dDiag, diag.data(), diag.size() * sizeof(uint64_t), cudaMemcpyHostToDevice), "copy diag");
    checkCuda(cudaMemcpy(dVals, vals.data(), vals.size() * sizeof(double), cudaMemcpyHostToDevice), "copy vals 1");
    checkCuda(cudaMemcpy(dRhs, rhs.data(), rhs.size() * sizeof(double), cudaMemcpyHostToDevice), "copy rhs");

    xianzhikrylov_handle_t handle = nullptr;
    xianzhikrylov_matrix_t matrix = nullptr;
    xianzhikrylov_vector_t xv = nullptr;
    xianzhikrylov_vector_t bv = nullptr;

    checkXk(xianzhikrylov_create(0, &handle), XIANZHIKRYLOV_SUCCESS, "create");
    checkXk(xianzhikrylov_config_set(handle, XIANZHIKRYLOV_CONFIG_FILE, configPath),
            XIANZHIKRYLOV_SUCCESS, "config_set", handle);

    xianzhikrylov_matrix_data_t data{};
    data.abi_version = XIANZHIKRYLOV_ABI_VERSION;
    data.struct_size = sizeof(data);
    data.format = XIANZHIKRYLOV_BSR;
    data.value_type = XIANZHIKRYLOV_FP64;
    data.index_type = XIANZHIKRYLOV_INDEX_I64;
    data.block_size = B;
    data.block_layout = XIANZHIKRYLOV_BLOCK_ROW_MAJOR;
    data.memory_space = XIANZHIKRYLOV_DEVICE;
    data.device_id = 0;
    data.rows = kRows;
    data.columns = kRows;
    data.entries = kBlocks;
    data.row_offsets = {dRows, rowOffsets.size() * sizeof(uint64_t)};
    data.column_indices = {dCols, cols.size() * sizeof(uint64_t)};
    data.diagonal_positions = {dDiag, diag.size() * sizeof(uint64_t)};
    data.values = {dVals, vals.size() * sizeof(double)};

    checkXk(xianzhikrylov_matrix_create(&matrix, &data), XIANZHIKRYLOV_SUCCESS, "matrix_create", handle);
    checkXk(xianzhikrylov_vector_create(&xv, kScalars, dSol, XIANZHIKRYLOV_FP64, XIANZHIKRYLOV_DEVICE),
            XIANZHIKRYLOV_SUCCESS, "vector_create x", handle);
    checkXk(xianzhikrylov_vector_create(&bv, kScalars, dRhs, XIANZHIKRYLOV_FP64, XIANZHIKRYLOV_DEVICE),
            XIANZHIKRYLOV_SUCCESS, "vector_create b", handle);

    checkXk(xianzhikrylov_execute(handle, XIANZHIKRYLOV_PHASE_ANALYSIS, matrix, nullptr, nullptr),
            XIANZHIKRYLOV_SUCCESS, "analysis", handle);
    checkXk(xianzhikrylov_execute(handle, XIANZHIKRYLOV_PHASE_NUMERIC_PREPARE, matrix, nullptr, nullptr),
            XIANZHIKRYLOV_SUCCESS, "numeric 1", handle);
    checkXk(xianzhikrylov_execute(handle, XIANZHIKRYLOV_PHASE_SOLVE, matrix, xv, bv),
            XIANZHIKRYLOV_CONVERGED, "solve 1", handle);
    checkXk(xianzhikrylov_print_info(handle), XIANZHIKRYLOV_SUCCESS, "print_info 1", handle);

    checkCuda(cudaMemcpy(sol.data(), dSol, sol.size() * sizeof(double), cudaMemcpyDeviceToHost), "download x 1");
    const double rel1 = computeTrueRelativeL2(kRows, rowOffsets, cols, vals, sol, rhs);
    if (!std::isfinite(rel1) || rel1 > 0.005) {
        std::fprintf(stderr, "Pass 1 residual check failed: rel1=%.12g\n", rel1);
        return 1;
    }

    // Pass 2: update numerical values while reusing structure analysis.
    for (double& v : vals) {
        v *= 1.5;
    }
    checkCuda(cudaMemcpy(dVals, vals.data(), vals.size() * sizeof(double), cudaMemcpyHostToDevice), "copy vals 2");
    checkXk(xianzhikrylov_execute(handle, XIANZHIKRYLOV_PHASE_NUMERIC_PREPARE, matrix, nullptr, nullptr),
            XIANZHIKRYLOV_SUCCESS, "numeric 2", handle);
    checkXk(xianzhikrylov_execute(handle, XIANZHIKRYLOV_PHASE_SOLVE, matrix, xv, bv),
            XIANZHIKRYLOV_CONVERGED, "solve 2", handle);
    checkCuda(cudaMemcpy(sol.data(), dSol, sol.size() * sizeof(double), cudaMemcpyDeviceToHost), "download x 2");
    const double rel2 = computeTrueRelativeL2(kRows, rowOffsets, cols, vals, sol, rhs);
    if (!std::isfinite(rel2) || rel2 > 0.005) {
        std::fprintf(stderr, "Pass 2 residual check failed: rel2=%.12g\n", rel2);
        return 1;
    }

    assertNoLinglongLoaded();

    xianzhikrylov_destroy(handle);
    xianzhikrylov_vector_destroy(bv);
    xianzhikrylov_vector_destroy(xv);
    xianzhikrylov_matrix_destroy(matrix);
    cudaFree(dSol);
    cudaFree(dRhs);
    cudaFree(dVals);
    cudaFree(dDiag);
    cudaFree(dCols);
    cudaFree(dRows);

    std::printf("XIANZHIKRYLOV_BSR3_SMOKE_PASS rel1=%.12g rel2=%.12g\n", rel1, rel2);
    return 0;
}
