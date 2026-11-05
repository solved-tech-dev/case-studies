#include "../config/config.hpp"

#ifdef OCP_USE_LINGLONGAMG

#include "LinglongAMGSolver.hpp"

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>

#include <cuda_runtime.h>

#include "GpuSolverCommon.hpp"
#include "Spe10NativeScaling.hpp"
#include "UtilError.hpp"

namespace {

using Clock = std::chrono::steady_clock;

double SecondsSince(const Clock::time_point& start)
{
    return std::chrono::duration<double>(Clock::now() - start).count();
}

void CheckCuda(const cudaError_t status, const char* phase)
{
    if (status != cudaSuccess) {
        OCP_ABORT(string("LinglongAMG CUDA ") + phase + " failed: " +
                  cudaGetErrorString(status));
    }
}

void CheckLinglong(const linglongamg_status_t status, const char* phase)
{
    if (status != LINGLONGAMG_STATUS_SUCCESS) {
        OCP_ABORT(string("LinglongAMG ") + phase + " failed: " +
                  linglongamg_status_string(status));
    }
}

bool IsAbsolutePath(const string& path)
{
    return !path.empty() && path.front() == '/';
}

} // namespace

LinglongAMGSolver::LinglongAMGSolver(const string& dir, const string& file,
                                     const OCPMatrix&)
{
    const char* envConfig = std::getenv("OCP_LINGLONGAMG_CONFIG");
    if (envConfig != nullptr && *envConfig != '\0') {
        configPath = envConfig;
    } else if (file != "linglongamg") {
        configPath = IsAbsolutePath(file) ? file : dir + file;
    }
    if (configPath.empty() || !std::ifstream(configPath).good()) {
        OCP_ABORT("LinglongAMG config is missing; set OCP_LINGLONGAMG_CONFIG");
    }
}

void LinglongAMGSolver::ReleaseDeviceStorage()
{
    if (deviceRowPtr != nullptr) CheckCuda(cudaFree(deviceRowPtr), "free-rowptr");
    if (deviceColIdx != nullptr) CheckCuda(cudaFree(deviceColIdx), "free-colidx");
    if (deviceValues != nullptr) CheckCuda(cudaFree(deviceValues), "free-values");
    if (deviceRhs != nullptr) CheckCuda(cudaFree(deviceRhs), "free-rhs");
    if (deviceSolution != nullptr) CheckCuda(cudaFree(deviceSolution), "free-solution");
    deviceRowPtr = nullptr;
    deviceColIdx = nullptr;
    deviceValues = nullptr;
    deviceRhs = nullptr;
    deviceSolution = nullptr;
}

void LinglongAMGSolver::AssembleMat(OCPMatrix& mat, const Domain* domain)
{
    if (domain == nullptr || domain->cs_numproc != 1) {
        OCP_ABORT("LinglongAMG integration currently requires one OpenCAEPoroX MPI rank");
    }

    const auto start = Clock::now();
    ReleaseDeviceStorage();
    const int factor = ocp_spe10_scaling::RefineFactor();
    const char* singleEnv = std::getenv("OCP_SINGLE_SOLVE_BENCHMARK");
    scalingBenchmark = factor > 1 ||
        (singleEnv != nullptr && *singleEnv != '\0' && *singleEnv != '0');
    if (scalingBenchmark) {
        auto csr = ocp_spe10_scaling::Build(mat, *domain, factor);
        rowPtr.swap(csr.rowPtr);
        colIdx.swap(csr.colIdx);
        values.swap(csr.values);
        rhs.swap(csr.rhs);
        solution.swap(csr.solution);
    } else {
        ocp_gpu_solver::BuildSortedScalarCsr(mat, rowPtr, colIdx, values);
        rhs.assign(mat.b.begin(), mat.b.end());
        solution.assign(mat.u.begin(), mat.u.end());
    }
    if (values.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        OCP_ABORT("LinglongAMG integration requires NNZ <= INT_MAX");
    }

    matrix = &mat;
    const std::size_t n = rowPtr.size() - 1;
    CheckCuda(cudaMalloc(reinterpret_cast<void**>(&deviceRowPtr),
                         rowPtr.size() * sizeof(int)), "malloc-rowptr");
    CheckCuda(cudaMalloc(reinterpret_cast<void**>(&deviceColIdx),
                         colIdx.size() * sizeof(int)), "malloc-colidx");
    CheckCuda(cudaMalloc(reinterpret_cast<void**>(&deviceValues),
                         values.size() * sizeof(double)), "malloc-values");
    CheckCuda(cudaMalloc(reinterpret_cast<void**>(&deviceRhs),
                         n * sizeof(double)), "malloc-rhs");
    CheckCuda(cudaMalloc(reinterpret_cast<void**>(&deviceSolution),
                         n * sizeof(double)), "malloc-solution");
    CheckCuda(cudaMemcpy(deviceRowPtr, rowPtr.data(),
                         rowPtr.size() * sizeof(int), cudaMemcpyHostToDevice),
              "copy-rowptr");
    CheckCuda(cudaMemcpy(deviceColIdx, colIdx.data(),
                         colIdx.size() * sizeof(int), cudaMemcpyHostToDevice),
              "copy-colidx");
    CheckCuda(cudaMemcpy(deviceValues, values.data(),
                         values.size() * sizeof(double), cudaMemcpyHostToDevice),
              "copy-values");
    CheckCuda(cudaMemcpy(deviceRhs, rhs.data(),
                         n * sizeof(double), cudaMemcpyHostToDevice), "copy-rhs");
    CheckCuda(cudaMemcpy(deviceSolution, solution.data(),
                         n * sizeof(double), cudaMemcpyHostToDevice),
              "copy-solution");

    cout << fixed << setprecision(6)
         << "[OCP_LINGLONGAMG] phase=preprocess rows=" << n
         << " nnz=" << values.size()
         << " seconds=" << SecondsSince(start) << endl;
}

OCP_INT LinglongAMGSolver::Solve()
{
    if (matrix == nullptr || deviceRowPtr == nullptr) {
        OCP_ABORT("LinglongAMG Solve called before matrix assembly");
    }

    const auto totalStart = Clock::now();
    linglongamg_handle_t solver = nullptr;
    CheckLinglong(linglongamg_create_from_config(&solver, configPath.c_str()),
                  "create");

    // Move CUDA context creation out of the measured AMG setup region.
    CheckCuda(cudaFree(nullptr), "warmup");

    const int nrows = static_cast<int>(rowPtr.size() - 1);
#ifdef OCP_LINGLONGAMG_OPAQUE_CSR_API
    linglongamg_csr_t* A = nullptr;
    CheckLinglong(linglongamg_read_csr(&A, nrows, nrows,
                                      static_cast<int>(values.size()),
                                      deviceRowPtr, deviceColIdx, deviceValues),
                  "read-csr");
#else
    linglongamg_csr_t AStorage{};
    AStorage.nrows = nrows;
    AStorage.ncols = nrows;
    AStorage.nnz = static_cast<int>(values.size());
    AStorage.row_ptr = deviceRowPtr;
    AStorage.col_idx = deviceColIdx;
    AStorage.val = deviceValues;
    linglongamg_csr_t* A = &AStorage;
#endif
    linglongamg_vec_t b{nrows, deviceRhs};
    linglongamg_vec_t x{nrows, deviceSolution};

    auto start = Clock::now();
    CheckLinglong(linglongamg_setup(solver, A), "setup");
    CheckCuda(cudaDeviceSynchronize(), "setup-sync");
    const double setupSeconds = SecondsSince(start);

    start = Clock::now();
    CheckLinglong(linglongamg_solve(solver, &b, &x), "solve");
    CheckCuda(cudaDeviceSynchronize(), "solve-sync");
    const double solveSeconds = SecondsSince(start);

    linglongamg_solve_info_t info{};
    CheckLinglong(linglongamg_get_solve_info(solver, &info), "get-solve-info");
    iterations = static_cast<USI>(std::max(info.iterations, 0));

    const auto verificationStart = Clock::now();
    CheckCuda(cudaMemcpy(solution.data(), deviceSolution,
                         (rowPtr.size() - 1) * sizeof(double),
                         cudaMemcpyDeviceToHost),
              "copy-solution-to-host");
    const auto residual = ocp_gpu_solver::TrueResidualNorms(
        rowPtr, colIdx, values, rhs.data(), solution.data());
    if (!scalingBenchmark) matrix->u = solution;

    linglongamg_destroy(solver);
#ifdef OCP_LINGLONGAMG_OPAQUE_CSR_API
    linglongamg_csr_destroy(A);
#endif
    ReleaseDeviceStorage();
    const double verificationSeconds = SecondsSince(verificationStart);

    cout << fixed << setprecision(6)
         << "[OCP_LINGLONGAMG] phase=linear-system-total setup_seconds="
         << setupSeconds << " solve_seconds=" << solveSeconds
         << " total_seconds=" << SecondsSince(totalStart)
         << " iterations=" << iterations
         << " converged=" << info.converged
         << " library_absolute_residual=" << scientific << info.final_residual
         << " library_relative_residual=" << info.relative_residual
         << fixed << endl;
    cout << scientific << setprecision(16)
         << "[OCP_LINGLONGAMG] phase=final-residual iterations=" << iterations
         << " absolute_residual=" << residual.absolute
         << " relative_residual=" << residual.relative
         << " rhs_l2_norm=" << residual.rhs << endl;
    cout << fixed << "[OCP_LINGLONGAMG] phase=verification-postprocess seconds="
         << verificationSeconds << endl;

    if (!info.converged) {
        return -static_cast<OCP_INT>(iterations == 0 ? 1 : iterations);
    }
    return static_cast<OCP_INT>(iterations);
}

#endif
