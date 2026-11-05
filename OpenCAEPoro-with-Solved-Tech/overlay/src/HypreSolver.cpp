#include "../config/config.hpp"

#ifdef OCP_USE_HYPRE

#include "HypreSolver.hpp"
#if defined(HYPRE_USING_CUDA)
#include <cuda_runtime.h>
#endif

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mutex>

#include <mpi.h>

#include "GpuSolverCommon.hpp"
#include "Spe10NativeScaling.hpp"
#include "UtilError.hpp"

namespace {

using Clock = std::chrono::steady_clock;
constexpr double kTolerance = 1.0e-6;

double SecondsSince(const Clock::time_point& start)
{
    return std::chrono::duration<double>(Clock::now() - start).count();
}

void CheckHypre(const HYPRE_Int status, const char* phase)
{
    if (status != 0) {
        OCP_ABORT(string("HYPRE ") + phase + " failed with status " +
                  std::to_string(status));
    }
}

void EnsureHypreInitialized()
{
    static std::once_flag flag;
    std::call_once(flag, []() { CheckHypre(HYPRE_Init(), "init"); });
}

} // namespace

HypreSolver::HypreSolver(const string&, const string&, const OCPMatrix&)
{
    EnsureHypreInitialized();
}

HypreSolver::~HypreSolver()
{
    ReleaseObjects();
}

void HypreSolver::ReleaseObjects()
{
    if (ijSolution != nullptr) HYPRE_IJVectorDestroy(ijSolution);
    if (ijRhs != nullptr) HYPRE_IJVectorDestroy(ijRhs);
    if (ijMatrix != nullptr) HYPRE_IJMatrixDestroy(ijMatrix);
    ijSolution = nullptr;
    ijRhs = nullptr;
    ijMatrix = nullptr;
    parSolution = nullptr;
    parRhs = nullptr;
    parMatrix = nullptr;
}

void HypreSolver::AssembleMat(OCPMatrix& mat, const Domain* domain)
{
    if (domain == nullptr || domain->cs_numproc != 1) {
        OCP_ABORT("HYPRE integration currently requires exactly one OpenCAEPoroX MPI rank");
    }

    const auto start = Clock::now();
    ReleaseObjects();
#if defined(HYPRE_USING_CUDA)
    CheckHypre(HYPRE_SetMemoryLocation(HYPRE_MEMORY_HOST), "assembly-host-memory");
    CheckHypre(HYPRE_SetExecutionPolicy(HYPRE_EXEC_HOST), "assembly-host-policy");
#endif

    const int factor = ocp_spe10_scaling::RefineFactor();
    const char* singleEnv = std::getenv("OCP_SINGLE_SOLVE_BENCHMARK");
    scalingBenchmark = factor > 1 ||
        (singleEnv != nullptr && *singleEnv != '\0' && *singleEnv != '0');
    if (scalingBenchmark) {
        auto csr = ocp_spe10_scaling::Build(mat, *domain, factor);
        rowPtr.assign(csr.rowPtr.begin(), csr.rowPtr.end());
        colIdx.swap(csr.colIdx);
        values.swap(csr.values);
        rhs.swap(csr.rhs);
        solution.swap(csr.solution);
    } else {
        ocp_gpu_solver::BuildSortedScalarCsr(mat, rowPtr, colIdx, values);
        rhs.assign(mat.b.begin(), mat.b.end());
        solution.assign(mat.u.begin(), mat.u.end());
    }

    matrix = &mat;
    const std::size_t n = rowPtr.size() - 1;
    if (n == 0 || n > static_cast<std::size_t>(std::numeric_limits<HYPRE_Int>::max())) {
        OCP_ABORT("HYPRE integration requires 0 < rows <= HYPRE_Int max");
    }

    rows.resize(n);
    rowSizes.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        rows[i] = static_cast<HYPRE_BigInt>(i);
        rowSizes[i] = rowPtr[i + 1] - rowPtr[i];
    }
    hypreColumns.assign(colIdx.begin(), colIdx.end());

    const HYPRE_BigInt last = static_cast<HYPRE_BigInt>(n - 1);
    CheckHypre(HYPRE_IJMatrixCreate(MPI_COMM_SELF, 0, last, 0, last, &ijMatrix),
               "matrix-create");
    CheckHypre(HYPRE_IJMatrixSetObjectType(ijMatrix, HYPRE_PARCSR),
               "matrix-object-type");
    // One MPI rank owns every column. Exact block sizes avoid millions of
    // per-row auxiliary allocations in the generic IJ insertion path.
    std::vector<HYPRE_Int> offdiagSizes(n, 0);
    CheckHypre(HYPRE_IJMatrixSetDiagOffdSizes(ijMatrix, rowSizes.data(), offdiagSizes.data()),
               "matrix-diag-offd-sizes");
    std::cout << "[OCP_HYPRE] phase=ij-preallocated seconds=" << SecondsSince(start) << std::endl;
    CheckHypre(HYPRE_IJMatrixInitialize(ijMatrix), "matrix-initialize");
    std::cout << "[OCP_HYPRE] phase=ij-initialized seconds=" << SecondsSince(start) << std::endl;
    CheckHypre(HYPRE_IJMatrixSetValues(
                   ijMatrix, static_cast<HYPRE_Int>(n), rowSizes.data(), rows.data(),
                   hypreColumns.data(), values.data()),
               "matrix-set-values");
    CheckHypre(HYPRE_IJMatrixAssemble(ijMatrix), "matrix-assemble");
    std::cout << "[OCP_HYPRE] phase=ij-assembled seconds=" << SecondsSince(start) << std::endl;
    CheckHypre(HYPRE_IJMatrixGetObject(ijMatrix,
               reinterpret_cast<void**>(&parMatrix)), "matrix-get-object");

    CheckHypre(HYPRE_IJVectorCreate(MPI_COMM_SELF, 0, last, &ijRhs), "rhs-create");
    CheckHypre(HYPRE_IJVectorSetObjectType(ijRhs, HYPRE_PARCSR), "rhs-object-type");
    CheckHypre(HYPRE_IJVectorInitialize(ijRhs), "rhs-initialize");
    CheckHypre(HYPRE_IJVectorSetValues(ijRhs, static_cast<HYPRE_Int>(n),
                                      rows.data(), rhs.data()), "rhs-set-values");
    CheckHypre(HYPRE_IJVectorAssemble(ijRhs), "rhs-assemble");
    CheckHypre(HYPRE_IJVectorGetObject(ijRhs,
               reinterpret_cast<void**>(&parRhs)), "rhs-get-object");

    CheckHypre(HYPRE_IJVectorCreate(MPI_COMM_SELF, 0, last, &ijSolution),
               "solution-create");
    CheckHypre(HYPRE_IJVectorSetObjectType(ijSolution, HYPRE_PARCSR),
               "solution-object-type");
    CheckHypre(HYPRE_IJVectorInitialize(ijSolution), "solution-initialize");
    CheckHypre(HYPRE_IJVectorSetValues(ijSolution, static_cast<HYPRE_Int>(n),
                                      rows.data(), solution.data()),
               "solution-set-values");
    CheckHypre(HYPRE_IJVectorAssemble(ijSolution), "solution-assemble");
    CheckHypre(HYPRE_IJVectorGetObject(ijSolution,
               reinterpret_cast<void**>(&parSolution)), "solution-get-object");

#if defined(HYPRE_USING_CUDA)
    CheckHypre(HYPRE_IJMatrixMigrate(ijMatrix, HYPRE_MEMORY_DEVICE), "matrix-to-gpu");
    CheckHypre(HYPRE_IJVectorMigrate(ijRhs, HYPRE_MEMORY_DEVICE), "rhs-to-gpu");
    CheckHypre(HYPRE_IJVectorMigrate(ijSolution, HYPRE_MEMORY_DEVICE), "solution-to-gpu");
    CheckHypre(HYPRE_SetMemoryLocation(HYPRE_MEMORY_DEVICE), "device-memory");
    CheckHypre(HYPRE_SetExecutionPolicy(HYPRE_EXEC_DEVICE), "device-policy");
    CheckHypre(HYPRE_SetSpGemmUseVendor(0), "device-native-spgemm");
    if (cudaDeviceSynchronize() != cudaSuccess) OCP_ABORT("HYPRE GPU transfer failed");
    std::cout << "[OCP_HYPRE] execution=CUDA memory=DEVICE coarsening=PMIS relax=l1-Jacobi spgemm=hypre" << std::endl;
#else
    std::cout << "[OCP_HYPRE] execution=CPU memory=HOST coarsening=PMIS relax=l1-Jacobi" << std::endl;
#endif
    HYPRE_Int major = 0;
    HYPRE_Int minor = 0;
    HYPRE_Int patch = 0;
    CheckHypre(HYPRE_VersionNumber(&major, &minor, &patch, nullptr), "version");
    std::cout << std::fixed << std::setprecision(6)
              << "[OCP_HYPRE] phase=preprocess version=" << major << "." << minor
              << "." << patch << " rows=" << n << " nnz=" << values.size()
              << " seconds=" << SecondsSince(start) << std::endl;
}

OCP_INT HypreSolver::Solve()
{
    if (parMatrix == nullptr || parRhs == nullptr || parSolution == nullptr) {
        OCP_ABORT("HYPRE Solve called before matrix assembly");
    }

    const auto totalStart = Clock::now();
    HYPRE_Solver solver = nullptr;
    HYPRE_Solver preconditioner = nullptr;

    CheckHypre(HYPRE_ParCSRGMRESCreate(MPI_COMM_SELF, &solver), "gmres-create");
    int restart = 10;
    if (const char* configured = std::getenv("OCP_HYPRE_GMRES_RESTART")) {
        char* end = nullptr;
        const long parsed = std::strtol(configured, &end, 10);
        if (end == configured || *end != '\0' || parsed < 1 || parsed > 100)
            OCP_ABORT("OCP_HYPRE_GMRES_RESTART must be between 1 and 100");
        restart = static_cast<int>(parsed);
    }
    std::cout << "[OCP_HYPRE] gmres_restart=" << restart
              << " relative_tolerance=" << kTolerance << std::endl;
    CheckHypre(HYPRE_ParCSRGMRESSetKDim(solver, restart), "gmres-restart");
    CheckHypre(HYPRE_ParCSRGMRESSetMaxIter(solver, 100), "gmres-max-iters");
    CheckHypre(HYPRE_ParCSRGMRESSetTol(solver, kTolerance), "gmres-tolerance");
    CheckHypre(HYPRE_ParCSRGMRESSetAbsoluteTol(solver, 0.0), "gmres-absolute-tolerance");
    CheckHypre(HYPRE_ParCSRGMRESSetLogging(solver, 1), "gmres-logging");
    CheckHypre(HYPRE_ParCSRGMRESSetPrintLevel(solver, 2), "gmres-print-level");

    CheckHypre(HYPRE_BoomerAMGCreate(&preconditioner), "amg-create");
    CheckHypre(HYPRE_BoomerAMGSetMaxLevels(preconditioner, 50), "amg-max-levels");
    CheckHypre(HYPRE_BoomerAMGSetMinCoarseSize(preconditioner, 2), "amg-min-coarse-size");
    CheckHypre(HYPRE_BoomerAMGSetCoarsenType(preconditioner, 8), "amg-coarsen-type");
    CheckHypre(HYPRE_BoomerAMGSetRelaxType(preconditioner, 18), "amg-relax-type");
    CheckHypre(HYPRE_BoomerAMGSetInterpType(preconditioner, 6), "amg-interp-type");
    int interpMaxElements = 4;
    if (const char* configured = std::getenv("OCP_HYPRE_INTERP_MAX_ELEMENTS")) {
        char* end = nullptr;
        const long parsed = std::strtol(configured, &end, 10);
        if (end == configured || *end != '\0' || parsed < 1 || parsed > 100)
            OCP_ABORT("OCP_HYPRE_INTERP_MAX_ELEMENTS must be between 1 and 100");
        interpMaxElements = static_cast<int>(parsed);
    }
    CheckHypre(HYPRE_BoomerAMGSetPMaxElmts(preconditioner,
               static_cast<HYPRE_Int>(interpMaxElements)), "amg-interp-max-elements");
    std::cout << "[OCP_HYPRE] interp_max_elements=" << interpMaxElements << std::endl;
    double strongThreshold = 0.25;
    if (const char* configured = std::getenv("OCP_HYPRE_STRONG_THRESHOLD")) {
        char* end = nullptr;
        strongThreshold = std::strtod(configured, &end);
        if (end == configured || *end != '\0' || strongThreshold <= 0.0 ||
            strongThreshold >= 1.0)
            OCP_ABORT("OCP_HYPRE_STRONG_THRESHOLD must be between 0 and 1");
    }
    CheckHypre(HYPRE_BoomerAMGSetStrongThreshold(preconditioner, strongThreshold),
               "amg-strong-threshold");
    std::cout << "[OCP_HYPRE] strong_threshold=" << strongThreshold << std::endl;
    int aggressiveLevels = 0;
    if (const char* configured = std::getenv("OCP_HYPRE_AGGRESSIVE_LEVELS")) {
        char* end = nullptr;
        const long parsed = std::strtol(configured, &end, 10);
        if (end == configured || *end != '\0' || parsed < 0 || parsed > 4)
            OCP_ABORT("OCP_HYPRE_AGGRESSIVE_LEVELS must be between 0 and 4");
        aggressiveLevels = static_cast<int>(parsed);
    }
    CheckHypre(HYPRE_BoomerAMGSetAggNumLevels(
                   preconditioner, static_cast<HYPRE_Int>(aggressiveLevels)),
               "amg-aggressive-levels");
    if (aggressiveLevels > 0) {
        CheckHypre(HYPRE_BoomerAMGSetNumPaths(preconditioner, 1),
                   "amg-aggressive-paths");
    }
    std::cout << "[OCP_HYPRE] aggressive_levels=" << aggressiveLevels << std::endl;
    int useRap2 = 1;
    if (const char* configured = std::getenv("OCP_HYPRE_RAP2")) {
        char* end = nullptr;
        const long parsed = std::strtol(configured, &end, 10);
        if (end == configured || *end != '\0' || parsed < 0 || parsed > 1)
            OCP_ABORT("OCP_HYPRE_RAP2 must be 0 or 1");
        useRap2 = static_cast<int>(parsed);
    }
    CheckHypre(HYPRE_BoomerAMGSetRAP2(preconditioner, useRap2), "amg-rap2");
    std::cout << "[OCP_HYPRE] rap2=" << useRap2 << std::endl;
    CheckHypre(HYPRE_BoomerAMGSetCycleType(preconditioner, 1), "amg-cycle-type");
    CheckHypre(HYPRE_BoomerAMGSetCycleRelaxType(preconditioner, 18, 1),
               "amg-down-relax-type");
    CheckHypre(HYPRE_BoomerAMGSetCycleRelaxType(preconditioner, 18, 2),
               "amg-up-relax-type");
    CheckHypre(HYPRE_BoomerAMGSetCycleRelaxType(preconditioner, 18, 3),
               "amg-coarse-relax-type");
    CheckHypre(HYPRE_BoomerAMGSetCycleNumSweeps(preconditioner, 2, 1),
               "amg-presweeps");
    CheckHypre(HYPRE_BoomerAMGSetCycleNumSweeps(preconditioner, 2, 2),
               "amg-postsweeps");
    CheckHypre(HYPRE_BoomerAMGSetCycleNumSweeps(preconditioner, 2, 3),
               "amg-coarsest-sweeps");
    std::cout << "[OCP_HYPRE] cycle=V sweeps=2/2/2" << std::endl;
    CheckHypre(HYPRE_BoomerAMGSetMaxRowSum(preconditioner, 0.9), "amg-max-row-sum");
    CheckHypre(HYPRE_BoomerAMGSetTol(preconditioner, 0.0), "amg-tolerance");
    CheckHypre(HYPRE_BoomerAMGSetMaxIter(preconditioner, 1), "amg-max-iters");
    CheckHypre(HYPRE_BoomerAMGSetPrintLevel(preconditioner, 1), "amg-print-level");
    CheckHypre(HYPRE_ParCSRGMRESSetPrecond(
                   solver,
                   reinterpret_cast<HYPRE_PtrToParSolverFcn>(HYPRE_BoomerAMGSolve),
                   reinterpret_cast<HYPRE_PtrToParSolverFcn>(HYPRE_BoomerAMGSetup),
                   preconditioner),
               "gmres-set-preconditioner");

    auto start = Clock::now();
    CheckHypre(HYPRE_ParCSRGMRESSetup(solver, parMatrix, parRhs, parSolution),
               "gmres-setup");
#if defined(HYPRE_USING_CUDA)
    if (cudaDeviceSynchronize() != cudaSuccess) OCP_ABORT("HYPRE GPU setup failed");
#endif
    const double setupSeconds = SecondsSince(start);

    start = Clock::now();
    CheckHypre(HYPRE_ParCSRGMRESSolve(solver, parMatrix, parRhs, parSolution),
               "gmres-solve");
#if defined(HYPRE_USING_CUDA)
    if (cudaDeviceSynchronize() != cudaSuccess) OCP_ABORT("HYPRE GPU solve failed");
#endif
    const double solveSeconds = SecondsSince(start);

    HYPRE_Int hypreIterations = 0;
    HYPRE_Real reportedRelative = 0.0;
    CheckHypre(HYPRE_ParCSRGMRESGetNumIterations(solver, &hypreIterations),
               "gmres-get-iterations");
    CheckHypre(HYPRE_ParCSRGMRESGetFinalRelativeResidualNorm(
                   solver, &reportedRelative), "gmres-get-relative-residual");
    iterations = static_cast<USI>(std::max<HYPRE_Int>(hypreIterations, 0));

#if defined(HYPRE_USING_CUDA)
    CheckHypre(HYPRE_IJVectorMigrate(ijSolution, HYPRE_MEMORY_HOST), "solution-to-host");
#endif
    CheckHypre(HYPRE_IJVectorGetValues(ijSolution,
               static_cast<HYPRE_Int>(rows.size()), rows.data(), solution.data()),
               "solution-get-values");
    const auto norms = ocp_gpu_solver::TrueResidualNorms(
        rowPtr, colIdx, values, rhs.data(), solution.data());
    if (!scalingBenchmark) matrix->u = solution;

    std::cout << std::fixed << std::setprecision(6)
              << "[OCP_HYPRE] phase=amg-runtime setup_seconds=" << setupSeconds
              << " solve_seconds=" << solveSeconds
              << " total_seconds=" << SecondsSince(totalStart) << std::endl;
    std::cout << std::scientific << std::setprecision(16)
              << "[OCP_HYPRE] phase=final-residual iterations=" << iterations
              << " absolute_residual=" << norms.absolute
              << " relative_residual=" << norms.relative
              << " rhs_l2_norm=" << norms.rhs
              << " hypre_reported_relative_residual=" << reportedRelative
              << std::fixed << std::endl;

    HYPRE_ParCSRGMRESDestroy(solver);
    HYPRE_BoomerAMGDestroy(preconditioner);

    const bool converged = std::isfinite(norms.relative) &&
                           norms.relative <= kTolerance;
    return converged
        ? static_cast<OCP_INT>(iterations)
        : -static_cast<OCP_INT>(iterations == 0 ? 1 : iterations);
}

#endif
