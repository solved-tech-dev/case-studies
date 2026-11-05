/*! \file    LinearSystem.cpp
 *  \brief   Contains Internal Matrix structure, ptr to external linearsolver and
 *interface \author  Shizhe Li \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "LinearSystem.hpp"
#include "Spe10NativeScaling.hpp"

#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <mpi.h>

namespace {

using BenchmarkClock = std::chrono::steady_clock;
const BenchmarkClock::time_point g_openCaeStart = BenchmarkClock::now();

double BenchmarkSecondsSince(const BenchmarkClock::time_point& start)
{
    return std::chrono::duration<double>(BenchmarkClock::now() - start).count();
}

bool EnvFlag(const char* name, const bool defaultValue)
{
    const char* value = std::getenv(name);
    if (value == nullptr) return defaultValue;
    const string text(value);
    return !(text.empty() || text == "0" || text == "false" || text == "FALSE");
}

void FinalizeMpiAndExit(const int status = 0)
{
    int finalized = 0;
    MPI_Finalized(&finalized);
    if (!finalized) MPI_Finalize();
    std::exit(status);
}

} // namespace

USI LinearSystem::Setup(const OCPModel& model, const string& dir, const string& file,
                        const Domain& d, const USI& nb)
{
    solveDir = dir;
    domain = &d;
    nbs.push_back(nb);

    mat.Allocate(*domain, nb);
    SetupLinearSolver(model, file);

    return LS.size() - 1;
}

void LinearSystem::SetWorkLS(const USI& i)
{
    wIndex = i;
    mat.SetBlockDim(nbs[wIndex]);
}

void LinearSystem::SetupLinearSolver(const OCPModel& model, const string& lsFile)
{
    string selectedSolver = lsFile;
    const char* solverOverride = std::getenv("OCP_LINEAR_SOLVER");
    if (solverOverride != nullptr && *solverOverride != '\0') {
        selectedSolver = solverOverride;
    }
    string lsMethod = selectedSolver;
    const auto pos = lsMethod.find_last_of('.');
    if (pos != string::npos) lsMethod = lsMethod.substr(pos + 1);
    transform(lsMethod.begin(), lsMethod.end(), lsMethod.begin(), ::tolower);

    if (false) {}
#ifdef OCP_USE_HYPRE
    else if (lsMethod == "hypre") {
        LS.push_back(new HypreSolver(solveDir, selectedSolver, mat));
        LStype.push_back(OCPLStype::hypre);
    }
#endif
#ifdef OCP_USE_LINGLONGAMG
    else if (lsMethod == "linglongamg") {
        LS.push_back(new LinglongAMGSolver(solveDir, selectedSolver, mat));
        LStype.push_back(OCPLStype::linglongamg);
    }
#endif
#ifdef OCP_USE_PARDISO
#if OCPFLOATTYPEWIDTH == 64
    else if (lsMethod == "pardiso") {
        if (mat.nb > 1) LS.push_back(new VectorPardisoSolver(solveDir, lsFile, mat));
        else            LS.push_back(new PardisoSolver(solveDir, lsFile, mat));
        LStype.push_back(OCPLStype::pardiso);
    }
#endif
#endif
#ifdef OCP_USE_SAMG
    else if (lsMethod == "samg") {
        if (mat.nb > 1) LS.push_back(new VectorSamgSolver(solveDir, lsFile, mat, model));
        else            LS.push_back(new ScalarSamgSolver(solveDir, lsFile, mat, model));
        LStype.push_back(OCPLStype::samg);
    }
#endif
#ifdef OCP_USE_FASP
    else if (lsMethod == "fasp") {
        if (mat.nb > 1) LS.push_back(new VectorFaspSolver(solveDir, lsFile, mat));
        else            LS.push_back(new ScalarFaspSolver(solveDir, lsFile, mat));
        LStype.push_back(OCPLStype::fasp);
    }
#endif
#ifdef OCP_USE_PETSCSOLVER
    else if (lsMethod == "petsc") {
        if (mat.nb > 1) LS.push_back(new VectorPetscSolver(solveDir, lsFile, mat, domain));
        else            LS.push_back(new ScalarPetscSolver(solveDir, lsFile, mat, domain));
        LStype.push_back(OCPLStype::petsc);
    }
#endif
    else {
        OCP_ABORT("Unsupported linear solver: " + lsMethod);
    }
}

void LinearSystem::AssembleMatLinearSolver()
{
#ifdef DEBUG
    mat.CheckLinearSystem();
#endif

    const auto baseReady = BenchmarkClock::now();
    const int refineFactor = ocp_spe10_scaling::RefineFactor();
    const bool singleSolveBenchmark = EnvFlag("OCP_SINGLE_SOLVE_BENCHMARK", false);
    double refinementSeconds = 0.0;
    if (refineFactor > 1 || singleSolveBenchmark) {
        const auto refinementStart = BenchmarkClock::now();
        ocp_spe10_scaling::Prepare(mat, *domain, refineFactor);
        refinementSeconds = BenchmarkSecondsSince(refinementStart);
    }
    if (domain->cs_rank == 0) {
        cout << fixed << setprecision(6)
             << "[OCP_COMMON] init_and_base_assembly_seconds="
             << std::chrono::duration<double>(baseReady - g_openCaeStart).count()
             << " matrix_generation_seconds=" << refinementSeconds
             << " common_total_seconds=" << BenchmarkSecondsSince(g_openCaeStart)
             << " factor=" << refineFactor << endl;
    }

    benchmarkCommonSeconds = BenchmarkSecondsSince(g_openCaeStart);
    const auto preprocessStart = BenchmarkClock::now();
    LS[wIndex]->AssembleMat(mat, domain);
    benchmarkPreprocessSeconds = BenchmarkSecondsSince(preprocessStart);
}

OCP_INT LinearSystem::Solve()
{
    const auto solverCallStart = BenchmarkClock::now();
    OCP_INT iters = LS[wIndex]->Solve();
    const double solverCallSeconds = BenchmarkSecondsSince(solverCallStart);
    if (domain->cs_rank == 0) {
        cout << fixed << setprecision(6)
             << "[OCP_TIMING] common_total_seconds=" << benchmarkCommonSeconds
             << " solver_preprocess_seconds=" << benchmarkPreprocessSeconds
             << " solver_call_seconds=" << solverCallSeconds
             << " solver_total_seconds="
             << benchmarkPreprocessSeconds + solverCallSeconds
             << " through_verification_seconds="
             << BenchmarkSecondsSince(g_openCaeStart) << endl;
    }

    const bool singleSolveBenchmark = EnvFlag("OCP_SINGLE_SOLVE_BENCHMARK", false);
    const char* refineFactor = std::getenv("OCP_SPE10_REFINE_FACTOR");
    if (singleSolveBenchmark ||
        (refineFactor != nullptr && std::strtol(refineFactor, nullptr, 10) > 1)) {
        const bool succeeded = iters >= 0;
        if (domain->cs_rank == 0) {
            cout << "[OCP_SPE10_SCALE] benchmark_complete=" << (succeeded ? 1 : 0)
                 << " simulator_exit_after_native_solve=1"
                 << " single_solve_benchmark=" << (singleSolveBenchmark ? 1 : 0)
                 << endl;
        }
        FinalizeMpiAndExit(succeeded ? 0 : 2);
    }

    if (iters < 0) {
        if (LStype[wIndex] == OCPLStype::fasp) {
            iters = -LS[wIndex]->GetNumIters();
        }
        if (domain->cs_rank == 0) {
            OCP_WARNING(to_string(CURRENT_RANK) + " : " +
                        to_string(domain->cs_numproc) +
                        "   linear solver failed! -- " + to_string(iters));
        }
    }

#ifdef DEBUG
    mat.CheckSolution();
#endif
    return iters;
}
