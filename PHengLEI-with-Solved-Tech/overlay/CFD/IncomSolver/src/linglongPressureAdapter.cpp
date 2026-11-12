#include "linglongPressureAdapter.h"
#include "linglongSession.h"
#include "hypreCalculator.h"
#include "IncomPhaseTiming.h"
#include "GlobalDataBase.h"
#include "PHMpi.h"
#include "TK_Exit.h"
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <type_traits>

namespace PHSPACE
{
namespace
{
std::string StringParameter(const char* name, const char* fallback)
{
    return GlobalDataBase::IsExist(name, PHSTRING, 1)
        ? GlobalDataBase::GetStrParaFromDB(name) : fallback;
}
}
struct linglongPressureAdapter::Impl
{
    Grid* grid = nullptr;
    std::unique_ptr<hypreCalculator> fallback;
    std::unique_ptr<LinglongSession> session;
    std::ofstream log;
};
linglongPressureAdapter::linglongPressureAdapter() : impl(new Impl) {}
linglongPressureAdapter::~linglongPressureAdapter() = default;

void linglongPressureAdapter::mathLibSolveInit(Grid* grid)
{
    static_assert(std::is_same<RDouble, double>::value, "FP64 PHengLEI required");
    try
    {
        if (!impl) throw std::runtime_error("adapter has been finalized");
        if (impl->grid == grid) return; // Repeated init on the same grid keeps the session.
        if (impl->grid) throw std::runtime_error("one adapter cannot serve multiple grid zones");
        if (PHMPI::GetNumberOfProcessor() != 1 || PHMPI::GetNumberofGlobalZones() != 1)
            throw std::runtime_error("Linglong pressure adapter requires one MPI rank and one grid zone");
        UnstructGrid* mesh = UnstructGridCast(grid);
        bool pressureOutlet = false;
        for (int face = 0; face < mesh->GetNBoundFace(); ++face)
        {
            const int type = mesh->GetBCRecord()[face]->GetKey();
            if (type == PHENGLEI::INTERFACE)
                throw std::runtime_error("coupled/processor boundary requires distributed assembly, which is not supported");
            pressureOutlet = pressureOutlet || type == PHENGLEI::OUTFLOW || type == PHENGLEI::PRESSURE_OUTLET;
        }
        if (!pressureOutlet && !GlobalDataBase::IsExist("incomPressureReferenceCell", PHINT, 1))
        {
            int reference = 0;
            GlobalDataBase::UpdateData("incomPressureReferenceCell", &reference, PHINT, 1);
        }
        if (!pressureOutlet && GlobalDataBase::GetIntParaFromDB("incomPressureReferenceCell") < 0)
            throw std::runtime_error("closed pressure systems require a nonnegative incomPressureReferenceCell");
        const int device = GlobalDataBase::IsExist("linglongDevice", PHINT, 1)
            ? GlobalDataBase::GetIntParaFromDB("linglongDevice") : 0;
        impl->log.open(StringParameter("linglongLogFile", "./results/linglong_pressure.csv").c_str());
        if (!impl->log) throw std::runtime_error("cannot open Linglong pressure CSV; create results/ first");
        impl->log << "call,rows,nnz,iterations,converged,initial_l2,final_l2,rhs_l2,relative_rhs,relative_initial,sdk_relative,reset_ms,upload_ms,setup_ms,solve_ms,query_ms,download_ms,verify_ms,total_ms,sdk_setup_ms,sdk_solve_ms\n"
                  << std::setprecision(17);
        impl->session.reset(new LinglongSession(StringParameter("linglongConfigFile", "./bin/linglongamg_config.json"), device));
        IncomPhaseTiming::Add(IDX::S_IP, INCOM_BACKEND_CREATE, impl->session->CreateMilliseconds()/1000.0);
        std::cout << "LINGLONG_SESSION create=1 device=" << impl->session->DeviceDescription()
                  << " create_ms=" << impl->session->CreateMilliseconds() << std::endl;
        // Construct only after the SDK session is ready. HYPRE owns other equations.
        impl->fallback.reset(new hypreCalculator);
        impl->fallback->mathLibSolveInit(grid);
        impl->grid = grid;
    }
    catch (const std::exception& error)
    {
        if (impl) impl->session.reset(); // ExceptionExit uses MPI_Abort, without stack unwinding.
        TK_Exit::ExceptionExit(std::string("Linglong pressure initialization: ") + error.what());
    }
}

void linglongPressureAdapter::mathLibSolve(Grid* grid, int equation, RDouble* x,
    RDouble* diagonal, RDouble* rhs, RDouble* upper, RDouble* lower)
{
    try
    {
        if (!impl || impl->grid != grid || !impl->session)
            throw std::runtime_error("pressure adapter was not initialized for this grid");
        if (equation != IDX::S_IP)
        {
            impl->fallback->mathLibSolve(grid, equation, x, diagonal, rhs, upper, lower);
            return;
        }
        UnstructGrid* mesh = UnstructGridCast(grid);
        const auto packStart = std::chrono::steady_clock::now();
        const AMGBenchmarkCSR matrix = AssembleAMGBenchmarkCSR(mesh->GetNTotalCell(), mesh->GetNTotalFace(),
            mesh->GetLeftCellOfFace(), mesh->GetRightCellOfFace(), diagonal, upper, lower);
        IncomPhaseTiming::Add(equation, INCOM_BACKEND_PACK,
            std::chrono::duration<double>(std::chrono::steady_clock::now()-packStart).count());
        const LinglongSessionReport r = impl->session->Solve(matrix, rhs, x);
        IncomPhaseTiming::Add(equation, INCOM_BACKEND_ASSEMBLE, (r.resetMs+r.uploadMs)/1000.0);
        IncomPhaseTiming::Add(equation, INCOM_LINEAR_SETUP, r.setupMs/1000.0);
        IncomPhaseTiming::Add(equation, INCOM_LINEAR_SOLVE, r.solveMs/1000.0);
        IncomPhaseTiming::Add(equation, INCOM_SOLUTION_GET, (r.queryMs+r.downloadMs+r.verifyMs)/1000.0);
        impl->log << r.call << ',' << mesh->GetNTotalCell() << ',' << matrix.values.size() << ','
            << r.iterations << ',' << r.converged << ',' << r.initialL2 << ',' << r.finalL2 << ','
            << r.rhsL2 << ',' << (r.rhsL2 ? r.finalL2/r.rhsL2 : r.finalL2) << ','
            << (r.initialL2 ? r.finalL2/r.initialL2 : r.finalL2) << ',' << r.sdkRelativeResidual << ','
            << r.resetMs << ',' << r.uploadMs << ',' << r.setupMs << ',' << r.solveMs << ','
            << r.queryMs << ',' << r.downloadMs << ',' << r.verifyMs << ',' << r.totalMs << ','
            << r.sdkSetupMs << ',' << r.sdkSolveMs << '\n';
        impl->log.flush();
        if (!impl->log) throw std::runtime_error("failed to write pressure CSV");
        if (!r.converged)
            std::cerr << "Linglong pressure call " << r.call << " reached its stopping limit; true residual=" << r.finalL2 << '\n';
    }
    catch (const std::exception& error)
    {
        if (impl) impl->session.reset();
        TK_Exit::ExceptionExit(std::string("Linglong pressure solve: ") + error.what());
    }
}
void linglongPressureAdapter::mathLibSolveFinalize() { impl.reset(); }
}
