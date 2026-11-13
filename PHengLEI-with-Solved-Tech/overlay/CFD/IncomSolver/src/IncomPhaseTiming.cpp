#include "IncomPhaseTiming.h"

#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include "PHMpi.h"
#include "TK_Log.h"
#include "TK_Time.h"

namespace PHSPACE
{
namespace
{
const int equationCount = INCOM_FLOW_EQUATION + 1;
const int parentCount = 6;
const int reportPhaseCount = INCOM_PHASE_COUNT + parentCount;

double phaseTotals[equationCount][INCOM_PHASE_COUNT] = {{0.0}};
double phaseCalls[equationCount][INCOM_PHASE_COUNT] = {{0.0}};
double equationSolveCalls[equationCount] = {0.0};
bool reportWritten = false;

const char *childNames[INCOM_PHASE_COUNT] = {
    "comm_bc", "gradient", "property", "face_flux",
    "lhs_init_bc", "lhs_convection", "lhs_diffusion", "lhs_transient", "lhs_other_relax",
    "rhs_init_bc", "rhs_convection", "rhs_diffusion", "rhs_transient", "rhs_source_relax",
    "backend_create", "backend_pack", "backend_assemble",
    "linear_setup", "linear_solve", "solution_get",
    "residual", "correction", "backend_wall"
};

struct ParentDefinition
{
    const char *name;
    int firstChild;
    int lastChild;
};

const ParentDefinition parents[parentCount] = {
    {"field_prep", INCOM_COMM_BC, INCOM_FACE_FLUX},
    {"fv_lhs", INCOM_LHS_INIT_BC, INCOM_LHS_OTHER_RELAX},
    {"fv_rhs", INCOM_RHS_INIT_BC, INCOM_RHS_SOURCE_RELAX},
    {"backend_assembly", INCOM_BACKEND_CREATE, INCOM_BACKEND_ASSEMBLE},
    {"linear_solver", INCOM_LINEAR_SETUP, INCOM_SOLUTION_GET},
    {"post_update", INCOM_RESIDUAL, INCOM_CORRECTION}
};

std::string EquationName(int equation)
{
    switch (equation)
    {
        case 0: return "U";
        case 1: return "V";
        case 2: return "W";
        case 3: return "P";
        case 5: return "H";
        case 6: return "K";
        case 7: return "EPSILON";
        case 8: return "OMEGA";
        case 9: return "SA";
        case 10: return "MUT";
        case 15: return "AIR";
        case 16: return "H2";
        case 17: return "H2_L";
        case 18: return "NH3_L";
        case 19: return "CH4";
        case 20: return "C2H4";
        case 21: return "C3H8";
        case 22: return "NG";
        case 23: return "LNG";
        case 24: return "H2S";
        case 25: return "CL2";
        case 26: return "NH3";
        case 27: return "C2H6";
        case INCOM_FLOW_EQUATION: return "FLOW";
        default:
        {
            std::ostringstream name;
            name << "E" << equation;
            return name.str();
        }
    }
}

void WritePhaseLine(const std::string &equationName, const char *phaseName,
                    double totalSeconds, double maxSeconds, double calls)
{
    const unsigned long long roundedCalls = static_cast<unsigned long long>(calls + 0.5);
    const double averageMicroseconds = calls > 0.0 ? totalSeconds * 1.0e6 / calls : 0.0;

    std::ostringstream line;
    line << std::fixed << std::setprecision(9)
         << "PERF_PHASE eq=" << equationName
         << " phase=" << phaseName
         << " total_s=" << totalSeconds
         << " max_s=" << maxSeconds
         << " calls=" << roundedCalls
         << std::setprecision(3)
         << " avg_us=" << averageMicroseconds;
    WriteLogFile(line.str());
}
}

void IncomPhaseTiming::Add(int equation, IncomPhase phase, double elapsedSeconds)
{
    if (equation < 0 || equation >= equationCount || phase < 0 || phase >= INCOM_PHASE_COUNT)
    {
        return;
    }
    phaseTotals[equation][phase] += elapsedSeconds;
    phaseCalls[equation][phase] += 1.0;
}

void IncomPhaseTiming::MarkEquationSolve(int equation)
{
    if (equation >= 0 && equation < equationCount)
    {
        equationSolveCalls[equation] += 1.0;
    }
}

void IncomPhaseTiming::Report()
{
    if (reportWritten)
    {
        return;
    }
    reportWritten = true;

    const int valueCount = equationCount * reportPhaseCount;
    std::vector<double> localTotals(valueCount, 0.0);
    std::vector<double> localCalls(valueCount, 0.0);

    for (int equation = 0; equation < equationCount; ++equation)
    {
        const int base = equation * reportPhaseCount;
        for (int phase = 0; phase < INCOM_PHASE_COUNT; ++phase)
        {
            localTotals[base + phase] = phaseTotals[equation][phase];
            localCalls[base + phase] = phaseCalls[equation][phase];
        }

        for (int parent = 0; parent < parentCount; ++parent)
        {
            const int reportIndex = base + INCOM_PHASE_COUNT + parent;
            localCalls[reportIndex] = equationSolveCalls[equation];
            for (int child = parents[parent].firstChild; child <= parents[parent].lastChild; ++child)
            {
                localTotals[reportIndex] += phaseTotals[equation][child];
            }
        }
    }

    std::vector<double> summedTotals(valueCount, 0.0);
    std::vector<double> maximumTotals(valueCount, 0.0);
    std::vector<double> summedCalls(valueCount, 0.0);
    PH_Reduce(&localTotals[0], &summedTotals[0], valueCount, PH_SUM);
    PH_Reduce(&localTotals[0], &maximumTotals[0], valueCount, PH_MAX);
    PH_Reduce(&localCalls[0], &summedCalls[0], valueCount, PH_SUM);

    if (PHMPI::GetCurrentProcessorID() != PHMPI::GetServerProcessorID())
    {
        return;
    }

    const double processorCount = static_cast<double>(PHMPI::GetNumberOfProcessor());
    for (int equation = 0; equation < equationCount; ++equation)
    {
        const int base = equation * reportPhaseCount;
        bool equationUsed = false;
        for (int phase = 0; phase < INCOM_PHASE_COUNT; ++phase)
        {
            if (summedCalls[base + phase] > 0.0)
            {
                equationUsed = true;
                break;
            }
        }
        if (!equationUsed)
        {
            continue;
        }

        const std::string equationName = EquationName(equation);
        for (int parent = 0; parent < parentCount; ++parent)
        {
            const int parentIndex = base + INCOM_PHASE_COUNT + parent;
            if (summedTotals[parentIndex] <= 0.0)
            {
                continue;
            }
            WritePhaseLine(equationName, parents[parent].name,
                           summedTotals[parentIndex] / processorCount,
                           maximumTotals[parentIndex],
                           summedCalls[parentIndex] / processorCount);

            for (int child = parents[parent].firstChild; child <= parents[parent].lastChild; ++child)
            {
                const int childIndex = base + child;
                if (summedCalls[childIndex] <= 0.0)
                {
                    continue;
                }
                WritePhaseLine(equationName, childNames[child],
                               summedTotals[childIndex] / processorCount,
                               maximumTotals[childIndex],
                               summedCalls[childIndex] / processorCount);
            }
        }
        const int backendIndex = base + INCOM_BACKEND_WALL;
        if (summedCalls[backendIndex] > 0.0)
        {
            WritePhaseLine(equationName, "backend_wall",
                           summedTotals[backendIndex] / processorCount,
                           maximumTotals[backendIndex],
                           summedCalls[backendIndex] / processorCount);
        }
    }
}

IncomPhaseScope::IncomPhaseScope(int equationIn, IncomPhase phaseIn)
    : equation(equationIn), phase(phaseIn), startTime(TIME_SPACE::GetWallTime())
{
}

IncomPhaseScope::~IncomPhaseScope()
{
    IncomPhaseTiming::Add(equation, phase, TIME_SPACE::GetWallTime() - startTime);
}

}
