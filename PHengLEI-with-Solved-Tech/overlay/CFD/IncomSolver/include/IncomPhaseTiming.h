#pragma once

namespace PHSPACE
{

enum IncomPhase
{
    INCOM_COMM_BC = 0,
    INCOM_GRADIENT,
    INCOM_PROPERTY,
    INCOM_FACE_FLUX,

    INCOM_LHS_INIT_BC,
    INCOM_LHS_CONVECTION,
    INCOM_LHS_DIFFUSION,
    INCOM_LHS_TRANSIENT,
    INCOM_LHS_OTHER_RELAX,

    INCOM_RHS_INIT_BC,
    INCOM_RHS_CONVECTION,
    INCOM_RHS_DIFFUSION,
    INCOM_RHS_TRANSIENT,
    INCOM_RHS_SOURCE_RELAX,

    INCOM_BACKEND_CREATE,
    INCOM_BACKEND_PACK,
    INCOM_BACKEND_ASSEMBLE,

    INCOM_LINEAR_SETUP,
    INCOM_LINEAR_SOLVE,
    INCOM_SOLUTION_GET,

    INCOM_RESIDUAL,
    INCOM_CORRECTION,

    INCOM_BACKEND_WALL,

    INCOM_PHASE_COUNT
};

const int INCOM_FLOW_EQUATION = 50;

class IncomPhaseTiming
{
public:
    static void Add(int equation, IncomPhase phase, double elapsedSeconds);
    static void MarkEquationSolve(int equation);
    static void Report();
};

class IncomPhaseScope
{
public:
    IncomPhaseScope(int equation, IncomPhase phase);
    ~IncomPhaseScope();

private:
    int equation;
    IncomPhase phase;
    double startTime;
};

}
