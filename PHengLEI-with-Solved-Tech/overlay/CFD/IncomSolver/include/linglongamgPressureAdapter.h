#pragma once

#include "LinearEquationSystemCalculator.h"
#include <memory>

namespace PHSPACE
{
// Pressure-only GPU adapter. Other equations retain the CPU HYPRE backend.
class linglongamgPressureAdapter : public LinearEquationSystemCalculator
{
public:
    linglongamgPressureAdapter();
    ~linglongamgPressureAdapter() override;
    void mathLibSolveInit(Grid* grid) override;
    void mathLibSolve(Grid* grid, int equation, RDouble* x, RDouble* diagonal,
                     RDouble* rhs, RDouble* upper, RDouble* lower) override;
    void mathLibSolveFinalize() override;
private:
    struct Impl;
    std::unique_ptr<Impl> impl;
};
}
