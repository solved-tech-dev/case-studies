// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef BenchLinear_H
#define BenchLinear_H
#include "lduMatrix.H"
extern "C" void benchFlush(const char*);
extern "C" double benchLinearSeconds();
extern "C" void benchSetOuter(double);
extern "C" void benchInitialize();
extern "C" void benchReleaseGPU();
namespace Foam {
class BenchLinear : public lduMatrix::solver {
public:
    TypeName("BenchLinear");
    BenchLinear(const word&, const lduMatrix&, const FieldField<Field, scalar>&,
        const FieldField<Field, scalar>&, const lduInterfaceFieldPtrsList&,
        const dictionary&);
    solverPerformance scalarSolve(solveScalarField&, const solveScalarField&,
        const direction = 0) const override;
    solverPerformance solve(scalarField&, const scalarField&,
        const direction = 0) const override;
};
}
#endif
