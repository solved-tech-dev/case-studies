// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef FoamPressure_H
#define FoamPressure_H
#include "lduMatrix.H"
extern "C" void fpOuter(double);
extern "C" void fpFlush(const char*);
extern "C" void fpShutdown();
extern "C" void fpNativePart(int, double);
namespace Foam {
class FoamPressure : public lduMatrix::solver {
public:
    TypeName("FoamPressure");
    FoamPressure(const word&, const lduMatrix&,
        const FieldField<Field, scalar>&, const FieldField<Field, scalar>&,
        const lduInterfaceFieldPtrsList&, const dictionary&);
    solverPerformance scalarSolve(solveScalarField&, const solveScalarField&,
        const direction = 0) const override;
    solverPerformance scalarSolveOnce(solveScalarField&, const solveScalarField&,
        const direction = 0) const;
    solverPerformance solve(scalarField&, const scalarField&,
        const direction = 0) const override;
};
}
#endif
