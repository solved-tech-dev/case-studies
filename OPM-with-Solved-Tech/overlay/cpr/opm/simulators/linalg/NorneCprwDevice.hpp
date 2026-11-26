#pragma once
#include <cstdint>
#include <memory>
#include <span>

namespace Opm::detail {
struct CprwDevice;
struct CprwDelete { void operator()(CprwDevice*) const; };
using CprwHandle = std::unique_ptr<CprwDevice, CprwDelete>;

struct CprwMatrixData {
    std::span<const uint64_t> rows, columns;
    std::span<const int> pressureRows, pressureColumns;
    int pressureIndex;
};
struct CprwNumericData {
    std::span<const double> matrix, weights, pressure;
};
struct CprwResult {
    bool converged;
    unsigned iterations;
    double residual, seconds;
};

CprwHandle createCprwDevice(const CprwMatrixData&, const char* amg, const char* amgConfig, double tolerance,
                           int maxIterations);
void prepareCprwDevice(CprwDevice&, const CprwNumericData&, void* well,
                       int (*wellApply)(void*, const double*, double*));
CprwResult solveCprwDevice(CprwDevice&, std::span<const double> rhs, std::span<double> solution);
void auditCprwOperator(CprwDevice&, std::span<const double> x, std::span<double> result);
void auditCprwCorrection(CprwDevice&, std::span<const double> rhs, std::span<double> coarseRhs,
                         std::span<double> coarseSolution, std::span<double> result);
}
