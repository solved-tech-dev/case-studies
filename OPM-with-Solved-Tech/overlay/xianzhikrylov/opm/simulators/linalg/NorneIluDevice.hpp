#pragma once

#include <cstdint>
#include <memory>
#include <span>

namespace Opm::detail {

struct IluDevice;
struct IluDeviceDelete
{
    void operator()(IluDevice*) const;
};
using IluDeviceHandle = std::unique_ptr<IluDevice, IluDeviceDelete>;

struct IluPrepareTiming
{
    double valuesUpload;
    double factor;
};

struct IluResult
{
    bool converged;
    unsigned iterations;
    double relativeResidual;
    double solve;
    double rhsUpload;
    double solutionDownload;
    double fineSpmv;
    double well;
    double iluApply;
    unsigned matvecs;
};

IluDeviceHandle createIluDevice(std::span<const uint64_t> rows,
                                std::span<const uint64_t> columns,
                                double tolerance,
                                int maxIterations);
IluPrepareTiming prepareIluDevice(IluDevice&,
                                  std::span<const double> values,
                                  void* well,
                                  int (*wellApply)(void*, const double*, double*));
IluResult solveIluDevice(IluDevice&, std::span<const double> rhs, std::span<double> solution);

} // namespace Opm::detail
