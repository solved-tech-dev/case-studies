#pragma once

#include "amgBenchmarkLinearSystem.h"
#include <memory>
#include <string>

namespace PHSPACE
{
struct LinglongAMGSessionReport
{
    unsigned long long call = 0;
    int iterations = 0;
    bool converged = false;
    double initialL2 = 0, finalL2 = 0, rhsL2 = 0, sdkRelativeResidual = 0;
    double resetMs = 0, uploadMs = 0, setupMs = 0, solveMs = 0;
    double queryMs = 0, downloadMs = 0, verifyMs = 0, totalMs = 0;
    double sdkSetupMs = 0, sdkSolveMs = 0;
};

// One instance owns one SDK solver handle. Matrix setup is rebuilt on each
// call; every result requires a successful SDK authorization/commit query.
// Not thread-safe. One selected GPU, host FP64 data and int32 CSR indices.
class LinglongAMGSession
{
public:
    LinglongAMGSession(const std::string& config, int device,
                    const std::string& library = std::string());
    ~LinglongAMGSession();
    LinglongAMGSessionReport Solve(const AMGBenchmarkCSR& matrix,
                               const double* rhs, double* solution);
    double CreateMilliseconds() const;
    const std::string& DeviceDescription() const;
    void Close();

private:
    struct Impl;
    std::unique_ptr<Impl> impl;
    LinglongAMGSession(const LinglongAMGSession&) = delete;
    LinglongAMGSession& operator=(const LinglongAMGSession&) = delete;
};
}
