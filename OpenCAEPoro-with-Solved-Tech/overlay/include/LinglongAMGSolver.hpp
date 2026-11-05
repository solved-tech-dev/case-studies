#include "../config/config.hpp"

#ifdef OCP_USE_LINGLONGAMG

#ifndef __LINGLONGAMGSOLVER_HEADER__
#define __LINGLONGAMGSOLVER_HEADER__

#include <cstdint>
#include <string>
#include <vector>

#include "LinearSolver.hpp"
#include "linglongamg.h"

class LinglongAMGSolver final : public LinearSolver
{
public:
    LinglongAMGSolver(const string& dir, const string& file, const OCPMatrix& mat);

    void AssembleMat(OCPMatrix& mat, const Domain* domain) override;
    OCP_INT Solve() override;
    USI GetNumIters() const override { return iterations; }

private:
    void ReleaseDeviceStorage();

    string configPath;
    std::vector<std::int32_t> rowPtr;
    std::vector<std::int32_t> colIdx;
    std::vector<double> values;
    std::vector<double> rhs;
    std::vector<double> solution;
    bool scalingBenchmark{false};
    OCPMatrix* matrix{nullptr};
    int* deviceRowPtr{nullptr};
    int* deviceColIdx{nullptr};
    double* deviceValues{nullptr};
    double* deviceRhs{nullptr};
    double* deviceSolution{nullptr};
    USI iterations{0};
};

#endif
#endif
