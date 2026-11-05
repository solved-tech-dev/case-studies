#include "../config/config.hpp"

#ifdef OCP_USE_HYPRE

#ifndef __HYPRESOLVER_HEADER__
#define __HYPRESOLVER_HEADER__

#include <cstdint>
#include <vector>

#include "HYPRE.h"
#include "HYPRE_IJ_mv.h"
#include "HYPRE_parcsr_ls.h"

#include "LinearSolver.hpp"

class HypreSolver : public LinearSolver
{
public:
    HypreSolver(const string& dir, const string& file, const OCPMatrix& mat);
    ~HypreSolver() override;

    void AssembleMat(OCPMatrix& mat, const Domain* domain) override;
    OCP_INT Solve() override;
    USI GetNumIters() const override { return iterations; }

private:
    void ReleaseObjects();

    std::vector<HYPRE_Int> rowPtr;
    std::vector<std::int32_t> colIdx;
    std::vector<double> values;
    std::vector<double> rhs;
    std::vector<double> solution;
    std::vector<HYPRE_BigInt> rows;
    std::vector<HYPRE_BigInt> hypreColumns;
    std::vector<HYPRE_Int> rowSizes;

    HYPRE_IJMatrix ijMatrix{nullptr};
    HYPRE_IJVector ijRhs{nullptr};
    HYPRE_IJVector ijSolution{nullptr};
    HYPRE_ParCSRMatrix parMatrix{nullptr};
    HYPRE_ParVector parRhs{nullptr};
    HYPRE_ParVector parSolution{nullptr};

    OCPMatrix* matrix{nullptr};
    bool scalingBenchmark{false};
    USI iterations{0};
};

#endif
#endif
