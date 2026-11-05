#ifndef __GPUSOLVERCOMMON_HEADER__
#define __GPUSOLVERCOMMON_HEADER__

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <vector>

#include "OCPMatrix.hpp"
#include "UtilError.hpp"

namespace ocp_gpu_solver {

struct ResidualNorms
{
    double absolute{0.0};
    double relative{0.0};
    double rhs{0.0};
};

template <typename RowPtr>
void BuildSortedScalarCsr(const OCPMatrix& mat,
                          std::vector<RowPtr>& rowPtr,
                          std::vector<std::int32_t>& colIdx,
                          std::vector<double>& values)
{
    const std::uint64_t n64 = static_cast<std::uint64_t>(mat.dim) * mat.nb;
    std::uint64_t nnz64 = 0;
    for (OCP_USI row = 0; row < mat.dim; ++row) {
        nnz64 += static_cast<std::uint64_t>(mat.colId[row].size()) * mat.nb2;
    }

    if (n64 > static_cast<std::uint64_t>(std::numeric_limits<std::int32_t>::max())) {
        OCP_ABORT("GPU solver scalar matrix dimension exceeds 32-bit column-index range");
    }
    if (nnz64 > static_cast<std::uint64_t>(std::numeric_limits<RowPtr>::max())) {
        OCP_ABORT("GPU solver scalar matrix NNZ exceeds row-pointer range");
    }

    rowPtr.resize(static_cast<std::size_t>(n64) + 1);
    colIdx.resize(static_cast<std::size_t>(nnz64));
    values.resize(static_cast<std::size_t>(nnz64));

    std::uint64_t cursor = 0;
    std::vector<std::size_t> order;
    for (OCP_USI blockRow = 0; blockRow < mat.dim; ++blockRow) {
        const std::size_t blockNnz = mat.colId[blockRow].size();
        order.resize(blockNnz);
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](std::size_t lhs, std::size_t rhs) {
            return mat.colId[blockRow][lhs] < mat.colId[blockRow][rhs];
        });

        for (USI r = 0; r < mat.nb; ++r) {
            const std::uint64_t scalarRow = static_cast<std::uint64_t>(blockRow) * mat.nb + r;
            rowPtr[static_cast<std::size_t>(scalarRow)] = static_cast<RowPtr>(cursor);
            for (const std::size_t entry : order) {
                const OCP_USI blockCol = mat.colId[blockRow][entry];
                if (blockCol >= mat.dim) {
                    OCP_ABORT("Single-rank GPU solver encountered a nonlocal/ghost matrix column");
                }
                const double* block = mat.val[blockRow].data() + entry * mat.nb2;
                for (USI c = 0; c < mat.nb; ++c) {
                    colIdx[static_cast<std::size_t>(cursor)] =
                        static_cast<std::int32_t>(static_cast<std::uint64_t>(blockCol) * mat.nb + c);
                    values[static_cast<std::size_t>(cursor)] = block[r * mat.nb + c];
                    ++cursor;
                }
            }
        }
    }
    rowPtr[static_cast<std::size_t>(n64)] = static_cast<RowPtr>(cursor);
}

template <typename RowPtr>
ResidualNorms TrueResidualNorms(const std::vector<RowPtr>& rowPtr,
                                const std::vector<std::int32_t>& colIdx,
                                const std::vector<double>& values,
                                const double* rhs,
                                const double* solution)
{
    long double numerator = 0.0L;
    long double denominator = 0.0L;
    const std::size_t n = rowPtr.size() - 1;
    for (std::size_t row = 0; row < n; ++row) {
        long double ax = 0.0L;
        for (std::uint64_t p = static_cast<std::uint64_t>(rowPtr[row]);
             p < static_cast<std::uint64_t>(rowPtr[row + 1]); ++p) {
            ax += static_cast<long double>(values[static_cast<std::size_t>(p)]) *
                  solution[static_cast<std::size_t>(colIdx[static_cast<std::size_t>(p)])];
        }
        const long double diff = ax - rhs[row];
        numerator += diff * diff;
        denominator += static_cast<long double>(rhs[row]) * rhs[row];
    }
    ResidualNorms result;
    result.absolute = static_cast<double>(std::sqrt(numerator));
    result.rhs = static_cast<double>(std::sqrt(denominator));
    result.relative = denominator > 0.0L
        ? static_cast<double>(std::sqrt(numerator / denominator))
        : result.absolute;
    return result;
}

template <typename RowPtr>
double TrueRelativeResidual(const std::vector<RowPtr>& rowPtr,
                            const std::vector<std::int32_t>& colIdx,
                            const std::vector<double>& values,
                            const double* rhs,
                            const double* solution)
{
    return TrueResidualNorms(rowPtr, colIdx, values, rhs, solution).relative;
}

} // namespace ocp_gpu_solver

#endif
