#ifndef __SPE10_NATIVE_SCALING_HEADER__
#define __SPE10_NATIVE_SCALING_HEADER__

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include "OCPMatrix.hpp"
#include "Spe10BenchmarkGrid.hpp"
#include "UtilError.hpp"

namespace ocp_spe10_scaling {

struct Csr
{
    std::vector<std::int32_t> rowPtr;
    std::vector<std::int32_t> colIdx;
    std::vector<double> values;
    std::vector<double> rhs;
    std::vector<double> solution;
};

inline double SixSignificant(const double value)
{
    if (value == 0.0 || !std::isfinite(value)) return value;
    char text[64];
    const int length = std::snprintf(text, sizeof(text), "%.6g", value);
    if (length <= 0 || static_cast<std::size_t>(length) >= sizeof(text)) {
        OCP_ABORT("failed to format SPE10 coefficient");
    }
    return std::strtod(text, nullptr);
}

inline int RefineFactor()
{
    const char* text = std::getenv("OCP_SPE10_REFINE_FACTOR");
    if (text == nullptr || *text == '\0') return 1;
    char* end = nullptr;
    const long factor = std::strtol(text, &end, 10);
    if (end == text || *end != '\0' || factor < 1 || factor > 5) {
        OCP_ABORT("OCP_SPE10_REFINE_FACTOR must be in [1,5]");
    }
    return static_cast<int>(factor);
}

inline double RefinedTimeRatio(const int factor)
{
    if (factor == 5) {
        const char* extension = std::getenv("OCP_SPE10_R5_USE_R4_TIME_RATIO");
        if (extension == nullptr || std::string(extension) != "1") {
            OCP_ABORT("R5 requires OCP_SPE10_R5_USE_R4_TIME_RATIO=1");
        }
        std::cout << "[OCP_SPE10_VARIANT] factor=5"
                  << " time_ratio=0.0309808932186387"
                  << " source=R4-fixed-time-ratio-extension"
                  << " native_r5_time_step=0" << std::endl;
        return 0.0309808932186387;
    }

    static constexpr double ratio[] = {
        0.0,
        1.0,
        0.181030563397388,
        0.0628794000204223,
        0.0309808932186387
    };
    return ratio[factor];
}

inline Csr BuildLightweightRefined(const OCPMatrix& mat, const Domain& domain,
                                   const int factor)
{
    const auto& source = ocp_spe10_benchmark_grid::Data();
    if (source.nx != 60 || source.ny != 220 || source.nz != 85 ||
        source.kx.empty() || source.ky.empty() || source.kz.empty()) {
        OCP_ABORT("SPE10 base-grid properties were not captured");
    }

    const std::uint64_t baseRows = domain.GetNumGridInterior();
    if (baseRows != 1094421ULL || baseRows > mat.dim) {
        OCP_ABORT("SPE10 refinement requires the assembled r=1 pressure matrix");
    }
    const auto& grid = domain.GetGridFull();
    if (grid.size() < baseRows) OCP_ABORT("SPE10 active-grid map is incomplete");

    std::vector<std::uint64_t> order(baseRows);
    std::iota(order.begin(), order.end(), 0ULL);
    std::sort(order.begin(), order.end(), [&](std::uint64_t a, std::uint64_t b) {
        return grid[static_cast<std::size_t>(a)] < grid[static_cast<std::size_t>(b)];
    });
    std::vector<std::uint64_t> sortedRank(baseRows);
    const std::uint64_t baseFullCount =
        static_cast<std::uint64_t>(source.nx) * source.ny * source.nz;
    std::vector<std::int32_t> baseFullToRank(baseFullCount, -1);
    for (std::uint64_t rank = 0; rank < baseRows; ++rank) {
        const auto local = order[static_cast<std::size_t>(rank)];
        sortedRank[static_cast<std::size_t>(local)] = rank;
        baseFullToRank[static_cast<std::size_t>(
            grid[static_cast<std::size_t>(local)])] = static_cast<std::int32_t>(rank);
    }

    std::vector<double> accumulation(baseRows, 0.0);
    std::vector<double> cellScale(baseRows, 0.0);
    std::vector<std::array<double, 6>> boundary(baseRows);
    std::uint64_t directedBaseEdges = 0;
    double globalScaleSum = 0.0;
    std::uint64_t globalScaleCount = 0;
    const std::int64_t nx = source.nx;
    const std::int64_t nxy = static_cast<std::int64_t>(source.nx) * source.ny;

    const auto direction = [&](const std::int64_t delta) {
        if (delta == 1) return 0;
        if (delta == -1) return 1;
        if (delta == nx) return 2;
        if (delta == -nx) return 3;
        if (delta == nxy) return 4;
        if (delta == -nxy) return 5;
        return -1;
    };

    const auto rockTrans = [&](const std::uint64_t a, const std::uint64_t b,
                               const int dir) {
        double areaA = 0.0;
        double areaB = 0.0;
        double permA = 0.0;
        double permB = 0.0;
        if (dir < 2) {
            areaA = 2.0 * source.dy[a] * source.dz[a] / source.dx[a];
            areaB = 2.0 * source.dy[b] * source.dz[b] / source.dx[b];
            permA = source.kx[a];
            permB = source.kx[b];
        } else if (dir < 4) {
            areaA = 2.0 * source.dz[a] * source.dx[a] / source.dy[a];
            areaB = 2.0 * source.dz[b] * source.dx[b] / source.dy[b];
            permA = source.ky[a];
            permB = source.ky[b];
        } else {
            areaA = 2.0 * source.dx[a] * source.dy[a] / source.dz[a];
            areaB = 2.0 * source.dx[b] * source.dy[b] / source.dz[b];
            permA = source.kz[a];
            permB = source.kz[b];
        }
        return 1.0 / (1.0 / (permA * areaA) + 1.0 / (permB * areaB));
    };

    for (std::uint64_t rank = 0; rank < baseRows; ++rank) {
        const auto local = order[static_cast<std::size_t>(rank)];
        const auto full = grid[static_cast<std::size_t>(local)];
        const auto& cols = mat.colId[static_cast<std::size_t>(local)];
        const auto& vals = mat.val[static_cast<std::size_t>(local)];
        double wellContribution = 0.0;
        for (std::size_t p = 0; p < cols.size(); ++p) {
            if (cols[p] >= baseRows) wellContribution += vals[p];
        }

        double rowSum = 0.0;
        double scaleSum = 0.0;
        std::uint64_t scaleCount = 0;
        for (std::size_t p = 0; p < cols.size(); ++p) {
            if (cols[p] >= baseRows) continue;
            double value = vals[p];
            if (cols[p] == local) value += wellContribution;
            value = SixSignificant(value);
            rowSum += value;
            if (cols[p] == local) continue;

            const auto neighborFull = grid[static_cast<std::size_t>(cols[p])];
            const int dir = direction(static_cast<std::int64_t>(neighborFull) -
                                      static_cast<std::int64_t>(full));
            if (dir < 0) OCP_ABORT("non-orthogonal edge in SPE10 pressure matrix");
            const double weight = std::abs(value);
            boundary[static_cast<std::size_t>(rank)][static_cast<std::size_t>(dir)] =
                weight;
            const double trans = rockTrans(full, neighborFull, dir);
            scaleSum += weight / trans;
            ++scaleCount;
            globalScaleSum += weight / trans;
            ++globalScaleCount;
            ++directedBaseEdges;
        }
        accumulation[static_cast<std::size_t>(rank)] = rowSum;
        if (scaleCount != 0) {
            cellScale[static_cast<std::size_t>(rank)] = scaleSum / scaleCount;
        }
    }
    if (globalScaleCount == 0) OCP_ABORT("SPE10 pressure matrix has no inter-cell edges");
    const double globalScale = globalScaleSum / globalScaleCount;
    for (auto& value : cellScale) if (value == 0.0) value = globalScale;

    const std::uint64_t r = static_cast<std::uint64_t>(factor);
    const std::uint64_t r2 = r * r;
    const std::uint64_t r3 = r2 * r;
    const std::uint64_t refinedNx = static_cast<std::uint64_t>(source.nx) * r;
    const std::uint64_t refinedNy = static_cast<std::uint64_t>(source.ny) * r;
    const std::uint64_t refinedNz = static_cast<std::uint64_t>(source.nz) * r;
    const std::uint64_t refinedNxy = refinedNx * refinedNy;
    const std::uint64_t refinedFullCount = refinedNxy * refinedNz;
    const std::uint64_t rows = baseRows * r3;
    const std::uint64_t internalDirected = baseRows * 6ULL * (r - 1ULL) * r2;
    const std::uint64_t nnz = rows + internalDirected + directedBaseEdges * r2;
    if (rows > static_cast<std::uint64_t>(std::numeric_limits<std::int32_t>::max()) ||
        nnz > static_cast<std::uint64_t>(std::numeric_limits<std::int32_t>::max())) {
        OCP_ABORT("refined SPE10 matrix exceeds 32-bit CSR limits");
    }

    std::vector<std::int32_t> refinedFullToRank(refinedFullCount, -1);
    std::uint64_t rankCursor = 0;
    for (std::uint64_t k = 0; k < refinedNz; ++k) {
        for (std::uint64_t j = 0; j < refinedNy; ++j) {
            for (std::uint64_t i = 0; i < refinedNx; ++i) {
                const std::uint64_t parent = (k / r) * source.nx * source.ny +
                    (j / r) * source.nx + i / r;
                if (baseFullToRank[static_cast<std::size_t>(parent)] >= 0) {
                    refinedFullToRank[static_cast<std::size_t>(
                        k * refinedNxy + j * refinedNx + i)] =
                        static_cast<std::int32_t>(rankCursor++);
                }
            }
        }
    }
    if (rankCursor != rows) OCP_ABORT("refined SPE10 active-row count mismatch");

    Csr out;
    out.rowPtr.resize(static_cast<std::size_t>(rows) + 1);
    out.colIdx.reserve(static_cast<std::size_t>(nnz));
    out.values.reserve(static_cast<std::size_t>(nnz));
    out.rhs.resize(static_cast<std::size_t>(rows));
    out.solution.assign(static_cast<std::size_t>(rows), 0.0);
    const double timeRatio = RefinedTimeRatio(factor);

    for (std::uint64_t full = 0; full < refinedFullCount; ++full) {
        const auto row = refinedFullToRank[static_cast<std::size_t>(full)];
        if (row < 0) continue;
        const std::uint64_t i = full % refinedNx;
        const std::uint64_t j = (full / refinedNx) % refinedNy;
        const std::uint64_t k = full / refinedNxy;
        const std::uint64_t parentFull = (k / r) * source.nx * source.ny +
            (j / r) * source.nx + i / r;
        const auto parentRank =
            baseFullToRank[static_cast<std::size_t>(parentFull)];

        std::vector<std::pair<std::int32_t, double>> entries;
        entries.reserve(7);
        double diagonal = accumulation[static_cast<std::size_t>(parentRank)] /
                          static_cast<double>(r3);
        const auto addNeighbor = [&](const std::uint64_t neighborFull,
                                     const int dir) {
            const auto col =
                refinedFullToRank[static_cast<std::size_t>(neighborFull)];
            if (col < 0) return;
            const std::uint64_t ni = neighborFull % refinedNx;
            const std::uint64_t nj = (neighborFull / refinedNx) % refinedNy;
            const std::uint64_t nk = neighborFull / refinedNxy;
            const std::uint64_t neighborParent =
                (nk / r) * source.nx * source.ny +
                (nj / r) * source.nx + ni / r;

            double weight = 0.0;
            if (neighborParent == parentFull) {
                double trans = 0.0;
                if (dir < 2) {
                    trans = source.kx[parentFull] * source.dy[parentFull] *
                            source.dz[parentFull] / (source.dx[parentFull] * r);
                } else if (dir < 4) {
                    trans = source.ky[parentFull] * source.dz[parentFull] *
                            source.dx[parentFull] / (source.dy[parentFull] * r);
                } else {
                    trans = source.kz[parentFull] * source.dx[parentFull] *
                            source.dy[parentFull] / (source.dz[parentFull] * r);
                }
                weight = trans *
                    cellScale[static_cast<std::size_t>(parentRank)] * timeRatio;
            } else {
                weight =
                    boundary[static_cast<std::size_t>(parentRank)]
                            [static_cast<std::size_t>(dir)] *
                    timeRatio / static_cast<double>(r);
            }
            entries.emplace_back(col, -weight);
            diagonal += weight;
        };

        if (i > 0) addNeighbor(full - 1, 1);
        if (i + 1 < refinedNx) addNeighbor(full + 1, 0);
        if (j > 0) addNeighbor(full - refinedNx, 3);
        if (j + 1 < refinedNy) addNeighbor(full + refinedNx, 2);
        if (k > 0) addNeighbor(full - refinedNxy, 5);
        if (k + 1 < refinedNz) addNeighbor(full + refinedNxy, 4);
        entries.emplace_back(row, diagonal);
        std::sort(entries.begin(), entries.end(),
                  [](const auto& a, const auto& b) {
                      return a.first < b.first;
                  });

        out.rowPtr[static_cast<std::size_t>(row)] =
            static_cast<std::int32_t>(out.colIdx.size());
        double rowSum = 0.0;
        for (const auto& entry : entries) {
            out.colIdx.push_back(entry.first);
            out.values.push_back(entry.second);
            rowSum += entry.second;
        }
        out.rhs[static_cast<std::size_t>(row)] = rowSum;
    }
    out.rowPtr.back() = static_cast<std::int32_t>(out.colIdx.size());
    if (out.colIdx.size() != static_cast<std::size_t>(nnz)) {
        OCP_ABORT("lightweight refined SPE10 NNZ accounting mismatch");
    }

    std::cout << "[OCP_SPE10_SCALE] factor=" << factor
              << " source=OpenCAEPoroX-lightweight-physical-refinement"
              << " rows=" << rows << " nnz=" << nnz
              << " known_solution=ones" << std::endl;
    return out;
}

inline Csr BuildFresh(const OCPMatrix& mat, const Domain& domain,
                      const int factor)
{
    if (mat.nb != 1) OCP_ABORT("SPE10 generation requires scalar IMPEC");
    const std::uint64_t baseRows = domain.GetNumGridInterior();
    const std::uint64_t r = static_cast<std::uint64_t>(factor);
    const std::uint64_t expectedRows = 1094421ULL * r * r * r;

    if (factor > 1 && baseRows == 1094421ULL) {
        return BuildLightweightRefined(mat, domain, factor);
    }
    if (baseRows != expectedRows || baseRows > mat.dim) {
        OCP_ABORT("SPE10 generation found an unexpected active-grid count");
    }
    if (baseRows >
        static_cast<std::uint64_t>(std::numeric_limits<std::int32_t>::max())) {
        OCP_ABORT("SPE10 row count exceeds the 32-bit solver index range");
    }

    const auto& grid = domain.GetGridFull();
    if (grid.size() < baseRows) OCP_ABORT("SPE10 active-grid map is incomplete");
    std::vector<std::uint64_t> order(baseRows);
    std::iota(order.begin(), order.end(), 0ULL);
    std::sort(order.begin(), order.end(), [&](std::uint64_t a, std::uint64_t b) {
        return grid[static_cast<std::size_t>(a)] < grid[static_cast<std::size_t>(b)];
    });
    std::vector<std::uint64_t> sortedRank(baseRows);
    for (std::uint64_t rank = 0; rank < baseRows; ++rank) {
        sortedRank[static_cast<std::size_t>(
            order[static_cast<std::size_t>(rank)])] = rank;
    }

    std::uint64_t nnz = 0;
    for (std::uint64_t sorted = 0; sorted < baseRows; ++sorted) {
        const auto baseRow = order[static_cast<std::size_t>(sorted)];
        for (const auto col : mat.colId[static_cast<std::size_t>(baseRow)]) {
            nnz += col < baseRows;
        }
    }
    if (nnz >
        static_cast<std::uint64_t>(std::numeric_limits<std::int32_t>::max())) {
        OCP_ABORT("SPE10 NNZ exceeds the 32-bit CSR pointer range");
    }

    Csr out;
    out.rowPtr.resize(static_cast<std::size_t>(baseRows) + 1);
    out.colIdx.reserve(static_cast<std::size_t>(nnz));
    out.values.reserve(static_cast<std::size_t>(nnz));
    out.rhs.resize(static_cast<std::size_t>(baseRows));
    out.solution.assign(static_cast<std::size_t>(baseRows), 0.0);

    for (std::uint64_t sorted = 0; sorted < baseRows; ++sorted) {
        const auto baseRow = order[static_cast<std::size_t>(sorted)];
        const auto& cols = mat.colId[static_cast<std::size_t>(baseRow)];
        const auto& vals = mat.val[static_cast<std::size_t>(baseRow)];
        double wellContribution = 0.0;
        for (std::size_t p = 0; p < cols.size(); ++p) {
            if (cols[p] >= baseRows) wellContribution += vals[p];
        }

        std::vector<std::pair<std::int32_t, double>> entries;
        entries.reserve(cols.size());
        for (std::size_t p = 0; p < cols.size(); ++p) {
            if (cols[p] >= baseRows) continue;
            double value = vals[p];
            if (cols[p] == baseRow) value += wellContribution;
            entries.emplace_back(
                static_cast<std::int32_t>(
                    sortedRank[static_cast<std::size_t>(cols[p])]),
                SixSignificant(value));
        }
        std::sort(entries.begin(), entries.end(),
                  [](const auto& a, const auto& b) {
                      return a.first < b.first;
                  });
        out.rowPtr[static_cast<std::size_t>(sorted)] =
            static_cast<std::int32_t>(out.colIdx.size());
        double rowSum = 0.0;
        for (const auto& entry : entries) {
            out.colIdx.push_back(entry.first);
            out.values.push_back(entry.second);
            rowSum += entry.second;
        }
        out.rhs[static_cast<std::size_t>(sorted)] = rowSum;
    }
    out.rowPtr.back() = static_cast<std::int32_t>(out.colIdx.size());
    if (out.colIdx.size() != static_cast<std::size_t>(nnz)) {
        OCP_ABORT("SPE10 NNZ accounting mismatch");
    }

    std::cout << "[OCP_SPE10_SCALE] factor=" << factor
              << " source=OpenCAEPoroX-assembled-pressure"
              << " rows=" << baseRows << " nnz=" << nnz
              << " known_solution=ones" << std::endl;
    return out;
}

inline std::unique_ptr<Csr>& PreparedCache()
{
    static std::unique_ptr<Csr> cache;
    return cache;
}

inline void Prepare(const OCPMatrix& mat, const Domain& domain,
                    const int factor)
{
    PreparedCache() = std::make_unique<Csr>(BuildFresh(mat, domain, factor));
}

inline Csr Build(const OCPMatrix& mat, const Domain& domain, const int factor)
{
    auto& cache = PreparedCache();
    if (cache) {
        Csr out = std::move(*cache);
        cache.reset();
        return out;
    }
    return BuildFresh(mat, domain, factor);
}

} // namespace ocp_spe10_scaling

#endif
