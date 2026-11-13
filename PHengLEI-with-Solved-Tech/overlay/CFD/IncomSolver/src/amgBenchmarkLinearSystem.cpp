#include "amgBenchmarkLinearSystem.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace PHSPACE
{
AMGBenchmarkCSR AssembleAMGBenchmarkCSR(int cells, int faces, const int* left,
    const int* right, const double* diagonal, const double* upper, const double* lower)
{
    if (cells <= 0 || faces < 0 || !diagonal || (faces && (!left || !right || !upper || !lower)))
        throw std::runtime_error("invalid finite-volume matrix arrays");
    // Upper bound before counting/allocating, including duplicate entries.
    if (static_cast<long long>(cells)+2LL*faces > std::numeric_limits<int>::max())
        throw std::runtime_error("matrix nonzero count exceeds benchmark int32 range");
    AMGBenchmarkCSR a;
    a.rows.assign(static_cast<std::size_t>(cells)+1, 1);
    a.rows[0] = 0;
    for (int face = 0; face < faces; ++face)
    {
        const int l = left[face], r = right[face];
        if (l < 0 || l >= cells || r < 0)
            throw std::runtime_error("invalid finite-volume face addressing");
        if (r < cells) { ++a.rows[l+1]; ++a.rows[r+1]; }
    }
    for (int row = 0; row < cells; ++row) a.rows[row+1] += a.rows[row];
    a.columns.resize(a.rows.back());
    a.values.resize(a.rows.back());
    std::vector<int> next = a.rows;
    const auto put = [&](int row, int column, double value)
    {
        if (!std::isfinite(value)) throw std::runtime_error("nonfinite matrix coefficient");
        const int pos = next[row]++;
        a.columns[pos] = column;
        a.values[pos] = value;
    };
    for (int row = 0; row < cells; ++row) put(row, row, diagonal[row]);
    for (int face = 0; face < faces; ++face)
        if (right[face] < cells)
        {
            put(left[face], right[face], upper[face]);
            put(right[face], left[face], lower[face]);
        }
    const std::vector<int> oldRows = a.rows;
    std::vector<std::pair<int, double> > entries;
    int write = 0;
    for (int row = 0; row < cells; ++row)
    {
        entries.clear();
        for (int k = oldRows[row]; k < oldRows[row+1]; ++k)
            entries.push_back(std::make_pair(a.columns[k], a.values[k]));
        std::sort(entries.begin(), entries.end());
        a.rows[row] = write;
        for (const auto& item : entries)
            if (write > a.rows[row] && a.columns[write-1] == item.first)
            {
                a.values[write-1] += item.second;
                if (!std::isfinite(a.values[write-1]))
                    throw std::runtime_error("nonfinite merged matrix coefficient");
            }
            else
            {
                a.columns[write] = item.first;
                a.values[write++] = item.second;
            }
    }
    a.rows[cells] = write;
    a.columns.resize(write);
    a.values.resize(write);
    return a;
}
}
