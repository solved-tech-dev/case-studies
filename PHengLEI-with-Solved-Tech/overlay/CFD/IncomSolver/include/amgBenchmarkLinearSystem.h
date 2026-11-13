#pragma once

#include <vector>

namespace PHSPACE
{
// Host matrix representation shared with the isolated benchmark worker.
struct AMGBenchmarkCSR
{
    std::vector<int> rows, columns;
    std::vector<double> values;
};

// Physical boundary contributions are already folded into diagonal/RHS by CFD.
// Coupled/processor boundaries must be rejected by the caller before assembly.
AMGBenchmarkCSR AssembleAMGBenchmarkCSR(int cells, int faces, const int* left,
    const int* right, const double* diagonal, const double* upper,
    const double* lower);
}
