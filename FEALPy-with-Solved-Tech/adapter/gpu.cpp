// SPDX-License-Identifier: AGPL-3.0-only
#include "gpu_backend.inc"

// All arrays belong to the caller. Indices are zero based. Values, RHS and
// solution use the scalar type selected when this shared library was built.
extern "C" int solve(int n, const int64_t* row_ptr, const int32_t* columns,
                     const Scalar* values, const Scalar* rhs, Scalar* solution) {
    try {
        if (n <= 0 || !row_ptr || !columns || !values || !rhs || !solution)
            throw std::invalid_argument("invalid solve arguments");
        solve_system(n, row_ptr, columns, values, rhs, solution);
        return 0;
    } catch (const std::exception& error) {
        std::fprintf(stderr, "FEALPy/PanguLU: %s\n", error.what());
        return 1;
    }
}
