/* SPDX-License-Identifier: AGPL-3.0-only */
#ifndef SOLVED_TECH_PANGULU_ADAPTER_H
#define SOLVED_TECH_PANGULU_ADAPTER_H
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif
int case_factor(int n, int64_t nnz, const int *rows, const int *cols,
                const double *real, const double *imag);
int case_solve(int n, const double *br, const double *bi, double *xr, double *xi);
int case_cleanup(void);
#ifdef __cplusplus
}
#endif
#endif
