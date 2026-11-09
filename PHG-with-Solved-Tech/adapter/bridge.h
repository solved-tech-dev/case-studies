/* SPDX-License-Identifier: LGPL-2.1-or-later */
#ifndef PHG_PANGULU_BRIDGE_H
#define PHG_PANGULU_BRIDGE_H
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif
/* CSR is zero-based, full unsymmetric storage. Inputs remain caller-owned.
 * Return zero on success, nonzero on failure. No implicit CPU fallback. */
int phg_pangulu_solve_csr(int64_t n, int64_t nnz, const int64_t *rp,
    const int32_t *ci, const double *a, const double *b, double *x);
const char *phg_pangulu_backend_identity(void);
#ifdef __cplusplus
}
#endif
#endif
