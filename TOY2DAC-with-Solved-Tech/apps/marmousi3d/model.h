// SPDX-License-Identifier: AGPL-3.0-only
// Copyright 2026 Solved Tech
#pragma once
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif
typedef struct { int *r, *c; double *ar, *ai; int64_t nnz; } Matrix;
int load_truth(const char *path, double *vp);
int build_matrix(double frequency, const double *vp, Matrix *out);
void release_matrix(Matrix *matrix);
#ifdef __cplusplus
}
#endif
