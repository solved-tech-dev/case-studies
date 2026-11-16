// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Solved Tech
#define _POSIX_C_SOURCE 200809L
#include "ngspice/ngspice.h"
#include "ngspice/cktdefs.h"
#include "ngspice/smpdefs.h"
#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

static unsigned int call_count;

static void fail(const char *message)
{
    fprintf(stderr, "CAPTURE_ERROR %s errno=%d\n", message, errno);
    exit(72);
}

static void save(const char *directory, const char *suffix, const void *values,
                 size_t count, size_t element_size)
{
    char path[4096];
    int length = snprintf(path, sizeof(path), "%s/call-%04u.%s", directory,
                          call_count, suffix);
    if (length < 0 || (size_t)length >= sizeof(path))
        fail("output path");
    FILE *file = fopen(path, "wbx");
    if (!file)
        fail("create original system");
    size_t written = fwrite(values, element_size, count, file);
    if (fclose(file) || written != count)
        fail("write original system");
}

void SolvedTechCaptureSolve(CKTcircuit *ckt, const char *directory)
{
    KLUmatrix *matrix = ckt->CKTmatrix->SMPkluMatrix;
    if (!ckt->CKTmatrix->CKTkluMODE || !matrix || !matrix->KLUmatrixN || matrix->KLUmatrixIsComplex)
        fail("real KLU matrix and capture directory required");
    size_t n = matrix->KLUmatrixN, nz = matrix->KLUmatrixNZ;
    if (!nz || matrix->KLUmatrixAp[0] != 0 || (size_t)matrix->KLUmatrixAp[n] != nz)
        fail("invalid assembled CSC");

    static int64_t *previous_ptr;
    static int32_t *previous_idx;
    static double *previous_a;
    static size_t previous_n, previous_nz;
    static unsigned int pattern_source, coefficient_source;
    const int same_dimensions = n == previous_n && nz == previous_nz;

    int64_t *ptr = calloc(n + 1, sizeof(*ptr));
    int64_t *next = malloc((n + 1) * sizeof(*next));
    int32_t *idx = malloc(nz * sizeof(*idx));
    double *a = malloc(nz * sizeof(*a));
    double *b = malloc(n * sizeof(*b));
    double *x = malloc(n * sizeof(*x));
    if (!ptr || !next || !idx || !a || !b || !x)
        fail("capture allocation");

    for (size_t k = 0; k < nz; ++k) {
        int row = matrix->KLUmatrixAi[k];
        if (row < 0 || (size_t)row >= n || !isfinite(matrix->KLUmatrixAx[k]))
            fail("invalid matrix entry");
        ++ptr[row + 1];
    }
    for (size_t i = 1; i <= n; ++i)
        ptr[i] += ptr[i - 1];
    memcpy(next, ptr, (n + 1) * sizeof(*ptr));
    for (size_t col = 0; col < n; ++col) {
        int lo = matrix->KLUmatrixAp[col], hi = matrix->KLUmatrixAp[col + 1];
        if (lo < 0 || hi < lo || (size_t)hi > nz)
            fail("invalid CSC column");
        for (int k = lo; k < hi; ++k) {
            int64_t pos = next[matrix->KLUmatrixAi[k]]++;
            idx[pos] = (int32_t)col;
            a[pos] = matrix->KLUmatrixAx[k];
        }
        unsigned int old = matrix->KLUmatrixNodeCollapsingNewToOld[col + 1];
        b[col] = old ? ckt->CKTrhs[old] : 0.0;
    }
    free(next);

    ++call_count;
    const int pattern_equal = same_dimensions &&
        !memcmp(ptr, previous_ptr, (n + 1) * sizeof(*ptr)) &&
        !memcmp(idx, previous_idx, nz * sizeof(*idx));
    const int coefficients_equal = pattern_equal &&
        !memcmp(a, previous_a, nz * sizeof(*a));
    if (!pattern_equal) {
        pattern_source = call_count;
        save(directory, "rp.i64", ptr, n + 1, sizeof(*ptr));
        save(directory, "ci.i32", idx, nz, sizeof(*idx));
    }
    if (!coefficients_equal) {
        coefficient_source = call_count;
        save(directory, "a.f64", a, nz, sizeof(*a));
    }
    free(previous_ptr);
    free(previous_idx);
    free(previous_a);
    previous_ptr = ptr;
    previous_idx = idx;
    previous_a = a;
    previous_n = n;
    previous_nz = nz;

    SMPsolve(ckt->CKTmatrix, ckt->CKTrhs, ckt->CKTrhsSpare);
    for (size_t i = 0; i < n; ++i) {
        unsigned int old = matrix->KLUmatrixNodeCollapsingNewToOld[i + 1];
        x[i] = old ? ckt->CKTrhs[old] : 0.0;
        if (!isfinite(x[i]) || !isfinite(b[i]))
            fail("nonfinite KLU reference vector");
    }
    long double sqerr = 0, sqrhs = 0;
    for (size_t i = 0; i < n; ++i) {
        long double r = -b[i];
        for (int64_t k = ptr[i]; k < ptr[i + 1]; ++k)
            r += (long double)a[k] * x[idx[k]];
        sqerr += r * r;
        sqrhs += (long double)b[i] * b[i];
    }
    double residual = (double)(sqrtl(sqerr) / (sqrhs > 0 ? sqrtl(sqrhs) : 1));
    if (!isfinite(residual))
        fail("nonfinite KLU reference residual");
    save(directory, "b.f64", b, n, sizeof(*b));
    save(directory, "x.f64", x, n, sizeof(*x));
    free(b);
    free(x);
    printf("NGSPICE_CAPTURE call=%u n=%zu nnz=%zu graph=%u values=%u residual=%.17g time=%.17g\n",
           call_count, n, nz, pattern_source, coefficient_source, residual, ckt->CKTtime);
    fflush(stdout);
}
