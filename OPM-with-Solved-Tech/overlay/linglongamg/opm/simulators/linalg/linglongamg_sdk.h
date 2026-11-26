#ifndef LINGLONGAMG_H
#define LINGLONGAMG_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#if defined(LINGLONGAMG_USE_FP32)
typedef float linglongamg_real_t;
#else
typedef double linglongamg_real_t;
#endif

typedef enum {
    LINGLONGAMG_STATUS_SUCCESS = 0,
    LINGLONGAMG_STATUS_INVALID_VALUE = 1,
    LINGLONGAMG_STATUS_ALLOC_FAILED = 2,
    LINGLONGAMG_STATUS_CUDA_ERROR = 3,
    LINGLONGAMG_STATUS_NOT_IMPLEMENTED = 4,
    LINGLONGAMG_STATUS_RUNTIME_ERROR = 5,
    LINGLONGAMG_STATUS_LICENSE_ERROR = 6
} linglongamg_status_t;

/* Vector values passed to solve reside on the currently selected CUDA device. */
typedef struct {
    int n;
    linglongamg_real_t *val;
} linglongamg_vec_t;

typedef struct {
    int converged;
    int iterations;
    linglongamg_real_t initial_residual;
    linglongamg_real_t final_residual;
    linglongamg_real_t relative_residual;
    double setup_ms;
    double solve_ms;
} linglongamg_solve_info_t;

struct linglongamg_solver_st;
typedef struct linglongamg_solver_st *linglongamg_handle_t;
struct linglongamg_csr;
typedef struct linglongamg_csr linglongamg_csr_t;

const char *linglongamg_status_string(linglongamg_status_t status);

/* Load all Setup/Solve parameters from config_file, validate the GPU-bound
 * license, and create an opaque solver. */
linglongamg_status_t linglongamg_create_from_config(
    linglongamg_handle_t *solver,
    const char *config_file);
/* Create an opaque CSR descriptor for arrays on the currently selected CUDA
 * device. row_ptr and col_idx use 32-bit indices. The arrays are not copied
 * or owned by the descriptor and must remain valid while the solver uses A. */
linglongamg_status_t linglongamg_read_csr(
    linglongamg_csr_t **A,
    int nrows,
    int ncols,
    int nnz,
    const int *device_row_ptr,
    const int *device_col_idx,
    const linglongamg_real_t *device_values);
void linglongamg_csr_destroy(linglongamg_csr_t *A);
linglongamg_status_t linglongamg_setup(
    linglongamg_handle_t solver,
    const linglongamg_csr_t *A);
linglongamg_status_t linglongamg_solve(
    linglongamg_handle_t solver,
    const linglongamg_vec_t *b,
    linglongamg_vec_t *x);
linglongamg_status_t linglongamg_get_solve_info(
    linglongamg_handle_t solver,
    linglongamg_solve_info_t *info);
linglongamg_status_t linglongamg_reset(linglongamg_handle_t solver);
linglongamg_status_t linglongamg_destroy(linglongamg_handle_t solver);

/* These two calls support the activation workflow.  machine_key intentionally
 * works before activation; create/setup/solve validate the signed license
 * inside the precompiled library. */
int linglongamg_machine_key(char *key, size_t key_size,
                           char *error, size_t error_size);
/* On success, error receives the signed licensee; on failure, it receives a
 * diagnostic message. The return value is zero only for a valid license. */
int linglongamg_license_check(char *error, size_t error_size);

#ifdef __cplusplus
}
#endif

#endif
