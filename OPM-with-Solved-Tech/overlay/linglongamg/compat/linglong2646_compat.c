#include "linglongamg.h"
#include <dlfcn.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>

static pthread_once_t once = PTHREAD_ONCE_INIT;
static void *sdk;
static int ready;
static struct {
    __typeof__(&linglongamg_create_from_config) create;
    __typeof__(&linglongamg_read_csr) read;
    __typeof__(&linglongamg_csr_destroy) csr_destroy;
    __typeof__(&linglongamg_setup) setup;
    __typeof__(&linglongamg_solve) solve;
    __typeof__(&linglongamg_get_solve_info) info;
    __typeof__(&linglongamg_reset) reset;
    __typeof__(&linglongamg_destroy) destroy;
} api;

static void load_sdk(void)
{
    const char *path = getenv("LINGLONGAMG_2646_LIBRARY");
    if (!path) { fputs("LinglongAMG 26.46 library path missing\n", stderr); return; }
    sdk = dlopen(path, RTLD_NOW | RTLD_LOCAL);
    if (!sdk) { fprintf(stderr, "LinglongAMG 26.46 load failed: %s\n", dlerror()); return; }
#define LOAD(field, name) do { \
    api.field = (__typeof__(api.field)) dlsym(sdk, name); \
    if (!api.field) { fprintf(stderr, "Missing SDK symbol: %s\n", name); return; } \
} while (0)
    LOAD(create, "linglongamg_create_from_config");
    LOAD(read, "linglongamg_read_csr");
    LOAD(csr_destroy, "linglongamg_csr_destroy");
    LOAD(setup, "linglongamg_setup");
    LOAD(solve, "linglongamg_solve");
    LOAD(info, "linglongamg_get_solve_info");
    LOAD(reset, "linglongamg_reset");
    LOAD(destroy, "linglongamg_destroy");
#undef LOAD
    fprintf(stderr, "[OPM LinglongAMG 26.46 compatibility] SDK=%s; solve -> get_solve_info -> consume x\n", path);
    ready = 1;
}

static int initialized(void) { pthread_once(&once, load_sdk); return ready; }

linglongamg_status_t linglongamg_create_from_config(linglongamg_handle_t *s, const char *c)
{ return initialized() ? api.create(s,c) : LINGLONGAMG_STATUS_RUNTIME_ERROR; }
linglongamg_status_t linglongamg_read_csr(linglongamg_csr_t **a, int nr, int nc, int nnz,
    const int *r, const int *c, const linglongamg_real_t *v)
{ return initialized() ? api.read(a,nr,nc,nnz,r,c,v) : LINGLONGAMG_STATUS_RUNTIME_ERROR; }
void linglongamg_csr_destroy(linglongamg_csr_t *a)
{ if (initialized()) api.csr_destroy(a); }
linglongamg_status_t linglongamg_setup(linglongamg_handle_t s, const linglongamg_csr_t *a)
{ return initialized() ? api.setup(s,a) : LINGLONGAMG_STATUS_RUNTIME_ERROR; }
linglongamg_status_t linglongamg_solve(linglongamg_handle_t s, const linglongamg_vec_t *b, linglongamg_vec_t *x)
{
    if (!initialized()) return LINGLONGAMG_STATUS_RUNTIME_ERROR;
    linglongamg_status_t status = api.solve(s,b,x);
    if (status != LINGLONGAMG_STATUS_SUCCESS) return status;
    linglongamg_solve_info_t info = {0};
    return api.info(s,&info);
}
linglongamg_status_t linglongamg_get_solve_info(linglongamg_handle_t s, linglongamg_solve_info_t *i)
{ return initialized() ? api.info(s,i) : LINGLONGAMG_STATUS_RUNTIME_ERROR; }
linglongamg_status_t linglongamg_reset(linglongamg_handle_t s)
{ return initialized() ? api.reset(s) : LINGLONGAMG_STATUS_RUNTIME_ERROR; }
linglongamg_status_t linglongamg_destroy(linglongamg_handle_t s)
{ return initialized() ? api.destroy(s) : LINGLONGAMG_STATUS_RUNTIME_ERROR; }
