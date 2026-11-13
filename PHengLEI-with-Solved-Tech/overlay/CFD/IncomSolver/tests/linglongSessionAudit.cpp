// Test-only forwarding shim. Every successful operation uses the REAL SDK.
// Fault injection can only deny a query; it cannot approve authorization.
#include "linglongamg.h"
#include <cstdio>
#include <cstdlib>
#include <dlfcn.h>

namespace
{
void* Real()
{
    static void* library = []() -> void*
    {
        const char* path = std::getenv("LINGLONG_AUDIT_REAL_LIBRARY");
        void* result = path ? dlopen(path, RTLD_NOW | RTLD_LOCAL) : nullptr;
        if (!result) { std::fprintf(stderr, "audit: real SDK load failed\n"); std::abort(); }
        return result;
    }();
    return library;
}
template<class T> T Function(const char* name)
{
    void* result = dlsym(Real(), name);
    if (!result) { std::fprintf(stderr, "audit: missing real symbol %s\n", name); std::abort(); }
    return reinterpret_cast<T>(result);
}
void Record(const char* event, const void* handle, int status)
{
    const char* path = std::getenv("LINGLONG_AUDIT_FILE");
    FILE* out = path ? std::fopen(path, "a") : nullptr;
    if (!out) { std::fprintf(stderr, "audit: LINGLONG_AUDIT_FILE required\n"); std::abort(); }
    std::fprintf(out, "%s %p %d\n", event, handle, status);
    std::fclose(out);
}
}

extern "C"
{
const char* linglongamg_status_string(linglongamg_status_t status)
{ return Function<decltype(&linglongamg_status_string)>("linglongamg_status_string")(status); }
linglongamg_status_t linglongamg_create_from_config(linglongamg_handle_t* solver, const char* config)
{
    const auto status = Function<decltype(&linglongamg_create_from_config)>("linglongamg_create_from_config")(solver, config);
    Record("create", *solver, status); return status;
}
linglongamg_status_t linglongamg_read_csr(linglongamg_csr_t** a, int n, int m, int nnz,
    const int* rows, const int* columns, const linglongamg_real_t* values)
{
    const auto status = Function<decltype(&linglongamg_read_csr)>("linglongamg_read_csr")(a, n, m, nnz, rows, columns, values);
    Record("read_csr", *a, status); return status;
}
void linglongamg_csr_destroy(linglongamg_csr_t* a)
{
    Function<decltype(&linglongamg_csr_destroy)>("linglongamg_csr_destroy")(a);
    Record("csr_destroy", a, 0);
}
linglongamg_status_t linglongamg_setup(linglongamg_handle_t solver, const linglongamg_csr_t* a)
{
    const auto status = Function<decltype(&linglongamg_setup)>("linglongamg_setup")(solver, a);
    Record("setup", solver, status); return status;
}
linglongamg_status_t linglongamg_solve(linglongamg_handle_t solver, const linglongamg_vec_t* b, linglongamg_vec_t* x)
{
    const auto status = Function<decltype(&linglongamg_solve)>("linglongamg_solve")(solver, b, x);
    Record("solve", solver, status); return status;
}
linglongamg_status_t linglongamg_get_solve_info(linglongamg_handle_t solver, linglongamg_solve_info_t* info)
{
    const auto status = std::getenv("LINGLONG_AUDIT_DENY_QUERY")
        ? LINGLONGAMG_STATUS_LICENSE_ERROR
        : Function<decltype(&linglongamg_get_solve_info)>("linglongamg_get_solve_info")(solver, info);
    Record("query", solver, status); return status;
}
linglongamg_status_t linglongamg_reset(linglongamg_handle_t solver)
{
    const auto status = Function<decltype(&linglongamg_reset)>("linglongamg_reset")(solver);
    Record("reset", solver, status); return status;
}
linglongamg_status_t linglongamg_destroy(linglongamg_handle_t solver)
{
    const auto status = Function<decltype(&linglongamg_destroy)>("linglongamg_destroy")(solver);
    Record("destroy", solver, status); return status;
}
int linglongamg_license_check(char* error, size_t size)
{
    const int status = Function<decltype(&linglongamg_license_check)>("linglongamg_license_check")(error, size);
    Record("license_check", nullptr, status); return status;
}
}
