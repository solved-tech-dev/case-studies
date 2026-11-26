#include <linglongamg.h>
#include <cuda_runtime.h>
#include <dlfcn.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <stdexcept>
#include <string>

template <class Function>
Function symbol(void* library, const char* name)
{
    const auto address = dlsym(library, name);
    if (!address) {
        throw std::runtime_error(name);
    }
    return reinterpret_cast<Function>(address);
}

void check(cudaError_t value)
{
    if (value != cudaSuccess) {
        throw std::runtime_error(cudaGetErrorString(value));
    }
}

void check(linglongamg_status_t value, const char* phase)
{
    if (value != LINGLONGAMG_STATUS_SUCCESS) {
        std::fprintf(stderr, "LinglongAMG %s returned %d\n", phase, int(value));
        throw std::runtime_error(phase);
    }
}

bool procMapsContains(const char* needle)
{
    std::ifstream maps("/proc/self/maps");
    std::string line;
    while (std::getline(maps, line)) {
        if (line.find(needle) != std::string::npos) {
            return true;
        }
    }
    return false;
}

template <class T, size_t N>
struct Device
{
    T* ptr{};
    explicit Device(const T (&host)[N])
    {
        check(cudaMalloc(&ptr, sizeof(host)));
        check(cudaMemcpy(ptr, host, sizeof(host), cudaMemcpyHostToDevice));
    }
    ~Device() { cudaFree(ptr); }
};

using ApplyCycleFn = linglongamg_status_t (*)(linglongamg_handle_t, const linglongamg_vec_t*, linglongamg_vec_t*);

int main(int argc, char** argv) try {
    if (argc != 3) {
        std::fprintf(stderr, "Usage: %s <liblinglongamg.so> <linglong-vcycle.json>\n", argv[0]);
        return 2;
    }
    void* library = dlopen(argv[1], RTLD_NOW | RTLD_LOCAL);
    if (!library) {
        throw std::runtime_error(dlerror());
    }
    if (procMapsContains("libXianzhiKrylov")) {
        throw std::runtime_error("Isolation violation: libXianzhiKrylov is loaded in LinglongAMG-only smoke test");
    }

    linglongamg_handle_t handle{};
    check(symbol<decltype(&linglongamg_create_from_config)>(library, "linglongamg_create_from_config")(
              &handle, argv[2]),
          "create");

    const int rows[]{0, 2, 5, 8, 10}, cols[]{0, 1, 0, 1, 2, 1, 2, 3, 2, 3};
    const double data[]{4, -1, -1, 4, -1, -1, 4, -1, -1, 4}, rhs[]{3, 2, 2, 3}, zeros[4]{};
    Device<int, 5> offsets(rows);
    Device<int, 10> columns(cols);
    Device<double, 10> values(data);
    Device<double, 4> b(rhs), x(zeros);
    linglongamg_csr_t* matrix{};
    check(symbol<decltype(&linglongamg_read_csr)>(library, "linglongamg_read_csr")(
              &matrix, 4, 4, 10, offsets.ptr, columns.ptr, values.ptr),
          "CSR");
    check(symbol<decltype(&linglongamg_setup)>(library, "linglongamg_setup")(handle, matrix), "setup");

    const linglongamg_vec_t input{4, b.ptr};
    linglongamg_vec_t output{4, x.ptr};
    if (auto cycle = reinterpret_cast<ApplyCycleFn>(dlsym(library, "linglongamg_apply_cycle"))) {
        check(cycle(handle, &input, &output), "apply cycle");
    } else {
        check(symbol<decltype(&linglongamg_solve)>(library, "linglongamg_solve")(handle, &input, &output),
              "solve single step");
        linglongamg_solve_info_t info{};
        check(symbol<decltype(&linglongamg_get_solve_info)>(library, "linglongamg_get_solve_info")(handle, &info),
              "solve info");
        if (info.iterations != 1) {
            throw std::runtime_error("expected one AMG step");
        }
    }

    double computed[4]{};
    check(cudaMemcpy(computed, x.ptr, sizeof computed, cudaMemcpyDeviceToHost));
    double squared = 0;
    for (int i = 0; i < 4; ++i) {
        double residual = rhs[i];
        for (int k = rows[i]; k < rows[i + 1]; ++k) {
            residual -= data[k] * computed[cols[k]];
        }
        squared += residual * residual;
    }
    const double relative = std::sqrt(squared / 26.0);
    if (!std::isfinite(relative) || relative >= 1.0) {
        throw std::runtime_error("AMG cycle did not reduce the residual");
    }
    check(symbol<decltype(&linglongamg_destroy)>(library, "linglongamg_destroy")(handle), "destroy");
    symbol<decltype(&linglongamg_csr_destroy)>(library, "linglongamg_csr_destroy")(matrix);
    dlclose(library);
    std::printf("LINGLONG_CYCLE_PASS relative=%.17g isolated_no_xianzhikrylov=1\n", relative);
    return 0;
} catch (const std::exception& error) {
    std::fprintf(stderr, "%s\n", error.what());
    return 1;
}
