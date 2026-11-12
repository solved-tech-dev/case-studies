#include "linglongSession.h"
#include "linglongamg.h"
#include <cuda_runtime_api.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <dlfcn.h>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <vector>

#ifndef PHENGLEI_LINGLONG_LIBRARY_PATH
#define PHENGLEI_LINGLONG_LIBRARY_PATH "liblinglongamg.so"
#endif

namespace PHSPACE
{
namespace
{
using Clock = std::chrono::steady_clock;
double Milliseconds(Clock::time_point start)
{
    return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}
void Cuda(cudaError_t status, const char* operation)
{
    if (status != cudaSuccess)
        throw std::runtime_error(std::string(operation) + ": " + cudaGetErrorString(status));
}
class DeviceScope
{
    int previous;
public:
    explicit DeviceScope(int device)
    {
        Cuda(cudaGetDevice(&previous), "cudaGetDevice");
        Cuda(cudaSetDevice(device), "cudaSetDevice");
    }
    ~DeviceScope() { cudaSetDevice(previous); }
};
template<class T> struct Buffer
{
    T* data = nullptr;
    std::size_t capacity = 0;
    void Ensure(std::size_t count)
    {
        if (count <= capacity) return;
        T* replacement = nullptr;
        Cuda(cudaMalloc(reinterpret_cast<void**>(&replacement), count * sizeof(T)), "cudaMalloc");
        if (data)
        {
            const cudaError_t status = cudaFree(data);
            if (status != cudaSuccess) { cudaFree(replacement); Cuda(status, "cudaFree"); }
        }
        data = replacement;
        capacity = count;
    }
    void Free()
    {
        if (data) Cuda(cudaFree(data), "cudaFree");
        data = nullptr;
        capacity = 0;
    }
    void Upload(const T* source, std::size_t count)
    {
        Cuda(cudaMemcpy(data, source, count * sizeof(T), cudaMemcpyHostToDevice), "cudaMemcpy H2D");
    }
};
template<class T> T Symbol(void* library, const char* name)
{
    dlerror();
    void* address = dlsym(library, name);
    const char* error = dlerror();
    if (error) throw std::runtime_error(std::string("LinglongAMG symbol ") + name + ": " + error);
    return reinterpret_cast<T>(address);
}
int Validate(const AMGBenchmarkCSR& a, const double* b, const double* x)
{
    if (a.rows.size() < 2 || a.rows.size()-1 > static_cast<std::size_t>(std::numeric_limits<int>::max())
        || a.values.empty() || a.values.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())
        || a.columns.size() != a.values.size() || !b || !x || a.rows.front() != 0
        || a.rows.back() != static_cast<int>(a.values.size()))
        throw std::invalid_argument("invalid FP64/int32 CSR matrix or vectors");
    const int n = static_cast<int>(a.rows.size()-1);
    for (int row = 0; row < n; ++row)
    {
        if (a.rows[row] < 0 || a.rows[row] > a.rows[row+1] || a.rows[row+1] > a.rows.back()
            || !std::isfinite(b[row]) || !std::isfinite(x[row]))
            throw std::invalid_argument("invalid CSR row offsets or nonfinite vectors");
        int previous = -1;
        for (int k = a.rows[row]; k < a.rows[row+1]; ++k)
        {
            if (a.columns[k] <= previous || a.columns[k] >= n || !std::isfinite(a.values[k]))
                throw std::invalid_argument("CSR columns must be sorted, unique and in range; values must be finite");
            previous = a.columns[k];
        }
    }
    return n;
}
double Residual(const AMGBenchmarkCSR& a, const double* b, const double* x)
{
    long double norm = 0;
    for (std::size_t row = 0; row+1 < a.rows.size(); ++row)
    {
        long double value = b[row];
        for (int k = a.rows[row]; k < a.rows[row+1]; ++k)
            value -= static_cast<long double>(a.values[k]) * x[a.columns[k]];
        norm = std::hypot(norm, value);
    }
    const double result = static_cast<double>(norm);
    if (!std::isfinite(result)) throw std::runtime_error("nonfinite independently computed residual");
    return result;
}
}

struct LinglongSession::Impl
{
    void* library = nullptr;
    linglongamg_handle_t solver = nullptr;
    linglongamg_csr_t* matrix = nullptr;
    decltype(&linglongamg_status_string) statusString = nullptr;
    decltype(&linglongamg_create_from_config) create = nullptr;
    decltype(&linglongamg_read_csr) read = nullptr;
    decltype(&linglongamg_csr_destroy) destroyMatrix = nullptr;
    decltype(&linglongamg_setup) setup = nullptr;
    decltype(&linglongamg_solve) solve = nullptr;
    decltype(&linglongamg_get_solve_info) query = nullptr;
    decltype(&linglongamg_reset) reset = nullptr;
    decltype(&linglongamg_destroy) destroy = nullptr;
    Buffer<int> rows, columns;
    Buffer<double> values, rhs, solution;
    int device = 0;
    bool hasSetup = false, failed = false;
    unsigned long long calls = 0;
    double createMs = 0;
    std::string description;

    void Check(linglongamg_status_t status, const char* operation)
    {
        if (status != LINGLONGAMG_STATUS_SUCCESS)
            throw std::runtime_error(std::string(operation) + ": " + statusString(status));
    }
    void Close()
    {
        if (!library) return;
        // The SDK may retain a handle when the original GPU cannot be verified.
        // In that case do not free arrays/descriptors or unload its library.
        if (solver || matrix || rows.data || columns.data || values.data || rhs.data || solution.data)
        {
            DeviceScope selected(device);
            if (solver) { Check(destroy(solver), "linglongamg_destroy"); solver = nullptr; }
            if (matrix) { destroyMatrix(matrix); matrix = nullptr; }
            rows.Free(); columns.Free(); values.Free(); rhs.Free(); solution.Free();
        }
        dlclose(library);
        library = nullptr;
    }
    ~Impl()
    {
        try { Close(); }
        catch (const std::exception& error)
        {
            std::fprintf(stderr, "LinglongAMG cleanup failed: %s\n", error.what());
        }
    }
};

LinglongSession::LinglongSession(const std::string& config, int device,
                                 const std::string& library) : impl(new Impl)
{
    static_assert(std::is_same<linglongamg_real_t, double>::value, "FP64 SDK required");
    if (config.empty() || device < 0) throw std::invalid_argument("config path and nonnegative device required");
    impl->device = device;
    DeviceScope selected(device);
    const char* overridePath = std::getenv("PHENGLEI_LINGLONG_LIBRARY");
    const std::string path = !library.empty() ? library
        : (overridePath && *overridePath ? overridePath : PHENGLEI_LINGLONG_LIBRARY_PATH);
    impl->library = dlopen(path.c_str(), RTLD_NOW | RTLD_LOCAL);
    if (!impl->library) throw std::runtime_error(std::string("load LinglongAMG: ") + dlerror());
#define LOAD(member, function) impl->member = Symbol<decltype(impl->member)>(impl->library, #function)
    LOAD(statusString, linglongamg_status_string);
    LOAD(create, linglongamg_create_from_config);
    LOAD(read, linglongamg_read_csr);
    LOAD(destroyMatrix, linglongamg_csr_destroy);
    LOAD(setup, linglongamg_setup);
    LOAD(solve, linglongamg_solve);
    LOAD(query, linglongamg_get_solve_info);
    LOAD(reset, linglongamg_reset);
    LOAD(destroy, linglongamg_destroy);
#undef LOAD
    cudaDeviceProp properties;
    Cuda(cudaGetDeviceProperties(&properties, device), "cudaGetDeviceProperties");
    impl->description = properties.name;
    const auto start = Clock::now();
    impl->Check(impl->create(&impl->solver, config.c_str()), "linglongamg_create_from_config");
    impl->createMs = Milliseconds(start);
    if (!impl->solver) throw std::runtime_error("SDK returned success without a solver handle");
}

LinglongSession::~LinglongSession() = default;
double LinglongSession::CreateMilliseconds() const { return impl->createMs; }
const std::string& LinglongSession::DeviceDescription() const { return impl->description; }
void LinglongSession::Close() { impl->Close(); }

LinglongSessionReport LinglongSession::Solve(const AMGBenchmarkCSR& a, const double* b, double* x)
{
    const auto totalStart = Clock::now();
    const int n = Validate(a, b, x);
    if (!impl->solver || impl->failed) throw std::runtime_error("LinglongAMG session is closed or failed; no automatic recreation");
    LinglongSessionReport report;
    report.initialL2 = Residual(a, b, x);
    long double rhsNorm = 0;
    for (int row = 0; row < n; ++row) rhsNorm = std::hypot(rhsNorm, static_cast<long double>(b[row]));
    report.rhsL2 = static_cast<double>(rhsNorm);
    if (!std::isfinite(report.rhsL2)) throw std::invalid_argument("RHS norm overflows FP64");
    try
    {
        DeviceScope selected(impl->device);
        auto start = Clock::now();
        // Reset setup state on the SAME handle before releasing old CSR storage.
        // No create, JSON reload, license_check or refresh_license in this loop.
        if (impl->hasSetup) impl->Check(impl->reset(impl->solver), "linglongamg_reset");
        impl->hasSetup = false;
        if (impl->matrix) { impl->destroyMatrix(impl->matrix); impl->matrix = nullptr; }
        report.resetMs = Milliseconds(start);
        start = Clock::now();
        impl->rows.Ensure(a.rows.size()); impl->columns.Ensure(a.columns.size());
        impl->values.Ensure(a.values.size()); impl->rhs.Ensure(n); impl->solution.Ensure(n);
        impl->rows.Upload(a.rows.data(), a.rows.size());
        impl->columns.Upload(a.columns.data(), a.columns.size());
        impl->values.Upload(a.values.data(), a.values.size());
        impl->rhs.Upload(b, n); impl->solution.Upload(x, n);
        impl->Check(impl->read(&impl->matrix, n, n, static_cast<int>(a.values.size()),
            impl->rows.data, impl->columns.data, impl->values.data), "linglongamg_read_csr");
        Cuda(cudaDeviceSynchronize(), "synchronize upload");
        report.uploadMs = Milliseconds(start);
        start = Clock::now();
        impl->Check(impl->setup(impl->solver, impl->matrix), "linglongamg_setup");
        impl->hasSetup = true;
        Cuda(cudaDeviceSynchronize(), "synchronize setup");
        report.setupMs = Milliseconds(start);
        linglongamg_vec_t rhs = {n, impl->rhs.data}, solution = {n, impl->solution.data};
        start = Clock::now();
        impl->Check(impl->solve(impl->solver, &rhs, &solution), "linglongamg_solve");
        Cuda(cudaDeviceSynchronize(), "synchronize solve");
        report.solveMs = Milliseconds(start);
        linglongamg_solve_info_t info = {};
        start = Clock::now();
        // Mandatory: query authorizes and COMMITS x. Never download before success.
        impl->Check(impl->query(impl->solver, &info), "linglongamg_get_solve_info");
        Cuda(cudaDeviceSynchronize(), "synchronize authorized commit");
        report.queryMs = Milliseconds(start);
        if (info.iterations < 0 || !std::isfinite(info.initial_residual)
            || !std::isfinite(info.final_residual) || !std::isfinite(info.relative_residual)
            || !std::isfinite(info.setup_ms) || !std::isfinite(info.solve_ms))
            throw std::runtime_error("nonfinite or invalid SDK solve metadata");
        start = Clock::now();
        std::vector<double> candidate(n);
        Cuda(cudaMemcpy(candidate.data(), impl->solution.data, n * sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy D2H");
        report.downloadMs = Milliseconds(start);
        start = Clock::now();
        for (double value : candidate)
            if (!std::isfinite(value)) throw std::runtime_error("nonfinite candidate solution");
        report.finalL2 = Residual(a, b, candidate.data());
        report.verifyMs = Milliseconds(start);
        report.iterations = info.iterations;
        report.converged = info.converged != 0;
        report.sdkRelativeResidual = info.relative_residual;
        report.sdkSetupMs = info.setup_ms; report.sdkSolveMs = info.solve_ms;
        report.call = ++impl->calls;
        // Publish only after authorization and all independent checks succeed.
        std::copy(candidate.begin(), candidate.end(), x);
        report.totalMs = Milliseconds(totalStart);
        return report;
    }
    catch (...) { impl->failed = true; throw; }
}
}
