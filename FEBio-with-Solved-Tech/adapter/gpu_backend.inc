// FEBio native solver interface. Each FEBio solver instance owns one factorization.
#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <dlfcn.h>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <omp.h>
#include <cuda_runtime_api.h>
#include <pangulu.h>
static_assert(PANGULU_ABI_VERSION == 2);
#include <FECore/LinearSolver.h>
#include <FECore/CompactUnSymmMatrix.h>
#include <FECore/FECoreKernel.h>
#include <FECore/sdk.h>
#include <FECore/version.h>
#include <unistd.h>
namespace {
using Clock = std::chrono::steady_clock;
using Moment = Clock::time_point;
thread_local int call = 0;
void require(bool ok, const char *message) {
    if (!ok) throw std::runtime_error(message);
}
void status(int value, const char *phase) {
    if (value) throw std::runtime_error(std::string(phase) + " status=" + std::to_string(value));
}
void event(const char *scope, const char *phase, Moment start, int code = 0) {
    double seconds = std::chrono::duration<double>(Clock::now() - start).count();
    std::printf("FEBIO_TIME call=%d scope=%s phase=%s seconds=%.9f status=%d\n", call, scope, phase, seconds, code);
    std::fflush(stdout);
}
void provider(const void *symbol) {
    Dl_info info{};
    require(dladdr(symbol, &info) && info.dli_fname, "native provider lookup");
    std::printf("FEBIO_PROVIDER %s\n", info.dli_fname);
}
struct Sparse {
    std::vector<int64_t> ptr;
    std::vector<int32_t> col;
    std::vector<double> values;
    std::vector<int32_t> from_csr;
    void refresh(const double *a) {
        auto start = Clock::now();
        for (size_t i=0; i<values.size(); ++i) {
            values[i] = a[from_csr[i]];
            require(std::isfinite(values[i]), "nonfinite updated coefficient");
        }
        event("service_part", "update_all_values", start);
    }
};
Sparse csr_to_csc(int n, int64_t nz, const int64_t *rp, const int32_t *ci, const double *a) {
    auto start = Clock::now();
    require(rp[0] == 0 && rp[n] == nz, "invalid CSR endpoints");
    Sparse out;
    out.ptr.assign(n + 1ULL, 0);
    out.col.resize(nz);
    out.values.resize(nz);
    out.from_csr.resize(nz);
    for (int row = 0; row < n; ++row) {
        require(rp[row] >= 0 && rp[row] <= rp[row+1] && rp[row+1] <= nz, "invalid CSR row pointer");
        for (auto k = rp[row]; k < rp[row+1]; ++k) {
            require(ci[k] >= 0 && ci[k] < n && std::isfinite(a[k]), "invalid original coefficient");
            ++out.ptr[ci[k]+1];
        }
    }
    for (int i = 0; i < n; ++i) out.ptr[i+1] += out.ptr[i];
    auto next = out.ptr;
    for (int row = 0; row < n; ++row)
        for (auto k = rp[row]; k < rp[row+1]; ++k) {
            auto pos = next[ci[k]]++;
            out.col[pos] = row;
            out.values[pos] = a[k];
            out.from_csr[pos] = int32_t(k);
        }
    for (int i = 0; i < n; ++i)
        for (auto k = out.ptr[i]+1; k < out.ptr[i+1]; ++k)
            require(out.col[k] > out.col[k-1], "duplicate original coefficient");
    event("service_part", "csr_to_csc", start);
    return out;
}
void device() {
    int ordinal = 0, runtime = 0;
    cudaDeviceProp properties{};
    status(cudaGetDevice(&ordinal), "CUDA device");
    status(cudaGetDeviceProperties(&properties, ordinal), "CUDA properties");
    status(cudaRuntimeGetVersion(&runtime), "CUDA version");
    std::printf("FEBIO_GPU name=%s; memory_bytes=%llu runtime=%d\n", properties.name,
                (unsigned long long)properties.totalGlobalMem, runtime);
}
struct DeviceBuffer {
    void *data = nullptr;
    explicit DeviceBuffer(size_t bytes) { status(cudaMalloc(&data, bytes), "cudaMalloc"); }
    ~DeviceBuffer() { if (data) cudaFree(data); }
    void upload(const void *values, size_t bytes) { status(cudaMemcpy(data, values, bytes, cudaMemcpyHostToDevice), "H2D"); }
    void download(void *values, size_t bytes) { status(cudaMemcpy(values, data, bytes, cudaMemcpyDeviceToHost), "D2H"); }
};
template<class Execute> void gpu_phase(const char *phase, Execute execute) {
    status(cudaDeviceSynchronize(), "pre phase synchronization");
    auto start = Clock::now();
    int code = execute();
    int sync = cudaDeviceSynchronize();
    event("native", phase, start, code ? code : sync);
    status(code, phase);
    status(sync, "post phase synchronization");
}
struct Context {
    Sparse input;
    DeviceBuffer rhs, sol;
    pangulu_handle_t handle = nullptr;
    pangulu_matrix_t matrix = nullptr;
    pangulu_vector_t b = nullptr, x = nullptr;
    Context(Sparse a, size_t bytes) : input(std::move(a)), rhs(bytes), sol(bytes) {}
    ~Context() {
        auto start = Clock::now();
        if (x) pangulu_vector_destroy(x);
        if (b) pangulu_vector_destroy(b);
        if (matrix) pangulu_matrix_destroy(matrix);
        if (handle) pangulu_destroy(handle);
        event("lifecycle", "release", start);
    }
};

std::unique_ptr<Context> create_context(int n, const std::vector<int64_t>& rp,
                                       const int* ci, const double* a) {
    device();
    auto nz = rp.back();
    provider(reinterpret_cast<void*>(pangulu_execute));
    auto state = std::make_unique<Context>(csr_to_csc(n,nz,rp.data(),ci,a), n*sizeof(double));
    auto &ctx = *state;
    status(pangulu_create(0, &ctx.handle), "PanguLU create");
    status(pangulu_config_set(ctx.handle, PANGULU_REORDERING_METHOD, L0), "PanguLU ordering");
    status(pangulu_config_set(ctx.handle, PANGULU_REUSE_LEVEL, L2), "PanguLU reuse");
    std::printf("FEBIO_CONFIG PanguLU version=26.46 abi=2 ordering=L0 reuse=L2\n");
    auto &matrix = ctx.input;
    status(pangulu_matrix_create(&ctx.matrix, PANGULU_MATRIX_LAYOUT_CSC, n, n, nz, 1,
                                  matrix.ptr.data(), matrix.col.data(), matrix.values.data(),
                                  PANGULU_R_64I, PANGULU_R_32I, PANGULU_R_64F,
                                  PANGULU_LU, PANGULU_MEMORY_HOST), "PanguLU matrix");
    status(pangulu_multiple_vector_create(&ctx.b, n, 1, n, ctx.rhs.data, PANGULU_R_64F, PANGULU_MEMORY_DEVICE), "PanguLU RHS");
    status(pangulu_multiple_vector_create(&ctx.x, n, 1, n, ctx.sol.data, PANGULU_R_64F, PANGULU_MEMORY_DEVICE), "PanguLU solution");
    gpu_phase("analysis", [&]{ return pangulu_execute(ctx.handle, PANGULU_PHASE_ANALYSIS, ctx.matrix, nullptr, nullptr); });
    return state;
}
void factor(Context& ctx, bool fresh, const double* values) {
    if (!fresh) ctx.input.refresh(values);
    auto phase = fresh ? PANGULU_PHASE_NUMERIC_FACTORISATION : PANGULU_PHASE_REFACTORISATION;
    gpu_phase(fresh ? "numeric" : "refactor_numeric", [&]{ return pangulu_execute(ctx.handle,phase,ctx.matrix,nullptr,nullptr); });
}
void backsolve(Context& ctx, int n, double* x, const double* b) {
    ctx.rhs.upload(b,n*sizeof(double));
    status(cudaMemset(ctx.sol.data,0,n*sizeof(double)), "reset solution");
    gpu_phase("solve",[&]{ return pangulu_execute(ctx.handle,PANGULU_PHASE_SOLVE,ctx.matrix,ctx.x,ctx.b); });
    ctx.sol.download(x,n*sizeof(double));
}
} // namespace
