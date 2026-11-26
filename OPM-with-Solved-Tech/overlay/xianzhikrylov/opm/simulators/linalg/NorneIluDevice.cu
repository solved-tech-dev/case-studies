#include "NorneIluDevice.hpp"
#include <XianzhiKrylov_CprComponents.h>
#include <XianzhiKrylov_MethodOptions.h>
#include <cuda_runtime.h>

#include <chrono>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <string>

namespace Opm::detail {
namespace {

using Clock = std::chrono::steady_clock;

double elapsed(Clock::time_point start)
{
    return std::chrono::duration<double>(Clock::now() - start).count();
}

void checked(cudaError_t code)
{
    if (code != cudaSuccess) {
        throw std::runtime_error(cudaGetErrorString(code));
    }
}

void ok(XianzhiKrylov_Status code, const char* action)
{
    if (code != XIANZHI_KRYLOV_SUCCESS) {
        throw std::runtime_error(std::string(action) + ": XianzhiKrylov status " + std::to_string(code));
    }
}

template <class T>
struct Buffer
{
    T* ptr = nullptr;
    std::size_t count = 0;

    Buffer() = default;
    Buffer(const Buffer&) = delete;
    ~Buffer() { cudaFree(ptr); }

    void allocate(std::size_t n)
    {
        count = n;
        checked(cudaMalloc(&ptr, n * sizeof(T)));
    }

    void put(std::span<const T> a)
    {
        if (a.size() != count) {
            throw std::logic_error("ILU device data changed size");
        }
        checked(cudaMemcpy(ptr, a.data(), count * sizeof(T), cudaMemcpyHostToDevice));
    }

    void get(std::span<T> a) const
    {
        if (a.size() != count) {
            throw std::logic_error("ILU result changed size");
        }
        checked(cudaMemcpy(a.data(), ptr, count * sizeof(T), cudaMemcpyDeviceToHost));
    }
};

template <class T>
void putNew(Buffer<T>& b, std::span<const T> s)
{
    b.allocate(s.size());
    b.put(s);
}

__global__ void fineSpmv(int n,
                         const uint64_t* row,
                         const uint64_t* col,
                         const double* a,
                         const double* x,
                         double* y)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) {
        return;
    }
    double v[3]{};
    for (auto k = row[i]; k < row[i + 1]; ++k) {
        const auto j = 3 * col[k];
        for (int r = 0; r < 3; ++r) {
            for (int c = 0; c < 3; ++c) {
                v[r] += a[9 * k + 3 * r + c] * x[j + c];
            }
        }
    }
    for (int r = 0; r < 3; ++r) {
        y[3 * i + r] = v[r];
    }
}

constexpr int threads = 256;
int blocks(int n)
{
    return (n + threads - 1) / threads;
}

template <class Run>
int callback(const Run& work) noexcept
{
    try {
        work();
        return 0;
    } catch (const std::exception& error) {
        std::fprintf(stderr, "NORNE_ILU_CALLBACK_FAILURE %s\n", error.what());
        return 1;
    }
}

} // namespace

struct IluDevice
{
    int n = 0;
    int device = 0;
    cudaEvent_t opStart{}, opFine{}, opEnd{};
    Buffer<uint64_t> row, col;
    Buffer<double> a, rhs, x;
    XkCprwBlockIlu* ilu = nullptr;
    XianzhiKrylov_Solver* solver = nullptr;
    Xianzhi_DeviceLinearSystem linearSystem{};
    void* well = nullptr;
    int (*wellApply)(void*, const double*, double*) = nullptr;
    double fineSpmvSeconds = 0.0;
    double wellSeconds = 0.0;
    double iluApplySeconds = 0.0;
    unsigned matvecs = 0;

    ~IluDevice()
    {
        if (solver) {
            xk_internal_cprw_krylov_destroy(solver);
        }
        if (ilu) {
            xk_internal_cprw_ilu_destroy(ilu);
        }
        if (opStart) {
            cudaEventDestroy(opStart);
        }
        if (opFine) {
            cudaEventDestroy(opFine);
        }
        if (opEnd) {
            cudaEventDestroy(opEnd);
        }
    }

    Xianzhi_DeviceLinearSystem system(Buffer<double>& b, Buffer<double>& solution)
    {
        Xianzhi_DeviceLinearSystem s{};
        s.abi_version = XIANZHI_DEVICE_LINEAR_SYSTEM_ABI_VERSION;
        s.struct_size = sizeof(s);
        s.matrix_format = XIANZHI_MATRIX_FORMAT_OPERATOR;
        s.scalar_type = XIANZHI_SCALAR_TYPE_FP64;
        s.block_rows = s.block_columns = s.vector_length = b.count;
        s.block_size = s.rhs_count = 1;
        s.device_id = device;
        s.device_rhs = b.ptr;
        s.device_solution = solution.ptr;
        s.topology_epoch = s.values_epoch = s.rhs_epoch = 1;
        s.lifetime_token = this;
        return s;
    }

    void fullOperator(const double* in, double* out)
    {
        checked(cudaEventRecord(opStart));
        fineSpmv<<<blocks(n), threads>>>(n, row.ptr, col.ptr, a.ptr, in, out);
        checked(cudaGetLastError());
        checked(cudaEventRecord(opFine));
        if (!wellApply || wellApply(well, in, out)) {
            throw std::runtime_error("OPM independent well operator failed");
        }
        checked(cudaEventRecord(opEnd));
        checked(cudaEventSynchronize(opEnd));
        float milliseconds = 0;
        checked(cudaEventElapsedTime(&milliseconds, opStart, opFine));
        fineSpmvSeconds += 0.001 * milliseconds;
        checked(cudaEventElapsedTime(&milliseconds, opFine, opEnd));
        wellSeconds += 0.001 * milliseconds;
        ++matvecs;
    }

    void configure(double tolerance, int maxIterations)
    {
        XianzhiKrylov_ExternalPreconditioner preconditioner{};
        preconditioner.abi_version = 1;
        preconditioner.struct_size = sizeof(preconditioner);
        preconditioner.context = this;
        preconditioner.properties = XIANZHI_PC_FIXED_LINEAR;
        preconditioner.apply = +[](void* context, uint32_t, uint64_t, const double* in, double* out, void*) {
            auto& self = *static_cast<IluDevice*>(context);
            return callback([&] {
                const auto start = Clock::now();
                ok(xk_internal_cprw_ilu_apply(self.ilu, in, out), "fine ILU0 apply");
                self.iluApplySeconds += elapsed(start);
            });
        };

        XianzhiKrylov_OperatorOptions options{};
        options.abi_version = 4;
        options.struct_size = sizeof(options);
        options.relative_to_initial = 1;
        options.context = this;
        options.external_preconditioner = &preconditioner;
        options.apply_operator = +[](void* context, uint64_t, const double* in, double* out, void*) {
            return callback([&] { static_cast<IluDevice*>(context)->fullOperator(in, out); });
        };

        XianzhiKrylov_SolverConfig config{};
        config.abi_version = 3;
        config.struct_size = sizeof(config);
        config.method = XIANZHI_KRYLOV_METHOD_BICGSTAB;
        config.preconditioner = XIANZHI_KRYLOV_PRECONDITIONER_EXTERNAL;
        config.relative_tolerance = tolerance;
        config.max_iterations = maxIterations;
        config.restart = 1;
        config.zero_initial_guess = 1;
        config.ordering_policy = XIANZHI_KRYLOV_ORDERING_NATIVE;
        config.method_options = reinterpret_cast<const XianzhiKrylov_MethodOptions*>(&options);
        ok(xk_internal_cprw_krylov_create(&config, &linearSystem, &solver), "BiCGStab create");
    }

    static XianzhiKrylov_Status run(XianzhiKrylov_Solver* solver,
                                    Xianzhi_DeviceLinearSystem& s,
                                    XianzhiKrylov_SolverResult& result)
    {
        XianzhiKrylov_SolveInput input{};
        input.abi_version = XIANZHI_KRYLOV_SOLVE_INPUT_ABI_VERSION;
        input.struct_size = sizeof(input);
        input.initial_guess = XIANZHI_KRYLOV_INITIAL_ZERO;
        input.flags = XIANZHI_KRYLOV_SOLVE_VERIFY_RESIDUAL;
        input.expected_topology_epoch = s.topology_epoch;
        input.expected_values_epoch = s.values_epoch;
        input.expected_rhs_epoch = ++s.rhs_epoch;
        result.abi_version = 3;
        result.struct_size = sizeof(result);
        return xk_internal_cprw_krylov_solve(solver, &s, &input, &result);
    }
};

void IluDeviceDelete::operator()(IluDevice* p) const
{
    delete p;
}

IluDeviceHandle createIluDevice(std::span<const uint64_t> rows,
                                std::span<const uint64_t> columns,
                                double tolerance,
                                int maxIterations)
{
    if (rows.size() < 2 || (rows.size() - 1) > std::numeric_limits<int>::max() / 3 ||
        columns.empty() || rows.back() != columns.size()) {
        throw std::invalid_argument("Invalid Norne ILU structure");
    }
    IluDeviceHandle p(new IluDevice);
    p->n = static_cast<int>(rows.size() - 1);
    checked(cudaGetDevice(&p->device));
    checked(cudaEventCreate(&p->opStart));
    checked(cudaEventCreate(&p->opFine));
    checked(cudaEventCreate(&p->opEnd));
    putNew(p->row, rows);
    putNew(p->col, columns);
    p->a.allocate(9 * columns.size());
    p->rhs.allocate(3 * p->n);
    p->x.allocate(3 * p->n);
    checked(cudaMemset(p->x.ptr, 0, p->x.count * sizeof(double)));
    p->linearSystem = p->system(p->rhs, p->x);
    const auto start = Clock::now();
    ok(xk_internal_cprw_ilu_create(p->n, p->col.count, p->row.ptr, p->col.ptr, nullptr, &p->ilu),
       "fine ILU0 analysis");
    std::fprintf(stderr, "NORNE_ILU_ANALYSIS rows=%d blocks=%zu analysis_s=%.9f\n",
                 p->n, p->col.count, elapsed(start));
    p->configure(tolerance, maxIterations);
    return p;
}

IluPrepareTiming prepareIluDevice(IluDevice& p,
                                  std::span<const double> values,
                                  void* well,
                                  int (*wellApply)(void*, const double*, double*))
{
    p.well = well;
    p.wellApply = wellApply;
    const auto start = Clock::now();
    p.a.put(values);
    const auto factorStart = Clock::now();
    ok(xk_internal_cprw_ilu_prepare(p.ilu, p.a.ptr), "fine ILU0 factor");
    ++p.linearSystem.values_epoch;
    return {std::chrono::duration<double>(factorStart - start).count(), elapsed(factorStart)};
}

IluResult solveIluDevice(IluDevice& p, std::span<const double> rhs, std::span<double> solution)
{
    p.fineSpmvSeconds = 0.0;
    p.wellSeconds = 0.0;
    p.iluApplySeconds = 0.0;
    p.matvecs = 0;
    const auto uploadStart = Clock::now();
    p.rhs.put(rhs);
    const auto solveStart = Clock::now();
    const double upload = std::chrono::duration<double>(solveStart - uploadStart).count();
    XianzhiKrylov_SolverResult result{};
    const auto status = p.run(p.solver, p.linearSystem, result);
    const double solve = elapsed(solveStart);
    if (status != XIANZHI_KRYLOV_CONVERGED && status != XIANZHI_KRYLOV_MAX_ITERATIONS) {
        ok(status, "BiCGStab solve");
    }
    const auto downloadStart = Clock::now();
    p.x.get(solution);
    const double download = elapsed(downloadStart);
    return {status == XIANZHI_KRYLOV_CONVERGED,
            result.iterations,
            result.residual_scale > 0 ? result.true_residual / result.residual_scale : result.true_residual,
            solve,
            upload,
            download,
            p.fineSpmvSeconds,
            p.wellSeconds,
            p.iluApplySeconds,
            p.matvecs};
}

} // namespace Opm::detail
