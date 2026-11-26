#include "NorneCprwDevice.hpp"
#include <XianzhiKrylov_CprComponents.h>
#include <XianzhiKrylov_MethodOptions.h>
#if __has_include(<linglongamg.h>)
#include <linglongamg.h>
#else
#include "linglongamg_sdk.h"
#endif
#include <cuda_runtime.h>
#include <dlfcn.h>

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <string>

namespace Opm::detail {
namespace {

using Clock = std::chrono::steady_clock;
using LinglongApplyCycleFn = linglongamg_status_t (*)(linglongamg_handle_t, const linglongamg_vec_t*, linglongamg_vec_t*);

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
            throw std::logic_error("CPRW device data changed size");
        }
        checked(cudaMemcpy(ptr, a.data(), count * sizeof(T), cudaMemcpyHostToDevice));
    }

    void get(std::span<T> a) const
    {
        if (a.size() != count) {
            throw std::logic_error("CPRW result changed size");
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

__global__ void pressureSpmv(int n, const int* row, const int* col, const double* a, const double* x, double* y)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) {
        return;
    }
    double v = 0.0;
    for (int k = row[i]; k < row[i + 1]; ++k) {
        v += a[k] * x[col[k]];
    }
    y[i] = v;
}

__global__ void restrictPressure(int n, int np, const double* weight, const double* residual, double* out)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= np) {
        return;
    }
    double value = 0.0;
    if (i < n) {
        for (int j = 0; j < 3; ++j) {
            value += weight[3 * i + j] * residual[3 * i + j];
        }
    }
    out[i] = value;
}

__global__ void prolongPressure(int n, int index, const double* p, double* out)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) {
        return;
    }
    for (int j = 0; j < 3; ++j) {
        out[3 * i + j] = j == index ? p[i] : 0.0;
    }
}

__global__ void combine(int n, double* y, const double* x, double sign)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) {
        y[i] = x[i] + sign * y[i];
    }
}

struct Amg
{
    void* library = nullptr;
    linglongamg_handle_t solver = nullptr;
    linglongamg_csr_t* matrix = nullptr;
    decltype(&linglongamg_create_from_config) create{};
    decltype(&linglongamg_read_csr) read{};
    decltype(&linglongamg_setup) setup{};
    decltype(&linglongamg_reset) reset{};
    LinglongApplyCycleFn cycle{};
    decltype(&linglongamg_solve) solve{};
    decltype(&linglongamg_get_solve_info) solveInfo{};
    decltype(&linglongamg_destroy) destroy{};
    decltype(&linglongamg_csr_destroy) destroyCsr{};

    ~Amg()
    {
        if (solver && destroy) {
            destroy(solver);
        }
        if (matrix && destroyCsr) {
            destroyCsr(matrix);
        }
        if (library) {
            dlclose(library);
        }
    }

    template <class Function>
    void get(Function& destination, const char* name)
    {
        destination = reinterpret_cast<Function>(dlsym(library, name));
        if (!destination) {
            throw std::runtime_error(std::string("LinglongAMG missing ") + name);
        }
    }

    static void check(linglongamg_status_t status, const char* action)
    {
        if (status != LINGLONGAMG_STATUS_SUCCESS) {
            throw std::runtime_error(std::string("LinglongAMG ") + action + ": status " + std::to_string(status));
        }
    }

    void open(const char* path, const char* cfg)
    {
        library = dlopen(path, RTLD_NOW | RTLD_LOCAL);
        if (!library) {
            throw std::runtime_error(dlerror());
        }
        get(create, "linglongamg_create_from_config");
        get(read, "linglongamg_read_csr");
        get(setup, "linglongamg_setup");
        get(reset, "linglongamg_reset");
        cycle = reinterpret_cast<LinglongApplyCycleFn>(dlsym(library, "linglongamg_apply_cycle"));
        if (!cycle) {
            if (!std::getenv("NORNE_AMG_ALLOW_SINGLE_STEP")) {
                throw std::runtime_error("LinglongAMG missing linglongamg_apply_cycle");
            }
            get(solve, "linglongamg_solve");
            get(solveInfo, "linglongamg_get_solve_info");
        }
        get(destroy, "linglongamg_destroy");
        get(destroyCsr, "linglongamg_csr_destroy");
        check(create(&solver, cfg), "create");
        std::fprintf(stderr, "NORNE_CPRW_LINGLONG api=%s\n", cycle ? "apply_cycle" : "amg_one_step");
    }
};

struct Timers
{
    double amgSetup{}, amgCycle{}, pressure{}, restriction{}, correction{}, iluFactor{}, iluApply{}, fineSpmv{}, well{},
        upload{};
    unsigned pressureCalls{}, pressureMax{}, amgCalls{}, outerMatvecs{};
};

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
        std::fprintf(stderr, "NORNE_CPRW_CALLBACK_FAILURE %s\n", error.what());
        return 1;
    }
}

} // namespace

struct CprwDevice
{
    int n = 0, np = 0, pressureIndex = 0, device = 0;
    cudaEvent_t opStart{}, opFine{}, opEnd{};
    Buffer<uint64_t> row, col;
    Buffer<int> pRow, pCol;
    Buffer<double> a, pa, weights, rhs, x, pRhs, pX, injected, corrected;
    XkCprwBlockIlu* ilu = nullptr;
    XianzhiKrylov_Solver *outer = nullptr, *inner = nullptr;
    Xianzhi_DeviceLinearSystem outerSystem{}, innerSystem{};
    Amg amg;
    void* well = nullptr;
    int (*wellApply)(void*, const double*, double*) = nullptr;
    Timers timers;
    bool prepared = false;

    ~CprwDevice()
    {
        if (outer) {
            xk_internal_cprw_krylov_destroy(outer);
        }
        if (inner) {
            xk_internal_cprw_krylov_destroy(inner);
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
        timers.fineSpmv += 0.001 * milliseconds;
        checked(cudaEventElapsedTime(&milliseconds, opFine, opEnd));
        timers.well += 0.001 * milliseconds;
        ++timers.outerMatvecs;
    }

    void pressureApply(const double* in, double* out)
    {
        const auto start = Clock::now();
        linglongamg_vec_t b{np, const_cast<double*>(in)}, xOut{np, out};
        if (amg.cycle) {
            Amg::check(amg.cycle(amg.solver, &b, &xOut), "apply cycle");
        } else {
            checked(cudaMemset(out, 0, np * sizeof(double)));
            Amg::check(amg.solve(amg.solver, &b, &xOut), "AMG one step");
            // Commit device solution and verify single-step V-cycle execution.
            linglongamg_solve_info_t info{};
            Amg::check(amg.solveInfo(amg.solver, &info), "AMG step count and solution");
            if (info.iterations != 1) {
                throw std::runtime_error("LinglongAMG expected exactly one V-cycle for nonzero GMRES input");
            }
        }
        checked(cudaDeviceSynchronize());
        timers.amgCycle += elapsed(start);
        ++timers.amgCalls;
    }

    XianzhiKrylov_Solver* createKrylov(bool pressure, double tolerance, int maxIterations)
    {
        XianzhiKrylov_ExternalPreconditioner pc{};
        pc.abi_version = 1;
        pc.struct_size = sizeof(pc);
        pc.context = this;
        pc.properties = pressure ? XIANZHI_PC_FIXED_LINEAR : 0;
        pc.apply = pressure
            ? +[](void* p, uint32_t, uint64_t, const double* in, double* out, void*) {
                  return callback([&] { static_cast<CprwDevice*>(p)->pressureApply(in, out); });
              }
            : +[](void* p, uint32_t, uint64_t, const double* in, double* out, void*) {
                  return callback([&] { static_cast<CprwDevice*>(p)->precondition(in, out); });
              };
        XianzhiKrylov_OperatorOptions options{};
        options.abi_version = 4;
        options.struct_size = sizeof(options);
        options.relative_to_initial = 1;
        options.context = this;
        options.external_preconditioner = &pc;
        options.apply_operator = pressure
            ? +[](void* p, uint64_t, const double* in, double* out, void*) {
                  const auto& q = *static_cast<CprwDevice*>(p);
                  pressureSpmv<<<blocks(q.np), threads>>>(q.np, q.pRow.ptr, q.pCol.ptr, q.pa.ptr, in, out);
                  return cudaGetLastError() == cudaSuccess ? 0 : 1;
              }
            : +[](void* p, uint64_t, const double* in, double* out, void*) {
                  return callback([&] { static_cast<CprwDevice*>(p)->fullOperator(in, out); });
              };
        XianzhiKrylov_SolverConfig c{};
        c.abi_version = 3;
        c.struct_size = sizeof(c);
        c.method = pressure ? XIANZHI_KRYLOV_METHOD_GMRES : XIANZHI_KRYLOV_METHOD_BICGSTAB;
        c.preconditioner = XIANZHI_KRYLOV_PRECONDITIONER_EXTERNAL;
        c.relative_tolerance = tolerance;
        c.max_iterations = maxIterations;
        c.restart = pressure ? 10 : 1;
        c.zero_initial_guess = 1;
        c.ordering_policy = XIANZHI_KRYLOV_ORDERING_NATIVE;
        c.method_options = reinterpret_cast<const XianzhiKrylov_MethodOptions*>(&options);
        XianzhiKrylov_Solver* output = nullptr;
        ok(xk_internal_cprw_krylov_create(&c, pressure ? &innerSystem : &outerSystem, &output),
           pressure ? "pressure GMRES create" : "outer BiCGStab create");
        return output;
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

    void precondition(const double* in, double* out)
    {
        auto start = Clock::now();
        restrictPressure<<<blocks(np), threads>>>(n, np, weights.ptr, in, pRhs.ptr);
        checked(cudaDeviceSynchronize());
        timers.restriction += elapsed(start);

        start = Clock::now();
        XianzhiKrylov_SolverResult pressure{};
        const auto status = run(inner, innerSystem, pressure);
        timers.pressure += elapsed(start);
        ++timers.pressureCalls;
        if (status == XIANZHI_KRYLOV_MAX_ITERATIONS) {
            ++timers.pressureMax;
        } else if (status != XIANZHI_KRYLOV_CONVERGED) {
            ok(status, "pressure GMRES");
        }

        start = Clock::now();
        prolongPressure<<<blocks(n), threads>>>(n, pressureIndex, pX.ptr, injected.ptr);
        checked(cudaGetLastError());
        fullOperator(injected.ptr, corrected.ptr);
        combine<<<blocks(3 * n), threads>>>(3 * n, corrected.ptr, in, -1.0);
        checked(cudaDeviceSynchronize());
        timers.correction += elapsed(start);

        start = Clock::now();
        ok(xk_internal_cprw_ilu_apply(ilu, corrected.ptr, out), "fine block ILU0 apply");
        timers.iluApply += elapsed(start);
        combine<<<blocks(3 * n), threads>>>(3 * n, out, injected.ptr, 1.0);
        checked(cudaDeviceSynchronize());
    }
};

void CprwDelete::operator()(CprwDevice* p) const
{
    delete p;
}

CprwHandle createCprwDevice(const CprwMatrixData& data,
                            const char* amgPath,
                            const char* amgCfg,
                            double tolerance,
                            int maxIterations)
{
    CprwHandle p(new CprwDevice);
    if (data.rows.size() < 2 || data.pressureRows.size() < data.rows.size() || data.pressureIndex < 0 ||
        data.pressureIndex > 2) {
        throw std::invalid_argument("Invalid CPRW structure");
    }
    checked(cudaGetDevice(&p->device));
    checked(cudaEventCreate(&p->opStart));
    checked(cudaEventCreate(&p->opFine));
    checked(cudaEventCreate(&p->opEnd));
    p->n = int(data.rows.size() - 1);
    p->np = int(data.pressureRows.size() - 1);
    p->pressureIndex = data.pressureIndex;
    putNew(p->row, data.rows);
    putNew(p->col, data.columns);
    putNew(p->pRow, data.pressureRows);
    putNew(p->pCol, data.pressureColumns);
    p->a.allocate(9 * data.columns.size());
    p->pa.allocate(data.pressureColumns.size());
    for (auto* b : {&p->weights, &p->rhs, &p->x, &p->injected, &p->corrected}) {
        b->allocate(3 * p->n);
    }
    p->pRhs.allocate(p->np);
    p->pX.allocate(p->np);
    checked(cudaMemset(p->x.ptr, 0, p->x.count * sizeof(double)));
    checked(cudaMemset(p->pX.ptr, 0, p->pX.count * sizeof(double)));
    p->outerSystem = p->system(p->rhs, p->x);
    p->innerSystem = p->system(p->pRhs, p->pX);
    p->amg.open(amgPath, amgCfg);
    ok(xk_internal_cprw_ilu_create(p->n, p->col.count, p->row.ptr, p->col.ptr, nullptr, &p->ilu),
       "fine block ILU0 analysis");
    p->inner = p->createKrylov(true, 1e-6, 100);
    p->outer = p->createKrylov(false, tolerance, maxIterations);
    return p;
}

void prepareCprwDevice(CprwDevice& p,
                       const CprwNumericData& data,
                       void* well,
                       int (*wellApply)(void*, const double*, double*))
{
    const bool reset = p.prepared;
    p.prepared = false;
    p.timers = {};
    p.well = well;
    p.wellApply = wellApply;
    auto start = Clock::now();
    if (reset) {
        Amg::check(p.amg.reset(p.amg.solver), "reset");
    }
    p.a.put(data.matrix);
    p.weights.put(data.weights);
    p.pa.put(data.pressure);
    p.timers.upload = elapsed(start);
    start = Clock::now();
    if (!p.amg.matrix) {
        Amg::check(p.amg.read(&p.amg.matrix, p.np, p.np, int(p.pa.count), p.pRow.ptr, p.pCol.ptr, p.pa.ptr),
                   "CSR bind");
    }
    Amg::check(p.amg.setup(p.amg.solver, p.amg.matrix), "setup");
    checked(cudaDeviceSynchronize());
    p.timers.amgSetup = elapsed(start);
    start = Clock::now();
    ok(xk_internal_cprw_ilu_prepare(p.ilu, p.a.ptr), "fine block ILU0 prepare");
    p.timers.iluFactor = elapsed(start);
    ++p.innerSystem.values_epoch;
    ++p.outerSystem.values_epoch;
    p.prepared = true;
}

void auditCprwOperator(CprwDevice& p, std::span<const double> x, std::span<double> result)
{
    const auto previous = p.timers;
    p.rhs.put(x);
    p.fullOperator(p.rhs.ptr, p.x.ptr);
    p.x.get(result);
    p.timers = previous;
}

void auditCprwCorrection(CprwDevice& p,
                         std::span<const double> rhs,
                         std::span<double> coarseRhs,
                         std::span<double> coarseSolution,
                         std::span<double> result)
{
    const auto previous = p.timers;
    p.rhs.put(rhs);
    p.precondition(p.rhs.ptr, p.x.ptr);
    p.pRhs.get(coarseRhs);
    p.pX.get(coarseSolution);
    p.x.get(result);
    p.timers = previous;
}

CprwResult solveCprwDevice(CprwDevice& p, std::span<const double> rhs, std::span<double> solution)
{
    if (!p.prepared) {
        throw std::logic_error("CPRW numeric setup missing");
    }
    auto start = Clock::now();
    p.rhs.put(rhs);
    const auto rhsUpload = elapsed(start);
    start = Clock::now();
    XianzhiKrylov_SolverResult result{};
    const auto status = p.run(p.outer, p.outerSystem, result);
    const auto seconds = elapsed(start);
    if (status != XIANZHI_KRYLOV_CONVERGED && status != XIANZHI_KRYLOV_MAX_ITERATIONS) {
        ok(status, "outer BiCGStab");
    }
    start = Clock::now();
    p.x.get(solution);
    const auto download = elapsed(start);
    const auto& t = p.timers;
    std::fprintf(stderr,
                 "NORNE_CPRW_STAGES outer_s=%.9f rhs_h2d_s=%.9f x_d2h_s=%.9f values_h2d_s=%.9f "
                 "amg_setup_s=%.9f amg_cycle_s=%.9f pressure_gmres_s=%.9f restriction_s=%.9f "
                 "correction_s=%.9f fine_spmv_s=%.9f well_s=%.9f ilu_factor_s=%.9f ilu_apply_s=%.9f "
                 "pressure_calls=%u pressure_max=%u amg_calls=%u matvecs=%u\n",
                 seconds, rhsUpload, download, t.upload, t.amgSetup, t.amgCycle, t.pressure, t.restriction,
                 t.correction, t.fineSpmv, t.well, t.iluFactor, t.iluApply, t.pressureCalls, t.pressureMax,
                 t.amgCalls, t.outerMatvecs);
    return {status == XIANZHI_KRYLOV_CONVERGED,
            result.iterations,
            result.residual_scale > 0 ? result.true_residual / result.residual_scale : result.true_residual,
            seconds};
}

} // namespace Opm::detail
