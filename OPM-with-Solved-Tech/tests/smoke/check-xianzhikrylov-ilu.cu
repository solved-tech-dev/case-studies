#include <XianzhiKrylov_CprComponents.h>
#include <XianzhiKrylov_MethodOptions.h>
#include <cuda_runtime.h>

#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <string>

void requireCuda(cudaError_t s)
{
    if (s != cudaSuccess) {
        throw std::runtime_error(cudaGetErrorString(s));
    }
}

void requireStatus(XianzhiKrylov_Status status, XianzhiKrylov_Status expected, const char* phase)
{
    if (status != expected) {
        std::fprintf(stderr, "%s: returned %d expected %d\n", phase, status, expected);
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
class Device
{
public:
    T* data = nullptr;
    explicit Device(const std::array<T, N>& host)
    {
        requireCuda(cudaMalloc(&data, N * sizeof(T)));
        requireCuda(cudaMemcpy(data, host.data(), N * sizeof(T), cudaMemcpyHostToDevice));
    }
    ~Device() { cudaFree(data); }
    std::array<T, N> read() const
    {
        std::array<T, N> output;
        requireCuda(cudaMemcpy(output.data(), data, N * sizeof(T), cudaMemcpyDeviceToHost));
        return output;
    }
};

__global__ void multiply(int n, const double* in, double* out)
{
    const int i = threadIdx.x;
    if (i < n) {
        out[i] = (2 + i) * in[i];
    }
}

int op(void*, uint64_t size, const double* in, double* out, void* stream)
{
    multiply<<<1, 32, 0, static_cast<cudaStream_t>(stream)>>>(size, in, out);
    return cudaGetLastError() == cudaSuccess ? 0 : 1;
}

int pc(void*, uint32_t, uint64_t n, const double* in, double* out, void* stream)
{
    return cudaMemcpyAsync(out, in, n * sizeof(double), cudaMemcpyDeviceToDevice, static_cast<cudaStream_t>(stream)) ==
                   cudaSuccess
               ? 0
               : 1;
}

int main()
{
    if (procMapsContains("liblinglongamg")) {
        std::fprintf(stderr, "Isolation violation: liblinglongamg is loaded in XianzhiKrylov-only smoke test\n");
        return 4;
    }

    Device<double, 6> b({2, 3, 4, 5, 6, 7}), x({});
    Xianzhi_DeviceLinearSystem s{};
    s.abi_version = 3;
    s.struct_size = sizeof(s);
    s.matrix_format = XIANZHI_MATRIX_FORMAT_OPERATOR;
    s.scalar_type = XIANZHI_SCALAR_TYPE_FP64;
    s.block_rows = s.block_columns = s.vector_length = 6;
    s.block_size = s.rhs_count = 1;
    s.lifetime_token = &b;
    s.device_rhs = b.data;
    s.device_solution = x.data;
    s.topology_epoch = s.values_epoch = s.rhs_epoch = 1;

    for (uint32_t method : {XIANZHI_KRYLOV_METHOD_BICGSTAB, XIANZHI_KRYLOV_METHOD_GMRES}) {
        XianzhiKrylov_ExternalPreconditioner pre{};
        pre.abi_version = 1;
        pre.struct_size = sizeof(pre);
        pre.apply = pc;
        pre.properties = method == XIANZHI_KRYLOV_METHOD_GMRES ? XIANZHI_PC_FIXED_LINEAR : 0;

        XianzhiKrylov_OperatorOptions o{};
        o.abi_version = 4;
        o.struct_size = sizeof(o);
        o.apply_operator = op;
        o.external_preconditioner = &pre;

        XianzhiKrylov_SolverConfig c{};
        c.abi_version = 3;
        c.struct_size = sizeof(c);
        c.method = method;
        c.preconditioner = XIANZHI_KRYLOV_PRECONDITIONER_EXTERNAL;
        c.zero_initial_guess = 1;
        c.max_iterations = 20;
        c.restart = 10;
        c.relative_tolerance = 1e-12;
        c.method_options = reinterpret_cast<const XianzhiKrylov_MethodOptions*>(&o);

        XianzhiKrylov_Solver* solver = nullptr;
        requireStatus(xk_internal_cprw_krylov_create(&c, &s, &solver), XIANZHI_KRYLOV_SUCCESS, "create");

        XianzhiKrylov_SolveInput input{};
        input.abi_version = 1;
        input.struct_size = sizeof(input);
        input.initial_guess = XIANZHI_KRYLOV_INITIAL_ZERO;
        input.flags = XIANZHI_KRYLOV_SOLVE_VERIFY_RESIDUAL;
        input.expected_topology_epoch = s.topology_epoch;
        input.expected_values_epoch = s.values_epoch;
        input.expected_rhs_epoch = s.rhs_epoch;

        XianzhiKrylov_SolverResult result{};
        result.abi_version = 3;
        result.struct_size = sizeof(result);

        auto status = xk_internal_cprw_krylov_solve(solver, &s, &input, &result);
        xk_internal_cprw_krylov_destroy(solver);
        requireStatus(status, XIANZHI_KRYLOV_CONVERGED, "solve");
        for (double value : x.read()) {
            if (std::abs(value - 1) > 1e-11) {
                return 2;
            }
        }
        std::printf("XIANZHIKRYLOV_COMPONENT_PASS method=%u iterations=%u true_residual=%.17g\n",
                    method, result.iterations, result.true_residual);
    }

    Device<uint64_t, 3> row({0, 2, 4});
    Device<uint64_t, 4> col({0, 1, 0, 1});
    Device<double, 36> values({4, 0, 0, 0, 5, 0, 0, 0, 6, 1, 0, 0, 0, 1, 0, 0, 0, 1,
                               2, 0, 0, 0, 2, 0, 0, 0, 2, 6, 0, 0, 0, 7, 0, 0, 0, 8});
    Device<double, 6> iluRhs({5, 6, 7, 8, 9, 10});
    XkCprwBlockIlu* ilu = nullptr;
    requireStatus(xk_internal_cprw_ilu_create(2, 4, row.data, col.data, nullptr, &ilu),
                  XIANZHI_KRYLOV_SUCCESS, "ilu create");
    for (int refresh = 0; refresh < 2; ++refresh) {
        requireStatus(xk_internal_cprw_ilu_prepare(ilu, values.data), XIANZHI_KRYLOV_SUCCESS, "ilu numeric");
        for (int apply = 0; apply < 2; ++apply) {
            requireStatus(xk_internal_cprw_ilu_apply(ilu, iluRhs.data, x.data), XIANZHI_KRYLOV_SUCCESS, "ilu apply");
            for (double value : x.read()) {
                if (std::abs(value - 1) > 1e-12) {
                    return 3;
                }
            }
        }
    }
    xk_internal_cprw_ilu_destroy(ilu);
    std::puts("XIANZHIKRYLOV_COMPONENT_PASS native_bsr3_ilu numeric_refresh=2 apply=4 isolated_no_linglongamg=1");
    return 0;
}
