#ifndef OPM_LINGLONGAMG_PRECONDITIONER_HPP
#define OPM_LINGLONGAMG_PRECONDITIONER_HPP

// OPM-side adapter for the unchanged LinglongAMG CUDA single-GPU SDK.
#include <opm/simulators/linalg/PreconditionerWithUpdate.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>

#include <opm/simulators/linalg/linglongamg_sdk.h>
#include <opm/simulators/linalg/AMGCallTiming.hpp>
#include <cuda_runtime.h>
#include <dlfcn.h>

#include <algorithm>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace Opm {

template <class Matrix, class Vector>
class LinglongAMGPreconditioner final : public Dune::PreconditionerWithUpdate<Vector, Vector>
{
public:
    LinglongAMGPreconditioner(const Matrix& matrix, const PropertyTree& options)
        : matrix_(matrix)
    {
        static_assert(Matrix::block_type::rows == 1 && Matrix::block_type::cols == 1);
        static_assert(std::is_same_v<typename Matrix::field_type, double>);
        static_assert(std::is_same_v<linglongamg_real_t, double>, "LinglongAMG OPM ABI is FP64");

        std::string defaultLibrary = "liblinglongamg.so";
        if (const char* envLib = std::getenv("LINGLONGAMG_OPM_LIBRARY"); envLib && *envLib) {
            defaultLibrary = envLib;
        } else if (const char* root = std::getenv("LINGLONGAMG_ROOT"); root && *root) {
            defaultLibrary = std::string(root) + "/lib/liblinglongamg.so";
        }
        const auto library = options.get<std::string>("library", defaultLibrary);

        std::string defaultConfig;
        if (const char* envCfg = std::getenv("LINGLONGAMG_OPM_CONFIG"); envCfg && *envCfg) {
            defaultConfig = envCfg;
        }
        const auto config = options.get<std::string>("config_file", defaultConfig);
        if (config.empty()) {
            throw std::runtime_error("LinglongAMG config_file is required");
        }

        library_ = dlopen(library.c_str(), RTLD_NOW | RTLD_LOCAL);
        if (!library_) {
            throw std::runtime_error(std::string("Cannot load LinglongAMG SDK: ") + dlerror());
        }
        try {
            create_ = symbol<decltype(&linglongamg_create_from_config)>("linglongamg_create_from_config");
            read_csr_ = symbol<decltype(&linglongamg_read_csr)>("linglongamg_read_csr");
            csr_destroy_ = symbol<decltype(&linglongamg_csr_destroy)>("linglongamg_csr_destroy");
            setup_ = symbol<decltype(&linglongamg_setup)>("linglongamg_setup");
            solve_ = symbol<decltype(&linglongamg_solve)>("linglongamg_solve");
            get_info_ = symbol<decltype(&linglongamg_get_solve_info)>("linglongamg_get_solve_info");
            reset_ = symbol<decltype(&linglongamg_reset)>("linglongamg_reset");
            destroy_ = symbol<decltype(&linglongamg_destroy)>("linglongamg_destroy");
            check(create_(&solver_, config.c_str()), "create_from_config");
            update();
        } catch (...) {
            release();
            throw;
        }
    }

    ~LinglongAMGPreconditioner() override
    {
        AMGCallTiming::linglong().report();
        release();
    }

    void update() override
    {
        const auto size = matrix_.N();
        const auto nonzeros = matrix_.nonzeroes();
        if (size == 0 || size > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
            nonzeros > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
            throw std::runtime_error("LinglongAMG needs nonempty 32-bit CSR indices");
        }
        const int n = static_cast<int>(size);
        std::vector<int> rows(size + 1, 0), columns;
        std::vector<double> values;
        columns.reserve(nonzeros);
        values.reserve(nonzeros);
        for (auto row = matrix_.begin(); row != matrix_.end(); ++row) {
            for (auto entry = row->begin(); entry != row->end(); ++entry) {
                columns.push_back(static_cast<int>(entry.index()));
                values.push_back((*entry)[0][0]);
            }
            rows[row.index() + 1] = static_cast<int>(columns.size());
        }

        // linglongamg_setup may alias the input CSR, so allocations survive
        // until the next setup has destroyed the previous hierarchy.
        if (capacity_n_ != n || capacity_nnz_ != static_cast<int>(nonzeros)) {
            if (solver_ready_) {
                check(reset_(solver_), "reset");
            }
            if (descriptor_) {
                csr_destroy_(descriptor_);
                descriptor_ = nullptr;
            }
            freeBuffers();
            capacity_n_ = n;
            capacity_nnz_ = static_cast<int>(nonzeros);
            cudaCheck(cudaMalloc(&d_rows_, rows.size() * sizeof(int)), "cudaMalloc(rows)");
            cudaCheck(cudaMalloc(&d_columns_, columns.size() * sizeof(int)), "cudaMalloc(columns)");
            cudaCheck(cudaMalloc(&d_values_, values.size() * sizeof(double)), "cudaMalloc(values)");
            cudaCheck(cudaMalloc(&d_rhs_, n * sizeof(double)), "cudaMalloc(rhs)");
            cudaCheck(cudaMalloc(&d_x_, n * sizeof(double)), "cudaMalloc(x)");
            check(read_csr_(&descriptor_, n, n, capacity_nnz_,
                            static_cast<int*>(d_rows_), static_cast<int*>(d_columns_),
                            static_cast<double*>(d_values_)),
                  "read_csr");
        }
        cudaCheck(cudaMemcpy(d_rows_, rows.data(), rows.size() * sizeof(int), cudaMemcpyHostToDevice), "copy rows");
        cudaCheck(cudaMemcpy(d_columns_, columns.data(), columns.size() * sizeof(int), cudaMemcpyHostToDevice),
                  "copy columns");
        cudaCheck(cudaMemcpy(d_values_, values.data(), values.size() * sizeof(double), cudaMemcpyHostToDevice),
                  "copy values");
        AMGCallTiming::linglong().measure(
            true, n, nonzeros,
            [] { cudaCheck(cudaDeviceSynchronize(), "setup timing synchronize"); },
            [&] { check(setup_(solver_, descriptor_), "setup"); });
        solver_ready_ = true;
    }

    void pre(Vector&, Vector&) override {}
    void post(Vector&) override {}

    void apply(Vector& result, const Vector& defect) override
    {
        const int n = capacity_n_;
        std::vector<double> input(n), output(n);
        for (int i = 0; i < n; ++i) {
            input[i] = defect[i][0];
        }
        cudaCheck(cudaMemcpy(d_rhs_, input.data(), n * sizeof(double), cudaMemcpyHostToDevice), "copy rhs");
        cudaCheck(cudaMemset(d_x_, 0, n * sizeof(double)), "zero x");
        linglongamg_vec_t rhs{n, static_cast<double*>(d_rhs_)};
        linglongamg_vec_t x{n, static_cast<double*>(d_x_)};
        auto& timing = AMGCallTiming::linglong();
        timing.measure(
            false, n, capacity_nnz_,
            [] { cudaCheck(cudaDeviceSynchronize(), "solve timing synchronize"); },
            [&] { check(solve_(solver_, &rhs, &x), "solve"); });
        // SDK inner-call residual; not an independently recomputed application residual.
        linglongamg_solve_info_t info{};
        if (get_info_(solver_, &info) == LINGLONGAMG_STATUS_SUCCESS) {
            ++timing.info_calls;
            timing.last_converged = info.converged;
            timing.last_iterations = info.iterations;
            timing.last_initial_residual = info.initial_residual;
            timing.last_final_residual = info.final_residual;
            timing.last_relative_residual = info.relative_residual;
        } else {
            ++timing.info_failures;
        }
        cudaCheck(cudaMemcpy(output.data(), d_x_, n * sizeof(double), cudaMemcpyDeviceToHost), "copy x");
        for (int i = 0; i < n; ++i) {
            result[i][0] = output[i];
        }
    }

    Dune::SolverCategory::Category category() const override { return Dune::SolverCategory::sequential; }
    bool hasPerfectUpdate() const override { return true; }

private:
    template <class Function>
    Function symbol(const char* name)
    {
        auto ptr = dlsym(library_, name);
        if (!ptr) {
            throw std::runtime_error(std::string("Missing LinglongAMG symbol: ") + name);
        }
        return reinterpret_cast<Function>(ptr);
    }
    static void check(linglongamg_status_t code, const char* operation)
    {
        if (code != LINGLONGAMG_STATUS_SUCCESS) {
            throw std::runtime_error(std::string("LinglongAMG ") + operation + " failed: status=" + std::to_string(code));
        }
    }
    static void cudaCheck(cudaError_t code, const char* operation)
    {
        if (code != cudaSuccess) {
            throw std::runtime_error(std::string("LinglongAMG ") + operation + ": " + cudaGetErrorString(code));
        }
    }
    void freeBuffers() noexcept
    {
        cudaFree(d_rows_);
        cudaFree(d_columns_);
        cudaFree(d_values_);
        cudaFree(d_rhs_);
        cudaFree(d_x_);
        d_rows_ = d_columns_ = d_values_ = d_rhs_ = d_x_ = nullptr;
    }
    void release() noexcept
    {
        if (solver_ && destroy_) {
            destroy_(solver_);
        }
        solver_ = nullptr;
        if (descriptor_ && csr_destroy_) {
            csr_destroy_(descriptor_);
        }
        descriptor_ = nullptr;
        freeBuffers();
        if (library_) {
            dlclose(library_);
        }
        library_ = nullptr;
    }

    const Matrix& matrix_;
    void* library_ = nullptr;
    linglongamg_handle_t solver_ = nullptr;
    linglongamg_csr_t* descriptor_ = nullptr;
    void *d_rows_ = nullptr, *d_columns_ = nullptr, *d_values_ = nullptr, *d_rhs_ = nullptr, *d_x_ = nullptr;
    int capacity_n_ = 0, capacity_nnz_ = 0;
    bool solver_ready_ = false;
    decltype(&linglongamg_create_from_config) create_ = nullptr;
    decltype(&linglongamg_read_csr) read_csr_ = nullptr;
    decltype(&linglongamg_csr_destroy) csr_destroy_ = nullptr;
    decltype(&linglongamg_setup) setup_ = nullptr;
    decltype(&linglongamg_solve) solve_ = nullptr;
    decltype(&linglongamg_get_solve_info) get_info_ = nullptr;
    decltype(&linglongamg_reset) reset_ = nullptr;
    decltype(&linglongamg_destroy) destroy_ = nullptr;
};

} // namespace Opm

#endif
