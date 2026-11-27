/*
  Copyright 2026 Equinor ASA / Solved-Tech Case Studies

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPM_NORNE_XIANZHI_HEADER_INCLUDED
#define OPM_NORNE_XIANZHI_HEADER_INCLUDED

#include <XianzhiKrylovAPI.h>

#include <cuda_runtime_api.h>
#include <dlfcn.h>

#include <dune/istl/operators.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>

#include <opm/simulators/linalg/PreconditionerWithUpdate.hpp>

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>

namespace Opm::detail
{

using NorneClock = std::chrono::steady_clock;

inline double norneElapsedSeconds(NorneClock::time_point start)
{
    return std::chrono::duration<double>(NorneClock::now() - start).count();
}

inline void norneCheckCuda(cudaError_t error, const char* what)
{
    if (error != cudaSuccess) {
        throw std::runtime_error(std::string("OPM XianzhiKrylov CUDA error (") + what + "): " +
                                 cudaGetErrorString(error));
    }
}

template <class T>
class NorneDeviceBuffer
{
public:
    NorneDeviceBuffer() = default;
    NorneDeviceBuffer(const NorneDeviceBuffer&) = delete;
    NorneDeviceBuffer& operator=(const NorneDeviceBuffer&) = delete;

    ~NorneDeviceBuffer()
    {
        reset();
    }

    void reset() noexcept
    {
        if (data_) {
            cudaFree(data_);
            data_ = nullptr;
            count_ = 0;
        }
    }

    void allocate(std::size_t count)
    {
        if (count == count_ && (data_ != nullptr || count == 0)) {
            return;
        }
        reset();
        count_ = count;
        if (count_ > 0) {
            norneCheckCuda(cudaMalloc(reinterpret_cast<void**>(&data_), bytes()), "cudaMalloc");
        }
    }

    void upload(const std::vector<T>& values)
    {
        if (values.size() != count_) {
            throw std::logic_error("NorneDeviceBuffer size mismatch on upload");
        }
        if (count_ > 0) {
            norneCheckCuda(cudaMemcpy(data_, values.data(), bytes(), cudaMemcpyHostToDevice),
                           "cudaMemcpyH2D");
        }
    }

    void download(std::vector<T>& values) const
    {
        if (values.size() != count_) {
            throw std::logic_error("NorneDeviceBuffer size mismatch on download");
        }
        if (count_ > 0) {
            norneCheckCuda(cudaMemcpy(values.data(), data_, bytes(), cudaMemcpyDeviceToHost),
                           "cudaMemcpyD2H");
        }
    }

    std::size_t bytes() const noexcept
    {
        return count_ * sizeof(T);
    }

    T* data() const noexcept
    {
        return data_;
    }

    xianzhikrylov_input_buffer_t input() const noexcept
    {
        return {data_, bytes()};
    }

private:
    T* data_ = nullptr;
    std::size_t count_ = 0;
};

template <class Matrix, class Vector>
class NorneXianzhi final : public Dune::InverseOperator<Vector, Vector>,
                           public Dune::PreconditionerWithUpdate<Vector, Vector>
{
public:
    static constexpr int BlockDim = 3;

    NorneXianzhi(const Matrix& A, const char* configPath, double tolerance)
        : A_(A), tolerance_(tolerance)
    {
        if (!configPath || *configPath == '\0') {
            throw std::invalid_argument("XianzhiKrylov configuration path must be non-empty");
        }
        using StatsFn = xianzhikrylov_status_t (*)(xianzhikrylov_handle_t, uint32_t*, double*, double*);
        statsFn_ = reinterpret_cast<StatsFn>(dlsym(RTLD_DEFAULT, "xk_internal_get_last_stats"));
        try {
            check(xianzhikrylov_create(0, &handle_), XIANZHIKRYLOV_SUCCESS, "create");
            check(xianzhikrylov_config_set(handle_, XIANZHIKRYLOV_CONFIG_FILE, configPath),
                  XIANZHIKRYLOV_SUCCESS, "config_set");
            rebuildStructure();
            updateNumericOnly();
        } catch (...) {
            release();
            throw;
        }
    }

    ~NorneXianzhi() override
    {
        release();
    }

    Dune::SolverCategory::Category category() const override
    {
        return Dune::SolverCategory::sequential;
    }

    bool hasPerfectUpdate() const override
    {
        return true;
    }

    void pre(Vector&, Vector&) override
    {
    }

    void post(Vector&) override
    {
    }

    void apply(Vector&, const Vector&) override
    {
        throw std::logic_error("NorneXianzhi implements the outer InverseOperator directly");
    }

    void update() override
    {
        if (structureChanged()) {
            rebuildStructure();
        }
        updateNumericOnly();
    }

    void apply(Vector& x, Vector& b, double tolerance, Dune::InverseOperatorResult& r) override
    {
        if (std::abs(tolerance - tolerance_) > 1e-12) {
            throw std::logic_error("XianzhiKrylov tolerance override differs from JSON configuration");
        }
        apply(x, b, r);
    }

    void apply(Vector& x, Vector& b, Dune::InverseOperatorResult& r) override
    {
        const auto start = NorneClock::now();
        for (std::size_t i = 0; i < A_.N(); ++i) {
            for (int j = 0; j < BlockDim; ++j) {
                rhs_[BlockDim * i + j] = b[i][j];
            }
        }
        auto stage = NorneClock::now();
        drhs_.upload(rhs_);
        const double rhsH2dSeconds = norneElapsedSeconds(stage);

        stage = NorneClock::now();
        const xianzhikrylov_status_t status =
            xianzhikrylov_execute(handle_, XIANZHIKRYLOV_PHASE_SOLVE, matrix_, xv_, bv_);
        const double solveSeconds = norneElapsedSeconds(stage);

        if (status != XIANZHIKRYLOV_CONVERGED && status != XIANZHIKRYLOV_MAX_ITERATIONS) {
            check(status, XIANZHIKRYLOV_CONVERGED, "solve");
        }

        stage = NorneClock::now();
        dx_.download(result_);
        const double downloadSeconds = norneElapsedSeconds(stage);

        for (std::size_t i = 0; i < A_.N(); ++i) {
            for (int j = 0; j < BlockDim; ++j) {
                x[i][j] = result_[BlockDim * i + j];
            }
        }

        uint32_t iterations = status == XIANZHIKRYLOV_CONVERGED ? 1u : 0u;
        double reduction = 0.0;
        if (statsFn_) {
            double residual = 0.0;
            double norm = 0.0;
            check(statsFn_(handle_, &iterations, &residual, &norm),
                  XIANZHIKRYLOV_SUCCESS, "read stats");
            reduction = norm > 0.0 ? residual / norm : residual;
        } else {
            Vector residualVec = b;
            const double rhsNorm = residualVec.two_norm();
            A_.mmv(x, residualVec);
            const double resNorm = residualVec.two_norm();
            reduction = rhsNorm > 0.0 ? resNorm / rhsNorm : resNorm;
        }

        r.clear();
        r.iterations = static_cast<int>(iterations);
        r.converged = (status == XIANZHIKRYLOV_CONVERGED);
        r.reduction = reduction;
        r.elapsed = norneElapsedSeconds(start);

        if (profileEnabled()) {
            std::fprintf(
                stderr,
                "OPM_XK_LINEAR rows=%zu seconds=%.9f sdk_solve_s=%.9f rhs_h2d_s=%.9f x_d2h_s=%.9f iterations=%u converged=%d relative=%.12g\n",
                A_.N(), r.elapsed, solveSeconds, rhsH2dSeconds, downloadSeconds,
                iterations, static_cast<int>(r.converged), r.reduction);
        }
    }

private:
    static bool profileEnabled() noexcept
    {
        return std::getenv("OPM_SOLVED_TECH_PROFILE") != nullptr ||
               std::getenv("NORNE_ILU_PROFILE") != nullptr;
    }

    bool structureChanged() const
    {
        if (offsets_.size() != A_.N() + 1 || columns_.size() != A_.nonzeroes()) {
            return true;
        }
        std::size_t block = 0;
        for (auto row = A_.begin(); row != A_.end(); ++row) {
            if (offsets_[row.index()] != block) {
                return true;
            }
            for (auto entry = row->begin(); entry != row->end(); ++entry) {
                if (block >= columns_.size() || columns_[block] != entry.index()) {
                    return true;
                }
                ++block;
            }
        }
        return block != columns_.size();
    }

    void rebuildStructure()
    {
        const auto start = NorneClock::now();
        destroyDescriptors();

        offsets_.assign(A_.N() + 1, 0);
        columns_.clear();
        columns_.reserve(A_.nonzeroes());
        std::vector<uint64_t> diagonals;
        diagonals.reserve(A_.N());

        for (auto row = A_.begin(); row != A_.end(); ++row) {
            const std::size_t before = diagonals.size();
            for (auto entry = row->begin(); entry != row->end(); ++entry) {
                if (entry.index() == row.index()) {
                    diagonals.push_back(columns_.size());
                }
                columns_.push_back(entry.index());
            }
            if (diagonals.size() != before + 1) {
                throw std::logic_error("XianzhiKrylov BSR3 requires exactly one diagonal block per row");
            }
            offsets_.at(row.index() + 1) = columns_.size();
        }

        values_.resize(static_cast<std::size_t>(BlockDim) * BlockDim * columns_.size());
        rhs_.resize(static_cast<std::size_t>(BlockDim) * A_.N());
        result_.resize(rhs_.size());

        doffsets_.allocate(offsets_.size());
        dcolumns_.allocate(columns_.size());
        ddiagonal_.allocate(diagonals.size());
        dvalues_.allocate(values_.size());
        drhs_.allocate(rhs_.size());
        dx_.allocate(result_.size());

        doffsets_.upload(offsets_);
        dcolumns_.upload(columns_);
        ddiagonal_.upload(diagonals);

        xianzhikrylov_matrix_data_t data{};
        data.abi_version = XIANZHIKRYLOV_ABI_VERSION;
        data.struct_size = sizeof(data);
        data.format = XIANZHIKRYLOV_BSR;
        data.value_type = XIANZHIKRYLOV_FP64;
        data.index_type = XIANZHIKRYLOV_INDEX_I64;
        data.block_size = BlockDim;
        data.block_layout = XIANZHIKRYLOV_BLOCK_ROW_MAJOR;
        data.memory_space = XIANZHIKRYLOV_DEVICE;
        data.device_id = 0;
        data.rows = static_cast<int64_t>(A_.N());
        data.columns = static_cast<int64_t>(A_.M());
        data.entries = static_cast<int64_t>(columns_.size());
        data.row_offsets = doffsets_.input();
        data.column_indices = dcolumns_.input();
        data.diagonal_positions = ddiagonal_.input();
        data.values = dvalues_.input();

        check(xianzhikrylov_matrix_create(&matrix_, &data), XIANZHIKRYLOV_SUCCESS, "matrix_create");
        check(xianzhikrylov_vector_create(
                  &xv_, static_cast<int64_t>(result_.size()), dx_.data(),
                  XIANZHIKRYLOV_FP64, XIANZHIKRYLOV_DEVICE),
              XIANZHIKRYLOV_SUCCESS, "solution_vector_create");
        check(xianzhikrylov_vector_create(
                  &bv_, static_cast<int64_t>(rhs_.size()), drhs_.data(),
                  XIANZHIKRYLOV_FP64, XIANZHIKRYLOV_DEVICE),
              XIANZHIKRYLOV_SUCCESS, "rhs_vector_create");
        check(xianzhikrylov_execute(handle_, XIANZHIKRYLOV_PHASE_ANALYSIS, matrix_, nullptr, nullptr),
              XIANZHIKRYLOV_SUCCESS, "analysis");

        if (profileEnabled()) {
            std::fprintf(stderr, "OPM_XK_ANALYSIS rows=%zu blocks=%zu seconds=%.9f\n",
                         A_.N(), columns_.size(), norneElapsedSeconds(start));
        }
    }

    void updateNumericOnly()
    {
        const auto start = NorneClock::now();
        std::size_t cursor = 0;
        for (auto row = A_.begin(); row != A_.end(); ++row) {
            for (auto entry = row->begin(); entry != row->end(); ++entry) {
                for (int i = 0; i < BlockDim; ++i) {
                    for (int j = 0; j < BlockDim; ++j) {
                        values_[cursor++] = (*entry)[i][j];
                    }
                }
            }
        }
        if (cursor != values_.size()) {
            throw std::logic_error("XianzhiKrylov matrix value count mismatch");
        }
        const double packSeconds = norneElapsedSeconds(start);

        auto stage = NorneClock::now();
        dvalues_.upload(values_);
        const double uploadSeconds = norneElapsedSeconds(stage);

        stage = NorneClock::now();
        check(xianzhikrylov_execute(handle_, XIANZHIKRYLOV_PHASE_NUMERIC_PREPARE, matrix_, nullptr, nullptr),
              XIANZHIKRYLOV_SUCCESS, "numeric_prepare");
        const double numericSeconds = norneElapsedSeconds(stage);

        if (profileEnabled()) {
            std::fprintf(
                stderr,
                "OPM_XK_PREPARE rows=%zu blocks=%zu pack_s=%.9f values_h2d_s=%.9f numeric_s=%.9f\n",
                A_.N(), columns_.size(), packSeconds, uploadSeconds, numericSeconds);
        }
    }

    void check(xianzhikrylov_status_t status, xianzhikrylov_status_t expected, const char* phase)
    {
        if (status == expected) {
            return;
        }
        if (handle_) {
            xianzhikrylov_print_info(handle_);
        }
        throw std::runtime_error(std::string("OPM XianzhiKrylov ") + phase +
                                 " failed with status " + std::to_string(status));
    }

    void destroyDescriptors() noexcept
    {
        if (bv_) {
            xianzhikrylov_vector_destroy(bv_);
            bv_ = nullptr;
        }
        if (xv_) {
            xianzhikrylov_vector_destroy(xv_);
            xv_ = nullptr;
        }
        if (matrix_) {
            xianzhikrylov_matrix_destroy(matrix_);
            matrix_ = nullptr;
        }
    }

    void release() noexcept
    {
        if (handle_) {
            xianzhikrylov_destroy(handle_);
            handle_ = nullptr;
        }
        destroyDescriptors();
        doffsets_.reset();
        dcolumns_.reset();
        ddiagonal_.reset();
        dvalues_.reset();
        drhs_.reset();
        dx_.reset();
    }

    const Matrix& A_;
    double tolerance_;
    std::vector<uint64_t> offsets_, columns_;
    std::vector<double> values_, rhs_, result_;
    NorneDeviceBuffer<uint64_t> doffsets_, dcolumns_, ddiagonal_;
    NorneDeviceBuffer<double> dvalues_, drhs_, dx_;
    xianzhikrylov_handle_t handle_ = nullptr;
    xianzhikrylov_matrix_t matrix_ = nullptr;
    xianzhikrylov_vector_t bv_ = nullptr;
    xianzhikrylov_vector_t xv_ = nullptr;
    xianzhikrylov_status_t (*statsFn_)(xianzhikrylov_handle_t, uint32_t*, double*, double*) = nullptr;
};

} // namespace Opm::detail

#endif
