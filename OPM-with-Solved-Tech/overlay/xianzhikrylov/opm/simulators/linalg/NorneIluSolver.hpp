#pragma once

#include "NorneIluDevice.hpp"
#include "NorneWellDevice.hpp"
#include <opm/simulators/linalg/PreconditionerWithUpdate.hpp>

#include <chrono>
#include <cstdio>
#include <stdexcept>
#include <vector>

namespace Opm::detail {

template <class Matrix, class Vector>
class NorneIluSolver final : public Dune::InverseOperator<Vector, Vector>,
                             public Dune::PreconditionerWithUpdate<Vector, Vector>
{
    using Operator = WellModelMatrixAdapter<Matrix, Vector, Vector>;
    using Clock = std::chrono::steady_clock;

public:
    NorneIluSolver(const Operator& op, const NorneWellSource& wells, const PropertyTree& config)
        : op_(op), wells_(wells), tolerance_(config.get<double>("tol")),
          maxIterations_(config.get<int>("maxiter"))
    {
        if (config.get<std::string>("solver") != "bicgstab" ||
            config.get<std::string>("preconditioner.type") != "paroverilu0" ||
            config.get<double>("preconditioner.relaxation") != 1.0 ||
            !(tolerance_ > 0.0 && tolerance_ < 1.0) || maxIterations_ < 1) {
            throw std::invalid_argument("Norne GPU ILU requires BiCGStab with unrelaxed ILU(0)");
        }
        update();
    }

    Dune::SolverCategory::Category category() const override
    {
        return Dune::SolverCategory::sequential;
    }

    bool hasPerfectUpdate() const override { return true; }
    void pre(Vector&, Vector&) override {}
    void post(Vector&) override {}

    void apply(Vector&, const Vector&) override
    {
        throw std::logic_error("Norne ILU preconditioning occurs inside the GPU BiCGStab solve");
    }

    void update() override
    {
        const auto start = Clock::now();
        const auto& matrix = op_.getmat();
        std::vector<uint64_t> rows(matrix.N() + 1);
        std::vector<uint64_t> columns;
        columns.reserve(matrix.nonzeroes());
        values_.resize(9 * matrix.nonzeroes());
        for (auto row = matrix.begin(); row != matrix.end(); ++row) {
            for (auto entry = row->begin(); entry != row->end(); ++entry) {
                const std::size_t block = columns.size();
                columns.push_back(entry.index());
                for (int i = 0; i < 3; ++i) {
                    for (int j = 0; j < 3; ++j) {
                        values_[9 * block + 3 * i + j] = (*entry)[i][j];
                    }
                }
            }
            rows[row.index() + 1] = columns.size();
        }
        const double pack = seconds(start);
        const auto wellStart = Clock::now();
        deviceWells_ = wells_.deviceWells(3 * matrix.N());
        const double wellPrepare = seconds(wellStart);
        const auto topologyStart = Clock::now();
        const bool rebuild = !device_ || row_ != rows || column_ != columns;
        if (rebuild) {
            device_.reset();
            device_ = createIluDevice(rows, columns, tolerance_, maxIterations_);
            row_ = std::move(rows);
            column_ = std::move(columns);
        }
        const double topology = seconds(topologyStart);
        const auto timing = prepareIluDevice(
            *device_, values_, deviceWells_.get(),
            +[](void* well, const double* x, double* y) {
                static_cast<NorneWellDevice*>(well)->apply(const_cast<double*>(x), y);
                return 0;
            });
        std::fprintf(stderr,
                     "NORNE_ILU_PREPARE setup_s=%.9f pack_s=%.9f well_prepare_s=%.9f "
                     "topology_s=%.9f values_h2d_s=%.9f ilu_factor_s=%.9f rebuild=%d\n",
                     seconds(start), pack, wellPrepare, topology, timing.valuesUpload, timing.factor, int(rebuild));
    }

    void apply(Vector& x, Vector& b, double tolerance, Dune::InverseOperatorResult& result) override
    {
        if (tolerance != tolerance_) {
            throw std::logic_error("Unexpected Norne BiCGStab tolerance override");
        }
        apply(x, b, result);
    }

    void apply(Vector& x, Vector& b, Dune::InverseOperatorResult& result) override
    {
        rhs_.resize(3 * b.size());
        solution_.resize(rhs_.size());
        for (std::size_t i = 0; i < b.size(); ++i) {
            for (int eq = 0; eq < 3; ++eq) {
                rhs_[3 * i + eq] = b[i][eq];
            }
        }

        const auto info = solveIluDevice(*device_, rhs_, solution_);
        for (std::size_t i = 0; i < x.size(); ++i) {
            for (int eq = 0; eq < 3; ++eq) {
                x[i][eq] = solution_[3 * i + eq];
            }
        }

        result.clear();
        result.iterations = info.iterations;
        result.converged = info.converged;
        result.reduction = info.relativeResidual;
        result.elapsed = info.solve;
        std::fprintf(stderr,
                     "NORNE_ILU_SOLVE solve_s=%.9f rhs_h2d_s=%.9f x_d2h_s=%.9f "
                     "reservoir_spmv_s=%.9f well_s=%.9f ilu_apply_s=%.9f matvecs=%u "
                     "iterations=%u converged=%d relative=%.17g\n",
                     info.solve, info.rhsUpload, info.solutionDownload, info.fineSpmv, info.well,
                     info.iluApply, info.matvecs, info.iterations, int(info.converged), info.relativeResidual);
    }

private:
    static double seconds(Clock::time_point start)
    {
        return std::chrono::duration<double>(Clock::now() - start).count();
    }

    const Operator& op_;
    const NorneWellSource& wells_;
    double tolerance_;
    int maxIterations_;
    std::vector<uint64_t> row_, column_;
    std::vector<double> values_, rhs_, solution_;
    std::unique_ptr<NorneWellDevice> deviceWells_;
    IluDeviceHandle device_;
};

} // namespace Opm::detail
