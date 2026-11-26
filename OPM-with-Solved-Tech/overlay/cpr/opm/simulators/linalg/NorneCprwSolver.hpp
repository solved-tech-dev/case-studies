#pragma once

#include "NorneCprwDevice.hpp"
#include "NorneWellDevice.hpp"
#include <opm/simulators/linalg/PressureBhpTransferPolicy.hpp>
#include <opm/simulators/linalg/ParallelOverlappingILU0.hpp>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace Opm::detail {

template <class Matrix, class Vector>
class NorneCprwSolver final : public Dune::InverseOperator<Vector, Vector>,
                              public Dune::PreconditionerWithUpdate<Vector, Vector>
{
    using Operator = WellModelMatrixAdapter<Matrix, Vector, Vector>;
    using Comm = Dune::Amg::SequentialInformation;
    using Policy = PressureBhpTransferPolicy<Operator, Comm, double, false>;
    using Clock = std::chrono::steady_clock;

public:
    NorneCprwSolver(const Operator& op,
                    const NorneWellSource& wells,
                    const PropertyTree& options,
                    const std::function<Vector()>& weights,
                    int pressureIndex)
        : op_(op), wells_(wells), settings_(options.get_child("preconditioner")), calculator_(weights),
          ip_(pressureIndex), maximum_(options.get<int>("maxiter", 20)), tol_(options.get<double>("tol", 0.005))
    {
        if (options.get<std::string>("solver") != "bicgstab" ||
            settings_.get<std::string>("type") != "cprw" ||
            settings_.get<std::string>("weight_type") != "trueimpes" ||
            settings_.get<std::string>("finesmoother.type") != "paroverilu0" ||
            settings_.get<double>("finesmoother.relaxation") != 1.0 ||
            settings_.get<int>("pre_smooth") != 0 || settings_.get<int>("post_smooth") != 1 ||
            !settings_.get<bool>("add_wells") || settings_.get<bool>("use_well_weights")) {
            throw std::runtime_error("The device CPRW route requires the Norne CPRW configuration");
        }
        if (const char* p = std::getenv("NORNE_AMG_LIBRARY"); p && *p) {
            amgPath_ = p;
        } else if (const char* p2 = std::getenv("LINGLONGAMG_OPM_LIBRARY"); p2 && *p2) {
            amgPath_ = p2;
        } else if (const char* root = std::getenv("LINGLONGAMG_ROOT"); root && *root) {
            amgPath_ = std::string(root) + "/lib/liblinglongamg.so";
        }

        if (const char* c = std::getenv("NORNE_AMG_CONFIG"); c && *c) {
            amgCfgPath_ = c;
        } else if (const char* c2 = std::getenv("LINGLONGAMG_OPM_CONFIG"); c2 && *c2) {
            amgCfgPath_ = c2;
        }
        if (amgPath_.empty() || amgCfgPath_.empty()) {
            throw std::runtime_error("Missing LinglongAMG library or configuration for CPR route");
        }
        update();
    }

    Dune::SolverCategory::Category category() const override { return Dune::SolverCategory::sequential; }
    bool hasPerfectUpdate() const override { return true; }
    void pre(Vector&, Vector&) override {}
    void post(Vector&) override {}

    void apply(Vector&, const Vector&) override
    {
        throw std::logic_error("Norne CPRW uses XianzhiKrylov as the outer inverse operator");
    }

    void update() override
    {
        const auto start = Clock::now();
        const auto& a = op_.getmat();
        const auto weights = calculator_();
        Policy pressure(comm_, weights, settings_, ip_);
        pressure.createCoarseLevelSystem(op_);
        const auto& ap = pressure.getCoarseLevelOperator()->getmat();
        if (ap.N() != a.N() + op_.getNumberOfExtraEquations()) {
            throw std::logic_error("Norne CPRW pressure build lost well equations");
        }
        if (ap.nonzeroes() > std::numeric_limits<int>::max() || a.N() > std::numeric_limits<int>::max() / 3) {
            throw std::overflow_error("Norne CPRW dimensions exceed the device route");
        }
        const double pressureBuild = seconds(start);
        const auto packStart = Clock::now();
        std::vector<uint64_t> rows(a.N() + 1), cols;
        std::vector<int> pRows(ap.N() + 1), pCols;
        cols.reserve(a.nonzeroes());
        pCols.reserve(ap.nonzeroes());
        fine_.resize(9 * a.nonzeroes());
        weight_.resize(3 * a.N());
        pressure_.resize(ap.nonzeroes());
        for (auto row = a.begin(); row != a.end(); ++row) {
            const auto i = row.index();
            for (int eq = 0; eq < 3; ++eq) {
                weight_[3 * i + eq] = weights[i][eq];
            }
            for (auto entry = row->begin(); entry != row->end(); ++entry) {
                const auto k = cols.size();
                cols.push_back(entry.index());
                for (int r = 0; r < 3; ++r) {
                    for (int c = 0; c < 3; ++c) {
                        fine_[9 * k + 3 * r + c] = (*entry)[r][c];
                    }
                }
            }
            rows[i + 1] = cols.size();
        }
        for (auto row = ap.begin(); row != ap.end(); ++row) {
            for (auto entry = row->begin(); entry != row->end(); ++entry) {
                pressure_[pCols.size()] = (*entry)[0][0];
                pCols.push_back(entry.index());
            }
            pRows[row.index() + 1] = pCols.size();
        }
        const double pack = seconds(packStart);
        const auto wellStart = Clock::now();
        deviceWell_ = wells_.deviceWells(3 * a.N());
        const double wellPrepare = seconds(wellStart);
        const auto topologyStart = Clock::now();
        const bool rebuild = !device_ || row_ != rows || col_ != cols || pRow_ != pRows || pCol_ != pCols;
        if (rebuild) {
            device_.reset();
            device_ = createCprwDevice({rows, cols, pRows, pCols, ip_},
                                       amgPath_.c_str(),
                                       amgCfgPath_.c_str(),
                                       tol_,
                                       maximum_);
            row_ = std::move(rows);
            col_ = std::move(cols);
            pRow_ = std::move(pRows);
            pCol_ = std::move(pCols);
        }
        const double topology = seconds(topologyStart);
        prepareCprwDevice(
            *device_, {fine_, weight_, pressure_}, deviceWell_.get(),
            +[](void* w, const double* x, double* y) {
                static_cast<NorneWellDevice*>(w)->apply(const_cast<double*>(x), y);
                return 0;
            });
        std::fprintf(stderr,
                     "NORNE_CPRW_SETUP rows=%zu blocks=%zu pressure_rows=%zu pressure_blocks=%zu "
                     "wells=%d pressure_build_s=%.9f pack_s=%.9f well_prepare_s=%.9f topology_s=%.9f rebuild=%d\n",
                     a.N(), a.nonzeroes(), ap.N(), ap.nonzeroes(), op_.getNumberOfExtraEquations(),
                     pressureBuild, pack, wellPrepare, topology, int(rebuild));
    }

    void apply(Vector& x, Vector& b, double tolerance, Dune::InverseOperatorResult& result) override
    {
        if (tolerance != tol_) {
            throw std::logic_error("Unexpected Norne outer tolerance override");
        }
        apply(x, b, result);
    }

    void apply(Vector& x, Vector& b, Dune::InverseOperatorResult& result) override
    {
        const bool audit = std::getenv("NORNE_GPU_CPRW_AUDIT") != nullptr;
        rhs_.resize(b.size() * 3);
        solution_.resize(rhs_.size());
        for (std::size_t i = 0; i != b.size(); ++i) {
            for (int eq = 0; eq < 3; ++eq) {
                rhs_[3 * i + eq] = b[i][eq];
            }
        }
        if (audit && !audited_) {
            Vector probe(b.size()), expected(b.size());
            for (std::size_t i = 0; i != b.size(); ++i) {
                for (int eq = 0; eq < 3; ++eq) {
                    solution_[3 * i + eq] = probe[i][eq] = double((i * 17 + eq * 7) % 113) / 113.0 - 0.5;
                }
            }
            std::vector<double> check(solution_.size());
            auditCprwOperator(*device_, solution_, check);
            op_.apply(probe, expected);
            const double reference = expected.two_norm();
            for (std::size_t i = 0; i != b.size(); ++i) {
                for (int eq = 0; eq < 3; ++eq) {
                    expected[i][eq] -= check[3 * i + eq];
                }
            }
            const double error = expected.two_norm() / std::max(reference, 1e-30);
            if (!std::isfinite(error) || error > 1e-10) {
                throw std::runtime_error("GPU well-inclusive operator failed independent comparison");
            }
            std::fprintf(stderr, "NORNE_CPRW_OPERATOR_AUDIT relative=%.17g\n", error);
            auditCorrection(probe, check);
            audited_ = true;
        }
        const auto info = solveCprwDevice(*device_, rhs_, solution_);
        for (std::size_t i = 0; i != x.size(); ++i) {
            for (int eq = 0; eq < 3; ++eq) {
                x[i][eq] = solution_[3 * i + eq];
            }
        }
        result.clear();
        result.iterations = info.iterations;
        result.converged = info.converged;
        result.reduction = info.residual;
        result.elapsed = info.seconds;
        std::fprintf(stderr, "NORNE_CPRW_GPU outer_s=%.9f iterations=%u converged=%d relative=%.17g\n",
                     info.seconds, info.iterations, int(info.converged), info.residual);
        if (audit) {
            Vector residual(b);
            const double norm = b.two_norm();
            op_.applyscaleadd(-1.0, x, residual);
            const double independent = residual.two_norm() / std::max(norm, 1e-30);
            std::fprintf(stderr, "NORNE_CPRW_TRUE_RESIDUAL iterations=%u converged=%d sdk=%.17g independent=%.17g\n",
                         info.iterations, int(info.converged), info.residual, independent);
            if (!std::isfinite(independent) || std::abs(independent - info.residual) > 1e-9 ||
                (info.converged && independent > tol_ + 1e-10)) {
                throw std::runtime_error("Norne original-system residual check failed");
            }
        }
    }

private:
    void auditCorrection(const Vector& probe, std::span<double> correction)
    {
        const auto weights = calculator_();
        Policy oracle(comm_, weights, settings_, ip_);
        oracle.createCoarseLevelSystem(op_);
        oracle.moveToCoarseLevel(probe);
        const auto& coarseRhs = oracle.getCoarseLevelRhs();
        auto coarseSolution = coarseRhs;
        std::vector<double> restricted(coarseRhs.size()), pressure(coarseRhs.size());
        auditCprwCorrection(*device_, solution_, restricted, pressure, correction);
        auto restrictionError = coarseRhs;
        for (std::size_t i = 0; i < coarseRhs.size(); ++i) {
            restrictionError[i][0] -= restricted[i];
            coarseSolution[i][0] = pressure[i];
        }
        auto residual = coarseRhs;
        oracle.getCoarseLevelOperator()->applyscaleadd(-1.0, coarseSolution, residual);
        const double rhsNorm = std::max(coarseRhs.two_norm(), 1e-30);
        const double restriction = restrictionError.two_norm() / rhsNorm;
        const double pressureResidual = residual.two_norm() / rhsNorm;

        Vector injection(probe.size()), fineRhs(probe), expected(probe.size());
        injection = 0.0;
        expected = 0.0;
        for (std::size_t i = 0; i < injection.size(); ++i) {
            injection[i][ip_] = coarseSolution[i][0];
        }
        op_.applyscaleadd(-1.0, injection, fineRhs);
        ParallelOverlappingILU0<Matrix, Vector, Vector, Comm> ilu(op_.getmat(), 1.0, MILU_VARIANT::ILU);
        ilu.pre(expected, fineRhs);
        ilu.apply(expected, fineRhs);
        ilu.post(expected);
        expected += injection;
        const double correctionNorm = std::max(expected.two_norm(), 1e-30);
        for (std::size_t i = 0; i < expected.size(); ++i) {
            for (int eq = 0; eq < 3; ++eq) {
                expected[i][eq] -= correction[3 * i + eq];
            }
        }
        const double correctionError = expected.two_norm() / correctionNorm;
        std::fprintf(stderr,
                     "NORNE_CPRW_CORRECTION_AUDIT restriction=%.17g pressure_residual=%.17g correction=%.17g\n",
                     restriction, pressureResidual, correctionError);
        if (!std::isfinite(restriction) || !std::isfinite(pressureResidual) || !std::isfinite(correctionError) ||
            restriction > 1e-12 || pressureResidual > 1.01e-6 || correctionError > 1e-9) {
            throw std::runtime_error("Norne GPU CPRW correction failed independent comparison");
        }
    }

    static double seconds(Clock::time_point t)
    {
        return std::chrono::duration<double>(Clock::now() - t).count();
    }

    const Operator& op_;
    const NorneWellSource& wells_;
    const PropertyTree settings_;
    const std::function<Vector()> calculator_;
    const int ip_, maximum_;
    const double tol_;
    const Comm comm_{};
    std::string amgPath_, amgCfgPath_;
    std::vector<uint64_t> row_, col_;
    std::vector<int> pRow_, pCol_;
    std::vector<double> fine_, weight_, pressure_, rhs_, solution_;
    std::unique_ptr<NorneWellDevice> deviceWell_;
    CprwHandle device_;
    bool audited_ = false;
};

} // namespace Opm::detail
