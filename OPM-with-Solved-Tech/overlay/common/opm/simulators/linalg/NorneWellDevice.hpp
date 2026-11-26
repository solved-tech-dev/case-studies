#pragma once

#include <opm/simulators/linalg/WellOperators.hpp>
#include <opm/simulators/linalg/gpubridge/cuda/cuWellContributions.hpp>
#include <opm/simulators/linalg/gpubridge/MultisegmentWellContribution.hpp>

#include <memory>
#include <stdexcept>

namespace Opm::detail {

class NorneWellDevice final : public WellContributionsCuda<double>
{
public:
    bool requiresHostWellSolve() const { return this->num_ms_wells != 0; }
};

struct NorneWellSource
{
    virtual ~NorneWellSource() = default;
    virtual std::unique_ptr<NorneWellDevice> deviceWells(unsigned n) const = 0;
};

template <class WellModel, class Vector>
class NorneWellOperator final : public WellModelAsLinearOperator<WellModel, Vector, Vector>,
                                public NorneWellSource
{
    using Parent = WellModelAsLinearOperator<WellModel, Vector, Vector>;

public:
    using Parent::Parent;

    std::unique_ptr<NorneWellDevice> deviceWells(unsigned n) const override
    {
        auto data = std::make_unique<NorneWellDevice>();
        data->setVectorSize(n);
        data->setCudaStream(nullptr);
        this->wellMod_.getWellContributions(*data);
        if (data->requiresHostWellSolve()) {
            throw std::runtime_error("OPM GPU device well path requires standard wells");
        }
        return data;
    }
};

} // namespace Opm::detail
