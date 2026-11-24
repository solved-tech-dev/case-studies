// SPDX-License-Identifier: GPL-3.0-or-later
#include "readRawCase.H"
#include "config_snapshot.hpp"
#include "XianzhiKrylovAPI.h"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>

using namespace Foam;

int main(int argc, char** argv)
{
    try {
        UPstream::initNull();
        argList::addBoolOption("check-case", "Validate input without creating GPU state");
        argList::addOption("expected-cells", "label", "Cell count from GS partition offsets");
        #include "setRootCase.H"
        #include "createTime.H"

        auto mesh = raw_resident_io::read_topology(runTime);
        auto parsed = raw_resident_io::read_input(mesh, runTime);
        const auto& control = runTime.controlDict();
        if (control.getOrDefault<word>("stopAt", "endTime") != "endTime"
            || control.getOrDefault<Switch>("adjustTimeStep", false))
            throw std::runtime_error("resident route requires fixed-step stopAt endTime");
        const double delta = runTime.deltaTValue();
        if (!std::isfinite(delta) || delta <= 0.)
            throw std::runtime_error("invalid fixed deltaT");
        const double span = (runTime.endTime().value() - runTime.value()) / delta;
        if (!std::isfinite(span) || span <= 0.
            || span > std::numeric_limits<int>::max()
            || std::abs(span - std::round(span)) > 1e-9)
            throw std::runtime_error("resident route requires a positive integral outer-step count");
        const int steps = static_cast<int>(std::round(span));
        if (const auto* residual = parsed.controls.solution.subDict("SIMPLE").findDict("residualControl"))
            if (!residual->empty())
                throw std::runtime_error("SIMPLE residualControl is not implemented");
        if (args.found("expected-cells")
            && args.get<label>("expected-cells") != mesh.nCells())
            throw std::runtime_error("GS partition offsets do not match the input cell count");
        if (!resident_configuration_snapshots())
            throw std::runtime_error("launch through scripts/run-resident.sh to provide configuration snapshots");

        auto data = raw_resident_io::marshal(mesh, parsed);
        parsed = raw_resident_io::RawParsedCase{};
        std::cout << "RAW_RESIDENT_INPUT cells=" << data.input.cells
                  << " faces=" << data.input.owner.size() << " steps=" << steps << std::endl;
        if (args.found("check-case")) {
            std::cout << "RAW_RESIDENT_INPUT_VALIDATED gpu_initialized=0" << std::endl;
            return 0;
        }

        xk_fv_resident::ResidentContext context(data.input, data.geometry);
        data = raw_resident_io::RawCase{};
        const char* output = std::getenv("RESIDENT_STATE_DUMPS");
        auto dump = [&](const std::string& state) {
            if (!output) return;
            const auto directory = std::filesystem::path(output) / state;
            std::filesystem::create_directories(directory);
            for (const auto* name : {"U", "p", "phi", "k", "omega", "nut",
                                    "U.boundary", "p.boundary", "k.boundary",
                                    "omega.boundary", "nut.boundary"}) {
                const auto values = context.output(name);
                for (double value : values)
                    if (!std::isfinite(value)) throw std::runtime_error("nonfinite output field");
                std::ofstream stream(directory / (std::string(name) + ".f64"), std::ios::binary);
                stream.write(reinterpret_cast<const char*>(values.data()), values.size() * sizeof(double));
                if (!stream) throw std::runtime_error("field output failed");
            }
        };
        dump("initial");
        for (int outer = 0; outer < steps; ++outer) {
            const auto result = context.advance();
            std::cout << "RAW_RESIDENT_STEP outer=" << outer << " stage=" << result.stage
                      << " cuda=" << int(result.device_status)
                      << " solver=" << int(result.solver_status) << std::endl;
            if (result.device_status != cudaSuccess
                || (result.solver_status != XIANZHIKRYLOV_SUCCESS
                    && result.solver_status != XIANZHIKRYLOV_CONVERGED))
                throw std::runtime_error(std::string("resident step failed: ") + result.stage);
            dump("step" + std::to_string(outer));
        }
        std::cout << "RAW_RESIDENT_COMPLETE steps=" << steps
                  << " requested_d2h_bytes=" << context.requested_download_bytes() << std::endl;
    } catch (const std::exception& error) {
        std::cerr << "RAW_RESIDENT_ERROR " << error.what() << std::endl;
        return 1;
    }
    return 0;
}
