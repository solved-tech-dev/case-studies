#ifndef __SPE10_BENCHMARK_GRID_HEADER__
#define __SPE10_BENCHMARK_GRID_HEADER__

#include <cstdint>
#include <vector>

namespace ocp_spe10_benchmark_grid {

struct GridData {
    std::uint32_t nx{0};
    std::uint32_t ny{0};
    std::uint32_t nz{0};
    std::vector<double> dx;
    std::vector<double> dy;
    std::vector<double> dz;
    std::vector<double> kx;
    std::vector<double> ky;
    std::vector<double> kz;
};

inline GridData& Data()
{
    static GridData data;
    return data;
}

inline void Capture(std::uint32_t nx, std::uint32_t ny, std::uint32_t nz,
                    const std::vector<double>& dx,
                    const std::vector<double>& dy,
                    const std::vector<double>& dz,
                    const std::vector<double>& kx,
                    const std::vector<double>& ky,
                    const std::vector<double>& kz)
{
    auto& data = Data();
    data.nx = nx;
    data.ny = ny;
    data.nz = nz;
    data.dx = dx;
    data.dy = dy;
    data.dz = dz;
    data.kx = kx;
    data.ky = ky;
    data.kz = kz;
}

} // namespace ocp_spe10_benchmark_grid

#endif
