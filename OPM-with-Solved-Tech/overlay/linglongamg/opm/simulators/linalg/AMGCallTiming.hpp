#ifndef OPM_AMG_CALL_TIMING_HPP
#define OPM_AMG_CALL_TIMING_HPP
// Only SDK calls and their completion synchronization; excludes transfers/CPR.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <unistd.h>
namespace Opm {
struct AMGCallTiming {
    const char* backend;
    std::uint64_t setup_calls = 0, solve_calls = 0, resetup_calls = 0;
    std::uint64_t setup_failures = 0, solve_failures = 0;
    std::uint64_t max_rows = 0, max_nnz = 0;
    double setup_seconds = 0, solve_seconds = 0;
    std::uint64_t info_calls = 0, info_failures = 0;
    int last_converged = 0, last_iterations = 0;
    double last_initial_residual = 0, last_final_residual = 0, last_relative_residual = 0;
    static AMGCallTiming& linglong() { static AMGCallTiming t{"linglongamg"}; return t; }
    static AMGCallTiming& hypreGpu() { static AMGCallTiming t{"hypre_gpu"}; return t; }
    static AMGCallTiming& hypreCpu() { static AMGCallTiming t{"hypre_cpu"}; return t; }
    static AMGCallTiming& amgx() { static AMGCallTiming t{"amgx"}; return t; }
    template<class Sync, class Operation>
    void measure(bool setup, std::uint64_t rows, std::uint64_t nnz,
                 Sync&& synchronize, Operation&& operation, bool resetup = false)
    {
        max_rows = std::max(max_rows, rows);
        max_nnz = std::max(max_nnz, nnz);
        synchronize(); // Drain earlier uploads/work outside the interval.
        const auto start = std::chrono::steady_clock::now();
        ++(setup ? setup_calls : solve_calls);
        if (resetup) ++resetup_calls;
        auto& seconds = setup ? setup_seconds : solve_seconds;
        try { operation(); synchronize(); }
        catch (...) {
            seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
            ++(setup ? setup_failures : solve_failures);
            throw;
        }
        seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    }
    static void number(std::ostream& out, double value)
    {
        if (std::isfinite(value)) out << value; else out << "null";
    }
    void report(bool final = false) const noexcept
    {
        if (!setup_calls && !solve_calls) return;
        try {
            const char* rank = std::getenv("OMPI_COMM_WORLD_RANK");
            if (!rank) rank = std::getenv("PMI_RANK");
            if (!rank) rank = std::getenv("SLURM_PROCID");
            std::ostringstream out;
            out << std::setprecision(17)
                << "{\"schema\":\"opm_amg_sdk_timing_v1\",\"backend\":\"" << backend
                << "\",\"pid\":" << getpid() << ",\"rank\":" << (rank ? std::atoi(rank) : 0)
                << ",\"final\":" << (final ? "true" : "false")
                << ",\"scope\":\"sdk_call_post_sync_excludes_transfers_cpr\""
                << ",\"setup_calls\":" << setup_calls << ",\"solve_calls\":" << solve_calls
                << ",\"resetup_calls\":" << resetup_calls
                << ",\"setup_seconds\":" << setup_seconds << ",\"solve_seconds\":" << solve_seconds
                << ",\"setup_failures\":" << setup_failures << ",\"solve_failures\":" << solve_failures
                << ",\"max_local_rows\":" << max_rows << ",\"max_local_nnz\":" << max_nnz
                << ",\"sdk_info_calls\":" << info_calls << ",\"sdk_info_failures\":" << info_failures;
            if (info_calls) {
                out << ",\"sdk_last_converged\":" << last_converged
                    << ",\"sdk_last_iterations\":" << last_iterations
                    << ",\"sdk_last_initial_residual\":"; number(out, last_initial_residual);
                out << ",\"sdk_last_final_residual\":"; number(out, last_final_residual);
                out << ",\"sdk_last_relative_residual\":"; number(out, last_relative_residual);
            }
            out << "}";
            std::cerr << "OPM_AMG_TIMING " << out.str() << '\n';
            if (const auto* path = std::getenv("OPM_AMG_TIMING_FILE")) {
                std::ofstream file(std::string(path) + "." + std::to_string(getpid()) + ".jsonl", std::ios::app);
                if (file) file << out.str() << '\n';
                else std::cerr << "OPM_AMG_TIMING_FILE_ERROR pid=" << getpid() << '\n';
            }
        } catch (...) {} // Reporting must not alter simulation completion.
    }
    ~AMGCallTiming() { report(true); }
};
} // namespace Opm
#endif
