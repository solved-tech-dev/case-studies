#include "linglongSession.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

using namespace PHSPACE;
namespace
{
void Require(bool condition, const char* reason)
{
    if (!condition) throw std::runtime_error(reason);
}
AMGBenchmarkCSR Matrix(int n, double shift, bool nonsymmetric)
{
    AMGBenchmarkCSR a;
    a.rows.push_back(0);
    for (int row = 0; row < n; ++row)
    {
        if (row) { a.columns.push_back(row-1); a.values.push_back(nonsymmetric ? -1.3 : -1.0); }
        a.columns.push_back(row); a.values.push_back(3.0+shift);
        if (row+1 < n) { a.columns.push_back(row+1); a.values.push_back(nonsymmetric ? -0.7 : -1.0); }
        a.rows.push_back(static_cast<int>(a.values.size()));
    }
    return a;
}
std::vector<double> Product(const AMGBenchmarkCSR& a, const std::vector<double>& x)
{
    std::vector<double> b(x.size(), 0);
    for (std::size_t i = 0; i < x.size(); ++i)
        for (int k = a.rows[i]; k < a.rows[i+1]; ++k) b[i] += a.values[k]*x[a.columns[k]];
    return b;
}
void Verify(LinglongSession& session, const AMGBenchmarkCSR& a, bool exactGuess, bool zero)
{
    const int n = static_cast<int>(a.rows.size()-1);
    std::vector<double> expected(n), x(n, 0);
    for (int i = 0; i < n; ++i) expected[i] = zero ? 0 : 0.5+std::sin(0.2*i);
    const std::vector<double> b = Product(a, expected);
    if (exactGuess) x = expected;
    if (zero && !exactGuess) std::fill(x.begin(), x.end(), 0.75);
    const LinglongSessionReport r = session.Solve(a, b.data(), x.data());
    if (!r.converged)
        std::cerr << "nonconverged call=" << r.call << " initial=" << r.initialL2
                  << " final=" << r.finalL2 << " sdk_relative=" << r.sdkRelativeResidual << '\n';
    // With rel_ini, roundoff-sized initial residuals may never improve by the
    // requested relative factor. Retain the SDK flag; verify the known exact
    // guess by independent absolute residual and solution error instead.
    Require(r.converged || (exactGuess && r.initialL2 < 1e-12 && r.finalL2 < 1e-12),
            "SDK reported nonconvergence for a nontrivial initial residual");
    Require(r.finalL2 <= 1e-8*std::max(1.0, r.rhsL2), "independent residual check failed");
    for (int i = 0; i < n; ++i) Require(std::abs(x[i]-expected[i]) < 1e-7, "known solution check failed");
    std::cout << "PASS solve=" << r.call << " rows=" << n << " exact_guess=" << exactGuess
              << " zero_rhs=" << zero << " sdk_converged=" << r.converged
              << " true_residual=" << r.finalL2 << '\n';
}
void Audit(unsigned calls, bool denied)
{
    const char* path = std::getenv("LINGLONG_AUDIT_FILE");
    Require(path && *path, "audit file required for actual API count verification");
    std::ifstream input(path);
    Require(static_cast<bool>(input), "cannot read audit file");
    std::map<std::string, unsigned> counts;
    std::string event, handle, created;
    int status;
    while (input >> event >> handle >> status)
    {
        ++counts[event];
        if (event == "create") created = handle;
        if (event == "setup" || event == "solve" || event == "query" || event == "reset" || event == "destroy")
            Require(!created.empty() && handle == created, "SDK handle changed during reuse");
        Require(status == 0 || (denied && event == "query" && status == 6), "unexpected real SDK error");
    }
    Require(counts["create"] == 1 && counts["destroy"] == 1, "expected exactly one create and destroy");
    Require(counts["setup"] == calls && counts["solve"] == calls && counts["query"] == calls, "setup/solve/query count mismatch");
    Require(counts["reset"] == calls-1, "reset count mismatch");
    Require(counts["license_check"] == 0, "adapter repeated an explicit authorization diagnostic");
    std::cout << "PASS lifecycle create=1 destroy=1 setup=" << calls << " solve=" << calls
              << " query=" << calls << " reset=" << calls-1 << " same_handle=1\n";
}
}

int main(int argc, char** argv)
{
    try
    {
        Require(argc == 2 || argc == 3, "usage: linglong_session_test config.json [--deny-query]");
        const bool denied = argc == 3 && std::string(argv[2]) == "--deny-query";
        unsigned calls = 0;
        {
            LinglongSession session(argv[1], 0);
            std::cout << "GPU " << session.DeviceDescription() << '\n';
            if (denied)
            {
                const AMGBenchmarkCSR a = Matrix(9, 0.2, true);
                std::vector<double> x(9, 0.375), expected(9, 1.0);
                const auto b = Product(a, expected);
                const auto unchanged = x;
                bool rejected = false;
                try { session.Solve(a, b.data(), x.data()); }
                catch (const std::runtime_error& error)
                { rejected = std::string(error.what()).find("linglongamg_get_solve_info") != std::string::npos; }
                Require(rejected && x == unchanged, "authorization denial published a candidate solution");
                rejected = false;
                try { session.Solve(a, b.data(), x.data()); }
                catch (const std::runtime_error&) { rejected = true; }
                Require(rejected && x == unchanged, "failed session was recreated or reused");
                calls = 1;
                std::cout << "PASS denied query leaves host x unchanged and session unusable\n";
            }
            else
            {
                // Includes duplicate faces and one physical boundary (ghost cell).
                const int left[] = {0, 0, 1, 0}, right[] = {1, 1, 2, 3};
                const double diagonal[] = {5, 6, 5}, upper[] = {-0.3, -0.4, -0.7, 0}, lower[] = {-0.6, -0.7, -1.3, 0};
                Verify(session, AssembleAMGBenchmarkCSR(3, 4, left, right, diagonal, upper, lower), false, false);
                ++calls;
                // Coefficients, topology and dimensions all change on one handle.
                for (int n : {3, 19, 128, 7})
                { Verify(session, Matrix(n, 0.1*calls, true), false, false); ++calls; }
                std::this_thread::sleep_for(std::chrono::milliseconds(1200)); // Exercise reuse after a delay.
                Verify(session, Matrix(64, 0, false), false, false); ++calls;
                Verify(session, Matrix(17, 0.2, false), false, true); ++calls;
                Verify(session, Matrix(17, 0.2, false), true, true); ++calls;
                Verify(session, Matrix(9, 0.5, true), true, false); ++calls;
                // Bad host input must not consume or replace the SDK session.
                auto bad = Matrix(3, 0, false); bad.rows[1] = -1;
                double b[] = {1, 2, 3}, x[] = {0, 0, 0};
                bool rejected = false;
                try { session.Solve(bad, b, x); } catch (const std::invalid_argument&) { rejected = true; }
                Require(rejected, "invalid CSR was accepted");
                Verify(session, Matrix(5, 0.3, true), false, false); ++calls;
            }
            session.Close();
            session.Close(); // Idempotent, without a second SDK destroy.
        }
        Audit(calls, denied);
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "FAIL " << error.what() << '\n';
        return 1;
    }
}
