// Single-process FrontISTR COO to PanguLU adapter. See NOTICE.md.
#include "gpu_backend.inc"
#include <charconv>

namespace {
class SolverThreads {
    const int previous = omp_get_max_threads();
public:
    SolverThreads() {
        if (const char* setting = std::getenv("FRONTISTR_PANGULU_THREADS")) {
            int count = 0;
            const char* end = setting + std::strlen(setting);
            const auto result = std::from_chars(setting, end, count);
            require(result.ec == std::errc{} && result.ptr == end && count > 0,
                    "FRONTISTR_PANGULU_THREADS must be a positive integer");
            omp_set_num_threads(count);
        }
        std::printf("SOLVER_THREADS entry=%d solver=%d\n", previous, omp_get_max_threads());
    }
    ~SolverThreads() {
        omp_set_num_threads(previous);
        std::printf("SOLVER_THREADS restored=%d\n", omp_get_max_threads());
    }
};

struct Matrix {
    std::vector<int64_t> offsets;
    std::vector<int32_t> columns;
    std::vector<double> values;
};

std::unique_ptr<Context> context;
Matrix previous;

// Convert one-based assembled entries to full CSR. Preserve the accumulation
// order of the benchmark adapter when several entries share a column.
Matrix assemble(int n, int64_t nz, const int32_t* rows, const int32_t* cols,
                const double* values, bool triangular) {
    Matrix full;
    full.offsets.resize(size_t(n) + 1);
    for (int64_t k = 0; k < nz; ++k) {
        require(rows[k] >= 1 && rows[k] <= n && cols[k] >= 1 && cols[k] <= n, "invalid COO index");
        require(std::isfinite(values[k]), "invalid coefficient");
        ++full.offsets[rows[k]];
        if (triangular && rows[k] != cols[k]) ++full.offsets[cols[k]];
    }
    for (int i = 0; i < n; ++i) full.offsets[i+1] += full.offsets[i];
    full.columns.resize(full.offsets.back());
    full.values.resize(full.offsets.back());
    auto next = full.offsets;
    auto append = [&](int row, int col, double value) {
        auto slot = next[row-1]++;
        full.columns[slot] = col-1;
        full.values[slot] = value;
    };
    for (int64_t k = 0; k < nz; ++k) {
        append(rows[k], cols[k], values[k]);
        if (triangular && rows[k] != cols[k]) append(cols[k], rows[k], values[k]);
    }
    int64_t write = 0;
    std::vector<std::pair<int32_t, double>> entries;
    for (int i = 0; i < n; ++i) {
        entries.clear();
        for (auto k = full.offsets[i]; k < full.offsets[i+1]; ++k)
            entries.emplace_back(full.columns[k], full.values[k]);
        std::sort(entries.begin(), entries.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
        full.offsets[i] = write;
        for (auto [col, value] : entries) {
            if (write > full.offsets[i] && col == full.columns[write-1]) {
                full.values[write-1] += value;
            } else {
                full.columns[write] = col;
                full.values[write++] = value;
            }
        }
    }
    full.offsets[n] = write;
    full.columns.resize(write);
    full.values.resize(write);
    return full;
}

double relative_residual(const Matrix& a, const double* b, const double* x, int n) {
    long double error = 0, norm = 0;
    for (int i = 0; i < n; ++i) {
        long double r = -b[i];
        for (auto k = a.offsets[i]; k < a.offsets[i+1]; ++k)
            r += static_cast<long double>(a.values[k]) * x[a.columns[k]];
        error += r*r;
        norm += static_cast<long double>(b[i])*b[i];
    }
    return std::sqrt(error / (norm > 0 ? norm : 1));
}
}

extern "C" int asb_solve_coo_i32(const int64_t* dimension, const int64_t* nz,
                                 const int32_t* row, const int32_t* col, const double* values,
                                 const double* rhs, double* x, const int32_t* triangle) {
    try {
        SolverThreads threads;
        require(*dimension > 0 && *dimension <= INT_MAX && *nz > 0, "unsupported system size");
        const int n = static_cast<int>(*dimension);
        ++call;
        auto a = assemble(n, *nz, row, col, values, *triangle != 0);
        const bool same_structure = context && a.offsets == previous.offsets && a.columns == previous.columns;
        const bool same_values = same_structure && a.values == previous.values;
        if (!same_structure) {
            context.reset();
            context = create_context(n, a.offsets, a.columns.data(), a.values.data());
        }
        if (!same_values) factor(*context, !same_structure, a.values.data());
        // Keep a copy: the benchmark version of the Fortran caller passes an in-place RHS.
        std::vector<double> b(rhs, rhs+n);
        backsolve(*context, n, x, b.data());
        double residual = relative_residual(a, b.data(), x, n);
        std::printf("FRONTISTR_SOLVE call=%d n=%d reuse_analysis=%d reuse_factor=%d residual=%.17g\n",
                    call, n, int(same_structure), int(same_values), residual);
        require(std::isfinite(residual) && residual <= 1e-8, "linear residual exceeds tolerance");
        previous = std::move(a);
        return 0;
    } catch (const std::exception& e) {
        context.reset();
        previous = {};
        std::fprintf(stderr, "FrontISTR PanguLU error: %s\n", e.what());
        return 1;
    }
}
