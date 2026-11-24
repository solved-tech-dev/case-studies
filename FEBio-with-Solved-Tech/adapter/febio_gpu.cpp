#include "gpu_backend.inc"

namespace {
class SolverThreads {
    int previous = omp_get_max_threads();
public:
    SolverThreads() { omp_set_num_threads(40); }
    ~SolverThreads() { omp_set_num_threads(previous); }
};
}

class FebioGpuLU : public LinearSolver {
public:
    explicit FebioGpuLU(FEModel* model) : LinearSolver(model) {}
    SparseMatrix* CreateSparseMatrix(Matrix_Type) override {
        Destroy();
        matrix = new CRSSparseMatrix(0);
        return matrix;
    }
    bool SetSparseMatrix(SparseMatrix* m) override {
        Destroy();
        matrix = dynamic_cast<CRSSparseMatrix*>(m);
        return matrix != nullptr;
    }
    // FENewtonSolver clears and rebuilds the same application-owned matrix
    // after Destroy(). Factor() compares its rebuilt graph before any reuse.
    void Destroy() override { ready=false; }
    bool PreProcess() override { ready=false; return matrix && !failed && LinearSolver::PreProcess(); }
    bool Factor() override {
        try {
            require(matrix && !failed, "missing matrix or prior failure");
            ready=false;
            const int n=matrix->Rows();
            if (n==0) return true;
            const auto nz=matrix->NonZeroes();
            require(n>0 && nz>0 && nz<=INT_MAX, "unsupported system size");
            const int *p=matrix->Pointers(), *i=matrix->Indices();
            const double *a=matrix->Values();
            require(p && i && a && p[0]==0 && p[n]==nz, "invalid original CSR");
            bool fresh=!ctx || ptr.size()!=size_t(n)+1 || col.size()!=size_t(nz)
                || !std::equal(ptr.begin(),ptr.end(),p) || !std::equal(col.begin(),col.end(),i);
            bool unchanged=!fresh && std::memcmp(original.data(),a,sizeof(double)*nz)==0;
            auto start=Clock::now();
            SolverThreads threads;
            if (fresh) {
                ctx.reset();
                ptr.assign(p,p+n+1);
                col.assign(i,i+nz);
                ctx=create_context(n,ptr,i,a);
            }
            if (!unchanged) {
                original.assign(a,a+nz);
                factor(*ctx,fresh,a);
            }
            ready=true;
            std::printf("FEBIO_FACTOR n=%d nz=%lld graph_new=%d identical_values=%d\n",n,(long long)nz,fresh,unchanged);
            event("service","factor",start);
            return true;
        } catch (const std::exception& e) { return error(e); }
    }
    bool BackSolve(double* x,double* b) override {
        try {
            require(matrix && !failed, "missing system");
            const int n=matrix->Rows();
            if (!n) return true;
            require(ready && ctx && b && x, "solve before successful factor");
            ++call;
            // Preserve the RHS even when FEBio calls through an in-place interface.
            saved_rhs.assign(b,b+n);
            SolverThreads threads;
            auto start=Clock::now();
            backsolve(*ctx,n,x,saved_rhs.data());
            event("service","solve",start);
            long double r2=0,b2=0;
            for (int row=0;row<n;++row) {
                require(std::isfinite(x[row]) && std::isfinite(saved_rhs[row]), "nonfinite solution or RHS");
                long double residual=-saved_rhs[row];
                for (auto k=ptr[row];k<ptr[row+1];++k)
                    residual+=(long double)original[k]*x[col[k]];
                r2+=residual*residual;
                b2+=(long double)saved_rhs[row]*saved_rhs[row];
            }
            const auto rel=b2>0 ? std::sqrt(r2/b2) : std::sqrt(r2);
            std::printf("FEBIO_RESIDUAL call=%d relative=%.17Lg\n",call,rel);
            require(std::isfinite(rel) && rel<=1e-8L,"original-CSR residual exceeds 1e-8");
            UpdateStats(1);
            return true;
        } catch (const std::exception& e) { return error(e); }
    }
private:
    bool error(const std::exception& e) {
        std::fprintf(stderr,"FEBIO_GPU_ERROR %s\n",e.what());
        failed=true; ready=false; return false;
    }
    CRSSparseMatrix* matrix=nullptr;
    std::unique_ptr<Context> ctx;
    std::vector<int64_t> ptr;
    std::vector<int> col;
    std::vector<double> original,saved_rhs;
    bool ready=false,failed=false;
    DECLARE_FECORE_CLASS();
};

BEGIN_FECORE_CLASS(FebioGpuLU,LinearSolver)
END_FECORE_CLASS();
FECORE_PLUGIN unsigned int GetSDKVersion() { return FE_SDK_VERSION; }
FECORE_PLUGIN void GetPluginVersion(int& major,int& minor,int& patch) { major=1; minor=0; patch=0; }
FECORE_PLUGIN void PluginInitialize(FECoreKernel& kernel) {
    FECoreKernel::SetInstance(&kernel);
    REGISTER_FECORE_CLASS(FebioGpuLU,"febio_gpu_lu");
}
