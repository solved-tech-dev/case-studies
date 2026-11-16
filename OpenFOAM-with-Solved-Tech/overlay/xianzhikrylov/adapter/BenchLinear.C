// SPDX-License-Identifier: GPL-3.0-or-later
#include "BenchLinear.H"
#include "PrecisionAdaptor.H"
#include "Pstream.H"
#include "processorLduInterface.H"
#include "processorLduInterfaceField.H"
#include "XianzhiKrylovAPI.h"
#include <mpi.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <unistd.h>
#include <sys/mman.h>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cstdlib>

// Private extension in the isolated benchmark core; forwards to the SDK's
// existing CSR per-call tolerance control without changing solver algorithms.
extern "C" xianzhikrylov_status_t xk_bench_set_tolerances(
    xianzhikrylov_handle_t,double,double);

static_assert(sizeof(Foam::scalar)==8 && sizeof(Foam::solveScalar)==8, "FP64 required");
static_assert(sizeof(Foam::label)==4 && sizeof(int)==4, "int32 required");
namespace {
double linearSeconds=0, outerTime=0;
std::ostringstream records;
double nativeParts[3]={0,0,0};
int rankNo=0, rankCount=1;
void barrier() { MPI_Barrier(MPI_COMM_WORLD); }
double finish(double begin) {
    barrier(); double local=MPI_Wtime()-begin, result=0;
    MPI_Allreduce(&local,&result,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    return result;
}
void check(xianzhikrylov_status_t status,const char* op,int expected=XIANZHIKRYLOV_SUCCESS) {
    if(status!=expected) throw std::runtime_error(std::string(op)+" returned "+std::to_string(status));
}
struct Objects {
    xianzhikrylov_handle_t handle=nullptr;
    xianzhikrylov_matrix_t matrix=nullptr;
    xianzhikrylov_vector_t b=nullptr,x=nullptr;
    ~Objects() {
        if(handle) xianzhikrylov_destroy(handle);
        if(x) xianzhikrylov_vector_destroy(x);
        if(b) xianzhikrylov_vector_destroy(b);
        if(matrix) xianzhikrylov_matrix_destroy(matrix);
    }
};
// Buffers are entirely overwritten by our matrix/MPI code. Avoid zero-filling
// another full matrix while all processor-local CSR copies are resident.
template<class T> class ReceiveBuffer {
    std::unique_ptr<T[]> data_;
    std::size_t size_=0,capacity_=0;
public:
    void resize(std::size_t n) {
        if(n>capacity_){data_.reset(new T[n]);capacity_=n;}
        size_=n;
    }
    void shrink(std::size_t n) { if(n<=size_) size_=n; }
    void clear() { data_.reset(); size_=capacity_=0; }
    T* data() { return data_.get(); }
    const T* data() const { return data_.get(); }
    T* begin() { return data(); }
    T* end() { return data()+size_; }
    std::size_t size() const { return size_; }
    T& operator[](std::size_t i) { return data_[i]; }
    const T& operator[](std::size_t i) const { return data_[i]; }
};
struct CSR {
    std::vector<int> rows;
    ReceiveBuffer<int> cols;
    ReceiveBuffer<double> vals;
    ReceiveBuffer<int> src;
};
struct LocalConversionTimings {
    double transfer=0,count=0,alloc=0,fill=0,sort=0;
};
struct ReceiveCSR { ReceiveBuffer<int> rows,cols; ReceiveBuffer<double> vals; };

struct PersistentGPU {
    std::unique_ptr<Objects> objects;
    ReceiveBuffer<double> values,b,x;
    std::string config;
    unsigned creates=0,analyses=0,numerics=0,solves=0;
    void reset() {
        objects.reset();values.clear();b.clear();x.clear();config.clear();
    }
};
PersistentGPU& persistentGPU(){static PersistentGPU p;return p;}
using TargetKey=std::pair<long long,std::string>;
const std::map<TargetKey,double>& residualTargets(){
    static const auto targets=[] {
        std::map<TargetKey,double> result;
        const char* path=std::getenv("XK_CPU_RESIDUAL_TARGETS");
        if(!path)return result;
        std::ifstream input(path);
        if(!input)throw std::runtime_error("Cannot read CPU residual targets");
        long long outer;std::string field;double tolerance;
        while(input>>outer>>field>>tolerance){
            if(!std::isfinite(tolerance)||tolerance<=0||
               !result.emplace(TargetKey{outer,field},tolerance).second)
                throw std::runtime_error("Invalid CPU residual target table");
        }
        if(!input.eof()||result.empty())throw std::runtime_error("Malformed CPU residual target table");
        return result;
    }();
    return targets;
}

// The sparsity pattern of a static mesh is identical for every field and
// component (lduAddr is mesh-owned). Cache one host-side structure for the
// whole process instead of one per field@component: at 65M rows this avoids
// ~6x duplication of rows/cols/src and keeps peak host memory acceptable.
// Every call refreshes matrix values. The GPU context and structure are
// persistent; each new values epoch still receives numeric preparation.
struct PersistentStructure {
    long long signature=0;
    bool analyzed=false;
    int localN=0,localNz=0,lowerCount=0;
    std::vector<int> sizes,offsets;
    ReceiveBuffer<int> rows,cols;                    // rank 0 global structure
    std::vector<int> localRows;                      // per-rank compact row offsets
    ReceiveBuffer<int> localSrc;                     // per-rank value source map
};
static PersistentStructure& persistentStructure(){
    static PersistentStructure s;
    return s;
}
long long structureSignature(const Foam::lduMatrix& m,const Foam::lduInterfaceFieldPtrsList& fields){
    long long s=static_cast<long long>(m.diag().size());
    s=s*1000003+static_cast<long long>(m.lduAddr().lowerAddr().size());
    forAll(fields,p) if(fields.set(p)) s=s*1000003+static_cast<long long>(m.lduAddr().patchAddr(p).size());
    return s;
}


// Includes processor-interface coefficients with the same sign as lduMatrix::Amul.
// Other coupled interfaces and transformed processor interfaces are rejected.
CSR localCSR(const Foam::lduMatrix& m,
    const Foam::FieldField<Foam::Field,Foam::scalar>& bou,
    const Foam::lduInterfaceFieldPtrsList& fields, int offset,
    LocalConversionTimings& timings) {
    using namespace Foam;
    const double stageStart=MPI_Wtime();
    const int n=m.diag().size();
    const auto& lower=m.lduAddr().lowerAddr();
    const auto& upper=m.lduAddr().upperAddr();
    const auto& interfaces=m.mesh().interfaces();
    labelField globalCells(n);
    for(int i=0;i<n;++i) globalCells[i]=offset+i;
    const label firstRequest=UPstream::nRequests();
    forAll(fields,p) if(fields.set(p)) {
        if(!interfaces.set(p) || !isA<processorLduInterface>(interfaces[p]))
            throw std::runtime_error("Only ordinary processor interfaces are supported");
        const auto* proc=dynamic_cast<const processorLduInterfaceField*>(&fields[p]);
        if(!proc || proc->doTransform())
            throw std::runtime_error("Transformed processor interfaces are unsupported");
        interfaces[p].initInternalFieldTransfer(Pstream::commsTypes::nonBlocking,globalCells);
    }
    UPstream::waitRequests(firstRequest);
    List<labelField> remote(fields.size());
    std::int64_t count=std::int64_t(n)+2*lower.size();
    forAll(fields,p) if(fields.set(p)) {
        remote[p]=interfaces[p].internalFieldTransfer(Pstream::commsTypes::nonBlocking,globalCells);
        if(remote[p].size()!=m.lduAddr().patchAddr(p).size())
            throw std::runtime_error("Processor interface addressing mismatch");
        count+=remote[p].size();
    }
    timings.transfer=finish(stageStart);
    if(count>std::numeric_limits<int>::max()) throw std::runtime_error("Local CSR too large");
    double stage=MPI_Wtime();
    CSR a; a.rows.assign(n+1,0);
    for(int i=0;i<n;++i) ++a.rows[i+1];
    forAll(lower,f) {++a.rows[lower[f]+1];++a.rows[upper[f]+1];}
    forAll(fields,p) if(fields.set(p))
        for(label cell:m.lduAddr().patchAddr(p)) ++a.rows[cell+1];
    for(int i=0;i<n;++i) a.rows[i+1]+=a.rows[i];
    timings.count=finish(stage);
    stage=MPI_Wtime();
    a.cols.resize(count); a.vals.resize(count); a.src.resize(count); auto next=a.rows;
    timings.alloc=finish(stage);
    stage=MPI_Wtime();
    // src encodes the value source in 32 bits: 0=diag, +f+1=upper[f],
    // -f-1=lower[f], and values below -lowerCount encode processor patches.
    // Sorted entries keep this mapping so a later call can refresh values
    // without rebuilding the sparsity pattern.
    if(lower.size()>=static_cast<label>(1U<<30))
        throw std::runtime_error("LDU source encoding overflow");
    auto put=[&](int row,int col,double value,int s) {
        if(!std::isfinite(value)) throw std::runtime_error("Nonfinite matrix value");
        const int k=next[row]++; a.cols[k]=col; a.vals[k]=value; a.src[k]=s;
    };
    for(int i=0;i<n;++i) put(i,offset+i,m.diag()[i],0);
    forAll(lower,f) {
        put(lower[f],offset+upper[f],m.upper()[f],static_cast<int>(f)+1);
        put(upper[f],offset+lower[f],m.lower()[f],-static_cast<int>(f)-1);
    }
    forAll(fields,p) if(fields.set(p)) {
        const auto& cells=m.lduAddr().patchAddr(p);
        if(p>=128 || cells.size()>=static_cast<label>(1U<<24))
            throw std::runtime_error("Processor patch source encoding overflow");
        const long long patchBase=-static_cast<long long>(lower.size())-1;
        forAll(cells,f) {
            const long long code=(static_cast<long long>(p)<<24)|f;
            const long long encoded=patchBase-code;
            if(encoded<std::numeric_limits<int>::min())
                throw std::runtime_error("Processor patch source encoding overflow");
            put(cells[f],remote[p][f],-bou[p][f],static_cast<int>(encoded));
        }
    }
    timings.fill=finish(stage);
    stage=MPI_Wtime();
    struct E { int col; int src; double val; };
    std::vector<E> row;
    int write=0;
    for(int i=0;i<n;++i) {
        const int read=a.rows[i], readEnd=a.rows[i+1];
        const int degree=readEnd-read;
        if(static_cast<int>(row.size())<degree) row.resize(degree);
        for(int k=0;k<degree;++k) {
            const int src=read+k;
            row[k]={a.cols[src],a.src[src],a.vals[src]};
        }
        std::sort(row.begin(),row.begin()+degree,
            [](const E& x,const E& y){return x.col<y.col;});
        a.rows[i]=write;
        for(int k=0;k<degree;++k) {
            const E& e=row[k];
            if(write>a.rows[i] && a.cols[write-1]==e.col) a.vals[write-1]+=e.val;
            else {a.cols[write]=e.col;a.src[write]=e.src;a.vals[write++]=e.val;}
        }
    }
    a.rows[n]=write;
    a.cols.shrink(write);
    a.src.shrink(write);
    a.vals.shrink(write);
    timings.sort=finish(stage);
    return a;
}

std::string statistics(xianzhikrylov_handle_t handle) {
    std::cout.flush(); std::fflush(stdout);
    int memoryFd=memfd_create("xk-statistics",MFD_CLOEXEC);
    FILE* file=memoryFd<0?nullptr:fdopen(memoryFd,"w+"); if(!file) throw std::runtime_error("Cannot capture SDK statistics");
    int saved=dup(STDOUT_FILENO);
    if(saved<0 || dup2(fileno(file),STDOUT_FILENO)<0)
        throw std::runtime_error("Cannot redirect statistics");
    auto status=xianzhikrylov_print_info(handle);
    std::cout.flush(); std::fflush(stdout);
    dup2(saved,STDOUT_FILENO);close(saved);std::rewind(file);
    std::string text;char buffer[8192];std::size_t got;
    while((got=std::fread(buffer,1,sizeof(buffer),file))) text.append(buffer,got);
    std::fclose(file);check(status,"print_info");
    return text;
}
}
extern "C" void benchNativePart(int part,double value){if(part>=0&&part<3)nativeParts[part]+=value;}
extern "C" void benchFlush(const char* directory){if(rankNo==0){std::ofstream f(std::string(directory)+"/linear.log");f<<records.str();}}
extern "C" double benchLinearSeconds(){return linearSeconds;}
extern "C" void benchSetOuter(double value){outerTime=value;}
extern "C" void benchInitialize(){(void)residualTargets();}
extern "C" void benchReleaseGPU(){
    const double begin=MPI_Wtime();
    auto& p=persistentGPU();p.reset();
    const double elapsed=finish(begin);
    if(rankNo==0)records<<"BENCH_CONTEXT_RELEASE wall_s="<<elapsed
        <<" handle_creates="<<p.creates<<" analysis_calls="<<p.analyses
        <<" numeric_calls="<<p.numerics<<" solve_calls="<<p.solves<<'\n';
}
namespace Foam {
defineTypeNameAndDebug(BenchLinear,0);
lduMatrix::solver::addsymMatrixConstructorToTable<BenchLinear> benchSym;
lduMatrix::solver::addasymMatrixConstructorToTable<BenchLinear> benchAsym;
BenchLinear::BenchLinear(const word& name,const lduMatrix& m,
    const FieldField<Field,scalar>& bou,const FieldField<Field,scalar>& internal,
    const lduInterfaceFieldPtrsList& interfaces,const dictionary& controls)
    : lduMatrix::solver(name,m,bou,internal,interfaces,controls) {}

solverPerformance BenchLinear::scalarSolve(solveScalarField& psi,
    const solveScalarField& source,const direction cmpt) const {
    try {
        MPI_Comm_rank(MPI_COMM_WORLD,&rankNo);MPI_Comm_size(MPI_COMM_WORLD,&rankCount);
        if(rankCount!=Pstream::nProcs(matrix_.mesh().comm()) || rankNo!=Pstream::myProcNo(matrix_.mesh().comm()))
            throw std::runtime_error("This adapter requires the OpenFOAM world communicator");
        const word backend=controlDict_.get<word>("backend");
        if(backend!="native" && backend!="xianzhi") throw std::runtime_error("Unknown backend");
        double targetRelativeL2=1e-8,effectiveRtol=1e-9;
        if(backend=="xianzhi"){
            const auto& targets=residualTargets();
            if(!targets.empty()){
                auto found=targets.find({std::llround(outerTime),fieldName_.c_str()});
                if(found==targets.end())throw std::runtime_error("Missing CPU target for current equation");
                targetRelativeL2=effectiveRtol=found->second;
            }
        }
        records.precision(17);
        barrier();const double started=MPI_Wtime();
        if(controlDict_.getOrDefault<bool>("zeroInitial",false)) psi=0;
        solverPerformance perf(typeName,fieldName_);
        solveScalarField product(psi.size()),work(psi.size());
        matrix_.Amul(product,psi,interfaceBouCoeffs_,interfaces_,cmpt);
        const solveScalar normalizer=normFactor(psi,source,product,work);
        if(!std::isfinite(normalizer) || normalizer<=0) throw std::runtime_error("Invalid normalizer");
        solveScalarField residual(source-product);
        perf.initialResidual()=gSumMag(residual,matrix_.mesh().comm())/normalizer;
        perf.finalResidual()=perf.initialResidual();
        const double initialCheck=finish(started);
        double nativeCtor=0; std::fill(nativeParts,nativeParts+3,0);
        double nativeWall=0,conversion=0,gather=0,scatter=0,analysis=0,prepare=0,solve=0,sdkLifecycle=0,rootPrepare=0;
        const bool skipped=minIter_==0 && perf.checkConvergence(tolerance_,relTol_,0);
        if(!skipped && backend=="native") {
            const double t=MPI_Wtime();
            {
                dictionary inner(controlDict_);
                word selected=controlDict_.get<word>("nativeSolver");
                if(selected=="smoothSolver")selected="BenchSmoothSolver";
                inner.set("solver",selected);
                const double ctorStart=MPI_Wtime();
                auto native=lduMatrix::solver::New(fieldName_,matrix_,interfaceBouCoeffs_,interfaceIntCoeffs_,interfaces_,inner);
                nativeCtor=MPI_Wtime()-ctorStart;
                perf=native->scalarSolve(psi,source,cmpt);
            }
            nativeWall=finish(t);
            double data[4]={nativeCtor,nativeParts[0],nativeParts[1],nativeParts[2]},maximum[4]={};
            MPI_Allreduce(data,maximum,4,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
            nativeCtor=maximum[0];for(int i=0;i<3;++i)nativeParts[i]=maximum[i+1];
        } else if(!skipped) {
            const int n=psi.size();std::vector<int> sizes(rankCount),offsets(rankCount),nnzs(rankCount),nzOffsets(rankCount);
            MPI_Allgather(&n,1,MPI_INT,sizes.data(),1,MPI_INT,MPI_COMM_WORLD);
            std::int64_t total=0;
            for(int p=0;p<rankCount;++p){offsets[p]=total;total+=sizes[p];}
            if(total>std::numeric_limits<int>::max()) throw std::runtime_error("Global row count exceeds int32");
            PersistentStructure& cache=persistentStructure();
            const long long sig=structureSignature(matrix_,interfaces_);
            const bool structureChanged=!cache.analyzed||cache.signature!=sig;
            double t=MPI_Wtime();
            CSR local;
            LocalConversionTimings conversionParts;
            if(structureChanged) {
                local=localCSR(matrix_,interfaceBouCoeffs_,interfaces_,offsets[rankNo],conversionParts);
                cache.localN=n;cache.localNz=local.vals.size();
                cache.lowerCount=matrix_.lduAddr().lowerAddr().size();
                cache.signature=sig;
                cache.sizes=sizes;cache.offsets=offsets;
                cache.analyzed=true;
            } else {
                local.vals.resize(cache.localNz);
                const auto& m=matrix_;const auto& bou=interfaceBouCoeffs_;
                for(int i=0;i<cache.localN;++i){
                    for(int k=cache.localRows[i];k<cache.localRows[i+1];++k){
                        const int s=cache.localSrc[k];
                        if(s==0) {
                            local.vals[k]=m.diag()[i];
                        } else if(s>0) {
                            local.vals[k]=m.upper()[static_cast<label>(s-1)];
                        } else if(s>=-cache.lowerCount) {
                            local.vals[k]=m.lower()[static_cast<label>(-s-1)];
                        } else {
                            const long long code=
                                -static_cast<long long>(s)-cache.lowerCount-1;
                            const int patch=static_cast<int>(code>>24);
                            const int face=static_cast<int>(code&0xFFFFFFLL);
                            local.vals[k]=-bou[patch][face];
                        }
                    }
                }
            }
            conversion=finish(t);
            if(structureChanged && rankNo==0) records
                << "BENCH_CONVERSION field="<<fieldName_
                <<" total_s="<<conversion
                <<" transfer_s="<<conversionParts.transfer
                <<" count_s="<<conversionParts.count
                <<" alloc_s="<<conversionParts.alloc
                <<" fill_s="<<conversionParts.fill
                <<" sort_s="<<conversionParts.sort<<nl;
            const int nz=cache.localNz;
            MPI_Allgather(&nz,1,MPI_INT,nnzs.data(),1,MPI_INT,MPI_COMM_WORLD);
            std::int64_t totalNz=0;
            for(int p=0;p<rankCount;++p){nzOffsets[p]=totalNz;totalNz+=nnzs[p];}
            if(totalNz>std::numeric_limits<int>::max()) throw std::runtime_error("Global nnz exceeds int32");
            auto& context=persistentGPU();
            if(structureChanged)context.reset();
            auto& gvals=context.values;auto& b=context.b;auto& x=context.x;
            t=MPI_Wtime();
            if(structureChanged){
                if(rankNo==0){cache.rows.resize(total+1);cache.cols.resize(totalNz);}
                if(rankNo==0)cache.rows[0]=0;
                MPI_Gatherv(local.rows.data()+1,n,MPI_INT,rankNo==0?cache.rows.data()+1:nullptr,sizes.data(),offsets.data(),MPI_INT,0,MPI_COMM_WORLD);
                MPI_Gatherv(local.cols.data(),nz,MPI_INT,rankNo==0?cache.cols.data():nullptr,nnzs.data(),nzOffsets.data(),MPI_INT,0,MPI_COMM_WORLD);
                if(rankNo==0)for(int p=0;p<rankCount;++p)
                    for(int i=offsets[p]+1;i<=offsets[p]+sizes[p];++i)cache.rows[i]+=nzOffsets[p];
            }
            if(rankNo==0)gvals.resize(totalNz);
            MPI_Gatherv(local.vals.data(),nz,MPI_DOUBLE,rankNo==0?gvals.data():nullptr,nnzs.data(),nzOffsets.data(),MPI_DOUBLE,0,MPI_COMM_WORLD);
            // Move the local structure into the process cache after the first gather.
            if(structureChanged) {
                cache.localRows=std::move(local.rows);
                cache.localSrc=std::move(local.src);
            } else {
                std::vector<int>().swap(local.rows);
                local.src.clear();
            }
            if(rankNo==0)b.resize(total);
            MPI_Gatherv(source.cdata(),n,MPI_DOUBLE,b.data(),sizes.data(),offsets.data(),MPI_DOUBLE,0,MPI_COMM_WORLD);
            if(rankNo==0)x.resize(total);
            MPI_Gatherv(psi.cdata(),n,MPI_DOUBLE,x.data(),sizes.data(),offsets.data(),MPI_DOUBLE,0,MPI_COMM_WORLD);
            gather=finish(t);
            int iterations=0;
            t=MPI_Wtime();
            if(rankNo==0) {
                if(controlDict_.getOrDefault<bool>("externalAMG",false)) {
                    bool allNeg=true,allPos=true;
                    for(int i=0;i<total;++i){
                        auto first=cache.cols.begin()+cache.rows[i],last=cache.cols.begin()+cache.rows[i+1];
                        auto d=std::lower_bound(first,last,i);
                        if(d==last||*d!=i) throw std::runtime_error("Missing global diagonal");
                        double v=gvals[d-cache.cols.begin()];allNeg&=v<0;allPos&=v>0;
                    }
                    if(!allNeg&&!allPos) throw std::runtime_error("AMG requires uniform diagonal sign");
                    if(allNeg){for(auto& v:gvals)v=-v;for(auto& v:b)v=-v;}
                }
                rootPrepare=MPI_Wtime()-t;
                const double sdkStart=MPI_Wtime();
                {
                    fileName config(controlDict_.get<fileName>("configFile"));config.expand();
                    if(context.objects && context.config!=config.c_str())
                        throw std::runtime_error("A persistent GPU context requires one solver configuration");
                    const bool first=!context.objects;
                    if(first){
                    context.objects.reset(new Objects);
                    context.config=config.c_str();++context.creates;
                    }
                    auto& d=*context.objects;
                    if(first){
                    check(xianzhikrylov_create(0,&d.handle),"create");
                    check(xianzhikrylov_config_set(d.handle,XIANZHIKRYLOV_CONFIG_FILE,config.c_str()),"config_set");
                    xianzhikrylov_matrix_data_t data{};data.abi_version=XIANZHIKRYLOV_ABI_VERSION;data.struct_size=sizeof(data);
                    data.format=XIANZHIKRYLOV_CSR;data.value_type=XIANZHIKRYLOV_FP64;data.index_type=XIANZHIKRYLOV_INDEX_I32;
                    data.block_size=1;data.block_layout=XIANZHIKRYLOV_BLOCK_ROW_MAJOR;data.memory_space=XIANZHIKRYLOV_HOST;data.device_id=-1;
                    data.rows=data.columns=total;data.entries=totalNz;
                    data.row_offsets={cache.rows.data(),cache.rows.size()*sizeof(int)};
                    data.column_indices={cache.cols.data(),cache.cols.size()*sizeof(int)};
                    data.values={gvals.data(),gvals.size()*sizeof(double)};
                    check(xianzhikrylov_matrix_create(&d.matrix,&data),"matrix_create");
                    check(xianzhikrylov_vector_create(&d.b,total,b.data(),XIANZHIKRYLOV_FP64,XIANZHIKRYLOV_HOST),"b_create");
                    check(xianzhikrylov_vector_create(&d.x,total,x.data(),XIANZHIKRYLOV_FP64,XIANZHIKRYLOV_HOST),"x_create");
                    double phase=MPI_Wtime();check(xianzhikrylov_execute(d.handle,XIANZHIKRYLOV_PHASE_ANALYSIS,d.matrix,nullptr,nullptr),"analysis");analysis=MPI_Wtime()-phase;
                    ++context.analyses;
                    }
                    check(xk_bench_set_tolerances(d.handle,effectiveRtol,1e-30),"per_call_tolerances");
                    double phase=MPI_Wtime();check(xianzhikrylov_execute(d.handle,XIANZHIKRYLOV_PHASE_NUMERIC_PREPARE,d.matrix,nullptr,nullptr),"prepare");prepare=MPI_Wtime()-phase;
                    ++context.numerics;
                    phase=MPI_Wtime();auto status=xianzhikrylov_execute(d.handle,XIANZHIKRYLOV_PHASE_SOLVE,d.matrix,d.x,d.b);solve=MPI_Wtime()-phase;
                    ++context.solves;
                    const std::string stats=statistics(d.handle);records<< "BENCH_SDK outer="<<outerTime<<" field="<<fieldName_<<"\n"<<stats;
                    check(status,"solve",XIANZHIKRYLOV_CONVERGED);
                    auto pos=stats.find(" iterations=");if(pos==std::string::npos) throw std::runtime_error("No SDK iteration statistic");
                    iterations=std::stoi(stats.substr(pos+12));
                }
                sdkLifecycle=MPI_Wtime()-sdkStart;
            }
            barrier();
            MPI_Bcast(&iterations,1,MPI_INT,0,MPI_COMM_WORLD);perf.nIterations()=iterations;
            t=MPI_Wtime();
            MPI_Scatterv(x.data(),sizes.data(),offsets.data(),MPI_DOUBLE,psi.data(),n,MPI_DOUBLE,0,MPI_COMM_WORLD);
            scatter=finish(t);
            if(rankNo==0) records<< "BENCH_GPU outer="<<outerTime<<" field="<<fieldName_
                <<" rows="<<total<<" nnz="<<totalNz<<" conversion_s="<<conversion
                <<" conversion_transfer_s="<<conversionParts.transfer
                <<" conversion_count_s="<<conversionParts.count
                <<" conversion_alloc_s="<<conversionParts.alloc
                <<" conversion_fill_s="<<conversionParts.fill
                <<" conversion_sort_s="<<conversionParts.sort
                <<" gather_s="<<gather<<" scatter_s="<<scatter<<" root_sign_prepare_s="<<rootPrepare<<" sdk_lifecycle_wall_s="<<sdkLifecycle
                <<" analysis_api_wall_s="<<analysis<<" prepare_api_wall_s="<<prepare<<" solve_api_wall_s="<<solve
                <<" handle_creates="<<context.creates<<" analysis_calls="<<context.analyses
                <<" numeric_calls="<<context.numerics<<" solve_calls="<<context.solves<<nl;
        }
        const double verifyStart=MPI_Wtime();
        matrix_.Amul(product,psi,interfaceBouCoeffs_,interfaces_,cmpt);
        residual=source-product;
        perf.finalResidual()=gSumMag(residual,matrix_.mesh().comm())/normalizer;
        long double localNorm[2]={0,0},globalNorm[2]={0,0};
        forAll(residual,i){
            if(!std::isfinite(residual[i]) || !std::isfinite(psi[i])) throw std::runtime_error("Nonfinite solution/residual");
            localNorm[0]+=static_cast<long double>(residual[i])*residual[i];localNorm[1]+=static_cast<long double>(source[i])*source[i];
        }
        MPI_Allreduce(localNorm,globalNorm,2,MPI_LONG_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        const double relativeL2=globalNorm[1]>0?std::sqrt(globalNorm[0]/globalNorm[1]):std::sqrt(globalNorm[0]);
        const bool foamConverged=perf.checkConvergence(tolerance_,relTol_,0);
        const bool accepted=backend=="native" || relativeL2<=targetRelativeL2*(1+1e-6);
        matrix_.setResidualField(residual,fieldName_,false);
        product.clear();work.clear();residual.clear();
        const double verification=finish(verifyStart),wall=finish(started);linearSeconds+=wall;
        if(rankNo==0) records<< "BENCH_LINEAR outer="<<outerTime<<" field="<<fieldName_<<" backend="<<backend
            <<" wall_s="<<wall<<" initial_check_s="<<initialCheck<<" verify_s="<<verification
            <<" native_call_s="<<nativeWall<<" native_ctor_s="<<nativeCtor
            <<" native_initial_s="<<nativeParts[0]<<" native_smoother_setup_s="<<nativeParts[1]<<" native_iterate_s="<<nativeParts[2]<<" sdk_lifecycle_wall_s="<<sdkLifecycle
            <<" wrapper_other_wall_s="<<wall-(backend=="native"?nativeWall:sdkLifecycle)<<" iterations="<<perf.nIterations()
            <<" initial_residual="<<perf.initialResidual()<<" final_residual="<<perf.finalResidual()
            <<" original_relative_l2="<<relativeL2<<" foam_converged="<<foamConverged<<" accepted="<<accepted<<" skipped="<<skipped<<nl;
        if(rankNo==0 && backend=="xianzhi")records<<"BENCH_TARGET outer="<<outerTime<<" field="<<fieldName_
            <<" target_relative_l2="<<targetRelativeL2<<" effective_rtol="<<effectiveRtol
            <<" actual_relative_l2="<<relativeL2<<" actual_to_cpu_target="<<relativeL2/targetRelativeL2<<nl;
        if(!accepted) throw std::runtime_error("Original OpenFOAM residual target was not met for "+std::string(fieldName_));
        return perf;
    } catch(const std::exception& e) {
        std::cerr<<"BENCH_FATAL rank="<<rankNo<<" "<<e.what()<<std::endl;
        MPI_Abort(MPI_COMM_WORLD,81);
    }
    return solverPerformance(typeName,fieldName_);
}
solverPerformance BenchLinear::solve(scalarField& psi,const scalarField& source,const direction cmpt) const {
    PrecisionAdaptor<solveScalar,scalar> p(psi);ConstPrecisionAdaptor<solveScalar,scalar> b(source);
    return scalarSolve(p.ref(),b(),cmpt);
}
}
