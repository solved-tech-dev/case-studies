// SPDX-License-Identifier: GPL-3.0-or-later
// MPI OpenFOAM pressure bridge for the unmodified LinglongAMG 26.46 C API.
#include "FoamPressure.H"
#include "PrecisionAdaptor.H"
#include "Pstream.H"
#include "processorLduInterface.H"
#include "processorLduInterfaceField.H"
#include "linglongamg.h"
#include "SHA1.H"
#include <cuda_runtime_api.h>
#include <mpi.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <numeric>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>

static_assert(sizeof(Foam::scalar)==8 && sizeof(Foam::solveScalar)==8, "FP64 required");
static_assert(sizeof(Foam::label)==4, "OpenFOAM DPInt32 required");
namespace {
struct Record {
    double outer=0,wall=0,initial=0,topology=0,conversion=0,gather=0;
    double create=0,alloc=0,h2d=0,reset=0,setup=0,solve=0,query=0,d2h=0,scatter=0,verify=0;
    double nativeCreate=0,nativeInitial=0,nativeSetup=0,nativeIterate=0,nativeFinalize=0,nativeDestroy=0;
    double sdkSetupMs=0,sdkSolveMs=0,initialFoam=0,finalFoam=0,initialL2=0,finalL2=0,rhsL2=0;
    double assemblyAuditError=0;
    double fingerprint=0;
    int repeat=0,hierarchyReuseRequested=0;
    std::string topologyDigest,matrixDigest,rhsDigest,x0Digest;
    int rows=0,nnz=0,iterations=0,sdkConverged=0,foamConverged=0,structureReused=0;
    std::string backend;
};
std::vector<Record> records;
Record* active=nullptr;
double currentOuter=0;
int currentRepeat=0;
const char* lastDeviceStage="not entered";
struct NativeResidual {double outer,l2,relative;int repeat,iteration;};
std::vector<NativeResidual> nativeResiduals;
double now(){return MPI_Wtime();}
double maxWall(double begin){double v=now()-begin,out=0;MPI_Allreduce(&v,&out,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);return out;}
void barrier(){MPI_Barrier(MPI_COMM_WORLD);}
void check(cudaError_t s){if(s!=cudaSuccess)throw std::runtime_error(cudaGetErrorString(s));}
void check(linglongamg_status_t s){if(s!=LINGLONGAMG_STATUS_SUCCESS)throw std::runtime_error(linglongamg_status_string(s));}
template<class Field> void hashField(Foam::SHA1& hash,const char* tag,const Field& field){
    hash.append(tag);
    const std::uint64_t count=field.size();
    hash.append(reinterpret_cast<const char*>(&count),sizeof(count));
    hash.append(reinterpret_cast<const char*>(field.cdata()),std::size_t(count)*sizeof(field[0]));
}
std::string distributedDigest(const Foam::SHA1& local){
    int ranks=1;MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    const std::string hex=local.str();std::vector<char> all(std::size_t(ranks)*40);
    MPI_Allgather(hex.data(),40,MPI_CHAR,all.data(),40,MPI_CHAR,MPI_COMM_WORLD);
    Foam::SHA1 global;global.append(all.data(),all.size());return global.str();
}
void fingerprintSystem(const Foam::lduMatrix& matrix,
    const Foam::FieldField<Foam::Field,Foam::scalar>& bou,
    const Foam::FieldField<Foam::Field,Foam::scalar>& internal,
    const Foam::lduInterfaceFieldPtrsList& fields,
    const Foam::solveScalarField& rhs,const Foam::solveScalarField& x0,Record& r){
    using namespace Foam;
    SHA1 topologyHash,matrixHash,rhsHash,xHash;
    const std::uint64_t n=matrix.diag().size();
    topologyHash.append(reinterpret_cast<const char*>(&n),sizeof(n));
    hashField(topologyHash,"lowerAddr",matrix.lduAddr().lowerAddr());
    hashField(topologyHash,"upperAddr",matrix.lduAddr().upperAddr());
    hashField(matrixHash,"diag",matrix.diag());
    hashField(matrixHash,"upper",matrix.upper());
    hashField(matrixHash,"lower",matrix.lower());
    const auto& ifaces=matrix.mesh().interfaces();
    forAll(fields,p){
        const int present=fields.set(p)?1:0,index=p;
        topologyHash.append(reinterpret_cast<const char*>(&index),sizeof(index));
        topologyHash.append(reinterpret_cast<const char*>(&present),sizeof(present));
        if(present){
            hashField(topologyHash,"patchAddr",matrix.lduAddr().patchAddr(p));
            if(!ifaces.set(p)||!isA<processorLduInterface>(ifaces[p]))
                throw std::runtime_error("Fingerprint requires processor interfaces");
            const int neighbour=refCast<const processorLduInterface>(ifaces[p]).neighbProcNo();
            topologyHash.append(reinterpret_cast<const char*>(&neighbour),sizeof(neighbour));
            hashField(matrixHash,"boundaryCoeffs",bou[p]);
            hashField(matrixHash,"internalCoeffs",internal[p]);
        }
    }
    hashField(rhsHash,"rhs",rhs);hashField(xHash,"x0",x0);
    r.topologyDigest=distributedDigest(topologyHash);r.matrixDigest=distributedDigest(matrixHash);
    r.rhsDigest=distributedDigest(rhsHash);r.x0Digest=distributedDigest(xHash);
}
template<class T> class Buffer {
    std::unique_ptr<T[]> ptr_;
public:
    void resize(std::size_t n){ptr_.reset(n?new T[n]:nullptr);}
    T* data(){return ptr_.get();}
    const T* data()const{return ptr_.get();}
    T& operator[](std::size_t i){return ptr_[i];}
};
struct Topology {
    bool ready=false;
    int n=0,nz=0,total=0,totalNz=0,rank=0,ranks=1;
    std::vector<int> sizes,offsets,nnzs,nzOffsets,rows,cols,diagPos,upPos,loPos,rowOrder;
    std::vector<std::vector<int>> patchPos;
    std::vector<double> values,localRhs,localX;
    Buffer<int> globalRows,globalCols;
    Buffer<double> globalValues,b,x;
    void build(const Foam::lduMatrix& m,
        const Foam::lduInterfaceFieldPtrsList& fields,const Foam::dictionary& controls) {
        using namespace Foam;
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);
        n=m.diag().size(); sizes.resize(ranks); offsets.resize(ranks);
        MPI_Allgather(&n,1,MPI_INT,sizes.data(),1,MPI_INT,MPI_COMM_WORLD);
        std::int64_t count=0;
        for(int p=0;p<ranks;++p){offsets[p]=count;count+=sizes[p];}
        if(count>std::numeric_limits<int>::max())throw std::runtime_error("Global cell count exceeds SDK int32");
        total=count;
        rowOrder.resize(n);std::iota(rowOrder.begin(),rowOrder.end(),0);
        if(controls.found("permutationFile")){
            std::string path=controls.get<fileName>("permutationFile");
            const auto marker=path.find("%r");
            if(marker!=std::string::npos)path.replace(marker,2,std::to_string(rank));
            std::ifstream input(path,std::ios::binary);
            if(!input.read(reinterpret_cast<char*>(rowOrder.data()),std::size_t(n)*sizeof(int)))
                throw std::runtime_error("Cannot read GPU row permutation: "+path);
            std::vector<unsigned char> seen(n,0);
            for(int value:rowOrder){
                if(value<0||value>=n||seen[value])throw std::runtime_error("GPU row permutation is not a bijection");
                seen[value]=1;
            }
        }
        localRhs.resize(n);localX.resize(n);
        labelField ids(n);for(int i=0;i<n;++i)ids[i]=offsets[rank]+rowOrder[i];
        const auto& ifaces=m.mesh().interfaces();
        const label first=UPstream::nRequests();
        forAll(fields,p)if(fields.set(p)) {
            if(!ifaces.set(p)||!isA<processorLduInterface>(ifaces[p]))
                throw std::runtime_error("Only untransformed processor coupled patches are supported");
            const auto* f=dynamic_cast<const processorLduInterfaceField*>(&fields[p]);
            if(!f||f->doTransform())throw std::runtime_error("Transformed processor patch is unsupported");
            ifaces[p].initInternalFieldTransfer(Pstream::commsTypes::nonBlocking,ids);
        }
        UPstream::waitRequests(first);
        List<labelField> remote(fields.size());
        const auto& lo=m.lduAddr().lowerAddr(); const auto& up=m.lduAddr().upperAddr();
        std::int64_t nz64=std::int64_t(n)+2*lo.size();
        forAll(fields,p)if(fields.set(p)) {
            remote[p]=ifaces[p].internalFieldTransfer(Pstream::commsTypes::nonBlocking,ids);
            if(remote[p].size()!=m.lduAddr().patchAddr(p).size())throw std::runtime_error("Processor addressing mismatch");
            nz64+=remote[p].size();
        }
        if(nz64>std::numeric_limits<int>::max())throw std::runtime_error("Local nnz exceeds SDK int32");
        rows.assign(n+1,0);
        for(int i=0;i<n;++i)++rows[i+1];
        forAll(lo,f){++rows[rowOrder[lo[f]]+1];++rows[rowOrder[up[f]]+1];}
        forAll(fields,p)if(fields.set(p))for(label c:m.lduAddr().patchAddr(p))++rows[rowOrder[c]+1];
        for(int i=0;i<n;++i)rows[i+1]+=rows[i];
        cols.resize(nz64);auto next=rows;
        auto put=[&](int r,int c){cols[next[r]++]=c;};
        for(int i=0;i<n;++i)put(rowOrder[i],offsets[rank]+rowOrder[i]);
        forAll(lo,f){put(rowOrder[lo[f]],offsets[rank]+rowOrder[up[f]]);put(rowOrder[up[f]],offsets[rank]+rowOrder[lo[f]]);}
        forAll(fields,p)if(fields.set(p)){
            const auto& cells=m.lduAddr().patchAddr(p);
            forAll(cells,f)put(rowOrder[cells[f]],remote[p][f]);
        }
        // Compact duplicate columns once, then retain exact coefficient-to-CSR maps.
        const auto oldRows=rows;int w=0;
        for(int i=0;i<n;++i){
            std::sort(cols.begin()+oldRows[i],cols.begin()+oldRows[i+1]);
            rows[i]=w;int prev=-1;
            for(int k=oldRows[i];k<oldRows[i+1];++k)if(cols[k]!=prev){prev=cols[k];cols[w++]=prev;}
        }
        rows[n]=w;cols.resize(w);nz=w;values.resize(nz);
        auto pos=[&](int r,int c){
            auto first=cols.begin()+rows[r],last=cols.begin()+rows[r+1];
            auto it=std::lower_bound(first,last,c);
            if(it==last||*it!=c)throw std::runtime_error("CSR map failure");
            return int(it-cols.begin());
        };
        diagPos.resize(n);for(int i=0;i<n;++i)diagPos[i]=pos(rowOrder[i],offsets[rank]+rowOrder[i]);
        upPos.resize(lo.size());loPos.resize(lo.size());
        forAll(lo,f){upPos[f]=pos(rowOrder[lo[f]],offsets[rank]+rowOrder[up[f]]);loPos[f]=pos(rowOrder[up[f]],offsets[rank]+rowOrder[lo[f]]);}
        patchPos.resize(fields.size());
        forAll(fields,p)if(fields.set(p)){
            const auto& cells=m.lduAddr().patchAddr(p);patchPos[p].resize(cells.size());
            forAll(cells,f)patchPos[p][f]=pos(rowOrder[cells[f]],remote[p][f]);
        }
        nnzs.resize(ranks);nzOffsets.resize(ranks);
        MPI_Allgather(&nz,1,MPI_INT,nnzs.data(),1,MPI_INT,MPI_COMM_WORLD);
        count=0;for(int p=0;p<ranks;++p){nzOffsets[p]=count;count+=nnzs[p];}
        if(count>std::numeric_limits<int>::max())throw std::runtime_error("Global nnz exceeds SDK int32");
        totalNz=count;
        if(rank==0){globalRows.resize(std::size_t(total)+1);globalRows[0]=0;globalCols.resize(totalNz);}
        MPI_Gatherv(rows.data()+1,n,MPI_INT,rank==0?globalRows.data()+1:nullptr,sizes.data(),offsets.data(),MPI_INT,0,MPI_COMM_WORLD);
        MPI_Gatherv(cols.data(),nz,MPI_INT,globalCols.data(),nnzs.data(),nzOffsets.data(),MPI_INT,0,MPI_COMM_WORLD);
        if(rank==0){
            for(int p=0;p<ranks;++p)for(int i=offsets[p]+1;i<=offsets[p]+sizes[p];++i)globalRows[i]+=nzOffsets[p];
            globalValues.resize(totalNz);b.resize(total);x.resize(total);
        }
        ready=true;
    }
    void update(const Foam::lduMatrix& m,const Foam::FieldField<Foam::Field,Foam::scalar>& bou,
        const Foam::lduInterfaceFieldPtrsList& fields,int sign){
        if(m.diag().size()!=n||m.upper().size()!=Foam::label(upPos.size()))throw std::runtime_error("Topology changed; rebuild required");
        std::fill(values.begin(),values.end(),0.0);
        for(int i=0;i<n;++i)values[diagPos[i]]+=sign*m.diag()[i];
        forAll(m.upper(),f){values[upPos[f]]+=sign*m.upper()[f];values[loPos[f]]+=sign*m.lower()[f];}
        forAll(fields,p)if(fields.set(p))forAll(bou[p],f)values[patchPos[p][f]]-=sign*bou[p][f];
    }
} topology;

struct Device {
    linglongamg_handle_t handle=nullptr;
    linglongamg_csr_t* csr=nullptr;
    int *rows=nullptr,*cols=nullptr;
    double *values=nullptr,*b=nullptr,*x=nullptr;
    bool used=false;
    void close(){
        if(handle){check(linglongamg_destroy(handle));handle=nullptr;}
        if(csr){linglongamg_csr_destroy(csr);csr=nullptr;}
        if(rows){check(cudaFree(rows));rows=nullptr;}if(cols){check(cudaFree(cols));cols=nullptr;}
        if(values){check(cudaFree(values));values=nullptr;}if(b){check(cudaFree(b));b=nullptr;}if(x){check(cudaFree(x));x=nullptr;}
    }
    void run(const Foam::dictionary& dict,Record& r){
        auto trace=[](const char* stage){lastDeviceStage=stage;if(std::getenv("FOAM_PRESSURE_TRACE")){std::fprintf(stderr,"FP_STAGE %s\n",stage);std::fflush(stderr);}};
        trace("cuda_device");
        check(cudaSetDevice(dict.getOrDefault<Foam::label>("device",0)));
        if(!handle){
            double t=now();Foam::fileName config(dict.get<Foam::fileName>("configFile"));config.expand();
            check(linglongamg_create_from_config(&handle,config.c_str()));r.create=now()-t;
            t=now();
            check(cudaMalloc(reinterpret_cast<void**>(&rows),(std::size_t(topology.total)+1)*sizeof(int)));
            check(cudaMalloc(reinterpret_cast<void**>(&cols),std::size_t(topology.totalNz)*sizeof(int)));
            check(cudaMalloc(reinterpret_cast<void**>(&values),std::size_t(topology.totalNz)*sizeof(double)));
            check(cudaMalloc(reinterpret_cast<void**>(&b),std::size_t(topology.total)*sizeof(double)));
            check(cudaMalloc(reinterpret_cast<void**>(&x),std::size_t(topology.total)*sizeof(double)));
            r.alloc=now()-t;
            t=now();
            check(cudaMemcpy(rows,topology.globalRows.data(),(std::size_t(topology.total)+1)*sizeof(int),cudaMemcpyHostToDevice));
            check(cudaMemcpy(cols,topology.globalCols.data(),std::size_t(topology.totalNz)*sizeof(int),cudaMemcpyHostToDevice));
            r.h2d=now()-t;
            check(linglongamg_read_csr(&csr,topology.total,topology.total,topology.totalNz,rows,cols,values));
            topology.globalRows.resize(0);topology.globalCols.resize(0);
        }
        const bool reuseHierarchy=std::getenv("FOAM_PRESSURE_HIERARCHY_REUSE") && std::string(std::getenv("FOAM_PRESSURE_HIERARCHY_REUSE"))=="1";
        r.hierarchyReuseRequested=used&&reuseHierarchy;
        if(used&&!reuseHierarchy){double t=now();check(linglongamg_reset(handle));check(cudaDeviceSynchronize());r.reset=now()-t;}
        double t=now();
        check(cudaMemcpy(values,topology.globalValues.data(),std::size_t(topology.totalNz)*sizeof(double),cudaMemcpyHostToDevice));
        check(cudaMemcpy(b,topology.b.data(),std::size_t(topology.total)*sizeof(double),cudaMemcpyHostToDevice));
        check(cudaMemcpy(x,topology.x.data(),std::size_t(topology.total)*sizeof(double),cudaMemcpyHostToDevice));
        r.h2d+=now()-t;
        trace("setup_begin");
        t=now();check(linglongamg_setup(handle,csr));check(cudaDeviceSynchronize());r.setup=now()-t;
        trace("setup_end");
        linglongamg_vec_t rhs{topology.total,b},sol{topology.total,x};
        trace("solve_begin");
        t=now();check(linglongamg_solve(handle,&rhs,&sol));check(cudaDeviceSynchronize());r.solve=now()-t;
        trace("solve_end");
        // SDK 26.46 only publishes x after this independently timed authorization query.
        linglongamg_solve_info_t info{};
        t=now();check(linglongamg_get_solve_info(handle,&info));check(cudaDeviceSynchronize());r.query=now()-t;
        r.sdkSetupMs=info.setup_ms;r.sdkSolveMs=info.solve_ms;r.iterations=info.iterations;r.sdkConverged=info.converged;
        if(!std::isfinite(info.final_residual)||!std::isfinite(info.setup_ms)||!std::isfinite(info.solve_ms))
            throw std::runtime_error("SDK returned invalid solve information");
        t=now();check(cudaMemcpy(topology.x.data(),x,std::size_t(topology.total)*sizeof(double),cudaMemcpyDeviceToHost));r.d2h=now()-t;
        used=true;
    }
} device;
}
extern "C" void fpOuter(double v){currentOuter=v;}
extern "C" void fpNativePart(int part,double seconds){
    if(!active)return;
    double maximum=0;MPI_Allreduce(&seconds,&maximum,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);seconds=maximum;
    if(part==0)active->nativeInitial=seconds;
    if(part==1)active->nativeSetup=seconds;
    if(part==2)active->nativeIterate=seconds;
    if(part==3)active->nativeFinalize=seconds;
}
extern "C" void fpNativeInitialResidual(long double q){
    long double sum=0;MPI_Allreduce(&q,&sum,1,MPI_LONG_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    if(active)active->initialL2=std::sqrt(sum);
}
extern "C" double fpNativeTrueRelativeL2(long double q,int iteration){
    long double sum=0;MPI_Allreduce(&q,&sum,1,MPI_LONG_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    const double l2=std::sqrt(sum);
    const double relative=active&&active->initialL2>0 ? l2/active->initialL2 : (l2==0?0:std::numeric_limits<double>::infinity());
    if(active)nativeResiduals.push_back({currentOuter,l2,relative,currentRepeat,iteration});
    return relative;
}
extern "C" void fpShutdown(){int rank=0;MPI_Comm_rank(MPI_COMM_WORLD,&rank);if(rank==0)device.close();}
extern "C" void fpFlush(const char* directory){
    int rank=0;MPI_Comm_rank(MPI_COMM_WORLD,&rank);if(rank)return;
    std::ofstream f(std::string(directory)+"/pressure.csv");f<<std::setprecision(17);
    f<<"outer,backend,rows,nnz,iterations,sdk_converged,foam_converged,structure_reused,wall_s,initial_check_s,topology_s,conversion_s,gather_s,create_s,allocate_s,h2d_s,reset_s,setup_api_s,solve_api_s,query_s,d2h_s,scatter_s,verify_s,native_constructor_s,native_initial_s,native_vcycle_setup_s,native_iterate_s,native_finalize_s,native_destroy_s,sdk_setup_ms,sdk_solve_ms,initial_foam,final_foam,initial_l2,final_l2,rhs_l2,assembly_audit_relative_error\n";
    for(const auto& r:records)f<<r.outer<<','<<r.backend<<','<<r.rows<<','<<r.nnz<<','<<r.iterations<<','<<r.sdkConverged<<','<<r.foamConverged<<','<<r.structureReused<<','<<r.wall<<','<<r.initial<<','<<r.topology<<','<<r.conversion<<','<<r.gather<<','<<r.create<<','<<r.alloc<<','<<r.h2d<<','<<r.reset<<','<<r.setup<<','<<r.solve<<','<<r.query<<','<<r.d2h<<','<<r.scatter<<','<<r.verify<<','<<r.nativeCreate<<','<<r.nativeInitial<<','<<r.nativeSetup<<','<<r.nativeIterate<<','<<r.nativeFinalize<<','<<r.nativeDestroy<<','<<r.sdkSetupMs<<','<<r.sdkSolveMs<<','<<r.initialFoam<<','<<r.finalFoam<<','<<r.initialL2<<','<<r.finalL2<<','<<r.rhsL2<<','<<r.assemblyAuditError<<'\n';
    std::ofstream hashes(std::string(directory)+"/system-fingerprints.csv");hashes<<std::setprecision(17);
    hashes<<"outer,linear_repeat,backend,topology_sha1,matrix_sha1,rhs_sha1,x0_sha1,fingerprint_s,hierarchy_reuse_requested\n";
    for(const auto& r:records)hashes<<r.outer<<','<<r.repeat<<','<<r.backend<<','<<r.topologyDigest<<','<<r.matrixDigest<<','<<r.rhsDigest<<','<<r.x0Digest<<','<<r.fingerprint<<','<<r.hierarchyReuseRequested<<'\n';
    std::ofstream curve(std::string(directory)+"/native-true-residual.csv");curve<<std::setprecision(17);
    curve<<"outer,linear_repeat,iteration,true_l2,true_relative_initial_l2\n";
    for(const auto& x:nativeResiduals)curve<<x.outer<<','<<x.repeat<<','<<x.iteration<<','<<x.l2<<','<<x.relative<<'\n';
}
namespace Foam {
defineTypeNameAndDebug(FoamPressure,0);
lduMatrix::solver::addsymMatrixConstructorToTable<FoamPressure> foamPressureSym;
lduMatrix::solver::addasymMatrixConstructorToTable<FoamPressure> foamPressureAsym;
FoamPressure::FoamPressure(const word& name,const lduMatrix& m,
    const FieldField<Field,scalar>& bou,const FieldField<Field,scalar>& internal,
    const lduInterfaceFieldPtrsList& fields,const dictionary& controls)
    : lduMatrix::solver(name,m,bou,internal,fields,controls){}
solverPerformance FoamPressure::scalarSolveOnce(solveScalarField& psi,const solveScalarField& source,const direction cmpt)const{
    try{
        int rank=0;MPI_Comm_rank(MPI_COMM_WORLD,&rank);
        if(Pstream::nProcs(matrix_.mesh().comm())!=Pstream::nProcs())throw std::runtime_error("World communicator required");
        barrier();const double begin=now();Record r;r.outer=currentOuter;r.repeat=currentRepeat;double initialNorm=0;
        r.backend=controlDict_.get<word>("backend");active=&r;
        const double digestBegin=now();fingerprintSystem(matrix_,interfaceBouCoeffs_,interfaceIntCoeffs_,interfaces_,source,psi,r);r.fingerprint=maxWall(digestBegin);
        solverPerformance perf(typeName,fieldName_);
        if(r.backend=="native"){
            long long localSize[2]={static_cast<long long>(psi.size()),static_cast<long long>(psi.size())+2*static_cast<long long>(matrix_.upper().size())},globalSize[2]{};
            forAll(interfaces_,p)if(interfaces_.set(p))localSize[1]+=interfaceBouCoeffs_[p].size();
            MPI_Allreduce(localSize,globalSize,2,MPI_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);
            r.rows=globalSize[0];r.nnz=globalSize[1];
            dictionary inner(controlDict_);inner.set("solver",word("GAMG"));inner.set("log",0);
            double t=now();
            auto native=lduMatrix::solver::New(fieldName_,matrix_,interfaceBouCoeffs_,interfaceIntCoeffs_,interfaces_,inner);
            r.nativeCreate=maxWall(t);
            perf=native->solve(psi,source,cmpt);
            t=now();native.clear();r.nativeDestroy=maxWall(t);
            r.iterations=perf.nIterations();r.initialFoam=perf.initialResidual();
        }else if(r.backend=="linglongamg"){
            double t=now();
            solveScalarField prod(psi.size()),work(psi.size());
            matrix_.Amul(prod,psi,interfaceBouCoeffs_,interfaces_,cmpt);
            const double norm=normFactor(psi,source,prod,work);initialNorm=norm;
            solveScalarField residual(source-prod);
            r.initialFoam=gSumMag(residual,matrix_.mesh().comm())/norm;
            long double q=0;forAll(residual,i)q+=static_cast<long double>(residual[i])*residual[i];
            long double qGlobal=0;MPI_Allreduce(&q,&qGlobal,1,MPI_LONG_DOUBLE,MPI_SUM,MPI_COMM_WORLD);r.initialL2=std::sqrt(qGlobal);
            prod.clear();work.clear();residual.clear();
            r.initial=maxWall(t);
            r.structureReused=topology.ready;
            if(!topology.ready){t=now();topology.build(matrix_,interfaces_,controlDict_);r.topology=maxWall(t);}
            r.rows=topology.total;r.nnz=topology.totalNz;
            int signs[2]={1,1},globalSigns[2]{};
            for(double v:matrix_.diag()){signs[0]&=v<0;signs[1]&=v>0;}
            MPI_Allreduce(signs,globalSigns,2,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
            if(!globalSigns[0]&&!globalSigns[1])throw std::runtime_error("AMG requires a uniform diagonal sign");
            const int sign=globalSigns[0]?-1:1;
            t=now();topology.update(matrix_,interfaceBouCoeffs_,interfaces_,sign);
            forAll(psi,i){topology.localRhs[topology.rowOrder[i]]=source[i];topology.localX[topology.rowOrder[i]]=psi[i];}
            r.conversion=maxWall(t);
            if(!r.structureReused && controlDict_.getOrDefault<bool>("auditAssembly",false)){
                if(topology.total>100000)throw std::runtime_error("Assembly audit is restricted to small validation cases");
                solveScalarField probe(psi.size()),expected(psi.size());
                forAll(probe,i)probe[i]=std::sin(0.37*(topology.offsets[rank]+i+1));
                matrix_.Amul(expected,probe,interfaceBouCoeffs_,interfaces_,cmpt);
                std::vector<double> globalProbe(topology.total),packedProbe(topology.n);
                forAll(probe,i)packedProbe[topology.rowOrder[i]]=probe[i];
                MPI_Allgatherv(packedProbe.data(),topology.n,MPI_DOUBLE,globalProbe.data(),topology.sizes.data(),topology.offsets.data(),MPI_DOUBLE,MPI_COMM_WORLD);
                long double q[2]={0,0},sum[2]{};
                forAll(probe,i){
                    long double value=0;
                    const int row=topology.rowOrder[i];
                    for(int k=topology.rows[row];k<topology.rows[row+1];++k)value+=topology.values[k]*globalProbe[topology.cols[k]];
                    const long double error=sign*value-expected[i];q[0]+=error*error;q[1]+=static_cast<long double>(expected[i])*expected[i];
                }
                MPI_Allreduce(q,sum,2,MPI_LONG_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
                r.assemblyAuditError=std::sqrt(sum[0]/std::max(sum[1],1e-300L));
                if(r.assemblyAuditError>1e-12)throw std::runtime_error("Global CSR matrix does not match OpenFOAM Amul");
            }
            t=now();
            MPI_Gatherv(topology.values.data(),topology.nz,MPI_DOUBLE,topology.globalValues.data(),topology.nnzs.data(),topology.nzOffsets.data(),MPI_DOUBLE,0,MPI_COMM_WORLD);
            MPI_Gatherv(topology.localRhs.data(),topology.n,MPI_DOUBLE,topology.b.data(),topology.sizes.data(),topology.offsets.data(),MPI_DOUBLE,0,MPI_COMM_WORLD);
            MPI_Gatherv(topology.localX.data(),topology.n,MPI_DOUBLE,topology.x.data(),topology.sizes.data(),topology.offsets.data(),MPI_DOUBLE,0,MPI_COMM_WORLD);
            if(rank==0&&sign<0)for(int i=0;i<topology.total;++i)topology.b[i]=-topology.b[i];
            r.gather=maxWall(t);
            if(rank==0)device.run(controlDict_,r);
            MPI_Bcast(&r.iterations,1,MPI_INT,0,MPI_COMM_WORLD);
            perf.nIterations()=r.iterations;perf.initialResidual()=r.initialFoam;
            t=now();MPI_Scatterv(topology.x.data(),topology.sizes.data(),topology.offsets.data(),MPI_DOUBLE,topology.localX.data(),topology.n,MPI_DOUBLE,0,MPI_COMM_WORLD);
            forAll(psi,i)psi[i]=topology.localX[topology.rowOrder[i]];
            r.scatter=maxWall(t);
        }else throw std::runtime_error("Unknown pressure backend");
        double t=now();
        solveScalarField prod(psi.size()),work(psi.size());
        matrix_.Amul(prod,psi,interfaceBouCoeffs_,interfaces_,cmpt);
        solveScalarField residual(source-prod);
        r.finalFoam=perf.finalResidual();
        long double local[3]={0,0,0},global[3]{};
        forAll(residual,i){
            if(!std::isfinite(residual[i])||!std::isfinite(psi[i]))throw std::runtime_error("Nonfinite result");
            local[0]+=static_cast<long double>(residual[i])*residual[i];
            local[1]+=static_cast<long double>(source[i])*source[i];local[2]+=std::abs(residual[i]);
        }
        MPI_Allreduce(local,global,3,MPI_LONG_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        r.finalL2=std::sqrt(global[0]);r.rhsL2=std::sqrt(global[1]);
        if(r.backend=="linglongamg"){
            r.finalFoam=double(global[2])/initialNorm;perf.finalResidual()=r.finalFoam;
            matrix_.setResidualField(residual,fieldName_,false);
        }
        r.foamConverged=perf.checkConvergence(tolerance_,relTol_,0);
        r.verify=maxWall(t);r.wall=maxWall(begin);active=nullptr;
        if(rank==0)records.push_back(r);
        return perf;
    }catch(const std::exception& e){active=nullptr;FatalErrorInFunction<<e.what()<<"; last SDK stage: "<<lastDeviceStage<<exit(FatalError);}
    return solverPerformance(typeName,fieldName_);
}
solverPerformance FoamPressure::scalarSolve(solveScalarField& psi,const solveScalarField& source,const direction cmpt)const{
    const char* requested=std::getenv("FOAM_PRESSURE_FROZEN_REPEATS");
    const int count=requested?std::atoi(requested):1;
    if(count<1||count>5)FatalErrorInFunction<<"Frozen replay count must be 1..5"<<exit(FatalError);
    const solveScalarField originalX0(psi);
    solverPerformance result(typeName,fieldName_);
    for(int i=1;i<=count;++i){
        psi=originalX0;currentRepeat=i;
        result=scalarSolveOnce(psi,source,cmpt);
    }
    return result;
}
solverPerformance FoamPressure::solve(scalarField& psi,const scalarField& source,const direction cmpt)const{
    PrecisionAdaptor<solveScalar,scalar> p(psi);
    return scalarSolve(p.ref(),ConstPrecisionAdaptor<solveScalar,scalar>(source)(),cmpt);
}
}
