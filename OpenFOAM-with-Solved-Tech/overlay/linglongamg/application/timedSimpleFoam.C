// SPDX-License-Identifier: GPL-3.0-or-later
// Upstream OpenFOAM v2412 simpleFoam equations with memory-buffered timing.
#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "fvOptions.H"
#include "FoamPressure.H"
#include <mpi.h>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <cstdio>
#include <iomanip>
#include <string>
#include <vector>
#include <streambuf>
extern "C" void benchIoStart();
extern "C" unsigned long benchIoStop();
extern "C" const char* benchIoFirstPath();
// The signed SDK prints its setup grid even at print_level=0. Keep C stdout
// in an anonymous memory stream and discard C++ cout during the compute loop.
// No descriptor redirection, /dev/null writes, or regular output files are used.
class QuietBuffer : public std::streambuf {
public:
    std::size_t bytes=0;
    std::streamsize xsputn(const char*,std::streamsize n) override {bytes+=n;return n;}
    int overflow(int c) override {if(c!=traits_type::eof())++bytes;return traits_type::not_eof(c);}
};
struct Phase {double outer,wall;std::string name;};
static std::vector<Phase> phases;
static double step=0;
static double stamp(){return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();}
static void phase(const char* name,double begin){
    const double local=stamp()-begin;double maximum=0;
    MPI_Allreduce(&local,&maximum,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    if(Foam::Pstream::master())phases.push_back({step,maximum,name});
    if(Foam::Pstream::master()&&std::getenv("FOAM_PRESSURE_TRACE")){
        std::fprintf(stderr,"FP_PHASE %.0f %s %.9f\n",step,name,maximum);std::fflush(stderr);
    }
}
template<class Field> static void summary(std::ofstream& out,const char* name,const Field& field){
    long double q=0,sum=0;double low=Foam::GREAT,high=-Foam::GREAT;
    for(const auto& v:field){double m=Foam::mag(v);q+=static_cast<long double>(m)*m;sum+=m;low=std::min(low,m);high=std::max(high,m);}
    long double data[2]={q,sum},total[2]{};double lo=0,hi=0;
    MPI_Allreduce(data,total,2,MPI_LONG_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&low,&lo,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    MPI_Allreduce(&high,&hi,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    if(Foam::Pstream::master())out<<name<<','<<std::sqrt(total[0])<<','<<total[1]<<','<<lo<<','<<hi<<'\n';
}
// Small validation cases only: retain all internal and patch values for a
// pointwise restart comparison. Formal runs leave this diagnostic disabled.
template<class Field> static void dump(const char* destination,const char* name,const Field& field){
    const std::string path=std::string(destination)+"/"+name+"-rank"+std::to_string(Foam::Pstream::myProcNo())+".f64";
    std::ofstream out(path,std::ios::binary);
    const auto& internal=field.primitiveField();
    if(internal.size())out.write(reinterpret_cast<const char*>(internal.cdata()),internal.size()*sizeof(internal[0]));
    for(const auto& patch:field.boundaryField())
        if(patch.size())out.write(reinterpret_cast<const char*>(patch.cdata()),patch.size()*sizeof(patch[0]));
    if(!out)FatalErrorInFunction<<"Cannot write diagnostic field values"<<Foam::exit(Foam::FatalError);
}
int main(int argc,char* argv[]){
    const double programBegin=stamp();
    const char* destination=std::getenv("FOAM_PRESSURE_RESULTS");
    if(!destination){std::cerr<<"FOAM_PRESSURE_RESULTS is required\n";return 2;}
    bool master=false;double beforeCleanup=0;unsigned executed=0;double startTime=0;
    const char* checkpointText=std::getenv("FOAM_BENCH_CHECKPOINT_STEP");
    const double checkpoint=checkpointText?std::stod(checkpointText):-1;
    const bool stopAfterCheckpoint=std::getenv("FOAM_BENCH_STOP_AFTER_CHECKPOINT");
    {
    // Suppress ordinary OpenFOAM log output; fatal errors retain their channel.
    Foam::messageStream::level=0;
    argList::addNote("OpenFOAM simpleFoam with pressure-only profiling");
    #include "postProcess.H"
    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    master=Pstream::master();
    #include "createTime.H"
    startTime=runTime.value();
    #include "createDynamicFvMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"
    turbulence->validate();
    phase("initialization",programBegin);
    char* quietData=nullptr;std::size_t quietBytes=0;
    FILE* quietStream=open_memstream(&quietData,&quietBytes);
    if(!quietStream)FatalErrorInFunction<<"Cannot create memory stdout sink"<<exit(FatalError);
    std::fflush(stdout);FILE* savedStdout=stdout;stdout=quietStream;
    QuietBuffer quietCout;std::streambuf* savedCout=std::cout.rdbuf(&quietCout);
    benchIoStart();
    const double loopBegin=stamp();
    while(simple.loop()){
        step=runTime.value();++executed;fpOuter(step);const double iterationBegin=stamp();
        double t=stamp();mesh.controlledUpdate();if(mesh.changing())MRF.update();phase("mesh_update",t);
        {
            t=stamp();
            #include "UEqn.H"
            phase("momentum",t);
            t=stamp();
            #include "pEqn.H"
            phase("pressure_total",t);
        }
        t=stamp();laminarTransport.correct();turbulence->correct();phase("turbulence",t);
        // Field writes and all function objects are disabled by construction.
        phase("outer_iteration",iterationBegin);
        // An explicit environment setting enables this only for approved
        // continuation runs or small restart validation cases.
        if(step==checkpoint){
            t=stamp();
            if(!runTime.writeNow())FatalErrorInFunction<<"Checkpoint write failed"<<exit(FatalError);
            phase("checkpoint_write",t);
            if(stopAfterCheckpoint)break;
        }
    }
    phase("iteration_loop",loopBegin);
    std::fflush(stdout);std::fclose(quietStream);stdout=savedStdout;
    std::cout.rdbuf(savedCout);
    const auto rejectedWrites=benchIoStop();
    unsigned long globalRejected=0;
    MPI_Allreduce(&rejectedWrites,&globalRejected,1,MPI_UNSIGNED_LONG,MPI_SUM,MPI_COMM_WORLD);
    if(globalRejected){std::cerr<<"BENCH_IO_REJECTED count="<<globalRejected<<" first="<<benchIoFirstPath()<<std::endl;MPI_Abort(MPI_COMM_WORLD,82);}
    if(master){std::ofstream audit(std::string(destination)+"/compute-io.json");audit<<"{\"status\":\"PASS\",\"rejected_regular_file_opens\":"<<globalRejected<<",\"CUDA_CACHE_DISABLE\":true,\"scope\":\"outer iteration loop including first AMG setup\"}\n";}
    if(master){std::ofstream audit(std::string(destination)+"/compute-stdout.json");audit<<"{\"status\":\"PASS\",\"c_stdout_memory_bytes\":"<<quietBytes<<",\"cpp_stdout_discarded_bytes\":"<<quietCout.bytes<<",\"output_descriptor_writes\":false}\n";}
    if(master){std::ofstream sdk(std::string(destination)+"/sdk-memory-output.log",std::ios::binary);sdk.write(quietData,quietBytes);}
    std::free(quietData);
    const double check=stamp();
    std::ofstream f;
    if(master){f.open(std::string(destination)+"/field-summary.csv");f<<std::setprecision(17)<<"field,l2,sum_magnitude,min_magnitude,max_magnitude\n";}
    summary(f,"p",p.primitiveField());summary(f,"U",U.primitiveField());
    summary(f,"k",turbulence->k()().primitiveField());summary(f,"omega",turbulence->omega()().primitiveField());
    if(std::getenv("FOAM_BENCH_DIAGNOSTIC_FIELD_DUMP")){
        long long local=mesh.nCells(),total=0;MPI_Allreduce(&local,&total,1,MPI_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);
        if(total>100000)FatalErrorInFunction<<"Pointwise dump is restricted to small validation cases"<<exit(FatalError);
        dump(destination,"p",p);dump(destination,"U",U);dump(destination,"phi",phi);
        dump(destination,"k",turbulence->k()());dump(destination,"omega",turbulence->omega()());dump(destination,"nut",turbulence->nut()());
    }
    phase("final_field_summary",check);
    double t=stamp();fpShutdown();phase("sdk_shutdown",t);
    phase("main_before_output_and_cleanup",programBegin);
    fpFlush(destination);
    if(master){
        std::ofstream out(std::string(destination)+"/phases.csv");out<<std::setprecision(17)<<"outer,stage,wall_s\n";
        for(const auto& p:phases)out<<p.outer<<','<<p.name<<','<<p.wall<<'\n';
    }
    beforeCleanup=stamp();
    }
    if(master){std::ofstream f(std::string(destination)+"/main.json");f<<std::setprecision(17)<<"{\"completed_steps\":"<<executed<<",\"last_outer\":"<<step<<",\"start_time\":"<<startTime<<",\"main_wall_s\":"<<stamp()-programBegin<<",\"cleanup_wall_s\":"<<stamp()-beforeCleanup<<"}\n";}
    return 0;
}
