/* SPDX-License-Identifier: LGPL-2.1-or-later */
#include "bridge.h"
#include "csc.hpp"
#include <pangulu.h>
#include <cuda_runtime_api.h>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <stdexcept>
#include <string>
#if PANGULU_ABI_VERSION != 2
#error "This adapter targets PanguLU 26.44 public ABI 2"
#endif
#ifdef _OPENMP
#include <omp.h>
#endif
namespace {
using Clock=std::chrono::steady_clock;
struct Context {
    pangulu_handle_t h=nullptr; pangulu_matrix_t A=nullptr;
    pangulu_vector_t b=nullptr,x=nullptr; double *db=nullptr,*dx=nullptr;
    ~Context(){if(x)pangulu_vector_destroy(x);if(b)pangulu_vector_destroy(b);
      if(A)pangulu_matrix_destroy(A);if(h)pangulu_destroy(h);
      if(dx)cudaFree(dx);if(db)cudaFree(db);}
};
struct HostOMP {
#ifdef _OPENMP
 int threads=omp_get_max_threads(),dynamic=omp_get_dynamic();
 ~HostOMP(){omp_set_num_threads(threads);omp_set_dynamic(dynamic);}
#endif
};
void check(pangulu_status_t s,const char*api){if(s!=PANGULU_STATUS_SUCCESS)throw std::runtime_error(std::string(api)+" status="+std::to_string(int(s)));}
void cuda_check(cudaError_t s){if(s!=cudaSuccess)throw std::runtime_error(cudaGetErrorString(s));}
void capture(unsigned long long call,int64_t n,int64_t nz,const int64_t*rp,const int32_t*ci,const double*a,const double*b,const double*x){
 const char*dir=std::getenv("PHG_REPRO_LINEAR_DIR");if(!dir||!*dir)return;
 const std::string prefix=std::string(dir)+"/call-"+std::to_string(call);
 auto write=[&](const char*s,const void*data,size_t width,size_t count){
   auto*f=std::fopen((prefix+s).c_str(),"wbx");if(!f)throw std::runtime_error("linear capture create failed");
   bool ok=std::fwrite(data,width,count,f)==count;ok=(std::fclose(f)==0)&&ok;
   if(!ok)throw std::runtime_error("linear capture write failed");
 };
 if(!x){write(".rp.i64",rp,sizeof(*rp),n+1);write(".ci.i32",ci,sizeof(*ci),nz);
 write(".a.f64",a,sizeof(*a),nz);write(".b.f64",b,sizeof(*b),n);}
 else write(".x.f64",x,sizeof(*x),n);
}
double phase(Context&c,pangulu_phase_t p,const char*name,unsigned long long call){
 cuda_check(cudaDeviceSynchronize());auto t=Clock::now();
 auto status=pangulu_execute(c.h,p,c.A,p==PANGULU_PHASE_SOLVE?c.x:nullptr,p==PANGULU_PHASE_SOLVE?c.b:nullptr);
 auto sync=cudaDeviceSynchronize();double dt=std::chrono::duration<double>(Clock::now()-t).count();
 std::fprintf(stderr,"PHG_PANGULU_PHASE call=%llu phase=%s seconds=%.9f status=%d cuda_status=%d\n",call,name,dt,int(status),int(sync));
 check(status,name);cuda_check(sync);return dt;
}
}
extern "C" const char* phg_pangulu_backend_identity(){
#ifdef PHG_PANGULU_TEST_API
return "TEST_ONLY_PANGULU_API_MOCK_NOT_REAL_CUDA";
#else
return "PanguLU public API " PANGULU_VERSION_STRING " / CUDA";
#endif
}
extern "C" int phg_pangulu_solve_csr(int64_t n,int64_t nz,const int64_t*rp,const int32_t*ci,const double*a,const double*b,double*x){
 HostOMP restore; static unsigned long long ordinal=0;const auto call=++ordinal;
 try {
    if(std::strcmp(PANGULU_VERSION_STRING,"26.44"))throw std::runtime_error("PanguLU 26.44 headers required");
#ifdef _OPENMP
    if(const char *threads=std::getenv("PHG_PANGULU_THREADS")){
      int count=std::stoi(threads);
      if(count<=0)throw std::invalid_argument("PHG_PANGULU_THREADS must be positive");
      omp_set_dynamic(0);omp_set_num_threads(count);
    }
#endif
    if(!b||!x)throw std::invalid_argument("null vector");
    for(int64_t i=0;i<n;++i)if(!std::isfinite(b[i]))throw std::invalid_argument("nonfinite RHS");
    auto t=Clock::now();auto csc=repro::convert(n,nz,rp,ci,a);
    double conversion=std::chrono::duration<double>(Clock::now()-t).count();
    capture(call,n,nz,rp,ci,a,b,nullptr);
    Context c;check(pangulu_create(0,&c.h),"create");
    check(pangulu_config_set(c.h,PANGULU_REORDERING_METHOD,L1),"ordering");
    const char *reuse=std::getenv("PHG_PANGULU_REUSE_LEVEL");
    check(pangulu_config_set(c.h,PANGULU_REUSE_LEVEL,reuse && !std::strcmp(reuse,"L2") ? L2 : L0),"reuse");
    cuda_check(cudaMalloc((void**)&c.db,n*sizeof(double)));cuda_check(cudaMalloc((void**)&c.dx,n*sizeof(double)));
    // Match the vendor example: initialize device buffers before binding vectors.
    cuda_check(cudaMemcpy(c.db,b,n*sizeof(double),cudaMemcpyHostToDevice));cuda_check(cudaMemset(c.dx,0,n*sizeof(double)));
    check(pangulu_matrix_create(&c.A,PANGULU_MATRIX_LAYOUT_CSC,n,n,csc.value.size(),1,csc.ptr.data(),csc.row.data(),csc.value.data(),PANGULU_R_64I,PANGULU_R_32I,PANGULU_R_64F,PANGULU_LU,PANGULU_MEMORY_HOST),"matrix_create");
    check(pangulu_vector_create(&c.b,n,c.db,PANGULU_R_64F,PANGULU_MEMORY_DEVICE),"rhs_create");
    check(pangulu_vector_create(&c.x,n,c.dx,PANGULU_R_64F,PANGULU_MEMORY_DEVICE),"solution_create");
    double core=phase(c,PANGULU_PHASE_ANALYSIS,"analysis",call);
    core+=phase(c,PANGULU_PHASE_NUMERIC_FACTORISATION,"numeric_factorisation",call);
    core+=phase(c,PANGULU_PHASE_SOLVE,"solve",call);
    cuda_check(cudaMemcpy(x,c.dx,n*sizeof(double),cudaMemcpyDeviceToHost));
    double r=repro::residual(n,rp,ci,a,b,x);
    std::fprintf(stderr,"PHG_PANGULU_RESULT call=%llu n=%lld nnz=%lld conversion_s=%.9f core_s=%.9f relative_residual=%.17g\n",call,(long long)n,(long long)nz,conversion,core,r);
    if(!std::isfinite(r)||r>1e-8)throw std::runtime_error("original CSR residual exceeds 1e-8");
    capture(call,n,nz,rp,ci,a,b,x);
    return 0;
 }catch(const std::exception&e){std::fprintf(stderr,"PHG_PANGULU_ERROR call=%llu %s\n",call,e.what());return 1;}
}
