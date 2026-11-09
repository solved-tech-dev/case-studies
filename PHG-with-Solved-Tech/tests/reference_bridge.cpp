/* SPDX-License-Identifier: LGPL-2.1-or-later
 * TEST ONLY. Tiny dense pivoted reference for validating the PHG connection.
 * This does not invoke, emulate, or benchmark the PanguLU runtime. */
#include "bridge.h"
#include "csc.hpp"
#include <cmath>
#include <cstdio>
#include <exception>
#include <utility>
#include <vector>
extern "C" const char*phg_pangulu_backend_identity(){return "TEST_ONLY_CPU_REFERENCE_NOT_PANGULU";}
extern "C" int phg_pangulu_solve_csr(int64_t n,int64_t nz,const int64_t*rp,const int32_t*ci,const double*a,const double*b,double*x){
 try{if(n>2000)throw std::runtime_error("test reference order limit is 2000");
    auto c=repro::convert(n,nz,rp,ci,a);std::vector<double>A(n*n,0),rhs(b,b+n);
    for(int64_t j=0;j<n;++j)for(int64_t k=c.ptr[j];k<c.ptr[j+1];++k)A[c.row[k]*n+j]+=c.value[k];
    for(int64_t j=0;j<n;++j){int64_t p=j;for(int64_t i=j+1;i<n;++i)if(std::abs(A[i*n+j])>std::abs(A[p*n+j]))p=i;
        if(A[p*n+j]==0)throw std::runtime_error("singular test system");
        if(p!=j){for(int64_t k=j;k<n;++k)std::swap(A[p*n+k],A[j*n+k]);std::swap(rhs[p],rhs[j]);}
        for(int64_t i=j+1;i<n;++i){double v=A[i*n+j]/A[j*n+j];for(int64_t k=j+1;k<n;++k)A[i*n+k]-=v*A[j*n+k];rhs[i]-=v*rhs[j];}}
    for(int64_t i=n;i-->0;){double v=rhs[i];for(int64_t j=i+1;j<n;++j)v-=A[i*n+j]*x[j];x[i]=v/A[i*n+i];}
    double r=repro::residual(n,rp,ci,a,b,x);std::fprintf(stderr,"TEST_ONLY_REFERENCE n=%lld relative_residual=%.17g\n",(long long)n,r);
    return !std::isfinite(r)||r>1e-8;
 }catch(const std::exception&e){std::fprintf(stderr,"TEST_ONLY_REFERENCE_ERROR %s\n",e.what());return 1;}
}
