/* SPDX-License-Identifier: LGPL-2.1-or-later */
#pragma once
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>
namespace repro {
struct CSC { std::vector<int64_t> ptr; std::vector<int32_t> row; std::vector<double> value; };
inline CSC convert(int64_t n, int64_t nnz, const int64_t *rp,
                   const int32_t *ci, const double *a) {
    if(n<=0 || n>INT32_MAX || nnz<0 || !rp || (nnz && (!ci || !a)) || rp[0]!=0 || rp[n]!=nnz)
        throw std::invalid_argument("invalid CSR dimensions/pointers");
    CSC c;
    c.ptr.assign(n+1,0);
    c.row.resize(nnz);
    c.value.resize(nnz);
    for(int64_t i=0;i<n;++i) {
        if(rp[i]<0 || rp[i+1]<rp[i] || rp[i+1]>nnz) throw std::invalid_argument("invalid CSR row pointer");
        for(int64_t k=rp[i];k<rp[i+1];++k) {
            if(ci[k]<0 || ci[k]>=n || !std::isfinite(a[k])) throw std::invalid_argument("invalid CSR coefficient/index");
            ++c.ptr[ci[k]+1];
        }
    }
    for(int64_t j=0;j<n;++j)c.ptr[j+1]+=c.ptr[j];
    auto cursor=c.ptr;
    for(int64_t i=0;i<n;++i)
        for(int64_t k=rp[i];k<rp[i+1];++k) {
            auto pos=cursor[ci[k]]++;
            c.row[pos]=static_cast<int32_t>(i);
            c.value[pos]=a[k];
        }
    return c;
}
inline double residual(int64_t n,const int64_t*rp,const int32_t*ci,const double*a,const double*b,const double*x) {
    long double r2=0,b2=0;
    for(int64_t i=0;i<n;++i) {
        if(!std::isfinite(b[i]) || !std::isfinite(x[i])) return std::numeric_limits<double>::infinity();
        long double r=-b[i];
        for(int64_t k=rp[i];k<rp[i+1];++k)r+=(long double)a[k]*x[ci[k]];
        r2+=r*r;b2+=(long double)b[i]*b[i];
    }
    return double(std::sqrt(r2)/(b2>0?std::sqrt(b2):1.0L));
}
}
