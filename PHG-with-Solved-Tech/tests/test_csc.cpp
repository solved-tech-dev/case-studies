/* SPDX-License-Identifier: LGPL-2.1-or-later */
#include "csc.hpp"
#include <cassert>
#include <iostream>
int main(){
 int64_t rp[]={0,4,6};int32_t ci[]={1,0,1,1,1,0};double a[]={3,2,-1,-2,4,1};
 auto c=repro::convert(2,6,rp,ci,a);
 assert((c.ptr==std::vector<int64_t>{0,2,6}));assert((c.row==std::vector<int32_t>{0,1,0,0,0,1}));assert((c.value==std::vector<double>{2,1,3,-1,-2,4}));
 a[0]=0.;assert(repro::convert(2,6,rp,ci,a).value[2]==0.);a[0]=3.;
 double b[]={2,5},x[]={1,1};assert(repro::residual(2,rp,ci,a,b,x)<1e-15);
 bool caught=false;try{int32_t bad[]={2};int64_t p[]={0,1,1};double v[]={1};repro::convert(2,1,p,bad,v);}catch(const std::invalid_argument&){caught=true;}assert(caught);
 std::cout<<"PASS CSR/CSC, stored-entry preservation, residual, invalid indices\n";
}
