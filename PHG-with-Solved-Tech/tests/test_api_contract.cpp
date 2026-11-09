/* TEST ONLY: compile the production bridge against the real external public
 * header, with in-process API doubles. No real GPU or license is accessed. */
#include <pangulu.h>
#include <cuda_runtime_api.h>
#include "bridge.h"
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
struct pangulu_handle { int stage; };
struct pangulu_matrix {int64_t n;int64_t*cp;int32_t*ri;double*values;};
struct pangulu_vector {double*data;};
static int objects=0,allocations=0,phase_calls=0,fail=0;
extern "C" {
cudaError_t cudaDeviceSynchronize(){return 0;}
cudaError_t cudaMalloc(void**p,size_t n){*p=std::malloc(n);++allocations;return 0;}
cudaError_t cudaFree(void*p){if(p)--allocations;std::free(p);return 0;}
cudaError_t cudaMemcpy(void*d,const void*s,size_t n,cudaMemcpyKind){std::memcpy(d,s,n);return 0;}
cudaError_t cudaMemset(void*d,int v,size_t n){std::memset(d,v,n);return 0;}
const char*cudaGetErrorString(cudaError_t){return "test CUDA status";}
pangulu_status_t pangulu_create(int32_t gpu,pangulu_handle_t*h){assert(gpu==0);if(fail==4)return PANGULU_STATUS_LICENSE_ERROR;*h=new pangulu_handle{0};++objects;return PANGULU_STATUS_SUCCESS;}
pangulu_status_t pangulu_config_set(pangulu_handle_t,pangulu_config_option_t k,pangulu_config_value_t v){assert((k==PANGULU_REORDERING_METHOD&&v==L1)||(k==PANGULU_REUSE_LEVEL&&v==L0));return PANGULU_STATUS_SUCCESS;}
pangulu_status_t pangulu_matrix_create(pangulu_matrix_t*A,pangulu_matrix_layout_t layout,int64_t m,int64_t n,int64_t nz,int64_t bs,void*ptr,void*index,void*values,pangulu_index_type_t pt,pangulu_index_type_t it,pangulu_value_type_t vt,pangulu_operation_t op,pangulu_memory_space_t mem){
 assert(layout==PANGULU_MATRIX_LAYOUT_CSC&&m==n&&nz==4&&bs==1&&pt==PANGULU_R_64I&&it==PANGULU_R_32I&&vt==PANGULU_R_64F&&op==PANGULU_LU&&mem==PANGULU_MEMORY_HOST);
 *A=new pangulu_matrix{n,(int64_t*)ptr,(int32_t*)index,(double*)values};++objects;return PANGULU_STATUS_SUCCESS;
}
pangulu_status_t pangulu_vector_create(pangulu_vector_t*v,int64_t n,void*data,pangulu_value_type_t t,pangulu_memory_space_t mem){assert(n==2&&t==PANGULU_R_64F&&mem==PANGULU_MEMORY_DEVICE);*v=new pangulu_vector{(double*)data};++objects;return PANGULU_STATUS_SUCCESS;}
pangulu_status_t pangulu_execute(pangulu_handle_t h,pangulu_phase_t p,pangulu_matrix_t A,pangulu_vector_t x,pangulu_vector_t b){
 ++phase_calls;
 if(fail==2)return PANGULU_STATUS_ALLOC_FAILED;
 if(p==PANGULU_PHASE_ANALYSIS){assert(h->stage==0&&!x&&!b);h->stage=1;}
 else if(p==PANGULU_PHASE_NUMERIC_FACTORISATION){assert(h->stage==1&&!x&&!b);h->stage=2;}
 else {assert(p==PANGULU_PHASE_SOLVE&&h->stage==2&&x&&b);double a[2][2]={{0,0},{0,0}};
  for(int j=0;j<2;++j)for(int64_t k=A->cp[j];k<A->cp[j+1];++k)a[A->ri[k]][j]+=A->values[k];
  double det=a[0][0]*a[1][1]-a[0][1]*a[1][0];x->data[0]=(b->data[0]*a[1][1]-a[0][1]*b->data[1])/det;x->data[1]=(a[0][0]*b->data[1]-b->data[0]*a[1][0])/det;if(fail==5)x->data[0]+=1;h->stage=3;}
 return PANGULU_STATUS_SUCCESS;
}
pangulu_status_t pangulu_matrix_destroy(pangulu_matrix_t p){delete p;--objects;return PANGULU_STATUS_SUCCESS;}
pangulu_status_t pangulu_vector_destroy(pangulu_vector_t p){delete p;--objects;return PANGULU_STATUS_SUCCESS;}
pangulu_status_t pangulu_destroy(pangulu_handle_t p){delete p;--objects;return PANGULU_STATUS_SUCCESS;}
}
int main(){
 assert(std::strstr(phg_pangulu_backend_identity(),"TEST_ONLY"));
 int64_t rp[]={0,2,4};int32_t ci[]={0,1,0,1};double a[]={5,0,2,3},b[]={5,5},x[2]={0,0};
 assert(phg_pangulu_solve_csr(2,4,rp,ci,a,b,x)==0);assert(x[0]==1&&x[1]==1&&phase_calls==3&&objects==0&&allocations==0);
 fail=2;phase_calls=0;assert(phg_pangulu_solve_csr(2,4,rp,ci,a,b,x)!=0);assert(phase_calls==1&&objects==0&&allocations==0);
 fail=4;phase_calls=0;assert(phg_pangulu_solve_csr(2,4,rp,ci,a,b,x)!=0);assert(phase_calls==0&&objects==0&&allocations==0);
 fail=5;phase_calls=0;assert(phg_pangulu_solve_csr(2,4,rp,ci,a,b,x)!=0);assert(phase_calls==3&&objects==0&&allocations==0);
 std::puts("PASS production bridge with TEST_ONLY API doubles: phase order, CSR conversion, solve transfer, allocation failure cleanup, license failure propagation, wrong-solution rejection");
}
