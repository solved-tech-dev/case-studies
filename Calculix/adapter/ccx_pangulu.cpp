// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 解得科技
// PanguLU-only CalculiX 2.23 spooles-entry replacement.
// Scope: the five documented real symmetric steady scalar examples.
// Original assembly, RHS, constraints, fields and application driver are retained.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <vector>
#include <time.h>
#include <dlfcn.h>
#include <sched.h>
#include <omp.h>
#include <zlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <cuda_runtime_api.h>
#include <pangulu.h>

static double now(clockid_t c=CLOCK_MONOTONIC){timespec t;clock_gettime(c,&t);return t.tv_sec+1e-9*t.tv_nsec;}
static bool failed=false;
static void fail(const char *s,int code=1){failed=true;fprintf(stderr,"HULL_FAILURE %s code=%d\n",s,code);fflush(nullptr);exit(80);}
static void ck(int s,const char *name){if(s)fail(name,s);}
static void provider(const char *name,void *f){Dl_info d{};if(dladdr(f,&d))printf("HULL_PROVIDER %s %s\n",name,d.dli_fname);}
static bool gpu=false;
static void syncgpu(){
  if(gpu)ck(cudaDeviceSynchronize(),"device synchronize");
}
template<class F> static void timed(const char *layer,const char *phase,F f){syncgpu();double w=now(),c=now(CLOCK_PROCESS_CPUTIME_ID);f();syncgpu();printf("HULL_TIME %s %s wall=%.9f cpu=%.9f\n",layer,phase,now()-w,now(CLOCK_PROCESS_CPUTIME_ID)-c);fflush(stdout);}
static int n=0, calls=0;
static bool live=false,reuse_pattern=false;
static std::vector<int32_t> rp,ci;
static std::vector<double> av,original_b;
static double service_start,service_cpu;

// Evidence is outside service/native timing. It is captured for each call or,
// when explicitly requested, compared with an existing external reference.
// Capture of this run alone is not a cross-solver input-equivalence check.
template<class T> static void evidence_array(gzFile f,const std::vector<T>& a,bool write,const char *name){
  size_t bytes=a.size()*sizeof(T);if(bytes>INT32_MAX)fail("compressed evidence array exceeds bounded size");
  if(write){if(gzwrite(f,a.data(),bytes)!=(int)bytes)fail("canonical write");return;}
  std::vector<T> ref(a.size());if(gzread(f,ref.data(),bytes)!=(int)bytes)fail("canonical truncated");
  size_t unequal=0;long double maxdiff=0,ds=0,rs=0;
  for(size_t i=0;i<a.size();++i){if(a[i]!=ref[i])++unequal;long double d=(long double)a[i]-ref[i];maxdiff=std::max(maxdiff,fabsl(d));ds+=d*d;rs+=(long double)ref[i]*ref[i];}
  printf("HULL_INPUT %s count=%zu unequal=%zu max_abs=%.12Le relative_l2=%.12Le\n",name,a.size(),unequal,maxdiff,sqrtl(ds/(rs?rs:1)));
  if(unequal&&!(calls>1&&strcmp(name,"rhs")==0))fail("actual matrix/first RHS differs from canonical");
  if(unequal)printf("HULL_INPUT_LIMITATION later_Newton_RHS_differs_from_reference=1 no_identical_input_claim_for_this_call=1\n");
}
static void evidence(double *b){
  const char *base=getenv("HULL_CANONICAL");if(!base)fail("canonical path required");
  char pathname[4096];int len=snprintf(pathname,sizeof(pathname),"%s.call%03d",base,calls);if(len<0||len>=(int)sizeof(pathname))fail("canonical path length");const char *p=pathname;
  bool write=getenv("HULL_CREATE_CANONICAL")!=nullptr;
  int fd=write?open(p,O_WRONLY|O_CREAT|O_EXCL,0660):-1;
  gzFile f=write?(fd>=0?gzdopen(fd,"wb1"):nullptr):gzopen(p,"rb");if(!f)fail("canonical open (never overwrites)");
  int64_t header[3]={0x484c4c4353523031LL,n,(int64_t)ci.size()};
  if(write){if(gzwrite(f,header,sizeof(header))!=sizeof(header))fail("header write");}
  else {int64_t h[3];if(gzread(f,h,sizeof(h))!=sizeof(h)||memcmp(h,header,sizeof(h)))fail("canonical dimensions");}
  evidence_array(f,rp,write,"rowptr");evidence_array(f,ci,write,"columns");evidence_array(f,av,write,"coefficients");
  std::vector<double> rhs(b,b+n);evidence_array(f,rhs,write,"rhs");
  if(!write&&gzgetc(f)!=-1)fail("canonical trailing bytes");
  if(gzclose(f)!=Z_OK)fail("canonical close");
  printf("HULL_CANONICAL %s n=%d nnz=%zu full_CSR_zeros_retained=1\n",write?"created":"entire_arrays_compared",n,ci.size());
  const char *actual=getenv("HULL_ACTUAL_PREFIX");if(!actual)fail("actual input retention prefix required");
  len=snprintf(pathname,sizeof(pathname),"%s.call%03d.csr.gz",actual,calls);if(len<0||len>=(int)sizeof(pathname))fail("actual retention path length");
  fd=open(pathname,O_WRONLY|O_CREAT|O_EXCL,0660);f=fd>=0?gzdopen(fd,"wb1"):nullptr;if(!f)fail("actual input retention open");
  if(gzwrite(f,header,sizeof(header))!=sizeof(header))fail("actual header write");
  evidence_array(f,rp,true,"rowptr");evidence_array(f,ci,true,"columns");evidence_array(f,av,true,"coefficients");evidence_array(f,rhs,true,"rhs");
  if(gzclose(f)!=Z_OK)fail("actual input retention close");
  printf("HULL_ACTUAL_INPUT_RETAINED call=%d full_CSR_RHS=1\n",calls);
}
static void convert(double *ad,double *au,int *icol,int *irow,int neq,int nz){
  n=neq;rp.assign(n+1,0);int64_t p=0;
  for(int j=0;j<n;++j){++rp[j+1];for(int k=0;k<icol[j];++k,++p){int r=irow[p]-1;if(r<=j||r>=n)fail("not strict lower column storage");++rp[r+1];++rp[j+1];}}
  if(p!=nz||2*p+n>INT32_MAX)fail("native nnz/count range");
  for(int i=0;i<n;++i)rp[i+1]+=rp[i];
  ci.resize(rp[n]);av.resize(rp[n]);auto cursor=rp;
  auto add=[&](int r,int c,double x){if(!std::isfinite(x))fail("nonfinite coefficient");int k=cursor[r]++;ci[k]=c;av[k]=x;};
  for(int j=0;j<n;++j)add(j,j,ad[j]);
  p=0;
  for(int j=0;j<n;++j)for(int k=0;k<icol[j];++k,++p){int r=irow[p]-1;add(r,j,au[p]);add(j,r,au[p]);}
  std::vector<std::pair<int,double>> row;
  for(int i=0;i<n;++i){row.clear();for(int k=rp[i];k<rp[i+1];++k)row.emplace_back(ci[k],av[k]);std::sort(row.begin(),row.end());for(size_t k=0;k<row.size();++k){if(k&&row[k].first==row[k-1].first)fail("duplicate graph entry");ci[rp[i]+k]=row[k].first;av[rp[i]+k]=row[k].second;}}
}

static double *db=nullptr,*dx=nullptr,*dv=nullptr;static int32_t *dp=nullptr,*di=nullptr;
static void gpuinit(){gpu=true;ck(cudaSetDevice(0),"cudaSetDevice");cudaDeviceProp p{};ck(cudaGetDeviceProperties(&p,0),"device properties");printf("HULL_GPU name=%s bytes=%zu\n",p.name,p.totalGlobalMem);provider("cudaMalloc",dlsym(RTLD_DEFAULT,"cudaMalloc"));ck(cudaMalloc((void**)&db,n*sizeof(double)),"malloc b");ck(cudaMalloc((void**)&dx,n*sizeof(double)),"malloc x");ck(cudaMemset(db,0,n*sizeof(double)),"zero b");ck(cudaMemset(dx,0,n*sizeof(double)),"zero x");}
static pangulu_handle_t h=nullptr;static pangulu_matrix_t a=nullptr;static pangulu_vector_t bobj=nullptr,xobj=nullptr;
static std::vector<int64_t> cp;static std::vector<int32_t> ri;static std::vector<double> vv;
static std::vector<int32_t> csr_to_csc;
static void setup(){gpuinit();provider("pangulu_execute",(void*)pangulu_execute);ck(pangulu_create(0,&h),"pangulu_create license/native");ck(pangulu_config_set(h,PANGULU_REORDERING_METHOD,L0),"pangulu reorder L0");ck(pangulu_config_set(h,PANGULU_REUSE_LEVEL,L2),"pangulu reuse L2");cp.assign(n+1,0);for(auto j:ci)++cp[j+1];for(int i=0;i<n;++i)cp[i+1]+=cp[i];auto cur=cp;ri.resize(av.size());vv.resize(av.size());csr_to_csc.resize(av.size());for(int i=0;i<n;++i)for(int k=rp[i];k<rp[i+1];++k){int64_t q=cur[ci[k]]++;ri[q]=i;vv[q]=av[k];csr_to_csc[k]=q;}ck(pangulu_matrix_create(&a,PANGULU_MATRIX_LAYOUT_CSC,n,n,av.size(),1,cp.data(),ri.data(),vv.data(),PANGULU_R_64I,PANGULU_R_32I,PANGULU_R_64F,PANGULU_LU,PANGULU_MEMORY_HOST),"pangulu matrix");ck(pangulu_vector_create(&bobj,n,db,PANGULU_R_64F,PANGULU_MEMORY_DEVICE),"pangulu b");ck(pangulu_vector_create(&xobj,n,dx,PANGULU_R_64F,PANGULU_MEMORY_DEVICE),"pangulu x");printf("HULL_CONFIG PanguLU version=%s ABI=%d ordering=L0 reuse=L2 persistent_object_and_value_arrays=1\n",PANGULU_VERSION_STRING,PANGULU_ABI_VERSION);}
static void update(){for(size_t k=0;k<av.size();++k)vv[csr_to_csc[k]]=av[k];}
static void factor(){if(!reuse_pattern)timed("native","analysis_order_symbolic_preprocess",[]{ck(pangulu_execute(h,PANGULU_PHASE_ANALYSIS,a,nullptr,nullptr),"pangulu analysis");});timed("native",reuse_pattern?"refactor_numeric":"numeric",[]{ck(pangulu_execute(h,reuse_pattern?PANGULU_PHASE_REFACTORISATION:PANGULU_PHASE_NUMERIC_FACTORISATION,a,nullptr,nullptr),"pangulu numeric/refactor");});}
static void solve(double *b){timed("service_part","rhs_h2d",[&]{ck(cudaMemcpy(db,b,n*sizeof(double),cudaMemcpyHostToDevice),"copy rhs");});timed("native","solve",[]{ck(pangulu_execute(h,PANGULU_PHASE_SOLVE,a,xobj,bobj),"pangulu solve");});timed("service_part","solution_d2h",[&]{ck(cudaMemcpy(b,dx,n*sizeof(double),cudaMemcpyDeviceToHost),"copy solution");});}
static void destroy(){ck(pangulu_vector_destroy(xobj),"destroy x");ck(pangulu_vector_destroy(bobj),"destroy b");ck(pangulu_matrix_destroy(a),"destroy matrix");ck(pangulu_destroy(h),"destroy handle");ck(cudaFree(db),"free b");ck(cudaFree(dx),"free x");}

static void final_cleanup(){if(live&&!failed){timed("lifecycle","destroy",[]{destroy();});live=false;}}

extern "C" {
void spooles_factor(double *ad,double *au,double *,double *,double *sigma,int *icol,int *irow,int *neq,int *nzs,int *sym,int *fmt,int *){
  if(++calls>16||*sigma!=0||*sym!=0||*fmt!=0)fail("only bounded symmetric unshifted static/steady systems are in scope");
  printf("HULL_ENTRY call=%d n=%d lower_nnz=%d symmetryflag=%d inputformat=%d sigma=%.9g\n",calls,*neq,*nzs,*sym,*fmt,*sigma);
  cpu_set_t mask;CPU_ZERO(&mask);sched_getaffinity(0,sizeof(mask),&mask);printf("HULL_AFFINITY cpus=%d OMP=%s\n",CPU_COUNT(&mask),getenv("OMP_NUM_THREADS"));
  service_start=now();service_cpu=now(CLOCK_PROCESS_CPUTIME_ID);
  omp_set_num_threads(40);
  printf("HULL_GPU_HOST_OMP max_threads=%d\n",omp_get_max_threads());
  auto previous_rp=std::move(rp);auto previous_ci=std::move(ci);auto previous_av=std::move(av);int previous_n=n;
  timed("service_part","full_csr_conversion",[&]{convert(ad,au,icol,irow,*neq,*nzs);});
  timed("service_part","graph_and_coefficient_audit",[&]{reuse_pattern=live&&previous_n==n&&previous_rp==rp&&previous_ci==ci;size_t changed=0;if(reuse_pattern)for(size_t k=0;k<av.size();++k)changed+=(av[k]!=previous_av[k]);printf("HULL_REUSE call=%d same_entire_graph=%d native_analysis_skipped=%d coefficient_unequal=%zu all_coefficients_refreshed=1 numeric_redone=1\n",calls,reuse_pattern,reuse_pattern,changed);});
  timed("service_part","previous_capture_release",[&]{std::vector<int32_t>().swap(previous_rp);std::vector<int32_t>().swap(previous_ci);std::vector<double>().swap(previous_av);});
  if(live&&!reuse_pattern)final_cleanup();
  if(reuse_pattern)timed("service_part","all_native_values_update",[]{update();});
  else {timed("service_part","backend_setup_format_transfer",[]{setup();});live=true;if(calls==1)atexit(final_cleanup);}
  factor();
  printf("HULL_TIME service factor wall=%.9f cpu=%.9f\n",now()-service_start,now(CLOCK_PROCESS_CPUTIME_ID)-service_cpu);
}
void spooles_solve(double *b,int *neq){
  if(*neq!=n||!live)fail("solve lifecycle");
  if(original_b.size()==(size_t)n){size_t changed=0;for(int i=0;i<n;++i)changed+=(original_b[i]!=b[i]);printf("HULL_RHS_UPDATE call=%d rhs_unequal=%zu\n",calls,changed);}
  original_b.assign(b,b+n);
  double w=now(),c=now(CLOCK_PROCESS_CPUTIME_ID);solve(b);
  std::vector<double> correction(n);int extra_solves=0;
  for(int attempt=0;attempt<=2;++attempt){
    long double rel=0;
    timed("service_part","refinement_residual_long_double",[&]{long double rr=0,bb=0;for(int i=0;i<n;++i){long double r=original_b[i];for(int k=rp[i];k<rp[i+1];++k)r-=(long double)av[k]*b[ci[k]];correction[i]=(double)r;rr+=r*r;bb+=(long double)original_b[i]*original_b[i];}rel=sqrtl(rr/(bb?bb:1));});
    printf("HULL_REFINEMENT check=%d relative_l2=%.12Le\n",attempt,rel);
    if(!std::isfinite(rel))fail("nonfinite refinement residual");
    if(rel<=1e-8L||attempt==2)break;
    solve(correction.data());++extra_solves;
    for(int i=0;i<n;++i)b[i]+=correction[i];
  }
  syncgpu();printf("HULL_SOLVE_CALLS initial=1 extra_same_factor=%d\n",extra_solves);
  printf("HULL_TIME service solve wall=%.9f cpu=%.9f\n",now()-w,now(CLOCK_PROCESS_CPUTIME_ID)-c);
  timed("evidence","full_input_comparison",[&]{evidence(original_b.data());});
  timed("evidence","independent_original_CSR_residual",[&]{long double rr=0,bb=0;for(int i=0;i<n;++i){if(!std::isfinite(b[i]))fail("nonfinite x");long double r=-original_b[i];for(int k=rp[i];k<rp[i+1];++k)r+=(long double)av[k]*b[ci[k]];rr+=r*r;bb+=(long double)original_b[i]*original_b[i];}long double rel=sqrtl(rr/(bb?bb:1));printf("HULL_RESIDUAL relative_l2=%.12Le\n",rel);if(!std::isfinite(rel)||rel>1e-8L)fail("original CSR residual exceeds 1e-8");});
}
void spooles_cleanup(){printf("HULL_LIFECYCLE call=%d cleanup_deferred_until_app_exit_or_graph_change=1\n",calls);
  omp_set_num_threads(4);}
void spooles(double *ad,double *au,double *adb,double *aub,double *sigma,double *b,int *icol,int *irow,int *neq,int *nzs,int *sym,int *fmt,int *nzs3){spooles_factor(ad,au,adb,aub,sigma,icol,irow,neq,nzs,sym,fmt,nzs3);spooles_solve(b,neq);spooles_cleanup();}
void spooles_factor_rad(double*,double*,double*,double*,double*,int*,int*,int*,int*,int*,int*,int*){fail("radiation not in scope");}
void spooles_solve_rad(double*,int*){fail("radiation not in scope");}
void spooles_cleanup_rad(){fail("radiation not in scope");}
}
