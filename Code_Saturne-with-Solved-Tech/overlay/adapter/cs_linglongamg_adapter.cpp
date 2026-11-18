/* Code_Saturne 9.2: library-native GPU GMRES + persistent LinglongAMG.
 * Same graph: upload current CSR values, without setup/reset/recreation.
 * Verified against the 2026-11-16 SDK on A100; all solves check current A.
 */
#include "cs_headers.h"
#include "alge/cs_matrix_hypre_priv.h"
#include "linglongamg.h"
#include <cuda_runtime_api.h>
#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace {
using Clock=std::chrono::steady_clock;
double ms(Clock::time_point t) {return std::chrono::duration<double,std::milli>(Clock::now()-t).count();}
struct Solver {
  linglongamg_handle_t handle=nullptr;
  linglongamg_csr_t *matrix=nullptr;
  int *dr=nullptr,*dc=nullptr;
  double *da=nullptr,*db=nullptr,*dx=nullptr;
  std::vector<int> row,col;
  std::vector<double> values,reference_row_norm,row_scale,gpu_values,scaled_rhs;
  std::string config;
  int n=0,device=-1;
  unsigned long creates=0,builds=0,updates=0,solves=0,free_calls=0;
  bool ready=false,changed=false,build_pending=false,row_scaling=false;
  double setup_pending=0,convert_pending=0,copy_pending=0,create_pending=0,build_ms_pending=0,scaling_pending=0;
};
static_assert(std::is_same<cs_real_t,double>::value,"FP64 required");
static_assert(std::is_same<linglongamg_real_t,double>::value,"FP64 SDK required");
void cuda_ok(cudaError_t e,const char *where) {
  if(e!=cudaSuccess)bft_error(__FILE__,__LINE__,0,"LinglongAMG %s: %s",where,cudaGetErrorString(e));
}
void ll_ok(linglongamg_status_t e,const char *where) {
  if(e!=LINGLONGAMG_STATUS_SUCCESS)bft_error(__FILE__,__LINE__,0,"LinglongAMG %s: %s",where,linglongamg_status_string(e));
}
void hypre_ok(HYPRE_Int e,const char *where) {
  if(e)bft_error(__FILE__,__LINE__,0,"LinglongAMG HYPRE CSR extraction %s: %d",where,int(e));
}
void device_check(Solver *s) {
  int d=-1;cuda_ok(cudaGetDevice(&d),"get device");
  if(d!=s->device)bft_error(__FILE__,__LINE__,0,"Selected GPU changed while LinglongAMG handle is alive");
}
void release(Solver *s) {
  if(!s->handle)return;
  device_check(s);ll_ok(linglongamg_destroy(s->handle),"destroy");s->handle=nullptr;
  linglongamg_csr_destroy(s->matrix);s->matrix=nullptr;
  for(void *p:{static_cast<void*>(s->dr),static_cast<void*>(s->dc),static_cast<void*>(s->da),
               static_cast<void*>(s->db),static_cast<void*>(s->dx)})cuda_ok(cudaFree(p),"free");
  s->dr=s->dc=nullptr;s->da=s->db=s->dx=nullptr;
}
void setup(void *context,const char *name,const cs_matrix_t *a,int) {
  auto *s=static_cast<Solver*>(context);const auto start=Clock::now();
  if(cs_glob_n_ranks>1 || cs_matrix_get_diag_block_size(a)!=1)
    bft_error(__FILE__,__LINE__,0,"LinglongAMG GPU adapter requires one MPI rank and a scalar pressure matrix");
  if(std::string(cs_matrix_get_type_name(a))!=cs_matrix_hypre_ij_type_name)
    bft_error(__FILE__,__LINE__,0,"LinglongAMG GPU adapter expects host HYPRE ParCSR, including periodic identifications");
  const auto convert=Clock::now();
  const cs_lnum_t n=cs_matrix_get_n_rows(a);
  if(n<=0 || n>=INT_MAX)bft_error(__FILE__,__LINE__,0,"Unsupported matrix size");
  auto *coeffs=cs_matrix_hypre_get_coeffs(a);
  HYPRE_ParCSRMatrix par=nullptr;
  // HYPRE APIs return the process-wide error flag, including older errors.
  // Start this independent extraction with a clean flag and check each call.
  const int previous_error=HYPRE_GetError();
  if(previous_error) {
    bft_printf("CS_HYPRE_EXTRACTION previous_error=%d cleared_before_extraction=1\n",previous_error);
    HYPRE_ClearAllErrors();
  }
  hypre_ok(HYPRE_IJMatrixGetObject(coeffs->hm,reinterpret_cast<void**>(&par)),"get matrix");
  if(!par)bft_error(__FILE__,__LINE__,0,"Cannot access assembled pressure ParCSR matrix");
  std::vector<int> row(static_cast<size_t>(n)+1,0),col;
  std::vector<double> values;
  std::vector<std::pair<int,double>> sorted;
  const HYPRE_BigInt first=coeffs->l_range[0];
  for(int i=0;i<n;++i) {
    HYPRE_Int count=0;HYPRE_BigInt *columns=nullptr;HYPRE_Complex *v=nullptr;
    hypre_ok(HYPRE_ParCSRMatrixGetRow(par,first+i,&count,&columns,&v),"get row");
    sorted.clear();
    for(int j=0;j<count;++j) {
      const auto c=columns[j]-first;
      if(c<0 || c>=n || !std::isfinite(v[j]))bft_error(__FILE__,__LINE__,0,"Invalid pressure CSR column/value");
      sorted.emplace_back(static_cast<int>(c),v[j]);
    }
    hypre_ok(HYPRE_ParCSRMatrixRestoreRow(par,first+i,&count,&columns,&v),"restore row");
    std::sort(sorted.begin(),sorted.end());
    for(size_t j=0;j<sorted.size();) {
      const int c=sorted[j].first;double value=0;
      do {value+=sorted[j++].second;}while(j<sorted.size() && sorted[j].first==c);
      if(col.size()>=INT_MAX)bft_error(__FILE__,__LINE__,0,"CSR exceeds int32 API limit");
      col.push_back(c);values.push_back(value);
    }
    row[i+1]=static_cast<int>(col.size());
  }
  const double conversion=ms(convert);
  const bool rebuild=!s->handle || row!=s->row || col!=s->col;
  s->changed=s->handle && values!=s->values;
  const auto scaling_start=Clock::now();
  const std::vector<double> *upload_values=&values;
  double minimum_scale=1,maximum_scale=1;
  if(s->row_scaling) {
    minimum_scale=HUGE_VAL;maximum_scale=0;
    if(rebuild) {s->reference_row_norm.resize(n);s->row_scale.resize(n);s->scaled_rhs.resize(n);}
    s->gpu_values.resize(values.size());
    for(int i=0;i<n;++i) {
      double norm=0;for(int j=row[i];j<row[i+1];++j)norm+=std::abs(values[j]);
      if(!std::isfinite(norm) || norm<=0)bft_error(__FILE__,__LINE__,0,"Invalid pressure row norm");
      if(rebuild)s->reference_row_norm[i]=norm;
      s->row_scale[i]=s->reference_row_norm[i]/norm;
      minimum_scale=std::min(minimum_scale,s->row_scale[i]);maximum_scale=std::max(maximum_scale,s->row_scale[i]);
      for(int j=row[i];j<row[i+1];++j)s->gpu_values[j]=values[j]*s->row_scale[i];
    }
    upload_values=&s->gpu_values;
  }
  const double scaling_ms=ms(scaling_start);
  double creation=0,build=0;
  const auto copy=Clock::now();
  if(rebuild) {
    release(s);s->n=static_cast<int>(n);s->row=std::move(row);s->col=std::move(col);
    cuda_ok(cudaGetDevice(&s->device),"get device");
    cudaDeviceProp properties{};cuda_ok(cudaGetDeviceProperties(&properties,s->device),"device properties");
    bft_printf("CS_LINGLONGAMG_DEVICE ordinal=%d name=\"%s\" compute=%d.%d\n",s->device,properties.name,properties.major,properties.minor);
    const auto creation_start=Clock::now();
    ll_ok(linglongamg_create_from_config(&s->handle,s->config.c_str()),"create_from_config");
    creation=ms(creation_start);++s->creates;
    cuda_ok(cudaMalloc(reinterpret_cast<void**>(&s->dr),s->row.size()*sizeof(int)),"allocate rows");
    cuda_ok(cudaMalloc(reinterpret_cast<void**>(&s->dc),s->col.size()*sizeof(int)),"allocate columns");
    cuda_ok(cudaMalloc(reinterpret_cast<void**>(&s->da),values.size()*sizeof(double)),"allocate values");
    cuda_ok(cudaMalloc(reinterpret_cast<void**>(&s->db),size_t(s->n)*sizeof(double)),"allocate rhs");
    cuda_ok(cudaMalloc(reinterpret_cast<void**>(&s->dx),size_t(s->n)*sizeof(double)),"allocate x");
    cuda_ok(cudaMemcpy(s->dr,s->row.data(),s->row.size()*sizeof(int),cudaMemcpyHostToDevice),"upload rows");
    cuda_ok(cudaMemcpy(s->dc,s->col.data(),s->col.size()*sizeof(int),cudaMemcpyHostToDevice),"upload columns");
    cuda_ok(cudaMemcpy(s->da,upload_values->data(),values.size()*sizeof(double),cudaMemcpyHostToDevice),"upload values");
    ll_ok(linglongamg_read_csr(&s->matrix,s->n,s->n,int(s->col.size()),s->dr,s->dc,s->da),"read_csr");
  } else {
    device_check(s);
    // The actual SDK uses these borrowed device values as GMRES' current A.
    // No setup/reset: keep AMG hierarchy and persistent work buffers.
    cuda_ok(cudaMemcpy(s->da,upload_values->data(),values.size()*sizeof(double),cudaMemcpyHostToDevice),"update current matrix");
    ++s->updates;
  }
  cuda_ok(cudaDeviceSynchronize(),"matrix transfer synchronize");
  const double matrix_copy=ms(copy)-creation;
  if(rebuild) {
    const auto b=Clock::now();ll_ok(linglongamg_setup(s->handle,s->matrix),"setup");
    cuda_ok(cudaDeviceSynchronize(),"setup synchronize");build=ms(b);++s->builds;s->build_pending=true;
  }
  s->values=std::move(values);s->ready=true;
  const double wall=ms(start);
  s->setup_pending+=wall;s->convert_pending+=conversion;s->copy_pending+=matrix_copy;
  s->create_pending+=creation;s->build_ms_pending+=build;
  s->scaling_pending+=scaling_ms;
  bft_printf("CS_LINGLONGAMG_SETUP equation=%s handle=%p creates=%lu hierarchy_builds=%lu matrix_updates=%lu "
             "rebuilt=%d values_changed=%d n=%d nnz=%zu setup_wall_ms=%.6f conversion_ms=%.6f "
             "matrix_upload_ms=%.6f create_ms=%.6f hierarchy_setup_ms=%.6f row_scaled=%d row_scaling_ms=%.6f "
             "minimum_row_scale=%.9e maximum_row_scale=%.9e\n",name,static_cast<void*>(s->handle),
             s->creates,s->builds,s->updates,int(rebuild),int(s->changed),s->n,s->col.size(),wall,
             conversion,matrix_copy,creation,build,int(s->row_scaling),scaling_ms,minimum_scale,maximum_scale);
}
double true_residual(Solver *s,const double *rhs,const double *x) {
  double sum=0;
  // Use a stable norm and the CURRENT host coefficients, not setup's matrix.
  for(int i=0;i<s->n;++i) {
    double ax=0;for(int j=s->row[i];j<s->row[i+1];++j)ax+=s->values[j]*x[s->col[j]];
    sum=std::hypot(sum,ax-rhs[i]);
  }
  return sum;
}
cs_sles_convergence_state_t solve(void *context,const char *name,const cs_matrix_t *a,int verbosity,
 double precision,double rnorm,int *iterations,double *residual,const cs_real_t *rhs,
 cs_real_t *initial,cs_real_t *x,size_t,void*) {
  auto *s=static_cast<Solver*>(context);const auto start=Clock::now();
  if(!s->ready)setup(s,name,a,verbosity);
  device_check(s);
  const size_t bytes=size_t(s->n)*sizeof(double);
  const auto upload=Clock::now();
  const double *upload_rhs=rhs;
  if(s->row_scaling) {
    for(int i=0;i<s->n;++i)s->scaled_rhs[i]=rhs[i]*s->row_scale[i];
    upload_rhs=s->scaled_rhs.data();
  }
  cuda_ok(cudaMemcpy(s->db,upload_rhs,bytes,cudaMemcpyHostToDevice),"upload rhs");
  if(initial)cuda_ok(cudaMemcpy(s->dx,initial,bytes,cudaMemcpyHostToDevice),"upload initial x");
  else cuda_ok(cudaMemset(s->dx,0,bytes),"zero initial x");
  cuda_ok(cudaDeviceSynchronize(),"vector upload synchronize");const double upload_ms=ms(upload);
  linglongamg_vec_t b{s->n,s->db},v{s->n,s->dx};
  const auto kernel=Clock::now();ll_ok(linglongamg_solve(s->handle,&b,&v),"solve");
  cuda_ok(cudaDeviceSynchronize(),"solve synchronize");const double solve_ms=ms(kernel);
  const auto query=Clock::now();linglongamg_solve_info_t info{};
  ll_ok(linglongamg_get_solve_info(s->handle,&info),"get_solve_info");
  cuda_ok(cudaDeviceSynchronize(),"solution commit synchronize");const double query_ms=ms(query);
  const auto download=Clock::now();cuda_ok(cudaMemcpy(x,s->dx,bytes,cudaMemcpyDeviceToHost),"download x");
  const double download_ms=ms(download);
  const auto verify=Clock::now();*residual=true_residual(s,rhs,x);const double verify_ms=ms(verify);
  *iterations=info.iterations;
  const double target=precision*rnorm;
  const auto state=std::isfinite(*residual) && *residual<=target ? CS_SLES_CONVERGED :
    (info.converged ? CS_SLES_DIVERGED : CS_SLES_MAX_ITERATION);
  bft_printf("CS_GPU_LINEAR route=linglongamg equation=%s step=%d solve_id=%lu handle=%p creates=%lu "
             "hierarchy_builds=%lu matrix_updates=%lu setup_rebuilt=%d n=%d nnz=%zu "
             "setup_wall_ms=%.6f conversion_ms=%.6f matrix_upload_ms=%.6f create_ms=%.6f "
             "hierarchy_setup_ms=%.6f row_scaling_ms=%.6f row_scaled=%d vector_upload_ms=%.6f solve_ms=%.6f "
             "sdk_setup_ms=%.6f sdk_solve_ms=%.6f query_commit_ms=%.6f download_ms=%.6f "
             "verify_ms=%.6f call_ms=%.6f iterations=%d status=%d sdk_converged=%d "
             "true_abs_residual=%.12e true_normalized_residual=%.12e tolerance=%.12e "
             "absolute_tolerance=%.12e\n",name,cs_glob_time_step?cs_glob_time_step->nt_cur:-1,++s->solves,
             static_cast<void*>(s->handle),s->creates,s->builds,s->updates,int(s->build_pending),s->n,s->col.size(),
             s->setup_pending,s->convert_pending,s->copy_pending,s->create_pending,s->build_ms_pending,
             s->scaling_pending,int(s->row_scaling),upload_ms,solve_ms,s->build_pending?info.setup_ms:0.,info.solve_ms,query_ms,download_ms,
             verify_ms,ms(start),*iterations,int(state),info.converged,*residual,
             rnorm>0?*residual/rnorm:NAN,precision,target);
  s->setup_pending=s->convert_pending=s->copy_pending=s->create_pending=s->build_ms_pending=0;
  s->scaling_pending=0;
  s->build_pending=false;s->changed=false;
  if(state!=CS_SLES_CONVERGED)
    bft_error(__FILE__,__LINE__,0,"LinglongAMG current-matrix residual failed: %.12e > %.12e; do not advance CFD with an invalid pressure solve",*residual,target);
  return state;
}
void free_context(void *context) {auto *s=static_cast<Solver*>(context);s->ready=false;++s->free_calls;}
void log(const void *,cs_log_t type) {
  cs_log_printf(type,"  Solver: LinglongAMG library-native GPU GMRES + AMG, persistent hierarchy\n");
}
void destroy(void **context) {
  auto *s=static_cast<Solver*>(*context);if(!s)return;
  bft_printf("CS_LINGLONGAMG_SUMMARY solves=%lu creates=%lu hierarchy_builds=%lu matrix_updates=%lu free_calls=%lu\n",
             s->solves,s->creates,s->builds,s->updates,s->free_calls);
  release(s);delete s;*context=nullptr;
}
}
extern "C" void cs_linglongamg_define(int field_id) {
  const char *config=std::getenv("CS_LINGLONGAMG_CONFIG");
  if(!config || !*config)bft_error(__FILE__,__LINE__,0,"CS_LINGLONGAMG_CONFIG is required");
  auto *s=new Solver;s->config=config;
  const char *scaling=std::getenv("CS_LINGLONGAMG_ROW_SCALE");
  s->row_scaling=scaling && std::string(scaling)=="1";
  #ifndef CS_SLES_PRESSURE_BACKEND_TYPE
#define CS_SLES_PRESSURE_BACKEND_TYPE "cs_sles_linglongamg_t"
#endif
  // Match the host registration type that selects ParCSR.
  cs_sles_define(field_id,nullptr,s,CS_SLES_PRESSURE_BACKEND_TYPE,setup,solve,free_context,log,nullptr,destroy);
}
