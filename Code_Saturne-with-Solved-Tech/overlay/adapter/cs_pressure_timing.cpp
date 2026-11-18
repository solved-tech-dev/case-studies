/* Measure the unchanged native CPU pressure backend through ELF interposition.
 * Only pressure is recorded; every real implementation is called unchanged.
 */
#include "base/cs_defs.h"
#include "base/cs_time_step.h"
#include "bft/bft_error.h"
#include "bft/bft_printf.h"
#include "alge/cs_matrix.h"
#include "alge/cs_sles_it.h"
#include <dlfcn.h>
#include <chrono>
#include <cmath>
#include <cstring>
#include <unordered_map>
#include <vector>
namespace {
using Clock=std::chrono::steady_clock;
double milliseconds(Clock::time_point t) {return std::chrono::duration<double,std::milli>(Clock::now()-t).count();}
struct Trace {double total_setup=0,pending=0;unsigned long setups=0,solves=0;};
std::unordered_map<void*,Trace> traces;
bool pressure(const char *name) {return name && std::strcmp(name,"pressure")==0;}
template<typename F> F real(const char *name) {
  auto f=reinterpret_cast<F>(dlsym(RTLD_NEXT,name));
  if(!f)bft_error(__FILE__,__LINE__,0,"Cannot locate original %s",name);
  return f;
}
}
extern "C" void cs_sles_it_setup(void *context,const char *name,const cs_matrix_t *a,int verbosity) {
  static const auto original=real<decltype(&cs_sles_it_setup)>("cs_sles_it_setup");
  if(!pressure(name)){original(context,name,a,verbosity);return;}
  auto start=Clock::now();original(context,name,a,verbosity);
  auto &t=traces[context];double wall=milliseconds(start);t.total_setup+=wall;t.pending+=wall;++t.setups;
}
extern "C" cs_sles_convergence_state_t cs_sles_it_solve(void *context,const char *name,
 const cs_matrix_t *a,int verbosity,double precision,double rnorm,int *iterations,double *residual,
 const cs_real_t *rhs,cs_real_t *initial,cs_real_t *x,size_t aux_size,void *aux) {
  static const auto original=real<decltype(&cs_sles_it_solve)>("cs_sles_it_solve");
  if(!pressure(name))return original(context,name,a,verbosity,precision,rnorm,iterations,residual,rhs,initial,x,aux_size,aux);
  const double setup_before=traces[context].total_setup;
  const auto start=Clock::now();
  auto state=original(context,name,a,verbosity,precision,rnorm,iterations,residual,rhs,initial,x,aux_size,aux);
  const double call_ms=milliseconds(start);
  auto &t=traces[context];const double solve_ms=call_ms-(t.total_setup-setup_before);
  const auto verify=Clock::now();
  const auto n=cs_matrix_get_n_rows(a);
  std::vector<double> ax(cs_matrix_get_n_columns(a));
  cs_matrix_vector_multiply(a,x,ax.data());
  double norm=0;for(cs_lnum_t i=0;i<n;++i)norm=std::hypot(norm,ax[i]-rhs[i]);
  const double verify_ms=milliseconds(verify);
  const bool passed=std::isfinite(norm) && norm<=precision*rnorm;
  if(!passed && state==CS_SLES_CONVERGED)state=CS_SLES_DIVERGED;
  *residual=norm;
  const cs_lnum_t *row=nullptr,*col=nullptr;const double *d=nullptr,*v=nullptr;
  long long nnz=-1;
  if(cs_matrix_get_type(a)==CS_MATRIX_MSR) {
    cs_matrix_get_msr_arrays(a,&row,&col,&d,&v);nnz=static_cast<long long>(row[n])+n;
  } else if(cs_matrix_get_type(a)==CS_MATRIX_CSR) {
    cs_matrix_get_csr_arrays(a,&row,&col,&v);nnz=row[n];
  }
  bft_printf("CS_CPU_LINEAR route=native equation=%s step=%d solve_id=%lu n=%lld nnz=%lld "
             "setup_count=%lu setup_ms=%.6f solve_ms=%.6f call_ms=%.6f verify_ms=%.6f "
             "iterations=%d status=%d true_abs_residual=%.12e true_normalized_residual=%.12e "
             "tolerance=%.12e absolute_tolerance=%.12e true_residual_pass=%d\n",name,
             cs_glob_time_step?cs_glob_time_step->nt_cur:-1,++t.solves,static_cast<long long>(n),nnz,
             t.setups,t.pending,solve_ms,call_ms,verify_ms,*iterations,int(state),norm,
             rnorm>0?norm/rnorm:NAN,precision,precision*rnorm,int(passed));
  t.pending=0;return state;
}
