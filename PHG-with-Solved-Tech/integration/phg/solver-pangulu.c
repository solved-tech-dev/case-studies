/* SPDX-License-Identifier: LGPL-2.1-or-later
 * New PHG OEM backend, single-process real double only. No MUMPS impersonation. */
#include "phg.h"
#include "bridge.h"
#include <dlfcn.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
static int Create(SOLVER*s){
    if(s->mat->rmap->nprocs!=1 || s->mat->cmap->nprocs!=1 || sizeof(FLOAT)!=sizeof(double))
        phgError(1,"PanguLU backend supports one PHG process and double precision only.\n");
    return 0;
}
static int Destroy(SOLVER*s){phgMatDestroy(&s->mat);s->sent_to_oem=TRUE;return 0;}
static int Solve(SOLVER*s,VEC*x,BOOLEAN destroy){
    MAT*A=s->mat;INT i,k,n=A->rmap->nlocal;int64_t nz=0,*rp;int32_t*ci;double*av;
    const char*name=getenv("PHG_PANGULU_BRIDGE");static void*lib=NULL;
    int(*solve)(int64_t,int64_t,const int64_t*,const int32_t*,const double*,const double*,double*);
    const char*(*identity)(void);
    if(n<=0 || n!=A->rmap->nglobal || A->cmap->nglobal!=n || n>INT32_MAX)
        phgError(1,"PanguLU: unsupported matrix layout.\n");
    if(!name || !*name)phgError(1,"Set PHG_PANGULU_BRIDGE to the explicit adapter library path.\n");
    if(!lib)lib=dlopen(name,RTLD_NOW|RTLD_LOCAL);
    if(!lib)phgError(1,"PanguLU adapter load failed: %s\n",dlerror());
    *(void**)(&solve)=dlsym(lib,"phg_pangulu_solve_csr");
    *(void**)(&identity)=dlsym(lib,"phg_pangulu_backend_identity");
    if(!solve || !identity)phgError(1,"PanguLU adapter ABI symbols missing.\n");
    phgPrintf("PHG_REPRO_BACKEND %s\n",identity());
    rp=phgAlloc(((size_t)n+1)*sizeof(*rp));rp[0]=0;
    for(i=0;i<n;++i){const MAT_ROW*r=phgMatGetRow(A,i);nz+=r->ncols;rp[i+1]=nz;}
    ci=phgAlloc((size_t)nz*sizeof(*ci));av=phgAlloc((size_t)nz*sizeof(*av));
    for(i=0;i<n;++i){const MAT_ROW*r=phgMatGetRow(A,i);
        for(k=0;k<r->ncols;++k){if(r->cols[k]<0 || r->cols[k]>=n)phgError(1,"PanguLU: invalid column.\n");
            ci[rp[i]+k]=(int32_t)r->cols[k];av[rp[i]+k]=(double)r->data[k];}}
    if(solve(n,nz,rp,ci,av,s->rhs->data,x->data)!=0)phgError(1,"PanguLU solve failed. No fallback was attempted.\n");
    phgFree(rp);phgFree(ci);phgFree(av);
    if(destroy){phgMatDestroy(&s->mat);s->sent_to_oem=TRUE;}
    s->residual=0.; /* Actual original-system residual is checked by the bridge. */
    return s->nits=1;
}
OEM_SOLVER phgSolverPanguLU_={"PanguLU",NULL,NULL,NULL,NULL,Create,Destroy,
    NULL,NULL,NULL,NULL,Solve,NULL,NULL,NULL,M_UNSYM,FALSE,FALSE,FALSE,FALSE};
