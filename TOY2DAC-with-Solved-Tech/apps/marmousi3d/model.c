// SPDX-License-Identifier: AGPL-3.0-only
// Copyright 2026 Solved Tech
// Marmousi-derived 3D elastic model generation and matrix assembly.
#define _POSIX_C_SOURCE 200809L
#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "model.h"
#include "grid.h"
typedef struct { int nx,ny,nz; double h,freq; const char *model; } Opt;
typedef struct { int *r,*c; double *ar,*ai; size_t nnz,cap; } Coo;
typedef struct { int nx,ny,nz; size_t nodes; double *vp,*vs,*rho,*la,*mu; } Model;

static void *xm(size_t n){void*p=malloc(n);if(!p){fprintf(stderr,"alloc %zu failed\n",n);exit(2);}return p;}
static size_t node(int i,int j,int k,int nx,int ny){return ((size_t)k*ny+j)*nx+i;}
static int dof(size_t q,int c){return (int)(3*q+c);}
static void add(Coo*a,int r,int c,double x,double y){if(a->nnz>=a->cap){fprintf(stderr,"COO overflow\n");exit(2);}a->r[a->nnz]=r;a->c[a->nnz]=c;a->ar[a->nnz]=x;a->ai[a->nnz]=y;a->nnz++;}

static float interp(const float*a,double x,double z){
 const int nx=681,nz=141;if(x<0)x=0;if(z<0)z=0;if(x>nx-1)x=nx-1;if(z>nz-1)z=nz-1;
 int x0=(int)floor(x),z0=(int)floor(z),x1=x0+1<nx?x0+1:x0,z1=z0+1<nz?z0+1:z0;
 double tx=x-x0,tz=z-z0,v00=a[(size_t)x0*nz+z0],v10=a[(size_t)x1*nz+z0],v01=a[(size_t)x0*nz+z1],v11=a[(size_t)x1*nz+z1];
 return (float)((1-tz)*((1-tx)*v00+tx*v10)+tz*((1-tx)*v01+tx*v11));
}

static int model_load(const Opt*o,Model*m){
 FILE*f=fopen(o->model,"rb");if(!f){fprintf(stderr,"open %s: %s\n",o->model,strerror(errno));return -1;}
 float*raw=xm(681U*141U*sizeof(float));size_t got=fread(raw,sizeof(float),681U*141U,f);fclose(f);
 if(got!=681U*141U){fprintf(stderr,"bad model size %zu\n",got);free(raw);return -1;}
 m->nx=o->nx+1;m->ny=o->ny+1;m->nz=o->nz+1;m->nodes=(size_t)m->nx*m->ny*m->nz;
 m->vp=xm(m->nodes*sizeof(double));m->vs=xm(m->nodes*sizeof(double));m->rho=xm(m->nodes*sizeof(double));m->la=xm(m->nodes*sizeof(double));m->mu=xm(m->nodes*sizeof(double));
 double lx=o->nx*o->h,ly=o->ny*o->h,x0=.5*(17000-lx),pi=3.141592653589793,vmin=1e30,vmax=-1e30;
 for(int k=0;k<m->nz;k++)for(int j=0;j<m->ny;j++)for(int i=0;i<m->nx;i++){
  double x=i*o->h,y=j*o->h,z=k*o->h,yc=y-.5*ly;
  /* Paper-guided 45-degree crossline extension plus a smooth y-z warp. */
  double warp=yc+150*sin(2*pi*y/(ly>0?ly:o->h))*sin(pi*z/(o->nz*o->h));
  double vp=interp(raw,(x0+x+warp)/25,z/25),vs=.67*vp,rho=2000;size_t q=node(i,j,k,m->nx,m->ny);
  m->vp[q]=vp;m->vs[q]=vs;m->rho[q]=rho;m->mu[q]=rho*vs*vs;m->la[q]=rho*(vp*vp-2*vs*vs);
  if(vp<vmin)vmin=vp;if(vp>vmax)vmax=vp;
 }
 free(raw);printf("[FREQ3D] model_nodes=%zu vp_min=%.6g vp_max=%.6g vs_over_vp=0.67 rho=2000\n",m->nodes,vmin,vmax);return 0;
}
static void model_free(Model*m){free(m->vp);free(m->vs);free(m->rho);free(m->la);free(m->mu);memset(m,0,sizeof(*m));}
static int bdist(int i,int j,int k,int nx,int ny,int nz){int d=i;if(nx-i<d)d=nx-i;if(j<d)d=j;if(ny-j<d)d=ny-j;if(k<d)d=k;if(nz-k<d)d=nz-k;return d;}
static void shift(int a,int s,int*di,int*dj,int*dk){*di=*dj=*dk=0;if(a==0)*di=s;else if(a==1)*dj=s;else*dk=s;}

static int assemble(const Opt*o,const Model*m,Coo*a){
 size_t nd=3*m->nodes;if(nd>INT32_MAX)return -1;int n=(int)nd;a->cap=15*nd;
 a->r=xm(a->cap*sizeof(int));a->c=xm(a->cap*sizeof(int));a->ar=xm(a->cap*sizeof(double));a->ai=xm(a->cap*sizeof(double));
 double w=2*3.141592653589793*o->freq,w2=w*w,h2=o->h*o->h;
 for(int k=0;k<m->nz;k++)for(int j=0;j<m->ny;j++)for(int i=0;i<m->nx;i++){
  size_t q=node(i,j,k,m->nx,m->ny);int d=bdist(i,j,k,o->nx,o->ny,o->nz);
  for(int comp=0;comp<3;comp++){int row=dof(q,comp);if(d==0){add(a,row,row,1,0);continue;}
   double scale=m->rho[q]*w2,diag=-1,sponge=0;if(d<=4){double t=(5-d)/4.0;sponge=1.5*t*t;}double diagi=-(.01+sponge);
   for(int ax=0;ax<3;ax++){double cq=comp==ax?m->la[q]+2*m->mu[q]:m->mu[q];
    for(int s=-1;s<=1;s+=2){int di,dj,dk;shift(ax,s,&di,&dj,&dk);size_t qn=node(i+di,j+dj,k+dk,m->nx,m->ny);
     double cn=comp==ax?m->la[qn]+2*m->mu[qn]:m->mu[qn],c=.5*(cq+cn)/(h2*scale);add(a,row,dof(qn,comp),-c,0);diag+=c;}}
   for(int other=0;other<3;other++)if(other!=comp){double c=(m->la[q]+m->mu[q])/(4*h2*scale);
    for(int s1=-1;s1<=1;s1+=2)for(int s2=-1;s2<=1;s2+=2){int a1,b1,c1,a2,b2,c2;shift(comp,s1,&a1,&b1,&c1);shift(other,s2,&a2,&b2,&c2);
     size_t qn=node(i+a1+a2,j+b1+b2,k+c1+c2,m->nx,m->ny);add(a,row,dof(qn,other),-(double)(s1*s2)*c,0);}}
   add(a,row,row,diag,diagi);
  }
 }
 printf("[FREQ3D] matrix_order=%d matrix_nnz=%zu frequency_hz=%.6g\n",n,a->nnz,o->freq);
 return n;
}

/* Reuse the released C model generator and assembly, with a supplied Vp volume. */
static Opt options(double frequency) {
    Opt value = {FWI_EX, FWI_EY, FWI_EZ, FWI_SPACING, frequency, NULL};
    return value;
}
int load_truth(const char *path, double *vp) {
    Opt opt = options(5);
    opt.model = path;
    Model model = {0};
    if (model_load(&opt, &model)) return 1;
    memcpy(vp, model.vp, model.nodes * sizeof(double));
    model_free(&model);
    return 0;
}
int build_matrix(double frequency, const double *vp, Matrix *out) {
    Opt opt = options(frequency);
    Model model = {0};
    model.nx = opt.nx + 1;
    model.ny = opt.ny + 1;
    model.nz = opt.nz + 1;
    model.nodes = (size_t)model.nx * model.ny * model.nz;
    model.rho = xm(model.nodes * sizeof(double));
    model.la = xm(model.nodes * sizeof(double));
    model.mu = xm(model.nodes * sizeof(double));
    for (size_t i = 0; i < model.nodes; ++i) {
        const double p = vp[i], s = .67*p, rho = 2000;
        model.rho[i] = rho;
        model.mu[i] = rho*s*s;
        model.la[i] = rho*(p*p - 2*s*s);
    }
    Coo matrix = {0};
    int n = assemble(&opt, &model, &matrix);
    model_free(&model);
    *out = (Matrix){matrix.r, matrix.c, matrix.ar, matrix.ai, (int64_t)matrix.nnz};
    return n;
}
void release_matrix(Matrix *matrix) {
    free(matrix->r); free(matrix->c); free(matrix->ar); free(matrix->ai);
    *matrix = (Matrix){0};
}
