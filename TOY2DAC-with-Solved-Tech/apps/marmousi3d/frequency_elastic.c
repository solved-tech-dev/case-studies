// SPDX-License-Identifier: AGPL-3.0-only
// Copyright 2026 Solved Tech
// Standalone Marmousi-derived 3D elastic forward model.
#define _POSIX_C_SOURCE 200809L
#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "pangulu_adapter.h"
typedef struct { int nx,ny,nz; double h,freq; const char *model,*out; } Opt;
typedef struct { int *r,*c; double *ar,*ai; size_t nnz,cap; } Coo;
typedef struct { int nx,ny,nz; size_t nodes; double *vp,*vs,*rho,*la,*mu; } Model;

static double now(void){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+1e-9*t.tv_nsec;}
static void *xm(size_t n){void*p=malloc(n);if(!p){fprintf(stderr,"alloc %zu failed\n",n);exit(2);}return p;}
static void *xc(size_t n,size_t s){void*p=calloc(n,s);if(!p){fprintf(stderr,"alloc failed\n");exit(2);}return p;}
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

static int assemble(const Opt*o,const Model*m,Coo*a,double**bro,double**bio){
 size_t nd=3*m->nodes;if(nd>INT32_MAX)return -1;int n=(int)nd;a->cap=15*nd;
 a->r=xm(a->cap*sizeof(int));a->c=xm(a->cap*sizeof(int));a->ar=xm(a->cap*sizeof(double));a->ai=xm(a->cap*sizeof(double));
 double*br=xc(nd,sizeof(double)),*bi=xc(nd,sizeof(double)),w=2*3.141592653589793*o->freq,w2=w*w,h2=o->h*o->h;
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
 int si=o->nx/2,sj=o->ny/2,sk=o->nz>2?2:1;size_t q=node(si,sj,sk,m->nx,m->ny);br[dof(q,2)]=1;
 printf("[FREQ3D] source_node=(%d,%d,%d) source_dof=%d component=z\n",si,sj,sk,dof(q,2));
 printf("[FREQ3D] matrix_order=%d matrix_nnz=%zu frequency_hz=%.6g\n",n,a->nnz,o->freq);*bro=br;*bio=bi;return n;
}
static void coo_free(Coo*a){free(a->r);free(a->c);free(a->ar);free(a->ai);memset(a,0,sizeof(*a));}
static double residual(int n,const Coo*a,const double*br,const double*bi,const double*xr,const double*xi){
 double*yr=xc(n,sizeof(double)),*yi=xc(n,sizeof(double));for(size_t p=0;p<a->nnz;p++){int r=a->r[p],c=a->c[p];yr[r]+=a->ar[p]*xr[c]-a->ai[p]*xi[c];yi[r]+=a->ar[p]*xi[c]+a->ai[p]*xr[c];}
 long double rn=0,bn=0;for(int i=0;i<n;i++){long double x=yr[i]-br[i],y=yi[i]-bi[i];rn+=x*x+y*y;bn+=(long double)br[i]*br[i]+(long double)bi[i]*bi[i];}free(yr);free(yi);return sqrt((double)(rn/(bn>0?bn:1)));
}

static int output(const Opt*o,const Model*m,const double*xr,const double*xi,double rr,double rr32){
 char p[4096];FILE*f;snprintf(p,sizeof(p),"%s/model.vtk",o->out);f=fopen(p,"w");if(!f)return -1;
 fprintf(f,"# vtk DataFile Version 3.0\nMarmousi-derived 3D model\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %d %d %d\nORIGIN 0 0 0\nSPACING %.17g %.17g %.17g\nPOINT_DATA %zu\nSCALARS Vp_m_per_s double 1\nLOOKUP_TABLE default\n",m->nx,m->ny,m->nz,o->h,o->h,o->h,m->nodes);
 for(size_t q=0;q<m->nodes;q++)fprintf(f,"%.17g\n",m->vp[q]);fprintf(f,"SCALARS Vs_m_per_s double 1\nLOOKUP_TABLE default\n");for(size_t q=0;q<m->nodes;q++)fprintf(f,"%.17g\n",m->vs[q]);fclose(f);
 snprintf(p,sizeof(p),"%s/field.vtk",o->out);f=fopen(p,"w");if(!f)return -1;
 fprintf(f,"# vtk DataFile Version 3.0\n3D elastic frequency field\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %d %d %d\nORIGIN 0 0 0\nSPACING %.17g %.17g %.17g\nPOINT_DATA %zu\nVECTORS displacement_real double\n",m->nx,m->ny,m->nz,o->h,o->h,o->h,m->nodes);
 for(size_t q=0;q<m->nodes;q++)fprintf(f,"%.17g %.17g %.17g\n",xr[3*q],xr[3*q+1],xr[3*q+2]);fprintf(f,"VECTORS displacement_imag double\n");for(size_t q=0;q<m->nodes;q++)fprintf(f,"%.17g %.17g %.17g\n",xi[3*q],xi[3*q+1],xi[3*q+2]);
 fprintf(f,"SCALARS displacement_amplitude double 1\nLOOKUP_TABLE default\n");for(size_t q=0;q<m->nodes;q++){double s=0;for(int c=0;c<3;c++)s+=xr[3*q+c]*xr[3*q+c]+xi[3*q+c]*xi[3*q+c];fprintf(f,"%.17g\n",sqrt(s));}fclose(f);
 snprintf(p,sizeof(p),"%s/field_complex.bin",o->out);f=fopen(p,"wb");if(!f)return -1;int32_t d[3]={m->nx,m->ny,m->nz};fwrite(d,4,3,f);for(size_t q=0;q<m->nodes;q++){double v[6]={xr[3*q],xi[3*q],xr[3*q+1],xi[3*q+1],xr[3*q+2],xi[3*q+2]};fwrite(v,8,6,f);}fclose(f);
 snprintf(p,sizeof(p),"%s/midplane_y.csv",o->out);f=fopen(p,"w");if(!f)return -1;fprintf(f,"x_m,z_m,vp,ux_real,ux_imag,uy_real,uy_imag,uz_real,uz_imag,amplitude\n");int j=o->ny/2;
 for(int k=0;k<m->nz;k++)for(int i=0;i<m->nx;i++){size_t q=node(i,j,k,m->nx,m->ny);double s=0;for(int c=0;c<3;c++)s+=xr[3*q+c]*xr[3*q+c]+xi[3*q+c]*xi[3*q+c];fprintf(f,"%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",i*o->h,k*o->h,m->vp[q],xr[3*q],xi[3*q],xr[3*q+1],xi[3*q+1],xr[3*q+2],xi[3*q+2],sqrt(s));}fclose(f);
 snprintf(p,sizeof(p),"%s/metadata.txt",o->out);f=fopen(p,"w");if(!f)return -1;fprintf(f,"element_grid=%dx%dx%d\nnode_grid=%dx%dx%d\nspacing_m=%.17g\nfrequency_hz=%.17g\nunknowns=%zu\nrelative_residual_against_fp64_system=%.17g\nrelative_residual_against_fp32_system=%.17g\n",o->nx,o->ny,o->nz,m->nx,m->ny,m->nz,o->h,o->freq,3*m->nodes,rr,rr32);fclose(f);return 0;
}
static int opts(int ac,char**av,Opt*o){o->nx=112;o->ny=16;o->nz=40;o->h=75;o->freq=5;o->model=o->out=NULL;for(int i=1;i<ac;i++){if(i+1>=ac)return -1;if(!strcmp(av[i],"--nx"))o->nx=atoi(av[++i]);else if(!strcmp(av[i],"--ny"))o->ny=atoi(av[++i]);else if(!strcmp(av[i],"--nz"))o->nz=atoi(av[++i]);else if(!strcmp(av[i],"--spacing"))o->h=atof(av[++i]);else if(!strcmp(av[i],"--frequency"))o->freq=atof(av[++i]);else if(!strcmp(av[i],"--marmousi"))o->model=av[++i];else if(!strcmp(av[i],"--output"))o->out=av[++i];else return -1;}return(!o->model||!o->out)?-1:0;}


int main(int argc, char **argv) {
    setvbuf(stdout, NULL, _IOLBF, 0);
    Opt options;
    if (opts(argc, argv, &options) || options.nx < 3 || options.ny < 3 || options.nz < 3 ||
        !(options.h > 0) || !(options.freq > 0)) {
        fprintf(stderr, "Invalid Marmousi-3D arguments\n");
        return 2;
    }
    Model model = {0};
    Coo matrix = {0};
    double *br = NULL, *bi = NULL, *xr = NULL, *xi = NULL;
    if (model_load(&options, &model)) return 2;
    int n = assemble(&options, &model, &matrix, &br, &bi);
    int status = n < 0 ? 2 : case_factor(n, (int64_t)matrix.nnz, matrix.r, matrix.c, matrix.ar, matrix.ai);
    if (!status) {
        xr = xm((size_t)n * sizeof(*xr));
        xi = xm((size_t)n * sizeof(*xi));
        status = case_solve(n, br, bi, xr, xi);
    }
    if (!status) {
        double rr = residual(n, &matrix, br, bi, xr, xi);
        printf("MARMOUSI3D_RESIDUAL n=%d nnz=%zu relative_l2=%.17g\n", n, matrix.nnz, rr);
        if (!isfinite(rr) || rr > 1e-8) status = 71;
        else if (output(&options, &model, xr, xi, rr, NAN)) status = 74;
    }
    if (case_cleanup() && !status) status = 70;
    free(xr);
    free(xi);
    free(br);
    free(bi);
    coo_free(&matrix);
    model_free(&model);
    return status;
}
