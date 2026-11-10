!Copyright 2008-2011 SEISCOPE project, All rights reserved.
!Copyright 2013-2021 SEISCOPEII project, All rights reserved.

!Computation of the gradient term with respect to the density assuming
!a standard 5 pts stencil on a staggered-grid 

subroutine sub_drho5(i1,i2,nb_wfld,&
     nsrc,omegac,pbdir,inv,inc_wfld,adj_wfld,grad_term)
  
  implicit none
#include "common.h"
#include "pbdirect.h"
  include "inversion.h"
    
  !PBDIRECT
  TYPE (pbdirect) :: pbdir
  !INVERSION
  TYPE (inversion) :: inv
  !IN
  integer :: i1,i2,nb_wfld,nsrc
  COMPLEX :: omegac
  complex,dimension(pbdir%n1e,pbdir%n2e,nb_wfld,nsrc) :: inc_wfld,adj_wfld
  !IN/OUT
  complex :: grad_term
  !Local variables
  integer :: isrc
  complex :: xi,c1,c2,c3,u1,u2,u3,u4,u5  
  parameter(xi=CMPLX(0.,1.))
  
  !COEFFICIENTS COMPUTATION
  c1=-(omegac**2)/((pbdir%vp(i1,i2)**2)*((1.-0.5*xi/pbdir%qp(i1,i2))**2)&
       *(pbdir%rho(i1,i2)**2))&
       +2./((pbdir%h**2)*pbdir%rho(i1,i2)**2)
  c2=-0.5/((pbdir%h**2)*pbdir%rho(i1,i2)**2)
  c3=-c2
  
  !COMPUTATION OF THE GRADIENT TERM
  grad_term=CMPLX(0.,0.)
  do isrc=1,nsrc
     u1=&
          c1*inc_wfld(pbdir%npml+i1,pbdir%npml+i2,1,isrc)+&
          c2*inc_wfld(pbdir%npml+i1-1,pbdir%npml+i2,1,isrc)+&
          c2*inc_wfld(pbdir%npml+i1+1,pbdir%npml+i2,1,isrc)+&
          c2*inc_wfld(pbdir%npml+i1,pbdir%npml+i2-1,1,isrc)+&
          c2*inc_wfld(pbdir%npml+i1,pbdir%npml+i2+1,1,isrc)
     grad_term=grad_term+u1*adj_wfld(pbdir%npml+i1,pbdir%npml+i2,1,isrc)
     u2=&          
          c3*inc_wfld(pbdir%npml+i1-1,pbdir%npml+i2,1,isrc)+&
          c2*inc_wfld(pbdir%npml+i1,pbdir%npml+i2,1,isrc)
     grad_term=grad_term+u2*adj_wfld(pbdir%npml+i1-1,pbdir%npml+i2,1,isrc)
     u3=&          
          c3*inc_wfld(pbdir%npml+i1+1,pbdir%npml+i2,1,isrc)+&
          c2*inc_wfld(pbdir%npml+i1,pbdir%npml+i2,1,isrc)
     grad_term=grad_term+u3*adj_wfld(pbdir%npml+i1+1,pbdir%npml+i2,1,isrc)
     u4=&          
          c3*inc_wfld(pbdir%npml+i1,pbdir%npml+i2-1,1,isrc)+&
          c2*inc_wfld(pbdir%npml+i1,pbdir%npml+i2,1,isrc)
     grad_term=grad_term+u4*adj_wfld(pbdir%npml+i1,pbdir%npml+i2-1,1,isrc)
     u5=&          
          c3*inc_wfld(pbdir%npml+i1,pbdir%npml+i2+1,1,isrc)+&
          c2*inc_wfld(pbdir%npml+i1,pbdir%npml+i2,1,isrc)
     grad_term=grad_term+u5*adj_wfld(pbdir%npml+i1,pbdir%npml+i2+1,1,isrc)
  enddo
  
end subroutine sub_drho5

subroutine sub_drho5_bis(i1,i2,ifld,nb_wfld,&
     nsrc,omegac,pbdir,inv,inc_wfld,adj_wfld,grad_term)
  
  implicit none
#include "common.h"
#include "pbdirect.h"
  include "inversion.h"
    
  !PBDIRECT
  TYPE (pbdirect) :: pbdir
  !INVERSION
  TYPE (inversion) :: inv
  
  !IN
  integer :: i1,i2,ifld,nb_wfld,nsrc
  complex :: omegac
  complex,dimension(pbdir%n1e,pbdir%n2e,nb_wfld,nsrc) :: inc_wfld,adj_wfld
  !IN/OUT
  complex :: grad_term
  complex,dimension(:),allocatable :: fld
  !Local variables
  integer :: isrc
  complex :: xi,c1,c2,c3
  complex,dimension(:,:),allocatable :: u
  parameter(xi=CMPLX(0.,1.))
  
  !COMPUTATION OF THE GRADIENT TERM
  grad_term=CMPLX(0.,0.)
  allocate(fld(5))
  do isrc=1,nsrc
     call sub_drho5_fld(nb_wfld,nsrc,isrc,i1,i2,ifld,omegac,pbdir,inc_wfld,fld)
     grad_term=grad_term+&
          fld(1)*adj_wfld(pbdir%npml+i1,pbdir%npml+i2,ifld,isrc)+&
          fld(2)*adj_wfld(pbdir%npml+i1-1,pbdir%npml+i2,ifld,isrc)+&
          fld(3)*adj_wfld(pbdir%npml+i1+1,pbdir%npml+i2,ifld,isrc)+&
          fld(4)*adj_wfld(pbdir%npml+i1,pbdir%npml+i2-1,ifld,isrc)+&
          fld(5)*adj_wfld(pbdir%npml+i1,pbdir%npml+i2+1,ifld,isrc)
  enddo
  deallocate(fld)
  
end subroutine sub_drho5_bis

subroutine sub_drho5_fld(nb_wfld,nsrc,isrc,i1,i2,ifld,omegac,pbdir,inc_wfld,fld)
  
  implicit none

#include "pbdirect.h"
  !PBDIRECT
  TYPE (pbdirect) :: pbdir
  !IN 
  integer :: nsrc,isrc,nb_wfld,i1,i2,ifld
  complex :: omegac
  complex,dimension(pbdir%n1e,pbdir%n2e,nb_wfld,nsrc) :: inc_wfld
  !IN/OUT
  complex,dimension(5) :: fld
  !LOCAL   
  complex :: xi,c1,c2,c3
  parameter(xi=CMPLX(0.,1.))

  !COEFFICIENTS COMPUTATION
  c1=-(omegac**2)/((pbdir%vp(i1,i2)**2)*((1.-0.5*xi/pbdir%qp(i1,i2))**2)&
       *(pbdir%rho(i1,i2)**2))&
       +2./((pbdir%h**2)*pbdir%rho(i1,i2)**2)
  c2=-0.5/((pbdir%h**2)*pbdir%rho(i1,i2)**2)
  c3=-c2
  
  fld(1)=&
       c1*inc_wfld(pbdir%npml+i1,pbdir%npml+i2,ifld,isrc)+&
       c2*inc_wfld(pbdir%npml+i1-1,pbdir%npml+i2,ifld,isrc)+&
       c2*inc_wfld(pbdir%npml+i1+1,pbdir%npml+i2,ifld,isrc)+&
       c2*inc_wfld(pbdir%npml+i1,pbdir%npml+i2-1,ifld,isrc)+&
       c2*inc_wfld(pbdir%npml+i1,pbdir%npml+i2+1,ifld,isrc)  
  fld(2)=&
       c3*inc_wfld(pbdir%npml+i1-1,pbdir%npml+i2,ifld,isrc)+&
       c2*inc_wfld(pbdir%npml+i1,pbdir%npml+i2,ifld,isrc)
  fld(3)=&       
       c3*inc_wfld(pbdir%npml+i1+1,pbdir%npml+i2,ifld,isrc)+&
       c2*inc_wfld(pbdir%npml+i1,pbdir%npml+i2,ifld,isrc)
  fld(4)=&
       c3*inc_wfld(pbdir%npml+i1,pbdir%npml+i2-1,ifld,isrc)+&
       c2*inc_wfld(pbdir%npml+i1,pbdir%npml+i2,ifld,isrc)
  fld(5)=&
       c3*inc_wfld(pbdir%npml+i1,pbdir%npml+i2+1,ifld,isrc)+&
       c2*inc_wfld(pbdir%npml+i1,pbdir%npml+i2,ifld,isrc)
  
end subroutine sub_drho5_fld

subroutine sub_drho5_fld_preco(nb_wfld,nsrc,isrc,i1,i2,ifld,omegac,pbdir,inc_wfld,fld)
  
  implicit none

#include "pbdirect.h"
  !PBDIRECT
  TYPE (pbdirect) :: pbdir
  !IN 
  integer :: nsrc,isrc,nb_wfld,i1,i2,ifld
  complex :: omegac
  complex,dimension(pbdir%n1e,pbdir%n2e,nb_wfld,nsrc) :: inc_wfld
  !IN/OUT
  complex,dimension(5) :: fld
  !LOCAL   
  complex :: xi,c1,c2,c3
  parameter(xi=CMPLX(0.,1.))

  !COEFFICIENTS COMPUTATION
  c1=-(omegac**2)/((pbdir%vp(i1,i2)**2)*((1.-0.5*xi/pbdir%qp(i1,i2))**2)&
       *(pbdir%rho(i1,i2)**2))&
       +2./((pbdir%h**2)*pbdir%rho(i1,i2)**2)
  c2=-0.5/((pbdir%h**2)*pbdir%rho(i1,i2)**2)
  c3=-c2
  
  fld(1)=&
       abs(c1)**2*abs(inc_wfld(pbdir%npml+i1,pbdir%npml+i2,ifld,isrc)**4)+&
       abs(c2)**2*abs(inc_wfld(pbdir%npml+i1-1,pbdir%npml+i2,ifld,isrc)**4)+&
       abs(c2)**2*abs(inc_wfld(pbdir%npml+i1+1,pbdir%npml+i2,ifld,isrc)**4)+&
       abs(c2)**2*abs(inc_wfld(pbdir%npml+i1,pbdir%npml+i2-1,ifld,isrc)**4)+&
       abs(c2)**2*abs(inc_wfld(pbdir%npml+i1,pbdir%npml+i2+1,ifld,isrc)**4)
  fld(2)=&
       abs(c3)**2*abs(inc_wfld(pbdir%npml+i1-1,pbdir%npml+i2,ifld,isrc)**4)+&
       abs(c2)**2*abs(inc_wfld(pbdir%npml+i1,pbdir%npml+i2,ifld,isrc)**4)
  fld(3)=&       
       abs(c3)**2*abs(inc_wfld(pbdir%npml+i1+1,pbdir%npml+i2,ifld,isrc)**4)+&
       abs(c2)**2*abs(inc_wfld(pbdir%npml+i1,pbdir%npml+i2,ifld,isrc)**4)
  fld(4)=&
       abs(c3)**2*abs(inc_wfld(pbdir%npml+i1,pbdir%npml+i2-1,ifld,isrc)**4)+&
       abs(c2)**2*abs(inc_wfld(pbdir%npml+i1,pbdir%npml+i2,ifld,isrc)**4)
  fld(5)=&
       abs(c3)**2*abs(inc_wfld(pbdir%npml+i1,pbdir%npml+i2+1,ifld,isrc)**4)+&
       abs(c2)**2*abs(inc_wfld(pbdir%npml+i1,pbdir%npml+i2,ifld,isrc)**4)
  
end subroutine sub_drho5_fld_preco



