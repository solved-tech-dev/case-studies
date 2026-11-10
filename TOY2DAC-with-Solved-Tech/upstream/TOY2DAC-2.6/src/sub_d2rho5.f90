!Copyright 2008-2011 SEISCOPE project, All rights reserved.
!Copyright 2013-2021 SEISCOPEII project, All rights reserved.
!Computation of the second-order derivative term with respect to the density assuming
!a standard 5 pts stencil on a staggered-grid 
subroutine sub_d2rho5_fld(nb_wfld,nsrc,isrc,i1,i2,ifld,omegac,pbdir,inc_wfld,fld)
  
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
  c1=(2.*omegac**2)/((pbdir%vp(i1,i2)**2)*((1.-0.5*xi/pbdir%qp(i1,i2))**2)&
       *(pbdir%rho(i1,i2)**3))&
       -4./((pbdir%h**2)*pbdir%rho(i1,i2)**3)
  c2=1./((pbdir%h**2)*pbdir%rho(i1,i2)**3)
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
  
end subroutine sub_d2rho5_fld



