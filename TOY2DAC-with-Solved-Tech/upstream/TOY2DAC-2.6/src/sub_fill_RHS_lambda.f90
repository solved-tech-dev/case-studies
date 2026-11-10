!Copyright 2008-2011 SEISCOPE project, All rights reserved.
!Copyright 2013-2021 SEISCOPEII project, All rights reserved.
!! SUBROUTINE TO FILL THE DENSE RHS OF MUMPS FOR COMPUTATION OF RHS FOR HESSIAN-VECTOR 
!! PRODUCT

subroutine sub_fill_RHS_lambda(pbdir,acqui,inv,adj_wfld,v,omega)

  implicit none
#include "common.h"
#include "pbdirect.h"
  include "acqui.h"
  include "inversion.h"

  !PBDIRECT
  TYPE (pbdirect) :: pbdir
  !ACQUI
  TYPE(acquisition) :: acqui
  !INVERSION
  TYPE (inversion) :: inv  
  !IN
  real :: omega
  complex,dimension(pbdir%n1e,pbdir%n2e,pbdir%nb_wfld,acqui%nsrc) :: adj_wfld
  !IN/OUT 
  real,dimension(pbdir%n1*pbdir%n2*inv%npar) :: v
  !Local variables
  integer :: i,i1,i2,isrc,ipar,k
  integer :: n1,n2,n1e,n2e,nsrc,npml,npar
  real :: h
  complex :: xi
  complex,dimension(:),allocatable :: fld
  complex :: omegac

  !The i complex
  xi=CMPLX(0.,1.)

  omegac=CMPLX(omega,pbdir%slaplace)
  !Size parameters
  n1=pbdir%n1
  n2=pbdir%n2
  n1e=pbdir%n1e
  n2e=pbdir%n2e
  npml=pbdir%npml
  h=pbdir%h
  nsrc=acqui%nsrc
  npar=inv%npar

  !Allocation for local variables
  allocate(fld(5))

  !Restriction of J_s*v products (computation of R_s^{T}R_sJ_s*v)
  !This routine modify directly pbdir%id%RHS
  if(pbdir%hicks.eq.0) then
     call sub_restriction(pbdir,acqui)        
  elseif(pbdir%hicks.eq.1) then
     call sub_restriction_Hicks(pbdir,acqui)        
  endif

  !Fill in RHS   
  do isrc = 1, nsrc
     do i=1,pbdir%nb_wfld        
        do ipar=1,npar
           do i2=1,n2
              do i1=1,n1                         
                 k=(isrc-1)*(n1e*n2e*pbdir%nb_wfld)+&
                      (i-1)*n1e*n2e+(i2+npml-1)*n1e+i1+npml       
                 if(inv%invpar(ipar).eq.1) then !VP
                    pbdir%id%RHS(k)=pbdir%id%RHS(k)+&
                         v((ipar-1)*n1*n2+(i2-1)*n1+i1)*&                
                         CONJG(-2.*omegac**2/(pbdir%rho(i1,i2)*&
                         pbdir%vp(i1,i2)**3*(1.-0.5*xi/pbdir%qp(i1,i2))**2))*&
                         CONJG(adj_wfld(i1+npml,i2+npml,i,isrc))
                 elseif (inv%invpar(ipar).eq.2) then !RHO
                    call sub_drho5_fld(pbdir%nb_wfld,acqui%nsrc,isrc,i1,i2,i,&
                         omegac,pbdir,adj_wfld,fld)
                    pbdir%id%RHS(k)=pbdir%id%RHS(k)+&
                         v((ipar-1)*n1*n2+(i2-1)*n1+i1)*CONJG(fld(1))
                    pbdir%id%RHS(k-1)=pbdir%id%RHS(k-1)-&
                         v((ipar-1)*n1*n2+(i2-1)*n1+i1)*CONJG(fld(2))
                    pbdir%id%RHS(k+1)=pbdir%id%RHS(k+1)-&
                         v((ipar-1)*n1*n2+(i2-1)*n1+i1)*CONJG(fld(3))
                    pbdir%id%RHS(k-n1e)=pbdir%id%RHS(k-n1e)-&
                         v((ipar-1)*n1*n2+(i2-1)*n1+i1)*CONJG(fld(4))
                    pbdir%id%RHS(k+n1e)=pbdir%id%RHS(k+n1e)-&
                         v((ipar-1)*n1*n2+(i2-1)*n1+i1)*CONJG(fld(5))
                 elseif (inv%invpar(ipar).eq.3) then !QP
                    pbdir%id%RHS(k)=pbdir%id%RHS(k)+&
                         v((ipar-1)*n1*n2+(i2-1)*n1+i1)*&                
                         CONJG(-xi*omegac**2/(pbdir%rho(i1,i2)*&
                         pbdir%vp(i1,i2)**2*(1.-0.5*xi/pbdir%qp(i1,i2))**3&
                         *pbdir%qp(i1,i2)**2))*&
                         CONJG(adj_wfld(i1+npml,i2+npml,i,isrc))
                 endif
              enddo
           enddo
        enddo
     enddo
  enddo

  !Fill in RHS   
!!$  do isrc = 1, nsrc
!!$     do i=1,pbdir%nb_wfld        
!!$        do i2=1,n2
!!$           do i1=1,n1                         
!!$              k=(isrc-1)*(n1e*n2e*pbdir%nb_wfld)+&
!!$                   (i-1)*n1e*n2e+(i2+npml-1)*n1e+i1+npml       
!!$              pbdir%id%RHS(k)=pbdir%id%RHS(k)+v((i2-1)*n1+i1)&                
!!$                   *CONJG(-2.*omegac**2/(pbdir%rho(i1,i2)*&
!!$                   pbdir%vp(i1,i2)**3*(1.-0.5*xi/pbdir%qp(i1,i2))**2))&
!!$                   *CONJG(adj_wfld(i1+npml,i2+npml,i,isrc))                           
!!$           enddo
!!$        enddo
!!$     enddo
!!$  enddo

  do k=1,pbdir%n1e*pbdir%n2e*acqui%nsrc*pbdir%nb_wfld
     pbdir%id%RHS(k)=-CONJG(pbdir%id%RHS(k))/(pbdir%h**2)
  enddo

end subroutine sub_fill_RHS_lambda

