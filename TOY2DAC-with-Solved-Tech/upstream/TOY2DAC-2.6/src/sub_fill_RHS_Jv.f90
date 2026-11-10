!Copyright 2008-2011 SEISCOPE project, All rights reserved.
!Copyright 2013-2021 SEISCOPEII project, All rights reserved.
!!$!SUBROUTINE TO FILL THE DENSE RHS OF MUMPS FOR COMPUTATION OF J_s*v
subroutine sub_fill_RHS_Jv(pbdir,acqui,inv,inc_wfld,v,omegac)
  
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
  
  integer :: npar,n1,n2,n1e,n2e,nsrc,npml
  real,dimension(pbdir%n1*pbdir%n2*inv%npar) :: v
  complex,dimension(pbdir%n1e,pbdir%n2e,pbdir%nb_wfld,acqui%nsrc) :: inc_wfld
  complex :: xi
  complex,dimension(:),allocatable :: fld
  COMPLEX :: omegac
  REAL :: h,xr,zr
  INTEGER :: isrc,irec,irec_true,ixre,izre,ixr,izr
  INTEGER :: i,ix,iz,i1,i2,k,ipar
  
  !The i complex
  xi=CMPLX(0.,1.)
  
  !Initialization of RHS to 0
  pbdir%id%RHS(:)=0.
  
  !Size parameters
  n1=pbdir%n1
  n2=pbdir%n2
  n1e=pbdir%n1e
  n2e=pbdir%n2e
  npml=pbdir%npml
  h=pbdir%h
  nsrc=acqui%nsrc
  npar=inv%npar
  
  allocate(fld(5))
  do isrc = 1, nsrc
     do i=1,pbdir%nb_wfld        
        do ipar=1,npar
           do i2=1,n2
              do i1=1,n1              
                 k=(isrc-1)*pbdir%nb_wfld*(n1e*n2e)+&
                      (i-1)*n1e*n2e+&                      
                      (i2+npml-1)*n1e+&
                      i1+npml
                 if(inv%invpar(ipar).eq.1) then !VP
                    pbdir%id%RHS(k)=pbdir%id%RHS(k)+&
                         v((ipar-1)*n1*n2+(i2-1)*n1+i1)&
                         *(2.*omegac**2/(pbdir%rho(i1,i2)*&
                         pbdir%vp(i1,i2)**3*(1.-0.5*xi/pbdir%qp(i1,i2))**2)*&
                         inc_wfld(i1+npml,i2+npml,i,isrc))/(pbdir%h**2)
                 elseif(inv%invpar(ipar).eq.2) then !RHO
                    call sub_drho5_fld(pbdir%nb_wfld,acqui%nsrc,isrc,i1,i2,i,&
                         omegac,pbdir,inc_wfld,fld)
                    pbdir%id%RHS(k)=pbdir%id%RHS(k)-&
                         v((ipar-1)*n1*n2+(i2-1)*n1+i1)*fld(1)/(pbdir%h**2)        
                    pbdir%id%RHS(k-1)=pbdir%id%RHS(k-1)-&
                         v((ipar-1)*n1*n2+(i2-1)*n1+i1)*fld(2)/(pbdir%h**2)        
                    pbdir%id%RHS(k+1)=pbdir%id%RHS(k+1)-&
                         v((ipar-1)*n1*n2+(i2-1)*n1+i1)*fld(3)/(pbdir%h**2)        
                    pbdir%id%RHS(k-n1e)=pbdir%id%RHS(k-n1e)-&
                         v((ipar-1)*n1*n2+(i2-1)*n1+i1)*fld(4)/(pbdir%h**2)        
                    pbdir%id%RHS(k+n1e)=pbdir%id%RHS(k+n1e)-&
                         v((ipar-1)*n1*n2+(i2-1)*n1+i1)*fld(5)/(pbdir%h**2)        
                 elseif(inv%invpar(ipar).eq.3) then !QP
                    pbdir%id%RHS(k)=pbdir%id%RHS(k)+&
                         v((ipar-1)*n1*n2+(i2-1)*n1+i1)&
                         *xi*omegac**2/(pbdir%rho(i1,i2)*&
                         pbdir%vp(i1,i2)**2*&
                         (1.-0.5*xi/pbdir%qp(i1,i2))**3*&
                         pbdir%qp(i1,i2)**2)*&
                         inc_wfld(i1+npml,i2+npml,i,isrc)/(pbdir%h**2)
                 endif
              enddo
           enddo
        enddo
     enddo
  enddo
  deallocate(fld)
  
  !!$  do isrc = 1, nsrc
!!$     do i=1,pbdir%nb_wfld
!!$        !Loop on discrete model space             
!!$        do i2=1,n2
!!$           do i1=1,n1   
!!$              k=(isrc-1)*pbdir%nb_wfld*(n1e*n2e)+&
!!$                   (i-1)*n1e*n2e+&
!!$                   (i2+npml-1)*n1e+&
!!$                   i1+npml
!!$              pbdir%id%RHS(k)=v((i2-1)*n1+i1)&
!!$                   *(2.*omegac**2/(pbdir%rho(i1,i2)*pbdir%vp(i1,i2)**3*(1.-0.5*xi/pbdir%qp(i1,i2))**2)*&
!!$                   inc_wfld(i1+npml,i2+npml,i,isrc))/(pbdir%h**2)
!!$           enddo
!!$        enddo
!!$     enddo
!!$  enddo
  
      
end subroutine sub_fill_RHS_Jv

