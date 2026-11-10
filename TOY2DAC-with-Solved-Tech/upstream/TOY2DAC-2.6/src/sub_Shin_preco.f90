!Copyright 2008-2011 SEISCOPE project, All rights reserved.
!Copyright 2013-2021 SEISCOPEII project, All rights reserved.
!****************************************************************************!
! Computation of the Shin pseudo-Hessian: the diagonal D of the reduced      !
! Hessian is approximated using virtual sources.                             !
!****************************************************************************!
subroutine Shin_preconditioning(nn,pbdir,acqui,inv,inc_wfld,v,v_preco)

  implicit none

  include "mpif.h"
#include "common.h"
#include "pbdirect.h"
  include "acqui.h"
  include "inversion.h"

  !PBDIRECT
  TYPE (pbdirect) :: pbdir
  !ACQUI
  TYPE(acquisition) :: acqui
  !INVERSION
  TYPE(inversion) :: inv

  integer :: nn,iw,is,ipar,i,i1,i2,j
  real,dimension(nn) :: v,v_preco
  complex,dimension(pbdir%n1e,pbdir%n2e,pbdir%nb_wfld,acqui%nsrc,pbdir%nw) :: inc_wfld  
  real,dimension(:,:,:),allocatable :: preco
  real :: threshold,max,norm_grad,norm_grad_preco
  complex :: xi,omegac
  complex,dimension(:),allocatable :: fld

  !COMPUTE PRECONDITIONER (MULTI-PARAMETER)
  xi=CMPLX(0.,1.)
  allocate(fld(5))
  allocate(preco(pbdir%n1,pbdir%n2,inv%npar))  
  preco(:,:,:)=0.
  do ipar=1,inv%npar
     if(inv%invpar(ipar).eq.1) then !VP
        do iw=1,pbdir%nw
           omegac=CMPLX(pbdir%freqw(iw),0)
           do is=1,acqui%nsrc
              do i=1,pbdir%nb_wfld
                 do i2=1,pbdir%n2
                    do i1=1,pbdir%n1
                       preco(i1,i2,ipar)=&
                            preco(i1,i2,ipar)+&
                            4.*REAL(omegac**4)/&                           
                            (pbdir%rho(i1,i2)**2)/&
                            (pbdir%vp(i1,i2)**6)/&
                            (abs(1.-0.5*xi/pbdir%qp(i1,i2))**4)*&
                            (abs(inc_wfld(i1+pbdir%npml,i2+pbdir%npml,i,is,iw))**4)
!WRITE(*,*)preco(i1,i2,ipar),REAL(omegac**4),(pbdir%rho(i1,i2)**2), (pbdir%vp(i1,i2)**6),abs(inc_wfld(i1+pbdir%npml,i2+pbdir%npml,i,is,iw))**4
                    enddo
                 enddo
              enddo
           enddo
        enddo
     elseif(inv%invpar(ipar).eq.2) then !RHO
        do iw=1,pbdir%nw
           omegac=CMPLX(pbdir%freqw(iw),0)
           do is=1,acqui%nsrc
              do i=1,pbdir%nb_wfld
                 do i2=1,pbdir%n2
                    do i1=1,pbdir%n1
                       call sub_drho5_fld_preco(pbdir%nb_wfld,acqui%nsrc,is,i1,i2,i,&
                            omegac,pbdir,inc_wfld,fld)
                       do j=1,5
                          preco(i1,i2,ipar)=&
                               preco(i1,i2,ipar)+fld(j)                               
                       enddo
                    enddo
                 enddo
              enddo
           enddo
        enddo
     elseif(inv%invpar(ipar).eq.3) then !QP
        do iw=1,pbdir%nw
           omegac=CMPLX(pbdir%freqw(iw),0)
           do is=1,acqui%nsrc
              do i=1,pbdir%nb_wfld
                 do i2=1,pbdir%n2
                    do i1=1,pbdir%n1
                       preco(i1,i2,ipar)=&
                            preco(i1,i2,ipar)+&
                            REAL(omegac**4)/&
                            (pbdir%rho(i1,i2)**2)/&
                            (pbdir%vp(i1,i2)**4)/&
                            (abs(1.-0.5*xi/pbdir%qp(i1,i2))**6)/&
                            (pbdir%qp(i1,i2)**4)*&
                            (abs(inc_wfld(i1+pbdir%npml,i2+pbdir%npml,i,is,iw))**4)
                    enddo
                 enddo
              enddo
           enddo
        enddo
     endif
  enddo

  !APPLY THRESHOLD
  threshold=inv%coeff_damping_preco
  max=maxval(preco(:,:,:))
  write(15,*) 'IN PCGNL max vaut : ',max
  write(15,*) 'IN PCGNL threshold vaut : ',threshold
  do ipar=1,inv%npar
     do i2=1,pbdir%n2
        do i1=1,pbdir%n1
           preco(i1,i2,ipar)=1./(preco(i1,i2,ipar)+threshold*max)
        enddo
     enddo
  enddo

  !APPLY PRECONDITIONER
  do ipar=1,inv%npar
     do i2=1,pbdir%n2
        do i1=1,pbdir%n1
           v_preco(i1+(i2-1)*pbdir%n1+(ipar-1)*pbdir%n1*pbdir%n2)=&
                preco(i1,i2,ipar)*v(i1+(i2-1)*pbdir%n1+(ipar-1)*pbdir%n1*pbdir%n2)
        enddo
     enddo
  enddo

  !NORMALIZE PRECONDITIONED GRADIENT
  v_preco(:)=v_preco(:)*inv%scal_preco

end subroutine Shin_preconditioning

!****************************************************************************!
! Computation of the Shin pseudo-Hessian: the diagonal D of the reduced      !
! Hessian is approximated using virtual sources.                             !
! Adaptation to the frozen coefficient case                                  !
!****************************************************************************!
subroutine Shin_preconditioning_frozen(nn,pbdir,acqui,inv,inc_wfld,v,v_preco)

  implicit none

  include "mpif.h"
#include "common.h"
#include "pbdirect.h"
  include "acqui.h"
  include "inversion.h"

  !PBDIRECT
  TYPE (pbdirect) :: pbdir
  !ACQUI
  TYPE(acquisition) :: acqui
  !INVERSION
  TYPE(inversion) :: inv

  integer :: nn,iw,is,i,i1,i2,k,j,ipar
  real,dimension(nn) :: v,v_preco
  complex,dimension(pbdir%n1e,pbdir%n2e,pbdir%nb_wfld,acqui%nsrc,pbdir%nw) :: inc_wfld
  real,dimension(:,:,:),allocatable :: preco
  real :: threshold,max,norm_grad,norm_grad_preco
  complex,dimension(:),allocatable :: fld
  complex :: xi,omegac

  !COMPUTE PRECONDITIONER (MULTI-PARAMETER)
  xi=CMPLX(0.,1.)
  allocate(fld(5))
  allocate(preco(pbdir%n1,pbdir%n2,inv%npar))  
  preco(:,:,:)=0.
  do ipar=1,inv%npar
     if(inv%invpar(ipar).eq.1) then !VP
        do iw=1,pbdir%nw
           omegac=CMPLX(pbdir%freqw(iw),0)
           do is=1,acqui%nsrc
              do i=1,pbdir%nb_wfld
                 do i2=1,pbdir%n2
                    do i1=1,pbdir%n1
                       preco(i1,i2,ipar)=&
                            preco(i1,i2,ipar)+&
                            4.*REAL(omegac**4)/&                           
                            (pbdir%rho(i1,i2)**2)/&
                            (pbdir%vp(i1,i2)**6)/&
                            (abs(1.-0.5*xi/pbdir%qp(i1,i2))**4)*&
                            (abs(inc_wfld(i1+pbdir%npml,i2+pbdir%npml,i,is,iw))**4)
                    enddo
                 enddo
              enddo
           enddo
        enddo
     elseif(inv%invpar(ipar).eq.2) then !RHO
        do iw=1,pbdir%nw
           omegac=CMPLX(pbdir%freqw(iw),0)
           do is=1,acqui%nsrc
              do i=1,pbdir%nb_wfld
                 do i2=1,pbdir%n2
                    do i1=1,pbdir%n1
                       call sub_drho5_fld_preco(pbdir%nb_wfld,acqui%nsrc,is,i1,i2,i,&
                            omegac,pbdir,inc_wfld,fld)
                       do j=1,5
                          preco(i1,i2,ipar)=&
                               preco(i1,i2,ipar)+fld(j)                               
                       enddo
                    enddo
                 enddo
              enddo
           enddo
        enddo
     elseif(inv%invpar(ipar).eq.3) then !QP
        do iw=1,pbdir%nw
           omegac=CMPLX(pbdir%freqw(iw),0)
           do is=1,acqui%nsrc
              do i=1,pbdir%nb_wfld
                 do i2=1,pbdir%n2
                    do i1=1,pbdir%n1
                       preco(i1,i2,ipar)=&
                            preco(i1,i2,ipar)+&
                            REAL(omegac**4)/&
                            (pbdir%rho(i1,i2)**2)/&
                            (pbdir%vp(i1,i2)**4)/&
                            (abs(1.-0.5*xi/pbdir%qp(i1,i2))**6)/&
                            (pbdir%qp(i1,i2)**4)*&
                            (abs(inc_wfld(i1+pbdir%npml,i2+pbdir%npml,i,is,iw))**4)
                    enddo
                 enddo
              enddo
           enddo
        enddo
     endif
  enddo

  !APPLY THRESHOLD
  threshold=inv%coeff_damping_preco
  max=maxval(preco(:,:,:))
 ! write(15,*) 'IN PCGNL max vaut : ',max
 ! write(15,*) 'IN PCGNL threshold vaut : ',threshold
  do ipar=1,inv%npar
     do i2=1,pbdir%n2
        do i1=1,pbdir%n1
           preco(i1,i2,ipar)=1./(preco(i1,i2,ipar)+threshold*max)
        enddo
     enddo
  enddo
  !APPLY PRECONDITIONER
  k=1
  do ipar=1,inv%npar
     do i2=1,pbdir%n2
        do i1=inv%ibathy(i2)+inv%ideadzone,pbdir%n1        
           v_preco(k)=preco(i1,i2,ipar)*v(k)
           k=k+1        
        enddo
     enddo
  enddo
  v_preco(:)=v_preco(:)*inv%scal_preco

end subroutine Shin_preconditioning_frozen

