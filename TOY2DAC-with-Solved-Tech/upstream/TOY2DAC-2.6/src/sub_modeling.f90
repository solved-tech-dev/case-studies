!Copyright 2008-2011 SEISCOPE project, All rights reserved.
!Copyright 2013-2021 SEISCOPEII project, All rights reserved.
subroutine sub_modeling(pbdir,acqui)
  
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
  TYPE(INVERSION) :: inv
  !INTERN
  INTEGER :: iw
  
  IF (pbdir%tool == 1 .OR. pbdir%tool == 2) THEN
     ! ACOUSTIC (AN)ISOTROPIC MODELING IN THE FREQUENCY DOMAIN
     DO iw=1,pbdir%nw
        IF (pbdir%tool.eq.1)  CALL sub_pbdirect_FDFD_MATRIX(pbdir,pbdir%freqw(iw))
        IF (pbdir%tool.eq.2)  CALL sub_pbdirect_FDFD_MATRIXaniso(pbdir,pbdir%freqw(iw))
        pbdir%type_RHS=1
        IF (pbdir%tool.eq.1) CALL sub_pbdirect_FDFD_SOLVE_ACQUI(pbdir,acqui,inv)
        IF (pbdir%tool.eq.2) CALL sub_pbdirect_FDFD_SOLVE_ACQUIaniso(pbdir,acqui,inv)
        CALL sub_pbdirect_FDFD_WRITE_IO_DATA(pbdir,acqui,iw)
        CALL sub_pbdirect_FDFD_FINISH(pbdir)
     END DO
  END IF
  
end subroutine sub_modeling

!-------------------------------------------------------------------!
! SUB_MODELING_FWI_fcost:                                           !
! 1. Computation of the cost function                               !
!-------------------------------------------------------------------!
subroutine sub_modeling_FWI_fcost(pbdir,acqui,inv)
  
  implicit none
  include "mpif.h"
#include "common.h"
#include "pbdirect.h"
  include "acqui.h"
  include "inversion.h"
  
  !PBDIRECT
  TYPE (pbdirect) :: pbdir
  !INVERSION
  TYPE (inversion) :: inv
  !ACQUI
  TYPE(acquisition) :: acqui
  
  INTEGER :: i1,i2,isrc,k,iw,i
  real :: fcost_bathy
  REAL(KIND=8) ::time1,time2
  
  if(mype.eq.0) then     
     inv%gradient(:,:,:)=0.
     inv%fcost=0.
     inv%residual(:)=0.     
  end if
  
  if ( (pbdir%tool .eq. 1) .OR. (pbdir%tool .eq. 2)) THEN     
     ! ACOUSTIC (AN)ISOTROPIC MODELING IN THE FREQUENCY DOMAIN
     do iw=1,pbdir%nw !Loop over the frequencies
        time1= MPI_WTime()
        !STEP 1 : LU FACTORIZATION OF S(p)
        if (pbdir%tool == 1) then !Isotropic
           CALL sub_pbdirect_FDFD_MATRIX(pbdir,pbdir%freqw(iw))           
        elseif (pbdir%tool == 2) then !Anisotropic
           CALL sub_pbdirect_FDFD_MATRIXaniso(pbdir,pbdir%freqw(iw))
        endif
        time2= MPI_WTime()
        if(mype==0) then
           WRITE(*,*)"TIME FOR MATRIX BUILDING+FACTO",time2-time1,'sec'
        endif
        !STEP 2 : COMPUTATION OF THE INCIDENT WAVEFIELD
        time1= MPI_WTime()
        pbdir%type_RHS=1        
        if (pbdir%tool == 1) then
           CALL sub_pbdirect_FDFD_SOLVE_ACQUI(pbdir,acqui,inv)
        elseif (pbdir%tool == 2) then
           CALL sub_pbdirect_FDFD_SOLVE_ACQUIaniso(pbdir,acqui,inv)
        endif
        !STEP 3 : COST COMPUTATION 
        if(mype==0)  then
           !Compute cost function 
           CALL sub_pbdirect_FDFD_FCOST(pbdir,inv,acqui,iw)           
        endif
        !Deallocation       
        CALL sub_pbdirect_FDFD_FINISH(pbdir) ! deallocation of RHS     
     enddo !end loop on the frequencies
     
     if(mype.eq.0) then        
        if( inv%regul.eq.1) then           
           call sub_Tikhonov_fcost(pbdir,inv)
        endif
        !Prior regularization
        if(inv%lambda_prior0>0.) then
           call sub_prior_fcost(pbdir,inv)
        endif
        !Scaling
        inv%fcost=inv%fcost/inv%scalingfactor                 
     endif
     
  endif
end subroutine sub_modeling_FWI_fcost

!-------------------------------------------------------------------!
! SUB_MODELING_FWI_grad :                                           !
! 1.Factorization of the modeling matrix and                        !
! 2. Computation of the incident wavefield for all sources          !
! 3. Computation of the adjoint wavefield for all sources           !
! 4. Computation of the cost function                               !
! 5. Computation of the gradient                                    !
!-------------------------------------------------------------------!
subroutine sub_modeling_FWI_grad(pbdir,acqui,inv,inc_wfld,adj_wfld)

  implicit none
  include "mpif.h"
#include "common.h"
#include "pbdirect.h"
  include "acqui.h"
  include "inversion.h"

  !PBDIRECT
  TYPE (pbdirect) :: pbdir
  !INVERSION
  TYPE (inversion) :: inv
  !ACQUI
  TYPE(acquisition) :: acqui

  INTEGER :: i1,i2,isrc,k,iw,i
  complex,dimension(pbdir%n1e_s,pbdir%n2e_s,pbdir%nb_wfld_s,acqui%nsrc,pbdir%nw_s) :: inc_wfld
  complex,dimension(pbdir%n1e_s,pbdir%n2e_s,pbdir%nb_wfld_s,acqui%nsrc,pbdir%nw_s) :: adj_wfld
  complex,dimension(:,:,:,:,:),allocatable :: ax,az  
  double precision:: t1,t2
  real :: fcost_bathy

  if(mype.eq.0) then
     inv%gradient(:,:,:)=0.
     inv%fcost=0.   
     inv%residual(:)=0.     
  end if

  IF ( (pbdir%tool.eq.1) .OR. (pbdir%tool.eq.2) ) THEN
     ! ACOUSTIC (AN)ISOTROPIC MODELING IN THE FREQUENCY DOMAIN
     DO iw=1,pbdir%nw                
        !STEP 1 : LU FACTORIZATION OF S(p)
        if (pbdir%tool == 1) then !Isotropic
           CALL sub_pbdirect_FDFD_MATRIX(pbdir,pbdir%freqw(iw))           
        elseif (pbdir%tool == 2)  then !Anisotropic
           CALL sub_pbdirect_FDFD_MATRIXaniso(pbdir,pbdir%freqw(iw))
        endif
        !STEP 2 : COMPUTATION OF THE INCIDENT WAVEFIELD
        pbdir%type_RHS=1  
        if (pbdir%tool == 1) then           
           CALL sub_pbdirect_FDFD_SOLVE_ACQUI(pbdir,acqui,inv)           
        elseif (pbdir%tool == 2) then           
           CALL sub_pbdirect_FDFD_SOLVE_ACQUIaniso(pbdir,acqui,inv)
        endif
        if(mype.eq.0)  then
           !Compute cost function
           CALL sub_pbdirect_FDFD_FCOST(pbdir,inv,acqui,iw)           
           !Store incident wavefield           
           k=1
           do isrc=1,acqui%nsrc
              do i=1,pbdir%nb_wfld
                 do i2=1,pbdir%n2e
                    do i1=1,pbdir%n1e
                       inc_wfld(i1,i2,i,isrc,iw)=pbdir%id%RHS(k)
                       k=k+1                    
                    enddo
                 enddo
              enddo
           enddo
        end if
        CALL sub_pbdirect_FDFD_FINISH(pbdir) ! deallocation of RHS

        !STEP 3 : COMPUTATION OF THE ADJOINT WAVEFIELD        
        pbdir%type_RHS=2
        if (pbdir%tool == 1) then
           CALL sub_pbdirect_FDFD_SOLVE_ACQUI(pbdir,acqui,inv)
        elseif (pbdir%tool == 2) then
           CALL sub_pbdirect_FDFD_SOLVE_ACQUIaniso(pbdir,acqui,inv)
        endif
        !STEP 4 : GRADIENT BUILDING
        if(mype.eq.0)  then
           !Store adjoint wavefield
           k=1
           do isrc=1,acqui%nsrc
              do i=1,pbdir%nb_wfld
                 do i2=1,pbdir%n2e
                    do i1=1,pbdir%n1e
                       adj_wfld(i1,i2,i,isrc,iw)=pbdir%id%RHS(k)
                       k=k+1
                    enddo
                 enddo
              enddo
           enddo
           !Gradient
           CALL sub_FDFD_gradient(pbdir,inv,pbdir%nb_wfld,acqui,&
                inc_wfld(:,:,:,:,iw),adj_wfld(:,:,:,:,iw),pbdir%freqw(iw),iw)           
           CALL sub_pbdirect_FDFD_FINISH(pbdir) ! deallocation of RHS           
        endif
     enddo

     !Tikhonov regularization
     if(mype.eq.0) then        
        if( inv%regul.eq.1) then
           call sub_Tikhonov_fcost(pbdir,inv)
           call sub_Tikhonov_fgrad(pbdir,inv)
        endif
        !Prior regularization
        if(inv%lambda_prior0>0.) then
           call sub_prior_fcost(pbdir,inv)
           call sub_prior_fgrad(pbdir,inv)
        endif

        if(inv%firstgrad==0) THEN
           !compute the scaling if first call
           inv%scalingfactor=1.
           if(inv%invpar(1)/=4 .AND. inv%invpar(1)/=34 .AND. inv%invpar(1)/=5 .AND. inv%invpar(1)/=35  ) THEN ! case  no epsilon and delta
              inv%scalingfactor = 10e-5 * &
                   SQRT(MAX(SUM(inv%model(:,:,1)*(inv%model(:,:,1))),&
                   1.e3*pbdir%n1*pbdir%n2)/(SUM(inv%gradient(:,:,1)*inv%gradient(:,:,1))))
           ELSE !case epsilon or delta
              inv%scalingfactor = 5e-3 * &
                   SQRT(MAX(SUM(inv%model(:,:,1)*(inv%model(:,:,1))),&
                   0.1*pbdir%n1*pbdir%n2)/(SUM(inv%gradient(:,:,1)*inv%gradient(:,:,1))))
           END IF
           inv%scalingfactor=1./inv%scalingfactor     
           inv%firstgrad=1
        END IF

        !Final scaling
        inv%fcost=inv%fcost/inv%scalingfactor        
        inv%gradient(:,:,:)=inv%gradient(:,:,:)/inv%scalingfactor
        if(inv%regul.eq.2) then  
           !apply smoothing
           CALL sub_precond2(inv,pbdir%vp,inv%gradient,pbdir%n1,pbdir%n2,pbdir%h,pbdir%freqw(pbdir%nw),inv%npar,-1) 
           !force point above bathy to keep 0 gradient       
           do i2=1,pbdir%n2
              do i1=1,inv%ibathy(i2)-1
                 inv%gradient(i1,i2,:)=0.
              enddo
           enddo
        endif
        !Final: not first grad anymore
        if(inv%firstgrad.eq.0) then
           inv%firstgrad=1
        endif
     endif
  endif


  OPEN(12,file='gradient',access='direct',recl=pbdir%n1*pbdir%n2*4*inv%npar)
  WRITE(12,rec=1)   inv%gradient
  CLOSE(12)

end subroutine sub_modeling_FWI_grad


!-------------------------------------------------------------------!
! SUB_MODELING_FWI_GN : Computation of  Hessian-vector products in  !
!  the Gauss-Newton approximation                                   !
! (isotropic and anisotropic acoustic modeling)                     !
!-------------------------------------------------------------------!
subroutine sub_modeling_FWI_GN(pbdir,acqui,inv,inc_wfld,v,Hv)
  
  implicit none
  include "mpif.h"
#include "common.h"
#include "pbdirect.h"
  include "acqui.h"
  include "inversion.h"
  
  !PBDIRECT
  TYPE (pbdirect) :: pbdir
  !INVERSION
  TYPE (inversion) :: inv
  !ACQUI
  TYPE(acquisition) :: acqui
  
  !IN
  REAL,dimension(pbdir%n1*pbdir%n2*inv%npar) :: v
  complex,dimension(pbdir%n1e_s,pbdir%n2e_s,pbdir%nb_wfld_s,acqui%nsrc,pbdir%nw_s) :: inc_wfld
  !IN/OUT
  REAL,dimension(pbdir%n1*pbdir%n2*inv%npar) :: Hv
  !Local variables
  integer :: i,i1,i2,isrc,k,iw,ipar
  COMPLEX :: omegac
  real(KIND=8) ::time1,time2
  complex :: xi,H_term
  complex,dimension(:),allocatable :: fld
  complex,allocatable,dimension(:,:,:,:) :: lambda
  
  !Allocate RHS, initialize Hv
  if(mype.eq.0) then
     allocate (pbdir%id%RHS(pbdir%n1e*pbdir%n2e*acqui%nsrc*pbdir%nb_wfld))
     Hv(:)=0.
  endif
  
  !Complex i
  xi=CMPLX(0.,1.)  
  allocate(fld(5))
  do iw=1,pbdir%nw !frequency loop     
     !If nw>1 then LU FACTORIZATION OF S(p)
     if(pbdir%nw>1) then
        if (pbdir%tool == 1) then !Isotropic
           CALL sub_pbdirect_FDFD_MATRIX(pbdir,pbdir%freqw(iw))           
        elseif (pbdir%tool == 2)  then !Anisotropic
           CALL sub_pbdirect_FDFD_MATRIXaniso(pbdir,pbdir%freqw(iw))
        endif
     endif
     
     !Computation of J_s*v for all s     
     pbdir%type_RHS=1
     call sub_pbdirect_FDFD_SOLVE_Jv(pbdir,acqui,inv,&
          inc_wfld(:,:,:,:,iw),v,pbdir%freqw(iw))
     
     !Computation of R_s^{T}R_sJ_s*v for all s
     if(mype.eq.0) then
        if(pbdir%hicks.eq.0) then
           call sub_restriction(pbdir,acqui)        
        elseif(pbdir%hicks.eq.1) then
           call sub_restriction_Hicks(pbdir,acqui)        
        endif
     endif
     
     !Computation of lambda (solve S(p)^{T}lambda=R^{T}RJ*v)
     pbdir%id%ICNTL(9)  = 0 ! adjoint
     pbdir%id%JOB = 3
     if(mype.eq.0) then        
        do k=1,pbdir%n1e*pbdir%n2e*acqui%nsrc*pbdir%nb_wfld
           pbdir%id%RHS(k)=-CONJG(pbdir%id%RHS(k))/(pbdir%h**2)
        enddo
        pbdir%id%NRHS = acqui%nsrc
        pbdir%id%LRHS = pbdir%id%N
     endif     
#ifdef DOUBLEPRECISION
     CALL zmumps( pbdir%id )     
#else
     CALL cmumps( pbdir%id ) 
#endif
     
     !Computation of Hv
     if(mype.eq.0)  then                      
        allocate(lambda(pbdir%n1e,pbdir%n2e,pbdir%nb_wfld,acqui%nsrc))
        omegac=CMPLX(pbdir%freqw(iw),pbdir%slaplace)
        k=1
        do isrc=1,acqui%nsrc
           do i=1,pbdir%nb_wfld
              do i2=1,pbdir%n2e
                 do i1=1,pbdir%n1e                       
                    lambda(i1,i2,i,isrc)=pbdir%id%RHS(k)
                    k=k+1
                 enddo
              enddo
           enddo
        enddo        
        do isrc=1,acqui%nsrc
           do i=1,pbdir%nb_wfld
              do ipar=1,inv%npar
                 do i2=1,pbdir%n2
                    do i1=1,pbdir%n1
                       if(inv%invpar(ipar).eq.1) then !VP
                          Hv((ipar-1)*pbdir%n1*pbdir%n2+(i2-1)*pbdir%n1+i1)=&
                               Hv((ipar-1)*pbdir%n1*pbdir%n2+(i2-1)*pbdir%n1+i1)&
                               -REAL(2.*omegac**2/(pbdir%rho(i1,i2)*pbdir%vp(i1,i2)**3*&
                               (1.-0.5*xi/pbdir%qp(i1,i2))**2)*&
                               inc_wfld(pbdir%npml+i1,pbdir%npml+i2,i,isrc,iw)*&
                               lambda(pbdir%npml+i1,pbdir%npml+i2,i,isrc))                               
                       elseif(inv%invpar(ipar).eq.2) then !RHO
                          call sub_drho5_fld(pbdir%nb_wfld,acqui%nsrc,isrc,i1,i2,i,omegac,&
                               pbdir,inc_wfld(:,:,:,:,iw),fld)
                          H_term=&
                               fld(1)*lambda(pbdir%npml+i1,pbdir%npml+i2,1,isrc)+&
                               fld(2)*lambda(pbdir%npml+i1-1,pbdir%npml+i2,1,isrc)+&
                               fld(3)*lambda(pbdir%npml+i1+1,pbdir%npml+i2,1,isrc)+&
                               fld(4)*lambda(pbdir%npml+i1,pbdir%npml+i2-1,1,isrc)+&
                               fld(5)*lambda(pbdir%npml+i1,pbdir%npml+i2+1,1,isrc)
                          Hv((ipar-1)*pbdir%n1*pbdir%n2+(i2-1)*pbdir%n1+i1)=&
                               Hv((ipar-1)*pbdir%n1*pbdir%n2+(i2-1)*pbdir%n1+i1)&
                               +REAL(H_term)                               
                       elseif(inv%invpar(ipar).eq.3) then !QP
                          Hv((ipar-1)*pbdir%n1*pbdir%n2+(i2-1)*pbdir%n1+i1)=&
                               Hv((ipar-1)*pbdir%n1*pbdir%n2+(i2-1)*pbdir%n1+i1)&
                               -REAL(xi*omegac**2/(pbdir%rho(i1,i2)*pbdir%vp(i1,i2)**2*&
                               (1.-0.5*xi/pbdir%qp(i1,i2))**3*&
                               pbdir%qp(i1,i2)**2)*&
                               inc_wfld(pbdir%npml+i1,pbdir%npml+i2,i,isrc,iw)*&
                               lambda(pbdir%npml+i1,pbdir%npml+i2,i,isrc))                               
                       endif
                    enddo
                 enddo
              enddo
           enddo
           Hv(:)=Hv(:)*inv%cc1(isrc,iw)/(ABS(inv%cc1(isrc,iw)))                                            
        enddo
        deallocate(lambda)        
     endif
  enddo
  deallocate(fld)

  !REGULARIZATION
  if(mype.eq.0) then
     !Tikhonov regularization
     if(inv%regul.eq.1) then
        call sub_Tikhonov_Hv(pbdir,inv,Hv,v)
     endif
     !Prior regularization 
     if(inv%lambda_prior0>0.) then
        call sub_prior_Hv(pbdir,inv,Hv,v)           
     endif
     !Scaling
     Hv(:)=Hv(:)/inv%scalingfactor                   
     deallocate (pbdir%id%RHS)
  end if
  
end subroutine sub_modeling_FWI_GN

!-------------------------------------------------------------------!
! SUB_MODELING_FWI_EN : Computation of  Hessian-vector products in  !
!  the Exact-Newton case (acoustic isotropic and anisotropic case)  !                       
!-------------------------------------------------------------------!
subroutine sub_modeling_FWI_EN(pbdir,acqui,inv,inc_wfld,adj_wfld,v,Hv)
  
  implicit none
  include "mpif.h"
#include "common.h"
#include "pbdirect.h"
  include "acqui.h"
  include "inversion.h"
  
  !PBDIRECT
  TYPE (pbdirect) :: pbdir
  !INVERSION
  TYPE (inversion) :: inv
  !ACQUI
  TYPE(acquisition) :: acqui
  
  !IN
  real,dimension(pbdir%n1*pbdir%n2*inv%npar) :: v
  complex,dimension(pbdir%n1e_s,pbdir%n2e_s,pbdir%nb_wfld_s,acqui%nsrc,pbdir%nw_s) :: inc_wfld
  complex,dimension(pbdir%n1e_s,pbdir%n2e_s,pbdir%nb_wfld_s,acqui%nsrc,pbdir%nw_s) :: adj_wfld
  !IN/OUT
  real,dimension(pbdir%n1*pbdir%n2*inv%npar) :: Hv
  !Local Variables
  integer :: i,i1,i2,isrc,k,iw,ipar,ipar2
  COMPLEX :: omegac
  real(KIND=8) ::time1,time2
  complex :: xi,H_term
  complex,dimension(:),allocatable :: fld
  complex,allocatable,dimension(:,:,:,:) :: lambda
  complex,allocatable,dimension(:,:,:,:) :: alpha
  
  !Allocate RHS
  if(mype.eq.0) then
     allocate (pbdir%id%RHS(pbdir%n1e*pbdir%n2e*acqui%nsrc*pbdir%nb_wfld))
     Hv(:)=0.
  endif
  
  !Complex i
  xi=CMPLX(0.,1.)
  
  allocate(fld(5))
  do iw=1,pbdir%nw !frequency loop
     
     !If nw>1 then LU FACTORIZATION OF S(p)
     if(pbdir%nw>1) then
        if (pbdir%tool == 1) then !Isotropic
           CALL sub_pbdirect_FDFD_MATRIX(pbdir,pbdir%freqw(iw))           
        elseif (pbdir%tool == 2)  then !Anisotropic
           CALL sub_pbdirect_FDFD_MATRIXaniso(pbdir,pbdir%freqw(iw))
        endif
     endif
     
     !Computation of J_s*v for all s     
     pbdir%type_RHS=1
     call sub_pbdirect_FDFD_SOLVE_Jv(pbdir,acqui,inv,&
          inc_wfld(:,:,:,:,iw),v,pbdir%freqw(iw))
     
     !Storing of J_s*v (ie alpha_s)
     if(mype.eq.0) then
        allocate(alpha(pbdir%n1e,pbdir%n2e,pbdir%nb_wfld,acqui%nsrc))
        k=1
        do isrc=1,acqui%nsrc
           do i=1,pbdir%nb_wfld
              do i2=1,pbdir%n2e
                 do i1=1,pbdir%n1e
                    alpha(i1,i2,i,isrc)=pbdir%id%RHS(k)
                    k=k+1
                 enddo
              enddo
           enddo
        enddo
     endif     
     !Computation of RHS for lambda_s     
     if(mype.eq.0) then
        call sub_fill_RHS_lambda(pbdir,acqui,inv,adj_wfld(:,:,:,:,iw),v,pbdir%freqw(iw))
     endif
     
     !Computation of lambda_s for all s
     pbdir%id%ICNTL(9)  = 0 ! solve A^{T}x=b
     pbdir%id%JOB = 3     
#ifdef DOUBLEPRECISION
     CALL zmumps( pbdir%id )     
#else
     CALL cmumps( pbdir%id ) 
#endif
     
     !Lambda_s copy from RHS
     if(mype.eq.0)  then              
        allocate(lambda(pbdir%n1e,pbdir%n2e,pbdir%nb_wfld,acqui%nsrc))
        omegac=CMPLX(pbdir%freqw(iw),pbdir%slaplace)
        k=1
        do isrc=1,acqui%nsrc
           do i=1,pbdir%nb_wfld
              do i2=1,pbdir%n2e
                 do i1=1,pbdir%n1e
                    lambda(i1,i2,i,isrc)=pbdir%id%RHS(k)
                    k=k+1
                 enddo
              enddo
           enddo
        enddo
     endif
     
     !Computation of Hv
     if(mype.eq.0)  then              
        do isrc=1,acqui%nsrc
           do i=1,pbdir%nb_wfld
              do ipar=1,inv%npar
                 do i2=1,pbdir%n2
                    do i1=1,pbdir%n1
                       if(inv%invpar(ipar).eq.1) then
                          Hv((ipar-1)*pbdir%n1*pbdir%n2+(i2-1)*pbdir%n1+i1)=&
                               Hv((ipar-1)*pbdir%n1*pbdir%n2+(i2-1)*pbdir%n1+i1)&
                               -REAL(2.*omegac**2/(pbdir%rho(i1,i2)*pbdir%vp(i1,i2)**3*&
                               (1.-0.5*xi/pbdir%qp(i1,i2))**2)*&
                               inc_wfld(pbdir%npml+i1,pbdir%npml+i2,i,isrc,iw)*&
                               lambda(pbdir%npml+i1,pbdir%npml+i2,i,isrc))&
                               -REAL(2.*omegac**2/(pbdir%rho(i1,i2)*pbdir%vp(i1,i2)**3*&
                               (1.-0.5*xi/pbdir%qp(i1,i2))**2)*&
                               alpha(pbdir%npml+i1,pbdir%npml+i2,i,isrc)*&
                               adj_wfld(pbdir%npml+i1,pbdir%npml+i2,i,isrc,iw))                               
                          do ipar2=1,inv%npar
                             if(inv%invpar(ipar2).eq.1) then !VP_VP
                                Hv((ipar-1)*pbdir%n1*pbdir%n2+(i2-1)*pbdir%n1+i1)=&
                                     Hv((ipar-1)*pbdir%n1*pbdir%n2+(i2-1)*pbdir%n1+i1)&
                                     +v((ipar2-1)*pbdir%n1*pbdir%n2+(i2-1)*pbdir%n1+i1)&
                                     *REAL(6.*omegac**2/(pbdir%rho(i1,i2)*&
                                     pbdir%vp(i1,i2)**4*&
                                     (1.-0.5*xi/pbdir%qp(i1,i2))**2)*&
                                     inc_wfld(pbdir%npml+i1,pbdir%npml+i2,i,isrc,iw)*&
                                     adj_wfld(pbdir%npml+i1,pbdir%npml+i2,i,isrc,iw))/&
                                     (pbdir%h**2)                                
                             elseif(inv%invpar(ipar2).eq.2) then !VP_RHO
                                Hv((ipar-1)*pbdir%n1*pbdir%n2+(i2-1)*pbdir%n1+i1)=&
                                     Hv((ipar-1)*pbdir%n1*pbdir%n2+(i2-1)*pbdir%n1+i1)&
                                     +v((ipar2-1)*pbdir%n1*pbdir%n2+(i2-1)*pbdir%n1+i1)&
                                     *REAL(2.*omegac**2/(pbdir%rho(i1,i2)**2*&
                                     pbdir%vp(i1,i2)**3*&
                                     (1.-0.5*xi/pbdir%qp(i1,i2))**2)*&
                                     inc_wfld(pbdir%npml+i1,pbdir%npml+i2,i,isrc,iw)*&
                                     adj_wfld(pbdir%npml+i1,pbdir%npml+i2,i,isrc,iw))/&                                     
                                     (pbdir%h**2)
                             elseif(inv%invpar(ipar2).eq.3) then !VP_QP
                                Hv((ipar-1)*pbdir%n1*pbdir%n2+(i2-1)*pbdir%n1+i1)=&
                                     Hv((ipar-1)*pbdir%n1*pbdir%n2+(i2-1)*pbdir%n1+i1)&
                                     +v((ipar2-1)*pbdir%n1*pbdir%n2+(i2-1)*pbdir%n1+i1)&
                                     *REAL(2.*xi*omegac**2/(pbdir%rho(i1,i2)**2*&
                                     pbdir%vp(i1,i2)**3*&
                                     (1.-0.5*xi/pbdir%qp(i1,i2))**2*&
                                     pbdir%qp(i1,i2)**2)*&
                                     inc_wfld(pbdir%npml+i1,pbdir%npml+i2,i,isrc,iw)*&
                                     adj_wfld(pbdir%npml+i1,pbdir%npml+i2,i,isrc,iw))/&                                     
                                     (pbdir%h**2)
                             endif
                          enddo
                       elseif(inv%invpar(ipar).eq.2) then                       
                          call sub_drho5_fld(pbdir%nb_wfld,acqui%nsrc,isrc,i1,i2,i,omegac,&
                               pbdir,inc_wfld(:,:,:,:,iw),fld)
                          H_term=&
                               fld(1)*lambda(pbdir%npml+i1,pbdir%npml+i2,i,isrc)+&
                               fld(2)*lambda(pbdir%npml+i1-1,pbdir%npml+i2,i,isrc)+&
                               fld(3)*lambda(pbdir%npml+i1+1,pbdir%npml+i2,i,isrc)+&
                               fld(4)*lambda(pbdir%npml+i1,pbdir%npml+i2-1,i,isrc)+&
                               fld(5)*lambda(pbdir%npml+i1,pbdir%npml+i2+1,i,isrc)
                          call sub_drho5_fld(pbdir%nb_wfld,acqui%nsrc,isrc,i1,i2,i,omegac,&
                               pbdir,alpha,fld)
                          H_term=H_term+&
                               fld(1)*adj_wfld(pbdir%npml+i1,pbdir%npml+i2,i,isrc,iw)+&
                               fld(2)*adj_wfld(pbdir%npml+i1-1,pbdir%npml+i2,i,isrc,iw)+&
                               fld(3)*adj_wfld(pbdir%npml+i1+1,pbdir%npml+i2,i,isrc,iw)+&
                               fld(4)*adj_wfld(pbdir%npml+i1,pbdir%npml+i2-1,i,isrc,iw)+&
                               fld(5)*adj_wfld(pbdir%npml+i1,pbdir%npml+i2+1,i,isrc,iw)
                          Hv((ipar-1)*pbdir%n1*pbdir%n2+(i2-1)*pbdir%n1+i1)=&
                               Hv((ipar-1)*pbdir%n1*pbdir%n2+(i2-1)*pbdir%n1+i1)&
                               +REAL(H_term)  
                          do ipar2=1,inv%npar
                             if(inv%invpar(ipar2).eq.1) then !RHO_VP
                                Hv((ipar-1)*pbdir%n1*pbdir%n2+(i2-1)*pbdir%n1+i1)=&
                                     Hv((ipar-1)*pbdir%n1*pbdir%n2+(i2-1)*pbdir%n1+i1)&
                                     +v((ipar2-1)*pbdir%n1*pbdir%n2+(i2-1)*pbdir%n1+i1)&
                                     *REAL(2.*omegac**2/(pbdir%rho(i1,i2)**2*&
                                     pbdir%vp(i1,i2)**3*&
                                     (1.-0.5*xi/pbdir%qp(i1,i2))**2)*&
                                     inc_wfld(pbdir%npml+i1,pbdir%npml+i2,i,isrc,iw)*&
                                     adj_wfld(pbdir%npml+i1,pbdir%npml+i2,i,isrc,iw))/&                                     
                                     (pbdir%h**2)
                             elseif(inv%invpar(ipar2).eq.2) then !RHO_RHO
                                call sub_d2rho5_fld(pbdir%nb_wfld,acqui%nsrc,isrc,&
                                     i1,i2,i,omegac,pbdir,inc_wfld(:,:,:,:,iw),fld)
                                H_term=&
                                     fld(1)*&
                                     adj_wfld(pbdir%npml+i1,pbdir%npml+i2,i,isrc,iw)+&
                                     fld(2)*&
                                     adj_wfld(pbdir%npml+i1-1,pbdir%npml+i2,i,isrc,iw)+&
                                     fld(3)*&
                                     adj_wfld(pbdir%npml+i1+1,pbdir%npml+i2,i,isrc,iw)+&
                                     fld(4)*&
                                     adj_wfld(pbdir%npml+i1,pbdir%npml+i2-1,i,isrc,iw)+&
                                     fld(5)*&
                                     adj_wfld(pbdir%npml+i1,pbdir%npml+i2+1,i,isrc,iw)
                                Hv((ipar-1)*pbdir%n1*pbdir%n2+(i2-1)*pbdir%n1+i1)=&
                                     Hv((ipar-1)*pbdir%n1*pbdir%n2+(i2-1)*pbdir%n1+i1)&
                                     +REAL(H_term)*&
                                     v((ipar2-1)*pbdir%n1*pbdir%n2+(i2-1)*pbdir%n1+i1)/&                                     
                                     (pbdir%h**2)
                             elseif(inv%invpar(ipar2).eq.3) then !RHO_QP
                                Hv((ipar-1)*pbdir%n1*pbdir%n2+(i2-1)*pbdir%n1+i1)=&
                                     Hv((ipar-1)*pbdir%n1*pbdir%n2+(i2-1)*pbdir%n1+i1)&
                                     +v((ipar2-1)*pbdir%n1*pbdir%n2+(i2-1)*pbdir%n1+i1)&
                                     *REAL(xi*omegac**2/(pbdir%rho(i1,i2)**2*&
                                     pbdir%vp(i1,i2)**2*&
                                     (1.-0.5*xi/pbdir%qp(i1,i2))**2*&
                                     pbdir%qp(i1,i2)**2)*&
                                     inc_wfld(pbdir%npml+i1,pbdir%npml+i2,i,isrc,iw)*&
                                     adj_wfld(pbdir%npml+i1,pbdir%npml+i2,i,isrc,iw))/&                                     
                                     (pbdir%h**2)
                             endif                             
                          enddo
                       elseif(inv%invpar(ipar).eq.3) then !QP
                          Hv((ipar-1)*pbdir%n1*pbdir%n2+(i2-1)*pbdir%n1+i1)=&
                               Hv((ipar-1)*pbdir%n1*pbdir%n2+(i2-1)*pbdir%n1+i1)&
                               -REAL(xi*omegac**2/(pbdir%rho(i1,i2)*pbdir%vp(i1,i2)**2*&
                               (1.-0.5*xi/pbdir%qp(i1,i2))**2*&
                               pbdir%qp(i1,i2)**2)*&
                               inc_wfld(pbdir%npml+i1,pbdir%npml+i2,i,isrc,iw)*&
                               lambda(pbdir%npml+i1,pbdir%npml+i2,i,isrc))&                               
                               -REAL(xi*omegac**2/(pbdir%rho(i1,i2)*pbdir%vp(i1,i2)**2*&
                               (1.-0.5*xi/pbdir%qp(i1,i2))**2*&
                               pbdir%qp(i1,i2)**2)*&
                               alpha(pbdir%npml+i1,pbdir%npml+i2,i,isrc)*&
                               adj_wfld(pbdir%npml+i1,pbdir%npml+i2,i,isrc,iw))                               
                          do ipar2=1,inv%npar
                             if(inv%invpar(ipar2).eq.1) then ! QP_VP
                                Hv((ipar-1)*pbdir%n1*pbdir%n2+(i2-1)*pbdir%n1+i1)=&
                                     Hv((ipar-1)*pbdir%n1*pbdir%n2+(i2-1)*pbdir%n1+i1)&
                                     +v((ipar2-1)*pbdir%n1*pbdir%n2+(i2-1)*pbdir%n1+i1)&
                                     *REAL(2.*xi*omegac**2/(pbdir%rho(i1,i2)**2*&
                                     pbdir%vp(i1,i2)**3*&
                                     (1.-0.5*xi/pbdir%qp(i1,i2))**2*&
                                     pbdir%qp(i1,i2)**2)*&
                                     inc_wfld(pbdir%npml+i1,pbdir%npml+i2,i,isrc,iw)*&
                                     adj_wfld(pbdir%npml+i1,pbdir%npml+i2,i,isrc,iw))/&                                     
                                     (pbdir%h**2)
                             elseif(inv%invpar(ipar2).eq.2) then ! QP_RHO
                                Hv((ipar-1)*pbdir%n1*pbdir%n2+(i2-1)*pbdir%n1+i1)=&
                                     Hv((ipar-1)*pbdir%n1*pbdir%n2+(i2-1)*pbdir%n1+i1)&
                                     +v((ipar2-1)*pbdir%n1*pbdir%n2+(i2-1)*pbdir%n1+i1)&
                                     *REAL(xi*omegac**2/(pbdir%rho(i1,i2)**2*&
                                     pbdir%vp(i1,i2)**2*&
                                     (1.-0.5*xi/pbdir%qp(i1,i2))**2*&
                                     pbdir%qp(i1,i2)**2)*&
                                     inc_wfld(pbdir%npml+i1,pbdir%npml+i2,i,isrc,iw)*&
                                     adj_wfld(pbdir%npml+i1,pbdir%npml+i2,i,isrc,iw))/&                                     
                                     (pbdir%h**2)
                             elseif(inv%invpar(ipar2).eq.3) then ! QP_QP
                                Hv((ipar-1)*pbdir%n1*pbdir%n2+(i2-1)*pbdir%n1+i1)=&
                                     Hv((ipar-1)*pbdir%n1*pbdir%n2+(i2-1)*pbdir%n1+i1)&
                                     +v((ipar2-1)*pbdir%n1*pbdir%n2+(i2-1)*pbdir%n1+i1)&
                                     *REAL(&
                                     (2.*xi*omegac**2*pbdir%qp(i1,i2)*&
                                     (1-0.5*xi/pbdir%qp(i1,i2))-1.5*xi)/&
                                     (pbdir%rho(i1,i2)*pbdir%vp(i1,i2)**2*&
                                     (1.-0.5*xi/pbdir%qp(i1,i2))**4*&
                                     pbdir%qp(i1,i2)**4)*&
                                     inc_wfld(pbdir%npml+i1,pbdir%npml+i2,i,isrc,iw)*&
                                     adj_wfld(pbdir%npml+i1,pbdir%npml+i2,i,isrc,iw))/&                                     
                                     (pbdir%h**2)
                             endif
                          enddo
                       endif
                    enddo
                 enddo
              enddo
           enddo
           Hv(:)=Hv(:)*inv%cc1(isrc,iw)/ABS(inv%cc1(isrc,iw))           
        enddo
        deallocate(lambda)
        deallocate(alpha)        
     end if
  enddo
  
  !Regularization and scaling
  if(mype.eq.0) then
     if(inv%regul.eq.1) then
        call sub_Tikhonov_Hv(pbdir,inv,Hv,v)
     endif
     !Prior regularization 
     if(inv%lambda_prior0>0.) then
        call sub_prior_Hv(pbdir,inv,Hv,v)           
     endif
     !Scaling
     Hv(:)=Hv(:)/inv%scalingfactor
     !Deallocation of RHS
     deallocate(pbdir%id%RHS)
  endif
  
end subroutine sub_modeling_FWI_EN
