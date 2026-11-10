!Copyright 2008-2011 SEISCOPE project, All rights reserved.
!Copyright 2013-2021 SEISCOPEII project, All rights reserved.

SUBROUTINE sub_precond2(inv,vp_table,table,n1,n2,h,omega,npar,mode)


  IMPLICIT NONE
#include "common.h"
  include "inversion.h"

  !INVERSION
  TYPE (inversion) :: inv

  INTEGER :: n1,n2,npar
  REAL ::h,omega
  REAL,DIMENSION(n1,n2,npar) ::table 
  REAL,DIMENSION(n1,n2) ::vp_table
  REAL,ALLOCATABLE,DIMENSION(:,:):: table_tmp 

  INTEGER :: mode

  REAL,DIMENSION(:,:),ALLOCATABLE ::table_unit 
  INTEGER :: iipar

  ALLOCATE(table_unit(n1,n2))
  ALLOCATE(table_tmp(n1,n2))
  DO iipar=1,npar
     table_unit(:,:) = 1.
     !PRECOND table
     table_tmp(:,:)=table(:,:,iipar)
! OPEN(12,file='gradient_1',access='direct',recl=n1*n2*4)
 !       WRITE(12,rec=1)table_tmp
 !       close(12)
     CALL smooth_gauss_adaptiv ( table_tmp, table_unit, vp_table,  n2,n1, h, inv%lambda_x,inv%lambda_z,omega )
write(*,*)'in gaussian smoothing',n1,n2,h, inv%lambda_x,inv%lambda_z,omega   
 !   CALL applyprecond(inv,table(:,:,iipar),n1,n2,mode)
     !PRECOND table_unit
     !   CALL applyprecond(inv,table_unit,n1,n2,mode)
     !NORMALIZE table by table_unit for an UNIT operator
     !    WRITE(*,*)'MAXVAL(table_unit)',MAXVAL(table_unit),'MINVAL(table_unit)',MINVAL(table_unit)
     table(:,:,iipar)=table_tmp(:,:)
 !OPEN(12,file='gradient_2',access='direct',recl=n1*n2*4)
 !       WRITE(12,rec=1)table_tmp
 !       close(12)
  END DO
  DEALLOCATE(table_tmp)
END SUBROUTINE sub_precond2


 SUBROUTINE smooth_gauss_adaptiv ( table_cart, sig_par_cart, table_vp_cart,  nx, nz, h_cart, xleng,zleng,omega )

    
    IMPLICIT NONE

    !IN
    INTEGER  :: nx, nz
    REAL     :: h_cart, sigma_smooth_gradient 
    REAL            :: hx,hz,xleng,zleng,omega
    !IN/OUT
    REAL     :: table_cart ( nz, nx ), sig_par_cart ( nz , nx ), table_vp_cart( nz , nx )

    REAL,ALLOCATABLE,DIMENSION(:,:) ::tmp_scaling_table,tmp_sig_par_cart


!!$    open(60,file='valid_smooth',access='direct',recl=4*nx*nz)
!!$    write(60,rec=1)table_cart

    CALL subsmooth ( table_cart, sig_par_cart, table_vp_cart,  nx, nz, h_cart, xleng,zleng,omega )
    ALLOCATE(tmp_scaling_table(nz,nx),tmp_sig_par_cart(nz,nx))
    tmp_scaling_table(:,:)=1.
    tmp_sig_par_cart(:,:)=1.
    CALL subsmooth ( tmp_scaling_table, tmp_sig_par_cart, table_vp_cart,  nx, nz, h_cart, xleng,zleng,omega )

    !rescale the table_cart because the operator of the linear system is not normed...
    table_cart(:,:)=table_cart(:,:)/tmp_scaling_table(:,:)
    deallocate(tmp_scaling_table,tmp_sig_par_cart) 


  END SUBROUTINE smooth_gauss_adaptiv


  SUBROUTINE subsmooth ( table_cart, sig_par_cart, table_vp_cart,  nx, nz, h_cart, xleng,zleng,omega )

    IMPLICIT NONE
#    include "common.h"

    !IN
    INTEGER  :: nx, nz
    REAL     :: h_cart, sigma_smooth_gradient 
    REAL            :: hx,hz,xleng,zleng,omega
    !IN/OUT
    REAL     :: table_cart ( nz, nx ), sig_par_cart ( nz , nx ), table_vp_cart( nz , nx )

    !INTERN
    REAL,ALLOCATABLE,DIMENSION(:) :: dl,d,du,rhs
    INTEGER         :: i1, i2
    INTEGER         :: ier
    REAL            :: sigma

    !      sigma2 = sigma_smooth_gradient * sigma_smooth_gradient

    ! WRITE(*,*) 'SMOOTHING BASED ON GAUSSIAN IN CARTESIAN PROJECTED MODEL'



    !DIRECTION 1 SMOOTHING   NZ

    allocate(dl(nz-1))
    allocate(d(nz))
    allocate(du(nz-1))
    allocate(rhs(nz))
    hx=h_cart
    hz=h_cart

    DO i2 = 1 , nx
       !different value on edge to impose neumann conditions to avoid edge effect onthe smoothing

       i1=1
       sigma=table_vp_cart(i1,i2)*pi/omega*zleng*hz**2 !*h**2 because  zleng is normalmized at the beginning of code...
       !write(*,*)'sigma',sigma
       ! construction of tridiagonal system over iz
       !     d(i1)=1.0+2.*exp(-hz/sigma)/(exp(hz/sigma)-exp(-hz/sigma))-1./(exp(hz/sigma)-exp(-hz/sigma))
       !     e(i1)=-1./(exp(hz/sigma)-exp(-hz/sigma))
       !     rhs(i1)=table_cart(i1,i2)
       d(i1)=1.0+2.*exp(-hz/sigma)/(exp(hz/sigma)-exp(-hz/sigma))-1./(exp(hz/sigma)-exp(-hz/sigma))
       du(i1)=-1./(exp(hz/sigma)-exp(-hz/sigma))
       rhs(i1)=table_cart(i1,i2)

       DO i1 = 2 , nz-1

          sigma=table_vp_cart(i1,i2)*pi/omega*zleng*hz**2 !*h**2 because  zleng is normalmized at the beginning of code...

          ! construction of tridiagonal system over iz
          ! d(i1)=1.0+2.*exp(-hz/sigma)/(exp(hz/sigma)-exp(-hz/sigma))
          ! e(i1)=-1./(exp(hz/sigma)-exp(-hz/sigma))
          ! rhs(i1)=table_cart(i1,i2)
          d(i1)=1.0+2.*exp(-hz/sigma)/(exp(hz/sigma)-exp(-hz/sigma))
          du(i1)=-1./(exp(hz/sigma)-exp(-hz/sigma))
          dl(i1-1)=-1./(exp(hz/sigma)-exp(-hz/sigma))
          rhs(i1)=table_cart(i1,i2)
       ENDDO  !nz
       i1=nz
       sigma=table_vp_cart(i1,i2)*pi/omega*zleng*hz**2 !*h**2 because  zleng is normalmized at the beginning of code...

       ! construction of tridiagonal system over iz
       !  d(i1)=1.0+2.*exp(-hz/sigma)/(exp(hz/sigma)-exp(-hz/sigma))-1./(exp(hz/sigma)-exp(-hz/sigma))
       !   e(i1)=-1./(exp(hz/sigma)-exp(-hz/sigma))
       ! rhs(i1)=table_cart(i1,i2)

       d(i1)=1.0+2.*exp(-hz/sigma)/(exp(hz/sigma)-exp(-hz/sigma))-1./(exp(hz/sigma)-exp(-hz/sigma))
       dl(i1-1)=-1./(exp(hz/sigma)-exp(-hz/sigma))
       rhs(i1)=table_cart(i1,i2)


       CALL SGTSV( nz, 1, dl, d, du, rhs, nz, ier )

       !   call trisim(d,e,rhs,nz,ier)    ! solving the tridiagonal system
       if(ier /= 0 ) then
          WRITE (*,*) 'sub_cartesian_smooth: PROBLEM WITH DATA smoothing '
          WRITE(*,*) ' sigma_smooth_gradient and h_cart along NZ direction'
          WRITE (*,*) 'REVIEW PARAMETERS FOR SMOOTHING BECAUSE THERE IS NO POSSIBLITY TO SMOOTH'
          WRITE (*,*) 'sigma_smooth_gradient',sigma_smooth_gradient,'h_cart',h_cart
          WRITE(*,*)'PROGRAM STOPPED NOW'
          CALL flushc()
          CALL sub_MPI_ABORT()
       endif
       DO i1=1,nz
          sigma=table_vp_cart(i1,i2)*pi/omega*zleng*hz**2 !*h**2 because  zleng is normalmized at the beginning of code...
          table_cart(i1,i2)=rhs(i1)/sigma*0.39894228
       ENDDO ! nz
    ENDDO ! nx


    deallocate(d)
    deallocate(dl)
    deallocate(du)
    deallocate(rhs)

    !DIRECTION 2 SMOOTHING   NX

    allocate(d(nx))
    allocate(dl(nx-1))
    allocate(du(nx-1))
    allocate(rhs(nx))
    DO i1 = 1 , nz
       !different value on edge to impose neumann conditions to avoid edge effect onthe smoothing
       i2=1
       sigma=table_vp_cart(i1,i2)*pi/omega*xleng*hx**2 !*h**2 because  xleng is normalmized at the beginning of code...

       ! construction of tridiagonal system over i2
!!!$       d(i2)=1.0+2.*exp(-hx/sigma)/(exp(hx/sigma)-exp(-hx/sigma))-1./(exp(hx/sigma)-exp(-hx/sigma))
!!!$       e(i2)=-1./(exp(hx/sigma)-exp(-hx/sigma))
!!!$       rhs(i2)=table_cart(i1,i2)
       d(i2)=1.0+2.*exp(-hx/sigma)/(exp(hx/sigma)-exp(-hx/sigma))-1./(exp(hx/sigma)-exp(-hx/sigma))
       du(i2)=-1./(exp(hx/sigma)-exp(-hx/sigma))
       rhs(i2)=table_cart(i1,i2)
       DO i2 = 2 , nx-1

          sigma=table_vp_cart(i1,i2)*pi/omega*xleng*hx**2 !*h**2 because  xleng is normalmized at the beginning of code...

          ! construction of tridiagonal system over i2
!!!$          d(i2)=1.0+2.*exp(-hx/sigma)/(exp(hx/sigma)-exp(-hx/sigma))
!!!$          e(i2)=-1./(exp(hx/sigma)-exp(-hx/sigma))
!!!$          rhs(i2)=table_cart(i1,i2)
          d(i2)=1.0+2.*exp(-hx/sigma)/(exp(hx/sigma)-exp(-hx/sigma))
          du(i2)=-1./(exp(hx/sigma)-exp(-hx/sigma))
          dl(i2-1)=-1./(exp(hx/sigma)-exp(-hx/sigma))
          rhs(i2)=table_cart(i1,i2)
       ENDDO  !nx
       i2=nx
       sigma=table_vp_cart(i1,i2)*pi/omega*xleng*hx**2 !*h**2 because  xleng is normalmized at the beginning of code...

       ! construction of tridiagonal system over i2
!!!$       d(i2)=1.0+2.*exp(-hx/sigma)/(exp(hx/sigma)-exp(-hx/sigma))-1./(exp(hx/sigma)-exp(-hx/sigma))
!!!$       e(i2)=-1./(exp(hx/sigma)-exp(-hx/sigma))
!!!$       rhs(i2)=table_cart(i1,i2)
       d(i2)=1.0+2.*exp(-hx/sigma)/(exp(hx/sigma)-exp(-hx/sigma))-1./(exp(hx/sigma)-exp(-hx/sigma))
       dl(i2-1)=-1./(exp(hx/sigma)-exp(-hx/sigma))
       rhs(i2)=table_cart(i1,i2)

       CALL SGTSV( nx, 1, dl, d, du, rhs, nx, ier )

       !   call trisim(d,e,rhs,nx,ier)    ! solving the tridiagonal system
       if(ier /= 0 ) then
          WRITE (*,*) 'sub_cartesian_smooth: PROBLEM WITH DATA smoothing '
          WRITE(*,*) ' sigma_smooth_gradient and h_cart along NX direction'
          WRITE (*,*) 'REVIEW PARAMETERS FOR SMOOTHING BECAUSE THERE IS NO POSSIBLITY TO SMOOTH'
          WRITE (*,*) 'sigma_smooth_gradient',sigma_smooth_gradient,'h_cart',h_cart
          WRITE(*,*)'PROGRAM STOPPED NOW'
          CALL flushc()
          CALL sub_MPI_ABORT()
       endif
       DO i2=1,nx
          sigma=table_vp_cart(i1,i2)*pi/omega*xleng*hx**2 !*h**2 because  xleng is normalmized at the beginning of code...

          table_cart( i1 , i2)=sig_par_cart( i1 , i2 )*sig_par_cart( i1 , i2 )*rhs(i2)/sigma*0.39894228

       ENDDO ! nx
    ENDDO ! nz

    deallocate(d)
    deallocate(dl)
    deallocate(du)
    deallocate(rhs)


  END SUBROUTINE subsmooth

! !BELLOW ARE OLD IMPLEMENTATION

! !*************************************************************************
! !ACCORDING TO GUITTON 2010 SEG ABSTRACT, THE PRECONDTIONNING IS DEFINED BY
! ! grad(p)= S^t grad(m)
! ! m=Sp
! ! p=S-1 m
! ! where m is the classic domain
! !       p is the smart domain
! !       grad(p) is the gradient in the 'smart' domain
! !       grad(m) is the gradient in the classic domain
! !       S is the transformation matrix 
! !-------------------------------------------------------------------------
! ! this routine applies S or S-1 depending on the value mode
! ! mode = 1 -> S
! ! mode = -1 -> S-1
! !*************************************************************************
! SUBROUTINE sub_precond(inv,vp_table,table,n1,n2,h,omega,npar,mode)
  
  
!   IMPLICIT NONE
!   include "inversion.h"

!   !INVERSION
!   TYPE (inversion) :: inv

!   INTEGER :: n1,n2,npar
!   REAL ::h,omega
!   REAL,DIMENSION(n1,n2,npar) ::table 
!   REAL,DIMENSION(n1,n2) ::vp_table
!   REAL,ALLOCATABLE,DIMENSION(:,:):: table_tmp 

!   INTEGER :: mode

!   REAL,DIMENSION(:,:),ALLOCATABLE ::table_unit 
!   INTEGER :: iipar

!   ALLOCATE(table_unit(n1,n2))
!   ALLOCATE(table_tmp(n1,n2))
!   DO iipar=1,npar
!      table_unit(:,:) = 1.
!      !PRECOND table
!      table_tmp(:,:)=table(:,:,iipar)     
!      CALL smooth_gauss_adaptiv ( table_tmp, table_unit, vp_table,  n2,n1, h, inv%lambda_x,inv%lambda_z,omega )
!      write(*,*)'in gaussian smoothing',n1,n2,h, inv%lambda_x,inv%lambda_z,omega        
!      table(:,:,iipar)=table_tmp(:,:)
!   END DO
!   DEALLOCATE(table_tmp)
! END SUBROUTINE sub_precond


! SUBROUTINE smooth_gauss_adaptiv ( table_cart, sig_par_cart, table_vp_cart,  nx, nz, h_cart, xleng,zleng,omega )

!   IMPLICIT NONE

!   !IN
!   INTEGER  :: nx, nz
!   REAL     :: h_cart, sigma_smooth_gradient 
!   REAL            :: hx,hz,xleng,zleng,omega
!   !IN/OUT
!   REAL     :: table_cart ( nz, nx ), sig_par_cart ( nz , nx ), table_vp_cart( nz , nx )
  
!   REAL,ALLOCATABLE,DIMENSION(:,:) ::tmp_scaling_table,tmp_sig_par_cart
   
!   CALL subsmooth ( table_cart, sig_par_cart, table_vp_cart,  nx, nz, h_cart, xleng,zleng,omega )
!   ALLOCATE(tmp_scaling_table(nz,nx),tmp_sig_par_cart(nz,nx))
!   tmp_scaling_table(:,:)=1.
!   tmp_sig_par_cart(:,:)=1.
!   CALL subsmooth ( tmp_scaling_table, tmp_sig_par_cart, table_vp_cart,  nx, nz, h_cart, xleng,zleng,omega )  
!   !rescale the table_cart because the operator of the linear system is not normed...
!   table_cart(:,:)=table_cart(:,:)/tmp_scaling_table(:,:)
!   deallocate(tmp_scaling_table,tmp_sig_par_cart) 
  
! END SUBROUTINE smooth_gauss_adaptiv


! SUBROUTINE subsmooth ( table_cart, sig_par_cart, table_vp_cart,  nx, nz, h_cart, xleng,zleng,omega )
  
!   IMPLICIT NONE
!   include "common.h"
  
!   !IN
!   INTEGER  :: nx, nz
!   REAL     :: h_cart, sigma_smooth_gradient 
!   REAL            :: hx,hz,xleng,zleng,omega
!   !IN/OUT
!   REAL     :: table_cart ( nz, nx ), sig_par_cart ( nz , nx ), table_vp_cart( nz , nx )
  
!   !INTERN
!   REAL,ALLOCATABLE,DIMENSION(:) :: dl,d,du,rhs
!   INTEGER         :: i1, i2
!   INTEGER         :: ier
!   REAL            :: sigma
  
!   !DIRECTION 1 SMOOTHING   NZ
  
!   allocate(dl(nz-1))
!   allocate(d(nz))
!   allocate(du(nz-1))
!   allocate(rhs(nz))
!   hx=h_cart
!   hz=h_cart
  
!   DO i2 = 1 , nx
!      !different value on edge to impose neumann conditions to avoid edge effect onthe smoothing
     
!      i1=1
!      sigma=table_vp_cart(i1,i2)*pi/omega*zleng*hz**2
!      d(i1)=1.0+2.*exp(-hz/sigma)/(exp(hz/sigma)-exp(-hz/sigma))-1./(exp(hz/sigma)-exp(-hz/sigma))
!      du(i1)=-1./(exp(hz/sigma)-exp(-hz/sigma))
!      rhs(i1)=table_cart(i1,i2)
     
!      DO i1 = 2 , nz-1        
!         sigma=table_vp_cart(i1,i2)*pi/omega*zleng*hz**2 !
!         d(i1)=1.0+2.*exp(-hz/sigma)/(exp(hz/sigma)-exp(-hz/sigma))
!         du(i1)=-1./(exp(hz/sigma)-exp(-hz/sigma))
!         dl(i1-1)=-1./(exp(hz/sigma)-exp(-hz/sigma))
!         rhs(i1)=table_cart(i1,i2)
!      ENDDO  !nz
!      i1=nz
!      sigma=table_vp_cart(i1,i2)*pi/omega*zleng*hz**2 
!      d(i1)=1.0+2.*exp(-hz/sigma)/(exp(hz/sigma)-exp(-hz/sigma))-1./(exp(hz/sigma)-exp(-hz/sigma))
!      dl(i1-1)=-1./(exp(hz/sigma)-exp(-hz/sigma))
!      rhs(i1)=table_cart(i1,i2)
!      CALL SGTSV( nz, 1, dl, d, du, rhs, nz, ier )
!      if(ier /= 0 ) then
!         WRITE (*,*) 'sub_cartesian_smooth: PROBLEM WITH DATA smoothing '
!         WRITE(*,*) ' sigma_smooth_gradient and h_cart along NZ direction'
!         WRITE (*,*) 'REVIEW PARAMETERS FOR SMOOTHING BECAUSE THERE IS NO POSSIBLITY TO SMOOTH'
!         WRITE (*,*) 'sigma_smooth_gradient',sigma_smooth_gradient,'h_cart',h_cart
!         WRITE(*,*)'PROGRAM STOPPED NOW'
!         CALL flushc()
!         CALL sub_MPI_ABORT()
!      endif
!      DO i1=1,nz
!         sigma=table_vp_cart(i1,i2)*pi/omega*zleng*hz**2 
!         table_cart(i1,i2)=rhs(i1)/sigma*0.39894228
!      ENDDO ! nz
!   ENDDO ! nx
  
  
!   deallocate(d)
!   deallocate(dl)
!   deallocate(du)
!   deallocate(rhs)
  
!   !DIRECTION 2 SMOOTHING   NX
  
!   allocate(d(nx))
!   allocate(dl(nx-1))
!   allocate(du(nx-1))
!   allocate(rhs(nx))
!   DO i1 = 1 , nz
!      i2=1
!      sigma=table_vp_cart(i1,i2)*pi/omega*xleng*hx**2 
!      d(i2)=1.0+2.*exp(-hx/sigma)/(exp(hx/sigma)-exp(-hx/sigma))-1./(exp(hx/sigma)-exp(-hx/sigma))
!      du(i2)=-1./(exp(hx/sigma)-exp(-hx/sigma))
!      rhs(i2)=table_cart(i1,i2)
!      DO i2 = 2 , nx-1        
!         sigma=table_vp_cart(i1,i2)*pi/omega*xleng*hx**2
!         d(i2)=1.0+2.*exp(-hx/sigma)/(exp(hx/sigma)-exp(-hx/sigma))
!         du(i2)=-1./(exp(hx/sigma)-exp(-hx/sigma))
!         dl(i2-1)=-1./(exp(hx/sigma)-exp(-hx/sigma))
!         rhs(i2)=table_cart(i1,i2)
!      ENDDO  !nx
!      i2=nx
!      sigma=table_vp_cart(i1,i2)*pi/omega*xleng*hx**2 
!      d(i2)=1.0+2.*exp(-hx/sigma)/(exp(hx/sigma)-exp(-hx/sigma))-1./(exp(hx/sigma)-exp(-hx/sigma))
!      dl(i2-1)=-1./(exp(hx/sigma)-exp(-hx/sigma))
!      rhs(i2)=table_cart(i1,i2)     
!      CALL SGTSV( nx, 1, dl, d, du, rhs, nx, ier )
!      if(ier /= 0 ) then
!         WRITE (*,*) 'sub_cartesian_smooth: PROBLEM WITH DATA smoothing '
!         WRITE(*,*) ' sigma_smooth_gradient and h_cart along NX direction'
!         WRITE (*,*) 'REVIEW PARAMETERS FOR SMOOTHING BECAUSE THERE IS NO POSSIBLITY TO SMOOTH'
!         WRITE (*,*) 'sigma_smooth_gradient',sigma_smooth_gradient,'h_cart',h_cart
!         WRITE(*,*)'PROGRAM STOPPED NOW'
!         CALL flushc()
!         CALL sub_MPI_ABORT()
!      endif
!      DO i2=1,nx
!         sigma=table_vp_cart(i1,i2)*pi/omega*xleng*hx**2         
!         table_cart( i1 , i2)=sig_par_cart( i1 , i2 )*sig_par_cart( i1 , i2 )*rhs(i2)/sigma*0.39894228
        
!      ENDDO ! nx
!   ENDDO ! nz
  
!   deallocate(d)
!   deallocate(dl)
!   deallocate(du)
!   deallocate(rhs)
  
  
! END SUBROUTINE subsmooth
