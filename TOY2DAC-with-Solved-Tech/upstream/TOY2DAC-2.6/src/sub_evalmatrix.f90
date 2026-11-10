!Copyright 2008-2011 SEISCOPE project, All rights reserved.
!Copyright 2013-2021 SEISCOPEII project, All rights reserved.

!
! subroutine subevalmatrix
!
! b: buoyancy
! mu: complex bulk modulus
! d1,d1b,d2,d2b: damping functions for PML
! mat,irn,icn: arrays for sparse impedance matrix storage
! n1 x n2: grid dimensions
! h: grid spacing
! ne: number of non zero coefficients in the matrix
! omegac: complex frequency

SUBROUTINE sub_evalmatrix(ve,b,mu,d1,d1b,d2,d2b,ne, &
     irn,icn,mat,n1e,n2e,h,omegac,nnz)

  IMPLICIT NONE
#include "common.h"
  integer ne,i2,i1,k,l,iparaxial,npml1,npml2,k0,i
  real h1,h,h2,h3
  real bpm,bmp,bpp,bmm,b00,bp0,bm0,b0p,b0m
  integer n1e,n2e,xi0,nnz,l0
  integer irn(ne),icn(ne)

  real ax,bx,cx,dx,ex,fx,gx,hx
  real az,bz,cz,dz,ez,fz,gz,hz,damp
  real wm1,wm2,wm3,w1,w2
  real b(0:n1e+1,0:n2e+1),ve(n1e,n2e)

  complex d1p,d1m,d1qp,d1qm,d2p,d2m,d2qp,d2qm
  complex omegac,omega2,xc0,m,ci
  complex d1(0:n1e+1),d1b(0:n1e+1),d2(0:n2e+1),d2b(0:n2e+1)
  complex mu(0:n1e+1,0:n2e+1)
  COMPLEX(KIND=creal) :: mat(ne)

  damp=0.
  wm1=0.6287326
  wm2=0.3712667
  wm3=1.-wm1-wm2
  wm2=0.25*wm2
  wm3=0.25*wm3
  w1=0.4382634
  w2=1.-w1

  ! set constants
  xi0=0
  xc0=cmplx(0.,0.)

  omega2=omegac*omegac
  h1=1./h
  h2=1./(h*h)
  h3=1./(h*h*h)
  ci=cmplx(0.,1.)
 
  ! Initialize sparse impedance matrix
  irn(:)=0.
  icn(:)=0.
  mat(:)=cmplx(0.,0.)
  
  nnz=0
  write(*,*) 'OMEGA2 = ',omega2
  write(*,*) 'N1 N2 = ',n1e,n2e
  k=0


  ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ! LOOP OVER ROWS OF THE IMPEDANCE MATRIX
  ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ax=1.
  bx=0.
  cx=0.
  dx=0.
  ex=0.
  fx=0.
  gx=0.
  hx=0.

  az=0.
  bz=1.
  cz=0.
  dz=0.
  ez=0.
  fz=0.
  gz=0.
  hz=0.

  do i2=1,n2e
     do i1=1,n1e
        call subbaverage(b,n1e,n2e,i1,i2,bpm,bmp,bpp,bmm,b00,bp0,bm0,b0p,b0m)
        d2p=d2b(i2)
        d2m=d2b(i2-1)
        d1p=d1b(i1)
        d1m=d1b(i1-1)
        
        ! Node 00
        k=(i2-1)*n1e+i1
        l=k
        nnz=nnz+1
        irn(nnz)=k
        icn(nnz)=l
        mat(nnz)= wm1*omega2/mu(i1,i2) &
             +w1*( &
             -0.25*h2*d2(i2)*ax*(bmp*d2p+bpm*d2m+bpp*d2p+bmm*d2m)  &
             +0.25*h2*d2(i2)*bx*(bmp*d1m+bpm*d1p-bpp*d1p-bmm*d1m)  &
             +0.25*h2*d1(i1)*az*(bmp*d2p+bpm*d2m-bpp*d2p-bmm*d2m)  &
             -0.25*h2*d1(i1)*bz*(bmp*d1m+bpm*d1p+bpp*d1p+bmm*d1m)) &
             +w2*(-d2(i2)*h2*ax*(d2p*b0p+d2m*b0m)-d1(i1)*h2*bz*(d1p*bp0+d1m*bm0) )
        mat(nnz)=mat(nnz)+damp*4*mat(nnz)

        ! Node 01  
        l=0
        IF (i2<n2e) then
           l=(i2-1+1)*n1e+i1
           nnz=nnz+1
           irn(nnz)=k
           icn(nnz)=l
           mat(nnz)= &
                wm2*omega2/mu(i1,i2+1) &
                +w1*( &
                0.25*h2*d2(i2)*ax*(bmp*d2p+bpp*d2p)  &
                +0.25*h2*d2(i2)*bx*(bmp*d1m-bpp*d1p) &
                +0.25*h2*d1(i1)*az*(-bmp*d2p+bpp*d2p) &
                +0.25*h2*d1(i1)*bz*(-bmp*d1m-bpp*d1p) ) &
                +w2*(d2(i2)*h2*ax*b0p*d2p + 0.25*d1(i1)*az*h2*(bp0*d2(i2)-bm0*d2(i2)) )
           mat(nnz)=mat(nnz)-damp*mat(nnz)
        end if
        
        ! Node 0-1 
        IF (i2>1) THEN
           l=(i2-1-1)*n1e+i1
           nnz=nnz+1
           irn(nnz)=k
           icn(nnz)=l
           mat(nnz)= &
                wm2*omega2/mu(i1,i2-1) &
                +w1*( &
                0.25*h2*d2(i2)*ax*(bpm*d2m+bmm*d2m) &
                +0.25*h2*d2(i2)*bx*(bpm*d1p-bmm*d1m) &
                +0.25*h2*d1(i1)*az*(-bpm*d2m+bmm*d2m) &
                +0.25*h2*d1(i1)*bz*(-bpm*d1p-bmm*d1m)) &
                +w2*(d2(i2)*h2*ax*b0m*d2m + 0.25*d1(i1)*az*h2*(-bp0*d2(i2)+bm0*d2(i2)) )
           mat(nnz)=mat(nnz)-damp*mat(nnz)

        end if

        ! Node -10 
        if (i1 > 1) then
           l=(i2-1)*n1e+i1-1
           nnz=nnz+1
           irn(nnz)=k
           icn(nnz)=l
           mat(nnz)= &
                wm2*omega2/mu(i1-1,i2) &
                +w1*( &
                0.25*h2*d2(i2)*ax*(-bmp*d2p-bmm*d2m) &
                +0.25*h2*d2(i2)*bx*(-bmp*d1m+bmm*d1m) &
                +0.25*h2*d1(i1)*az*(bmp*d2p-bmm*d2m) &
                +0.25*h2*d1(i1)*bz*(bmp*d1m+bmm*d1m)) &
                +w2*(d1(i1)*h2*bz*bm0*d1m + 0.25*d2(i2)*bx*h2*(-b0p*d1(i1)+b0m*d1(i1)) )
           mat(nnz)=mat(nnz)-damp*mat(nnz)
        end if

        ! Node 10        
        if (i1 < n1e) then
           l=(i2-1)*n1e+i1+1
           nnz=nnz+1
           irn(nnz)=k
           icn(nnz)=l
           mat(nnz)=   &
                wm2*omega2/mu(i1+1,i2) &
                +w1*( &
                0.25*h2*d2(i2)*ax*(-bpm*d2m-bpp*d2p) &
                +0.25*h2*d2(i2)*bx*(-bpm*d1p+bpp*d1p) &
                +0.25*h2*d1(i1)*az*(bpm*d2m-bpp*d2p) &
                +0.25*h2*d1(i1)*bz*(bpm*d1p+bpp*d1p)) &
                +w2*(d1(i1)*h2*bz*bp0*d1p + 0.25*d2(i2)*bx*h2*(-b0m*d1(i1)+b0p*d1(i1)) )
           mat(nnz)=mat(nnz)-damp*mat(nnz)
        end if

        ! Node -11 
        if (i1 > 1 .AND. i2 < n2e) then
           l=(i2-1+1)*n1e+i1-1
           nnz=nnz+1
           irn(nnz)=k
           icn(nnz)=l
           mat(nnz)= &
                wm3*omega2/mu(i1-1,i2+1) &
                +w1*( &
                0.25*h2*d2(i2)*ax*bmp*d2p &
                -0.25*h2*d2(i2)*bx*bmp*d1m &
                -0.25*h2*d1(i1)*az*bmp*d2p &
                +0.25*h2*d1(i1)*bz*bmp*d1m) &
                +w2*(-0.25*d2(i2)*h2*bx*b0p*d1(i1) - 0.25*h2*d1(i1)*az*d2(i2)*bm0 )
        end if
        
        ! Node +1-1         
        if (i1 < n1e .AND. i2 >1) then
           l=(i2-1-1)*n1e+i1+1
           nnz=nnz+1
           irn(nnz)=k
           icn(nnz)=l
           mat(nnz)= &
                wm3*omega2/mu(i1+1,i2-1) &
                +w1*( &
                0.25*h2*d2(i2)*ax*bpm*d2m &
                -0.25*h2*d2(i2)*bx*bpm*d1p &
                -0.25*h2*d1(i1)*az*bpm*d2m &
                +0.25*h2*d1(i1)*bz*bpm*d1p ) &
                +w2*(-0.25*d2(i2)*h2*bx*b0m*d1(i1) - 0.25*h2*d1(i1)*az*d2(i2)*bp0)
        end if
        
        ! Node +1+1         
        if (i1 < n1e .AND. i2 < n2e) then
           l=(i2-1+1)*n1e+i1+1
           nnz=nnz+1
           irn(nnz)=k
           icn(nnz)=l
           mat(nnz)= &
                wm3*omega2/mu(i1+1,i2+1) &
                +w1*( &
                0.25*h2*d2(i2)*ax*bpp*d2p &
                +0.25*h2*d2(i2)*bx*bpp*d1p &
                +0.25*h2*d1(i1)*az*bpp*d2p &
                +0.25*h2*d1(i1)*bz*bpp*d1p ) &
                +w2*(0.25*d2(i2)*h2*bx*b0p*d1(i1) +0.25*h2*d1(i1)*az*d2(i2)*bp0 )
        end if
        
        ! Node -1-1 
        if (i1 > 1 .AND. i2 > 1) then
           l=(i2-1-1)*n1e+i1-1
           nnz=nnz+1
           irn(nnz)=k
           icn(nnz)=l
           mat(nnz)= &
                wm3*omega2/mu(i1-1,i2-1) &
                +w1*( &
                0.25*h2*d2(i2)*ax*bmm*d2m &
                +0.25*h2*d2(i2)*bx*bmm*d1m &
                +0.25*h2*d1(i1)*az*bmm*d2m &
                +0.25*h2*d1(i1)*bz*bmm*d1m ) &
                +w2*(0.25*d2(i2)*h2*bx*b0m*d1(i1) + 0.25*h2*d1(i1)*az*d2(i2)*bm0 )           
        end if
     end do            !end do i1=1,n1
  end do            !end do i2=1,n2
  
  write(*,*) 'NNZ NE = ',nnz,ne
  
end SUBROUTINE sub_evalmatrix
