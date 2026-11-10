!Copyright 2008-2011 SEISCOPE project, All rights reserved.
!Copyright 2013-2021 SEISCOPEII project, All rights reserved.
! ------------------------------------------------------------------------------------------------------
! subroutine subevalmatrix
! Build impedance matrix for various anisotropic media using the equation
! system of Zhou, SEG, 2006.
! Rotated staggered grid stencil is implemented.
! 
! version from S. Operto Modified by R. Brossier 10/2011
!
! ianiso=1: isotropic medium
! ianiso=2: elliptic medium
! ianiso=3: tilted elliptic medium
! ianiso=4: VTI medium
! ianiso=5: TTI medium
!
!            -----------
!           |  A  |  B  | | p |     | s_p |
!           | -  -  -   | |   |  =  |     |
!           |  C  |  D  | | q |     | s_q |
!            -----------
!
! b: buoyancy
! mu: complex bulk modulus
! delta: Thomsen parameter
! epsilon: Thomsen parameter
! teta: anisotropy tilt
! ianiso: type of anisotropy
! d1,d1b,d2,d2b: damping functions for PML
! mat,irn,icn: arrays for sparse impedance matrix storage
! n1e x n2e: grid dimensions
! h: grid spacing
! ne: number of non zero coefficients in the matrix
! omegac: complex frequency
!
! Copyright 2008-2011 SEISCOPE project, All rights reserved.
! 
! ---------------------------------------------------------------------------------------------------------


SUBROUTINE sub_evalmatrix_TTI(ve,b,mu,delta,epsilon,teta,ianiso,&
     d1,d1b,d2,d2b,d1q,d1qb,d2q,d2qb,ne, &
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
  real b(0:n1e+1,0:n2e+1),ve(n1e,n2e),delta(n1e,n2e),epsilon(n1e,n2e),teta(n1e,n2e)

  complex d1p,d1m,d1qp,d1qm,d2p,d2m,d2qp,d2qm
  complex omegac,omega2,xc0,m,ci
  complex d1(0:n1e+1),d1b(0:n1e+1),d2(0:n2e+1),d2b(0:n2e+1)
  complex d1q(0:n1e+1),d1qb(0:n1e+1),d2q(0:n2e+1),d2qb(0:n2e+1)
  complex mu(0:n1e+1,0:n2e+1)
  complex(kind=creal) mat(ne)
  integer :: nne,ianiso

  damp=0.

  if (ianiso.eq.1) then
     write(*,*) 'SIMULATION IN ISOTROPIC MEDIA'
  else if (ianiso.eq.2) then
     write(*,*) 'SIMULATION IN ELLIPTIC MEDIA'
  else if (ianiso.eq.3) then
     write(*,*) 'SIMULATION IN TILTED ELLIPTIC MEDIA'
  else if (ianiso.eq.4) then
     write(*,*) 'SIMULATION IN VTI MEDIA'
  else if (ianiso.eq.5) then
     write(*,*) 'SIMULATION IN TTI MEDIA'
  else
     write(*,*) 'WARNING: CASE NOT TAKEN INTO ACCOUNT'
  end if

  ! Build anisotropic parameters

  if (ianiso.eq.1) then
     delta(:,:)=0.
     epsilon(:,:)=0.
     teta(:,:)=0.
  end if

  if (ianiso.eq.2) then
     epsilon(:,:)=delta(:,:)
     teta(:,:)=0.
  end if

  if (ianiso.eq.3) then
     epsilon(:,:)=delta(:,:)
  end if

  if (ianiso.eq.4) then
     teta(:,:)=0.
  end if

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
  nne=n1e*n2e
  ci=cmplx(0.,1.)

  ! Initialize sparse impedance matrix

  irn(:)=0.
  icn(:)=0.
  mat(:)=cmplx(0.,0.)


  nnz=0
  write(*,*) 'OMEGA2 = ',omega2
  write(*,*) 'N1 N2 = ',n1e,n2e
  k=0


  ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ! LOOP OVER ROWS OF THE IMPEDANCE MATRIX
  ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  do i2=1,n2e
     do i1=1,n1e

        ax=1.+2.*delta(i1,i2)*cos(teta(i1,i2))*cos(teta(i1,i2))
        bx=-delta(i1,i2)*sin(2.*teta(i1,i2))
        cx=(1.+2.*delta(i1,i2))*cos(teta(i1,i2))*cos(teta(i1,i2))
        dx=-0.5*(1.+2.*delta(i1,i2))*sin(2.*teta(i1,i2))
        ex=2.*(epsilon(i1,i2)-delta(i1,i2))*cos(teta(i1,i2))*cos(teta(i1,i2))
        fx=-(epsilon(i1,i2)-delta(i1,i2))*sin(2.*teta(i1,i2))
        gx=2.*(epsilon(i1,i2)-delta(i1,i2))*cos(teta(i1,i2))*cos(teta(i1,i2))
        hx=-(epsilon(i1,i2)-delta(i1,i2))*sin(2.*teta(i1,i2))

        az=-delta(i1,i2)*sin(2.*teta(i1,i2))
        bz=1.+2.*delta(i1,i2)*sin(teta(i1,i2))*sin(teta(i1,i2))
        cz=-0.5*(1.+2.*delta(i1,i2))*sin(2.*teta(i1,i2))
        dz=(1.+2.*delta(i1,i2))*sin(teta(i1,i2))*sin(teta(i1,i2))
        ez=-(epsilon(i1,i2)-delta(i1,i2))*sin(2.*teta(i1,i2))
        fz=2.*(epsilon(i1,i2)-delta(i1,i2))*sin(teta(i1,i2))*sin(teta(i1,i2))
        gz=-(epsilon(i1,i2)-delta(i1,i2))*sin(2.*teta(i1,i2))
        hz=2.*(epsilon(i1,i2)-delta(i1,i2))*sin(teta(i1,i2))*sin(teta(i1,i2))

      !  write(*,*)'ax,bx,cx,dx,ex,fx,gx,hx',ax,bx,cx,dx,ex,fx,gx,hx
      !  write(*,*)'az,bz,cz,dz,ez,fz,gz,hz',az,bz,cz,dz,ez,fz,gz,hz


        call subbaverage(b,n1e,n2e,i1,i2,bpm,bmp,bpp,bmm,b00,bp0,bm0,b0p,b0m)

        d2p=d2b(i2)
        d2m=d2b(i2-1)
        d1p=d1b(i1)
        d1m=d1b(i1-1)
        d2qp=d2qb(i2)
        d2qm=d2qb(i2-1)
        d1qp=d1qb(i1)
        d1qm=d1qb(i1-1)


        ! ##################  BLOCK MATRIX A ########################################################"
        ! Block matrix A
        !
        ! Node 00 -----------------------------------

        k=(i2-1)*n1e+i1
        !
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
        !
        ! Node 01  ----------------------------------
        !

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
        !
        ! Node 0-1 --------------------------------------------------------------
        !

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
        !
        ! Node -10 ---------------------------------------------------------------
        !
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
        !
        ! Node 10
        !
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

        !
        ! Node -11 -----------------------------------------------------------------
        !
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
        !
        ! Node +1-1 ---------------------------------------------------------------
        !
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
        !
        ! Node +1+1 ---------------------------------------------------------------
        !
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

        !
        ! Node -1-1 ---------------------------------------------------------------
        !
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


        ! Block matrices B, C and D computed only for VTI/TTI media

        IF (ianiso.eq.4.or.ianiso.eq.5) THEN    !VTI OR TTI
           ! ##################  BLOCK MATRIX B ###################################

           ! Block matrix B
           ! Node 00 -----------------------------------

           !
           l=(i2-1)*n1e+i1+nne
           nnz=nnz+1
           irn(nnz)=k
           icn(nnz)=l
           mat(nnz)= &
                w1*( &
                -0.25*h2*d2(i2)*cx*(bmp*d2qp+bpm*d2qm+bpp*d2qp+bmm*d2qm)  &
                +0.25*h2*d2(i2)*dx*(bmp*d1qm+bpm*d1qp-bpp*d1qp-bmm*d1qm)  &
                +0.25*h2*d1(i1)*cz*(bmp*d2qp+bpm*d2qm-bpp*d2qp-bmm*d2qm)  &
                -0.25*h2*d1(i1)*dz*(bmp*d1qm+bpm*d1qp+bpp*d1qp+bmm*d1qm)) &
                +w2*(-d2(i2)*h2*cx*(d2qp*b0p+d2qm*b0m)-d1(i1)*h2*dz*(d1qp*bp0+d1qm*bm0) )

           !
           ! Node 01  ----------------------------------
           !

           l=0

           IF (i2<n2e) then

              l=(i2-1+1)*n1e+i1+nne

              nnz=nnz+1

              irn(nnz)=k
              icn(nnz)=l
              mat(nnz)= &
                   w1*( &
                   0.25*h2*d2(i2)*cx*(bmp*d2qp+bpp*d2qp)  &
                   +0.25*h2*d2(i2)*dx*(bmp*d1qm-bpp*d1qp) &
                   +0.25*h2*d1(i1)*cz*(-bmp*d2qp+bpp*d2qp) &
                   +0.25*h2*d1(i1)*dz*(-bmp*d1qm-bpp*d1qp) ) &
                   +w2*(d2(i2)*h2*cx*b0p*d2qp + 0.25*d1(i1)*cz*h2*(bp0*d2q(i2)-bm0*d2q(i2)) )

           end if
           !
           ! Node 0-1 --------------------------------------------------------------
           !

           IF (i2>1) THEN

              l=(i2-1-1)*n1e+i1+nne
              nnz=nnz+1

              irn(nnz)=k
              icn(nnz)=l
              mat(nnz)= &
                   w1*( &
                   0.25*h2*d2(i2)*cx*(bpm*d2qm+bmm*d2qm) &
                   +0.25*h2*d2(i2)*dx*(bpm*d1qp-bmm*d1qm) &
                   +0.25*h2*d1(i1)*cz*(-bpm*d2qm+bmm*d2qm) &
                   +0.25*h2*d1(i1)*dz*(-bpm*d1qp-bmm*d1qm)) &
                   +w2*(d2(i2)*h2*cx*b0m*d2qm + 0.25*d1(i1)*cz*h2*(-bp0*d2q(i2)+bm0*d2q(i2)) )



           end if
           !
           ! Node -10 ---------------------------------------------------------------
           !
           if (i1 > 1) then

              l=(i2-1)*n1e+i1-1+nne

              nnz=nnz+1

              irn(nnz)=k
              icn(nnz)=l
              mat(nnz)= &
                   w1*( &
                   0.25*h2*d2(i2)*cx*(-bmp*d2qp-bmm*d2qm) &
                   +0.25*h2*d2(i2)*dx*(-bmp*d1qm+bmm*d1qm) &
                   +0.25*h2*d1(i1)*cz*(bmp*d2qp-bmm*d2qm) &
                   +0.25*h2*d1(i1)*dz*(bmp*d1qm+bmm*d1qm)) &
                   +w2*(d1(i1)*h2*dz*bm0*d1qm + 0.25*d2(i2)*dx*h2*(-b0p*d1q(i1)+b0m*d1q(i1)) )



           endif
           !
           ! Node 10
           !
           if (i1 < n1e) then

              l=(i2-1)*n1e+i1+1+nne

              nnz=nnz+1

              irn(nnz)=k
              icn(nnz)=l
              mat(nnz)=   &
                   w1*( &
                   0.25*h2*d2(i2)*cx*(-bpm*d2qm-bpp*d2qp) &
                   +0.25*h2*d2(i2)*dx*(-bpm*d1qp+bpp*d1qp) &
                   +0.25*h2*d1(i1)*cz*(bpm*d2qm-bpp*d2qp) &
                   +0.25*h2*d1(i1)*dz*(bpm*d1qp+bpp*d1qp)) &
                   +w2*(d1(i1)*h2*dz*bp0*d1qp + 0.25*d2(i2)*dx*h2*(-b0m*d1q(i1)+b0p*d1q(i1)) )


           end if

           !
           ! Node -11 -----------------------------------------------------------------
           !
           if (i1 > 1 .AND. i2 < n2e) then

              l=(i2-1+1)*n1e+i1-1+nne

              nnz=nnz+1

              irn(nnz)=k
              icn(nnz)=l
              mat(nnz)= &
                   w1*( &
                   0.25*h2*d2(i2)*cx*bmp*d2qp &
                   -0.25*h2*d2(i2)*dx*bmp*d1qm &
                   -0.25*h2*d1(i1)*cz*bmp*d2qp &
                   +0.25*h2*d1(i1)*dz*bmp*d1qm) &
                   +w2*(-0.25*d2(i2)*h2*dx*b0p*d1q(i1) - 0.25*h2*d1(i1)*cz*d2q(i2)*bm0 )


           end if
           !
           ! Node +1-1 ---------------------------------------------------------------
           !
           if (i1 < n1e .AND. i2 >1) then

              l=(i2-1-1)*n1e+i1+1+nne

              nnz=nnz+1

              irn(nnz)=k
              icn(nnz)=l
              mat(nnz)= &
                   w1*( &
                   0.25*h2*d2(i2)*cx*bpm*d2qm &
                   -0.25*h2*d2(i2)*dx*bpm*d1qp &
                   -0.25*h2*d1(i1)*cz*bpm*d2qm &
                   +0.25*h2*d1(i1)*dz*bpm*d1qp ) &
                   +w2*(-0.25*d2(i2)*h2*dx*b0m*d1q(i1) - 0.25*h2*d1(i1)*cz*d2q(i2)*bp0 )

           end if
           !
           ! Node +1+1 ---------------------------------------------------------------
           !
           if (i1 < n1e .AND. i2 < n2e) then

              l=(i2-1+1)*n1e+i1+1+nne


              nnz=nnz+1

              irn(nnz)=k
              icn(nnz)=l
              mat(nnz)= &
                   w1*( &
                   0.25*h2*d2(i2)*cx*bpp*d2qp &
                   +0.25*h2*d2(i2)*dx*bpp*d1qp &
                   +0.25*h2*d1(i1)*cz*bpp*d2qp &
                   +0.25*h2*d1(i1)*dz*bpp*d1qp ) &
                   +w2*(0.25*d2(i2)*h2*dx*b0p*d1q(i1) + 0.25*h2*d1(i1)*cz*d2q(i2)*bp0 )


           end if

           !
           ! Node -1-1 ---------------------------------------------------------------
           !
           if (i1 > 1 .AND. i2 > 1) then

              l=(i2-1-1)*n1e+i1-1+nne

              nnz=nnz+1

              irn(nnz)=k
              icn(nnz)=l
              mat(nnz)= &
                   w1*( &
                   0.25*h2*d2(i2)*cx*bmm*d2qm &
                   +0.25*h2*d2(i2)*dx*bmm*d1qm &
                   +0.25*h2*d1(i1)*cz*bmm*d2qm &
                   +0.25*h2*d1(i1)*dz*bmm*d1qm ) &
                   +w2*(0.25*d2(i2)*h2*dx*b0m*d1q(i1) + 0.25*h2*d1(i1)*cz*d2q(i2)*bm0 )



           end if

           ! ##################  BLOCK MATRIX C ####################################
           !
           !
           ! Block matrix C
           !
           ! COEFFICIENTS E, F

           ! Node 00 -----------------------------------

           k=(i2-1)*n1e+i1+nne
           !
           l=(i2-1)*n1e+i1
           nnz=nnz+1
           irn(nnz)=k
           icn(nnz)=l
           mat(nnz)=  &
                w1*( &
                -0.25*h2*d2q(i2)*ex*(bmp*d2p+bpm*d2m+bpp*d2p+bmm*d2m)  &
                +0.25*h2*d2q(i2)*fx*(bmp*d1m+bpm*d1p-bpp*d1p-bmm*d1m)  &
                +0.25*h2*d1q(i1)*ez*(bmp*d2p+bpm*d2m-bpp*d2p-bmm*d2m)  &
                -0.25*h2*d1q(i1)*fz*(bmp*d1m+bpm*d1p+bpp*d1p+bmm*d1m)) &
                +w2*(-d2q(i2)*h2*ex*(d2p*b0p+d2m*b0m)-d1q(i1)*h2*fz*(d1p*bp0+d1m*bm0) )

           !
           ! Node 01  ----------------------------------
           !

           l=0

           IF (i2<n2e) then

              l=(i2-1+1)*n1e+i1

              nnz=nnz+1

              irn(nnz)=k
              icn(nnz)=l
              mat(nnz)= &
                   w1*( &
                   0.25*h2*d2q(i2)*ex*(bmp*d2p+bpp*d2p)  &
                   +0.25*h2*d2q(i2)*fx*(bmp*d1m-bpp*d1p) &
                   +0.25*h2*d1q(i1)*ez*(-bmp*d2p+bpp*d2p) &
                   +0.25*h2*d1q(i1)*fz*(-bmp*d1m-bpp*d1p) ) &
                   +w2*(d2q(i2)*h2*ex*b0p*d2p + 0.25*d1q(i1)*ez*h2*(bp0*d2(i2)-bm0*d2(i2)) )


           end if
           !
           ! Node 0-1 --------------------------------------------------------------
           !

           IF (i2>1) THEN

              l=(i2-1-1)*n1e+i1
              nnz=nnz+1

              irn(nnz)=k
              icn(nnz)=l
              mat(nnz)= &
                   w1*( &
                   0.25*h2*d2q(i2)*ex*(bpm*d2m+bmm*d2m) &
                   +0.25*h2*d2q(i2)*fx*(bpm*d1p-bmm*d1m) &
                   +0.25*h2*d1q(i1)*ez*(-bpm*d2m+bmm*d2m) &
                   +0.25*h2*d1q(i1)*fz*(-bpm*d1p-bmm*d1m)) &
                   +w2*(d2q(i2)*h2*ex*b0m*d2m + 0.25*d1q(i1)*ez*h2*(-bp0*d2(i2)+bm0*d2(i2)) )


           end if
           !
           ! Node -10 ---------------------------------------------------------------
           !
           if (i1 > 1) then

              l=(i2-1)*n1e+i1-1

              nnz=nnz+1

              irn(nnz)=k
              icn(nnz)=l
              mat(nnz)= &
                   w1*( &
                   0.25*h2*d2q(i2)*ex*(-bmp*d2p-bmm*d2m) &
                   +0.25*h2*d2q(i2)*fx*(-bmp*d1m+bmm*d1m) &
                   +0.25*h2*d1q(i1)*ez*(bmp*d2p-bmm*d2m) &
                   +0.25*h2*d1q(i1)*fz*(bmp*d1m+bmm*d1m)) &
                   +w2*(d1q(i1)*h2*fz*bm0*d1m + 0.25*d2q(i2)*fx*h2*(-b0p*d1(i1)+b0m*d1(i1)) )


           end if
           !
           ! Node 10
           !
           if (i1 < n1e) then

              l=(i2-1)*n1e+i1+1

              nnz=nnz+1

              irn(nnz)=k
              icn(nnz)=l
              mat(nnz)=   &
                   w1*( &
                   0.25*h2*d2q(i2)*ex*(-bpm*d2m-bpp*d2p) &
                   +0.25*h2*d2q(i2)*fx*(-bpm*d1p+bpp*d1p) &
                   +0.25*h2*d1q(i1)*ez*(bpm*d2m-bpp*d2p) &
                   +0.25*h2*d1q(i1)*fz*(bpm*d1p+bpp*d1p)) &
                   +w2*(d1q(i1)*h2*fz*bp0*d1p + 0.25*d2q(i2)*fx*h2*(-b0m*d1(i1)+b0p*d1(i1)) )


           end if

           !
           ! Node -11 -----------------------------------------------------------------
           !
           if (i1 > 1 .AND. i2 < n2e) then

              l=(i2-1+1)*n1e+i1-1

              nnz=nnz+1

              irn(nnz)=k
              icn(nnz)=l
              mat(nnz)= &
                   w1*( &
                   0.25*h2*d2q(i2)*ex*bmp*d2p &
                   -0.25*h2*d2q(i2)*fx*bmp*d1m &
                   -0.25*h2*d1q(i1)*ez*bmp*d2p &
                   +0.25*h2*d1q(i1)*fz*bmp*d1m) &
                   +w2*(-0.25*d2q(i2)*h2*fx*b0p*d1(i1) - 0.25*h2*d1q(i1)*ez*d2(i2)*bm0 )


           end if
           !
           ! Node +1-1 ---------------------------------------------------------------
           !
           if (i1 < n1e .AND. i2 >1) then

              l=(i2-1-1)*n1e+i1+1

              nnz=nnz+1

              irn(nnz)=k
              icn(nnz)=l
              mat(nnz)= &
                   w1*( &
                   0.25*h2*d2q(i2)*ex*bpm*d2m &
                   -0.25*h2*d2q(i2)*fx*bpm*d1p &
                   -0.25*h2*d1q(i1)*ez*bpm*d2m &
                   +0.25*h2*d1q(i1)*fz*bpm*d1p ) &
                   +w2*(-0.25*d2q(i2)*h2*fx*b0m*d1(i1) - 0.25*h2*d1q(i1)*ez*d2(i2)*bp0 )


           end if
           !
           ! Node +1+1 ---------------------------------------------------------------
           !
           if (i1 < n1e .AND. i2 < n2e) then

              l=(i2-1+1)*n1e+i1+1


              nnz=nnz+1

              irn(nnz)=k
              icn(nnz)=l
              mat(nnz)= &
                   w1*( &
                   0.25*h2*d2q(i2)*ex*bpp*d2p &
                   +0.25*h2*d2q(i2)*fx*bpp*d1p &
                   +0.25*h2*d1q(i1)*ez*bpp*d2p &
                   +0.25*h2*d1q(i1)*fz*bpp*d1p ) &
                   +w2*(0.25*d2q(i2)*h2*fx*b0p*d1(i1) + 0.25*h2*d1q(i1)*ez*d2(i2)*bp0 )


           end if

           !
           ! Node -1-1 ---------------------------------------------------------------
           !
           if (i1 > 1 .AND. i2 > 1) then

              l=(i2-1-1)*n1e+i1-1

              nnz=nnz+1

              irn(nnz)=k
              icn(nnz)=l
              mat(nnz)= &
                   w1*( &
                   0.25*h2*d2q(i2)*ex*bmm*d2m &
                   +0.25*h2*d2q(i2)*fx*bmm*d1m &
                   +0.25*h2*d1q(i1)*ez*bmm*d2m &
                   +0.25*h2*d1q(i1)*fz*bmm*d1m ) &
                   +w2*(0.25*d2q(i2)*h2*fx*b0m*d1(i1) + 0.25*h2*d1q(i1)*ez*d2(i2)*bm0 )


           end if

           ! ##################  BLOCK MATRIX D ####################################

           ! Block matrix D
           !
           ! COEFFICIENTS G, H

           ! Node 00 -----------------------------------

           k=(i2-1)*n1e+i1+nne
           !
           l=(i2-1)*n1e+i1+nne
           nnz=nnz+1
           irn(nnz)=k
           icn(nnz)=l
           mat(nnz)= wm1*omega2/mu(i1,i2) &
                +w1*( &
                -0.25*h2*d2q(i2)*gx*(bmp*d2qp+bpm*d2qm+bpp*d2qp+bmm*d2qm)  &
                +0.25*h2*d2q(i2)*hx*(bmp*d1qm+bpm*d1qp-bpp*d1qp-bmm*d1qm)  &
                +0.25*h2*d1q(i1)*gz*(bmp*d2qp+bpm*d2qm-bpp*d2qp-bmm*d2qm)  &
                -0.25*h2*d1q(i1)*hz*(bmp*d1qm+bpm*d1qp+bpp*d1qp+bmm*d1qm)) &
                +w2*(-d2q(i2)*h2*gx*(d2qp*b0p+d2qm*b0m)-d1q(i1)*h2*hz*(d1qp*bp0+d1qm*bm0) )
           mat(nnz)=mat(nnz)+4*damp*mat(nnz)


           !
           ! Node 01  ----------------------------------
           !

           l=0

           IF (i2<n2e) then

              l=(i2-1+1)*n1e+i1+nne

              nnz=nnz+1

              irn(nnz)=k
              icn(nnz)=l
              mat(nnz)= &
                   wm2*omega2/mu(i1,i2+1) &
                   +w1*( &
                   0.25*h2*d2q(i2)*gx*(bmp*d2qp+bpp*d2qp)  &
                   +0.25*h2*d2q(i2)*hx*(bmp*d1qm-bpp*d1qp) &
                   +0.25*h2*d1q(i1)*gz*(-bmp*d2qp+bpp*d2qp) &
                   +0.25*h2*d1q(i1)*hz*(-bmp*d1qm-bpp*d1qp) ) &
                   +w2*(d2q(i2)*h2*gx*b0p*d2qp + 0.25*d1q(i1)*gz*h2*(bp0*d2q(i2)-bm0*d2q(i2)) )
              mat(nnz)=mat(nnz)-damp*mat(nnz)

           end if
           !
           ! Node 0-1 --------------------------------------------------------------
           !

           IF (i2>1) THEN

              l=(i2-1-1)*n1e+i1+nne
              nnz=nnz+1

              irn(nnz)=k
              icn(nnz)=l
              mat(nnz)= &
                   wm2*omega2/mu(i1,i2-1) &
                   +w1*( &
                   0.25*h2*d2q(i2)*gx*(bpm*d2qm+bmm*d2qm) &
                   +0.25*h2*d2q(i2)*hx*(bpm*d1qp-bmm*d1qm) &
                   +0.25*h2*d1q(i1)*gz*(-bpm*d2qm+bmm*d2qm) &
                   +0.25*h2*d1q(i1)*hz*(-bpm*d1qp-bmm*d1qm)) &
                   +w2*(d2q(i2)*h2*gx*b0m*d2qm + 0.25*d1q(i1)*gz*h2*(-bp0*d2q(i2)+bm0*d2q(i2)) )
              mat(nnz)=mat(nnz)-damp*mat(nnz)


           end if
           !
           ! Node -10 ---------------------------------------------------------------
           !
           if (i1 > 1) then

              l=(i2-1)*n1e+i1-1+nne

              nnz=nnz+1

              irn(nnz)=k
              icn(nnz)=l
              mat(nnz)= &
                   wm2*omega2/mu(i1-1,i2) &
                   +w1*( &
                   0.25*h2*d2q(i2)*gx*(-bmp*d2qp-bmm*d2qm) &
                   +0.25*h2*d2q(i2)*hx*(-bmp*d1qm+bmm*d1qm) &
                   +0.25*h2*d1q(i1)*gz*(bmp*d2qp-bmm*d2qm) &
                   +0.25*h2*d1q(i1)*hz*(bmp*d1qm+bmm*d1qm)) &
                   +w2*(d1q(i1)*h2*hz*bm0*d1qm + 0.25*d2q(i2)*hx*h2*(-b0p*d1q(i1)+b0m*d1q(i1)) )
              mat(nnz)=mat(nnz)-damp*mat(nnz)


           end if
           !
           ! Node 10
           !
           if (i1 < n1e) then

              l=(i2-1)*n1e+i1+1+nne

              nnz=nnz+1

              irn(nnz)=k
              icn(nnz)=l
              mat(nnz)=   &
                   wm2*omega2/mu(i1+1,i2) &
                   +w1*( &
                   0.25*h2*d2q(i2)*gx*(-bpm*d2qm-bpp*d2qp) &
                   +0.25*h2*d2q(i2)*hx*(-bpm*d1qp+bpp*d1qp) &
                   +0.25*h2*d1q(i1)*gz*(bpm*d2qm-bpp*d2qp) &
                   +0.25*h2*d1q(i1)*hz*(bpm*d1qp+bpp*d1qp)) &
                   +w2*(d1q(i1)*h2*hz*bp0*d1qp + 0.25*d2q(i2)*hx*h2*(-b0m*d1q(i1)+b0p*d1q(i1)) )
              mat(nnz)=mat(nnz)-damp*mat(nnz)



           end if

           !
           ! Node -11 -----------------------------------------------------------------
           !
           if (i1 > 1 .AND. i2 < n2e) then

              l=(i2-1+1)*n1e+i1-1+nne

              nnz=nnz+1

              irn(nnz)=k
              icn(nnz)=l
              mat(nnz)= &
                   wm3*omega2/mu(i1-1,i2+1) &
                   +w1*( &
                   0.25*h2*d2q(i2)*gx*bmp*d2qp &
                   -0.25*h2*d2q(i2)*hx*bmp*d1qm &
                   -0.25*h2*d1q(i1)*gz*bmp*d2qp &
                   +0.25*h2*d1q(i1)*hz*bmp*d1qm) &
                   +w2*(-0.25*d2q(i2)*h2*hx*b0p*d1q(i1) - 0.25*h2*d1q(i1)*gz*d2q(i2)*bm0 )


           end if
           !
           ! Node +1-1 ---------------------------------------------------------------
           !
           if (i1 < n1e .AND. i2 >1) then

              l=(i2-1-1)*n1e+i1+1+nne

              nnz=nnz+1

              irn(nnz)=k
              icn(nnz)=l
              mat(nnz)= &
                   wm3*omega2/mu(i1+1,i2-1) &
                   +w1*( &
                   0.25*h2*d2q(i2)*gx*bpm*d2qm &
                   -0.25*h2*d2q(i2)*hx*bpm*d1qp &
                   -0.25*h2*d1q(i1)*gz*bpm*d2qm &
                   +0.25*h2*d1q(i1)*hz*bpm*d1qp ) &
                   +w2*(-0.25*d2q(i2)*h2*hx*b0m*d1q(i1) - 0.25*h2*d1q(i1)*gz*d2q(i2)*bp0 )




           end if
           !
           ! Node +1+1 ---------------------------------------------------------------
           !
           if (i1 < n1e .AND. i2 < n2e) then

              l=(i2-1+1)*n1e+i1+1+nne


              nnz=nnz+1

              irn(nnz)=k
              icn(nnz)=l
              mat(nnz)= &
                   wm3*omega2/mu(i1+1,i2+1) &
                   +w1*( &
                   0.25*h2*d2q(i2)*gx*bpp*d2qp &
                   +0.25*h2*d2q(i2)*hx*bpp*d1qp &
                   +0.25*h2*d1q(i1)*gz*bpp*d2qp &
                   +0.25*h2*d1q(i1)*hz*bpp*d1qp ) &
                   +w2*(0.25*d2q(i2)*h2*hx*b0p*d1q(i1) + 0.25*h2*d1q(i1)*gz*d2q(i2)*bp0 )


           end if

           !
           ! Node -1-1 ---------------------------------------------------------------
           !
           if (i1 > 1 .AND. i2 > 1) then

              l=(i2-1-1)*n1e+i1-1+nne

              nnz=nnz+1

              irn(nnz)=k
              icn(nnz)=l
              mat(nnz)= &
                   wm3*omega2/mu(i1-1,i2-1) &
                   +w1*( &
                   0.25*h2*d2q(i2)*gx*bmm*d2qm &
                   +0.25*h2*d2q(i2)*hx*bmm*d1qm &
                   +0.25*h2*d1q(i1)*gz*bmm*d2qm &
                   +0.25*h2*d1q(i1)*hz*bmm*d1qm ) &
                   +w2*(0.25*d2q(i2)*h2*hx*b0m*d1q(i1) + 0.25*h2*d1q(i1)*gz*d2q(i2)*bm0 )


           end if

        END IF

     end do            !end do i1=1,n1
  end do            !end do i2=1,n2

  write(*,*) 'NNZ NE = ',nnz,ne
  write(*,*)'MAXVAL(irn)',MAXVAL(IRN)
  write(*,*)'MAXVAL(JCN)',MAXVAL(ICN)

end SUBROUTINE sub_evalmatrix_TTI
