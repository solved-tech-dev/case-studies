!=======================================================================
! Copyright 2008-2011 SEISCOPE project, All rights reserved.
!Copyright 2013-2021 SEISCOPEII project, All rights reserved.
!=======================================================================

! Modification from an original subroutine developed by V. Etienne

! Reference: Arbitrary source and receiver positioning in finite-difference schemes
! using Kaiser windowed sinc functions, G. Hicks, Geophysics, 67(1), p. 156-166, 2002.
! Inputs:
! ndim: dimension of the problem (1: 1D, 2: 2D, 3: 3D)
! ifreesurf [0/1]: mirror coefficients when free surface boundary conditions must be taken into account (1)
! irad: radius in number of grid points of the windowed sinc function
! nalpha: number of tabulated values
! Outputs:
! bhicks_monopole,  bhicks_dipole: optimal Hicks factors for Kaiser functions (equation 6)
! irad1, irad2, irad3: radius in number of grid points of the windowed sinc function along the 3 directions.
! For 2D problems, irad3 = 0, for 1D problems, irad2 and irad3 = 0.
! nalpha1, nalpha2, nalpha3: number of tabulated values along 3 directions.
! For 2D problems, nalpha3 = 1, for 1D problems, nalpha2 and nalpha3 = 1.
! ishift1 (0/1), ishift2 (0/1), ishift3 (0/1): index equal to 0 or 1 to manage 1D, 2D, 3D casis
! coef_hicks_monopole(-irad1+ishift1:irad1, -irad2+ishift2:irad2, -irad3+ishift3:irad3, nalpha1, nalpha2, nalpha3, 0:irad1-1))
! coef_hicks_dipole(-irad1+ishift1:irad1, -irad2+ishift2:irad2, -irad3+ishift3:irad3, nalpha1, nalpha2, nalpha3, 0:irad1-1))
! weighting coefficients of the interpolation function for different values of alpha (see Hicks paper) and different positioning
! with respect to the free surface
!
! Version 6 March 2013 (S. Operto)
!
! Modified 29 September 2013 (F. Lavoue)
! to distinguish between horizontal and vertical dipoles
! + correct dipole interpolation (05/10/2013).
! ------------------------------------------------------------------------------------------------------------------------------------

SUBROUTINE sub_Hicks(pbdir,ndim, ifreesurf, irad, nalpha,&
     bhicks_monopole, bhicks_dipole, &
     irad1, irad2, irad3, &
     nalpha1, nalpha2, nalpha3, &
     ishift1, ishift2, ishift3)

  IMPLICIT NONE
#include "common.h"
# include "pbdirect.h"

  type(pbdirect):: pbdir

  INTEGER :: irad,irad1,irad2,irad3
  INTEGER :: nalpha,nalpha1,nalpha2,nalpha3
  INTEGER :: i_alpha,i1_alpha,i2_alpha,i3_alpha,i1_sym
  INTEGER :: iloop,iloop1,iloop2,iloop3
  INTEGER :: ndim,ifreesurf,ip,debug,ishift1, ishift2, ishift3 

  REAL :: bhicks_monopole, bhicks_dipole
  REAL :: i0d_monopole,i0d_dipole,i0n_monopole,i0n_dipole,kaiser_monopole,kaiser_dipole
  REAL :: alpha_courant, BESSI0, HICKS_EPS_TEST, HICKS_EPS_KAISER

  REAL, ALLOCATABLE :: coef_hicks1D_monopole(:,:),coef_hicks_temp_monopole(:,:,:,:,:,:)
  REAL, ALLOCATABLE :: coef_hicks1D_dipole(:,:),coef_hicks_temp_dipole_x(:,:,:,:,:,:),coef_hicks_temp_dipole_z(:,:,:,:,:,:)
  REAL, ALLOCATABLE :: alpha(:)

  HICKS_EPS_TEST = 1e-15
  HICKS_EPS_KAISER = 1e-15
  debug = 0

  ! compute Bessel func. for bb
  i0d_monopole = BESSI0(bhicks_monopole)   ! denominator, equation (6)
  i0d_dipole   = BESSI0(bhicks_dipole)

  ! Adaptation as a function of the dimension of propagation
  IF ( ndim == 1 ) THEN
     irad1 = irad
     irad2 = 0
     irad3 = 0
     ishift1=1
     ishift2=0
     ishift3=0
     nalpha1 = nalpha
     nalpha2 = 1
     nalpha3 = 1

  ELSE IF ( ndim == 2 ) THEN
     irad1 = irad
     irad2 = irad
     irad3 = 0
     ishift1=1
     ishift2=1
     ishift3=0
     nalpha1 = nalpha
     nalpha2 = nalpha
     nalpha3 = 1

  ELSE   !3D
     irad1 = irad
     irad2 = irad
     irad3 = irad
     ishift1=1
     ishift2=1
     ishift3=1
     nalpha1 = nalpha
     nalpha2 = nalpha
     nalpha3 = nalpha
  END IF

  ALLOCATE(coef_hicks1D_monopole(-irad+1:irad, nalpha)) 
  ALLOCATE(coef_hicks1D_dipole(-irad+1:irad, nalpha)) 
  coef_hicks1D_monopole(:,:) = 0.
  coef_hicks1D_dipole(:,:) = 0.

  ALLOCATE(coef_hicks_temp_monopole(-irad1+ishift1:irad1, -irad2+ishift2:irad2, -irad3+ishift3:irad3, nalpha1, nalpha2, nalpha3))
  ALLOCATE(coef_hicks_temp_dipole_x(-irad1+ishift1:irad1, -irad2+ishift2:irad2, -irad3+ishift3:irad3, nalpha1, nalpha2, nalpha3))
  ALLOCATE(coef_hicks_temp_dipole_z(-irad1+ishift1:irad1, -irad2+ishift2:irad2, -irad3+ishift3:irad3, nalpha1, nalpha2, nalpha3))
  coef_hicks_temp_monopole(:,:,:,:,:,:) = 0.
  coef_hicks_temp_dipole_x(:,:,:,:,:,:) = 0.
  coef_hicks_temp_dipole_z(:,:,:,:,:,:) = 0.

  ! Tabulation of Hicks coefficients: set the interval positions between two grid nodes ( 0 <= alpha < 1 )
  ! For alpha = 0, the source is on the grid node.

  ALLOCATE(alpha(nalpha))

  DO i_alpha = 1, nalpha
     alpha(i_alpha) = REAL(i_alpha - 1) / REAL(nalpha)
     IF ( debug .GT. 0 ) WRITE(*,*) "ALPHA = ",i_alpha,alpha(i_alpha)
  ENDDO

  ! Convention to number the Hicks coefficients coef_hicks(-irad+1:irad)

  ! -------------------------------------------------------------------------------------------> x,y,z
  !     x	x	x	x	x	x	x	x	x	x	x	x 	(grid points)
  !                                     * 								(Position of source)
  !            -3      -2      -1       0       1	2	3	4				(indexing in array(-irad+1;irad))
  !

  ! Outer loop over the tabulated positions in the interval [0 ; 1]
  DO i_alpha = 1, nalpha 

     !==================================================================================================================================================
     !FL, 05/10/2013: fill the array for ALL points in interpolation region and alpha values,
     !                then specify values for i_irad = 0 and alpha = 0 (case of src on grid node)
     !    With the previous "if/else" statement, coef_hicks1D_dipole(iloop,i_alpha) remained 0 for alpha=0 and iloop=/=0, whereas it should be non-zero.
     !==================================================================================================================================================
     !
     !! case of the position at the grid point location (commented 05/10/2013, FL: now after inner loop on interpolation region)
     !IF (alpha(i_alpha) <  HICKS_EPS_TEST) THEN      ! for i_alpha = 1 and i_irad = 0, the source is on the grid point of index 0 (see above; 4th value in the array)
     !   coef_hicks1D_monopole(0, i_alpha) = 1.   !coef_hicks_monopole(0,1) = 1.
     !   coef_hicks1D_dipole(0, i_alpha) = 0.     !coef_hicks_dipole(0,1) = 0.
     !
     !ELSE   ! coef_hicks1D(0, 2:nalpha) we fill the array for alpha > 0

     ! Inner loop over the spatial support of the interpolating function
     ! summation is between -irad and irad-1 to avoid negative value in the square root 
     ! in the argument of the Bessel function  (alpha_courant never bigger than irad * h)
     DO iloop = -irad+1, irad      

        ! Discrete position of the potential source x = n + alpha (see equation (3) in Hicks, Geophysics, 2002).
        ! Note that alpha_courant always smaller than irad because alpha >= 0
        alpha_courant = REAL(iloop)  - alpha(i_alpha)

        ! compute coef.
        IF (ABS(alpha_courant) > HICKS_EPS_TEST) THEN

           i0n_monopole = BESSI0(bhicks_monopole* SQRT(1. - alpha_courant*alpha_courant / (irad*irad)))     !Equation (6), numerator
           i0n_dipole   = BESSI0(bhicks_dipole  * SQRT(1. - alpha_courant*alpha_courant / (irad*irad)))                      
           kaiser_monopole = i0n_monopole / i0d_monopole                                                 !Equation 6
           kaiser_dipole   = i0n_dipole   / i0d_dipole
           ! index iloop+1 such that the support of the interpolation below the free surface is greater than irad*h
           ! this requires exploitation of the symmetry of the Sinc function for some kind of demultiplexing... (trick of VE)

           coef_hicks1D_monopole(iloop,i_alpha ) = kaiser_monopole * sin(pi*alpha_courant) / (pi*alpha_courant)      !Equations (4), (5)

           coef_hicks1D_dipole( iloop, i_alpha ) = kaiser_dipole * ( cos(pi*alpha_courant) - ( sin(pi*alpha_courant) / (pi*alpha_courant) ) ) &   ! Equation (9)
                / alpha_courant   !and not pi*alpha, corrected 05/10/2013, FL.

           ! check small values
           !IF (ABS(coef_hicks1D_monopole(iloop , nalpha - i_alpha + 2)) <  HICKS_EPS_KAISER ) coef_hicks1D_monopole(iloop , nalpha - i_alpha + 2) = 0.
           !IF (ABS(coef_hicks1D_dipole(iloop , nalpha - i_alpha + 2)) <  HICKS_EPS_KAISER ) coef_hicks1D_dipole(iloop , nalpha - i_alpha + 2) = 0.
           IF (ABS(coef_hicks1D_monopole(iloop , i_alpha)) <  HICKS_EPS_KAISER ) coef_hicks1D_monopole(iloop , i_alpha) = 0.
           IF (ABS(coef_hicks1D_dipole(iloop , i_alpha)) <  HICKS_EPS_KAISER ) coef_hicks1D_dipole(iloop , i_alpha) = 0.

        ENDIF

     ENDDO   !irad

     ! case of the position at the grid point location (FL, 05/10/2013: now after loop on interpolation region to force values for i_irad=0, alpha=0)
     IF (alpha(i_alpha) <  HICKS_EPS_TEST) THEN      ! for i_alpha = 1 and i_irad = 0, the source is on the grid point of index 0 (see above; 4th value in the array)
        coef_hicks1D_monopole(0, i_alpha) = 1.   !coef_hicks_monopole(0,1) = 1.
        coef_hicks1D_dipole(0, i_alpha) = 0.     !coef_hicks_dipole(0,1) = 0.
     END IF

  ENDDO   !i_alpha

  IF ( debug .GT. 0 ) THEN
     WRITE(*,*) 'alpha', alpha
     DO i_alpha = 1, nalpha 
        WRITE(*,*) 'i_alpha, coef_hicks1D_monopole(:, i_alpha)', i_alpha, coef_hicks1D_monopole(:, i_alpha)
        WRITE(*,*) 'i_alpha, coef_hicks1D_dipole(:, i_alpha)', i_alpha, coef_hicks1D_dipole(:, i_alpha)
     ENDDO
  ENDIF


  ! compute the coef. in 3D without interaction of free surface
  !==============================================================

  ! Tensorial construction of the 3D interpolation function
  !   (eq. 12 in Hicks, 2002)
  ! For 1D and 2D cases, note that iloop3 and/or iloop2 range between 0 and 0.

  ! loop on the spatial samples
  DO i3_alpha = 1, nalpha3      
     DO i2_alpha = 1, nalpha2  
        DO i1_alpha = 1, nalpha1 

           ! loop on the support of the interpolation
           DO iloop3 = -irad3 + ishift3,irad3
              DO iloop2 = -irad2 + ishift2,irad2
                 DO iloop1 = -irad1 + ishift1,irad1
                    !monopole
                    coef_hicks_temp_monopole(iloop1,iloop2,iloop3,i1_alpha,i2_alpha,i3_alpha) = coef_hicks1D_monopole(iloop3,i3_alpha) &
                         *coef_hicks1D_monopole(iloop2,i2_alpha) &
                         *coef_hicks1D_monopole(iloop1,i1_alpha)
                    !horizontal dipole
                    coef_hicks_temp_dipole_x(iloop1,iloop2,iloop3,i1_alpha,i2_alpha,i3_alpha) = coef_hicks1D_monopole(iloop3,i3_alpha) &
                         *coef_hicks1D_dipole(iloop2,i2_alpha)   &
                         *coef_hicks1D_monopole(iloop1,i1_alpha)
                    !vertical dipole
                    coef_hicks_temp_dipole_z(iloop1,iloop2,iloop3,i1_alpha,i2_alpha,i3_alpha) = coef_hicks1D_monopole(iloop3,i3_alpha) &
                         *coef_hicks1D_monopole(iloop2,i2_alpha) &
                         *coef_hicks1D_dipole(iloop1,i1_alpha)
                 ENDDO
              ENDDO
           ENDDO   !irad

        ENDDO
     ENDDO
  ENDDO   !i_alpha

  IF ( debug .GT. 0 ) THEN
     WRITE(*,*) 'coef_hicks_temp_monopole(:,:,:,1,1,1)', coef_hicks_temp_monopole(:,:,:,1,1,1)
     WRITE(*,*) 'coef_hicks_temp_dipole_x(:,:,:,1,1,1)', coef_hicks_temp_dipole_x(:,:,:,1,1,1)
     WRITE(*,*) 'coef_hicks_temp_dipole_z(:,:,:,1,1,1)', coef_hicks_temp_dipole_z(:,:,:,1,1,1)
  ENDIF

  ! Dynamic allocation for the output array coef_hicks; the seven^th field is used to manage free surface interaction
  ! at most, there are irad - 1 non zero coefficients above free surface (the coefficient on the free surface is 0)

  ! test if free surface

  ALLOCATE(pbdir%coef_hicks_monopole(-irad1+ishift1:irad1, -irad2+ishift2:irad2, -irad3+ishift3:irad3, nalpha1, nalpha2, nalpha3, 0:irad1-1))
  ALLOCATE(pbdir%coef_hicks_dipole_x(-irad1+ishift1:irad1, -irad2+ishift2:irad2, -irad3+ishift3:irad3, nalpha1, nalpha2, nalpha3, 0:irad1-1))
  ALLOCATE(pbdir%coef_hicks_dipole_z(-irad1+ishift1:irad1, -irad2+ishift2:irad2, -irad3+ishift3:irad3, nalpha1, nalpha2, nalpha3, 0:irad1-1))
  pbdir%coef_hicks_monopole(:,:,:,:,:,:,:) = 0.
  pbdir%coef_hicks_dipole_x(:,:,:,:,:,:,:) = 0.
  pbdir%coef_hicks_dipole_z(:,:,:,:,:,:,:) = 0.

  ! No interaction with free surface; field 7 is 0
  !==========================================
  pbdir%coef_hicks_monopole(:,:,:,:,:,:, 0) = coef_hicks_temp_monopole(:,:,:,:,:,:)
  pbdir%coef_hicks_dipole_x(:,:,:,:,:,:, 0) = coef_hicks_temp_dipole_x(:,:,:,:,:,:)
  pbdir%coef_hicks_dipole_z(:,:,:,:,:,:, 0) = coef_hicks_temp_dipole_z(:,:,:,:,:,:)


  ! Free-surface interaction; the field 7 is ip index, the number of points above the free surface
  !=======================================

  IF ( ifreesurf == 1 ) THEN

     ! loop over the possible number of nodes above free surface (nb interpolated points above free surface) 
     ! ip: number of nodes above free surface
     ! the "local coordinate system" for numbering indexes is [ - irad + 1; irad ]
     DO ip = 1, irad1 - 1

        ! Initialization
        pbdir%coef_hicks_monopole(:,:,:,:,:,:, ip) = coef_hicks_temp_monopole(:,:,:,:,:,:)
        pbdir%coef_hicks_dipole_x(:,:,:,:,:,:, ip) = coef_hicks_temp_dipole_x(:,:,:,:,:,:)
        pbdir%coef_hicks_dipole_z(:,:,:,:,:,:, ip) = coef_hicks_temp_dipole_z(:,:,:,:,:,:)

        ! Index of the node located at the free surface 
        i1_sym = -irad1 + ip + 1

        ! apply change of coef. below free surface
        DO iloop1 = 1, ip

           ! All the coefficients at and above the free surface at set to zeroes
           pbdir%coef_hicks_monopole(i1_sym - iloop1,:,:,:,:,:,ip) = 0.
           pbdir%coef_hicks_dipole_x(i1_sym - iloop1,:,:,:,:,:,ip) = 0.
           pbdir%coef_hicks_dipole_z(i1_sym - iloop1,:,:,:,:,:,ip) = 0.

           ! Mirroring of coefficients below free surface (anti-symmetry)
           pbdir%coef_hicks_monopole(i1_sym+iloop1,:,:,:,:,:,ip) =&
                coef_hicks_temp_monopole(i1_sym+iloop1,:,:,:,:,:) - &
                coef_hicks_temp_monopole(i1_sym-iloop1,:,:,:,:,:)

           pbdir%coef_hicks_dipole_x(i1_sym+iloop1,:,:,:,:,:,ip) =&
                coef_hicks_temp_dipole_x(i1_sym+iloop1,:,:,:,:,:) -&
                coef_hicks_temp_dipole_x(i1_sym-iloop1,:,:,:,:,:)

           pbdir%coef_hicks_dipole_z(i1_sym+iloop1,:,:,:,:,:,ip) =&
                coef_hicks_temp_dipole_z(i1_sym+iloop1,:,:,:,:,:) -&
                coef_hicks_temp_dipole_z(i1_sym-iloop1,:,:,:,:,:)
        ENDDO

        IF ( debug .GT. 0 ) THEN
           WRITE(*,*) 'ip, coef_hicks_monopole(:,:,:,1,1,1, ip)', ip, pbdir%coef_hicks_monopole(:,:,:,1,1,1, ip)
           WRITE(*,*) 'ip, coef_hicks_dipole_x(:,:,:,1,1,1, ip)', ip, pbdir%coef_hicks_dipole_x(:,:,:,1,1,1, ip)
           WRITE(*,*) 'ip, coef_hicks_dipole_z(:,:,:,1,1,1, ip)', ip, pbdir%coef_hicks_dipole_z(:,:,:,1,1,1, ip)
        ENDIF

     ENDDO

  ENDIF ! IF ( ifreesurf == 1) THEN


  ! deallocate local tables
  DEALLOCATE(coef_hicks1D_monopole)
  DEALLOCATE(coef_hicks_temp_monopole)
  DEALLOCATE(coef_hicks1D_dipole)
  DEALLOCATE(coef_hicks_temp_dipole_x)
  DEALLOCATE(coef_hicks_temp_dipole_z)
  DEALLOCATE(alpha)

END SUBROUTINE sub_Hicks



! ---------------------------------------------------------------------------------------------------------------
! FUNCTION BESSI0 (FROM NUMERICAL RECIPES)
! ---------------------------------------------------------------------------------------------------------------

FUNCTION BESSI0(X)

  !--------------------------------------- 
  ! Compute BESSEL function
  !---------------------------------------

  REAL BESSI0, X, AX
  REAL*8 Y,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9
  DATA P1,P2,P3,P4,P5,P6,P7/1.0D0,3.5156229D0,3.0899424D0,1.2067492D0,0.2659732D0,0.360768D-1,0.45813D-2/
  DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,0.1328592D-1,0.225319D-2,-0.157565D-2,0.916281D-2,-0.2057706D-1,0.2635537D-1,-0.1647633D-1,0.392377D-2/

  IF (ABS(X).LT.3.75) THEN
     Y=(X/3.75)**2
     BESSI0=P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
  ELSE
     AX=ABS(X)
     Y=3.75/AX
     BESSI0=(EXP(AX)/SQRT(AX))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9))))))))
  ENDIF

  RETURN 

END FUNCTION BESSI0


