!Copyright 2008-2011 SEISCOPE project, All rights reserved.
!Copyright 2013-2021 SEISCOPEII project, All rights reserved.
! SUBROUTINE SUBSAMPLEHICKS
! ix1, ix2, ix3: index of the shot / receiver position at the nearest point by default
! ialpha_1, ialpha_2, ialpha_3: indices in the table of dimension nalpha (pick the good tabulated value)
! nalpha_frees: number of grid points above free surface.
!
SUBROUTINE sub_sample_Hicks(x1,x2,x3,ix1,ix2,ix3,ialpha_1,ialpha_2,ialpha_3,nalpha_frees,&
     h,irad,nalpha,npml)
  
  IMPLICIT NONE
  
  integer :: ix1,ix2,ix3,ialpha_1,ialpha_2,ialpha_3,nalpha_frees
  integer :: nsp1,nsp2,nsp3
  integer :: irad,nalpha,npml
  real :: x1,x2,x3,x1_0,x2_0,x3_0,alpha_1,alpha_2,alpha_3,dalpha,xoversurf
  real :: h
  
  dalpha = 1. / float( nalpha )
  xoversurf = float(irad -1) * h
  
  ix1 = INT(x1 / h ) + 1 
  ix2 = INT(x2 / h ) + 1 
  ix3 = INT(x3 / h ) + 1 
  
  x1_0 = FLOAT(ix1 - 1)*h
  x2_0 = FLOAT(ix2 - 1)*h
  x3_0 = FLOAT(ix3 - 1)*h
    
  alpha_1 = ( x1 - x1_0 ) / h
  alpha_2 = ( x2 - x2_0 ) / h
  alpha_3 = ( x3 - x3_0 ) / h
  
  ialpha_1 = NINT ( alpha_1 / dalpha ) + 1
  ialpha_1 = MIN(nalpha,ialpha_1)
  ialpha_2 = NINT ( alpha_2 / dalpha ) + 1 
  ialpha_2 = MIN(nalpha,ialpha_2)
  ialpha_3 = NINT ( alpha_3 / dalpha ) + 1 
  ialpha_3 = MIN(nalpha,ialpha_3)

  nalpha_frees =  irad - 1 - INT(x1/h)
  nalpha_frees =  MAX ( nalpha_frees , 0 )
  
  ! ----------------------------------------------------------------------------
  ! WARNING - WARNING - WARNING - WARNING - WARNING - Account for PML
  ! ----------------------------------------------------------------------------
  nsp1=npml
  nsp2=npml
  nsp3=0
  ix1 = ix1 + nsp1
  ix2 = ix2 + nsp2
  ix3 = ix3 + nsp3
  
END SUBROUTINE sub_sample_Hicks


