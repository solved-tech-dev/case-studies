!Copyright 2008-2011 SEISCOPE project, All rights reserved.
!Copyright 2013-2021 SEISCOPEII project, All rights reserved.
!SUBROUTINE TO FILL THE DENSE RHS OF MUMPS FOR SOURCES
SUBROUTINE sub_fill_RHS_SOURCES(RHS,n1,n2,npml,h,nsrc,coord_src)

  IMPLICIT NONE
#include "common.h"
  INTEGER :: n1,n2,nsrc,npml
  COMPLEX(kind=creal) :: RHS(n1*n2*nsrc)
  REAL :: h
  REAL :: coord_src(nsrc,2),scaling_factor
  INTEGER :: isrc
  INTEGER :: ix,iz,k

  RHS(:)=0.  
  DO isrc = 1, nsrc
     ix = NINT(coord_src(isrc,1)/h)+1+npml
     iz = NINT(coord_src(isrc,2)/h)+1+npml
     !ix = INT(coord_src(isrc,1)/h)+1+npml
     !iz = INT(coord_src(isrc,2)/h)+1+npml
     !write(99,*) 'ix, iz are  : ',ix, iz
     k=(isrc-1)*(n1*n2)+(ix-1)*n1+iz
     RHS(k)=1./(h**2)          
  END DO
END SUBROUTINE SUB_FILL_RHS_SOURCES

SUBROUTINE sub_fill_RHS_SOURCESaniso(RHS,n1,n2,npml,h,nsrc,coord_src)
  
  IMPLICIT NONE
#include "common.h" 
  INTEGER :: n1,n2,nsrc,npml
  COMPLEX(kind=creal) :: RHS(n1*n2*nsrc*2)
  REAL :: h
  REAL :: coord_src(nsrc,2)
  INTEGER :: isrc
  INTEGER :: ix,iz,k

  RHS(:)=0.
  DO isrc = 1, nsrc
     ix = NINT(coord_src(isrc,1)/h)+1+npml
     iz = NINT(coord_src(isrc,2)/h)+1+npml
     k=(isrc-1)*(n1*n2*2)+(ix-1)*n1+iz
     RHS(k)=1./(h**2)     
  END DO  
END SUBROUTINE SUB_FILL_RHS_SOURCESANISO
