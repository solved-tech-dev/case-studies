!Copyright 2008-2011 SEISCOPE project, All rights reserved.
!Copyright 2013-2021 SEISCOPEII project, All rights reserved.
!SUBROUTINE TO MANAGE COMPLEX VELOCITY
SUBROUTINE sub_vcmplx(vp,qp,rho,mu,n1,n2)
  IMPLICIT NONE


  REAL,DIMENSION(n1,n2) :: vp,qp,rho
  COMPLEX,DIMENSION(0:n1+1,0:n2+1) :: mu
  INTEGER :: n1,n2,npml
  REAL :: omega,omegar
  INTEGER :: imq

  !hard code the classic attenuation model
  imq=1

  !parameters not used in subv2mu
  npml=0
  omega=1.
  omegar=1.

  CALL sub_v2mu(rho,vp,qp,mu,n1,n2,npml,npml,omega,omegar,imq)


END SUBROUTINE sub_vcmplx
