!Copyright 2008-2011 SEISCOPE project, All rights reserved.
!Copyright 2013-2021 SEISCOPEII project, All rights reserved.
!MANAGEMENT OF PML COEF

SUBROUTINE sub_PML(damp1p,damp1pb,damp2p,damp2pb,n1,n2,h,apml,npml,&
     omega,omegac)

  IMPLICIT NONE

  COMPLEX,DIMENSION(0:n1+1) :: damp1p,damp1pb
  COMPLEX,DIMENSION(0:n2+1) :: damp2p,damp2pb
  INTEGER :: n1,n2
  REAL :: h
  REAL :: apml,omega,betad,alphad,bpml
  COMPLEX :: omegac
  INTEGER :: idp1,idpb1,idp2,idpb2,itd,npml

  idp1=1
  idpb1=1
  idp2=1
  idpb2=1
  alphad=1.
  betad=0.
  itd=2
  bpml=0.

  call subdamp(damp1p,damp1pb,n1,npml,h,apml,bpml,omega,idp1,idpb1,alphad,betad,itd,omegac)
  call subdamp(damp2p,damp2pb,n2,npml,h,apml,bpml,omega,idp2,idpb2,alphad,betad,itd,omegac)

END SUBROUTINE sub_PML
