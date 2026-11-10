!Copyright 2008-2011 SEISCOPE project, All rights reserved.
!Copyright 2013-2021 SEISCOPEII project, All rights reserved.
!SUBROUTINE TO READ INPUT MODEL FILES

SUBROUTINE  sub_read_models(file_vp,file_qp,file_rho,vp,qp,rho,n1,n2)

  IMPLICIT NONE

  INTEGER :: n1,n2
  REAL,DIMENSION(n1,n2):: vp,rho,qp
  CHARACTER(LEN=80) :: file_vp,file_rho,file_qp

  OPEN(13,file=file_vp,access='direct',recl=4*n1*n2)
  READ(13,rec=1)vp(:,:)
  CLOSE(13)
  OPEN(14,file=file_qp,access='direct',recl=4*n1*n2)
  READ(14,rec=1)qp(:,:)
  CLOSE(14)
  OPEN(15,file=file_rho,access='direct',recl=4*n1*n2)
  READ(15,rec=1)rho(:,:)
  CLOSE(15)

END SUBROUTINE SUB_READ_MODELS
