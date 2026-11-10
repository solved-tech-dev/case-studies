! Copyright 2008-2011 SEISCOPE project, All rights reserved.

!############################################################################
!  SUBROUTINE TO WRITE IN OUTPUT FILE ALL BUFFERS OF ALL 
!             PROC IN A SEQUENTIAL ORDER
!############################################################################

SUBROUTINE flush()
  IMPLICIT NONE
  include "mpif.h"
#include "common.h"

  INTEGER  :: i,infompi

  DO i = 1, nproc
     CALL MPI_BARRIER(MPI_COMM_WORLD,infompi)
     IF ( mype == (i-1) ) THEN
        CALL flushc()
     END IF
  END DO

END SUBROUTINE flush
