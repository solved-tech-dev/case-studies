!Copyright 2008-2011 SEISCOPE project, All rights reserved.
!Copyright 2013-2021 SEISCOPEII project, All rights reserved.

!####################################################################
! SUBORUTINE TO KILL THE CODE THROUGH A MPI_ABORT CALL IN CASE OF 
! PROBLEM DETECTED ON A SINGLE PROCESS
!####################################################################

SUBROUTINE sub_MPI_ABORT()


  IMPLICIT NONE

  include "mpif.h"
#include "common.h"

  INTEGER :: errorcode,ierror
  WRITE(*,*)'######################################################'
  WRITE(*,*)' '
  WRITE(*,*)'CODE IS ASKED TO BE KILLED NOW BY PROCESS ',mype
  WRITE(*,*)' '
  WRITE(*,*)'******************************************************'
  WRITE(*,*)' '
  WRITE(*,*)'A MPI_ABORT CALL IS SENT NOW TO KILL ALL PROCESS'
  WRITE(*,*)' '
  WRITE(*,*)'               BYE BYE'
  WRITE(*,*)' '
  WRITE(*,*)'######################################################'
  CALL flushc()

  CALL MPI_ABORT(MPI_COMM_WORLD,errorcode,ierror)

END SUBROUTINE sub_MPI_ABORT
