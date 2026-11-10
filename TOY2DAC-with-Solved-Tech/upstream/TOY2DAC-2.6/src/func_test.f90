! Copyright 2008-2011 SEISCOPE project, All rights reserved.
! Copyright 2013-2015 SEISCOPEII project, All rights reserved.

!############################################################################
!  func_test include different function to test equality, superiority,
!  inferiority of  real without direct test (== , <= >=) which can give 
!  problems. 
! Test is computed as a comparison with average cell length
!############################################################################

INTEGER FUNCTION func_test_equality(xa,xb)
  IMPLICIT NONE
#include "common.h"
  REAL  :: xa,xb  
  IF (ABS(xa-xb)/lmoy < epsilondist) THEN
     func_test_equality=0
  ELSE
     func_test_equality=1
  END IF
  RETURN
END FUNCTION func_test_equality
