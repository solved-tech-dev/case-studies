!Copyright 2008-2011 SEISCOPE project, All rights reserved.
!Copyright 2013-2021 SEISCOPEII project, All rights reserved.
!***************************************!
! OPTIMIZATION ROUINTE INITIALIZATION   !
!***************************************!
subroutine sub_init_optim(pbdir,acqui,inv,optim)

  implicit none

  include "mpif.h"
#include "common.h"
#include "pbdirect.h"
  include "acqui.h"
  include "inversion.h"
  include "optim_type.h"

  !IN
  TYPE(acquisition) :: acqui   
  TYPE (pbdirect) :: pbdir  
  TYPE(inversion) :: inv
  !IN/OUT
  TYPE(optim_type) :: optim

  !Local variables
  integer :: nn,iipar

  !Compute the problem size nn
  nn=pbdir%n1*pbdir%n2*inv%npar  

  !GLOBAL OPTIONS  
  optim%conv=inv%convergence_criterion
  optim%niter_max=inv%niter_nonlin

  !PRINTING OPTIONS
  if(mype.eq.0) then  
     optim%print_flag=1
     if(inv%debug_option.eq.1) then
        optim%debug=.true.
     else
        optim%debug=.false.
     endif
  endif

  !BOUND CONSTRAINTS
  optim%bound=inv%bound 
  if(optim%bound.eq.1) then
     allocate(optim%ub(nn),optim%lb(nn))
     DO iipar=1,inv%npar
        optim%lb(pbdir%n1*pbdir%n2*(iipar-1)+1:pbdir%n1*pbdir%n2*(iipar))=inv%lb(iipar)
        optim%ub(pbdir%n1*pbdir%n2*(iipar-1)+1:pbdir%n1*pbdir%n2*(iipar))=inv%ub(iipar)

        optim%threshold=inv%threshold 
     END DO
  END IF

end subroutine sub_init_optim

