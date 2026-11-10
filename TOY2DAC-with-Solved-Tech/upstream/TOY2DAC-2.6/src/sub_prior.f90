!Copyright 2008-2011 SEISCOPE project, All rights reserved.
!Copyright 2013-2021 SEISCOPEII project, All rights reserved.
!IMPLEMENTATION OF PRIOR INFORMATION REGULARIZATION TERMS 
subroutine sub_prior_fcost(pbdir,inv)
  
  IMPLICIT NONE
#include "common.h"
#include "pbdirect.h"
  include "inversion.h"

  !PBDIRECT
  TYPE (pbdirect) :: pbdir
  !INVERSION
  TYPE (inversion) :: inv
  real :: fcost_prior
  integer :: i1,i2
  
  fcost_prior=0.
  do i2=1,pbdir%n2     
     do i1=1,pbdir%n1
        fcost_prior=fcost_prior+&
             (inv%model(i1,i2,1)-inv%prior0(i1,i2))**2
     enddo
  enddo
  inv%fcost=inv%fcost+inv%lambda_prior0*0.5*fcost_prior  
!  write(15,*) 'inv%lambda_prior0 vaut : ',inv%lambda_prior0     
!  write(15,*) 'inv%fcost vaut : ',inv%fcost
!  write(15,*) 'fcost_prior vaut : ',inv%lambda_prior0*0.5*fcost_prior  
  
end subroutine sub_prior_fcost


subroutine sub_prior_fgrad(pbdir,inv)
  
  IMPLICIT NONE
#include "common.h"
#include "pbdirect.h"
  include "inversion.h"

  !PBDIRECT
  TYPE (pbdirect) :: pbdir
  !INVERSION
  TYPE (inversion) :: inv
  real,dimension(:,:),allocatable :: fgrad_prior
  integer :: i1,i2
  
  allocate(fgrad_prior(pbdir%n1,pbdir%n2))
  fgrad_prior(:,:)=0.
  do i2=1,pbdir%n2
     do i1=1,pbdir%n1     
        fgrad_prior(i1,i2)=inv%model(i1,i2,1)-inv%prior0(i1,i2)
     enddo
  enddo
  inv%gradient(:,:,1)=inv%gradient(:,:,1)+inv%lambda_prior0*fgrad_prior(:,:)  
  deallocate(fgrad_prior)
  
end subroutine sub_prior_fgrad

subroutine sub_prior_Hv(pbdir,inv,Hv,v)
  
  IMPLICIT NONE
#include "common.h"
#include "pbdirect.h"
  include "inversion.h"

  !PBDIRECT
  TYPE (pbdirect) :: pbdir
  !INVERSION
  TYPE (inversion) :: inv
  real,dimension(pbdir%n1*pbdir%n2) :: v,Hv
  real,allocatable,dimension(:) :: vect_prior 
  integer :: i1,i2
  
  Hv(:)=Hv(:)+inv%lambda_prior0*v(:)  
      
end subroutine sub_prior_Hv

