!Copyright 2008-2011 SEISCOPE project, All rights reserved.
!Copyright 2013-2021 SEISCOPEII project, All rights reserved.
!NEW IMPLEMENTATION OF TIKHONOV REGULARIZATION 

subroutine sub_Tikhonov_fcost(pbdir,inv)
  
  IMPLICIT NONE
#include "common.h"
#include "pbdirect.h"
  include "inversion.h"

  !PBDIRECT
  TYPE (pbdirect) :: pbdir
  !INVERSION
  TYPE (inversion) :: inv
  
  real :: fcost_tikho
  integer :: i,i1,i2

  !STORE FCOST DATA
  inv%fcost_data=inv%fcost/inv%scalingfactor
  
  !COMPUTE REGUL IN Z DIRECTION
  fcost_tikho=0.
  do i=1,inv%npar
     do i2=1,pbdir%n2     
        do i1=inv%ibathy(i2),pbdir%n1-1        
           fcost_tikho=fcost_tikho+&
                inv%lambda(i)*(inv%model(i1+1,i2,i)-inv%model(i1,i2,i))**2
        enddo
     enddo
  enddo  
  inv%fcost=inv%fcost+1e0/(pbdir%h**2)*inv%lambda_z*0.5*fcost_tikho  

  !STORE FCOST REG Z
  inv%fcost_reg=1e0/(pbdir%h**2)*inv%lambda_z*0.5*fcost_tikho

  !COMPUTE REGUL IN X DIRECTION
  fcost_tikho=0.
  do i=1,inv%npar
     do i2=1,pbdir%n2-1     
        do i1=inv%ibathy(i2),pbdir%n1-1        
           fcost_tikho=fcost_tikho+&
                inv%lambda(i)*(inv%model(i1,i2+1,i)-inv%model(i1,i2,i))**2
        enddo
     enddo
  enddo
  inv%fcost=inv%fcost+1e0/(pbdir%h**2)*inv%lambda_x*0.5*fcost_tikho  
  
  !STORE FCOST REG X
  inv%fcost_reg=inv%fcost_reg+&
       1e0/(pbdir%h**2)*inv%lambda_x*0.5*fcost_tikho  

  !NORMALIZE FCOST REG
  inv%fcost_reg=inv%fcost_reg/inv%scalingfactor
  
end subroutine sub_Tikhonov_fcost

subroutine sub_Tikhonov_fgrad(pbdir,inv)
  
  IMPLICIT NONE
#include "common.h"
#include "pbdirect.h"
  include "inversion.h"

  !PBDIRECT
  TYPE (pbdirect) :: pbdir
  !INVERSION
  TYPE (inversion) :: inv

  real,allocatable,dimension(:,:,:) :: fgrad_tikho_x,fgrad_tikho_z
  integer :: i,i1,i2

  !Vertical term
  allocate(fgrad_tikho_z(pbdir%n1,pbdir%n2,inv%npar))
  fgrad_tikho_z(:,:,:)=0.
  do i=1,inv%npar
     do i2=1,pbdir%n2     
        do i1=inv%ibathy(i2)+1,pbdir%n1-1        
           fgrad_tikho_z(i1,i2,i)=-inv%model(i1+1,i2,i)+2.*inv%model(i1,i2,i)-inv%model(i1-1,i2,i)
        enddo
     enddo
  enddo
  do i=1,inv%npar
     do i2=1,pbdir%n2             
        fgrad_tikho_z(inv%ibathy(i2),i2,i)=-inv%model(inv%ibathy(i2)+1,i2,i)+inv%model(inv%ibathy(i2),i2,i)             
        fgrad_tikho_z(pbdir%n1,i2,i)=inv%model(pbdir%n1,i2,i)-inv%model(pbdir%n1-1,i2,i)
     enddo
  enddo
  do i=1,inv%npar
     inv%gradient(:,:,i)=inv%gradient(:,:,i)+inv%lambda(i)/(pbdir%h**2)*inv%lambda_z*fgrad_tikho_z(:,:,i)       
  enddo
  deallocate(fgrad_tikho_z)
  
  !Horizontal term
  allocate(fgrad_tikho_x(pbdir%n1,pbdir%n2,inv%npar))
  fgrad_tikho_x(:,:,:)=0.
  do i=1,inv%npar
     do i2=2,pbdir%n2-1
        do i1=inv%ibathy(i2),pbdir%n1-1        
           fgrad_tikho_x(i1,i2,i)=-inv%model(i1,i2+1,i)+2.*inv%model(i1,i2,i)-inv%model(i1,i2-1,i)
        enddo
     enddo
  enddo
  do i=1,inv%npar
     do i1=inv%ibathy(1),pbdir%n1-1     
        fgrad_tikho_x(i1,1,i)=-inv%model(i1,2,i)+inv%model(i1,1,i)
     enddo
     do i1=inv%ibathy(pbdir%n2),pbdir%n1-1
        fgrad_tikho_x(i1,pbdir%n2,i)=inv%model(i1,pbdir%n2,i)-inv%model(i1,pbdir%n2-1,i)
     enddo
  enddo
  do i=1,inv%npar
     inv%gradient(:,:,i)=inv%gradient(:,:,i)+inv%lambda(i)/(pbdir%h**2)*inv%lambda_x*fgrad_tikho_x(:,:,i)
  enddo
  deallocate(fgrad_tikho_x)
  
end subroutine sub_Tikhonov_fgrad

subroutine sub_Tikhonov_Hv(pbdir,inv,Hv,v)
  
  IMPLICIT NONE
#include "common.h"
#include "pbdirect.h"
  include "inversion.h"

  !PBDIRECT
  TYPE (pbdirect) :: pbdir
  !INVERSION
  TYPE (inversion) :: inv
  real,dimension(pbdir%n1*pbdir%n2*inv%npar) :: v,Hv
  real,allocatable,dimension(:) :: vect_x,vect_z
  integer :: i,i1,i2
  
   !Vertical term
  allocate(vect_z(pbdir%n1*pbdir%n2*inv%npar))
  vect_z(:)=0.
  do i=1,inv%npar     
     do i2=1,pbdir%n2     
        do i1=inv%ibathy(i2)+1,pbdir%n1-1
           vect_z(i1+(i2-1)*pbdir%n1+(i-1)*pbdir%n1*pbdir%n2)=&
                inv%lambda(i)*(&
                -v(i1+1+(i2-1)*pbdir%n1+(i-1)*pbdir%n1*pbdir%n2)+&
                2.*v(i1+(i2-1)*pbdir%n1+(i-1)*pbdir%n1*pbdir%n2)-&
                v(i1-1+(i2-1)*pbdir%n1+(i-1)*pbdir%n1*pbdir%n2)&
                )
        enddo
     enddo
     do i2=1,pbdir%n2     
        vect_z(inv%ibathy(i2)+(i2-1)*pbdir%n1+(i-1)*pbdir%n1*pbdir%n2)=&
             inv%lambda(i)*(&
             -v(inv%ibathy(i2)+1+(i2-1)*pbdir%n1+(i-1)*pbdir%n1*pbdir%n2)+&
             v(inv%ibathy(i2)+(i2-1)*pbdir%n1+(i-1)*pbdir%n1*pbdir%n2)&
             )
        vect_z(pbdir%n1+(i2-1)*pbdir%n1+(i-1)*pbdir%n1*pbdir%n2)=&
             inv%lambda(i)*(&
             v(pbdir%n1+(i2-1)*pbdir%n1+(i-1)*pbdir%n1*pbdir%n2)-&
             v(pbdir%n1-1+(i2-1)*pbdir%n1+(i-1)*pbdir%n1*pbdir%n2)&
             )
     enddo
  enddo
  Hv(:)=Hv(:)+1e0/(pbdir%h**2)*inv%lambda_z*vect_z(:)
  deallocate(vect_z)
  
  !Horizontal term  
  allocate(vect_x(pbdir%n1*pbdir%n2*inv%npar))
  vect_x(:)=0.
  do i=1,inv%npar  
     do i2=2,pbdir%n2-1
        do i1=inv%ibathy(i2),pbdir%n1  
           vect_x(i1+(i2-1)*pbdir%n1+(i-1)*pbdir%n1*pbdir%n2)=&
                inv%lambda(i)*(&
                -v(i1+(i2+1-1)*pbdir%n1+(i-1)*pbdir%n1*pbdir%n2)+&
                2.*v(i1+(i2-1)*pbdir%n1+(i-1)*pbdir%n1*pbdir%n2)-&
                v(i1+(i2-1-1)*pbdir%n1+(i-1)*pbdir%n1*pbdir%n2)&
                )          
        enddo
     enddo
     do i1=inv%ibathy(1),pbdir%n1  
        vect_x(i1+(i-1)*pbdir%n1*pbdir%n2)=&
             inv%lambda(i)*(&
             -v(i1+(2-1)*pbdir%n1+(i-1)*pbdir%n1*pbdir%n2)+&
             v(i1+(i-1)*pbdir%n1*pbdir%n2)&
             )       
     enddo
     do i1=inv%ibathy(pbdir%n2),pbdir%n1  
        vect_x(i1+(pbdir%n2-1)*pbdir%n1+(i-1)*pbdir%n1*pbdir%n2)=&
             inv%lambda(i)*(&
             v(i1+(pbdir%n2-1)*pbdir%n1+(i-1)*pbdir%n1*pbdir%n2)-&
             v(i1+(pbdir%n2-1-1)*pbdir%n1+(i-1)*pbdir%n1*pbdir%n2)&
             )
     enddo
  enddo
  Hv(:)=Hv(:)+1e0/(pbdir%h**2)*inv%lambda_x*vect_x(:)  
  deallocate(vect_x)
     
end subroutine sub_Tikhonov_Hv

