!Copyright 2008-2011 SEISCOPE project, All rights reserved.
!Copyright 2013-2021 SEISCOPEII project, All rights reserved.
!****************************************************************************************!
! THIS FILE CONTAINS THE SUBROUTINE ASSOCIATED WITH THE FREEZING STRATEGY:               !
! THE FILE fbathy INDICATES WHICH SUBSURFACE PARAMETER SHALL BE UNCHANDED,               !
! WHICH REDUCE THE SIZE OF THE LINEAR SYSTEM CORRESPONDING TO THE COMPUTATION OF         !
! THE NEWTON DESCENT DIRECTION                                                           !
!****************************************************************************************!


!********************************************************************************!
! This subroutine reduces computes the size of the reduced parameter space       !
! nn_tilde                                                                       !
!********************************************************************************!
subroutine sub_compute_nn_tilde (nn,pbdir,inv,nn_tilde)

  implicit none 
  
#include "common.h"
  include "inversion.h"
#include "pbdirect.h"
  !INVERSION
  TYPE (inversion) :: inv
  TYPE (pbdirect) :: pbdir
  
  integer :: nn,nn_tilde,ipar,i2
  
  nn_tilde=nn
  
  do ipar=1,inv%npar
     do i2=1,pbdir%n2
        !nn_tilde=nn_tilde-(inv%ibathy(i2)-1)
        nn_tilde=nn_tilde-(inv%ibathy(i2)+inv%ideadzone-1)
     enddo
  enddo
  
end subroutine sub_compute_nn_tilde


!********************************************************************************!
! This subroutine reduces the vector v in the parameter space                    !
! of size nn into the vector vtilde in the reduced parameter space of size       !
! nn_tilde.                                                                      !
!********************************************************************************!
subroutine sub_reduce_vector(nn,pbdir,acqui,inv,v,nn_tilde,vtilde)

  implicit none
  
#include "common.h"
#include "pbdirect.h"
  include "inversion.h"
  include "acqui.h"
  
  !PBDIRECT
  TYPE (pbdirect) :: pbdir
  !ACQUI
  TYPE(acquisition) :: acqui
  !INVERSION
  TYPE (inversion) :: inv
  
  integer :: nn,nn_tilde,ipar,i1,i2,k
  real,dimension(nn) :: v
  real,dimension(nn_tilde) :: vtilde
  
  k=1
  do ipar=1,inv%npar     
     do i2=1,pbdir%n2
        !do i1=inv%ibathy(i2),pbdir%n1
        do i1=inv%ibathy(i2)+inv%ideadzone,pbdir%n1
           vtilde(k)=v(i1+(i2-1)*(pbdir%n1)+(ipar-1)*pbdir%n1*pbdir%n2)
           k=k+1
        enddo
     enddo
  enddo
  
end subroutine sub_reduce_vector


!********************************************************************************!
! This subroutine is augments the vector vtilde in the reduced parameter space   !
! of size nn_tilde into the vector v in the parameter space of size nn by adding !
! corresponding 0 entries                                                        !
!********************************************************************************!
subroutine sub_augment_vector(nn,pbdir,acqui,inv,v,nn_tilde,vtilde)
  
  implicit none

#include "common.h"
#include "pbdirect.h"
  include "inversion.h"
  include "acqui.h"
  
  !PBDIRECT
  TYPE (pbdirect) :: pbdir
  !ACQUI
  TYPE(acquisition) :: acqui
  !INVERSION
  TYPE (inversion) :: inv
  
  integer :: nn,nn_tilde,ipar,i1,i2,k
  real,dimension(nn) :: v
  real,dimension(nn_tilde) :: vtilde
  
  k=1
  do ipar=1,inv%npar
     do i2=1,pbdir%n2
        !do i1=1,inv%ibathy(i2)-1
        do i1=1,inv%ibathy(i2)+inv%ideadzone-1
           v(i1+(i2-1)*(pbdir%n1)+(ipar-1)*pbdir%n1*pbdir%n2)=0.
        enddo
        !do i1=inv%ibathy(i2),pbdir%n1
        do i1=inv%ibathy(i2)+inv%ideadzone,pbdir%n1
           v(i1+(i2-1)*(pbdir%n1)+(ipar-1)*pbdir%n1*pbdir%n2)=vtilde(k)
           k=k+1
        enddo
     enddo
  enddo
end subroutine sub_augment_vector


subroutine sub_augment_model(nn,pbdir,acqui,inv,v,v_ref,nn_tilde,vtilde)
  
  implicit none

#include "common.h"
#include "pbdirect.h"
  include "inversion.h"
  include "acqui.h"
  
  !PBDIRECT
  TYPE (pbdirect) :: pbdir
  !ACQUI
  TYPE(acquisition) :: acqui
  !INVERSION
  TYPE (inversion) :: inv
  
  integer :: nn,nn_tilde,ipar,i1,i2,k
  real,dimension(nn) :: v,v_ref
  real,dimension(nn_tilde) :: vtilde
  
  k=1
  do ipar=1,inv%npar
     do i2=1,pbdir%n2
        !do i1=1,inv%ibathy(i2)-1
        do i1=1,inv%ibathy(i2)+inv%ideadzone-1
           v(i1+(i2-1)*(pbdir%n1)+(ipar-1)*pbdir%n1*pbdir%n2)=&
                v_ref(i1+(i2-1)*(pbdir%n1)+(ipar-1)*pbdir%n1*pbdir%n2)
        enddo
        !do i1=inv%ibathy(i2),pbdir%n1
        do i1=inv%ibathy(i2)+inv%ideadzone,pbdir%n1
           v(i1+(i2-1)*(pbdir%n1)+(ipar-1)*pbdir%n1*pbdir%n2)=vtilde(k)
           k=k+1
        enddo
     enddo
  enddo
end subroutine sub_augment_model

!**********************************************************!
! This subroutine converts a vector of the parameter space !
! into its 2D matrix representation                        !
!**********************************************************!
subroutine sub_vect2mat(nn,npar,pbdir,v,mat)
  
  implicit none 
  
#include "pbdirect.h"
  !INVERSION
  TYPE (pbdirect) :: pbdir
  real,dimension(nn) :: v
  real,dimension(pbdir%n1,pbdir%n2,npar) :: mat
  integer :: nn,npar,ipar,i1,i2

  do ipar=1,npar
     do i2=1,pbdir%n2
        do i1=1,pbdir%n1
           mat(i1,i2,ipar)=v(i1+(i2-1)*pbdir%n1+(ipar-1)*pbdir%n1*pbdir%n2)
        enddo
     enddo
  enddo
end subroutine sub_vect2mat
