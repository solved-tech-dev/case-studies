!Copyright 2008-2011 SEISCOPE project, All rights reserved.
!Copyright 2013-2021 SEISCOPEII project, All rights reserved.
!SUBOROUTINES TO COMPUTE GRADIENT FOR BASE PARAMETER -> VP, QP, DELTA, EPSILON
! OTHER PARAMETERIZATION WILL USE THEM 

REAL FUNCTION gradient_vp(omegac,rho,vp,qp,inc_wave_i,adj_wave_i,cc1)
  IMPLICIT NONE
  COMPLEX :: omegac
  REAL :: rho,vp,qp
  COMPLEX :: inc_wave_i,adj_wave_i,cc1

  COMPLEX :: xi

  xi=CMPLX(0.,1.)

  gradient_vp = - REAL(2.*omegac**2/(rho*vp**3 *(1.-0.5*xi/qp)**2 )*inc_wave_i*adj_wave_i*cc1/ABS(cc1))

END FUNCTION gradient_vp


REAL FUNCTION gradient_qp(omegac,rho,vp,qp,inc_wave_i,adj_wave_i,cc1)
  IMPLICIT NONE
  COMPLEX :: omegac
  REAL :: rho,vp,qp
  COMPLEX :: inc_wave_i,adj_wave_i,cc1

  COMPLEX :: xi

  xi=CMPLX(0.,1.)

  gradient_qp = - REAL(xi*omegac**2/(rho*vp**2 *(1.-0.5*xi/qp)**3 * qp**2)*inc_wave_i*adj_wave_i*cc1/ABS(cc1))

END FUNCTION gradient_qp


REAL FUNCTION gradient_epsilon(pbdir,acqui,inv,inc_wfld,adj_wfld,i1,i2,iw,nb_wfld)
  IMPLICIT NONE
#include "pbdirect.h"
#include "common.h"
  include "acqui.h"
  include "inversion.h"

  !PBDIRECT
  TYPE (pbdirect) :: pbdir
  !ACQUI
  TYPE(acquisition) :: acqui
  !INVERSION
  TYPE (inversion) :: inv
  COMPLEX,DIMENSION(pbdir%n1e,pbdir%n2e,nb_wfld,acqui%nsrc) :: inc_wfld,adj_wfld
  INTEGER :: i1,i2,isrc,iw,nb_wfld

  REAL ::   h2,gx,gz,hx,hz,b00,bp0,bm0,b0p,b0m

  REAL:: tmp

  gradient_epsilon=0.

  gx=2.*cos(pbdir%theta(i1,i2))*cos(pbdir%theta(i1,i2))
  hx=-sin(2.*pbdir%theta(i1,i2))
  gz=-sin(2.*pbdir%theta(i1,i2))
  hz=2.*sin(pbdir%theta(i1,i2))*sin(pbdir%theta(i1,i2))

  IF(i1<pbdir%n1)THEN
     bp0=0.5*(1./pbdir%rho(i1,i2)+1./pbdir%rho(i1+1,i2))
  ELSE
     bp0=1./pbdir%rho(i1,i2)
  END IF
  IF(i1>1)THEN                 
     bm0=0.5*(1./pbdir%rho(i1,i2)+1./pbdir%rho(i1-1,i2))
  ELSE
     bm0=1./pbdir%rho(i1,i2)
  END IF
  IF(i2<pbdir%n2)THEN
     b0p=0.5*(1./pbdir%rho(i1,i2)+1./pbdir%rho(i1,i2+1))
  ELSE
     b0p=1./pbdir%rho(i1,i2)
  END IF
  IF(i2>1)THEN        
     b0m=0.5*(1./pbdir%rho(i1,i2)+1./pbdir%rho(i1,i2-1))
  ELSE
     b0m=1./pbdir%rho(i1,i2)
  END IF

  h2=1./(pbdir%h*pbdir%h)

  DO isrc=1,acqui%nsrc
     !term sur ligne qi= -2Hq -2 Hp    -> H du papier de Stephane
     tmp=tmp + REAL(&
                                !term 00 sur q
          ((-h2*gx*(b0p+b0m)-h2*hz*(bp0+bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2,2,isrc) +&
                                !term 01 sur q
          (h2*gx*b0p + 0.25*gz*h2*(bp0-bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2+1,2,isrc) +&
                                !term 0-1 sur q
          (h2*gx*b0m + 0.25*gz*h2*(-bp0+bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2-1,2,isrc) +&
                                !term -10 sur q
          (h2*hz*bm0 + 0.25*hx*h2*(-b0p+b0m) )*inc_wfld(pbdir%npml+i1-1,pbdir%npml+i2,2,isrc) +&
                                !term 10 sur q
          (h2*hz*bp0 + 0.25*hx*h2*(-b0m+b0p) )*inc_wfld(pbdir%npml+i1+1,pbdir%npml+i2,2,isrc) +&
                                !term 00 sur p
          (-h2*gx*(b0p+b0m)-h2*hz*(bp0+bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2,1,isrc) +&
                                !term 01 sur p
          (h2*gx*b0p + 0.25*gz*h2*(bp0-bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2+1,1,isrc) +&
                                !term 0-1 sur p
          (h2*gx*b0m + 0.25*gz*h2*(-bp0+bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2-1,1,isrc) +&
                                !term -10 sur p
          (h2*hz*bm0 + 0.25*hx*h2*(-b0p+b0m) )*inc_wfld(pbdir%npml+i1-1,pbdir%npml+i2,1,isrc) +&
                                !term 10 sur p
          (h2*hz*bp0 + 0.25*hx*h2*(-b0m+b0p) )*inc_wfld(pbdir%npml+i1+1,pbdir%npml+i2,1,isrc) ) &
                                !multiply the qi term of the adjoint wavefield
          * adj_wfld(pbdir%npml+i1,pbdir%npml+i2,2,isrc)*inv%cc1(isrc,iw)/ABS(inv%cc1(isrc,iw)) )

  END DO

  gradient_epsilon=tmp

END FUNCTION gradient_epsilon


REAL FUNCTION gradient_delta(pbdir,acqui,inv,inc_wfld,adj_wfld,i1,i2,iw,nb_wfld)
  IMPLICIT NONE
#include "pbdirect.h"
#include "common.h"
  include "acqui.h"
  include "inversion.h"

  !PBDIRECT
  TYPE (pbdirect) :: pbdir
  !ACQUI
  TYPE(acquisition) :: acqui
  !INVERSION
  TYPE (inversion) :: inv

  COMPLEX,DIMENSION(pbdir%n1e,pbdir%n2e,nb_wfld,acqui%nsrc) :: inc_wfld,adj_wfld
  INTEGER :: i1,i2,isrc,iw,nb_wfld

  REAL ::   h2,gx,gz,hx,hz,b00,bp0,bm0,b0p,b0m

  REAL:: tmp

  gradient_delta=0.

  !bellow are the g? and h? derived with regard of pbdir%delta
  gx=-2.*cos(pbdir%theta(i1,i2))*cos(pbdir%theta(i1,i2))
  hx=sin(2.*pbdir%theta(i1,i2))
  gz=sin(2.*pbdir%theta(i1,i2))
  hz=-2.*sin(pbdir%theta(i1,i2))*sin(pbdir%theta(i1,i2))

  IF(i1<pbdir%n1)THEN
     bp0=0.5*(1./pbdir%rho(i1,i2)+1./pbdir%rho(i1+1,i2))
  ELSE
     bp0=1./pbdir%rho(i1,i2)
  END IF
  IF(i1>1)THEN                 
     bm0=0.5*(1./pbdir%rho(i1,i2)+1./pbdir%rho(i1-1,i2))
  ELSE
     bm0=1./pbdir%rho(i1,i2)
  END IF
  IF(i2<pbdir%n2)THEN
     b0p=0.5*(1./pbdir%rho(i1,i2)+1./pbdir%rho(i1,i2+1))
  ELSE
     b0p=1./pbdir%rho(i1,i2)
  END IF
  IF(i2>1)THEN        
     b0m=0.5*(1./pbdir%rho(i1,i2)+1./pbdir%rho(i1,i2-1))
  ELSE
     b0m=1./pbdir%rho(i1,i2)
  END IF

  h2=1./(pbdir%h*pbdir%h)

  DO isrc=1,acqui%nsrc
     !term sur ligne pi= -2Hq -2 Hp    -> H du papier de Stephane -> we take -gx -gz -hx and -hz inthe following lines
     tmp=tmp + REAL(&
                                !term 00 sur q
          ((-h2*(-gx)*(b0p+b0m)-h2*(-hz)*(bp0+bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2,2,isrc) +&
                                !term 01 sur q
          (h2*(-gx)*b0p + 0.25*(-gz)*h2*(bp0-bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2+1,2,isrc) +&
                                !term 0-1 sur q
          (h2*(-gx)*b0m + 0.25*(-gz)*h2*(-bp0+bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2-1,2,isrc) +&
                                !term -10 sur q
          (h2*(-hz)*bm0 + 0.25*(-hx)*h2*(-b0p+b0m) )*inc_wfld(pbdir%npml+i1-1,pbdir%npml+i2,2,isrc) +&
                                !term 10 sur q
          (h2*(-hz)*bp0 + 0.25*(-hx)*h2*(-b0m+b0p) )*inc_wfld(pbdir%npml+i1+1,pbdir%npml+i2,2,isrc) +&
                                !term 00 sur p
          (-h2*(-gx)*(b0p+b0m)-h2*(-hz)*(bp0+bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2,1,isrc) +&
                                !term 01 sur p
          (h2*(-gx)*b0p + 0.25*(-gz)*h2*(bp0-bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2+1,1,isrc) +&
                                !term 0-1 sur p
          (h2*(-gx)*b0m + 0.25*(-gz)*h2*(-bp0+bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2-1,1,isrc) +&
                                !term -10 sur p
          (h2*(-hz)*bm0 + 0.25*(-hx)*h2*(-b0p+b0m) )*inc_wfld(pbdir%npml+i1-1,pbdir%npml+i2,1,isrc) +&
                                !term 10 sur p
          (h2*(-hz)*bp0 + 0.25*(-hx)*h2*(-b0m+b0p) )*inc_wfld(pbdir%npml+i1+1,pbdir%npml+i2,1,isrc) ) &
                                !multiply the pi term of the adjoint wavefield
          * adj_wfld(pbdir%npml+i1,pbdir%npml+i2,1,isrc)*inv%cc1(isrc,iw)/ABS(inv%cc1(isrc,iw)) )
     !write(*,*)inv%gradient(i1,i2,iipar)

     !term sur ligne qi= 2Hq +2 Hp    -> H du papier de Stephane
     tmp=tmp + REAL(&
                                !term 00 sur q
          ((-h2*gx*(b0p+b0m)-h2*hz*(bp0+bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2,2,isrc) +&
                                !term 01 sur q
          (h2*gx*b0p + 0.25*gz*h2*(bp0-bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2+1,2,isrc) +&
                                !term 0-1 sur q
          (h2*gx*b0m + 0.25*gz*h2*(-bp0+bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2-1,2,isrc) +&
                                !term -10 sur q
          (h2*hz*bm0 + 0.25*hx*h2*(-b0p+b0m) )*inc_wfld(pbdir%npml+i1-1,pbdir%npml+i2,2,isrc) +&
                                !term 10 sur q
          (h2*hz*bp0 + 0.25*hx*h2*(-b0m+b0p) )*inc_wfld(pbdir%npml+i1+1,pbdir%npml+i2,2,isrc) +&
                                !term 00 sur p
          (-h2*gx*(b0p+b0m)-h2*hz*(bp0+bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2,1,isrc) +&
                                !term 01 sur p
          (h2*gx*b0p + 0.25*gz*h2*(bp0-bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2+1,1,isrc) +&
                                !term 0-1 sur p
          (h2*gx*b0m + 0.25*gz*h2*(-bp0+bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2-1,1,isrc) +&
                                !term -10 sur p
          (h2*hz*bm0 + 0.25*hx*h2*(-b0p+b0m) )*inc_wfld(pbdir%npml+i1-1,pbdir%npml+i2,1,isrc) +&
                                !term 10 sur p
          (h2*hz*bp0 + 0.25*hx*h2*(-b0m+b0p) )*inc_wfld(pbdir%npml+i1+1,pbdir%npml+i2,1,isrc) ) &
                                !multiply the qi term of the adjoint wavefield
          * adj_wfld(pbdir%npml+i1,pbdir%npml+i2,2,isrc)*inv%cc1(isrc,iw)/ABS(inv%cc1(isrc,iw)) )
     !write(*,*)inv%gradient(i1,i2,iipar)
  END DO

  gradient_delta=tmp

END FUNCTION gradient_delta


subroutine sub_gradient_rho(i1,i2,iw,iipar,nb_wfld,nsrc,omegac,pbdir,inv,inc_wfld,adj_wfld,grad_term)

  implicit none
  include "mpif.h"
#include "common.h"
#include "pbdirect.h"
  include "inversion.h"

  !PBDIRECT
  TYPE (pbdirect) :: pbdir
  !INVERSION
  TYPE (inversion) :: inv

  INTEGER :: i1,i2,nsrc,ianiso,iw,isrc,nb_wfld,iipar
  COMPLEX :: grad_term
  COMPLEX  :: omegac
  COMPLEX,DIMENSION(pbdir%n1e,pbdir%n2e,nb_wfld,nsrc) :: inc_wfld,adj_wfld
  REAL,DIMENSION(-2:2,-2:2) :: b
  REAL,DIMENSION(-1:1,-1:1) :: epsilo,theta,delta
  COMPLEX,DIMENSION(-2:2,-2:2) :: mu1
  complex :: vc,xi,val
  INTEGER :: i1e,i2e,ne,nnz,ne2,nnz2
  INTEGER,DIMENSION(100) :: icn,irn
  COMPLEX,DIMENSION(100) :: mat
  INTEGER,DIMENSION(10) :: ivect2
  COMPLEX,DIMENSION(10) :: vect2
  INTEGER :: mode
  ne=100
  ne2=10

  xi=CMPLX(0.,1.)


  i1e=i1+pbdir%npml
  i2e=i2+pbdir%npml

  if(pbdir%tool == 1) ianiso = 0
  if(pbdir%tool == 2) ianiso = 5

  mu1(:,:)=0.
  b(:,:)=0.
  IF(pbdir%tool == 1) THEN
     epsilo(:,:)=0.
     delta(:,:)=0.
     theta(:,:)=0. 
  ELSE
     epsilo(:,:)=pbdir%epsilon(i1,i2)
     delta(:,:)=pbdir%delta(i1,i2)
     theta(:,:)=pbdir%theta(i1,i2)
  END IF
  vc=pbdir%vp(i1,i2)*(1.-0.5*xi/pbdir%qp(i1,i2))
  mu1(0,0)=1./(pbdir%rho(i1,i2)**2*vc*vc) !=(1 /rho^2 vp^2)
  b(0,0)=1./pbdir%rho(i1,i2)**2 ! 1/rho^2
  mode=2

  ! write(*,*)'mu1(0,0)',mu1(0,0)
  grad_term=0.
  !step1 : dAdmi
  CALL sub_dAdmi(i1e,i2e,b,mu1,delta,epsilo,theta,ianiso,ne, irn,icn,mat,pbdir%h,omegac,nnz,pbdir%n1e,pbdir%n2e,mode)
  mat(:)=-mat(:) !minus sign for consistency with other dAdmi computation

  DO isrc=1,nsrc
     !step 2 : sparse matrix vector product -> sparse vector
     !write(*,*)'mat',mat
     !   write(*,*)'pbdir%n1e*pbdir%n2e*nb_wfld',pbdir%n1e,pbdir%n2e,nb_wfld,pbdir%n1e*pbdir%n2e*nb_wfld
     CALL sub_P_dAdmi(mat,irn,icn,nnz,ne, inc_wfld(:,:,:,isrc),pbdir%n1e*pbdir%n2e*nb_wfld,vect2,ivect2,ne2,nnz2)
     !step 3 : sparse vector/dense vector product
     !write(*,*)'vect2',vect2

     CALL sub_P_dAdmi_P(adj_wfld(:,:,:,isrc),pbdir%n1e*pbdir%n2e*nb_wfld,vect2,ivect2,ne2,nnz2,val)
     !step4 : take into account source term
     grad_term=grad_term+val*inv%cc1(isrc,iw)/ABS(inv%cc1(isrc,iw))
     !write(*,*)'grad_term',grad_term
  END DO

end subroutine sub_gradient_rho
