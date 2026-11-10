!Copyright 2008-2011 SEISCOPE project, All rights reserved.
!Copyright 2013-2021 SEISCOPEII project, All rights reserved.

subroutine sub_v2mu(rho,v,q,mu,n1e,n2e,npml1,npml2,&
     omega,omegar,imq)
  IMPLICIT NONE
  integer n1e,n2e,imq,i2,i1,npml1,npml2
  real v(n1e,n2e),rho(n1e,n2e),q(n1e,n2e),omega,omegar
  complex mu(0:n1e+1,0:n2e+1)
  complex vc,ci,vci

  real pi
  parameter (pi=3.141592654)

  ci=cmplx(0.,1.)


  write(*,*) 'SUBROUTINE SUBV2MU'


  if (imq.eq.1) then
     do i2=1,n2e
        do i1=1,n1e
           vc=v(i1,i2)*(1.-0.5*ci/q(i1,i2))
           mu(i1,i2)=rho(i1,i2)*vc*vc
        end do
     end do
  else if (imq.eq.2) then
     do i2=1,n2e
        do i1=1,n1e
           vci=cmplx((1./v(i1,i2))+(1./(pi*v(i1,i2)*q(i1,i2)))&
                *log(omegar/omega),0.5/(v(i1,i2)*q(i1,i2)))
           vc=1./vci
           mu(i1,i2)=rho(i1,i2)*vc*vc
        end do
     end do
  end if

  do i1=1,n1e
     mu(i1,0)=mu(i1,1)
     mu(i1,n2e+1)=mu(i1,n2e)
  end do

  do i2=1,n2e
     mu(0,i2)=mu(1,i2)
     mu(n1e+1,i2)=mu(n1e,i2)
  end do

  mu(0,0)=mu(1,1)
  mu(0,n2e+1)=mu(1,n2e)
  mu(n1e+1,0)=mu(n1e,1)
  mu(n1e+1,n2e+1)=mu(n1e,n2e)

!OPEN(88,file='tmp_vp',access='direct',recl=4*n1e*n2e)
!WRITE(88,rec=1)v
!CLOSE(88)


end SUBROUTINE sub_v2mu


