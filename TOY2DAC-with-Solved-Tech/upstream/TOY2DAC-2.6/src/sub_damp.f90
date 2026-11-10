!Copyright 2008-2011 SEISCOPE project, All rights reserved.
!Copyright 2013-2021 SEISCOPEII project, All rights reserved.
!
! subroutine subdamp.f
! set PML functions
! damp is discretized on the reference grid
! dampb is discretized on staggered grid.
! damp=1.+i . gamma/omega where gamma ia a smoothly
! decreasing function
!
! id: build (1) or not (0) damping function for the reference grid
! idb: build (1) or not (0) damping function for the staggered grid
! itd: type of damping function (1: polynomial else cosinus)
! alphad: add a constant to the xi function
! betad: add an imaginary part to the frequency


subroutine subdamp(damp,dampb,n,npml,dz,a,b,omega,id,idb,alphad,betamax,itd,omegac)
  IMPLICIT NONE
  complex damp(0:n+1),dampb(0:n+1),ci,omegac,test
  real omega,dz,a,b,pi,xpml,xmax,x,xb,eps,epsb,pi2,xb1
  real alphad,betamax,alpha,alphab,beta,betab,betam,betad,betadb
  integer n,npml,itd,i,id,idb
!WRITE(*,*)'itd',itd,omegac

  pi=3.14159265
  pi2=pi/2
  ci=cmplx(0.,1.)

  !       omegac=omegac+cmplx(omega,betad*omega)

  ! values of damp and dampn outside PML (equal to 1)

  do i=0,n+1
     damp(i)=cmplx(1.,0.)
     dampb(i)=cmplx(1.,0.)
  end do

  betam=betamax*omega
  !
  xpml=float(npml)*dz
  xmax=float(n-1)*dz

  !   x     x    x    x    x    x    x
  !      xb   xb   xb   xb   xb   xb

  do i=1,npml

     x=float(i-1)*dz
     xb=float(i-1)*dz+0.5*dz

     betad=x*betam/xpml
     betadb=xb*betam/xpml

     !  compute eps and alpha

     if (itd.eq.1) then
        eps=a*exp(b*log(((xpml-x)/xpml)+1e-10))
        alpha=(alphad-1.)*exp(b*log(((xpml-x)/xpml)+1e-10))+1.
     else
        eps=a*(1.-cos((xpml-x)*pi2/xpml))
        alpha=(alphad-1.)*(1.-cos((xpml-x)*pi2/xpml))+1.
     end if

     ! compute epsb and alphab

     if (itd.eq.1) then
        epsb=a*exp(b*log(((xpml-xb)/xpml)+1e-10))
        alphab=(alphad-1.)*exp(b*log(((xpml-xb)/xpml)+1e-10))+1.
     else
        epsb=a*(1.-cos((xpml-xb)*pi2/xpml))
        alphab=(alphad-1.)*(1.-cos((xpml-xb)*pi2/xpml))+1.
     end if

     if (abs(betamax).lt.1e-3) then
	betad=0.
	betadb=0.
     end if

     if (id.eq.1) then
        damp(i)=1./(alpha+ci*eps/(omegac+ci*betad))
     end if

     if (idb.eq.1) then
        dampb(i)=1./(alphab+ci*epsb/(omegac+ci*betadb))
     end if

     damp(n-i+1)=damp(i)
  end do
  damp(0)=damp(1)
  damp(n+1)=damp(n)

  do i=1,npml+1
     xb=xmax+0.5*dz-float(i-1)*dz
     xb1=float(i-1)*dz-0.5*dz

     betadb=xb1*betam/xpml

     if (itd.eq.1) then
        epsb=a*exp(b*log(((xb-(xmax-xpml))/xpml)+1e-10))
        alphab=(alphad-1.)*exp(b*log(((xb-(xmax-xpml))/xpml)+1e-10))+1.
     else
        epsb=a*(1.-cos((xb-(xmax-xpml))*pi2/xpml))
        alphab=(alphad-1.)*(1.-cos((xb-(xmax-xpml))*pi2/xpml))+1.
     end if

     if (abs(betamax).lt.1e-3) then
	betad=0.
	betadb=0.
     end if

     if (idb.eq.1) then
        dampb(n-i+1)=1./(alphab+ci*epsb/(omegac+ci*betadb))
     end if

  end do

  dampb(0)=dampb(1)
  dampb(n+1)=dampb(n)

 ! write(*,*) "******************"
 ! do i=0,n+1
 !    write(*,*) "damp(i)",i,damp(i)
 !    write(*,*) "dampb(i)",i,dampb(i)
!  end do

  return
end subroutine subdamp
