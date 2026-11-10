!Copyright 2008-2011 SEISCOPE project, All rights reserved.
!Copyright 2013-2021 SEISCOPEII project, All rights reserved.
!
! subroutine subrho2b.f
! compute buoyancy from density
! if density is 0, buoyancy is forced to 1e10

 subroutine sub_rho2b(rho,n1,n2,bu)
 IMPLICIT NONE
 integer n1,n2,i2,i1
 real rho(n1,n2),bu(0:n1+1,0:n2+1)


 do i2=1,n2
 do i1=1,n1
 if (rho(i1,i2).ne.0.) then
 bu(i1,i2)=1./rho(i1,i2)
 else
 bu(i1,i2)=1e10
 end if
 end do
 end do
 do i2=1,n2
 bu(0,i2)=bu(1,i2)
 bu(n1+1,i2)=bu(n1,i2)
 end do
 do i1=1,n1
 bu(i1,0)=bu(i1,1)
 bu(i1,n2+1)=bu(i1,n2)
 end do
 bu(0,0)=bu(1,1)
 bu(0,n2+1)=bu(1,n2)
 bu(n1+1,0)=bu(n1,1)
 bu(n1+1,n2+1)=bu(n1,n2)
 return
 end
