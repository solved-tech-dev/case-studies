!Copyright 2008-2011 SEISCOPE project, All rights reserved.
!Copyright 2013-2015 SEISCOPEII project, All rights reserved.
!=====================================================================
! Extension of the model defined by v(nz,nx) into ve(nz+2*nsp,nz+2*nsp)
! by copying values at edges into PML zones
!
!=====================================================================
	subroutine sub_modext(v,nx,nz,ve,nspx,nspz)
       
       IMPLICIT NONE
       
       integer nz,nx,nspz,nspx,ix,iz
       real v(nz,nx),ve(nz+2*nspz,nx+2*nspx)
      
!
! kernel
!
       do ix=1,nx
       do iz=1,nz
       ve(iz+nspz,ix+nspx)=v(iz,ix)
       end do
       end do
!
! right- and left- sides
!
       do ix=1,nspx
       do iz=1,nz
       ve(iz+nspz,ix)=v(iz,1)
       ve(iz+nspz,ix+nspx+nx)=v(iz,nx)
       end do
       end do
!
! upper- and  lower- side
!

       do ix=1,nx
       do iz=1,nspz
       ve(iz,ix+nspx)=v(1,ix)
       ve(nz+iz+nspz,ix+nspx)=v(nz,ix)
       end do
       end do
!
! corners
!
       do ix=1,nspx
       do iz=1,nspz
       ve(iz,ix)=v(1,1)
       ve(iz,ix+nspx+nx)=v(1,nx)
       ve(iz+nz+nspz,ix)=v(nz,1)
       ve(iz+nz+nspz,ix+nspx+nx)=v(nz,nx)
       end do
       end do
!
       return
       end
