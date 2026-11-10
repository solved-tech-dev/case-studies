!Copyright 2008-2011 SEISCOPE project, All rights reserved.
!Copyright 2013-2021 SEISCOPEII project, All rights reserved.

!SUBROUTINE TO FILL THE DENSE RHS OF MUMPS FOR RECEIVERS
SUBROUTINE sub_fill_RHS_RECEIVERS(RHS,n1,n2,npml,h,nsrc,nrec_tot,&
     ntrace,nrecmax,coord_rec,&
     nbr_rec, associated_rec,residual)

  IMPLICIT NONE
#include "common.h"

  INTEGER :: n1,n2,nsrc,nrec_tot,ntrace,nrecmax
  COMPLEX(kind=creal) :: RHS(n1*n2*nsrc)
  INTEGER :: nbr_rec(nsrc),associated_rec(nsrc,nrecmax)
  REAL    :: h
  REAL    :: coord_rec(nrec_tot,2)
  COMPLEX(kind=creal) :: residual(ntrace)

  INTEGER :: isrc,irec
  REAL    :: xr,zr
  INTEGER :: ixr,izr,k,irec_true,kmat,npml


  RHS(:)=0.
  k=1
  DO isrc = 1, nsrc
     DO irec = 1, nbr_rec(isrc)
        irec_true = associated_rec(isrc,irec)
        !position of the receiver in the grid
        xr=coord_rec(irec_true,1)
        zr=coord_rec(irec_true,2)
        ixr=NINT(xr/h)+1+npml
        izr=NINT(zr/h)+1+npml
        kmat=(isrc-1)*(n1*n2)+(ixr-1)*n1+izr        
        RHS(kmat)=CONJG(residual(k))/(h**2)        
        k=k+1
     END DO
  END DO

  !open(10,file='test',access='direct',recl=4*n1*n2)
  !write(10,rec=1) real(RHS(1:n1*n2))
  !close(10)

END SUBROUTINE SUB_FILL_RHS_RECEIVERS


SUBROUTINE sub_fill_RHS_RECEIVERSaniso(RHS,n1,n2,npml,h,nsrc,nrec_tot,&
     ntrace,nrecmax,coord_rec,&
     nbr_rec, associated_rec,residual)

  IMPLICIT NONE
#include "common.h"
  INTEGER :: n1,n2,nsrc,nrec_tot,ntrace,nrecmax
  COMPLEX(kind=creal) :: RHS(n1*n2*nsrc*2)
  INTEGER :: nbr_rec(nsrc),associated_rec(nsrc,nrecmax)
  REAL    :: h
  REAL    :: coord_rec(nrec_tot,2)
  COMPLEX(kind=creal) :: residual(ntrace)

  INTEGER :: isrc,irec
  REAL    :: xr,zr
  INTEGER :: ixr,izr,k,irec_true,kmat,npml


  RHS(:)=0.
  k=1
  DO isrc = 1, nsrc
     DO irec = 1, nbr_rec(isrc)
        irec_true = associated_rec(isrc,irec)
        !position of the source in the grid
        xr=coord_rec(irec_true,1)
        zr=coord_rec(irec_true,2)
        ixr=NINT(xr/h)+1+npml
        izr=NINT(zr/h)+1+npml

        kmat=(isrc-1)*(n1*n2*2)+(ixr-1)*n1+izr
        RHS(kmat)=CONJG(residual(k))/h**2
        k=k+1
     END DO
  END DO


END SUBROUTINE SUB_FILL_RHS_RECEIVERSANISO
