!Copyright 2008-2011 SEISCOPE project, All rights reserved.
!Copyright 2013-2021 SEISCOPEII project, All rights reserved.

!Multiplication by R_s^{T}R_s where R_s is the restriction operator associated with source s
!Not efficient in memory cost: RHS is duplicated, something more clever will have to be done

subroutine sub_restriction(pbdir,acqui)
  
  implicit none
#include "pbdirect.h"
  include "acqui.h"
  
  !PBDIRECT
  TYPE (pbdirect) :: pbdir
  !ACQUI
  TYPE(acquisition) :: acqui
  
  integer :: n1e,n2e,nsrc  
  complex,allocatable,dimension(:) :: RHS_temp
  integer :: isrc,irec
  real    :: xr,zr,h
  integer :: ixr,izr,k,irec_true,npml
  
  !Problem size
  n1e=pbdir%n1e
  n2e=pbdir%n2e
  npml=pbdir%npml  
  h=pbdir%h
  nsrc=acqui%nsrc

  !Duplicate RHS  
  allocate(RHS_temp(n1e*n2e*nsrc*pbdir%nb_wfld))
  do k=1,n1e*n2e*nsrc*pbdir%nb_wfld
     RHS_temp(k)=pbdir%id%RHS(k)
  enddo
  
  !Set RHS to 0
  pbdir%id%RHS(:)=0.
    
  !Loop on sources
  do isrc = 1, nsrc
     !Loop on receivers
     do irec = 1, acqui%nbr_rec(isrc)
        irec_true = acqui%associated_rec(isrc,irec)
        !Position of the receiver in the grid
        xr=acqui%coord_rec(irec_true,1)
        zr=acqui%coord_rec(irec_true,2)
        ixr=NINT(xr/h)+1+npml
        izr=NINT(zr/h)+1+npml
        k=(isrc-1)*(n1e*n2e*pbdir%nb_wfld)+(ixr-1)*n1e+izr
        pbdir%id%RHS(k)=RHS_temp(k)                                
     enddo
  enddo
  deallocate(RHS_temp)
  
end subroutine sub_restriction

!In case the HICKS interpolation is chosen, the restriction operator is changed.
!This subroutined accounts for this change.
subroutine sub_restriction_HICKS(pbdir,acqui)
  
  implicit none
#include "pbdirect.h"
  include "acqui.h"
  
  !IN
  TYPE(acquisition) :: acqui
  !IN/OUT
  TYPE (pbdirect) :: pbdir
  !Local variables  
  integer :: n1e,n2e,nsrc,npml,i1,i2,isrc,irec,irec_true,kmat
  integer :: ix1,ix2,ialpha1,ialpha2,nalpha_frees
  integer :: irad1min,irad1max,irad2min,irad2max
  real    :: coeff,h
  complex,allocatable,dimension(:) :: RHS_temp

  !Problem size
  n1e=pbdir%n1e
  n2e=pbdir%n2e
  npml=pbdir%npml  
  h=pbdir%h
  nsrc=acqui%nsrc

  !Duplicate RHS  
  allocate(RHS_temp(n1e*n2e*nsrc*pbdir%nb_wfld))
  RHS_temp(:)=pbdir%id%RHS(:)
  
  !Set RHS to 0
  pbdir%id%RHS(:)=0.
    
  !Loop on sources
  do isrc = 1, nsrc
     !Loop on receivers
     do irec = 1, acqui%nbr_rec(isrc)
        irec_true = acqui%associated_rec(isrc,irec)
        !Receiver Hicks parameter
        ix1=pbdir%rec_tab(irec_true,1)
        ix2=pbdir%rec_tab(irec_true,2)        
        ialpha1=pbdir%rec_tab(irec_true,4)
        ialpha2=pbdir%rec_tab(irec_true,5)
        nalpha_frees=pbdir%rec_tab(irec_true,7)
        !Index
        irad1min=MAX(1-ix1,-pbdir%irad1+1)
        irad1max=MIN(pbdir%n1e-ix1,pbdir%irad1)
        irad2min=MAX(1-ix2,-pbdir%irad2+1)
        irad2max=MIN(pbdir%n2e-ix2,pbdir%irad2)
        do i1=irad1min,irad1max
           do i2=irad2min,irad2max
              coeff = pbdir%coef_hicks_monopole(i1,i2,0,ialpha1,ialpha2,1,nalpha_frees)
              if(pbdir%tool.eq.1) then !ISOTROPIC MODELING
                 kmat=(isrc-1)*(pbdir%n1e*pbdir%n2e)+(ix2+i2-1)*pbdir%n1e+ix1+i1           
              elseif(pbdir%tool.eq.2) then !ANISOTROPIC MODELING
                 kmat=(isrc-1)*(2*pbdir%n1e*pbdir%n2e)+(ix2+i2-1)*pbdir%n1e+ix1+i1
              endif
              !MODIF 15/08/13 => SUMMATION FOR ACCOUNTING FOR CLOSE INTERACTIONS
              !pbdir%id%RHS(kmat)=coeff*RHS_temp(kmat)
              pbdir%id%RHS(kmat)=pbdir%id%RHS(kmat)+coeff*RHS_temp(kmat)
           enddo
        enddo
     enddo
  enddo
  deallocate(RHS_temp)
  
end subroutine sub_restriction_HICKS


