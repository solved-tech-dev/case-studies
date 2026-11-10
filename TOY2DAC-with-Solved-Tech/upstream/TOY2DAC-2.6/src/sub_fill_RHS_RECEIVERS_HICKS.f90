!Copyright 2008-2011 SEISCOPE project, All rights reserved.
!Copyright 2013-2021 SEISCOPEII project, All rights reserved.
!SUBROUTINE TO FILL THE DENSE RHS OF MUMPS FOR RECEIVERS USING HICKS INTERPOLATION
!THE RHS IS FILLED WITH THE RESIDUALS 
SUBROUTINE sub_fill_RHS_RECEIVERS_HICKS(pbdir,nsrc,ntrace,nrecmax,&
     nbr_rec,associated_rec,residual,RHS)
  
  implicit none
#include "common.h"
#include "pbdirect.h"
  
  !IN
  integer :: nsrc,ntrace,nrecmax
  integer,dimension(nsrc) :: nbr_rec
  integer,dimension(nsrc,nrecmax) :: associated_rec  
  complex(kind=creal),dimension(ntrace) :: residual
  TYPE(pbdirect) :: pbdir
  !IN/OUT
  complex(kind=creal),dimension(pbdir%n1e*pbdir%n2e*nsrc) :: RHS
  !Local variables
  integer :: i1,i2,k,kmat,isrc,irec,irec_true
  integer :: ix1,ix2,ialpha1,ialpha2,nalpha_frees
  integer :: irad1min,irad1max,irad2min,irad2max
  real :: coeff
  
  if (pbdir%ityperec.eq.0) then					!hydrophones
    
  RHS(:)=0.
  k=1
  do isrc = 1, nsrc
     do irec = 1, nbr_rec(isrc)
        irec_true = associated_rec(isrc,irec)
        !Receiver-Source Parameter
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
              if(pbdir%free_surface==0 .OR.(ix1+i1).ne.(pbdir%npml+1)) then !MODIF 17/04/13                   
                 if(coeff.ne.0.) then                                       
                    !CAREFUL 14/08/13=> account for close receivers and interaction points !SUMMATION
                    kmat=(isrc-1)*(pbdir%n1e*pbdir%n2e)+(ix2+i2-1)*pbdir%n1e+ix1+i1
                    RHS(kmat)=RHS(kmat)+coeff*CONJG(residual(k))/(pbdir%h**2)                                        
                 endif
              endif
           enddo
        enddo
        k=k+1
     enddo
  enddo
  
  else								!vertical geophones
  
  RHS(:)=0.
  k=1
  do isrc = 1, nsrc
     do irec = 1, nbr_rec(isrc)
        irec_true = associated_rec(isrc,irec)
        !Receiver-Source Parameter
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
              coeff = pbdir%coef_hicks_dipole_z(i1,i2,0,ialpha1,ialpha2,1,nalpha_frees)              
              if(pbdir%free_surface==0 .OR.(ix1+i1).ne.(pbdir%npml+1)) then !MODIF 17/04/13                   
                 if(coeff.ne.0.) then                                       
                    !CAREFUL 14/08/13=> account for close receivers and interaction points !SUMMATION
                    kmat=(isrc-1)*(pbdir%n1e*pbdir%n2e)+(ix2+i2-1)*pbdir%n1e+ix1+i1
                    RHS(kmat)=RHS(kmat)+coeff*CONJG(residual(k))/(pbdir%h**2)                                        
                 endif
              endif
           enddo
        enddo
        k=k+1
     enddo
  enddo
  
  end if
  
  
  !open(10,file='test',access='direct',recl=4*pbdir%n1e*pbdir%n2e)
  !write(10,rec=1) real(RHS(1:pbdir%n1e*pbdir%n2e))
  !close(10)
  
end subroutine sub_fill_RHS_RECEIVERS_HICKS


!WITH ANISOTROPY
subroutine sub_fill_RHS_RECEIVERSaniso_HICKS(pbdir,nsrc,ntrace,nrecmax,&
     nbr_rec,associated_rec,residual,RHS)
  
  implicit none
#include "common.h"  
#include "pbdirect.h"
  !IN
  integer :: nsrc,ntrace,nrecmax
  integer,dimension(nsrc) :: nbr_rec
  integer,dimension(nsrc,nrecmax) :: associated_rec
  complex(kind=creal),dimension(ntrace) :: residual
  type(pbdirect) :: pbdir
  !IN/OUT
  complex(kind=creal),dimension(2*pbdir%n1e*pbdir%n2e*nsrc) :: RHS
  !Local variables
  integer :: i1,i2,k,kmat,isrc,irec,irec_true
  integer :: ix1,ix2,ialpha1,ialpha2,nalpha_frees
  integer :: irad1min,irad1max,irad2min,irad2max
  real :: coeff
  
  if (pbdir%ityperec.eq.0) then							!hydrophones

  RHS(:)=0.
  k=1
  do isrc = 1, nsrc
     do irec = 1, nbr_rec(isrc)
        irec_true = associated_rec(isrc,irec)
        !Receiver-Source Parameter
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
              if(pbdir%free_surface==0 .OR.(ix1+i1).ne.(pbdir%npml+1)) then !MODIF 17/04/13
                 if(coeff.ne.0) then
                    !CAREFUL 14/08/13=> account for close receivers  and interaction points
                    kmat=(isrc-1)*(2*pbdir%n1e*pbdir%n2e)+(ix2+i2-1)*pbdir%n1e+ix1+i1
                    RHS(kmat)=RHS(kmat)+coeff*CONJG(residual(k))/(pbdir%h**2) 
                 endif
              endif
           enddo
        enddo
        k=k+1
     enddo
  enddo
  
  else										!vertical geophones
  
  RHS(:)=0.
  k=1
  do isrc = 1, nsrc
     do irec = 1, nbr_rec(isrc)
        irec_true = associated_rec(isrc,irec)
        !Receiver-Source Parameter
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
              coeff = pbdir%coef_hicks_dipole_z(i1,i2,0,ialpha1,ialpha2,1,nalpha_frees)
              if(pbdir%free_surface==0 .OR.(ix1+i1).ne.(pbdir%npml+1)) then !MODIF 17/04/13
                 if(coeff.ne.0) then
                    !CAREFUL 14/08/13=> account for close receivers  and interaction points
                    kmat=(isrc-1)*(2*pbdir%n1e*pbdir%n2e)+(ix2+i2-1)*pbdir%n1e+ix1+i1
                    RHS(kmat)=RHS(kmat)+coeff*CONJG(residual(k))/(pbdir%h**2) 
                 endif
              endif
           enddo
        enddo
        k=k+1
     enddo
  enddo
  
  end if
  
    
END SUBROUTINE sub_fill_RHS_RECEIVERSaniso_HICKS

