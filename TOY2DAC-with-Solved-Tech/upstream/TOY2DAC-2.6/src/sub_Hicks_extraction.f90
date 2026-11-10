!Copyright 2008-2011 SEISCOPE project, All rights reserved.
!Copyright 2013-2021 SEISCOPEII project, All rights reserved.
!HICKS EXTRACTION
subroutine Hicks_extraction(pbdir,acqui,isrc,irec_true,data_cal_k)
  
  implicit none
#include "common.h"
#include "pbdirect.h"
  include "acqui.h"
  
  !IN
  integer :: isrc,irec_true
  type(acquisition) :: acqui
  type(pbdirect) :: pbdir
  !IN/OUT
  complex(kind=creal) :: data_cal_k
  !Local variables
  integer :: i1,i2,ix1,ix2,ialpha1,ialpha2,nalpha_frees
  integer :: irad1min,irad1max,irad2min,irad2max
  real :: coeff
  
  !Index of the receiver irect_true
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

  data_cal_k=CMPLX(0.,0.)
  do i1=irad1min,irad1max
     do i2=irad2min,irad2max
        coeff = pbdir%coef_hicks_monopole(i1,i2,0,ialpha1,ialpha2,1,nalpha_frees)
        if(pbdir%tool.eq.1) then
           !DO NOT TAKE POINT ON THE !FREE SURFACE
           if((ix1+i1).ne.(pbdir%npml+1)) then               
              data_cal_k=data_cal_k+coeff*&
                   pbdir%id%RHS((isrc-1)*(pbdir%n1e*pbdir%n2e)+(ix2+i2-1)*pbdir%n1e+ix1+i1)
           endif
        elseif(pbdir%tool.eq.2) then
           !DO NOT TAKE POINT ON THE !FREE SURFACE
           if((ix1+i1).ne.(pbdir%npml+1)) then               
              data_cal_k=data_cal_k+coeff*pbdir%id%RHS((isrc-1)*&
                   (2*pbdir%n1e*pbdir%n2e)+(ix2+i2-1)*pbdir%n1e+ix1+i1)
           endif
        endif
     enddo
  enddo
    
end subroutine Hicks_extraction



