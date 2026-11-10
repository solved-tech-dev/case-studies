!Copyright 2008-2011 SEISCOPE project, All rights reserved.
!Copyright 2013-2021 SEISCOPEII project, All rights reserved.
!SUBROUTINE TO BUILD CONSTANT FOR HICKS INTERPOLATION
subroutine sub_Hicks_constant(pbdir,acqui,coord_src,coord_rec)
  
  implicit none
#include "pbdirect.h"
  include "acqui.h"
  include "hicks.h"

  !IN 
  TYPE(acquisition) :: acqui
  real,dimension(acqui%nsrc,2) :: coord_src
  real,dimension(acqui%nrec_tot,2) :: coord_rec
  
  !IN/OUT
  TYPE(pbdirect) :: pbdir
  !Local variables 
  integer :: is,ir,ix1,ix2,ix3,ialpha_1,ialpha_2,ialpha_3,nalpha_frees
  real :: xs1,xs2,xr1,xr2,bhicks_monopole,bhicks_dipole
  
  !GET TABULATED VALUES
  pbdir%bhicks0(:)=bhicks0_tab(:)
  pbdir%bhicks1(:)=bhicks1_tab(:)
  
  !FIXED CONSTANT
  pbdir%irad=4
  pbdir%nalpha=10  
  pbdir%nsp1=pbdir%npml
  pbdir%nsp2=pbdir%npml  
  pbdir%nsp3=1
  bhicks_monopole=pbdir%bhicks0(pbdir%irad)
  bhicks_dipole=pbdir%bhicks1(pbdir%irad)

  !BUILT HICKS COEFFICIENT !MODIF 03/05/2013
  CALL sub_Hicks(pbdir,2, pbdir%free_surface, pbdir%irad, pbdir%nalpha, &
       bhicks_monopole, bhicks_dipole, &
       pbdir%irad1, pbdir%irad2, pbdir%irad3,&
       pbdir%nalpha1, pbdir%nalpha2, pbdir%nalpha3,&
       pbdir%ishift1,pbdir%ishift2,pbdir%ishift3)
    
  !CONSTRUCT INTERPOLATION TABLES RELATED TO SOURCES AND RECEIVERS
  allocate(pbdir%src_tab(acqui%nsrc,7))  
  do is=1,acqui%nsrc   
     xs1=coord_src(is,2)
     xs2=coord_src(is,1)     
     write(59,*) 'src ',is,' : xs1 :', coord_src(is,2), ' xs2 :',coord_src(is,1)
     CALL sub_sample_Hicks(xs1,xs2,0.,ix1,ix2,ix3,ialpha_1,ialpha_2,ialpha_3,nalpha_frees,&
          pbdir%h,pbdir%irad,pbdir%nalpha,pbdir%npml)
     if(pbdir%free_surface.eq.0) then
        nalpha_frees=0
     endif
     pbdir%src_tab(is,1) = ix1
     pbdir%src_tab(is,2) = ix2
     pbdir%src_tab(is,3) = ix3
     pbdir%src_tab(is,4) = ialpha_1
     pbdir%src_tab(is,5) = ialpha_2
     pbdir%src_tab(is,6) = ialpha_3
     pbdir%src_tab(is,7) = nalpha_frees
  !   write(*,*)pbdir%src_tab(is,:)
  enddo
  allocate(pbdir%rec_tab(acqui%nrec_tot,7))
  do ir=1,acqui%nrec_tot
     write(58,*) 'rec ',ir,' : xr1 :', coord_rec(ir,2), ' xr2 :',coord_rec(ir,1)
  enddo  
  do ir=1,acqui%nrec_tot
     xr1=coord_rec(ir,2) ! careful with the index
     xr2=coord_rec(ir,1)     
     CALL sub_sample_Hicks(xr1,xr2,0.,ix1,ix2,ix3,ialpha_1,ialpha_2,ialpha_3,nalpha_frees,&
          pbdir%h,pbdir%irad,pbdir%nalpha,pbdir%npml)
     if(pbdir%free_surface.eq.0) then
        nalpha_frees=0
     endif
     pbdir%rec_tab(ir,1) = ix1
     pbdir%rec_tab(ir,2) = ix2
     pbdir%rec_tab(ir,3) = ix3
     pbdir%rec_tab(ir,4) = ialpha_1
     pbdir%rec_tab(ir,5) = ialpha_2
     pbdir%rec_tab(ir,6) = ialpha_3
     pbdir%rec_tab(ir,7) = nalpha_frees
  enddo
  
end subroutine sub_Hicks_constant

