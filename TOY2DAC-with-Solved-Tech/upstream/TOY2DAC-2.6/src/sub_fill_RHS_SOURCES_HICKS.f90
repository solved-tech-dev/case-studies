!Copyright 2008-2011 SEISCOPE project, All rights reserved.
!Copyright 2013-2021 SEISCOPEII project, All rights reserved.
!SUBROUTINE TO FILL THE DENSE RHS OF MUMPS FOR SOURCES USING HICKS INTERPOLATION
subroutine sub_fill_RHS_SOURCES_HICKS(pbdir,nsrc,RHS)
  
  implicit none
#include "common.h"
#include "pbdirect.h"
  
  !IN
  integer :: nsrc
  TYPE (pbdirect) :: pbdir
  !IN/OUT
  complex(kind=creal),dimension(pbdir%n1e*pbdir%n2e*nsrc) :: RHS  
  !Local variables
  integer :: isrc,ix1,ix2,ialpha1,ialpha2,nalpha_frees
  integer :: irad1min,irad1max,irad2min,irad2max
  integer :: i1,i2,k
  real :: coeff
  
  
  if (pbdir%itypeso.eq.0) then			!explosion
  
  RHS(:)=0.  
  do isrc = 1, nsrc
     !Source Parameter     
     ix1=pbdir%src_tab(isrc,1) 
     ix2=pbdir%src_tab(isrc,2)
     ialpha1=pbdir%src_tab(isrc,4)
     ialpha2=pbdir%src_tab(isrc,5)
     nalpha_frees=pbdir%src_tab(isrc,7)     
     !write(99,*) pbdir%src_tab(isrc,:)
     !Index     
     irad1min=MAX(1-ix1,-pbdir%irad1+1)
     irad1max=MIN(pbdir%n1e-ix1,pbdir%irad1)
     irad2min=MAX(1-ix2,-pbdir%irad2+1)
     irad2max=MIN(pbdir%n2e-ix2,pbdir%irad2)
     !write(99,*) 'irad1min : ', irad1min
     !write(99,*) 'irad1max : ', irad1max
     !write(99,*) 'irad2min : ', irad2min
     !write(99,*) 'irad2max : ', irad2max
     do i1=irad1min,irad1max
        do i2=irad2min,irad2max
           coeff = pbdir%coef_hicks_monopole(i1,i2,0,ialpha1,ialpha2,1,nalpha_frees)
           !if(coeff.ne.0.) then
           !   write(94,*) 'i1+ix1,i2+ix2,coeff : ',i1+ix1,i2+ix2,coeff
           !endif
           if(pbdir%free_surface==0 .OR. (ix1+i1).ne.(pbdir%npml+1)) then !DO NOT EXCITE THE FREE SURFACE
              if(coeff.ne.0.) then
                 k=(isrc-1)*(pbdir%n1e*pbdir%n2e)+(ix2+i2-1)*pbdir%n1e+ix1+i1
                 !RHS(k)=coeff/(pbdir%h**2)                                   
                 RHS(k)=RHS(k)+coeff/(pbdir%h**2)                                   
              endif
           endif
        enddo
     enddo
  enddo
  
  else						!vertical force
  
  RHS(:)=0.  
  do isrc = 1, nsrc
     !Source Parameter     
     ix1=pbdir%src_tab(isrc,1) 
     ix2=pbdir%src_tab(isrc,2)
     ialpha1=pbdir%src_tab(isrc,4)
     ialpha2=pbdir%src_tab(isrc,5)
     nalpha_frees=pbdir%src_tab(isrc,7)     
     !write(99,*) pbdir%src_tab(isrc,:)
     !Index     
     irad1min=MAX(1-ix1,-pbdir%irad1+1)
     irad1max=MIN(pbdir%n1e-ix1,pbdir%irad1)
     irad2min=MAX(1-ix2,-pbdir%irad2+1)
     irad2max=MIN(pbdir%n2e-ix2,pbdir%irad2)
     !write(99,*) 'irad1min : ', irad1min
     !write(99,*) 'irad1max : ', irad1max
     !write(99,*) 'irad2min : ', irad2min
     !write(99,*) 'irad2max : ', irad2max
     do i1=irad1min,irad1max
        do i2=irad2min,irad2max
!write(*,*)i1,i2,0,ialpha1,ialpha2,1,nalpha_frees,pbdir%coef_hicks_dipole_z(i1,i2,0,ialpha1,ialpha2,1,nalpha_frees)
           coeff = pbdir%coef_hicks_dipole_z(i1,i2,0,ialpha1,ialpha2,1,nalpha_frees)
           !if(coeff.ne.0.) then
           !   write(94,*) 'i1+ix1,i2+ix2,coeff : ',i1+ix1,i2+ix2,coeff
           !endif
           if(pbdir%free_surface==0 .OR. (ix1+i1).ne.(pbdir%npml+1)) then !DO NOT EXCITE THE FREE SURFACE
              if(coeff.ne.0.) then
                 k=(isrc-1)*(pbdir%n1e*pbdir%n2e)+(ix2+i2-1)*pbdir%n1e+ix1+i1
                 !RHS(k)=coeff/(pbdir%h**2)                                   
                 RHS(k)=RHS(k)+coeff/(pbdir%h**2)                                   
              endif
!else
!write(*,*)'toto'
           endif
        enddo
     enddo
  enddo
  
  end if
  
  !open(10,file='test',access='direct',recl=4*pbdir%n1e*pbdir%n2e)
  !write(10,rec=1) real(RHS(1:pbdir%n1e*pbdir%n2e))
  !close(10)

end subroutine sub_fill_RHS_SOURCES_HICKS


!ANISOTROPIC VERSION
subroutine sub_fill_RHS_SOURCESaniso_HICKS(pbdir,nsrc,RHS)
  
  implicit none
#include "common.h"
#include "pbdirect.h"
  
  !IN
  integer :: nsrc  
  TYPE (pbdirect) :: pbdir
  !IN/OUT
  complex(kind=creal),dimension(2*pbdir%n1e*pbdir%n2e*nsrc) :: RHS  
  !Local variables
  integer :: isrc,ix1,ix2,ialpha1,ialpha2,nalpha_frees
  integer :: irad1min,irad1max,irad2min,irad2max
  integer :: i1,i2,k
  real :: coeff
  
  if (pbdir%itypeso.eq.0) then   	!explosion
  
  RHS(:)=0.  
  do isrc = 1, nsrc
     !Source Parameter
     ix1=pbdir%src_tab(isrc,1)
     ix2=pbdir%src_tab(isrc,2)         
     ialpha1=pbdir%src_tab(isrc,4)
     ialpha2=pbdir%src_tab(isrc,5)
     nalpha_frees=pbdir%src_tab(isrc,7)     
     !Index     
     irad1min=MAX(1-ix1,-pbdir%irad1+1)
     irad1max=MIN(pbdir%n1e-ix1,pbdir%irad1)
     irad2min=MAX(1-ix2,-pbdir%irad2+1)
     irad2max=MIN(pbdir%n2e-ix2,pbdir%irad2)
     do i1=irad1min,irad1max
        do i2=irad2min,irad2max
           coeff = pbdir%coef_hicks_monopole(i1,i2,0,ialpha1,ialpha2,1,nalpha_frees)
           if(pbdir%free_surface==0 .OR.(ix1+i1).ne.(pbdir%npml+1)) then !DO NOT EXCITE THE FREE SURFACE              
              if(coeff.ne.0.) then
                 k=(isrc-1)*(2*pbdir%n1e*pbdir%n2e)+(ix2+i2-1)*pbdir%n1e+ix1+i1
                 !RHS(k)=coeff/(pbdir%h**2)                            
                 !CAREFUL 14/08/13=> account for close sources and interaction points
                 RHS(k)=RHS(k)+coeff/(pbdir%h**2)                            
              endif
           endif
        enddo
     enddo
  enddo
  
  else					!vertical force
  
  RHS(:)=0.  
  do isrc = 1, nsrc
     !Source Parameter
     ix1=pbdir%src_tab(isrc,1)
     ix2=pbdir%src_tab(isrc,2)         
     ialpha1=pbdir%src_tab(isrc,4)
     ialpha2=pbdir%src_tab(isrc,5)
     nalpha_frees=pbdir%src_tab(isrc,7)     
     !Index     
     irad1min=MAX(1-ix1,-pbdir%irad1+1)
     irad1max=MIN(pbdir%n1e-ix1,pbdir%irad1)
     irad2min=MAX(1-ix2,-pbdir%irad2+1)
     irad2max=MIN(pbdir%n2e-ix2,pbdir%irad2)
     do i1=irad1min,irad1max
        do i2=irad2min,irad2max
           coeff = pbdir%coef_hicks_dipole_z(i1,i2,0,ialpha1,ialpha2,1,nalpha_frees)
           if(pbdir%free_surface==0 .OR.(ix1+i1).ne.(pbdir%npml+1)) then !DO NOT EXCITE THE FREE SURFACE              
              if(coeff.ne.0.) then
                 k=(isrc-1)*(2*pbdir%n1e*pbdir%n2e)+(ix2+i2-1)*pbdir%n1e+ix1+i1
                 !RHS(k)=coeff/(pbdir%h**2)                            
                 !CAREFUL 14/08/13=> account for close sources and interaction points
                 RHS(k)=RHS(k)+coeff/(pbdir%h**2)                            
              endif
           endif
        enddo
     enddo
  enddo
  
  end if
  
end subroutine sub_fill_RHS_SOURCESaniso_HICKS

