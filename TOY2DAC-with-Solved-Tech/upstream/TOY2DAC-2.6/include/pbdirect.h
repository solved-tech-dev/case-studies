#ifdef DOUBLEPRECISION
include "zmumps_struc.h"
#else
include "cmumps_struc.h"
#endif

TYPE pbdirect
                                                                      
SEQUENCE

 !###MUMPS
#ifdef DOUBLEPRECISION
 TYPE (zmumps_struc) :: id
#else
TYPE (cmumps_struc) :: id
#endif

 INTEGER ::type_RHS
 !###PARAMETERS
 INTEGER :: mode,tool
 !###MODEL
 REAL,DIMENSION(:,:),ALLOCATABLE :: vp,rho,qp,epsilon,theta,delta
 !### TYPE OF SOURCE AND RECEIVERS
 INTEGER  :: itypeso,ityperec
 !###FREQUENCY
 INTEGER  :: nw
 !###IMAGINARY PART OF FREQUENCY (LAPLACE CONSTANT)
 REAL :: slaplace
 REAL,DIMENSION(:),ALLOCATABLE :: freqw,freq
 !###GRID
 INTEGER :: n1,n2,n1e,n2e
 REAL    :: h
 !###PML
 REAL    :: apml
 INTEGER :: npml
 !###FREE SURFACE
 INTEGER :: free_surface
 !###NUMBER OF WAVEFIELDS
 INTEGER :: nb_wfld
 !###HICKS INTERPOLATION
 INTEGER :: hicks 
 INTEGER :: irad,nalpha,irad1,irad2,irad3,nalpha1,nalpha2,nalpha3,ishift1,ishift2,ishift3
 INTEGER :: nsp1,nsp2,nsp3
 INTEGER,dimension(:,:),allocatable :: src_tab,rec_tab
 REAL,DIMENSION(10) :: bhicks0,bhicks1
 !REAL,DIMENSION(:,:,:,:,:,:,:),pointer :: coef_hicks_monopole,coef_hicks_dipole
 REAL,DIMENSION(:,:,:,:,:,:,:),allocatable :: coef_hicks_monopole,coef_hicks_dipole_z,coef_hicks_dipole_x

  !fake sizes
 INTEGER :: n1e_s,n2e_s,nb_wfld_s,nw_s  
  
  !MULTICOMPONENT
  integer :: multicomponent

END TYPE pbdirect
