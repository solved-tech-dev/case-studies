!Copyright 2008-2011 SEISCOPE project, All rights reserved.
!Copyright 2013-2021 SEISCOPEII project, All rights reserved.
subroutine sub_init_FDFD(pbdir)
  implicit none

#include "pbdirect.h"
  include "mpif.h"
#include "common.h"
  TYPE (pbdirect) :: pbdir
  CHARACTER(LEN=80) :: file_vp,file_rho,file_qp
  INTEGER           :: icntl_7,icntl_14,icntl_23,keep_84

  OPEN(11,file='fdfd_input')
  READ(11,*)pbdir%n1,pbdir%n2
  READ(11,*)pbdir%h
  READ(11,*)file_vp,file_qp,file_rho
  READ(11,*)pbdir%apml,pbdir%npml
  READ(11,*)pbdir%hicks
  READ(11,*)pbdir%free_surface
  READ(11,*)pbdir%itypeso,pbdir%ityperec
  READ(11,*)pbdir%slaplace
  CLOSE(11)
  lmoy=pbdir%h

  IF(pbdir%hicks==0 .AND. (pbdir%itypeso==1 .OR. pbdir%ityperec==1)) THEN
     WRITE(*,*) '**************************************************'
     WRITE(*,*) ' PROBLEM IN SETTING : DIPOLE SOURCE OR RECEIVER '
     WRITE(*,*) ' REQUIRES HICKS interpolator'
     WRITE(*,*) '--------------------------------------------------'
     WRITE(*,*) ' CODE STOPPED NOW'
     WRITE(*,*) '**************************************************'
     CALL MPI_ABORT()
     STOP
  END IF

  IF(pbdir%ityperec==1) THEN
     WRITE(*,*) '**************************************************'
     WRITE(*,*) ' PROBLEM IN SETTING : DIPOLE RECEIVER '
     WRITE(*,*) ' NOT YET IMPLEMENTED'
     WRITE(*,*) '--------------------------------------------------'
     WRITE(*,*) ' CODE STOPPED NOW'
     WRITE(*,*) '**************************************************'
     CALL MPI_ABORT()
     STOP

  END IF


  IF (mype ==0) THEN
     !read fdfd input file
     pbdir%n1e=pbdir%n1+2*pbdir%npml
     pbdir%n2e=pbdir%n2+2*pbdir%npml
     ALLOCATE(pbdir%vp(pbdir%n1,pbdir%n2))
     ALLOCATE(pbdir%rho(pbdir%n1,pbdir%n2))
     ALLOCATE(pbdir%qp(pbdir%n1,pbdir%n2))
     !READ MODELS
     CALL sub_read_models(file_vp,file_qp,file_rho,&
          pbdir%vp,pbdir%qp,pbdir%rho,pbdir%n1,pbdir%n2) 
  END IF
  ! READ Frequency_Domain frequency files management
  CALL sub_read_freqdom_freq_list(pbdir)

  !INIT MUMPS 
  ! Parallel session 
  pbdir%id%COMM = MPI_COMM_WORLD
  ! Unsymmetric matrix
  pbdir%id%SYM =  0
  ! Host working
  pbdir%id%PAR =  1
  ! Initilization of MUMPS
  pbdir%id%JOB = -1
  ! Initialize MUMPS
#ifdef DOUBLEPRECISION
  CALL zmumps( pbdir%id )
#else
  CALL cmumps( pbdir%id )
#endif


  IF (mype == 0) THEN
     OPEN(16,file='mumps_input')
     READ(16,*)icntl_7
     READ(16,*)icntl_14
     READ(16,*)icntl_23
     READ(16,*)keep_84
     CLOSE(16)
     pbdir%id%ICNTL(20)  =  0     !dense RHS
     pbdir%id%ICNTL(21)  =  0     !0->centralized solution 
     pbdir%id%ICNTL(7)   =  icntl_7     !ordering for analyse
     pbdir%id%ICNTL(14)  =  icntl_14  ! 
     pbdir%id%ICNTL(23)  =  icntl_23
     pbdir%id%KEEP(84)   =  keep_84   ! blocking factor for Multiple RHS
     write(*,*)'INIT OF MUMPS OK'
  END IF


  flag_first_call_to_mumps = 0
end subroutine sub_init_FDFD

subroutine sub_init_FDFDaniso(pbdir)
  implicit none

#include "pbdirect.h"
  include "mpif.h"
#include "common.h"

  TYPE (pbdirect) :: pbdir
  CHARACTER(LEN=80) :: file_vp,file_rho,file_qp,file_epsilon,file_delta,file_theta
  INTEGER           :: icntl_7,icntl_14,icntl_23,keep_84

  OPEN(11,file='fdfd_input')
  READ(11,*)pbdir%n1,pbdir%n2
  READ(11,*)pbdir%h
  READ(11,*)file_vp,file_qp,file_rho,file_epsilon,file_delta,file_theta
  READ(11,*)pbdir%apml,pbdir%npml
  READ(11,*)pbdir%hicks
  READ(11,*)pbdir%free_surface
  READ(11,*)pbdir%itypeso,pbdir%ityperec
  READ(11,*)pbdir%slaplace
  CLOSE(11)

  IF (mype ==0) THEN     
     pbdir%n1e=pbdir%n1+2*pbdir%npml
     pbdir%n2e=pbdir%n2+2*pbdir%npml
     ALLOCATE(pbdir%vp(pbdir%n1,pbdir%n2))
     ALLOCATE(pbdir%rho(pbdir%n1,pbdir%n2))
     ALLOCATE(pbdir%qp(pbdir%n1,pbdir%n2))
     ALLOCATE(pbdir%epsilon(pbdir%n1,pbdir%n2))
     ALLOCATE(pbdir%theta(pbdir%n1,pbdir%n2))
     ALLOCATE(pbdir%delta(pbdir%n1,pbdir%n2))
     !READ MODELS
     CALL sub_read_models(file_vp,file_qp,file_rho,&
          pbdir%vp,pbdir%qp,pbdir%rho,pbdir%n1,pbdir%n2) 
     CALL sub_read_models(file_epsilon,file_theta,file_delta,&
          pbdir%epsilon,pbdir%theta,pbdir%delta,pbdir%n1,pbdir%n2) 
  END IF

  ! READ Frequency_Domain frequency files management
  CALL sub_read_freqdom_freq_list(pbdir)

  !INIT MUMPS 
  ! Parallel session 
  pbdir%id%COMM = MPI_COMM_WORLD
  ! Unsymmetric matrix
  pbdir%id%SYM =  0
  ! Host working
  pbdir%id%PAR =  1
  ! Initilization of MUMPS
  pbdir%id%JOB = -1
  ! Initialize MUMPS
#ifdef DOUBLEPRECISION
  CALL zmumps( pbdir%id )
#else
  CALL cmumps( pbdir%id )
#endif
  
  IF (mype == 0) THEN
     OPEN(16,file='mumps_input')
     READ(16,*)icntl_7
     READ(16,*)icntl_14
     READ(16,*)icntl_23
     READ(16,*)keep_84
     CLOSE(16)
     pbdir%id%ICNTL(20)  =  0     !dense RHS
     pbdir%id%ICNTL(21)  =  0     !0->centralized solution 
     pbdir%id%ICNTL(7)   =  icntl_7     !ordering for analyse
     pbdir%id%ICNTL(14)  =  icntl_14  ! 
     pbdir%id%ICNTL(23)  =  icntl_23
     pbdir%id%KEEP(84)   =  keep_84   ! blocking factor for Multiple RHS
  END IF
  write(*,*)'INIT OF MUMPS OK'
  flag_first_call_to_mumps = 0
  
end subroutine sub_init_FDFDaniso

subroutine sub_pbdirect_FDFD_MATRIX(pbdir,omega)
  ! ACOUSTIC ISOTROPIC MODELING IN THE FREQUENCY DOMAIN

  implicit none
  include "mpif.h"
#include "common.h"
#include "pbdirect.h"

  !PBDIRECT
  TYPE (pbdirect) :: pbdir
  REAL :: omega
  COMPLEX :: omegac
  !PML
  COMPLEX, ALLOCATABLE :: damp1p(:),damp2p(:),damp1pb(:),damp2pb(:)
  COMPLEX,DIMENSION(:,:),ALLOCATABLE :: mu
  REAL,DIMENSION(:,:),ALLOCATABLE :: bu,vpe,qpe,rhoe
  !MATRIX
  INTEGER :: ne,nnz
  ! MPI 
  INTEGER :: infompi

  IF (mype ==0) THEN
     !PML
     ALLOCATE( damp1p(0:pbdir%n1e+1))
     ALLOCATE( damp2p(0:pbdir%n2e+1))
     ALLOCATE( damp1pb(0:pbdir%n1e+1))
     ALLOCATE( damp2pb(0:pbdir%n2e+1))
     ALLOCATE( vpe(pbdir%n1e,pbdir%n2e))
     ALLOCATE( rhoe(pbdir%n1e,pbdir%n2e))
     ALLOCATE( qpe(pbdir%n1e,pbdir%n2e))
     CALL sub_modext(pbdir%vp,pbdir%n2,pbdir%n1,vpe,pbdir%npml,pbdir%npml)
     CALL sub_modext(pbdir%rho,pbdir%n2,pbdir%n1,rhoe,pbdir%npml,pbdir%npml)
     CALL sub_modext(pbdir%qp,pbdir%n2,pbdir%n1,qpe,pbdir%npml,pbdir%npml)
     omegac=CMPLX(omega,pbdir%slaplace)
     !STEP 1.1 : PML
     CALL sub_PML(damp1p,damp1pb,damp2p,damp2pb,&
          pbdir%n1e,pbdir%n2e,pbdir%h,pbdir%apml,pbdir%npml,omega,omegac)
     WRITE(*,*)'sub_PML OK'

     !STEP 1.2 : COMPLEX VELOCITIES
     ALLOCATE(mu(0:pbdir%n1e+1,0:pbdir%n2e+1))
     CALL sub_vcmplx(vpe,qpe,rhoe,mu,pbdir%n1e,pbdir%n2e)
     WRITE(*,*)'sub_vcmplx OK'

     !STEP 1.3 : byuancy
     ALLOCATE(bu(0:pbdir%n1e+1,0:pbdir%n2e+1))
     CALL sub_rho2b(rhoe,pbdir%n1e,pbdir%n2e,bu)
     WRITE(*,*)'subrho2b OK'

     !STEP 1.4 : allocate MUMPS table
     !estimate maximum size of matrix
     ne=9*pbdir%n1e*pbdir%n2e
     
     !allocate
     ALLOCATE(pbdir%id%A(ne))
     ALLOCATE(pbdir%id%IRN(ne))
     ALLOCATE(pbdir%id%JCN(ne))
     CLOSE(99)
     
     !STEP 1.5 : build matrix
     if(pbdir%free_surface.eq.0) then
        CALL  sub_evalmatrix(vpe,bu,mu,damp1p,damp1pb,damp2p,damp2pb,ne, &
             pbdir%id%irn,pbdir%id%jcn,pbdir%id%A,pbdir%n1e,pbdir%n2e,pbdir%h,omegac,nnz)
     elseif(pbdir%free_surface.eq.1) then
        CALL  sub_evalmatrix_fs(pbdir%npml,vpe,bu,mu,damp1p,damp1pb,damp2p,damp2pb,ne, &
             pbdir%id%irn,pbdir%id%jcn,pbdir%id%A,pbdir%n1e,pbdir%n2e,pbdir%h,omegac,nnz)
     endif
     WRITE(*,*)'subevalmatrix OK'
  END IF
  !STAGE 2 : LU DECOMPOSITION
  
  !STEP 2.1 : set up MUMPS
  IF (mype == 0 ) THEN 
     pbdir%id%N = pbdir%n1e*pbdir%n2e
     pbdir%id%Nz = nnz     
  END IF
  IF ( flag_first_call_to_mumps == 0 ) THEN
     IF (mype == 0 )  WRITE(*,*)'ANALYSIS OF THE MATRIX'
     !ANALYSE OF MATRIX
     pbdir%id%JOB = 1
#ifdef DOUBLEPRECISION
  CALL zmumps( pbdir%id )
#else
  CALL cmumps( pbdir%id )
#endif
     flag_first_call_to_mumps=1
  END IF

  if (mype == 0 ) then 
     WRITE(*,*)'FACTORISATION OF THE MATRIX'
  endif
  pbdir%id%JOB = 2
#ifdef DOUBLEPRECISION
  CALL zmumps( pbdir%id )
#else
  CALL cmumps( pbdir%id )
#endif
  
  !DEALLOCATION
  if (mype==0)  then
     deallocate(damp1p,damp1pb,damp2p,damp2pb,bu,mu,vpe,rhoe,qpe,&
          pbdir%id%A,pbdir%id%IRN,pbdir%id%JCN)
  endif
  
end subroutine sub_pbdirect_FDFD_MATRIX

subroutine sub_pbdirect_FDFD_MATRIXaniso(pbdir,omega)
  ! ACOUSTIC ISOTROPIC MODELING IN THE FREQUENCY DOMAIN

  implicit none
  include "mpif.h"
#include "common.h"
#include "pbdirect.h"

  !PBDIRECT
  TYPE (pbdirect) :: pbdir
  REAL :: omega
  COMPLEX :: omegac
  !PML
  COMPLEX, ALLOCATABLE :: damp1p(:),damp2p(:),damp1pb(:),damp2pb(:)
  COMPLEX, ALLOCATABLE :: damp1q(:),damp2q(:),damp1qb(:),damp2qb(:)
  COMPLEX,DIMENSION(:,:),ALLOCATABLE :: mu
  REAL,DIMENSION(:,:),ALLOCATABLE :: bu,vpe,qpe,rhoe,epsilone,deltae,thetae
  !MATRIX
  INTEGER :: ne,nnz
  ! MPI 
  INTEGER :: infompi
  INTEGER :: ianiso

  IF (mype ==0) THEN
     !PML
     ALLOCATE( damp1p(0:pbdir%n1e+1))
     ALLOCATE( damp2p(0:pbdir%n2e+1))
     ALLOCATE( damp1pb(0:pbdir%n1e+1))
     ALLOCATE( damp2pb(0:pbdir%n2e+1))
     ALLOCATE( damp1q(0:pbdir%n1e+1))
     ALLOCATE( damp2q(0:pbdir%n2e+1))
     ALLOCATE( damp1qb(0:pbdir%n1e+1))
     ALLOCATE( damp2qb(0:pbdir%n2e+1))
     ALLOCATE( vpe(pbdir%n1e,pbdir%n2e))
     ALLOCATE( rhoe(pbdir%n1e,pbdir%n2e))
     ALLOCATE( qpe(pbdir%n1e,pbdir%n2e))
     ALLOCATE( epsilone(pbdir%n1e,pbdir%n2e))
     ALLOCATE( deltae(pbdir%n1e,pbdir%n2e))
     ALLOCATE( thetae(pbdir%n1e,pbdir%n2e))
     CALL sub_modext(pbdir%vp,pbdir%n2,pbdir%n1,vpe,pbdir%npml,pbdir%npml)
     CALL sub_modext(pbdir%rho,pbdir%n2,pbdir%n1,rhoe,pbdir%npml,pbdir%npml)
     CALL sub_modext(pbdir%qp,pbdir%n2,pbdir%n1,qpe,pbdir%npml,pbdir%npml)
     CALL sub_modext(pbdir%epsilon,pbdir%n2,pbdir%n1,epsilone,pbdir%npml,pbdir%npml)
     CALL sub_modext(pbdir%delta,pbdir%n2,pbdir%n1,deltae,pbdir%npml,pbdir%npml)
     CALL sub_modext(pbdir%theta,pbdir%n2,pbdir%n1,thetae,pbdir%npml,pbdir%npml)
     omegac=CMPLX(omega,pbdir%slaplace)

     !STEP 1.1 : PML
     CALL sub_PML(damp1p,damp1pb,damp2p,damp2pb,&
          pbdir%n1e,pbdir%n2e,pbdir%h,pbdir%apml,pbdir%npml,omega,omegac)
     CALL sub_PML(damp1q,damp1qb,damp2q,damp2qb,&
          pbdir%n1e,pbdir%n2e,pbdir%h,pbdir%apml,pbdir%npml,omega,omegac)
     WRITE(*,*)'sub_PML OK'

     !STEP 1.2 : COMPLEX VELOCITIES
     ALLOCATE(mu(0:pbdir%n1e+1,0:pbdir%n2e+1))
     CALL sub_vcmplx(vpe,qpe,rhoe,mu,pbdir%n1e,pbdir%n2e)
     WRITE(*,*)'sub_vcmplx OK'

     !STEP 1.3 : byuancy
     ALLOCATE(bu(0:pbdir%n1e+1,0:pbdir%n2e+1))
     CALL sub_rho2b(rhoe,pbdir%n1e,pbdir%n2e,bu)
     WRITE(*,*)'subrho2b OK'

     !STEP 1.4 : allocate MUMPS table, estimate maximum size of matrix
     ne=9*pbdir%n1e*pbdir%n2e*4 !*4 for anisotropy

     !allocate
     ALLOCATE(pbdir%id%A(ne))
     ALLOCATE(pbdir%id%IRN(ne))
     ALLOCATE(pbdir%id%JCN(ne))
     CLOSE(99)

     !STEP 1.5 : build matrix    
     ianiso=5
     if(pbdir%free_surface.eq.0) then
        CALL  sub_evalmatrix_TTI(vpe,bu,mu,deltae,epsilone,thetae,&
             ianiso,damp1p,damp1pb,damp2p,damp2pb,damp1q,damp1qb,damp2q,damp2qb,ne,&
             pbdir%id%irn,pbdir%id%jcn,pbdir%id%A,pbdir%n1e,pbdir%n2e,pbdir%h,omegac,nnz)
     elseif(pbdir%free_surface.eq.1) then
        CALL  sub_evalmatrix_TTI_fs(pbdir%npml,vpe,bu,mu,deltae,epsilone,thetae,&
             ianiso,damp1p,damp1pb,damp2p,damp2pb,damp1q,damp1qb,damp2q,damp2qb,ne,&
             pbdir%id%irn,pbdir%id%jcn,pbdir%id%A,pbdir%n1e,pbdir%n2e,pbdir%h,omegac,nnz)
     endif
     WRITE(*,*)'subevalmatrix OK'
  END IF

  !STAGE 2 : LU DECOMPOSITION
  !STEP 2.1 : set up MUMPS
  IF (mype == 0 ) THEN 
     pbdir%id%N = pbdir%n1e*pbdir%n2e*2
     pbdir%id%Nz = nnz
  END IF
  IF ( flag_first_call_to_mumps == 0 ) THEN
     IF (mype == 0 )  WRITE(*,*)'ANALYSIS OF THE MATRIX'
     !ANALYSE OF MATRIX
     pbdir%id%JOB = 1
#ifdef DOUBLEPRECISION
  CALL zmumps( pbdir%id )
#else
  CALL cmumps( pbdir%id )
#endif
     flag_first_call_to_mumps=1
  END IF
  IF (mype == 0 ) WRITE(*,*)'FACTORISATION OF THE MATRIX'
  pbdir%id%JOB = 2
#ifdef DOUBLEPRECISION
  CALL zmumps( pbdir%id )
#else
  CALL cmumps( pbdir%id )
#endif
  if (mype==0)  then
     deallocate(damp1p,damp1pb,damp2p,damp2pb,damp1q,damp1qb,damp2q,damp2qb,&
          pbdir%id%A,pbdir%id%IRN,pbdir%id%JCN,bu,mu,vpe,deltae,epsilone,thetae)
  endif

end subroutine sub_pbdirect_FDFD_MATRIXaniso

subroutine sub_pbdirect_FDFD_SOLVE_acqui(pbdir,acqui,inv)
  ! ACOUSTIC ISOTROPIC MODELING IN THE FREQUENCY DOMAIN
  
  implicit none
  include "mpif.h"
#include "common.h"
#include "pbdirect.h"
  include "acqui.h"
  include "inversion.h"

  !PBDIRECT
  TYPE (pbdirect) :: pbdir
  !ACQUI
  TYPE(acquisition) :: acqui
  !INVERSION
  TYPE(inversion)   :: inv

  !TIME 
  real(KIND=8) :: t1,t2

  !STAGE 1 : SOLVES

  !STEP 1.1 : fill the RHS
  t1=MPI_WTime()
  if (mype == 0 ) then
     allocate(pbdir%id%RHS(pbdir%n1e*pbdir%n2e*acqui%nsrc))
     if(pbdir%type_RHS.eq.1) then ! SRC SOURCES
        if(pbdir%hicks.eq.0) then
           call sub_fill_RHS_SOURCES(pbdir%id%RHS,pbdir%n1e,pbdir%n2e,pbdir%npml,pbdir%h,acqui%nsrc,acqui%coord_src)
        elseif(pbdir%hicks.eq.1) then
           call sub_fill_RHS_SOURCES_HICKS(pbdir,acqui%nsrc,pbdir%id%RHS)
        endif
     elseif(pbdir%type_RHS.eq.2) then ! REC SOURCES
        if(pbdir%hicks.eq.0) then           
           call sub_fill_RHS_RECEIVERS(pbdir%id%RHS,pbdir%n1e,pbdir%n2e,pbdir%npml,&
                pbdir%h,acqui%nsrc,acqui%nrec_tot,acqui%ntrace,acqui%nrecmax,&
                acqui%coord_rec,acqui%nbr_rec,acqui%associated_rec,inv%residual)
        elseif(pbdir%hicks.eq.1) then
           call sub_fill_RHS_RECEIVERS_HICKS(pbdir,acqui%nsrc,acqui%ntrace,acqui%nrecmax,&
                acqui%nbr_rec,acqui%associated_rec,inv%residual,pbdir%id%RHS)           
        endif
     endif
     !STEP 3.2 : Solve
     pbdir%id%NRHS = acqui%nsrc
     pbdir%id%LRHS = pbdir%id%N
  endif
  t2=MPI_WTime()
  if(mype.eq.0) then
     write(*,*) '1st part SOLVE ACQUI :',t2-t1, ' s' 
  endif
  IF(pbdir%type_RHS == 1) THEN ! SRC SOURCES -> solve Ax=b
     pbdir%id%ICNTL(9)  = 1
  ELSEIF(pbdir%type_RHS == 2) THEN ! REC SOURCES -> Adjoint system -> solve A^tx=b
     pbdir%id%ICNTL(9)  = 0
  END IF
  
  pbdir%id%JOB = 3
  !if(mype.eq.0) then
  !   call CPU_TIME(t1)
  !endif
   t1=MPI_WTime()
#ifdef DOUBLEPRECISION
  CALL zmumps( pbdir%id )
#else
  CALL cmumps( pbdir%id )
#endif
    t2=MPI_WTime()
  if(mype.eq.0) then
     !call CPU_TIME(t2)
     write(*,*) 'TIME FOR MUMPS RESOLUTION :',t2-t1, ' s' 
  endif
  
end subroutine sub_pbdirect_FDFD_SOLVE_acqui

subroutine sub_pbdirect_FDFD_SOLVE_acquianiso(pbdir,acqui,inv)
  ! ACOUSTIC ANISOTROPIC MODELING IN THE FREQUENCY DOMAIN
  
  implicit none
  include "mpif.h"
#include "common.h"
#include "pbdirect.h"
  include "acqui.h"
  include "inversion.h"
  
  !PBDIRECT
  TYPE (pbdirect) :: pbdir
  !ACQUI
  TYPE(acquisition) :: acqui
  !INVERSION
  TYPE(inversion)   :: inv

  !STAGE 1 : SOLVES

  !STEP 1.1 : fill the RHS
  if (mype == 0 ) then
     allocate (pbdir%id%RHS(pbdir%n1e*pbdir%n2e*acqui%nsrc*2))
     if(pbdir%type_RHS.eq.1) then ! SRC SOURCES
        if(pbdir%hicks.eq.0) then
           call sub_fill_RHS_SOURCESaniso(pbdir%id%RHS,pbdir%n1e,pbdir%n2e,&
                pbdir%npml,pbdir%h,acqui%nsrc,acqui%coord_src)
        elseif(pbdir%hicks.eq.1) then           
           call sub_fill_RHS_SOURCESaniso_HICKS(pbdir,acqui%nsrc,pbdir%id%RHS)
        endif
     elseif(pbdir%type_RHS.eq.2) then ! REC SOURCES
        if(pbdir%hicks.eq.0) then
           call sub_fill_RHS_RECEIVERSaniso(pbdir%id%RHS,pbdir%n1e,pbdir%n2e,pbdir%npml,&
                pbdir%h,acqui%nsrc,acqui%nrec_tot,acqui%ntrace,acqui%nrecmax,&
                acqui%coord_rec,acqui%nbr_rec,acqui%associated_rec,inv%residual)
        elseif(pbdir%hicks.eq.1) then
           call sub_fill_RHS_RECEIVERSaniso_HICKS(pbdir,acqui%nsrc,acqui%ntrace,&
                acqui%nrecmax,acqui%nbr_rec,acqui%associated_rec,inv%residual,pbdir%id%RHS)
        endif
     endif
     !STEP 3.2 : Solve
     pbdir%id%NRHS = acqui%nsrc
     pbdir%id%LRHS = pbdir%id%N
  END IF
  IF(pbdir%type_RHS == 1) THEN ! SRC SOURCES -> solve Ax=b
     pbdir%id%ICNTL(9)  = 1
  ELSEIF(pbdir%type_RHS == 2) THEN ! REC SOURCES -> Adjoint system -> solve A^tx=b
     pbdir%id%ICNTL(9)  = 0
  END IF
  
  pbdir%id%JOB = 3
#ifdef DOUBLEPRECISION
  CALL zmumps( pbdir%id )
#else
  CALL cmumps( pbdir%id )
#endif
   
end subroutine sub_pbdirect_FDFD_SOLVE_acquianiso

subroutine sub_pbdirect_FDFD_SOLVE_Jv(pbdir,acqui,inv,inc_wfld,v,omega)
  ! Computation of the matrix vector product J_sv

  implicit none
  include "mpif.h"
#include "common.h"
#include "pbdirect.h"
  include "acqui.h"
  include "inversion.h"

  !PBDIRECT
  TYPE (pbdirect) :: pbdir
  !ACQUI
  TYPE(acquisition) :: acqui
  !INVERSION
  TYPE(inversion) :: inv

  integer :: npar
  complex,dimension(pbdir%n1e_s,pbdir%n2e_s,pbdir%nb_wfld_s,acqui%nsrc) :: inc_wfld
  real,dimension(pbdir%n1*pbdir%n2*inv%npar):: v
  REAL :: omega
  COMPLEX :: omegac
  double precision :: time1,time2

  omegac=CMPLX(omega,pbdir%slaplace)

  !STEP 3.1 : fill the RHS
  IF (mype == 0 ) THEN      
     CALL sub_fill_RHS_Jv(pbdir,acqui,inv,inc_wfld,v,omegac)     
     pbdir%id%NRHS = acqui%nsrc
     pbdir%id%LRHS = pbdir%id%N
  END IF
  !STEP 3.2 : Solve
  pbdir%id%ICNTL(9)  = 1 ! solve Ax=b
  pbdir%id%JOB = 3
#ifdef DOUBLEPRECISION
  CALL zmumps( pbdir%id )
#else
  CALL cmumps( pbdir%id )
#endif

end subroutine sub_pbdirect_FDFD_SOLVE_Jv

subroutine sub_pbdirect_FDFD_WRITE_IO_DATA(pbdir,acqui,iw)
  implicit none
  include "mpif.h"
#include "common.h"
#include "pbdirect.h"
  include "acqui.h"
  
  !PBDIRECT
  TYPE (pbdirect) :: pbdir
  !ACQUI
  TYPE(acquisition) :: acqui
  INTEGER :: iw

  COMPLEX(kind=creal),DIMENSION(:),ALLOCATABLE ::data_cal
  INTEGER :: isrc,irec,irec_true,ixr,izr,k
  REAL :: xrnzr,xr,zr
  complex(kind=creal) :: data_cal_k

  IF (mype ==0) THEN
     ALLOCATE(data_cal(acqui%ntrace))
     data_cal(:)=0.
     k=1
     !loop over source
     DO isrc = 1, acqui%nsrc
        !loop over receivers of the sources
        DO irec = 1, acqui%nbr_rec(isrc)
           irec_true = acqui%associated_rec(isrc,irec)
           if(pbdir%hicks.eq.0) then
              !position of the source in the grid
              xr=acqui%coord_rec(irec_true,1)
              zr=acqui%coord_rec(irec_true,2)
              ixr=NINT(xr/pbdir%h)+1+pbdir%npml
              izr=NINT(zr/pbdir%h)+1+pbdir%npml
              if(pbdir%tool.eq.1)  then 
                 data_cal(k)=pbdir%id%RHS((isrc-1)*(pbdir%n1e*pbdir%n2e)&
                      +(ixr-1)*pbdir%n1e+izr)
              elseif(pbdir%tool.eq.2) then
                 data_cal(k)=pbdir%id%RHS((isrc-1)*(pbdir%n1e*pbdir%n2e*2)&
                      +(ixr-1)*pbdir%n1e+izr)
              endif
           elseif(pbdir%hicks.eq.1) then
              call Hicks_extraction(pbdir,acqui,isrc,irec_true,data_cal_k)
              data_cal(k)=data_cal_k
           endif
           k=k+1
        END DO!irec
     END DO!isrc

     
     
     OPEN(20,file='data_modeling',access='direct',recl=2*creal*acqui%ntrace)
     WRITE(20,rec=iw)data_cal(:)
     CLOSE(20)
     WRITE(*,*) 'DATA MODELING WRITING OK'
          
     if(pbdir%tool == 1) THEN
        OPEN(20,file='wavefield',access='direct',recl=4*pbdir%n1e*pbdir%n2e)
        WRITE(20,rec=iw)SNGL(REAL(pbdir%id%RHS(1:pbdir%n1e*pbdir%n2e)))
     END IF
     if(pbdir%tool == 2) THEN
        OPEN(20,file='wavefield',access='direct',recl=4*pbdir%n1e*pbdir%n2e*2)
        WRITE(20,rec=iw) SNGL(REAL(pbdir%id%RHS(1:pbdir%n1e*pbdir%n2e*2)))
     END IF
     CLOSE(20)
     WRITE(*,*) 'WAVEFIELD MODELING WRITING OK'
     
  END IF
  
end subroutine sub_pbdirect_FDFD_WRITE_IO_DATA

subroutine sub_pbdirect_FDFD_FCOST(pbdir,inv,acqui,iw)
  implicit none
#include "common.h"
#include "pbdirect.h"
  include "inversion.h"
  include "acqui.h"
 

  !PBDIRECT
  TYPE (pbdirect) :: pbdir
  !ACQUI
  TYPE(acquisition) :: acqui
  !INVERSION
  TYPE (inversion) :: inv
  INTEGER ::iw

  COMPLEX(kind=creal),DIMENSION(:),ALLOCATABLE ::data_cal,data_cal_weight
  INTEGER :: isrc,irec,irec_true,ixr,izr,k
  REAL :: xrnzr,xr,zr
  COMPLEX(kind=creal) :: valnum,valden
  complex(kind=creal) :: data_cal_k

  REAL(kind=creal) :: weight

  !extract computed data at receiver position
  ALLOCATE(data_cal(acqui%ntrace))
  ALLOCATE(data_cal_weight(acqui%ntrace))
  data_cal(:)=0.
  valnum=0.
  valden=0.
  k=1
  !loop over source
  DO isrc = 1, acqui%nsrc
     IF(inv%firstgrad==0 .and. (inv%src_estim==2 .OR. inv%src_estim==4))  THEN !REINIT SRC ESTIMATION
        valnum=0.
        valden=0.
     END IF
     !loop over receivers of the sources
     DO irec = 1, acqui%nbr_rec(isrc)
        irec_true = acqui%associated_rec(isrc,irec)
        if(inv%src_estim == 1 .OR. inv%src_estim == 2) THEN
           weight=exp(pbdir%slaplace*acqui%first_arrival_t0(isrc,irec))*acqui%weight_data(k)
        ELSEif(inv%src_estim == 3 .OR. inv%src_estim == 4) THEN
           weight=exp(pbdir%slaplace*acqui%first_arrival_t0(isrc,irec))
        end if
        if(pbdir%hicks.eq.0) then
           !position of the receiver in the grid
           xr=acqui%coord_rec(irec_true,1)
           zr=acqui%coord_rec(irec_true,2)
           ixr=NINT(xr/pbdir%h)+1+pbdir%npml
           izr=NINT(zr/pbdir%h)+1+pbdir%npml
           if(pbdir%tool == 1) then
              data_cal(k)=pbdir%id%RHS((isrc-1)*&
                   (pbdir%n1e*pbdir%n2e)+(ixr-1)*pbdir%n1e+izr)
           elseif(pbdir%tool == 2)  then 
              data_cal(k)=pbdir%id%RHS((isrc-1)*&
                   (pbdir%n1e*pbdir%n2e*2)+(ixr-1)*pbdir%n1e+izr)
           endif
        elseif(pbdir%hicks.eq.1) then
           call Hicks_extraction(pbdir,acqui,isrc,irec_true,data_cal_k)
           data_cal(k)=data_cal_k                     
        endif
        !update on the fly the src estimation terms
        IF(inv%firstgrad==0 .and. inv%src_estim>0 ) THEN
           valnum=valnum+CONJG(data_cal(k)*weight)*inv%obs_data(k,iw)*weight
           valden=valden+CONJG(data_cal(k)*weight)*data_cal(k)*weight
        END IF
        k=k+1
     END DO!irec
     !compute source estimation
     IF(inv%firstgrad==0 .and. ( inv%src_estim==2 .OR. inv%src_estim==4))  THEN
        inv%cc1(isrc,iw)=valnum/valden      
        WRITE(*,*)'source',inv%cc1(isrc,iw),valnum,valden
     END IF
  END DO!isrc

  ! !compute source if global source
  IF(inv%firstgrad==0 .and. ( inv%src_estim==1 .OR.  inv%src_estim==3)) THEN
     inv%cc1(:,iw)=valnum/valden
     WRITE(*,*)'source',inv%cc1(1,iw),valnum,valden
  END IF

  ! correct data from the source + compute fcost
  k=1
  DO isrc = 1, acqui%nsrc
     DO irec = 1, acqui%nbr_rec(isrc)
        weight=dexp(dble(pbdir%slaplace)*dble(acqui%first_arrival_t0(isrc,irec)))*acqui%weight_data(k)
        !WRITE(*,*)weight
        data_cal(k)=data_cal(k)*inv%cc1(isrc,iw)
        data_cal_weight(k)=data_cal(k)*weight
        inv%residual(k)=(inv%obs_data(k,iw)-data_cal(k))*weight/ABS(inv%cc1(isrc,iw))                
        inv%fcost=inv%fcost+0.5*(inv%residual(k)*CONJG(inv%residual(k)))
        inv%residual(k)=inv%residual(k)*weight
        k=k+1
     END DO
  END DO

  OPEN(20,file='data_cal',access='direct',recl=4*acqui%ntrace)
  WRITE(20,rec=iw)SNGL(REAL(data_cal(:)))
  CLOSE(20)
  OPEN(20,file='data_cal_weight',access='direct',recl=4*acqui%ntrace)
  WRITE(20,rec=iw)SNGL(REAL(data_cal_weight(:)))
  CLOSE(20)


  deallocate(data_cal)

end subroutine sub_pbdirect_FDFD_FCOST

subroutine sub_FDFD_gradient(pbdir,inv,nb_wfld,acqui,inc_wfld,adj_wfld,omega,iw)
  implicit none
#include "common.h"
#include "pbdirect.h"
  include "inversion.h"
  include "acqui.h"

  !PBDIRECT
  TYPE (pbdirect) :: pbdir
  !ACQUI
  TYPE(acquisition) :: acqui
  !INVERSION
  TYPE (inversion) :: inv
  REAL :: omega,gradient_vp,gradient_qp
  INTEGER :: isrc,i1,i2,iipar,iw,nb_thread,nb_wfld,i

  COMPLEX,DIMENSION(pbdir%n1e,pbdir%n2e,nb_wfld,acqui%nsrc) :: inc_wfld,adj_wfld
  COMPLEX :: grad_term,xi
  COMPLEX :: omegac
  REAL ::   h2,gx,gz,hx,hz,b00,bp0,bm0,b0p,b0m

 omegac=CMPLX(omega,pbdir%slaplace)
  
  xi=CMPLX(0.,1.)

  DO iipar=1,inv%npar
     DO i2=1,pbdir%n2
        DO i1=inv%ibathy(i2)+inv%ideadzone,pbdir%n1
           IF(inv%invfamily==1) THEN
              IF(inv%invpar(iipar)==1) THEN !VP in family 1
                 !  cas vitesse
                 !term pi = -2 omega^2 / (rho vp^3 (1 - 0.5 i / Q)^2)
                 !term qi = -2 omega^2 / (rho vp^3 (1 - 0.5 i / Q)^2)
                 DO isrc=1,acqui%nsrc
                    DO i=1,nb_wfld
                       inv%gradient(i1,i2,iipar)=inv%gradient(i1,i2,iipar)+ gradient_vp(omegac,pbdir%rho(i1,i2),pbdir%vp(i1,i2),pbdir%qp(i1,i2),inc_wfld(pbdir%npml+i1,pbdir%npml+i2,i,isrc),adj_wfld(pbdir%npml+i1,pbdir%npml+i2,i,isrc),inv%cc1(isrc,iw))     
!write(*,*)    inv%gradient(i1,i2,iipar)              
                    END DO
                 END DO
              ELSEIF(inv%invpar(iipar)==2) THEN !RHO in family 1                 
                 !CALL sub_graddAdmi(i1,i2,iw,iipar,nb_wfld,acqui%nsrc,&
                 !     omega,pbdir,inv,inc_wfld,adj_wfld,grad_term)                 
                 do i=1,nb_wfld
                    call sub_drho5_bis(i1,i2,i,nb_wfld,acqui%nsrc,&
                         omegac,pbdir,inv,inc_wfld,adj_wfld,grad_term)
                    inv%gradient(i1,i2,iipar)=inv%gradient(i1,i2,iipar)+REAL(grad_term)
                 enddo
              ELSEIF(inv%invpar(iipar)==3) THEN !Q in family 1
                 ! cas attenuation
                 !term pi = -2 i omega^2 / (rho vp^2 (1 - 0.5 i / Q)^2 Q^2)
                 !term qi = -2 i omega^2 / (rho vp^2 (1 - 0.5 i / Q)^2 Q^2) 
                 DO isrc=1,acqui%nsrc
                    DO i=1,nb_wfld
                       inv%gradient(i1,i2,iipar)=inv%gradient(i1,i2,iipar)+&
                            gradient_qp(omegac,pbdir%rho(i1,i2),pbdir%vp(i1,i2),&
                            pbdir%qp(i1,i2),inc_wfld(pbdir%npml+i1,pbdir%npml+i2,i,isrc)&
                            ,adj_wfld(pbdir%npml+i1,pbdir%npml+i2,i,isrc),inv%cc1(isrc,iw))
                    END DO
                 END DO
              ELSEIF(inv%invpar(iipar)==4) THEN !EPSILON in family 1
                 !bellow are the g? and h? derived with regard of pbdir%epsilon
                 gx=2.*cos(pbdir%theta(i1,i2))*cos(pbdir%theta(i1,i2))
                 hx=-sin(2.*pbdir%theta(i1,i2))
                 gz=-sin(2.*pbdir%theta(i1,i2))
                 hz=2.*sin(pbdir%theta(i1,i2))*sin(pbdir%theta(i1,i2))

                 IF(i1<pbdir%n1)THEN
                    bp0=0.5*(1./pbdir%rho(i1,i2)+1./pbdir%rho(i1+1,i2))
                 ELSE
                    bp0=1./pbdir%rho(i1,i2)
                 END IF
                 IF(i1>1)THEN                 
                    bm0=0.5*(1./pbdir%rho(i1,i2)+1./pbdir%rho(i1-1,i2))
                 ELSE
                    bm0=1./pbdir%rho(i1,i2)
                 END IF
                 IF(i2<pbdir%n2)THEN
                    b0p=0.5*(1./pbdir%rho(i1,i2)+1./pbdir%rho(i1,i2+1))
                 ELSE
                    b0p=1./pbdir%rho(i1,i2)
                 END IF
                 IF(i2>1)THEN        
                    b0m=0.5*(1./pbdir%rho(i1,i2)+1./pbdir%rho(i1,i2-1))
                 ELSE
                    b0m=1./pbdir%rho(i1,i2)
                 END IF

                 h2=1./(pbdir%h*pbdir%h)

                 DO isrc=1,acqui%nsrc
                    !term sur ligne qi= -2Hq -2 Hp    -> H du papier de Stephane
                    inv%gradient(i1,i2,iipar)=inv%gradient(i1,i2,iipar) + REAL(&
                                !term 00 sur q
                         ((-h2*gx*(b0p+b0m)-h2*hz*(bp0+bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2,2,isrc) +&
                                !term 01 sur q
                         (h2*gx*b0p + 0.25*gz*h2*(bp0-bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2+1,2,isrc) +&
                                !term 0-1 sur q
                         (h2*gx*b0m + 0.25*gz*h2*(-bp0+bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2-1,2,isrc) +&
                                !term -10 sur q
                         (h2*hz*bm0 + 0.25*hx*h2*(-b0p+b0m) )*inc_wfld(pbdir%npml+i1-1,pbdir%npml+i2,2,isrc) +&
                                !term 10 sur q
                         (h2*hz*bp0 + 0.25*hx*h2*(-b0m+b0p) )*inc_wfld(pbdir%npml+i1+1,pbdir%npml+i2,2,isrc) +&
                                !term 00 sur p
                         (-h2*gx*(b0p+b0m)-h2*hz*(bp0+bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2,1,isrc) +&
                                !term 01 sur p
                         (h2*gx*b0p + 0.25*gz*h2*(bp0-bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2+1,1,isrc) +&
                                !term 0-1 sur p
                         (h2*gx*b0m + 0.25*gz*h2*(-bp0+bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2-1,1,isrc) +&
                                !term -10 sur p
                         (h2*hz*bm0 + 0.25*hx*h2*(-b0p+b0m) )*inc_wfld(pbdir%npml+i1-1,pbdir%npml+i2,1,isrc) +&
                                !term 10 sur p
                         (h2*hz*bp0 + 0.25*hx*h2*(-b0m+b0p) )*inc_wfld(pbdir%npml+i1+1,pbdir%npml+i2,1,isrc) ) &
                                !multiply the qi term of the adjoint wavefield
                         * adj_wfld(pbdir%npml+i1,pbdir%npml+i2,2,isrc)*inv%cc1(isrc,iw) )

                 END DO

              ELSEIF(inv%invpar(iipar)==5) THEN !DELTA in family 1

                 !bellow are the g? and h? derived with regard of pbdir%delta
                 gx=-2.*cos(pbdir%theta(i1,i2))*cos(pbdir%theta(i1,i2))
                 hx=sin(2.*pbdir%theta(i1,i2))
                 gz=sin(2.*pbdir%theta(i1,i2))
                 hz=-2.*sin(pbdir%theta(i1,i2))*sin(pbdir%theta(i1,i2))

                 IF(i1<pbdir%n1)THEN
                    bp0=0.5*(1./pbdir%rho(i1,i2)+1./pbdir%rho(i1+1,i2))
                 ELSE
                    bp0=1./pbdir%rho(i1,i2)
                 END IF
                 IF(i1>1)THEN                 
                    bm0=0.5*(1./pbdir%rho(i1,i2)+1./pbdir%rho(i1-1,i2))
                 ELSE
                    bm0=1./pbdir%rho(i1,i2)
                 END IF
                 IF(i2<pbdir%n2)THEN
                    b0p=0.5*(1./pbdir%rho(i1,i2)+1./pbdir%rho(i1,i2+1))
                 ELSE
                    b0p=1./pbdir%rho(i1,i2)
                 END IF
                 IF(i2>1)THEN        
                    b0m=0.5*(1./pbdir%rho(i1,i2)+1./pbdir%rho(i1,i2-1))
                 ELSE
                    b0m=1./pbdir%rho(i1,i2)
                 END IF

                 h2=1./(pbdir%h*pbdir%h)

                 DO isrc=1,acqui%nsrc
                    !term sur ligne pi= -2Hq -2 Hp    -> H du papier de Stephane -> we take -gx -gz -hx and -hz inthe following lines
                    inv%gradient(i1,i2,iipar)=inv%gradient(i1,i2,iipar) + REAL(&
                                !term 00 sur q
                         ((-h2*(-gx)*(b0p+b0m)-h2*(-hz)*(bp0+bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2,2,isrc) +&
                                !term 01 sur q
                         (h2*(-gx)*b0p + 0.25*(-gz)*h2*(bp0-bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2+1,2,isrc) +&
                                !term 0-1 sur q
                         (h2*(-gx)*b0m + 0.25*(-gz)*h2*(-bp0+bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2-1,2,isrc) +&
                                !term -10 sur q
                         (h2*(-hz)*bm0 + 0.25*(-hx)*h2*(-b0p+b0m) )*inc_wfld(pbdir%npml+i1-1,pbdir%npml+i2,2,isrc) +&
                                !term 10 sur q
                         (h2*(-hz)*bp0 + 0.25*(-hx)*h2*(-b0m+b0p) )*inc_wfld(pbdir%npml+i1+1,pbdir%npml+i2,2,isrc) +&
                                !term 00 sur p
                         (-h2*(-gx)*(b0p+b0m)-h2*(-hz)*(bp0+bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2,1,isrc) +&
                                !term 01 sur p
                         (h2*(-gx)*b0p + 0.25*(-gz)*h2*(bp0-bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2+1,1,isrc) +&
                                !term 0-1 sur p
                         (h2*(-gx)*b0m + 0.25*(-gz)*h2*(-bp0+bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2-1,1,isrc) +&
                                !term -10 sur p
                         (h2*(-hz)*bm0 + 0.25*(-hx)*h2*(-b0p+b0m) )*inc_wfld(pbdir%npml+i1-1,pbdir%npml+i2,1,isrc) +&
                                !term 10 sur p
                         (h2*(-hz)*bp0 + 0.25*(-hx)*h2*(-b0m+b0p) )*inc_wfld(pbdir%npml+i1+1,pbdir%npml+i2,1,isrc) ) &
                                !multiply the pi term of the adjoint wavefield
                         * adj_wfld(pbdir%npml+i1,pbdir%npml+i2,1,isrc)*inv%cc1(isrc,iw) )
                    !write(*,*)inv%gradient(i1,i2,iipar)

                    !term sur ligne qi= 2Hq +2 Hp    -> H du papier de Stephane
                    inv%gradient(i1,i2,iipar)=inv%gradient(i1,i2,iipar) + REAL(&
                                !term 00 sur q
                         ((-h2*gx*(b0p+b0m)-h2*hz*(bp0+bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2,2,isrc) +&
                                !term 01 sur q
                         (h2*gx*b0p + 0.25*gz*h2*(bp0-bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2+1,2,isrc) +&
                                !term 0-1 sur q
                         (h2*gx*b0m + 0.25*gz*h2*(-bp0+bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2-1,2,isrc) +&
                                !term -10 sur q
                         (h2*hz*bm0 + 0.25*hx*h2*(-b0p+b0m) )*inc_wfld(pbdir%npml+i1-1,pbdir%npml+i2,2,isrc) +&
                                !term 10 sur q
                         (h2*hz*bp0 + 0.25*hx*h2*(-b0m+b0p) )*inc_wfld(pbdir%npml+i1+1,pbdir%npml+i2,2,isrc) +&
                                !term 00 sur p
                         (-h2*gx*(b0p+b0m)-h2*hz*(bp0+bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2,1,isrc) +&
                                !term 01 sur p
                         (h2*gx*b0p + 0.25*gz*h2*(bp0-bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2+1,1,isrc) +&
                                !term 0-1 sur p
                         (h2*gx*b0m + 0.25*gz*h2*(-bp0+bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2-1,1,isrc) +&
                                !term -10 sur p
                         (h2*hz*bm0 + 0.25*hx*h2*(-b0p+b0m) )*inc_wfld(pbdir%npml+i1-1,pbdir%npml+i2,1,isrc) +&
                                !term 10 sur p
                         (h2*hz*bp0 + 0.25*hx*h2*(-b0m+b0p) )*inc_wfld(pbdir%npml+i1+1,pbdir%npml+i2,1,isrc) ) &
                                !multiply the qi term of the adjoint wavefield
                         * adj_wfld(pbdir%npml+i1,pbdir%npml+i2,2,isrc)*inv%cc1(isrc,iw) )
                    !write(*,*)inv%gradient(i1,i2,iipar)
                 END DO


              ELSEIF(inv%invpar(iipar)==11.OR. inv%invpar(iipar)==10) THEN
                 !LOG(VP) in family 1
                 !multiply by pbdir%vp(i1,i2) because of log : d vp / d log(vp) = vp
                 !  cas vitesse
                 !term pi = -2 vp omega^2 / (rho vp^3 (1 - 0.5 i / Q)^2)
                 !term qi = -2 vp omega^2 / (rho vp^3 (1 - 0.5 i / Q)^2)
                 DO isrc=1,acqui%nsrc
                    DO i=1,nb_wfld
                       inv%gradient(i1,i2,iipar)=inv%gradient(i1,i2,iipar)-REAL(2.*omegac**2/(pbdir%rho(i1,i2)*pbdir%vp(i1,i2)**2 * (1.-0.5*xi/pbdir%qp(i1,i2))**2 )*&
                            inc_wfld(pbdir%npml+i1,pbdir%npml+i2,i,isrc)*adj_wfld(pbdir%npml+i1,pbdir%npml+i2,i,isrc)*inv%cc1(isrc,iw))
                    END DO
                 END DO
                 
              ELSEIF(inv%invpar(iipar)==12) THEN !LOG(RHO) in family 1

                 !gradient from dA/dmi
                 !step1 compute grad_term
                 CALL sub_graddAdmi(i1,i2,iw,iipar,nb_wfld,acqui%nsrc,omegac,pbdir,inv,inc_wfld,adj_wfld,grad_term)

                 !add value
                 inv%gradient(i1,i2,iipar)=inv%gradient(i1,i2,iipar)+REAL(grad_term)
                 
              ELSEIF(inv%invpar(iipar)==31) THEN !LOG(VP/VP_0) in family 1 (VP_0 disappears iun gradient due to log(vp/vp0) )
                 !term pi = -2 vp omega^2 / (rho vp^3 (1 - 0.5 i / Q)^2)
                 !term qi = -2 vp omega^2 / (rho vp^3 (1 - 0.5 i / Q)^2)
                 DO isrc=1,acqui%nsrc
                    DO i=1,nb_wfld
                       inv%gradient(i1,i2,iipar)=inv%gradient(i1,i2,iipar)&
                            -REAL(2.*omegac**2/(pbdir%rho(i1,i2)*pbdir%vp(i1,i2)**2 * (1.-0.5*xi/pbdir%qp(i1,i2))**2 )*&
                            inc_wfld(pbdir%npml+i1,pbdir%npml+i2,i,isrc)*adj_wfld(pbdir%npml+i1,pbdir%npml+i2,i,isrc)*inv%cc1(isrc,iw))
                    END DO
                 END DO
                 
              ELSEIF(inv%invpar(iipar)==32) THEN !LOG(RHO/RHO_0) in family 1

                 !gradient from dA/dmi
                 !step1 compute grad_term
                 CALL sub_graddAdmi(i1,i2,iw,iipar,nb_wfld,acqui%nsrc,omegac,pbdir,inv,inc_wfld,adj_wfld,grad_term)

                 !add value
                 inv%gradient(i1,i2,iipar)=inv%gradient(i1,i2,iipar)+REAL(grad_term)

              ELSEIF(inv%invpar(iipar)==33) THEN !!LOG(QP/QP_0) in family 1

                 ! cas attenuation
                 !term pi = -2 i Q omega^2 / (rho vp^2 (1 - 0.5 i / Q)^2 Q^2)
                 !term qi = -2 i Q omega^2 / (rho vp^2 (1 - 0.5 i / Q)^2 Q^2) 
                 DO isrc=1,acqui%nsrc
                    DO i=1,nb_wfld
                       inv%gradient(i1,i2,iipar)=inv%gradient(i1,i2,iipar)&
                            -REAL(2.*xi*omegac**2/(pbdir%rho(i1,i2)*pbdir%vp(i1,i2)**2*(1.-0.5*xi/pbdir%qp(i1,i2))**2*pbdir%qp(i1,i2))*&
                            inc_wfld(pbdir%npml+i1,pbdir%npml+i2,i,isrc)*adj_wfld(pbdir%npml+i1,pbdir%npml+i2,i,isrc)*inv%cc1(isrc,iw))
                    END DO
                 END DO
              ELSEIF(inv%invpar(iipar)==34) THEN !LOG(1+EPSILON) in family 1 GRAD(LOG(1+EPSILON))=GRAD(EPSILON)*(1+EPSILON) 
                 !bellow are the g? and h? derived with regard of pbdir%epsilon
                 gx=2.*cos(pbdir%theta(i1,i2))*cos(pbdir%theta(i1,i2))
                 hx=-sin(2.*pbdir%theta(i1,i2))
                 gz=-sin(2.*pbdir%theta(i1,i2))
                 hz=2.*sin(pbdir%theta(i1,i2))*sin(pbdir%theta(i1,i2))

                 IF(i1<pbdir%n1)THEN
                    bp0=0.5*(1./pbdir%rho(i1,i2)+1./pbdir%rho(i1+1,i2))
                 ELSE
                    bp0=1./pbdir%rho(i1,i2)
                 END IF
                 IF(i1>1)THEN                 
                    bm0=0.5*(1./pbdir%rho(i1,i2)+1./pbdir%rho(i1-1,i2))
                 ELSE
                    bm0=1./pbdir%rho(i1,i2)
                 END IF
                 IF(i2<pbdir%n2)THEN
                    b0p=0.5*(1./pbdir%rho(i1,i2)+1./pbdir%rho(i1,i2+1))
                 ELSE
                    b0p=1./pbdir%rho(i1,i2)
                 END IF
                 IF(i2>1)THEN        
                    b0m=0.5*(1./pbdir%rho(i1,i2)+1./pbdir%rho(i1,i2-1))
                 ELSE
                    b0m=1./pbdir%rho(i1,i2)
                 END IF

                 h2=1./(pbdir%h*pbdir%h)

                 DO isrc=1,acqui%nsrc
                    !term sur ligne qi= -2Hq -2 Hp    -> H du papier de Stephane
                    inv%gradient(i1,i2,iipar)=inv%gradient(i1,i2,iipar) + REAL( (1+pbdir%epsilon(i1,i2)) * &
                                !term 00 sur q
                         ((-h2*gx*(b0p+b0m)-h2*hz*(bp0+bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2,2,isrc) +&
                                !term 01 sur q
                         (h2*gx*b0p + 0.25*gz*h2*(bp0-bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2+1,2,isrc) +&
                                !term 0-1 sur q
                         (h2*gx*b0m + 0.25*gz*h2*(-bp0+bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2-1,2,isrc) +&
                                !term -10 sur q
                         (h2*hz*bm0 + 0.25*hx*h2*(-b0p+b0m) )*inc_wfld(pbdir%npml+i1-1,pbdir%npml+i2,2,isrc) +&
                                !term 10 sur q
                         (h2*hz*bp0 + 0.25*hx*h2*(-b0m+b0p) )*inc_wfld(pbdir%npml+i1+1,pbdir%npml+i2,2,isrc) +&
                                !term 00 sur p
                         (-h2*gx*(b0p+b0m)-h2*hz*(bp0+bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2,1,isrc) +&
                                !term 01 sur p
                         (h2*gx*b0p + 0.25*gz*h2*(bp0-bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2+1,1,isrc) +&
                                !term 0-1 sur p
                         (h2*gx*b0m + 0.25*gz*h2*(-bp0+bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2-1,1,isrc) +&
                                !term -10 sur p
                         (h2*hz*bm0 + 0.25*hx*h2*(-b0p+b0m) )*inc_wfld(pbdir%npml+i1-1,pbdir%npml+i2,1,isrc) +&
                                !term 10 sur p
                         (h2*hz*bp0 + 0.25*hx*h2*(-b0m+b0p) )*inc_wfld(pbdir%npml+i1+1,pbdir%npml+i2,1,isrc) ) &
                                !multiply the qi term of the adjoint wavefield
                         * adj_wfld(pbdir%npml+i1,pbdir%npml+i2,2,isrc)*inv%cc1(isrc,iw) )

                 END DO

              ELSEIF(inv%invpar(iipar)==35) THEN !!LOG(1+DELTA) in family 1 GRAD(LOG(1+DELTA))=GRAD(DELTA)*(1+DELTA)

                 !bellow are the g? and h? derived with regard of pbdir%delta
                 gx=-2.*cos(pbdir%theta(i1,i2))*cos(pbdir%theta(i1,i2))
                 hx=sin(2.*pbdir%theta(i1,i2))
                 gz=sin(2.*pbdir%theta(i1,i2))
                 hz=-2.*sin(pbdir%theta(i1,i2))*sin(pbdir%theta(i1,i2))

                 IF(i1<pbdir%n1)THEN
                    bp0=0.5*(1./pbdir%rho(i1,i2)+1./pbdir%rho(i1+1,i2))
                 ELSE
                    bp0=1./pbdir%rho(i1,i2)
                 END IF
                 IF(i1>1)THEN                 
                    bm0=0.5*(1./pbdir%rho(i1,i2)+1./pbdir%rho(i1-1,i2))
                 ELSE
                    bm0=1./pbdir%rho(i1,i2)
                 END IF
                 IF(i2<pbdir%n2)THEN
                    b0p=0.5*(1./pbdir%rho(i1,i2)+1./pbdir%rho(i1,i2+1))
                 ELSE
                    b0p=1./pbdir%rho(i1,i2)
                 END IF
                 IF(i2>1)THEN        
                    b0m=0.5*(1./pbdir%rho(i1,i2)+1./pbdir%rho(i1,i2-1))
                 ELSE
                    b0m=1./pbdir%rho(i1,i2)
                 END IF

                 h2=1./(pbdir%h*pbdir%h)

                 DO isrc=1,acqui%nsrc
                    !term sur ligne pi= -2Hq -2 Hp    -> H du papier de Stephane -> we take -gx -gz -hx and -hz inthe following lines
                    inv%gradient(i1,i2,iipar)=inv%gradient(i1,i2,iipar) + REAL((1+pbdir%delta(i1,i2)) * &
                                !term 00 sur q
                         ((-h2*(-gx)*(b0p+b0m)-h2*(-hz)*(bp0+bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2,2,isrc) +&
                                !term 01 sur q
                         (h2*(-gx)*b0p + 0.25*(-gz)*h2*(bp0-bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2+1,2,isrc) +&
                                !term 0-1 sur q
                         (h2*(-gx)*b0m + 0.25*(-gz)*h2*(-bp0+bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2-1,2,isrc) +&
                                !term -10 sur q
                         (h2*(-hz)*bm0 + 0.25*(-hx)*h2*(-b0p+b0m) )*inc_wfld(pbdir%npml+i1-1,pbdir%npml+i2,2,isrc) +&
                                !term 10 sur q
                         (h2*(-hz)*bp0 + 0.25*(-hx)*h2*(-b0m+b0p) )*inc_wfld(pbdir%npml+i1+1,pbdir%npml+i2,2,isrc) +&
                                !term 00 sur p
                         (-h2*(-gx)*(b0p+b0m)-h2*(-hz)*(bp0+bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2,1,isrc) +&
                                !term 01 sur p
                         (h2*(-gx)*b0p + 0.25*(-gz)*h2*(bp0-bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2+1,1,isrc) +&
                                !term 0-1 sur p
                         (h2*(-gx)*b0m + 0.25*(-gz)*h2*(-bp0+bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2-1,1,isrc) +&
                                !term -10 sur p
                         (h2*(-hz)*bm0 + 0.25*(-hx)*h2*(-b0p+b0m) )*inc_wfld(pbdir%npml+i1-1,pbdir%npml+i2,1,isrc) +&
                                !term 10 sur p
                         (h2*(-hz)*bp0 + 0.25*(-hx)*h2*(-b0m+b0p) )*inc_wfld(pbdir%npml+i1+1,pbdir%npml+i2,1,isrc) ) &
                                !multiply the pi term of the adjoint wavefield
                         * adj_wfld(pbdir%npml+i1,pbdir%npml+i2,1,isrc)*inv%cc1(isrc,iw) )
                    !write(*,*)inv%gradient(i1,i2,iipar)

                    !term sur ligne qi= 2Hq +2 Hp    -> H du papier de Stephane
                    inv%gradient(i1,i2,iipar)=inv%gradient(i1,i2,iipar) + REAL((1+pbdir%delta(i1,i2)) *&
                                !term 00 sur q
                         ((-h2*gx*(b0p+b0m)-h2*hz*(bp0+bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2,2,isrc) +&
                                !term 01 sur q
                         (h2*gx*b0p + 0.25*gz*h2*(bp0-bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2+1,2,isrc) +&
                                !term 0-1 sur q
                         (h2*gx*b0m + 0.25*gz*h2*(-bp0+bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2-1,2,isrc) +&
                                !term -10 sur q
                         (h2*hz*bm0 + 0.25*hx*h2*(-b0p+b0m) )*inc_wfld(pbdir%npml+i1-1,pbdir%npml+i2,2,isrc) +&
                                !term 10 sur q
                         (h2*hz*bp0 + 0.25*hx*h2*(-b0m+b0p) )*inc_wfld(pbdir%npml+i1+1,pbdir%npml+i2,2,isrc) +&
                                !term 00 sur p
                         (-h2*gx*(b0p+b0m)-h2*hz*(bp0+bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2,1,isrc) +&
                                !term 01 sur p
                         (h2*gx*b0p + 0.25*gz*h2*(bp0-bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2+1,1,isrc) +&
                                !term 0-1 sur p
                         (h2*gx*b0m + 0.25*gz*h2*(-bp0+bm0) )*inc_wfld(pbdir%npml+i1,pbdir%npml+i2-1,1,isrc) +&
                                !term -10 sur p
                         (h2*hz*bm0 + 0.25*hx*h2*(-b0p+b0m) )*inc_wfld(pbdir%npml+i1-1,pbdir%npml+i2,1,isrc) +&
                                !term 10 sur p
                         (h2*hz*bp0 + 0.25*hx*h2*(-b0m+b0p) )*inc_wfld(pbdir%npml+i1+1,pbdir%npml+i2,1,isrc) ) &
                                !multiply the qi term of the adjoint wavefield
                         * adj_wfld(pbdir%npml+i1,pbdir%npml+i2,2,isrc)*inv%cc1(isrc,iw) )
                    !write(*,*)inv%gradient(i1,i2,iipar)
                 END DO

              end if
           end if
        END DO
        !$OMP END DO 
     END DO

  END DO


  !$OMP END PARALLEL

end subroutine sub_FDFD_gradient


subroutine sub_pbdirect_FDFD_FINISH(pbdir)
  implicit none
  include "mpif.h"
#include "common.h"
#include "pbdirect.h"

  !PBDIRECT
  TYPE (pbdirect) :: pbdir

  IF(mype==0) DEALLOCATE(pbdir%id%RHS)

end subroutine sub_pbdirect_FDFD_FINISH

subroutine sub_finish_FDFD(pbdir)
  implicit none
  include "mpif.h"
#include "common.h"
#include "pbdirect.h"

  !PBDIRECT
  TYPE (pbdirect) :: pbdir
  pbdir%id%JOB = -2
  ! FINALIZE MUMPS
#ifdef DOUBLEPRECISION
  CALL zmumps( pbdir%id )
#else
  CALL cmumps( pbdir%id )
#endif

end subroutine sub_finish_FDFD


subroutine sub_graddAdmi(i1,i2,iw,iipar,nb_wfld,nsrc,omegac,pbdir,inv,inc_wfld,adj_wfld,grad_term)

  implicit none
  include "mpif.h"
#include "common.h"
#include "pbdirect.h"
  include "inversion.h"

  !PBDIRECT
  TYPE (pbdirect) :: pbdir
  !INVERSION
  TYPE (inversion) :: inv

  INTEGER :: i1,i2,nsrc,ianiso,iw,isrc,nb_wfld,iipar
  COMPLEX :: grad_term
  COMPLEX  :: omegac
  COMPLEX,DIMENSION(pbdir%n1e,pbdir%n2e,nb_wfld,nsrc) :: inc_wfld,adj_wfld
  REAL,DIMENSION(-2:2,-2:2) :: b
  REAL,DIMENSION(-1:1,-1:1) :: epsilo,theta,delta
  COMPLEX,DIMENSION(-2:2,-2:2) :: mu1
  complex :: vc,xi,val
  INTEGER :: i1e,i2e,ne,nnz,ne2,nnz2
  INTEGER,DIMENSION(100) :: icn,irn
  COMPLEX,DIMENSION(100) :: mat
  INTEGER,DIMENSION(10) :: ivect2
  COMPLEX,DIMENSION(10) :: vect2
  INTEGER :: mode
  ne=100
  ne2=10

  !omegac=omega
  xi=CMPLX(0.,1.)


  i1e=i1+pbdir%npml
  i2e=i2+pbdir%npml

  if(pbdir%tool == 1) ianiso = 0
  if(pbdir%tool == 2) ianiso = 5

  if(inv%invfamily==1) THEN
     if(inv%invpar(iipar)==1) THEN !VP in family 1
     elseif(inv%invpar(iipar)==2) THEN !RHO in family 1
        mu1(:,:)=0.
        b(:,:)=0.
        IF(pbdir%tool == 1) THEN
           epsilo(:,:)=0.
           delta(:,:)=0.
           theta(:,:)=0. 
        ELSE
           epsilo(:,:)=pbdir%epsilon(i1,i2)
           delta(:,:)=pbdir%delta(i1,i2)
           theta(:,:)=pbdir%theta(i1,i2)
        END IF
        vc=pbdir%vp(i1,i2)*(1.-0.5*xi/pbdir%qp(i1,i2))
        mu1(0,0)=1./(pbdir%rho(i1,i2)**2*vc*vc) !=(1 /rho^2 vp^2)
        b(0,0)=1./pbdir%rho(i1,i2)**2 ! 1/rho^2
        mode=2
     elseif(inv%invpar(iipar)==12) THEN !LOG(RHO) in family 1
        mu1(:,:)=0.
        b(:,:)=0.
        IF(pbdir%tool == 1) THEN
           epsilo(:,:)=0.
           delta(:,:)=0.
           theta(:,:)=0. 
        ELSE
           epsilo(:,:)=pbdir%epsilon(i1,i2)
           delta(:,:)=pbdir%delta(i1,i2)
           theta(:,:)=pbdir%theta(i1,i2)
        END IF
        vc=pbdir%vp(i1,i2)*(1.-0.5*xi/pbdir%qp(i1,i2))
        mu1(0,0)=1./(pbdir%rho(i1,i2)*vc*vc) !=(rho /rho^2 vp^2)
        b(0,0)=1./pbdir%rho(i1,i2) ! rho/rho^2
        mode=2
    elseif(inv%invpar(iipar)==32) THEN !LOG(RHO/RHO_0) in family 1
        mu1(:,:)=0.
        b(:,:)=0.
        IF(pbdir%tool == 1) THEN
           epsilo(:,:)=0.
           delta(:,:)=0.
           theta(:,:)=0. 
        ELSE
           epsilo(:,:)=pbdir%epsilon(i1,i2)
           delta(:,:)=pbdir%delta(i1,i2)
           theta(:,:)=pbdir%theta(i1,i2)
        END IF
        vc=pbdir%vp(i1,i2)*(1.-0.5*xi/pbdir%qp(i1,i2))
        mu1(0,0)=1./(pbdir%rho(i1,i2)*vc*vc) !=(rho /rho^2 vp^2)
        b(0,0)=1./pbdir%rho(i1,i2) ! rho/rho^2
        mode=2
     elseif(inv%invpar(iipar)==4) THEN !EPSILON in family 1
     elseif(inv%invpar(iipar)==10) THEN !LOG(VP) in family 1
        mu1(:,:)=0.
        b(:,:)=0.
        epsilo(:,:)=0.
        delta(:,:)=0.
        theta(:,:)=0.
        vc=pbdir%vp(i1,i2)*(1.-0.5*xi/pbdir%qp(i1,i2))
        mu1(0,0)=-2./(pbdir%rho(i1,i2)*vc*vc) !=(-2 vp /rho vp^3)
        mode=1
     end if
  end if

  grad_term=0.
  !step1 : dAdmi
  CALL sub_dAdmi(i1e,i2e,b,mu1,delta,epsilo,theta,ianiso,ne, irn,icn,mat,pbdir%h,omegac,nnz,pbdir%n1e,pbdir%n2e,mode)
  mat(:)=-mat(:) !minus sign for consistency with other dAdmi computation
  DO isrc=1,nsrc
     !step 2 : sparse matrix vector product -> sparse vector    
     CALL sub_P_dAdmi(mat,irn,icn,nnz,ne, inc_wfld(:,:,:,isrc),pbdir%n1e*pbdir%n2e*nb_wfld,vect2,ivect2,ne2,nnz2)
     !step 3 : sparse vector/dense vector product
     CALL sub_P_dAdmi_P(adj_wfld(:,:,:,isrc),pbdir%n1e*pbdir%n2e*nb_wfld,vect2,ivect2,ne2,nnz2,val)
     !step4 : take into account source term
     grad_term=grad_term+val*inv%cc1(isrc,iw)
  END DO
end subroutine sub_graddAdmi
