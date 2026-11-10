!Copyright 2008-2011 SEISCOPE project, All rights reserved.
!Copyright 2013-2021 SEISCOPEII project, All rights reserved.

!*********************************************************************
!                          TOY2DAC code                              !
!    a user friendly 2D acoustic Full Waveform Inversion package     !
!--------------------------------------------------------------------!
!This package have been developed for testing methodological concepts!
!of FWI in the framework of 2D acoustic approximation, in the        !
!frequency-domain. This package is developped basicaly in sequential !
!(with openMP threads) in order to keep things simple for quick tests!
!and developements... only the LU decomposition done by MUMPS is done!
!with MPI directives                                                 !
!--------------------------------------------------------------------!
!                    DEVELOPMENT HISTORY                             !
! V2.6  - 05/2019  - R. Brossier   & L. Metivier                     !
! correct bug                                                        !
!                                                                    !
! V2.5  - 07/2018  - R. Brossier                                     !
! few bugs + patch the code for compatibilty with gfortran           !
!                                                                    !
! V2.4  - 08/2016  - R. Brossier                                     !
! *correction of several bugs associated to erroneous allocations of !
! table and other problems feedbacked from users                     !
! *merge back the branche developed for the Nankai application of    !
! Andrew including now : complex frequencies, vertical source        !
! weight function...                                                 !
! *correct the documentation from errors                             !
! *implement the choice of single or double precision for MUMPS      !
! usefull for complex frequencies                                    !
!                                                                    !
! V2.3  - 05/2015  - R. Brossier - L. Metivier                       !
! correct few bugs                                                   !
!                                                                    !
! V2.2  - 10/2014  - R. Brossier - L. Metivier                       !
! few bugs associated to anisotropic inversion should be fixed       !
!                                                                    !
! V2.1  - 08/2014  - R. Brossier - L. Metivier                       !
! Reintegration of the V2 in the trunk as reference version          !
!                                                                    !
! V2.0  - 09/2013  - R. Brossier - L. Metivier                       !
!Development in a branche, in parallel of the simple V1.0 version    !
!This new version is connected to the SEISCOPE OPTIMIZATION TOOLBOX  !
!and implements second-order adjoint state formulas for the          !
!computation of Hessian-vector products. The use of the different    ! 
!optimization schemes proposed in the TOOLBOX is thus possible       !
!Three additional functionalities have aslo been implemented in this !
! version                                                            !
! - free surface condition                                           !
! - acoustic anisotropic modeling                                    !
! - Hicks interpolation for sources and receivers                    !
!                                                                    !
! V1.0  - 09/2011  - R. Brossier                                     !
!        version based on FWT2D V5.4 code of S. Operto for the       !
!        forward problem.                                            !
!        inversion based on LBFGSB routine of Nocedal                !
!--------------------------------------------------------------------!
!********************************************************************!

program toy2dac

  USE MPI_TIMES
  USE sub_read_acqui

  implicit none

  include "mpif.h"
  include "acqui.h"
#include "pbdirect.h"
#include "common.h"

  !GENERAL MODES
  INTEGER :: mode

  !PBDIRECT
  TYPE (pbdirect) :: pbdir

  !MPI
  INTEGER                               :: ierr,infompi, mynode_lenght, iproc
  character*(MPI_MAX_PROCESSOR_NAME)    ::  mynode

  !acquisition
  TYPE (acquisition) :: acqui
  REAL,   DIMENSION(:,:),ALLOCATABLE :: coord_rec_tmp
  INTEGER                               ::  isrc


  !*************************
  ! INIT MPI

  CALL MPI_INIT(infompi)

  ! Define Host and Slave Processors
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,mype,infompi)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,infompi)
  CALL MPI_GET_PROCESSOR_NAME(mynode, mynode_lenght, infompi)
  CALL MPI_TIME(0)


  if(mype==0) THEN
     WRITE(*,*)
     WRITE(*,*)'***********************************'
     WRITE(*,*)'  Welcome in TOY2DAC code     '
     WRITE(*,*)'  V2.6  '
     WRITE(*,*)'***********************************'
#ifdef DOUBLEPRECISION
     WRITE(*,*)'  TOY2DAC is compiled with double precision flag'
     WRITE(*,*)'  ATTENTION, this code will read and write data files'
     WRITE(*,*)'  IN DOUBLE PRECISION !!!!'
#else
     WRITE(*,*)'  TOY2DAC is compiled with single precision flag'
     WRITE(*,*)'  ATTENTION, this code will read and write data files'
     WRITE(*,*)'  IN SINGLE PRECISION !!!!'
#endif

     WRITE(*,*)
  END IF


  !****************

  !STEP 1 -> READ INPUT
  OPEN(10,file='toy2dac_input')
  READ(10,*)mode
  READ(10,*)pbdir%tool
  READ(10,*)acqui%file_acqui  
  CLOSE(10)
  
  !****************
  !STEP 2 -> PRE PROCESSING FOR FORWARD MODELING
  !Hicks interpolation: read and build coefficients  
  IF (pbdir%tool == 1) THEN
     ! ACOUSTIC ISOTROPIC MODELING IN THE FREQUENCY DOMAIN
     pbdir%nb_wfld=1
     CALL sub_init_FDFD(pbdir)
  ELSEIF(pbdir%tool == 2) THEN
     ! ACOUSTIC ANISOTROPIC MODELING IN THE FREQUENCY DOMAIN
     pbdir%nb_wfld=2
     CALL sub_init_FDFDaniso(pbdir)     
  END IF
  !****************
  !STEP 3 -> REAL WORK

  !STEP 3.1 -> READ ACQUISITION
  IF(mype==0) THEN
     CALL sub_read_acqui1(acqui%file_acqui,acqui%nline,acqui%ntrace,&
          acqui%nsrc,acqui%nrecmax)
     !ALLOCATE definitive table
     ALLOCATE(acqui%coord_src(acqui%nsrc,2))
     ALLOCATE(acqui%nbr_rec(acqui%nsrc))
     ALLOCATE(acqui%associated_rec(acqui%nsrc,acqui%nrecmax))
     ALLOCATE(acqui%first_arrival_t0(acqui%nsrc,acqui%nrecmax))
     acqui%coord_src(:,:)=0.
     acqui%nbr_rec(:)=0
     acqui%associated_rec(:,:)=0
     !ALLOCATE TMP table
     ALLOCATE(coord_rec_tmp(acqui%ntrace,2))
     coord_rec_tmp(:,:)=0.
     CALL sub_read_acqui2(acqui%file_acqui,acqui%nline,acqui%ntrace,&
          acqui%nsrc,acqui%nrecmax,acqui%coord_src,&
          acqui%nbr_rec,acqui%associated_rec,acqui%first_arrival_t0,&
          coord_rec_tmp,acqui%nrec_tot)
     !ALLOCATE definitive table
     ALLOCATE(acqui%coord_rec(acqui%nrec_tot,2))
     acqui%coord_rec(:,:)=coord_rec_tmp(1: acqui%nrec_tot,:)
     DEALLOCATE(coord_rec_tmp)
     !Hicks interpolation
     if(pbdir%hicks.eq.1) then        
        call sub_Hicks_constant(pbdir,acqui,acqui%coord_src,acqui%coord_rec)
     endif
  END IF

  !BCAST acqui%nsrc
  CALL MPI_BCAST(acqui%nsrc,1, MPI_INTEGER,0,MPI_COMM_WORLD,infompi)
  !BCAST pbdir%n1e,pbdir%n2e
  CALL MPI_BCAST(pbdir%n1e,1, MPI_INTEGER,0,MPI_COMM_WORLD,infompi)
  CALL MPI_BCAST(pbdir%n2e,1, MPI_INTEGER,0,MPI_COMM_WORLD,infompi)

  !STEP 3.2 WORK....
  IF (mode==0) THEN
     ! PB DIRECT
     pbdir%mode=1 ! fcost only
     CALL sub_modeling(pbdir,acqui)
  ELSEIF(mode==1) THEN
     !FWI
     CALL sub_FWI(pbdir,acqui)
  ELSEIF(mode==2) THEN
     !RTM
  ELSEIF(mode==3) THEN
     !MVA
  END IF
  
  !****************
  !STEP FINAL -> FINALIZE 
  IF (pbdir%tool == 1) THEN
     ! ACOUSTIC ISOTROPIC MODELING IN THE FREQUENCY DOMAIN
     CALL sub_finish_FDFD(pbdir)
  END IF
  CALL MPI_BARRIER( MPI_COMM_WORLD, infompi )
  CALL MPI_TIME(1)
  CALL flush()
  WRITE(*,*) '-----------END OF PROGRAM---------------'
  CALL flush()
  CALL MPI_FINALIZE(ierr)

  
end program toy2dac
