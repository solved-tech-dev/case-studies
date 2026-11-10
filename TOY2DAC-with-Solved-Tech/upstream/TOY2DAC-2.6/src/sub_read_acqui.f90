!Copyright 2008-2011 SEISCOPE project, All rights reserved.
!Copyright 2013-2021 SEISCOPEII project, All rights reserved.

!###########################################################################
! Subroutine to scan acquisition file to build all table relatives to 
! source-receiver positions and organisation
!
! INPUT : unit of file to scan
!
! OUTPUT : nline -> number of line in file
!          ntrace -> total number of traces
!          nrec_tot -> total number of independant receivers (i.e. if 1 receiver 
!                      is commun to 2 sources it is counted only once
!          nsrc -> number of sources
!          nrecmax -> maximum number of receivers for one source
!          coord_src(nsrc,2) -> coordinates of sources (x,z)
!          coord_rec(nrec_tot,2) -> coordinates of receivers (x,z)
!          associated_rec(nsrc,nrecmax) -> associate for each source the indice of all its recievers
!          nbr_rec(nsrc) -> number of receivers per source
!
! MODIFICATION : 10/2008 V1.2: Adding acquistruct%first_arrival_t0 table
!
!############################################################################          
MODULE sub_read_acqui 
CONTAINS  
  SUBROUTINE sub_read_acqui1(file_acqui,nline,ntrace,nsrc,nrecmax) 


    IMPLICIT NONE
#include "common.h"

    !INPUT
    CHARACTER(LEN=80) :: file_acqui

    !OUTPUT
    INTEGER :: nline,ntrace,nsrc,nrecmax

    !INTERNAL VARIABLE
    INTEGER                            :: i,j,nrec_tmp
    REAL                               :: x1,x2,t0,inc
    INTEGER                            :: icode

    !INIT VARIABLE
    nline=0
    ntrace=0
    nsrc=0
    nrecmax=0
    nrec_tmp=0

    OPEN(11,file=file_acqui)

    !READ FILE FOR COMPUTATION OF nrecmax,nline,nsrc,ntrace
    REWIND(11)
10  CONTINUE
    READ(11,*,end=99) x2,x1,t0,inc,icode    
    nline=nline+1
    IF (icode == 0) THEN !SOURCE
       nsrc=nsrc+1
       nrec_tmp=0
    ELSE IF (icode ==1) THEN ! RECEIVER
       ntrace=ntrace+1
       nrec_tmp=nrec_tmp+1
       IF (nrec_tmp > nrecmax) nrecmax=nrec_tmp
    END IF
    GOTO 10
    
99  CONTINUE  ! FILE READ ONCE, TEMPORARY TABLE CAN NOW BE ALLOCATE & INITIALIZED

    WRITE(*,*)'SCAN ACQUISITION '
    WRITE(*,*)'in file ', 11
    WRITE(*,*)'ntrace= ' ,ntrace
    WRITE(*,*)'nsrc= ' ,nsrc
    WRITE(*,*)'nrecmax= ' ,nrecmax

    IF(nsrc<=0 .OR. nrecmax<=0 .OR. ntrace<=0) THEN
       WRITE(*,*)'PROBLEM IN THE ACQUISITION FILE'
       WRITE(*,*)'CODE STOPPED NOW'
       CALL sub_MPI_ABORT()
    END IF
    
    
  END SUBROUTINE sub_read_acqui1
  
  SUBROUTINE sub_read_acqui2(file_acqui,nline,ntrace,nsrc,nrecmax,coord_src,nbr_rec,associated_rec,first_arrival_t0,coord_rec_tmp,nrec_tot) 
    
    IMPLICIT NONE
#include "common.h"

    CHARACTER(LEN=80) :: file_acqui
    !IN/OUT
    REAL :: coord_src(nsrc,2), first_arrival_t0(nsrc,nrecmax),coord_rec_tmp(ntrace,2)
    INTEGER :: nbr_rec(nsrc),associated_rec(nsrc,nrecmax)
    INTEGER :: nline,ntrace,nsrc,nrecmax
    INTEGER                               ::  nrec_tot 
    REAL                               :: x1,x2,t0,inc
    INTEGER                            :: icode,isrc,irec,i,nrec_tmp
    INTEGER                            :: func_test_equality

    REWIND(11)

    isrc=0
    irec=0
    nrec_tot=0

    DO i=1,nline
       READ(11,*) x2,x1,t0,inc,icode
       IF (icode == 0) THEN !SOURCE
          nrec_tmp=0
          isrc=isrc+1
          coord_src(isrc,1)=x1  !x coord
          coord_src(isrc,2)=x2  !z coord
       ELSE IF (icode == 1) THEN ! RECEIVER
          nbr_rec(isrc)=nbr_rec(isrc)+1
          DO irec=1,nrec_tot
             !TEST PREVIOUS COORDINATES
             IF (func_test_equality(x1,coord_rec_tmp(irec,1)) == 0) THEN
                IF (func_test_equality(x2,coord_rec_tmp(irec,2)) == 0) THEN
                   associated_rec( isrc, nbr_rec(isrc))= irec     
                   ! No of current receiver of source isrc in global table of receivers 
                   first_arrival_t0( isrc, nbr_rec(isrc))= t0                   
                   GOTO 20
                END IF
             END IF
          END DO
          ! WE ARRIVE HERE ONLY IF RECEIVER IS NOT YET FUND
          nrec_tot                     = nrec_tot + 1
          coord_rec_tmp( nrec_tot, 1 ) = x1
          coord_rec_tmp( nrec_tot, 2 ) = x2
          associated_rec( isrc , nbr_rec(isrc) ) = nrec_tot          
          first_arrival_t0( isrc, nbr_rec(isrc))= t0
20        CONTINUE
       END IF
    END DO
    REWIND(11)
    
  END SUBROUTINE sub_read_acqui2
END MODULE sub_read_acqui
