INTEGER   :: maxfreq
PARAMETER (maxfreq = 1000 )

!MPI
INTEGER :: mype,nproc
COMMON mype,nproc

REAL :: epsilondist,pi

PARAMETER (epsilondist =0.001) 
PARAMETER (pi = 3.14159265)
INTEGER flag_first_call_to_mumps
COMMON flag_first_call_to_mumps
REAL :: lmoy
COMMON lmoy
!PARAMETER(lmoy=1.)

INTEGER :: creal


#ifdef DOUBLEPRECISION
PARAMETER(creal=8)
#else
PARAMETER(creal=4)
#endif


