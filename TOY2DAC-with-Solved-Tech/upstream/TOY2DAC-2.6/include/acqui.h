TYPE acquisition
                                                                      
SEQUENCE

 INTEGER       :: nsrc,nrecmax,ntrace,nline
 CHARACTER(LEN=80) :: file_acqui
  REAL,   DIMENSION(:,:),ALLOCATABLE :: coord_src,coord_rec,first_arrival_t0
  INTEGER,DIMENSION(:,:),ALLOCATABLE :: associated_rec
  REAL,   DIMENSION(:),ALLOCATABLE ::  weight_data
  INTEGER,DIMENSION(:),ALLOCATABLE :: nbr_rec
  INTEGER                      ::  nrec_tot
END TYPE acquisition
