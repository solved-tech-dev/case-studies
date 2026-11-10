TYPE inversion
                                                                      
SEQUENCE

COMPLEX(kind=creal),ALLOCATABLE,DIMENSION(:,:) :: obs_data
COMPLEX(kind=creal),ALLOCATABLE,DIMENSION(:) ::residual

  !INVERSION
  INTEGER :: npar,npar0,invfamily,invpar(10),opt_meth,misfit_fct_type
  REAL,ALLOCATABLE,DIMENSION(:,:) :: limit_par,limit_pbdir_par
  REAL,DIMENSION(:,:,:),ALLOCATABLE :: gradient
  REAL,DIMENSION(:,:,:),ALLOCATABLE :: model
  REAL,DIMENSION(:,:),ALLOCATABLE :: previous_model
  REAL,DIMENSION(:),ALLOCATABLE     :: bathy
  INTEGER,DIMENSION(:),ALLOCATABLE     :: ibathy
  REAL :: scalingfactor,fcost_model
  REAL :: fcost
  REAL,dimension(:,:),ALLOCATABLE :: prior0
  INTEGER :: firstgrad,src_estim,zbathy,niter_nonlin,niter_CG,nls_max
  INTEGER :: ideadzone !zone bellow bathy where gradient is not computed
  ! but is allowed to change through regularization or preconditionning
  REAL :: coeff_damping_preco,scal_preco
  INTEGER :: lbfgs_m
  INTEGER :: debug_option
  COMPLEX,DIMENSION(:,:),ALLOCATABLE :: cc1
  INTEGER :: regul
  REAL,dimension(:),allocatable :: lambda
  REAL :: lambda_x,lambda_z,lambda_prior0
  REAL :: norm_grad_pure      
  REAL,DIMENSION(10) :: refvalue

  !BOUND CONSTRAINTS
  integer :: bound
  real ::  threshold
  REAL ::  lb(10),ub(10)
    
  !HESSIAN EIGENVECTORS COMPUTATION
  integer :: nev
  character*2 :: which,opt_type
  
  !REGULARIZATION INFORMATION
  real :: fcost_data,fcost_reg

  !STOPPING CRITERION
  REAL :: convergence_criterion
  REAL :: convergence_criterion_model

  !WEIGHTING
   CHARACTER(LEN=80) :: data_weight_file 
   REAL    :: dx_data_weight
   INTEGER :: n_data_weight

END TYPE inversion
