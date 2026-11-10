!Copyright 2008-2011 SEISCOPE project, All rights reserved.
!Copyright 2013-2021 SEISCOPEII project, All rights reserved.
subroutine sub_FWI(pbdir,acqui)

  implicit none

  include "mpif.h"
#include "common.h"
#include "pbdirect.h"
  include "acqui.h"
  include "inversion.h"
  include "optim_type.h"

  !PBDIRECT
  TYPE (pbdirect) :: pbdir
  !ACQUI
  TYPE(acquisition) :: acqui
  !INVERSION
  TYPE(inversion) :: inv
  !OPTIMIZATION
  TYPE(optim_type) :: optim

  !DATA
  CHARACTER(LEN=80) :: name_data_obs,name_bathy

  REAL :: xs,xr,zs,zr
  REAL :: dist,red_dist
  INTEGER :: irec,isrc,idist,k,irec_true
  REAL,DIMENSION(:),ALLOCATABLE :: weight_table

  INTEGER :: i,j,iipar,iiipar,ifound
  REAL depth_deadzone



  !STEP 1 -> READ DATA
  !Main data
  OPEN(15,file='fwi_input')
  READ(15,*)name_data_obs !data filename
  READ(15,*)inv%invfamily !parameter family 
  READ(15,*)inv%npar,(inv%invpar(i),i=1,inv%npar) !number of parameters and index in the family  
  READ(15,*) inv%src_estim    !type of source estimation
  READ(15,*) name_bathy       !bathymetry filename
  READ(15,*) depth_deadzone   !depth to start the inversion

  !Optimization setup
  READ(15,*) inv%opt_meth     !choice of the optimization method 
  READ(15,*) inv%convergence_criterion !misfit reduction for convergence 
  READ(15,*) inv%convergence_criterion_model !model change for convergence  
  READ(15,*) inv%niter_nonlin !maximum number of nonlinear iterations  
  READ(15,*) inv%lbfgs_m      !number of gradient stored for the l-BFGS strategy
  READ(15,*) inv%niter_CG     !maximum number of inner iterations for Truncated Newton 
  READ(15,*) inv%coeff_damping_preco   !threshold parameter for the Shin preconditioner   
  READ(15,*) inv%debug_option          !threshold parameter for the Shin preconditioner   

  !Regularization setup
  ALLOCATE(inv%lambda(inv%npar))
  READ(15,*) inv%lambda(:)              !Tikhonov regularization 
  READ(15,*) inv%lambda_x, inv%lambda_z !Direction weights for the Tikhonov regularization
  READ(15,*) inv%lambda_prior0          !regularization parameter for prior information  

  !Bound constraints   
  READ(15,*) inv%bound       !bound constraints activation
  !  if(inv%bound.eq.1) then
  READ(15,*) (inv%ub(i),i=1,inv%npar)        !upper bound scalar value
  READ(15,*) (inv%lb(i),i=1,inv%npar)        !lower bound scalar value
  READ(15,*) inv%threshold !tolerance for bound constraints
  !  endif

  READ(15,*)inv%data_weight_file                      !input weight file
  READ(15,*)inv%n_data_weight,inv%dx_data_weight      ! number of sample in the weight file and spatial sampling

  CLOSE(15) ! close fwi_input
  !END READ DATA
  !CHECK IF REGULARIZATION IS ACTIVATED
  if(inv%lambda(1).lt.0.) then
     inv%regul=2 ! Preconditioning by smoothing the gradient direction
     !NORMALIZE inv%lambda_z and inv%lambda_x by h**2 because of the finite difference scheme used in preconditioning
     inv%lambda_x=inv%lambda_x/pbdir%h**2
     inv%lambda_z=inv%lambda_z/pbdir%h**2
  else
     do i=1,inv%npar
        if(inv%lambda(i).gt.0.) then
           inv%regul=1
        endif
     enddo
  endif

  if(mype==0) then
     if(inv%lambda_prior0>0.) then
        allocate(inv%prior0(pbdir%n1,pbdir%n2))
        open(20,file='prior0',access='direct',recl=4*pbdir%n1*pbdir%n2)
        read(20,rec=1) inv%prior0(:,:)
        close(20)
     endif
  endif

  if(mype==0) then
     !READ DATA
     ALLOCATE(inv%obs_data(acqui%ntrace,pbdir%nw))
     OPEN(16,file=name_data_obs,access='direct',recl=2*creal*acqui%ntrace*pbdir%nw)
     READ(16,rec=1)inv%obs_data(:,:)
     CLOSE(16)          
     !READ bathy
     ALLOCATE(inv%bathy(pbdir%n2))
     ALLOCATE(inv%ibathy(pbdir%n2))
     OPEN(17,file=name_bathy,access='direct',recl=4*pbdir%n2)
     READ(17,rec=1)inv%bathy(:)
     CLOSE(17)
     !COMPUTE IBATHY
     write(*,*)'computed bathymetry indices'
     do i=1,pbdir%n2
        inv%ibathy(i)=NINT(inv%bathy(i)/pbdir%h)+1
     end do
     write(*,*) inv%ibathy
     inv%ideadzone=NINT(depth_deadzone/pbdir%h)
     !compute weighting on data
     ALLOCATE(acqui%weight_data(acqui%ntrace))
     ALLOCATE(weight_table(inv%n_data_weight))
     OPEN(45,file=inv%data_weight_file)
     READ(45,*) (weight_table(i),i=1,inv%n_data_weight)
     CLOSE(45)

     !compute weight
     k=1
     DO isrc = 1, acqui%nsrc
        xs=acqui%coord_src(isrc,1)
        zs=acqui%coord_src(isrc,2)
        DO irec = 1, acqui%nbr_rec(isrc)
           irec_true = acqui%associated_rec(isrc,irec)
           xr=acqui%coord_rec(irec_true,1)
           zr=acqui%coord_rec(irec_true,2)
           !compute distance between source and receiver
           dist=SQRT((xs-xr)**2)! +(zs-zr)**2) -> here offset
           idist = INT(dist/inv%dx_data_weight)+1
           red_dist=dist/inv%dx_data_weight-idist+1 !reduced coordinate between 0 and 1
           !linear interpolation +extrapolation of the last value if the tabel is not big enough
           acqui%weight_data(k) = weight_table(min(idist,inv%n_data_weight)) * (1-red_dist) + weight_table(min(idist+1,inv%n_data_weight)) * red_dist
           k=k+1
        END DO
     END DO
     DEALLOCATE(weight_table)

  end if

  ! STEP 2 -> define inverted model
  if(mype==0) then
     !define here the inverted parameters and the other fixed parameters in the current family
     if(pbdir%tool==1) then
        inv%npar0=3 ! 3 parameters in isotropic FDFD
     elseif(pbdir%tool==2) then
        inv%npar0=6 ! 6 parameters in aniso FDFD
     end if
     !fill the inv%invpar with non-inverted parameters of the family
     iipar=inv%npar
     if(pbdir%tool==1.OR.pbdir%tool==2) then
        !check for vp
        !check if already set in inverted para
        ifound=0
        do iiipar=1,inv%npar
           if(inv%invpar(iiipar)==1 .OR. inv%invpar(iiipar)==10 .OR. inv%invpar(iiipar)==11 .OR. inv%invpar(iiipar)==31) then
              ifound=1
              exit
           end if
        end do
        if(ifound==0) then
           iipar=iipar+1
           inv%invpar(iipar)=1
        end if
        !check for  rho
        !check if already set in inverted para
        ifound=0
        do iiipar=1,inv%npar
           if(inv%invpar(iiipar)==2 .OR. inv%invpar(iiipar)==12 .OR. inv%invpar(iiipar)==32) then
              ifound=1
              exit
           end if
        end do
        if(ifound==0) then
           iipar=iipar+1
           inv%invpar(iipar)=2
        end if
        !check for  qp
        !check if already set in inverted para
        ifound=0
        do iiipar=1,inv%npar
           if(inv%invpar(iiipar)==3 .OR. inv%invpar(iiipar)==33) then
              ifound=1
              exit
           end if
        end do
        if(ifound==0) then
           iipar=iipar+1
           inv%invpar(iipar)=3
        end if
     end if
     if(pbdir%tool==2) then
        !check for epsilon
        !check if already set in inverted para
        ifound=0
        do iiipar=1,inv%npar
           if(inv%invpar(iiipar)==4 .OR. inv%invpar(iiipar)==34) then
              ifound=1
              exit
           end if
        end do
        if(ifound==0) then
           iipar=iipar+1
           inv%invpar(iipar)=4
        end if
        !check for  delta
        !check if already set in inverted para
        ifound=0
        do iiipar=1,inv%npar
           if(inv%invpar(iiipar)==5 .OR. inv%invpar(iiipar)==35) then
              ifound=1
              exit
           end if
        end do
        if(ifound==0) then
           iipar=iipar+1
           inv%invpar(iipar)=5
        end if
        !check for  theta
        !check if already set in inverted para
        ifound=0
        do iiipar=1,inv%npar
           if(inv%invpar(iiipar)==6 .OR. inv%invpar(iiipar)==36) then
              ifound=1
              exit
           end if
        end do
        if(ifound==0) then
           iipar=iipar+1
           inv%invpar(iipar)=6
        end if
     end if
     write(*,*) inv%invpar(:)
     allocate(inv%model(pbdir%n1,pbdir%n2,inv%npar0))
     call sub_setmodelpbdir2inv(pbdir,inv)     
  end if


  !CHECK WHAT IS POSSIBLE !!! 
  !PSEUDO HESSIAN NOT IMPLEMENTED FOR PARAMETER OTHER THAT VP, QP AND RHO
  !STOP THE CODE IF BAD CHOICE OF OPTIMIZATION
  do iiipar=1,inv%npar
     IF((inv%invpar(iiipar)/=1 .AND.inv%invpar(iiipar)/=2 .AND. inv%invpar(iiipar)/=3) &
          .AND. ( inv%opt_meth ==2 .OR. inv%opt_meth ==4  .OR. inv%opt_meth ==6 .OR. inv%opt_meth ==9 .OR. inv%opt_meth ==10)) THEN
        WRITE(*,*)'**************************************'
        WRITE(*,*)'PROBLEM WITH THE CHOICE OF PARAMETERS'
        WRITE(*,*)'AND THE OPTIMIZATION...'
        WRITE(*,*)'A PRECONDITIONNED OPTIMIZATION IS ASKED'
        WRITE(*,*)'BUT PRECONDITIONNER DOES NOT EXIST'
        WRITE(*,*)'**************************************'
        WRITE(*,*)'     CODE STOPPED NOW, SORRY'
        WRITE(*,*)'**************************************'
        CALL sub_MPI_ABORT
     END IF

  END DO

  !STEP 3 -> OPTIMIZATION  
  inv%firstgrad=0
  allocate (inv%gradient(pbdir%n1,pbdir%n2,inv%npar))
  inv%gradient(:,:,:)=0.
  allocate(inv%residual(acqui%ntrace))
  inv%residual(:)=0.  
  allocate(inv%cc1(acqui%nsrc,pbdir%nw))
  inv%cc1(:,:)=1.    

  !common nitialization of the optimization routine
  call sub_init_optim(pbdir,acqui,inv,optim)

  call sub_FWI_OPTIM(pbdir,acqui,inv,optim)

!!$  if(inv%opt_meth.eq.1) then                          !1  PSTD TOOLBOX NO PRECONDITIONING
!!$     call sub_FWI_PSTD(pbdir,acqui,inv,optim,0)                
!!$  elseif(inv%opt_meth.eq.2) then                      !2  PSTD TOOLBOX WITH PRECONDITIONING
!!$     call sub_FWI_PSTD(pbdir,acqui,inv,optim,1)                
!!$  elseif(inv%opt_meth .eq. 3) then                    !3  PNLCG TOOLBOX NO PRECONDITIONING
!!$     call sub_FWI_PNLCG(pbdir,acqui,inv,optim,0)                
!!$  elseif(inv%opt_meth .eq. 4) then                    !4  PNLCG WITH PRECONDITIONING
!!$     call sub_FWI_PNLCG(pbdir,acqui,inv,optim,1)                
!!$  elseif(inv%opt_meth .eq. 5) then                    !5  LBFGS TOOLBOX
!!$     call sub_FWI_LBFGS(pbdir,acqui,inv,optim)                
!!$  elseif(inv%opt_meth .eq. 6) then                    !6  PLBFGS TOOLBOX
!!$     call sub_FWI_PLBFGS(pbdir,acqui,inv,optim)                
!!$  elseif(inv%opt_meth .eq. 7) then                    !7  TGN TOOLBOX/ GN approximation
!!$     call sub_FWI_TRN(pbdir,acqui,inv,optim,'GN')                
!!$  elseif(inv%opt_meth .eq. 8) then                    !8  TRN TOOLBOX
!!$     call sub_FWI_TRN(pbdir,acqui,inv,optim,'EN')                 
!!$  elseif(inv%opt_meth .eq. 9) then                    !9  PTGN TOOLBOX/ GN approximation
!!$     call sub_FWI_PTRN(pbdir,acqui,inv,optim,'GN')                
!!$  elseif(inv%opt_meth .eq. 10) then                   !10 PTRN TOOLBOX
!!$     call sub_FWI_PTRN(pbdir,acqui,inv,optim,'EN')                
!!$  endif

  if(mype.eq.0) then
     call sub_modelinv2pbdirforwriting(pbdir,inv,1,1)
  endif

end subroutine sub_FWI



