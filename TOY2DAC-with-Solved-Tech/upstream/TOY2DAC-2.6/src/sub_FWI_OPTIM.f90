!Copyright 2008-2011 SEISCOPE project, All rights reserved.
!Copyright 2013-2021 SEISCOPEII project, All rights reserved.
subroutine sub_FWI_OPTIM(pbdir,acqui,inv,optim)

  implicit none

  include "mpif.h"
#include "common.h"
#include "pbdirect.h"
  include "acqui.h"
  include "inversion.h"
  include "optim_type.h"

  !IN
  TYPE(acquisition) :: acqui
  !IN/OUT
  TYPE (pbdirect) :: pbdir  
  TYPE(inversion) :: inv
  TYPE(optim_type) :: optim
  !Local variables
  logical :: convergence_test
  character*4 :: FLAG
  integer :: nn,cpt_iter,nn_tilde
  real :: fcost_init,norm_grad,norm_grad_preco,fcost0,norm_diff,norm_init
  real,dimension(:),allocatable :: x,grad,grad_preco
  complex,allocatable,dimension(:,:,:,:,:) :: inc_wfld
  complex,allocatable,dimension(:,:,:,:,:) :: adj_wfld    
  integer :: opt_preco,opt_newton
  character(LEN=2) :: flag_Newton
  real,dimension(:),allocatable :: x_ref,x_tilde,grad_tilde,big_d,big_Hd,q_preco
  real,dimension(:),allocatable :: x_prev,diff
  real,dimension(:),allocatable :: v_unit,Hv_unit
  real,dimension(:,:),allocatable :: Hv_mat

  !MPI
  INTEGER :: infompi

  !******************************
  !       INITIALIZATION        !
  !******************************
  !Compute the problem size nn
  nn=pbdir%n1*pbdir%n2*inv%npar  
  nn_tilde=nn

  opt_preco=0
  opt_newton=0

  IF(inv%opt_meth.eq.2.OR.inv%opt_meth .eq. 4) opt_preco=1 !activate the preconditioning flag for PSTD or PNLCG or PLBFGS
  IF(inv%opt_meth .eq. 7 .OR. inv%opt_meth .eq. 9) flag_Newton='GN'
  IF(inv%opt_meth .eq. 8 .OR. inv%opt_meth .eq. 10) flag_Newton='EN'
  IF(inv%opt_meth.eq. 7 .OR. inv%opt_meth .eq. 9.OR.inv%opt_meth .eq. 8 .OR. inv%opt_meth .eq. 10) opt_newton=1 !newton like optimization, requiring second order adjoint

  !recompute here the good nn_tilde for hessian based optim, as it works on reduced space parameterization 
  IF(mype==0.and. opt_newton==1) call sub_compute_nn_tilde (nn,pbdir,inv,nn_tilde)
  CALL MPI_BCAST(nn_tilde,1, MPI_INTEGER,0,MPI_COMM_WORLD,infompi)

  !HERE MANAGE FAKE SIZE FOR SLAVES
  IF(mype==0)THEN
     pbdir%n1e_s=pbdir%n1e
     pbdir%n2e_s=pbdir%n2e
     pbdir%nb_wfld_s=pbdir%nb_wfld
     pbdir%nw_s=pbdir%nw
  ELSE   
     pbdir%n1e_s=1
     pbdir%n2e_s=1
     pbdir%nb_wfld_s=1
     pbdir%nw_s=pbdir%nw
  end if
  !Allocations
  allocate(inc_wfld(pbdir%n1e_s,pbdir%n2e_s,pbdir%nb_wfld_s,acqui%nsrc,pbdir%nw_s))
  allocate(adj_wfld(pbdir%n1e_s,pbdir%n2e_s,pbdir%nb_wfld_s,acqui%nsrc,pbdir%nw_s))
  IF(mype==0) THEN
     allocate(x(nn),grad(nn))
     if(opt_newton==0) ALLOCATE(grad_preco(nn))
     if(opt_newton==1) ALLOCATE(x_tilde(nn_tilde),grad_tilde(nn_tilde),x_ref(nn))
  END IF

  !  !Scaling

  !  call sub_scaling(nn,pbdir,acqui,inv,inc_wfld,adj_wfld,x,grad) 

  !Factorize the initial modeling matrix and recompute the initial gradient
  call sub_modeling_FWI_grad(pbdir,acqui,inv,inc_wfld,adj_wfld)   

  if(mype.eq.0) then     
     call sub_setmodel2x_real(nn,x,pbdir,inv,inv%model)
     call sub_setmodel2x_real(nn,grad,pbdir,inv,inv%gradient)
     if(opt_newton==0) THEN
        if(opt_preco.eq.1 .or.inv%opt_meth .eq. 6 ) then
           !Preconditioning
           !1. Computation of the scaling factor
           inv%scal_preco=1.
           call Shin_preconditioning(nn,pbdir,acqui,inv,inc_wfld,grad,grad_preco)
           call normL2(nn,grad_preco,norm_grad_preco) 
           call normL2(nn,grad,norm_grad) 
           inv%scal_preco=norm_grad/norm_grad_preco
           !2. Actual preconditioning
           call Shin_preconditioning(nn,pbdir,acqui,inv,inc_wfld,grad,grad_preco)        
        else
           grad_preco(:)=grad(:)
        endif
     elseif(opt_newton==1) THEN
        !Reduce vectors
        x_ref(:)=x(:)
        call sub_reduce_vector(nn,pbdir,acqui,inv,grad,nn_tilde,grad_tilde)
        call sub_reduce_vector(nn,pbdir,acqui,inv,x,nn_tilde,x_tilde)
        !Compute preconditioned gradient
        allocate(grad_preco(nn_tilde))
        inv%scal_preco=1.
        call Shin_preconditioning_frozen(nn_tilde,pbdir,acqui,inv,inc_wfld,grad_tilde,grad_preco)
        call normL2(nn_tilde,grad_preco,norm_grad_preco)       
        call normL2(nn_tilde,grad_tilde,norm_grad)       
        inv%scal_preco=norm_grad/norm_grad_preco
        call Shin_preconditioning_frozen(nn_tilde,pbdir,acqui,inv,inc_wfld,grad_tilde,grad_preco)
     end if
  endif

  !******************************
  !       OPTIMIZATION LOOP     !
  !******************************  

  if(opt_newton==1) THEN
     !INTERNAL CG MAX ITERATION NUMBER
     optim%niter_max_CG=inv%niter_CG
  end if

  !ADDITIONAL STOPPING CRITERION ON MODEL CHANGES
  if(mype.eq.0) then
     allocate(x_prev(nn))
     x_prev(:)=x(:)
     allocate(diff(nn))
     call norml2(nn,x_prev,norm_init)
     norm_diff=0e0
     open(14,file='norm_model')
     write(14,*) 'NORM MODEL INIT : ',norm_init
  endif

  if(inv%opt_meth.eq.5.OR.inv%opt_meth.eq.6)   optim%l=inv%lbfgs_m !LBFGS MEMORY OPTION

  !PRECONDITIONING temporary variable for PLBFGS
  if(mype==0 .and. inv%opt_meth.eq.6) allocate(q_preco(nn))

  !INFORMATION ON REGULARIZATION COST / DATA COST
  if(mype.eq.0) then
     fcost0=inv%fcost
     open(30,file='fcost_data_reg')
  endif

  FLAG='INIT'       
  cpt_iter=0

  if(opt_newton==1) THEN
     allocate(big_d(nn),big_Hd(nn))
  end if

  do while( (FLAG.ne.'CONV').and.(FLAG.ne.'FAIL'))  
     !OPTIMIZER
     if(mype.eq.0) then
        WRITE(*,*)'inv%fcost, tmp version as not yet accepted',inv%fcost
        !!WRITE(33,*)'inv%fcost',inv%fcost
        if(inv%opt_meth.eq.1.OR.inv%opt_meth.eq.2)   call PSTD(nn,x,inv%fcost,grad,grad_preco,optim,FLAG)
        if(inv%opt_meth.eq.3.OR.inv%opt_meth.eq.4)   call PNLCG(nn,x,inv%fcost,grad,grad_preco,optim,FLAG)
        if(inv%opt_meth.eq.5)                        call LBFGS(nn,x,inv%fcost,grad,optim,FLAG)
        if(inv%opt_meth.eq.6)                        call PLBFGS(nn,x,inv%fcost,grad,grad_preco,optim,FLAG)
        if(inv%opt_meth.eq.7.or.inv%opt_meth.eq.8)   call TRN(nn_tilde,x_tilde,inv%fcost,grad_tilde,optim,FLAG)
        if(inv%opt_meth.eq.9.or.inv%opt_meth.eq.10)  call PTRN(nn_tilde,x_tilde,inv%fcost,grad_tilde,grad_preco,optim,FLAG)
     endif

     CALL MPI_BCAST(FLAG,4, MPI_CHARACTER,0,MPI_COMM_WORLD,infompi)

     if(FLAG.eq.'GRAD') then
        if(mype.eq.0) then
           if(opt_newton==1) call sub_augment_model(nn,pbdir,acqui,inv,x,x_ref,nn_tilde,x_tilde)
           CALL sub_setx2model_real(nn,x,pbdir,inv,inv%model)
           CALL sub_modelinv2pbdir(pbdir,inv)
        endif
        call sub_modeling_FWI_grad(pbdir,acqui,inv,inc_wfld,adj_wfld)
        if(mype.eq.0.AND.opt_newton==0) then
           CALL sub_setmodel2x_real(nn,x,pbdir,inv,inv%model)
           CALL sub_setmodel2x_real(nn,grad,pbdir,inv,inv%gradient)
           if(opt_preco.eq.1) then
              !Preconditioning
              !1. Computation of the scaling factor
              inv%scal_preco=1.
              call Shin_preconditioning(nn,pbdir,acqui,inv,inc_wfld,grad,grad_preco)
              call normL2(nn,grad_preco,norm_grad_preco) 
              call normL2(nn,grad,norm_grad) 
              inv%scal_preco=norm_grad/norm_grad_preco
              !2. Actual preconditioning
              call Shin_preconditioning(nn,pbdir,acqui,inv,inc_wfld,grad,grad_preco)
           else
              grad_preco(:)=grad(:)
           endif
        elseif(mype.eq.0.AND.opt_newton==1) then
           call sub_setmodel2x_real(nn,grad,pbdir,inv,inv%gradient)
           call sub_reduce_vector(nn,pbdir,acqui,inv,grad,nn_tilde,grad_tilde)
           !Compute preconditioned gradient 
           inv%scal_preco=1.
           call Shin_preconditioning_frozen(nn_tilde,pbdir,acqui,inv,inc_wfld,&
                grad_tilde,grad_preco)
           call normL2(nn_tilde,grad_preco,norm_grad_preco)       
           call normL2(nn_tilde,grad_tilde,norm_grad)       
           inv%scal_preco=norm_grad/norm_grad_preco
           call Shin_preconditioning_frozen(nn_tilde,pbdir,acqui,inv,inc_wfld,&
                grad_tilde,grad_preco) 
        endif
     elseif(FLAG.eq.'HESS') then
        !Hessian-vector computation
        if(mype.eq.0) then
           call sub_augment_vector(nn,pbdir,acqui,inv,big_d,nn_tilde,optim%d)
        endif
        CALL MPI_BCAST(big_d,nn,MPI_REAL,0,MPI_COMM_WORLD,infompi)
        if(flag_Newton.eq.'GN') then
           call sub_modeling_FWI_GN(pbdir,acqui,inv,inc_wfld,big_d,big_Hd)
        elseif(flag_Newton.eq.'EN') then        
           call sub_modeling_FWI_EN(pbdir,acqui,inv,inc_wfld,adj_wfld,big_d,big_Hd)
        endif
        if(mype.eq.0) then
           call sub_reduce_vector(nn,pbdir,acqui,inv,big_Hd,nn_tilde,optim%Hd)
        endif
     elseif(FLAG.eq.'PREC'.AND. inv%opt_meth.eq.6) then
        if(mype.eq.0) then           
           CALL sub_setx2model_real(nn,x,pbdir,inv,inv%model) !ADDED 02/05/13
           CALL sub_modelinv2pbdir(pbdir,inv)
           inv%scal_preco=1.
           call Shin_preconditioning(nn,pbdir,acqui,inv,inc_wfld,grad,grad_preco)
           call normL2(nn,grad,norm_grad)
           call normL2(nn,grad_preco,norm_grad_preco)
           inv%scal_preco=norm_grad/norm_grad_preco           
           call Shin_preconditioning(nn,pbdir,acqui,inv,inc_wfld,optim%q_plb,q_preco)
           optim%q_plb(:)=q_preco(:)  
        endif
     elseif(FLAG.eq.'PREC'.and.opt_newton==1) then
        !Preconditioner applied to optim%residual !MODIFIED 02/05/13
        if(mype.eq.0) then           
           call Shin_preconditioning_frozen(nn_tilde,pbdir,acqui,inv,inc_wfld,&
                optim%residual,optim%residual_preco)           
        endif
     elseif(FLAG.eq.'NSTE') then
        if(mype.eq.0) then
           cpt_iter=cpt_iter+1
           CALL sub_setx2model_real(nn,x,pbdir,inv,inv%model)
           CALL sub_modelinv2pbdirforwriting(pbdir,inv,cpt_iter,0)
           call sub_modelinv2pbdirforwriting(pbdir,inv,1,1) !write the last temporary model in final files
           open(20,file='invparinter',access='direct',recl=4*nn)
           write(20,rec=cpt_iter)inv%model(:,:,1:inv%npar)
           close(20)
           !WRITE REGULARIZATION INFORMATION
           write(30,*) cpt_iter,&
                inv%fcost/fcost0,&
                inv%fcost_data/fcost0,&
                inv%fcost_reg/fcost0
           diff(:)=x(:)-x_prev(:)
           call normL2(nn,diff,norm_diff)           
           write(14,*) norm_diff,norm_diff/norm_init
           norm_diff=norm_diff/norm_init           
           if(norm_diff<inv%convergence_criterion_model) then
              FLAG='CONV'
               WRITE(*,*)'CONVERGENCE ASSOCIATED TO THE MODEL UPDATE CRITERION'
           endif
           x_prev(:)=x(:)
        endif
        CALL MPI_BCAST(FLAG,4, MPI_CHARACTER,0,MPI_COMM_WORLD,infompi)
     endif
  enddo
  if(mype.eq.0) then
     cpt_iter=cpt_iter+1
     write(30,*) cpt_iter,&
          inv%fcost/fcost0,&
          inv%fcost_data/fcost0,&
          inv%fcost_reg/fcost0
     close(30)
  endif

  WRITE(*,*)'EXIT THE OPTIMIZATION LOOP WITH FLAG VALUE ',FLAG

end subroutine sub_FWI_OPTIM

