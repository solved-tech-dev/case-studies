subroutine sub_FWI_LBFGS(pbdir,acqui,inv)
  implicit none


  include "mpif.h"
  include "common.h"
  include "pbdirect.h"
  include "acqui.h"
  include "inversion.h"
  include "lbfgsb.h"

  !PBDIRECT
  TYPE (pbdirect) :: pbdir
  !ACQUI
  TYPE(acquisition) :: acqui
  !INVERSION
  TYPE(inversion) :: inv

  INTEGER :: iteration_number
  ! REAL :: lbfgs_conv_crit

  !MPI
  INTEGER :: infompi


  !lbfgs_conv_crit=1e7
  !lbfgs_conv_crit=1e3
  m=inv%lbfgs_m

  !******************
  ! FWI with LBFGSB
  !******************

  iteration_number=0

  if(mype==0) then 
     iprint = 1

     !     We specify the tolerances in the stopping criteria.

     factr=inv%lbfgs_conv_crit!1.0d+12
     pgtol=1.0d-35


     !     We specify the dimension n of the sample problem and the number
     !        m of limited memory corrections stored.  (n and m should not
     !        exceed the limits nmax and mmax respectively.)
     n=pbdir%n1*pbdir%n2*inv%npar


     allocate(nbd(n), iwa(3*n))
     allocate(x(n), l(n), u(n), g(n), wa(2*m*n+5*n+11*m*m+8*m))

     !     We now provide nbd which defines the bounds on the variables:
     !                    l   specifies the lower bounds,
     !                    u   specifies the upper bounds. 

     !     First set bounds on the odd-numbered variables.

     CALL sub_setbound(n,l,u,nbd,pbdir,inv)

     !     We now define the starting point.
     CALL sub_setmodel2x(n,x,pbdir,inv,inv%model)

     ALLOCATE(inv%cc1(acqui%nsrc,pbdir%nw))
     inv%cc1(:,:)=1.

     write (6,16)
16   format(/,5x, 'Solving sample problem.',/,5x, ' (f = 0.0 at the optimal solution.)',/)               

     !     We start the iteration by initializing task.
     ! 
     task = 'START'

     !        ------- the beginning of the loop ----------
  END IF

111 continue

  !     This is the call to the L-BFGS-B code.
  if(mype==0) THEN
     write(*,*)'in lbfgs'
     call setulb(n,m,x,l,u,nbd,f,g,factr,pgtol,wa,iwa,task,iprint,&
          csave,lsave,isave,dsave)
     write(*,*)'out lbfgs'
  END IF

  !BCAST task
  CALL MPI_BCAST(task,60, MPI_CHARACTER,0,MPI_COMM_WORLD,infompi)
  if (task(1:2) .eq. 'FG') then
     !        the minimization routine has returned to request the
     !c        function f and gradient g values at the current x.

     !c        Compute function value f for the sample problem.

     if(mype==0)   THEN
        CALL sub_setx2model(n,x,pbdir,inv,inv%model)
        CALL sub_modelinv2pbdir(pbdir,inv)
     end if

     CALL sub_modeling_FWI(pbdir,acqui,inv)

     if(mype==0)  then
        write(*,*)'inv%fcost',inv%fcost
        !CALL MPI_BARRIER(MPI_COMM_WORLD,infompi)
        !stop
        f=inv%fcost
        CALL sub_setmodel2x(n,g,pbdir,inv,inv%gradient)

        !          go back to the minimization routine.

     endif

     CALL MPI_BARRIER(MPI_COMM_WORLD,infompi)
     !call sub_mpi_abort()
     goto 111

  elseif (task(1:5) .eq. 'NEW_X')  THEN
     if(mype==0) THEN
        iteration_number=iteration_number+1
        WRITE(*,*)'iteration number',iteration_number
        CALL sub_setx2model(n,x,pbdir,inv,inv%model)
        !        CALL sub_unprecond(inv,pbdir)
        CALL sub_modelinv2pbdirforwrinting(pbdir,inv,iteration_number,0)
        open(71,file='invparinter',access='direct',recl=4*n)
        write(71,rec=iteration_number)inv%model(:,:,1:inv%npar)
        close(71)

     endif

     CALL MPI_BARRIER(MPI_COMM_WORLD,infompi)
     goto 111

  END IF
  !        the minimization routine has returned with a new iterate,
  !         and we have opted to continue the iteration.

  !           ---------- the end of the loop -------------

  !     If task is neither FG nor NEW_X we terminate execution.


  !EXIT ->OUTPUT TH FINAL MODEL
  if(mype==0) THEN
     CALL sub_setx2model(n,x,pbdir,inv,inv%model)
     !     CALL sub_unprecond(inv,pbdir)
     CALL sub_modelinv2pbdirforwrinting(pbdir,inv,iteration_number,1)
     open(71,file='invparfinal',access='direct',recl=4*n)
     write(71,rec=1)inv%model(:,:,1:inv%npar)
     close(71)
  endif


end subroutine sub_FWI_LBFGS

subroutine sub_setbound(n,l,u,nbd,pbdir,inv)

  IMPLICIT NONE

  include "pbdirect.h"
  include "inversion.h"


  !PBDIRECT
  TYPE (pbdirect) :: pbdir
  !INVERSION
  TYPE (inversion) :: inv

  integer :: n
  real,dimension(n) :: l,u
  integer ,dimension(n) :: nbd

  INTEGER :: k,ipar,i1,i2

  CALL sub_setlimit_par(pbdir,inv)

  k=1
  DO ipar=1,inv%npar
     DO i2=1,pbdir%n2
        DO i1=1,pbdir%n1
           l(k)=inv%limit_par(ipar,1)
           u(k)=inv%limit_par(ipar,2)
           nbd(k)=2
           k=k+1     
        END DO
     END DO
  END DO

end subroutine sub_setbound

subroutine  sub_setmodel2x(n,x,pbdir,inv,model)

  IMPLICIT NONE

  include "pbdirect.h"
  include "inversion.h"


  !PBDIRECT
  TYPE (pbdirect) :: pbdir
  !INVERSION
  TYPE (inversion) :: inv

  REAL,DIMENSION(pbdir%n1,pbdir%n2,inv%npar) :: model

  integer :: n
  real,dimension(n) :: x

  integer :: k,i1,i2,ipar

  k=1
  DO ipar=1,inv%npar
     DO i2=1,pbdir%n2
        DO i1=1,pbdir%n1
           x(k)=model(i1,i2,ipar)
           k=k+1
        END DO
     END DO
  END DO

end subroutine sub_setmodel2x

subroutine  sub_setx2model(n,x,pbdir,inv,model)

  IMPLICIT NONE

  include "pbdirect.h"
  include "inversion.h"


  !PBDIRECT
  TYPE (pbdirect) :: pbdir
  !INVERSION
  TYPE (inversion) :: inv
  REAL,DIMENSION(pbdir%n1,pbdir%n2,inv%npar) :: model


  integer :: n
  real,dimension(n) :: x

  integer :: k,i1,i2,ipar

  k=1
  DO ipar=1,inv%npar
     DO i2=1,pbdir%n2
        DO i1=1,pbdir%n1
           model(i1,i2,ipar)=x(k)
           k=k+1
        END DO
     END DO
  END DO

end subroutine sub_setx2model

subroutine sub_modelinv2pbdir(pbdir,inv)

  IMPLICIT NONE

  include "pbdirect.h"
  include "inversion.h"


  !PBDIRECT
  TYPE (pbdirect) :: pbdir
  !INVERSION
  TYPE (inversion) :: inv


  integer :: i1,i2,iipar

  ! if(inv%lambda < 0. ) CALL sub_precond(inv,inv%model,pbdir%n1,pbdir%n2,inv%npar,1)

  DO iipar=1,inv%npar0!inv%npar
     if(inv%invfamily==1) THEN
        if(inv%invpar(iipar)==1) THEN
           !VP in family 1
           DO i2=1,pbdir%n2
              DO i1=1,pbdir%n1
                 pbdir%vp(i1,i2)=inv%model(i1,i2,iipar)
              END DO
           END DO
        ELSEIF(inv%invpar(iipar)==2) THEN
           !RHO in family 1
           DO i2=1,pbdir%n2
              DO i1=1,pbdir%n1
                 pbdir%rho(i1,i2)=inv%model(i1,i2,iipar)
              END DO
           END DO
        ELSEIF(inv%invpar(iipar)==11.OR.inv%invpar(iipar)==10 ) THEN
           !LOG(VP) in family 1
           DO i2=1,pbdir%n2
              DO i1=1,pbdir%n1
                 pbdir%vp(i1,i2)=EXP(inv%model(i1,i2,iipar))
              END DO
           END DO
        ELSEIF(inv%invpar(iipar)==12) THEN
           !LOG(RHO) in family 1
           DO i2=1,pbdir%n2
              DO i1=1,pbdir%n1
                 pbdir%rho(i1,i2)=EXP(inv%model(i1,i2,iipar))
              END DO
           END DO
        ELSEIF(inv%invpar(iipar)==31 ) THEN
           !LOG(VP/VP_0) in family 1
           DO i2=1,pbdir%n2
              DO i1=1,pbdir%n1
                 pbdir%vp(i1,i2)=inv%refvalue(iipar)*EXP(inv%model(i1,i2,iipar))
              END DO
           END DO
           WRITE(*,*)'pbdir%vp(1,1)',pbdir%vp(1,1)
        ELSEIF(inv%invpar(iipar)==32) THEN
           !LOG(RHO/RHO_0) in family 1
           DO i2=1,pbdir%n2
              DO i1=1,pbdir%n1
                 pbdir%rho(i1,i2)=inv%refvalue(iipar)*EXP(inv%model(i1,i2,iipar))
              END DO
           END DO
        END IF
     END IF
     if(inv%invfamily==2) THEN
        if(inv%invpar(iipar)==1) THEN
           !VP in family 2
           DO i2=1,pbdir%n2
              DO i1=1,pbdir%n1
                 pbdir%vp(i1,i2)=inv%model(i1,i2,iipar)
              END DO
           END DO
        ELSEIF(inv%invpar(iipar)==2) THEN
           !IP in family 2
           DO i2=1,pbdir%n2
              DO i1=1,pbdir%n1
                 pbdir%rho(i1,i2)=inv%model(i1,i2,iipar)/pbdir%vp(i1,i2)
              END DO
           END DO
        END IF
     END IF
     if(inv%invfamily==1 .OR. inv%invfamily==2) THEN
        IF(inv%invpar(iipar)==3) THEN
           !QP in family 1
           DO i2=1,pbdir%n2
              DO i1=1,pbdir%n1
                 pbdir%qp(i1,i2)=inv%model(i1,i2,iipar)
              END DO
           END DO
        ELSEIF(inv%invpar(iipar)==4) THEN
           !EPSILON in family 1
           DO i2=1,pbdir%n2
              DO i1=1,pbdir%n1
                 pbdir%epsilon(i1,i2)=inv%model(i1,i2,iipar)
              END DO
           END DO
        ELSEIF(inv%invpar(iipar)==5) THEN
           !DELTA in family 1
           DO i2=1,pbdir%n2
              DO i1=1,pbdir%n1
                 pbdir%delta(i1,i2)=inv%model(i1,i2,iipar)
              END DO
           END DO
        ELSEIF(inv%invpar(iipar)==6) THEN
           !DELTA in family 1
           DO i2=1,pbdir%n2
              DO i1=1,pbdir%n1
                 pbdir%theta(i1,i2)=inv%model(i1,i2,iipar)
              END DO
           END DO
!!$        ELSEIF(inv%invpar(iipar)==21) THEN
!!$           !VP/Vp_min in family 1
!!$           DO i2=1,pbdir%n2
!!$              DO i1=1,pbdir%n1
!!$                 pbdir%vp(i1,i2)=inv%model(i1,i2,iipar)*inv%limit_par(iipar,1)
!!$              END DO
!!$           END DO
!!$        ELSEIF(inv%invpar(iipar)==22) THEN
!!$           !RHO/rho_min in family 1
!!$           DO i2=1,pbdir%n2
!!$              DO i1=1,pbdir%n1
!!$                 pbdir%rho(i1,i2)=inv%model(i1,i2,iipar)*inv%limit_par(iipar,1)
!!$              END DO
!!$           END DO
        ELSEIF(inv%invpar(iipar)==33) THEN
           !LOG(QP/QP_0) in family 1
           DO i2=1,pbdir%n2
              DO i1=1,pbdir%n1
                 pbdir%qp(i1,i2)=inv%refvalue(iipar)*EXP(inv%model(i1,i2,iipar))
              END DO
           END DO
        ELSEIF(inv%invpar(iipar)==34) THEN
           !LOG(1+EPSILON) in family 1
           DO i2=1,pbdir%n2
              DO i1=1,pbdir%n1
                 pbdir%epsilon(i1,i2)=EXP(inv%model(i1,i2,iipar))-1.
              END DO
           END DO
        ELSEIF(inv%invpar(iipar)==35) THEN
           !LOG(1+DELTA) in family 1
           DO i2=1,pbdir%n2
              DO i1=1,pbdir%n1
                 pbdir%delta(i1,i2)=EXP(inv%model(i1,i2,iipar))-1.
              END DO
           END DO
        end if
     end if
  END DO

end subroutine sub_modelinv2pbdir


subroutine sub_modelinv2pbdirforwrinting(pbdir,inv,iteration_number,mode)

  IMPLICIT NONE

  include "pbdirect.h"
  include "inversion.h"


  !PBDIRECT
  TYPE (pbdirect) :: pbdir
  !INVERSION
  TYPE (inversion) :: inv
  INTEGER :: iteration_number,mode

  integer :: i1,i2,iipar

  CALL sub_modelinv2pbdir(pbdir,inv)

  IF(mode==0) THEN !intermediate writing

     OPEN(41,file='param_vp_inter',access='direct',recl=4*pbdir%n1*pbdir%n2)
     WRITE(41,rec=iteration_number)pbdir%vp(:,:)
     CLOSE(41)

     OPEN(42,file='param_qp_inter',access='direct',recl=4*pbdir%n1*pbdir%n2)
     WRITE(42,rec=iteration_number)pbdir%qp(:,:)
     CLOSE(42)

     OPEN(43,file='param_rho_inter',access='direct',recl=4*pbdir%n1*pbdir%n2)
     WRITE(43,rec=iteration_number)pbdir%rho(:,:)
     CLOSE(43)

     IF(pbdir%tool==2) THEN

        OPEN(44,file='param_epsilon_inter',access='direct',recl=4*pbdir%n1*pbdir%n2)
        WRITE(44,rec=iteration_number)pbdir%epsilon(:,:)
        CLOSE(44)

        OPEN(45,file='param_delta_inter',access='direct',recl=4*pbdir%n1*pbdir%n2)
        WRITE(45,rec=iteration_number)pbdir%delta(:,:)
        CLOSE(45)

        OPEN(46,file='param_theta_inter',access='direct',recl=4*pbdir%n1*pbdir%n2)
        WRITE(46,rec=iteration_number)pbdir%theta(:,:)
        CLOSE(46)

     END IF
  ELSE !final writing

     OPEN(41,file='param_vp_final',access='direct',recl=4*pbdir%n1*pbdir%n2)
     WRITE(41,rec=1)pbdir%vp(:,:)
     CLOSE(41)

     OPEN(42,file='param_qp_final',access='direct',recl=4*pbdir%n1*pbdir%n2)
     WRITE(42,rec=1)pbdir%qp(:,:)
     CLOSE(42)

     OPEN(43,file='param_rho_final',access='direct',recl=4*pbdir%n1*pbdir%n2)
     WRITE(43,rec=1)pbdir%rho(:,:)
     CLOSE(43)

     IF(pbdir%tool==2) THEN

        OPEN(44,file='param_epsilon_final',access='direct',recl=4*pbdir%n1*pbdir%n2)
        WRITE(44,rec=1)pbdir%epsilon(:,:)
        CLOSE(44)

        OPEN(45,file='param_delta_final',access='direct',recl=4*pbdir%n1*pbdir%n2)
        WRITE(45,rec=1)pbdir%delta(:,:)
        CLOSE(45)

        OPEN(46,file='param_theta_final',access='direct',recl=4*pbdir%n1*pbdir%n2)
        WRITE(46,rec=1)pbdir%theta(:,:)
        CLOSE(46)

     END IF

  END IF

end subroutine sub_modelinv2pbdirforwrinting
subroutine sub_setmodelpbdir2inv(pbdir,inv)
  IMPLICIT NONE

  include "pbdirect.h"
  include "inversion.h"


  !PBDIRECT
  TYPE (pbdirect) :: pbdir
  !INVERSION
  TYPE (inversion) :: inv
  integer :: i1,i2,iipar




  do iipar=1,inv%npar0
     if(inv%invfamily==1) THEN

        IF(inv%invpar(iipar)==1) THEN
           !VP
           inv%model(:,:,iipar)=pbdir%vp(:,:)
        ELSEIF(inv%invpar(iipar)==2) THEN
           !RHO
           inv%model(:,:,iipar)=pbdir%rho(:,:)
        ELSEIF(inv%invpar(iipar)==11 .OR. inv%invpar(iipar)==10) THEN
           !LOG(VP)
           inv%model(:,:,iipar)=LOG(pbdir%vp(:,:))
        ELSEIF(inv%invpar(iipar)==12) THEN
           !LOG(RHO)
           inv%model(:,:,iipar)=LOG(pbdir%rho(:,:))
        ELSEIF(inv%invpar(iipar)==31) THEN
           !LOG(VP/VP_0)
           !    inv%refvalue(iipar)=SUM(pbdir%vp(:,:))/(pbdir%n1*pbdir%n2)
           inv%refvalue(iipar)=inv%limit_pbdir_par(1,1)!MINVAL(pbdir%vp(:,:))
           inv%model(:,:,iipar)=LOG(pbdir%vp(:,:)/inv%limit_pbdir_par(1,1))
           WRITE(*,*)'inv%refvalue(iipar)',inv%refvalue(iipar)
           WRITE(*,*)'inv%model(1,1,iipar)',inv%model(1,1,iipar)
        ELSEIF(inv%invpar(iipar)==32) THEN
           !LOG(RHO/RHO_0)
           !  inv%refvalue(iipar)=SUM(pbdir%rho(:,:))/(pbdir%n1*pbdir%n2)
           inv%refvalue(iipar)=inv%limit_pbdir_par(3,1)!MINVAL(pbdir%rho(:,:))
           inv%model(:,:,iipar)=LOG(pbdir%rho(:,:)/inv%refvalue(iipar))
           WRITE(*,*)'inv%refvalue(iipar)',inv%refvalue(iipar)
           WRITE(*,*)'inv%model(1,1,iipar)',inv%model(1,1,iipar)
        END IF

     END IF
     if(inv%invfamily==2) THEN
        IF(inv%invpar(iipar)==1) THEN
           !VP
           inv%model(:,:,iipar)=pbdir%vp(:,:)
        ELSEIF(inv%invpar(iipar)==2) THEN
           !IP
           inv%model(:,:,iipar)=pbdir%rho(:,:)*pbdir%vp(:,:)
        END IF
     END IF
     if(inv%invfamily==1 .OR. inv%invfamily==2) THEN
        IF(inv%invpar(iipar)==3) THEN
           !QP
           inv%model(:,:,iipar)=pbdir%qp(:,:)
        ELSEIF(inv%invpar(iipar)==4) THEN
           !EPSILON
           inv%model(:,:,iipar)=pbdir%epsilon(:,:)
        ELSEIF(inv%invpar(iipar)==5) THEN
           !DELTA
           inv%model(:,:,iipar)=pbdir%delta(:,:)
        ELSEIF(inv%invpar(iipar)==6) THEN
           !THETA
           inv%model(:,:,iipar)=pbdir%theta(:,:)
!!$        ELSEIF(inv%invpar(iipar)==21) THEN
!!$           !VP/vp_min
!!$           inv%model(:,:,iipar)=pbdir%vp(:,:)/inv%limit_par(iipar,1)
!!$        ELSEIF(inv%invpar(iipar)==22) THEN
!!$           !RHO/rho_min
!!$           inv%model(:,:,iipar)=pbdir%rho(:,:)/inv%limit_par(iipar,1)
        ELSEIF(inv%invpar(iipar)==33) THEN
           !LOG(QP/QP_0)
           inv%refvalue(iipar)=inv%limit_pbdir_par(2,1)!MINVAL(pbdir%qp(:,:))
           ! inv%refvalue(iipar)=SUM(pbdir%qp(:,:))/(pbdir%n1*pbdir%n2)
           inv%model(:,:,iipar)=LOG(pbdir%qp(:,:)/inv%refvalue(iipar))
           WRITE(*,*)'inv%refvalue(iipar)',inv%refvalue(iipar)
           WRITE(*,*)'inv%model(1,1,iipar)',inv%model(1,1,iipar)
        ELSEIF(inv%invpar(iipar)==34) THEN
           !LOG(1+EPSILON)
           inv%model(:,:,iipar)=LOG(1+pbdir%epsilon(:,:))
           WRITE(*,*)'inv%model(1,1,iipar)',inv%model(1,1,iipar)
        ELSEIF(inv%invpar(iipar)==35) THEN
           !LOG(1+delta)
           inv%model(:,:,iipar)=LOG(1+pbdir%delta(:,:))
           WRITE(*,*)'inv%model(1,1,iipar)',inv%model(1,1,iipar)
        END IF
     END IF

  END DO



  ! if(inv%lambda < 0. ) CALL sub_precond(inv,inv%model,pbdir%n1,pbdir%n2,inv%npar,-1)

END SUBROUTINE sub_setmodelpbdir2inv

subroutine sub_setlimit_par(pbdir,inv)

  IMPLICIT NONE

  include "pbdirect.h"
  include "inversion.h"


  !PBDIRECT
  TYPE (pbdirect) :: pbdir
  !INVERSION
  TYPE (inversion) :: inv


  integer :: i1,i2,iipar

  ! if(inv%lambda < 0. ) CALL sub_precond(inv,inv%model,pbdir%n1,pbdir%n2,inv%npar,1)
  ALLOCATE(inv%limit_par(inv%npar,2))

  DO iipar=1,inv%npar
     if(inv%invfamily==1) THEN
        if(inv%invpar(iipar)==1) THEN
           !VP in family 1
           inv%limit_par(iipar,1)=inv%limit_pbdir_par(1,1)
           inv%limit_par(iipar,2)=inv%limit_pbdir_par(1,2)
        ELSEIF(inv%invpar(iipar)==2) THEN
           !RHO in family 1
           inv%limit_par(iipar,1)=inv%limit_pbdir_par(3,1)
           inv%limit_par(iipar,2)=inv%limit_pbdir_par(3,2) 
        ELSEIF(inv%invpar(iipar)==11.OR.inv%invpar(iipar)==10 ) THEN
           !LOG(VP) in family 1
           inv%limit_par(iipar,1)=LOG(inv%limit_pbdir_par(1,1))
           inv%limit_par(iipar,2)=LOG(inv%limit_pbdir_par(1,2)) 
        ELSEIF(inv%invpar(iipar)==12) THEN
           !LOG(RHO) in family 1
           inv%limit_par(iipar,1)=LOG(inv%limit_pbdir_par(3,1))
           inv%limit_par(iipar,2)=LOG(inv%limit_pbdir_par(3,2)) 
        ELSEIF(inv%invpar(iipar)==31 ) THEN
           !LOG(VP/VP_0) in family 1
           inv%limit_par(iipar,1)=LOG(inv%limit_pbdir_par(1,1)/inv%refvalue(iipar))
           inv%limit_par(iipar,2)=LOG(inv%limit_pbdir_par(1,2)/inv%refvalue(iipar))
        ELSEIF(inv%invpar(iipar)==32) THEN
           !LOG(RHO/RHO_0) in family 1
           inv%limit_par(iipar,1)=LOG(inv%limit_pbdir_par(3,1)/inv%refvalue(iipar))
           inv%limit_par(iipar,2)=LOG(inv%limit_pbdir_par(3,2)/inv%refvalue(iipar))
        end if
     end if
     if(inv%invfamily==2) THEN
        if(inv%invpar(iipar)==1) THEN
           !VP in family 2
           inv%limit_par(iipar,1)=inv%limit_pbdir_par(1,1)
           inv%limit_par(iipar,2)=inv%limit_pbdir_par(1,2)
        ELSEIF(inv%invpar(iipar)==2) THEN
           !IP in family 2
           inv%limit_par(iipar,1)=inv%limit_pbdir_par(1,1)*inv%limit_pbdir_par(3,1)
           inv%limit_par(iipar,2)=inv%limit_pbdir_par(1,2)*inv%limit_pbdir_par(3,2)
        end if
     end if 
     if(inv%invfamily==1 .OR. inv%invfamily==2) THEN
        IF(inv%invpar(iipar)==3) THEN
           !QP in family 1&2
           inv%limit_par(iipar,1)=inv%limit_pbdir_par(2,1)
           inv%limit_par(iipar,2)=inv%limit_pbdir_par(2,2) 
        ELSEIF(inv%invpar(iipar)==4) THEN
           !EPSILON in family 1&2
           inv%limit_par(iipar,1)=inv%limit_pbdir_par(4,1)
           inv%limit_par(iipar,2)=inv%limit_pbdir_par(4,2) 
        ELSEIF(inv%invpar(iipar)==5) THEN
           !DELTA in family 1&2
           inv%limit_par(iipar,1)=inv%limit_pbdir_par(5,1)
           inv%limit_par(iipar,2)=inv%limit_pbdir_par(5,2) 


!!$        ELSEIF(inv%invpar(iipar)==21) THEN
!!$           !VP/Vp_min in family 1
!!$           DO i2=1,pbdir%n2
!!$              DO i1=1,pbdir%n1
!!$                 pbdir%vp(i1,i2)=inv%model(i1,i2,iipar)*inv%limit_par(iipar,1)
!!$              END DO
!!$           END DO
!!$        ELSEIF(inv%invpar(iipar)==22) THEN
!!$           !RHO/rho_min in family 1
!!$           DO i2=1,pbdir%n2
!!$              DO i1=1,pbdir%n1
!!$                 pbdir%rho(i1,i2)=inv%model(i1,i2,iipar)*inv%limit_par(iipar,1)
!!$              END DO
!!$           END DO
        ELSEIF(inv%invpar(iipar)==33) THEN
           !LOG(QP/QP_0) in family 1&2
           inv%limit_par(iipar,1)=LOG(inv%limit_pbdir_par(2,1)/inv%refvalue(iipar))
           inv%limit_par(iipar,2)=LOG(inv%limit_pbdir_par(2,2)/inv%refvalue(iipar))
        ELSEIF(inv%invpar(iipar)==34) THEN
           !LOG(1+EPSILON) in family 1&2
           inv%limit_par(iipar,1)=LOG(1+inv%limit_pbdir_par(4,1))
           inv%limit_par(iipar,2)=LOG(1+inv%limit_pbdir_par(4,2))
        ELSEIF(inv%invpar(iipar)==35) THEN
           !LOG(1+DELTA) in family 1&2
           inv%limit_par(iipar,1)=LOG(1+inv%limit_pbdir_par(5,1))
           inv%limit_par(iipar,2)=LOG(1+inv%limit_pbdir_par(5,2))    
        end if
     end if
     WRITE(*,*)'SET LIMITS ON INVERTED PARA',iipar,inv%limit_par(iipar,1),inv%limit_par(iipar,2)
  END DO

end subroutine sub_setlimit_par
