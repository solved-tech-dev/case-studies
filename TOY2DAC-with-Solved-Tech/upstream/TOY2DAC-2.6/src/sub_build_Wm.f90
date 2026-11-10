!Copyright 2008-2011 SEISCOPE project, All rights reserved.
!Copyright 2013-2021 SEISCOPEII project, All rights reserved.

SUBROUTINE sub_build_Wm_d2dx2(Wm_val,Wm_irn,Wm_jcn,nn,nf,inv,pbdir)
  IMPLICIT NONE

  include "pbdirect.h"
  include "inversion.h"

  !PBDIRECT
  TYPE (pbdirect) :: pbdir
  !INVERSION
  TYPE (inversion) :: inv
  INTEGER :: nn,nf  !size of table, real number of non-zero terms
  REAL,DIMENSION(nn) :: Wm_val
  INTEGER,DIMENSION(nn) ::Wm_irn,Wm_jcn

  INTEGER :: i1,i2,kii,kij,ii1,ii2
  REAL    :: val_coloc,val

  !INIT
  Wm_irn(:)=0
  Wm_jcn(:)=0
  Wm_val(:)=0.
  nf=0

  ! BUILD MATRIX IN SPARSE COORDINATE FORMAT
  DO i2=1,pbdir%n2
     DO i1=1,pbdir%n1

        kii=i1+(i2-1)*pbdir%n1
        !init coloc term
        val_coloc=0.

        !horizontal direction
        ii1=i1
        DO ii2=i2-1,i2+1,2
           val=-inv%lambda_x
           !CHECK IF VALUE EXIST OR SHOULD BE PUT ON THE COLOC TERM
           IF(ii2==0 .OR. ii2==pbdir%n2+1) THEN
              val_coloc=val_coloc+val
           ELSE
              kij=ii1+(ii2-1)*pbdir%n1
              nf=nf+1
              Wm_irn(nf)=kii
              Wm_jcn(nf)=kij
              Wm_val(nf)=val
           END IF
        END DO !ii2

        !vertical direction
        ii2=i2
        DO ii1=i1-1,i1+1,2
           val=-inv%lambda_z
           !CHECK IF VALUE EXIST OR SHOULD BE PUT ON THE COLOC TERM
           IF(ii1==0 .OR. ii1==pbdir%n1+1) THEN
              val_coloc=val_coloc+val
           ELSE
              kij=ii1+(ii2-1)*pbdir%n1
              nf=nf+1
              Wm_irn(nf)=kii
              Wm_jcn(nf)=kij
              Wm_val(nf)=val
           END IF
        END DO !ii1
        !colocation term
        val=val_coloc+2.*inv%lambda_z+2.*inv%lambda_x
        nf=nf+1
        Wm_irn(nf)=kii
        Wm_jcn(nf)=kii
        Wm_val(nf)=val
     END DO !i1
  END DO !i2

write(*,*)'Wm_irn(1:nf/100)',Wm_irn(1:nf/100)
write(*,*)'Wm_jcn(1:nf/100)',Wm_jcn(1:nf/100)
write(*,*)'Wm_val(1:nf/100)',Wm_val(1:nf/100)
write(*,*)'MAXVAL(Wm_irn)',MAXVAL(Wm_irn)
write(*,*)'MAXVAL(Wm_jcn)',MAXVAL(Wm_jcn)
write(*,*)'MAXVAL(Wm_val)',MAXVAL(Wm_val)


 END SUBROUTINE sub_build_Wm_d2dx2
