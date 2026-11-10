!Copyright 2008-2011 SEISCOPE project, All rights reserved.
!Copyright 2013-2021 SEISCOPEII project, All rights reserved.
! ------------------------------------------------------------------------------------------------------
! subroutine subDadmi
! Build local derivative of impedance matrix for various 
! anisotropic media using the equation
! system of Zhou, SEG, 2006.
! Rotated staggered grid stencil is implemented.
! 
! version from S. Operto Modified by R. Brossier 10/2011

!
! b: buoyancy
! mu: complex bulk modulus
! delta: Thomsen parameter
! epsilon: Thomsen parameter
! teta: anisotropy tilt
! ianiso: type of anisotropy
! d1,d1b,d2,d2b: damping functions for PML
! mat,irn,icn: arrays for sparse impedance matrix storage
! n1e x n2e: grid dimensions
! h: grid spacing
! ne: number of non zero coefficients in the matrix
! omegac: complex frequency
!
! Copyright 2008-2011 SEISCOPE project, All rights reserved.
! 
! ---------------------------------------------------------------------------------------------------------


SUBROUTINE sub_dAdmi(i1e,i2e,b,mu1,delta,epsilon,teta,ianiso,&
     ne, irn,icn,mat,h,omegac,nnz,n1e,n2e,mode)

  IMPLICIT NONE

  integer ne,i2,i1,k,l,iparaxial,k0,i,i1e,i2e
  real h1,h,h2,h3
  real bpm,bmp,bpp,bmm,b00,bp0,bm0,b0p,b0m
  integer n1e,n2e,xi0,nnz,l0
  integer irn(ne),icn(ne)

  real ax,bx,cx,dx,ex,fx,gx,hx
  real az,bz,cz,dz,ez,fz,gz,hz,damp
  real wm1,wm2,wm3,w1,w2
  real b(i1e-2:i1e+2,i2e-2:i2e+2),pi,delta(i1e-1:i1e+1,i2e-1:i2e+1),epsilon(i1e-1:i1e+1,i2e-1:i2e+1),teta(i1e-1:i1e+1,i2e-1:i2e+1)

  complex d1p,d1m,d1qp,d1qm,d2p,d2m,d2qp,d2qm
  complex omegac,omega2,xc0,m,ci
  complex mu1(i1e-2:i1e+2,i2e-2:i2e+2),mat(ne)

  integer :: nne,ianiso,born,mode

  damp=0.
  ! Build anisotropic parameters
    wm1=1.
  w1=0.
  wm2=0.
  wm3=0.
  w2=1.-w1

  ! set constants
  xi0=0
  xc0=cmplx(0.,0.)
  pi=3.141592654
  omega2=omegac*omegac
  h1=1./h
  h2=1./(h*h)
  h3=1./(h*h*h)
  nne=n1e*n2e
  ci=cmplx(0.,1.)

  ! Initialize sparse impedance matrix
  irn(:)=0
  icn(:)=0
  mat(:)=cmplx(0.,0.)
  nnz=0
  k=0

  
  ! LOOP OVER ROWS OF THE IMPEDANCE MATRIX 
  if(mode==1) then
     !inversion of vp, epsilon, delta etc... build only at collocation point
     born=0
  elseif(mode==2) then
     !inversion of density -> build at 5 point arround collocation point
     born=1
  end if

  do i2=i2e-born,i2e+born
     do i1=i1e-born,i1e+born
        If((i2==i2e-1 .AND. i1==i1e) .OR. (i2==i2e+1 .AND. i1==i1e) .OR. (i2==i2e .AND. i1==i1e+1).OR. (i2==i2e .AND. i1==i1e-1).OR. (i2==i2e .AND. i1==i1e))THEN
           
           ax=1.+2.*delta(i1,i2)*cos(teta(i1,i2))*cos(teta(i1,i2))
           bx=-delta(i1,i2)*sin(2.*teta(i1,i2))
           cx=(1.+2.*delta(i1,i2))*cos(teta(i1,i2))*cos(teta(i1,i2))
           dx=-0.5*(1.+2.*delta(i1,i2))*sin(2.*teta(i1,i2))
           ex=2.*(epsilon(i1,i2)-delta(i1,i2))*cos(teta(i1,i2))*cos(teta(i1,i2))
           fx=-(epsilon(i1,i2)-delta(i1,i2))*sin(2.*teta(i1,i2))
           gx=2.*(epsilon(i1,i2)-delta(i1,i2))*cos(teta(i1,i2))*cos(teta(i1,i2))
           hx=-(epsilon(i1,i2)-delta(i1,i2))*sin(2.*teta(i1,i2))

           az=-delta(i1,i2)*sin(2.*teta(i1,i2))
           bz=1.+2.*delta(i1,i2)*sin(teta(i1,i2))*sin(teta(i1,i2))
           cz=-0.5*(1.+2.*delta(i1,i2))*sin(2.*teta(i1,i2))
           dz=(1.+2.*delta(i1,i2))*sin(teta(i1,i2))*sin(teta(i1,i2))
           ez=-(epsilon(i1,i2)-delta(i1,i2))*sin(2.*teta(i1,i2))
           fz=2.*(epsilon(i1,i2)-delta(i1,i2))*sin(teta(i1,i2))*sin(teta(i1,i2))
           gz=-(epsilon(i1,i2)-delta(i1,i2))*sin(2.*teta(i1,i2))
           hz=2.*(epsilon(i1,i2)-delta(i1,i2))*sin(teta(i1,i2))*sin(teta(i1,i2))
                      
           bp0=0.5*(b(i1,i2)+b(i1+1,i2))
           bm0=0.5*(b(i1,i2)+b(i1-1,i2))
           b0p=0.5*(b(i1,i2)+b(i1,i2+1))
           b0m=0.5*(b(i1,i2)+b(i1,i2-1))
                      
           ! Block matrix A          
           k=(i2-1)*n1e+i1
           l=k
           nnz=nnz+1
           irn(nnz)=k
           icn(nnz)=l
           mat(nnz)= wm1*omega2*mu1(i1,i2) &
                +w2*(-h2*ax*(b0p+b0m)-h2*bz*(bp0+bm0) )

           ! Node 01
           l=0
           IF (i2<n2e) then
              l=(i2-1+1)*n1e+i1
              nnz=nnz+1
              irn(nnz)=k
              icn(nnz)=l
              mat(nnz)= &
                   +w2*(h2*ax*b0p + 0.25*az*h2*(bp0-bm0) )              
           end if

           ! Node 0-1
           IF (i2>1) THEN
              l=(i2-1-1)*n1e+i1
              nnz=nnz+1
              irn(nnz)=k
              icn(nnz)=l
              mat(nnz)= &
                   +w2*(h2*ax*b0m + 0.25*az*h2*(-bp0+bm0) )
           end if
          
           ! Node -10 
           if (i1 > 1) then
              l=(i2-1)*n1e+i1-1
              nnz=nnz+1
              irn(nnz)=k
              icn(nnz)=l
              mat(nnz)= &
                   +w2*(h2*bz*bm0 + 0.25*bx*h2*(-b0p+b0m) )
           end if
                      
           ! Node 10           
           if (i1 < n1e) then
              l=(i2-1)*n1e+i1+1
              nnz=nnz+1
              irn(nnz)=k
              icn(nnz)=l
              mat(nnz)=   &
                   +w2*(h2*bz*bp0 + 0.25*bx*h2*(-b0m+b0p) )
           end if

           ! Block matrices B, C and D computed only for VTI/TTI media
           
           IF (ianiso.eq.4.or.ianiso.eq.5) THEN    !VTI OR TTI

              ! Block matrix B
              ! Node 00 -----------------------------------              
              l=(i2-1)*n1e+i1+nne
              nnz=nnz+1
              irn(nnz)=k
              icn(nnz)=l
              mat(nnz)= &
                   +w2*(-h2*cx*(b0p+b0m)-h2*dz*(bp0+bm0) )
              
              ! Node 01 
              l=0
              IF (i2<n2e) then
                 l=(i2-1+1)*n1e+i1+nne
                 nnz=nnz+1
                 irn(nnz)=k
                 icn(nnz)=l
                 mat(nnz)= &
                      +w2*(h2*cx*b0p + 0.25*cz*h2*(bp0-bm0) )
              end if


              ! Node 0-1 
              IF (i2>1) THEN
                 l=(i2-1-1)*n1e+i1+nne
                 nnz=nnz+1
                 irn(nnz)=k
                 icn(nnz)=l
                 mat(nnz)= &
                      +w2*(h2*cx*b0m + 0.25*cz*h2*(-bp0+bm0) )
              end if

              ! Node -10 
              if (i1 > 1) then
                 l=(i2-1)*n1e+i1-1+nne
                 nnz=nnz+1
                 irn(nnz)=k
                 icn(nnz)=l
                 mat(nnz)= &
                      +w2*(h2*dz*bm0 + 0.25*dx*h2*(-b0p+b0m) )                
              endif
                           
              ! Node 10              
              if (i1 < n1e) then
                 l=(i2-1)*n1e+i1+1+nne
                 nnz=nnz+1
                 irn(nnz)=k
                 icn(nnz)=l
                 mat(nnz)=   &
                      +w2*(h2*dz*bp0 + 0.25*dx*h2*(-b0m+b0p) )
              end if
              
              ! Block matrix C              
              ! COEFFICIENTS E, F
              ! Node 00 

              k=(i2-1)*n1e+i1+nne
              l=(i2-1)*n1e+i1
              nnz=nnz+1
              irn(nnz)=k
              icn(nnz)=l
              mat(nnz)=  &
                   +w2*(-h2*ex*(b0p+b0m)-h2*fz*(bp0+bm0) )


              ! Node 01 
              l=0
              IF (i2<n2e) then
                 l=(i2-1+1)*n1e+i1
                 nnz=nnz+1
                 irn(nnz)=k
                 icn(nnz)=l
                 mat(nnz)= &
                      +w2*(h2*ex*b0p + 0.25*ez*h2*(bp0-bm0) )
              end if
              
              ! Node 0-1 
              IF (i2>1) THEN
                 l=(i2-1-1)*n1e+i1
                 nnz=nnz+1
                 irn(nnz)=k
                 icn(nnz)=l
                 mat(nnz)= &
                      +w2*(h2*ex*b0m + 0.25*ez*h2*(-bp0+bm0) )
              end if
              
              ! Node -10               
              if (i1 > 1) then
                 l=(i2-1)*n1e+i1-1
                 nnz=nnz+1
                 irn(nnz)=k
                 icn(nnz)=l
                 mat(nnz)= &
                      +w2*(h2*fz*bm0 + 0.25*fx*h2*(-b0p+b0m) )
              end if

              ! Node 10
              if (i1 < n1e) then
                 l=(i2-1)*n1e+i1+1
                 nnz=nnz+1
                 irn(nnz)=k
                 icn(nnz)=l
                 mat(nnz)=   &
                      +w2*(h2*fz*bp0 + 0.25*fx*h2*(-b0m+b0p) )
              end if

              
              ! Block matrix D
              ! COEFFICIENTS G, H
              ! Node 00 
              
              k=(i2-1)*n1e+i1+nne
              l=(i2-1)*n1e+i1+nne
              nnz=nnz+1
              irn(nnz)=k
              icn(nnz)=l
              mat(nnz)= wm1*omega2*mu1(i1,i2) &
                   +w2*(-h2*gx*(b0p+b0m)-h2*hz*(bp0+bm0) )
              
              ! Node 01 
              l=0
              IF (i2<n2e) then
                 l=(i2-1+1)*n1e+i1+nne
                 nnz=nnz+1
                 irn(nnz)=k
                 icn(nnz)=l
                 mat(nnz)= &
                      +w2*(h2*gx*b0p + 0.25*gz*h2*(bp0-bm0) )
              end if
              
              ! Node 0-1              
              IF (i2>1) THEN
                 l=(i2-1-1)*n1e+i1+nne
                 nnz=nnz+1
                 irn(nnz)=k
                 icn(nnz)=l
                 mat(nnz)= &
                      +w2*(h2*gx*b0m + 0.25*gz*h2*(-bp0+bm0) )
              end if

              ! Node -10 
              if (i1 > 1) then
                 l=(i2-1)*n1e+i1-1+nne
                 nnz=nnz+1
                 irn(nnz)=k
                 icn(nnz)=l
                 mat(nnz)= &
                      +w2*(h2*hz*bm0 + 0.25*hx*h2*(-b0p+b0m) )                                  
              end if
                            
              ! Node 10              
              if (i1 < n1e) then
                 l=(i2-1)*n1e+i1+1+nne
                 nnz=nnz+1
                 irn(nnz)=k
                 icn(nnz)=l
                 mat(nnz)=   &
                      +w2*(h2*hz*bp0 + 0.25*hx*h2*(-b0m+b0p) )
              end if
           END IF
           
           
        END IF
     end do            !end do i1=1,n1
  end do            !end do i2=1,n2
end SUBROUTINE sub_dAdmi

subroutine sub_P_dAdmi(mat,irn,jcn,nnz,ne,vect1,size1,vect2,ivect2,ne2,nnz2)
  implicit none

  INTEGER :: nnz,ne,size1,ne2,nnz2
  COMPLEX :: mat(ne),vect1(size1),vect2(ne2)
  INTEGER :: irn(ne),jcn(ne),ivect2(ne2)
  INTEGER :: i,flag_store,icheck
  COMPLEX :: val
!write(*,*)'size1',size1
  nnz2=0
  DO i=1,nnz
     val=mat(i)*vect1(jcn(i))

     !CHECK IF A VALUE IS ALREADY STORED IN EMPLACEMENT vect2(irn(i))
     flag_store = 0
     DO icheck = nnz2 , 1 , -1
        IF (  ivect2 ( icheck ) == irn(i)) THEN ! VALUE IS STORED, WE ADD THE VALUE coef to the old one
           vect2( icheck ) = vect2( icheck )+ val
           flag_store = 1
           EXIT
        END IF
     END DO  !icheck

     IF ( flag_store == 0 ) THEN ! the value is not already stored, we store it now
        nnz2 = nnz2 + 1
        vect2 ( nnz2 ) = val
        ivect2( nnz2 ) = irn(i)
     END IF

  END DO

end subroutine sub_P_dAdmi

subroutine sub_P_dAdmi_P(vect1,size1,vect2,ivect2,ne2,nnz2,sol)
  implicit none

  INTEGER :: ne,size1,ne2,nnz2
  COMPLEX :: vect1(size1),vect2(ne2)
  INTEGER :: ivect2(ne2)
  INTEGER :: i
  COMPLEX :: sol
  sol=0.

  DO i=1,nnz2
     sol = sol + vect1(ivect2(i))*vect2(i)
  END DO

end subroutine sub_P_dAdmi_P
