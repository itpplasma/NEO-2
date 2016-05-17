!
  !SUBROUTINE magdata_for_particles(phi,y,bhat,geodcu,h_phi,dlogbdphi)
  SUBROUTINE magdata_for_particles(phi,bhat,geodcu,h_phi,dlogbdphi)
!
!   y(1)                     - $R$
!   y(2)                     - $Z$
!   y(3)-y(5)                - $\nabla \psi$ 
!   y(6)                     - $\int \rd s / B$
!   y(7)                     - $\int \rd s |\nabla \psi| / B$
!   y(8)                     - $\int \rd s \br \cdot \nabla \psi / B$
!   y(9)                     - $\int \rd s B$
!   y(10)                    - $\int \rd s B^2$
!   y(11)                    - $\int \rd s$ y(7)
!   y(12)                    - $\int \rd s$ y(8)
!   y(13)                    - $\int \rd s$ y(9)
!   y(14)                    - $\int \rd s$ y(10)
!
  USE size_mod
  USE partpa_mod
  USE rk4_kin_mod, ONLY : y
  USE mag_interface_mod, ONLY: mag_coordinates,boozer_s

!
  IMPLICIT NONE
!
  DOUBLE PRECISION                           :: phi
  !DOUBLE PRECISION, DIMENSION(ndim0)         :: y,dery
  DOUBLE PRECISION, DIMENSION(ndim0)         :: dery
!
  ! Winny: ipmin moved to module
  INTEGER :: i,j,npassing ! ,ipmin
  DOUBLE PRECISION                 :: bmod,sqrtg
  DOUBLE PRECISION, DIMENSION(3)   :: x,bder,hcovar,hctrvr
  DOUBLE PRECISION, DIMENSION(3,3) :: hcoder,hctder
  DOUBLE PRECISION :: geodcu,bhat,h_phi,dlogbdphi
!
  IF (mag_coordinates .EQ. 0) THEN
     ! cylindrical
     x(1)=y(1)
     x(2)=phi
     x(3)=y(2)
  ELSE
     ! Boozer
     x(1)=boozer_s
     x(2)=phi
     x(3)=y(1)
  END IF
!
  CALL mag(x,bmod,sqrtg,bder,hcovar,hctrvr,hcoder,hctder)
  pardeb0=(hctrvr(1)*bder(1)+hctrvr(2)*bder(2)+hctrvr(3)*bder(3))/hctrvr(2)
!
  bhat=bmod/bmod0
!
! geodcu = $k_G |\nabla\psi|$
  geodcu=(hcovar(1)*bder(2)*y(5)        &
        + hcovar(2)*bder(3)*y(3)        &
        + hcovar(3)*bder(1)*y(4)        &
        - hcovar(1)*bder(3)*y(4)        &
        - hcovar(2)*bder(1)*y(5)        &
        - hcovar(3)*bder(2)*y(3))/sqrtg
!
  h_phi=hctrvr(2)
!
  dlogbdphi=(hctrvr(1)*bder(1)+hctrvr(2)*bder(2)+hctrvr(3)*bder(3))/hctrvr(2)
!
  !print *, y
  !print *, x,bmod,sqrtg,bder,hcovar,hctrvr,hcoder,hctder
  !print *, pardeb0,bhat,geodcu,h_phi,dlogbdphi 
  !pause 
  RETURN
  END
!
