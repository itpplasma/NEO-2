!
!SUBROUTINE magdata_for_particles(phi,y,bhat,geodcu,h_phi,dlogbdphi)
SUBROUTINE magdata_for_particles(phi,bhat,geodcu,h_phi,dlogbdphi,&
     bcovar_s_hat,dlogbds,dbcovar_s_hat_dphi,bnoverb0)
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
  !! Modifications by Andreas F. Martitsch (09.03.2014)
  ! Collection of subroutines (mag.f90) converted to a module.
  ! This allows one to make use of generic interfaces 
  ! for module procedures (e.g., for a different number of
  ! arguments).
  ! Note: This requires changes in "flint_prepare" and 
  ! "write_volume_data" (both flint.f90), and "rhs_kin"
  ! and "magdata_for_particles". 
  USE mag_sub, ONLY: mag
  !! End Modifications by Andreas F. Martitsch (09.03.2014)
  !! Modifications by Andreas F. Martitsch (11.06.2014)
  ! Addtional output needed for modeling the plasma rotation - 
  ! Routine for the direct computation of the perturbation
  ! field from the Boozer file
  USE neo_magfie_perturbation, ONLY : neo_magfie_pert
  !! End Modifications by Andreas F. Martitsch (11.06.2014)
  !
  IMPLICIT NONE
  !
  !! Modifications by Andreas F. Martitsch (11.03.2014)
  ! Optional output (necessary for modeling the magnetic rotation)
  DOUBLE PRECISION, INTENT(in)            :: phi
  DOUBLE PRECISION, INTENT(out)           :: geodcu,bhat,h_phi,dlogbdphi
  DOUBLE PRECISION, OPTIONAL, INTENT(out) :: bcovar_s_hat, &
       dlogbds, dbcovar_s_hat_dphi
  ! Variable "y" already specified by rk4_kin_mod
  !DOUBLE PRECISION, DIMENSION(ndim0)         :: y,dery
  ! Variable "dery" seems to be unused at the moment
  ! (also not specified through any used module)
  !DOUBLE PRECISION, DIMENSION(ndim0)         :: dery
  !! End Modifications by Andreas F. Martitsch (11.03.2014)
  !! Modifications by Andreas F. Martitsch (11.06.2014)
  ! Addtional output needed for modeling the plasma rotation -
  DOUBLE COMPLEX, OPTIONAL, INTENT(out) :: bnoverb0
  !! End Modifications by Andreas F. Martitsch (11.06.2014)
  !
  ! Local definitions:
  ! Winny: ipmin moved to module
  INTEGER :: i,j,npassing ! ,ipmin
  DOUBLE PRECISION                 :: bmod,sqrtg
  DOUBLE PRECISION, DIMENSION(3)   :: x,bder,hcovar,hctrvr,bcovar_s_hat_der
  DOUBLE PRECISION, DIMENSION(3,3) :: hcoder,hctder
  DOUBLE COMPLEX                   :: bn_hat_pert
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
  !! Modifications by Andreas F. Martitsch (11.03.2014)
  ! Optional output (necessary for modeling the magnetic rotation)
  !
  ! Extra output from "mag" needed for the computation of dbcovar_s_dphi
  IF (PRESENT(dbcovar_s_hat_dphi)) THEN
     CALL mag(x,bmod,sqrtg,bder,hcovar,hctrvr,hcoder,hctder,bcovar_s_hat_der)
     ! Compute the derivative along the field line (phi_mfl) of the covariant radial
     ! B-field component. Note: In case of Boozer coordinates and for an
     ! axisymmetric B-field, this derivative reduces to a derivative over phi_mfl. 
     ! [In ripple_solver we work with field-aligned Boozer coordinates -
     ! ($\varphi_0$ = field line label,$\varphi_s$ = field line parameter)]
     dbcovar_s_hat_dphi=(hctrvr(1)*bcovar_s_hat_der(1) + &
          hctrvr(2)*bcovar_s_hat_der(2)  + &
          hctrvr(3)*bcovar_s_hat_der(3)) / &
          hctrvr(2)
     !PRINT *,'------------'
     !PRINT *,'dbcovar_s_hat_dphi:        ',dbcovar_s_hat_dphi
     !PRINT *,'dbcovar_s_hat_dphi(check): ',hctrvr(3)*bcovar_s_hat_der(3)/hctrvr(2)
     !PRINT *,'------------'
  ELSE ! otherwise use the old "mag"-routine
     IF (PRESENT(dlogbds) .OR. PRESENT(bcovar_s_hat)) THEN
        ! All three optional quantities are required for the 
        ! computation of the magnetic rotation (-> There
        ! might be an error, if one component is missing)
        STOP "Exit magdata_for_particles: optional argument dbcovar_s_dphi missing"
     END IF
     CALL mag(x,bmod,sqrtg,bder,hcovar,hctrvr,hcoder,hctder)
  END IF
  !
  ! Compute the covariant radial B-field component
  IF (PRESENT(bcovar_s_hat)) THEN
     bcovar_s_hat=hcovar(1)*(bmod/bmod0)
  ELSE
     IF (PRESENT(dbcovar_s_hat_dphi) .OR. PRESENT(dlogbds)) THEN
        ! All three optional quantities are required for the 
        ! computation of the magnetic rotation (-> There
        ! might be an error, if one component is missing)
        STOP "Exit magdata_for_particles: optional argument bcovar_s missing"
     END IF
  END IF
  !
  ! Compute $\frac{\partial\log{B}}{\partial s}$
  IF (PRESENT(dlogbds)) THEN
     dlogbds=bder(1)
  ELSE
     IF (PRESENT(dbcovar_s_hat_dphi) .OR. PRESENT(bcovar_s_hat)) THEN
        ! All three optional quantities are required for the 
        ! computation of the magnetic rotation (-> There
        ! might be an error, if one component is missing)
        STOP "Exit magdata_for_particles: optional dlogbds argument  missing"
     END IF
  END IF
  !
  !
  !! End Modifications by Andreas F. Martitsch (11.03.2014)
  !
  pardeb0=(hctrvr(1)*bder(1)+hctrvr(2)*bder(2)+hctrvr(3)*bder(3))/hctrvr(2)
  !
  bhat=bmod/bmod0
  !
  !! Modifications by Andreas F. Martitsch (11.06.2014)
  ! Addtional output needed for modeling the plasma rotation -
  IF (PRESENT(bnoverb0) .AND. (mag_coordinates .NE. 0)) THEN
     CALL neo_magfie_pert(x,bn_hat_pert)
     bnoverb0=bn_hat_pert/bhat
     !PRINT *,'bnoverb0: ',bnoverb0
     !PRINT *,x(2),x(3)
  ELSE ! for testing artificial perturbation field can be created within ripple_solver
     IF (mag_coordinates .EQ. 0) THEN
        ! At the moment the quasilinear approach for modeling the plasma rotation
        ! is only implemented for B-field spectra in Boozer coordinates.
        ! Make sure that your Boozer file matches one of the specified
        ! formats (see neo_read_pert() within neo_magfie_perturbation.f90).
        STOP "Exit magdata_for_particles: Extraction of the perturbation&
             & field only implemented for Boozer files (at the moment)"
     END IF
  END IF
  !! End Modifications by Andreas F. Martitsch (11.06.2014)  
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
END SUBROUTINE magdata_for_particles
!
