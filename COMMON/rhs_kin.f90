!
SUBROUTINE rhs_kin(phi,y,dery)
  USE mag_interface_mod, ONLY: mag_coordinates,boozer_s
  USE size_mod
  USE partpa_mod
  !! Modifications by Andreas F. Martitsch (12.03.2014)
  ! boozer_curr_tor, boozer_curr_pol, boozer_psi_pr,
  ! boozer_sqrtg11 and boozer_isqrg are now converted
  ! to cgs-units within neo_magfie.
  ! This step requires changes within rhs_kin.f90 and
  ! ripple_solver.f90!
  USE neo_magfie_mod, ONLY: boozer_iota,boozer_curr_pol_hat,&
       boozer_curr_tor_hat,boozer_sqrtg11,boozer_psi_pr_hat,&
       boozer_isqrg
  !! End Modifications by Andreas F. Martitsch (12.03.2014)
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

  ! New Boozer
  !   y(1) = boozer_theta
  !   y(2) not used
  !   y(3) - y(14) unchanged at the moment has to be adapted


  !

  IMPLICIT NONE
  !
  DOUBLE PRECISION                           :: phi
  DOUBLE PRECISION, DIMENSION(ndim0)         :: y,dery
  !
  INTEGER :: i,j,npassing
  DOUBLE PRECISION                 :: bmod,sqrtg
  DOUBLE PRECISION, DIMENSION(3)   :: x,bder,hcovar,hctrvr
  DOUBLE PRECISION, DIMENSION(3,3) :: hcoder,hctder
  DOUBLE PRECISION :: pardeb,bhat,subsq,drive,derg
  DOUBLE PRECISION :: tdery,g11,sqrg11_tok
  DOUBLE PRECISION :: dery7_sergei
  DOUBLE PRECISION :: dery6_sergei
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
  !
  bhat=bmod/bmod0
  !
  pardeb=(hctrvr(1)*bder(1)+hctrvr(2)*bder(2)+hctrvr(3)*bder(3))/hctrvr(2)
  !
!!$  ! Winny 
!!$  IF(ipmin.EQ.0.AND.pardeb*pardeb0.LE.0.d0.AND.pardeb.GT.0.d0) THEN
!!$     ipmin = 1
!!$  ELSE
!!$     !ipmin = 0
!!$  END IF
!!$  IF(ipmax.EQ.0.AND.pardeb*pardeb0.LE.0.d0.AND.pardeb.LT.0.d0) THEN 
!!$     ipmax = 1
!!$  ELSE
!!$     !ipmax = 0
!!$  END IF
  pardeb0=pardeb
  !
  IF (mag_coordinates .EQ. 0) THEN
     ! cylindrical
     dery(1)=hctrvr(1)/hctrvr(2)
     dery(2)=hctrvr(3)/hctrvr(2)
     !
     DO j=3,5
        dery(j)=0.d0
        DO i=3,5
           dery(j)=dery(j)-hctder(j-2,i-2)*y(i)
        ENDDO
        dery(j)=dery(j)/hctrvr(2)
     ENDDO
     !
     dery(6)=1.d0/(bhat*hctrvr(2))
     dery(7)=dsqrt(y(3)**2+(y(4)/x(1))**2+y(5)**2)*dery(6)
     dery(8)=(y(1)*y(3)+y(2)*y(5))*dery(6)
     dery(9)=bhat**2*dery(6)
     dery(10)=bhat**3*dery(6)
     dery(11)=y(7)
     dery(12)=y(8)
     dery(13)=y(9)
     dery(14)=y(10)
  ELSE
     ! Boozer 
     dery(1)=boozer_iota
     dery(2)=0.0d0
     dery(3)=0.0d0
     dery(4)=0.0d0
     dery(5)=0.0d0
     !
     !dery(6)=1.d0/(bhat*hctrvr(2))
     !dery6_sergei = 1.d0/(bhat*hctrvr(2))
     !tdery  =(boozer_curr_pol + boozer_iota * boozer_curr_tor) / bhat**2

     !dery(6) = 1.d0 / boozer_isqrg
     
     !! Modifications by Andreas F. Martitsch (12.03.2014)
     ! boozer_curr_tor, boozer_curr_pol, boozer_psi_pr,
     ! boozer_sqrtg11 and boozer_isqrg are now converted
     ! to cgs-units within neo_magfie.
     ! This step requires changes within rhs_kin.f90 and
     ! ripple_solver.f90!
     !This was the old version (conversion to cgs is done here):
     !dery(6) = 1.d0 / boozer_isqrg * bmod0 * 1.0d2
     !Now the quantities are already converted within neo_magfie:
     dery(6) = 1.0d0/(boozer_isqrg*boozer_psi_pr_hat)
     !! End Modifications by Andreas F. Martitsch (12.03.2014)
     
     !tdery  =(boozer_curr_pol + boozer_iota * boozer_curr_tor) / bmod0 / bhat**2 * 1.0d2
     !tdery  =(boozer_curr_pol + boozer_iota * boozer_curr_tor) / bmod**2
     !PRINT *, 'bmod0            = ', bmod0 
     !PRINT *, 'bhat             = ', bhat
     !PRINT *, 'dery6            = ', dery(6)
     !PRINT *, 'tdery6           = ', tdery
     !PRINT *, 'dery6_sergei     = ', dery6_sergei
     !PRINT *, 'fac6             = ', dery(6) / tdery
     !PRINT *, 'boozer_curr_pol  = ', boozer_curr_pol
     !PRINT *, 'boozer_curr_tor  = ', boozer_curr_tor
     !PRINT *, 'boozer_iota      = ', boozer_iota
     !PRINT *, 'boozer_s         = ', boozer_s
     !PRINT *, 'boozer_sqrtg11   = ', boozer_sqrtg11
     !PRINT *, 'boozer_psi_pr    = ', boozer_psi_pr
     !PRINT *, 'boozer_isqrg     = ', boozer_isqrg

     !dery(7)=dery(6) * boozer_sqrtg11 / boozer_psi_pr * 1.d-2 / bmod0


     !! Modifications by Andreas F. Martitsch (12.03.2014)
     ! boozer_curr_tor, boozer_curr_pol, boozer_psi_pr,
     ! boozer_sqrtg11 and boozer_isqrg are now converted
     ! to cgs-units within neo_magfie.
     ! This step requires changes within rhs_kin.f90 and
     ! ripple_solver.f90!
     !This was the old version (conversion to cgs is done here):
     !dery(7)= dery(6) * boozer_sqrtg11 / boozer_psi_pr / 1.d2
     !Now the quantities are already converted within neo_magfie:
     dery(7) = dery(6)*boozer_sqrtg11
     !PRINT *,'dery(6),dery(7) :',dery(6),dery(7)
     !! End Modifications by Andreas F. Martitsch (12.03.2014)
     
     !g11=4.d0*boozer_s/30.d0**2                                           !<-in
     !dery7_sergei = dsqrt(y(3)**2*g11)*dery6_sergei                                !<-in
     !sqrg11_tok = 2.0d0 * SQRT(boozer_s) / 3.d-1 ! * boozer_psi_pr
     !tdery =   dery(6) * sqrg11_tok
     !PRINT *, 'sqrg11_tok       = ', sqrg11_tok
     !PRINT *, 'dery7            = ', dery(7)
     !PRINT *, 'tdery7           = ', tdery
     !PRINT *, 'dery7_sergei     = ', dery7_sergei
     !PRINT *, 'fac7             = ', dery(7) / tdery
     !PRINT *, ' '
     !PRINT *, 'd6/d7            = ', dery(6) / dery(7)
     !PRINT *, '(d6/d7)_sergei   = ', dery6_sergei / dery7_sergei
     !PRINT *, 'd6/d6_sergei     = ', dery(6) / dery6_sergei
     !PRINT *, 'd7/d7_sergei     = ', dery(7) / dery7_sergei
     !PRINT *, '(d6/d7)/(d6/d7)s = ', dery(6) / dery(7) / (dery6_sergei / dery7_sergei)





     dery(8)=0.0d0
     dery(9)=bhat**2*dery(6)
     dery(10)=bhat**3*dery(6)
     dery(11)=y(7)
     dery(12)=y(8)
     dery(13)=y(9)
     dery(14)=y(10)    
     !print *, dery
  END IF
  !
  RETURN
END SUBROUTINE rhs_kin
!
