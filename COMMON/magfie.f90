! Usage of program
!
!   USE magfie_mod, ONLY: magfie, magfie_deallocate, stevvo
!   USE neo_magfie_mod, ONLY: magfie_spline, magfie_sarray, magfie_result
!
! ATTENTION!!!!!!
!
! Standard outut of magfie produces arrays for SMT and BMC
!
! For other usage put magfie_result=1 (NEO)   
!
! EXAMPLE with splines (3 surfaces, can be less or more)
!
!   magfie_spline = 1
!   ALLOCATE(magfie_sarray(3))
!   magfie_sarray = (/0.4_dp,0.5_dp,0.6_dp/)
!   CALL stevvo( RT0, R0i, L1i, cbfi, bz0i, bf0 )
!   CALL magfie( x, bmoda, sqrtg, bder, hcovar, hctrvr, hcurl )
!
! EXAMPLE with splines (1 surface)
!
!   magfie_spline = 1
!   ALLOCATE(magfie_sarray(1))
!   magfie_sarray = (/0.4_dp/)
!   CALL stevvo( RT0, R0i, L1i, cbfi, bz0i, bf0 )
!   CALL magfie( x, bmoda, sqrtg, bder, hcovar, hctrvr, hcurl )
!
! EXAMPLE without splines
!
!   magfie_spline = 0
!   CALL stevvo( RT0, R0i, L1i, cbfi, bz0i, bf0 )
!   CALL magfie( x, bmoda, sqrtg, bder, hcovar, hctrvr, hcurl )
!
! EXAMPLE with loop over different s-values, spline, 3 surfaces
!
!   CALL stevvo( RT0, R0i, L1i, cbfi, bz0i, bf0 )
!   magfie_spline = 1
!   ALLOCATE( magfie_sarray(3) )
!   DO ...... ! different Fluxsurfaces
!      magfie_sarray = ......
!      DO ..... ! theta, phi
!         CALL magfie( x, bmoda, sqrtg, bder, hcovar, hctrvr, hcurl )
!      END DO
!      CALL magfie_deallocate
!   END DO
!   DEALLOCATE( magfie_sarray )

MODULE magfie_mod

  USE neo_precision,  ONLY: dp 
  IMPLICIT NONE

  !! Modifications by Andreas F. Martitsch (09.03.2014)
  ! Optional output (necessary for modeling the magnetic rotation)
  INTERFACE magfie
     MODULE PROCEDURE magfie_x, magfie_x1, magfie_x2
  END INTERFACE magfie
  !! End Modifications by Andreas F. Martitsch (09.03.2014)

  INTERFACE magfie_deallocate
     MODULE PROCEDURE magfie_deallocate_x
  END INTERFACE magfie_deallocate

  INTERFACE stevvo
     MODULE PROCEDURE stevvo_x
  END INTERFACE stevvo

CONTAINS

  SUBROUTINE magfie_x( x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
    USE neo_magfie_mod, ONLY: neo_magfie
    REAL(dp), DIMENSION(:),       INTENT(in)         :: x
    REAL(dp),                     INTENT(out)        :: bmod
    REAL(dp),                     INTENT(out)        :: sqrtg
    REAL(dp), DIMENSION(SIZE(x)), INTENT(out)        :: bder
    REAL(dp), DIMENSION(SIZE(x)), INTENT(out)        :: hcovar
    REAL(dp), DIMENSION(SIZE(x)), INTENT(out)        :: hctrvr
    REAL(dp), DIMENSION(SIZE(x)), INTENT(out)        :: hcurl
    CALL neo_magfie( x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
  END SUBROUTINE magfie_x

  !! Modifications by Andreas F. Martitsch (09.03.2014)
  ! Optional output (necessary for modeling the magnetic rotation)
  SUBROUTINE magfie_x1( x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl, bcovar_s_hat_der)
    USE neo_magfie_mod, ONLY: neo_magfie
    REAL(dp), DIMENSION(:),       INTENT(in)         :: x
    REAL(dp),                     INTENT(out)        :: bmod
    REAL(dp),                     INTENT(out)        :: sqrtg
    REAL(dp), DIMENSION(SIZE(x)), INTENT(out)        :: bder
    REAL(dp), DIMENSION(SIZE(x)), INTENT(out)        :: hcovar
    REAL(dp), DIMENSION(SIZE(x)), INTENT(out)        :: hctrvr
    REAL(dp), DIMENSION(SIZE(x)), INTENT(out)        :: hcurl
    REAL(dp), DIMENSION(SIZE(x)), INTENT(out)        :: bcovar_s_hat_der
    CALL neo_magfie( x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl, bcovar_s_hat_der)
  END SUBROUTINE magfie_x1
  !! End Modifications by Andreas F. Martitsch (09.03.2014)

  !! Modifications by Andreas F. Martitsch (13.11.2014)
  ! Optional output for NTV output
  SUBROUTINE magfie_x2( x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl, bcovar_s_hat_der, R, Z)
    USE neo_magfie_mod, ONLY: neo_magfie
    REAL(dp), DIMENSION(:),       INTENT(in)         :: x
    REAL(dp),                     INTENT(out)        :: bmod
    REAL(dp),                     INTENT(out)        :: sqrtg
    REAL(dp), DIMENSION(SIZE(x)), INTENT(out)        :: bder
    REAL(dp), DIMENSION(SIZE(x)), INTENT(out)        :: hcovar
    REAL(dp), DIMENSION(SIZE(x)), INTENT(out)        :: hctrvr
    REAL(dp), DIMENSION(SIZE(x)), INTENT(out)        :: hcurl
    REAL(dp), DIMENSION(SIZE(x)), INTENT(out)        :: bcovar_s_hat_der
    REAL(dp),                     INTENT(out)        :: R
    REAL(dp),                     INTENT(out)        :: Z
    CALL neo_magfie( x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl, bcovar_s_hat_der, R, Z)
  END SUBROUTINE magfie_x2
  !! End Modifications by Andreas F. Martitsch (13.11.2014)
  
  !! Modifications by Andreas F. Martitsch (11.03.2014)
  ! Deallocate the storage arrays for the 2d spline interpolation
  ! (over the flux-surface) of the additionally needed metric
  ! tensor elements
  SUBROUTINE magfie_deallocate_x
    USE neo_magfie_mod, ONLY: curr_tor_array, curr_tor_s_array, &
         curr_pol_array, curr_pol_s_array, iota_array,          &
         bmod_spl, bb_s_spl, bb_tb_spl, bb_pb_spl,              &
         gstb_spl, gspb_spl, gstb_tb_spl, gspb_tb_spl,          &
         gstb_pb_spl, gspb_pb_spl, R_spl,                       &
         magfie_newspline
    
    IF ( ALLOCATED( curr_tor_array ) )   DEALLOCATE( curr_tor_array )
    IF ( ALLOCATED( curr_tor_s_array ) ) DEALLOCATE( curr_tor_s_array )
    IF ( ALLOCATED( curr_pol_array ) )   DEALLOCATE( curr_pol_array )
    IF ( ALLOCATED( curr_pol_s_array ) ) DEALLOCATE( curr_pol_s_array )
    IF ( ALLOCATED( iota_array ) )       DEALLOCATE( iota_array )
    IF ( ALLOCATED( bmod_spl ) )         DEALLOCATE( bmod_spl )
    IF ( ALLOCATED( bb_s_spl ) )         DEALLOCATE( bb_s_spl )
    IF ( ALLOCATED( bb_tb_spl ) )        DEALLOCATE( bb_tb_spl )
    IF ( ALLOCATED( bb_pb_spl ) )        DEALLOCATE( bb_pb_spl )
    !! Modifications by Andreas F. Martitsch (11.03.2014)
    IF ( ALLOCATED( gstb_spl ) )         DEALLOCATE( gstb_spl )
    IF ( ALLOCATED( gspb_spl ) )         DEALLOCATE( gspb_spl )
    IF ( ALLOCATED( gstb_tb_spl ) )      DEALLOCATE( gstb_tb_spl )
    IF ( ALLOCATED( gspb_tb_spl ) )      DEALLOCATE( gspb_tb_spl )
    IF ( ALLOCATED( gstb_pb_spl ) )      DEALLOCATE( gstb_pb_spl )
    IF ( ALLOCATED( gspb_pb_spl ) )      DEALLOCATE( gspb_pb_spl )
    IF ( ALLOCATED( R_spl ) )            DEALLOCATE( R_spl )
    !! End Modifications by Andreas F. Martitsch (11.03.2014)
    magfie_newspline = 1
    
  END SUBROUTINE magfie_deallocate_x
  !! End Modifications by Andreas F. Martitsch (11.03.2014)

  SUBROUTINE stevvo_x( bigR, R0i, L1i, cbfi, bz0i, bf0 )
    USE neo_input,  ONLY: nfp 
    USE neo_exchange,  ONLY: rt0, bmref_g
    USE neo_magfie_mod, ONLY: magfie_spline, magfie_sarray
    REAL(dp),                     INTENT(out)        :: bigR
    REAL(dp),                     INTENT(out)        :: R0i
    INTEGER,                      INTENT(out)        :: L1i
    REAL(dp),                     INTENT(out)        :: cbfi
    REAL(dp),                     INTENT(out)        :: bz0i
    REAL(dp),                     INTENT(out)        :: bf0

    REAL(dp), DIMENSION(3) :: x
    REAL(dp)               :: bmod
    REAL(dp)               :: sqrtg
    REAL(dp), DIMENSION(3) :: bder
    REAL(dp), DIMENSION(3) :: hcovar
    REAL(dp), DIMENSION(3) :: hctrvr
    REAL(dp), DIMENSION(3) :: hcurl

    ! just to initialize
    x    = 0.0_dp
    IF (magfie_spline .EQ. 1) THEN
       x(1) = magfie_sarray(1)
    ELSE
       x(1) = 0.5
    END IF

    CALL magfie( x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl )

!!!    bigR = rt0      ! big radius
    bigR = rt0 * 1d2 ! big radius (cm)                    !!!
    R0i  = 0.3_dp    ! small radius - not really used
    L1i  = nfp       ! number of field periods
    cbfi = 0.0_dp    ! unused
    bz0i = 0.0_dp    ! unused
    bf0  = bmref_g   ! reference magnetic field

  END SUBROUTINE stevvo_x
    
END MODULE magfie_mod







