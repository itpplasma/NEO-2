MODULE neo_magfie_mod

  USE neo_precision
  USE neo_input,                                                       &
       ONLY: es, ixm, ixn, mnmax, psi_pr, pixm, pixn, nfp
  USE neo_control,                                                     &
       ONLY: fluxs_interp, write_progress, phi_n, theta_n, lab_swi
  USE neo_sub_mod,                                                     &
       ONLY: neo_read_control, neo_init, neo_init_spline
  USE neo_spline_data,                                                 &
       ONLY: r_mhalf,                                                  &
       a_bmnc, b_bmnc, c_bmnc, d_bmnc,                                 &
       a_rmnc, b_rmnc, c_rmnc, d_rmnc,                                 &
       a_zmnc, b_zmnc, c_zmnc, d_zmnc,                                 &
       a_lmnc, b_lmnc, c_lmnc, d_lmnc,                                 &
       a_iota, b_iota, c_iota, d_iota,                                 &
       a_curr_tor, b_curr_tor, c_curr_tor, d_curr_tor,                 &
       a_curr_pol, b_curr_pol, c_curr_pol, d_curr_pol,                 &
       a_pprime,   b_pprime,   c_pprime,   d_pprime,                   &
       a_sqrtg00,  b_sqrtg00,  c_sqrtg00,  d_sqrtg00
  USE inter_interfaces,                                                &
       ONLY: splint_horner3,                                           &
       tf, tfp, tfpp, tfppp,                                           &
       tfone, tfzero
  USE neo_work,                                                        &
       ONLY: cosmth, cosnph, sinmth, sinnph, theta_int, phi_int,       &
       theta_start, theta_end, phi_start, phi_end
  USE neo_actual_fluxs, ONLY : s_sqrtg00
  USE spline_mod, ONLY: spl2d, poi2d, eva2d

  !---------------------------------------------------------------------------
  !USE var_sub_misc, ONLY: fac_c,iota_m ! fac_m
  !---------------------------------------------------------------------------

  IMPLICIT NONE
  REAL(dp), DIMENSION(:), ALLOCATABLE              :: magfie_sarray
  INTEGER                                          :: magfie_spline    = 0
  INTEGER                                          :: magfie_newspline = 1
  ! switch for different output:
  !  0:  output for SMT and BMC
  !  1:  output for NEO 
  INTEGER                                          :: magfie_result    = 0
  INTEGER                                          :: magfie_sarray_len

  REAL(dp), DIMENSION(:), ALLOCATABLE              :: curr_tor_array
  REAL(dp), DIMENSION(:), ALLOCATABLE              :: curr_tor_s_array
  REAL(dp), DIMENSION(:), ALLOCATABLE              :: curr_pol_array
  REAL(dp), DIMENSION(:), ALLOCATABLE              :: curr_pol_s_array
  REAL(dp), DIMENSION(:), ALLOCATABLE              :: iota_array
  REAL(dp), DIMENSION(:), ALLOCATABLE              :: pprime_array
  REAL(dp), DIMENSION(:), ALLOCATABLE              :: sqrtg00_array

  REAL(dp)                                         :: s_pprime

  REAL(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE, TARGET      :: bmod_spl
  REAL(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE, TARGET      :: bb_s_spl
  REAL(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE, TARGET      :: bb_tb_spl
  REAL(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE, TARGET      :: bb_pb_spl
  REAL(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE, TARGET      :: gval_spl

  REAL(dp) :: boozer_iota
  REAL(dp) :: boozer_sqrtg00
  REAL(dp) :: boozer_curr_pol
  REAL(dp) :: boozer_curr_tor
  REAL(dp) :: boozer_psi_pr
  REAL(dp) :: boozer_sqrtg11 ! Test
  REAL(dp) :: boozer_isqrg
 
  REAL(dp), PRIVATE :: av_b2_m ! Klaus

  INTERFACE neo_magfie
     MODULE PROCEDURE neo_magfie_a
  END INTERFACE

CONTAINS

  SUBROUTINE neo_magfie_a( x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl )
    ! input / output
    REAL(dp), DIMENSION(:),       INTENT(in)         :: x
    REAL(dp),                     INTENT(out)        :: bmod
    REAL(dp),                     INTENT(out)        :: sqrtg
    REAL(dp), DIMENSION(SIZE(x)), INTENT(out)        :: bder
    REAL(dp), DIMENSION(SIZE(x)), INTENT(out)        :: hcovar
    REAL(dp), DIMENSION(SIZE(x)), INTENT(out)        :: hctrvr
    REAL(dp), DIMENSION(SIZE(x)), INTENT(out)        :: hcurl
    ! local definitions
    INTEGER(i4b)                                     :: swd = 1
    INTEGER                                          :: i, m, n
    INTEGER                                          :: npsi
    REAL(dp)                                         :: m0  = 0.0_dp
    REAL(dp)                                         :: yp, ypp, yppp

    REAL(dp)                                         :: bmnc, bmnc_s
    REAL(dp)                                         :: sinv, cosv
    REAL(dp)                                         :: iota
    REAL(dp)                                         :: curr_tor, curr_tor_s
    REAL(dp)                                         :: curr_pol, curr_pol_s
    REAL(dp)                                         :: bb_s, bb_tb, bb_pb
    REAL(dp)                                         :: fac, fac1

    INTEGER                                          :: k_es
    INTEGER                                          :: s_detected
    INTEGER                                          :: imn    
    INTEGER                                          :: it, ip, im, in
    INTEGER                                          :: mt = 1
    INTEGER                                          :: mp = 1
    INTEGER                                          :: theta_ind, phi_ind
    INTEGER                                          :: ierr
    REAL(dp)                                         :: s
    REAL(dp)                                         :: magfie_epsi = 1.e-9
    REAL(dp)                                         :: bi, bi_s, ri, zi, li
    REAL(dp)                                         :: theta_d, phi_d

    REAL(dp), DIMENSION(:), ALLOCATABLE              :: s_bmnc, s_bmnc_s
    REAL(dp), DIMENSION(:), ALLOCATABLE              :: s_rmnc, s_zmnc, s_lmnc
    REAL(dp), DIMENSION(:,:), ALLOCATABLE            :: bmod_a
    REAL(dp), DIMENSION(:,:), ALLOCATABLE            :: bb_s_a
    REAL(dp), DIMENSION(:,:), ALLOCATABLE            :: bb_tb_a
    REAL(dp), DIMENSION(:,:), ALLOCATABLE            :: bb_pb_a
    REAL(dp), DIMENSION(:,:), ALLOCATABLE            :: r,z,l
    REAL(dp), DIMENSION(:,:), ALLOCATABLE            :: r_tb,z_tb,p_tb
    REAL(dp), DIMENSION(:,:), ALLOCATABLE            :: r_pb,z_pb,p_pb
    REAL(dp), DIMENSION(:,:), ALLOCATABLE            :: gtbtb,gpbpb,gtbpb
    REAL(dp), DIMENSION(:,:), ALLOCATABLE            :: sqrg11_met

    REAL(dp), DIMENSION(:,:,:,:), POINTER            :: p_spl
    
    REAL(dp) :: isqrg, sqrg11

    !REAL(dp) :: s_sqrtg00_m Klaus
    
    !*******************************************************************
    ! Initialisation if necessary
    !*******************************************************************
    IF ( .NOT. ALLOCATED(es) ) THEN
       CALL neo_read_control()
       fluxs_interp = 1
       CALL neo_init(npsi)
       PRINT *, 'theta_start,theta_end,phi_start,phi_end'
       PRINT *, theta_start,theta_end,phi_start,phi_end
    END IF
    !*******************************************************************
    ! Spline of surfaces in magfie_sarray
    !*******************************************************************
    IF (magfie_spline .EQ. 1 .AND. magfie_newspline .EQ. 1) THEN
       magfie_sarray_len =  SIZE(magfie_sarray)
       !****************************************************************
       ! Allocation
       !****************************************************************
       ALLOCATE( bmod_spl(4,4,theta_n,phi_n,magfie_sarray_len) )
       ALLOCATE( bb_s_spl(4,4,theta_n,phi_n,magfie_sarray_len) )
       ALLOCATE( bb_tb_spl(4,4,theta_n,phi_n,magfie_sarray_len) )
       ALLOCATE( bb_pb_spl(4,4,theta_n,phi_n,magfie_sarray_len) )
       ALLOCATE( gval_spl(4,4,theta_n,phi_n,magfie_sarray_len) )

       ALLOCATE( curr_tor_array(magfie_sarray_len) )
       ALLOCATE( curr_tor_s_array(magfie_sarray_len) )
       ALLOCATE( curr_pol_array(magfie_sarray_len) )
       ALLOCATE( curr_pol_s_array(magfie_sarray_len) )
       ALLOCATE( iota_array(magfie_sarray_len) )
       ALLOCATE( pprime_array(magfie_sarray_len) )
       ALLOCATE( sqrtg00_array(magfie_sarray_len) )
       !****************************************************************
       ! Loop over predefined s-values 
       !****************************************************************
       DO k_es = 1, magfie_sarray_len
          s = magfie_sarray(k_es)
          !*************************************************************
          ! Surface
          !*************************************************************
          IF (write_progress .EQ. 1) THEN
             PRINT *, 'Initialize Surface, k_es = ',k_es
          END IF
          ALLOCATE ( s_bmnc(mnmax) )
          ALLOCATE ( s_bmnc_s(mnmax) )
          ALLOCATE ( s_rmnc(mnmax) )
          ALLOCATE ( s_zmnc(mnmax) )
          ALLOCATE ( s_lmnc(mnmax) )
          DO imn = 1, mnmax
             swd = 1
             CALL splint_horner3(es,                                   &
                  a_bmnc(:,imn), b_bmnc(:,imn),                        &
                  c_bmnc(:,imn), d_bmnc(:,imn),                        &
                  swd, r_mhalf(imn),                                   &
                  s, tf, tfp, tfpp, tfppp,                             &
                  s_bmnc(imn), s_bmnc_s(imn), ypp, yppp)
             swd = 0
             CALL splint_horner3(es,                                   &
                  a_rmnc(:,imn), b_rmnc(:,imn),                        &
                  c_rmnc(:,imn), d_rmnc(:,imn),                        &
                  swd, r_mhalf(imn),                                   &
                  s, tf, tfp, tfpp, tfppp,                             &
                  s_rmnc(imn), yp, ypp, yppp)
             CALL splint_horner3(es,                                   &
                  a_zmnc(:,imn), b_zmnc(:,imn),                        &
                  c_zmnc(:,imn), d_zmnc(:,imn),                        &
                  swd, r_mhalf(imn),                                   &
                  s, tf, tfp, tfpp, tfppp,                             &
                  s_zmnc(imn), yp, ypp, yppp)
             CALL splint_horner3(es,                                   &
                  a_lmnc(:,imn), b_lmnc(:,imn),                        &
                  c_lmnc(:,imn), d_lmnc(:,imn),                        &
                  swd, r_mhalf(imn),                                   &
                  s, tf, tfp, tfpp, tfppp,                             &
                  s_lmnc(imn), yp, ypp, yppp)
          END DO
          !*************************************************************
          ! Fourier summation for the full theta-phi array
          !*************************************************************
          IF (write_progress .EQ. 1) THEN
             PRINT *, 'Do Fourier'
          END IF
          ALLOCATE( bmod_a(theta_n,phi_n) ) 
          ALLOCATE( bb_s_a(theta_n,phi_n) ) 
          ALLOCATE( bb_tb_a(theta_n,phi_n) ) 
          ALLOCATE( bb_pb_a(theta_n,phi_n) )
          bmod_a  = 0.0_dp
          bb_s_a  = 0.0_dp
          bb_tb_a = 0.0_dp
          bb_pb_a = 0.0_dp

          ALLOCATE( r(theta_n,phi_n) )  ! NEW
          ALLOCATE( z(theta_n,phi_n) ) 
          ALLOCATE( l(theta_n,phi_n) ) 
          ALLOCATE( r_tb(theta_n,phi_n) ) 
          ALLOCATE( z_tb(theta_n,phi_n) ) 
          ALLOCATE( p_tb(theta_n,phi_n) ) 
          ALLOCATE( r_pb(theta_n,phi_n) ) 
          ALLOCATE( z_pb(theta_n,phi_n) ) 
          ALLOCATE( p_pb(theta_n,phi_n) ) 
          r = 0.0d0
          z = 0.0d0
          l = 0.0d0
          r_tb = 0.0d0
          z_tb = 0.0d0
          p_tb = 0.0d0
          r_pb = 0.0d0
          z_pb = 0.0d0
          p_pb = 0.0d0

          DO imn=1,mnmax
             ri = s_rmnc(imn) ! NEW
             zi = s_zmnc(imn)
             li = s_lmnc(imn)

             bi   = s_bmnc(imn)
             bi_s = s_bmnc_s(imn)
             m = ixm(imn)
             n = ixn(imn)
             im = pixm(imn)
             in = pixn(imn)
             DO ip=1,phi_n
                DO it=1,theta_n
                   cosv = cosmth(it,im) * cosnph(ip,in) + sinmth(it,im) * sinnph(ip,in)
                   sinv = sinmth(it,im) * cosnph(ip,in) - cosmth(it,im) * sinnph(ip,in)
                   
                   bmod_a(it,ip)   = bmod_a(it,ip)   +     bi   * cosv
                   bb_s_a(it,ip)   = bb_s_a(it,ip)   +     bi_s * cosv
                   bb_tb_a(it,ip)  = bb_tb_a(it,ip)  - m * bi   * sinv
                   bb_pb_a(it,ip)  = bb_pb_a(it,ip)  + n * bi   * sinv

                   r(it,ip) = r(it,ip) + ri*cosv
                   z(it,ip) = z(it,ip) + zi*sinv
                   l(it,ip) = l(it,ip) + li*sinv

                   r_tb(it,ip) = r_tb(it,ip) - m*ri*sinv
                   r_pb(it,ip) = r_pb(it,ip) + n*ri*sinv
                   z_tb(it,ip) = z_tb(it,ip) + m*zi*cosv
                   z_pb(it,ip) = z_pb(it,ip) - n*zi*cosv
                   p_tb(it,ip) = p_tb(it,ip) - m*li*cosv ! -l_tb
                   p_pb(it,ip) = p_pb(it,ip) + n*li*cosv ! -l_pb

                END DO
             END DO
          END DO
          DEALLOCATE( s_bmnc )
          DEALLOCATE( s_bmnc_s )
          DEALLOCATE( s_rmnc )
          DEALLOCATE( s_zmnc )
          DEALLOCATE( s_lmnc )

          IF (lab_swi .EQ. 5 .OR. lab_swi .EQ. 3) THEN ! CHS, LHD
             p_tb = - p_tb
             p_pb = 1 - p_pb
          ELSE
             p_tb = p_tb * twopi / nfp
             p_pb = 1.0_dp + p_pb * twopi / nfp
          END IF
          
          ! **********************************************************************
          ! Ensure periodicity boundaries to be the same
          ! **********************************************************************
          r(theta_n,:) = r(1,:)
          r(:,phi_n)   = r(:,1)
          z(theta_n,:) = z(1,:)
          z(:,phi_n)   = z(:,1)
          l(theta_n,:) = l(1,:)
          l(:,phi_n)   = l(:,1)
          bmod_a(theta_n,:) = bmod_a(1,:)
          bmod_a(:,phi_n)   = bmod_a(:,1)
          r_tb(theta_n,:) = r_tb(1,:)
          r_tb(:,phi_n)   = r_tb(:,1)
          r_pb(theta_n,:) = r_pb(1,:)
          r_pb(:,phi_n)   = r_pb(:,1)
          z_tb(theta_n,:) = z_tb(1,:)
          z_tb(:,phi_n)   = z_tb(:,1)
          z_pb(theta_n,:) = z_pb(1,:)
          z_pb(:,phi_n)   = z_pb(:,1)
          p_tb(theta_n,:) = p_tb(1,:)
          p_tb(:,phi_n)   = p_tb(:,1)
          p_pb(theta_n,:) = p_pb(1,:)
          p_pb(:,phi_n)   = p_pb(:,1)
          bb_tb_a(theta_n,:) = bb_tb_a(1,:)
          bb_tb_a(:,phi_n)   = bb_tb_a(:,1)
          bb_s_a(theta_n,:)  = bb_tb_a(1,:)
          bb_s_a(:,phi_n)    = bb_tb_a(:,1)
          bb_pb_a(theta_n,:) = bb_pb_a(1,:)
          bb_pb_a(:,phi_n)   = bb_pb_a(:,1)
          
          ! **********************************************************************
          ! Derived quantities
          ! **********************************************************************
          ALLOCATE( gtbtb(theta_n,phi_n) ) 
          ALLOCATE( gpbpb(theta_n,phi_n) ) 
          ALLOCATE( gtbpb(theta_n,phi_n) ) 
          ALLOCATE( sqrg11_met(theta_n,phi_n) ) 
          ! metric tensor
          gtbtb = r_tb*r_tb + z_tb*z_tb + r*r*p_tb*p_tb  
          gpbpb = r_pb*r_pb + z_pb*z_pb + r*r*p_pb*p_pb  
          gtbpb = r_tb*r_pb + z_tb*z_pb + r*r*p_tb*p_pb

          ! Winny for Klaus
          av_b2_m = theta_n * phi_n / SUM(1 / (bmod_a*bmod_a))
          !PRINT *, 'theta_n,phi_n ',theta_n,phi_n 
          !PRINT *, 'av_b2_m ',av_b2_m
          ! Winny for Klaus - Ende

  
          ! $1/sqrt(g)$
          !! fac = s_curr_pol + s_iota * s_curr_tor  ! (J + iota I)
          !! isqrg  = b*b / fac 
          ! $sqrt(g^{11})$
          !IF (g11_swi .EQ. 0) THEN
          !   sqrg11 = 1.0_dp
          !ELSE
             !sqrg11 = SQRT( (gtbtb*gpbpb - gtbpb*gtbpb ) * isqrg**2)
             sqrg11_met = SQRT( (gtbtb*gpbpb - gtbpb*gtbpb ) )
          !END IF

             !PRINT *, 'max_gtbtb = ',maxval(gtbtb)
             !PRINT *, 'min_gtbtb = ',MINVAL(gtbtb)
             !PRINT *, 'max_gpbpb = ',maxval(gpbpb)
             !PRINT *, 'min_gpbpb = ',MINVAL(gpbpb)
             !PRINT *, 'max_gtbpb = ',MAXVAL(gtbpb)
             !PRINT *, 'min_gtbpb = ',MINVAL(gtbpb)
            
          !*************************************************************
          ! Do the 2-D periodic spline
          !*************************************************************
          IF (write_progress .EQ. 1) THEN
             PRINT *, 'Do 2-D spline of bmod'
          END IF
          p_spl => bmod_spl(:,:,:,:,k_es)
          CALL spl2d(theta_n,phi_n,theta_int,phi_int,mt,mp,            &
               bmod_a,p_spl)
          IF (write_progress .EQ. 1) THEN
             PRINT *, 'Do 2-D spline of bb_s'
          END IF
          p_spl => bb_s_spl(:,:,:,:,k_es)
          CALL spl2d(theta_n,phi_n,theta_int,phi_int,mt,mp,            &
               bb_s_a,p_spl)
          IF (write_progress .EQ. 1) THEN
             PRINT *, 'Do 2-D spline of bb_tb'
          END IF
          p_spl => bb_tb_spl(:,:,:,:,k_es)
          CALL spl2d(theta_n,phi_n,theta_int,phi_int,mt,mp,            &
               bb_tb_a,p_spl)
          IF (write_progress .EQ. 1) THEN
             PRINT *, 'Do 2-D spline of bb_pb'
          END IF
          p_spl => bb_pb_spl(:,:,:,:,k_es)
          CALL spl2d(theta_n,phi_n,theta_int,phi_int,mt,mp,            &
               bb_pb_a,p_spl)
          IF (write_progress .EQ. 1) THEN
             PRINT *, 'Do 2-D spline of sqrg11'
          END IF
          p_spl => gval_spl(:,:,:,:,k_es)
          CALL spl2d(theta_n,phi_n,theta_int,phi_int,mt,mp,            &
               sqrg11_met,p_spl)

          !PRINT *, 'max_g11 = ', maxval(sqrg11_met)
          !PRINT *, 'min_g11 = ', minval(sqrg11_met)

          DEALLOCATE( bmod_a )
          DEALLOCATE( bb_s_a )
          DEALLOCATE( bb_tb_a )
          DEALLOCATE( bb_pb_a )

          DEALLOCATE( r )  ! NEW
          DEALLOCATE( z ) 
          DEALLOCATE( l ) 
          DEALLOCATE( r_tb ) 
          DEALLOCATE( z_tb ) 
          DEALLOCATE( p_tb ) 
          DEALLOCATE( r_pb ) 
          DEALLOCATE( z_pb ) 
          DEALLOCATE( p_pb ) 

          DEALLOCATE( gtbtb ) 
          DEALLOCATE( gpbpb ) 
          DEALLOCATE( gtbpb ) 
          DEALLOCATE( sqrg11_met ) 

          !*************************************************************
          ! Provide curr_tor, curr_tor_s, curr_pol, curr_pol_s, iota
          !*************************************************************
          IF (write_progress .EQ. 1) THEN
             PRINT *, 'Prep of currents: ',s
          END IF
          swd = 1 ! derivative
          CALL splint_horner3(es,                                      &
               a_curr_tor, b_curr_tor, c_curr_tor, d_curr_tor,         &
               swd, m0,                                                &
               s, tfone, tfzero, tfzero, tfzero,                       &
               curr_tor_array(k_es), curr_tor_s_array(k_es), ypp, yppp)
          swd = 1 ! derivative
          CALL splint_horner3(es,                                      &
               a_curr_pol, b_curr_pol, c_curr_pol, d_curr_pol,         &
               swd, m0,                                                &
               s, tfone, tfzero, tfzero, tfzero,                       &
               curr_pol_array(k_es), curr_pol_s_array(k_es) ,ypp, yppp)    
          swd = 0 ! no derivative
          CALL splint_horner3(es,                                      &
               a_iota, b_iota, c_iota, d_iota, swd, m0,                &
               s, tfone, tfzero, tfzero, tfzero,                       &
               iota_array(k_es), yp, ypp, yppp)
          CALL splint_horner3(es,                                      &
               a_pprime, b_pprime, c_pprime, d_pprime, swd, m0,                &
               s, tfone, tfzero, tfzero, tfzero,                       &
               pprime_array(k_es), yp, ypp, yppp)
          CALL splint_horner3(es,                                      &
               a_sqrtg00, b_sqrtg00, c_sqrtg00, d_sqrtg00, swd, m0,                &
               s, tfone, tfzero, tfzero, tfzero,                       &
               sqrtg00_array(k_es), yp, ypp, yppp)
       END DO
       magfie_newspline = 0
    END IF

    s_detected = 0
    IF (magfie_spline .EQ. 1) THEN
       s = x(1)
       !****************************************************************
       ! Detection of index
       !****************************************************************
       DO k_es = 1, magfie_sarray_len
          IF ( ABS(s-magfie_sarray(k_es)) .LT. magfie_epsi) THEN
             s_detected = 1
             EXIT
          END IF
       END DO
       IF (s_detected .EQ. 1) THEN
          curr_tor   = curr_tor_array(k_es)
          curr_tor_s = curr_tor_s_array(k_es)
          curr_pol   = curr_pol_array(k_es)
          curr_pol_s = curr_pol_s_array(k_es)
          iota       = iota_array(k_es)
          s_pprime   = pprime_array(k_es) ! only local
          s_sqrtg00  = sqrtg00_array(k_es)
          ! ************************************************************
          ! Evaluation of 2d-splines
          ! ************************************************************
          CALL poi2d(theta_int,phi_int,mt,mp,                          &
               theta_start,theta_end,phi_start,phi_end,                &
               x(3),x(2),theta_ind,phi_ind,theta_d,phi_d,ierr)
          p_spl => bmod_spl(:,:,:,:,k_es)
          CALL eva2d(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,    &
               p_spl,bmod)
          p_spl => bb_s_spl(:,:,:,:,k_es)
          CALL eva2d(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,    &
               p_spl,bb_s)
          p_spl => bb_tb_spl(:,:,:,:,k_es)
          CALL eva2d(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,    &
               p_spl,bb_tb)
          p_spl => bb_pb_spl(:,:,:,:,k_es)
          CALL eva2d(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,    &
               p_spl,bb_pb)
          p_spl => gval_spl(:,:,:,:,k_es)
          CALL eva2d(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,    &
               p_spl,sqrg11)
          ! $1/sqrt(g)$
          fac = curr_pol + iota * curr_tor  ! (J + iota I)
          isqrg  = bmod*bmod / fac

          ! Winny for Klaus
          !s_sqrtg00_m = fac / av_b2_m
          s_sqrtg00 = fac / av_b2_m
          !PRINT *, 's_sqrtg00, s_sqrtg00_m ',s_sqrtg00, s_sqrtg00_m
          !PRINT *, 'fac, av_b2_m ',fac, av_b2_m
          !PAUSE
          ! Winny for Klaus - Ende

          !PRINT *, ' '
          !PRINT *, 'curr_pol = ',curr_pol
          !PRINT *, 'curr_tor = ',curr_tor
          !PRINT *, 'iota     = ',iota
          !PRINT *, 'fac      = ',fac
          !PRINT *, 'bmod     = ',bmod
          !PRINT *, 'isqrg    = ',isqrg
          !PRINT *, 'sqrg     = ',1.d0 / isqrg

          !PRINT *, 'sqrg11_n = ',sqrg11
          sqrg11 = sqrg11 * ABS(isqrg)
          !PRINT *, 'sqrg11   = ',sqrg11
       ELSE
          PRINT *, 'neo_magfie: s not detected!'
          STOP
       END IF
    END IF

    IF (magfie_spline .EQ. 0 .OR. s_detected .EQ. 0) THEN
       IF (magfie_spline .EQ. 1 .AND. s_detected .EQ. 0) THEN
          PRINT *, 'WARNING from neo_magfie - s out of range: ',s
          PRINT *, ' Using Fourier Summation directly'
       END IF
       
       PRINT *, 'magfie_spline .EQ. 0 not implemented'
       STOP

       !****************************************************************
       ! Direct summation of Fourier components
       !****************************************************************
       bmod   = 0.0_dp
       bb_s   = 0.0_dp
       bb_tb  = 0.0_dp
       bb_pb  = 0.0_dp

       DO i = 1, mnmax
          swd = 1
          CALL splint_horner3(es,                                    &
               a_bmnc(:,i), b_bmnc(:,i), c_bmnc(:,i), d_bmnc(:,i),   &
               swd, r_mhalf(i),                                      &
               x(1), tf, tfp, tfpp, tfppp,                           &
               bmnc, bmnc_s, ypp, yppp)

          m = ixm(i)
          n = ixn(i)
          sinv = SIN(m*x(3) - n*x(2))
          cosv = COS(m*x(3) - n*x(2))
          
          bmod   = bmod   +     bmnc   * cosv
          bb_s   = bb_s   +     bmnc_s * cosv
          bb_tb  = bb_tb  - m * bmnc   * sinv
          bb_pb  = bb_pb  + n * bmnc   * sinv
       END DO

       swd = 1
       CALL splint_horner3(es,                                       &
            a_curr_tor, b_curr_tor, c_curr_tor, d_curr_tor,          &
            swd, m0,                                                 &
            x(1), tfone, tfzero, tfzero, tfzero,                     &
            curr_tor, curr_tor_s, ypp, yppp)
       CALL splint_horner3(es,                                       &
            a_curr_pol, b_curr_pol, c_curr_pol, d_curr_pol,          &
            swd, m0,                                                 &
            x(1), tfone, tfzero, tfzero, tfzero,                     &
            curr_pol, curr_pol_s ,ypp, yppp)    
       swd = 0 ! no derivative
       CALL splint_horner3(es,                                       &
            a_iota, b_iota, c_iota, d_iota, swd, m0,                 &
            x(1), tfone, tfzero, tfzero, tfzero,                     &
            iota, yp, ypp, yppp)       
    END IF

    IF (magfie_result .EQ. 1) THEN
       ! This was the original version:     
       ! derived quantities in (s,theta_b,phi_b)-system
       fac   = (curr_pol + iota * curr_tor) * psi_pr
       fac1  = fac  / bmod                 ! sqrtg*bmod
       sqrtg = fac1 / bmod 

       bder(1) = bb_s
       bder(2) = bb_tb
       bder(3) = bb_pb

       hcovar(1) = 0.0_dp
       hcovar(2) = curr_tor / bmod
       hcovar(3) = curr_pol / bmod

       hctrvr(1) = 0.0_dp
       hctrvr(2) = iota / fac1
       hctrvr(3) = 1.0_dp / fac1

       hcurl(1)  = (curr_pol * bb_pb      - curr_tor * bb_tb     ) / fac 
       hcurl(2)  = (curr_pol * bb_s       - bmod     * curr_pol_s) / fac 
       hcurl(3)  = (bmod     * curr_tor_s - curr_tor * bb_s      ) / fac 
       ! Remark by Winny:
       ! The consisteny check for curr_pol shows a problem in all
       ! Greifswald (standard) input files
       ! According to the consistency check, 
       ! curr_pol has to be changed to -curr_pol

    ELSEIF ( magfie_result .EQ. 0 ) THEN
       ! Modifications made by Sergie for use in SMT
       ! derived quantities in (s,theta_b,phi_b)-system
       !fac   = (curr_pol + iota * curr_tor) * psi_pr
       ! This is used in NEO2
       fac   =  curr_pol + iota * curr_tor                       !!!
       fac1  = fac  / bmod                 ! sqrtg*bmod
       fac = fac * psi_pr                                        !!!
       !    sqrtg = fac1 / bmod 
       sqrtg = - fac1 / bmod * psi_pr * 1d6                      !!!
       !---------------------------------------------------------------------------
       !  iota_m = iota
       ! fac_m  =  (curr_pol + iota * curr_tor) * 1d6 * psi_pr
       !  fac_c  =  (curr_pol + iota * curr_tor) * 1d6 
       !---------------------------------------------------------------------------

       bder(1) = bb_s
       bder(3) = bb_tb
       bder(2) = bb_pb
       bder=bder / bmod                                          !!!

       hcovar(1) = 0.0_dp
       hcovar(3) = curr_tor / bmod
       hcovar(2) = curr_pol / bmod
       hcovar=hcovar * 1.d2                                      !!!

       hctrvr(1) = 0.0_dp
       hctrvr(3) = iota / fac1
       hctrvr(2) = 1.0_dp / fac1
       hctrvr=hctrvr * 1d-2                                      !!!

       !    hcurl(1)  = (curr_pol * bb_pb      - curr_tor * bb_tb     ) / fac 
       hcurl(1)  = (curr_tor * bb_pb      - curr_pol * bb_tb     ) / fac  !!!
       hcurl(3)  = (curr_pol * bb_s       - bmod     * curr_pol_s) / fac 
       hcurl(2)  = (bmod     * curr_tor_s - curr_tor * bb_s      ) / fac 
       hcurl=hcurl * 1d-4                                                 !!!
    END IF
    
    boozer_iota = iota
    boozer_sqrtg00 = s_sqrtg00
    boozer_curr_tor = curr_tor
    boozer_curr_pol = curr_pol
    boozer_psi_pr = psi_pr 
    boozer_sqrtg11 = sqrg11
    boozer_isqrg = isqrg

  END SUBROUTINE neo_magfie_a

END MODULE neo_magfie_mod

