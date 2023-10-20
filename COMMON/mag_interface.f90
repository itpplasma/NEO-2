MODULE gfactor_mod

  IMPLICIT NONE
  INTEGER :: ienter=1,npoia,npoib
  DOUBLE PRECISION :: ha,hb
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: garr

CONTAINS

  DOUBLE PRECISION FUNCTION gfactor(a,b)

    INTEGER :: nistep,i,k,j
    DOUBLE PRECISION :: a,b,atmp,btmp,wa,wb,x,hint
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: arr

    IF(ienter.EQ.1) THEN
       ienter=0
       npoia=100
       npoib=100
       nistep=100
       ha=1.d0/npoia
       hb=1.d0/(npoib+1)
       hint=1.d0/nistep
       ALLOCATE(garr(0:npoia,0:npoib),arr(0:nistep))
       DO j=0,npoia
          atmp=j*ha
          DO k=0,npoib
             btmp=k*hb
             DO i=0,nistep
                x=i*hint
                arr(i)=SQRT(ABS(1.d0-atmp*(3.d0*x**2-2.d0*x**3)))                 &
                     /(1.d0-btmp*(3.d0*(x-1.d0)**2+2.d0*(x-1.)**3))
             ENDDO
             garr(j,k)=0.5d0*hint*(SUM(arr(0:nistep-1))+SUM(arr(1:nistep)))
          ENDDO
       ENDDO
       DEALLOCATE(arr)
    ENDIF

    wa=a/ha
    j=INT(wa)
    j=MIN(npoia-1,MAX(0,j))
    wa=wa-j

    wb=b/hb
    k=INT(wb)
    k=MIN(npoib-1,MAX(0,k))
    wb=MIN(1.d0,wb-k)

    gfactor=(1.d0-wa)*(garr(j,k)*(1.d0-wb)+garr(j,k+1)*wb)                    &
         +       wa*(garr(j+1,k)*(1.d0-wb)+garr(j+1,k+1)*wb)

  END FUNCTION gfactor
END MODULE gfactor_mod

MODULE extremum_mod

  USE magnetics_mod
  USE plagrange_mod

  IMPLICIT NONE
  INTEGER, PARAMETER, PRIVATE :: dp = KIND(1.0d0)
  INTEGER, PARAMETER :: nlagrange = 5

  PUBLIC find_extremum
  PRIVATE find_ext
  INTERFACE find_extremum
     MODULE PROCEDURE find_ext
  END INTERFACE

CONTAINS
  !---------------------------------------------------------------------

  SUBROUTINE find_ext(fieldperiod,x1i,x2i,dxi,x,y)
    TYPE(fieldperiod_struct), POINTER :: fieldperiod
    REAL(kind=dp), INTENT(in)   :: x1i,x2i,dxi 
    REAL(kind=dp), INTENT(out)  :: x,y 

    REAL(kind=dp) :: A,d,x_hlp,d2,dummy,dx,x1_p,x2_p
    REAL(kind=dp) :: x1_in,x2_in,x1,x2,x3,x4
    REAL(kind=dp) :: fx1,fx2,fxm,fx3,fx4
    INTEGER       :: n,k
    LOGICAL       :: L

    x1 = x1i
    x2 = x2i

    A = 2.0d0 / (1.0d0 + SQRT(5.0d0))
    n = 0
    dx = ABS(dxi)

    IF (x2 .LT. x1) THEN
       x_hlp = x2
       x2 = x1
       x1 = x_hlp
    END IF

    x1_in = x1
    x2_in = x2
    
    ! second derivative
    CALL plagrange_interp(fieldperiod,x1,nlagrange,fx1,dummy)
    CALL plagrange_interp(fieldperiod,x2,nlagrange,fx2,dummy)
    CALL plagrange_interp(fieldperiod,(x1+x2)/2.0d0,nlagrange,fxm,dummy)
    d2 = (fx1+fx2-2.0d0*fxm) / ((x2-x1)/2.0d0)**2

    DO k = 1,2

       d = x2 - x1
       x3 = x1 + A * d
       x4 = x2 - A * d
    
       CALL plagrange_interp(fieldperiod,x3,nlagrange,fx3,dummy)
       CALL plagrange_interp(fieldperiod,x4,nlagrange,fx4,dummy)

       DO WHILE (ABS(d)/ABS(MAX(x1_in,x2_in)) .GT. dx)
          n = n + 1
          IF (d2 < 0) THEN
             L = fx4 < fx3
          ELSE
             L = fx3 < fx4
          END IF
          x1_p = x1
          x2_p = x2
          IF (L) THEN
             x1 = x4
             x4 = x3
             fx4 = fx3
             d = x2 - x1
             x3 = x1 + A * d
             CALL plagrange_interp(fieldperiod,x3,nlagrange,fx3,dummy)
          ELSE
             x2 = x3
             x3 = x4
             fx3 = fx4 
             d = x2 - x1
             x4 = x2 - A * d
             CALL plagrange_interp(fieldperiod,x4,nlagrange,fx4,dummy)
          END IF

          IF (x1_p .EQ. x1 .AND. x2_p .EQ. x2) THEN
             EXIT
          END IF
       END DO

       CALL plagrange_interp(fieldperiod,x1,nlagrange,fx1,dummy)
       CALL plagrange_interp(fieldperiod,x2,nlagrange,fx2,dummy)
    
       IF (d2 .LT. 0) THEN
          L = fx2 < fx1
       ELSE
          L = fx1 < fx2
       END IF
    
       IF (L) THEN
          x = x1 
          y = fx1
       ELSE
          x = x2
          y = fx2
       END IF

       IF (x .GT. x1_in .AND. x .LT. x2_in) THEN
          EXIT
       ELSEIF (x .LE. x1_in) THEN
          x2 = x1_in
          x1 = x1_in - (x2_in - x1_in)
       ELSE
          x1 = x2_in
          x2 = x2_in + (x2_in - x1_in)
       END IF
    END DO

  END SUBROUTINE find_ext

END MODULE extremum_mod

MODULE mag_interface_mod
  ! This module provides public routines for:
  !
  !  make_magnetics
  !   creates a device
  !   creates a surface
  !   creates a fieldline
  !    (two private routines at the moment make_mag_fieldline and make_mag_fieldline_new)
  !    must be resolved
  !
  !  ripple_eta_magnetics
  USE magnetics_mod
  USE device_mod

  IMPLICIT NONE

  ! double precision
  INTEGER, PARAMETER, PRIVATE :: dp = KIND(1.0d0)

  ! public 
  INTEGER, PUBLIC       :: mag_local_sigma = 0
  LOGICAL, PUBLIC       :: mag_symmetric = .FALSE.
  LOGICAL, PUBLIC       :: mag_symmetric_shorten = .FALSE.
  LOGICAL, PUBLIC       :: split_inflection_points = .TRUE.
  LOGICAL, PUBLIC       :: split_at_period_boundary = .FALSE.
  REAL(kind=dp), PUBLIC :: hphi_lim = 1.0d-10
  REAL(kind=dp), PUBLIC :: sigma_shield_factor = 5.0d0  
  INTEGER, PUBLIC       :: mag_magfield = 1
  INTEGER, PUBLIC       :: mag_coordinates = 0
  INTEGER, PUBLIC       :: mag_nperiod_min
  INTEGER, PUBLIC       :: mag_save_memory
  INTEGER, PUBLIC       :: mag_cycle_ripples
  INTEGER, PUBLIC       :: mag_start_special
  INTEGER, PUBLIC       :: mag_max_prop = 30
  INTEGER, PUBLIC       :: magnetic_device = 1
  INTEGER, PUBLIC       :: mag_close_fieldline = 1
  INTEGER, PUBLIC       :: mag_ripple_contribution = 1


  REAL(kind=dp), PUBLIC :: aiota_tokamak
  REAL(kind=dp), PUBLIC :: efit_raxis
  REAL(kind=dp), PUBLIC :: efit_zaxis
  REAL(kind=dp), PUBLIC :: boozer_s
  REAL(kind=dp), PUBLIC :: boozer_theta_beg
  REAL(kind=dp), PUBLIC :: boozer_phi_beg
  REAL(kind=dp), PUBLIC :: boozer_bmod0
  REAL(kind=dp), PUBLIC :: mag_dbhat_min
  REAL(kind=dp), PUBLIC :: mag_dphi_inf_min
  REAL(kind=dp), PUBLIC :: mag_inflection_mult
  REAL(kind=dp), PUBLIC :: average_bhat
  REAL(kind=dp), PUBLIC :: average_bhat2
  REAL(kind=dp), PUBLIC :: average_one_over_bhat
  REAL(kind=dp), PUBLIC :: surface_boozer_B00
  REAL(kind=dp), PUBLIC :: travis_convfac
  

  ! default values
  CHARACTER(len=100), PRIVATE :: name_def    = 'unnamed'
  INTEGER,            PRIVATE :: nperiod_def = 50
  INTEGER,            PRIVATE :: nstep_def   = 480
  INTEGER,            PRIVATE :: ndim_def    = 14

  ! private
  LOGICAL,                    PRIVATE :: first_ripple
  INTEGER,                    PRIVATE :: imin_ripple
  INTEGER,                    PRIVATE :: n_prop
  INTEGER,                    PRIVATE :: n_per
  INTEGER,                    PRIVATE :: n_ripple_tot
  REAL(kind=dp),              PRIVATE :: period_length
  REAL(kind=dp),              PRIVATE :: period_boundary

  REAL(kind=dp), ALLOCATABLE, PRIVATE :: phi_arr(:)
  INTEGER,       ALLOCATABLE, PRIVATE :: ripple_prop_bound(:)
  INTEGER,       ALLOCATABLE, PRIVATE :: ripple_period_bound(:)
  REAL(kind=dp), ALLOCATABLE, PRIVATE :: x1(:)
  REAL(kind=dp), ALLOCATABLE, PRIVATE :: x2(:)
  REAL(kind=dp), ALLOCATABLE, PRIVATE :: x3(:)
  REAL(kind=dp), ALLOCATABLE, PRIVATE :: bhat(:)
  REAL(kind=dp), ALLOCATABLE, PRIVATE :: geodcu(:)
  REAL(kind=dp), ALLOCATABLE, PRIVATE :: h_phi(:)
  REAL(kind=dp), ALLOCATABLE, PRIVATE :: dlogbdphi(:)
  !! Modifications by Andreas F. Martitsch (13.03.2014)
  ! Optional output (necessary for modeling the magnetic rotation)
  REAL(kind=dp), ALLOCATABLE, PRIVATE :: dbcovar_s_hat_dphi(:)
  REAL(kind=dp), ALLOCATABLE, PRIVATE :: bcovar_s_hat(:)
  REAL(kind=dp), ALLOCATABLE, PRIVATE :: dlogbds(:)
  !! End Modifications by Andreas F. Martitsch (13.03.2014)
  REAL(kind=dp), ALLOCATABLE, PRIVATE :: hlp_arr(:)

  ! internal constants
  REAL(kind=dp), PARAMETER, PRIVATE :: pi=3.14159265358979_dp

  ! ---------------------------------------------------------------------------
  ! make the real stuff
  PUBLIC make_magnetics
  PRIVATE                        &
       make_mag_device,          &
       make_mag_surface,         &
       make_mag_fieldline_newperiod
  INTERFACE make_magnetics
     MODULE PROCEDURE                   &
          make_mag_device,              &
          make_mag_surface,             &
          make_mag_fieldline_newperiod
  END INTERFACE
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  ! fixes eta_values for the ripples
  PUBLIC ripple_eta_magnetics
  PRIVATE ripple_eta_mag
  INTERFACE ripple_eta_magnetics
     MODULE PROCEDURE ripple_eta_mag
  END INTERFACE

  ! joins all propagators in a fieldripple
  PUBLIC ripple_prop_joiner
  PRIVATE ripple_prop_join
  INTERFACE ripple_prop_joiner
     MODULE PROCEDURE ripple_prop_join
  END INTERFACE


  ! ---------------------------------------------------------------------------
  ! private routines
  ! ---------------------------------------------------------------------------

  ! find_extrema_period in periods
  PRIVATE find_extrema
  INTERFACE find_extrema
     MODULE PROCEDURE find_extrema
  END INTERFACE

  PRIVATE find_next_extremum
  INTERFACE find_next_extremum
     MODULE PROCEDURE find_next_extremum
  END INTERFACE

  PRIVATE find_width
  INTERFACE find_width
     MODULE PROCEDURE find_width
  END INTERFACE

  PRIVATE compute_width
  INTERFACE compute_width
     MODULE PROCEDURE compute_width
  END INTERFACE

  PRIVATE find_poly
  INTERFACE find_poly
     MODULE PROCEDURE find_poly
  END INTERFACE

  PRIVATE eval_poly
  INTERFACE eval_poly
     MODULE PROCEDURE eval_poly
  END INTERFACE

  ! find_extrema_period in periods
  PRIVATE setup_fieldpropagators
  INTERFACE setup_fieldpropagators
     MODULE PROCEDURE setup_fieldpropagators
  END INTERFACE

CONTAINS

  ! ---------------------------------------------------------------------------
  ! make for device_struct
  ! 
  SUBROUTINE make_mag_device(name)
    use neo_magfie, only: magfie_result,magfie_spline,magfie_sarray
    USE magfie_mod, ONLY : stevvo
    USE field_eq_mod, ONLY : rtf
    !! Modifications by Andreas F. Martitsch (18.09.2015)
    ! Used within neo_get_b00 (neo_sub.f90/Boozer coordinates)
    ! to obtain the normalization of the magnetic field (Bref=B_00(s))
    USE neo_actual_fluxs, ONLY: s_es
    !! End Modifications by Andreas F. Martitsch (18.09.2015)
    CHARACTER(len=*),    INTENT(in), OPTIONAL :: name
    ! stevvo related stuff
    REAL(kind=dp) :: r0i,cbfi,bz0i,bf0
    ! name
    CALL construct_magnetics(device)
    IF (      PRESENT(name)) device%name = name
    IF (.NOT. PRESENT(name)) device%name = name_def
    ! this is still in a form I do not like at all
    IF (mag_coordinates .EQ. 0) THEN
       ! cylindrical coordinates
       IF (mag_magfield .EQ. 0) THEN ! homogeneous case
          device%r0  = 200.0_dp
          device%z0  = 0.0_dp
          device%nfp = 1
       ELSEIF (mag_magfield .EQ. 1) THEN ! Biot-Savart
          IF (magnetic_device .EQ. 0) THEN ! Tokamak
             CALL stevvo_0(device%r0,r0i,device%nfp,cbfi,bz0i,bf0)
             device%z0 = 0.0_dp
          ELSEIF (magnetic_device .EQ. 1) THEN ! W7-AS
             CALL stevvo_1(device%r0,r0i,device%nfp,cbfi,bz0i,bf0)
             device%z0 = 0.0_dp
          ELSEIF (magnetic_device .EQ. 2) THEN ! W7-AS
             CALL stevvo_l(device%r0,r0i,device%nfp,cbfi,bz0i,bf0)
             device%z0 = 0.0_dp
          ELSE
             PRINT *, 'Magnetic Device not implemented'
             STOP
          END IF
       ELSEIF (mag_magfield .EQ. 2) THEN ! Legendre
          CALL stevvo_l(device%r0,r0i,device%nfp,cbfi,bz0i,bf0)
          device%z0 = 0.0_dp
       ELSEIF (mag_magfield .EQ. 3) THEN ! EFIT
          magnetic_device = 0 ! Tokamak
          device%r0 = efit_raxis
          device%z0 = efit_zaxis
          device%nfp = 1
       ELSE ! does not exist
          PRINT *, 'Not implemented: mag_magfield = ',mag_magfield
          STOP          
       END IF
    ELSE
       ! Boozer
       magfie_result = 0
       magfie_spline = 1
       ALLOCATE(magfie_sarray(1))
       magfie_sarray = boozer_s
       !**********************************************************
       ! For neo_fourier() consitency check
       !**********************************************************
       s_es = boozer_s
       write (*,*) "Flux surface: ", s_es
       !**********************************************************
       CALL stevvo(device%r0,r0i,device%nfp,cbfi,bz0i,bf0)
       device%z0  = 0.0_dp
       boozer_bmod0 = bf0
    END IF
  END SUBROUTINE make_mag_device
  ! ---------------------------------------------------------------------------
  
  ! ---------------------------------------------------------------------------
  ! make for surface_struct
  SUBROUTINE make_mag_surface( &
       bmod0,nperiod_in,nstep_in,ndim_in)
    REAL(kind=dp),           INTENT(in) :: bmod0
    INTEGER,       OPTIONAL, INTENT(in) :: nperiod_in
    INTEGER,       OPTIONAL, INTENT(in) :: nstep_in
    INTEGER,       OPTIONAL, INTENT(in) :: ndim_in

    CALL construct_magnetics(device,surface)
    surface%bmod0 = bmod0
    IF (.NOT. PRESENT(nperiod_in)) surface%nperiod = nperiod_def
    IF (      PRESENT(nperiod_in)) surface%nperiod = nperiod_in
    IF (.NOT. PRESENT(nstep_in)  ) surface%nstep   = nstep_def
    IF (      PRESENT(nstep_in)  ) surface%nstep   = nstep_in
    IF (.NOT. PRESENT(ndim_in)   ) surface%ndim    = ndim_def
    IF (      PRESENT(ndim_in)   ) surface%ndim    = ndim_in
  END SUBROUTINE make_mag_surface
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  ! make for fieldline_struct
  ! 
  SUBROUTINE make_mag_fieldline_newperiod(xstart)

    USE rk4_kin_mod, ONLY : y
    USE plagrange_mod
    use neo_magfie, only : boozer_iota
    use collisionality_mod, only : isw_axisymm
    use field_eq_mod, only : dpsidr, dpsidz
    !! Modifications by Andreas F. Martitsch (11.06.2014)
    ! Optional output (necessary for modeling the magnetic rotation)
    ! Note: This requires changes in "modify_propagator"
    ! (flint.f90; block using the routine commented out) and
    ! "mag_interface_mod" (mag_interface.f90).
    INTERFACE
       SUBROUTINE magdata_for_particles(phi,bhat,geodcu,h_phi,dlogbdphi,&
            bcovar_s_hat,dlogbds,dbcovar_s_hat_dphi)
         DOUBLE PRECISION, INTENT(in)            :: phi
         DOUBLE PRECISION, INTENT(out)           :: geodcu,bhat,h_phi,dlogbdphi
         DOUBLE PRECISION, OPTIONAL, INTENT(out) :: bcovar_s_hat, &
              dlogbds, dbcovar_s_hat_dphi
       END SUBROUTINE magdata_for_particles
    END INTERFACE
    !! End Modifications by Andreas F. Martitsch (11.06.2014)

    REAL(kind=dp),                INTENT(in) :: xstart(3)
    
    ! local
    INTEGER       :: nperiod,nstep,ndim,nfp
    INTEGER       :: i_period,i

    INTEGER       :: u1 = 117

    REAL(kind=dp) :: phibeg,phi
    REAL(kind=dp) :: aiota
    REAL(kind=dp) :: dist_min_req,dist
    REAL(kind=dp) :: phimi,phima,h
    REAL(kind=dp) :: r_start,z_start,theta_start,theta_b,theta_e

    REAL(kind=dp) :: x1_start,x2_start,x3_start,bhat_start,geodcu_start,h_phi_start,dlogbdphi_start
    REAL(kind=dp) :: x1_end,x2_end,x3_end,bhat_end,geodcu_end,h_phi_end,dlogbdphi_end
    REAL(kind=dp) :: x1_mid,x3_mid,bhat_mid,geodcu_mid,h_phi_mid,dlogbdphi_mid
    !! Modifications by Andreas F. Martitsch (11.06.2014)
    ! Optional output (necessary for modeling the magnetic rotation)
    REAL(kind=dp) :: dbcovar_s_hat_dphi_start, dbcovar_s_hat_dphi_end, dbcovar_s_hat_dphi_mid
    REAL(kind=dp) :: bcovar_s_hat_start, bcovar_s_hat_end, bcovar_s_hat_mid
    REAL(kind=dp) :: dlogbds_start, dlogbds_end, dlogbds_mid
    !! End Modifications by Andreas F. Martitsch (11.06.2014)
    REAL(kind=dp) :: phi_start,phi_end,phi_span,phi_mult
    REAL(kind=dp) :: one_over_h_phi_int,b_over_h_phi_int,b2_over_h_phi_int

    INTEGER       :: count_sum_bhat,uplim
    INTEGER       :: back_period,back_period_end,i_start
    TYPE(fieldperiod_struct),     POINTER :: fieldperiod_delete

    r_start = 1.234e+5
    z_start = 1.234e+5
    x1_start = 1.234e+5
    x2_start = 1.234e+5
    x3_start = 1.234e+5
    bhat_start = 1.234e+5
    geodcu_start = 1.234e+5
    h_phi_start = 1.234e+5
    dbcovar_s_hat_dphi_start = 1.234e+5
    bcovar_s_hat_start = 1.234e+5
    dlogbdphi_start = 1.234e+5
    dlogbds_start = 1.234e+5

    ! local copies
    nperiod = surface%nperiod
    nstep   = surface%nstep  
    ndim    = surface%ndim   
    nfp     = surface%parent%nfp

    ! y-vector at the beginning
    IF (mag_coordinates .EQ. 0) THEN
       ! cylindrical coordinates
       y(1) = xstart(1)
       phibeg  = xstart(2)
       y(2) = xstart(3)
       y(3) = 1.0_dp
       y(4:ndim) = 0.0_dp
       if (mag_magfield .eq. 3) then  ! ASDEX EFIT
          call rk4_kin(phibeg, 0.d0)
          y(3) = dpsidr
          y(5) = dpsidz
       end if
    ELSE
       ! Boozer
       y(1) = xstart(3) ! boozer_theta_beg
       phibeg  = xstart(2)
       y(2) = 1.0_dp
       y(3) = 1.0_dp
       y(4:ndim) = 0.0_dp
    END IF
       
    ! construct a fieldline
    CALL construct_magnetics(surface,fieldline)
    fieldline%xstart = xstart
    IF (mag_coordinates .EQ. 0) THEN
       ! cylindrical
       IF (mag_magfield .NE. 0 .AND. magnetic_device .EQ. 0) THEN ! Tokamak
          surface%aiota = aiota_tokamak
          aiota = aiota_tokamak
       ELSE
          surface%aiota       =  0.0_dp
          aiota               =  0.0_dp
       END IF
       dist_min_req        =  1.0_dp
    ELSE
       ! Boozer
       aiota = boozer_iota
       surface%aiota = aiota
       dist_min_req        =  0.1_dp
    END IF

    fieldline%b_abs_max = -1.0d100
    fieldline%b_abs_min =  1.0d100

    surface%r_max       = -1.0d100
    surface%r_min       =  1.0d100
    surface%z_max       = -1.0d100
    surface%z_min       =  1.0d100


    ! construct periods
    ! period and steplength (private variables)
    IF (mag_coordinates .EQ. 0) THEN
       ! cylindrical
       IF (mag_magfield .NE. 0 .AND. magnetic_device .EQ. 0) THEN ! Tokamak
          IF (isw_axisymm .EQ. 0) THEN
             period_length   = 2.0_dp*pi/aiota_tokamak/CEILING(1.0_dp/aiota_tokamak)
          ELSEIF (isw_axisymm .EQ. 1) THEN
             ! the number 5 ensures that there are at least two periods
             ! otherwise there is a problem with logics regarding extra 
             ! children of a period (there is only one possible)
             period_length   = 5.0_dp*pi/aiota_tokamak
             split_inflection_points = .FALSE. 
          ELSE
             PRINT *, 'isw_axisymm = ',isw_axisymm,' not implemented!'
             STOP
          END IF
       ELSE ! Stellarator, homogeneous
          period_length   = 2.0_dp*pi/DBLE(ABS(nfp))
       END IF
    ELSE
       ! Boozer
       IF (magnetic_device .EQ. 0) THEN ! Tokamak
          IF (isw_axisymm .EQ. 0) THEN
             period_length   = 2.0_dp*pi/boozer_iota/CEILING(1.0_dp/boozer_iota)
          ELSEIF (isw_axisymm .EQ. 1) THEN
             ! the number 5 ensures that there are at least two periods
             ! otherwise there is a problem with logics regarding extra 
             ! children of a period (there is only one possible)
             period_length   = 5.0_dp*pi/abs(boozer_iota)
             split_inflection_points = .FALSE.
          ELSE
             PRINT *, 'isw_axisymm = ',isw_axisymm,' not implemented!'
             STOP
          END IF
       ELSE
          period_length   = 2.0_dp*pi/DBLE(ABS(nfp))
       END IF
    END IF
    phimi = phibeg
    phima = phimi + period_length
    h=(phima-phimi)/nstep

    i_period = 0
    phi = phibeg
    IF (mag_coordinates .EQ. 0) THEN
       ! cylindrical
       r_start = y(1)
       z_start = y(2)
       theta_start = ATAN2(y(2),y(1)-device%r0)
    ELSE
       ! Boozer
       theta_start = y(1)
    END IF

    construct_periods: DO
       i_period = i_period + 1
       ! preparation
       CALL construct_magnetics(fieldline,fieldperiod) ! construction
       ALLOCATE(fieldperiod%coords)
       ALLOCATE(fieldperiod%coords%x1(0:nstep))
       ALLOCATE(fieldperiod%coords%x2(0:nstep))
       ALLOCATE(fieldperiod%coords%x3(0:nstep))
       ALLOCATE(fieldperiod%mdata)
       ALLOCATE(fieldperiod%mdata%bhat(0:nstep))
       ALLOCATE(fieldperiod%mdata%geodcu(0:nstep))
       ALLOCATE(fieldperiod%mdata%h_phi(0:nstep))
       ALLOCATE(fieldperiod%mdata%dlogbdphi(0:nstep))
       !! Modifications by Andreas F. Martitsch (11.06.2014)
       ! Optional output (necessary for modeling the magnetic rotation) ->
       ! Allocate the additional entries
       ALLOCATE(fieldperiod%mdata%dbcovar_s_hat_dphi(0:nstep))
       ALLOCATE(fieldperiod%mdata%bcovar_s_hat(0:nstep))
       ALLOCATE(fieldperiod%mdata%dlogbds(0:nstep))
       !! End Modifications by Andreas F. Martitsch (11.06.2014)
       ALLOCATE(fieldperiod%mdata%ybeg(1:ndim))
       ALLOCATE(fieldperiod%mdata%yend(1:ndim))
       ! first storage
       fieldperiod%mdata%ybeg = y
       fieldperiod%phi_l = phi
       IF (mag_coordinates .EQ. 0) THEN
          ! cylindrical       
          theta_b = ATAN2(y(2),y(1)-device%r0)
       ELSE
          ! Boozer
          theta_b = y(1)
       END IF
       fieldperiod%theta_b = theta_b
       !! Modifications by Andreas F. Martitsch (11.06.2014)
       ! Optional output (necessary for modeling the magnetic rotation) ->
       ! Compute the additional entries
       CALL magdata_for_particles(                  &
            phi,                                    & 
            fieldperiod%mdata%bhat(0),              &
            fieldperiod%mdata%geodcu(0),            &
            fieldperiod%mdata%h_phi(0),             &
            fieldperiod%mdata%dlogbdphi(0),         &
            fieldperiod%mdata%bcovar_s_hat(0),      &
            fieldperiod%mdata%dlogbds(0),           &
            fieldperiod%mdata%dbcovar_s_hat_dphi(0) &
            )

       !! End Modifications by Andreas F. Martitsch (11.06.2014)
       IF (mag_coordinates .EQ. 0) THEN
          ! cylindrical
          fieldperiod%coords%x1(0) = y(1)
          fieldperiod%coords%x2(0) = phi
          fieldperiod%coords%x3(0) = y(2)
       ELSE
          ! Boozer%
          fieldperiod%coords%x1(0) = boozer_s
          fieldperiod%coords%x2(0) = phi
          fieldperiod%coords%x3(0) = y(1) ! boozer_theta
       END IF

       IF (i_period .EQ. 1) THEN
          x1_start                 = fieldperiod%coords%x1(0)
          x2_start                 = fieldperiod%coords%x2(0)
          x3_start                 = fieldperiod%coords%x3(0)
          bhat_start               = fieldperiod%mdata%bhat(0)
          geodcu_start             = fieldperiod%mdata%geodcu(0)
          h_phi_start              = fieldperiod%mdata%h_phi(0)
          dlogbdphi_start          = fieldperiod%mdata%dlogbdphi(0)
          !! Modifications by Andreas F. Martitsch (11.06.2014)
          ! Optional output (necessary for modeling the magnetic rotation)
          bcovar_s_hat_start       = fieldperiod%mdata%bcovar_s_hat(0)
          dlogbds_start            = fieldperiod%mdata%dlogbds(0)
          dbcovar_s_hat_dphi_start = fieldperiod%mdata%dbcovar_s_hat_dphi(0)
          !! End Modifications by Andreas F. Martitsch (11.06.2014)
       END IF

       IF (mag_coordinates .EQ. 0) THEN
          ! cylindrical
          surface%r_max = MAX(surface%r_max,y(1))
          surface%r_min = MIN(surface%r_min,y(1))
          surface%z_max = MAX(surface%z_max,y(2))
          surface%z_min = MIN(surface%z_min,y(2))
       END IF

       ! all steps within period
       all_steps: DO i = 1, nstep
          CALL rk4_kin(phi,h)
          !! Modifications by Andreas F. Martitsch (11.06.2014)
          ! Optional output (necessary for modeling the magnetic rotation) ->
          ! Compute the additional entries
          CALL magdata_for_particles(                  &
               phi,                                    & 
               fieldperiod%mdata%bhat(i),              &
               fieldperiod%mdata%geodcu(i),            &
               fieldperiod%mdata%h_phi(i),             &
               fieldperiod%mdata%dlogbdphi(i),         &
               fieldperiod%mdata%bcovar_s_hat(i),      &
               fieldperiod%mdata%dlogbds(i),           &
               fieldperiod%mdata%dbcovar_s_hat_dphi(i) &
               )

          !! End Modifications by Andreas F. Martitsch (11.06.2014)
          if (mag_coordinates .eq. 0) then
             ! cylindrical
             fieldperiod%coords%x1(i) = y(1)
             fieldperiod%coords%x2(i) = phi
             fieldperiod%coords%x3(i) = y(2)
             surface%r_max = MAX(surface%r_max,y(1))
             surface%r_min = MIN(surface%r_min,y(1))
             surface%z_max = MAX(surface%z_max,y(2))
             surface%z_min = MIN(surface%z_min,y(2))
          ELSE
             ! Boozer
             IF (y(1) .GT. 2.0_dp*pi) y(1) = y(1) - 2.0_dp*pi
             IF (y(1) < 0.d0) y(1) = y(1) + 2.0_dp*pi
             fieldperiod%coords%x1(i) = boozer_s
             fieldperiod%coords%x2(i) = phi
             fieldperiod%coords%x3(i) = y(1)
          end if

       END DO all_steps
       ! final storage 
       fieldperiod%mdata%yend = y
       fieldperiod%phi_r = phi
       
       IF (mag_coordinates .EQ. 0) THEN
          ! cylindrical
          ! distance and aiota
          dist = SQRT( (y(1)-r_start)**2 + (y(2)-z_start)**2 )
          theta_e = ATAN2(y(2),y(1)-device%r0)
          IF (.NOT. (mag_magfield .NE. 0 .AND. magnetic_device .EQ. 0)) THEN ! no Tokamak       
             aiota = aiota + theta_e - theta_b
             IF(theta_e .LT. theta_b) aiota = aiota + 2.0_dp*pi
          END IF
       ELSE
          ! Boozer
          theta_e = y(1)
          IF (theta_e .GE. 2.0_dp*pi) theta_e = theta_e - 2.0_dp*pi
          dist = ABS(theta_e - theta_start)
          IF (dist .GE. 2.0_dp*pi) dist = dist - 2.0_dp*pi
          dist = MIN(dist,2.0_dp*pi-dist)
       END IF

       ! final decision
       IF (mag_coordinates .EQ. 0) THEN
          ! cylindrical
          IF (mag_magfield .EQ. 0) THEN ! homogeneous
             IF (i_period .EQ. 3) THEN ! we need at least three periods
                EXIT construct_periods
             END IF
          ELSEIF (mag_magfield .NE. 0 .AND. magnetic_device .EQ. 0) THEN ! Tokamak
             IF ( ABS(theta_start - theta_e) .LT. 1.0d-3 .OR. &
                  ABS(theta_start - theta_e + 2.0_dp*pi) .LT. 1.0d-3 .OR. &
                  ABS(theta_start - theta_e - 2.0_dp*pi) .LT. 1.0d-3 ) THEN
                EXIT construct_periods
             END IF
          ELSE ! no Tokamak
             IF (i_period .LE. mag_nperiod_min) THEN 
                dist_min_req = MIN(dist_min_req,dist)
             ELSE
                IF (dist .LT. dist_min_req) THEN
                   EXIT construct_periods
                END IF
             END IF
          END IF ! end tok - no tok
       ELSE
          ! Boozer
          IF (magnetic_device .EQ. 0) THEN ! Tokamak
             IF ( ABS(theta_start - theta_e) .LT. 1.0d-3 .OR. &
                  ABS(theta_start - theta_e + 2.0_dp*pi) .LT. 1.0d-3 .OR. &
                  ABS(theta_start - theta_e - 2.0_dp*pi) .LT. 1.0d-3 ) THEN
                EXIT construct_periods
             END IF
          ELSE ! no Tokamak
             IF (i_period .LE. mag_nperiod_min) THEN 
                dist_min_req = MIN(dist_min_req,dist)
             ELSE
                IF (dist .LT. dist_min_req) THEN
                   EXIT construct_periods
                END IF
             END IF
          END IF ! end tok - no tok          
       END IF
    END DO construct_periods

    ! end points for fixing periodicity
    x1_end                 = fieldperiod%coords%x1(nstep)
    x2_end                 = fieldperiod%coords%x2(nstep)
    x3_end                 = fieldperiod%coords%x3(nstep)
    bhat_end               = fieldperiod%mdata%bhat(nstep)
    geodcu_end             = fieldperiod%mdata%geodcu(nstep)
    h_phi_end              = fieldperiod%mdata%h_phi(nstep)
    dlogbdphi_end          = fieldperiod%mdata%dlogbdphi(nstep)
    !! Modifications by Andreas F. Martitsch (11.06.2014)
    ! Optional output (necessary for modeling the magnetic rotation)
    bcovar_s_hat_end       = fieldperiod%mdata%bcovar_s_hat(nstep)
    dlogbds_end            = fieldperiod%mdata%dlogbds(nstep)
    dbcovar_s_hat_dphi_end = fieldperiod%mdata%dbcovar_s_hat_dphi(nstep)
    !! End Modifications by Andreas F. Martitsch (11.06.2014)
    ! output
    IF (mag_coordinates .EQ. 0) THEN
       ! cylindrical
       IF (.NOT. (mag_magfield .NE. 0 .AND. magnetic_device .EQ. 0)) THEN ! no Tokamak       
          surface%aiota = aiota / (phi - phibeg)
       END IF
    END IF
    PRINT *, '-------------------------------------------------'
    PRINT *, ' Final Joining Point Found!'
    PRINT *, '  period:                  ', i_period
    PRINT *, '  aiota:                   ', surface%aiota
    IF (mag_magfield .NE. 0 .AND. magnetic_device .EQ. 0) THEN ! Tokamak
       PRINT *, '  theta_start, theta_end:  ',theta_start,theta_e, &
            ' (',theta_e + 2.0_dp*pi,theta_e - 2.0_dp*pi,')'
    ELSE ! no Tokamak
       PRINT *, '  dist, required dist:     ', dist,dist_min_req
    END IF

    ! This is the new stuff with going back in fieldperiods
    IF (mag_symmetric) THEN
       fieldline%ch_act => fieldline%ch_fir ! go back from the first
       fieldperiod => fieldline%ch_act
       y = fieldperiod%mdata%ybeg
       phi = fieldperiod%phi_l
       IF (mag_symmetric_shorten) THEN
          back_period_end = INT(DBLE(i_period)/2.0d0)
       ELSE
          back_period_end = i_period
       END IF
       construct_periods_back: DO back_period = 1,back_period_end
          CALL construct_magnetics(fieldline,fieldperiod,-1) ! construction backward
          ALLOCATE(fieldperiod%coords)
          ALLOCATE(fieldperiod%coords%x1(0:nstep))
          ALLOCATE(fieldperiod%coords%x2(0:nstep))
          ALLOCATE(fieldperiod%coords%x3(0:nstep))
          ALLOCATE(fieldperiod%mdata)
          ALLOCATE(fieldperiod%mdata%bhat(0:nstep))
          ALLOCATE(fieldperiod%mdata%geodcu(0:nstep))
          ALLOCATE(fieldperiod%mdata%h_phi(0:nstep))
          ALLOCATE(fieldperiod%mdata%dlogbdphi(0:nstep))
          !! Modifications by Andreas F. Martitsch (11.06.2014)
          ! Optional output (necessary for modeling the magnetic rotation) ->
          ! Allocate the additional entries
          ALLOCATE(fieldperiod%mdata%dbcovar_s_hat_dphi(0:nstep))
          ALLOCATE(fieldperiod%mdata%bcovar_s_hat(0:nstep))
          ALLOCATE(fieldperiod%mdata%dlogbds(0:nstep))
          !! End Modifications by Andreas F. Martitsch (11.06.2014)
          ALLOCATE(fieldperiod%mdata%ybeg(1:ndim))
          ALLOCATE(fieldperiod%mdata%yend(1:ndim))
          ! first storage
          fieldperiod%mdata%yend = y
          fieldperiod%phi_r = phi
          IF (mag_coordinates .EQ. 0) THEN
             ! cylindrical       
             theta_b = ATAN2(y(2),y(1)-device%r0)
          ELSE
             ! Boozer
             theta_b = y(1)
          END IF
          fieldperiod%theta_b = theta_b ! theta_b is now on the right side

          !! Modifications by Andreas F. Martitsch (11.06.2014)
          ! Optional output (necessary for modeling the magnetic rotation) ->
          ! Compute the additional entries
          CALL magdata_for_particles(                  &
               phi,                                    & 
               fieldperiod%mdata%bhat(nstep),              &
               fieldperiod%mdata%geodcu(nstep),            &
               fieldperiod%mdata%h_phi(nstep),             &
               fieldperiod%mdata%dlogbdphi(nstep),         &
               fieldperiod%mdata%bcovar_s_hat(nstep),      &
               fieldperiod%mdata%dlogbds(nstep),           &
               fieldperiod%mdata%dbcovar_s_hat_dphi(nstep) &
               )

          !! End Modifications by Andreas F. Martitsch (11.06.2014)
          IF (mag_coordinates .EQ. 0) THEN
             ! cylindrical
             fieldperiod%coords%x1(nstep) = y(1)
             fieldperiod%coords%x2(nstep) = phi
             fieldperiod%coords%x3(nstep) = y(2)
          ELSE
             ! Boozer
             fieldperiod%coords%x1(nstep) = boozer_s
             fieldperiod%coords%x2(nstep) = phi
             fieldperiod%coords%x3(nstep) = y(1) ! boozer_theta
          END IF

          IF (back_period .EQ. i_period) THEN
             x1_start                 = fieldperiod%coords%x1(0)
             x2_start                 = fieldperiod%coords%x2(0)
             x3_start                 = fieldperiod%coords%x3(0)
             bhat_start               = fieldperiod%mdata%bhat(0)
             geodcu_start             = fieldperiod%mdata%geodcu(0)
             h_phi_start              = fieldperiod%mdata%h_phi(0)
             dlogbdphi_start          = fieldperiod%mdata%dlogbdphi(0)
             !! Modifications by Andreas F. Martitsch (11.06.2014)
             ! Optional output (necessary for modeling the magnetic rotation)
             bcovar_s_hat_start       = fieldperiod%mdata%bcovar_s_hat(0)
             dlogbds_start            = fieldperiod%mdata%dlogbds(0)
             dbcovar_s_hat_dphi_start = fieldperiod%mdata%dbcovar_s_hat_dphi(0)
             !! End Modifications by Andreas F. Martitsch (11.06.2014)
          END IF

          IF (mag_coordinates .EQ. 0) THEN
             ! cylindrical
             surface%r_max = MAX(surface%r_max,y(1))
             surface%r_min = MIN(surface%r_min,y(1))
             surface%z_max = MAX(surface%z_max,y(2))
             surface%z_min = MIN(surface%z_min,y(2))
          END IF

          ! all steps within period - back
          all_steps_back: DO i = nstep-1,0,-1
             CALL rk4_kin(phi,-h)
             !! Modifications by Andreas F. Martitsch (11.06.2014)
             ! Optional output (necessary for modeling the magnetic rotation) ->
             ! Compute the additional entries
             CALL magdata_for_particles(                  &
                  phi,                                    & 
                  fieldperiod%mdata%bhat(i),              &
                  fieldperiod%mdata%geodcu(i),            &
                  fieldperiod%mdata%h_phi(i),             &
                  fieldperiod%mdata%dlogbdphi(i),         &
                  fieldperiod%mdata%bcovar_s_hat(i),      &
                  fieldperiod%mdata%dlogbds(i),           &
                  fieldperiod%mdata%dbcovar_s_hat_dphi(i) &
                  )

             !! End Modifications by Andreas F. Martitsch (11.06.2014)
             IF (mag_coordinates .EQ. 0) THEN
                ! cylindrical
                fieldperiod%coords%x1(i) = y(1)
                fieldperiod%coords%x2(i) = phi
                fieldperiod%coords%x3(i) = y(2)
                surface%r_max = MAX(surface%r_max,y(1))
                surface%r_min = MIN(surface%r_min,y(1))
                surface%z_max = MAX(surface%z_max,y(2))
                surface%z_min = MIN(surface%z_min,y(2))
             ELSE
                ! Boozer
                IF (y(1) .GT. 2.0_dp*pi) y(1) = y(1) - 2.0_dp*pi
                fieldperiod%coords%x1(i) = boozer_s
                fieldperiod%coords%x2(i) = phi
                fieldperiod%coords%x3(i) = y(1)
             END IF
          END DO all_steps_back
          ! final storage 
          fieldperiod%mdata%ybeg = y
          fieldperiod%phi_l = phi

       END DO construct_periods_back

       IF (mag_symmetric_shorten) THEN
          fieldperiod => fieldline%ch_fir
          DO i = 1,i_period-1
             fieldperiod => fieldperiod%next
          END DO
          destruct_mag: DO
             IF (ASSOCIATED(fieldperiod%next)) THEN
                fieldperiod_delete => fieldperiod%next
                CALL destruct_magnetics(fieldperiod_delete)
             ELSE
                EXIT destruct_mag
             END IF
          END DO destruct_mag
       ELSE   
          i_period = 2 * i_period
       END IF
       ! fix period_tag
       fieldperiod => fieldline%ch_fir
       fix_period_tag: DO
          IF ( ASSOCIATED(fieldperiod%prev) ) THEN
             fieldperiod%tag = fieldperiod%prev%tag + 1
          ELSE
             fieldperiod%tag = 1
          END IF
          IF ( ASSOCIATED(fieldperiod%next) ) THEN
             fieldperiod => fieldperiod%next
          ELSE
             EXIT fix_period_tag
          END IF
       END DO fix_period_tag
       ! now they are numbered from 1

    END IF
    ! This is the new stuff with going back in fieldperiods - End

    ! fix periodicity
    IF (mag_symmetric) THEN ! new symmetric case
       fieldperiod => fieldline%ch_fir ! first period
       x1_start        = fieldperiod%coords%x1(0)
       x2_start        = fieldperiod%coords%x2(0)
       x3_start        = fieldperiod%coords%x3(0)
       bhat_start      = fieldperiod%mdata%bhat(0)
       geodcu_start    = fieldperiod%mdata%geodcu(0)
       h_phi_start     = fieldperiod%mdata%h_phi(0)
       dlogbdphi_start = fieldperiod%mdata%dlogbdphi(0)
       !! Modifications by Andreas F. Martitsch (11.06.2014)
       ! Optional output (necessary for modeling the magnetic rotation)
       bcovar_s_hat_start       = fieldperiod%mdata%bcovar_s_hat(0)
       dlogbds_start            = fieldperiod%mdata%dlogbds(0)
       dbcovar_s_hat_dphi_start = fieldperiod%mdata%dbcovar_s_hat_dphi(0)
       !! End Modifications by Andreas F. Martitsch (11.06.2014)
       fieldperiod => fieldline%ch_las ! last period
       x1_end        = fieldperiod%coords%x1(nstep)
       x2_end        = fieldperiod%coords%x2(nstep)
       x3_end        = fieldperiod%coords%x3(nstep)
       bhat_end      = fieldperiod%mdata%bhat(nstep)
       geodcu_end    = fieldperiod%mdata%geodcu(nstep)
       h_phi_end     = fieldperiod%mdata%h_phi(nstep)
       dlogbdphi_end = fieldperiod%mdata%dlogbdphi(nstep)
       !! Modifications by Andreas F. Martitsch (11.06.2014)
       ! Optional output (necessary for modeling the magnetic rotation)
       bcovar_s_hat_end       = fieldperiod%mdata%bcovar_s_hat(nstep)
       dlogbds_end            = fieldperiod%mdata%dlogbds(nstep)
       dbcovar_s_hat_dphi_end = fieldperiod%mdata%dbcovar_s_hat_dphi(nstep)
       !! End Modifications by Andreas F. Martitsch (11.06.2014)
       
       ! new end points for both sides
       x1_mid = (x1_start + x1_end) / 2.0d0
       x3_mid = (x3_start + x3_end) / 2.0d0
       bhat_mid = (bhat_start + bhat_end) / 2.0d0
       geodcu_mid = (geodcu_start + geodcu_end) / 2.0d0
       h_phi_mid = (h_phi_start + h_phi_end) / 2.0d0
       dlogbdphi_mid = (dlogbdphi_start + dlogbdphi_end) / 2.0d0
       !! Modifications by Andreas F. Martitsch (11.06.2014)
       ! Optional output (necessary for modeling the magnetic rotation)
       bcovar_s_hat_mid       = &
            ( bcovar_s_hat_start + bcovar_s_hat_end ) / 2.0d0
       dlogbds_mid            = &
            ( dlogbds_start + dlogbds_end ) / 2.0d0
       dbcovar_s_hat_dphi_mid = &
            ( dbcovar_s_hat_dphi_start  + dbcovar_s_hat_dphi_end ) / 2.0d0
       !! End Modifications by Andreas F. Martitsch (11.06.2014)

       ! fix last period
       fieldperiod => fieldline%ch_las ! last period
       phi_start = fieldperiod%coords%x2(0)
       phi_end   = x2_end
       phi_span  = phi_end - phi_start
       DO i = 1, nstep
          phi_mult = (fieldperiod%coords%x2(i) - phi_start) / phi_span
          fieldperiod%coords%x1(i) = fieldperiod%coords%x1(i) + &
               (x1_mid - x1_end) * phi_mult
          fieldperiod%coords%x3(i) = fieldperiod%coords%x3(i) + &
               (x3_mid - x3_end) * phi_mult
          fieldperiod%mdata%bhat(i) = fieldperiod%mdata%bhat(i) + &
               (bhat_mid - bhat_end) * phi_mult
          fieldperiod%mdata%geodcu(i) = fieldperiod%mdata%geodcu(i) + &
               (geodcu_mid - geodcu_end) * phi_mult
          fieldperiod%mdata%h_phi(i) = fieldperiod%mdata%h_phi(i) + &
               (h_phi_mid - h_phi_end) * phi_mult
          fieldperiod%mdata%dlogbdphi(i) = fieldperiod%mdata%dlogbdphi(i) + &
               (dlogbdphi_mid - dlogbdphi_end) * phi_mult
          !! Modifications by Andreas F. Martitsch (11.06.2014)
          ! Optional output (necessary for modeling the magnetic rotation)
          fieldperiod%mdata%bcovar_s_hat(i) = &
               fieldperiod%mdata%bcovar_s_hat(i) + &
               (bcovar_s_hat_mid - bcovar_s_hat_end) * phi_mult
          fieldperiod%mdata%dlogbds(i) = &
               fieldperiod%mdata%dlogbds(i) + &
               (dlogbds_mid - dlogbds_end) * phi_mult
          fieldperiod%mdata%dbcovar_s_hat_dphi(i) = &
               fieldperiod%mdata%dbcovar_s_hat_dphi(i) + &
               (dbcovar_s_hat_dphi_mid - dbcovar_s_hat_dphi_end) * phi_mult
          !! End Modifications by Andreas F. Martitsch (11.06.2014)
       END DO
       ! fix first period
       fieldperiod => fieldline%ch_fir ! first period
       phi_start = x2_start
       phi_end   = fieldperiod%coords%x2(nstep)
       phi_span  = phi_end - phi_start
       DO i = 0, nstep-1
          phi_mult = (phi_end - fieldperiod%coords%x2(i)) / phi_span
          fieldperiod%coords%x1(i) = fieldperiod%coords%x1(i) + &
               (x1_mid - x1_start) * phi_mult
          fieldperiod%coords%x3(i) = fieldperiod%coords%x3(i) + &
               (x3_mid - x3_start) * phi_mult
          fieldperiod%mdata%bhat(i) = fieldperiod%mdata%bhat(i) + &
               (bhat_mid - bhat_start) * phi_mult
          fieldperiod%mdata%geodcu(i) = fieldperiod%mdata%geodcu(i) + &
               (geodcu_mid - geodcu_start) * phi_mult
          fieldperiod%mdata%h_phi(i) = fieldperiod%mdata%h_phi(i) + &
               (h_phi_mid - h_phi_start) * phi_mult
          fieldperiod%mdata%dlogbdphi(i) = fieldperiod%mdata%dlogbdphi(i) + &
               (dlogbdphi_mid - dlogbdphi_start) * phi_mult
          !! Modifications by Andreas F. Martitsch (11.06.2014)
          ! Optional output (necessary for modeling the magnetic rotation)
          fieldperiod%mdata%bcovar_s_hat(i) = &
               fieldperiod%mdata%bcovar_s_hat(i) + &
               (bcovar_s_hat_mid - bcovar_s_hat_start) * phi_mult
          fieldperiod%mdata%dlogbds(i) = &
               fieldperiod%mdata%dlogbds(i) + &
               (dlogbds_mid - dlogbds_start) * phi_mult
          fieldperiod%mdata%dbcovar_s_hat_dphi(i) = &
               fieldperiod%mdata%dbcovar_s_hat_dphi(i) + &
               (dbcovar_s_hat_dphi_mid - dbcovar_s_hat_dphi_start) * phi_mult
          !! End Modifications by Andreas F. Martitsch (11.06.2014)
       END DO

    ELSE ! normal not symmetric case
       IF (mag_close_fieldline .EQ. 1) THEN 
          ! modify all
          phi_start = x2_start
          phi_end   = x2_end
          phi_span  = phi_end - phi_start
          fieldperiod => fieldline%ch_fir ! first period
          DO ! go through all periods
             DO i = 0, nstep
                phi_mult = (fieldperiod%coords%x2(i) - phi_start) / phi_span
                fieldperiod%coords%x1(i) = fieldperiod%coords%x1(i) + &
                     (x1_start - x1_end) * phi_mult
                fieldperiod%coords%x3(i) = fieldperiod%coords%x3(i) + &
                     (x3_start - x3_end) * phi_mult
                fieldperiod%mdata%bhat(i) = fieldperiod%mdata%bhat(i) + &
                     (bhat_start - bhat_end) * phi_mult
                fieldperiod%mdata%geodcu(i) = fieldperiod%mdata%geodcu(i) + &
                     (geodcu_start - geodcu_end) * phi_mult
                fieldperiod%mdata%h_phi(i) = fieldperiod%mdata%h_phi(i) + &
                     (h_phi_start - h_phi_end) * phi_mult
                fieldperiod%mdata%dlogbdphi(i) = fieldperiod%mdata%dlogbdphi(i) + &
                     (dlogbdphi_start - dlogbdphi_end) * phi_mult
                !! Modifications by Andreas F. Martitsch (11.06.2014)
                ! Optional output (necessary for modeling the magnetic rotation)
                fieldperiod%mdata%bcovar_s_hat(i) = &
                     fieldperiod%mdata%bcovar_s_hat(i) + &
                     (bcovar_s_hat_start - bcovar_s_hat_end) * phi_mult
                fieldperiod%mdata%dlogbds(i) = &
                     fieldperiod%mdata%dlogbds(i) + &
                     (dlogbds_start - dlogbds_end) * phi_mult
                fieldperiod%mdata%dbcovar_s_hat_dphi(i) = &
                     fieldperiod%mdata%dbcovar_s_hat_dphi(i) + &
                     (dbcovar_s_hat_dphi_start - dbcovar_s_hat_dphi_end) * phi_mult               
                !! End Modifications by Andreas F. Martitsch (11.06.2014)
             END DO
             IF (.NOT.(ASSOCIATED(fieldperiod%next))) EXIT
             fieldperiod => fieldperiod%next
          END DO
       ELSEIF (mag_close_fieldline .EQ. 2) THEN 
          ! modify only last period
          fieldperiod => fieldline%ch_las ! last period

          phi_start = fieldperiod%coords%x2(0)
          phi_end   = x2_end
          phi_span  = phi_end - phi_start
          DO i = 1, nstep
             phi_mult = (fieldperiod%coords%x2(i) - phi_start) / phi_span
             fieldperiod%coords%x1(i) = fieldperiod%coords%x1(i) + &
                  (x1_start - x1_end) * phi_mult
             fieldperiod%coords%x3(i) = fieldperiod%coords%x3(i) + &
                  (x3_start - x3_end) * phi_mult
             fieldperiod%mdata%bhat(i) = fieldperiod%mdata%bhat(i) + &
                  (bhat_start - bhat_end) * phi_mult
             fieldperiod%mdata%geodcu(i) = fieldperiod%mdata%geodcu(i) + &
                  (geodcu_start - geodcu_end) * phi_mult
             fieldperiod%mdata%h_phi(i) = fieldperiod%mdata%h_phi(i) + &
                  (h_phi_start - h_phi_end) * phi_mult
             fieldperiod%mdata%dlogbdphi(i) = fieldperiod%mdata%dlogbdphi(i) + &
                  (dlogbdphi_start - dlogbdphi_end) * phi_mult
             !! Modifications by Andreas F. Martitsch (11.06.2014)
             ! Optional output (necessary for modeling the magnetic rotation)
             fieldperiod%mdata%bcovar_s_hat(i) = &
                  fieldperiod%mdata%bcovar_s_hat(i) + &
                  (bcovar_s_hat_start - bcovar_s_hat_end) * phi_mult
             fieldperiod%mdata%dlogbds(i) = &
                  fieldperiod%mdata%dlogbds(i) + &
                  (dlogbds_start - dlogbds_end) * phi_mult
             fieldperiod%mdata%dbcovar_s_hat_dphi(i) = &
                  fieldperiod%mdata%dbcovar_s_hat_dphi(i) + &
                  (dbcovar_s_hat_dphi_start - dbcovar_s_hat_dphi_end) * phi_mult
             !! End Modifications by Andreas F. Martitsch (11.06.2014)
          END DO
       END IF
    END IF

    ! make new summation
    fieldperiod => fieldline%ch_fir ! first field_period
    average_bhat  = 0.0d0
    average_bhat2 = 0.0d0
    average_one_over_bhat = 0.0d0
    count_sum_bhat = 0
    one_over_h_phi_int=0.0d0
    b_over_h_phi_int=0.0d0 
    b2_over_h_phi_int=0.0d0
    DO ! go through all periods
       uplim = UBOUND(fieldperiod%mdata%bhat,1)
       average_bhat  = average_bhat  + SUM(fieldperiod%mdata%bhat(1:uplim))
       average_bhat2 = average_bhat2 + SUM(fieldperiod%mdata%bhat(1:uplim)**2)
       average_one_over_bhat = average_one_over_bhat + SUM(1.0d0 / fieldperiod%mdata%bhat(1:uplim))
       count_sum_bhat = count_sum_bhat + uplim
       one_over_h_phi_int=one_over_h_phi_int                                &
                         +SUM(1.d0/fieldperiod%mdata%h_phi(1:uplim))
       b_over_h_phi_int  =b_over_h_phi_int                                  &
                         +SUM(fieldperiod%mdata%bhat(1:uplim)               &
                         /    fieldperiod%mdata%h_phi(1:uplim))
       b2_over_h_phi_int=b2_over_h_phi_int                                  &
                        +SUM(fieldperiod%mdata%bhat(1:uplim)**2             &
                        /    fieldperiod%mdata%h_phi(1:uplim))

       IF (.NOT.(ASSOCIATED(fieldperiod%next))) EXIT
       fieldperiod => fieldperiod%next
    END DO
    average_bhat  = average_bhat  / count_sum_bhat
    average_bhat2 = average_bhat2 / count_sum_bhat
    average_one_over_bhat = average_one_over_bhat / count_sum_bhat
    surface_boozer_B00=b2_over_h_phi_int/b_over_h_phi_int
    travis_convfac=one_over_h_phi_int/b_over_h_phi_int

    PRINT *, '  average_bhat:           ',average_bhat,' points: ',count_sum_bhat
    PRINT *, '  average_bhat2:          ',average_bhat2
    PRINT *, '  average_one_over_bhat:  ',average_one_over_bhat
    PRINT *, '  surface_boozer_B00:     ',surface_boozer_B00
    PRINT *, '  travis_convfac:         ',travis_convfac

    ! make new summation - end

    ! redo fieldline for ybeg and yend
    IF (mag_symmetric) THEN
       fieldperiod => fieldline%ch_fir
       correct_ybeg_yend: DO
          IF (fieldperiod%tag .EQ. 1) THEN
             y = fieldperiod%mdata%ybeg ! new starting point
             y(6:ndim) = 0.0_dp ! reset of quantities
             fieldperiod%mdata%ybeg = y
          ELSE
             fieldperiod%mdata%ybeg = fieldperiod%prev%mdata%yend
          END IF

          phi = fieldperiod%coords%x2(0)
          all_steps_redo: DO i = 1, nstep
             CALL rk4_kin(phi,h)
          END DO all_steps_redo
          fieldperiod%mdata%yend = y

          IF (ASSOCIATED(fieldperiod%next)) THEN
             fieldperiod => fieldperiod%next
          ELSE
             EXIT correct_ybeg_yend
          END IF
       END DO correct_ybeg_yend

    END IF
    ! end redo fieldline for ybeg and yend

    ! search for extrema
    IF (mag_magfield .NE. 0) THEN
       CALL find_extrema
    END IF

    CALL setup_fieldpropagators
    surface%b_abs_max = fieldline%b_abs_max
    surface%b_abs_min = fieldline%b_abs_min


    PRINT *, '  Absolute Maximum:        ',fieldline%b_abs_max,' Prop-Tag ',fieldline%abs_max_ptag
    PRINT *, '  Absolute Minimum:        ',fieldline%b_abs_min,' Prop-Tag ',fieldline%abs_min_ptag
    PRINT *, '  Periods:                 ',fieldline%ch_fir%tag,fieldline%ch_las%tag
    PRINT *, '  Propagators:             ',fieldline%ch_fir%ch_fir%tag,fieldline%ch_las%ch_las%tag
    ! in the homogeneous case this is not associated
    IF ( ASSOCIATED(fieldline%ch_fir%ch_ext) .AND. ASSOCIATED(fieldline%ch_las%ch_ext) ) THEN
       PRINT *, '  Propagators-Extra:       ',fieldline%ch_fir%ch_ext%tag,fieldline%ch_las%ch_ext%tag
    END IF
    PRINT *, '  Ripples:                 ',fieldline%ch_fir%ch_fir%ch_act%tag,fieldline%ch_las%ch_las%ch_act%tag
    PRINT *, '-------------------------------------------------'


    RETURN

  END SUBROUTINE make_mag_fieldline_newperiod

  ! ---------------------------------------------------------------------------
  ! joins all props in a fieldripple on a whole fieldline
  ! practical only for tokamaks
  ! 
  SUBROUTINE ripple_prop_join(fieldline)
    TYPE(fieldline_struct),         POINTER :: fieldline
 
    TYPE(fieldperiod_struct),       POINTER :: fieldperiod
    TYPE(fieldpropagator_struct),   POINTER :: fieldpropagator
    TYPE(fieldpropagator_struct),   POINTER :: fieldpropagator_add
    TYPE(fieldripple_struct),       POINTER :: fieldripple

    INTEGER :: proptag,proptag_first,proptag_last

    fieldperiod => fieldline%ch_fir
    fieldpropagator => fieldperiod%ch_fir
    fieldripple => fieldpropagator%ch_act
       
    fieldpropagator => fieldripple%pa_fir
    proptag_first = fieldpropagator%tag
    proptag_last  = fieldripple%pa_las%tag
    PRINT *, '-----------'
    PRINT *, 'fieldperiod%phi_l          ',fieldperiod%phi_l
    PRINT *, 'fieldperiod%phi_r          ',fieldperiod%phi_r
    PRINT *, 'fieldperiod%next%phi_l     ',fieldperiod%next%phi_l
    PRINT *, 'fieldperiod%next%phi_r     ',fieldperiod%next%phi_r
    PRINT *, '-----------'
    PRINT *, 'fieldripple%b_max_l        ',fieldripple%b_max_l
    PRINT *, 'fieldripple%b_max_r        ',fieldripple%b_max_r
    PRINT *, 'fieldripple%b_min          ',fieldripple%b_min
    PRINT *, '-----------'
    PRINT *, 'fieldripple%pa_fir%tag     ',fieldripple%pa_fir%tag
    PRINT *, 'fieldripple%pa_las%tag     ',fieldripple%pa_las%tag
    PRINT *, 'fieldripple%pa_fir%b_l     ',fieldripple%pa_fir%b_l
    PRINT *, 'fieldripple%pa_fir%b_r     ',fieldripple%pa_fir%b_r
    PRINT *, 'fieldripple%pa_las%b_l     ',fieldripple%pa_las%b_l
    PRINT *, 'fieldripple%pa_las%b_r     ',fieldripple%pa_las%b_r
    PRINT *, '-----------'
    PRINT *, 'fieldripple%pa_fir%phi_l   ',fieldripple%pa_fir%phi_l
    PRINT *, 'fieldripple%pa_fir%phi_r   ',fieldripple%pa_fir%phi_r
    PRINT *, 'fieldripple%pa_las%phi_l   ',fieldripple%pa_las%phi_l
    PRINT *, 'fieldripple%pa_las%phi_r   ',fieldripple%pa_las%phi_r
    PRINT *, '-----------'
    
    fieldpropagator_add => fieldpropagator
    IF (proptag_first .NE. proptag_last) THEN
       allpropsinripple: DO proptag = proptag_first,proptag_last-1
          fieldpropagator_add => fieldpropagator%next ! this one should be added to the previous one
          fieldpropagator%next => fieldpropagator%next%next
          PRINT *, 'fieldpropagator%tag        ',fieldpropagator%tag
          PRINT *, 'fieldpropagator%parent%tag ',fieldpropagator%parent%tag
          ! physical content
          fieldpropagator%phi_r = fieldpropagator_add%phi_r
          fieldpropagator%b_r = fieldpropagator_add%b_r
          IF (fieldpropagator_add%has_min .NE. 0) THEN
             fieldpropagator%has_min = fieldpropagator_add%has_min
             fieldpropagator%phi_min = fieldpropagator_add%phi_min
             fieldpropagator%b_min = fieldpropagator_add%b_min
             fieldpropagator%i_min = fieldpropagator_add%i_min ! ATTENTION NOT CORRECT
          END IF
          ! HERE
          ! destruct the propagator without destructing the ripple
          CALL set_magnetics_data(fieldpropagator_add)
          DEALLOCATE(fieldpropagator_add%coords)
          DEALLOCATE(fieldpropagator_add%mdata)
          NULLIFY(fieldpropagator_add%coords)
          NULLIFY(fieldpropagator_add%mdata)
          IF (ALLOCATED(fieldpropagator_add%phi_eta_ind)) DEALLOCATE(fieldpropagator_add%phi_eta_ind)
          NULLIFY(fieldpropagator_add%ch_act)
          NULLIFY(fieldpropagator_add%prev)
          NULLIFY(fieldpropagator_add%next)
          NULLIFY(fieldpropagator_add%parent)
          DEALLOCATE(fieldpropagator_add)
          NULLIFY(fieldpropagator_add)

       END DO allpropsinripple
    END IF

  END SUBROUTINE ripple_prop_join
  ! ---------------------------------------------------------------------------
  ! looks for all extrema in all periods and stores the information within
  ! the fieldperiod structure
  ! 
  SUBROUTINE find_extrema

    USE plagrange_mod
    USE extremum_mod

    INTEGER :: nstep
    INTEGER :: count
    INTEGER :: i,j,ii
    INTEGER :: niter = 32

    LOGICAL :: reg_ext

    REAL(kind=dp) :: phi,phi_start,phi_end
    REAL(kind=dp) :: bhat,der0,der1,der1i
    REAL(kind=dp) :: d2er0,d2er1,d2er1i
    REAL(kind=dp) :: h,hh
    REAL(kind=dp) :: h_mul = 1.0_dp
    REAL(kind=dp), DIMENSION(200)     :: phi_ext,bhat_ext,dbhat_ext,d2bhat_ext
    INTEGER,       DIMENSION(200)     :: minmax_ext
    REAL(kind=dp) :: delta_phi = 1.0d-4
    REAL(kind=dp) :: dummy,bhat_m,bhat_p
    REAL(kind=dp) :: nstep_phi
    REAL(kind=dp) :: phi_1,phi_2,phi_extrem,bhat_extrem
    REAL(kind=dp) :: dphi = 1.0d-16

    fieldperiod => fieldline%ch_fir
    nstep = fieldperiod%parent%parent%nstep

    search_all_periods: DO
       count = 0
       phi_ext = 0.0_dp
       bhat_ext = 0.0_dp
       dbhat_ext = 0.0_dp
       d2bhat_ext = 0.0_dp
       minmax_ext = 0
       phi = fieldperiod%coords%x2(0)
       nstep_phi = fieldperiod%coords%x2(0) - phi
       CALL plagrange_interp(fieldperiod,phi,nlagrange,bhat,der0,d2er0)
       h = fieldperiod%coords%x2(1) - fieldperiod%coords%x2(0)

       all_steps: DO i = 1,nstep
          phi = fieldperiod%coords%x2(i)
          CALL plagrange_interp(fieldperiod,phi,nlagrange,bhat,der1,d2er1)

          IF (der0*der1 .LE. 0.0_dp .AND. der1 .NE. 0.0_dp) THEN
             der1i = der1
             hh = h
             iteration_ext: DO j = 1,niter
                hh = hh*dsign(0.5_dp,der0*der1i)
                der0 = der1i
                phi = phi + hh
                CALL plagrange_interp(fieldperiod,phi,nlagrange,bhat,der1i,d2er1i)
             END DO iteration_ext

             IF (count .EQ. 0) THEN
                count = count + 1
             ELSE
                IF (minmax_ext(count) .LE. 1) THEN
                   count = count + 1
                ELSEIF (minmax_ext(count) .GT. 1 .AND. &
                     phi-phi_ext(count) .GT. mag_dphi_inf_min) THEN
                   count = count + 1
                END IF
             END IF
             phi_ext(count)    = phi
             bhat_ext(count)   = bhat
             dbhat_ext(count)  = der1i
             d2bhat_ext(count) = d2er1i
             IF (d2er1i .LT. 0.0d0) THEN
                minmax_ext(count) = 1 ! Maximum
             ELSE
                minmax_ext(count) = 0 ! Minimum
             END IF
             
          ELSEIF (d2er0*d2er1 .LE. 0.0_dp .AND. d2er1 .NE. 0.0_dp .AND. split_inflection_points) THEN
             d2er1i = d2er1
             hh = h
             iteration_inf: DO j = 1,niter
                hh = hh*dsign(0.5_dp,d2er0*d2er1i)
                d2er0 = d2er1i
                phi = phi + hh
                CALL plagrange_interp(fieldperiod,phi,nlagrange,bhat,der1i,d2er1i)
             END DO iteration_inf

             IF (ABS(der1i) .LT.  mag_dbhat_min) THEN
                reg_ext = .FALSE.
                IF (count .EQ. 0) THEN
                   reg_ext = .TRUE.
                ELSEIF (minmax_ext(count) .LE. 1 .AND. &
                     phi-phi_ext(count) .GT. mag_dphi_inf_min) THEN
                   reg_ext = .TRUE.
                ELSEIF (minmax_ext(count) .GT. 1 .AND. &
                     phi-phi_ext(count) .GT. mag_dphi_inf_min) THEN
                   reg_ext = .TRUE.
                END IF
                IF (reg_ext) THEN
                   count = count + 1
                   phi_ext(count)    = phi
                   bhat_ext(count)   = bhat
                   dbhat_ext(count)  = der1i
                   d2bhat_ext(count) = d2er1i
                   IF (der1i .LT. 0.0d0) THEN
                      minmax_ext(count) = 3 ! Left
                   ELSE
                      minmax_ext(count) = 4 ! Right
                   END IF

                END IF
             END IF
          END IF
          der0  = der1
          d2er0 = d2er1
       END DO all_steps
       
       IF (count .GT. 0) THEN
          ALLOCATE(fieldperiod%phi_ext(count))
          fieldperiod%phi_ext = phi_ext(1:count)
          ALLOCATE(fieldperiod%bhat_ext(count))
          fieldperiod%bhat_ext = bhat_ext(1:count)
          ALLOCATE(fieldperiod%dbp_ext(count))
          fieldperiod%dbp_ext = dbhat_ext(1:count)
          ALLOCATE(fieldperiod%d2bp_ext(count))
          fieldperiod%d2bp_ext = d2bhat_ext(1:count)
          ALLOCATE(fieldperiod%minmax(count))
          fieldperiod%minmax = minmax_ext(1:count)
          ALLOCATE(fieldperiod%width_left(count))
          fieldperiod%width_left = 0.0_dp
          ALLOCATE(fieldperiod%width_right(count))
          fieldperiod%width_right = 0.0_dp
       END IF

       IF (ASSOCIATED(fieldperiod%next)) THEN
          fieldperiod => fieldperiod%next
       ELSE
          EXIT search_all_periods
       END IF
    END DO search_all_periods

    ! find the appropriate width for exact treatment around maxima
    ! and inflection points   
    fieldperiod => fieldline%ch_fir
    nstep = fieldperiod%parent%parent%nstep
    search_all_periods_width: DO
       IF (ALLOCATED(fieldperiod%minmax)) THEN
          count = UBOUND(fieldperiod%minmax,1)
       ELSE
          count = 0
       END IF
       cycle_through_minmax: DO i = 1,count
          IF (fieldperiod%minmax(i) .NE. 0) THEN
             ! maximum or inflection point found
             CALL find_width(i,count)
          END IF
       END DO cycle_through_minmax
       IF (ASSOCIATED(fieldperiod%next)) THEN
          fieldperiod => fieldperiod%next
       ELSE
          EXIT search_all_periods_width
       END IF
    END DO search_all_periods_width
       
  END SUBROUTINE find_extrema

  SUBROUTINE find_width(start,count)
    INTEGER, INTENT(in)         :: start,count
    INTEGER :: dir,minmax
    REAL(kind=dp), DIMENSION(1:2) :: phi,bhat,dbhat,d2bhat
    REAL(kind=dp), DIMENSION(1:6) :: poly
    REAL(kind=dp) :: width
    phi(1) = fieldperiod%phi_ext(start)
    bhat(1) = fieldperiod%bhat_ext(start)
    dbhat(1) = fieldperiod%dbp_ext(start)
    d2bhat(1) = fieldperiod%d2bp_ext(start)
    DO dir = -1,1,2
       CALL find_next_extremum(start,count,dir,phi(2),bhat(2),dbhat(2),d2bhat(2),minmax)
       phi(2) = phi(2) - phi(1)
       phi(1) = 0.0_dp
       CALL find_poly(phi,bhat,dbhat,d2bhat,poly)
       CALL compute_width(phi,poly,width)
       IF (dir .EQ. -1) fieldperiod%width_left(start)  = width
       IF (dir .EQ.  1) fieldperiod%width_right(start) = width
    END DO
  END SUBROUTINE find_width

  SUBROUTINE find_next_extremum(start,count,dir,phi,bhat,dbhat,d2bhat,minmax)
    INTEGER, INTENT(in)         :: start,count,dir
    REAL(kind=dp), INTENT(out)  :: phi,bhat,dbhat,d2bhat
    INTEGER, INTENT(out)        :: minmax
    INTEGER :: i,s,c
    REAL(kind=dp) :: end_phi

    TYPE(fieldperiod_struct),     POINTER :: periodwalker
    
    periodwalker => fieldperiod
    s = start
    c = count

    IF (dir .EQ. -1) THEN ! look to the left
       end_phi = 0.0d0
       walker_left: DO
          i = s
          count_left: DO
             i = i - 1
             IF (i .LT. 1) EXIT count_left
             minmax = periodwalker%minmax(i)
             phi = periodwalker%phi_ext(i) - end_phi
             bhat = periodwalker%bhat_ext(i)
             dbhat = periodwalker%dbp_ext(i)
             d2bhat = periodwalker%d2bp_ext(i)
             EXIT walker_left
          END DO count_left
          IF (ASSOCIATED(periodwalker%prev)) THEN
             periodwalker => periodwalker%prev
          ELSE
             periodwalker => fieldline%ch_las
             end_phi = fieldline%ch_las%coords%x2(UBOUND(periodwalker%coords%x2,1))
          END IF
          IF (ALLOCATED(periodwalker%minmax)) THEN
             c = UBOUND(periodwalker%minmax,1)
          ELSE
             c = 0
          END IF
          s = c + 1
       END DO walker_left
    ELSE
       end_phi = 0.0d0
       walker_right: DO
          i = s
          count_right: DO
             i = i + 1
             IF (i .GT. c) EXIT count_right
             minmax = periodwalker%minmax(i)
             phi = periodwalker%phi_ext(i) + end_phi
             bhat = periodwalker%bhat_ext(i)
             dbhat = periodwalker%dbp_ext(i)
             d2bhat = periodwalker%d2bp_ext(i)
             EXIT walker_right
          END DO count_right
          IF (ASSOCIATED(periodwalker%next)) THEN
             periodwalker => periodwalker%next
          ELSE
             periodwalker => fieldline%ch_fir
             end_phi = fieldline%ch_las%coords%x2(UBOUND(periodwalker%coords%x2,1))
          END IF
          IF (ALLOCATED(periodwalker%minmax)) THEN
             c = UBOUND(periodwalker%minmax,1)
          ELSE
             c = 0
          END IF
          s = 0
       END DO walker_right
    END IF
  END SUBROUTINE find_next_extremum

  SUBROUTINE find_poly(x,y,yp,ypp,p)
    REAL(kind=dp), DIMENSION(1:2), INTENT(in)  :: x,y,yp,ypp
    REAL(kind=dp), DIMENSION(1:6), INTENT(out) :: p
    REAL(kind=dp) :: d
    
    d = -x(2)**5 - 10*x(1)**2*x(2)**3 + 10*x(2)**2*x(1)**3 - &
         5*x(1)**4*x(2) + 5*x(1)*x(2)**4 + x(1)**5
    p(1) = -(12*y(2) - 12*y(1) - 2*x(1)*x(2)*ypp(2) + 6*x(1)*yp(2) - &
         x(1)**2*ypp(1) + x(2)**2*ypp(2) + x(1)**2*ypp(2) + 2*x(1)*x(2)*ypp(1) - &
         6*x(2)*yp(2) - 6*x(2)*yp(1) - x(2)**2*ypp(1) + 6*x(1)*yp(1))
    p(2) = (3*x(1)**3*ypp(2) - 2*x(1)**3*ypp(1) + x(1)**2*x(2)*ypp(1) + &
         16*x(1)**2*yp(2) + 14*x(1)**2*yp(1) - 4*x(1)**2*x(2)*ypp(2) + &
         4*x(1)*x(2)**2*ypp(1) + 30*x(1)*y(2) - 30*x(1)*y(1) + 2*x(1)*x(2)*yp(1) - &
         x(1)*x(2)**2*ypp(2) - 2*x(1)*x(2)*yp(2) + 30*x(2)*y(2) + 2*x(2)**3*ypp(2) - &
         3*x(2)**3*ypp(1) - 14*x(2)**2*yp(2) - 30*x(2)*y(1) - 16*x(2)**2*yp(1))
    p(3) = -(3*x(1)**4*ypp(2) - x(1)**4*ypp(1) + 8*x(1)**3*yp(1) + &
         12*x(1)**3*yp(2) - 4*x(1)**3*x(2)*ypp(1) + 8*x(1)**2*x(2)**2*ypp(1) - &
         20*x(1)**2*y(1) + 32*x(1)**2*x(2)*yp(1) - 8*x(1)**2*x(2)**2*ypp(2) + &
         20*x(1)**2*y(2) + 28*x(1)**2*x(2)*yp(2) - 80*x(1)*x(2)*y(1) - 28*x(1)*x(2)**2*yp(1) - &
         32*x(1)*x(2)**2*yp(2) + 4*x(1)*x(2)**3*ypp(2) + 80*x(1)*x(2)*y(2) + 20*x(2)**2*y(2) + &
         x(2)**4*ypp(2) - 12*x(2)**3*yp(1) - 20*x(2)**2*y(1) - 3*x(2)**4*ypp(1) - &
         8*x(2)**3*yp(2))
    p(4) = (3*x(2)**4*x(1)*ypp(2) + x(1)**5*ypp(2) - 24*x(2)**3*x(1)*yp(2) + &
         8*x(2)**3*x(1)**2*ypp(1) - 4*x(2)**4*x(1)*ypp(1) - 36*x(2)**3*x(1)*yp(1) + &
         60*x(1)**2*x(2)*y(2) - 12*x(1)**2*x(2)**2*yp(2) - 60*x(1)**2*x(2)*y(1) + &
         12*x(1)**2*x(2)**2*yp(1) + 60*x(2)**2*x(1)*y(2) - 60*x(2)**2*x(1)*y(1) - &
         x(2)**5*ypp(1) + 4*x(1)**4*x(2)*ypp(2) - 8*x(1)**3*x(2)**2*ypp(2) - 3*x(1)**4*x(2)*ypp(1) + &
         36*x(1)**3*x(2)*yp(2) + 24*x(1)**3*x(2)*yp(1))
    p(5) = -(-3*x(1)**4*x(2)**2*ypp(1) + 10*x(1)**4*x(2)*yp(2) - x(1)**4*x(2)**2*ypp(2) + &
         2*x(1)**5*x(2)*ypp(2) - 24*x(1)**2*x(2)**3*yp(2) - 10*x(2)**4*x(1)*yp(1) - &
         2*x(2)**5*x(1)*ypp(1) + 3*x(2)**4*x(1)**2*ypp(2) + x(2)**4*x(1)**2*ypp(1) - &
         16*x(1)**2*x(2)**3*yp(1) - 4*x(1)**3*x(2)**3*ypp(2) + 4*x(1)**3*x(2)**3*ypp(1) + &
         16*x(2)**2*x(1)**3*yp(2) + 24*x(2)**2*x(1)**3*yp(1) - 2*x(1)**5*yp(2) + 2*x(2)**5*yp(1) - &
         60*x(2)**2*x(1)**2*y(1) + 60*x(2)**2*x(1)**2*y(2))
    p(6) = (-x(2)**5*x(1)**2*ypp(1) - 2*x(2)**5*y(1) + 2*x(2)**5*x(1)*yp(1) + &
         2*x(2)**4*x(1)**3*ypp(1) - 10*x(2)**4*x(1)**2*yp(1) + x(2)**4*x(1)**3*ypp(2) + &
         10*x(2)**4*x(1)*y(1) - 2*x(2)**3*x(1)**4*ypp(2) - x(2)**3*x(1)**4*ypp(1) + &
         8*x(2)**3*x(1)**3*yp(1) - 8*x(2)**3*x(1)**3*yp(2) - 20*x(2)**3*x(1)**2*y(1) + &
         10*x(2)**2*yp(2)*x(1)**4 - 10*y(2)*x(1)**4*x(2) + 20*y(2)*x(2)**2*x(1)**3 + &
         2*y(2)*x(1)**5 - 2*x(2)*yp(2)*x(1)**5 + x(2)**2*ypp(2)*x(1)**5)
    p = p / d / 2.0_dp
  END SUBROUTINE find_poly

  SUBROUTINE eval_poly(p,x,y)
    REAL(kind=dp), DIMENSION(1:6), INTENT(in) :: p
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: x
    REAL(kind=dp), DIMENSION(:), INTENT(inout) :: y
    y = p(6) + x * (p(5) + x * (p(4) + x * (p(3) + x * (p(2) + x * p(1)))))
  END SUBROUTINE eval_poly

  SUBROUTINE compute_width(phi,p,width)
    REAL(kind=dp), DIMENSION(1:2), INTENT(in)  :: phi
    REAL(kind=dp), DIMENSION(1:6), INTENT(in)  :: p
    REAL(kind=dp), INTENT(out) :: width
    
    ! test
    REAL(kind=dp), DIMENSION(100) :: x,y
    INTEGER :: i

    width = ABS(phi(2)-phi(1)) / 3.0_dp
    
    ! test
    x(1) = phi(1)
    DO i = 2,99
       x(i) = x(i-1) + (phi(2) - phi(1))/99.0_dp
    END DO
    x(100) = phi(2)
    CALL eval_poly(p,x,y)

  END SUBROUTINE compute_width

  ! ---------------------------------------------------------------------------
  ! with the information in fieldperiods now the fieldpropagators are set up
  ! 
  SUBROUTINE setup_fieldpropagators

    USE plagrange_mod
    
    INTEGER, PARAMETER :: nlagrange = 5

    INTEGER :: count,i,k
    INTEGER :: nstep
    REAL(kind=dp) :: phi_span

    LOGICAL :: last_min,last_max
    INTEGER :: last_period_count
    REAL(kind=dp) :: last_phi_min,last_bhat_min,last_d2bp_min
    REAL(kind=dp) :: last_phi_max,last_bhat_max,last_d2bp_max

    LOGICAL :: first_min,first_max
    INTEGER :: first_period_count
    REAL(kind=dp) :: first_phi_min,first_bhat_min,first_d2bp_min
    REAL(kind=dp) :: first_phi_max,first_bhat_max,first_d2bp_max

    LOGICAL :: first_period,last_period
    LOGICAL :: first_extra,last_extra
    LOGICAL :: ripple_start,ripple_end
    LOGICAL :: all_props_end
    LOGICAL :: has_min,has_max,has_inf

    INTEGER :: props_num,props_count
    INTEGER :: min_num,max_num,inf_num
    INTEGER :: max_count,min_count,inf_count
    INTEGER :: ripple_count

    INTEGER,DIMENSION(:), ALLOCATABLE :: prop_borders
    INTEGER,DIMENSION(:), ALLOCATABLE :: prop_hasmin
    INTEGER,DIMENSION(:), ALLOCATABLE :: prop_hasmax
    INTEGER,DIMENSION(:), ALLOCATABLE :: prop_hasinf

    REAL(kind=dp) :: ripple_phi_l,ripple_phi_r,ripple_phi_min,ripple_b_min
    REAL(kind=dp) :: ripple_d2bp_max_r,ripple_d2bp_max_l,ripple_d2bp_min

    REAL(kind=dp),DIMENSION(:), ALLOCATABLE :: phi_borders
    REAL(kind=dp) :: h1
    INTEGER,DIMENSION(:), ALLOCATABLE :: phi_num
    INTEGER :: phi_ub
    INTEGER :: n_addlim = 4

    REAL(kind=dp), DIMENSION(100) :: phi_inflection
    REAL(kind=dp), DIMENSION(100) :: b_inflection
    REAL(kind=dp), DIMENSION(100) :: dbdp_inflection
    INTEGER :: icount

    ! To avoid warning maby used uninitialized (should be false positive).
    min_count = -1

    ripple_start = .false.
    ripple_end = .false.

    first_phi_min = 1.234e+5
    first_bhat_min = 1.234e+5
    first_d2bp_min = 1.234e+5

    last_phi_min = 1.234e+5
    last_bhat_min = 1.234e+5

    ripple_phi_min = 1.234e+5
    ripple_b_min = 1.234e+5
    ripple_phi_l = 1.234e+5

    ripple_d2bp_min = 0.0
    ripple_d2bp_max_l = 0.0

    IF (mag_magfield .EQ. 0) THEN ! homogeneous case
       fieldline => fieldperiod%parent
       fieldperiod => fieldline%ch_fir
       periods_hom: DO ! all fieldperiods
          !propagator
          CALL construct_magnetics(fieldperiod,fieldpropagator)
          fieldpropagator%phi_l = fieldperiod%phi_l
          fieldpropagator%phi_r = fieldperiod%phi_r
          ! coordinates
          ALLOCATE(fieldpropagator%coords)
          ALLOCATE(fieldpropagator%coords%x1(LBOUND(fieldperiod%coords%x1,1):UBOUND(fieldperiod%coords%x1,1)))
          fieldpropagator%coords%x1 = fieldperiod%coords%x1
          ALLOCATE(fieldpropagator%coords%x2(LBOUND(fieldperiod%coords%x2,1):UBOUND(fieldperiod%coords%x2,1)))
          fieldpropagator%coords%x2 = fieldperiod%coords%x2
          ALLOCATE(fieldpropagator%coords%x3(LBOUND(fieldperiod%coords%x3,1):UBOUND(fieldperiod%coords%x3,1)))
          fieldpropagator%coords%x3 = fieldperiod%coords%x3
          ! mdata
          ALLOCATE(fieldpropagator%mdata)
          ALLOCATE(fieldpropagator%mdata%bhat(LBOUND(fieldperiod%mdata%bhat,1):UBOUND(fieldperiod%mdata%bhat,1)))
          fieldpropagator%mdata%bhat = fieldperiod%mdata%bhat
          ALLOCATE(fieldpropagator%mdata%geodcu(LBOUND(fieldperiod%mdata%geodcu,1):UBOUND(fieldperiod%mdata%geodcu,1)))
          fieldpropagator%mdata%geodcu = fieldperiod%mdata%geodcu
          ALLOCATE(fieldpropagator%mdata%h_phi(LBOUND(fieldperiod%mdata%h_phi,1):UBOUND(fieldperiod%mdata%h_phi,1)))
          fieldpropagator%mdata%h_phi = fieldperiod%mdata%h_phi
          ALLOCATE(fieldpropagator%mdata%dlogbdphi(LBOUND(fieldperiod%mdata%dlogbdphi,1):UBOUND(fieldperiod%mdata%dlogbdphi,1)))
          fieldpropagator%mdata%dlogbdphi = fieldperiod%mdata%dlogbdphi
          !! Modifications by Andreas F. Martitsch (11.06.2014)
          ! Optional output (necessary for modeling the magnetic rotation) ->
          ! Allocate the additional entries
          ALLOCATE( fieldpropagator%mdata%dbcovar_s_hat_dphi( &
                     LBOUND(fieldperiod%mdata%dbcovar_s_hat_dphi,1):&
                     UBOUND(fieldperiod%mdata%dbcovar_s_hat_dphi,1) &
                    ) &
                  )
          fieldpropagator%mdata%dbcovar_s_hat_dphi=fieldperiod%mdata%dbcovar_s_hat_dphi
          ALLOCATE( fieldpropagator%mdata%bcovar_s_hat( &
                     LBOUND(fieldperiod%mdata%bcovar_s_hat,1):&
                     UBOUND(fieldperiod%mdata%bcovar_s_hat,1) &
                    ) &
                  )
          fieldpropagator%mdata%bcovar_s_hat=fieldperiod%mdata%bcovar_s_hat
          ALLOCATE( fieldpropagator%mdata%dlogbds( &
                     LBOUND(fieldperiod%mdata%dlogbds,1):&
                     UBOUND(fieldperiod%mdata%dlogbds,1) &
                    ) &
                  )
          fieldpropagator%mdata%dlogbds=fieldperiod%mdata%dlogbds       
          !! End Modifications by Andreas F. Martitsch (11.06.2014)
          ALLOCATE(fieldpropagator%mdata%ybeg(LBOUND(fieldperiod%mdata%ybeg,1):UBOUND(fieldperiod%mdata%ybeg,1)))
          fieldpropagator%mdata%ybeg = fieldperiod%mdata%ybeg
          ALLOCATE(fieldpropagator%mdata%yend(LBOUND(fieldperiod%mdata%yend,1):UBOUND(fieldperiod%mdata%yend,1)))
          fieldpropagator%mdata%yend = fieldperiod%mdata%yend
          fieldpropagator%b_l = fieldpropagator%mdata%bhat(LBOUND(fieldpropagator%mdata%bhat,1))
          fieldpropagator%b_r = fieldpropagator%mdata%bhat(UBOUND(fieldpropagator%mdata%bhat,1))
          ! ripple
          CALL construct_magnetics(fieldpropagator,fieldripple)
          fieldripple%b_max_l = fieldpropagator%b_l
          fieldripple%b_max_r = fieldpropagator%b_r
          fieldripple%b_min   = fieldpropagator%b_r
          fieldripple%width   = fieldpropagator%phi_r - fieldpropagator%phi_l
          fieldripple%parent => fieldpropagator
          fieldripple%pa_fir => fieldpropagator
          fieldripple%pa_las => fieldpropagator
          fieldpropagator%ch_act => fieldripple
          fieldpropagator%ch_tag =  fieldripple%tag
          
          IF (ASSOCIATED(fieldperiod%next)) THEN
             fieldperiod => fieldperiod%next
          ELSE
             EXIT periods_hom
          END IF
       END DO periods_hom
       fieldline%abs_max_ptag = fieldpropagator%tag
       fieldline%abs_min_ptag = fieldpropagator%tag
       fieldline%b_abs_min = fieldripple%b_min
       fieldline%b_abs_max = fieldripple%b_max_l
       
       RETURN
    END IF

    icount=0
    
    nstep = fieldperiod%parent%parent%nstep
    phi_span = fieldline%ch_las%coords%x2(nstep) - fieldline%ch_fir%coords%x2(0)

    ! find last maximum at the end
    last_min = .FALSE.
    last_max = .FALSE.
    fieldperiod => fieldline%ch_las
    last_period_count = 0
    last_periods: DO
       last_period_count = last_period_count + 1
       IF (ALLOCATED(fieldperiod%minmax)) THEN
          count = UBOUND(fieldperiod%minmax,1)
          search_last_max: DO i = count,1,-1
             IF (fieldperiod%minmax(i) .EQ. 0) THEN
                last_min = .TRUE.
                last_phi_min  = fieldperiod%phi_ext(i)
                last_bhat_min = fieldperiod%bhat_ext(i)
                last_d2bp_min = fieldperiod%d2bp_ext(i)
             ELSEIF (fieldperiod%minmax(i) .EQ. 1) THEN
                last_max = .TRUE.
                last_phi_max  = fieldperiod%phi_ext(i)
                last_bhat_max = fieldperiod%bhat_ext(i)
                last_d2bp_max = fieldperiod%d2bp_ext(i)
                EXIT search_last_max
             END IF
          END DO search_last_max
       END IF
       IF (last_max) THEN
          EXIT last_periods
       ELSE
          IF (ASSOCIATED(fieldperiod%prev)) THEN
             fieldperiod => fieldperiod%prev
          ELSE
             PRINT *, 'ERROR in setup_fieldpropagators: no last maximum'
          END IF
       END IF
    END DO last_periods

    ! find first maximum at the beginning
    first_min = .FALSE.
    first_max = .FALSE.
    fieldperiod => fieldline%ch_fir
    first_period_count = 0
    first_periods: DO
       first_period_count = first_period_count + 1
       IF (ALLOCATED(fieldperiod%minmax)) THEN
          count = UBOUND(fieldperiod%minmax,1)
          search_first_max: DO i = 1,count
             IF (fieldperiod%minmax(i) .EQ. 0) THEN
                first_min = .TRUE.
                first_phi_min  = fieldperiod%phi_ext(i)
                first_bhat_min = fieldperiod%bhat_ext(i)
                first_d2bp_min = fieldperiod%d2bp_ext(i)
             ELSEIF (fieldperiod%minmax(i) .EQ. 1) THEN
                first_max = .TRUE.
                first_phi_max  = fieldperiod%phi_ext(i)
                first_bhat_max = fieldperiod%bhat_ext(i)
                first_d2bp_max = fieldperiod%d2bp_ext(i)
                EXIT search_first_max
             END IF
          END DO search_first_max
       END IF
       IF (first_max) THEN
          EXIT first_periods
       ELSE
          IF (ASSOCIATED(fieldperiod%next)) THEN
             fieldperiod => fieldperiod%next
          ELSE
             PRINT *, 'ERROR in setup_fieldpropagators: no first maximum'
          END IF
       END IF
    END DO first_periods

    ! setup the fieldpropagators and fieldripples
    fieldperiod => fieldline%ch_fir
    ripple_count = 0

    fieldline%abs_max_ptag = 0
    fieldline%abs_min_ptag = 0
    fieldline%b_abs_min = 1.d+100
    fieldline%b_abs_max = 0.0_dp

    setup_propagators: DO
       IF (ASSOCIATED(fieldperiod%prev)) THEN
          first_period = .FALSE.
       ELSE
          first_period = .TRUE.
       END IF
       IF (ASSOCIATED(fieldperiod%next)) THEN
          last_period = .FALSE.
       ELSE
          last_period = .TRUE.
       END IF

       IF (first_period) first_extra = .TRUE.
       last_extra = .FALSE.


       ! Here a new fieldperiod is analyzed
       IF (ALLOCATED(fieldperiod%minmax)) THEN
          min_num = 0
          max_num = 0
          inf_num = 0
          props_numfind: DO i = 1,UBOUND(fieldperiod%minmax,1)
             IF (fieldperiod%minmax(i) .EQ. 0) min_num = min_num + 1
             IF (fieldperiod%minmax(i) .EQ. 1) max_num = max_num + 1
             IF (fieldperiod%minmax(i) .GT. 1) inf_num = inf_num + 1
          END DO props_numfind
       ELSE
          min_num = 0
          max_num = 0
          inf_num = 0
       END IF
       props_num = max_num + inf_num + 1
       ALLOCATE(prop_borders(props_num))
       ALLOCATE(prop_hasmin(props_num))
       ALLOCATE(prop_hasmax(props_num))
       ALLOCATE(prop_hasinf(props_num))
       prop_borders = 0
       prop_hasmax  = 0
       prop_hasinf  = 0
       prop_hasmin  = 0
       IF (props_num .EQ. 1) THEN
          IF (min_num .GT. 0) prop_hasmin = 1
       ELSE
          i = 0
          DO k = 1,props_num
             i = i + 1
             IF (i .LE. UBOUND(fieldperiod%minmax,1)) THEN
                IF (fieldperiod%minmax(i) .EQ. 0) THEN
                   prop_hasmin(k) = i
                   i = i + 1
                END IF
             END IF
             IF (i .LE. UBOUND(fieldperiod%minmax,1)) THEN
                IF (fieldperiod%minmax(i) .EQ. 1) prop_hasmax(k) = i
                IF (fieldperiod%minmax(i) .GT. 1) prop_hasinf(k) = i
             END IF
          END DO
       END IF
       ! Initialization
       props_count = 0

       all_props: DO
          ! make all necessary propagators in one period
          CALL construct_magnetics(fieldperiod,fieldpropagator)
          
          IF (first_extra) THEN
             ! the first extra propagator
             fieldperiod%ch_ext => fieldpropagator ! extra propagator
             fieldpropagator%phi_l = last_phi_max - phi_span ! last maximum
             fieldpropagator%b_l   = last_bhat_max
             fieldpropagator%phi_r = fieldline%ch_fir%coords%x2(0) ! begining of first period
             fieldpropagator%b_r   = fieldline%ch_fir%mdata%bhat(0)
             ripple_d2bp_max_l = last_d2bp_max
             IF (last_min) THEN
                fieldpropagator%phi_min = last_phi_min - phi_span ! last minimum
                fieldpropagator%b_min = last_bhat_min
                ripple_d2bp_min = last_d2bp_min
                fieldpropagator%has_min = 1
             ELSE
                fieldpropagator%phi_min = fieldpropagator%phi_r ! right side
                fieldpropagator%b_min = fieldpropagator%b_r
                fieldpropagator%has_min = 0             
             END IF

             ripple_start = .TRUE.
             ripple_end   = .FALSE.
             first_extra  = .FALSE. 
             all_props_end = .FALSE.
          ELSEIF (last_extra) THEN
             ! the last extra propagator
             fieldperiod%ch_ext => fieldpropagator ! extra propagator
             fieldpropagator%phi_l = fieldline%ch_las%coords%x2(nstep) ! end of last period
             fieldpropagator%b_l   = fieldline%ch_las%mdata%bhat(nstep)
             fieldpropagator%phi_r = first_phi_max + phi_span ! first maximum shifted
             fieldpropagator%b_r   = first_bhat_max
             ripple_d2bp_max_r = first_d2bp_max
             IF (first_min) THEN
                fieldpropagator%phi_min = first_phi_min + phi_span ! first minimum shifted
                fieldpropagator%b_min = first_bhat_min
                ripple_d2bp_min = first_d2bp_min
                fieldpropagator%has_min = 1
             ELSE
                fieldpropagator%phi_min = fieldpropagator%phi_l ! left side
                fieldpropagator%b_min = fieldpropagator%b_l
                fieldpropagator%has_min = 0             
             END IF

             ripple_end = .TRUE.
             last_extra  = .FALSE. 
             all_props_end = .TRUE.
          ELSE
             ! normal propagators
             props_count = props_count + 1
             ! find relevant minimum
             IF (prop_hasmin(props_count) .EQ. 0) THEN
                has_min = .FALSE.
             ELSE
                has_min = .TRUE.
                min_count = prop_hasmin(props_count)
                ripple_d2bp_min = fieldperiod%d2bp_ext(min_count)
             END IF
             ! find relevant maximum
             IF (prop_hasmax(props_count) .EQ. 0) THEN
                has_max = .FALSE.
             ELSE
                has_max = .TRUE.
                max_count = prop_hasmax(props_count)
                ripple_d2bp_max_r = fieldperiod%d2bp_ext(max_count)
             END IF
             ! find relevant inflection
             IF (ripple_start) THEN
                phi_inflection  = 0.0_dp
                b_inflection    = 0.0_dp
                dbdp_inflection = 0.0_dp
                icount = 0
             END IF
             IF (prop_hasinf(props_count) .EQ. 0) THEN
                has_inf = .FALSE.
             ELSE
                has_inf = .TRUE.
                inf_count = prop_hasinf(props_count)
                icount = icount + 1
                phi_inflection(icount)  = fieldperiod%phi_ext(inf_count)
                b_inflection(icount)    = fieldperiod%bhat_ext(inf_count)
                dbdp_inflection(icount) = fieldperiod%dbp_ext(inf_count)
             END IF

             ! left side
             fieldpropagator%phi_l  = fieldpropagator%prev%phi_r 
             fieldpropagator%b_l    = fieldpropagator%prev%b_r 
             ! right side
             IF (has_max) THEN ! maximum
                fieldpropagator%phi_r = fieldperiod%phi_ext(max_count)
                fieldpropagator%b_r   = fieldperiod%bhat_ext(max_count)
                ripple_end = .TRUE.                
             ELSEIF (has_inf) THEN ! inflection
                fieldpropagator%phi_r = fieldperiod%phi_ext(inf_count)
                fieldpropagator%b_r   = fieldperiod%bhat_ext(inf_count)
             ELSE ! end of period
                fieldpropagator%phi_r = fieldperiod%coords%x2(nstep)
                fieldpropagator%b_r   = fieldperiod%mdata%bhat(nstep)
             END IF
             ! minimum
             IF (has_min) THEN
                fieldpropagator%phi_min = fieldperiod%phi_ext(min_count)
                fieldpropagator%b_min   = fieldperiod%bhat_ext(min_count)
                fieldpropagator%has_min = 1
             ELSE
                fieldpropagator%has_min = 0
                IF (fieldpropagator%b_l .LT. fieldpropagator%b_r) THEN
                   fieldpropagator%phi_min = fieldpropagator%phi_l
                   fieldpropagator%b_min   = fieldpropagator%b_l
                ELSE
                   fieldpropagator%phi_min = fieldpropagator%phi_r
                   fieldpropagator%b_min   = fieldpropagator%b_r
                END IF
             END IF

             IF (props_count .LT. props_num) THEN
                all_props_end = .FALSE.
                ! exit if maximum is at the end of period (otherwise zerowidth propagator)
                IF (fieldpropagator%phi_r .EQ. fieldperiod%coords%x2(nstep) ) THEN
                   all_props_end = .TRUE.
                   props_count = props_count + 1 
                END IF
             ELSE
                all_props_end = .TRUE.
             END IF
             

             IF (last_period .AND. props_count .EQ. props_num) THEN
                ! handles the last extra propagator
                last_extra = .TRUE.
                all_props_end = .FALSE.
             END IF
          END IF

          ! absolute maxima and minima
          IF (MAX(fieldpropagator%b_l,fieldpropagator%b_r) .GE. fieldline%b_abs_max) THEN
             fieldline%b_abs_max = MAX(fieldpropagator%b_l,fieldpropagator%b_r)
             fieldline%abs_max_ptag = fieldpropagator%tag
          END IF
          IF (fieldpropagator%b_min .LE. fieldline%b_abs_min) THEN
             fieldline%b_abs_min = fieldpropagator%b_min
             fieldline%abs_min_ptag = fieldpropagator%tag
          END IF

          ! ripples
          IF (ripple_start) THEN
             ripple_count = ripple_count + 1
             ripple_start = .FALSE.
             CALL construct_magnetics(fieldpropagator,fieldripple)
             fieldripple%pa_fir => fieldpropagator
             fieldripple%b_max_l = fieldpropagator%b_l
             ripple_phi_l = fieldpropagator%phi_l
          END IF
          IF (fieldpropagator%has_min .EQ. 1) THEN 
             ripple_phi_min = fieldpropagator%phi_min
             ripple_b_min = fieldpropagator%b_min
          END IF
          IF (ripple_end) THEN
             ripple_end = .FALSE.
             ripple_start = .TRUE.
             fieldripple%b_max_r = fieldpropagator%b_r
             fieldripple%b_min   = ripple_b_min
             ripple_phi_r = fieldpropagator%phi_r
             fieldripple%width = (ripple_phi_r - ripple_phi_l) * device%r0
             fieldripple%width_l = (ripple_phi_min - ripple_phi_l) * device%r0
             fieldripple%width_r = (ripple_phi_r - ripple_phi_min) * device%r0
             IF (ripple_count .EQ. 1) THEN
                fieldripple%d2bp_max_l = ripple_d2bp_max_l
             ELSE
                fieldripple%d2bp_max_l = fieldripple%prev%d2bp_max_r
             END IF
             fieldripple%d2bp_max_r = ripple_d2bp_max_r
             fieldripple%d2bp_min = ripple_d2bp_min
             ripple_d2bp_max_l = 0.0_dp
             ripple_d2bp_max_r = 0.0_dp
             ripple_d2bp_min = 0.0_dp
             IF (icount .GT. 0) THEN
                ALLOCATE(fieldripple%phi_inflection(icount))
                fieldripple%phi_inflection = phi_inflection(1:icount)
                ALLOCATE(fieldripple%b_inflection(icount))
                fieldripple%b_inflection = b_inflection(1:icount)
                ALLOCATE(fieldripple%dbdp_inflection(icount))
                fieldripple%dbdp_inflection = dbdp_inflection(1:icount)
             END IF
             IF (fieldripple%width_l .LT. 0.0_dp .OR. fieldripple%width_r .LT. 0.0_dp) THEN
                PRINT *, 'There is no minimum in ripple ',fieldripple%tag
                PRINT *, 'I better stop!'
                STOP
             END IF
          END IF
          fieldpropagator%ch_act => fieldripple
          fieldpropagator%ch_tag = fieldripple%tag
          fieldripple%pa_las => fieldpropagator

          IF (all_props_end) EXIT all_props

       END DO all_props

       ! final
       IF (ALLOCATED(prop_borders)) DEALLOCATE(prop_borders)
       IF (ALLOCATED(prop_hasmin)) DEALLOCATE(prop_hasmin)
       IF (ALLOCATED(prop_hasmax)) DEALLOCATE(prop_hasmax)
       IF (ALLOCATED(prop_hasinf)) DEALLOCATE(prop_hasinf)

       IF (last_period) THEN
          EXIT setup_propagators
       ELSE
          fieldperiod => fieldperiod%next
       END IF
    END DO setup_propagators
    ! remove the extra fieldpropagators from first and last child
    ! here is the problem if there is only one period
    fieldperiod => fieldline%ch_fir
    fieldperiod%ch_fir => fieldperiod%ch_ext%next
    fieldperiod => fieldline%ch_las
    fieldperiod%ch_las => fieldperiod%ch_ext%prev

    ! now go through all fielpropagators and fill in the
    ! physical quantities
    fieldperiod => fieldline%ch_fir
    fieldpropagator => fieldperiod%ch_ext ! %ch_fir
    all_props_fill: DO
       ALLOCATE(fieldpropagator%coords)
       ALLOCATE(fieldpropagator%mdata)
       fieldperiod => fieldpropagator%parent
       period_length = fieldperiod%phi_r - fieldperiod%phi_l
       IF (fieldpropagator%has_min .EQ. 1) THEN
          ALLOCATE(phi_borders(3))
          ALLOCATE(phi_num(2))
          phi_borders(1) = fieldpropagator%phi_l
          phi_borders(2) = fieldpropagator%phi_min
          phi_borders(3) = fieldpropagator%phi_r
          phi_num(1) = CEILING((phi_borders(2)-phi_borders(1)) / period_length * DBLE(nstep))
          phi_num(1) = MAX(n_addlim,phi_num(1)-MOD(phi_num(1),2))
          phi_num(2) = CEILING((phi_borders(3)-phi_borders(2)) / period_length * DBLE(nstep))      
          phi_num(2) = MAX(n_addlim,phi_num(2)-MOD(phi_num(2),2))
       ELSE
          ALLOCATE(phi_borders(2))
          ALLOCATE(phi_num(1))
          phi_borders(1) = fieldpropagator%phi_l
          phi_borders(2) = fieldpropagator%phi_r
          phi_num(1) = CEILING((phi_borders(2)-phi_borders(1)) / period_length * DBLE(nstep))
          phi_num(1) = MAX(n_addlim,phi_num(1)-MOD(phi_num(1),2))
       END IF

       ! this is even to have 0:SUM(phi_num) odd
       phi_ub = SUM(phi_num)
       ALLOCATE(fieldpropagator%coords%x2(0:phi_ub))
       h1 = (phi_borders(2)-phi_borders(1)) / DBLE(phi_num(1))
       DO i = 0,phi_num(1)-1
          fieldpropagator%coords%x2(i) = phi_borders(1) + DBLE(i) * h1
       END DO
       fieldpropagator%coords%x2(phi_num(1)) = phi_borders(2)
       IF (fieldpropagator%has_min .EQ. 1) THEN
          h1 = (phi_borders(3)-phi_borders(2)) / DBLE(phi_num(2))
          DO i = 1,phi_num(2)-1
             fieldpropagator%coords%x2(phi_num(1)+i) = phi_borders(2) + DBLE(i) * h1
          END DO
          fieldpropagator%coords%x2(phi_ub) = phi_borders(3)
       END IF
       fieldpropagator%i_min = phi_num(1)
       IF (fieldpropagator%has_min .EQ. 0 .AND. & 
            fieldpropagator%b_l .LT. fieldpropagator%b_r) fieldpropagator%i_min = 0
       DEALLOCATE(phi_borders)
       DEALLOCATE(phi_num)

       ALLOCATE(fieldpropagator%coords%x1(0:phi_ub))
       ALLOCATE(fieldpropagator%coords%x3(0:phi_ub))
       ALLOCATE(fieldpropagator%mdata%bhat(0:phi_ub))
       ALLOCATE(fieldpropagator%mdata%geodcu(0:phi_ub))
       ALLOCATE(fieldpropagator%mdata%h_phi(0:phi_ub))
       ALLOCATE(fieldpropagator%mdata%dlogbdphi(0:phi_ub))
       !! Modifications by Andreas F. Martitsch (11.06.2014)
       ! Optional output (necessary for modeling the magnetic rotation) ->
       ! Allocate the additional entries
       ALLOCATE(fieldpropagator%mdata%dbcovar_s_hat_dphi(0:phi_ub))
       ALLOCATE(fieldpropagator%mdata%bcovar_s_hat(0:phi_ub))
       ALLOCATE(fieldpropagator%mdata%dlogbds(0:phi_ub))      
       !! End Modifications by Andreas F. Martitsch (11.06.2014)

       DO i = 0,phi_ub
          !! Modifications by Andreas F. Martitsch (11.06.2014)
          ! Optional output (necessary for modeling the magnetic rotation)
          CALL plagrange_interp(fieldperiod,        &
               fieldpropagator%coords%x2(i),        &
               nlagrange,                           &
               fieldpropagator%coords%x1(i),        &
               fieldpropagator%coords%x3(i),        &
               fieldpropagator%mdata%bhat(i),       &
               fieldpropagator%mdata%geodcu(i),     &
               fieldpropagator%mdata%h_phi(i),      &
               fieldpropagator%mdata%dlogbdphi(i),  &
               fieldpropagator%mdata%dbcovar_s_hat_dphi(i), &
               fieldpropagator%mdata%bcovar_s_hat(i),       &
               fieldpropagator%mdata%dlogbds(i)             &
               )
          !! End Modifications by Andreas F. Martitsch (11.06.2014)
       END DO
       
       IF (ASSOCIATED(fieldpropagator%next)) THEN
          fieldpropagator => fieldpropagator%next
       ELSE
          EXIT all_props_fill
       END IF

    END DO all_props_fill

  END SUBROUTINE setup_fieldpropagators


  ! ---------------------------------------------------------------------------
  ! place eta in ripples 
  ! SUBROUTINE ripple_eta_mag(fieldline,collpar,eta_s_lim)
  SUBROUTINE ripple_eta_mag(collpar,eta_s_lim)

    USE gfactor_mod, ONLY : ienter,garr,gfactor

    REAL(kind=dp), INTENT(in) :: collpar
    REAL(kind=dp), INTENT(in) :: eta_s_lim

    TYPE(fieldripple_struct),     POINTER :: ripplewalker
    TYPE(fieldripple_struct),     POINTER :: rippleback

    REAL(kind=dp) :: b_shade,b_max_loc,b_max_locb,b_min_loc,b_min_locb,gamma1,gamma2
    REAL(kind=dp) :: eta_loc,lam_loc,col_loc
    REAL(kind=dp) :: etas2,etas
    REAL(kind=dp) :: phi_b,r0

    REAL(kind=dp), ALLOCATABLE :: eta_s_left(:),eta_s_right(:)
    REAL(kind=dp), ALLOCATABLE :: eta_x0_left(:),eta_x0_right(:)
    REAL(kind=dp), ALLOCATABLE :: eta_cl_left(:),eta_cl_right(:)
    REAL(kind=dp), ALLOCATABLE :: eta_shield_left(:)
    REAL(kind=dp), ALLOCATABLE :: eta_type_left(:)

    REAL(kind=dp) :: ripple_phi_left,ripple_phi_right
    REAL(kind=dp) :: period_phi_right

    INTEGER :: inf_start,inf_end,inf_ind
    REAL(kind=dp) :: phi_inf_loc,b_inf_loc,dbdp_inf_loc
    REAL(kind=dp) :: b_max_locb_left,b_max_locb_right,width_left,width_right

    REAL(kind=dp) :: eta,eta_b_abs_max,eta_loc_shield  
    REAL(kind=dp) :: eta_for_sigma,etas2_contrib,count_shielding

    INTEGER :: ripple_tag,walker_tag,back_tag
    INTEGER :: iwfirst,ilr,i_s
    INTEGER :: add_extra_contr
    INTEGER :: count
    INTEGER :: last_ripple_tag
    LOGICAL :: fin_ripple

    ! go to the first fieldripple
    fieldperiod     => fieldline%ch_fir
    fieldpropagator => fieldperiod%ch_fir
    fieldripple     => fieldpropagator%ch_act

    r0 = fieldline%parent%parent%r0
    eta_b_abs_max = 1.0d0 / fieldline%b_abs_max

    count = 0

    IF (mag_ripple_contribution .EQ. 1) THEN
       PRINT *, 'This is a dead branch'
       STOP
       ! now go through all fieldripples
       allripples : DO 
          ripple_tag = fieldripple%tag
          ! make proper allocation
          IF (ASSOCIATED(fieldripple%eta_x0)) THEN
             CALL delete_all(fieldripple%eta_x0)
             NULLIFY(fieldripple%eta_x0)
          END IF
          IF (ASSOCIATED(fieldripple%eta_s)) THEN
             CALL delete_all(fieldripple%eta_s)
             NULLIFY(fieldripple%eta_s)
          END IF
          count = 0
          leftright: DO ilr = 1,2 ! first left then right
             iwfirst = 1
             ripplewalker => fieldripple
             IF (ilr .EQ. 1) THEN ! left
                b_shade = ripplewalker%b_max_l
             ELSE ! right
                b_shade = ripplewalker%b_max_r
             END IF
             ! now walk away
             walk: DO
                walker_tag = ripplewalker%tag
                IF (ilr .EQ. 1) THEN ! left
                   b_max_loc = ripplewalker%b_max_l
                ELSE ! right
                   b_max_loc = ripplewalker%b_max_r
                END IF
                eta_loc = 1.0_dp /  b_max_loc
                IF (iwfirst .EQ. 1 .OR. b_max_loc .GT. b_shade) THEN ! contributes
                   IF (iwfirst .EQ. 1) THEN
                      add_extra_contr = 1
                   ELSE
                      add_extra_contr = 0
                   END IF

                   b_shade = b_max_loc
                   rippleback => ripplewalker
                   etas2 = 0.0_dp
                   walkback: DO
                      back_tag = rippleback%tag
                      ! here the contributions are computed
                      IF (ilr .EQ. 1) THEN ! left
                         b_max_locb = rippleback%b_max_l
                         phi_b = SQRT(2.0_dp * b_max_locb / ABS(rippleback%d2bp_max_l))
                      ELSE ! right
                         b_max_locb = rippleback%b_max_r
                         phi_b = SQRT(2.0_dp * b_max_locb / ABS(rippleback%d2bp_max_r))
                      END IF
                      b_min_loc = rippleback%b_min
                      eta_loc = 1.0d0 / b_max_locb ! local values

                      ! Correction by Winny because of tiny but negativ value
                      lam_loc = 1.0d0 - b_min_loc/b_max_locb
                      IF (lam_loc .GT. 0.0d0) THEN
                         lam_loc = SQRT(lam_loc)
                      ELSE
                         lam_loc = 0.0d0
                      END IF
                      col_loc = rippleback%width * collpar
                      etas2 = etas2 + eta_loc * lam_loc * col_loc

                      ! end of contributions
                      ! make the next step backwards or exit 
                      IF (back_tag .EQ. ripple_tag) THEN
                         EXIT walkback
                      ELSE
                         IF (ilr .EQ. 1) THEN ! left now back to the right
                            IF (ASSOCIATED(rippleback%next)) THEN
                               rippleback => rippleback%next
                            ELSE
                               rippleback => rippleback%parent%parent%parent%ch_fir%ch_fir%ch_act
                            END IF
                         ELSE ! right now back to the left
                            IF (ASSOCIATED(rippleback%prev)) THEN
                               rippleback => rippleback%prev
                            ELSE
                               rippleback => rippleback%parent%parent%parent%ch_las%ch_las%ch_act
                            END IF
                         END IF
                      END IF
                   END DO walkback
                   ! set the values
                   etas = SQRT(etas2)

                      count = count + 1

                      CALL set_new(fieldripple%eta_x0,1.0_dp/b_max_loc)
                      CALL set_new(fieldripple%eta_s,etas)
                      IF (add_extra_contr .EQ. 1) CALL set_new(fieldripple%eta_cl,DBLE(count))
                      IF (add_extra_contr .EQ. 1 .AND. mag_local_sigma .EQ. 1) THEN
                         etas = collpar * r0 * phi_b
                      ELSEIF (add_extra_contr .EQ. 1 .AND. mag_local_sigma .EQ. 2) THEN
                         etas = SQRT(collpar * r0 * phi_b)
                         etas = etas * (etas + rippleback%width/r0/ABS(phi_b))
                      END IF
                      IF (add_extra_contr .EQ. 1 .AND. &
                           (mag_local_sigma .EQ. 1 .OR. mag_local_sigma .EQ. 2) ) THEN
                         count = count + 1
                         CALL set_new(fieldripple%eta_x0,1.0_dp/b_max_loc)
                         CALL set_new(fieldripple%eta_s,etas)
                         CALL set_new(fieldripple%eta_cl,DBLE(count))
                      END IF
                      
                END IF
                ! next ripple when walking away or exit
                IF (mag_cycle_ripples .EQ. 0) THEN
                   IF (ilr .EQ. 1) THEN ! left
                      IF (ASSOCIATED(ripplewalker%prev)) THEN
                         ripplewalker => ripplewalker%prev
                      ELSE
                         EXIT walk
                      END IF
                   ELSE ! right
                      IF (ASSOCIATED(ripplewalker%next)) THEN
                         ripplewalker => ripplewalker%next
                      ELSE
                         EXIT walk
                      END IF
                   END IF
                ELSE
                   IF (ripplewalker%tag .EQ. ripple_tag .AND. iwfirst .NE. 1) THEN
                      EXIT walk
                   ELSE
                      IF (ilr .EQ. 1) THEN ! left
                         IF (ASSOCIATED(ripplewalker%prev)) THEN
                            ripplewalker => ripplewalker%prev
                         ELSE
                            ripplewalker => ripplewalker%parent%parent%parent%ch_las%ch_las%ch_act
                         END IF
                      ELSE ! right
                         IF (ASSOCIATED(ripplewalker%next)) THEN
                            ripplewalker => ripplewalker%next
                         ELSE
                            ripplewalker => ripplewalker%parent%parent%parent%ch_fir%ch_fir%ch_act
                         END IF
                      END IF
                   END IF
                END IF

                iwfirst = 0
             END DO walk

          END DO leftright

          ! go to the next fieldripple or exit
          IF (ASSOCIATED(fieldripple%next)) THEN
             fieldripple => fieldripple%next
          ELSE
             EXIT allripples
          END IF
       END DO allripples

    ELSEIF (mag_ripple_contribution .EQ. 2) THEN
       ! new stuff to look for shielding with respect to global maximum (bootstrap)
       ! here things are stored at the left and right maxima
       count_shielding = 0
       allripples2pre : DO
          ripple_tag = fieldripple%tag
          fieldripple%shielding_ll = .FALSE.
          fieldripple%shielding_lr = .FALSE.
          fieldripple%shielding_rl = .FALSE.
          fieldripple%shielding_rr = .FALSE.
          leftright2pre: DO ilr = 1,2 ! first left then right
             ripplewalker => fieldripple
             IF (ilr .EQ. 1) THEN ! go left
                eta = 1.0d0 / fieldripple%b_max_r
             ELSE ! go right
                eta = 1.0d0 / fieldripple%b_max_l
             END IF
             eta_for_sigma = eta
             ! now walk away
             etas2 = 0.0d0
             iwfirst = 1
             walk2pre: DO
                walker_tag = ripplewalker%tag
                IF (walker_tag .EQ. ripple_tag .AND. iwfirst .NE. 1) THEN ! I am home
                   EXIT walk2pre
                END IF
                IF (ilr .EQ. 1) THEN ! go left
                   eta_loc_shield = 1.0d0 / ripplewalker%b_max_l ! other side
                   eta_for_sigma = MIN(eta,1.0d0 / ripplewalker%b_max_r)
                ELSE ! go right
                   eta_loc_shield = 1.0d0 / ripplewalker%b_max_r ! other side
                   eta_for_sigma = MIN(eta,1.0d0 / ripplewalker%b_max_l)
                END IF

                b_min_loc = ripplewalker%b_min
                lam_loc = 1.0d0 - b_min_loc * eta_for_sigma
                IF (lam_loc .GT. 0.0d0) THEN
                   lam_loc = SQRT(lam_loc)
                ELSE
                   lam_loc = 0.0d0
                END IF
                !left half of ripplewalker
                b_max_loc = ripplewalker%b_max_l
                col_loc = 2.0_dp * ripplewalker%width_l * collpar
                gamma1 = eta_for_sigma * (b_max_loc - b_min_loc) / (1.0_dp - eta_for_sigma * b_min_loc)
                gamma1 = MIN(1.0d0,gamma1) ! safety messure
                gamma2 = (b_max_loc - b_min_loc) / b_max_loc
                etas2_contrib = col_loc * eta_for_sigma * lam_loc / b_max_loc * gfactor(gamma1,gamma2)
                etas2 = etas2 + etas2_contrib
                IF (etas2_contrib .LT. 0.0d0) PRINT *,'neg 2a ',walker_tag,ilr,gamma1,gfactor(gamma1,gamma2)
                !right half of ripplewalker
                b_max_loc = ripplewalker%b_max_r
                col_loc = 2.0_dp * ripplewalker%width_r * collpar
                gamma1 = eta_for_sigma * (b_max_loc - b_min_loc) / (1.0_dp - eta_for_sigma * b_min_loc)
                gamma1 = MIN(1.0d0,gamma1) ! safety messure
                gamma2 = (b_max_loc - b_min_loc) / b_max_loc
                etas2_contrib = col_loc * eta_for_sigma * lam_loc / b_max_loc * gfactor(gamma1,gamma2)
                etas2 = etas2 + etas2_contrib
                !two halfs finished
                etas = SQRT(etas2) ! total sigma up to this point
                iwfirst = 0
                IF (eta - sigma_shield_factor*etas .GE. eta_loc_shield) THEN ! shielded
                   IF (ilr .EQ. 1) THEN ! go left
                      fieldripple%shielding_rl = .TRUE.
                   ELSE ! go right
                      fieldripple%shielding_lr = .TRUE.
                   END IF
                   count_shielding = count_shielding + 1
                   EXIT walk2pre
                END IF
                IF (eta - sigma_shield_factor*etas .LE. eta_b_abs_max) THEN ! eta_b_abs_max reached
                   EXIT walk2pre
                END IF

                IF (ilr .EQ. 1) THEN ! left
                   IF (ASSOCIATED(ripplewalker%prev)) THEN
                      ripplewalker => ripplewalker%prev
                      ! Winny avoid double counting of last and first ripple which are identical
                      IF (.NOT. ASSOCIATED(ripplewalker%prev) .AND. &
                           ripplewalker%parent%parent%parent%ch_las%ch_las%ch_act%tag .EQ. ripple_tag) &
                           ripplewalker => ripplewalker%parent%parent%parent%ch_las%ch_las%ch_act
                   ELSE
                      ripplewalker => ripplewalker%parent%parent%parent%ch_las%ch_las%ch_act
                      ! if the last is not connected to a propagator (max at the end)
                      IF (ASSOCIATED(ripplewalker%next)) ripplewalker => ripplewalker%next
                      IF ( ripplewalker%tag .NE. ripple_tag) ripplewalker => ripplewalker%prev
                   END IF
                ELSE ! right
                   IF (ASSOCIATED(ripplewalker%next)) THEN
                      ripplewalker => ripplewalker%next
                      ! Winny avoid double counting of last and first ripple which are identical
                      IF (.NOT. ASSOCIATED(ripplewalker%next) .AND. &
                           ripplewalker%parent%parent%parent%ch_fir%ch_fir%ch_act%tag .EQ. ripple_tag) &
                           ripplewalker => ripplewalker%parent%parent%parent%ch_fir%ch_fir%ch_act
                   ELSE
                      ripplewalker => ripplewalker%parent%parent%parent%ch_fir%ch_fir%ch_act
                      ! if the first is not connected to a propagator (max at the beginning)
                      IF (ASSOCIATED(ripplewalker%prev)) ripplewalker => ripplewalker%prev
                      ! Winny avoid double counting of last and first ripple which are identical
                      IF ( ripplewalker%tag .NE. ripple_tag) ripplewalker => ripplewalker%next
                   END IF
                END IF
             END DO walk2pre
          END DO leftright2pre

          ! go to the next ripple or exit
          IF (ASSOCIATED(fieldripple%next)) THEN
             fieldripple => fieldripple%next
          ELSE
             EXIT allripples2pre
          END IF
       END DO allripples2pre

       ! now store also the shielding information about the other direction
       ! from previous or next ripple
       ! go to the first fieldripple
       fieldperiod     => fieldline%ch_fir
       fieldpropagator => fieldperiod%ch_fir
       fieldripple     => fieldpropagator%ch_act
       allripples2prefix : DO
          ripple_tag = fieldripple%tag
          IF (ASSOCIATED(fieldripple%prev)) THEN
             fieldripple%shielding_ll = fieldripple%prev%shielding_rl
          ELSE
             ripplewalker => ripplewalker%parent%parent%parent%ch_las%ch_las%ch_act
             IF (ASSOCIATED(ripplewalker%next)) ripplewalker => ripplewalker%next ! fix
             ripplewalker => ripplewalker%prev ! to fix doubling of first and last ripple
             fieldripple%shielding_ll = ripplewalker%shielding_rl
          END IF
          IF (ASSOCIATED(fieldripple%next)) THEN
             fieldripple%shielding_rr = fieldripple%next%shielding_lr
          ELSE
             ripplewalker => ripplewalker%parent%parent%parent%ch_fir%ch_fir%ch_act
             IF (ASSOCIATED(ripplewalker%prev)) ripplewalker => ripplewalker%prev ! fix
             ripplewalker => ripplewalker%next ! to fix doubling of first and last ripple
             fieldripple%shielding_rr = ripplewalker%shielding_lr
          END IF

          ! go to the next ripple or exit
          IF (ASSOCIATED(fieldripple%next)) THEN
             fieldripple => fieldripple%next
          ELSE
             EXIT allripples2prefix
          END IF
       END DO allripples2prefix
       PRINT *, 'after all ripples pre: shielding / not shielding ',count_shielding,' / ', &
            2*ripple_tag-count_shielding
       ! end - new stuff to look for shielding with respect to global maximum (bootstrap)
       ! each side now knows whether it is shielded or not

       ! now go through all fieldripples
       ! go to the first fieldripple
       fieldperiod     => fieldline%ch_fir
       fieldpropagator => fieldperiod%ch_fir
       fieldripple     => fieldpropagator%ch_act
       allripples2 : DO 
          ripple_tag = fieldripple%tag
          ! make proper allocation
          IF (ASSOCIATED(fieldripple%eta_x0)) THEN
             CALL delete_all(fieldripple%eta_x0)
             NULLIFY(fieldripple%eta_x0)
          END IF
          IF (ASSOCIATED(fieldripple%eta_s)) THEN
             CALL delete_all(fieldripple%eta_s)
             NULLIFY(fieldripple%eta_s)
          END IF
          IF (ASSOCIATED(fieldripple%eta_cl)) THEN
             CALL delete_all(fieldripple%eta_cl)
             NULLIFY(fieldripple%eta_cl)
          END IF
          IF (ASSOCIATED(fieldripple%eta_shield)) THEN
             CALL delete_all(fieldripple%eta_shield)
             NULLIFY(fieldripple%eta_shield)
          END IF
          IF (ASSOCIATED(fieldripple%eta_type)) THEN
             CALL delete_all(fieldripple%eta_type)
             NULLIFY(fieldripple%eta_type)
          END IF
          count = 0
          leftright2: DO ilr = 1,2 ! first left then right
             iwfirst = 1
             ripplewalker => fieldripple
             IF (ilr .EQ. 1) THEN ! go left
                b_shade = ripplewalker%b_max_l
             ELSE ! go right
                b_shade = ripplewalker%b_max_r
             END IF
             ! now walk away
             walk2: DO
                walker_tag = ripplewalker%tag

                IF (ilr .EQ. 1) THEN ! left
                   b_max_loc = ripplewalker%b_max_l
                ELSE ! right
                   b_max_loc = ripplewalker%b_max_r
                END IF
                eta_loc = 1.0_dp /  b_max_loc
                ! WINNY_NEU
                ! Winny shift shade a little to be able to use almost shielded maxima
                IF (iwfirst .EQ. 1 .OR. b_max_loc .GT. b_shade - 1.0d-11) THEN ! contributes
                   IF (iwfirst .EQ. 1) THEN
                      add_extra_contr = 1
                   ELSE
                      add_extra_contr = 0
                   END IF
                   b_shade = b_max_loc
                   rippleback => ripplewalker
                   etas2 = 0.0_dp
                   walkback2: DO
                      back_tag = rippleback%tag
                      ! This is just used for extra local contribution, see below
                      IF (ilr .EQ. 1) THEN ! left
                         b_max_locb = rippleback%b_max_l
                         phi_b = SQRT(2.0_dp * b_max_locb / ABS(rippleback%d2bp_max_l))
                      ELSE ! right
                         b_max_locb = rippleback%b_max_r
                         phi_b = SQRT(2.0_dp * b_max_locb / ABS(rippleback%d2bp_max_r))
                      END IF

                      ! local b_min
                      b_min_locb = rippleback%b_min
                      ! local lambda
                      lam_loc = 1.0d0 - b_min_locb * eta_loc
                      IF (lam_loc .GT. 0.0d0) THEN
                         lam_loc = SQRT(lam_loc)
                      ELSE
                         lam_loc = 0.0d0
                      END IF
                      ! always two contributions
                      ! the left half of the ripple
                      b_max_locb = rippleback%b_max_l
                      col_loc = 2.0_dp * rippleback%width_l * collpar
                      gamma1 = eta_loc * (b_max_locb - b_min_locb) / (1.0_dp - eta_loc * b_min_locb)
                      gamma1 = MIN(1.0d0,gamma1)
                      gamma2 = (b_max_locb - b_min_locb) / b_max_locb
                      etas2_contrib = col_loc * eta_loc * lam_loc / b_max_locb * gfactor(gamma1,gamma2)
                      etas2 = etas2 + etas2_contrib                      
                      IF (etas2_contrib .LT. 0.0_dp) THEN
                         PRINT *, 'negative etas2 contribution left',ripple_tag
                         STOP
                      END IF
                      ! the right half of the ripple
                      b_max_locb = rippleback%b_max_r
                      col_loc = 2.0_dp * rippleback%width_r * collpar
                      gamma1 = eta_loc * (b_max_locb - b_min_locb) / (1.0_dp - eta_loc * b_min_locb)
                      gamma1 = MIN(1.0d0,gamma1)
                      gamma2 = (b_max_locb - b_min_locb) / b_max_locb
                      etas2_contrib = col_loc * eta_loc * lam_loc / b_max_locb * gfactor(gamma1,gamma2)
                      etas2 = etas2 + etas2_contrib                      
                      IF (etas2_contrib .LT. 0.0_dp) THEN
                         PRINT *, 'negative etas2 contribution right',ripple_tag
                         STOP
                      END IF

                      IF (back_tag .EQ. ripple_tag) THEN
                         EXIT walkback2
                      ELSE
                         IF (ilr .EQ. 1) THEN ! left now back to the right
                            IF (ASSOCIATED(rippleback%next)) THEN
                               rippleback => rippleback%next
                               ! Winny avoid double counting of last and first ripple which are identical
                               IF (.NOT. ASSOCIATED(rippleback%next) .AND. &
                                    rippleback%parent%parent%parent%ch_fir%ch_fir%ch_act%tag .EQ. ripple_tag) &
                                    rippleback => rippleback%parent%parent%parent%ch_fir%ch_fir%ch_act
                            ELSE
                               rippleback => rippleback%parent%parent%parent%ch_fir%ch_fir%ch_act 
                               ! if the first is not connected to a propagator (max at the beginning)
                               IF (ASSOCIATED(rippleback%prev)) rippleback => rippleback%prev
                               ! Winny avoid double counting of last and first ripple which are identical
                               IF ( rippleback%tag .NE. ripple_tag) rippleback => rippleback%next
                            END IF
                         ELSE ! right now back to the left
                            IF (ASSOCIATED(rippleback%prev)) THEN
                               rippleback => rippleback%prev
                               ! Winny avoid double counting of last and first ripple which are identical
                               IF (.NOT. ASSOCIATED(rippleback%prev) .AND. &
                                    rippleback%parent%parent%parent%ch_las%ch_las%ch_act%tag .EQ. ripple_tag) &
                                    rippleback => rippleback%parent%parent%parent%ch_las%ch_las%ch_act
                            ELSE
                               rippleback => rippleback%parent%parent%parent%ch_las%ch_las%ch_act
                               ! if the last is not connected to a propagator (max at the end)
                               IF (ASSOCIATED(rippleback%next)) rippleback => rippleback%next
                               ! Winny avoid double counting of last and first ripple which are identical
                               IF ( rippleback%tag .NE. ripple_tag) rippleback => rippleback%prev
                            END IF
                         END IF
                      END IF
                   END DO walkback2

                   ! set the values
                   etas = SQRT(etas2)

                   count = count + 1
                   CALL set_new(fieldripple%eta_x0,1.0_dp/b_max_loc)
                   CALL set_new(fieldripple%eta_s,etas)
                   IF (add_extra_contr .EQ. 1) THEN
                      CALL set_new(fieldripple%eta_type,2.0d0) ! local level
                      CALL set_new(fieldripple%eta_cl,DBLE(count))
                      IF (ilr .EQ. 1) THEN ! go left
                         IF (ripplewalker%shielding_lr .AND. ripplewalker%shielding_ll) THEN
                            CALL set_new(fieldripple%eta_shield,2.0d0) ! shielded
                         ELSE
                            CALL set_new(fieldripple%eta_shield,0.0d0) ! not shielded
                         END IF
                      ELSE ! go right
                         IF (ripplewalker%shielding_rl .AND. ripplewalker%shielding_rr) THEN
                            CALL set_new(fieldripple%eta_shield,2.0d0) ! shielded
                         ELSE  
                            CALL set_new(fieldripple%eta_shield,0.0d0) ! not shielded
                         END IF
                      END IF
                   ELSE
                      CALL set_new(fieldripple%eta_type,0.0d0) ! normal level
                      IF (ilr .EQ. 1) THEN ! go left
                         IF (ripplewalker%shielding_lr .AND. ripplewalker%shielding_ll) THEN
                            CALL set_new(fieldripple%eta_shield,1.0d0) ! shielded
                         ELSE
                            CALL set_new(fieldripple%eta_shield,0.0d0) ! not shielded
                         END IF
                      ELSE ! go right
                         IF (ripplewalker%shielding_rl .AND. ripplewalker%shielding_rr) THEN
                            CALL set_new(fieldripple%eta_shield,1.0d0) ! shielded
                         ELSE  
                            CALL set_new(fieldripple%eta_shield,0.0d0) ! not shielded
                         END IF
                      END IF
                   END IF

                   IF (add_extra_contr .EQ. 1 .AND. mag_local_sigma .EQ. 1) THEN
                      etas = collpar * r0 * phi_b
                      etas = etas / 3.0d0 ! Sergei ?
                   ELSEIF (add_extra_contr .EQ. 1 .AND. mag_local_sigma .EQ. 2) THEN
                      etas = SQRT(collpar * r0 * phi_b)
                      etas = etas * (etas + rippleback%width/r0/ABS(phi_b))
                   END IF
                   IF (add_extra_contr .EQ. 1 .AND. &
                        (mag_local_sigma .EQ. 1 .OR. mag_local_sigma .EQ. 2) ) THEN
                      count = count + 1
                      CALL set_new(fieldripple%eta_x0,1.0_dp/b_max_loc)
                      CALL set_new(fieldripple%eta_s,etas)
                      CALL set_new(fieldripple%eta_cl,DBLE(count))
                      CALL set_new(fieldripple%eta_shield,0.0d0) ! not shielded
                      CALL set_new(fieldripple%eta_type,3.0d0)                       
                   END IF

                END IF
                ! next ripple when walking away or exit
                IF (mag_cycle_ripples .EQ. 0) THEN
                   IF (ilr .EQ. 1) THEN ! left
                      IF (ASSOCIATED(ripplewalker%prev)) THEN
                         ripplewalker => ripplewalker%prev
                      ELSE
                         EXIT walk2
                      END IF
                   ELSE ! right
                      IF (ASSOCIATED(ripplewalker%next)) THEN
                         ripplewalker => ripplewalker%next
                      ELSE
                         EXIT walk2
                      END IF
                   END IF
                ELSE
                   IF (ripplewalker%tag .EQ. ripple_tag .AND. iwfirst .NE. 1) THEN
                      EXIT walk2
                   ELSE
                      IF (ilr .EQ. 1) THEN ! left
                         IF (ASSOCIATED(ripplewalker%prev)) THEN
                            ripplewalker => ripplewalker%prev
                         ELSE
                            ripplewalker => ripplewalker%parent%parent%parent%ch_las%ch_las%ch_act
                            ! if the last is not connected to a propagator (max at the end)
                            IF (ASSOCIATED(ripplewalker%next)) ripplewalker => ripplewalker%next
                            IF (ripplewalker%tag .NE. ripple_tag) THEN
                               ripplewalker => ripplewalker%prev ! there is one extra
                            END IF
                         END IF
                      ELSE ! right
                         IF (ASSOCIATED(ripplewalker%next)) THEN
                            ripplewalker => ripplewalker%next
                         ELSE
                            ripplewalker => ripplewalker%parent%parent%parent%ch_fir%ch_fir%ch_act
                            ! if the first is not connected to a propagator (max at the beginning)
                            IF (ASSOCIATED(ripplewalker%prev)) ripplewalker => ripplewalker%prev
                         END IF
                      END IF
                   END IF
                END IF

                iwfirst = 0
             END DO walk2

          END DO leftright2

          ! Here now inflection points are handled
          if (split_inflection_points) then
             leftright2inf: DO ilr = 1,2 ! first left then right
                iwfirst = 1
                ripplewalker => fieldripple
                IF (ilr .EQ. 1) THEN ! go left
                   b_shade = ripplewalker%b_max_l
                ELSE ! go right
                   b_shade = ripplewalker%b_max_r
                END IF
                ! now walk away
                walk2inf: DO
                   walker_tag = ripplewalker%tag

                   IF ( ALLOCATED(ripplewalker%b_inflection) ) THEN
                      inf_start = LBOUND(ripplewalker%b_inflection,1)
                      inf_end   = UBOUND(ripplewalker%b_inflection,1)
                   ELSE ! nothing to do
                      inf_start = 0
                      inf_end   = -1
                   END IF

                   infloop: DO inf_ind = inf_start,inf_end
                      phi_inf_loc = ripplewalker%phi_inflection(inf_ind)
                      ! find out if inflection point is left or right of minimum
                      IF (phi_inf_loc .LT. ripplewalker%pa_fir%phi_l + ripplewalker%width_l) THEN
                         ! I am left of min
                         IF (ilr .EQ. 2) CYCLE infloop 
                      ELSE
                         ! I am right of min
                         IF (ilr .EQ. 1) CYCLE infloop
                      END IF
                      b_inf_loc = ripplewalker%b_inflection(inf_ind)
                      eta_loc = 1.0d0 / b_inf_loc
                      dbdp_inf_loc = ripplewalker%dbdp_inflection(inf_ind)

                      IF (ilr .EQ. 1) THEN ! left
                         b_max_loc = ripplewalker%b_max_l
                      ELSE ! right
                         b_max_loc = ripplewalker%b_max_r
                      END IF

                      IF (iwfirst .EQ. 1 .OR. b_inf_loc .GT. b_shade) THEN ! contributes

                         IF (iwfirst .EQ. 1) THEN
                            add_extra_contr = 1
                         ELSE
                            add_extra_contr = 0
                         END IF
                         b_shade = b_max_loc
                         rippleback => ripplewalker
                         etas2 = 0.0_dp
                         walkback2inf: DO
                            back_tag = rippleback%tag

                            ! local b_min
                            b_min_locb = rippleback%b_min
                            ! local lambda
                            lam_loc = 1.0d0 - b_min_locb * eta_loc
                            IF (lam_loc .GT. 0.0d0) THEN
                               lam_loc = SQRT(lam_loc)
                            ELSE
                               lam_loc = 0.0d0
                            END IF
                            ! always two contributions
                            IF (rippleback%tag .EQ. ripplewalker%tag) THEN
                               IF (ilr .EQ. 1) THEN
                                  b_max_locb_left  = b_inf_loc
                                  b_max_locb_right = rippleback%b_max_r
                                  width_left  = rippleback%width_l - (phi_inf_loc - rippleback%pa_fir%phi_l)
                                  width_right = rippleback%width_r
                               ELSE
                                  b_max_locb_left  = rippleback%b_max_l
                                  b_max_locb_right = b_inf_loc
                                  width_left  = rippleback%width_l
                                  width_right = rippleback%width_r - (rippleback%pa_las%phi_r-phi_inf_loc)
                               END IF
                            ELSE
                               b_max_locb_left  = rippleback%b_max_l
                               b_max_locb_right = rippleback%b_max_r
                               width_left  = rippleback%width_l
                               width_right = rippleback%width_r
                            END IF
                            ! the left side of the ripple
                            col_loc = 2.0_dp * width_left * collpar
                            gamma1 = eta_loc * (b_max_locb_left - b_min_locb) / (1.0_dp - eta_loc * b_min_locb)
                            gamma1 = MIN(1.0d0,gamma1)
                            gamma2 = (b_max_locb_left - b_min_locb) / b_max_locb_left
                            etas2_contrib = col_loc * eta_loc * lam_loc / b_max_locb_left * gfactor(gamma1,gamma2)
                            etas2 = etas2 + etas2_contrib                      
                            IF (etas2_contrib .LT. 0.0_dp) THEN
                               PRINT *, 'negative etas2 contribution left',ripple_tag,walker_tag
                               PRINT *, 'width ',width_left
                               PRINT *, 'col_loc ',col_loc,eta_loc,lam_loc
                               PRINT *, 'gamma ',gamma1,gamma2,gfactor(gamma1,gamma2)
                               STOP
                            END IF
                            ! the right half of the ripple
                            col_loc = 2.0_dp * width_right * collpar
                            gamma1 = eta_loc * (b_max_locb_right - b_min_locb) / (1.0_dp - eta_loc * b_min_locb)
                            gamma1 = MIN(1.0d0,gamma1)
                            gamma2 = (b_max_locb_right - b_min_locb) / b_max_locb_right
                            etas2_contrib = col_loc * eta_loc * lam_loc / b_max_locb_right * gfactor(gamma1,gamma2)
                            etas2 = etas2 + etas2_contrib                      
                            IF (etas2_contrib .LT. 0.0_dp) THEN
                               PRINT *, 'negative etas2 contribution right',ripple_tag,walker_tag
                               STOP
                            END IF

                            IF (back_tag .EQ. ripple_tag) THEN
                               EXIT walkback2inf
                            ELSE
                               IF (ilr .EQ. 1) THEN ! left now back to the right
                                  IF (ASSOCIATED(rippleback%next)) THEN
                                     rippleback => rippleback%next
                                     ! Winny avoid double counting of last and first ripple which are identical
                                     IF (.NOT. ASSOCIATED(rippleback%next) .AND. &
                                          rippleback%parent%parent%parent%ch_fir%ch_fir%ch_act%tag .EQ. ripple_tag) &
                                          rippleback => rippleback%parent%parent%parent%ch_fir%ch_fir%ch_act
                                  ELSE
                                     rippleback => rippleback%parent%parent%parent%ch_fir%ch_fir%ch_act 
                                     ! if the first is not connected to a propagator (max at the beginning)
                                     IF (ASSOCIATED(rippleback%prev)) rippleback => rippleback%prev
                                     ! Winny avoid double counting of last and first ripple which are identical
                                     IF ( rippleback%tag .NE. ripple_tag) rippleback => rippleback%next
                                  END IF
                               ELSE ! right now back to the left
                                  IF (ASSOCIATED(rippleback%prev)) THEN
                                     rippleback => rippleback%prev
                                     ! Winny avoid double counting of last and first ripple which are identical
                                     IF (.NOT. ASSOCIATED(rippleback%prev) .AND. &
                                          rippleback%parent%parent%parent%ch_las%ch_las%ch_act%tag .EQ. ripple_tag) &
                                          rippleback => rippleback%parent%parent%parent%ch_las%ch_las%ch_act
                                  ELSE
                                     rippleback => rippleback%parent%parent%parent%ch_las%ch_las%ch_act
                                     ! if the last is not connected to a propagator (max at the end)
                                     IF (ASSOCIATED(rippleback%next)) rippleback => rippleback%next
                                     ! Winny avoid double counting of last and first ripple which are identical
                                     IF ( rippleback%tag .NE. ripple_tag) rippleback => rippleback%prev
                                  END IF
                               END IF
                            END IF
                         END DO walkback2inf

                         ! set the values
                         etas = SQRT(etas2)
                         count = count + 1
                         CALL set_new(fieldripple%eta_x0,eta_loc)
                         CALL set_new(fieldripple%eta_s,etas)
                         IF (add_extra_contr .EQ. 1) THEN
                            CALL set_new(fieldripple%eta_type,5.0d0) ! local level
                            IF (ilr .EQ. 1) THEN ! go left
                               IF (ripplewalker%shielding_lr .AND. ripplewalker%shielding_ll) THEN
                                  CALL set_new(fieldripple%eta_shield,2.0d0) ! shielded
                               ELSE
                                  CALL set_new(fieldripple%eta_shield,0.0d0) ! not shielded
                               END IF
                            ELSE ! go right
                               IF (ripplewalker%shielding_rl .AND. ripplewalker%shielding_rr) THEN
                                  CALL set_new(fieldripple%eta_shield,2.0d0) ! shielded
                               ELSE  
                                  CALL set_new(fieldripple%eta_shield,0.0d0) ! not shielded
                               END IF
                            END IF
                         ELSE
                            CALL set_new(fieldripple%eta_type,4.0d0) ! normal level
                            IF (ilr .EQ. 1) THEN ! go left
                               IF (ripplewalker%shielding_lr .AND. ripplewalker%shielding_ll) THEN
                                  CALL set_new(fieldripple%eta_shield,1.0d0) ! shielded
                               ELSE
                                  CALL set_new(fieldripple%eta_shield,0.0d0) ! not shielded
                               END IF
                            ELSE ! go right
                               IF (ripplewalker%shielding_rl .AND. ripplewalker%shielding_rr) THEN
                                  CALL set_new(fieldripple%eta_shield,1.0d0) ! shielded
                               ELSE  
                                  CALL set_new(fieldripple%eta_shield,0.0d0) ! not shielded
                               END IF
                            END IF
                         END IF

                      END IF
                   END DO infloop
                   ! next ripple when walking away or exit
                   IF (mag_cycle_ripples .EQ. 0) THEN
                      IF (ilr .EQ. 1) THEN ! left
                         IF (ASSOCIATED(ripplewalker%prev)) THEN
                            ripplewalker => ripplewalker%prev
                         ELSE
                            EXIT walk2inf
                         END IF
                      ELSE ! right
                         IF (ASSOCIATED(ripplewalker%next)) THEN
                            ripplewalker => ripplewalker%next
                         ELSE
                            EXIT walk2inf
                         END IF
                      END IF
                   ELSE
                      IF (ripplewalker%tag .EQ. ripple_tag .AND. iwfirst .NE. 1) THEN
                         EXIT walk2inf
                      ELSE
                         IF (ilr .EQ. 1) THEN ! left
                            IF (ASSOCIATED(ripplewalker%prev)) THEN
                               ripplewalker => ripplewalker%prev
                            ELSE
                               ripplewalker => ripplewalker%parent%parent%parent%ch_las%ch_las%ch_act
                               ! if the last is not connected to a propagator (max at the end)
                               IF (ASSOCIATED(ripplewalker%next)) ripplewalker => ripplewalker%next
                               IF (ripplewalker%tag .NE. ripple_tag) THEN
                                  ripplewalker => ripplewalker%prev ! there is one extra
                               END IF
                            END IF
                         ELSE ! right
                            IF (ASSOCIATED(ripplewalker%next)) THEN
                               ripplewalker => ripplewalker%next
                            ELSE
                               ripplewalker => ripplewalker%parent%parent%parent%ch_fir%ch_fir%ch_act
                               ! if the first is not connected to a propagator (max at the beginning)
                               IF (ASSOCIATED(ripplewalker%prev)) ripplewalker => ripplewalker%prev
                            END IF
                         END IF
                      END IF
                   END IF

                   iwfirst = 0
                END DO walk2inf
             END DO leftright2inf
          end if
          ! End of inflection handling

          if (split_at_period_boundary) then
             ! Handling of intersections with period boundaries
             ripple_phi_left  = fieldripple%pa_fir%phi_l
             ripple_phi_right = fieldripple%pa_las%phi_r
             fieldperiod => fieldripple%pa_fir%parent
             DO
                period_phi_right = fieldperiod%phi_r
                IF (period_phi_right .LT. ripple_phi_right) THEN
                   etas = 1.0d-4
                   CALL set_new(fieldripple%eta_x0,1.0_dp/fieldperiod%ch_las%b_r)
                   CALL set_new(fieldripple%eta_s,etas)
                   CALL set_new(fieldripple%eta_type,6.0d0)   ! intersection level
                   CALL set_new(fieldripple%eta_shield,0.0d0) ! not shielded
                END IF
                IF (fieldperiod%ch_las%tag .GE. fieldripple%pa_las%tag) EXIT
                IF (.NOT. ASSOCIATED(fieldperiod%next)) EXIT
                fieldperiod => fieldperiod%next
             END DO
             ! End of handling of intersections with period boundaries
          end if

          ! sort the eta_x0 from small to large (together with eta_s)
          ! CALL sort(fieldripple%eta_s)
          ! go to the next fieldripple or exit
          IF (ASSOCIATED(fieldripple%next)) THEN
             fieldripple => fieldripple%next
          ELSE
             EXIT allripples2
          END IF
       END DO allripples2
       DEALLOCATE( garr )
       ienter = 1

       ! There is a Tokamak-problem with this matching of left and right - Winny
       IF (magnetic_device .NE. 0) THEN

          ! modification of local sigmas to have a match at the common boundary
          ! last ripple
          fieldperiod     => fieldline%ch_las
          fieldpropagator => fieldperiod%ch_las
          fieldripple     => fieldpropagator%ch_act
          IF (ASSOCIATED(fieldripple%next)) fieldripple => fieldripple%next
          last_ripple_tag = fieldripple%tag

          ! first ripple
          fieldperiod     => fieldline%ch_fir
          fieldpropagator => fieldperiod%ch_fir
          fieldripple     => fieldpropagator%ch_act
          fin_ripple = .FALSE.
          allripplesmod: DO
             CALL extract_array(fieldripple%eta_x0,eta_x0_left,1)   ! left ripple
             CALL extract_array(fieldripple%eta_s,eta_s_left,1)   ! left ripple
             CALL extract_array(fieldripple%eta_cl,eta_cl_left,1)
             CALL delete_all(fieldripple%eta_s)
             NULLIFY(fieldripple%eta_s)

             IF (ASSOCIATED(fieldripple%next) .AND. fieldripple%tag .LT. last_ripple_tag-1) THEN
                CALL extract_array(fieldripple%next%eta_x0,eta_x0_right,1) ! right ripple
                CALL extract_array(fieldripple%next%eta_s,eta_s_right,1) ! right ripple
                CALL extract_array(fieldripple%next%eta_cl,eta_cl_right,1)
                CALL delete_all(fieldripple%next%eta_s)
                NULLIFY(fieldripple%next%eta_s)
             ELSE ! final
                CALL extract_array(fieldperiod%ch_fir%ch_act%eta_x0,eta_x0_right,1)   ! first
                CALL extract_array(fieldperiod%ch_fir%ch_act%eta_s,eta_s_right,1)   ! first
                CALL extract_array(fieldperiod%ch_fir%ch_act%eta_cl,eta_cl_right,1)
                CALL delete_all(fieldperiod%ch_fir%ch_act%eta_s)
                NULLIFY(fieldperiod%ch_fir%ch_act%eta_s)
                CALL delete_all(fieldperiod%ch_fir%ch_act%eta_x0)
                NULLIFY(fieldperiod%ch_fir%ch_act%eta_x0)
                fin_ripple = .TRUE.
             END IF

             IF (mag_local_sigma .EQ. 0) THEN
                IF ( ALLOCATED(eta_s_left) .AND. ALLOCATED(eta_s_right) ) THEN ! WINNY-TOK
                   eta_s_left(INT(eta_cl_left(2))) = &
                        MAX(eta_s_right(INT(eta_cl_right(1))),eta_s_left(INT(eta_cl_left(2))))
                   eta_s_right(INT(eta_cl_right(1))) = &
                        MAX(eta_s_right(INT(eta_cl_right(1))),eta_s_left(INT(eta_cl_left(2))))
                END IF
             ELSE
                IF ( ALLOCATED(eta_s_left) .AND. ALLOCATED(eta_s_right) ) THEN ! WINNY-TOK
                   eta_s_left(INT(eta_cl_left(3))) = &
                        MAX(eta_s_right(INT(eta_cl_right(1))),eta_s_left(INT(eta_cl_left(3))))
                   eta_s_left(INT(eta_cl_left(4))) = &
                        MAX(eta_s_right(INT(eta_cl_right(4))),eta_s_left(INT(eta_cl_left(4))))
                   eta_s_right(INT(eta_cl_right(1))) = &
                        MAX(eta_s_right(INT(eta_cl_right(1))),eta_s_left(INT(eta_cl_left(3))))
                   eta_s_right(INT(eta_cl_right(2))) = &
                        MAX(eta_s_right(INT(eta_cl_right(2))),eta_s_left(INT(eta_cl_left(4))))
                END IF
             END IF

             IF (fin_ripple) THEN
                IF (mag_local_sigma .EQ. 0) THEN
                   eta_x0_right(INT(eta_cl_right(1))) = eta_x0_left(INT(eta_cl_left(2)))
                ELSE
                   eta_x0_right(INT(eta_cl_right(1))) = eta_x0_left(INT(eta_cl_left(3)))
                   eta_x0_right(INT(eta_cl_right(2))) = eta_x0_left(INT(eta_cl_left(4)))
                END IF
             END IF

             set_left: DO i_s = LBOUND(eta_s_left,1),UBOUND(eta_s_left,1)
                CALL set_new(fieldripple%eta_s,eta_s_left(i_s))
             END DO set_left
             IF (.NOT. fin_ripple) THEN
                set_right: DO i_s = LBOUND(eta_s_right,1),UBOUND(eta_s_right,1)
                   CALL set_new(fieldripple%next%eta_s,eta_s_right(i_s))
                END DO set_right
             ELSE
                IF ( ALLOCATED(eta_s_left) .AND. ALLOCATED(eta_s_right) ) THEN ! WINNY-TOK
                   set_final: DO i_s = LBOUND(eta_s_right,1),UBOUND(eta_s_right,1)
                      CALL set_new(fieldperiod%ch_fir%ch_act%eta_s,eta_s_right(i_s))
                      CALL set_new(fieldperiod%ch_fir%ch_act%eta_x0,eta_x0_right(i_s))
                   END DO set_final
                END IF
             END IF
             IF (fin_ripple) EXIT allripplesmod
             fieldripple => fieldripple%next
          END DO allripplesmod

          ! WINNY-TOK - still a problem
          ! make the last ripple an exact copy of the first ripple 
          ! first ripple
          fieldperiod     => fieldline%ch_fir
          fieldpropagator => fieldperiod%ch_fir
          fieldripple     => fieldpropagator%ch_act
          CALL extract_array(fieldripple%eta_x0,eta_x0_left,1)   ! left ripple
          CALL extract_array(fieldripple%eta_s,eta_s_left,1)   ! left ripple
          CALL extract_array(fieldripple%eta_cl,eta_cl_left,1)
          CALL extract_array(fieldripple%eta_shield,eta_shield_left,1)
          CALL extract_array(fieldripple%eta_type,eta_type_left,1)

          ! last ripple
          fieldperiod     => fieldline%ch_las
          fieldpropagator => fieldperiod%ch_las
          fieldripple     => fieldpropagator%ch_act
          IF (ASSOCIATED(fieldripple%next)) fieldripple => fieldripple%next
          CALL delete_all(fieldripple%eta_x0)
          NULLIFY(fieldripple%eta_x0)
          CALL delete_all(fieldripple%eta_s)
          NULLIFY(fieldripple%eta_s)
          CALL delete_all(fieldripple%eta_cl)
          NULLIFY(fieldripple%eta_cl)
          CALL delete_all(fieldripple%eta_shield)
          NULLIFY(fieldripple%eta_shield)
          CALL delete_all(fieldripple%eta_type)
          NULLIFY(fieldripple%eta_type)

          set_first_last: DO i_s = LBOUND(eta_s_left,1),UBOUND(eta_s_left,1)
             CALL set_new(fieldripple%eta_x0,eta_x0_left(i_s))
             CALL set_new(fieldripple%eta_s,eta_s_left(i_s))
             CALL set_new(fieldripple%eta_shield,eta_shield_left(i_s))
             CALL set_new(fieldripple%eta_type,eta_type_left(i_s))
          END DO set_first_last

          set_first_last_cl: DO i_s = LBOUND(eta_cl_left,1),UBOUND(eta_cl_left,1)
             CALL set_new(fieldripple%eta_cl,eta_cl_left(i_s))
          END DO set_first_last_cl

          IF (ALLOCATED(eta_s_left)) DEALLOCATE(eta_s_left)
          IF (ALLOCATED(eta_s_right)) DEALLOCATE(eta_s_right)
          IF (ALLOCATED(eta_x0_left)) DEALLOCATE(eta_x0_left)
          IF (ALLOCATED(eta_x0_right)) DEALLOCATE(eta_x0_right)
          IF (ALLOCATED(eta_cl_left)) DEALLOCATE(eta_cl_left)
          IF (ALLOCATED(eta_cl_right)) DEALLOCATE(eta_cl_right)
          IF (ALLOCATED(eta_shield_left)) DEALLOCATE(eta_shield_left)
          IF (ALLOCATED(eta_type_left)) DEALLOCATE(eta_type_left)
       END IF ! This is this Tokamak if

    END IF

    NULLIFY(fieldripple)
    NULLIFY(ripplewalker)
    NULLIFY(rippleback)
    NULLIFY(fieldpropagator)
    NULLIFY(fieldperiod)

  END SUBROUTINE ripple_eta_mag
  ! end ripple place eta
  ! ---------------------------------------------------------------------------

END MODULE mag_interface_mod
