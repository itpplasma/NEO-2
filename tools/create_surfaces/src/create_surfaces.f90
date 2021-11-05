!> \brief Program to create input profile from given profile file.
!>
!> This program creates a profile file from a given input profile, for
!> a number of radial points in a range (boozer s).
!>
!> Radial points are created in two steps. First an equidistant grid is
!> created with the appropriate bounds. In a second step it is made sure
!> that none of the surfaces corresponds to a low-number rational
!> surface, by shifting those surfaces slightly where it is the case.
!> If shifting the surface does not work, the the surface is skipped.
program create_surfaces
  use neo_sub_mod,     only : neo_read_control, neo_filenames, neo_read, neo_dealloc
  use neo_input,       only : es, nfp
  use neo_exchange,    only : iota
  use nrtype
  use neo_spline_data, only :  a_iota, b_iota, c_iota, d_iota, sp_index
  use inter_interfaces, only : splinecof3, splint_horner3, &
       tf, tfzero, tfone

  implicit none

  !*************************************
  ! Estimation of number of periods
  !*************************************
  double precision   :: s
  double precision   :: s_step
  integer            :: ns, k, l
  integer            :: i_period
  double precision   :: dist_min_req
  double precision   :: dist, theta_start
  double precision   :: theta_e
  double precision   :: s_corr
  double precision, parameter :: s_corr_shift = 0.001_dp
  integer            :: i_counter     = 0
  integer            :: i_counter_max = 100
  integer            :: surf_counter  = 0
  double precision, dimension(:), allocatable :: s_surf
  double precision, dimension(:), allocatable :: n_surf
  double precision, dimension(:), allocatable :: T_surf
  double precision, dimension(:), allocatable :: Z_surf

  !*************************************
  ! Spline of iota
  !*************************************
  integer            :: sw1, sw2, swd
  double precision   :: m0, c1, cn, s_iota
  real(dp), dimension(:), allocatable :: lambda
  real(dp)                            :: yp, ypp, yppp

  !*************************************
  ! Reading of namelist
  !*************************************
  integer            :: u1, ios

  !*************************************
  ! Settings from namelist
  !*************************************
  double precision   :: s_beg                = 0.10d0
  double precision   :: s_end                = 0.99d0
  integer            :: s_steps              = 10
  integer            :: mag_nperiod_min      = 200
  integer            :: mag_nperiod_max      = 300
  character(len=500) :: output_file          = 'surfaces.dat'
  logical            :: isw_create_surfaces  = .false.
  character(len=500) :: profiles_file        = 'profiles.in'
  namelist /settings/                                                        &
       s_beg, s_end, s_steps, mag_nperiod_min, mag_nperiod_max, output_file, &
       isw_create_surfaces, profiles_file

  !*************************************
  ! Profile
  !*************************************
  integer            :: n_p
  double precision, dimension(:), allocatable :: s_prof, n_prof, T_prof, Z_prof
  double precision, dimension(:), allocatable :: a_n_prof, b_n_prof, c_n_prof, d_n_prof
  double precision, dimension(:), allocatable :: a_T_prof, b_T_prof, c_T_prof, d_T_prof
  double precision, dimension(:), allocatable :: a_Z_prof, b_Z_prof, c_Z_prof, d_Z_prof

  !*****************************************
  ! Constants for collisionality calculation
  !*****************************************
  double precision,parameter  :: c        = 2.9979d10
  double precision,parameter  :: e_charge = 4.8032d-10
  double precision,parameter  :: e_mass   = 9.1094d-28
  double precision,parameter  :: p_mass   = 1.6726d-24
  double precision,parameter  :: ev       = 1.6022d-12
  double precision :: collog, v_te, tau_ee, collpar

  integer :: mag_nperiod_max_old


  !*************************************
  ! Read input arguments
  !*************************************
  u1 = 100
  open(unit=u1,file='create_surfaces.in',status='old',iostat=ios)
  if (ios .ne. 0) then
    print *, 'WARNING: Settings file cannot be opened, using defaults!'
    print *, ''
  else
    read(u1,nml=settings,iostat=ios)
    if (ios .ne. 0) then
      print *, 'WARNING: group settings cannot be READ!'
      print *, ''
    end if
  end if

  !*************************************
  ! Extract s and iota from boozer file
  !*************************************
  write (*,*) "Reading s, iota and nfp from boozer file."
  call neo_read_control
  call neo_filenames
  call neo_read
  write (*,*) "Number of field periods: ", nfp

  !*************************************
  !Preparing spline for iota
  !*************************************
  ns = ubound(es,1) - lbound(es,1) + 1
  sw1 = 2
  sw2 = 4
  m0  = 0.0_dp
  c1  = 0.0_dp
  cn  = 0.0_dp
  allocate (lambda(ns))
  lambda = 1.0D0
  allocate (a_iota(ns), b_iota(ns))
  allocate (c_iota(ns), d_iota(ns))
  allocate (sp_index(ns))
  sp_index = (/ (k, k=1,ns) /)

  splinecof_compatibility = .false.
  call splinecof3(es, iota, c1, cn, lambda, sp_index, sw1, sw2, &
       a_iota, b_iota, c_iota, d_iota, m0, tf)

  !*************************************
  ! Allocate output variables
  !*************************************
  allocate(s_surf(s_steps), n_surf(s_steps))
  allocate(T_surf(s_steps), Z_surf(s_steps))

  i_counter_max = floor(((s_end - s_beg)/(s_steps-1) / (4*s_corr_shift)))
  write (*,*) "Maximum number of tries per surfaces: ", i_counter_max

  !*************************************
  ! Loop over desired number of surfaces
  !*************************************
  write (*,*) ""
  write (*,*) "                          s                        iota         periods"
  s_step = (s_end - s_beg) / (s_steps-1)

  mag_nperiod_max_old = mag_nperiod_max

  i_counter_max = s_step / (2*s_corr_shift)
  write (*,*) i_counter_max

  create_surfs: do l = 1, s_steps

    s = s_beg + s_step*(l-1)
    s_corr = s_corr_shift
    mag_nperiod_max = mag_nperiod_max_old

    correct_s: do

      swd = 0 ! no derivative
      call splint_horner3(es,                                   &
         & a_iota, b_iota, c_iota, d_iota, swd, m0,             &
         & s, tfone, tfzero, tfzero, tfzero,                    &
         & s_iota, yp, ypp, yppp)

      !**********************************************************
      ! Estimate number of periods
      ! This is the same part as in NEO-2
      !**********************************************************
      dist_min_req = 0.1_dp
      dist = dist_min_req
      theta_start  = 3.1415926d0
      theta_e      = theta_start
      i_period = 0
      do
        i_period = i_period + 1

        ! Changed NEO-2 behavior for theta_e, but with same results
        theta_e = theta_e + (2.0_dp*pi)/nfp * s_iota

        if (theta_e .ge. 2.0_dp*pi) theta_e = theta_e - 2.0_dp*pi
        dist = abs(theta_e - theta_start)
        if (dist .ge. 2.0_dp*pi) dist = dist - 2.0_dp*pi
        dist = min(dist,2.0_dp*pi-dist)
        if (i_period .le. mag_nperiod_min) then
          dist_min_req = min(dist_min_req,dist)
        else
          if (dist .lt. dist_min_req) then
            exit
          end if
        end if
      end do

      !write (*,*) s, s_iota, i_period

      !******************************************************
      ! Check if condition for number of periods is fulfilled
      !******************************************************
      if (i_period > mag_nperiod_max) then
        !write (*,*) "** Warning, too many periods, trying again with new s **"

        !****************************************************************
        ! Emergency exit, if no suitable s is found in the current region
        !****************************************************************
        i_counter = i_counter + 1
        if (i_counter > i_counter_max) then
          i_counter = 0
          write (*,*) "Resonance at s=", (s_beg + s_step*(l-1))

          mag_nperiod_max = mag_nperiod_max + 200
          write (*,*) "Increasing mag_nperiod_max to ", mag_nperiod_max
          s = s_beg + s_step*(l-1)
          s_corr = s_corr_shift

          cycle correct_s

          !cycle create_surfs
        end if
        !s = s + (s_beg + s_step*(l-1)) * s_corr
        s = s + s_corr

        if (s > s_end) then
          s = s_beg + s_step*(l-1)
          s_corr = -s_corr*0.1_dp
        end if
      else
        !***************************************************
        ! Found good s, write results
        !***************************************************
        write (*,'(E28.17, E28.17, I16)') s, s_iota, i_period
        write (123,'(E28.17, I16)') s, i_period

        i_counter = 0
        surf_counter = surf_counter + 1
        s_surf(surf_counter) = s
        exit correct_s
      end if

    !write (*,*) s, s_corr, s_beg + s_step*(l-1)

    end do correct_s
  end do create_surfs

  deallocate(a_iota, b_iota, c_iota, d_iota)

  !*****************************************
  ! Close results file
  !*****************************************
  write (*,*) ""
  write (*,*) "Found   ", surf_counter, " suitable surfaces."
  write (*,*) "Skipped ", s_steps - surf_counter, " surfaces."

  !********************************************
  ! Read profile and interpolate n, T and Z_eff
  !********************************************
  if (isw_create_surfaces) then
    open(unit=u1,file=profiles_file,status='old',iostat=ios)
    read(u1, *) n_p
    allocate(s_prof(n_p), n_prof(n_p))
    allocate(T_prof(n_p), Z_prof(n_p))
    do k = 1, n_p
      read(u1, *) s_prof(k), n_prof(k), T_prof(k), Z_prof(k)
    end do
    close(u1)

    !*************************************
    ! Open results file
    !*************************************
    open(unit=u1, file=output_file, status='replace', iostat=ios)
    write (*,*) "s, kappa, Z_eff, T (eV)"

    !**********************************************
    ! Prepare for spline
    !**********************************************
    if (allocated(lambda))   deallocate(lambda)
    if (allocated(sp_index)) deallocate(sp_index)
    allocate (lambda(n_p))
    lambda = 1.0d0
    allocate (sp_index(n_p))
    sp_index = (/ (k, k=1,n_p) /)

    splinecof_compatibility = .false.
    allocate(a_n_prof(n_p), b_n_prof(n_p), c_n_prof(n_p), d_n_prof(n_p))
    call splinecof3(s_prof, n_prof, c1, cn, lambda, sp_index, sw1, sw2, &
       & a_n_prof, b_n_prof, c_n_prof, d_n_prof, m0, tf)

    allocate(a_T_prof(n_p), b_T_prof(n_p), c_T_prof(n_p), d_T_prof(n_p))
    call splinecof3(s_prof, T_prof, c1, cn, lambda, sp_index, sw1, sw2, &
       & a_T_prof, b_T_prof, c_T_prof, d_T_prof, m0, tf)

    allocate(a_Z_prof(n_p), b_Z_prof(n_p), c_Z_prof(n_p), d_Z_prof(n_p))
    call splinecof3(s_prof, Z_prof, c1, cn, lambda, sp_index, sw1, sw2, &
       & a_Z_prof, b_Z_prof, c_Z_prof, d_Z_prof, m0, tf)

    !**********************************************************
    ! Loop over the best surfaces to compute the collisionality
    !**********************************************************
    do k = 1, surf_counter

      swd = 0 ! no derivative
      call splint_horner3(s_prof,                               &
         & a_n_prof, b_n_prof, c_n_prof, d_n_prof, swd, m0,     &
         & s_surf(k), tfone, tfzero, tfzero, tfzero,            &
         & n_surf(k), yp, ypp, yppp)

      call splint_horner3(s_prof,                               &
         & a_T_prof, b_T_prof, c_T_prof, d_T_prof, swd, m0,     &
         & s_surf(k), tfone, tfzero, tfzero, tfzero,            &
         & T_surf(k), yp, ypp, yppp)

      call splint_horner3(s_prof,                               &
         & a_Z_prof, b_Z_prof, c_Z_prof, d_Z_prof, swd, m0,     &
         & s_surf(k), tfone, tfzero, tfzero, tfzero,            &
         & Z_surf(k), yp, ypp, yppp)

      !*************************************************************
      ! Computation of the collisionality parameter
      ! from /proj/plasma/Neo2/HGW_2011/TEST
      !*************************************************************
      collog=39.1d0-1.15d0*log10(n_surf(k)*1.d6)+2.3d0*log10(T_surf(k)*1.d-3)
      v_te=sqrt(2.d0*T_surf(k)*ev/e_mass)
      tau_ee=3.d0*e_mass**2*v_te**3/(16.d0*sqrt(pi)*n_surf(k)*e_charge**4*collog)

      !**********************************************************
      ! Change by Gernot Kapper - 11.02.2016
      ! Internal consistency check of collisionality shows that in NEO-2
      ! $\kappa = 2/l_c$, where $l_c = 1/(v_{te} T_{ee})$
      !**********************************************************
      ! collpar=-4.d0/(v_te*tau_ee)    ! Before patch
      collpar=-2.d0/(v_te*tau_ee)      ! After patch

      write (*,'(F10.5, E20.8, E20.8, E20.8, ES20.8)') s_surf(k), n_surf(k),T_surf(k), Z_surf(k)
      !write (*,*) collog, v_te, tau_ee, collpar

      !**********************************************
      ! Write to output file
      !**********************************************
      write (u1, *) s_surf(k), collpar, Z_surf(k), T_surf(k)

    end do

    !*************************************************
    ! Close results file and dellocate space
    !*************************************************
    close(u1)
    deallocate(s_prof, T_prof, n_prof, Z_prof)
    deallocate(a_n_prof, b_n_prof, c_n_prof, d_n_prof)
    deallocate(a_T_prof, b_T_prof, c_T_prof, d_T_prof)
    deallocate(a_Z_prof, b_Z_prof, c_Z_prof, d_Z_prof)

  end if

  !****************************************************
  ! Deallocate memory
  !****************************************************
  if (allocated(lambda))   deallocate(lambda)
  if (allocated(sp_index)) deallocate(sp_index)
  deallocate(s_surf, n_surf, T_surf, Z_surf)

end program create_surfaces
