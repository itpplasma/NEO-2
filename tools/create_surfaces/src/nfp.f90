program nfp1
  use neo_sub_mod,     only: neo_read_control, neo_filenames, neo_read, neo_dealloc
  use neo_input,       only: es, nfp
  use neo_exchange,    only: iota
  use nrtype
  use neo_spline_data,                        &
       only:  a_iota, b_iota, c_iota, d_iota, &
       sp_index
  use inter_interfaces, only: splinecof3, splint_horner3, &
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
  double precision, parameter :: s_corr_shift = 0.00001_dp
  integer            :: i_counter_max = 100

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

  !*****************************************
  ! Constants for collisionality calculation
  !*****************************************
  double precision,parameter  :: c        = 2.9979d10
  double precision,parameter  :: e_charge = 4.8032d-10
  double precision,parameter  :: e_mass   = 9.1094d-28
  double precision,parameter  :: p_mass   = 1.6726d-24
  double precision,parameter  :: ev       = 1.6022d-12

  splinecof_compatibility = .true.

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
  write (*,*) "Reading s,  iota and nfp from boozer file."
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
  call splinecof3(es, iota, c1, cn, lambda, sp_index, sw1, sw2, &
       a_iota, b_iota, c_iota, d_iota, m0, tf)

  i_counter_max = floor((s_end - s_beg)/(s_steps-1) / s_corr_shift)
  write (*,*) "Maximum number of tries per surfaces: ", i_counter_max

  !*************************************
  !Loop over desired number of surfaces
  !*************************************
  write (*,*) ""
  write (*,*) "                          s                        iota         periods"
  s_step = (s_end - s_beg) / (s_steps-1)
  create_surfs: do l = 1, s_steps

    s = s_beg + s_step*(l-1)
    s_corr = s_corr_shift

    swd = 0 ! no derivative
    call splint_horner3(es,                           &
           & a_iota, b_iota, c_iota, d_iota, swd, m0, &
           & s, tfone, tfzero, tfzero, tfzero,        &
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

    write (234,*) s, s_iota
    write (235,*) s, i_period

  end do create_surfs

  deallocate(a_iota, b_iota, c_iota, d_iota)

  !call neo_dealloc
end program nfp1
