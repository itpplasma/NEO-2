! A single-species precomputed collision operator written with z_eff = 0
! must be readable again and give the same matrices as a direct computation,
! also when it is read with z_eff /= 0.
program test_precom_collop_zeff
  use collop
  use rkstep_mod, only: lag, leg, legmax
  use collisionality_mod, only: collpar, v_min_resolution, v_max_resolution, &
    phi_x_max, isw_lorentz, isw_relativistic, lsw_multispecies, num_spec, &
    m_spec, T_spec, collop_bspline_order, collop_bspline_dist
  implicit none

  real(kind=dp), allocatable :: direct0(:,:), direct1(:,:)
  logical :: ok

  lag = 3
  leg = 3
  legmax = 5
  collop_base_prj = 11
  collop_base_exp = 11
  collop_bspline_order = 3
  collop_bspline_dist = 1.0_dp
  phi_x_max = 4.0_dp
  isw_lorentz = 0
  isw_relativistic = 0
  lsw_multispecies = .false.
  num_spec = 1
  allocate(m_spec(0:0), T_spec(0:0))
  m_spec = 9.1093836d-28
  T_spec = 1.6022d-9
  collpar = 1.0d-3
  v_min_resolution = 0.1_dp
  v_max_resolution = 3.0_dp

  z_eff = 1.0_dp
  lsw_read_precom = .false.
  lsw_write_precom = .false.
  call collop_load()
  direct1 = anumm_a(:,:,0)
  call collop_unload()

  z_eff = 0.0_dp
  lsw_write_precom = .true.
  call collop_load()
  direct0 = anumm_a(:,:,0)
  call collop_unload()

  ok = .true.
  lsw_write_precom = .false.
  lsw_read_precom = .true.
  call collop_load()
  ok = ok .and. same(anumm_a(:,:,0), direct0, 'read with z_eff = 0')
  call collop_unload()

  z_eff = 1.0_dp
  call collop_load()
  ok = ok .and. same(anumm_a(:,:,0), direct1, 'read with z_eff = 1')
  call collop_unload()

  if (ok) then
    print *, 'All tests passed!'
  else
    print *, 'FAIL'
    error stop 1
  end if

contains

  logical function same(a, b, label)
    real(kind=dp), intent(in) :: a(:,:), b(:,:)
    character(len=*), intent(in) :: label

    same = maxval(abs(a - b)) <= 1d-12 * maxval(abs(b))
    if (.not. same) print *, 'FAIL: ', label, maxval(abs(a - b))
  end function same

end program test_precom_collop_zeff
