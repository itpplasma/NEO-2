program test_magnetics_pitch_support
  use hdf5, only: HID_T
  use hdf5_tools, only: h5_init, h5_open, h5_open_group, h5_close_group, h5_close, h5_get
  use magnetics_mod, only: device_struct, surface_struct, fieldline_struct, &
       fieldperiod_struct, fieldpropagator_struct, fieldripple_struct, &
       coordinates_struct, magneticdata_struct, construct_magnetics, &
       h5_magnetics, mag_talk, mag_write_hdf5, h5_magnetics_file_name
  implicit none

  integer, parameter :: dp = kind(1.0d0)
  real(kind=dp), parameter :: pi = acos(-1.0_dp)
  real(kind=dp), parameter :: b_ref = 5.0_dp
  real(kind=dp), parameter :: ripple_depth = 0.2_dp
  real(kind=dp), parameter :: tol = 100.0_dp * epsilon(1.0_dp)

  type(device_struct), pointer :: device
  type(surface_struct), pointer :: surface
  type(fieldline_struct), pointer :: fieldline
  type(fieldperiod_struct), pointer :: fieldperiod
  type(fieldpropagator_struct), pointer :: fieldpropagator
  type(fieldripple_struct), pointer :: fieldripple
  integer(HID_T) :: file_id, category_id, group_id
  integer :: status, owner_tag, ios
  logical :: output_exists
  real(kind=dp) :: stored_b_ref, stored_b_min, stored_eta_left
  real(kind=dp) :: stored_bhat(3), stored_eta(3), expected_b(3)

  status = 0
  h5_magnetics_file_name = 'test_pitch_support_magnetics.h5'
  open(unit=99, file=h5_magnetics_file_name, status='old', iostat=ios)
  if (ios == 0) close(unit=99, status='delete')

  call h5_init()
  mag_talk = .false.
  mag_write_hdf5 = .false.
  call make_analytic_ripple()
  call h5_magnetics(device)
  inquire (file=h5_magnetics_file_name, exist=output_exists)
  if (output_exists) then
     print *, 'FAIL: default-off producer created an HDF5 file'
     status = status + 1
     open(unit=99, file=h5_magnetics_file_name, status='old', iostat=ios)
     if (ios == 0) close(unit=99, status='delete')
  end if

  mag_write_hdf5 = .true.
  call h5_magnetics(device)

  call h5_open(h5_magnetics_file_name, file_id)
  call h5_open_group(file_id, 'surface', category_id)
  call h5_open_group(category_id, '1', group_id)
  call h5_get(group_id, 'bmod0', stored_b_ref)
  call h5_close_group(group_id)
  call h5_close_group(category_id)

  call h5_open_group(file_id, 'fieldpropagator', category_id)
  call h5_open_group(category_id, '1', group_id)
  call h5_get(group_id, 'bhat', stored_bhat)
  call h5_close_group(group_id)
  call h5_close_group(category_id)

  call h5_open_group(file_id, 'fieldripple', category_id)
  call h5_open_group(category_id, '1', group_id)
  call h5_get(group_id, 'b_min', stored_b_min)
  call h5_get(group_id, 'eta_boundary_left', stored_eta_left)
  call h5_get(group_id, 'eta', stored_eta)
  call h5_get(group_id, 'pa_fir', owner_tag)
  call h5_close_group(group_id)
  call h5_close_group(category_id)
  call h5_close(file_id)

  ! Independent analytic oracle: B(phi)=B_ref*(1-ripple_depth*cos(phi)).
  expected_b = b_ref * [1.0_dp + ripple_depth, 1.0_dp - ripple_depth, &
       1.0_dp + ripple_depth]
  call assert_close(stored_b_ref * stored_bhat, expected_b, 'physical B samples', status)
  call assert_scalar(stored_b_min, 1.0_dp - ripple_depth, 'well minimum', status)
  call assert_scalar(stored_eta_left, 1.0_dp / (1.0_dp + ripple_depth), &
       'left trapped/passing boundary', status)
  call assert_scalar(stored_eta(2), stored_eta_left, 'critical pitch node', status)
  if (owner_tag /= fieldpropagator%tag) then
     print *, 'FAIL: left endpoint owner tag', owner_tag, fieldpropagator%tag
     status = status + 1
  end if

  open(unit=99, file=h5_magnetics_file_name, status='old', iostat=ios)
  if (ios == 0) close(unit=99, status='delete')

  if (status == 0) then
     print *, 'All tests passed!'
  else
     print *, 'FAIL: magnetics pitch-support tests:', status
     error stop
  end if

contains

  subroutine make_analytic_ripple()
    real(kind=dp) :: phi(0:2), bhat(0:2), zeros(0:2), state(1:2)

    call construct_magnetics(device)
    call construct_magnetics(device, surface)
    call construct_magnetics(surface, fieldline)
    call construct_magnetics(fieldline, fieldperiod)
    call construct_magnetics(fieldperiod, fieldpropagator)
    call construct_magnetics(fieldpropagator, fieldripple)

    device%name = 'analytic ripple'
    device%r0 = 1.0_dp
    device%nfp = 1
    surface%bmod0 = b_ref
    surface%b_abs_min = 1.0_dp - ripple_depth
    surface%b_abs_max = 1.0_dp + ripple_depth
    surface%nperiod = 1
    surface%nstep = 2
    surface%ndim = 2
    fieldline%b_abs_min = surface%b_abs_min
    fieldline%b_abs_max = surface%b_abs_max
    fieldline%abs_min_ptag = fieldpropagator%tag
    fieldline%abs_max_ptag = fieldpropagator%tag

    phi = [-pi, 0.0_dp, pi]
    bhat = [1.0_dp + ripple_depth, 1.0_dp - ripple_depth, &
         1.0_dp + ripple_depth]
    zeros = 0.0_dp
    state = 0.0_dp

    fieldperiod%phi_l = phi(0)
    fieldperiod%phi_r = phi(2)
    allocate(fieldperiod%phi_ext(1), fieldperiod%bhat_ext(1), &
         fieldperiod%dbp_ext(1), fieldperiod%d2bp_ext(1), &
         fieldperiod%minmax(1), fieldperiod%width_left(1), &
         fieldperiod%width_right(1))
    fieldperiod%phi_ext = 0.0_dp
    fieldperiod%bhat_ext = 1.0_dp - ripple_depth
    fieldperiod%dbp_ext = 0.0_dp
    fieldperiod%d2bp_ext = ripple_depth
    fieldperiod%minmax = -1
    fieldperiod%width_left = pi
    fieldperiod%width_right = pi
    allocate(fieldperiod%coords, fieldperiod%mdata)
    call allocate_sample_data(fieldperiod%coords, fieldperiod%mdata, phi, bhat, zeros, state)

    fieldpropagator%phi_l = phi(0)
    fieldpropagator%phi_r = phi(2)
    fieldpropagator%b_l = 1.0_dp + ripple_depth
    fieldpropagator%b_r = 1.0_dp + ripple_depth
    fieldpropagator%phi_min = 0.0_dp
    fieldpropagator%b_min = 1.0_dp - ripple_depth
    fieldpropagator%i_min = 1
    fieldpropagator%has_min = 1
    allocate(fieldpropagator%phi_eta_ind(0:2, 2))
    fieldpropagator%phi_eta_ind = 0
    allocate(fieldpropagator%coords, fieldpropagator%mdata)
    call allocate_sample_data(fieldpropagator%coords, fieldpropagator%mdata, &
         phi, bhat, zeros, state)

    fieldripple%pa_fir => fieldpropagator
    fieldripple%pa_las => fieldpropagator
    fieldripple%b_max_l = 1.0_dp + ripple_depth
    fieldripple%b_max_r = 1.0_dp + ripple_depth
    fieldripple%b_min = 1.0_dp - ripple_depth
    fieldripple%width = 2.0_dp * pi
    fieldripple%width_l = pi
    fieldripple%width_r = pi
    allocate(fieldripple%phi_inflection(1), fieldripple%b_inflection(1), &
         fieldripple%dbdp_inflection(1), fieldripple%eta(3), &
         fieldripple%eta_loc(3))
    fieldripple%phi_inflection = 0.5_dp * pi
    fieldripple%b_inflection = 1.0_dp
    fieldripple%dbdp_inflection = ripple_depth
    fieldripple%eta = [0.0_dp, 1.0_dp / (1.0_dp + ripple_depth), &
         1.0_dp / (1.0_dp - ripple_depth)]
    fieldripple%eta_loc = fieldripple%eta
    fieldripple%eta_boundary_left = fieldripple%eta(2)
    fieldripple%eta_boundary_right = fieldripple%eta(2)
    allocate(fieldripple%eta_x0, fieldripple%eta_s, fieldripple%eta_cl, &
         fieldripple%eta_shield, fieldripple%eta_type)
    fieldripple%eta_x0%d = fieldripple%eta(2)
    fieldripple%eta_s%d = fieldripple%eta(2)
    fieldripple%eta_cl%d = fieldripple%eta(2)
    fieldripple%eta_shield%d = fieldripple%eta(2)
    fieldripple%eta_type%d = 1.0_dp
  end subroutine make_analytic_ripple

  subroutine allocate_sample_data(coords, mdata, phi, bhat, zeros, state)
    type(coordinates_struct), pointer, intent(inout) :: coords
    type(magneticdata_struct), pointer, intent(inout) :: mdata
    real(kind=dp), intent(in) :: phi(0:2), bhat(0:2), zeros(0:2), state(1:2)

    allocate(coords%x1(0:2), coords%x2(0:2), coords%x3(0:2))
    coords%x1 = zeros
    coords%x2 = phi
    coords%x3 = zeros
    allocate(mdata%bhat(0:2), mdata%geodcu(0:2), mdata%h_phi(0:2), &
         mdata%dlogbdphi(0:2), mdata%ybeg(1:2), mdata%yend(1:2))
    mdata%bhat = bhat
    mdata%geodcu = zeros
    mdata%h_phi = zeros
    mdata%dlogbdphi = zeros
    mdata%ybeg = state
    mdata%yend = state
  end subroutine allocate_sample_data

  subroutine assert_close(actual, expected, label, status)
    real(kind=dp), intent(in) :: actual(:), expected(:)
    character(len=*), intent(in) :: label
    integer, intent(inout) :: status

    if (maxval(abs(actual - expected)) > tol * max(1.0_dp, maxval(abs(expected)))) then
       print *, 'FAIL: ', label, actual, expected
       status = status + 1
    end if
  end subroutine assert_close

  subroutine assert_scalar(actual, expected, label, status)
    real(kind=dp), intent(in) :: actual, expected
    character(len=*), intent(in) :: label
    integer, intent(inout) :: status

    if (abs(actual - expected) > tol * max(1.0_dp, abs(expected))) then
       print *, 'FAIL: ', label, actual, expected
       status = status + 1
    end if
  end subroutine assert_scalar

end program test_magnetics_pitch_support
