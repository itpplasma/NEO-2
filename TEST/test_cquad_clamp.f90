program test_cquad_clamp
  ! fint1d_cquad must never sample outside [x_low, x_up]. The segment
  ! [8/3, 3] reproduces the knot-bounded segment traced in the ql golden
  ! case: the first Clenshaw-Curtis node mid - half rounds to
  ! 2.66666666666666607, one ulp below x_low. The integrands below are large
  ! outside the segment, as the collision-operator integrands are in the
  ! neighbouring knot cell. They must never be sampled there, and they must
  ! integrate to their exact values over the segment.
  use, intrinsic :: iso_fortran_env, only: wp => real64
  use gsl_integration_routines_mod, only: fint1d_cquad

  implicit none

  real(wp), parameter :: x_low = 8.0_wp/3.0_wp
  real(wp), parameter :: x_up = 3.0_wp
  real(wp), parameter :: outside = 5.4e5_wp
  real(wp), parameter :: epsabs = 1.0e-13_wp
  real(wp), parameter :: epsrel = 1.0e-13_wp
  real(wp) :: a, b, exact, res(2)
  logical :: sampled_outside = .false.

  a = x_low
  b = x_up
  if (0.5_wp*(a + b) - 0.5_wp*(b - a) >= a) then
    print *, "FAIL: segment no longer has a node below x_low"
    stop 1
  end if

  res = fint1d_cquad(step_zero, a, b, epsabs, epsrel)
  call assert_close(res(1), 0.0_wp, 0.0_wp, "cquad zero integrand on segment")

  exact = 0.5_wp*(b - a)**2
  res = fint1d_cquad(step_ramp, a, b, epsabs, epsrel)
  call assert_close(res(1), exact, 1.0e-14_wp*exact, "cquad ramp on segment")

  res = fint1d_cquad(step_scaled, 2.0_wp, a, b, epsabs, epsrel)
  call assert_close(res(1), 2.0_wp*exact, 2.0e-14_wp*exact, &
       "cquad parameter ramp on segment")

  if (sampled_outside) then
    print *, "FAIL: cquad sampled the integrand outside [x_low, x_up]"
    stop 1
  end if

  print *, "All tests passed!"

contains

  subroutine assert_close(actual, expected, tol, label)
    real(wp), intent(in) :: actual, expected, tol
    character(*), intent(in) :: label

    if (abs(actual - expected) > tol) then
      print *, "FAIL: ", label, " got ", actual, " expected ", expected
      stop 1
    end if
  end subroutine assert_close

  function step_zero(x)
    real(wp) :: step_zero
    real(wp) :: x

    if (x < x_low .or. x > x_up) then
      sampled_outside = .true.
      step_zero = outside
    else
      step_zero = 0.0_wp
    end if
  end function step_zero

  function step_ramp(x)
    real(wp) :: step_ramp
    real(wp) :: x

    if (x < x_low .or. x > x_up) then
      sampled_outside = .true.
      step_ramp = outside
    else
      step_ramp = x - x_low
    end if
  end function step_ramp

  function step_scaled(x, param1)
    real(wp) :: step_scaled
    real(wp) :: x, param1

    if (x < x_low .or. x > x_up) then
      sampled_outside = .true.
      step_scaled = outside
    else
      step_scaled = param1*(x - x_low)
    end if
  end function step_scaled
end program test_cquad_clamp
