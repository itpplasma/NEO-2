program test_integration_routines
  use, intrinsic :: iso_fortran_env, only: wp => real64
  use gsl_integration_routines_mod, only: fint1d_qag, fint1d_qags, &
       fint1d_qagp, fint1d_cquad, fint1d_qagiu

  implicit none

  real(wp), parameter :: pi = acos(-1.0_wp)
  real(wp), parameter :: epsabs = 1.0e-10_wp
  real(wp), parameter :: epsrel = 1.0e-10_wp
  real(wp) :: res(2)

  res = fint1d_qag(cosine, 0.0_wp, 0.5_wp*pi, epsabs, epsrel, 3)
  call assert_close(res(1), 1.0_wp, 1.0e-9_wp, "qag cos over [0,pi/2]")

  res = fint1d_qags(inv_sqrt, 0.0_wp, 1.0_wp, epsabs, epsrel)
  call assert_close(res(1), 2.0_wp, 1.0e-6_wp, "qags 1/sqrt(x) over [0,1]")

  res = fint1d_qagp(abs_kink, [0.0_wp, 1.0_wp, 2.0_wp], epsabs, epsrel)
  call assert_close(res(1), 1.0_wp, 1.0e-9_wp, "qagp |x-1| over [0,2]")

  res = fint1d_cquad(cosine, 0.0_wp, 0.5_wp*pi, epsabs, epsrel)
  call assert_close(res(1), 1.0_wp, 1.0e-9_wp, "cquad cos over [0,pi/2]")

  res = fint1d_qagiu(decaying, 0.0_wp, epsabs, epsrel)
  call assert_close(res(1), 1.0_wp, 1.0e-9_wp, "qagiu exp(-x) over [0,inf)")

  res = fint1d_cquad(decaying_param, 2.0_wp, 0.0_wp, 1.0_wp, epsabs, epsrel)
  call assert_close(res(1), 0.5_wp*(1.0_wp - exp(-2.0_wp)), 1.0e-9_wp, &
       "cquad exp(-p*x) over [0,1]")

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

  function cosine(x)
    real(wp) :: cosine
    real(wp) :: x

    cosine = cos(x)
  end function cosine

  function inv_sqrt(x)
    real(wp) :: inv_sqrt
    real(wp) :: x

    inv_sqrt = 1.0_wp/sqrt(x)
  end function inv_sqrt

  function abs_kink(x)
    real(wp) :: abs_kink
    real(wp) :: x

    abs_kink = abs(x - 1.0_wp)
  end function abs_kink

  function decaying(x)
    real(wp) :: decaying
    real(wp) :: x

    decaying = exp(-x)
  end function decaying

  function decaying_param(x, param1)
    real(wp) :: decaying_param
    real(wp) :: x, param1

    decaying_param = exp(-param1*x)
  end function decaying_param
end program test_integration_routines
