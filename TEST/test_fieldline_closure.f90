module test_fieldline_closure_rhs
  implicit none
  double precision, parameter :: x_on = 0.25d0
contains
  ! y1' = cos(x), y2' = 1 for x > x_on, else 0.
  ! Started at x = x_on with y2 = 0, the second component is exactly zero
  ! with zero derivative, like the running integrals of the field line.
  subroutine rhs_switch_on(x, y, dydx)
    double precision, intent(in)  :: x
    double precision, intent(in)  :: y(:)
    double precision, intent(out) :: dydx(:)

    dydx(1) = cos(x)
    dydx(2) = merge(1.0d0, 0.0d0, x > x_on)
  end subroutine rhs_switch_on

  subroutine rhs_smooth(x, y, dydx)
    double precision, intent(in)  :: x
    double precision, intent(in)  :: y(:)
    double precision, intent(out) :: dydx(:)

    dydx(1) = cos(x)
    dydx(2) = -y(2)
  end subroutine rhs_smooth
end module test_fieldline_closure_rhs

program test_fieldline_closure
  use odeint_allroutines_sub, only: odeint_allroutines, odeint_has_failed
  use odeint_robust_mod, only: odeint_robust
  use mag_interface_mod, only: tokamak_fieldline_closed
  use test_fieldline_closure_rhs
  implicit none

  double precision, parameter :: pi = 3.14159265358979d0
  double precision, parameter :: eps = 1.0d-10, x_end = 1.25d0
  double precision :: y(2), y_ref(2)
  integer :: ierr, nfail

  nfail = 0

  ! 1. Zero component with zero start derivative: plain odeint gives up,
  !    odeint_robust must integrate to the analytic solution.
  y = [sin(x_on), 0.0d0]
  call odeint_allroutines(y, 2, x_on, x_end, eps, rhs_switch_on)
  print *, 'plain odeint failed on zero-start component: ', odeint_has_failed()

  y = [sin(x_on), 0.0d0]
  call odeint_robust(y, x_on, x_end, eps, rhs_switch_on, ierr)
  call check(ierr == 0, 'odeint_robust reports success')
  call check(abs(y(1) - sin(x_end)) < 1.0d-8, 'smooth component y1 = sin(x)')
  call check(abs(y(2) - (x_end - x_on)) < 1.0d-8, 'zero-start component y2 = x - x_on')

  ! 2. Regular problem: odeint_robust must reproduce odeint bit for bit.
  y_ref = [0.5d0, 2.0d0]
  call odeint_allroutines(y_ref, 2, 0.0d0, x_end, eps, rhs_smooth)
  y = [0.5d0, 2.0d0]
  call odeint_robust(y, 0.0d0, x_end, eps, rhs_smooth, ierr)
  call check(ierr == 0 .and. all(y == y_ref), 'bitwise equal to odeint when it succeeds')
  call check(abs(y(1) - (0.5d0 + sin(x_end))) < 1.0d-8 .and. &
       abs(y(2) - 2.0d0*exp(-x_end)) < 1.0d-8, 'regular problem matches analytic solution')

  ! 3. Tokamak closure test: wrap-around of theta and minimum period count.
  call check(tokamak_fieldline_closed(2, 2, 0.2d0, 0.2d0 + 1.0d-8), 'closes at expected period')
  call check(tokamak_fieldline_closed(2, 2, pi - 1.0d-4, -pi + 1.0d-4), 'closes across theta = pi')
  call check(.not. tokamak_fieldline_closed(1, 2, 0.2d0, 0.2d0), 'no closure before expected period')
  call check(.not. tokamak_fieldline_closed(3, 2, 0.2d0, 0.2d0 + 2.0d-3), 'open line is not closed')

  if (nfail == 0) then
    print *, 'All tests passed!'
  else
    print *, 'FAIL: ', nfail, ' checks failed'
    error stop 1
  end if

contains

  subroutine check(ok, name)
    logical, intent(in) :: ok
    character(*), intent(in) :: name
    if (ok) then
      print *, 'PASS: ', name
    else
      print *, 'FAIL: ', name
      nfail = nfail + 1
    end if
  end subroutine check

end program test_fieldline_closure
