module test_rk4_kin_retry_rhs
  implicit none
  double precision, parameter :: x_on = 0.25d0
contains
  ! y(1:10)' = cos(x); y(10+k)' = k for x > x_on, else 0. Started at x_on,
  ! y(11:14) are exactly zero with zero derivative, like the running
  ! integrals of the field line.
  subroutine rhs_switch_on(x, y, dydx)
    double precision, intent(in)  :: x
    double precision, intent(in)  :: y(:)
    double precision, intent(out) :: dydx(:)
    integer :: k

    dydx(1:10) = cos(x)
    do k = 1, 4
      dydx(10 + k) = merge(dble(k), 0.0d0, x > x_on)
    end do
  end subroutine rhs_switch_on
end module test_rk4_kin_retry_rhs

program test_rk4_kin_retry
  use odeint_allroutines_sub, only: odeint_allroutines, odeint_has_failed
  use rk4_kin_mod, only: odeint_kin_retry
  use test_rk4_kin_retry_rhs
  implicit none

  double precision, parameter :: eps = 1.0d-10, x_end = 1.25d0
  double precision :: y(14), y_ref(14)
  integer :: k
  logical :: failed

  y_ref(1:10) = sin(x_end)
  y_ref(11:14) = [(k*(x_end - x_on), k = 1, 4)]

  y(1:10) = sin(x_on)
  y(11:14) = 0.0d0
  call odeint_allroutines(y, 14, x_on, x_end, eps, rhs_switch_on)
  if (.not. odeint_has_failed()) error stop 'FAIL: plain odeint did not fail'

  y(1:10) = sin(x_on)
  y(11:14) = 0.0d0
  call odeint_kin_retry(y, x_on, x_end, eps, rhs_switch_on, failed)
  if (failed) error stop 'FAIL: retry reported failure'
  if (maxval(abs(y - y_ref)) > 1.0d-8) then
    print *, 'FAIL: max error ', maxval(abs(y - y_ref))
    error stop 1
  end if
  print *, 'All tests passed!'
end program test_rk4_kin_retry
