program test_spline_sparse_assembly
  use inter_interfaces, only: splinecof3, splint_horner3
  use nrtype, only: DP, I4B

  implicit none

  integer(I4B), parameter :: n = 8
  real(DP) :: x(n), y(n), lambda(n)
  real(DP) :: a(n), b(n), c(n), d(n)
  integer(I4B) :: indx(n)
  real(DP) :: c1, cn, y_eval, yp, ypp, yppp
  real(DP) :: max_err
  integer(I4B) :: i

  x = [0.0_DP, 0.08_DP, 0.21_DP, 0.37_DP, 0.52_DP, 0.70_DP, 0.86_DP, 1.0_DP]
  y = sin(2.0_DP*x) + 0.25_DP*x*x
  lambda = 1.0_DP
  indx = [(i, i = 1, n)]
  c1 = 0.0_DP
  cn = 0.0_DP

  call splinecof3(x, y, c1, cn, lambda, indx, 1_I4B, 2_I4B, &
       a, b, c, d, 1.0_DP, unit_weight)

  max_err = 0.0_DP
  do i = 1, n
    call splint_horner3(x, a, b, c, d, 0_I4B, 1.0_DP, x(i), &
         unit_weight, zero_weight, zero_weight, zero_weight, y_eval, yp, ypp, yppp)
    max_err = max(max_err, abs(y_eval - y(i)))
  end do

  if (max_err > 1.0e-10_DP) then
    print *, "FAIL: spline reconstruction mismatch", max_err
    stop 1
  end if

  print *, "All tests passed!"

contains

  real(DP) function unit_weight(x_value, m_value)
    real(DP), intent(in) :: x_value, m_value

    unit_weight = 1.0_DP + 0.0_DP*x_value + 0.0_DP*m_value
  end function unit_weight

  real(DP) function zero_weight(x_value, m_value)
    real(DP), intent(in) :: x_value, m_value

    zero_weight = 0.0_DP*x_value + 0.0_DP*m_value
  end function zero_weight
end program test_spline_sparse_assembly
