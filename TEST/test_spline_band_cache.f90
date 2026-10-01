program test_spline_band_cache
  ! splinecof3 reuses the band LU factorization when consecutive calls share
  ! grid, test function and weights. Oracles: a result obtained through the
  ! cache must be bit-identical to a fresh factorization of the same system,
  ! and every call must interpolate its own data (a stale factorization after
  ! a change of grid or test function would not).
  use inter_interfaces, only: splinecof3, splint_horner3
  use nrtype, only: DP, I4B

  implicit none

  integer(I4B), parameter :: n = 9
  real(DP) :: x1(n), x2(n), y1(n), y2(n), lambda(n)
  real(DP), dimension(n) :: a, b, c, d, a_hit, b_hit, c_hit, d_hit
  integer(I4B) :: indx(n), i

  x1 = [0.0_DP, 0.08_DP, 0.21_DP, 0.37_DP, 0.52_DP, 0.70_DP, 0.81_DP, 0.93_DP, 1.0_DP]
  x2 = [0.0_DP, 0.11_DP, 0.19_DP, 0.33_DP, 0.55_DP, 0.64_DP, 0.80_DP, 0.90_DP, 1.0_DP]
  y1 = sin(2.0_DP*x1) + 0.25_DP*x1*x1
  y2 = exp(-x1) * cos(3.0_DP*x1)
  lambda = 1.0_DP
  indx = [(i, i = 1, n)]

  ! fresh factorization, then a cache hit for different data
  call spline(x1, y1, unit_weight, a, b, c, d)
  call assert_interpolates(x1, y1, a, b, c, d, unit_weight, 1.0_DP)
  call spline(x1, y2, unit_weight, a_hit, b_hit, c_hit, d_hit)
  call assert_interpolates(x1, y2, a_hit, b_hit, c_hit, d_hit, unit_weight, 1.0_DP)

  ! new grid of the same size: must refactorize
  call spline(x2, y1, unit_weight, a, b, c, d)
  call assert_interpolates(x2, y1, a, b, c, d, unit_weight, 1.0_DP)

  ! new test function on the same grid: must refactorize
  call spline(x2, y1, two_weight, a, b, c, d)
  call assert_interpolates(x2, y1, a, b, c, d, two_weight, 2.0_DP)

  ! fresh factorization of the earlier system: bit-identical to the cache hit
  call spline(x1, y2, unit_weight, a, b, c, d)
  if (any(a /= a_hit) .or. any(b /= b_hit) .or. any(c /= c_hit) &
      .or. any(d /= d_hit)) then
    print *, "FAIL: cached solve differs from fresh solve"
    stop 1
  end if

  print *, "All tests passed!"

contains

  subroutine spline(x, y, f, a_coef, b_coef, c_coef, d_coef)
    real(DP), intent(in) :: x(:), y(:)
    real(DP), intent(out) :: a_coef(:), b_coef(:), c_coef(:), d_coef(:)
    interface
      function f(x_value, m_value)
        import :: DP
        real(DP), intent(in) :: x_value, m_value
        real(DP) :: f
      end function f
    end interface
    real(DP) :: c1, cn

    c1 = 0.0_DP
    cn = 0.0_DP
    call splinecof3(x, y, c1, cn, lambda, indx, 2_I4B, 4_I4B, &
         a_coef, b_coef, c_coef, d_coef, 0.0_DP, f)
  end subroutine spline

  subroutine assert_interpolates(x_values, y_values, a_coef, b_coef, c_coef, &
       d_coef, f, f_value)
    real(DP), intent(in) :: x_values(:), y_values(:)
    real(DP), intent(in) :: a_coef(:), b_coef(:), c_coef(:), d_coef(:)
    real(DP), intent(in) :: f_value
    interface
      function f(x_value, m_value)
        import :: DP
        real(DP), intent(in) :: x_value, m_value
        real(DP) :: f
      end function f
    end interface
    real(DP) :: y_eval, yp, ypp, yppp, max_err

    max_err = 0.0_DP
    do i = 1, size(x_values)
      call splint_horner3(x_values, a_coef, b_coef, c_coef, d_coef, 0_I4B, &
           0.0_DP, x_values(i), f, zero_weight, zero_weight, &
           zero_weight, y_eval, yp, ypp, yppp)
      max_err = max(max_err, abs(y_eval - y_values(i)))
    end do
    ! s(x) = y/f at the nodes: the coefficients scale with 1/f
    if (max_err > 1.0e-12_DP .or. abs(a_coef(2)*f_value - y_values(2)) > 1.0e-12_DP) then
      print *, "FAIL: spline does not interpolate its data", max_err
      stop 1
    end if
  end subroutine assert_interpolates

  real(DP) function unit_weight(x_value, m_value)
    real(DP), intent(in) :: x_value, m_value

    unit_weight = 1.0_DP + 0.0_DP*x_value + 0.0_DP*m_value
  end function unit_weight

  real(DP) function two_weight(x_value, m_value)
    real(DP), intent(in) :: x_value, m_value

    two_weight = 2.0_DP + 0.0_DP*x_value + 0.0_DP*m_value
  end function two_weight

  real(DP) function zero_weight(x_value, m_value)
    real(DP), intent(in) :: x_value, m_value

    zero_weight = 0.0_DP*x_value + 0.0_DP*m_value
  end function zero_weight
end program test_spline_band_cache
