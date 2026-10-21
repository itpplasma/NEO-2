program test_radial_lagrange
  ! Local Lagrange radial coefficients evaluated with splint_horner3 against
  ! y(s) = s^m g(s): exact for cubic g, order 4 (value) and 3 (d/ds) for smooth
  ! g, and the default switch reproduces splinecof3_hi_driv bitwise.
  use nrtype, only : I4B, DP
  use inter_interfaces, only : splinecof3_hi_driv, splint_horner3, &
    tf, tfp, tfpp, tfppp
  use neo_spline_data, only : lsw_lagrange_boozer
  use radial_lagrange_cof, only : lagrange_cof3, radial_cof3_hi_driv
  implicit none

  logical :: ok = .true.
  real(DP) :: e1(2), e2(2)

  call expect(cubic_err(0.0_DP, 0.0_DP) < 1.0e-11_DP, 'cubic m=0')
  call expect(cubic_err(1.5_DP, 0.0_DP) < 1.0e-11_DP, 'cubic m=1.5')
  call expect(cubic_err(2.0_DP, 0.05_DP) < 1.0e-11_DP, 'cubic m=2')
  e1 = max_err2(100)
  e2 = max_err2(200)
  print '(a,2f6.2)', 'observed orders value, d/ds: ', log(e1 / e2) / log(2.0_DP)
  call expect(e1(1) / e2(1) > 12.0_DP, 'value order 4')
  call expect(e1(2) / e2(2) > 6.0_DP, 'derivative order 3')
  call check_switch()

  if (.not. ok) error stop 'FAIL'
  print *, 'All tests passed!'

contains

  subroutine expect(cond, msg)
    logical, intent(in) :: cond
    character(*), intent(in) :: msg
    if (.not. cond) then
      ok = .false.
      print *, 'FAIL: ', msg
    end if
  end subroutine expect

  pure function g(s, cubic) result(r)
    ! Value and d/ds of the profile factor.
    real(DP), intent(in) :: s
    logical, intent(in) :: cubic
    real(DP) :: r(2)
    if (cubic) then
      r = [0.7_DP - 1.3_DP*s + 2.1_DP*s**2 - 0.9_DP*s**3, &
           -1.3_DP + 4.2_DP*s - 2.7_DP*s**2]
    else
      r = [exp(sin(3.0_DP*s)), 3.0_DP*cos(3.0_DP*s)*exp(sin(3.0_DP*s))]
    end if
  end function g

  function cubic_err(m, s0) result(err)
    ! Max error of value and d/ds for cubic g on 9 non-uniform nodes from s0
    ! to 1, sampled in every interval including the first one.
    real(DP), intent(in) :: m, s0
    integer(I4B), parameter :: n = 9
    real(DP) :: err, x(n), y(n), a(n), b(n), c(n), d(n), s, yv, ypv, ypp, yppp
    real(DP) :: gs(2)
    integer(I4B) :: j, indx(n)

    x = s0 + (1.0_DP - s0) * [((real(j - 1, DP) / (n - 1))**1.3_DP, j = 1, n)]
    indx = [(j, j = 1, n)]
    do j = 1, n
      gs = g(x(j), .true.)
      y(j) = tf(x(j), m) * gs(1)
    end do
    call lagrange_cof3(x, y, m, a, b, c, d, indx, tf)
    err = 0.0_DP
    do j = 1, n - 1
      s = x(j) + 0.37_DP * (x(j + 1) - x(j))
      call splint_horner3(x, a, b, c, d, 1_I4B, m, s, tf, tfp, tfpp, tfppp, &
        yv, ypv, ypp, yppp)
      gs = g(s, .true.)
      err = max(err, abs(yv - s**m * gs(1)), &
        abs(ypv - s**m * gs(2) - m * s**(m - 1.0_DP) * gs(1)))
    end do
  end function cubic_err

  function max_err2(n) result(err)
    ! Value and derivative errors (m = 0) for the smooth profile.
    integer(I4B), intent(in) :: n
    real(DP) :: err(2), x(n), y(n), a(n), b(n), c(n), d(n), s, yv, ypv, ypp, yppp
    real(DP) :: gs(2)
    integer(I4B) :: j, indx(n)

    x = [((real(j - 1, DP) / (n - 1))**1.3_DP, j = 1, n)]
    indx = [(j, j = 1, n)]
    do j = 1, n
      gs = g(x(j), .false.)
      y(j) = gs(1)
    end do
    call lagrange_cof3(x, y, 0.0_DP, a, b, c, d, indx, tf)
    err = 0.0_DP
    do j = n / 10, n - 1
      s = 0.5_DP * (x(j) + x(j + 1))
      call splint_horner3(x, a, b, c, d, 1_I4B, 0.0_DP, s, tf, tfp, tfpp, &
        tfppp, yv, ypv, ypp, yppp)
      gs = g(s, .false.)
      err = max(err, abs([yv, ypv] - gs))
    end do
  end function max_err2

  subroutine check_switch()
    ! Off: bitwise splinecof3_hi_driv. On: column-wise lagrange_cof3.
    integer(I4B), parameter :: n = 40
    real(DP) :: x(n), y(n, 2), mh(2), gs(2)
    real(DP), dimension(n, 2) :: a, b, c, d, a0, b0, c0, d0
    integer(I4B) :: i, j, indx(n)

    x = [((real(j - 1, DP) / (n - 1))**1.3_DP, j = 1, n)]
    indx = [(j, j = 1, n)]
    mh = [0.0_DP, 0.5_DP]
    do i = 1, 2
      do j = 1, n
        gs = g(x(j), .false.)
        y(j, i) = tf(x(j), mh(i)) * gs(1)
      end do
    end do

    lsw_lagrange_boozer = .false.
    call radial_cof3_hi_driv(x, y, mh, a, b, c, d, indx, tf)
    call splinecof3_hi_driv(x, y, mh, a0, b0, c0, d0, indx, tf)
    call expect(all([a == a0, b == b0, c == c0, d == d0]), 'switch off')

    lsw_lagrange_boozer = .true.
    call radial_cof3_hi_driv(x, y, mh, a, b, c, d, indx, tf)
    do i = 1, 2
      call lagrange_cof3(x, y(:,i), mh(i), a0(:,i), b0(:,i), c0(:,i), &
        d0(:,i), indx, tf)
    end do
    call expect(all([a == a0, b == b0, c == c0, d == d0]), 'switch on')
  end subroutine check_switch

end program test_radial_lagrange
