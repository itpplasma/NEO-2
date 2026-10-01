program test_radial_lagrange
  ! Local Lagrange radial interpolation (radial_lagrange_cof) against
  ! analytic profiles y(s) = s^m g(s), evaluated with splint_horner3 as in
  ! neo_magfie:
  !  1. exact for cubic g, including the axis interval and m /= 0,
  !  2. fourth-order convergence of y and third-order of dy/ds for smooth g,
  !  3. the switch selects Lagrange or spline in radial_cof3_hi_driv.
  use nrtype, only : I4B, DP
  use inter_interfaces, only : splint_horner3, tf, tfp, tfpp, tfppp
  use radial_lagrange_cof, only : lsw_lagrange_boozer, lagrange_cof3, &
    radial_cof3_hi_driv
  implicit none

  logical :: ok = .true.
  real(DP) :: e1(2), e2(2), e3(2)

  call check_cubic(0.0_DP, .true.)
  call check_cubic(1.5_DP, .true.)
  call check_cubic(2.0_DP, .false.)

  call max_err(50, 1.0_DP, e1)
  call max_err(100, 1.0_DP, e2)
  call max_err(200, 1.0_DP, e3)
  write (*,'(a,3es10.2)') 'value errors n=50,100,200: ', e1(1), e2(1), e3(1)
  write (*,'(a,3es10.2)') 'deriv errors n=50,100,200: ', e1(2), e2(2), e3(2)
  call expect(e2(1) / e3(1) > 12.0_DP, 'value converges at order 4')
  call expect(e2(2) / e3(2) > 6.0_DP, 'derivative converges at order 3')
  call expect(e3(1) < 1.0e-7_DP, 'value accurate at n=200')

  call check_switch()

  if (ok) then
    print *, 'All tests passed!'
  else
    print *, 'FAIL'
    error stop 1
  end if

contains

  subroutine expect(cond, msg)
    logical, intent(in) :: cond
    character(*), intent(in) :: msg
    if (.not. cond) then
      ok = .false.
      print *, 'FAIL: ', msg
    end if
  end subroutine expect

  pure real(DP) function gcub(s)
    real(DP), intent(in) :: s
    gcub = 0.7_DP - 1.3_DP * s + 2.1_DP * s**2 - 0.9_DP * s**3
  end function gcub

  pure real(DP) function gcub_s(s)
    real(DP), intent(in) :: s
    gcub_s = -1.3_DP + 4.2_DP * s - 2.7_DP * s**2
  end function gcub_s

  pure real(DP) function gsm(s)
    real(DP), intent(in) :: s
    gsm = exp(sin(3.0_DP * s))
  end function gsm

  pure real(DP) function gsm_s(s)
    real(DP), intent(in) :: s
    gsm_s = 3.0_DP * cos(3.0_DP * s) * exp(sin(3.0_DP * s))
  end function gsm_s

  subroutine nodes(n, x)
    ! Non-uniform grid, starting at the axis.
    integer(I4B), intent(in) :: n
    real(DP), intent(out) :: x(n)
    integer(I4B) :: j
    x = [((real(j - 1, DP) / real(n - 1, DP))**1.3_DP, j = 1, n)]
  end subroutine nodes

  subroutine eval(x, a, b, c, d, m, s, y, yp)
    real(DP), intent(in) :: x(:), a(:), b(:), c(:), d(:), m, s
    real(DP), intent(out) :: y, yp
    real(DP) :: ypp, yppp
    call splint_horner3(x, a, b, c, d, 1_I4B, m, s, tf, tfp, tfpp, tfppp, &
      y, yp, ypp, yppp)
  end subroutine eval

  subroutine check_cubic(m, axis)
    real(DP), intent(in) :: m
    logical, intent(in) :: axis
    integer(I4B), parameter :: n = 9, ne = 101
    real(DP) :: x(n), y(n), a(n), b(n), c(n), d(n)
    real(DP) :: s, yv, ypv, yex, ypex, err
    integer(I4B) :: j, indx(n)
    character(64) :: msg

    call nodes(n, x)
    if (.not. axis) x = 0.05_DP + 0.9_DP * x
    indx = [(j, j = 1, n)]
    do j = 1, n
      y(j) = tf(x(j), m) * gcub(x(j))
    end do
    call lagrange_cof3(x, y, m, a, b, c, d, indx, tf)
    err = 0.0_DP
    do j = 1, ne
      s = x(1) + 1.0e-3_DP + (x(n) - x(1) - 2.0e-3_DP) * real(j - 1, DP) / (ne - 1)
      call eval(x, a, b, c, d, m, s, yv, ypv)
      yex = s**m * gcub(s)
      ypex = s**m * gcub_s(s)
      if (m /= 0.0_DP) ypex = ypex + m * s**(m - 1.0_DP) * gcub(s)
      err = max(err, abs(yv - yex), abs(ypv - ypex))
    end do
    write (msg, '(a,f4.1)') 'cubic profile reproduced exactly, m=', m
    call expect(err < 1.0e-11_DP, trim(msg))
  end subroutine check_cubic

  subroutine max_err(n, m, err)
    integer(I4B), intent(in) :: n
    real(DP), intent(in) :: m
    real(DP), intent(out) :: err(2)
    real(DP) :: x(n), y(n), a(n), b(n), c(n), d(n), s, yv, ypv
    integer(I4B) :: j, indx(n)

    call nodes(n, x)
    indx = [(j, j = 1, n)]
    do j = 1, n
      y(j) = tf(x(j), m) * gsm(x(j))
    end do
    call lagrange_cof3(x, y, m, a, b, c, d, indx, tf)
    err = 0.0_DP
    ! Midpoints away from the axis, where the grid is quasi-uniform
    do j = n / 10, n - 1
      s = 0.5_DP * (x(j) + x(j + 1))
      call eval(x, a, b, c, d, m, s, yv, ypv)
      err(1) = max(err(1), abs(yv - s**m * gsm(s)))
      err(2) = max(err(2), abs(ypv - (s**m * gsm_s(s) &
        + m * s**(m - 1.0_DP) * gsm(s))))
    end do
  end subroutine max_err

  subroutine check_switch()
    ! Lagrange is exact for the cubic column, the natural spline is not;
    ! both agree with the smooth column to interpolation accuracy.
    integer(I4B), parameter :: n = 40, nc = 2
    real(DP) :: x(n), y(n, nc), mh(nc), s, yv, ypv
    real(DP), dimension(n, nc) :: a, b, c, d
    real(DP) :: err_lag(nc), err_spl(nc)
    integer(I4B) :: j, indx(n)

    call nodes(n, x)
    indx = [(j, j = 1, n)]
    mh = [1.0_DP, 0.5_DP]
    do j = 1, n
      y(j, 1) = tf(x(j), mh(1)) * gcub(x(j))
      y(j, 2) = tf(x(j), mh(2)) * gsm(x(j))
    end do
    s = 0.5_DP * (x(n - 1) + x(n))

    lsw_lagrange_boozer = .true.
    call radial_cof3_hi_driv(x, y, mh, a, b, c, d, indx, tf)
    call eval(x, a(:,1), b(:,1), c(:,1), d(:,1), mh(1), s, yv, ypv)
    err_lag(1) = abs(yv - s * gcub(s))
    call eval(x, a(:,2), b(:,2), c(:,2), d(:,2), mh(2), s, yv, ypv)
    err_lag(2) = abs(yv - sqrt(s) * gsm(s))

    lsw_lagrange_boozer = .false.
    call radial_cof3_hi_driv(x, y, mh, a, b, c, d, indx, tf)
    call eval(x, a(:,1), b(:,1), c(:,1), d(:,1), mh(1), s, yv, ypv)
    err_spl(1) = abs(yv - s * gcub(s))
    call eval(x, a(:,2), b(:,2), c(:,2), d(:,2), mh(2), s, yv, ypv)
    err_spl(2) = abs(yv - sqrt(s) * gsm(s))

    write (*,'(a,2es10.2,a,2es10.2)') 'last interval errors: lagrange', &
      err_lag, ' spline', err_spl
    call expect(err_lag(1) < 1.0e-12_DP, 'switch on selects Lagrange')
    call expect(err_spl(1) > 1.0e-8_DP, 'switch off selects natural spline')
    call expect(maxval(err_lag(2:)) < 1.0e-4_DP .and. &
      maxval(err_spl(2:)) < 1.0e-3_DP, 'both accurate for smooth data')
  end subroutine check_switch

end program test_radial_lagrange
