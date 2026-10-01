!> Regression test for theta_rz_solver_mod (used by calc_thetaB_RZloc).
!> Oracle: points generated from a known theta on an analytic, up-down
!> asymmetric, shaped curve; the solver must recover that theta. Includes the
!> three Boozer angles (24, 10 and 90 times 2 pi/100) at which the former
!> separate R/Z Newton iterations stopped on the golden-record ql surface,
!> lower-half-plane points and points next to theta = 0.
program test_theta_rz_solver
  use nrtype, only: dp, twopi
  use theta_rz_solver_mod, only: find_theta_of_rz
  implicit none

  integer :: nfail, i
  real(dp) :: th0, R0, Z0, R_tb, Z_tb, th, dist, extent
  real(dp), parameter :: special(5) = [ &
       24.0_dp * twopi / 100.0_dp, &   ! theta_B = 1.508: Newton gave 0.734 / 1.571
       10.0_dp * twopi / 100.0_dp, &   ! theta_B = 0.628: Newton gave 5.822 / 0.628
       90.0_dp * twopi / 100.0_dp, &   ! theta_B = 5.655: Newton gave 5.655 / 0.793
       1.0e-9_dp, twopi - 1.0e-9_dp ]

  nfail = 0

  ! 1) recover theta on a dense set of angles, both half planes
  do i = 0, 240
     th0 = real(i, dp) * twopi / 241.0_dp + 1.0e-3_dp
     call check_on_curve(th0)
  end do
  do i = 1, size(special)
     call check_on_curve(special(i))
  end do

  ! 2) coarse scan still selects the right branch
  call shaped(1.508_dp, R0, R_tb, Z0, Z_tb)
  call find_theta_of_rz(shaped, R0, Z0, th, dist, extent, nscan=12)
  call expect(angle_diff(th, 1.508_dp) < 1.0e-9_dp, 'coarse scan, theta = 1.508')

  ! 3) a point 10 cm off the curve along the outward normal at theta = 1:
  !    nearest point is theta = 1 and the distance is reported as 10 cm
  call shaped(1.0_dp, R0, R_tb, Z0, Z_tb)
  R0 = R0 + 10.0_dp * Z_tb / hypot(R_tb, Z_tb)
  Z0 = Z0 - 10.0_dp * R_tb / hypot(R_tb, Z_tb)
  call find_theta_of_rz(shaped, R0, Z0, th, dist, extent)
  call expect(abs(dist - 10.0_dp) < 1.0e-6_dp .and. dist > 1.0e-6_dp * extent, &
       'off-curve point: distance reported')

  if (nfail == 0) then
     print *, 'All tests passed!'
  else
     print *, 'FAIL: ', nfail, ' check(s)'
     error stop 1
  end if

contains

  !> Shaped, up-down asymmetric curve (cm), similar in size to the golden ql
  !> surface: R = 170 + 35 cos(t + 0.3 sin t), Z = 2 + 55 sin t + 3 cos 2t
  subroutine shaped(t, R, R_t, Z, Z_t)
    real(dp), intent(in)  :: t
    real(dp), intent(out) :: R, R_t, Z, Z_t
    R = 170.0_dp + 35.0_dp * cos(t + 0.3_dp * sin(t))
    R_t = -35.0_dp * sin(t + 0.3_dp * sin(t)) * (1.0_dp + 0.3_dp * cos(t))
    Z = 2.0_dp + 55.0_dp * sin(t) + 3.0_dp * cos(2.0_dp * t)
    Z_t = 55.0_dp * cos(t) - 6.0_dp * sin(2.0_dp * t)
  end subroutine shaped

  real(dp) function angle_diff(a, b)
    real(dp), intent(in) :: a, b
    angle_diff = abs(modulo(a - b + 0.5_dp * twopi, twopi) - 0.5_dp * twopi)
  end function angle_diff

  subroutine check_on_curve(t)
    real(dp), intent(in) :: t
    real(dp) :: Rt, Zt, Rt_tb, Zt_tb, tt, d, ext
    character(len=64) :: label
    call shaped(t, Rt, Rt_tb, Zt, Zt_tb)
    call find_theta_of_rz(shaped, Rt, Zt, tt, d, ext)
    write (label, '(a,es12.5)') 'recover theta = ', t
    call expect(angle_diff(tt, t) < 1.0e-9_dp .and. d < 1.0e-9_dp &
         .and. tt >= 0.0_dp .and. tt < twopi, trim(label))
  end subroutine check_on_curve

  subroutine expect(ok, label)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: label
    if (.not. ok) then
       nfail = nfail + 1
       print *, 'FAIL ', label
    end if
  end subroutine expect

end program test_theta_rz_solver
