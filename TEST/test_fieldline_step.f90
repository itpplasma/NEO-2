module test_fieldline_step_model
    ! Analytic model with the layout of rhs_kin in Boozer form:
    !   theta' = iota, y(2:5) constant,
    !   f_6 = 1 + eps*cos(theta), f_7 = a f_6, f_8 = d H(phi), f_9 = b f_6,
    !   f_10 = c f_6, y(11:14)' = y(7:10).
    ! All running integrals start at zero; y(11:14) also start with zero
    ! derivative, as at the start of a field line. In rhs_kin the purely
    ! relative error test then fails at every step size until roundoff
    ! decides, which is platform dependent. The switch-on integrand f_8
    ! (H the Heaviside step, zero at phi = 0 but d at every later RK stage)
    ! makes that failure deterministic: its error estimate is a fixed
    ! fraction of h*d, while the relative tolerance at phi = 0 is 1e-30*rtol.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: &
        ieee_value, ieee_quiet_nan, ieee_positive_inf
    implicit none
    real(dp), parameter :: iota = 0.37_dp, eps = 0.3_dp, theta0 = 0.4_dp
    real(dp), parameter :: a = 2.5e3_dp, b = 0.9_dp, c = 0.8_dp, d = 0.7_dp
    real(dp) :: phi_nan = huge(1.0_dp)
contains
    subroutine rhs_model(phi, y, dydphi)
        real(dp), intent(in)  :: phi
        real(dp), intent(in)  :: y(:)
        real(dp), intent(out) :: dydphi(:)
        real(dp) :: f6

        f6 = 1.0_dp + eps*cos(y(1))
        dydphi = 0.0_dp
        dydphi(1) = iota
        dydphi(6) = f6
        dydphi(7) = a*f6
        if (phi > 0.0_dp) dydphi(8) = d
        dydphi(9) = b*f6
        dydphi(10) = c*f6
        dydphi(11:14) = y(7:10)
        if (phi > phi_nan) dydphi(6) = ieee_value(1.0_dp, ieee_quiet_nan)
    end subroutine rhs_model

    function y_exact(phi) result(y)
        real(dp), intent(in) :: phi
        real(dp) :: y(14), i1, i2, theta

        theta = theta0 + iota*phi
        i1 = phi + eps/iota*(sin(theta) - sin(theta0))
        i2 = 0.5_dp*phi**2 &
             - eps/iota*((cos(theta) - cos(theta0))/iota + phi*sin(theta0))
        y = [theta, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
             i1, a*i1, d*max(phi, 0.0_dp), b*i1, c*i1, &
             a*i2, 0.5_dp*d*max(phi, 0.0_dp)**2, b*i2, c*i2]
    end function y_exact
    subroutine rhs_cylindrical(phi, y, dydphi)
        real(dp), intent(in) :: phi, y(:)
        real(dp), intent(out) :: dydphi(:)
        dydphi = 0.0_dp
        dydphi(2) = 10.0_dp*cos(phi)
        dydphi(3) = -y(5)
        dydphi(5) = y(3)
        dydphi(6) = 1.0_dp
        dydphi(7) = sqrt(200.0_dp**2 + 50.0_dp**2)
        dydphi(8) = y(1)*y(3) + y(2)*y(5)
        dydphi(9:10) = 1.0_dp
        dydphi(11:14) = y(7:10)
    end subroutine rhs_cylindrical

    function y_cylindrical_exact(phi) result(y)
        real(dp), intent(in) :: phi
        real(dp) :: y(14), i8, j8, grad
        grad = sqrt(200.0_dp**2 + 50.0_dp**2)
        i8 = 200.0_dp*(300.0_dp*sin(phi) + 20.0_dp*(1.0_dp - cos(phi)) &
                       + 10.0_dp*(0.5_dp*phi - 0.25_dp*sin(2.0_dp*phi)))
        j8 = 200.0_dp*(300.0_dp*(1.0_dp - cos(phi)) &
                       + 20.0_dp*(phi - sin(phi)) &
                       + 10.0_dp*(0.25_dp*phi**2 + 0.125_dp*(cos(2.0_dp*phi) - 1.0_dp)))
        y(1) = 300.0_dp
        y(2) = 20.0_dp + 10.0_dp*sin(phi)
        y(3) = 200.0_dp*cos(phi)
        y(4) = 15000.0_dp
        y(5) = 200.0_dp*sin(phi)
        y(6) = phi
        y(7) = grad*phi
        y(8) = i8
        y(9:10) = phi
        y(11) = 0.5_dp*grad*phi**2
        y(12) = j8
        y(13:14) = 0.5_dp*phi**2
    end function y_cylindrical_exact
end module test_fieldline_step_model

program test_fieldline_step
    use test_fieldline_step_model
    use fieldline_step_mod, only: fieldline_step, fieldline_rtol
    use odeint_allroutines_sub, only: odeint_allroutines_checked
    use mag_interface_mod, only: tokamak_fieldline_closed
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp) :: y(14), y_ref(14), phi, dphi, err
    integer :: ierr, i, nfail
    logical :: closed
    character(len=16) :: mode

    call get_command_argument(1, mode)
    if (trim(mode) == 'closure_cap') then
        ! Must stop at the expected period when the line did not close.
        closed = tokamak_fieldline_closed(2, 2, 0.2_dp, 1.2_dp)
        print *, 'FAIL: closure cap did not stop the run'
        stop
    end if

    nfail = 0
    dphi = 2.0_dp*pi/50.0_dp

    ! 0. The previous call in rk4_kin (purely relative error control) fails
    !    on the first interval of a line and must report it.
    y = y_exact(0.0_dp)
    call odeint_allroutines_checked(y, size(y), 0.0_dp, dphi, fieldline_rtol, &
                                    rhs_model, ierr)
    call check(ierr /= 0, 'relative error control reports failure')

    ! 1. Field-line model over 20 toroidal turns against the analytic solution.
    y = y_exact(0.0_dp)
    phi = 0.0_dp
    ierr = 0
    do i = 1, 1000
        call fieldline_step(phi, dphi, y, rhs_model, 1.0_dp, ierr)
        if (ierr /= 0) exit
    end do
    y_ref = y_exact(phi)
    err = maxval(abs(y - y_ref)/max(abs(y_ref), 1.0_dp))
    print *, 'max relative error after ', i - 1, ' steps: ', err
    call check(ierr == 0, 'all steps succeed from zero running integrals')
    call check(abs(phi - 1000*dphi) < 1.0e-9_dp, 'phi advanced by the steps')
    call check(err < 1.0e-10_dp, 'matches analytic solution')

    ! Backward output intervals and a zero interval have analytic endpoints.
    phi = 0.0_dp
    y = y_exact(phi)
    call fieldline_step(phi, 0.0_dp, y, rhs_model, 1.0_dp, ierr)
    y_ref = y_exact(phi)
    call check(ierr == 0 .and. phi == 0.0_dp .and. all(y == y_ref), &
               'zero interval preserves state and angle')
    do i = 1, 1000
        call fieldline_step(phi, -dphi, y, rhs_model, 1.0_dp, ierr)
        if (ierr /= 0) exit
    end do
    y_ref = y_exact(phi)
    err = maxval(abs(y - y_ref)/max(abs(y_ref), 1.0_dp))
    print *, 'backward max relative error: ', err
    call check(ierr == 0 .and. abs(phi + 1000*dphi) < 1.0e-9_dp, &
               'backward steps reach endpoint')
    call check(err < 1.0e-10_dp, 'backward steps match analytic solution')

    ! Cylindrical cm units, rotating gradients and a sign-changing quadrature.
    phi = 0.0_dp
    y = y_cylindrical_exact(phi)
    do i = 1, 200
        call fieldline_step(phi, 0.1_dp, y, rhs_cylindrical, 300.0_dp, ierr)
        if (ierr /= 0) exit
    end do
    y_ref = y_cylindrical_exact(phi)
    err = maxval(abs(y - y_ref)/max(abs(y_ref), 1.0_dp))
    print *, 'cylindrical max relative error: ', err
    call check(ierr == 0 .and. abs(phi - 20.0_dp) < 1.0e-12_dp, &
               'cylindrical steps reach endpoint')
    call check(err < 1.0e-9_dp, 'cylindrical units match analytic solution')

    ! 2. A failing step returns ierr /= 0 and leaves phi and y unchanged.
    phi_nan = 0.5_dp*dphi
    y = y_exact(0.0_dp)
    phi = 0.0_dp
    call fieldline_step(phi, dphi, y, rhs_model, 1.0_dp, ierr)
    call check(ierr /= 0, 'non-finite right-hand side reports failure')
    call check(phi == 0.0_dp .and. all(y == y_exact(0.0_dp)), &
               'failed step does not advance phi or y')
    phi_nan = huge(1.0_dp)

    y = y_exact(0.0_dp)
    y(6) = ieee_value(1.0_dp, ieee_positive_inf)
    phi = 0.0_dp
    call fieldline_step(phi, dphi, y, rhs_model, 1.0_dp, ierr)
    call check(ierr /= 0 .and. phi == 0.0_dp, 'nonfinite input is rejected')
    y = y_exact(0.0_dp)
    call fieldline_step(phi, dphi, y, rhs_model, 0.0_dp, ierr)
    call check(ierr /= 0 .and. phi == 0.0_dp, 'zero physical scale is rejected')

    ! 3. Tokamak closure: only at the expected period.
    call check(tokamak_fieldline_closed(2, 2, 0.2_dp, 0.2_dp + 1.0e-8_dp), &
               'closes at expected period')
    call check(tokamak_fieldline_closed(2, 2, pi - 1.0e-4_dp, -pi + 1.0e-4_dp), &
               'closes across theta = pi')
    call check(.not. tokamak_fieldline_closed(1, 2, 0.2_dp, 0.2_dp), &
               'no closure before expected period')
    call check(.not. tokamak_fieldline_closed(1, 2, 0.2_dp, 0.2_dp + 2.0e-3_dp), &
               'open line is not closed')

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

end program test_fieldline_step
