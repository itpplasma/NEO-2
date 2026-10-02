module fieldline_step_mod
    ! One output step of the field-line ODE of rhs_kin with mixed
    ! absolute/relative error control.
    !
    ! State vector (see rhs_kin):
    !   y(1:2)   trajectory: R, Z (cylindrical) or theta, const (Boozer)
    !   y(3:5)   grad psi, transported along the line
    !   y(6:10)  running integrals of f_i along the line (f_6 = 1/(B^phi B))
    !   y(11:14) running integrals of y(7:10)
    ! Only y(1:5) enter the right-hand side; y(6:14) are quadratures that do
    ! not feed back into the trajectory.
    !
    ! A purely relative error test, |err_i| <= rtol*(|y_i| + |h*f_i|), cannot
    ! be met by y(11:14) at the start of a line, where y(11:14) = 0 and their
    ! derivatives y(7:10) = 0 as well, so the step size underflows. Each
    ! component therefore gets an absolute tolerance atol_i = rtol*S_i with S_i
    ! its physical scale over the output interval dphi:
    !   y(1:2)   S = length_scale (R in cm, or 1 rad for the Boozer angle);
    !            Z and theta cross zero, R sets the scale of both.
    !   y(3:5)   covariant grad psi (dpsi/dR, dpsi/dphi, dpsi/dZ); the
    !            components rotate and cross zero. S = |grad psi| for y(3),
    !            y(5) and R |grad psi| for y(4) (constant in Boozer form).
    !   y(6:10)  S = |dphi|*|f_i|, the increment over the interval.
    !            f_8 = (r . grad psi) f_6 changes sign, so its scale is the
    !            Cauchy-Schwarz bound |r| |grad psi| |f_6|.
    !   y(11:14) S = |dphi|*|y_{i-4}| + dphi**2/2*S'_{i-4}, the increment of
    !            the second integral (S' the integrand scale of y_{i-4}).
    ! These characteristic magnitudes control the local error estimator.
    ! They do not bound accumulated error relative to a sign-cancelling integral.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use odeint_allroutines_sub, only: odeint_allroutines_checked
    implicit none
    private

    public :: fieldline_step, fieldline_abs_tol

    integer, parameter, public :: fieldline_ndim = 14
    real(dp), parameter, public :: fieldline_rtol = 1.0e-10_dp

    abstract interface
        subroutine rhs_interface(x, y, dydx)
            import :: dp
            real(dp), intent(in)  :: x
            real(dp), intent(in)  :: y(:)
            real(dp), intent(out) :: dydx(:)
        end subroutine rhs_interface
    end interface

contains

    ! Integrate y from phi to phi + dphi. On success ierr = 0 and phi is
    ! advanced. On failure (step size underflow, step limit, non-finite state)
    ! ierr = 1 and phi and y are left at the start of the interval.
    subroutine fieldline_step(phi, dphi, y, rhs, length_scale, ierr)
        real(dp), intent(inout) :: phi
        real(dp), intent(in)    :: dphi, length_scale
        real(dp), intent(inout), contiguous :: y(:)
        procedure(rhs_interface) :: rhs
        integer, intent(out)    :: ierr

        real(dp) :: dydphi(fieldline_ndim), atol(fieldline_ndim), &
                    y_start(fieldline_ndim)

        if (size(y) /= fieldline_ndim) error stop 'fieldline_step: size(y) /= 14'
        ierr = 1
        if (.not. all(ieee_is_finite(y))) return
        if (.not. ieee_is_finite(phi)) return
        if (.not. ieee_is_finite(dphi)) return
        if (.not. ieee_is_finite(length_scale)) return
        if (length_scale == 0.0_dp) return
        call rhs(phi, y, dydphi)
        if (.not. all(ieee_is_finite(dydphi))) return
        atol = fieldline_abs_tol(dphi, y, dydphi, length_scale, fieldline_rtol)

        y_start = y
        call odeint_allroutines_checked(y, size(y), phi, phi + dphi, fieldline_rtol, &
                                        rhs, ierr, atol=atol)
        if (ierr /= 0) then
            y = y_start
            return
        end if
        ierr = 0
        phi = phi + dphi
    end subroutine fieldline_step

    pure function fieldline_abs_tol(dphi, y, dydphi, length_scale, rtol) result(atol)
        real(dp), intent(in) :: dphi, y(:), dydphi(:), length_scale, rtol
        real(dp) :: atol(fieldline_ndim)

        real(dp) :: f_scale(6:10), grad_psi, h, r_scale

        h = abs(dphi)
        r_scale = abs(length_scale)
        ! |grad psi| from the covariant components (dpsi/dR, dpsi/dphi, dpsi/dZ)
        grad_psi = sqrt(y(3)**2 + (y(4)/r_scale)**2 + y(5)**2)

        f_scale = abs(dydphi(6:10))
        f_scale(8) = max(f_scale(8), norm2(y(1:2))*grad_psi*abs(dydphi(6)))

        atol(1:2) = rtol*r_scale
        atol(3) = rtol*grad_psi
        atol(4) = rtol*r_scale*grad_psi
        atol(5) = rtol*grad_psi
        atol(6:10) = rtol*h*f_scale
        atol(11:14) = rtol*(h*abs(y(7:10)) + 0.5_dp*h**2*f_scale(7:10))
    end function fieldline_abs_tol

end module fieldline_step_mod
