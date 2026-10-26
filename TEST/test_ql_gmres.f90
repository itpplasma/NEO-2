program test_ql_gmres
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use ql_gmres_mod, only: gmres_solve
    implicit none
    integer, parameter :: n = 24
    complex(dp) :: a(n, n), exact(n), b(n), x(n), w(n)
    real(dp) :: residual
    integer :: i, info, steps, ncall

    ! B=I-A with a non-normal, complex fixed-point map whose diagonal has
    ! |lambda(A)|>1. The independently constructed exact solution is known.
    a = 0
    do i = 1, n
        a(i, i) = cmplx(3._dp+i/real(n, dp), 0.2_dp, dp)
        if (i < n) a(i, i+1) = (0.4_dp, -0.3_dp)
        exact(i) = cmplx(1._dp/i, (-1._dp)**i, dp)
    end do
    b = matmul(a, exact)
    call gmres_solve(n, b, x, apply, dot, 1.e-12_dp, 200, 8, info, steps, residual)
    if (info /= 0) error stop 'FAIL: non-normal solve'
    if (maxval(abs(x-exact)) > 1.e-10_dp) error stop 'FAIL: solution'
    w = b - matmul(a, x)
    if (sqrt(real(dot(n, w, w), dp)/real(dot(n, b, b), dp)) > 1.e-12_dp) &
        error stop 'FAIL: independent residual'

    call gmres_solve(n, b, x, apply, dot, 1.e-14_dp, 1, 1, info, steps, residual)
    if (info /= 1 .or. residual <= 1.e-14_dp) error stop 'FAIL: exhaustion'

    b = 0
    call gmres_solve(n, b, x, apply, dot, 1.e-12_dp, 10, 4, info, steps, residual)
    if (info /= 0 .or. any(x /= 0) .or. steps /= 0) error stop 'FAIL: zero source'

    ! Rounding-like noise of 1e-9 relative that changes with every application:
    ! the solve must stop at the noise level with the best iterate (info=3).
    b = matmul(a, exact)
    ncall = 0
    call gmres_solve(n, b, x, apply_noisy, dot, 1.e-14_dp, 2000, 8, info, steps, &
                     residual)
    if (info /= 3) error stop 'FAIL: stagnation not detected'
    if (steps >= 2000) error stop 'FAIL: stagnation detected too late'
    if (residual > 1.e-7_dp) error stop 'FAIL: stagnated residual'
    if (maxval(abs(x-exact)) > 1.e-6_dp) error stop 'FAIL: stagnated solution'

    b(1) = ieee_value(0._dp, ieee_quiet_nan)
    call gmres_solve(n, b, x, apply, dot, 1.e-12_dp, 10, 4, info, steps, residual)
    if (info == 0) error stop 'FAIL: nonfinite accepted'
    print *, 'All tests passed!'
contains
    subroutine apply(n, x, y)
        integer, intent(in) :: n
        complex(dp), intent(in) :: x(n)
        complex(dp), intent(out) :: y(n)
        y = matmul(a, x)
    end subroutine
    subroutine apply_noisy(n, x, y)
        integer, intent(in) :: n
        complex(dp), intent(in) :: x(n)
        complex(dp), intent(out) :: y(n)
        integer :: k
        ncall = ncall + 1
        y = matmul(a, x)
        do k = 1, n
            y(k) = y(k) + 1.e-9_dp*norm2(abs(x))*cmplx(sin(1.7_dp*k*ncall), &
                                                       cos(2.3_dp*k+ncall), dp)
        end do
    end subroutine
    function dot(n, x, y) result(z)
        integer, intent(in) :: n
        complex(dp), intent(in) :: x(n), y(n)
        complex(dp) :: z
        z = sum(conjg(x)*y)
    end function
end program
