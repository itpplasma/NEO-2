! fixed_point_gmres_mod on x = c + M x where the plain iteration diverges
! (an eigenvalue of M is 1.5) but I - M is regular; checked against the
! LAPACK solution of (I - M) x = c.
program test_fixed_point_gmres
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fixed_point_gmres_mod, only: gmres_solve_real
    implicit none
    integer, parameter :: n = 60
    real(dp) :: m(n, n), a(n, n), c(n), x(n), xref(n), relres
    integer :: i, j, info, its, ipiv(n)
    do j = 1, n
        do i = 1, n
            m(i, j) = 0.3_dp*sin(real(3*i + 7*j, dp))/sqrt(real(n, dp))
        end do
        c(j) = cos(real(j, dp))
    end do
    m(1, 1) = 1.5_dp
    a = -m
    do i = 1, n
        a(i, i) = a(i, i) + 1
    end do
    xref = c
    call dgesv(n, 1, a, n, ipiv, xref, n, info)
    if (info /= 0) error stop 'FAIL: reference solve'
    call gmres_solve_real(n, c, x, apply, 1e-12_dp, 500, 20, info, its, relres)
    print *, 'info, iterations, relres, error:', info, its, relres, maxval(abs(x - xref))
    if (info /= 0 .or. maxval(abs(x - xref)) > 1e-9_dp*maxval(abs(xref))) then
        print *, 'FAIL'
        error stop 1
    end if
    print *, 'All tests passed!'
contains
    subroutine apply(nn, v, w)
        integer, intent(in) :: nn
        real(dp), intent(in) :: v(nn)
        real(dp), intent(out) :: w(nn)
        w = v - matmul(m, v)
    end subroutine apply
end program test_fixed_point_gmres
