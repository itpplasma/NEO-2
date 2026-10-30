!> Restarted GMRES on x = b + T x with a slowly contracting, non-normal T.
!> Oracle: LAPACK solve of (I - T) x = b.
program test_fixed_point_gmres
    use fixed_point_gmres_mod, only: fixed_point_gmres_t
    implicit none

    integer, parameter :: dp = kind(1.0d0), n = 300
    real(dp) :: t(n, n), a(n, n), b(n), xref(n)
    integer :: ipiv(n), info, i, j
    integer(8) :: seed
    type(fixed_point_gmres_t) :: gm

    seed = 12345_8
    t = 0.0_dp
    do i = 1, n
        t(i, i) = -0.4_dp + 1.37_dp*real(i - 1, dp)/real(n - 1, dp)
        do j = i + 1, n
            t(i, j) = 0.3_dp*(lcg(seed) - 0.5_dp)/sqrt(real(n, dp))
        end do
        b(i) = 1.0_dp + lcg(seed)
    end do
    a = -t
    do i = 1, n
        a(i, i) = a(i, i) + 1.0_dp
    end do
    xref = b
    call dgesv(n, 1, a, n, ipiv, xref, n, info)
    if (info /= 0) error stop 'FAIL: LAPACK reference solve'

    call run(1.0e-10_dp, 1000, 30)
    if (.not. gm%converged) error stop 'FAIL: GMRES(30) did not converge'
    if (maxval(abs(gm%x - xref)) > 1.0e-8_dp*maxval(abs(xref))) &
        error stop 'FAIL: GMRES(30) disagrees with LAPACK'
    print '(a,i0)', ' GMRES(30) applications: ', gm%napply

    call run(1.0e-10_dp, 5000, 5)
    if (.not. gm%converged) error stop 'FAIL: GMRES(5) did not converge'
    if (maxval(abs(gm%x - xref)) > 1.0e-8_dp*maxval(abs(xref))) &
        error stop 'FAIL: GMRES(5) disagrees with LAPACK'
    print '(a,i0)', ' GMRES(5) applications: ', gm%napply

    call run(1.0e-14_dp, 7, 30)
    if (gm%converged .or. gm%napply /= 7) &
        error stop 'FAIL: budget exhaustion not reported exactly'

    print *, 'All tests passed!'

contains

    subroutine run(eps, maxapply, m)
        real(dp), intent(in) :: eps
        integer, intent(in) :: maxapply, m

        call gm%start(b, eps, maxapply, m)
        do while (gm%needs_apply())
            call gm%put(matmul(t, gm%vin))
        end do
    end subroutine run

    real(dp) function lcg(s)
        integer(8), intent(inout) :: s
        s = modulo(16807_8*s, 2147483647_8)
        lcg = real(s, dp)/2147483647.0_dp
    end function lcg

end program test_fixed_point_gmres
