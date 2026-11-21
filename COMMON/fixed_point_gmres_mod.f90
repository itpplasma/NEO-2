! Restarted real GMRES for a fixed-point equation x = c + M x written as
! (I - M) x = c, matrix-free. Used when the integral-part iteration of
! NEO-2-PAR stalls; the real counterpart of NEO-2-QL's ql_gmres_mod.
! Acceptance requires a fresh true residual.
module fixed_point_gmres_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: gmres_solve_real

    abstract interface
        subroutine matvec(n, x, y)
            import dp
            integer, intent(in) :: n
            real(dp), intent(in) :: x(n)
            real(dp), intent(out) :: y(n)
        end subroutine
    end interface
contains
    ! Solve B x = b from x = 0. info = 0: true relative 2-norm residual <= eps;
    ! 1: maxiter exhausted; 2: breakdown or non-finite; 3: stagnation (true
    ! residual did not halve over five restarts); x is then the iterate with
    ! the smallest true residual relres.
    subroutine gmres_solve_real(n, b, x, apply, eps, maxiter, restart, info, &
            iterations, relres)
        integer, intent(in) :: n, maxiter, restart
        real(dp), intent(in) :: b(n), eps
        real(dp), intent(out) :: x(n)
        procedure(matvec) :: apply
        integer, intent(out) :: info, iterations
        real(dp), intent(out) :: relres
        real(dp), allocatable :: v(:, :), h(:, :), w(:), g(:), y(:), sn(:), cs(:), xbest(:)
        real(dp) :: bnorm, beta, wnorm, hnorm, best, z, rot, rdiag
        integer :: m, i, j, pass, stalled
        logical :: breakdown

        x = 0
        iterations = 0
        info = 2
        relres = huge(1._dp)
        if (restart < 1 .or. maxiter < 1 .or. eps <= 0) return
        m = min(restart, maxiter)
        allocate (v(n, m + 1), h(m + 1, m), w(n), g(m + 1), y(m), sn(m), cs(m), xbest(n))
        bnorm = norm2(b)
        if (.not. ieee_is_finite(bnorm)) return
        if (bnorm == 0) then
            relres = 0
            info = 0
            return
        end if
        w = b
        best = huge(1._dp)
        stalled = 0
        do
            beta = norm2(w)
            relres = beta/bnorm
            if (.not. ieee_is_finite(relres)) return
            if (relres <= eps) then
                info = 0
                return
            end if
            if (relres < 0.5_dp*best) then
                stalled = 0
            else
                stalled = stalled + 1
            end if
            if (relres < best) then
                best = relres
                xbest = x
            end if
            if (stalled >= 5) then
                x = xbest
                relres = best
                info = 3
                return
            end if
            if (iterations >= maxiter) then
                info = 1
                return
            end if
            v(:, 1) = w/beta
            h = 0
            g = 0
            g(1) = beta
            do j = 1, min(m, maxiter - iterations)
                call apply(n, v(:, j), w)
                iterations = iterations + 1
                wnorm = norm2(w)
                if (.not. ieee_is_finite(wnorm)) return
                do pass = 1, 2
                    do i = 1, j
                        z = dot_product(v(:, i), w)
                        h(i, j) = h(i, j) + z
                        w = w - z*v(:, i)
                    end do
                end do
                hnorm = norm2(w)
                h(j + 1, j) = hnorm
                breakdown = hnorm <= epsilon(1._dp)*wnorm
                if (.not. breakdown) v(:, j + 1) = w/hnorm
                do i = 1, j - 1
                    rot = cs(i)*h(i, j) + sn(i)*h(i + 1, j)
                    h(i + 1, j) = -sn(i)*h(i, j) + cs(i)*h(i + 1, j)
                    h(i, j) = rot
                end do
                call dlartg(h(j, j), h(j + 1, j), cs(j), sn(j), rdiag)
                h(j, j) = rdiag
                h(j + 1, j) = 0
                g(j + 1) = -sn(j)*g(j)
                g(j) = cs(j)*g(j)
                if (abs(g(j + 1)) <= eps*bnorm .or. breakdown .or. &
                    j == m .or. iterations == maxiter) then
                    y(1:j) = g(1:j)
                    do i = j, 1, -1
                        if (abs(h(i, i)) <= tiny(1._dp)) return
                        y(i) = (y(i) - sum(h(i, i + 1:j)*y(i + 1:j)))/h(i, i)
                    end do
                    x = x + matmul(v(:, 1:j), y(1:j))
                    call apply(n, x, w)
                    w = b - w
                    exit
                end if
            end do
        end do
    end subroutine gmres_solve_real
end module fixed_point_gmres_mod
