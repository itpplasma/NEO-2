! Matrix-free complex GMRES for the LU-preconditioned QL equation.
! As in PAR's fixed_point_gmres_mod, acceptance requires a fresh residual.
module ql_gmres_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: gmres_solve, ql_gmres, gmres_restart, gmres_maxiter, gmres_tol

    logical :: ql_gmres = .false.
    integer :: gmres_restart = 100, gmres_maxiter = 3000
    real(dp) :: gmres_tol = 1.e-7_dp

    abstract interface
        subroutine matvec(n, x, y)
            import dp
            integer, intent(in) :: n
            complex(dp), intent(in) :: x(n)
            complex(dp), intent(out) :: y(n)
        end subroutine
        function inner_product(n, x, y) result(z)
            import dp
            integer, intent(in) :: n
            complex(dp), intent(in) :: x(n), y(n)
            complex(dp) :: z
        end function
    end interface
contains
    ! Solve B x = b from x=0. dot must sum over every participating species.
    ! info=0: true relative 2-norm residual <= eps; 1: exhausted; 2: breakdown.
    subroutine gmres_solve(n, b, x, apply, dot, eps, maxiter, restart, &
                           info, iterations, relres)
        integer, intent(in) :: n, maxiter, restart
        complex(dp), intent(in) :: b(n)
        complex(dp), intent(out) :: x(n)
        procedure(matvec) :: apply
        procedure(inner_product) :: dot
        real(dp), intent(in) :: eps
        integer, intent(out) :: info, iterations
        real(dp), intent(out) :: relres
        complex(dp), allocatable :: v(:, :), h(:, :), w(:), g(:), y(:), sn(:)
        real(dp), allocatable :: cs(:)
        real(dp) :: bnorm, beta, wnorm, hnorm
        complex(dp) :: z, rot, rdiag
        integer :: m, i, j, pass
        logical :: breakdown

        x = 0
        iterations = 0
        info = 2
        relres = huge(1._dp)
        if (restart < 1 .or. maxiter < 1 .or. eps <= 0) return
        m = min(restart, maxiter)
        allocate(v(n, m+1), h(m+1, m), w(n), g(m+1), y(m), sn(m), cs(m))
        bnorm = sqrt(real(dot(n, b, b), dp))
        if (.not. ieee_is_finite(bnorm)) return
        if (bnorm == 0) then
            relres = 0
            info = 0
            return
        end if
        w = b
        do
            beta = sqrt(real(dot(n, w, w), dp))
            relres = beta / bnorm
            print *, 'QL GMRES true residual:', iterations, relres
            if (.not. ieee_is_finite(relres)) return
            if (relres <= eps) then
                info = 0
                return
            end if
            if (iterations >= maxiter) then
                info = 1
                return
            end if
            v(:, 1) = w / beta
            h = 0
            g = 0
            g(1) = beta
            do j = 1, min(m, maxiter-iterations)
                call apply(n, v(:, j), w)
                iterations = iterations + 1
                wnorm = sqrt(real(dot(n, w, w), dp))
                if (.not. ieee_is_finite(wnorm)) return
                ! Reorthogonalize: the high-collisionality map is strongly non-normal.
                do pass = 1, 2
                    do i = 1, j
                        z = dot(n, v(:, i), w)
                        h(i, j) = h(i, j) + z
                        w = w - z*v(:, i)
                    end do
                end do
                hnorm = sqrt(real(dot(n, w, w), dp))
                if (.not. ieee_is_finite(hnorm)) return
                h(j+1, j) = hnorm
                breakdown = hnorm <= epsilon(1._dp)*wnorm
                if (.not. breakdown) v(:, j+1) = w / hnorm
                do i = 1, j-1
                    rot = cs(i)*h(i, j) + sn(i)*h(i+1, j)
                    h(i+1, j) = -conjg(sn(i))*h(i, j) + cs(i)*h(i+1, j)
                    h(i, j) = rot
                end do
                call zlartg(h(j, j), h(j+1, j), cs(j), sn(j), rdiag)
                h(j, j) = rdiag
                h(j+1, j) = 0
                g(j+1) = -conjg(sn(j))*g(j)
                g(j) = cs(j)*g(j)
                if (abs(g(j+1)) <= eps*bnorm .or. breakdown .or. &
                    j == m .or. iterations == maxiter) then
                    y(1:j) = g(1:j)
                    do i = j, 1, -1
                        if (abs(h(i, i)) <= tiny(1._dp)) return
                        y(i) = (y(i)-sum(h(i, i+1:j)*y(i+1:j))) / h(i, i)
                    end do
                    x = x + matmul(v(:, 1:j), y(1:j))
                    call apply(n, x, w)
                    w = b - w
                    exit
                end if
            end do
        end do
    end subroutine gmres_solve
end module ql_gmres_mod
