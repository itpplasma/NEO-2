! Matrix-free reference solve for L = P - D W, with optional exact constraints.
! Vectors are global/concatenated; this is NOT an MPI-distributed solver.
module moment_schur_solver
    use, intrinsic :: iso_fortran_env, only: real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    integer, parameter :: dp = real64
    integer, parameter, public :: schur_ok = 0, schur_invalid = 1, schur_callback = 2
    integer, parameter, public :: schur_nonfinite = 3, schur_breakdown = 4, schur_limit = 5
    integer, parameter, public :: schur_original_residual = 6, schur_constraint = 7
    integer, parameter, public :: schur_allocation = 8, schur_reconstruction = 9

    abstract interface
        subroutine action(x, y, info)
            import dp
            complex(dp), intent(in) :: x(:)
            complex(dp), intent(out) :: y(:)
            integer, intent(out) :: info
        end subroutine
    end interface

    type, public :: moment_schur_operator_t
        integer :: nfield = 0, nmoment = 0, nconstraint = 0
        procedure(action), pointer, nopass :: inverse_p => null()
        procedure(action), pointer, nopass :: feedback => null() ! D: moments -> field
        procedure(action), pointer, nopass :: moments => null()  ! W: field -> moments
        procedure(action), pointer, nopass :: original => null() ! ORIGINAL L, not P
        procedure(action), pointer, nopass :: force => null()    ! K: multipliers -> field
        procedure(action), pointer, nopass :: constraints => null() ! B: field -> constraints
    end type

    type, public :: moment_schur_options_t
        integer :: restart = 30, maxiter = 300
        real(dp) :: rtol = 1.0e-10_dp, atol = 0.0_dp
        ! Absolute tolerance for normalized constraint rows; set explicitly for physical units.
        real(dp) :: constraint_atol = 1.0e-12_dp
    end type

    type, public :: moment_schur_report_t
        integer :: status = schur_invalid, iterations = 0, applications = 0, inverse_calls = 0
        real(dp) :: reduced_residual = huge(1.0_dp), original_residual = huge(1.0_dp)
        real(dp) :: bordered_residual = huge(1.0_dp), constraint_residual = huge(1.0_dp)
    end type
    public :: solve_moment_schur

    interface
        subroutine zlartg(f, g, c, s, r)
            import dp
            complex(dp), intent(in) :: f, g
            real(dp), intent(out) :: c
            complex(dp), intent(out) :: s, r
        end subroutine
    end interface

contains
    ! Solve [L K; B 0] [f;lambda] = [q;g] by eliminating the full field.
    ! Success ALSO requires L f = q: a nonzero K lambda is never silently
    ! accepted as a solution of the original kinetic problem.
    subroutine solve_moment_schur(op, q, g, f, lambda, options, report)
        type(moment_schur_operator_t), intent(in) :: op
        complex(dp), intent(in) :: q(:), g(:)
        complex(dp), intent(out) :: f(:), lambda(:)
        type(moment_schur_options_t), intent(in) :: options
        type(moment_schur_report_t), intent(out) :: report
        complex(dp), allocatable :: f0(:), rhs(:), z(:), source(:), correction(:), work(:), bg(:)
        integer :: n, r, c, k, info, stat, solve_status
        real(dp) :: target

        report = moment_schur_report_t()
        f = (0.0_dp, 0.0_dp); lambda = (0.0_dp, 0.0_dp)
        n = op%nfield; r = op%nmoment; c = op%nconstraint
        if (n < 1 .or. r < 0 .or. c < 0) return
        if (r > huge(0) - c) return
        k = r + c
        if (size(q) /= n .or. size(f) /= n .or. size(g) /= c .or. size(lambda) /= c) return
        if (.not. associated(op%inverse_p) .or. .not. associated(op%original)) return
        if (r > 0) then
            if (.not. associated(op%feedback) .or. .not. associated(op%moments)) return
        end if
        if (c > 0) then
            if (.not. associated(op%force) .or. .not. associated(op%constraints)) return
        end if
        if (options%restart < 1 .or. options%maxiter < 1) return
        if (.not. ieee_is_finite(options%rtol) .or. .not. ieee_is_finite(options%atol)) return
        if (.not. ieee_is_finite(options%constraint_atol)) return
        if (options%rtol < 0 .or. options%atol < 0 .or. options%constraint_atol < 0) return
        if (options%rtol <= 0 .and. options%atol <= 0) return
        if (.not. finite(q) .or. .not. finite(g)) then
            report%status = schur_nonfinite
            return
        end if
        allocate (f0(n), rhs(k), z(k), source(n), correction(n), work(n), bg(c), stat=stat)
        if (stat /= 0) then
            report%status = schur_allocation
            return
        end if
        call invoke(op%inverse_p, q, f0, info)
        report%inverse_calls = report%inverse_calls + 1
        if (info /= 0) then
            report%status = info
            return
        end if
        if (r > 0) then
            call invoke(op%moments, f0, rhs(1:r), info)
            if (info /= 0) then
                report%status = info
                return
            end if
        end if
        if (c > 0) then
            call invoke(op%constraints, f0, bg, info)
            if (info /= 0) then
                report%status = info
                return
            end if
            rhs(r + 1:k) = g - bg
        end if
        z = (0.0_dp, 0.0_dp)
        if (k > 0) then
            call gmres(op, rhs, z, options, report, solve_status)
        else
            solve_status = schur_ok
            report%reduced_residual = 0.0_dp
        end if
        if (solve_status == schur_callback .or. solve_status == schur_nonfinite &
            .or. solve_status == schur_allocation) then
            report%status = solve_status
            return
        end if
        ! Recover all fields even after a finite nonconverged GMRES iteration.
        call field_correction(op, z, correction, report, info)
        if (info /= 0) then
            report%status = info
            return
        end if
        f = f0 + correction
        lambda = z(r + 1:k)
        if (.not. finite(f)) then
            report%status = schur_nonfinite
            return
        end if
        call invoke(op%original, f, work, info)
        if (info /= 0) then
            report%status = info
            return
        end if
        work = work - q
        report%original_residual = vnorm(work)
        if (c > 0) then
            call invoke(op%force, lambda, source, info)
            if (info /= 0) then
                report%status = info
                return
            end if
            work = work + source
            call invoke(op%constraints, f, bg, info)
            if (info /= 0) then
                report%status = info
                return
            end if
        end if
        report%bordered_residual = vnorm(work)
        report%constraint_residual = vnorm(bg - g)
        target = options%atol + options%rtol*vnorm(q)
        report%status = solve_status
        if (solve_status /= schur_ok) return
        if (.not. ieee_is_finite(target) .or. .not. ieee_is_finite(report%original_residual) &
            .or. .not. ieee_is_finite(report%bordered_residual) &
            .or. .not. ieee_is_finite(report%constraint_residual)) then
            report%status = schur_nonfinite
        else if (report%bordered_residual > target) then
            report%status = schur_reconstruction
        else if (report%constraint_residual > options%constraint_atol) then
            report%status = schur_constraint
        else if (report%original_residual > target) then
            report%status = schur_original_residual
        end if

    end subroutine solve_moment_schur

    subroutine field_correction(op, x, y, report, info)
        type(moment_schur_operator_t), intent(in) :: op
        complex(dp), intent(in) :: x(:)
        complex(dp), intent(out) :: y(:)
        type(moment_schur_report_t), intent(inout) :: report
        integer, intent(out) :: info
        complex(dp) :: source(op%nfield), work(op%nfield)
        integer :: r, k
        r = op%nmoment; k = r + op%nconstraint
        source = (0.0_dp, 0.0_dp)
        info = schur_ok
        if (r > 0) then
            call invoke(op%feedback, x(1:r), source, info)
            if (info /= 0) return
        end if
        if (op%nconstraint > 0) then
            call invoke(op%force, x(r + 1:k), work, info)
            if (info /= 0) return
            source = source - work
        end if
        call invoke(op%inverse_p, source, y, info)
        report%inverse_calls = report%inverse_calls + 1
    end subroutine field_correction

    subroutine reduced_action(op, x, y, report, info)
        type(moment_schur_operator_t), intent(in) :: op
        complex(dp), intent(in) :: x(:)
        complex(dp), intent(out) :: y(:)
        type(moment_schur_report_t), intent(inout) :: report
        integer, intent(out) :: info
        complex(dp) :: field(op%nfield)
        integer :: r, k
        r = op%nmoment; k = r + op%nconstraint
        call field_correction(op, x, field, report, info)
        report%applications = report%applications + 1
        if (info /= 0) return
        if (r > 0) then
            call invoke(op%moments, field, y(1:r), info)
            if (info /= 0) return
            y(1:r) = x(1:r) - y(1:r)
        end if
        if (op%nconstraint > 0) call invoke(op%constraints, field, y(r + 1:k), info)
        if (info == schur_ok) then
            if (.not. finite(y)) info = schur_nonfinite
        end if
    end subroutine reduced_action

    ! Restarted complex GMRES; two-pass MGS, LAPACK stable complex Givens,
    ! explicit residual after every cycle, and no division by a unit Ritz gap.
    subroutine gmres(op, b, x, opts, report, status)
        type(moment_schur_operator_t), intent(in) :: op
        complex(dp), intent(in) :: b(:)
        complex(dp), intent(inout) :: x(:)
        type(moment_schur_options_t), intent(in) :: opts
        type(moment_schur_report_t), intent(inout) :: report
        integer, intent(out) :: status
        complex(dp), allocatable :: v(:, :), h(:, :), s(:), gg(:), yy(:), w(:), residual(:)
        real(dp), allocatable :: cs(:)
        complex(dp) :: a, t, rot
        integer :: n, m, j, i, pass, info, stat, count, used
        real(dp) :: beta, initial_norm, remaining_norm, target
        logical :: closed

        n = size(b); m = min(opts%restart, n, opts%maxiter)
        allocate (v(n, m + 1), h(m + 1, m), s(m), gg(m + 1), yy(m), w(n), residual(n), cs(m), stat=stat)
        if (stat /= 0) then
            status = schur_allocation
            return
        end if
        target = opts%atol + opts%rtol*vnorm(b)
        if (.not. ieee_is_finite(target)) then
            status = schur_nonfinite
            return
        end if
        count = 0
        call actual_residual(info)
        if (info /= 0) then
            status = info
            return
        end if
        do while (report%reduced_residual > target .and. count < opts%maxiter)
            beta = report%reduced_residual
            v = (0.0_dp, 0.0_dp); h = (0.0_dp, 0.0_dp); gg = (0.0_dp, 0.0_dp)
            v(:, 1) = residual/beta
            gg(1) = cmplx(beta, 0.0_dp, dp)
            closed = .false.
            used = min(m, opts%maxiter - count)
            do j = 1, used
                call reduced_action(op, v(:, j), w, report, info)
                if (info /= 0) then
                    status = info
                    return
                end if
                count = count + 1
                report%iterations = count
                initial_norm = vnorm(w)
                do pass = 1, 2
                    do i = 1, j
                        a = dot_product(v(:, i), w)
                        h(i, j) = h(i, j) + a
                        w = w - a*v(:, i)
                    end do
                end do
                remaining_norm = vnorm(w)
                if (.not. ieee_is_finite(initial_norm) .or. .not. ieee_is_finite(remaining_norm)) then
                    status = schur_nonfinite
                    return
                end if
                closed = remaining_norm <= 64*epsilon(1.0_dp)*max(initial_norm, tiny(1.0_dp))
                h(j + 1, j) = cmplx(remaining_norm, 0.0_dp, dp)
                if (.not. closed) v(:, j + 1) = w/remaining_norm
                do i = 1, j - 1
                    t = cs(i)*h(i, j) + s(i)*h(i + 1, j)
                    h(i + 1, j) = -conjg(s(i))*h(i, j) + cs(i)*h(i + 1, j)
                    h(i, j) = t
                end do
                call zlartg(h(j, j), h(j + 1, j), cs(j), s(j), rot)
                h(j, j) = rot
                h(j + 1, j) = (0.0_dp, 0.0_dp)
                gg(j + 1) = -conjg(s(j))*gg(j)
                gg(j) = cs(j)*gg(j)
                if (abs(gg(j + 1)) <= target .or. closed) then
                    used = j
                    exit
                end if
            end do
            yy(1:used) = gg(1:used)
            do i = used, 1, -1
                if (abs(h(i, i)) <= tiny(1.0_dp)) then
                    status = schur_breakdown
                    return
                end if
                yy(i) = (yy(i) - sum(h(i, i + 1:used)*yy(i + 1:used)))/h(i, i)
            end do
            x = x + matmul(v(:, 1:used), yy(1:used))
            if (.not. finite(x)) then
                status = schur_nonfinite
                return
            end if
            call actual_residual(info)
            if (info /= 0) then
                status = info
                return
            end if
            if (closed .and. report%reduced_residual > target) then
                status = schur_breakdown
                return
            end if
        end do
        status = schur_limit
        if (report%reduced_residual <= target) status = schur_ok

    contains
        subroutine actual_residual(ierr)
            integer, intent(out) :: ierr
            call reduced_action(op, x, w, report, ierr)
            if (ierr /= 0) return
            residual = b - w
            report%reduced_residual = vnorm(residual)
            if (.not. ieee_is_finite(report%reduced_residual)) ierr = schur_nonfinite
        end subroutine
    end subroutine gmres

    subroutine invoke(fn, x, y, info)
        procedure(action) :: fn
        complex(dp), intent(in) :: x(:)
        complex(dp), intent(out) :: y(:)
        integer, intent(out) :: info
        integer :: ierr
        if (.not. finite(x)) then
            info = schur_nonfinite
            return
        end if
        call fn(x, y, ierr)
        info = schur_ok
        if (ierr /= 0) then
            info = schur_callback
        else if (.not. finite(y)) then
            info = schur_nonfinite
        end if
    end subroutine

    logical function finite(x)
        complex(dp), intent(in) :: x(:)
        finite = all(ieee_is_finite(real(x, dp))) .and. all(ieee_is_finite(aimag(x)))
    end function

    real(dp) function vnorm(x) result(value)
        complex(dp), intent(in) :: x(:)
        real(dp) :: scale
        value = 0.0_dp
        if (size(x) == 0) return
        scale = maxval(abs(x))
        if (scale > 0) value = scale*sqrt(sum(abs(x/scale)**2))
    end function
end module moment_schur_solver
