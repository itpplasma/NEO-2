! B-spline basis evaluation on a breakpoint grid, replacing the former
! FGSL/GSL wrapper (issue #88). Basis values follow the Cox-de Boor
! recurrence and derivatives the standard knot-difference formula, both
! as given in DLMF section "Splines" / de Boor's published algorithm
! descriptions (clean-room implementation, no GSL code).
!
! The breakpoint grid xd(1:nbreak) is augmented to a knot vector with
! order-fold end knots, giving nbreak + order - 2 basis functions, the
! same convention as gsl_bspline_knots. Beyond the last breakpoint the
! basis is continued by a Taylor expansion around xd(nbreak), as the old
! wrapper did.
module bspline_routines_mod

    use nrtype, only : dp

    implicit none

    public init_bspline, set_bspline_knots, bspline_eval, bspline_deriv_eval
    public taylorExpansion

    integer, private :: spline_order = 0
    integer, private :: n_break = 0
    integer, private :: n_basis = 0
    real(dp), dimension(:), allocatable, private :: knots_full
    real(dp), private :: xd_end = 0.0d0
    real(dp), dimension(:), allocatable, private :: by_end
    real(dp), dimension(:, :), allocatable, private :: dby_end

    logical :: taylorExpansion = .true.

contains

    subroutine init_bspline(order, nbreak)
        integer, intent(in) :: order, nbreak

        spline_order = order
        n_break = nbreak
        n_basis = nbreak + order - 2

        if (allocated(knots_full)) deallocate(knots_full)
        allocate(knots_full(n_basis + order))
    end subroutine init_bspline

    subroutine set_bspline_knots(xd)
        real(dp), dimension(:), intent(in) :: xd

        integer :: i, k

        knots_full(1:spline_order) = xd(1)
        do i = 2, n_break - 1
            knots_full(spline_order + i - 1) = xd(i)
        end do
        knots_full(n_basis + 1:n_basis + spline_order) = xd(n_break)
        xd_end = xd(n_break)

        if (allocated(by_end)) deallocate(by_end)
        allocate(by_end(n_basis))
        call eval_basis_deriv(xd_end, 0, by_end)

        if (allocated(dby_end)) deallocate(dby_end)
        allocate(dby_end(0:spline_order, n_basis))
        do k = 0, spline_order
            call eval_basis_deriv(xd_end, k, dby_end(k, :))
        end do
    end subroutine set_bspline_knots

    recursive subroutine bspline_eval(x, by_loc)
        real(dp), intent(in) :: x
        real(dp), dimension(:), intent(out) :: by_loc

        integer :: k, m
        real(dp) :: gam

        if (x .le. xd_end) then
            call eval_basis_deriv(x, 0, by_loc)
        elseif (taylorExpansion) then
            by_loc = by_end
            do m = 1, n_basis
                gam = 1d0
                do k = 1, spline_order - 1
                    gam = k*gam
                    by_loc(m) = by_loc(m) + dby_end(k, m)/gam*(x - xd_end)**k
                end do
            end do
        else
            by_loc = 0.0d0
        end if
    end subroutine bspline_eval

    recursive subroutine bspline_deriv_eval(x, nder, dby_loc)
        real(dp), intent(in) :: x
        integer, intent(in) :: nder
        real(dp), dimension(:), intent(out) :: dby_loc

        integer :: k, m
        real(dp) :: gam

        if (x .le. xd_end) then
            call eval_basis_deriv(x, nder, dby_loc)
        elseif (taylorExpansion) then
            dby_loc = dby_end(nder, :)
            do m = 1, n_basis
                gam = 1d0
                do k = 1, spline_order - nder - 1
                    gam = k*gam
                    dby_loc(m) = dby_loc(m) + dby_end(k + nder, m)/gam &
                                 *(x - xd_end)**k
                end do
            end do
        else
            dby_loc = 0.0d0
        end if
    end subroutine bspline_deriv_eval

    ! nder-th derivative of all n_basis basis functions at x. The order-1
    ! indicator of the knot span containing x is raised to the full order
    ! by the Cox-de Boor recurrence; the last nder raising steps use the
    ! derivative recurrence instead.
    pure subroutine eval_basis_deriv(x, nder, b)
        real(dp), intent(in) :: x
        integer, intent(in) :: nder
        real(dp), dimension(:), intent(out) :: b

        real(dp), dimension(n_basis + spline_order - 1) :: work
        integer :: i, j, span
        real(dp) :: dl, dr, term

        b = 0.0d0
        if (nder .ge. spline_order) return

        span = find_span(x)
        work = 0.0d0
        work(span) = 1.0d0

        do j = 2, spline_order - nder
            do i = max(1, span - j + 1), span
                term = 0.0d0
                dl = knots_full(i + j - 1) - knots_full(i)
                if (dl .gt. 0.0d0) term = (x - knots_full(i))/dl*work(i)
                dr = knots_full(i + j) - knots_full(i + 1)
                if (dr .gt. 0.0d0) term = term &
                    + (knots_full(i + j) - x)/dr*work(i + 1)
                work(i) = term
            end do
        end do

        do j = spline_order - nder + 1, spline_order
            do i = max(1, span - j + 1), span
                term = 0.0d0
                dl = knots_full(i + j - 1) - knots_full(i)
                if (dl .gt. 0.0d0) term = work(i)/dl
                dr = knots_full(i + j) - knots_full(i + 1)
                if (dr .gt. 0.0d0) term = term - work(i + 1)/dr
                work(i) = real(j - 1, dp)*term
            end do
        end do

        b(1:n_basis) = work(1:n_basis)
    end subroutine eval_basis_deriv

    ! Index of the knot interval [t_span, t_span+1) containing x, clamped
    ! so that x at or beyond the last breakpoint lands in the last
    ! nonempty interval.
    pure function find_span(x) result(span)
        real(dp), intent(in) :: x
        integer :: span

        span = spline_order
        do while (span .lt. n_basis .and. x .ge. knots_full(span + 1))
            span = span + 1
        end do
    end function find_span

end module bspline_routines_mod
