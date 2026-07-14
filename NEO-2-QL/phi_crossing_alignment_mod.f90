module phi_crossing_alignment_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use plagrange_mod, only: plagrange_coeff

    implicit none

    private

    integer, parameter, public :: crossing_ok = 0
    integer, parameter, public :: crossing_no_bracket = 1
    integer, parameter, public :: crossing_non_simple = 2
    integer, parameter, public :: crossing_invalid_stencil = 3
    integer, parameter :: crossing_stencil_size = 6

    public :: bracketed_phi_crossing
    public :: claim_crossing_target
    public :: nearest_phi_crossing_offset

contains

    pure integer function nearest_phi_crossing_offset(bhat, eta_cross) &
            result(offset)
        real(dp), intent(in) :: bhat(-1:1)
        real(dp), intent(in) :: eta_cross

        real(dp) :: residual(-1:1)

        residual = abs(bhat*eta_cross - 1.0_dp)
        offset = 0
        if (residual(-1) < residual(offset)) offset = -1
        if (residual(1) < residual(offset)) offset = 1
    end function nearest_phi_crossing_offset

    subroutine bracketed_phi_crossing(ibeg, iend, phi, bhat, bracket_left, &
            eta_cross, phi_cross, stencil_left, &
            weights, derivative_weights, target, status)
        integer, intent(in) :: ibeg, iend, bracket_left
        real(dp), intent(in) :: phi(ibeg:iend), bhat(ibeg:iend)
        real(dp), intent(in) :: eta_cross
        real(dp), intent(out) :: phi_cross
        integer, intent(out) :: stencil_left, target, status
        real(dp), intent(out) :: weights(crossing_stencil_size)
        real(dp), intent(out) :: derivative_weights(crossing_stencil_size)

        integer, parameter :: maximum_iterations = 100
        integer, parameter :: monotonic_samples = 9
        real(dp) :: bhat_cross, bhat_derivative, bhat_scale, derivative_scale
        real(dp) :: derivative_reference, f_cross, f_left, f_right
        real(dp) :: left, phi_scale, phi_tolerance, right, sample_phi
        real(dp) :: value_tolerance
        integer :: iteration, sample

        phi_cross = 0.0_dp
        stencil_left = ibeg
        target = ibeg
        status = crossing_invalid_stencil
        weights = 0.0_dp
        derivative_weights = 0.0_dp

        if (iend - ibeg + 1 < crossing_stencil_size) return
        if (bracket_left < ibeg .or. bracket_left >= iend) return
        if (eta_cross <= 0.0_dp) return
        if (any(phi(ibeg + 1:iend) <= phi(ibeg:iend - 1))) return

        stencil_left = max(ibeg, min(iend - crossing_stencil_size + 1, &
            bracket_left - crossing_stencil_size/2 + 1))
        left = phi(bracket_left)
        right = phi(bracket_left + 1)
        phi_scale = max(1.0_dp, abs(left), abs(right))
        phi_tolerance = 128.0_dp*epsilon(1.0_dp)*phi_scale
        bhat_scale = max(1.0_dp, &
            maxval(abs(bhat(stencil_left: &
            stencil_left + crossing_stencil_size - 1))))
        value_tolerance = 128.0_dp*epsilon(1.0_dp)*bhat_scale*eta_cross
        derivative_scale = 128.0_dp*epsilon(1.0_dp)*bhat_scale/(right - left)

        call evaluate_profile(left, &
            phi(stencil_left:stencil_left + crossing_stencil_size - 1), &
            bhat(stencil_left:stencil_left + crossing_stencil_size - 1), &
            bhat_cross, bhat_derivative, weights, &
            derivative_weights)
        f_left = 1.0_dp - eta_cross*bhat_cross
        derivative_reference = bhat_derivative
        if (abs(derivative_reference) <= derivative_scale) then
            status = crossing_non_simple
            return
        end if

        call evaluate_profile(right, &
            phi(stencil_left:stencil_left + crossing_stencil_size - 1), &
            bhat(stencil_left:stencil_left + crossing_stencil_size - 1), &
            bhat_cross, bhat_derivative, weights, &
            derivative_weights)
        f_right = 1.0_dp - eta_cross*bhat_cross
        if (abs(bhat_derivative) <= derivative_scale .or. &
            bhat_derivative*derivative_reference <= 0.0_dp) then
            status = crossing_non_simple
            return
        end if

        do sample = 1, monotonic_samples - 2
            sample_phi = left + (right - left)*real(sample, dp)/ &
                real(monotonic_samples - 1, dp)
            call evaluate_profile(sample_phi, &
                phi(stencil_left:stencil_left + crossing_stencil_size - 1), &
                bhat(stencil_left:stencil_left + crossing_stencil_size - 1), &
                bhat_cross, bhat_derivative, weights, &
                derivative_weights)
            if (abs(bhat_derivative) <= derivative_scale .or. &
                bhat_derivative*derivative_reference <= 0.0_dp) then
                status = crossing_non_simple
                return
            end if
        end do

        if (abs(f_left) <= value_tolerance) then
            phi_cross = left
        else if (abs(f_right) <= value_tolerance) then
            phi_cross = right
        else
            if (f_left*f_right >= 0.0_dp) then
                status = crossing_no_bracket
                return
            end if
            do iteration = 1, maximum_iterations
                phi_cross = 0.5_dp*(left + right)
                call evaluate_profile(phi_cross, &
                    phi(stencil_left:stencil_left + crossing_stencil_size - 1), &
                    bhat(stencil_left:stencil_left + crossing_stencil_size - 1), &
                    bhat_cross, bhat_derivative, weights, &
                    derivative_weights)
                f_cross = 1.0_dp - eta_cross*bhat_cross
                if (abs(f_cross) <= value_tolerance .or. &
                    right - left <= phi_tolerance) exit
                if (f_left*f_cross < 0.0_dp) then
                    right = phi_cross
                    f_right = f_cross
                else
                    left = phi_cross
                    f_left = f_cross
                end if
            end do
        end if

        call evaluate_profile(phi_cross, &
            phi(stencil_left:stencil_left + crossing_stencil_size - 1), &
            bhat(stencil_left:stencil_left + crossing_stencil_size - 1), &
            bhat_cross, bhat_derivative, weights, &
            derivative_weights)
        if (abs(bhat_derivative) <= derivative_scale) then
            status = crossing_non_simple
            return
        end if

        if (abs(phi_cross - phi(bracket_left)) <= phi_tolerance) then
            target = bracket_left
        else if (abs(phi_cross - phi(bracket_left + 1)) <= phi_tolerance) then
            target = bracket_left + 1
        else if (phi_cross - phi(bracket_left) <= &
                phi(bracket_left + 1) - phi_cross) then
            target = bracket_left
        else
            target = bracket_left + 1
        end if
        if (target == ibeg .and. phi_cross > phi(ibeg) + phi_tolerance) &
            target = ibeg + 1
        if (target == iend .and. phi_cross < phi(iend) - phi_tolerance) &
            target = iend - 1
        status = crossing_ok
    end subroutine bracketed_phi_crossing

    logical function claim_crossing_target(claimed, target) result(accepted)
        logical, intent(inout) :: claimed(0:)
        integer, intent(in) :: target

        accepted = .false.
        if (target < lbound(claimed, 1) .or. target > ubound(claimed, 1)) return
        if (claimed(target)) return
        claimed(target) = .true.
        accepted = .true.
    end function claim_crossing_target

    subroutine evaluate_profile(phi_value, phi, bhat, &
            bhat_value, bhat_derivative, weights, &
            derivative_weights)
        real(dp), intent(in) :: phi_value
        real(dp), intent(in) :: phi(crossing_stencil_size)
        real(dp), intent(in) :: bhat(crossing_stencil_size)
        real(dp), intent(out) :: bhat_value, bhat_derivative
        real(dp), intent(out) :: weights(crossing_stencil_size)
        real(dp), intent(out) :: derivative_weights(crossing_stencil_size)

        real(dp) :: coefficients(0:1, crossing_stencil_size)
        call plagrange_coeff(crossing_stencil_size, 1, phi_value, &
            phi, coefficients)
        weights = coefficients(0, :)
        derivative_weights = coefficients(1, :)
        bhat_value = sum(weights*bhat)
        bhat_derivative = sum(derivative_weights*bhat)
    end subroutine evaluate_profile

end module phi_crossing_alignment_mod
