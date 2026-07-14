program test_phi_crossing_alignment
    use, intrinsic :: iso_fortran_env, only: real64
    use phi_crossing_alignment_mod, only: bracketed_phi_crossing, &
        claim_crossing_target, crossing_no_bracket, crossing_non_simple, &
        crossing_ok, nearest_phi_crossing_offset

    implicit none

    real(real64), parameter :: eta_cross = 1.0076495897856284_real64
    real(real64), parameter :: tolerance = 1.0e-12_real64
    real(real64) :: bhat(-1:1)
    real(real64) :: bhat_profile(0:7), derivative_weights(6)
    real(real64) :: phi_cross, phi_profile(0:7), weights(6)
    real(real64) :: eta_root
    logical :: claimed(0:7)
    integer :: offset, status, stencil_left, target

    bhat = [0.99203602638228305_real64, 0.99274003698367863_real64, &
        0.99340_real64]
    offset = nearest_phi_crossing_offset(bhat, eta_cross)
    if (offset /= 0) &
        error stop 'FAIL: closest current point was not selected'

    bhat = [1.0_real64/1.00766_real64, 1.0_real64/1.00790_real64, &
        1.0_real64/1.00820_real64]
    offset = nearest_phi_crossing_offset(bhat, eta_cross)
    if (offset /= -1) &
        error stop 'FAIL: closest preceding point was not selected'

    bhat = [1.0_real64/1.00820_real64, 1.0_real64/1.00790_real64, &
        1.0_real64/1.00766_real64]
    offset = nearest_phi_crossing_offset(bhat, eta_cross)
    if (offset /= 1) &
        error stop 'FAIL: closest following point was not selected'

    bhat = [1.0_real64/1.00770_real64, 1.0_real64/1.00820_real64, &
        1.0_real64/1.00766_real64]
    offset = nearest_phi_crossing_offset(bhat, eta_cross)
    if (offset /= 1) &
        error stop 'FAIL: nearest neighbor was not selected'

    phi_profile = [(real(offset, real64), offset = 0, 7)]
    bhat_profile = 1.0_real64 + 0.2_real64*phi_profile
    eta_root = 1.0_real64/(1.0_real64 + 0.2_real64*2.3_real64)
    call bracketed_phi_crossing(0, 7, phi_profile, bhat_profile, 2, &
        eta_root, phi_cross, stencil_left, weights, derivative_weights, &
        target, status)
    if (status /= crossing_ok) error stop 'FAIL: linear root was rejected'
    if (abs(phi_cross - 2.3_real64) > tolerance) &
        error stop 'FAIL: linear root was not recovered'
    if (target /= 2) error stop 'FAIL: nearest bracket node was not selected'
    if (abs(sum(weights) - 1.0_real64) > tolerance) &
        error stop 'FAIL: interpolation does not preserve constants'
    if (abs(sum(weights*bhat_profile(stencil_left:stencil_left + 5)) &
        *eta_root - 1.0_real64) > tolerance) &
        error stop 'FAIL: interpolated field misses the pitch boundary'
    if (abs(sum(derivative_weights* &
        bhat_profile(stencil_left:stencil_left + 5)) - 0.2_real64) &
        > tolerance) &
        error stop 'FAIL: field derivative is inconsistent'

    eta_root = 1.0_real64/(1.0_real64 + 0.2_real64*0.1_real64)
    call bracketed_phi_crossing(0, 7, phi_profile, bhat_profile, 0, &
        eta_root, phi_cross, stencil_left, weights, derivative_weights, &
        target, status)
    if (status /= crossing_ok .or. target /= 1) &
        error stop 'FAIL: periodic endpoint was selected for movement'

    eta_root = 1.0_real64/4.0_real64
    call bracketed_phi_crossing(0, 7, phi_profile, bhat_profile, 2, &
        eta_root, phi_cross, stencil_left, weights, derivative_weights, &
        target, status)
    if (status /= crossing_no_bracket) &
        error stop 'FAIL: missing bracket was not rejected'

    bhat_profile = 1.0_real64 + (phi_profile - 2.5_real64)**3
    call bracketed_phi_crossing(0, 7, phi_profile, bhat_profile, 2, &
        1.0_real64, phi_cross, stencil_left, weights, derivative_weights, &
        target, status)
    if (status /= crossing_non_simple) &
        error stop 'FAIL: non-simple root was not rejected'

    claimed = .false.
    if (.not. claim_crossing_target(claimed, 3)) &
        error stop 'FAIL: first crossing claim was rejected'
    if (claim_crossing_target(claimed, 3)) &
        error stop 'FAIL: competing crossing claim was accepted'

    print '(a)', 'All tests passed!'
end program test_phi_crossing_alignment
