program test_phi_crossing_alignment
    use, intrinsic :: iso_fortran_env, only: real64
    use phi_crossing_alignment_mod, only: nearest_phi_crossing_offset

    implicit none

    real(real64), parameter :: eta_cross = 1.0076495897856284_real64
    real(real64) :: bhat(-1:1)
    integer :: offset

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

    print '(a)', 'All tests passed!'
end program test_phi_crossing_alignment
