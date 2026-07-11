program test_helical_source
    use helical_source_mod, only: add_helical_source, &
        apply_reconstructed_incoming_rows
    implicit none

    integer, parameter :: dp = kind(1.0d0), ibeg = 0, iend = 4
    integer, parameter :: replay_ibeg = 0, replay_iend = 2, replay_lag = 1
    complex(dp) :: even_source(24, 3), odd_source(24, 3), profile(3, ibeg:iend)
    complex(dp) :: forward_even(24, 3), forward_odd(24, 3)
    complex(dp) :: backward_even(24, 3), backward_odd(24, 3)
    complex(dp) :: weighted_sum
    real(dp) :: factors(ibeg:iend), moments(0:0), zero_factors(ibeg:iend)
    integer :: i, ind_start(ibeg:iend), npl(ibeg:iend)
    real(dp) :: replay_source(24, 3), original_source(24, 3)
    real(dp) :: flux_left(4, 3), flux_right(4, 3)
    integer :: replay_ind_start(replay_ibeg:replay_iend)
    integer :: replay_npl(replay_ibeg:replay_iend)
    logical :: incoming_row(24)

    even_source = (0.0d0, 0.0d0)
    odd_source = (0.0d0, 0.0d0)
    forward_even = (0.0d0, 0.0d0)
    forward_odd = (0.0d0, 0.0d0)
    backward_even = (0.0d0, 0.0d0)
    backward_odd = (0.0d0, 0.0d0)
    zero_factors = 0.0d0
    factors = [1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0]
    moments = 2.0d0
    npl = [1, 2, 1, 2, 1]
    ind_start = [0, 4, 10, 14, 20]
    do i = ibeg, iend
        profile(:, i) = cmplx([real(i + 1, dp), real(2*i + 1, dp), real(3*i + 2, dp)], &
            [real(3*i + 1, dp), real(i + 2, dp), real(2*i + 4, dp)], kind=dp)
    end do

    call add_helical_source(even_source, profile, moments, 1, 1, ibeg, iend, &
        0, npl, ind_start, factors, factors, factors, factors)
    call add_helical_source(odd_source, profile, moments, 1, -1, ibeg, iend, &
        0, npl, ind_start, factors, factors, factors, factors)
    call add_helical_source(forward_even, profile, moments, 1, 1, ibeg, iend, &
        0, npl, ind_start, factors, factors, zero_factors, zero_factors)
    call add_helical_source(forward_odd, profile, moments, 1, -1, ibeg, iend, &
        0, npl, ind_start, factors, factors, zero_factors, zero_factors)
    call add_helical_source(backward_even, profile, moments, 1, 1, ibeg, iend, &
        0, npl, ind_start, zero_factors, zero_factors, factors, factors)
    call add_helical_source(backward_odd, profile, moments, 1, -1, ibeg, iend, &
        0, npl, ind_start, zero_factors, zero_factors, factors, factors)

    if (any(even_source(:, 2:3) /= (0.0d0, 0.0d0))) &
        error stop 'FAIL: helical source changed an unselected force column'
    if (maxval(abs(forward_even - forward_odd)) > 1.0d-14) &
        error stop 'FAIL: sigma sign changed the forward source'
    if (maxval(abs(backward_even + backward_odd)) > 1.0d-14) &
        error stop 'FAIL: sigma-odd backward source has the wrong sign'
    if (maxval(abs(aimag(even_source))) == 0.0d0) &
        error stop 'FAIL: complex source quadrature was discarded'
    weighted_sum = sum([(even_source(i, 1)*real(i, dp), i = 1, 24)])
    if (abs(weighted_sum - cmplx(11874.666666666666d0, 14272.0d0, dp)) > 1.0d-10) &
        error stop 'FAIL: shared source stencil changed'
    if (abs(sum(abs(even_source(:, 1))**2) - 111813.33333333333d0) > 1.0d-9) &
        error stop 'FAIL: shared source stencil norm changed'

    replay_source = reshape([(real(i, dp), i = 1, size(replay_source))], &
        shape(replay_source))
    original_source = replay_source
    flux_left = reshape([(100.0_dp + real(i, dp), i = 1, size(flux_left))], &
        shape(flux_left))
    flux_right = reshape([(200.0_dp + real(i, dp), i = 1, size(flux_right))], &
        shape(flux_right))
    replay_npl = 1
    replay_ind_start = [0, 8, 16]
    call apply_reconstructed_incoming_rows(replay_source, flux_left, flux_right, &
        replay_ibeg, replay_iend, replay_lag, replay_npl, replay_ind_start)

    if (any(replay_source(1:2, :) /= flux_left(1:2, :))) &
        error stop 'FAIL: first left incoming Laguerre block was not restored'
    if (any(replay_source(5:6, :) /= flux_left(3:4, :))) &
        error stop 'FAIL: second left incoming Laguerre block was not restored'
    if (any(replay_source(19:20, :) /= flux_right(2:1:-1, :))) &
        error stop 'FAIL: first right incoming Laguerre block was not restored'
    if (any(replay_source(23:24, :) /= flux_right(4:3:-1, :))) &
        error stop 'FAIL: second right incoming Laguerre block was not restored'
    incoming_row = .false.
    incoming_row([1, 2, 5, 6, 19, 20, 23, 24]) = .true.
    do i = 1, size(replay_source, 1)
        if (incoming_row(i)) cycle
        if (any(replay_source(i, :) /= original_source(i, :))) &
            error stop 'FAIL: reconstruction changed an interior or outgoing row'
    end do

    print *, 'All tests passed!'
end program test_helical_source
