program test_qflux_profile
    use qflux_profile_mod, only: qflux_contributions_from_flux_vector, &
        qflux_point_components_from_flux_vector
    implicit none

    integer, parameter :: dp = kind(1.0d0)
    integer, parameter :: lag = 1
    integer, parameter :: npoint = 3
    integer :: block_npassing(npoint), block_base(npoint)
    integer :: i, m, nblk, base, col, first, ntot
    real(dp), allocatable :: flux_row(:), source_col(:)
    real(dp) :: contribution(npoint)
    real(dp) :: raw_plus(npoint), raw_minus(npoint), split_contribution(npoint)
    real(dp) :: step_plus(npoint), step_minus(npoint)
    real(dp) :: expected_plus, expected_minus, split_channel
    real(dp) :: channel, reference, sum_contrib

    block_npassing = [2, 1, 3]
    base = 0
    do i = 1, npoint
        block_base(i) = base
        nblk = 2*(lag + 1)*(block_npassing(i) + 1)
        base = base + nblk
    end do
    ntot = base ! 12 + 8 + 16 = 36
    allocate (flux_row(ntot), source_col(ntot))
    do col = 1, ntot
        flux_row(col) = real(col, dp)*0.5_dp - 3.0_dp ! spans negative and positive
        source_col(col) = 1.0_dp + 0.25_dp*real(mod(col, 5), dp)
    end do
    ! 1. Partition identity: the block decomposition reproduces the full
    !    flux_vector . source_vector contraction (an off-by-one in block base
    !    or size breaks this).
    reference = dot_product(flux_row, source_col)
    call qflux_contributions_from_flux_vector(flux_row, source_col, block_base, &
        block_npassing, lag, contribution, channel)
    if (abs(channel - reference) > 1.0e-11_dp*max(1.0_dp, abs(reference))) &
        error stop 'FAIL: per-point contributions do not sum to the qflux current channel'

    ! 2. Channel equals the sum of per-point contributions.
    sum_contrib = sum(contribution)
    if (abs(channel - sum_contrib) > 1.0e-12_dp*max(1.0_dp, abs(channel))) &
        error stop 'FAIL: channel is not the sum of per-point contributions'

    ! 3. Known closed form: unit flux and source give the block size per point.
    flux_row = 1.0_dp
    source_col = 1.0_dp
    call qflux_contributions_from_flux_vector(flux_row, source_col, block_base, &
        block_npassing, lag, contribution, channel)
    do i = 1, npoint
        nblk = 2*(lag + 1)*(block_npassing(i) + 1)
        if (abs(contribution(i) - real(nblk, dp)) > 1.0e-12_dp) &
            error stop 'FAIL: unit contraction does not equal the block size'
    end do
    if (abs(channel - real(ntot, dp)) > 1.0e-12_dp) &
        error stop 'FAIL: unit contraction total does not equal the vector length'

    ! 4. Current-row sign propagates: a negative flux_vector entry (the
    !    counter-passing current weight) yields a negative contribution.
    flux_row = 1.0_dp
    flux_row(block_base(2) + 1) = -50.0_dp
    source_col = 1.0_dp
    call qflux_contributions_from_flux_vector(flux_row, source_col, block_base, &
        block_npassing, lag, contribution, channel)
    if (contribution(2) >= 0.0_dp) &
        error stop 'FAIL: negative current-row weight does not produce a negative contribution'

    step_plus = [0.25_dp, 0.5_dp, 0.75_dp]
    step_minus = [0.4_dp, 0.6_dp, 0.8_dp]
    source_col = 1.0_dp
    do i = 1, npoint
        nblk = block_npassing(i) + 1
        do m = 0, lag
            first = block_base(i) + 2*m*nblk + 1
            flux_row(first:first + nblk - 1) = &
                step_plus(i)*real(i + m + 1, dp)
            first = first + nblk
            flux_row(first:first + nblk - 1) = &
                -step_minus(i)*real(2*i + m + 1, dp)
        end do
    end do
    call qflux_point_components_from_flux_vector(flux_row, source_col, block_base, &
        block_npassing, lag, step_plus, step_minus, raw_plus, raw_minus, &
        split_contribution, split_channel)
    do i = 1, npoint
        nblk = block_npassing(i) + 1
        expected_plus = real(nblk, dp)*sum([(real(i + m + 1, dp), m = 0, lag)])
        expected_minus = -real(nblk, dp)*sum([(real(2*i + m + 1, dp), m = 0, lag)])
        if (abs(raw_plus(i) - expected_plus) > 1.0e-12_dp) &
            error stop 'FAIL: co-passing component retains the spatial weight'
        if (abs(raw_minus(i) - expected_minus) > 1.0e-12_dp) &
            error stop 'FAIL: counter-passing component loses its current sign'
        if (abs(split_contribution(i) - step_plus(i)*expected_plus &
            - step_minus(i)*expected_minus) > 1.0e-12_dp) &
            error stop 'FAIL: directional components do not reconstruct the contribution'
    end do
    if (abs(split_channel - sum(split_contribution)) > 1.0e-12_dp) &
        error stop 'FAIL: split channel does not equal the contribution sum'

    print *, 'All tests passed!'
end program test_qflux_profile
