module qflux_profile_mod
    use nrtype, only: dp
    implicit none
    private

    public :: qflux_contributions_from_flux_vector

contains

    !> Per-point partition of the raw qflux current channel.
    !>
    !> The phi-integrated current channel is
    !>   qflux(2, k) = sum_col flux_vector(2, col) * source_vector(col, k),
    !> and flux_vector's columns partition into disjoint per-field-line-point
    !> blocks: point i owns columns block_base(i)+1 .. block_base(i)+nblk with
    !> nblk = 2*(lag+1)*(npassing(i)+1). The block dot product is that point's
    !> contribution to the integral. Because flux_vector already carries the
    !> current-row sign (the counter-passing block enters with an explicit
    !> minus, unlike density and energy) and the step_factor phi-quadrature
    !> weight, each block dot product is one contribution to the normalized
    !> qflux current channel. Summing all contributions reproduces qflux(2, k)
    !> exactly; that partition identity is the validation contract.
    !>
    !> This routine only decomposes the already-assembled flux_vector /
    !> source_vector contraction; it introduces no new sign, weight, or
    !> normalization. It does not produce a pointwise physical current density.
    subroutine qflux_contributions_from_flux_vector(flux_row, source_col, &
            block_base, block_npassing, lag, contribution, channel)
        real(dp), intent(in) :: flux_row(:) !< flux_vector(2, :)
        real(dp), intent(in) :: source_col(:) !< source_vector(:, k)
        integer, intent(in) :: block_base(:) !< 0-based column offset per point
        integer, intent(in) :: block_npassing(:) !< npassing per point
        integer, intent(in) :: lag !< Laguerre degree
        real(dp), intent(out) :: contribution(:) !< per-point integral contribution
        real(dp), intent(out) :: channel !< sum(contribution) = qflux(2, k)
        integer :: i, npoint, base, nblk, col
        real(dp) :: acc

        npoint = size(block_base)
        if (size(block_npassing) /= npoint .or. size(contribution) /= npoint) &
            error stop 'qflux_contributions_from_flux_vector: point arrays disagree in size'

        channel = 0.0_dp
        do i = 1, npoint
            base = block_base(i)
            nblk = 2*(lag + 1)*(block_npassing(i) + 1)
            if (base < 0 .or. base + nblk > size(flux_row) &
                .or. base + nblk > size(source_col)) &
                error stop 'qflux_contributions_from_flux_vector: block exceeds vector bounds'
            acc = 0.0_dp
            do col = base + 1, base + nblk
                acc = acc + flux_row(col)*source_col(col)
            end do
            contribution(i) = acc
            channel = channel + acc
        end do
    end subroutine qflux_contributions_from_flux_vector
end module qflux_profile_mod
