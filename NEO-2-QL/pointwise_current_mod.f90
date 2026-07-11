module pointwise_current_mod
    use nrtype, only: dp
    implicit none
    private

    public :: current_profile_from_flux_vector

contains

    !> Pointwise parallel-current profile along the field line.
    !>
    !> The phi-integrated current channel is
    !>   qflux(2, k) = sum_col flux_vector(2, col) * source_vector(col, k),
    !> and flux_vector's columns partition into disjoint per-field-line-point
    !> blocks: point i owns columns block_base(i)+1 .. block_base(i)+nblk with
    !> nblk = 2*(lag+1)*(npassing(i)+1). The block dot product is that point's
    !> contribution to the integral. Because flux_vector already carries the
    !> current-row sign (the counter-passing block enters with an explicit
    !> minus, unlike density and energy) and the step_factor phi-quadrature
    !> weight, the contribution is the current integrand times its quadrature
    !> weight; dividing by phi_weight(i) recovers the pointwise current density
    !> j_par at that point. Summing all contributions reproduces qflux(2, k)
    !> exactly; that identity is the validation contract for the profile.
    !>
    !> This routine only decomposes the already-assembled flux_vector /
    !> source_vector contraction; it introduces no new sign, weight, or
    !> normalization, so the pointwise current inherits the verified
    !> convention of the harmonic response.
    subroutine current_profile_from_flux_vector(flux_row, source_col, &
            block_base, block_npassing, lag, phi_weight, contribution, &
            current_density, channel)
        real(dp), intent(in) :: flux_row(:)         !< flux_vector(2, :)
        real(dp), intent(in) :: source_col(:)       !< source_vector(:, k)
        integer, intent(in) :: block_base(:)        !< 0-based column offset per point
        integer, intent(in) :: block_npassing(:)    !< npassing per point
        integer, intent(in) :: lag                  !< Laguerre degree
        real(dp), intent(in) :: phi_weight(:)       !< step_factor phi measure per point
        real(dp), intent(out) :: contribution(:)    !< per-point integral contribution
        real(dp), intent(out) :: current_density(:) !< contribution / phi_weight (pointwise j_par)
        real(dp), intent(out) :: channel            !< sum(contribution) = qflux(2, k)
        integer :: i, npoint, base, nblk, col
        real(dp) :: acc

        npoint = size(block_base)
        if (size(block_npassing) /= npoint .or. size(phi_weight) /= npoint &
                .or. size(contribution) /= npoint &
                .or. size(current_density) /= npoint) &
            error stop 'current_profile_from_flux_vector: point arrays disagree in size'

        channel = 0.0_dp
        do i = 1, npoint
            base = block_base(i)
            nblk = 2*(lag + 1)*(block_npassing(i) + 1)
            if (base < 0 .or. base + nblk > size(flux_row) &
                    .or. base + nblk > size(source_col)) &
                error stop 'current_profile_from_flux_vector: block exceeds vector bounds'
            acc = 0.0_dp
            do col = base + 1, base + nblk
                acc = acc + flux_row(col)*source_col(col)
            end do
            contribution(i) = acc
            if (phi_weight(i) /= 0.0_dp) then
                current_density(i) = acc/phi_weight(i)
            else
                current_density(i) = 0.0_dp
            end if
            channel = channel + acc
        end do
    end subroutine current_profile_from_flux_vector
end module pointwise_current_mod
