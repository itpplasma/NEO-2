module qflux_profile_mod
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use nrtype, only: dp
    implicit none
    private

    character(len=:), allocatable, save :: interface_trace_filename
    integer, save :: interface_trace_sequence = 0
    logical, save :: interface_trace_initialized = .false.

    public :: qflux_contributions_from_flux_vector
    public :: qflux_point_components_from_flux_vector
    public :: record_qflux_interface_traces

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

    subroutine qflux_point_components_from_flux_vector(flux_row, source_col, &
            block_base, block_npassing, lag, step_plus, step_minus, raw_plus, &
            raw_minus, contribution, channel)
        real(dp), intent(in) :: flux_row(:), source_col(:)
        integer, intent(in) :: block_base(:), block_npassing(:), lag
        real(dp), intent(in) :: step_plus(:), step_minus(:)
        real(dp), intent(out) :: raw_plus(:), raw_minus(:), contribution(:)
        real(dp), intent(out) :: channel
        integer :: i, m, npoint, nband, base, first
        real(dp) :: weighted_plus, weighted_minus

        npoint = size(block_base)
        if (size(block_npassing) /= npoint .or. size(step_plus) /= npoint &
            .or. size(step_minus) /= npoint .or. size(raw_plus) /= npoint &
            .or. size(raw_minus) /= npoint .or. size(contribution) /= npoint) &
            error stop 'qflux_point_components_from_flux_vector: point arrays disagree in size'

        channel = 0.0_dp
        do i = 1, npoint
            if (step_plus(i) == 0.0_dp .or. step_minus(i) == 0.0_dp) &
                error stop 'qflux_point_components_from_flux_vector: zero spatial weight'
            nband = block_npassing(i) + 1
            base = block_base(i)
            if (base < 0 .or. base + 2*(lag + 1)*nband > size(flux_row) &
                .or. base + 2*(lag + 1)*nband > size(source_col)) &
                error stop 'qflux_point_components_from_flux_vector: block exceeds vector bounds'

            weighted_plus = 0.0_dp
            weighted_minus = 0.0_dp
            do m = 0, lag
                first = base + 2*m*nband + 1
                weighted_plus = weighted_plus + dot_product( &
                    flux_row(first:first + nband - 1), &
                    source_col(first:first + nband - 1))
                first = first + nband
                weighted_minus = weighted_minus + dot_product( &
                    flux_row(first:first + nband - 1), &
                    source_col(first:first + nband - 1))
            end do
            raw_plus(i) = weighted_plus/step_plus(i)
            raw_minus(i) = weighted_minus/step_minus(i)
            contribution(i) = weighted_plus + weighted_minus
            channel = channel + contribution(i)
        end do
    end subroutine qflux_point_components_from_flux_vector

    subroutine record_qflux_interface_traces(tag, phi, eta_grid, block_base, &
            block_npassing, lag, step_plus, step_minus, flux_row, source_rhs, &
            source_solution, ierr)
        integer, intent(in) :: tag, block_base(:), block_npassing(:), lag
        real(dp), intent(in) :: phi(:), eta_grid(0:), step_plus(:), step_minus(:)
        real(dp), intent(in) :: flux_row(:), source_rhs(:, :)
        real(dp), intent(in) :: source_solution(:, :)
        integer, intent(out) :: ierr
        character(len=5) :: side
        integer :: band, base, column, endpoint, force, iunit, m, nband, point
        integer :: status
        real(dp) :: constant_coordinate(3), source_solution_orthogonal

        ierr = 0
        if (.not. interface_trace_initialized) call initialize_interface_trace(ierr)
        if (ierr /= 0 .or. .not. allocated(interface_trace_filename)) return
        if (size(phi) /= size(block_base) &
            .or. size(block_npassing) /= size(block_base) &
            .or. size(step_plus) /= size(block_base) &
            .or. size(step_minus) /= size(block_base) &
            .or. size(source_rhs, 1) /= size(flux_row) &
            .or. any(shape(source_solution) /= shape(source_rhs)) &
            .or. size(source_rhs, 2) /= 3 &
            .or. lag < 0) then
            ierr = 1
            return
        end if
        if (size(block_base) < 2 &
            .or. maxval(block_npassing) >= size(eta_grid) &
            .or. any(block_npassing < 0) &
            .or. any(step_plus == 0.0_dp) &
            .or. any(step_minus == 0.0_dp)) then
            ierr = 1
            return
        end if
        if (.not. all(ieee_is_finite(phi)) &
            .or. .not. all(ieee_is_finite(eta_grid)) &
            .or. .not. all(ieee_is_finite(step_plus)) &
            .or. .not. all(ieee_is_finite(step_minus)) &
            .or. .not. all(ieee_is_finite(flux_row)) &
            .or. .not. all(ieee_is_finite(source_rhs)) &
            .or. .not. all(ieee_is_finite(source_solution))) then
            ierr = 1
            return
        end if

        open(newunit=iunit, file=interface_trace_filename, status='old', &
            position='append', action='write', iostat=status)
        if (status /= 0) then
            ierr = status
            return
        end if
        do endpoint = 1, 2
            if (endpoint == 1) then
                point = 1
                side = 'left'
            else
                point = size(block_base)
                side = 'right'
            end if
            nband = block_npassing(point) + 1
            base = block_base(point)
            ! The Lorentz collision operator's discrete right-null state is
            ! one constant shared by the co- and counter-passing lag=0 rows.
            ! Export the Euclidean-orthogonal component as a diagnostic; the
            ! raw solution and all solver inputs remain unchanged.
            do force = 1, 3
                constant_coordinate(force) = sum( &
                    source_solution(base + 1:base + nband, force))
                constant_coordinate(force) = constant_coordinate(force) + sum( &
                    source_solution(base + nband + 1:base + 2*nband, force))
            end do
            constant_coordinate = constant_coordinate/real(2*nband, dp)
            do m = 0, lag
                base = block_base(point) + 2*m*nband
                do band = 1, nband
                    column = base + band
                    do force = 1, 3
                        source_solution_orthogonal = &
                            source_solution(column, force)
                        if (m == 0) source_solution_orthogonal = &
                            source_solution_orthogonal - &
                            constant_coordinate(force)
                        call write_interface_trace(iunit, tag, side, 'p', m, &
                            force, band - 1, eta_grid(band - 1), phi(point), &
                            flux_row(column)/step_plus(point), &
                            source_rhs(column, force), &
                            source_solution(column, force), &
                            source_solution_orthogonal, status)
                        if (status /= 0) exit
                    end do
                    if (status /= 0) exit
                    column = base + 2*nband - band + 1
                    do force = 1, 3
                        source_solution_orthogonal = &
                            source_solution(column, force)
                        if (m == 0) source_solution_orthogonal = &
                            source_solution_orthogonal - &
                            constant_coordinate(force)
                        call write_interface_trace(iunit, tag, side, 'm', m, &
                            force, band - 1, eta_grid(band - 1), phi(point), &
                            flux_row(column)/step_minus(point), &
                            source_rhs(column, force), &
                            source_solution(column, force), &
                            source_solution_orthogonal, status)
                        if (status /= 0) exit
                    end do
                    if (status /= 0) exit
                end do
                if (status /= 0) exit
            end do
            if (status /= 0) exit
        end do
        close(iunit, iostat=ierr)
        if (status /= 0) ierr = status
    end subroutine record_qflux_interface_traces

    subroutine write_interface_trace(iunit, tag, side, direction, laguerre, &
            force, eta_index, eta_value, phi, flux_kernel, source_rhs, &
            source_solution, source_solution_orthogonal, status)
        integer, intent(in) :: iunit, tag, laguerre, force, eta_index
        character(len=*), intent(in) :: side, direction
        real(dp), intent(in) :: eta_value, phi, flux_kernel, source_rhs
        real(dp), intent(in) :: source_solution, source_solution_orthogonal
        integer, intent(out) :: status

        interface_trace_sequence = interface_trace_sequence + 1
        write(iunit, '(2(i0,","),a,",",a,",",3(i0,","),5(es25.16e3,","),' // &
            'es25.16e3)', iostat=status) interface_trace_sequence, tag, &
            trim(side), direction, laguerre, force, eta_index, eta_value, phi, &
            flux_kernel, source_rhs, source_solution, &
            source_solution_orthogonal
    end subroutine write_interface_trace

    subroutine initialize_interface_trace(ierr)
        integer, intent(out) :: ierr
        integer :: iunit, length, status

        ierr = 0
        interface_trace_initialized = .true.
        call get_environment_variable('NEO2_INTERFACE_TRACE_FILE', &
            length=length, status=status)
        if (status /= 0 .or. length == 0) return
        allocate(character(len=length) :: interface_trace_filename)
        call get_environment_variable('NEO2_INTERFACE_TRACE_FILE', &
            value=interface_trace_filename, status=status)
        if (status == 0) then
            open(newunit=iunit, file=interface_trace_filename, status='replace', &
                action='write', iostat=status)
        end if
        if (status == 0) then
            write(iunit, '(a)', iostat=status) &
                'sequence,tag,side,direction,laguerre,force,eta_index,eta,phi,' // &
                'flux_kernel,source_rhs,source_solution,' // &
                'source_solution_orthogonal'
            close(iunit, iostat=ierr)
        end if
        if (status /= 0) ierr = status
        if (ierr /= 0 .and. allocated(interface_trace_filename)) &
            deallocate(interface_trace_filename)
    end subroutine initialize_interface_trace
end module qflux_profile_mod
