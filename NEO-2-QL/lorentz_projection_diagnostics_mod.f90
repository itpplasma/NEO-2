module lorentz_projection_diagnostics_mod
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
    private

    type, public :: local_projection_residuals
        real(real64) :: source_defect(3)
        real(real64) :: source_scale(3)
        real(real64) :: flux_null(3)
        real(real64) :: flux_scale(3)
        real(real64) :: measure
        real(real64) :: left_transport
        real(real64) :: left_scale
        real(real64) :: right_transport
        real(real64) :: right_scale
        real(real64) :: intertwining
        real(real64) :: intertwining_scale
    end type local_projection_residuals

    character(len=:), allocatable, save :: output_filename
    integer, save :: record_sequence = 0
    logical, save :: output_initialized = .false.

    public :: compute_local_projection_residuals
    public :: record_local_projection_residuals

contains

    subroutine compute_local_projection_residuals(source_p, source_m, flux_p, &
            flux_m, amat_p_p, amat_m_p, amat_p_m, amat_m_m, eta_l, eta_r, &
            eta_boundary_l, eta_boundary_r, residuals, ierr)
        real(real64), intent(in) :: source_p(:, :), source_m(:, :)
        real(real64), intent(in) :: flux_p(:, :), flux_m(:, :)
        real(real64), intent(in) :: amat_p_p(:, :), amat_m_p(:, :)
        real(real64), intent(in) :: amat_p_m(:, :), amat_m_m(:, :)
        real(real64), intent(in) :: eta_l(0:), eta_r(0:)
        real(real64), intent(in) :: eta_boundary_l, eta_boundary_r
        type(local_projection_residuals), intent(out) :: residuals
        integer, intent(out) :: ierr

        real(real64), allocatable :: matrix(:, :), weights_in(:), weights_out(:)
        real(real64), allocatable :: applied_weights(:), column_sum(:)
        real(real64), allocatable :: intertwining_matrix(:, :)
        integer :: nl, nr, ndim

        nl = size(amat_p_p, 2)
        nr = size(amat_p_p, 1)
        call validate_shapes(nl, nr, source_p, source_m, flux_p, flux_m, &
            amat_m_p, amat_p_m, amat_m_m, eta_l, eta_r, ierr)
        if (ierr /= 0) return
        ndim = nl + nr
        allocate(matrix(ndim, ndim), weights_in(ndim), weights_out(ndim))
        allocate(applied_weights(ndim), column_sum(ndim))
        call assemble_map(matrix, amat_p_p, amat_m_p, amat_p_m, amat_m_m, nl, nr)
        call assemble_weights(weights_in, weights_out, eta_l, eta_r, &
            eta_boundary_l, eta_boundary_r, nl, nr, ierr)
        if (ierr /= 0) return
        residuals%measure = sum(weights_out)
        applied_weights = matmul(matrix, weights_in)
        column_sum = sum(matrix, dim=1)
        residuals%source_defect = sum(source_p, dim=1) + sum(source_m, dim=1)
        residuals%source_scale = max(sum(abs(source_p), dim=1) + &
            sum(abs(source_m), dim=1), tiny(1.0_real64))
        residuals%flux_null = matmul(flux_p, weights_in(1:nl)) + &
            matmul(flux_m, weights_in(nl + 1:ndim))
        residuals%flux_scale = max(matmul(abs(flux_p), abs(weights_in(1:nl))) + &
            matmul(abs(flux_m), abs(weights_in(nl + 1:ndim))), tiny(1.0_real64))
        residuals%left_transport = maxval(abs(column_sum - 1.0_real64))
        residuals%left_scale = max(maxval(sum(abs(matrix), dim=1) + 1.0_real64), &
            tiny(1.0_real64))
        residuals%right_transport = maxval(abs(applied_weights - weights_out))
        residuals%right_scale = max(maxval(matmul(abs(matrix), abs(weights_in)) + &
            abs(weights_out)), tiny(1.0_real64))
        allocate(intertwining_matrix(ndim, ndim))
        intertwining_matrix = (spread(applied_weights, 2, ndim) - &
            spread(weights_out, 2, ndim)*spread(column_sum, 1, ndim)) &
            /residuals%measure
        residuals%intertwining = maxval(abs(intertwining_matrix))
        residuals%intertwining_scale = max(maxval((spread(abs(applied_weights), &
            2, ndim) + spread(abs(weights_out), 2, ndim)* &
            spread(abs(column_sum), 1, ndim))/residuals%measure), &
            tiny(1.0_real64))
        ierr = 0
    end subroutine compute_local_projection_residuals

    subroutine validate_shapes(nl, nr, source_p, source_m, flux_p, flux_m, &
            amat_m_p, amat_p_m, amat_m_m, eta_l, eta_r, ierr)
        integer, intent(in) :: nl, nr
        real(real64), intent(in) :: source_p(:, :), source_m(:, :)
        real(real64), intent(in) :: flux_p(:, :), flux_m(:, :)
        real(real64), intent(in) :: amat_m_p(:, :), amat_p_m(:, :)
        real(real64), intent(in) :: amat_m_m(:, :), eta_l(0:), eta_r(0:)
        integer, intent(out) :: ierr

        ierr = 1
        if (nl < 1 .or. nr < 1) return
        if (size(source_p, 1) /= nr .or. size(source_m, 1) /= nl) return
        if (size(source_p, 2) /= 3 .or. size(source_m, 2) /= 3) return
        if (size(flux_p, 1) /= 3 .or. size(flux_p, 2) /= nl) return
        if (size(flux_m, 1) /= 3 .or. size(flux_m, 2) /= nr) return
        if (any(shape(amat_m_p) /= [nr, nr])) return
        if (any(shape(amat_p_m) /= [nl, nl])) return
        if (any(shape(amat_m_m) /= [nl, nr])) return
        if (size(eta_l) < nl .or. size(eta_r) < nr) return
        ierr = 0
    end subroutine validate_shapes

    subroutine assemble_map(matrix, amat_p_p, amat_m_p, amat_p_m, amat_m_m, &
            nl, nr)
        real(real64), intent(out) :: matrix(:, :)
        real(real64), intent(in) :: amat_p_p(:, :), amat_m_p(:, :)
        real(real64), intent(in) :: amat_p_m(:, :), amat_m_m(:, :)
        integer, intent(in) :: nl, nr

        matrix(1:nr, 1:nl) = amat_p_p
        matrix(1:nr, nl + 1:nl + nr) = amat_m_p
        matrix(nr + 1:nr + nl, 1:nl) = amat_p_m
        matrix(nr + 1:nr + nl, nl + 1:nl + nr) = amat_m_m
    end subroutine assemble_map

    subroutine assemble_weights(weights_in, weights_out, eta_l, eta_r, &
            eta_boundary_l, eta_boundary_r, nl, nr, ierr)
        real(real64), intent(out) :: weights_in(:), weights_out(:)
        real(real64), intent(in) :: eta_l(0:), eta_r(0:)
        real(real64), intent(in) :: eta_boundary_l, eta_boundary_r
        integer, intent(in) :: nl, nr
        integer, intent(out) :: ierr
        real(real64), allocatable :: weights_l(:), weights_r(:)

        allocate(weights_l(nl), weights_r(nr))
        call boundary_weights(eta_l, eta_boundary_l, weights_l, ierr)
        if (ierr /= 0) return
        call boundary_weights(eta_r, eta_boundary_r, weights_r, ierr)
        if (ierr /= 0) return
        weights_in = [weights_l, weights_r]
        weights_out = [weights_r, weights_l]
    end subroutine assemble_weights

    subroutine boundary_weights(eta, eta_boundary, weights, ierr)
        real(real64), intent(in) :: eta(0:), eta_boundary
        real(real64), intent(out) :: weights(:)
        integer, intent(out) :: ierr
        integer :: n

        n = size(weights)
        if (n > 1) weights(1:n - 1) = eta(1:n - 1) - eta(0:n - 2)
        weights(n) = eta_boundary - eta(n - 1)
        ierr = 0
        if (.not. all(ieee_is_finite(weights)) .or. &
            .not. ieee_is_finite(eta_boundary)) ierr = 2
        if (ierr == 0 .and. any(weights <= 0.0_real64)) ierr = 2
    end subroutine boundary_weights

    subroutine record_local_projection_residuals(tag, source_p, source_m, &
            flux_p, flux_m, amat_p_p, amat_m_p, amat_p_m, amat_m_m, eta_l, &
            eta_r, eta_boundary_l, eta_boundary_r, ierr)
        integer, intent(in) :: tag
        real(real64), intent(in) :: source_p(:, :), source_m(:, :)
        real(real64), intent(in) :: flux_p(:, :), flux_m(:, :)
        real(real64), intent(in) :: amat_p_p(:, :), amat_m_p(:, :)
        real(real64), intent(in) :: amat_p_m(:, :), amat_m_m(:, :)
        real(real64), intent(in) :: eta_l(0:), eta_r(0:)
        real(real64), intent(in) :: eta_boundary_l, eta_boundary_r
        integer, intent(out) :: ierr
        type(local_projection_residuals) :: residuals
        integer :: force, iunit, status

        call initialize_output(ierr)
        if (ierr /= 0 .or. .not. allocated(output_filename)) return
        call compute_local_projection_residuals(source_p, source_m, flux_p, &
            flux_m, amat_p_p, amat_m_p, amat_p_m, amat_m_m, eta_l, eta_r, &
            eta_boundary_l, eta_boundary_r, residuals, ierr)
        if (ierr /= 0) return
        open(newunit=iunit, file=output_filename, status='old', position='append', &
            action='write', iostat=status)
        if (status /= 0) then
            ierr = 3
            return
        end if
        do force = 1, 3
            call write_value(iunit, tag, 'source_defect', force, &
                residuals%source_defect(force), residuals%source_scale(force), status)
            call write_value(iunit, tag, 'flux_null', force, &
                residuals%flux_null(force), residuals%flux_scale(force), status)
        end do
        call write_value(iunit, tag, 'left_transport', -1, &
            residuals%left_transport, residuals%left_scale, status)
        call write_value(iunit, tag, 'right_transport', -1, &
            residuals%right_transport, residuals%right_scale, status)
        call write_value(iunit, tag, 'intertwining', -1, residuals%intertwining, &
            residuals%intertwining_scale, status)
        call write_value(iunit, tag, 'measure', -1, residuals%measure, &
            max(abs(residuals%measure), tiny(1.0_real64)), status)
        close(iunit, iostat=ierr)
        if (status /= 0 .or. ierr /= 0) ierr = 3
    end subroutine record_local_projection_residuals

    subroutine write_value(iunit, tag, kind, index, value, scale, status)
        integer, intent(in) :: iunit, tag, index
        character(len=*), intent(in) :: kind
        real(real64), intent(in) :: value, scale
        integer, intent(inout) :: status

        if (status /= 0) return
        if (.not. ieee_is_finite(value) .or. .not. ieee_is_finite(scale)) then
            status = 1
            return
        end if
        if (scale <= 0.0_real64) then
            status = 1
            return
        end if
        record_sequence = record_sequence + 1
        write(iunit, '(i0,",",i0,",",a,",",i0,2(",",es25.16e3))', &
            iostat=status) record_sequence, tag, trim(kind), index, value, scale
    end subroutine write_value

    subroutine initialize_output(ierr)
        integer, intent(out) :: ierr
        integer :: iunit, length, status

        ierr = 0
        if (output_initialized) return
        output_initialized = .true.
        call get_environment_variable('NEO2_LOCAL_PROJECTION_TRACE_FILE', &
            length=length, status=status)
        if (status /= 0 .or. length == 0) return
        allocate(character(len=length) :: output_filename)
        call get_environment_variable('NEO2_LOCAL_PROJECTION_TRACE_FILE', &
            value=output_filename, status=status)
        if (status == 0) open(newunit=iunit, file=output_filename, &
            status='replace', action='write', iostat=status)
        if (status == 0) then
            write(iunit, '(a)', iostat=status) &
                'sequence,propagator,kind,index,value,scale'
            close(iunit, iostat=ierr)
        end if
        if (status /= 0 .or. ierr /= 0) ierr = 3
        if (ierr /= 0 .and. allocated(output_filename)) deallocate(output_filename)
    end subroutine initialize_output
end module lorentz_projection_diagnostics_mod
