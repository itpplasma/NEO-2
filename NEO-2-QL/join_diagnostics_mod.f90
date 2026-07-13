module join_diagnostics_mod
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
    private

    character(len=:), allocatable, save :: normalization_output_filename
    character(len=:), allocatable, save :: join_end_output_filename
    integer, save :: normalization_sequence = 0
    integer, save :: join_end_record = 0
    integer, save :: join_end_sequence = 0
    logical, save :: normalization_output_initialized = .false.
    logical, save :: join_end_output_initialized = .false.

    public :: join_end_trace_enabled, record_join_end_compatibility
    public :: report_join_failure, validate_join_normalization

contains

    subroutine report_join_failure(stage, info, old_tags, new_tags, ndim, ndim1, &
            matrix, rhs, ierr)
        integer, intent(in) :: stage, info, old_tags(2), new_tags(2), ndim, ndim1
        real(real64), intent(in) :: matrix(:, :), rhs(:, :)
        integer, intent(out) :: ierr

        write(*, '(a,i0)') 'join_ripples: DGBSV stage=', stage
        write(*, '(a,i0)') 'join_ripples: DGBSV info=', info
        write(*, '(a,2(1x,i0))') 'join_ripples: old propagator tags=', old_tags
        write(*, '(a,2(1x,i0))') 'join_ripples: new propagator tags=', new_tags
        write(*, '(a,2(1x,i0))') 'join_ripples: dimensions=', ndim, ndim1
        write(*, '(a,i0,a,i0)') 'join_ripples: finite matrix entries=', &
            count(ieee_is_finite(matrix)), '/', size(matrix)
        write(*, '(a,i0,a,i0)') 'join_ripples: finite rhs entries=', &
            count(ieee_is_finite(rhs)), '/', size(rhs)
        ierr = stage
    end subroutine report_join_failure

    subroutine validate_join_normalization(facnorm, direction, column, old_tags, &
            new_tags, npass_l, npass_r, ierr)
        real(real64), intent(in) :: facnorm
        character(len=1), intent(in) :: direction
        integer, intent(in) :: column, old_tags(2), new_tags(2), npass_l, npass_r
        integer, intent(out) :: ierr

        if (ieee_is_finite(facnorm) .and. facnorm /= 0.0d0) then
            call record_join_normalization(facnorm, direction, column, &
                old_tags, new_tags, npass_l, npass_r, ierr)
            return
        end if

        write(*, '(a,a)') 'join_ripples: normalization direction=', direction
        write(*, '(a,i0)') 'join_ripples: normalization column=', column
        write(*, '(a,es23.15)') 'join_ripples: normalization factor=', facnorm
        write(*, '(a,2(1x,i0))') 'join_ripples: old propagator tags=', old_tags
        write(*, '(a,2(1x,i0))') 'join_ripples: new propagator tags=', new_tags
        write(*, '(a,2(1x,i0))') 'join_ripples: pass dimensions=', npass_l, npass_r
        ierr = 5
    end subroutine validate_join_normalization

    subroutine record_join_normalization(facnorm, direction, column, old_tags, &
            new_tags, npass_l, npass_r, ierr)
        real(real64), intent(in) :: facnorm
        character(len=1), intent(in) :: direction
        integer, intent(in) :: column, old_tags(2), new_tags(2), npass_l, npass_r
        integer, intent(out) :: ierr

        integer :: iunit, status

        ierr = 0
        if (.not. normalization_output_initialized) &
            call initialize_normalization_output(ierr)
        if (ierr /= 0 .or. .not. allocated(normalization_output_filename)) return
        open(newunit=iunit, file=normalization_output_filename, status='old', &
            position='append', action='write', iostat=status)
        if (status /= 0) then
            write(*, '(a,i0)') 'join_ripples: normalization output error=', status
            ierr = 6
            return
        end if
        normalization_sequence = normalization_sequence + 1
        write(iunit, &
            '(i0,",",a,",",7(i0,","),es25.16e3)', iostat=status) &
            normalization_sequence, direction, column, old_tags, new_tags, &
            npass_l, npass_r, facnorm
        close(iunit)
        if (status /= 0) then
            write(*, '(a,i0)') 'join_ripples: normalization output error=', status
            ierr = 6
        end if
    end subroutine record_join_normalization

    subroutine initialize_normalization_output(ierr)
        integer, intent(out) :: ierr

        integer :: iunit, length, status

        ierr = 0
        normalization_output_initialized = .true.
        call get_environment_variable('NEO2_JOIN_NORMALIZATION_FILE', &
            length=length, status=status)
        if (status /= 0 .or. length == 0) return
        allocate(character(len=length) :: normalization_output_filename)
        call get_environment_variable('NEO2_JOIN_NORMALIZATION_FILE', &
            value=normalization_output_filename, status=status)
        if (status == 0) then
            open(newunit=iunit, file=normalization_output_filename, &
                status='replace', action='write', iostat=status)
        end if
        if (status == 0) then
            write(iunit, '(a)', iostat=status) &
                'sequence,direction,column,old_start,old_end,new_start,' // &
                'new_end,npass_l,npass_r,factor'
            close(iunit)
        end if
        if (status /= 0) then
            write(*, '(a,i0)') 'join_ripples: normalization output error=', status
            if (allocated(normalization_output_filename)) &
                deallocate(normalization_output_filename)
            ierr = 6
        end if
    end subroutine initialize_normalization_output

    logical function join_end_trace_enabled(ierr)
        integer, intent(out) :: ierr

        ierr = 0
        if (.not. join_end_output_initialized) call initialize_join_end_output(ierr)
        join_end_trace_enabled = allocated(join_end_output_filename)
    end function join_end_trace_enabled

    subroutine record_join_end_compatibility(source_factor, source_before, &
            source_after, source_scale, measure_sum, compatibility, &
            dropped_p, dropped_m, &
            compatibility_scale, dropped_p_scale, dropped_m_scale, left_null, &
            right_null, left_residual, right_residual, left_scale, right_scale, &
            transfer_error, ierr)
        real(real64), intent(in) :: source_factor(3), source_before(3)
        real(real64), intent(in) :: source_after(3), source_scale(3), measure_sum
        real(real64), intent(in) :: compatibility(:, :), dropped_p(:, :)
        real(real64), intent(in) :: dropped_m(:, :), left_null(:, :)
        real(real64), intent(in) :: compatibility_scale(:, :)
        real(real64), intent(in) :: dropped_p_scale(:, :), dropped_m_scale(:, :)
        real(real64), intent(in) :: right_null(:, :), left_residual(:)
        real(real64), intent(in) :: right_residual(:), left_scale(:)
        real(real64), intent(in) :: right_scale(:), transfer_error(2)
        integer, intent(out) :: ierr

        integer :: force, index, iunit, laguerre, status

        ierr = 0
        if (.not. join_end_trace_enabled(ierr)) return
        if (ierr /= 0) return
        if (size(compatibility, 2) /= 3 &
            .or. any(shape(dropped_p) /= shape(compatibility)) &
            .or. any(shape(dropped_m) /= shape(compatibility)) &
            .or. any(shape(compatibility_scale) /= shape(compatibility)) &
            .or. any(shape(dropped_p_scale) /= shape(compatibility)) &
            .or. any(shape(dropped_m_scale) /= shape(compatibility)) &
            .or. size(left_null, 2) /= size(compatibility, 1) &
            .or. any(shape(right_null) /= shape(left_null)) &
            .or. size(left_residual) /= size(compatibility, 1) &
            .or. size(right_residual) /= size(compatibility, 1) &
            .or. size(left_scale) /= size(compatibility, 1) &
            .or. size(right_scale) /= size(compatibility, 1)) then
            ierr = 7
            return
        end if
        if (.not. all(ieee_is_finite(source_factor)) &
            .or. .not. all(ieee_is_finite(source_before)) &
            .or. .not. all(ieee_is_finite(source_after)) &
            .or. .not. all(ieee_is_finite(source_scale)) &
            .or. .not. ieee_is_finite(measure_sum) &
            .or. .not. all(ieee_is_finite(compatibility)) &
            .or. .not. all(ieee_is_finite(dropped_p)) &
            .or. .not. all(ieee_is_finite(dropped_m)) &
            .or. .not. all(ieee_is_finite(compatibility_scale)) &
            .or. .not. all(ieee_is_finite(dropped_p_scale)) &
            .or. .not. all(ieee_is_finite(dropped_m_scale)) &
            .or. .not. all(ieee_is_finite(left_null)) &
            .or. .not. all(ieee_is_finite(right_null)) &
            .or. .not. all(ieee_is_finite(left_residual)) &
            .or. .not. all(ieee_is_finite(right_residual)) &
            .or. .not. all(ieee_is_finite(left_scale)) &
            .or. .not. all(ieee_is_finite(right_scale)) &
            .or. any(compatibility_scale <= 0.0_real64) &
            .or. any(dropped_p_scale <= 0.0_real64) &
            .or. any(dropped_m_scale <= 0.0_real64) &
            .or. any(left_scale <= 0.0_real64) &
            .or. any(right_scale <= 0.0_real64) &
            .or. any(source_scale <= 0.0_real64) &
            .or. .not. all(ieee_is_finite(transfer_error))) then
            ierr = 7
            return
        end if

        open(newunit=iunit, file=join_end_output_filename, status='old', &
            position='append', action='write', iostat=status)
        if (status /= 0) then
            ierr = 8
            return
        end if
        join_end_record = join_end_record + 1
        do force = 1, 3
            call write_join_end_value(iunit, 'source_factor', -1, force, -1, &
                source_factor(force), status)
            call write_join_end_value(iunit, 'source_before', -1, force, -1, &
                source_before(force), status)
            call write_join_end_value(iunit, 'source_after', -1, force, -1, &
                source_after(force), status)
            call write_join_end_value(iunit, 'source_scale', -1, force, -1, &
                source_scale(force), status)
        end do
        call write_join_end_value(iunit, 'measure_sum', -1, 0, -1, &
            measure_sum, status)
        call write_join_end_value(iunit, 'forward_identity', -1, 0, -1, &
            transfer_error(1), status)
        call write_join_end_value(iunit, 'backward_identity', -1, 0, -1, &
            transfer_error(2), status)
        do laguerre = 1, size(compatibility, 1)
            do force = 1, 3
                call write_join_end_value(iunit, 'compatibility', &
                    laguerre - 1, force, -1, &
                    compatibility(laguerre, force), status)
                call write_join_end_value(iunit, 'compatibility_scale', &
                    laguerre - 1, force, -1, &
                    compatibility_scale(laguerre, force), status)
                call write_join_end_value(iunit, 'dropped_p', laguerre - 1, &
                    force, -1, dropped_p(laguerre, force), status)
                call write_join_end_value(iunit, 'dropped_p_scale', &
                    laguerre - 1, force, -1, &
                    dropped_p_scale(laguerre, force), status)
                call write_join_end_value(iunit, 'dropped_m', laguerre - 1, &
                    force, -1, dropped_m(laguerre, force), status)
                call write_join_end_value(iunit, 'dropped_m_scale', &
                    laguerre - 1, force, -1, &
                    dropped_m_scale(laguerre, force), status)
                call write_join_end_value(iunit, 'dropped_difference', &
                    laguerre - 1, force, -1, &
                    dropped_p(laguerre, force) &
                    - dropped_m(laguerre, force), status)
            end do
            call write_join_end_value(iunit, 'left_residual', laguerre - 1, &
                0, -1, left_residual(laguerre), status)
            call write_join_end_value(iunit, 'left_scale', laguerre - 1, &
                0, -1, left_scale(laguerre), status)
            call write_join_end_value(iunit, 'right_residual', laguerre - 1, &
                0, -1, right_residual(laguerre), status)
            call write_join_end_value(iunit, 'right_scale', laguerre - 1, &
                0, -1, right_scale(laguerre), status)
            do index = 1, size(left_null, 1)
                call write_join_end_value(iunit, 'left_null', laguerre - 1, &
                    0, index, left_null(index, laguerre), status)
                call write_join_end_value(iunit, 'right_null', laguerre - 1, &
                    0, index, right_null(index, laguerre), status)
            end do
        end do
        close(iunit, iostat=ierr)
        if (status /= 0 .or. ierr /= 0) ierr = 8
    end subroutine record_join_end_compatibility

    subroutine write_join_end_value(iunit, kind, laguerre, force, index, value, &
            status)
        integer, intent(in) :: iunit, laguerre, force, index
        character(len=*), intent(in) :: kind
        real(real64), intent(in) :: value
        integer, intent(inout) :: status

        if (status /= 0) return
        join_end_sequence = join_end_sequence + 1
        write(iunit, '(2(i0,","),a,",",3(i0,","),es25.16e3)', iostat=status) &
            join_end_sequence, join_end_record, trim(kind), laguerre, force, &
            index, value
    end subroutine write_join_end_value

    subroutine initialize_join_end_output(ierr)
        integer, intent(out) :: ierr

        integer :: iunit, length, status

        ierr = 0
        join_end_output_initialized = .true.
        call get_environment_variable('NEO2_JOIN_END_TRACE_FILE', &
            length=length, status=status)
        if (status /= 0 .or. length == 0) return
        allocate(character(len=length) :: join_end_output_filename)
        call get_environment_variable('NEO2_JOIN_END_TRACE_FILE', &
            value=join_end_output_filename, status=status)
        if (status == 0) then
            open(newunit=iunit, file=join_end_output_filename, &
                status='replace', action='write', iostat=status)
        end if
        if (status == 0) then
            write(iunit, '(a)', iostat=status) &
                'sequence,join,kind,laguerre,force,index,value'
            close(iunit, iostat=ierr)
        end if
        if (status /= 0) ierr = 8
        if (ierr /= 0 .and. allocated(join_end_output_filename)) &
            deallocate(join_end_output_filename)
    end subroutine initialize_join_end_output
end module join_diagnostics_mod
