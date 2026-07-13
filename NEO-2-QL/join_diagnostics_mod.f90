module join_diagnostics_mod
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private

    character(len=:), allocatable, save :: normalization_output_filename
    integer, save :: normalization_sequence = 0
    logical, save :: normalization_output_initialized = .false.

    public :: report_join_failure, validate_join_normalization

contains

    subroutine report_join_failure(stage, info, old_tags, new_tags, ndim, ndim1, &
            matrix, rhs, ierr)
        integer, intent(in) :: stage, info, old_tags(2), new_tags(2), ndim, ndim1
        real(kind(1.0d0)), intent(in) :: matrix(:, :), rhs(:, :)
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
        real(kind(1.0d0)), intent(in) :: facnorm
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
        real(kind(1.0d0)), intent(in) :: facnorm
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
end module join_diagnostics_mod
