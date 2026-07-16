module stage2_attribution_mod
    !> Env-guarded stage-2 attribution experiment for the periodic-closure
    !> interface jump.
    !>
    !> When NEO2_STAGE2_ATTRIBUTION_FILE names a footprint file, the stage-2
    !> local re-solve subtracts the listed per-band source corrections from
    !> the assembled right-hand side before sparse_solve.  The footprint rows
    !> are the source_correction_p/source_correction_m values exported by the
    !> stage-0 join-end trace (join_diagnostics_mod): the per-band, per-force
    !> Lorentz solvability corrections delta_eta*facnorm that join_ends
    !> subtracts from the joined boundary sources.  The band vectors are
    !> mapped onto the local propagator phase-space grid, matching the
    !> source_p/source_m extraction order of ripple_solver_axi_test.f90:
    !>
    !>   side p: co-passing block at the last spatial step,
    !>           row = ind_start(iend) + band,          band = 1..npl(iend)+1
    !>   side m: counter-passing block at the first spatial step,
    !>           row = ind_start(ibeg) + 2*(npl(ibeg)+1) + 1 - band,
    !>                                                   band = 1..npl(ibeg)+1
    !>
    !> Only lag = 0 (Lorentz) states are accepted, matching the scope of the
    !> join_ends correction.  This mode measures the remainder of the jump:
    !> it does not assume the footprint identity holds.  The subtracted rows
    !> are appended to <file>.applied so the experiment records exactly what
    !> was removed; the remaining jump is measured by the unchanged interface
    !> diagnostics downstream.
    !>
    !> File format: one header line, then CSV rows
    !>   tag,side,force,band,value
    !> with side one of p|m, force in 1..3, and band 1-based.  With the
    !> environment variable unset this module leaves the solve untouched.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
    private

    !> Distinct hard-failure codes: never 3, which propagator.f90 treats as
    !> a retryable phi-refinement request.
    integer, parameter, public :: stage2_attribution_input_error = 11
    integer, parameter, public :: stage2_attribution_apply_error = 12

    character(len=:), allocatable, save :: footprint_filename
    integer, allocatable, save :: row_tag(:), row_force(:), row_band(:)
    character(len=1), allocatable, save :: row_side(:)
    real(real64), allocatable, save :: row_value(:)
    logical, save :: footprint_initialized = .false.
    logical, save :: applied_output_initialized = .false.

    public :: stage2_attribution_enabled, apply_stage2_attribution

contains

    subroutine stage2_attribution_enabled(enabled, ierr)
        logical, intent(out) :: enabled
        integer, intent(out) :: ierr

        ierr = 0
        if (.not. footprint_initialized) call initialize_footprint(ierr)
        enabled = ierr == 0 .and. allocated(footprint_filename)
    end subroutine stage2_attribution_enabled

    !> Subtract the footprint rows recorded for this propagator tag from the
    !> assembled right-hand side.  Rows for other tags are ignored; a tag
    !> without rows leaves the solve untouched.
    subroutine apply_stage2_attribution(tag, rhs, npl, ind_start, ibeg, iend, &
            lag, ierr)
        integer, intent(in) :: tag
        real(real64), intent(inout) :: rhs(:, :)
        integer, intent(in) :: npl(ibeg:), ind_start(ibeg:)
        integer, intent(in) :: ibeg, iend, lag
        integer, intent(out) :: ierr

        integer :: applied_m, applied_p, entry, iunit, row, status
        real(real64) :: l1_m, l1_p

        ierr = 0
        if (.not. allocated(footprint_filename)) return
        if (size(row_tag) == 0) return
        if (all(row_tag /= tag)) return
        if (lag /= 0 .or. ibeg > iend .or. size(rhs, 2) /= 3) then
            write (*, '(a,i0)') &
                'stage2 attribution: unsupported state layout for tag=', tag
            ierr = stage2_attribution_apply_error
            return
        end if

        call initialize_applied_output(status)
        if (status /= 0) then
            ierr = stage2_attribution_input_error
            return
        end if
        open(newunit=iunit, file=footprint_filename//'.applied', &
            status='old', position='append', action='write', iostat=status)
        if (status /= 0) then
            ierr = stage2_attribution_input_error
            return
        end if

        applied_p = 0
        applied_m = 0
        l1_p = 0.0_real64
        l1_m = 0.0_real64
        do entry = 1, size(row_tag)
            if (row_tag(entry) /= tag) cycle
            if (row_side(entry) == 'p') then
                if (row_band(entry) > npl(iend) + 1) then
                    call report_band_mismatch(tag, entry, npl(iend) + 1)
                    ierr = stage2_attribution_apply_error
                    exit
                end if
                row = ind_start(iend) + row_band(entry)
                applied_p = applied_p + 1
                l1_p = l1_p + abs(row_value(entry))
            else
                if (row_band(entry) > npl(ibeg) + 1) then
                    call report_band_mismatch(tag, entry, npl(ibeg) + 1)
                    ierr = stage2_attribution_apply_error
                    exit
                end if
                row = ind_start(ibeg) + 2*(npl(ibeg) + 1) + 1 - row_band(entry)
                applied_m = applied_m + 1
                l1_m = l1_m + abs(row_value(entry))
            end if
            if (row < 1 .or. row > size(rhs, 1)) then
                call report_band_mismatch(tag, entry, size(rhs, 1))
                ierr = stage2_attribution_apply_error
                exit
            end if
            rhs(row, row_force(entry)) = rhs(row, row_force(entry)) &
                - row_value(entry)
            write(iunit, '(i0,",",a,",",3(i0,","),es25.16e3)', iostat=status) &
                tag, row_side(entry), row_force(entry), row_band(entry), row, &
                row_value(entry)
            if (status /= 0) then
                ierr = stage2_attribution_input_error
                exit
            end if
        end do
        close(iunit, iostat=status)
        if (ierr == 0 .and. status /= 0) ierr = stage2_attribution_input_error
        if (ierr /= 0) return
        write (*, '(a,i0,a,i0,a,es13.6,a,i0,a,es13.6)') &
            'stage2 attribution: tag=', tag, ' rows_p=', applied_p, &
            ' l1_p=', l1_p, ' rows_m=', applied_m, ' l1_m=', l1_m
    end subroutine apply_stage2_attribution

    subroutine report_band_mismatch(tag, entry, bound)
        integer, intent(in) :: tag, entry, bound

        write (*, '(3(a,i0),a,i0)') &
            'stage2 attribution: band outside local grid for tag=', tag, &
            ' entry=', entry, ' band=', row_band(entry), ' bound=', bound
    end subroutine report_band_mismatch

    subroutine initialize_footprint(ierr)
        integer, intent(out) :: ierr

        character(len=1024) :: line
        character(len=1) :: side
        integer :: band, force, iunit, length, n_rows, status, tag
        real(real64) :: value

        ierr = 0
        footprint_initialized = .true.
        call get_environment_variable('NEO2_STAGE2_ATTRIBUTION_FILE', &
            length=length, status=status)
        if (status /= 0 .or. length == 0) return
        allocate(character(len=length) :: footprint_filename)
        call get_environment_variable('NEO2_STAGE2_ATTRIBUTION_FILE', &
            value=footprint_filename, status=status)
        if (status == 0) then
            open(newunit=iunit, file=footprint_filename, status='old', &
                action='read', iostat=status)
        end if
        if (status /= 0) then
            call reject_footprint(0)
            ierr = stage2_attribution_input_error
            return
        end if
        read(iunit, '(a)', iostat=status) line
        if (status /= 0) then
            close(iunit)
            call reject_footprint(1)
            ierr = stage2_attribution_input_error
            return
        end if
        n_rows = 0
        do
            read(iunit, '(a)', iostat=status) line
            if (status /= 0) exit
            if (len_trim(line) == 0) cycle
            n_rows = n_rows + 1
        end do
        allocate(row_tag(n_rows), row_side(n_rows), row_force(n_rows))
        allocate(row_band(n_rows), row_value(n_rows))
        rewind(iunit)
        read(iunit, '(a)', iostat=status) line
        n_rows = 0
        do
            read(iunit, '(a)', iostat=status) line
            if (status < 0) exit
            if (status > 0) then
                close(iunit)
                call reject_footprint(n_rows + 2)
                ierr = stage2_attribution_input_error
                return
            end if
            if (len_trim(line) == 0) cycle
            read(line, *, iostat=status) tag, side, force, band, value
            if (status /= 0 .or. (side /= 'p' .and. side /= 'm') &
                .or. force < 1 .or. force > 3 .or. band < 1 &
                .or. .not. ieee_is_finite(value)) then
                close(iunit)
                call reject_footprint(n_rows + 2)
                ierr = stage2_attribution_input_error
                return
            end if
            n_rows = n_rows + 1
            row_tag(n_rows) = tag
            row_side(n_rows) = side
            row_force(n_rows) = force
            row_band(n_rows) = band
            row_value(n_rows) = value
        end do
        close(iunit)
    end subroutine initialize_footprint

    subroutine reject_footprint(line_number)
        integer, intent(in) :: line_number

        write (*, '(a,a,a,i0)') 'stage2 attribution: footprint file ', &
            'rejected, line ', 'number=', line_number
        if (allocated(footprint_filename)) deallocate(footprint_filename)
        if (allocated(row_tag)) deallocate(row_tag)
        if (allocated(row_side)) deallocate(row_side)
        if (allocated(row_force)) deallocate(row_force)
        if (allocated(row_band)) deallocate(row_band)
        if (allocated(row_value)) deallocate(row_value)
    end subroutine reject_footprint

    subroutine initialize_applied_output(status)
        integer, intent(out) :: status

        integer :: iunit

        status = 0
        if (applied_output_initialized) return
        open(newunit=iunit, file=footprint_filename//'.applied', &
            status='replace', action='write', iostat=status)
        if (status == 0) then
            write(iunit, '(a)', iostat=status) 'tag,side,force,band,row,value'
            close(iunit)
        end if
        if (status == 0) applied_output_initialized = .true.
    end subroutine initialize_applied_output
end module stage2_attribution_mod
