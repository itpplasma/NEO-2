module peak_rss_mod
    !> Per-process peak resident-set-size accounting.
    !>
    !> current_peak_rss_kib reads VmHWM (the resident-set high-water mark of
    !> the calling process, in KiB) from /proc/self/status.  record_peak_rss
    !> writes it to the CSV file named by NEO2_PEAK_RSS_FILE, following the
    !> opt-in NEO2_*_FILE pattern of qflux_profile_mod; with the variable
    !> unset nothing is written.  The value is per process: under MPI each
    !> rank sees its own high-water mark, and the caller decides which ranks
    !> report.
    use, intrinsic :: iso_fortran_env, only: int64
    implicit none
    private

    public :: current_peak_rss_kib, record_peak_rss

contains

    subroutine current_peak_rss_kib(peak_kib, ierr)
        integer(int64), intent(out) :: peak_kib
        integer, intent(out) :: ierr

        character(len=256) :: line
        character(len=8) :: unit_token
        integer :: iunit, status

        peak_kib = -1_int64
        ierr = 1
        open(newunit=iunit, file='/proc/self/status', status='old', &
            action='read', iostat=status)
        if (status /= 0) return
        do
            read(iunit, '(a)', iostat=status) line
            if (status /= 0) exit
            if (line(1:6) /= 'VmHWM:') cycle
            read(line(7:), *, iostat=status) peak_kib, unit_token
            if (status == 0 .and. trim(unit_token) == 'kB' &
                .and. peak_kib > 0_int64) ierr = 0
            exit
        end do
        close(iunit)
        if (ierr /= 0) peak_kib = -1_int64
    end subroutine current_peak_rss_kib

    subroutine record_peak_rss(rank, ierr)
        integer, intent(in) :: rank
        integer, intent(out) :: ierr

        character(len=:), allocatable :: filename
        integer(int64) :: peak_kib
        integer :: iunit, length, status

        ierr = 0
        call get_environment_variable('NEO2_PEAK_RSS_FILE', length=length, &
            status=status)
        if (status /= 0 .or. length == 0) return
        allocate(character(len=length) :: filename)
        call get_environment_variable('NEO2_PEAK_RSS_FILE', value=filename, &
            status=status)
        if (status /= 0) then
            ierr = 2
            return
        end if
        call current_peak_rss_kib(peak_kib, ierr)
        if (ierr /= 0) return
        open(newunit=iunit, file=filename, status='replace', action='write', &
            iostat=status)
        if (status /= 0) then
            ierr = 2
            return
        end if
        write(iunit, '(a)', iostat=status) 'rank,peak_rss_kib'
        if (status == 0) write(iunit, '(i0,",",i0)', iostat=status) rank, &
            peak_kib
        close(iunit, iostat=ierr)
        if (status /= 0 .or. ierr /= 0) ierr = 2
    end subroutine record_peak_rss
end module peak_rss_mod
