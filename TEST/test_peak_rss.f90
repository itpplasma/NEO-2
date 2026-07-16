program test_peak_rss
    use, intrinsic :: iso_fortran_env, only: int64
    use peak_rss_mod, only: current_peak_rss_kib, record_peak_rss
    implicit none

    real(kind(1.0d0)), allocatable :: ballast(:)
    character(len=1024) :: filename
    character(len=64) :: header
    integer(int64) :: peak_after, peak_before, recorded
    integer :: ierr, iunit, rank, status

    call current_peak_rss_kib(peak_before, ierr)
    if (ierr /= 0 .or. peak_before <= 0_int64) &
        error stop 'FAIL: VmHWM was not read from /proc/self/status'

    ! Touch 64 MiB so the high-water mark provably moves.
    allocate(ballast(8*1024*1024))
    call random_number(ballast)
    call current_peak_rss_kib(peak_after, ierr)
    if (ierr /= 0) error stop 'FAIL: VmHWM was not re-read'
    if (peak_after < peak_before + 32768_int64) &
        error stop 'FAIL: peak RSS did not track the 64 MiB allocation'
    if (ballast(1) < 0.0d0) error stop 'FAIL: ballast values are invalid'

    call record_peak_rss(3, ierr)
    if (ierr /= 0) error stop 'FAIL: peak RSS record was not written'

    call get_environment_variable('NEO2_PEAK_RSS_FILE', filename)
    open(newunit=iunit, file=trim(filename), status='old', action='read', &
        iostat=status)
    if (status /= 0) error stop 'FAIL: peak RSS output is absent'
    read(iunit, '(a)', iostat=status) header
    if (status /= 0 .or. trim(header) /= 'rank,peak_rss_kib') &
        error stop 'FAIL: peak RSS header is wrong'
    read(iunit, *, iostat=status) rank, recorded
    close(iunit)
    if (status /= 0 .or. rank /= 3) error stop 'FAIL: peak RSS record is wrong'
    if (recorded < peak_after) &
        error stop 'FAIL: recorded peak is below the observed peak'

    write(*, '(a)') 'All tests passed!'
end program test_peak_rss
