program test_stage2_attribution
    use, intrinsic :: iso_fortran_env, only: real64
    use stage2_attribution_mod, only: apply_stage2_attribution, &
        stage2_attribution_apply_error, stage2_attribution_enabled
    implicit none

    character(len=1024) :: filename
    character(len=64) :: header
    character(len=1) :: side
    real(real64) :: rhs(14, 3), value
    integer :: band, force, ierr, ind_start(3), iunit, npl(3), row, status, tag
    integer :: line_count
    logical :: enabled

    ! The footprint file is written before the module reads it lazily.
    call get_environment_variable('NEO2_STAGE2_ATTRIBUTION_FILE', filename)
    if (len_trim(filename) == 0) &
        error stop 'FAIL: NEO2_STAGE2_ATTRIBUTION_FILE is not set'
    open(newunit=iunit, file=trim(filename), status='replace', action='write')
    write(iunit, '(a)') 'tag,side,force,band,value'
    write(iunit, '(a)') '7,p,1,1,0.25'
    write(iunit, '(a)') '7,p,1,2,0.5'
    write(iunit, '(a)') '7,m,2,1,0.125'
    write(iunit, '(a)') '7,m,2,2,0.0625'
    write(iunit, '(a)') '9,p,3,1,1.0'
    write(iunit, '(a)') '9,m,1,5,0.5'
    close(iunit)

    call stage2_attribution_enabled(enabled, ierr)
    if (ierr /= 0 .or. .not. enabled) &
        error stop 'FAIL: valid footprint file was not accepted'

    ! Local grid: three spatial steps, npl = [1, 2, 1], lag = 0, so the
    ! per-step block sizes are [4, 6, 4] and ind_start = [0, 4, 10].
    npl = [1, 2, 1]
    ind_start = [0, 4, 10]

    ! A tag without footprint rows must leave the solve untouched.
    rhs = 1.0_real64
    call apply_stage2_attribution(5, rhs, npl, ind_start, 1, 3, 0, ierr)
    if (ierr /= 0 .or. any(rhs /= 1.0_real64)) &
        error stop 'FAIL: tag without rows modified the right-hand side'

    call apply_stage2_attribution(7, rhs, npl, ind_start, 1, 3, 0, ierr)
    if (ierr /= 0) error stop 'FAIL: valid footprint application failed'
    ! side p, band b: row = ind_start(iend) + b
    if (rhs(11, 1) /= 0.75_real64 .or. rhs(12, 1) /= 0.5_real64) &
        error stop 'FAIL: co-passing correction landed on the wrong rows'
    ! side m, band b: row = ind_start(ibeg) + 2*(npl(ibeg)+1) + 1 - b
    if (rhs(4, 2) /= 0.875_real64 .or. rhs(3, 2) /= 0.9375_real64) &
        error stop 'FAIL: counter-passing correction landed on the wrong rows'
    rhs(11, 1) = 1.0_real64
    rhs(12, 1) = 1.0_real64
    rhs(4, 2) = 1.0_real64
    rhs(3, 2) = 1.0_real64
    if (any(rhs /= 1.0_real64)) &
        error stop 'FAIL: correction touched rows outside the footprint'

    ! Non-Lorentz state layouts are outside the join_ends correction scope.
    rhs = 1.0_real64
    call apply_stage2_attribution(7, rhs, npl, ind_start, 1, 3, 1, ierr)
    if (ierr /= stage2_attribution_apply_error) &
        error stop 'FAIL: lag > 0 state was accepted'

    ! A band beyond the local boundary grid must be a hard failure.
    rhs = 1.0_real64
    call apply_stage2_attribution(9, rhs, npl, ind_start, 1, 3, 0, ierr)
    if (ierr /= stage2_attribution_apply_error) &
        error stop 'FAIL: band outside the local grid was accepted'

    ! The applied trace records exactly what was subtracted.
    open(newunit=iunit, file=trim(filename)//'.applied', status='old', &
        action='read', iostat=status)
    if (status /= 0) error stop 'FAIL: applied trace is absent'
    read(iunit, '(a)', iostat=status) header
    if (status /= 0 .or. trim(header) /= 'tag,side,force,band,row,value') &
        error stop 'FAIL: applied trace header is wrong'
    line_count = 0
    do
        read(iunit, *, iostat=status) tag, side, force, band, row, value
        if (status < 0) exit
        if (status > 0) error stop 'FAIL: malformed applied trace row'
        line_count = line_count + 1
        if (line_count == 1) then
            if (tag /= 7 .or. side /= 'p' .or. force /= 1 .or. band /= 1 &
                .or. row /= 11 .or. value /= 0.25_real64) &
                error stop 'FAIL: applied trace record is wrong'
        end if
    end do
    close(iunit)
    ! Four rows from the accepted tag-7 application plus the valid tag-9
    ! co-passing row written before its counter-passing band was rejected.
    if (line_count /= 5) error stop 'FAIL: applied trace row count differs'

    write(*, '(a)') 'All tests passed!'
end program test_stage2_attribution
