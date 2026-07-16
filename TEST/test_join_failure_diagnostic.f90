program test_join_failure_diagnostic
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use join_diagnostics_mod, only: record_join_end_compatibility, &
        report_join_failure, &
        validate_join_normalization
    use lapack_band, only: gbsv
    implicit none

    real(kind(1.0d0)) :: bad_dropped(2, 3), compatibility(1, 3)
    real(kind(1.0d0)) :: compatibility_scale(1, 3)
    real(kind(1.0d0)) :: delta_eta_l(2), delta_eta_r(2)
    real(kind(1.0d0)) :: dropped_m(1, 3), dropped_p(1, 3), factor
    real(kind(1.0d0)) :: dropped_m_scale(1, 3), dropped_p_scale(1, 3)
    real(kind(1.0d0)) :: left_null(2, 1), left_residual(1), left_scale(1)
    real(kind(1.0d0)) :: matrix(2, 2), measure_sum, rhs(2, 1)
    real(kind(1.0d0)) :: right_null(2, 1), right_residual(1), right_scale(1)
    real(kind(1.0d0)) :: source_after(3), source_before(3)
    real(kind(1.0d0)) :: source_factor(3), source_scale(3), transfer_error(2)
    real(kind(1.0d0)) :: value
    character(len=1024) :: filename
    character(len=128) :: header
    character(len=32) :: record_kind
    character(len=1) :: direction
    integer :: column, force, index, info, ierr, iunit, join, laguerre, line_count
    integer :: new_tags(2), old_tags(2), npass(2), pivot(2), sequence, status
    logical :: found_correction_p, found_delta_eta_l, found_right_null
    logical :: found_source_factor

    matrix = reshape([1.0d0, 2.0d0, 2.0d0, 4.0d0], shape(matrix))
    rhs(:, 1) = [1.0d0, 2.0d0]
    call gbsv(2, 2, matrix, pivot, rhs, info)

    if (info /= 2) error stop 'FAIL: singular solve did not return pivot 2'
    call report_join_failure(2, info, [10, 20], [21, 30], 2, 2, &
        matrix, rhs, ierr)
    if (ierr /= 2) error stop 'FAIL: join stage was not retained'

    call validate_join_normalization(0.0d0, 'p', 2, [12, 13], [13, 13], &
        24, 24, ierr)
    if (ierr /= 5) error stop 'FAIL: zero normalization was not rejected'

    call validate_join_normalization(1.0d0, 'p', 2, [12, 13], [13, 13], &
        24, 24, ierr)
    if (ierr /= 0) error stop 'FAIL: finite normalization was rejected'

    call get_environment_variable('NEO2_JOIN_NORMALIZATION_FILE', filename)
    open(newunit=iunit, file=trim(filename), status='old', action='read', &
        iostat=status)
    if (status /= 0) error stop 'FAIL: join normalization output is absent'
    read(iunit, '(a)', iostat=status) header
    if (status /= 0) error stop 'FAIL: join normalization header is absent'
    read(iunit, *, iostat=status) sequence, direction, column, old_tags, &
        new_tags, npass, factor
    close(iunit)
    if (status /= 0 .or. sequence /= 1 .or. direction /= 'p' &
        .or. column /= 2 .or. any(old_tags /= [12, 13]) &
        .or. any(new_tags /= [13, 13]) .or. any(npass /= [24, 24]) &
        .or. factor /= 1.0d0) &
        error stop 'FAIL: join normalization record is wrong'

    source_factor = [1.0d0, 2.0d0, 3.0d0]
    source_before = [4.0d0, 5.0d0, 6.0d0]
    source_after = [0.0d0, 0.0d0, 0.0d0]
    source_scale = [25.0d0, 26.0d0, 27.0d0]
    measure_sum = 2.0d0
    delta_eta_l = [0.5d0, 1.5d0]
    delta_eta_r = [2.5d0, 3.5d0]
    compatibility = reshape([7.0d0, 8.0d0, 9.0d0], [1, 3])
    dropped_p = reshape([10.0d0, 11.0d0, 12.0d0], [1, 3])
    dropped_m = reshape([13.0d0, 14.0d0, 15.0d0], [1, 3])
    compatibility_scale = 20.0d0
    dropped_p_scale = 21.0d0
    dropped_m_scale = 22.0d0
    left_null(:, 1) = [1.0d0, 1.0d0]
    right_null(:, 1) = [3.0d0, 4.0d0]
    left_residual = [16.0d0]
    right_residual = [17.0d0]
    left_scale = [23.0d0]
    right_scale = [24.0d0]
    transfer_error = [18.0d0, 19.0d0]
    call record_join_end_compatibility(source_factor, source_before, &
        source_after, source_scale, measure_sum, delta_eta_l, delta_eta_r, &
        compatibility, dropped_p, dropped_m, &
        compatibility_scale, dropped_p_scale, dropped_m_scale, left_null, &
        right_null, left_residual, right_residual, left_scale, right_scale, &
        transfer_error, ierr)
    if (ierr /= 0) error stop 'FAIL: valid join-end trace was rejected'

    bad_dropped = 0.0d0
    call record_join_end_compatibility(source_factor, source_before, &
        source_after, source_scale, measure_sum, delta_eta_l, delta_eta_r, &
        compatibility, bad_dropped, dropped_m, &
        compatibility_scale, dropped_p_scale, dropped_m_scale, left_null, &
        right_null, left_residual, right_residual, left_scale, right_scale, &
        transfer_error, ierr)
    if (ierr /= 7) error stop 'FAIL: invalid join-end shape was accepted'

    source_scale(1) = ieee_value(0.0d0, ieee_quiet_nan)
    call record_join_end_compatibility(source_factor, source_before, &
        source_after, source_scale, measure_sum, delta_eta_l, delta_eta_r, &
        compatibility, dropped_p, &
        dropped_m, compatibility_scale, dropped_p_scale, dropped_m_scale, &
        left_null, right_null, left_residual, right_residual, left_scale, &
        right_scale, transfer_error, ierr)
    if (ierr /= 7) error stop 'FAIL: non-finite join-end value was accepted'
    source_scale(1) = 25.0d0

    delta_eta_r(2) = ieee_value(0.0d0, ieee_quiet_nan)
    call record_join_end_compatibility(source_factor, source_before, &
        source_after, source_scale, measure_sum, delta_eta_l, delta_eta_r, &
        compatibility, dropped_p, &
        dropped_m, compatibility_scale, dropped_p_scale, dropped_m_scale, &
        left_null, right_null, left_residual, right_residual, left_scale, &
        right_scale, transfer_error, ierr)
    if (ierr /= 7) error stop 'FAIL: non-finite band correction was accepted'
    delta_eta_r(2) = 3.5d0

    call get_environment_variable('NEO2_JOIN_END_TRACE_FILE', filename)
    open(newunit=iunit, file=trim(filename), status='old', action='read', &
        iostat=status)
    if (status /= 0) error stop 'FAIL: join-end trace is absent'
    read(iunit, '(a)', iostat=status) header
    if (status /= 0 .or. trim(header) /= &
        'sequence,join,kind,laguerre,force,index,value') &
        error stop 'FAIL: join-end trace header is wrong'
    line_count = 0
    found_correction_p = .false.
    found_delta_eta_l = .false.
    found_right_null = .false.
    found_source_factor = .false.
    do
        read(iunit, *, iostat=status) sequence, join, record_kind, laguerre, &
            force, index, value
        if (status < 0) exit
        if (status > 0) error stop 'FAIL: malformed join-end trace row'
        line_count = line_count + 1
        if (sequence /= line_count) &
            error stop 'FAIL: join-end trace sequence is wrong'
        if (join /= 1) error stop 'FAIL: join-end trace join is wrong'
        if (trim(record_kind) == 'source_factor' .and. force == 2) then
            found_source_factor = value == 2.0d0
        end if
        if (trim(record_kind) == 'right_null' .and. index == 2) then
            found_right_null = value == 4.0d0
        end if
        if (trim(record_kind) == 'delta_eta_l' .and. index == 2) then
            found_delta_eta_l = value == 1.5d0
        end if
        if (trim(record_kind) == 'source_correction_p' .and. force == 2 &
            .and. index == 1) then
            ! delta_eta_r(1)*source_factor(2) = 2.5*2.0
            found_correction_p = value == 5.0d0
        end if
    end do
    close(iunit)
    if (line_count /= 60 .or. .not. found_source_factor &
        .or. .not. found_right_null .or. .not. found_delta_eta_l &
        .or. .not. found_correction_p) &
        error stop 'FAIL: join-end trace content is wrong'
end program test_join_failure_diagnostic
