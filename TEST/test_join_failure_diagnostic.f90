program test_join_failure_diagnostic
    use join_diagnostics_mod, only: report_join_failure, &
        validate_join_normalization
    use lapack_band, only: gbsv
    implicit none

    real(kind(1.0d0)) :: matrix(2, 2), rhs(2, 1)
    real(kind(1.0d0)) :: factor
    character(len=1024) :: filename
    character(len=128) :: header
    character(len=1) :: direction
    integer :: column, info, ierr, iunit, new_tags(2), old_tags(2), &
        npass(2), pivot(2), sequence, status

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
end program test_join_failure_diagnostic
