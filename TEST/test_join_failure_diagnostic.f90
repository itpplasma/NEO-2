program test_join_failure_diagnostic
    use join_diagnostics_mod, only: report_join_failure
    use lapack_band, only: gbsv
    implicit none

    real(kind(1.0d0)) :: matrix(2, 2), rhs(2, 1)
    integer :: info, ierr, pivot(2)

    matrix = reshape([1.0d0, 2.0d0, 2.0d0, 4.0d0], shape(matrix))
    rhs(:, 1) = [1.0d0, 2.0d0]
    call gbsv(2, 2, matrix, pivot, rhs, info)

    if (info /= 2) error stop 'FAIL: singular solve did not return pivot 2'
    call report_join_failure(2, info, [10, 20], [21, 30], 2, 2, &
        matrix, rhs, ierr)
    if (ierr /= 2) error stop 'FAIL: join stage was not retained'
end program test_join_failure_diagnostic
