program test_k_coefficient_only
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use hdf5_tools, only: HID_T, h5_init, h5_deinit, h5_open, h5_get, h5_exists, &
        h5_close
    implicit none

    integer(HID_T) :: file_id
    character(len=1024) :: work_dir
    character(len=16) :: suffix, full_suffix
    real(dp) :: expected(3), expected_b(3), actual(3), qflux(3, 3)
    logical :: k_only

    call get_command_argument(1, work_dir)
    call get_command_argument(2, suffix)
    call get_command_argument(3, full_suffix)
    call h5_init()
    call h5_open(trim(work_dir)//'/full'//trim(full_suffix)//'/fulltransp.h5', file_id)
    call h5_get(file_id, 'k_cof', expected(1))
    call h5_get(file_id, 'qflux', qflux)
    expected(2:3) = [qflux(2, 1), qflux(2, 3)]
    if (h5_exists(file_id, 'k_coefficient_only')) &
        error stop 'Default transport schema changed.'
    call h5_close(file_id)

    call h5_open(trim(work_dir)//'/only'//trim(suffix)//'/fulltransp.h5', file_id)
    call h5_get(file_id, 'k_cof', actual(1))
    call h5_get(file_id, 'qcurr_g', actual(2))
    call h5_get(file_id, 'qcurr_t', actual(3))
    call h5_get(file_id, 'k_coefficient_only', k_only)
    if (.not. k_only) error stop 'Missing k-only output marker.'
    if (h5_exists(file_id, 'qflux') .or. h5_exists(file_id, 'gamma_out') .or. &
        h5_exists(file_id, 'full_version')) &
        error stop 'Incomplete transport matrix was published.'
    call h5_close(file_id)

    if (.not. all(ieee_is_finite(expected)) .or. &
        .not. all(ieee_is_finite(actual))) error stop 'Nonfinite transport result.'
    if (abs(expected(2)) <= tiny(1.0_dp)) error stop 'Zero current response.'
    if (any(actual /= expected)) then
        print *, 'Full result: ', expected
        print *, 'K-only result: ', actual
        error stop 'K-only changed the current responses or k coefficient.'
    end if

    call h5_open(trim(work_dir)//'/second/fulltransp.h5', file_id)
    call h5_get(file_id, 'k_cof', expected_b(1))
    call h5_get(file_id, 'qcurr_g', expected_b(2))
    call h5_get(file_id, 'qcurr_t', expected_b(3))
    call h5_close(file_id)
    if (.not. all(ieee_is_finite(expected_b))) error stop 'Nonfinite second result.'
    if (abs(expected_b(2)) <= tiny(1.0_dp)) error stop 'Zero second current response.'
    call check_batch('aba'//trim(suffix)//'/k_coefficients.h5', [1, 2, 1], 1)
    call check_batch('bab'//trim(suffix)//'/k_coefficients.h5', [2, 1, 2], 1)
    call check_batch('chunk'//trim(suffix)//'/chunk.h5', [2], 2)
    call check_batch('laststep'//trim(suffix)//'/k_coefficients.h5', [1], 1)
    call check_batch('exhausted'//trim(suffix)//'/k_coefficients.h5', [1], 1, 0)
    call h5_deinit()
    print *, 'All tests passed!'

contains

    subroutine check_batch(name, identities, first_expected, convergence_expected)
        character(*), intent(in) :: name
        integer, intent(in) :: identities(:), first_expected
        integer, optional, intent(in) :: convergence_expected
        real(dp) :: results(3, size(identities)), s(size(identities))
        real(dp) :: r(size(identities)), z(size(identities)), nu(size(identities))
        real(dp), parameter :: expected_s(2) = [0.5_dp, 0.49_dp]
        real(dp), parameter :: expected_nu(2) = [0.01_dp, 0.02_dp]
        integer :: converged(size(identities)), first, last, row, expected_flag

        expected_flag = 1
        if (present(convergence_expected)) expected_flag = convergence_expected
        call h5_open(trim(work_dir)//'/'//name, file_id)
        call h5_get(file_id, 'k_coefficient_only', k_only)
        call h5_get(file_id, 'surface_first', first)
        call h5_get(file_id, 'surface_last', last)
        call h5_get(file_id, 'k_cof', results(1, :))
        call h5_get(file_id, 'qcurr_g', results(2, :))
        call h5_get(file_id, 'qcurr_t', results(3, :))
        call h5_get(file_id, 'boozer_s', s)
        call h5_get(file_id, 'rbeg', r)
        call h5_get(file_id, 'zbeg', z)
        call h5_get(file_id, 'conl_over_mfp', nu)
        call h5_get(file_id, 'integral_converged', converged)
        call h5_close(file_id)
        if (.not. k_only) error stop 'Missing batch output marker.'
        if (first /= first_expected .or. last /= first + size(identities) - 1) &
            error stop 'Incorrect batch range.'
        if (.not. all(ieee_is_finite(results))) error stop 'Nonfinite batch result.'
        if (any(s /= expected_s(identities)) .or. any(nu /= expected_nu(identities)) &
            .or. any(r /= 0.0_dp) .or. any(z /= 0.0_dp)) &
            error stop 'Incorrect batch surface identities.'
        if (any(converged /= expected_flag)) error stop 'Incorrect convergence status.'
        ! An exhausted iterate is finite but differs from the converged reference.
        if (expected_flag == 0) return
        do row = 1, size(identities)
            if (identities(row) == 1) then
                actual = expected
            else
                actual = expected_b
            end if
            if (any(results(:, row) /= actual)) then
                print *, name, ' row ', row, ': ', results(:, row), ' expected ', actual
                error stop 'Batch changed the current responses or k coefficient.'
            end if
        end do
    end subroutine check_batch
end program test_k_coefficient_only
