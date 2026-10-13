program test_k_coefficient_only
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use hdf5_tools, only: HID_T, h5_init, h5_deinit, h5_open, h5_get, h5_exists, &
        h5_close
    implicit none

    integer(HID_T) :: file_id
    character(len=1024) :: work_dir
    character(len=16) :: only_case
    real(dp) :: expected(3), expected_b(3), actual(3), qflux(3, 3)
    real(dp) :: expected_gmres(3), expected_gmres_b(3)
    logical :: k_only

    call get_command_argument(1, work_dir)
    call get_command_argument(2, only_case)
    if (len_trim(only_case) == 0) only_case = 'only'
    call h5_init()
    call h5_open(trim(work_dir)//'/full/fulltransp.h5', file_id)
    call h5_get(file_id, 'k_cof', expected(1))
    call h5_get(file_id, 'qflux', qflux)
    expected(2:3) = [qflux(2, 1), qflux(2, 3)]
    if (h5_exists(file_id, 'k_coefficient_only')) &
        error stop 'Default transport schema changed.'
    call h5_close(file_id)

    call h5_open(trim(work_dir)//'/'//trim(only_case)//'/fulltransp.h5', file_id)
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
    call check_batch('aba/k_coefficients.h5', [1, 2, 1], 1)
    call check_batch('bab/k_coefficients.h5', [2, 1, 2], 1)
    call check_batch('chunk/chunk.h5', [2], 2)
    call check_batch('batch_laststep/k_coefficients.h5', [1], 1)
    call check_batch('batch_exhausted/k_coefficients.h5', [1], 1, &
                     expected_convergence=0)
    call check_gmres('gmres', 3, expected_gmres, final_application=3)
    call check_gmres('gmres_second', 100, expected_gmres_b)
    call check_batch('gmres_aba/k_coefficients.h5', [1, 2, 1], 1, &
                     expected_gmres, expected_gmres_b)
    call check_batch('gmres_laststep/k_coefficients.h5', [1], 1, &
                     expected_gmres, expected_gmres_b)
    call h5_deinit()
    print *, 'All tests passed!'

contains

    subroutine check_gmres(name, max_applications, results, final_application)
        character(*), intent(in) :: name
        integer, intent(in) :: max_applications
        real(dp), intent(out) :: results(3)
        integer, intent(in), optional :: final_application
        real(dp) :: residual(2)
        integer :: solver, source_ids(2), converged(2), applications(2)

        call h5_open(trim(work_dir)//'/'//name//'/fulltransp.h5', file_id)
        call h5_get(file_id, 'k_cof', results(1))
        call h5_get(file_id, 'qcurr_g', results(2))
        call h5_get(file_id, 'qcurr_t', results(3))
        call h5_get(file_id, 'k_iter_solver', solver)
        call h5_get(file_id, 'k_source_ids', source_ids)
        call h5_get(file_id, 'k_source_converged', converged)
        call h5_get(file_id, 'k_source_applications', applications)
        call h5_get(file_id, 'k_source_residual', residual)
        call h5_close(file_id)
        if (.not. all(ieee_is_finite(results)) .or. &
            .not. all(ieee_is_finite(residual))) error stop 'Nonfinite GMRES result.'
        if (solver /= 1 .or. any(source_ids /= [1,3]) .or. any(converged /= 1)) &
            error stop 'Missing GMRES source convergence.'
        if (any(residual > 5.d-2) .or. any(residual < 0.d0)) &
            error stop 'GMRES residual exceeds the requested tolerance.'
        if (any(applications < 2) .or. any(applications > max_applications)) &
            error stop 'GMRES exceeded the application budget.'
        if (present(final_application)) then
            if (any(applications /= final_application)) &
                error stop 'GMRES did not converge on the final allowed application.'
        end if
        if (abs(results(2)) <= tiny(1.0_dp)) error stop 'Zero GMRES current response.'
        if (results(1) /= 2.5_dp-results(3)/results(2)) &
            error stop 'GMRES K coefficient disagrees with current responses.'
    end subroutine check_gmres

    subroutine check_batch(name, identities, first_expected, reference_a, &
                           reference_b, expected_convergence)
        character(*), intent(in) :: name
        integer, intent(in) :: identities(:), first_expected
        real(dp), intent(in), optional :: reference_a(3), reference_b(3)
        integer, intent(in), optional :: expected_convergence
        real(dp) :: results(3, size(identities)), s(size(identities))
        real(dp) :: r(size(identities)), z(size(identities)), nu(size(identities))
        real(dp), parameter :: expected_s(2) = [0.5_dp, 0.49_dp]
        real(dp), parameter :: expected_nu(2) = [0.01_dp, 0.02_dp]
        integer :: converged(size(identities)), first, last, row, expected_flag, solver

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
        if (present(reference_a)) then
            call h5_get(file_id, 'k_iter_solver', solver)
            if (solver /= 1) error stop 'Missing batch GMRES solver marker.'
        else if (h5_exists(file_id, 'k_iter_solver')) then
            error stop 'Default batch schema changed.'
        end if
        call h5_close(file_id)
        if (.not. k_only) error stop 'Missing batch output marker.'
        if (first /= first_expected .or. last /= first + size(identities) - 1) &
            error stop 'Incorrect batch range.'
        if (.not. all(ieee_is_finite(results))) error stop 'Nonfinite batch result.'
        if (any(s /= expected_s(identities)) .or. any(nu /= expected_nu(identities)) &
            .or. any(r /= 0.0_dp) .or. any(z /= 0.0_dp)) &
            error stop 'Incorrect batch surface identities.'
        expected_flag = 1
        if (present(expected_convergence)) expected_flag = expected_convergence
        if (any(converged /= expected_flag)) &
            error stop 'Incorrect batch convergence status.'
        if (expected_flag == 0) return
        do row = 1, size(identities)
            if (identities(row) == 1) then
                actual = expected
                if (present(reference_a)) actual = reference_a
            else
                actual = expected_b
                if (present(reference_b)) actual = reference_b
            end if
            if (any(results(:, row) /= actual)) then
                print *, name, ' row ', row, ': ', results(:, row), ' expected ', actual
                error stop 'Batch changed the current responses or k coefficient.'
            end if
        end do
    end subroutine check_batch
end program test_k_coefficient_only
