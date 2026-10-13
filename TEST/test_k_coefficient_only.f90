program test_k_coefficient_only
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use hdf5_tools, only: HID_T, h5_init, h5_deinit, h5_open, h5_get, h5_exists, &
        h5_close, h5_open_group, h5_close_group
    implicit none

    integer(HID_T) :: file_id, group_id
    character(len=1024) :: work_dir
    character(len=16) :: only_case
    real(dp) :: expected(3), actual(3), qflux(3, 3)
    real(dp) :: residual(2)
    integer :: solver, applications(2), budget
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

    call h5_open(trim(work_dir)//'/gmres/fulltransp.h5', file_id)
    call h5_get(file_id, 'k_cof', actual(1))
    call h5_get(file_id, 'qcurr_g', actual(2))
    call h5_get(file_id, 'qcurr_t', actual(3))
    call h5_get(file_id, 'k_iter_solver', solver)
    call h5_get(file_id, 'k_source_applications', applications)
    call h5_get(file_id, 'k_source_residual', residual)
    call h5_close(file_id)
    call h5_open(trim(work_dir)//'/gmres/neo2_config.h5', file_id)
    call h5_open_group(file_id, 'settings', group_id)
    call h5_get(group_id, 'niter', budget)
    call h5_close_group(group_id)
    call h5_close(file_id)
    call h5_deinit()
    if (.not. all(ieee_is_finite(actual)) .or. &
        .not. all(ieee_is_finite(residual))) error stop 'Nonfinite GMRES result.'
    if (solver /= 1) error stop 'Missing GMRES solver flag.'
    if (any(residual > 5.d-2) .or. any(residual < 0.d0)) &
        error stop 'GMRES residual exceeds the requested tolerance.'
    if (maxval(applications) /= budget) &
        error stop 'GMRES did not converge on the final allowed application.'
    if (abs(actual(2)) <= tiny(1.0_dp)) error stop 'Zero GMRES current response.'
    if (actual(1) /= 2.5_dp-actual(3)/actual(2)) &
        error stop 'GMRES K coefficient disagrees with current responses.'
    print *, 'All tests passed!'
end program test_k_coefficient_only
