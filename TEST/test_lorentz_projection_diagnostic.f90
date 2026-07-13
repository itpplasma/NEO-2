program test_lorentz_projection_diagnostic
    use, intrinsic :: iso_fortran_env, only: real64
    use lorentz_projection_diagnostics_mod, only: &
        assemble_local_constant_state, compute_local_projection_residuals, &
        compute_sparse_constant_residual, local_projection_residuals, &
        record_local_constant_row, record_local_constant_stage_residuals, &
        record_local_projection_residuals
    implicit none

    real(real64) :: source_p(1, 3), source_m(1, 3)
    real(real64) :: flux_p(3, 1), flux_m(3, 1)
    real(real64) :: amat_p_p(1, 1), amat_m_p(1, 1)
    real(real64) :: amat_p_m(1, 1), amat_m_m(1, 1)
    real(real64) :: eta_l(0:0), eta_r(0:0)
    real(real64) :: constant_state(8), constant_rhs(8), matrix_values(8)
    real(real64) :: residual, scale
    real(real64) :: eta(0:1), bhat(1:2)
    integer :: npl(1:2), ind_start(1:2), matrix_index(8), residual_index
    type(local_projection_residuals) :: residuals
    character(len=1024) :: filename, line, row_filename
    integer :: ierr, iunit, line_count, status

    source_p = reshape([1.0_real64, 2.0_real64, 3.0_real64], [1, 3])
    source_m = reshape([-0.5_real64, -1.0_real64, -1.5_real64], [1, 3])
    flux_p = reshape([2.0_real64, 3.0_real64, 4.0_real64], [3, 1])
    flux_m = reshape([-1.0_real64, -2.0_real64, -3.0_real64], [3, 1])
    amat_p_p = 0.0_real64
    amat_m_m = 0.0_real64
    amat_m_p = 1.0_real64
    amat_p_m = 1.0_real64
    eta_l = 0.0_real64
    eta_r = 0.0_real64

    call compute_local_projection_residuals(source_p, source_m, flux_p, flux_m, &
        amat_p_p, amat_m_p, amat_p_m, amat_m_m, eta_l, eta_r, 0.4_real64, &
        0.6_real64, residuals, ierr)
    if (ierr /= 0) error stop 'FAIL: valid projection residual input rejected'
    if (maxval(abs(residuals%source_defect - [0.5_real64, 1.0_real64, &
        1.5_real64])) > 1.0e-15_real64) &
        error stop 'FAIL: source defect is incorrect'
    if (maxval(abs(residuals%flux_null - [0.2_real64, 0.0_real64, &
        -0.2_real64])) > 1.0e-15_real64) &
        error stop 'FAIL: flux null contraction is incorrect'
    if (residuals%measure /= 1.0_real64) &
        error stop 'FAIL: boundary measure is incorrect'
    if (residuals%left_transport /= 0.0_real64) &
        error stop 'FAIL: conservative left transport was not recognized'
    if (residuals%right_transport /= 0.0_real64) &
        error stop 'FAIL: null-vector transport was not recognized'
    if (residuals%intertwining /= 0.0_real64) &
        error stop 'FAIL: projector intertwining was not recognized'

    npl = [1, 1]
    ind_start = [0, 4]
    eta = [0.0_real64, 0.2_real64]
    bhat = [2.0_real64, 2.5_real64]
    call assemble_local_constant_state(constant_state, constant_rhs, npl, &
        ind_start, eta, bhat, 1, 2, 0, ierr)
    if (ierr /= 0) error stop 'FAIL: valid constant state input rejected'
    if (maxval(abs(constant_state - [0.2_real64, 0.3_real64, 0.3_real64, &
        0.2_real64, 0.2_real64, 0.2_real64, 0.2_real64, 0.2_real64])) &
        > 1.0e-15_real64) error stop 'FAIL: constant state is incorrect'
    if (maxval(abs(constant_rhs - [0.2_real64, 0.3_real64, 0.0_real64, &
        0.0_real64, 0.0_real64, 0.0_real64, 0.2_real64, 0.2_real64])) &
        > 1.0e-15_real64) error stop 'FAIL: constant boundary rhs is incorrect'
    matrix_index = [(status, status=1, 8)]
    matrix_values = 1.0_real64
    call compute_sparse_constant_residual(matrix_index, matrix_index, &
        matrix_values, constant_state, constant_state, residual, scale, &
        residual_index, ierr)
    if (ierr /= 0 .or. residual /= 0.0_real64) &
        error stop 'FAIL: exact sparse constant state was not recognized'

    call record_local_projection_residuals(7, source_p, source_m, flux_p, &
        flux_m, amat_p_p, amat_m_p, amat_p_m, amat_m_m, eta_l, eta_r, &
        0.4_real64, 0.6_real64, ierr)
    if (ierr /= 0) error stop 'FAIL: projection trace was not written'
    call record_local_constant_stage_residuals(7, residual, scale, &
        residual_index, residual, scale, residual_index, 2, 0, -1, 1, ierr)
    if (ierr /= 0) error stop 'FAIL: constant-stage trace was not written'
    call record_local_constant_row(7, 1, matrix_index, matrix_index, &
        matrix_values, constant_state, constant_state, ierr)
    if (ierr /= 0) error stop 'FAIL: constant-row trace was not written'
    call get_environment_variable('NEO2_LOCAL_PROJECTION_TRACE_FILE', &
        value=filename, status=status)
    if (status /= 0) error stop 'FAIL: projection trace path is absent'
    open (newunit=iunit, file=trim(filename), status='old', action='read')
    line_count = 0
    do
        read (iunit, '(a)', iostat=status) line
        if (status /= 0) exit
        line_count = line_count + 1
    end do
    close (iunit)
    if (line_count /= 17) error stop 'FAIL: projection trace row count differs'
    row_filename = trim(filename)//'.rows'
    open (newunit=iunit, file=trim(row_filename), status='old', action='read')
    line_count = 0
    do
        read (iunit, '(a)', iostat=status) line
        if (status /= 0) exit
        line_count = line_count + 1
    end do
    close (iunit)
    if (line_count /= 2) error stop 'FAIL: constant-row trace row count differs'

    call compute_local_projection_residuals(source_p, source_m, flux_p, flux_m, &
        amat_p_p, amat_m_p, amat_p_m, amat_m_m, eta_l, eta_r, -0.4_real64, &
        0.6_real64, residuals, ierr)
    if (ierr == 0) error stop 'FAIL: nonpositive boundary weight accepted'

    print *, 'All tests passed!'
end program test_lorentz_projection_diagnostic
