program test_qflux_interface_diagnostic
    use, intrinsic :: iso_fortran_env, only: real64
    use qflux_profile_mod, only: record_qflux_interface_traces
    implicit none

    character(len=512) :: filename, header
    character(len=5) :: side
    character(len=1) :: direction
    integer :: band, force, i, ierr, iunit, laguerre, sequence, status, tag
    integer :: block_base(2), block_npassing(2), rows
    real(real64) :: eta, eta_grid(0:1), flux_kernel, flux_row(8), phi
    real(real64) :: phi_grid(2), solution(8, 3), source_rhs, rhs(8, 3)
    real(real64) :: step_minus(2), step_plus(2), value
    logical :: found_counter_record

    block_base = [0, 4]
    block_npassing = 1
    eta_grid = [0.25_real64, 0.75_real64]
    phi_grid = [1.0_real64, 2.0_real64]
    step_plus = [2.0_real64, 4.0_real64]
    step_minus = [5.0_real64, 10.0_real64]
    flux_row = [(real(i, real64), i = 1, size(flux_row))]
    rhs = reshape([(100.0_real64 + real(i, real64), i = 1, size(rhs))], &
        shape(rhs))
    solution = rhs + 1000.0_real64

    call record_qflux_interface_traces(12, phi_grid, eta_grid, block_base, &
        block_npassing, 0, step_plus, step_minus, flux_row, rhs, solution, ierr)
    if (ierr /= 0) error stop 'FAIL: interface trace record returned an error'
    call get_environment_variable('NEO2_INTERFACE_TRACE_FILE', filename)
    open(newunit=iunit, file=trim(filename), status='old', action='read', &
        iostat=status)
    if (status /= 0) error stop 'FAIL: interface trace output is absent'
    read(iunit, '(a)', iostat=status) header
    if (status /= 0 .or. trim(header) /= &
        'sequence,tag,side,direction,laguerre,force,eta_index,eta,phi,' // &
        'flux_kernel,source_rhs,source_solution') &
        error stop 'FAIL: interface trace header is wrong'

    rows = 0
    found_counter_record = .false.
    do
        read(iunit, *, iostat=status) sequence, tag, side, direction, laguerre, &
            force, band, eta, phi, flux_kernel, source_rhs, value
        if (status < 0) exit
        if (status > 0) error stop 'FAIL: interface trace row cannot be read'
        rows = rows + 1
        if (sequence /= rows .or. tag /= 12 .or. laguerre /= 0) &
            error stop 'FAIL: interface trace identity is wrong'
        if (trim(side) == 'right' .and. direction == 'm' .and. band == 1 &
            .and. force == 3) then
            found_counter_record = .true.
            if (abs(eta - 0.75_real64) > 1.0e-15_real64 &
                .or. abs(phi - 2.0_real64) > 1.0e-15_real64 &
                .or. abs(flux_kernel - 0.7_real64) > 1.0e-15_real64 &
                .or. abs(source_rhs - rhs(7, 3)) > 1.0e-15_real64 &
                .or. abs(value - solution(7, 3)) > 1.0e-15_real64) &
                error stop 'FAIL: counter-passing trace mapping is wrong'
        end if
    end do
    close(iunit)
    if (rows /= 24) error stop 'FAIL: interface trace row count is wrong'
    if (.not. found_counter_record) &
        error stop 'FAIL: expected counter-passing trace record is absent'

    print '(a)', 'All tests passed!'
end program test_qflux_interface_diagnostic
