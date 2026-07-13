program test_phi_placement_diagnostic
    use, intrinsic :: iso_fortran_env, only: real64
    use flint_mod, only: record_phi_placement

    implicit none

    character(len=512) :: filename, header
    character(len=16) :: outcome
    integer :: decision_count, eta_index, interval, ierr, iunit, part
    integer :: result_count, sequence, status, tag
    real(real64) :: hphi, phi_end, phi_start

    call record_phi_placement(12, 'placed', 3, 2, 7, 1.0_real64, 1.5_real64, &
        0.1_real64, 5, 5, ierr)
    if (ierr /= 0) error stop 'FAIL: phi placement record returned an error'
    call get_environment_variable('NEO2_PHI_PLACEMENT_FILE', filename)
    open(newunit=iunit, file=trim(filename), status='old', action='read', &
        iostat=status)
    if (status /= 0) error stop 'FAIL: phi placement output is absent'
    read(iunit, '(a)', iostat=status) header
    if (status /= 0 .or. trim(header) /= &
        'sequence,tag,outcome,interval,part,eta_index,phi_start,phi_end,hphi,'// &
        'decision_count,result_count') &
        error stop 'FAIL: phi placement header is wrong'
    read(iunit, *, iostat=status) sequence, tag, outcome, interval, part, &
        eta_index, phi_start, phi_end, hphi, decision_count, result_count
    close(iunit)
    if (status /= 0 .or. sequence /= 1 .or. tag /= 12 &
        .or. trim(outcome) /= 'placed' .or. interval /= 3 &
        .or. part /= 2 .or. eta_index /= 7 .or. decision_count /= 5 &
        .or. result_count /= 5 .or. abs(phi_start - 1.0_real64) &
        > 1.0e-15_real64 .or. abs(phi_end - 1.5_real64) &
        > 1.0e-15_real64 .or. abs(hphi - 0.1_real64) &
        > 1.0e-15_real64) &
        error stop 'FAIL: phi placement record is wrong'

    print '(a)', 'All tests passed!'
end program test_phi_placement_diagnostic
