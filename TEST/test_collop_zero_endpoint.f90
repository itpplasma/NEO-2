program test_collop_zero_endpoint
    use, intrinsic :: ieee_arithmetic, only: ieee_get_flag, ieee_invalid, &
        ieee_set_flag, ieee_set_halting_mode
    use collop_compute, only: compute_energyscattering, gamma_ab, &
        init_collop, isw_relativistic, make_ortho, T_a, T_b
    use rkstep_mod, only: lag, leg
    implicit none

    logical :: invalid_raised
    real(kind(1.0d0)) :: denmm(1, 1)

    lag = 0
    leg = 1
    call ieee_set_flag(ieee_invalid, .false.)
    call ieee_set_halting_mode(ieee_invalid, .true.)
    call init_collop(0, 0, 0.0d0, 0.0d0)
    gamma_ab = 1.0d0
    isw_relativistic = 0
    make_ortho = .false.
    T_a = 1.0d0
    T_b = 1.0d0
    call compute_energyscattering(denmm)
    call ieee_get_flag(ieee_invalid, invalid_raised)

    if (invalid_raised) error stop 'FAIL: zero endpoint raised IEEE invalid'
    if (denmm(1, 1) /= 0.0d0) error stop 'FAIL: constant basis energy scattering is not zero'

    gamma_ab = sqrt(2.0d0)
    T_a = 2.0d0
    call compute_energyscattering(denmm)
    call ieee_get_flag(ieee_invalid, invalid_raised)

    if (invalid_raised) error stop 'FAIL: unequal-temperature energy scattering raised IEEE invalid'
    if (denmm(1, 1) == 0.0d0) error stop 'FAIL: unequal-temperature energy scattering was skipped'
    print *, 'All tests passed!'
end program test_collop_zero_endpoint
