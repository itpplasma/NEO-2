program test_collop_zero_endpoint
    use, intrinsic :: ieee_arithmetic, only: ieee_get_flag, ieee_invalid, &
        ieee_set_flag, ieee_set_halting_mode
    use collop_compute, only: init_collop
    use rkstep_mod, only: lag, leg
    implicit none

    logical :: invalid_raised

    ! Polynomial bases keep every integrand nonzero, so any IEEE invalid
    ! comes from the inner-kernel endpoint at x = 0.
    lag = 0
    leg = 1
    call ieee_set_halting_mode(ieee_invalid, .false.)
    call ieee_set_flag(ieee_invalid, .false.)
    call init_collop(0, 0, 0.0d0, 0.0d0)
    call ieee_get_flag(ieee_invalid, invalid_raised)

    if (invalid_raised) error stop 'FAIL: inner-kernel precomputation raised IEEE invalid'
    print *, 'All tests passed!'
end program test_collop_zero_endpoint
