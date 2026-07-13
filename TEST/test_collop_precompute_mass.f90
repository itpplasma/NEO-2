program test_collop_precompute_mass
    use, intrinsic :: ieee_arithmetic, only: ieee_get_flag, ieee_invalid, &
        ieee_set_flag, ieee_set_halting_mode
    use collop_compute, only: scaled_num_sub_intervals
    implicit none

    integer :: n_sub
    logical :: invalid_raised

    call ieee_set_flag(ieee_invalid, .false.)
    call ieee_set_halting_mode(ieee_invalid, .true.)
    n_sub = scaled_num_sub_intervals()
    call ieee_get_flag(ieee_invalid, invalid_raised)

    if (invalid_raised) error stop 'FAIL: default subdivision used an undefined mass ratio'
    if (n_sub /= 5) error stop 'FAIL: default subdivision is not the equal-mass value'
    print *, 'All tests passed!'
end program test_collop_precompute_mass
