program test_collop_bspline_precompute
    use, intrinsic :: ieee_arithmetic, only: ieee_get_flag, ieee_invalid, &
        ieee_set_flag, ieee_set_halting_mode
    use collisionality_mod, only: collop_bspline_dist, collop_bspline_order, &
        collop_bspline_taylor, phi_x_max
    use collop_compute, only: init_collop
    use rkstep_mod, only: lag, leg
    implicit none

    logical :: invalid_raised

    lag = 3
    leg = 3
    phi_x_max = 4.0d0
    collop_bspline_order = 2
    collop_bspline_dist = 1.0d0
    collop_bspline_taylor = .true.
    call ieee_set_flag(ieee_invalid, .false.)
    call ieee_set_halting_mode(ieee_invalid, .true.)
    call init_collop(11, 11, 0.0d0, 0.0d0)
    call ieee_get_flag(ieee_invalid, invalid_raised)

    if (invalid_raised) error stop 'FAIL: B-spline precompute raised IEEE invalid'
    print *, 'All tests passed!'
end program test_collop_bspline_precompute
