program test_ripple_solver_route
    use ripple_solver_route_mod, only: uses_single_propagator_route
    implicit none

    if (.not. uses_single_propagator_route(0, 1, 1)) &
        error stop 'FAIL: axisymmetric tokamak route was rejected'
    if (.not. uses_single_propagator_route(1, 0, 0)) &
        error stop 'FAIL: homogeneous-field route was rejected'
    if (uses_single_propagator_route(0, 0, 1)) &
        error stop 'FAIL: non-axisymmetric tokamak route was accepted'
    if (uses_single_propagator_route(1, 0, 1)) &
        error stop 'FAIL: general multi-propagator route was accepted'

    print *, 'All tests passed!'
end program test_ripple_solver_route
