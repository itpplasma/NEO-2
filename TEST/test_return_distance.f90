program test_return_distance
    use mag_interface_mod, only: return_distance_improves
    implicit none

    integer :: failures
    real(kind(1.0d0)), parameter :: first_distance = &
        6.2831853072124133d-2
    real(kind(1.0d0)), parameter :: numerical_tie = &
        6.2831849552580366d-2
    real(kind(1.0d0)), parameter :: closed_return = &
        3.8516034805979871d-10

    failures = 0
    if (return_distance_improves(numerical_tie, first_distance)) then
        print *, 'FAIL: integration noise replaced a tied return distance'
        failures = failures + 1
    end if
    if (.not. return_distance_improves(closed_return, first_distance)) then
        print *, 'FAIL: closed field-line return was not accepted'
        failures = failures + 1
    end if
    if (return_distance_improves(first_distance, first_distance)) then
        print *, 'FAIL: equal return distance was accepted'
        failures = failures + 1
    end if

    if (failures == 0) then
        print *, 'All tests passed!'
    else
        error stop 'return-distance comparison contract failed'
    end if
end program test_return_distance
