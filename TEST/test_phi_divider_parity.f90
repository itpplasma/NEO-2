program test_phi_divider_parity
    use device_mod, only: fieldpropagator
    use flint_mod, only: phi_divide, phi_divider, phiarr
    implicit none

    integer, parameter :: ninterval = 4
    integer :: forward_map(0:ninterval), reverse_map(0:ninterval)
    real(kind(1.0d0)) :: forward_nodes(0:ninterval), reverse_nodes(0:ninterval)
    integer :: i, failures

    failures = 0
    call refined_node_map([2, 1, 2, 2], forward_map, forward_nodes)
    call refined_node_map([2, 2, 1, 2], reverse_map, reverse_nodes)

    do i = 0, ninterval
        if (forward_nodes(i) /= real(i, kind(1.0d0)) .or. &
            reverse_nodes(i) /= real(i, kind(1.0d0))) then
            print *, 'FAIL: refinement dropped an original grid point'
            failures = failures + 1
        end if
        if (reverse_map(i) /= forward_map(ninterval) - forward_map(ninterval - i)) then
            print *, 'FAIL: reversed refinement changed the mapped grid'
            failures = failures + 1
        end if
    end do

    if (any(forward_map(1:) <= forward_map(:ninterval - 1)) .or. &
        any(reverse_map(1:) <= reverse_map(:ninterval - 1))) then
        print *, 'FAIL: refinement did not preserve grid ordering'
        failures = failures + 1
    end if

    if (mod(forward_map(ninterval), 2) /= 0) then
        print *, 'FAIL: refinement made the final Runge-Kutta point a half step'
        failures = failures + 1
    end if

    if (failures == 0) then
        print *, 'All tests passed!'
    else
        error stop 'phi-divider parity contract failed'
    end if

contains

    subroutine refined_node_map(divisions, node_map, node_values)
        integer, intent(in) :: divisions(ninterval)
        integer, intent(out) :: node_map(0:ninterval)
        real(kind(1.0d0)), intent(out) :: node_values(0:ninterval)

        integer :: i, phi_eta_ind(0:ninterval, 2)

        allocate(fieldpropagator)
        allocate(fieldpropagator%coords)
        allocate(fieldpropagator%coords%x2(0:ninterval))
        fieldpropagator%coords%x2 = [(real(i, kind(1.0d0)), i = 0, ninterval)]

        allocate(phi_divide(ninterval))
        phi_divide = divisions
        do i = 0, ninterval
            phi_eta_ind(i, :) = i
        end do

        call phi_divider(ninterval, phi_eta_ind)
        node_map = phi_eta_ind(:, 1)
        node_values = phiarr(node_map)

        deallocate(phiarr, phi_divide)
        deallocate(fieldpropagator%coords%x2)
        deallocate(fieldpropagator%coords)
        deallocate(fieldpropagator)
    end subroutine refined_node_map
end program test_phi_divider_parity
