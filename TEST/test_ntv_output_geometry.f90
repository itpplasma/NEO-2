program test_ntv_output_geometry
    use ntv_mod, only: y_ntv_mod, resolve_ntv_output_geometry
    implicit none

    real(kind(1.0d0)), allocatable :: output_y(:)

    if (allocated(y_ntv_mod)) deallocate(y_ntv_mod)
    call resolve_ntv_output_geometry(output_y, [1.0d0, 2.0d0, 3.0d0])
    if (any(output_y /= [1.0d0, 2.0d0, 3.0d0])) &
        error stop 'FAIL: final geometry was not selected'

    allocate(y_ntv_mod(2))
    y_ntv_mod = [4.0d0, 5.0d0]
    call resolve_ntv_output_geometry(output_y)
    if (any(output_y /= y_ntv_mod)) &
        error stop 'FAIL: legacy exchange geometry was not selected'

    print *, 'All tests passed!'
end program test_ntv_output_geometry
