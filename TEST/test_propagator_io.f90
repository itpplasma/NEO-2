program test_propagator_io
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use propagator_mod, only: propagator, read_propagator_content, &
        write_propagator_content
    implicit none

    type(propagator), pointer :: loaded, written
    integer :: i, unit
    real(dp) :: expected(3, 3)

    allocate (loaded, written)
    written%nr_joined = 0
    written%fieldpropagator_tag_s = 41
    written%fieldpropagator_tag_e = 41
    written%fieldperiod_tag_s = 7
    written%fieldperiod_tag_e = 7
    written%phi_l = 0.25_dp
    written%phi_r = 0.5_dp
    written%p%npart = 1
    written%p%npass_l = 1
    written%p%npass_r = 1
    written%p%nvelocity = 0
    written%p%eta_boundary_l = 1.0_dp
    written%p%eta_boundary_r = 1.0_dp
    expected = reshape([(real(i, dp), i = 1, 9)], shape(expected))
    allocate (written%p%qflux(3, 3))
    written%p%qflux = expected

    call write_propagator_content(written, 3)
    call read_propagator_content(loaded, 3, 41, 41, 1)

    if (.not. allocated(loaded%p%qflux)) &
        error stop 'FAIL: qflux was not allocated while reading propagator content'
    if (any(shape(loaded%p%qflux) /= [3, 3])) &
        error stop 'FAIL: reconstructed qflux has the wrong shape'
    if (any(loaded%p%qflux /= expected)) &
        error stop 'FAIL: reconstructed qflux values changed during serialization'

    open (newunit=unit, file='propagator_41_41.prop', status='old')
    close (unit, status='delete')
    deallocate (loaded, written)
    print *, 'All tests passed!'
end program test_propagator_io
