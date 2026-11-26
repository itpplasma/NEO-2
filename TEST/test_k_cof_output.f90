!> k_cof written to the multispecies output must be 5/2 - D32/D31 for every
!> species pair, and 0 where D31 vanishes (axisymmetric solution not computed).
program test_k_cof_output
    use nrtype, only: dp
    use hdf5_tools, only: HID_T, h5_init, h5_deinit, h5_create, h5_open, h5_close, h5_get
    use ntv_mod, only: write_k_cof
    implicit none

    character(len=*), parameter :: fname = 'k_cof_test.h5'
    real(dp), parameter :: d31(0:3) = [-0.8_dp, 0.1_dp, 0.0_dp, 1.25_dp]
    real(dp), parameter :: d32(0:3) = [0.4_dp, -0.3_dp, 0.7_dp, 2.0_dp]
    real(dp), parameter :: expected(0:3) = [3.0_dp, 5.5_dp, 0.0_dp, 0.9_dp]
    real(dp) :: k_read(0:3)
    integer(HID_T) :: h5id

    call h5_init()
    call h5_create(fname, h5id)
    call write_k_cof(h5id, d31, d32)
    call h5_close(h5id)

    k_read = -1.0_dp
    call h5_open(fname, h5id)
    call h5_get(h5id, 'k_cof', k_read)
    call h5_close(h5id)
    call h5_deinit()

    if (any(abs(k_read - expected) > 1.0e-14_dp)) then
        print *, "FAIL k_cof read back:", k_read, " expected:", expected
        error stop 1
    end if
    print *, "All tests passed!"
end program test_k_cof_output
