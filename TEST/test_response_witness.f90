program test_response_witness
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use response_witness_mod
    use hdf5_tools, only: hid_t, h5_open, h5_get, h5_close, h5_stream_write
    implicit none
    integer(hid_t) :: file_id
    integer :: i, j, k, complete, unit
    character(len=160) :: first, second, disabled
    character(len=32) :: mode
    logical :: exists
    complex(dp) :: values(2, 3, 3), original(2, 3, 3)
    real(dp) :: saved(2, 3, 3), saved2(2, 3), saved1(2)

    call get_command_argument(1, mode)
    if (trim(mode) == 'overwrite') then
        first = 'response_witness_p17_a2_c1.h5'
        open(newunit=unit, file=trim(first), status='replace')
        close(unit)
        lsw_response_witness = .true.
        call response_witness_begin(17, 2, file_id, first)
        stop 0 ! Wrong acceptance must fail the WILL_FAIL test.
    end if
    if (lsw_response_witness) error stop 'FAIL witness enabled by default'
    call response_witness_begin(17, 2, file_id, disabled)
    if (file_id /= -1 .or. len_trim(disabled) /= 0) &
        error stop 'FAIL default path opened a file'

    do k = 1, 3
        do j = 1, 3
            do i = 1, 2
                ! Distinct cross-species/drive coefficients, including phase.
                values(i,j,k) = cmplx(i+10*j+100*k, 2*i-3*j+7*k, dp)
            end do
        end do
    end do
    original = values
    lsw_response_witness = .true.
    call response_witness_begin(17, 2, file_id, first)
    if (h5_stream_write) error stop 'FAIL changed global writer default'
    call response_witness_complex(file_id, 'y0', values)
    call response_witness_complex(file_id, 'rhs', values(:,:,2))
    call response_witness_complex(file_id, 'integral', values(:,3,1))
    call response_witness_end(file_id)
    if (file_id /= -1) error stop 'FAIL witness left open'
    if (any(values /= original)) error stop 'FAIL observer changed coefficients'
    call h5_open(trim(first), file_id)
    call h5_get(file_id, 'y0_real', saved)
    if (any(saved /= real(original, dp))) error stop 'FAIL real reconstruction'
    call h5_get(file_id, 'y0_imag', saved)
    if (any(saved /= aimag(original))) error stop 'FAIL phase reconstruction'
    call h5_get(file_id, 'rhs_imag', saved2)
    if (any(saved2 /= aimag(original(:,:,2)))) &
        error stop 'FAIL matrix phase reconstruction'
    call h5_get(file_id, 'integral_real', saved1)
    if (any(saved1 /= real(original(:,3,1), dp))) &
        error stop 'FAIL vector reconstruction'
    call h5_get(file_id, 'complete', complete)
    if (complete /= 1) error stop 'FAIL premature completion'
    call h5_close(file_id)

    h5_stream_write = .true.
    call response_witness_begin(117, 12, file_id, second)
    if (.not. h5_stream_write) error stop 'FAIL changed enabled writer default'
    h5_stream_write = .false.
    if (first == second) error stop 'FAIL witness overwrote another propagator'
    call response_witness_end(file_id)
    lsw_response_witness = .false.
    call response_witness_begin(999, 999, file_id, disabled)
    call response_witness_complex(file_id, 'ignored', values)
    call response_witness_end(file_id)
    inquire(file='response_witness_p999_a999_c3.h5', exist=exists)
    if (exists) error stop 'FAIL disabled observer wrote a file'
    if (any(values /= original)) error stop 'FAIL disabled observer changed input'
    open(newunit=unit, file=trim(first), status='old')
    close(unit, status='delete')
    open(newunit=unit, file=trim(second), status='old')
    close(unit, status='delete')
    print *, 'All tests passed!'
end program test_response_witness
