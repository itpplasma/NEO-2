! Opt-in, raw native Arnoldi coefficients. No equation or normalization change.
module response_witness_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use hdf5_tools, only: hid_t, h5_create, h5_add, h5_close, h5_stream_write
    implicit none
    private
    public :: lsw_response_witness, response_witness_begin
    public :: response_witness_complex, response_witness_end

    logical :: lsw_response_witness = .false.
    integer :: invocation = 0
    interface response_witness_complex
        module procedure witness_complex_1, witness_complex_2, witness_complex_3
    end interface
contains
    subroutine response_witness_begin(tag, species, file_id, filename)
        integer, intent(in) :: tag, species
        integer(hid_t), intent(out) :: file_id
        character(len=*), intent(out) :: filename
        logical :: exists, previous_stream_write

        file_id = -1
        filename = ''
        if (.not. lsw_response_witness) return
        invocation = invocation + 1
        write(filename, '(a,i0,a,i0,a,i0,a)') 'response_witness_p', tag, &
            '_a', species, '_c', invocation, '.h5'
        inquire(file=trim(filename), exist=exists)
        if (exists) error stop 'Response witness would overwrite an existing file'
        ! Streaming is a property of this new handle. Restore the global
        ! creation default immediately, before any native solver resumes.
        previous_stream_write = h5_stream_write
        h5_stream_write = .true.
        call h5_create(trim(filename), file_id)
        h5_stream_write = previous_stream_write
        call h5_add(file_id, 'schema_version', 1)
        call h5_add(file_id, 'propagator_tag', tag)
        call h5_add(file_id, 'receiving_species', species)
        call h5_add(file_id, 'invocation', invocation)
        call h5_add(file_id, 'physical_collision_action_qualified', 0)
        call h5_add(file_id, 'residual_response_qualified', 0)
    end subroutine response_witness_begin

    subroutine witness_complex_1(file_id, name, value)
        integer(hid_t), intent(in) :: file_id
        character(len=*), intent(in) :: name
        complex(dp), intent(in) :: value(:)
        if (file_id < 0) return
        call h5_add(file_id, name//'_real', real(value, dp), &
            lbound(value), ubound(value))
        call h5_add(file_id, name//'_imag', aimag(value), &
            lbound(value), ubound(value))
    end subroutine witness_complex_1

    subroutine witness_complex_2(file_id, name, value)
        integer(hid_t), intent(in) :: file_id
        character(len=*), intent(in) :: name
        complex(dp), intent(in) :: value(:, :)
        if (file_id < 0) return
        call h5_add(file_id, name//'_real', real(value, dp), &
            lbound(value), ubound(value))
        call h5_add(file_id, name//'_imag', aimag(value), &
            lbound(value), ubound(value))
    end subroutine witness_complex_2

    subroutine witness_complex_3(file_id, name, value)
        integer(hid_t), intent(in) :: file_id
        character(len=*), intent(in) :: name
        complex(dp), intent(in) :: value(:, :, :)
        if (file_id < 0) return
        call h5_add(file_id, name//'_real', real(value, dp), &
            lbound(value), ubound(value))
        call h5_add(file_id, name//'_imag', aimag(value), &
            lbound(value), ubound(value))
    end subroutine witness_complex_3

    subroutine response_witness_end(file_id)
        integer(hid_t), intent(inout) :: file_id
        if (file_id < 0) return
        call h5_add(file_id, 'complete', 1)
        call h5_close(file_id)
        file_id = -1
    end subroutine response_witness_end
end module response_witness_mod
