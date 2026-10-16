module neo2_sparse_config
    ! Per-process application policy. FortSparse owns solver payload accounting.
    use, intrinsic :: iso_fortran_env, only: int64, error_unit
    use mpi, only: MPI_INITIALIZED, MPI_FINALIZED, MPI_COMM_RANK, &
                   MPI_COMM_WORLD, MPI_ABORT, MPI_SUCCESS
    implicit none
    private
    public :: sparse_host_memory_limit_bytes, sparse_low_memory
    public :: sparse_validate_memory_options, sparse_abort
    integer(int64) :: sparse_host_memory_limit_bytes = 0_int64
    logical :: sparse_low_memory = .false.
contains
    subroutine sparse_validate_memory_options()
        if (sparse_host_memory_limit_bytes < 0_int64) then
            call sparse_abort('sparse_host_memory_limit_bytes must be nonnegative')
        end if
#ifndef NEO2_NATIVE_SPARSE
        if (sparse_host_memory_limit_bytes /= 0_int64 .or. sparse_low_memory) then
            call sparse_abort('memory controls require NEO2_NATIVE_SPARSE=ON')
        end if
#endif
    end subroutine sparse_validate_memory_options

    subroutine sparse_abort(message, code)
        ! Solves need not occur on all ranks: never insert a collective here.
        ! MPI_ABORT terminates peers blocked in the existing scheduler/collectives.
        character(*), intent(in) :: message
        integer, optional, intent(in) :: code
        integer :: exit_code, ierr, rank
        logical :: initialized, finalized
        exit_code = 1
        if (present(code)) exit_code = code
        if (len_trim(message) > 0) &
            write (error_unit, '(a)') 'NEO-2 sparse: '//trim(message)
        initialized = .false.
        finalized = .true.
        call MPI_INITIALIZED(initialized, ierr)
        if (ierr == MPI_SUCCESS .and. initialized) then
            call MPI_FINALIZED(finalized, ierr)
            if (ierr == MPI_SUCCESS .and. .not. finalized) then
                call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
                if (ierr == MPI_SUCCESS) &
                    write (error_unit, '(a,i0)') 'NEO-2 sparse failure rank=', rank
                flush (error_unit)
                call MPI_ABORT(MPI_COMM_WORLD, exit_code, ierr)
            end if
        end if
        flush (error_unit)
        error stop exit_code
    end subroutine sparse_abort
end module neo2_sparse_config
