program test_k_sparse_backend
    ! Known nonsymmetric solution; the second solve reuses the factor.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use mpi, only: MPI_Init, MPI_Finalize
    use k_sparse_backend_mod, only: k_sparse_backend_t
    implicit none

    integer, parameter :: pcol(4) = [1, 3, 6, 8]
    integer, parameter :: row(7) = [1, 2, 1, 2, 3, 2, 3]
    real(dp), parameter :: val(7) = [4._dp, -2._dp, 1._dp, 3._dp, &
                                     1._dp, 2._dp, 5._dp]
    real(dp), parameter :: expected(3) = [1._dp, 2._dp, 3._dp]
    real(dp), parameter :: forcing(3) = [6._dp, 10._dp, 17._dp]
    type(k_sparse_backend_t) :: solver
    real(dp) :: rhs(3)
    integer :: multiple, ierr

    call MPI_Init(ierr)
    call solver%factor(3, 3, 7, pcol, row, val)
    do multiple = 1, 2
        rhs = multiple*forcing
        call solver%solve(rhs)
        if (.not. maxval(abs(rhs - multiple*expected)) <= 64*epsilon(1._dp)*3*multiple) &
            error stop 'Sparse solution differs from the known solution.'
    end do
    call solver%free()
    call MPI_Finalize(ierr)
    print *, 'All tests passed!'
end program test_k_sparse_backend
