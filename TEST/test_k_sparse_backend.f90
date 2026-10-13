program test_k_sparse_backend
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use mpi, only: MPI_Init, MPI_Finalize, MPI_SUCCESS
    use k_sparse_backend_mod, only: k_sparse_backend_t, require_k_sparse_backend
    implicit none

    integer, parameter :: pcol(4) = [1, 3, 6, 8]
    integer, parameter :: row(7) = [1, 2, 1, 2, 3, 2, 3]
    real(dp), parameter :: val(7) = [4._dp, -2._dp, 1._dp, 3._dp, &
                                    1._dp, 2._dp, 5._dp]
    real(dp), parameter :: matrix(3, 3) = reshape([4._dp, -2._dp, 0._dp, &
        1._dp, 3._dp, 1._dp, 0._dp, 2._dp, 5._dp], [3, 3])
    real(dp), parameter :: expected(3) = [1._dp, 2._dp, 3._dp]
    real(dp), parameter :: forcing(3) = [6._dp, 10._dp, 17._dp]
    real(dp), parameter :: tolerance = 64._dp*epsilon(1._dp)
    type(k_sparse_backend_t) :: solver
    real(dp) :: rhs(3), residual(3), scale
    integer :: backend, multiple, ierr, ios
    character(16) :: argument

    call get_command_argument(1, argument)
    read(argument, *, iostat=ios) backend
    if (ios /= 0) error stop 'Expected an integer backend argument.'
    call MPI_Init(ierr)
    if (ierr /= MPI_SUCCESS) error stop 'MPI initialization failed.'
    call require_k_sparse_backend(backend)
    call solver%factor(backend, 3, 3, 7, pcol, row, val)
    do multiple = 1, 2
        rhs = real(multiple, dp)*forcing
        call solver%solve(rhs)
        if (.not. all(ieee_is_finite(rhs))) error stop 'Nonfinite solution.'
        if (maxval(abs(rhs-real(multiple, dp)*expected)) > &
            tolerance*maxval(abs(real(multiple, dp)*expected))) &
            error stop 'Sparse solution differs from the known solution.'
        residual = matmul(matrix, rhs)-real(multiple, dp)*forcing
        scale = maxval(sum(abs(matrix), dim=2))*maxval(abs(rhs)) + &
            real(multiple, dp)*maxval(abs(forcing))
        if (maxval(abs(residual)) > tolerance*scale) &
            error stop 'Sparse solution has excessive backward error.'
    end do
    call solver%free()
    call MPI_Finalize(ierr)
    if (ierr /= MPI_SUCCESS) error stop 'MPI finalization failed.'
    print *, 'All tests passed!'
end program test_k_sparse_backend
