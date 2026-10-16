program test_native_sparse_memory
    use, intrinsic :: iso_fortran_env, only: int64, output_unit
    use mpi, only: MPI_INIT, MPI_COMM_RANK, MPI_COMM_WORLD, MPI_BARRIER
    use sparse_mod, only: sparse_solve, sparse_solve_method, sparse_talk, &
                          sparse_host_memory_limit_bytes, sparse_low_memory
    use fortsparse, only: FORTSPARSE_RESOURCE_LIMIT
    implicit none
    integer, parameter :: dp = kind(1d0), n = 4
    real(dp) :: a(n, n), exact(n), rhs(n), b(n), block(n, 2)
    complex(dp) :: z(n, n), truth(n), zr(n), zb(n), zblock(n, 2)
    integer :: method, j, low, ierr, rank
    character(32) :: argument
    call get_command_argument(1, argument)
    sparse_talk = .false.
    a = reshape([0d0, 3d0, 0d0, 1d0, 2d0, 4d0, 1d0, 0d0, &
                 0d0, 1d0, 5d0, 2d0, 1d0, 0d0, 2d0, 6d0], [n, n])
    exact = [0.5d0, -1d0, 2d0, 0.25d0]
    rhs = matmul(a, exact)
    z = cmplx(a, 0d0, dp)
    z(1, 2) = z(1, 2) + (0d0, 0.5d0)
    z(3, 1) = (0d0, -0.25d0)
    truth = cmplx(exact, exact(n:1:-1), dp)
    zr = matmul(z, truth)
    sparse_solve_method = 3
    if (trim(argument) == '--mpi-cap') then
        call MPI_INIT(ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
        if (rank == 0) then
            ! A different rank fails during an asynchronous sparse operation.
            call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            error stop 'peer failure did not terminate the MPI job'
        end if
    end if
    if (len_trim(argument) > 0) then
        print '(a,i0)', 'Expected resource code=', FORTSPARSE_RESOURCE_LIMIT
        flush (output_unit)
        b = rhs
        if (trim(argument) == '--retained-cap') call sparse_solve(a, b, 1)
        sparse_host_memory_limit_bytes = 1_int64
        if (trim(argument) == '--cap-complex') then
            zb = zr
            call sparse_solve(z, zb)
        else if (trim(argument) == '--retained-cap') then
            call sparse_solve(a, b, 2)
        else
            call sparse_solve(a, b)
        end if
        error stop 'tiny payload cap unexpectedly accepted'
    end if
    do method = 2, 3
        sparse_solve_method = method
        do low = 0, 1
            sparse_low_memory = low == 1
            sparse_host_memory_limit_bytes = 1048576_int64
            b = rhs
            call sparse_solve(a, b, 1)
            call close (b, rhs, 'factor-only RHS')
            call sparse_solve(a, b, 2)
            call close (b, exact, 'bounded real solve')
            call close (matmul(a, b), rhs, 'original real residual')
            ! Only LOW_MEMORY changes: the supplied matrix must be refactored.
            ! A stale factor would return 2*exact and fail this original-A oracle.
            ! A tighter-cap change is separately checked as a terminal error.
            sparse_low_memory = .not. sparse_low_memory
            do j = 1, 2
                block(:, j) = 2*j*rhs
            end do
            call sparse_solve(2*a, block, 2)
            do j = 1, 2
                call close (matmul(2*a, block(:, j)), 2*j*rhs, &
                            'LOW-only real refactor residual')
            end do
            call sparse_solve(a, b, 3)
            zb = zr
            call sparse_solve(z, zb, 1)
            if (maxval(abs(zb - zr)) > 1d-12) error stop 'complex factor-only RHS'
            call sparse_solve(z, zb, 2)
            if (maxval(abs(zb - truth)) > 1d-12) error stop 'bounded complex solve'
            sparse_low_memory = .not. sparse_low_memory
            sparse_host_memory_limit_bytes = 1048576_int64
            do j = 1, 2
                zblock(:, j) = 2*j*zr
            end do
            call sparse_solve(2*z, zblock, 2)
            do j = 1, 2
                if (maxval(abs(matmul(2*z, zblock(:, j)) - 2*j*zr)) > 1d-12) &
                    error stop 'LOW-only complex refactor residual'
            end do
            call sparse_solve(z, zb, 3)
        end do
    end do
    sparse_host_memory_limit_bytes = 0_int64
    sparse_low_memory = .false.
    b = rhs
    call sparse_solve(a, b)
    call close (b, exact, 'unlimited recovery after free')
    print '(a)', 'All tests passed!'
contains
    subroutine close (actual, expected, description)
        real(dp), intent(in) :: actual(:), expected(:)
        character(*), intent(in) :: description
        if (maxval(abs(actual - expected)) > 1d-12) then
            print *, description
            error stop 'bounded native residual oracle'
        end if
    end subroutine close
end program test_native_sparse_memory
