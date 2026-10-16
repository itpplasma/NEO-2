program test_sparse_mumps
    ! sparse_solve_method = 4 against a manufactured solution: real and complex,
    ! one-shot, factor reuse and multi-RHS block solves through sparse_mod.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use mpi, only: MPI_Init, MPI_Finalize
    use sparse_mod, only: sparse_solve, sparse_solve_method
    implicit none

    integer, parameter :: n = 5, nz = 13, nrhs = 3
    ! Nonsymmetric tridiagonal matrix plus one corner entry, CSC, 1-based.
    integer, parameter :: pcol(n + 1) = [1, 3, 6, 9, 12, 14]
    integer, parameter :: irow(nz) = [1, 2, 1, 2, 3, 2, 3, 4, 3, 4, 5, 1, 5]
    real(dp), parameter :: aval(nz) = [4._dp, -1._dp, 2._dp, 5._dp, -2._dp, &
                                       1._dp, 6._dp, 3._dp, -1._dp, 7._dp, 2._dp, &
                                       0.5_dp, 8._dp]
    real(dp), parameter :: tol = 1.0e-12_dp
    complex(dp), parameter :: zshift = (0._dp, 1._dp)
    real(dp) :: x(n, nrhs), b(n, nrhs), b1(n)
    complex(dp) :: zval(nz), zx(n, nrhs), zb(n, nrhs), zb1(n)
    integer :: j, k, ierr

    call MPI_Init(ierr)
    sparse_solve_method = 4

    do k = 1, nrhs
        do j = 1, n
            x(j, k) = real(j + 2*k, dp)/real(n, dp)
        end do
    end do
    zval = aval*(1._dp + 0.25_dp*zshift)
    zval(nz) = zval(nz) + zshift
    zx = cmplx(x, -2._dp*x, kind=dp)
    b = matmul_csc(aval, x)
    zb = zmatmul_csc(zval, zx)

    ! Real one-shot
    b1 = b(:, 1)
    call sparse_solve(n, n, nz, irow, pcol, aval, b1, 0)
    call check(maxval(abs(b1 - x(:, 1))), 'real one-shot')

    ! Real factor once, reuse for a vector and a block
    b1 = b(:, 2)
    call sparse_solve(n, n, nz, irow, pcol, aval, b1, 1)
    call sparse_solve(n, n, nz, irow, pcol, aval, b1, 2)
    call check(maxval(abs(b1 - x(:, 2))), 'real reuse vector')
    call sparse_solve(n, n, nz, irow, pcol, aval, b, 2)
    call check(maxval(abs(b - x)), 'real reuse block')
    call sparse_solve(n, n, nz, irow, pcol, aval, b1, 3)

    ! Complex one-shot
    zb1 = zb(:, 3)
    call sparse_solve(n, n, nz, irow, pcol, zval, zb1, 0)
    call check(maxval(abs(zb1 - zx(:, 3))), 'complex one-shot')

    ! Complex factor once, reuse for a block and a vector
    call sparse_solve(n, n, nz, irow, pcol, zval, zb1, 1)
    call sparse_solve(n, n, nz, irow, pcol, zval, zb, 2)
    call check(maxval(abs(zb - zx)), 'complex reuse block')
    zb = zmatmul_csc(zval, zx)
    zb1 = zb(:, 1)
    call sparse_solve(n, n, nz, irow, pcol, zval, zb1, 2)
    call check(maxval(abs(zb1 - zx(:, 1))), 'complex reuse vector')
    call sparse_solve(n, n, nz, irow, pcol, zval, zb1, 3)

    call MPI_Finalize(ierr)
    print *, 'All tests passed!'

contains

    function matmul_csc(val, xin) result(y)
        real(dp), intent(in) :: val(:), xin(:, :)
        real(dp) :: y(n, size(xin, 2))
        integer :: col, p

        y = 0._dp
        do col = 1, n
            do p = pcol(col), pcol(col + 1) - 1
                y(irow(p), :) = y(irow(p), :) + val(p)*xin(col, :)
            end do
        end do
    end function matmul_csc

    function zmatmul_csc(val, xin) result(y)
        complex(dp), intent(in) :: val(:), xin(:, :)
        complex(dp) :: y(n, size(xin, 2))
        integer :: col, p

        y = (0._dp, 0._dp)
        do col = 1, n
            do p = pcol(col), pcol(col + 1) - 1
                y(irow(p), :) = y(irow(p), :) + val(p)*xin(col, :)
            end do
        end do
    end function zmatmul_csc

    subroutine check(err, label)
        real(dp), intent(in) :: err
        character(*), intent(in) :: label

        if (.not. err <= tol) then
            print *, 'FAIL: ', label, ' error ', err
            error stop 1
        end if
    end subroutine check

end program test_sparse_mumps
