program test_sparse_solve_reuse
    ! Factorize once (iopt=1), solve with the stored factors (iopt=2), free
    ! (iopt=3). Oracles: manufactured solution and bitwise equality with the
    ! one-shot path (iopt=0), for real/complex and 1-D/2-D right-hand sides.
    use sparse_mod, only: sparse_solve, sparse_solve_method
    implicit none

    integer, parameter :: dp = kind(1.0d0), n = 30, nrhs = 2
    real(dp), parameter :: tol = 1.0e-12_dp
    complex(dp) :: a(n, n), xt(n, nrhs)
    integer, allocatable :: irow(:), pcol(:)
    real(dp), allocatable :: val(:)
    complex(dp), allocatable :: valc(:)
    integer :: i, j, k, nz, method
    logical :: ok

    ! Nonsymmetric banded matrix, diagonally dominant, no zero real parts.
    a = (0.0_dp, 0.0_dp)
    do j = 1, n
        do i = max(1, j - 1), min(n, j + 3)
            a(i, j) = cmplx(-1.0_dp + 0.3_dp*sin(real(i + 2*j, dp)), &
                            0.2_dp*cos(real(i - j, dp)), dp)
        end do
        a(j, j) = cmplx(8.0_dp + 0.1_dp*j, 0.5_dp, dp)
        xt(j, 1) = cmplx(cos(0.3_dp*j), sin(0.2_dp*j), dp)
        xt(j, 2) = cmplx(1.0_dp/j, -0.1_dp*j, dp)
    end do

    nz = count(a /= (0.0_dp, 0.0_dp))
    allocate (irow(nz), pcol(n + 1), valc(nz))
    k = 0
    do j = 1, n
        pcol(j) = k + 1
        do i = 1, n
            if (a(i, j) /= (0.0_dp, 0.0_dp)) then
                k = k + 1
                irow(k) = i
                valc(k) = a(i, j)
            end if
        end do
    end do
    pcol(n + 1) = k + 1
    val = real(valc, dp)

    ok = .true.
    do method = 2, 3
        sparse_solve_method = method
        call check_real(ok)
        call check_complex(ok)
    end do

    if (ok) then
        print *, 'All tests passed!'
    else
        print *, 'FAIL'
        stop 1
    end if

contains

    subroutine check_real(ok)
        logical, intent(inout) :: ok
        real(dp) :: b(n, nrhs), x0(n, nrhs), x(n, nrhs), y(n)
        integer :: r

        b = matmul(real(a, dp), real(xt, dp))
        do r = 1, nrhs
            y = b(:, r)
            call sparse_solve(n, n, nz, irow, pcol, val, y, 0)
            x0(:, r) = y
        end do

        y = b(:, 1)
        call sparse_solve(n, n, nz, irow, pcol, val, y, 1)
        do r = 1, nrhs
            y = b(:, r)
            call sparse_solve(n, n, nz, irow, pcol, val, y, 2)
            x(:, r) = y
        end do
        ok = ok .and. all(x == x0)
        ok = ok .and. maxval(abs(x - real(xt, dp))) < tol

        x = b
        call sparse_solve(n, n, nz, irow, pcol, val, x, 2)
        ok = ok .and. all(x == x0)
        call sparse_solve(n, n, nz, irow, pcol, val, y, 3)
    end subroutine check_real

    subroutine check_complex(ok)
        logical, intent(inout) :: ok
        complex(dp) :: b(n, nrhs), x0(n, nrhs), x(n, nrhs), y(n)
        integer :: r

        b = matmul(a, xt)
        do r = 1, nrhs
            y = b(:, r)
            call sparse_solve(n, n, nz, irow, pcol, valc, y, 0)
            x0(:, r) = y
        end do

        y = b(:, 1)
        call sparse_solve(n, n, nz, irow, pcol, valc, y, 1)
        do r = 1, nrhs
            y = b(:, r)
            call sparse_solve(n, n, nz, irow, pcol, valc, y, 2)
            x(:, r) = y
        end do
        ok = ok .and. all(x == x0)
        ok = ok .and. maxval(abs(x - xt)) < tol

        x = b
        call sparse_solve(n, n, nz, irow, pcol, valc, x, 2)
        ok = ok .and. all(x == x0)
        call sparse_solve(n, n, nz, irow, pcol, valc, y, 3)
    end subroutine check_complex

end program test_sparse_solve_reuse
