program test_sparse_solve_reuse
  ! Factorize once (iopt=1), solve repeatedly (iopt=2), free (iopt=3) with
  ! the SuiteSparse back end of sparse_mod.  Oracles: a manufactured solution
  ! x_true with b = A*x_true computed here from the CSC data, and bitwise
  ! equality with the one-shot path (iopt=0) for every right-hand side.
  use sparse_mod, only: sparse_solve, sparse_solve_method
  implicit none

  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: n = 40, nrhs = 3
  real(dp), parameter :: tol = 1.0d-12
  integer, allocatable :: irow(:), pcol(:)
  real(dp), allocatable :: val(:)
  complex(dp), allocatable :: valc(:)
  integer :: nz, method, imeth, nfail

  nfail = 0
  call build_matrix(irow, pcol, val, valc, nz)

  do imeth = 1, 2
     method = merge(3, 2, imeth == 1)
     sparse_solve_method = method
     call check_real(method)
     call check_complex(method)
        call check_real_lifecycle(method)
        call check_complex_lifecycle(method)
  end do

  if (nfail == 0) then
     print *, 'All tests passed!'
  else
     print *, 'FAIL: ', nfail, ' checks failed'
     stop 1
  end if

contains

  subroutine build_matrix(irow, pcol, val, valc, nz)
    ! Nonsymmetric, diagonally dominant: diagonal plus neighbours at
    ! distance 1 and 7 (wrapped), stored column by column.
    integer, allocatable, intent(out) :: irow(:), pcol(:)
    real(dp), allocatable, intent(out) :: val(:)
    complex(dp), allocatable, intent(out) :: valc(:)
    integer, intent(out) :: nz
    integer :: j, k, i, rows(5)

    nz = 5*n
    allocate(irow(nz), pcol(n+1), val(nz), valc(nz))
    k = 0
    do j = 1, n
       pcol(j) = k + 1
       rows = [modulo(j-8, n)+1, modulo(j-2, n)+1, j, modulo(j, n)+1, modulo(j+6, n)+1]
       call sort5(rows)
       do i = 1, 5
          k = k + 1
          irow(k) = rows(i)
          if (rows(i) == j) then
             val(k) = 10.0_dp + 0.1_dp*j
          else
             val(k) = -1.0_dp + 0.37_dp*sin(real(rows(i) + 3*j, dp))
          end if
          valc(k) = cmplx(val(k), 0.5_dp*cos(real(2*rows(i) - j, dp)), dp)
       end do
    end do
    pcol(n+1) = k + 1
  end subroutine build_matrix

  subroutine sort5(a)
    integer, intent(inout) :: a(5)
    integer :: i, j, t
    do i = 2, 5
       t = a(i)
       j = i - 1
       do while (j >= 1)
          if (a(j) <= t) exit
          a(j+1) = a(j)
          j = j - 1
       end do
       a(j+1) = t
    end do
  end subroutine sort5

  function matvec_real(x) result(y)
    real(dp), intent(in) :: x(:)
    real(dp) :: y(n)
    integer :: j, k
    y = 0.0_dp
    do j = 1, n
       do k = pcol(j), pcol(j+1) - 1
          y(irow(k)) = y(irow(k)) + val(k)*x(j)
       end do
    end do
  end function matvec_real

  function matvec_cmplx(x) result(y)
    complex(dp), intent(in) :: x(:)
    complex(dp) :: y(n)
    integer :: j, k
    y = (0.0_dp, 0.0_dp)
    do j = 1, n
       do k = pcol(j), pcol(j+1) - 1
          y(irow(k)) = y(irow(k)) + valc(k)*x(j)
       end do
    end do
  end function matvec_cmplx

  subroutine expect(ok, what, method)
    logical, intent(in) :: ok
    character(*), intent(in) :: what
    integer, intent(in) :: method
    if (.not. ok) then
       print *, 'FAIL: ', what, ' (sparse_solve_method=', method, ')'
       nfail = nfail + 1
    end if
  end subroutine expect

  subroutine check_real(method)
    integer, intent(in) :: method
    real(dp) :: xt(n, nrhs), b(n, nrhs), x1(n), x0(n), x2(n, nrhs), x20(n, nrhs)
    integer :: r, j

    do r = 1, nrhs
       do j = 1, n
          xt(j, r) = cos(0.3_dp*j*r) + 0.01_dp*r
       end do
       b(:, r) = matvec_real(xt(:, r))
    end do

    ! one-shot reference results
    do r = 1, nrhs
       x0 = b(:, r)
       call sparse_solve(n, n, nz, irow, pcol, val, x0, 0)
       x20(:, r) = x0
    end do

    ! factorize once, then repeated 1-D solves
    x1 = b(:, 1)
    call sparse_solve(n, n, nz, irow, pcol, val, x1, 1)
    do r = 1, nrhs
       x1 = b(:, r)
       call sparse_solve(n, n, nz, irow, pcol, val, x1, 2)
       call expect(maxval(abs(x1 - xt(:, r))) < tol*maxval(abs(xt(:, r))), &
            'real b1 manufactured solution', method)
       call expect(all(x1 == x20(:, r)), 'real b1 bitwise iopt=2 vs iopt=0', method)
    end do

    ! repeated 2-D solve with the same factors
    x2 = b
    call sparse_solve(n, n, nz, irow, pcol, val, x2, 2)
    call expect(all(x2 == x20), 'real b2 bitwise iopt=2 vs iopt=0', method)
    call sparse_solve(n, n, nz, irow, pcol, val, x1, 3)
  end subroutine check_real

  subroutine check_complex(method)
    integer, intent(in) :: method
    complex(dp) :: xt(n, nrhs), b(n, nrhs), x1(n), x0(n), x2(n, nrhs), x20(n, nrhs)
    integer :: r, j

    do r = 1, nrhs
       do j = 1, n
          xt(j, r) = cmplx(cos(0.3_dp*j*r), sin(0.2_dp*j + r), dp)
       end do
       b(:, r) = matvec_cmplx(xt(:, r))
    end do

    do r = 1, nrhs
       x0 = b(:, r)
       call sparse_solve(n, n, nz, irow, pcol, valc, x0, 0)
       x20(:, r) = x0
    end do

    x1 = b(:, 1)
    call sparse_solve(n, n, nz, irow, pcol, valc, x1, 1)
    do r = 1, nrhs
       x1 = b(:, r)
       call sparse_solve(n, n, nz, irow, pcol, valc, x1, 2)
       call expect(maxval(abs(x1 - xt(:, r))) < tol*maxval(abs(xt(:, r))), &
            'complex b1 manufactured solution', method)
       call expect(all(x1 == x20(:, r)), 'complex b1 bitwise iopt=2 vs iopt=0', method)
    end do

    x2 = b
    call sparse_solve(n, n, nz, irow, pcol, valc, x2, 2)
    call expect(all(x2 == x20), 'complex b2 bitwise iopt=2 vs iopt=0', method)
    call sparse_solve(n, n, nz, irow, pcol, valc, x1, 3)
  end subroutine check_complex

    subroutine lifecycle_matrix(m, rows, pointers, values, valuesc, dense, densec)
        integer, intent(in) :: m
        integer, allocatable, intent(out) :: rows(:), pointers(:)
        real(dp), allocatable, intent(out) :: values(:), dense(:, :)
        complex(dp), allocatable, intent(out) :: valuesc(:), densec(:, :)
        integer :: i, j, k

        allocate(rows(3*m - 2), pointers(m + 1), values(3*m - 2), &
                 valuesc(3*m - 2), dense(m, m), densec(m, m))
        dense = 0.0_dp
        densec = (0.0_dp, 0.0_dp)
        k = 0
        do j = 1, m
            pointers(j) = k + 1
            do i = max(1, j - 1), min(m, j + 1)
                k = k + 1
                rows(k) = i
                if (i == j) then
                    values(k) = 4.0_dp + 0.1_dp*j
                else if (i < j) then
                    values(k) = -0.6_dp
                else
                    values(k) = -1.1_dp
                end if
                valuesc(k) = cmplx(values(k), 0.05_dp*(2*i - j), dp)
                dense(i, j) = values(k)
                densec(i, j) = valuesc(k)
            end do
        end do
        pointers(m + 1) = k + 1
    end subroutine lifecycle_matrix

    subroutine check_real_lifecycle(method)
        integer, intent(in) :: method
        integer, parameter :: sizes(4) = [5, 12, 4, 9]
        integer, parameter :: counts(4) = [2, 4, 1, 3]
        integer :: stage, m, nb, nzloc, i, r, solve_method
        integer, allocatable :: rows(:), pointers(:)
        real(dp), allocatable :: values(:), dense(:, :), xt(:, :), b(:, :)
        real(dp), allocatable :: x1(:), x2(:, :)
        complex(dp), allocatable :: valuesc(:), densec(:, :)

        do stage = 1, size(sizes)
            m = sizes(stage)
            nb = counts(stage)
            call lifecycle_matrix(m, rows, pointers, values, valuesc, dense, densec)
            nzloc = size(values)
            if (allocated(xt)) deallocate(xt, b, x1, x2)
            allocate(xt(m, nb), b(m, nb), x1(m), x2(m, nb))
            do r = 1, nb
                do i = 1, m
                    xt(i, r) = cos(0.23_dp*i*r) + 0.01_dp*m
                end do
            end do
            ! Dense multiplication supplies the manufactured oracle independently
            ! of the sparse solve. Refactorization changes both matrix and shape.
            b = matmul(dense, xt)
            sparse_solve_method = method
            x2 = b
            call sparse_solve(m, m, nzloc, rows, pointers, values, x2, 1)
            do solve_method = 2, 3
                sparse_solve_method = solve_method
                x1 = b(:, 1)
                call sparse_solve(m, m, nzloc, rows, pointers, values, x1, 2)
                call expect(all(abs(x1 - xt(:, 1)) < &
                                tol*maxval(abs(xt(:, 1)))), &
                            'real lifecycle b1 manufactured solution', solve_method)
                x2 = b
                call sparse_solve(m, m, nzloc, rows, pointers, values, x2(:, 1:1), 2)
                call expect(all(abs(x2(:, 1) - xt(:, 1)) < &
                                tol*maxval(abs(xt(:, 1)))), &
                            'real lifecycle single-column b2', solve_method)
                x2 = b
                call sparse_solve(m, m, nzloc, rows, pointers, values, x2, 2)
                call expect(all(abs(x2 - xt) < tol*maxval(abs(xt))), &
                            'real lifecycle multiple-column b2', solve_method)
            end do
        end do

        ! Free b2-created factors through b1; then create factors automatically
        ! through each rank of RHS and free them through the opposite rank.
        call sparse_solve(m, m, nzloc, rows, pointers, values, x1, 3)
        sparse_solve_method = method
        x1 = b(:, 1)
        call sparse_solve(m, m, nzloc, rows, pointers, values, x1, 2)
        call expect(all(abs(x1 - xt(:, 1)) < tol*maxval(abs(xt(:, 1)))), &
                    'real lifecycle automatic b1 factorization', method)
        call sparse_solve(m, m, nzloc, rows, pointers, values, x2, 3)
        x2 = b
        call sparse_solve(m, m, nzloc, rows, pointers, values, x2, 2)
        call expect(all(abs(x2 - xt) < tol*maxval(abs(xt))), &
                    'real lifecycle automatic b2 factorization', method)
        call sparse_solve(m, m, nzloc, rows, pointers, values, x1, 3)
        x2 = b
        call sparse_solve(m, m, nzloc, rows, pointers, values, x2, 0)
        call expect(all(abs(x2 - xt) < tol*maxval(abs(xt))), &
                    'real lifecycle one-shot b2', method)
    end subroutine check_real_lifecycle

    subroutine check_complex_lifecycle(method)
        integer, intent(in) :: method
        integer, parameter :: sizes(4) = [5, 12, 4, 9]
        integer, parameter :: counts(4) = [2, 4, 1, 3]
        integer :: stage, m, nb, nzloc, i, r, solve_method
        integer, allocatable :: rows(:), pointers(:)
        real(dp), allocatable :: values(:), dense(:, :)
        complex(dp), allocatable :: valuesc(:), densec(:, :), xt(:, :), b(:, :)
        complex(dp), allocatable :: x1(:), x2(:, :)

        do stage = 1, size(sizes)
            m = sizes(stage)
            nb = counts(stage)
            call lifecycle_matrix(m, rows, pointers, values, valuesc, dense, densec)
            nzloc = size(valuesc)
            if (allocated(xt)) deallocate(xt, b, x1, x2)
            allocate(xt(m, nb), b(m, nb), x1(m), x2(m, nb))
            do r = 1, nb
                do i = 1, m
                    xt(i, r) = cmplx(cos(0.23_dp*i*r), sin(0.17_dp*i + r), dp)
                end do
            end do
            b = matmul(densec, xt)
            sparse_solve_method = method
            x2 = b
            call sparse_solve(m, m, nzloc, rows, pointers, valuesc, x2, 1)
            do solve_method = 2, 3
                sparse_solve_method = solve_method
                x1 = b(:, 1)
                call sparse_solve(m, m, nzloc, rows, pointers, valuesc, x1, 2)
                call expect(all(abs(x1 - xt(:, 1)) < &
                                tol*maxval(abs(xt(:, 1)))), &
                            'complex lifecycle b1 manufactured solution', solve_method)
                x2 = b
                call sparse_solve(m, m, nzloc, rows, pointers, valuesc, x2(:, 1:1), 2)
                call expect(all(abs(x2(:, 1) - xt(:, 1)) < &
                                tol*maxval(abs(xt(:, 1)))), &
                            'complex lifecycle single-column b2', solve_method)
                x2 = b
                call sparse_solve(m, m, nzloc, rows, pointers, valuesc, x2, 2)
                call expect(all(abs(x2 - xt) < tol*maxval(abs(xt))), &
                            'complex lifecycle multiple-column b2', solve_method)
            end do
        end do

        call sparse_solve(m, m, nzloc, rows, pointers, valuesc, x1, 3)
        sparse_solve_method = method
        x1 = b(:, 1)
        call sparse_solve(m, m, nzloc, rows, pointers, valuesc, x1, 2)
        call expect(all(abs(x1 - xt(:, 1)) < tol*maxval(abs(xt(:, 1)))), &
                    'complex lifecycle automatic b1 factorization', method)
        call sparse_solve(m, m, nzloc, rows, pointers, valuesc, x2, 3)
        x2 = b
        call sparse_solve(m, m, nzloc, rows, pointers, valuesc, x2, 2)
        call expect(all(abs(x2 - xt) < tol*maxval(abs(xt))), &
                    'complex lifecycle automatic b2 factorization', method)
        call sparse_solve(m, m, nzloc, rows, pointers, valuesc, x1, 3)
        x2 = b
        call sparse_solve(m, m, nzloc, rows, pointers, valuesc, x2, 0)
        call expect(all(abs(x2 - xt) < tol*maxval(abs(xt))), &
                    'complex lifecycle one-shot b2', method)
    end subroutine check_complex_lifecycle

end program test_sparse_solve_reuse
