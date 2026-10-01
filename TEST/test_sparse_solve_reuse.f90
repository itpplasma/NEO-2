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

end program test_sparse_solve_reuse
