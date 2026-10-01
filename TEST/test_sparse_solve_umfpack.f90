program test_sparse_solve_umfpack
  ! Solve nonsymmetric sparse systems with a manufactured solution through
  ! the UMFPACK path of sparse_mod (factorize once, solve repeatedly, free),
  ! for real and complex matrices. The matrix is a 2D upwinded
  ! convection-diffusion stencil, periodic in one direction like a closed
  ! field line, so the fill-reducing ordering has real work to do.
  use sparse_mod, only: sparse_solve, sparse_solve_method, column_full2pointer
  implicit none

  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: nx = 60, ny = 40, n = nx*ny
  integer, allocatable :: irow(:), icol(:), pcol(:)
  real(dp), allocatable :: aval(:), x_true(:), b(:)
  complex(dp), allocatable :: zval(:), z_true(:), zb(:)
  integer :: nz, k
  real(dp) :: err
  integer :: m, rows(5)
  real(dp) :: vals(5)

  sparse_solve_method = 3
  call build_matrix(irow, icol, aval, nz)
  call column_full2pointer(icol, pcol)

  allocate(x_true(n), b(n))
  do k = 1, n
     x_true(k) = sin(0.37_dp*k) + 0.1_dp*k/n
  end do

  ! factorize, then two solves with the same factor
  b = matvec(aval, x_true)
  call sparse_solve(n, n, nz, irow, pcol, aval, b, 1)
  call sparse_solve(n, n, nz, irow, pcol, aval, b, 2)
  err = maxval(abs(b - x_true))/maxval(abs(x_true))
  call check('real solve 1', err)
  b = matvec(aval, 2.0_dp*x_true)
  call sparse_solve(n, n, nz, irow, pcol, aval, b, 2)
  err = maxval(abs(b - 2.0_dp*x_true))/maxval(abs(2.0_dp*x_true))
  call check('real solve 2', err)
  call sparse_solve(n, n, nz, irow, pcol, aval, b, 3)

  allocate(zval(nz), z_true(n), zb(n))
  zval = cmplx(aval, 0.0_dp, dp)
  do k = 1, nz
     if (irow(k) == icol(k)) zval(k) = zval(k) + cmplx(0.0_dp, 0.5_dp, dp)
  end do
  z_true = cmplx(x_true, cos(0.11_dp*[(k, k = 1, n)]), dp)
  zb = zmatvec(zval, z_true)
  call sparse_solve(n, n, nz, irow, pcol, zval, zb, 1)
  call sparse_solve(n, n, nz, irow, pcol, zval, zb, 2)
  err = maxval(abs(zb - z_true))/maxval(abs(z_true))
  call check('complex solve', err)
  call sparse_solve(n, n, nz, irow, pcol, zval, zb, 3)

  print *, 'All tests passed!'

contains

  integer function idx(i, j)
    integer, intent(in) :: i, j
    idx = modulo(i - 1, nx) + 1 + (j - 1)*nx
  end function idx

  subroutine build_matrix(irow, icol, aval, nz)
    ! Columns sorted, rows sorted within each column (CSC order).
    integer, allocatable, intent(out) :: irow(:), icol(:)
    real(dp), allocatable, intent(out) :: aval(:)
    integer, intent(out) :: nz
    integer :: i, j, c, r, p(5)

    allocate(irow(5*n), icol(5*n), aval(5*n))
    nz = 0
    do j = 1, ny
       do i = 1, nx
          c = idx(i, j)
          ! entries A(r, c): column c of the transpose stencil
          m = 0
          call add(idx(i - 1, j), -1.0_dp)
          call add(idx(i + 1, j), -1.0_dp - 0.8_dp)
          if (j > 1) call add(idx(i, j - 1), -1.0_dp)
          if (j < ny) call add(idx(i, j + 1), -1.0_dp)
          call add(c, 4.0_dp + 0.8_dp + 0.05_dp)
          p(1:m) = sort_index(rows(1:m))
          do r = 1, m
             nz = nz + 1
             irow(nz) = rows(p(r)); icol(nz) = c; aval(nz) = vals(p(r))
          end do
       end do
    end do
    irow = irow(1:nz); icol = icol(1:nz); aval = aval(1:nz)
  end subroutine build_matrix

  subroutine add(row, v)
    integer, intent(in) :: row
    real(dp), intent(in) :: v
    m = m + 1; rows(m) = row; vals(m) = v
  end subroutine add

  function sort_index(a) result(p)
    integer, intent(in) :: a(:)
    integer :: p(size(a)), i, j, t
    p = [(i, i = 1, size(a))]
    do i = 2, size(a)
       j = i
       do while (j > 1)
          if (a(p(j - 1)) <= a(p(j))) exit
          t = p(j); p(j) = p(j - 1); p(j - 1) = t
          j = j - 1
       end do
    end do
  end function sort_index

  function matvec(v, x) result(y)
    real(dp), intent(in) :: v(:), x(:)
    real(dp) :: y(n)
    integer :: k
    y = 0.0_dp
    do k = 1, nz
       y(irow(k)) = y(irow(k)) + v(k)*x(icol(k))
    end do
  end function matvec

  function zmatvec(v, x) result(y)
    complex(dp), intent(in) :: v(:), x(:)
    complex(dp) :: y(n)
    integer :: k
    y = (0.0_dp, 0.0_dp)
    do k = 1, nz
       y(irow(k)) = y(irow(k)) + v(k)*x(icol(k))
    end do
  end function zmatvec

  subroutine check(label, err)
    character(*), intent(in) :: label
    real(dp), intent(in) :: err
    print '(A,": max rel error ",ES10.3)', label, err
    if (.not. (err < 1.0e-11_dp)) then
       print *, 'FAIL: ', label
       error stop 1
    end if
  end subroutine check

end program test_sparse_solve_umfpack
