program test_sparse_solve_factorized
  ! Solves several right-hand sides concurrently with sparse_solve_factorized.
  ! Oracles: the manufactured solution xtrue of A*x = A*xtrue, and bitwise
  ! equality to the serial sparse_solve(..., iopt=2) path.
  use sparse_mod, only: sparse_solve, sparse_solve_method, sparse_solve_factorized
  implicit none

  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: n = 4000, nrhs = 8
  integer :: nz, i, j, k, nfail
  integer, allocatable :: irow(:), pcol(:)
  real(dp), allocatable :: val(:), xtrue(:,:), b(:,:), xser(:,:), xpar(:,:)
  real(dp) :: dummy(n), err

  ! Nonsymmetric, diagonally dominant tridiagonal matrix, compressed columns.
  nz = 3*n - 2
  allocate(irow(nz), pcol(n+1), val(nz))
  k = 0
  do j = 1, n
    pcol(j) = k + 1
    do i = max(1, j - 1), min(n, j + 1)
      k = k + 1
      irow(k) = i
      if (i == j) then
        val(k) = 4.0_dp + 0.5_dp*sin(0.01_dp*j)
      else if (i < j) then
        val(k) = -1.3_dp
      else
        val(k) = -0.7_dp
      end if
    end do
  end do
  pcol(n+1) = nz + 1

  allocate(xtrue(n,nrhs), b(n,nrhs), xser(n,nrhs), xpar(n,nrhs))
  do j = 1, nrhs
    do i = 1, n
      xtrue(i,j) = cos(0.003_dp*i*j) + real(j, dp)
    end do
  end do
  b = 0.0_dp
  do j = 1, n
    do k = pcol(j), pcol(j+1) - 1
      b(irow(k),:) = b(irow(k),:) + val(k)*xtrue(j,:)
    end do
  end do

  sparse_solve_method = 3
  dummy = 0.0_dp
  call sparse_solve(n, n, nz, irow, pcol, val, dummy, 1)

  xser = b
  do j = 1, nrhs
    call sparse_solve(n, n, nz, irow, pcol, val, xser(:,j), 2)
  end do

  xpar = b
  !$omp parallel do num_threads(4) schedule(static,1)
  do j = 1, nrhs
    call sparse_solve_factorized(xpar(:,j))
  end do
  !$omp end parallel do

  call sparse_solve(n, n, nz, irow, pcol, val, dummy, 3)

  nfail = 0
  do j = 1, nrhs
    err = maxval(abs(xpar(:,j) - xtrue(:,j)))/maxval(abs(xtrue(:,j)))
    if (err > 1.0e-12_dp) then
      print *, 'FAIL: rhs ', j, ' relative error to exact solution ', err
      nfail = nfail + 1
    end if
    if (any(xpar(:,j) /= xser(:,j))) then
      print *, 'FAIL: rhs ', j, ' concurrent solve differs from serial iopt=2 solve'
      nfail = nfail + 1
    end if
  end do

  if (nfail == 0) print *, 'All tests passed!'

end program test_sparse_solve_factorized
