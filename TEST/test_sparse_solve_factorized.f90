program test_sparse_solve_factorized
  ! Checks sparse_solve_factorized, the state-free solve with an existing
  ! SuiteSparse factorization, when several right-hand sides are solved
  ! concurrently from an OpenMP loop.
  ! Oracles: (1) the manufactured exact solution x_true of A*x = A*x_true,
  ! (2) bitwise equality to the serial sparse_solve(..., iopt=2) path.
  use sparse_mod, only: sparse_solve, sparse_solve_method, sparse_solve_factorized
  implicit none

  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: n = 4000, nrhs = 8
  integer :: nz, i, j, k, nfail
  integer, allocatable :: irow(:), pcol(:)
  real(dp), allocatable :: val(:), xtrue(:,:), b(:,:), xser(:,:), xpar(:,:)
  real(dp) :: dummy(n), err

  ! Nonsymmetric, periodic, diagonally dominant "advection-diffusion" matrix
  ! in compressed-column form: columns j hold rows j-1, j, j+1 (cyclic).
  nz = 3*n
  allocate(irow(nz), pcol(n+1), val(nz))
  k = 0
  do j = 1, n
    pcol(j) = k + 1
    do i = j - 1, j + 1
      k = k + 1
      irow(k) = modulo(i - 1, n) + 1
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
  call sort_columns()

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

contains

  subroutine sort_columns()
    ! Row indices within each column must be increasing for UMFPACK.
    integer :: c, p, q, itmp
    real(dp) :: vtmp
    do c = 1, n
      do p = pcol(c), pcol(c+1) - 1
        do q = p + 1, pcol(c+1) - 1
          if (irow(q) < irow(p)) then
            itmp = irow(p); irow(p) = irow(q); irow(q) = itmp
            vtmp = val(p); val(p) = val(q); val(q) = vtmp
          end if
        end do
      end do
    end do
  end subroutine sort_columns

end program test_sparse_solve_factorized
