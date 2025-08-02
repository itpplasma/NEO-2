program test_amg_simple_matrix
  use amg_precond_mod
  use amg_types_mod
  use sparse_types_mod, only: dp
  implicit none
  
  ! Test AMG on a simple tridiagonal matrix
  integer, parameter :: n = 25  ! Matrix size
  integer, parameter :: nnz = 73  ! 3*n-2 for tridiagonal
  integer :: row_ptr(n+1), col_idx(nnz)
  real(dp) :: values(nnz)
  type(amg_hierarchy) :: hier
  real(dp) :: b(n), x(n), r(n)
  integer :: i, j, k, idx
  
  write(*,*) '========================================='
  write(*,*) 'Simple AMG Test on 2D Poisson Matrix'
  write(*,*) '========================================='
  
  ! Build a simple tridiagonal matrix for testing
  ! 2 on diagonal, -1 on off-diagonals
  idx = 1
  do i = 1, n
    row_ptr(i) = idx
    
    ! Lower diagonal
    if (i > 1) then
      col_idx(idx) = i-1
      values(idx) = -1.0_dp
      idx = idx + 1
    end if
    
    ! Main diagonal
    col_idx(idx) = i
    values(idx) = 2.0_dp
    idx = idx + 1
    
    ! Upper diagonal
    if (i < n) then
      col_idx(idx) = i+1
      values(idx) = -1.0_dp
      idx = idx + 1
    end if
  end do
  row_ptr(n+1) = idx
  
  write(*,'(A,I0,A,I0,A,I0)') 'Matrix: ', n, 'x', n, ' with ', idx-1, ' nonzeros'
  
  ! Set up AMG
  call amg_precond_setup(hier, n, idx-1, row_ptr, col_idx, values)
  
  write(*,'(A,I0)') 'AMG levels built: ', hier%n_levels
  
  ! Test AMG application
  b = 1.0_dp  ! RHS = ones
  x = 0.0_dp  ! Initial guess = zeros
  
  write(*,*) 'Testing AMG preconditioner application...'
  call amg_precond_apply(hier, x, b)
  
  write(*,'(A,ES12.5)') 'Max |x|: ', maxval(abs(x))
  write(*,'(A,L1)') 'All finite: ', all(x == x .and. abs(x) < huge(1.0_dp))
  
  if (all(x == x .and. abs(x) < huge(1.0_dp))) then
    write(*,*) 'SUCCESS: AMG preconditioner works on simple matrix'
  else
    write(*,*) 'FAILURE: AMG preconditioner produces non-finite values'
  end if
  
  call amg_precond_destroy(hier)
  
end program test_amg_simple_matrix