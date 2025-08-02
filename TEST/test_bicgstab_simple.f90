program test_bicgstab_simple
  use nrtype, only: I4B, DP
  use sparse_solvers_mod, only: bicgstab_abs_tolerance, bicgstab_rel_tolerance, &
                                bicgstab_max_iter, sparse_solve_method, SOLVER_BICGSTAB, &
                                SOLVER_UMFPACK
  use sparse_mod, only: sparse_solve
  implicit none
  
  ! Simple 3x3 test matrix: A = [2 -1  0]
  !                             [-1  2 -1] 
  !                             [ 0 -1  2]
  ! With RHS b = [1, 0, 1], exact solution x = [1, 1, 1]
  
  integer(I4B), parameter :: n = 3
  integer(I4B), parameter :: nz = 7
  
  ! CSC format arrays
  integer(I4B) :: irow(nz) = [1, 2, 1, 2, 3, 2, 3]  ! Row indices
  integer(I4B) :: pcol(n+1) = [1, 3, 6, 8]          ! Column pointers
  real(DP) :: amat(nz) = [2.0_DP, -1.0_DP, -1.0_DP, 2.0_DP, -1.0_DP, -1.0_DP, 2.0_DP]
  
  real(DP) :: b(n) = [1.0_DP, 0.0_DP, 1.0_DP]       ! RHS
  real(DP) :: x_exact(n) = [1.0_DP, 1.0_DP, 1.0_DP] ! Exact solution
  real(DP) :: x_bicg(n), x_umf(n)
  real(DP) :: error_bicg, error_umf
  integer :: i
  
  ! Save current settings
  real(DP) :: saved_abs_tol, saved_rel_tol
  integer :: saved_max_iter, saved_method
  
  write(*,'(A)') '=== Simple BiCGSTAB Test ===  '
  write(*,'(A)') 'Matrix A = [2 -1  0]  (well-conditioned tridiagonal)'
  write(*,'(A)') '           [-1  2 -1]'
  write(*,'(A)') '           [ 0 -1  2]'
  write(*,'(A)') 'RHS b = [1, 0, 1]'
  write(*,'(A)') 'Exact solution = [1, 1, 1]'
  write(*,'(A)') ''
  
  saved_abs_tol = bicgstab_abs_tolerance
  saved_rel_tol = bicgstab_rel_tolerance
  saved_max_iter = bicgstab_max_iter
  saved_method = sparse_solve_method
  
  ! Test with UMFPACK first (reference solution)
  sparse_solve_method = SOLVER_UMFPACK
  x_umf = b  ! Copy RHS
  call sparse_solve(n, n, nz, irow, pcol, amat, x_umf, 1)
  error_umf = maxval(abs(x_umf - x_exact))
  
  write(*,'(A)') 'UMFPACK solution:'
  write(*,'(A,3F10.6)') '  x = [', x_umf
  write(*,'(A,E12.4)') '  Max error = ', error_umf
  write(*,'(A)') ''
  
  ! Test with BiCGSTAB
  sparse_solve_method = SOLVER_BICGSTAB
  bicgstab_abs_tolerance = 1.0e-12_DP
  bicgstab_rel_tolerance = 1.0e-12_DP
  bicgstab_max_iter = 100
  
  x_bicg = b  ! Copy RHS
  call sparse_solve(n, n, nz, irow, pcol, amat, x_bicg, 1)
  error_bicg = maxval(abs(x_bicg - x_exact))
  
  write(*,'(A)') 'BiCGSTAB solution:'
  write(*,'(A,3F10.6)') '  x = [', x_bicg
  write(*,'(A,E12.4)') '  Max error = ', error_bicg
  write(*,'(A)') ''
  
  ! Restore settings
  bicgstab_abs_tolerance = saved_abs_tol
  bicgstab_rel_tolerance = saved_rel_tol
  bicgstab_max_iter = saved_max_iter
  sparse_solve_method = saved_method
  
  if (error_bicg < 1.0e-10_DP) then
    write(*,'(A)') 'BiCGSTAB test PASSED!'
  else
    write(*,'(A)') 'BiCGSTAB test FAILED!'
    stop 1
  end if
  
end program test_bicgstab_simple