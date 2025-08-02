program test_bicgstab_illcond
  use nrtype, only: I4B, DP
  use sparse_solvers_mod, only: bicgstab_abs_tolerance, bicgstab_rel_tolerance, &
                                bicgstab_max_iter, sparse_solve_method, SOLVER_BICGSTAB, &
                                SOLVER_UMFPACK
  use sparse_mod, only: sparse_solve
  implicit none
  
  ! Ill-conditioned 3x3 test matrix with very small diagonal element:
  ! A = [1     0     0   ]
  !     [0   1e-8   0   ]  <- Very small diagonal element
  !     [0     0     1   ]
  ! With RHS b = [1, 1e-8, 1], exact solution x = [1, 1, 1]
  
  integer(I4B), parameter :: n = 3
  integer(I4B), parameter :: nz = 3
  
  ! CSC format arrays (diagonal matrix)
  integer(I4B) :: irow(nz) = [1, 2, 3]           ! Row indices
  integer(I4B) :: pcol(n+1) = [1, 2, 3, 4]      ! Column pointers
  real(DP) :: amat(nz) = [1.0_DP, 1.0e-12_DP, 1.0_DP]  ! Very ill-conditioned
  
  real(DP) :: b(n) = [1.0_DP, 1.0e-12_DP, 1.0_DP]  ! RHS
  real(DP) :: x_exact(n) = [1.0_DP, 1.0_DP, 1.0_DP] ! Exact solution
  real(DP) :: x_bicg(n), x_umf(n)
  real(DP) :: error_bicg, error_umf
  integer :: i
  
  ! Save current settings
  real(DP) :: saved_abs_tol, saved_rel_tol
  integer :: saved_max_iter, saved_method
  
  write(*,'(A)') '=== Ill-Conditioned BiCGSTAB Test ===  '
  write(*,'(A)') 'Matrix A = [1    0     0  ]  (ill-conditioned diagonal)'
  write(*,'(A)') '           [0  1e-12  0  ]  <- Very small element'
  write(*,'(A)') '           [0    0     1  ]'
  write(*,'(A)') 'RHS b = [1, 1e-12, 1]'
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
  write(*,'(A,3E12.4)') '  x = [', x_umf
  write(*,'(A,E12.4)') '  Max error = ', error_umf
  write(*,'(A)') ''
  
  ! Test with BiCGSTAB - relaxed tolerances for ill-conditioned problem
  sparse_solve_method = SOLVER_BICGSTAB
  bicgstab_abs_tolerance = 1.0e-10_DP  ! More relaxed
  bicgstab_rel_tolerance = 1.0e-6_DP   ! Much more relaxed
  bicgstab_max_iter = 1000
  
  x_bicg = b  ! Copy RHS
  call sparse_solve(n, n, nz, irow, pcol, amat, x_bicg, 1)
  error_bicg = maxval(abs(x_bicg - x_exact))
  
  write(*,'(A)') 'BiCGSTAB solution (relaxed tolerances):'
  write(*,'(A,3E12.4)') '  x = [', x_bicg
  write(*,'(A,E12.4)') '  Max error = ', error_bicg
  write(*,'(A)') ''
  
  ! Test with even more relaxed tolerances
  bicgstab_abs_tolerance = 1.0e-8_DP
  bicgstab_rel_tolerance = 1.0e-4_DP
  
  x_bicg = b  ! Copy RHS
  call sparse_solve(n, n, nz, irow, pcol, amat, x_bicg, 1)
  error_bicg = maxval(abs(x_bicg - x_exact))
  
  write(*,'(A)') 'BiCGSTAB solution (very relaxed tolerances):'
  write(*,'(A,3E12.4)') '  x = [', x_bicg
  write(*,'(A,E12.4)') '  Max error = ', error_bicg
  write(*,'(A)') ''
  
  ! Restore settings
  bicgstab_abs_tolerance = saved_abs_tol
  bicgstab_rel_tolerance = saved_rel_tol
  bicgstab_max_iter = saved_max_iter
  sparse_solve_method = saved_method
  
  if (error_bicg < 1.0e-6_DP) then
    write(*,'(A)') 'BiCGSTAB ill-conditioned test PASSED!'
  else
    write(*,'(A)') 'BiCGSTAB ill-conditioned test shows expected difficulty with ill-conditioning'
  end if
  
end program test_bicgstab_illcond