program test_ilu_spline_direct
  ! Direct test of ILU-preconditioned BiCGSTAB on splines
  
  use nrtype, only: I4B, DP
  use sparse_solvers_mod, only: ilu_fill_level, bicgstab_verbose, &
                                sparse_solve_method, SOLVER_BICGSTAB
  use sparse_mod, only: sparse_solve
  implicit none

  write(*,'(A)') '=== Direct ILU+BiCGSTAB Test for Splines ==='
  write(*,'(A)') ''
  
  ! Show current configuration
  write(*,'(A,I0)') 'Current ilu_fill_level = ', ilu_fill_level
  write(*,'(A,L1)') 'Current bicgstab_verbose = ', bicgstab_verbose  
  write(*,'(A,I0)') 'Current sparse_solve_method = ', sparse_solve_method
  
  ! Force BiCGSTAB method and enable verbose output
  sparse_solve_method = SOLVER_BICGSTAB
  bicgstab_verbose = .TRUE.
  
  write(*,'(A)') ''
  write(*,'(A)') 'Testing with forced BiCGSTAB method and verbose output:'
  write(*,'(A,I0)') 'Set sparse_solve_method = ', sparse_solve_method
  write(*,'(A,L1)') 'Set bicgstab_verbose = ', bicgstab_verbose
  write(*,'(A)') ''
  
  write(*,'(A)') 'If ILU is working, you should see ILU factorization messages.'
  write(*,'(A)') 'Test program ready - actual spline testing would require spline_cof integration.'
  
end program test_ilu_spline_direct