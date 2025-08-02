PROGRAM test_dual_tolerance
  USE sparse_mod
  USE sparse_types_mod, ONLY: dp
  USE sparse_conversion_mod
  USE sparse_solvers_mod, ONLY: bicgstab_abs_tolerance, bicgstab_rel_tolerance, &
                                bicgstab_tolerance, bicgstab_max_iter, bicgstab_verbose
  IMPLICIT NONE
  
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: A_full
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: b_large, b_small, x, b_orig
  REAL(kind=dp) :: max_abs_err, max_rel_err
  INTEGER :: n, i
  
  PRINT *, "=== Testing Dual Tolerance BiCGSTAB ==="
  
  ! Create a simple test matrix (diagonal dominant)
  n = 100
  ALLOCATE(A_full(n, n))
  ALLOCATE(b_large(n), b_small(n), x(n), b_orig(n))
  
  ! Create diagonal dominant matrix
  A_full = 0.0_dp
  DO i = 1, n
    A_full(i,i) = 4.0_dp
    IF (i > 1) A_full(i,i-1) = -1.0_dp
    IF (i < n) A_full(i,i+1) = -1.0_dp
  END DO
  
  ! Test 1: Large RHS - relative tolerance should dominate
  PRINT *, ""
  PRINT *, "Test 1: Large RHS (||b|| = 1000)"
  b_large = 1000.0_dp
  b_orig = b_large
  
  sparse_solve_method = 4  ! BiCGSTAB
  bicgstab_abs_tolerance = 1.0e-14_dp
  bicgstab_rel_tolerance = 1.0e-8_dp
  bicgstab_verbose = .TRUE.
  
  x = b_large
  CALL sparse_solve(A_full, x)
  CALL sparse_solver_test(A_full, x, b_orig, max_abs_err, max_rel_err)
  PRINT '(A,E12.3,A,E12.3)', "Accuracy: abs=", max_abs_err, ", rel=", max_rel_err
  PRINT '(A,E12.3)', "Expected tolerance (rel*||b||): ", bicgstab_rel_tolerance * 1000.0_dp
  
  ! Test 2: Small RHS - absolute tolerance should dominate
  PRINT *, ""
  PRINT *, "Test 2: Small RHS (||b|| = 1e-10)"
  b_small = 1.0e-10_dp
  b_orig = b_small
  
  x = b_small
  CALL sparse_solve(A_full, x)
  CALL sparse_solver_test(A_full, x, b_orig, max_abs_err, max_rel_err)
  PRINT '(A,E12.3,A,E12.3)', "Accuracy: abs=", max_abs_err, ", rel=", max_rel_err
  PRINT '(A,E12.3)', "Expected tolerance (abs): ", bicgstab_abs_tolerance
  
  ! Test 3: Legacy parameter support
  PRINT *, ""
  PRINT *, "Test 3: Legacy tolerance parameter"
  bicgstab_tolerance = 1.0e-6_dp  ! Set legacy parameter
  bicgstab_verbose = .FALSE.
  
  b_large = 100.0_dp
  b_orig = b_large
  x = b_large
  CALL sparse_solve(A_full, x)
  CALL sparse_solver_test(A_full, x, b_orig, max_abs_err, max_rel_err)
  PRINT '(A,E12.3,A,E12.3)', "Accuracy: abs=", max_abs_err, ", rel=", max_rel_err
  PRINT '(A,E12.3)', "Used legacy tolerance: ", bicgstab_tolerance
  
  DEALLOCATE(A_full, b_large, b_small, x, b_orig)
  
  PRINT *, ""
  PRINT *, "=== Dual Tolerance Tests Complete ==="
  
END PROGRAM test_dual_tolerance