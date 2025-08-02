PROGRAM test_adaptive
  USE sparse_mod
  USE sparse_types_mod, ONLY: dp
  IMPLICIT NONE
  
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: A_full
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: b, x, b_orig
  REAL(kind=dp) :: max_abs_err, max_rel_err
  INTEGER :: i
  
  ! Test adaptive selection mode
  sparse_solve_method = 0  ! Auto-select
  sparse_talk = .TRUE.
  
  PRINT *, "=== Testing Adaptive Solver Selection ==="
  
  ! Test 1: Small matrix (should use UMFPACK)
  CALL load_mini_example(A_full)
  ALLOCATE(b(5), x(5), b_orig(5))
  b = 1.0_dp
  b_orig = b
  x = b
  
  PRINT *, ""
  PRINT *, "Test 1: Small 5x5 matrix (should auto-select UMFPACK)"
  CALL sparse_solve(A_full, x)
  CALL sparse_solver_test(A_full, x, b_orig, max_abs_err, max_rel_err)
  PRINT '(A,E12.3)', "Max error:", max_abs_err
  
  DEALLOCATE(A_full, b, x, b_orig)
  
  ! Test 2: Create a larger matrix (should use BiCGSTAB)
  ALLOCATE(A_full(150,150), b(150), x(150), b_orig(150))
  A_full = 0.0_dp
  ! Create a simple diagonal matrix
  DO i = 1, 150
    A_full(i,i) = REAL(i, dp)
  END DO
  b = 1.0_dp
  b_orig = b
  x = b
  
  PRINT *, ""
  PRINT *, "Test 2: Large 150x150 matrix (should auto-select BiCGSTAB)"
  CALL sparse_solve(A_full, x)
  CALL sparse_solver_test(A_full, x, b_orig, max_abs_err, max_rel_err)
  PRINT '(A,E12.3)', "Max error:", max_abs_err
  
  PRINT *, ""
  PRINT *, "=== Adaptive Selection Test Complete ==="
  
END PROGRAM test_adaptive