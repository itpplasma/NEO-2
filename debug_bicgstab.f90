PROGRAM debug_bicgstab
  USE sparse_mod
  USE sparse_types_mod, ONLY: dp
  IMPLICIT NONE
  
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: A_full
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: b, x_bicgstab, x_umfpack, b_orig
  REAL(kind=dp) :: max_abs_err, max_rel_err, diff_norm
  INTEGER :: i
  
  ! Load the problematic mini example
  CALL load_mini_example(A_full)
  
  ! Create RHS (all ones)
  ALLOCATE(b(5), x_bicgstab(5), x_umfpack(5), b_orig(5))
  b = 1.0_dp
  b_orig = b
  
  PRINT *, "=== Debug BiCGSTAB vs UMFPACK ==="
  PRINT *, "Matrix A:"
  DO i = 1, 5
    PRINT '(5F12.1)', A_full(i,:)
  END DO
  PRINT *, ""
  
  ! Test with UMFPACK (method 3)
  sparse_solve_method = 3
  sparse_talk = .TRUE.
  x_umfpack = b
  PRINT *, "Solving with UMFPACK (method 3)..."
  CALL sparse_solve(A_full, x_umfpack)
  CALL sparse_solver_test(A_full, x_umfpack, b_orig, max_abs_err, max_rel_err)
  PRINT '(A,E12.3,A,E12.3)', "UMFPACK - Max abs error:", max_abs_err, ", Max rel error:", max_rel_err
  PRINT '(A,5F12.6)', "UMFPACK solution:", x_umfpack
  PRINT *, ""
  
  ! Test with BiCGSTAB (method 4)
  sparse_solve_method = 4
  x_bicgstab = b
  PRINT *, "Solving with BiCGSTAB (method 4)..."
  CALL sparse_solve(A_full, x_bicgstab)
  CALL sparse_solver_test(A_full, x_bicgstab, b_orig, max_abs_err, max_rel_err)
  PRINT '(A,E12.3,A,E12.3)', "BiCGSTAB - Max abs error:", max_abs_err, ", Max rel error:", max_rel_err
  PRINT '(A,5F12.6)', "BiCGSTAB solution:", x_bicgstab
  PRINT *, ""
  
  ! Compare solutions
  diff_norm = SQRT(SUM((x_bicgstab - x_umfpack)**2))
  PRINT '(A,E12.3)', "Solution difference norm:", diff_norm
  PRINT '(A,5E12.3)', "Difference:", x_bicgstab - x_umfpack
  
  ! Check matrix condition
  PRINT *, ""
  PRINT '(A,F12.1)', "Matrix element range: ", MINVAL(A_full), " to ", MAXVAL(A_full)
  
END PROGRAM debug_bicgstab