PROGRAM debug_bicgstab_small
  USE sparse_mod
  USE sparse_types_mod, ONLY: dp
  USE sparse_conversion_mod
  USE sparse_solvers_mod, ONLY: bicgstab_tolerance, bicgstab_max_iter, bicgstab_verbose
  IMPLICIT NONE
  
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: A_full
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: b, x_bicg, x_umf, b_orig
  REAL(kind=dp) :: max_abs_err, max_rel_err, diff_norm
  INTEGER :: nrow, ncol, i, j, test_size
  
  PRINT *, "=== Systematic BiCGSTAB Debugging ==="
  
  ! Test progressively larger spline-like matrices
  DO test_size = 5, 50, 5
    PRINT *, ""
    PRINT '(A,I0,A,I0,A)', "=== Testing ", test_size, "x", test_size, " spline matrix ==="
    
    nrow = test_size
    ncol = test_size
    ALLOCATE(A_full(nrow, ncol))
    A_full = 0.0_dp
    
    ! Create spline coefficient matrix pattern (tridiagonal + boundary conditions)
    DO i = 1, nrow
      A_full(i,i) = 4.0_dp  ! Main diagonal
      
      IF (i > 1) A_full(i,i-1) = 1.0_dp      ! Sub-diagonal
      IF (i < nrow) A_full(i,i+1) = 1.0_dp   ! Super-diagonal
    END DO
    
    ! Add boundary condition modifications (typical in spline problems)
    A_full(1,1) = 2.0_dp      ! Natural boundary at start
    A_full(nrow,nrow) = 2.0_dp ! Natural boundary at end
    
    ! Create RHS
    ALLOCATE(b(nrow), x_bicg(nrow), x_umf(nrow), b_orig(nrow))
    b = 1.0_dp
    b_orig = b
    
    ! Test UMFPACK
    sparse_solve_method = 3
    sparse_talk = .FALSE.
    x_umf = b
    CALL sparse_solve(A_full, x_umf)
    CALL sparse_solver_test(A_full, x_umf, b_orig, max_abs_err, max_rel_err)
    PRINT '(A,E10.2)', "UMFPACK accuracy: ", max_abs_err
    
    ! Test BiCGSTAB with verbose output
    sparse_solve_method = 4
    sparse_talk = .FALSE.
    bicgstab_verbose = .TRUE.
    bicgstab_tolerance = 1.0e-8_dp  ! Use practical tolerance from NEO-2
    bicgstab_max_iter = 1000
    x_bicg = b
    CALL sparse_solve(A_full, x_bicg)
    CALL sparse_solver_test(A_full, x_bicg, b_orig, max_abs_err, max_rel_err)
    PRINT '(A,E10.2)', "BiCGSTAB accuracy: ", max_abs_err
    
    ! Compare solutions
    diff_norm = SQRT(SUM((x_bicg - x_umf)**2))
    PRINT '(A,E10.2)', "Solution difference: ", diff_norm
    
    ! Check if BiCGSTAB failed catastrophically
    IF (max_abs_err > 1.0e-6_dp) THEN
      PRINT *, "*** BiCGSTAB FAILED at size ", test_size, " ***"
      PRINT *, "Matrix condition analysis:"
      PRINT '(A,2F10.3)', "Matrix element range: ", MINVAL(A_full), " to ", MAXVAL(A_full)
      
      ! Exit at first failure to analyze
      EXIT
    END IF
    
    DEALLOCATE(A_full, b, x_bicg, x_umf, b_orig)
  END DO
  
  PRINT *, ""
  PRINT *, "=== Testing even simpler matrices ==="
  
  ! Test 1: Simple 3x3 diagonal matrix
  PRINT *, ""
  PRINT *, "Test 1: 3x3 diagonal matrix"
  ALLOCATE(A_full(3,3), b(3), x_bicg(3), x_umf(3), b_orig(3))
  A_full = 0.0_dp
  A_full(1,1) = 1.0_dp
  A_full(2,2) = 2.0_dp  
  A_full(3,3) = 3.0_dp
  b = [1.0_dp, 2.0_dp, 3.0_dp]
  b_orig = b
  
  ! Expected solution: [1, 1, 1]
  sparse_solve_method = 4
  bicgstab_verbose = .TRUE.
  x_bicg = b
  CALL sparse_solve(A_full, x_bicg)
  PRINT '(A,3F8.4)', "BiCGSTAB solution: ", x_bicg
  PRINT '(A,3F8.4)', "Expected solution: ", [1.0_dp, 1.0_dp, 1.0_dp]
  
  DEALLOCATE(A_full, b, x_bicg, x_umf, b_orig)
  
  ! Test 2: Simple 3x3 tridiagonal matrix
  PRINT *, ""
  PRINT *, "Test 2: 3x3 tridiagonal matrix"
  ALLOCATE(A_full(3,3), b(3), x_bicg(3), x_umf(3), b_orig(3))
  A_full = 0.0_dp
  A_full(1,1) = 2.0_dp; A_full(1,2) = -1.0_dp
  A_full(2,1) = -1.0_dp; A_full(2,2) = 2.0_dp; A_full(2,3) = -1.0_dp
  A_full(3,2) = -1.0_dp; A_full(3,3) = 2.0_dp
  b = [1.0_dp, 0.0_dp, 1.0_dp]
  b_orig = b
  
  sparse_solve_method = 4
  bicgstab_verbose = .TRUE.
  x_bicg = b
  CALL sparse_solve(A_full, x_bicg)
  CALL sparse_solver_test(A_full, x_bicg, b_orig, max_abs_err, max_rel_err)
  PRINT '(A,E10.2)', "BiCGSTAB accuracy: ", max_abs_err
  
END PROGRAM debug_bicgstab_small