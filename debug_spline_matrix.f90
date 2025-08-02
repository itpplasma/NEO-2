PROGRAM debug_spline_matrix
  USE sparse_mod
  USE sparse_types_mod, ONLY: dp
  USE sparse_conversion_mod
  USE sparse_solvers_mod, ONLY: bicgstab_tolerance, bicgstab_max_iter, bicgstab_verbose
  IMPLICIT NONE
  
  ! Test the problematic spline matrix
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: A_full
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: b, x_bicg, x_umf, b_orig
  INTEGER, DIMENSION(:), ALLOCATABLE :: irow, pcol
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: val
  REAL(kind=dp) :: max_abs_err, max_rel_err, diff_norm, cond_est
  INTEGER :: nrow, ncol, nz, i, j
  
  PRINT *, "=== Debugging Spline Matrix BiCGSTAB Failure ==="
  
  ! Create a 348x348 spline-like matrix (tridiagonal with off-diagonals)
  nrow = 348
  ncol = 348
  ALLOCATE(A_full(nrow, ncol))
  A_full = 0.0_dp
  
  ! Create a spline coefficient matrix pattern (pentadiagonal)
  DO i = 1, nrow
    ! Main diagonal
    A_full(i,i) = 4.0_dp
    
    ! Super and sub diagonals
    IF (i > 1) A_full(i,i-1) = 1.0_dp
    IF (i < nrow) A_full(i,i+1) = 1.0_dp
    
    ! Second super and sub diagonals (typical in spline problems)
    IF (i > 2) A_full(i,i-2) = -0.1_dp
    IF (i < nrow-1) A_full(i,i+2) = -0.1_dp
  END DO
  
  ! Add some ill-conditioning
  DO i = 1, MIN(10, nrow)
    A_full(i,i) = A_full(i,i) * (1.0_dp + 1.0e-8_dp * i)
  END DO
  
  ! Create RHS
  ALLOCATE(b(nrow), x_bicg(nrow), x_umf(nrow), b_orig(nrow))
  b = 1.0_dp
  b_orig = b
  
  PRINT '(A,I0,A,I0)', "Matrix size: ", nrow, " x ", ncol
  
  ! Test 1: UMFPACK solution
  sparse_solve_method = 3
  sparse_talk = .FALSE.
  x_umf = b
  PRINT *, ""
  PRINT *, "=== UMFPACK Solution ==="
  CALL sparse_solve(A_full, x_umf)
  CALL sparse_solver_test(A_full, x_umf, b_orig, max_abs_err, max_rel_err)
  PRINT '(A,E12.3,A,E12.3)', "UMFPACK accuracy: abs=", max_abs_err, ", rel=", max_rel_err
  
  ! Test 2: BiCGSTAB with default parameters
  sparse_solve_method = 4
  sparse_talk = .TRUE.
  bicgstab_verbose = .TRUE.
  bicgstab_tolerance = 1.0e-8_dp  ! Use practical default tolerance
  bicgstab_max_iter = 1000
  x_bicg = b
  PRINT *, ""
  PRINT *, "=== BiCGSTAB Solution (default parameters) ==="
  CALL sparse_solve(A_full, x_bicg)
  CALL sparse_solver_test(A_full, x_bicg, b_orig, max_abs_err, max_rel_err)
  PRINT '(A,E12.3,A,E12.3)', "BiCGSTAB accuracy: abs=", max_abs_err, ", rel=", max_rel_err
  
  ! Compare solutions
  diff_norm = SQRT(SUM((x_bicg - x_umf)**2))
  PRINT '(A,E12.3)', "Solution difference norm: ", diff_norm
  
  ! Test 3: BiCGSTAB with relaxed tolerance
  PRINT *, ""
  PRINT *, "=== BiCGSTAB with relaxed tolerance ==="
  bicgstab_tolerance = 1.0e-6_dp
  bicgstab_max_iter = 2000
  x_bicg = b
  CALL sparse_solve(A_full, x_bicg)
  CALL sparse_solver_test(A_full, x_bicg, b_orig, max_abs_err, max_rel_err)
  PRINT '(A,E12.3,A,E12.3)', "BiCGSTAB (relaxed) accuracy: abs=", max_abs_err, ", rel=", max_rel_err
  
  ! Test 4: Check matrix properties
  PRINT *, ""
  PRINT *, "=== Matrix Analysis ==="
  CALL full2sparse(A_full, irow, pcol, val, nrow, ncol, nz)
  PRINT '(A,I0)', "Non-zeros: ", nz
  PRINT '(A,F12.3)', "Sparsity: ", REAL(nz) / REAL(nrow*ncol) * 100.0_dp
  PRINT '(A,2F12.3)', "Matrix element range: ", MINVAL(val), MAXVAL(val)
  
  ! Simple condition number estimate
  cond_est = MAXVAL(ABS(val)) / MINVAL(ABS(val), MASK=ABS(val) > 1.0e-15_dp)
  PRINT '(A,E12.3)', "Element condition estimate: ", cond_est
  
END PROGRAM debug_spline_matrix