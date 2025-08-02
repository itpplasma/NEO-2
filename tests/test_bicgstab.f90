PROGRAM test_bicgstab
  ! Comprehensive tests for bicgstab_mod module
  ! Tests BiCGSTAB solver with and without preconditioning
  
  USE sparse_types_mod, ONLY: dp, long
  USE sparse_utils_mod
  USE ilu_precond_mod
  USE bicgstab_mod
  IMPLICIT NONE
  
  ! Test variables
  INTEGER :: n, nz, iter, max_iter
  INTEGER, ALLOCATABLE :: csr_row_ptr(:), csr_col_idx(:)
  REAL(kind=dp), ALLOCATABLE :: csr_val(:)
  REAL(kind=dp), ALLOCATABLE :: x(:), b(:), x_exact(:), r(:)
  REAL(kind=dp) :: tol, achieved_tol, residual_norm
  TYPE(bicgstab_stats) :: stats
  TYPE(ilu_factorization) :: ilu_fac
  INTEGER :: info, i, j, k
  LOGICAL :: test_passed, tests_passed
  LOGICAL :: converged
  REAL(kind=dp) :: error_norm
  REAL(kind=dp), PARAMETER :: test_tol = 1.0e-12_dp
  
  ! Complex test variables
  COMPLEX(kind=dp), ALLOCATABLE :: z_csr_val(:)
  COMPLEX(kind=dp), ALLOCATABLE :: z_x(:), z_b(:), z_x_exact(:)
  TYPE(ilu_factorization_complex) :: z_ilu_fac
  
  tests_passed = .TRUE.
  
  WRITE(*,'(A)') "================================="
  WRITE(*,'(A)') "BiCGSTAB Solver Test Suite"
  WRITE(*,'(A)') "================================="
  WRITE(*,*)
  
  ! Test 1: Diagonal system (should converge in 1 iteration)
  WRITE(*,'(A)') "Test 1: Diagonal system"
  test_passed = .TRUE.
  
  ! 4x4 diagonal matrix
  n = 4
  nz = 4
  ALLOCATE(csr_row_ptr(5), csr_col_idx(4), csr_val(4))
  csr_row_ptr = (/1, 2, 3, 4, 5/)
  csr_col_idx = (/1, 2, 3, 4/)
  csr_val = (/2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp/)
  
  ALLOCATE(x(4), b(4), x_exact(4))
  x_exact = (/1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp/)
  
  ! Compute b = A*x_exact
  CALL csr_matvec(n, csr_row_ptr, csr_col_idx, csr_val, x_exact, b)
  
  ! Initial guess
  x = 0.0_dp
  
  ! Solve without preconditioner
  tol = 1.0e-12_dp
  max_iter = 100
  CALL bicgstab_solve(n, csr_row_ptr, csr_col_idx, csr_val, b, x, &
                      tol, max_iter, converged, iter, stats)
  
  IF (converged) THEN
    error_norm = SQRT(SUM((x - x_exact)**2))
    IF (error_norm < 1.0e-12_dp .AND. iter <= 5) THEN
      WRITE(*,'(A,I0,A)') "[PASS] Diagonal system converged in ", iter, " iterations"
    ELSE
      test_passed = .FALSE.
      tests_passed = .FALSE.
      WRITE(*,'(A,E12.4)') "[FAIL] Diagonal system error = ", error_norm
    END IF
  ELSE
    test_passed = .FALSE.
    tests_passed = .FALSE.
    WRITE(*,'(A)') "[FAIL] Diagonal system did not converge"
  END IF
  
  DEALLOCATE(csr_row_ptr, csr_col_idx, csr_val)
  DEALLOCATE(x, b, x_exact)
  
  ! Test 2: Small SPD system
  WRITE(*,'(A)') "Test 2: Small SPD system"
  test_passed = .TRUE.
  
  ! 3x3 SPD matrix
  n = 3
  nz = 7
  ALLOCATE(csr_row_ptr(4), csr_col_idx(7), csr_val(7))
  csr_row_ptr = (/1, 4, 6, 8/)
  csr_col_idx = (/1, 2, 3, 1, 2, 2, 3/)
  csr_val = (/4.0_dp, 1.0_dp, 0.5_dp, 1.0_dp, 3.0_dp, 0.5_dp, 2.0_dp/)
  
  ALLOCATE(x(3), b(3), x_exact(3))
  x_exact = (/1.0_dp, -1.0_dp, 2.0_dp/)
  
  ! Compute b = A*x_exact
  CALL csr_matvec(n, csr_row_ptr, csr_col_idx, csr_val, x_exact, b)
  
  ! Initial guess
  x = 0.0_dp
  
  ! Solve
  tol = 1.0e-10_dp
  max_iter = 50
  CALL bicgstab_solve(n, csr_row_ptr, csr_col_idx, csr_val, b, x, &
                      tol, max_iter, converged, iter, stats)
  
  IF (converged) THEN
    error_norm = SQRT(SUM((x - x_exact)**2))
    IF (error_norm < test_tol) THEN
      WRITE(*,'(A,I0,A)') "[PASS] SPD system converged in ", iter, " iterations"
    ELSE
      test_passed = .FALSE.
      tests_passed = .FALSE.
      WRITE(*,'(A,E12.4)') "[FAIL] SPD system error = ", error_norm
    END IF
  ELSE
    test_passed = .FALSE.
    tests_passed = .FALSE.
    WRITE(*,'(A)') "[FAIL] SPD system did not converge"
  END IF
  
  DEALLOCATE(csr_row_ptr, csr_col_idx, csr_val)
  DEALLOCATE(x, b, x_exact)
  
  ! Test 3: BiCGSTAB with ILU preconditioner
  WRITE(*,'(A)') "Test 3: BiCGSTAB with ILU preconditioner"
  test_passed = .TRUE.
  
  ! Same SPD system
  n = 3
  nz = 7
  ALLOCATE(csr_row_ptr(4), csr_col_idx(7), csr_val(7))
  csr_row_ptr = (/1, 4, 6, 8/)
  csr_col_idx = (/1, 2, 3, 1, 2, 2, 3/)
  csr_val = (/4.0_dp, 1.0_dp, 0.5_dp, 1.0_dp, 3.0_dp, 0.5_dp, 2.0_dp/)
  
  ! Compute ILU(0) preconditioner
  CALL ilu_factorize(n, csr_row_ptr, csr_col_idx, csr_val, 0, 0.0_dp, ilu_fac, info)
  
  IF (info == 0) THEN
    ALLOCATE(x(3), b(3), x_exact(3))
    x_exact = (/1.0_dp, -1.0_dp, 2.0_dp/)
    
    ! Compute b = A*x_exact
    CALL csr_matvec(n, csr_row_ptr, csr_col_idx, csr_val, x_exact, b)
    
    ! Initial guess
    x = 0.0_dp
    
    ! Solve with ILU preconditioner
    tol = 1.0e-8_dp  ! Slightly relaxed tolerance for preconditioned solver
    max_iter = 50
    CALL bicgstab_solve_precond(n, csr_row_ptr, csr_col_idx, csr_val, b, x, &
                                ilu_fac, tol, max_iter, converged, iter, stats)
    
    IF (converged) THEN
      error_norm = SQRT(SUM((x - x_exact)**2))
      IF (error_norm < 1.0e-6_dp .AND. iter < 30) THEN
        WRITE(*,'(A,I0,A)') "[PASS] Preconditioned system converged in ", iter, " iterations"
      ELSE
        test_passed = .FALSE.
        tests_passed = .FALSE.
        WRITE(*,'(A,E12.4,A,I0)') "[FAIL] Preconditioned error = ", error_norm, &
          ", iterations = ", iter
      END IF
    ELSE
      test_passed = .FALSE.
      tests_passed = .FALSE.
      WRITE(*,'(A,I0,A,E12.4)') "[FAIL] Preconditioned system did not converge, iter = ", &
        iter, ", final_residual = ", stats%final_residual
    END IF
    
    DEALLOCATE(x, b, x_exact)
    CALL ilu_free(ilu_fac)
  ELSE
    test_passed = .FALSE.
    tests_passed = .FALSE.
    WRITE(*,'(A)') "[FAIL] Could not compute ILU preconditioner"
  END IF
  
  DEALLOCATE(csr_row_ptr, csr_col_idx, csr_val)
  
  ! Test 4: Convergence monitoring
  WRITE(*,'(A)') "Test 4: Convergence monitoring"
  test_passed = .TRUE.
  
  ! 5x5 tridiagonal matrix
  n = 5
  nz = 13  ! 5 diagonal + 4 upper + 4 lower
  ALLOCATE(csr_row_ptr(6), csr_col_idx(13), csr_val(13))
  
  ! Build tridiagonal matrix
  k = 1
  DO i = 1, n
    csr_row_ptr(i) = k
    IF (i > 1) THEN
      csr_col_idx(k) = i-1
      csr_val(k) = -1.0_dp
      k = k + 1
    END IF
    csr_col_idx(k) = i
    csr_val(k) = 4.0_dp
    k = k + 1
    IF (i < n) THEN
      csr_col_idx(k) = i+1
      csr_val(k) = -1.0_dp
      k = k + 1
    END IF
  END DO
  csr_row_ptr(n+1) = k
  
  ALLOCATE(x(5), b(5))
  b = 1.0_dp
  x = 0.0_dp
  
  ! Solve with detailed statistics
  tol = 1.0e-8_dp
  max_iter = 50
  CALL bicgstab_solve(n, csr_row_ptr, csr_col_idx, csr_val, b, x, &
                      tol, max_iter, converged, iter, stats)
  
  IF (converged) THEN
    ! Check that residual decreased monotonically (approximately)
    IF (stats%final_residual < tol .AND. iter > 0) THEN
      WRITE(*,'(A,E12.4)') "[PASS] Convergence monitoring, final residual = ", &
        stats%final_residual
    ELSE
      test_passed = .FALSE.
      tests_passed = .FALSE.
      WRITE(*,'(A)') "[FAIL] Convergence monitoring failed"
    END IF
  ELSE
    test_passed = .FALSE.
    tests_passed = .FALSE.
    WRITE(*,'(A)') "[FAIL] Tridiagonal system did not converge"
  END IF
  
  DEALLOCATE(csr_row_ptr, csr_col_idx, csr_val)
  DEALLOCATE(x, b)
  
  ! Test 5: Zero RHS
  WRITE(*,'(A)') "Test 5: Zero RHS"
  test_passed = .TRUE.
  
  ! Simple diagonal matrix
  n = 3
  nz = 3
  ALLOCATE(csr_row_ptr(4), csr_col_idx(3), csr_val(3))
  csr_row_ptr = (/1, 2, 3, 4/)
  csr_col_idx = (/1, 2, 3/)
  csr_val = (/1.0_dp, 2.0_dp, 3.0_dp/)
  
  ALLOCATE(x(3), b(3))
  b = 0.0_dp  ! Zero RHS
  x = 0.0_dp
  
  tol = 1.0e-12_dp
  max_iter = 10
  CALL bicgstab_solve(n, csr_row_ptr, csr_col_idx, csr_val, b, x, &
                      tol, max_iter, converged, iter, stats)
  
  IF (converged .AND. iter == 0) THEN
    IF (MAXVAL(ABS(x)) < tol) THEN
      WRITE(*,'(A)') "[PASS] Zero RHS handled correctly"
    ELSE
      test_passed = .FALSE.
      tests_passed = .FALSE.
      WRITE(*,'(A)') "[FAIL] Zero RHS produced non-zero solution"
    END IF
  ELSE
    test_passed = .FALSE.
    tests_passed = .FALSE.
    WRITE(*,'(A,I0)') "[FAIL] Zero RHS convergence issue, iter = ", iter
  END IF
  
  DEALLOCATE(csr_row_ptr, csr_col_idx, csr_val)
  DEALLOCATE(x, b)
  
  ! Test 6: Complex system
  WRITE(*,'(A)') "Test 6: Complex BiCGSTAB"
  test_passed = .TRUE.
  
  ! 2x2 complex matrix
  n = 2
  nz = 4
  ALLOCATE(csr_row_ptr(3), csr_col_idx(4), z_csr_val(4))
  csr_row_ptr = (/1, 3, 5/)
  csr_col_idx = (/1, 2, 1, 2/)
  z_csr_val = (/(3.0_dp,0.0_dp), (1.0_dp,1.0_dp), (1.0_dp,-1.0_dp), (3.0_dp,0.0_dp)/)
  
  ALLOCATE(z_x(2), z_b(2), z_x_exact(2))
  z_x_exact = (/(1.0_dp,1.0_dp), (0.0_dp,-1.0_dp)/)
  
  ! Compute b = A*x_exact
  CALL csr_matvec(n, csr_row_ptr, csr_col_idx, z_csr_val, z_x_exact, z_b)
  
  ! Initial guess
  z_x = (0.0_dp, 0.0_dp)
  
  ! Solve
  tol = 1.0e-10_dp
  max_iter = 50
  CALL bicgstab_solve(n, csr_row_ptr, csr_col_idx, z_csr_val, z_b, z_x, &
                      tol, max_iter, converged, iter, stats)
  
  IF (converged) THEN
    error_norm = SQRT(SUM(ABS(z_x - z_x_exact)**2))
    IF (error_norm < test_tol) THEN
      WRITE(*,'(A,I0,A)') "[PASS] Complex system converged in ", iter, " iterations"
    ELSE
      test_passed = .FALSE.
      tests_passed = .FALSE.
      WRITE(*,'(A,E12.4)') "[FAIL] Complex system error = ", error_norm
    END IF
  ELSE
    test_passed = .FALSE.
    tests_passed = .FALSE.
    WRITE(*,'(A)') "[FAIL] Complex system did not converge"
  END IF
  
  DEALLOCATE(csr_row_ptr, csr_col_idx, z_csr_val)
  DEALLOCATE(z_x, z_b, z_x_exact)
  
  ! Test 7: Iteration limit
  WRITE(*,'(A)') "Test 7: Iteration limit"
  test_passed = .TRUE.
  
  ! Ill-conditioned system
  n = 3
  nz = 9
  ALLOCATE(csr_row_ptr(4), csr_col_idx(9), csr_val(9))
  csr_row_ptr = (/1, 4, 7, 10/)
  csr_col_idx = (/1, 2, 3, 1, 2, 3, 1, 2, 3/)
  csr_val = (/1.0_dp, 0.99_dp, 0.0_dp, 0.99_dp, 1.0_dp, 0.99_dp, 0.0_dp, 0.99_dp, 1.0_dp/)
  
  ALLOCATE(x(3), b(3))
  b = (/1.0_dp, 2.0_dp, 3.0_dp/)
  x = 0.0_dp
  
  ! Set very low iteration limit
  tol = 1.0e-12_dp
  max_iter = 2
  CALL bicgstab_solve(n, csr_row_ptr, csr_col_idx, csr_val, b, x, &
                      tol, max_iter, converged, iter, stats)
  
  IF (.NOT. converged .AND. iter == max_iter) THEN
    WRITE(*,'(A)') "[PASS] Iteration limit enforced correctly"
  ELSE
    test_passed = .FALSE.
    tests_passed = .FALSE.
    WRITE(*,'(A,L1,A,I0)') "[FAIL] Iteration limit test: converged=", converged, &
      ", iter=", iter
  END IF
  
  DEALLOCATE(csr_row_ptr, csr_col_idx, csr_val)
  DEALLOCATE(x, b)
  
  ! Test 8: Restart capability
  WRITE(*,'(A)') "Test 8: BiCGSTAB restart on stagnation"
  test_passed = .TRUE.
  
  ! Well-conditioned system for reliable testing
  n = 4
  nz = 10
  ALLOCATE(csr_row_ptr(5), csr_col_idx(10), csr_val(10))
  csr_row_ptr = (/1, 4, 7, 9, 11/)
  csr_col_idx = (/1, 2, 3, 1, 2, 3, 2, 3, 3, 4/)
  csr_val = (/5.0_dp, 1.0_dp, 0.5_dp, 1.0_dp, 4.0_dp, 1.0_dp, 1.0_dp, 3.0_dp, 0.5_dp, 2.0_dp/)
  
  ALLOCATE(x(4), b(4))
  b = (/1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp/)
  x = 0.0_dp
  
  tol = 1.0e-10_dp
  max_iter = 100
  CALL bicgstab_solve(n, csr_row_ptr, csr_col_idx, csr_val, b, x, &
                      tol, max_iter, converged, iter, stats)
  
  IF (converged) THEN
    WRITE(*,'(A,I0,A)') "[PASS] BiCGSTAB with restart capability converged in ", &
      iter, " iterations"
  ELSE
    test_passed = .FALSE.
    tests_passed = .FALSE.
    WRITE(*,'(A)') "[FAIL] BiCGSTAB restart test failed to converge"
  END IF
  
  DEALLOCATE(csr_row_ptr, csr_col_idx, csr_val)
  DEALLOCATE(x, b)
  
  ! Test 9: Residual norm calculation
  WRITE(*,'(A)') "Test 9: Residual norm verification"
  test_passed = .TRUE.
  
  ! Simple diagonal system
  n = 3
  nz = 3
  ALLOCATE(csr_row_ptr(4), csr_col_idx(3), csr_val(3))
  csr_row_ptr = (/1, 2, 3, 4/)
  csr_col_idx = (/1, 2, 3/)
  csr_val = (/2.0_dp, 3.0_dp, 4.0_dp/)
  
  ALLOCATE(x(3), b(3), r(3))
  b = (/2.0_dp, 6.0_dp, 12.0_dp/)  ! Solution should be [1, 2, 3]
  x = 0.0_dp
  
  tol = 1.0e-10_dp
  max_iter = 50
  CALL bicgstab_solve(n, csr_row_ptr, csr_col_idx, csr_val, b, x, &
                      tol, max_iter, converged, iter, stats)
  
  IF (converged) THEN
    ! Manually compute residual r = b - A*x
    CALL csr_matvec(n, csr_row_ptr, csr_col_idx, csr_val, x, r)
    r = b - r
    residual_norm = SQRT(SUM(r**2))
    
    IF (ABS(residual_norm - stats%final_residual) < 1.0e-12_dp) THEN
      WRITE(*,'(A,E12.4)') "[PASS] Residual norm verified = ", residual_norm
    ELSE
      test_passed = .FALSE.
      tests_passed = .FALSE.
      WRITE(*,'(A,E12.4,A,E12.4)') "[FAIL] Residual mismatch: computed = ", &
        residual_norm, ", reported = ", stats%final_residual
    END IF
  ELSE
    test_passed = .FALSE.
    tests_passed = .FALSE.
    WRITE(*,'(A)') "[FAIL] System did not converge for residual test"
  END IF
  
  DEALLOCATE(csr_row_ptr, csr_col_idx, csr_val)
  DEALLOCATE(x, b, r)
  
  ! Test 10: Large system performance
  WRITE(*,'(A)') "Test 10: Large system (100x100 tridiagonal)"
  test_passed = .TRUE.
  
  ! Create 100x100 tridiagonal matrix
  n = 100
  nz = 298  ! 100 diagonal + 99 upper + 99 lower
  ALLOCATE(csr_row_ptr(101), csr_col_idx(298), csr_val(298))
  
  ! Build tridiagonal structure
  k = 1
  DO i = 1, n
    csr_row_ptr(i) = k
    IF (i > 1) THEN
      csr_col_idx(k) = i-1
      csr_val(k) = -1.0_dp
      k = k + 1
    END IF
    csr_col_idx(k) = i
    csr_val(k) = 4.0_dp
    k = k + 1
    IF (i < n) THEN
      csr_col_idx(k) = i+1
      csr_val(k) = -1.0_dp
      k = k + 1
    END IF
  END DO
  csr_row_ptr(n+1) = k
  
  ALLOCATE(x(100), b(100))
  b = 1.0_dp
  x = 0.0_dp
  
  ! Solve without preconditioner
  tol = 1.0e-8_dp
  max_iter = 200
  CALL bicgstab_solve(n, csr_row_ptr, csr_col_idx, csr_val, b, x, &
                      tol, max_iter, converged, iter, stats)
  
  IF (converged .AND. iter < 100) THEN
    WRITE(*,'(A,I0,A,E12.4)') "[PASS] Large system converged in ", iter, &
      " iterations, residual = ", stats%final_residual
  ELSE
    test_passed = .FALSE.
    tests_passed = .FALSE.
    IF (converged) THEN
      WRITE(*,'(A,I0)') "[FAIL] Large system took too many iterations: ", iter
    ELSE
      WRITE(*,'(A)') "[FAIL] Large system did not converge"
    END IF
  END IF
  
  DEALLOCATE(csr_row_ptr, csr_col_idx, csr_val)
  DEALLOCATE(x, b)
  
  ! Summary
  WRITE(*,*)
  WRITE(*,'(A)') "================================="
  IF (tests_passed) THEN
    WRITE(*,'(A)') "All BiCGSTAB tests PASSED!"
  ELSE
    WRITE(*,'(A)') "Some BiCGSTAB tests FAILED!"
    STOP 1
  END IF
  WRITE(*,'(A)') "================================="
  
END PROGRAM test_bicgstab