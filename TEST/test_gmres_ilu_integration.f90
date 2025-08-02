PROGRAM test_gmres_ilu_integration
  ! TDD RED phase: Test GMRES with ILU preconditioning integration
  ! Verify that GMRES can use existing ILU preconditioner effectively
  
  USE nrtype, ONLY: I4B, DP
  USE gmres_mod, ONLY: gmres_workspace, create_gmres_workspace, destroy_gmres_workspace, &
                       gmres_solve_structured_preconditioned
  USE ilu_precond_mod, ONLY: ilu_factorization, ilu_factorize, ilu_free
  USE sparse_utils_mod, ONLY: csc_to_csr_real
  IMPLICIT NONE
  
  ! Test parameters
  INTEGER(I4B), PARAMETER :: n = 50         ! Matrix size
  INTEGER(I4B), PARAMETER :: restart = 10   ! GMRES restart parameter
  REAL(DP), PARAMETER :: tol = 1.0e-10_DP   ! Convergence tolerance
  
  TYPE(gmres_workspace) :: workspace
  TYPE(ilu_factorization) :: ilu_fac
  
  ! Sparse matrix in CSC format (for compatibility)
  INTEGER(I4B), ALLOCATABLE :: irow(:), pcol(:)   ! CSC indices
  REAL(DP), ALLOCATABLE :: amat(:)                ! CSC values
  INTEGER(I4B) :: nnz                             ! Number of nonzeros
  
  ! CSR format for ILU
  INTEGER(I4B), ALLOCATABLE :: row_ptr(:), col_idx(:)
  REAL(DP), ALLOCATABLE :: val_csr(:)
  
  ! Vectors
  REAL(DP), ALLOCATABLE :: b(:), x(:), x_exact(:), x_initial(:)
  REAL(DP) :: residual_norm, error_norm
  INTEGER(I4B) :: max_iter, iterations
  LOGICAL :: converged, test_passed
  INTEGER(I4B) :: i, j, k, info
  INTEGER(I4B) :: total_tests, passed_tests
  
  WRITE(*,'(A)') '============================================='
  WRITE(*,'(A)') 'GMRES ILU Preconditioning Test Suite (TDD RED)'
  WRITE(*,'(A)') 'Testing GMRES with ILU integration before implementation'
  WRITE(*,'(A)') '============================================='
  WRITE(*,*)
  
  total_tests = 0
  passed_tests = 0
  
  ! Setup workspace
  CALL create_gmres_workspace(workspace, n, restart)
  
  ! Test 1: GMRES with ILU(0) on tridiagonal system
  WRITE(*,'(A)', ADVANCE='NO') 'Test 1: GMRES + ILU(0) on tridiagonal system... '
  total_tests = total_tests + 1
  
  ! Create tridiagonal system in CSC format
  nnz = 3*n - 2
  ALLOCATE(irow(nnz), pcol(n+1), amat(nnz))
  ALLOCATE(row_ptr(n+1), col_idx(nnz), val_csr(nnz))
  ALLOCATE(b(n), x(n), x_exact(n), x_initial(n))
  
  ! Build tridiagonal matrix: -1, 2, -1
  k = 0
  pcol(1) = 1
  DO j = 1, n
    IF (j > 1) THEN
      k = k + 1
      irow(k) = j - 1
      amat(k) = -1.0_DP
    END IF
    k = k + 1
    irow(k) = j
    amat(k) = 2.0_DP
    IF (j < n) THEN
      k = k + 1
      irow(k) = j + 1
      amat(k) = -1.0_DP
    END IF
    pcol(j+1) = k + 1
  END DO
  
  ! Convert to CSR for ILU
  CALL csc_to_csr_real(n, n, nnz, pcol, irow, amat, row_ptr, col_idx, val_csr)
  
  ! Setup ILU(0) factorization
  CALL ilu_factorize(n, row_ptr, col_idx, val_csr, 0, 0.0_DP, ilu_fac, info)
  
  ! Setup test problem
  x_exact = 1.0_DP
  b = 0.0_DP
  DO i = 1, n
    DO j = MAX(1, i-1), MIN(n, i+1)
      IF (i == j) THEN
        b(i) = b(i) + 2.0_DP * x_exact(j)
      ELSE
        b(i) = b(i) - 1.0_DP * x_exact(j)
      END IF
    END DO
  END DO
  x_initial = 0.0_DP
  
  ! Solve with preconditioned GMRES
  max_iter = n
  CALL gmres_solve_structured_preconditioned(workspace, n, row_ptr, col_idx, val_csr, &
                                             b, x_initial, max_iter, tol, ilu_fac, &
                                             x, iterations, residual_norm, converged, info)
  
  ! Check convergence and accuracy
  error_norm = SQRT(SUM((x - x_exact)**2))
  test_passed = converged .AND. error_norm < tol * 100.0_DP
  
  IF (test_passed) THEN
    WRITE(*,'(A)') '[PASS]'
    passed_tests = passed_tests + 1
  ELSE
    WRITE(*,'(A)') '[FAIL]'
  END IF
  
  ! Clean up ILU
  CALL ilu_free(ilu_fac)
  
  ! Test 2: Compare convergence with and without preconditioning
  WRITE(*,'(A)', ADVANCE='NO') 'Test 2: Convergence improvement with ILU preconditioning... '
  total_tests = total_tests + 1
  
  ! Solve without preconditioning first
  CALL gmres_solve_structured_preconditioned(workspace, n, row_ptr, col_idx, val_csr, &
                                             b, x_initial, 20, tol, ilu_fac, &
                                             x, iterations, residual_norm, converged, info, &
                                             use_preconditioner=.FALSE.)
  k = iterations  ! Store unpreconditioned iteration count
  
  ! Setup ILU(1) for better preconditioning
  CALL ilu_factorize(n, row_ptr, col_idx, val_csr, 1, 0.0_DP, ilu_fac, info)
  
  ! Solve with preconditioning
  CALL gmres_solve_structured_preconditioned(workspace, n, row_ptr, col_idx, val_csr, &
                                             b, x_initial, 20, tol, ilu_fac, &
                                             x, iterations, residual_norm, converged, info, &
                                             use_preconditioner=.TRUE.)
  
  ! Preconditioning should reduce iterations significantly
  test_passed = iterations < k / 2  ! At least 2x speedup expected
  
  IF (test_passed) THEN
    WRITE(*,'(A)') '[PASS]'
    WRITE(*,'(A,I0,A,I0)') '    Without ILU: ', k, ' iterations, With ILU(1): ', iterations
    passed_tests = passed_tests + 1
  ELSE
    WRITE(*,'(A)') '[FAIL]'
  END IF
  
  CALL ilu_free(ilu_fac)
  
  ! Test 3: ILU(k) with different fill levels
  WRITE(*,'(A)', ADVANCE='NO') 'Test 3: GMRES with different ILU fill levels... '
  total_tests = total_tests + 1
  
  test_passed = .TRUE.
  
  ! Test ILU(0), ILU(1), ILU(2)
  DO k = 0, 2
    CALL ilu_factorize(n, row_ptr, col_idx, val_csr, k, 0.0_DP, ilu_fac, info)
    
    CALL gmres_solve_structured_preconditioned(workspace, n, row_ptr, col_idx, val_csr, &
                                               b, x_initial, 50, tol, ilu_fac, &
                                               x, iterations, residual_norm, converged, info)
    
    IF (.NOT. converged) THEN
      test_passed = .FALSE.
      EXIT
    END IF
    
    CALL ilu_free(ilu_fac)
  END DO
  
  IF (test_passed) THEN
    WRITE(*,'(A)') '[PASS]'
    passed_tests = passed_tests + 1
  ELSE
    WRITE(*,'(A)') '[FAIL]'
  END IF
  
  ! Test 4: Ill-conditioned system
  WRITE(*,'(A)', ADVANCE='NO') 'Test 4: GMRES + ILU on ill-conditioned system... '
  total_tests = total_tests + 1
  
  ! Modify diagonal to create ill-conditioning
  DO i = 1, n
    j = row_ptr(i)
    DO WHILE (j < row_ptr(i+1))
      IF (col_idx(j) == i) THEN
        val_csr(j) = 2.0_DP + 0.1_DP * SIN(REAL(i, DP))  ! Variable diagonal
        EXIT
      END IF
      j = j + 1
    END DO
  END DO
  
  ! Setup ILU(2) for ill-conditioned system
  CALL ilu_factorize(n, row_ptr, col_idx, val_csr, 2, 1.0e-4_DP, ilu_fac, info)
  
  ! Recompute RHS for modified matrix
  b = 0.0_DP
  DO i = 1, n
    DO j = row_ptr(i), row_ptr(i+1) - 1
      b(i) = b(i) + val_csr(j) * x_exact(col_idx(j))
    END DO
  END DO
  
  ! Solve with preconditioned GMRES
  max_iter = 100
  CALL gmres_solve_structured_preconditioned(workspace, n, row_ptr, col_idx, val_csr, &
                                             b, x_initial, max_iter, tol, ilu_fac, &
                                             x, iterations, residual_norm, converged, info)
  
  error_norm = SQRT(SUM((x - x_exact)**2))
  test_passed = converged .AND. error_norm < tol * 1000.0_DP
  
  IF (test_passed) THEN
    WRITE(*,'(A)') '[PASS]'
    passed_tests = passed_tests + 1
  ELSE
    WRITE(*,'(A)') '[FAIL]'
  END IF
  
  CALL ilu_free(ilu_fac)
  
  ! Test 5: Left vs right preconditioning
  WRITE(*,'(A)', ADVANCE='NO') 'Test 5: Left preconditioning mode... '
  total_tests = total_tests + 1
  
  ! Reset to standard tridiagonal
  CALL csc_to_csr_real(n, n, nnz, pcol, irow, amat, row_ptr, col_idx, val_csr)
  CALL ilu_factorize(n, row_ptr, col_idx, val_csr, 1, 0.0_DP, ilu_fac, info)
  
  ! Recompute RHS for standard tridiagonal matrix
  b = 0.0_DP
  DO i = 1, n
    DO j = MAX(1, i-1), MIN(n, i+1)
      IF (i == j) THEN
        b(i) = b(i) + 2.0_DP * x_exact(j)
      ELSE
        b(i) = b(i) - 1.0_DP * x_exact(j)
      END IF
    END DO
  END DO
  
  ! Solve with left preconditioning (standard)
  CALL gmres_solve_structured_preconditioned(workspace, n, row_ptr, col_idx, val_csr, &
                                             b, x_initial, max_iter, tol, ilu_fac, &
                                             x, iterations, residual_norm, converged, info, &
                                             precond_side='left')
  
  error_norm = SQRT(SUM((x - x_exact)**2))
  test_passed = converged .AND. error_norm < tol * 100.0_DP
  
  IF (test_passed) THEN
    WRITE(*,'(A)') '[PASS]'
    passed_tests = passed_tests + 1
  ELSE
    WRITE(*,'(A)') '[FAIL]'
  END IF
  
  CALL ilu_free(ilu_fac)
  
  WRITE(*,*)
  WRITE(*,'(A)') '============================================='
  WRITE(*,'(A)') 'GMRES ILU Integration Test Results'
  WRITE(*,'(A)') '============================================='
  WRITE(*,'(A,I0,A,I0,A)') 'Tests passed: ', passed_tests, ' / ', total_tests, ' total'
  WRITE(*,'(A,F6.1,A)') 'Success rate: ', REAL(passed_tests) / REAL(total_tests) * 100.0, '%'
  
  IF (passed_tests == total_tests) THEN
    WRITE(*,'(A)') 'All tests PASSED - GMRES ILU integration working correctly!'
  ELSE
    WRITE(*,'(A)') 'Some tests FAILED - Implementation needs work'
    WRITE(*,'(A)') 'This is expected in TDD RED phase before implementation'
  END IF
  
  WRITE(*,'(A)') '============================================='
  
  ! Cleanup
  CALL destroy_gmres_workspace(workspace)
  DEALLOCATE(irow, pcol, amat, row_ptr, col_idx, val_csr)
  DEALLOCATE(b, x, x_exact, x_initial)
  
END PROGRAM test_gmres_ilu_integration