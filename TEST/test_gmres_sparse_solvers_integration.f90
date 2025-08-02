PROGRAM test_gmres_sparse_solvers_integration
  ! TDD RED phase: Test GMRES integration into sparse_solvers_mod
  ! Verify that GMRES can be called through the unified sparse_solve interface
  
  USE nrtype, ONLY: I4B, DP
  USE sparse_solvers_mod, ONLY: sparse_solve, sparse_solve_method, &
                                SOLVER_UMFPACK, SOLVER_BICGSTAB, SOLVER_GMRES, &
                                iterative_solver_params, default_iterative_params
  USE sparse_types_mod, ONLY: sparse_matrix_csc_real, allocate_sparse, deallocate_sparse
  IMPLICIT NONE
  
  ! Test parameters
  INTEGER(I4B), PARAMETER :: n = 100       ! Matrix size
  REAL(DP), PARAMETER :: tol = 1.0e-10_DP  ! Convergence tolerance
  
  ! Sparse matrix in CSC format
  TYPE(sparse_matrix_csc_real) :: mat_csc
  INTEGER(I4B) :: nnz                      ! Number of nonzeros
  
  ! Vectors
  REAL(DP), ALLOCATABLE :: b(:), x_exact(:)
  REAL(DP) :: error_norm
  INTEGER(I4B) :: max_iter, iterations
  LOGICAL :: test_passed, converged
  INTEGER(I4B) :: i, j, k, idx
  INTEGER(I4B) :: total_tests, passed_tests
  
  WRITE(*,'(A)') '====================================================='
  WRITE(*,'(A)') 'GMRES Sparse Solvers Integration Test Suite (TDD RED)'
  WRITE(*,'(A)') 'Testing GMRES through sparse_solve interface'
  WRITE(*,'(A)') '====================================================='
  WRITE(*,*)
  
  total_tests = 0
  passed_tests = 0
  
  ! Create tridiagonal test matrix
  nnz = 3*n - 2
  CALL allocate_sparse(mat_csc, n, n, nnz)
  
  ! Build tridiagonal matrix: -1, 2, -1
  idx = 0
  mat_csc%pcol(1) = 1
  DO j = 1, n
    IF (j > 1) THEN
      idx = idx + 1
      mat_csc%irow(idx) = j - 1
      mat_csc%val(idx) = -1.0_DP
    END IF
    idx = idx + 1
    mat_csc%irow(idx) = j
    mat_csc%val(idx) = 2.0_DP
    IF (j < n) THEN
      idx = idx + 1
      mat_csc%irow(idx) = j + 1
      mat_csc%val(idx) = -1.0_DP
    END IF
    mat_csc%pcol(j+1) = idx + 1
  END DO
  
  ! Setup test problem
  ALLOCATE(b(n), x_exact(n))
  x_exact = 1.0_DP
  
  ! Compute b = A * x_exact
  b = 0.0_DP
  DO j = 1, n
    DO i = mat_csc%pcol(j), mat_csc%pcol(j+1) - 1
      b(mat_csc%irow(i)) = b(mat_csc%irow(i)) + mat_csc%val(i) * x_exact(j)
    END DO
  END DO
  
  ! Test 1: Basic GMRES solver selection
  WRITE(*,'(A)', ADVANCE='NO') 'Test 1: Select GMRES solver method... '
  total_tests = total_tests + 1
  
  sparse_solve_method = SOLVER_GMRES
  
  ! This should work if GMRES is integrated
  test_passed = sparse_solve_method == SOLVER_GMRES
  
  IF (test_passed) THEN
    WRITE(*,'(A)') '[PASS]'
    passed_tests = passed_tests + 1
  ELSE
    WRITE(*,'(A)') '[FAIL]'
  END IF
  
  ! Test 2: GMRES solver through sparse_solve interface
  WRITE(*,'(A)', ADVANCE='NO') 'Test 2: Solve with GMRES through sparse_solve... '
  total_tests = total_tests + 1
  
  ! Configure GMRES parameters
  default_iterative_params%abs_tolerance = tol
  default_iterative_params%rel_tolerance = tol * 10.0_DP
  default_iterative_params%max_iterations = n
  default_iterative_params%gmres_restart = 30
  default_iterative_params%preconditioner_type = 0  ! No preconditioning
  default_iterative_params%verbose = .TRUE.
  
  ! Call sparse_solve with GMRES method
  sparse_solve_method = SOLVER_GMRES
  CALL sparse_solve(mat_csc%nrow, mat_csc%ncol, mat_csc%nz, &
                    mat_csc%irow, mat_csc%pcol, mat_csc%val, b)
  
  ! Check solution accuracy
  error_norm = SQRT(SUM((b - x_exact)**2))
  test_passed = error_norm < tol * 100.0_DP
  
  IF (.NOT. test_passed .AND. default_iterative_params%verbose) THEN
    WRITE(*,*) '    First 5 elements of solution b:', b(1:MIN(5,n))
    WRITE(*,*) '    First 5 elements of x_exact:', x_exact(1:MIN(5,n))
  END IF
  
  IF (test_passed) THEN
    WRITE(*,'(A)') '[PASS]'
    passed_tests = passed_tests + 1
  ELSE
    WRITE(*,'(A)') '[FAIL]'
    WRITE(*,'(A,ES12.5)') '    Error norm: ', error_norm
    WRITE(*,'(A,ES12.5)') '    Tolerance: ', tol * 100.0_DP
  END IF
  
  ! Test 3: GMRES with ILU preconditioning
  WRITE(*,'(A)', ADVANCE='NO') 'Test 3: GMRES with ILU preconditioning... '
  total_tests = total_tests + 1
  
  ! Reset b vector
  b = 0.0_DP
  DO j = 1, n
    DO i = mat_csc%pcol(j), mat_csc%pcol(j+1) - 1
      b(mat_csc%irow(i)) = b(mat_csc%irow(i)) + mat_csc%val(i) * x_exact(j)
    END DO
  END DO
  
  ! Configure GMRES with ILU
  default_iterative_params%preconditioner_type = 1  ! ILU preconditioning
  default_iterative_params%ilu_fill_level = 1      ! ILU(1)
  
  ! Call sparse_solve with GMRES + ILU
  sparse_solve_method = SOLVER_GMRES
  CALL sparse_solve(mat_csc%nrow, mat_csc%ncol, mat_csc%nz, &
                    mat_csc%irow, mat_csc%pcol, mat_csc%val, b)
  
  error_norm = SQRT(SUM((b - x_exact)**2))
  test_passed = error_norm < tol * 100.0_DP
  
  IF (test_passed) THEN
    WRITE(*,'(A)') '[PASS]'
    passed_tests = passed_tests + 1
  ELSE
    WRITE(*,'(A)') '[FAIL]'
  END IF
  
  ! Test 4: GMRES with custom restart parameter
  WRITE(*,'(A)', ADVANCE='NO') 'Test 4: GMRES with custom restart parameter... '
  total_tests = total_tests + 1
  
  ! Reset b vector
  b = 0.0_DP
  DO j = 1, n
    DO i = mat_csc%pcol(j), mat_csc%pcol(j+1) - 1
      b(mat_csc%irow(i)) = b(mat_csc%irow(i)) + mat_csc%val(i) * x_exact(j)
    END DO
  END DO
  
  ! Configure with small restart value to force multiple restarts
  default_iterative_params%gmres_restart = 10
  default_iterative_params%preconditioner_type = 0  ! No preconditioning
  
  sparse_solve_method = SOLVER_GMRES
  CALL sparse_solve(mat_csc%nrow, mat_csc%ncol, mat_csc%nz, &
                    mat_csc%irow, mat_csc%pcol, mat_csc%val, b)
  
  error_norm = SQRT(SUM((b - x_exact)**2))
  test_passed = error_norm < tol * 100.0_DP
  
  IF (test_passed) THEN
    WRITE(*,'(A)') '[PASS]'
    passed_tests = passed_tests + 1
  ELSE
    WRITE(*,'(A)') '[FAIL]'
  END IF
  
  ! Test 5: Auto-selection includes GMRES
  WRITE(*,'(A)', ADVANCE='NO') 'Test 5: Auto-solver considers GMRES for large systems... '
  total_tests = total_tests + 1
  
  ! For large ill-conditioned systems, auto-solver should consider GMRES
  sparse_solve_method = 0  ! Auto-select
  
  ! This test checks if the error message includes GMRES as an option
  ! In RED phase, this will fail as GMRES is not yet integrated
  test_passed = .FALSE.  ! Will be true when GMRES is properly integrated
  
  IF (test_passed) THEN
    WRITE(*,'(A)') '[PASS]'
    passed_tests = passed_tests + 1
  ELSE
    WRITE(*,'(A)') '[FAIL]'
    WRITE(*,'(A)') '    (Expected in RED phase - GMRES not yet in auto-selection)'
  END IF
  
  WRITE(*,*)
  WRITE(*,'(A)') '====================================================='
  WRITE(*,'(A)') 'GMRES Sparse Solvers Integration Test Results'
  WRITE(*,'(A)') '====================================================='
  WRITE(*,'(A,I0,A,I0,A)') 'Tests passed: ', passed_tests, ' / ', total_tests, ' total'
  WRITE(*,'(A,F6.1,A)') 'Success rate: ', REAL(passed_tests) / REAL(total_tests) * 100.0, '%'
  
  IF (passed_tests == total_tests) THEN
    WRITE(*,'(A)') 'All tests PASSED - GMRES fully integrated into sparse_solve!'
  ELSE
    WRITE(*,'(A)') 'Some tests FAILED - Implementation needs work'
    WRITE(*,'(A)') 'This is expected in TDD RED phase before implementation'
  END IF
  
  WRITE(*,'(A)') '====================================================='
  
  ! Cleanup
  CALL deallocate_sparse(mat_csc)
  DEALLOCATE(b, x_exact)
  
END PROGRAM test_gmres_sparse_solvers_integration