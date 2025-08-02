PROGRAM test_gmres_complete
  ! TDD RED phase: Test complete restarted GMRES algorithm
  ! Based on IterativeSolvers.jl template - comprehensive GMRES implementation
  
  USE nrtype, ONLY: I4B, DP
  USE gmres_mod, ONLY: gmres_workspace, create_gmres_workspace, destroy_gmres_workspace, &
                       gmres_solve_structured, gmres_iterate, &
                       initialize_gmres, finalize_gmres_iteration
  IMPLICIT NONE
  
  ! Test parameters
  INTEGER(I4B), PARAMETER :: n = 50        ! Matrix size
  INTEGER(I4B), PARAMETER :: restart = 10  ! GMRES restart parameter  
  REAL(DP), PARAMETER :: tol = 1.0e-10_DP  ! Convergence tolerance
  
  TYPE(gmres_workspace) :: workspace
  REAL(DP), ALLOCATABLE :: A(:,:)          ! Test matrix (dense)
  REAL(DP), ALLOCATABLE :: b(:), x(:)      ! RHS and solution vectors
  REAL(DP), ALLOCATABLE :: x_exact(:)      ! Exact solution
  REAL(DP), ALLOCATABLE :: x_initial(:)    ! Initial guess
  REAL(DP) :: residual_norm, error_norm
  INTEGER(I4B) :: max_iter, iterations
  LOGICAL :: converged, test_passed
  INTEGER(I4B) :: i, j, info
  INTEGER(I4B) :: total_tests, passed_tests
  
  WRITE(*,'(A)') '============================================='
  WRITE(*,'(A)') 'Complete GMRES Algorithm Test Suite (TDD RED)'
  WRITE(*,'(A)') 'Testing full restarted GMRES before implementation'
  WRITE(*,'(A)') '============================================='
  WRITE(*,*)
  
  total_tests = 0
  passed_tests = 0
  
  ! Setup workspace and test problem
  CALL create_gmres_workspace(workspace, n, restart)
  ALLOCATE(A(n, n), b(n), x(n), x_exact(n), x_initial(n))
  
  ! Test 1: Diagonal system (should converge in 1 iteration)
  WRITE(*,'(A)', ADVANCE='NO') 'Test 1: Diagonal system solve... '
  total_tests = total_tests + 1
  
  ! Create diagonal system A = diag(1, 2, 3, ..., n)
  A = 0.0_DP
  DO i = 1, n
    A(i, i) = REAL(i, DP)
  END DO
  
  ! Exact solution
  x_exact = 1.0_DP
  b = MATMUL(A, x_exact)
  
  ! Initial guess
  x_initial = 0.0_DP
  
  ! Solve with structured GMRES
  max_iter = restart * 2
  CALL gmres_solve_structured(workspace, A, b, x_initial, max_iter, tol, &
                              x, iterations, residual_norm, converged, info)
  
  ! Check convergence and accuracy
  error_norm = SQRT(SUM((x - x_exact)**2))
  test_passed = converged .AND. error_norm < tol * 10.0_DP
  
  IF (test_passed) THEN
    WRITE(*,'(A)') '[PASS]'
    passed_tests = passed_tests + 1
  ELSE
    WRITE(*,'(A)') '[FAIL]'
  END IF
  
  ! Test 2: Tridiagonal system 
  WRITE(*,'(A)', ADVANCE='NO') 'Test 2: Tridiagonal system solve... '
  total_tests = total_tests + 1
  
  ! Create tridiagonal system: -1, 2, -1
  A = 0.0_DP
  DO i = 1, n
    A(i, i) = 2.0_DP
    IF (i > 1) A(i, i-1) = -1.0_DP
    IF (i < n) A(i, i+1) = -1.0_DP
  END DO
  
  ! Setup known solution
  DO i = 1, n
    x_exact(i) = SIN(REAL(i, DP) * 3.14159_DP / REAL(n+1, DP))
  END DO
  b = MATMUL(A, x_exact)
  x_initial = 0.0_DP
  
  ! Solve with GMRES
  max_iter = restart * 5
  CALL gmres_solve_structured(workspace, A, b, x_initial, max_iter, tol, &
                              x, iterations, residual_norm, converged, info)
  
  error_norm = SQRT(SUM((x - x_exact)**2))
  test_passed = converged .AND. error_norm < tol * 100.0_DP
  
  IF (test_passed) THEN
    WRITE(*,'(A)') '[PASS]'
    passed_tests = passed_tests + 1
  ELSE
    WRITE(*,'(A)') '[FAIL]'
  END IF
  
  ! Test 3: Restart behavior
  WRITE(*,'(A)', ADVANCE='NO') 'Test 3: GMRES restart behavior... '
  total_tests = total_tests + 1
  
  ! Use same tridiagonal system but force restarts
  max_iter = restart * 3  ! Allow multiple restarts
  CALL gmres_solve_structured(workspace, A, b, x_initial, max_iter, tol, &
                              x, iterations, residual_norm, converged, info)
  
  ! Should still converge even with restarts
  error_norm = SQRT(SUM((x - x_exact)**2))
  test_passed = converged .AND. error_norm < tol * 100.0_DP .AND. iterations > restart
  
  IF (test_passed) THEN
    WRITE(*,'(A)') '[PASS]'
    passed_tests = passed_tests + 1
  ELSE
    WRITE(*,'(A)') '[FAIL]'
  END IF
  
  ! Test 4: GMRES iteration interface
  WRITE(*,'(A)', ADVANCE='NO') 'Test 4: GMRES iteration interface... '
  total_tests = total_tests + 1
  
  ! Initialize for iterative solve
  x = x_initial
  CALL initialize_gmres(workspace, A, b, x, residual_norm)
  
  converged = .FALSE.
  iterations = 0
  
  ! Manual iteration loop
  DO WHILE (.NOT. converged .AND. iterations < max_iter)
    CALL gmres_iterate(workspace, A, converged, residual_norm, info)
    iterations = iterations + 1
    
    IF (residual_norm < tol) converged = .TRUE.
    
    ! Check for restart
    IF (workspace%k > restart) THEN
      CALL finalize_gmres_iteration(workspace, x)
      IF (.NOT. converged) THEN
        CALL initialize_gmres(workspace, A, b, x, residual_norm)
      END IF
    END IF
  END DO
  
  IF (converged) THEN
    CALL finalize_gmres_iteration(workspace, x)
  END IF
  
  error_norm = SQRT(SUM((x - x_exact)**2))
  test_passed = converged .AND. error_norm < tol * 100.0_DP
  
  IF (test_passed) THEN
    WRITE(*,'(A)') '[PASS]'
    passed_tests = passed_tests + 1
  ELSE
    WRITE(*,'(A)') '[FAIL]'
  END IF
  
  ! Test 5: Non-symmetric system
  WRITE(*,'(A)', ADVANCE='NO') 'Test 5: Non-symmetric system solve... '
  total_tests = total_tests + 1
  
  ! Create non-symmetric system
  A = 0.0_DP
  DO i = 1, n
    A(i, i) = 3.0_DP
    IF (i > 1) A(i, i-1) = -1.0_DP
    IF (i < n) A(i, i+1) = -2.0_DP  ! Different from lower diagonal
  END DO
  
  x_exact = 1.0_DP
  b = MATMUL(A, x_exact)
  x_initial = 0.0_DP
  
  max_iter = n  ! Allow full Krylov space
  CALL gmres_solve_structured(workspace, A, b, x_initial, max_iter, tol, &
                              x, iterations, residual_norm, converged, info)
  
  error_norm = SQRT(SUM((x - x_exact)**2))
  test_passed = converged .AND. error_norm < tol * 100.0_DP
  
  IF (test_passed) THEN
    WRITE(*,'(A)') '[PASS]'
    passed_tests = passed_tests + 1
  ELSE
    WRITE(*,'(A)') '[FAIL]'
  END IF
  
  ! Test 6: Convergence with initial guess
  WRITE(*,'(A)', ADVANCE='NO') 'Test 6: Convergence with good initial guess... '
  total_tests = total_tests + 1
  
  ! Use diagonal system with good initial guess
  A = 0.0_DP
  DO i = 1, n
    A(i, i) = REAL(i, DP)
  END DO
  
  x_exact = 1.0_DP
  b = MATMUL(A, x_exact)
  x_initial = 0.9_DP  ! Close to exact solution
  
  max_iter = 10
  CALL gmres_solve_structured(workspace, A, b, x_initial, max_iter, tol, &
                              x, iterations, residual_norm, converged, info)
  
  ! Should converge quickly with good initial guess
  error_norm = SQRT(SUM((x - x_exact)**2))
  test_passed = converged .AND. error_norm < tol * 10.0_DP .AND. iterations <= 5
  
  IF (test_passed) THEN
    WRITE(*,'(A)') '[PASS]'
    passed_tests = passed_tests + 1
  ELSE
    WRITE(*,'(A)') '[FAIL]'
  END IF
  
  WRITE(*,*)
  WRITE(*,'(A)') '============================================='
  WRITE(*,'(A)') 'Complete GMRES Algorithm Test Results'
  WRITE(*,'(A)') '============================================='
  WRITE(*,'(A,I0,A,I0,A)') 'Tests passed: ', passed_tests, ' / ', total_tests, ' total'
  WRITE(*,'(A,F6.1,A)') 'Success rate: ', REAL(passed_tests) / REAL(total_tests) * 100.0, '%'
  
  IF (passed_tests == total_tests) THEN
    WRITE(*,'(A)') 'All tests PASSED - Complete GMRES working correctly!'
  ELSE
    WRITE(*,'(A)') 'Some tests FAILED - Implementation needs work'
    WRITE(*,'(A)') 'This is expected in TDD RED phase before implementation'
  END IF
  
  WRITE(*,'(A)') '============================================='
  
  ! Cleanup
  CALL destroy_gmres_workspace(workspace)
  DEALLOCATE(A, b, x, x_exact, x_initial)
  
END PROGRAM test_gmres_complete