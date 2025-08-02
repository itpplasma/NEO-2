PROGRAM test_gmres_data_structures
  ! TDD RED phase: Test GMRES core data structures
  ! Based on IterativeSolvers.jl template structure
  
  USE nrtype, ONLY: I4B, DP
  USE gmres_mod, ONLY: arnoldi_decomp, gmres_workspace, &
                       create_arnoldi_decomp, destroy_arnoldi_decomp, &
                       create_gmres_workspace, destroy_gmres_workspace
  IMPLICIT NONE
  
  ! Test parameters
  INTEGER(I4B), PARAMETER :: n = 100      ! Matrix size
  INTEGER(I4B), PARAMETER :: restart = 30 ! Restart dimension
  
  TYPE(arnoldi_decomp) :: arnoldi
  TYPE(gmres_workspace) :: workspace
  INTEGER(I4B) :: i, j
  LOGICAL :: test_passed
  INTEGER(I4B) :: total_tests, passed_tests
  
  WRITE(*,'(A)') '============================================='
  WRITE(*,'(A)') 'GMRES Data Structures Test Suite (TDD RED)'
  WRITE(*,'(A)') 'Testing core data structures before implementation'
  WRITE(*,'(A)') '============================================='
  WRITE(*,*)
  
  total_tests = 0
  passed_tests = 0
  
  ! Test 1: Arnoldi decomposition creation
  WRITE(*,'(A)', ADVANCE='NO') 'Test 1: Create Arnoldi decomposition... '
  total_tests = total_tests + 1
  
  CALL create_arnoldi_decomp(arnoldi, n, restart)
  
  ! Check V matrix dimensions (n x restart+1)
  test_passed = .TRUE.
  IF (.NOT. ALLOCATED(arnoldi%V)) test_passed = .FALSE.
  IF (test_passed .AND. SIZE(arnoldi%V, 1) /= n) test_passed = .FALSE.
  IF (test_passed .AND. SIZE(arnoldi%V, 2) /= restart + 1) test_passed = .FALSE.
  
  ! Check H matrix dimensions (restart+1 x restart)  
  IF (.NOT. ALLOCATED(arnoldi%H)) test_passed = .FALSE.
  IF (test_passed .AND. SIZE(arnoldi%H, 1) /= restart + 1) test_passed = .FALSE.
  IF (test_passed .AND. SIZE(arnoldi%H, 2) /= restart) test_passed = .FALSE.
  
  ! Check order parameter
  IF (test_passed .AND. arnoldi%order /= restart) test_passed = .FALSE.
  
  IF (test_passed) THEN
    WRITE(*,'(A)') '[PASS]'
    passed_tests = passed_tests + 1
  ELSE
    WRITE(*,'(A)') '[FAIL]'
  END IF
  
  ! Test 2: Arnoldi matrices initialized to zero
  WRITE(*,'(A)', ADVANCE='NO') 'Test 2: Arnoldi matrices zero-initialized... '
  total_tests = total_tests + 1
  
  test_passed = .TRUE.
  DO i = 1, n
    DO j = 1, restart + 1
      IF (ABS(arnoldi%V(i, j)) > 1.0e-15_DP) THEN
        test_passed = .FALSE.
        EXIT
      END IF
    END DO
    IF (.NOT. test_passed) EXIT
  END DO
  
  IF (test_passed) THEN
    DO i = 1, restart + 1
      DO j = 1, restart
        IF (ABS(arnoldi%H(i, j)) > 1.0e-15_DP) THEN
          test_passed = .FALSE.
          EXIT
        END IF
      END DO
      IF (.NOT. test_passed) EXIT
    END DO
  END IF
  
  IF (test_passed) THEN
    WRITE(*,'(A)') '[PASS]'
    passed_tests = passed_tests + 1
  ELSE
    WRITE(*,'(A)') '[FAIL]'
  END IF
  
  ! Test 3: GMRES workspace creation
  WRITE(*,'(A)', ADVANCE='NO') 'Test 3: Create GMRES workspace... '
  total_tests = total_tests + 1
  
  CALL create_gmres_workspace(workspace, n, restart)
  
  test_passed = .TRUE.
  
  ! Check that Arnoldi is properly embedded
  IF (.NOT. ALLOCATED(workspace%arnoldi%V)) test_passed = .FALSE.
  IF (.NOT. ALLOCATED(workspace%arnoldi%H)) test_passed = .FALSE.
  
  ! Check Givens rotation arrays
  IF (.NOT. ALLOCATED(workspace%givens_c)) test_passed = .FALSE.
  IF (.NOT. ALLOCATED(workspace%givens_s)) test_passed = .FALSE.
  IF (test_passed .AND. SIZE(workspace%givens_c) /= restart) test_passed = .FALSE.
  IF (test_passed .AND. SIZE(workspace%givens_s) /= restart) test_passed = .FALSE.
  
  ! Check QR RHS vector  
  IF (.NOT. ALLOCATED(workspace%rhs_qr)) test_passed = .FALSE.
  IF (test_passed .AND. SIZE(workspace%rhs_qr) /= restart + 1) test_passed = .FALSE.
  
  ! Check initialization values
  IF (test_passed .AND. ABS(workspace%residual_norm) > 1.0e-15_DP) test_passed = .FALSE.
  IF (test_passed .AND. workspace%k /= 0) test_passed = .FALSE.
  
  IF (test_passed) THEN
    WRITE(*,'(A)') '[PASS]'
    passed_tests = passed_tests + 1
  ELSE
    WRITE(*,'(A)') '[FAIL]'
  END IF
  
  ! Test 4: Memory cleanup - Arnoldi
  WRITE(*,'(A)', ADVANCE='NO') 'Test 4: Destroy Arnoldi decomposition... '
  total_tests = total_tests + 1
  
  CALL destroy_arnoldi_decomp(arnoldi)
  
  test_passed = .TRUE.
  IF (ALLOCATED(arnoldi%V)) test_passed = .FALSE.
  IF (ALLOCATED(arnoldi%H)) test_passed = .FALSE.
  
  IF (test_passed) THEN
    WRITE(*,'(A)') '[PASS]'
    passed_tests = passed_tests + 1
  ELSE
    WRITE(*,'(A)') '[FAIL]'
  END IF
  
  ! Test 5: Memory cleanup - GMRES workspace
  WRITE(*,'(A)', ADVANCE='NO') 'Test 5: Destroy GMRES workspace... '
  total_tests = total_tests + 1
  
  CALL destroy_gmres_workspace(workspace)
  
  test_passed = .TRUE.
  IF (ALLOCATED(workspace%arnoldi%V)) test_passed = .FALSE.
  IF (ALLOCATED(workspace%arnoldi%H)) test_passed = .FALSE.
  IF (ALLOCATED(workspace%givens_c)) test_passed = .FALSE.
  IF (ALLOCATED(workspace%givens_s)) test_passed = .FALSE.
  IF (ALLOCATED(workspace%rhs_qr)) test_passed = .FALSE.
  
  IF (test_passed) THEN
    WRITE(*,'(A)') '[PASS]'
    passed_tests = passed_tests + 1
  ELSE
    WRITE(*,'(A)') '[FAIL]'
  END IF
  
  ! Test 6: Real and complex variants
  WRITE(*,'(A)', ADVANCE='NO') 'Test 6: Complex data structure support... '
  total_tests = total_tests + 1
  
  ! This will test complex versions when implemented
  test_passed = .FALSE.  ! Will fail until complex support added
  
  IF (test_passed) THEN
    WRITE(*,'(A)') '[PASS]'
    passed_tests = passed_tests + 1
  ELSE
    WRITE(*,'(A)') '[FAIL - Complex support not yet implemented]'
  END IF
  
  WRITE(*,*)
  WRITE(*,'(A)') '============================================='
  WRITE(*,'(A)') 'GMRES Data Structures Test Results'
  WRITE(*,'(A)') '============================================='
  WRITE(*,'(A,I0,A,I0,A)') 'Tests passed: ', passed_tests, ' / ', total_tests, ' total'
  WRITE(*,'(A,F6.1,A)') 'Success rate: ', REAL(passed_tests) / REAL(total_tests) * 100.0, '%'
  
  IF (passed_tests == total_tests) THEN
    WRITE(*,'(A)') 'All tests PASSED - Data structures working correctly!'
  ELSE
    WRITE(*,'(A)') 'Some tests FAILED - Implementation needs work'
    WRITE(*,'(A)') 'This is expected in TDD RED phase before implementation'
  END IF
  
  WRITE(*,'(A)') '============================================='
  
  ! Expected to fail with compilation errors since gmres_mod doesn't exist yet
  
END PROGRAM test_gmres_data_structures