PROGRAM test_gmres_arnoldi
  ! TDD RED phase: Test Arnoldi orthogonalization process
  ! Based on IterativeSolvers.jl template - Modified Gram-Schmidt
  
  USE nrtype, ONLY: I4B, DP
  USE gmres_mod, ONLY: arnoldi_decomp, create_arnoldi_decomp, destroy_arnoldi_decomp, &
                       arnoldi_expand, modified_gram_schmidt, check_orthogonality
  IMPLICIT NONE
  
  ! Test parameters
  INTEGER(I4B), PARAMETER :: n = 20       ! Small matrix size for testing
  INTEGER(I4B), PARAMETER :: restart = 5  ! Small restart for focused testing
  REAL(DP), PARAMETER :: tol = 1.0e-12_DP ! Orthogonality tolerance
  
  TYPE(arnoldi_decomp) :: arnoldi
  REAL(DP), ALLOCATABLE :: A(:,:)          ! Test matrix (dense for verification)
  REAL(DP), ALLOCATABLE :: v_input(:)      ! Input vector for expansion
  REAL(DP), ALLOCATABLE :: v_output(:)     ! Output vector from expansion
  REAL(DP) :: h_value                      ! Hessenberg entry
  REAL(DP) :: orthogonality_error         ! Measure of orthogonality
  REAL(DP), ALLOCATABLE :: h_coeffs(:)    ! Array for orthogonalization coefficients
  INTEGER(I4B) :: i, j, k
  LOGICAL :: test_passed
  INTEGER(I4B) :: total_tests, passed_tests
  
  WRITE(*,'(A)') '============================================='
  WRITE(*,'(A)') 'GMRES Arnoldi Orthogonalization Test Suite (TDD RED)'
  WRITE(*,'(A)') 'Testing Arnoldi process before implementation'
  WRITE(*,'(A)') '============================================='
  WRITE(*,*)
  
  total_tests = 0
  passed_tests = 0
  
  ! Setup test matrix and workspace
  CALL create_arnoldi_decomp(arnoldi, n, restart)
  ALLOCATE(A(n, n), v_input(n), v_output(n), h_coeffs(restart))
  
  ! Create a simple test matrix (tridiagonal)
  A = 0.0_DP
  DO i = 1, n
    A(i, i) = 2.0_DP                        ! Diagonal
    IF (i > 1) A(i, i-1) = -1.0_DP         ! Lower diagonal
    IF (i < n) A(i, i+1) = -1.0_DP         ! Upper diagonal
  END DO
  
  ! Initialize first Arnoldi vector (normalized random)
  DO i = 1, n
    arnoldi%V(i, 1) = REAL(i, DP)
  END DO
  h_value = SQRT(DOT_PRODUCT(arnoldi%V(:, 1), arnoldi%V(:, 1)))
  arnoldi%V(:, 1) = arnoldi%V(:, 1) / h_value
  
  ! Test 1: Modified Gram-Schmidt orthogonalization
  WRITE(*,'(A)', ADVANCE='NO') 'Test 1: Modified Gram-Schmidt orthogonalization... '
  total_tests = total_tests + 1
  
  ! Create a vector to orthogonalize against V(:, 1)
  v_input = 1.0_DP  ! All ones vector
  v_output = v_input
  
  ! Call modified Gram-Schmidt to orthogonalize v_output against V(:, 1)
  CALL modified_gram_schmidt(arnoldi%V(:, 1:1), v_output, h_coeffs(1:1))
  
  ! Check that v_output is orthogonal to V(:, 1)
  h_value = DOT_PRODUCT(v_output, arnoldi%V(:, 1))
  test_passed = ABS(h_value) < tol
  
  IF (test_passed) THEN
    WRITE(*,'(A)') '[PASS]'
    passed_tests = passed_tests + 1
  ELSE
    WRITE(*,'(A)') '[FAIL]'
  END IF
  
  ! Test 2: Arnoldi expansion step
  WRITE(*,'(A)', ADVANCE='NO') 'Test 2: Arnoldi expansion (A * v_k)... '
  total_tests = total_tests + 1
  
  ! Expand Arnoldi subspace: v_{k+1} = A * v_k
  k = 1
  CALL arnoldi_expand(A, arnoldi%V(:, k), v_output)
  
  ! Verify this is actually A * v_k
  v_input = MATMUL(A, arnoldi%V(:, k))  ! Reference result
  h_value = SQRT(SUM((v_output - v_input)**2))
  test_passed = h_value < tol
  
  IF (test_passed) THEN
    WRITE(*,'(A)') '[PASS]'
    passed_tests = passed_tests + 1
  ELSE
    WRITE(*,'(A)') '[FAIL]'
  END IF
  
  ! Test 3: Full Arnoldi step (expand + orthogonalize + normalize)
  WRITE(*,'(A)', ADVANCE='NO') 'Test 3: Complete Arnoldi step... '
  total_tests = total_tests + 1
  
  ! Build second Arnoldi vector
  k = 1
  CALL arnoldi_step(A, arnoldi, k)
  
  ! Check that H(1,1) and H(2,1) are computed correctly
  ! H(1,1) should be <A*v_1, v_1>
  h_value = DOT_PRODUCT(MATMUL(A, arnoldi%V(:, 1)), arnoldi%V(:, 1))
  test_passed = ABS(arnoldi%H(1, 1) - h_value) < tol
  
  ! Check that V(:, 2) is normalized
  h_value = SQRT(DOT_PRODUCT(arnoldi%V(:, 2), arnoldi%V(:, 2)))
  test_passed = test_passed .AND. ABS(h_value - 1.0_DP) < tol
  
  IF (test_passed) THEN
    WRITE(*,'(A)') '[PASS]'
    passed_tests = passed_tests + 1
  ELSE
    WRITE(*,'(A)') '[FAIL]'
  END IF
  
  ! Test 4: Multi-step Arnoldi process
  WRITE(*,'(A)', ADVANCE='NO') 'Test 4: Multi-step Arnoldi process... '
  total_tests = total_tests + 1
  
  ! Build full Arnoldi basis
  DO k = 1, restart
    CALL arnoldi_step(A, arnoldi, k)
  END DO
  
  ! Verify Arnoldi relation: A * V_k = V_{k+1} * H_k
  test_passed = .TRUE.
  DO j = 1, restart
    v_output = MATMUL(A, arnoldi%V(:, j))
    v_input = 0.0_DP
    DO i = 1, j + 1
      v_input = v_input + arnoldi%H(i, j) * arnoldi%V(:, i)
    END DO
    h_value = SQRT(SUM((v_output - v_input)**2))
    IF (h_value > tol) THEN
      test_passed = .FALSE.
      EXIT
    END IF
  END DO
  
  IF (test_passed) THEN
    WRITE(*,'(A)') '[PASS]'
    passed_tests = passed_tests + 1
  ELSE
    WRITE(*,'(A)') '[FAIL]'
  END IF
  
  ! Test 5: Orthogonality verification
  WRITE(*,'(A)', ADVANCE='NO') 'Test 5: Orthogonality of Arnoldi basis... '
  total_tests = total_tests + 1
  
  CALL check_orthogonality(arnoldi%V(:, 1:restart), orthogonality_error)
  test_passed = orthogonality_error < tol * 10.0_DP  ! Allow some numerical error
  
  IF (test_passed) THEN
    WRITE(*,'(A)') '[PASS]'
    passed_tests = passed_tests + 1
  ELSE
    WRITE(*,'(A)') '[FAIL]'
  END IF
  
  ! Test 6: Sparse matrix compatibility
  WRITE(*,'(A)', ADVANCE='NO') 'Test 6: Sparse matrix Arnoldi expansion... '
  total_tests = total_tests + 1
  
  ! Test with CSR format (simplified)
  CALL arnoldi_expand_csr(n, [1, 3, 5, 7, 9, 11], & ! Row pointers for tridiagonal
                          [1, 2, 1, 2, 3, 2, 3, 4, 3, 4], & ! Column indices  
                          [2.0_DP, -1.0_DP, -1.0_DP, 2.0_DP, -1.0_DP, &
                           -1.0_DP, 2.0_DP, -1.0_DP, -1.0_DP, 2.0_DP], & ! Values
                          arnoldi%V(:, 1), v_output)
  
  ! Should match dense matrix result
  v_input = MATMUL(A(1:4, 1:4), arnoldi%V(1:4, 1))  ! Only first 4x4 block
  h_value = SQRT(SUM((v_output(1:4) - v_input)**2))
  test_passed = h_value < tol
  
  IF (test_passed) THEN
    WRITE(*,'(A)') '[PASS]'
    passed_tests = passed_tests + 1
  ELSE
    WRITE(*,'(A)') '[FAIL]'
  END IF
  
  WRITE(*,*)
  WRITE(*,'(A)') '============================================='
  WRITE(*,'(A)') 'GMRES Arnoldi Test Results'
  WRITE(*,'(A)') '============================================='
  WRITE(*,'(A,I0,A,I0,A)') 'Tests passed: ', passed_tests, ' / ', total_tests, ' total'
  WRITE(*,'(A,F6.1,A)') 'Success rate: ', REAL(passed_tests) / REAL(total_tests) * 100.0, '%'
  
  IF (passed_tests == total_tests) THEN
    WRITE(*,'(A)') 'All tests PASSED - Arnoldi process working correctly!'
  ELSE
    WRITE(*,'(A)') 'Some tests FAILED - Implementation needs work'
    WRITE(*,'(A)') 'This is expected in TDD RED phase before implementation'
  END IF
  
  WRITE(*,'(A)') '============================================='
  
  ! Cleanup
  CALL destroy_arnoldi_decomp(arnoldi)
  DEALLOCATE(A, v_input, v_output, h_coeffs)
  
END PROGRAM test_gmres_arnoldi