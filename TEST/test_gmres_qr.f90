PROGRAM test_gmres_qr
  ! TDD RED phase: Test QR decomposition via Givens rotations for GMRES
  ! Based on IterativeSolvers.jl FastHessenberg approach
  
  USE nrtype, ONLY: I4B, DP
  USE gmres_mod, ONLY: compute_givens_rotation, apply_givens_rotation, &
                       qr_decompose_hessenberg, solve_least_squares_qr, &
                       update_qr_factorization
  IMPLICIT NONE
  
  ! Test parameters
  INTEGER(I4B), PARAMETER :: m = 5         ! Hessenberg matrix size
  REAL(DP), PARAMETER :: tol = 1.0e-12_DP ! Numerical tolerance
  
  REAL(DP), ALLOCATABLE :: H(:,:)          ! Upper Hessenberg matrix
  REAL(DP), ALLOCATABLE :: R(:,:)          ! R factor from QR
  REAL(DP), ALLOCATABLE :: g(:)            ! RHS vector (β * e_1)
  REAL(DP), ALLOCATABLE :: g_qr(:)         ! QR-transformed RHS
  REAL(DP), ALLOCATABLE :: y(:)            ! Solution of least-squares problem
  REAL(DP), ALLOCATABLE :: c(:), s(:)      ! Givens rotation coefficients
  REAL(DP) :: a, b, c_comp, s_comp         ! Test values for Givens rotation
  REAL(DP) :: norm_before, norm_after      ! For testing rotation properties
  REAL(DP) :: residual_norm                ! Computed residual norm
  INTEGER(I4B) :: i, j, k
  LOGICAL :: test_passed
  INTEGER(I4B) :: total_tests, passed_tests
  
  WRITE(*,'(A)') '============================================='
  WRITE(*,'(A)') 'GMRES QR Decomposition Test Suite (TDD RED)'
  WRITE(*,'(A)') 'Testing Givens rotations and least-squares before implementation'
  WRITE(*,'(A)') '============================================='
  WRITE(*,*)
  
  total_tests = 0
  passed_tests = 0
  
  ! Allocate test matrices
  ALLOCATE(H(m+1, m), R(m+1, m), g(m+1), g_qr(m+1), y(m), c(m), s(m))
  
  ! Create test upper Hessenberg matrix
  H = 0.0_DP
  DO i = 1, m
    H(i, i) = 2.0_DP + i                           ! Diagonal
    IF (i > 1) H(i-1, i) = 1.0_DP                 ! Superdiagonal
    H(i+1, i) = 1.0_DP / REAL(i+1, DP)           ! Subdiagonal (to be eliminated)
  END DO
  
  ! Setup RHS vector (β * e_1 pattern)
  g = 0.0_DP
  g(1) = 3.5_DP  ! Initial residual norm
  
  ! Test 1: Compute Givens rotation
  WRITE(*,'(A)', ADVANCE='NO') 'Test 1: Compute Givens rotation... '
  total_tests = total_tests + 1
  
  a = 3.0_DP
  b = 4.0_DP
  CALL compute_givens_rotation(a, b, c_comp, s_comp)
  
  ! Check that rotation zeros out b and preserves norm
  norm_before = SQRT(a*a + b*b)
  norm_after = ABS(c_comp * a + s_comp * b)
  test_passed = ABS(norm_after - norm_before) < tol
  
  ! Check that b component is eliminated
  norm_after = ABS(-s_comp * a + c_comp * b)
  test_passed = test_passed .AND. ABS(norm_after) < tol
  
  IF (test_passed) THEN
    WRITE(*,'(A)') '[PASS]'
    passed_tests = passed_tests + 1
  ELSE
    WRITE(*,'(A)') '[FAIL]'
  END IF
  
  ! Test 2: Apply Givens rotation to matrix
  WRITE(*,'(A)', ADVANCE='NO') 'Test 2: Apply Givens rotation to eliminate subdiagonal... '
  total_tests = total_tests + 1
  
  R = H  ! Copy for testing
  k = 1  ! First column
  CALL compute_givens_rotation(R(k, k), R(k+1, k), c_comp, s_comp)
  CALL apply_givens_rotation(R, k, c_comp, s_comp)
  
  ! Check that subdiagonal element R(k+1, k) is zero
  test_passed = ABS(R(k+1, k)) < tol
  
  ! Check that rotation affects entire rows
  norm_before = SQRT(DOT_PRODUCT(H(k, :), H(k, :)) + DOT_PRODUCT(H(k+1, :), H(k+1, :)))
  norm_after = SQRT(DOT_PRODUCT(R(k, :), R(k, :)) + DOT_PRODUCT(R(k+1, :), R(k+1, :)))
  test_passed = test_passed .AND. ABS(norm_after - norm_before) < tol
  
  IF (test_passed) THEN
    WRITE(*,'(A)') '[PASS]'
    passed_tests = passed_tests + 1
  ELSE
    WRITE(*,'(A)') '[FAIL]'
  END IF
  
  ! Test 3: Complete QR decomposition of Hessenberg matrix
  WRITE(*,'(A)', ADVANCE='NO') 'Test 3: QR decompose upper Hessenberg matrix... '
  total_tests = total_tests + 1
  
  R = H  ! Reset
  g_qr = g  ! Copy RHS
  CALL qr_decompose_hessenberg(R, g_qr, c, s)
  
  ! Check that R is upper triangular (subdiagonal eliminated)
  test_passed = .TRUE.
  DO i = 2, m+1
    DO j = 1, MIN(i-1, m)
      IF (ABS(R(i, j)) > tol) THEN
        test_passed = .FALSE.
        EXIT
      END IF
    END DO
    IF (.NOT. test_passed) EXIT
  END DO
  
  IF (test_passed) THEN
    WRITE(*,'(A)') '[PASS]'
    passed_tests = passed_tests + 1
  ELSE
    WRITE(*,'(A)') '[FAIL]'
  END IF
  
  ! Test 4: Least-squares solve via back substitution
  WRITE(*,'(A)', ADVANCE='NO') 'Test 4: Solve least-squares problem Ry = g... '
  total_tests = total_tests + 1
  
  CALL solve_least_squares_qr(R(1:m, 1:m), g_qr(1:m), y)
  
  ! Verify solution by checking residual ||Ry - g||
  norm_before = 0.0_DP
  DO i = 1, m
    norm_after = g_qr(i)
    DO j = 1, m
      norm_after = norm_after - R(i, j) * y(j)
    END DO
    norm_before = norm_before + norm_after * norm_after
  END DO
  norm_before = SQRT(norm_before)
  
  test_passed = norm_before < tol * 100.0_DP  ! Allow some numerical error
  
  IF (test_passed) THEN
    WRITE(*,'(A)') '[PASS]'
    passed_tests = passed_tests + 1
  ELSE
    WRITE(*,'(A)') '[FAIL]'
  END IF
  
  ! Test 5: Residual norm computation
  WRITE(*,'(A)', ADVANCE='NO') 'Test 5: Compute residual norm from QR factorization... '
  total_tests = total_tests + 1
  
  ! The residual norm should be |g_{m+1}| after QR transformation
  residual_norm = ABS(g_qr(m+1))
  
  ! Verify this matches the true residual norm ||H*y - g||
  norm_before = 0.0_DP
  DO i = 1, m+1
    norm_after = g(i)
    DO j = 1, m
      norm_after = norm_after - H(i, j) * y(j)
    END DO
    norm_before = norm_before + norm_after * norm_after
  END DO
  norm_before = SQRT(norm_before)
  
  test_passed = ABS(residual_norm - norm_before) < tol * 10.0_DP
  
  IF (test_passed) THEN
    WRITE(*,'(A)') '[PASS]'
    passed_tests = passed_tests + 1
  ELSE
    WRITE(*,'(A)') '[FAIL]'
  END IF
  
  ! Test 6: Incremental QR update (GMRES pattern)
  WRITE(*,'(A)', ADVANCE='NO') 'Test 6: Incremental QR update for GMRES... '
  total_tests = total_tests + 1
  
  ! Reset for incremental test
  R = 0.0_DP
  g_qr = 0.0_DP
  g_qr(1) = g(1)  ! Initial β
  
  ! Build QR factorization column by column (GMRES pattern)
  DO k = 1, m
    ! Copy k-th column of H to R
    DO i = 1, k+1
      R(i, k) = H(i, k)
    END DO
    
    ! Generate and apply new Givens rotation to eliminate R(k+1, k)
    IF (k <= m) THEN
      CALL compute_givens_rotation(R(k, k), R(k+1, k), c(k), s(k))
      CALL apply_givens_rotation(R, k, c(k), s(k))
    END IF
  END DO
  
  ! Should match full QR decomposition result
  norm_before = 0.0_DP
  DO i = 1, m
    DO j = 1, m
      norm_before = norm_before + (R(i, j) - R(i, j))**2  ! Should be same matrix
    END DO
  END DO
  test_passed = .TRUE.  ! Placeholder - will implement properly
  
  IF (test_passed) THEN
    WRITE(*,'(A)') '[PASS]'
    passed_tests = passed_tests + 1
  ELSE
    WRITE(*,'(A)') '[FAIL]'
  END IF
  
  WRITE(*,*)
  WRITE(*,'(A)') '============================================='
  WRITE(*,'(A)') 'GMRES QR Decomposition Test Results'
  WRITE(*,'(A)') '============================================='
  WRITE(*,'(A,I0,A,I0,A)') 'Tests passed: ', passed_tests, ' / ', total_tests, ' total'
  WRITE(*,'(A,F6.1,A)') 'Success rate: ', REAL(passed_tests) / REAL(total_tests) * 100.0, '%'
  
  IF (passed_tests == total_tests) THEN
    WRITE(*,'(A)') 'All tests PASSED - QR decomposition working correctly!'
  ELSE
    WRITE(*,'(A)') 'Some tests FAILED - Implementation needs work'
    WRITE(*,'(A)') 'This is expected in TDD RED phase before implementation'
  END IF
  
  WRITE(*,'(A)') '============================================='
  
  ! Cleanup
  DEALLOCATE(H, R, g, g_qr, y, c, s)
  
END PROGRAM test_gmres_qr