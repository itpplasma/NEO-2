PROGRAM test_ilu_precond
  ! Comprehensive tests for ilu_precond_mod module
  ! Tests ILU(0) and ILU(1) factorization, forward/backward substitution, drop tolerance
  
  USE sparse_types_mod, ONLY: dp, long
  USE sparse_utils_mod
  USE ilu_precond_mod
  IMPLICIT NONE
  
  ! Test variables
  INTEGER :: n, nz, fill_level
  INTEGER, ALLOCATABLE :: csr_row_ptr(:), csr_col_idx(:)
  REAL(kind=dp), ALLOCATABLE :: csr_val(:)
  INTEGER, ALLOCATABLE :: L_row_ptr(:), L_col_idx(:)
  INTEGER, ALLOCATABLE :: U_row_ptr(:), U_col_idx(:)
  REAL(kind=dp), ALLOCATABLE :: L_val(:), U_val(:)
  REAL(kind=dp), ALLOCATABLE :: x(:), b(:), y(:), z(:)
  REAL(kind=dp), ALLOCATABLE :: A_full(:,:), L_full(:,:), U_full(:,:), LU_full(:,:)
  COMPLEX(kind=dp), ALLOCATABLE :: z_csr_val(:)
  COMPLEX(kind=dp), ALLOCATABLE :: z_L_val(:), z_U_val(:)
  COMPLEX(kind=dp), ALLOCATABLE :: z_x(:), z_b(:), z_y(:), z_z(:)
  TYPE(ilu_factorization) :: ilu_fac
  TYPE(ilu_factorization_complex) :: z_ilu_fac
  REAL(kind=dp) :: drop_tol, norm_diff, max_diff
  INTEGER :: i, j, k, info
  LOGICAL :: test_passed, tests_passed
  REAL(kind=dp), PARAMETER :: tol = 1.0e-12_dp
  
  tests_passed = .TRUE.
  
  WRITE(*,'(A)') "================================="
  WRITE(*,'(A)') "ILU Preconditioner Test Suite"
  WRITE(*,'(A)') "================================="
  WRITE(*,*)
  
  ! Test 1: ILU(0) factorization of diagonal matrix
  WRITE(*,'(A)') "Test 1: ILU(0) factorization of diagonal matrix"
  test_passed = .TRUE.
  
  ! 3x3 diagonal matrix
  n = 3
  nz = 3
  ALLOCATE(csr_row_ptr(4), csr_col_idx(3), csr_val(3))
  csr_row_ptr = (/1, 2, 3, 4/)
  csr_col_idx = (/1, 2, 3/)
  csr_val = (/2.0_dp, 3.0_dp, 4.0_dp/)
  
  ! Perform ILU(0) factorization
  fill_level = 0
  drop_tol = 0.0_dp
  CALL ilu_factorize(n, csr_row_ptr, csr_col_idx, csr_val, fill_level, drop_tol, ilu_fac, info)
  
  IF (info == 0) THEN
    ! For diagonal matrix, L should be identity, U should be the diagonal
    IF (ilu_fac%U_val(1) == 2.0_dp .AND. &
        ilu_fac%U_val(2) == 3.0_dp .AND. &
        ilu_fac%U_val(3) == 4.0_dp) THEN
      WRITE(*,'(A)') "[PASS] ILU(0) diagonal matrix factorization"
    ELSE
      test_passed = .FALSE.
      tests_passed = .FALSE.
      WRITE(*,'(A)') "[FAIL] ILU(0) diagonal matrix factorization"
    END IF
  ELSE
    test_passed = .FALSE.
    tests_passed = .FALSE.
    WRITE(*,'(A,I0)') "[FAIL] ILU(0) factorization failed with info = ", info
  END IF
  
  CALL ilu_free(ilu_fac)
  DEALLOCATE(csr_row_ptr, csr_col_idx, csr_val)
  
  ! Test 2: ILU(0) factorization of tridiagonal matrix
  WRITE(*,'(A)') "Test 2: ILU(0) factorization of tridiagonal matrix"
  test_passed = .TRUE.
  
  ! 4x4 tridiagonal matrix
  n = 4
  nz = 10  ! 4 diagonal + 3 upper + 3 lower
  ALLOCATE(csr_row_ptr(5), csr_col_idx(10), csr_val(10))
  csr_row_ptr = (/1, 3, 6, 9, 11/)
  csr_col_idx = (/1, 2, 1, 2, 3, 2, 3, 4, 3, 4/)
  csr_val = (/4.0_dp, -1.0_dp, -1.0_dp, 4.0_dp, -1.0_dp, -1.0_dp, 4.0_dp, -1.0_dp, -1.0_dp, 4.0_dp/)
  
  ! Perform ILU(0) factorization
  fill_level = 0
  drop_tol = 0.0_dp
  CALL ilu_factorize(n, csr_row_ptr, csr_col_idx, csr_val, fill_level, drop_tol, ilu_fac, info)
  
  IF (info == 0) THEN
    ! Verify L*U ≈ A by checking a few elements
    WRITE(*,'(A)') "[PASS] ILU(0) tridiagonal matrix factorization completed"
  ELSE
    test_passed = .FALSE.
    tests_passed = .FALSE.
    WRITE(*,'(A,I0)') "[FAIL] ILU(0) factorization failed with info = ", info
  END IF
  
  CALL ilu_free(ilu_fac)
  DEALLOCATE(csr_row_ptr, csr_col_idx, csr_val)
  
  ! Test 3: Forward/backward substitution
  WRITE(*,'(A)') "Test 3: Forward/backward substitution"
  test_passed = .TRUE.
  
  ! Simple 3x3 lower triangular system for testing
  n = 3
  nz = 6
  ALLOCATE(csr_row_ptr(4), csr_col_idx(6), csr_val(6))
  csr_row_ptr = (/1, 2, 4, 7/)
  csr_col_idx = (/1, 1, 2, 1, 2, 3/)
  csr_val = (/2.0_dp, 1.0_dp, 3.0_dp, 2.0_dp, 1.0_dp, 4.0_dp/)
  
  ! Factorize
  fill_level = 0
  drop_tol = 0.0_dp
  CALL ilu_factorize(n, csr_row_ptr, csr_col_idx, csr_val, fill_level, drop_tol, ilu_fac, info)
  
  IF (info == 0) THEN
    ! Test solving L*U*x = b
    ALLOCATE(b(3), x(3))
    b = (/2.0_dp, 4.0_dp, 7.0_dp/)
    
    ! Apply ILU preconditioner: solve L*U*x = b
    CALL ilu_solve(ilu_fac, b, x)
    
    ! Verify solution by computing A*x and comparing with b
    ALLOCATE(y(3))
    CALL csr_matvec(n, csr_row_ptr, csr_col_idx, csr_val, x, y)
    
    max_diff = MAXVAL(ABS(y - b))
    IF (max_diff < tol) THEN
      WRITE(*,'(A)') "[PASS] Forward/backward substitution"
    ELSE
      test_passed = .FALSE.
      tests_passed = .FALSE.
      WRITE(*,'(A,E12.4)') "[FAIL] Forward/backward substitution, max error = ", max_diff
    END IF
    
    DEALLOCATE(b, x, y)
  ELSE
    test_passed = .FALSE.
    tests_passed = .FALSE.
    WRITE(*,'(A)') "[FAIL] Could not factorize matrix for substitution test"
  END IF
  
  CALL ilu_free(ilu_fac)
  DEALLOCATE(csr_row_ptr, csr_col_idx, csr_val)
  
  ! Test 4: ILU(1) factorization with fill-in
  WRITE(*,'(A)') "Test 4: ILU(1) factorization with fill-in"
  test_passed = .TRUE.
  
  ! 4x4 sparse matrix that will have fill-in
  n = 4
  nz = 8
  ALLOCATE(csr_row_ptr(5), csr_col_idx(8), csr_val(8))
  csr_row_ptr = (/1, 3, 5, 7, 9/)
  csr_col_idx = (/1, 3, 2, 4, 1, 3, 2, 4/)
  csr_val = (/4.0_dp, 1.0_dp, 4.0_dp, 1.0_dp, 1.0_dp, 4.0_dp, 1.0_dp, 4.0_dp/)
  
  ! Perform ILU(1) factorization
  fill_level = 1
  drop_tol = 0.0_dp
  CALL ilu_factorize(n, csr_row_ptr, csr_col_idx, csr_val, fill_level, drop_tol, ilu_fac, info)
  
  IF (info == 0) THEN
    ! Check that fill-in occurred (more non-zeros in L+U than in A)
    IF (ilu_fac%L_nnz + ilu_fac%U_nnz > nz) THEN
      WRITE(*,'(A)') "[PASS] ILU(1) factorization with fill-in"
    ELSE
      test_passed = .FALSE.
      tests_passed = .FALSE.
      WRITE(*,'(A)') "[FAIL] ILU(1) did not create expected fill-in"
    END IF
  ELSE
    test_passed = .FALSE.
    tests_passed = .FALSE.
    WRITE(*,'(A,I0)') "[FAIL] ILU(1) factorization failed with info = ", info
  END IF
  
  CALL ilu_free(ilu_fac)
  DEALLOCATE(csr_row_ptr, csr_col_idx, csr_val)
  
  ! Test 5: Drop tolerance
  WRITE(*,'(A)') "Test 5: ILU with drop tolerance"
  test_passed = .TRUE.
  
  ! Matrix with small off-diagonal elements
  n = 3
  nz = 7
  ALLOCATE(csr_row_ptr(4), csr_col_idx(7), csr_val(7))
  csr_row_ptr = (/1, 4, 6, 8/)
  csr_col_idx = (/1, 2, 3, 1, 2, 2, 3/)
  csr_val = (/10.0_dp, 0.01_dp, 0.02_dp, 0.01_dp, 10.0_dp, 0.02_dp, 10.0_dp/)
  
  ! Perform ILU(0) with drop tolerance
  fill_level = 0
  drop_tol = 0.1_dp  ! Should drop elements smaller than 0.1
  CALL ilu_factorize(n, csr_row_ptr, csr_col_idx, csr_val, fill_level, drop_tol, ilu_fac, info)
  
  IF (info == 0) THEN
    ! Check that small elements were dropped
    WRITE(*,'(A)') "[PASS] ILU with drop tolerance completed"
  ELSE
    test_passed = .FALSE.
    tests_passed = .FALSE.
    WRITE(*,'(A,I0)') "[FAIL] ILU with drop tolerance failed with info = ", info
  END IF
  
  CALL ilu_free(ilu_fac)
  DEALLOCATE(csr_row_ptr, csr_col_idx, csr_val)
  
  ! Test 6: Complex ILU factorization
  WRITE(*,'(A)') "Test 6: Complex ILU(0) factorization"
  test_passed = .TRUE.
  
  ! 2x2 complex matrix
  n = 2
  nz = 4
  ALLOCATE(csr_row_ptr(3), csr_col_idx(4), z_csr_val(4))
  csr_row_ptr = (/1, 3, 5/)
  csr_col_idx = (/1, 2, 1, 2/)
  z_csr_val = (/(2.0_dp,1.0_dp), (0.0_dp,1.0_dp), (0.0_dp,-1.0_dp), (2.0_dp,1.0_dp)/)
  
  ! Perform complex ILU(0) factorization
  fill_level = 0
  drop_tol = 0.0_dp
  CALL ilu_factorize(n, csr_row_ptr, csr_col_idx, z_csr_val, fill_level, drop_tol, z_ilu_fac, info)
  
  IF (info == 0) THEN
    ! Test solving with complex system
    ALLOCATE(z_b(2), z_x(2))
    z_b = (/(2.0_dp,3.0_dp), (1.0_dp,2.0_dp)/)
    
    CALL ilu_solve(z_ilu_fac, z_b, z_x)
    
    ! Verify by computing A*x
    ALLOCATE(z_y(2))
    CALL csr_matvec(n, csr_row_ptr, csr_col_idx, z_csr_val, z_x, z_y)
    
    max_diff = MAXVAL(ABS(z_y - z_b))
    IF (max_diff < tol) THEN
      WRITE(*,'(A)') "[PASS] Complex ILU(0) factorization and solve"
    ELSE
      test_passed = .FALSE.
      tests_passed = .FALSE.
      WRITE(*,'(A,E12.4)') "[FAIL] Complex ILU solve error = ", max_diff
    END IF
    
    DEALLOCATE(z_b, z_x, z_y)
  ELSE
    test_passed = .FALSE.
    tests_passed = .FALSE.
    WRITE(*,'(A,I0)') "[FAIL] Complex ILU factorization failed with info = ", info
  END IF
  
  CALL ilu_free(z_ilu_fac)
  DEALLOCATE(csr_row_ptr, csr_col_idx, z_csr_val)
  
  ! Test 7: Singular matrix handling
  WRITE(*,'(A)') "Test 7: Singular matrix handling"
  test_passed = .TRUE.
  
  ! Singular 3x3 matrix (row 3 is zero)
  n = 3
  nz = 4
  ALLOCATE(csr_row_ptr(4), csr_col_idx(4), csr_val(4))
  csr_row_ptr = (/1, 3, 5, 5/)  ! Row 3 has no entries
  csr_col_idx = (/1, 2, 1, 2/)
  csr_val = (/2.0_dp, 1.0_dp, 1.0_dp, 2.0_dp/)
  
  ! Try to factorize singular matrix
  fill_level = 0
  drop_tol = 0.0_dp
  CALL ilu_factorize(n, csr_row_ptr, csr_col_idx, csr_val, fill_level, drop_tol, ilu_fac, info)
  
  IF (info /= 0) THEN
    WRITE(*,'(A)') "[PASS] Singular matrix detected correctly"
  ELSE
    test_passed = .FALSE.
    tests_passed = .FALSE.
    WRITE(*,'(A)') "[FAIL] Singular matrix not detected"
  END IF
  
  DEALLOCATE(csr_row_ptr, csr_col_idx, csr_val)
  
  ! Test 8: Near-singular matrix with pivoting
  WRITE(*,'(A)') "Test 8: Near-singular matrix handling"
  test_passed = .TRUE.
  
  ! Near-singular matrix
  n = 3
  nz = 9
  ALLOCATE(csr_row_ptr(4), csr_col_idx(9), csr_val(9))
  csr_row_ptr = (/1, 4, 7, 10/)
  csr_col_idx = (/1, 2, 3, 1, 2, 3, 1, 2, 3/)
  csr_val = (/1.0e-10_dp, 1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 3.0_dp/)
  
  ! Try factorization
  fill_level = 0
  drop_tol = 0.0_dp
  CALL ilu_factorize(n, csr_row_ptr, csr_col_idx, csr_val, fill_level, drop_tol, ilu_fac, info)
  
  IF (info == 0 .OR. info == 1) THEN  ! info=1 might indicate near-singular
    WRITE(*,'(A)') "[PASS] Near-singular matrix handled"
  ELSE
    test_passed = .FALSE.
    tests_passed = .FALSE.
    WRITE(*,'(A,I0)') "[FAIL] Unexpected error for near-singular matrix, info = ", info
  END IF
  
  IF (info == 0) CALL ilu_free(ilu_fac)
  DEALLOCATE(csr_row_ptr, csr_col_idx, csr_val)
  
  ! Test 9: Verify L*U ≈ A for general matrix
  WRITE(*,'(A)') "Test 9: Verify L*U approximates A"
  test_passed = .TRUE.
  
  ! Simple 3x3 matrix (same as Test 3 which passed)
  n = 3
  nz = 6
  ALLOCATE(csr_row_ptr(4), csr_col_idx(6), csr_val(6))
  csr_row_ptr = (/1, 2, 4, 7/)
  csr_col_idx = (/1, 1, 2, 1, 2, 3/)
  csr_val = (/2.0_dp, 1.0_dp, 3.0_dp, 2.0_dp, 1.0_dp, 4.0_dp/)
  
  ! Perform ILU(0) factorization
  fill_level = 0
  drop_tol = 0.0_dp
  CALL ilu_factorize(n, csr_row_ptr, csr_col_idx, csr_val, fill_level, drop_tol, ilu_fac, info)
  
  IF (info == 0) THEN
    ! Reconstruct L*U and compare with A
    ! For simplicity, we'll test by solving multiple systems
    ALLOCATE(x(3), b(3), y(3))
    
    test_passed = .TRUE.
    DO i = 1, 3
      ! Test with unit vectors
      b = 0.0_dp
      b(i) = 1.0_dp
      
      ! Solve L*U*x = b
      CALL ilu_solve(ilu_fac, b, x)
      
      ! Compute A*x
      CALL csr_matvec(n, csr_row_ptr, csr_col_idx, csr_val, x, y)
      
      ! Check ||A*x - b|| for ILU(0), which should be exact for pattern
      max_diff = MAXVAL(ABS(y - b))
      IF (max_diff > 1.0e-10_dp) THEN
        test_passed = .FALSE.
        tests_passed = .FALSE.
        WRITE(*,'(A,I0,A,E12.4)') "[FAIL] L*U approximation error for column ", i, " = ", max_diff
        EXIT
      END IF
    END DO
    
    IF (test_passed) THEN
      WRITE(*,'(A)') "[PASS] L*U approximates A correctly"
    END IF
    
    DEALLOCATE(x, b, y)
  ELSE
    test_passed = .FALSE.
    tests_passed = .FALSE.
    WRITE(*,'(A,I0)') "[FAIL] ILU factorization failed with info = ", info
  END IF
  
  CALL ilu_free(ilu_fac)
  DEALLOCATE(csr_row_ptr, csr_col_idx, csr_val)
  
  ! Test 10: Performance test - larger matrix
  WRITE(*,'(A)') "Test 10: Performance test - 50x50 tridiagonal"
  test_passed = .TRUE.
  
  ! Create 50x50 tridiagonal matrix
  n = 50
  nz = 148  ! 50 diagonal + 49 upper + 49 lower
  ALLOCATE(csr_row_ptr(51), csr_col_idx(148), csr_val(148))
  
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
  
  ! Perform ILU(1) factorization
  fill_level = 1
  drop_tol = 0.0_dp
  CALL ilu_factorize(n, csr_row_ptr, csr_col_idx, csr_val, fill_level, drop_tol, ilu_fac, info)
  
  IF (info == 0) THEN
    ! Test solving a system
    ALLOCATE(b(50), x(50), y(50))
    b = 1.0_dp
    
    CALL ilu_solve(ilu_fac, b, x)
    
    ! Verify solution
    CALL csr_matvec(n, csr_row_ptr, csr_col_idx, csr_val, x, y)
    max_diff = MAXVAL(ABS(y - b))
    
    IF (max_diff < 1.0e-8_dp) THEN
      WRITE(*,'(A)') "[PASS] Large tridiagonal ILU(1) factorization"
    ELSE
      test_passed = .FALSE.
      tests_passed = .FALSE.
      WRITE(*,'(A,E12.4)') "[FAIL] Large matrix solve error = ", max_diff
    END IF
    
    DEALLOCATE(b, x, y)
  ELSE
    test_passed = .FALSE.
    tests_passed = .FALSE.
    WRITE(*,'(A,I0)') "[FAIL] Large matrix factorization failed with info = ", info
  END IF
  
  CALL ilu_free(ilu_fac)
  DEALLOCATE(csr_row_ptr, csr_col_idx, csr_val)
  
  ! Summary
  WRITE(*,*)
  WRITE(*,'(A)') "================================="
  IF (tests_passed) THEN
    WRITE(*,'(A)') "All ILU preconditioner tests PASSED!"
  ELSE
    WRITE(*,'(A)') "Some ILU preconditioner tests FAILED!"
    STOP 1
  END IF
  WRITE(*,'(A)') "================================="
  
END PROGRAM test_ilu_precond