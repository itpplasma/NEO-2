PROGRAM test_sparse_utils
  ! Comprehensive tests for sparse_utils_mod module
  ! Tests CSC<->CSR conversions, matrix-vector multiplication, diagonal extraction
  
  USE sparse_types_mod, ONLY: dp, long
  USE sparse_utils_mod
  IMPLICIT NONE
  
  ! Test variables
  INTEGER :: nrow, ncol, nz
  INTEGER, ALLOCATABLE :: csc_col_ptr(:), csc_row_idx(:)
  INTEGER, ALLOCATABLE :: csr_row_ptr(:), csr_col_idx(:)
  REAL(kind=dp), ALLOCATABLE :: csc_val(:), csr_val(:)
  REAL(kind=dp), ALLOCATABLE :: x(:), y(:), y_expected(:)
  REAL(kind=dp), ALLOCATABLE :: diag(:), diag_expected(:)
  COMPLEX(kind=dp), ALLOCATABLE :: z_csc_val(:), z_csr_val(:)
  COMPLEX(kind=dp), ALLOCATABLE :: z_x(:), z_y(:), z_y_expected(:)
  COMPLEX(kind=dp), ALLOCATABLE :: z_diag(:), z_diag_expected(:)
  INTEGER, ALLOCATABLE :: csr2_row_ptr(:), csr2_col_idx(:)
  REAL(kind=dp), ALLOCATABLE :: csr2_val(:)
  INTEGER :: i, j
  LOGICAL :: test_passed, tests_passed
  REAL(kind=dp), PARAMETER :: tol = 1.0e-14_dp
  
  tests_passed = .TRUE.
  
  WRITE(*,'(A)') "================================="
  WRITE(*,'(A)') "Sparse Utils Module Test Suite"
  WRITE(*,'(A)') "================================="
  WRITE(*,*)
  
  ! Test 1: Basic CSC to CSR conversion (3x3 matrix)
  WRITE(*,'(A)') "Test 1: Basic CSC to CSR conversion"
  test_passed = .TRUE.
  
  ! Example matrix:
  ! [4 1 0]
  ! [1 3 0]
  ! [0 0 2]
  nrow = 3
  ncol = 3
  nz = 5
  
  ! CSC format (column-major)
  ALLOCATE(csc_col_ptr(4), csc_row_idx(5), csc_val(5))
  csc_col_ptr = (/1, 3, 5, 6/)  ! Column pointers (1-based)
  csc_row_idx = (/1, 2, 1, 2, 3/)  ! Row indices for each value
  csc_val = (/4.0_dp, 1.0_dp, 1.0_dp, 3.0_dp, 2.0_dp/)
  
  ! Expected CSR format (row-major)
  ALLOCATE(csr_row_ptr(4), csr_col_idx(5), csr_val(5))
  
  ! Convert CSC to CSR
  CALL csc_to_csr(nrow, ncol, nz, csc_col_ptr, csc_row_idx, csc_val, &
                  csr_row_ptr, csr_col_idx, csr_val)
  
  ! Check CSR row pointers (should be [1, 3, 5, 6])
  IF (csr_row_ptr(1) /= 1 .OR. csr_row_ptr(2) /= 3 .OR. &
      csr_row_ptr(3) /= 5 .OR. csr_row_ptr(4) /= 6) THEN
    test_passed = .FALSE.
    WRITE(*,'(A)') "[FAIL] CSR row pointers incorrect"
  END IF
  
  ! Check values and column indices
  ! Row 1: (1,4), (2,1)
  ! Row 2: (1,1), (2,3)
  ! Row 3: (3,2)
  
  IF (test_passed) THEN
    WRITE(*,'(A)') "[PASS] Basic CSC to CSR conversion"
  ELSE
    tests_passed = .FALSE.
  END IF
  
  DEALLOCATE(csc_col_ptr, csc_row_idx, csc_val)
  DEALLOCATE(csr_row_ptr, csr_col_idx, csr_val)
  
  ! Test 2: CSR to CSC conversion (round-trip)
  WRITE(*,'(A)') "Test 2: CSR to CSC conversion (round-trip)"
  test_passed = .TRUE.
  
  ! Start with CSR format
  nrow = 4
  ncol = 4
  nz = 7
  
  ALLOCATE(csr_row_ptr(5), csr_col_idx(7), csr_val(7))
  csr_row_ptr = (/1, 3, 5, 6, 8/)  ! Row pointers
  csr_col_idx = (/1, 3, 2, 4, 3, 1, 4/)  ! Column indices
  csr_val = (/2.0_dp, -1.0_dp, 3.0_dp, 1.0_dp, 4.0_dp, -2.0_dp, 5.0_dp/)
  
  ! Convert to CSC
  ALLOCATE(csc_col_ptr(5), csc_row_idx(7), csc_val(7))
  CALL csr_to_csc(nrow, ncol, nz, csr_row_ptr, csr_col_idx, csr_val, &
                  csc_col_ptr, csc_row_idx, csc_val)
  
  ! Convert back to CSR
  ALLOCATE(csr2_row_ptr(5), csr2_col_idx(7), csr2_val(7))
  
  CALL csc_to_csr(nrow, ncol, nz, csc_col_ptr, csc_row_idx, csc_val, &
                  csr2_row_ptr, csr2_col_idx, csr2_val)
  
  ! Check if we get back the original
  DO i = 1, 5
    IF (csr_row_ptr(i) /= csr2_row_ptr(i)) THEN
      test_passed = .FALSE.
      WRITE(*,'(A,I0,A)') "[FAIL] Row pointer mismatch at index ", i
    END IF
  END DO
  
  IF (test_passed) THEN
    WRITE(*,'(A)') "[PASS] CSR to CSC round-trip conversion"
  ELSE
    tests_passed = .FALSE.
  END IF
  
  DEALLOCATE(csr_row_ptr, csr_col_idx, csr_val)
  DEALLOCATE(csc_col_ptr, csc_row_idx, csc_val)
  DEALLOCATE(csr2_row_ptr, csr2_col_idx, csr2_val)
  
  ! Test 3: CSR matrix-vector multiplication
  WRITE(*,'(A)') "Test 3: CSR matrix-vector multiplication"
  test_passed = .TRUE.
  
  ! Matrix A (3x3):
  ! [4 1 0]
  ! [1 3 0]
  ! [0 0 2]
  nrow = 3
  ncol = 3
  nz = 5
  
  ALLOCATE(csr_row_ptr(4), csr_col_idx(5), csr_val(5))
  csr_row_ptr = (/1, 3, 5, 6/)
  csr_col_idx = (/1, 2, 1, 2, 3/)
  csr_val = (/4.0_dp, 1.0_dp, 1.0_dp, 3.0_dp, 2.0_dp/)
  
  ! Vector x
  ALLOCATE(x(3))
  x = (/1.0_dp, 2.0_dp, 3.0_dp/)
  
  ! Expected y = A*x = [6, 7, 6]
  ALLOCATE(y(3), y_expected(3))
  y_expected = (/6.0_dp, 7.0_dp, 6.0_dp/)
  
  ! Compute y = A*x
  CALL csr_matvec(nrow, csr_row_ptr, csr_col_idx, csr_val, x, y)
  
  ! Check result
  DO i = 1, nrow
    IF (ABS(y(i) - y_expected(i)) > tol) THEN
      test_passed = .FALSE.
      WRITE(*,'(A,I0,A,F10.6,A,F10.6)') "[FAIL] y(", i, ") = ", y(i), &
        " expected ", y_expected(i)
    END IF
  END DO
  
  IF (test_passed) THEN
    WRITE(*,'(A)') "[PASS] CSR matrix-vector multiplication"
  ELSE
    tests_passed = .FALSE.
  END IF
  
  DEALLOCATE(csr_row_ptr, csr_col_idx, csr_val)
  DEALLOCATE(x, y, y_expected)
  
  ! Test 4: Diagonal extraction from CSR
  WRITE(*,'(A)') "Test 4: Diagonal extraction from CSR"
  test_passed = .TRUE.
  
  ! Matrix with diagonal [4, 3, 2]
  nrow = 3
  ncol = 3
  nz = 5
  
  ALLOCATE(csr_row_ptr(4), csr_col_idx(5), csr_val(5))
  csr_row_ptr = (/1, 3, 5, 6/)
  csr_col_idx = (/1, 2, 1, 2, 3/)
  csr_val = (/4.0_dp, 1.0_dp, 1.0_dp, 3.0_dp, 2.0_dp/)
  
  ALLOCATE(diag(3), diag_expected(3))
  diag_expected = (/4.0_dp, 3.0_dp, 2.0_dp/)
  
  ! Extract diagonal
  CALL csr_extract_diagonal(nrow, csr_row_ptr, csr_col_idx, csr_val, diag)
  
  ! Check diagonal
  DO i = 1, nrow
    IF (ABS(diag(i) - diag_expected(i)) > tol) THEN
      test_passed = .FALSE.
      WRITE(*,'(A,I0,A,F10.6,A,F10.6)') "[FAIL] diag(", i, ") = ", diag(i), &
        " expected ", diag_expected(i)
    END IF
  END DO
  
  IF (test_passed) THEN
    WRITE(*,'(A)') "[PASS] Diagonal extraction from CSR"
  ELSE
    tests_passed = .FALSE.
  END IF
  
  DEALLOCATE(csr_row_ptr, csr_col_idx, csr_val)
  DEALLOCATE(diag, diag_expected)
  
  ! Test 5: Complex CSC to CSR conversion
  WRITE(*,'(A)') "Test 5: Complex CSC to CSR conversion"
  test_passed = .TRUE.
  
  ! Complex 2x2 matrix:
  ! [(2,1)  (0,3)]
  ! [(0,-2) (4,0)]
  nrow = 2
  ncol = 2
  nz = 4
  
  ALLOCATE(csc_col_ptr(3), csc_row_idx(4), z_csc_val(4))
  csc_col_ptr = (/1, 3, 5/)
  csc_row_idx = (/1, 2, 1, 2/)
  z_csc_val = (/(2.0_dp,1.0_dp), (0.0_dp,-2.0_dp), (0.0_dp,3.0_dp), (4.0_dp,0.0_dp)/)
  
  ALLOCATE(csr_row_ptr(3), csr_col_idx(4), z_csr_val(4))
  
  ! Convert complex CSC to CSR
  CALL csc_to_csr(nrow, ncol, nz, csc_col_ptr, csc_row_idx, z_csc_val, &
                  csr_row_ptr, csr_col_idx, z_csr_val)
  
  ! Verify structure
  IF (csr_row_ptr(1) == 1 .AND. csr_row_ptr(2) == 3 .AND. csr_row_ptr(3) == 5) THEN
    WRITE(*,'(A)') "[PASS] Complex CSC to CSR conversion"
  ELSE
    test_passed = .FALSE.
    tests_passed = .FALSE.
    WRITE(*,'(A)') "[FAIL] Complex CSC to CSR conversion"
  END IF
  
  DEALLOCATE(csc_col_ptr, csc_row_idx, z_csc_val)
  DEALLOCATE(csr_row_ptr, csr_col_idx, z_csr_val)
  
  ! Test 6: Complex CSR matrix-vector multiplication
  WRITE(*,'(A)') "Test 6: Complex CSR matrix-vector multiplication"
  test_passed = .TRUE.
  
  ! Complex matrix-vector product
  nrow = 2
  ncol = 2
  nz = 4
  
  ALLOCATE(csr_row_ptr(3), csr_col_idx(4), z_csr_val(4))
  csr_row_ptr = (/1, 3, 5/)
  csr_col_idx = (/1, 2, 1, 2/)
  z_csr_val = (/(2.0_dp,1.0_dp), (0.0_dp,3.0_dp), (0.0_dp,-2.0_dp), (4.0_dp,0.0_dp)/)
  
  ALLOCATE(z_x(2), z_y(2), z_y_expected(2))
  z_x = (/(1.0_dp,0.0_dp), (0.0_dp,1.0_dp)/)
  ! Expected: Row 1: (2,1)*(1,0) + (0,3)*(0,i) = (2,1) + (0,3)*i = (2,1) + (-3,0) = (-1,1)
  !           Row 2: (0,-2)*(1,0) + (4,0)*(0,i) = (0,-2) + (4,0)*i = (0,-2) + (0,4) = (0,2)
  z_y_expected = (/(-1.0_dp,1.0_dp), (0.0_dp,2.0_dp)/)
  
  CALL csr_matvec(nrow, csr_row_ptr, csr_col_idx, z_csr_val, z_x, z_y)
  
  DO i = 1, nrow
    IF (ABS(z_y(i) - z_y_expected(i)) > tol) THEN
      test_passed = .FALSE.
      WRITE(*,'(A,I0,A)') "[FAIL] Complex matrix-vector product at row ", i
    END IF
  END DO
  
  IF (test_passed) THEN
    WRITE(*,'(A)') "[PASS] Complex CSR matrix-vector multiplication"
  ELSE
    tests_passed = .FALSE.
  END IF
  
  DEALLOCATE(csr_row_ptr, csr_col_idx, z_csr_val)
  DEALLOCATE(z_x, z_y, z_y_expected)
  
  ! Test 7: Edge case - empty matrix
  WRITE(*,'(A)') "Test 7: Edge case - empty matrix"
  test_passed = .TRUE.
  
  nrow = 0
  ncol = 0
  nz = 0
  
  ALLOCATE(csr_row_ptr(1), csr_col_idx(0), csr_val(0))
  csr_row_ptr(1) = 1
  
  ALLOCATE(csc_col_ptr(1), csc_row_idx(0), csc_val(0))
  
  CALL csr_to_csc(nrow, ncol, nz, csr_row_ptr, csr_col_idx, csr_val, &
                  csc_col_ptr, csc_row_idx, csc_val)
  
  IF (csc_col_ptr(1) == 1) THEN
    WRITE(*,'(A)') "[PASS] Empty matrix handling"
  ELSE
    test_passed = .FALSE.
    tests_passed = .FALSE.
    WRITE(*,'(A)') "[FAIL] Empty matrix handling"
  END IF
  
  DEALLOCATE(csr_row_ptr, csr_col_idx, csr_val)
  DEALLOCATE(csc_col_ptr, csc_row_idx, csc_val)
  
  ! Test 8: Rectangular matrix CSC to CSR
  WRITE(*,'(A)') "Test 8: Rectangular matrix CSC to CSR"
  test_passed = .TRUE.
  
  ! 2x3 matrix:
  ! [1 0 2]
  ! [0 3 0]
  nrow = 2
  ncol = 3
  nz = 3
  
  ALLOCATE(csc_col_ptr(4), csc_row_idx(3), csc_val(3))
  csc_col_ptr = (/1, 2, 3, 4/)
  csc_row_idx = (/1, 2, 1/)
  csc_val = (/1.0_dp, 3.0_dp, 2.0_dp/)
  
  ALLOCATE(csr_row_ptr(3), csr_col_idx(3), csr_val(3))
  
  CALL csc_to_csr(nrow, ncol, nz, csc_col_ptr, csc_row_idx, csc_val, &
                  csr_row_ptr, csr_col_idx, csr_val)
  
  ! Expected: row_ptr = [1, 3, 4], col_idx = [1, 3, 2], val = [1, 2, 3]
  IF (csr_row_ptr(1) == 1 .AND. csr_row_ptr(2) == 3 .AND. csr_row_ptr(3) == 4) THEN
    WRITE(*,'(A)') "[PASS] Rectangular matrix conversion"
  ELSE
    test_passed = .FALSE.
    tests_passed = .FALSE.
    WRITE(*,'(A)') "[FAIL] Rectangular matrix conversion"
  END IF
  
  DEALLOCATE(csc_col_ptr, csc_row_idx, csc_val)
  DEALLOCATE(csr_row_ptr, csr_col_idx, csr_val)
  
  ! Test 9: Diagonal extraction from non-square matrix
  WRITE(*,'(A)') "Test 9: Diagonal extraction from non-square matrix"
  test_passed = .TRUE.
  
  ! 3x4 matrix - diagonal should be min(nrow,ncol) = 3
  nrow = 3
  ncol = 4
  nz = 5
  
  ALLOCATE(csr_row_ptr(4), csr_col_idx(5), csr_val(5))
  csr_row_ptr = (/1, 3, 4, 6/)
  csr_col_idx = (/1, 3, 2, 3, 4/)
  csr_val = (/2.0_dp, 1.0_dp, 3.0_dp, 4.0_dp, 5.0_dp/)
  
  ALLOCATE(diag(3), diag_expected(3))
  diag_expected = (/2.0_dp, 3.0_dp, 4.0_dp/)
  
  CALL csr_extract_diagonal(nrow, csr_row_ptr, csr_col_idx, csr_val, diag)
  
  DO i = 1, MIN(nrow, ncol)
    IF (ABS(diag(i) - diag_expected(i)) > tol) THEN
      test_passed = .FALSE.
      WRITE(*,'(A,I0,A)') "[FAIL] Non-square diagonal at position ", i
    END IF
  END DO
  
  IF (test_passed) THEN
    WRITE(*,'(A)') "[PASS] Non-square matrix diagonal extraction"
  ELSE
    tests_passed = .FALSE.
  END IF
  
  DEALLOCATE(csr_row_ptr, csr_col_idx, csr_val)
  DEALLOCATE(diag, diag_expected)
  
  ! Test 10: Performance test - larger matrix
  WRITE(*,'(A)') "Test 10: Performance test - 100x100 tridiagonal"
  test_passed = .TRUE.
  
  ! Create 100x100 tridiagonal matrix
  nrow = 100
  ncol = 100
  nz = 298  ! 100 diagonal + 99 upper + 99 lower
  
  ALLOCATE(csr_row_ptr(101), csr_col_idx(298), csr_val(298))
  
  ! Build tridiagonal structure
  j = 1
  DO i = 1, nrow
    csr_row_ptr(i) = j
    IF (i > 1) THEN
      csr_col_idx(j) = i-1
      csr_val(j) = -1.0_dp
      j = j + 1
    END IF
    csr_col_idx(j) = i
    csr_val(j) = 2.0_dp
    j = j + 1
    IF (i < ncol) THEN
      csr_col_idx(j) = i+1
      csr_val(j) = -1.0_dp
      j = j + 1
    END IF
  END DO
  csr_row_ptr(nrow+1) = j
  
  ! Test matrix-vector product
  ALLOCATE(x(100), y(100))
  x = 1.0_dp
  
  CALL csr_matvec(nrow, csr_row_ptr, csr_col_idx, csr_val, x, y)
  
  ! Check first, middle and last elements
  IF (ABS(y(1) - 1.0_dp) < tol .AND. &
      ABS(y(50) - 0.0_dp) < tol .AND. &
      ABS(y(100) - 1.0_dp) < tol) THEN
    WRITE(*,'(A)') "[PASS] Large tridiagonal matrix operations"
  ELSE
    test_passed = .FALSE.
    tests_passed = .FALSE.
    WRITE(*,'(A)') "[FAIL] Large tridiagonal matrix operations"
  END IF
  
  DEALLOCATE(csr_row_ptr, csr_col_idx, csr_val)
  DEALLOCATE(x, y)
  
  ! Summary
  WRITE(*,*)
  WRITE(*,'(A)') "================================="
  IF (tests_passed) THEN
    WRITE(*,'(A)') "All sparse utils tests PASSED!"
  ELSE
    WRITE(*,'(A)') "Some sparse utils tests FAILED!"
    STOP 1
  END IF
  WRITE(*,'(A)') "================================="
  
END PROGRAM test_sparse_utils