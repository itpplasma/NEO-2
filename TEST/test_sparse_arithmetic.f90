PROGRAM test_sparse_arithmetic
  ! Test for sparse_arithmetic_mod module
  
  USE sparse_types_mod, ONLY: dp
  USE sparse_arithmetic_mod
  USE sparse_conversion_mod, ONLY: full2sparse, column_pointer2full
  IMPLICIT NONE
  
  ! Test variables
  INTEGER :: nrow, ncol, nz
  INTEGER, ALLOCATABLE :: irow(:), pcol(:), icol(:)
  REAL(kind=dp), ALLOCATABLE :: val(:), x(:), r(:), b(:)
  REAL(kind=dp), ALLOCATABLE :: x_2d(:,:), r_2d(:,:)
  REAL(kind=dp), ALLOCATABLE :: A_full(:,:)
  COMPLEX(kind=dp), ALLOCATABLE :: z_val(:), z_x(:), z_r(:)
  COMPLEX(kind=dp), ALLOCATABLE :: z_x_2d(:,:), z_r_2d(:,:)
  COMPLEX(kind=dp), ALLOCATABLE :: z_A_full(:,:)
  REAL(kind=dp) :: max_abs_err, max_rel_err
  INTEGER :: i, j
  LOGICAL :: test_passed
  REAL(kind=dp), PARAMETER :: tol = 1.0e-14_dp
  
  test_passed = .TRUE.
  
  WRITE(*,'(A)') "================================="
  WRITE(*,'(A)') "Sparse Arithmetic Module Test"
  WRITE(*,'(A)') "================================="
  WRITE(*,*)
  
  ! Test 1: Sparse matrix-vector multiplication (Real, 1D)
  WRITE(*,'(A)') "Test 1: Sparse matrix-vector multiplication (Real, 1D)"
  
  ! Create test matrix in CSC format
  nrow = 4
  ncol = 4
  nz = 7
  ALLOCATE(irow(nz), pcol(ncol+1), val(nz))
  
  pcol = (/1, 3, 4, 6, 8/)
  irow = (/1, 3, 2, 1, 3, 2, 4/)
  val = (/1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp, 7.0_dp/)
  
  ! Test vector
  ALLOCATE(x(ncol))
  x = (/1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp/)
  
  ! Compute matrix-vector product
  CALL sparse_matmul(nrow, ncol, irow, pcol, val, x, r)
  
  ! Check result (computed manually)
  IF (SIZE(r) == nrow .AND. &
      ABS(r(1) - 13.0_dp) < tol .AND. &  ! 1*1 + 4*3
      ABS(r(2) - 30.0_dp) < tol .AND. &  ! 3*2 + 6*4
      ABS(r(3) - 17.0_dp) < tol .AND. &  ! 2*1 + 5*3
      ABS(r(4) - 28.0_dp) < tol) THEN    ! 7*4
    WRITE(*,'(A)') "[PASS] Sparse matmul (Real, 1D)"
  ELSE
    WRITE(*,'(A)') "[FAIL] Sparse matmul (Real, 1D)"
    test_passed = .FALSE.
  END IF
  
  DEALLOCATE(x, r)
  
  ! Test 2: Sparse matrix-vector multiplication (Real, 2D)
  WRITE(*,'(A)') "Test 2: Sparse matrix-vector multiplication (Real, 2D)"
  
  ! Test with multiple vectors
  ALLOCATE(x_2d(ncol, 2))
  x_2d(:,1) = (/1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp/)
  x_2d(:,2) = (/0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp/)
  
  CALL sparse_matmul(nrow, ncol, irow, pcol, val, x_2d, r_2d)
  
  IF (SIZE(r_2d,1) == nrow .AND. SIZE(r_2d,2) == 2 .AND. &
      ABS(r_2d(1,1) - 13.0_dp) < tol .AND. &
      ABS(r_2d(2,2) - 9.0_dp) < tol) THEN  ! 3*1 + 6*1
    WRITE(*,'(A)') "[PASS] Sparse matmul (Real, 2D)"
  ELSE
    WRITE(*,'(A)') "[FAIL] Sparse matmul (Real, 2D)"
    test_passed = .FALSE.
  END IF
  
  DEALLOCATE(x_2d, r_2d, irow, pcol, val)
  
  ! Test 3: Sparse matrix-vector multiplication (Complex, 1D)
  WRITE(*,'(A)') "Test 3: Sparse matrix-vector multiplication (Complex, 1D)"
  
  nrow = 3
  ncol = 3
  nz = 5
  ALLOCATE(irow(nz), pcol(ncol+1), z_val(nz))
  
  pcol = (/1, 3, 4, 6/)
  irow = (/1, 2, 2, 1, 3/)
  z_val = (/(1.0_dp, 0.0_dp), (0.0_dp, 1.0_dp), &
            (2.0_dp, -1.0_dp), (3.0_dp, 2.0_dp), (4.0_dp, 0.0_dp)/)
  
  ALLOCATE(z_x(ncol))
  z_x = (/(1.0_dp, 0.0_dp), (0.0_dp, 1.0_dp), (1.0_dp, 1.0_dp)/)
  
  CALL sparse_matmul(nrow, ncol, irow, pcol, z_val, z_x, z_r)
  
  IF (SIZE(z_r) == nrow) THEN
    WRITE(*,'(A)') "[PASS] Sparse matmul (Complex, 1D)"
  ELSE
    WRITE(*,'(A)') "[FAIL] Sparse matmul (Complex, 1D)"
    test_passed = .FALSE.
  END IF
  
  DEALLOCATE(z_x, z_r, irow, pcol, z_val)
  
  ! Test 4: Full matrix interface (Real)
  WRITE(*,'(A)') "Test 4: Full matrix interface (Real)"
  
  ALLOCATE(A_full(3,3))
  A_full = 0.0_dp
  A_full(1,1) = 2.0_dp
  A_full(2,2) = 3.0_dp
  A_full(3,3) = 4.0_dp
  A_full(1,2) = 1.0_dp
  
  ALLOCATE(x(3))
  x = (/1.0_dp, 2.0_dp, 3.0_dp/)
  
  ! Convert to sparse for Test 5 to use
  CALL full2sparse(A_full, irow, pcol, val, nrow, ncol, nz)
  
  CALL sparse_matmul(A_full, x, r)
  
  IF (SIZE(r) == 3 .AND. &
      ABS(r(1) - 4.0_dp) < tol .AND. &  ! 2*1 + 1*2
      ABS(r(2) - 6.0_dp) < tol .AND. &  ! 3*2
      ABS(r(3) - 12.0_dp) < tol) THEN   ! 4*3
    WRITE(*,'(A)') "[PASS] Full matrix interface (Real)"
  ELSE
    WRITE(*,'(A)') "[FAIL] Full matrix interface (Real)"
    test_passed = .FALSE.
  END IF
  
  DEALLOCATE(A_full, x, r)
  
  ! Test 5: Sparse solver test
  WRITE(*,'(A)') "Test 5: Sparse solver test"
  
  ! Use the same matrix data from Test 4
  
  ALLOCATE(x(ncol), b(nrow))
  x = (/1.0_dp, 2.0_dp, 3.0_dp/)
  b = (/4.0_dp, 6.0_dp, 12.0_dp/)  ! A*x with known values
  
  CALL sparse_solver_test(nrow, ncol, irow, pcol, val, x, b, max_abs_err, max_rel_err)
  
  IF (max_abs_err >= 0.0_dp .AND. max_rel_err >= 0.0_dp) THEN
    WRITE(*,'(A)') "[PASS] Sparse solver test"
  ELSE
    WRITE(*,'(A)') "[FAIL] Sparse solver test"
    test_passed = .FALSE.
  END IF
  
  DEALLOCATE(x, b, irow, pcol, val)
  
  ! Test 6: Edge case - empty result
  WRITE(*,'(A)') "Test 6: Edge case - empty matrix"
  
  nrow = 2
  ncol = 2
  nz = 0
  ALLOCATE(irow(1), pcol(ncol+1), val(1))  ! Allocate minimal size
  pcol = (/1, 1, 1/)  ! No elements
  
  ALLOCATE(x(ncol))
  x = (/1.0_dp, 2.0_dp/)
  
  CALL sparse_matmul(nrow, ncol, irow(1:nz), pcol, val(1:nz), x, r)
  
  IF (SIZE(r) == nrow .AND. ALL(r == 0.0_dp)) THEN
    WRITE(*,'(A)') "[PASS] Empty matrix multiplication"
  ELSE
    WRITE(*,'(A)') "[FAIL] Empty matrix multiplication"
    test_passed = .FALSE.
  END IF
  
  DEALLOCATE(irow, pcol, val, x, r)
  
  ! Test 7: Complex solver test
  WRITE(*,'(A)') "Test 7: Complex solver test"
  
  nrow = 3
  ncol = 3
  nz = 5
  ALLOCATE(irow(nz), pcol(ncol+1), z_val(nz))
  
  pcol = (/1, 3, 4, 6/)
  irow = (/1, 2, 2, 1, 3/)
  z_val = (/(2.0_dp, 0.0_dp), (0.0_dp, 1.0_dp), &
            (3.0_dp, 0.0_dp), (0.0_dp, -1.0_dp), (4.0_dp, 0.0_dp)/)
  
  ALLOCATE(z_x(ncol), b(nrow))
  z_x = (/(1.0_dp, 0.0_dp), (0.0_dp, 1.0_dp), (1.0_dp, 1.0_dp)/)
  
  ! Compute expected b = A*x
  ALLOCATE(z_r(nrow))
  CALL sparse_matmul(nrow, ncol, irow, pcol, z_val, z_x, z_r)
  b = REAL(z_r)  ! For test purposes
  
  CALL sparse_solver_test(nrow, ncol, irow, pcol, z_val, z_x, CMPLX(b, 0.0_dp, dp), max_abs_err, max_rel_err)
  
  IF (max_abs_err >= 0.0_dp .AND. max_rel_err >= 0.0_dp) THEN
    WRITE(*,'(A)') "[PASS] Complex solver test"
  ELSE
    WRITE(*,'(A)') "[FAIL] Complex solver test"
    test_passed = .FALSE.
  END IF
  
  DEALLOCATE(z_x, z_r, b, irow, pcol, z_val)
  
  ! Test 8: 2D solver test
  WRITE(*,'(A)') "Test 8: 2D solver test"
  
  nrow = 3
  ncol = 3
  nz = 5
  ALLOCATE(irow(nz), pcol(ncol+1), val(nz))
  
  pcol = (/1, 3, 4, 6/)
  irow = (/1, 2, 2, 1, 3/)
  val = (/2.0_dp, 1.0_dp, 3.0_dp, 1.0_dp, 4.0_dp/)
  
  ALLOCATE(x_2d(ncol, 2))
  x_2d(:,1) = (/1.0_dp, 1.0_dp, 1.0_dp/)
  x_2d(:,2) = (/0.0_dp, 1.0_dp, 0.0_dp/)
  
  ! Compute expected b = A*x for both columns
  CALL sparse_matmul(nrow, ncol, irow, pcol, val, x_2d, r_2d)
  
  CALL sparse_solver_test(nrow, ncol, irow, pcol, val, x_2d, r_2d, max_abs_err, max_rel_err)
  
  IF (max_abs_err >= 0.0_dp .AND. max_rel_err >= 0.0_dp) THEN
    WRITE(*,'(A)') "[PASS] 2D solver test"
  ELSE
    WRITE(*,'(A)') "[FAIL] 2D solver test"
    test_passed = .FALSE.
  END IF
  
  DEALLOCATE(x_2d, r_2d, irow, pcol, val)
  
  ! Test 9: Full matrix solver test
  WRITE(*,'(A)') "Test 9: Full matrix solver test"
  
  ALLOCATE(A_full(3,3))
  A_full = 0.0_dp
  A_full(1,1) = 2.0_dp
  A_full(2,2) = 3.0_dp
  A_full(3,3) = 4.0_dp
  A_full(1,2) = 1.0_dp
  
  ALLOCATE(x(3), b(3))
  x = (/1.0_dp, 2.0_dp, 3.0_dp/)
  
  ! Compute b = A*x
  CALL sparse_matmul(A_full, x, r)
  b = r
  
  CALL sparse_solver_test(A_full, x, b, max_abs_err, max_rel_err)
  
  IF (max_abs_err >= 0.0_dp .AND. max_rel_err >= 0.0_dp) THEN
    WRITE(*,'(A)') "[PASS] Full matrix solver test"
  ELSE
    WRITE(*,'(A)') "[FAIL] Full matrix solver test"
    test_passed = .FALSE.
  END IF
  
  DEALLOCATE(A_full, x, b, r)
  
  ! Test 10: Complex 2D multiplication
  WRITE(*,'(A)') "Test 10: Complex 2D multiplication"
  
  nrow = 2
  ncol = 2
  nz = 3
  ALLOCATE(irow(nz), pcol(ncol+1), z_val(nz))
  
  pcol = (/1, 2, 4/)
  irow = (/1, 1, 2/)
  z_val = (/(1.0_dp, 1.0_dp), (2.0_dp, 0.0_dp), (3.0_dp, -1.0_dp)/)
  
  ALLOCATE(z_x_2d(ncol, 2))
  z_x_2d(:,1) = (/(1.0_dp, 0.0_dp), (0.0_dp, 1.0_dp)/)
  z_x_2d(:,2) = (/(1.0_dp, 1.0_dp), (1.0_dp, -1.0_dp)/)
  
  CALL sparse_matmul(nrow, ncol, irow, pcol, z_val, z_x_2d, z_r_2d)
  
  IF (SIZE(z_r_2d,1) == nrow .AND. SIZE(z_r_2d,2) == 2) THEN
    WRITE(*,'(A)') "[PASS] Complex 2D multiplication"
  ELSE
    WRITE(*,'(A)') "[FAIL] Complex 2D multiplication"
    test_passed = .FALSE.
  END IF
  
  DEALLOCATE(z_x_2d, z_r_2d, irow, pcol, z_val)
  
  ! Test 11: Complex full matrix interface
  WRITE(*,'(A)') "Test 11: Complex full matrix interface"
  
  ALLOCATE(z_A_full(2,2))
  z_A_full(1,1) = (1.0_dp, 0.0_dp)
  z_A_full(1,2) = (0.0_dp, 1.0_dp)
  z_A_full(2,1) = (0.0_dp, -1.0_dp)
  z_A_full(2,2) = (2.0_dp, 0.0_dp)
  
  ALLOCATE(z_x(2))
  z_x = (/(1.0_dp, 0.0_dp), (0.0_dp, 1.0_dp)/)
  
  CALL sparse_matmul(z_A_full, z_x, z_r)
  
  IF (SIZE(z_r) == 2) THEN
    WRITE(*,'(A)') "[PASS] Complex full matrix interface"
  ELSE
    WRITE(*,'(A)') "[FAIL] Complex full matrix interface"
    test_passed = .FALSE.
  END IF
  
  DEALLOCATE(z_A_full, z_x, z_r)
  
  ! Test 12: Complex 2D full matrix interface
  WRITE(*,'(A)') "Test 12: Complex 2D full matrix interface"
  
  ALLOCATE(z_A_full(2,2))
  z_A_full(1,1) = (2.0_dp, 0.0_dp)
  z_A_full(1,2) = (1.0_dp, 1.0_dp)
  z_A_full(2,1) = (1.0_dp, -1.0_dp)
  z_A_full(2,2) = (2.0_dp, 0.0_dp)
  
  ALLOCATE(z_x_2d(2, 2))
  z_x_2d(:,1) = (/(1.0_dp, 0.0_dp), (1.0_dp, 0.0_dp)/)
  z_x_2d(:,2) = (/(0.0_dp, 1.0_dp), (0.0_dp, -1.0_dp)/)
  
  CALL sparse_matmul(z_A_full, z_x_2d, z_r_2d)
  
  IF (SIZE(z_r_2d,1) == 2 .AND. SIZE(z_r_2d,2) == 2) THEN
    WRITE(*,'(A)') "[PASS] Complex 2D full matrix interface"
  ELSE
    WRITE(*,'(A)') "[FAIL] Complex 2D full matrix interface"
    test_passed = .FALSE.
  END IF
  
  DEALLOCATE(z_A_full, z_x_2d, z_r_2d)
  
  ! Test 13: Complex 2D solver test
  WRITE(*,'(A)') "Test 13: Complex 2D solver test"
  
  nrow = 2
  ncol = 2
  nz = 4
  ALLOCATE(irow(nz), pcol(ncol+1), z_val(nz))
  
  pcol = (/1, 3, 5/)
  irow = (/1, 2, 1, 2/)
  z_val = (/(2.0_dp, 0.0_dp), (0.0_dp, 1.0_dp), (1.0_dp, 0.0_dp), (2.0_dp, 0.0_dp)/)
  
  ALLOCATE(z_x_2d(ncol, 2))
  z_x_2d(:,1) = (/(1.0_dp, 0.0_dp), (1.0_dp, 0.0_dp)/)
  z_x_2d(:,2) = (/(0.0_dp, 1.0_dp), (0.0_dp, 1.0_dp)/)
  
  ! Compute expected b = A*x for both columns
  ALLOCATE(z_r_2d(nrow, 2))
  CALL sparse_matmul(nrow, ncol, irow, pcol, z_val, z_x_2d, z_r_2d)
  
  CALL sparse_solver_test(nrow, ncol, irow, pcol, z_val, z_x_2d, z_r_2d, max_abs_err, max_rel_err)
  
  IF (max_abs_err >= 0.0_dp .AND. max_rel_err >= 0.0_dp) THEN
    WRITE(*,'(A)') "[PASS] Complex 2D solver test"
  ELSE
    WRITE(*,'(A)') "[FAIL] Complex 2D solver test"
    test_passed = .FALSE.
  END IF
  
  DEALLOCATE(z_x_2d, z_r_2d, irow, pcol, z_val)
  
  ! Test 14: Complex full matrix solver test
  WRITE(*,'(A)') "Test 14: Complex full matrix solver test"
  
  ALLOCATE(z_A_full(2,2))
  z_A_full(1,1) = (3.0_dp, 0.0_dp)
  z_A_full(1,2) = (1.0_dp, 1.0_dp)
  z_A_full(2,1) = (1.0_dp, -1.0_dp)
  z_A_full(2,2) = (3.0_dp, 0.0_dp)
  
  ALLOCATE(z_x(2))
  z_x = (/(1.0_dp, 0.0_dp), (1.0_dp, 0.0_dp)/)
  
  ! Compute b = A*x
  ALLOCATE(z_r(2))
  CALL sparse_matmul(z_A_full, z_x, z_r)
  
  CALL sparse_solver_test(z_A_full, z_x, z_r, max_abs_err, max_rel_err)
  
  IF (max_abs_err >= 0.0_dp .AND. max_rel_err >= 0.0_dp) THEN
    WRITE(*,'(A)') "[PASS] Complex full matrix solver test"
  ELSE
    WRITE(*,'(A)') "[FAIL] Complex full matrix solver test"
    test_passed = .FALSE.
  END IF
  
  DEALLOCATE(z_A_full, z_x, z_r)
  
  ! Test 15: Full matrix 2D solver test
  WRITE(*,'(A)') "Test 15: Full matrix 2D solver test"
  
  ALLOCATE(A_full(3,3))
  A_full = 0.0_dp
  A_full(1,1) = 5.0_dp
  A_full(2,2) = 5.0_dp
  A_full(3,3) = 5.0_dp
  A_full(1,2) = 1.0_dp
  A_full(2,3) = 1.0_dp
  
  ALLOCATE(x_2d(3, 2))
  x_2d(:,1) = (/1.0_dp, 1.0_dp, 1.0_dp/)
  x_2d(:,2) = (/2.0_dp, 2.0_dp, 2.0_dp/)
  
  ! Compute b = A*x
  CALL sparse_matmul(A_full, x_2d, r_2d)
  
  CALL sparse_solver_test(A_full, x_2d, r_2d, max_abs_err, max_rel_err)
  
  IF (max_abs_err >= 0.0_dp .AND. max_rel_err >= 0.0_dp) THEN
    WRITE(*,'(A)') "[PASS] Full matrix 2D solver test"
  ELSE
    WRITE(*,'(A)') "[FAIL] Full matrix 2D solver test"
    test_passed = .FALSE.
  END IF
  
  DEALLOCATE(A_full, x_2d, r_2d)
  
  ! Test 16: Complex full matrix 2D solver test
  WRITE(*,'(A)') "Test 16: Complex full matrix 2D solver test"
  
  ALLOCATE(z_A_full(2,2))
  z_A_full(1,1) = (4.0_dp, 0.0_dp)
  z_A_full(1,2) = (0.0_dp, 2.0_dp)
  z_A_full(2,1) = (0.0_dp, -2.0_dp)
  z_A_full(2,2) = (4.0_dp, 0.0_dp)
  
  ALLOCATE(z_x_2d(2, 2))
  z_x_2d(:,1) = (/(1.0_dp, 0.0_dp), (1.0_dp, 0.0_dp)/)
  z_x_2d(:,2) = (/(0.0_dp, 1.0_dp), (0.0_dp, 1.0_dp)/)
  
  ! Compute b = A*x
  ALLOCATE(z_r_2d(2, 2))
  CALL sparse_matmul(z_A_full, z_x_2d, z_r_2d)
  
  CALL sparse_solver_test(z_A_full, z_x_2d, z_r_2d, max_abs_err, max_rel_err)
  
  IF (max_abs_err >= 0.0_dp .AND. max_rel_err >= 0.0_dp) THEN
    WRITE(*,'(A)') "[PASS] Complex full matrix 2D solver test"
  ELSE
    WRITE(*,'(A)') "[FAIL] Complex full matrix 2D solver test"
    test_passed = .FALSE.
  END IF
  
  DEALLOCATE(z_A_full, z_x_2d, z_r_2d)
  
  ! Summary
  WRITE(*,*)
  WRITE(*,'(A)') "================================="
  IF (test_passed) THEN
    WRITE(*,'(A)') "All tests PASSED!"
  ELSE
    WRITE(*,'(A)') "Some tests FAILED!"
    STOP 1
  END IF
  WRITE(*,'(A)') "================================="
  
END PROGRAM test_sparse_arithmetic