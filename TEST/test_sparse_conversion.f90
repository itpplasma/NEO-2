PROGRAM test_sparse_conversion
  ! Test for sparse_conversion_mod module
  
  USE sparse_types_mod, ONLY: dp
  USE sparse_conversion_mod
  IMPLICIT NONE
  
  ! Test variables
  INTEGER :: nrow, ncol, nz
  INTEGER, ALLOCATABLE :: irow(:), pcol(:), icol(:), pcol_reconstructed(:)
  REAL(kind=dp), ALLOCATABLE :: val(:), A_full(:,:), A_reconstructed(:,:)
  COMPLEX(kind=dp), ALLOCATABLE :: z_val(:), z_A_full(:,:), z_A_reconstructed(:,:)
  INTEGER :: i, j, nz_out
  LOGICAL :: test_passed
  REAL(kind=dp), PARAMETER :: tol = 1.0e-14_dp
  
  test_passed = .TRUE.
  
  WRITE(*,'(A)') "================================="
  WRITE(*,'(A)') "Sparse Conversion Module Test"
  WRITE(*,'(A)') "================================="
  WRITE(*,*)
  
  ! Test 1: Column pointer to full conversion
  WRITE(*,'(A)') "Test 1: Column pointer to full conversion"
  ALLOCATE(pcol(5))
  pcol = (/1, 3, 4, 6, 8/)
  CALL column_pointer2full(pcol, icol)
  
  IF (SIZE(icol) == 7 .AND. &
      icol(1) == 1 .AND. icol(2) == 1 .AND. &
      icol(3) == 2 .AND. icol(4) == 3 .AND. &
      icol(5) == 3 .AND. icol(6) == 4 .AND. &
      icol(7) == 4) THEN
    WRITE(*,'(A)') "[PASS] Column pointer to full"
  ELSE
    WRITE(*,'(A)') "[FAIL] Column pointer to full"
    test_passed = .FALSE.
  END IF
  
  ! Test 2: Column full to pointer conversion
  WRITE(*,'(A)') "Test 2: Column full to pointer conversion"
  CALL column_full2pointer(icol, pcol_reconstructed)
  
  IF (ALL(pcol == pcol_reconstructed)) THEN
    WRITE(*,'(A)') "[PASS] Column full to pointer"
  ELSE
    WRITE(*,'(A)') "[FAIL] Column full to pointer"
    test_passed = .FALSE.
  END IF
  
  DEALLOCATE(pcol, icol, pcol_reconstructed)
  
  ! Test 3: Sparse to full conversion (Real)
  WRITE(*,'(A)') "Test 3: Sparse to full conversion (Real)"
  nrow = 4
  ncol = 4
  nz = 7
  ALLOCATE(irow(nz), pcol(ncol+1), val(nz))
  
  pcol = (/1, 3, 4, 6, 8/)
  irow = (/1, 3, 2, 1, 3, 2, 4/)
  val = (/1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp, 7.0_dp/)
  
  CALL sparse2full(irow, pcol, val, nrow, ncol, A_full)
  
  ! Check specific values
  IF (ABS(A_full(1,1) - 1.0_dp) < tol .AND. &
      ABS(A_full(3,1) - 2.0_dp) < tol .AND. &
      ABS(A_full(2,2) - 3.0_dp) < tol .AND. &
      ABS(A_full(1,3) - 4.0_dp) < tol .AND. &
      ABS(A_full(3,3) - 5.0_dp) < tol .AND. &
      ABS(A_full(2,4) - 6.0_dp) < tol .AND. &
      ABS(A_full(4,4) - 7.0_dp) < tol .AND. &
      ABS(A_full(4,1)) < tol) THEN
    WRITE(*,'(A)') "[PASS] Sparse to full (Real)"
  ELSE
    WRITE(*,'(A)') "[FAIL] Sparse to full (Real)"
    test_passed = .FALSE.
  END IF
  
  DEALLOCATE(irow, pcol, val)
  
  ! Test 4: Full to sparse conversion (Real)
  WRITE(*,'(A)') "Test 4: Full to sparse conversion (Real)"
  CALL full2sparse(A_full, irow, pcol, val, nrow, ncol, nz_out)
  
  IF (nz_out == 7) THEN
    WRITE(*,'(A)') "[PASS] Full to sparse nonzero count"
  ELSE
    WRITE(*,'(A)') "[FAIL] Full to sparse nonzero count"
    test_passed = .FALSE.
  END IF
  
  ! Convert back and check
  CALL sparse2full(irow, pcol, val, nrow, ncol, A_reconstructed)
  
  IF (MAXVAL(ABS(A_full - A_reconstructed)) < tol) THEN
    WRITE(*,'(A)') "[PASS] Full to sparse roundtrip"
  ELSE
    WRITE(*,'(A)') "[FAIL] Full to sparse roundtrip"
    test_passed = .FALSE.
  END IF
  
  DEALLOCATE(A_full, A_reconstructed, irow, pcol, val)
  
  ! Test 5: Sparse to full conversion (Complex)
  WRITE(*,'(A)') "Test 5: Sparse to full conversion (Complex)"
  nrow = 3
  ncol = 3
  nz = 5
  ALLOCATE(irow(nz), pcol(ncol+1), z_val(nz))
  
  pcol = (/1, 3, 4, 6/)
  irow = (/1, 2, 2, 1, 3/)
  z_val = (/(1.0_dp, 0.0_dp), (0.0_dp, 1.0_dp), &
            (2.0_dp, -1.0_dp), (3.0_dp, 2.0_dp), (4.0_dp, 0.0_dp)/)
  
  CALL sparse2full(irow, pcol, z_val, nrow, ncol, z_A_full)
  
  IF (ABS(z_A_full(1,1) - CMPLX(1.0_dp, 0.0_dp, dp)) < tol .AND. &
      ABS(z_A_full(2,1) - CMPLX(0.0_dp, 1.0_dp, dp)) < tol .AND. &
      ABS(z_A_full(2,2) - CMPLX(2.0_dp, -1.0_dp, dp)) < tol .AND. &
      ABS(z_A_full(1,3) - CMPLX(3.0_dp, 2.0_dp, dp)) < tol .AND. &
      ABS(z_A_full(3,3) - CMPLX(4.0_dp, 0.0_dp, dp)) < tol) THEN
    WRITE(*,'(A)') "[PASS] Sparse to full (Complex)"
  ELSE
    WRITE(*,'(A)') "[FAIL] Sparse to full (Complex)"
    test_passed = .FALSE.
  END IF
  
  DEALLOCATE(irow, pcol, z_val)
  
  ! Test 6: Full to sparse conversion (Complex)
  WRITE(*,'(A)') "Test 6: Full to sparse conversion (Complex)"
  CALL full2sparse(z_A_full, irow, pcol, z_val, nrow, ncol, nz_out)
  
  IF (nz_out == 5) THEN
    WRITE(*,'(A)') "[PASS] Full to sparse (Complex) nonzero count"
  ELSE
    WRITE(*,'(A)') "[FAIL] Full to sparse (Complex) nonzero count"
    test_passed = .FALSE.
  END IF
  
  ! Convert back and check
  CALL sparse2full(irow, pcol, z_val, nrow, ncol, z_A_reconstructed)
  
  IF (MAXVAL(ABS(z_A_full - z_A_reconstructed)) < tol) THEN
    WRITE(*,'(A)') "[PASS] Full to sparse (Complex) roundtrip"
  ELSE
    WRITE(*,'(A)') "[FAIL] Full to sparse (Complex) roundtrip"
    test_passed = .FALSE.
  END IF
  
  DEALLOCATE(z_A_full, z_A_reconstructed, irow, pcol, z_val)
  
  ! Test 7: Edge case - empty matrix
  WRITE(*,'(A)') "Test 7: Edge case - empty matrix"
  ALLOCATE(A_full(3,3))
  A_full = 0.0_dp
  
  CALL full2sparse(A_full, irow, pcol, val, nrow, ncol, nz_out)
  
  IF (nz_out == 0 .AND. SIZE(irow) == 0 .AND. SIZE(val) == 0) THEN
    WRITE(*,'(A)') "[PASS] Empty matrix conversion"
  ELSE
    WRITE(*,'(A)') "[FAIL] Empty matrix conversion"
    test_passed = .FALSE.
  END IF
  
  DEALLOCATE(A_full)
  
  ! Test 8: Column indices as input
  WRITE(*,'(A)') "Test 8: Sparse to full with column indices"
  nrow = 3
  ncol = 3
  nz = 4
  IF (ALLOCATED(irow)) DEALLOCATE(irow)
  IF (ALLOCATED(icol)) DEALLOCATE(icol)
  IF (ALLOCATED(val)) DEALLOCATE(val)
  ALLOCATE(irow(nz), icol(nz), val(nz))
  
  irow = (/1, 2, 2, 3/)
  icol = (/1, 2, 3, 3/)
  val = (/1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp/)
  
  CALL sparse2full(irow, icol, val, nrow, ncol, A_full)
  
  IF (ABS(A_full(1,1) - 1.0_dp) < tol .AND. &
      ABS(A_full(2,2) - 2.0_dp) < tol .AND. &
      ABS(A_full(2,3) - 3.0_dp) < tol .AND. &
      ABS(A_full(3,3) - 4.0_dp) < tol) THEN
    WRITE(*,'(A)') "[PASS] Sparse to full with column indices"
  ELSE
    WRITE(*,'(A)') "[FAIL] Sparse to full with column indices"
    test_passed = .FALSE.
  END IF
  
  DEALLOCATE(irow, icol, val, A_full)
  
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
  
END PROGRAM test_sparse_conversion