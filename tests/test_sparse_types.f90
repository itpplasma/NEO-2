PROGRAM test_sparse_types
  ! Test for sparse_types_mod module
  ! NOTE: This test is currently minimal as the generic interfaces
  ! sparse_get_dimensions and sparse_deallocate are not yet implemented
  
  USE sparse_types_mod
  IMPLICIT NONE
  
  ! Test variables
  TYPE(sparse_matrix_csc_real) :: mat_csc_r
  TYPE(sparse_matrix_csc_complex) :: mat_csc_c
  TYPE(sparse_matrix_csr_real) :: mat_csr_r
  TYPE(sparse_matrix_csr_complex) :: mat_csr_c
  
  INTEGER :: nrow, ncol, nz
  LOGICAL :: test_passed
  
  test_passed = .TRUE.
  
  WRITE(*,'(A)') "================================="
  WRITE(*,'(A)') "Sparse Types Module Test"
  WRITE(*,'(A)') "================================="
  WRITE(*,*)
  
  ! Test 1: CSC real matrix type creation and basic operations
  WRITE(*,'(A)') "Test 1: CSC Real Matrix Type"
  CALL allocate_sparse(mat_csc_r, 4, 4, 7)
  
  IF (mat_csc_r%nrow == 4 .AND. mat_csc_r%ncol == 4 .AND. mat_csc_r%nz == 7) THEN
    WRITE(*,'(A)') "[PASS] CSC real allocation"
  ELSE
    WRITE(*,'(A)') "[FAIL] CSC real allocation"
    test_passed = .FALSE.
  END IF
  
  IF (ALLOCATED(mat_csc_r%irow) .AND. SIZE(mat_csc_r%irow) == 7) THEN
    WRITE(*,'(A)') "[PASS] CSC real irow allocation"
  ELSE
    WRITE(*,'(A)') "[FAIL] CSC real irow allocation"
    test_passed = .FALSE.
  END IF
  
  IF (ALLOCATED(mat_csc_r%pcol) .AND. SIZE(mat_csc_r%pcol) == 5) THEN
    WRITE(*,'(A)') "[PASS] CSC real pcol allocation"
  ELSE
    WRITE(*,'(A)') "[FAIL] CSC real pcol allocation"
    test_passed = .FALSE.
  END IF
  
  CALL deallocate_sparse(mat_csc_r)
  WRITE(*,*)
  
  ! Test 2: CSC complex matrix type creation
  WRITE(*,'(A)') "Test 2: CSC Complex Matrix Type"
  CALL allocate_sparse(mat_csc_c, 3, 3, 5)
  
  IF (mat_csc_c%nrow == 3 .AND. mat_csc_c%ncol == 3 .AND. mat_csc_c%nz == 5) THEN
    WRITE(*,'(A)') "[PASS] CSC complex allocation"
  ELSE
    WRITE(*,'(A)') "[FAIL] CSC complex allocation"
    test_passed = .FALSE.
  END IF
  
  CALL deallocate_sparse(mat_csc_c)
  WRITE(*,*)
  
  ! Test 3: CSR real matrix type creation
  WRITE(*,'(A)') "Test 3: CSR Real Matrix Type"
  CALL allocate_sparse(mat_csr_r, 3, 4, 6)
  
  IF (mat_csr_r%nrow == 3 .AND. mat_csr_r%ncol == 4 .AND. mat_csr_r%nz == 6) THEN
    WRITE(*,'(A)') "[PASS] CSR real allocation"
  ELSE
    WRITE(*,'(A)') "[FAIL] CSR real allocation"
    test_passed = .FALSE.
  END IF
  
  IF (ALLOCATED(mat_csr_r%jcol) .AND. SIZE(mat_csr_r%jcol) == 6) THEN
    WRITE(*,'(A)') "[PASS] CSR real jcol allocation"
  ELSE
    WRITE(*,'(A)') "[FAIL] CSR real jcol allocation"
    test_passed = .FALSE.
  END IF
  
  IF (ALLOCATED(mat_csr_r%prow) .AND. SIZE(mat_csr_r%prow) == 4) THEN
    WRITE(*,'(A)') "[PASS] CSR real prow allocation"
  ELSE
    WRITE(*,'(A)') "[FAIL] CSR real prow allocation"
    test_passed = .FALSE.
  END IF
  
  CALL deallocate_sparse(mat_csr_r)
  WRITE(*,*)
  
  ! Test 4: CSR complex matrix type creation
  WRITE(*,'(A)') "Test 4: CSR Complex Matrix Type"
  CALL allocate_sparse(mat_csr_c, 2, 2, 3)
  
  IF (mat_csr_c%nrow == 2 .AND. mat_csr_c%ncol == 2 .AND. mat_csr_c%nz == 3) THEN
    WRITE(*,'(A)') "[PASS] CSR complex allocation"
  ELSE
    WRITE(*,'(A)') "[FAIL] CSR complex allocation"
    test_passed = .FALSE.
  END IF
  
  CALL deallocate_sparse(mat_csr_c)
  
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
  
END PROGRAM test_sparse_types