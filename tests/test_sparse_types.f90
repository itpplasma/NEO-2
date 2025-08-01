PROGRAM test_sparse_types
  ! Test for sparse_types_mod module
  
  USE sparse_types_mod
  IMPLICIT NONE
  
  ! Test variables
  TYPE(sparse_matrix_csc_real) :: mat_csc_r
  TYPE(sparse_matrix_csc_complex) :: mat_csc_c
  TYPE(sparse_matrix_coo_real) :: mat_coo_r
  TYPE(sparse_matrix_coo_complex) :: mat_coo_c
  TYPE(sparse_matrix_csr_real) :: mat_csr_r
  TYPE(sparse_matrix_csr_complex) :: mat_csr_c
  
  INTEGER :: nrow, ncol, nz
  INTEGER :: i
  LOGICAL :: test_passed
  
  test_passed = .TRUE.
  
  WRITE(*,'(A)') "================================="
  WRITE(*,'(A)') "Sparse Types Module Test"
  WRITE(*,'(A)') "================================="
  WRITE(*,*)
  
  ! Test 1: Create and initialize CSC real matrix
  WRITE(*,'(A)') "Test 1: CSC Real Matrix"
  mat_csc_r%nrow = 4
  mat_csc_r%ncol = 4
  mat_csc_r%nz = 7
  ALLOCATE(mat_csc_r%irow(mat_csc_r%nz))
  ALLOCATE(mat_csc_r%pcol(mat_csc_r%ncol+1))
  ALLOCATE(mat_csc_r%val(mat_csc_r%nz))
  
  mat_csc_r%pcol = (/1, 3, 4, 6, 8/)
  mat_csc_r%irow = (/1, 3, 2, 1, 3, 2, 4/)
  mat_csc_r%val = (/1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp, 7.0_dp/)
  
  ! Test get_dimensions
  CALL sparse_get_dimensions(mat_csc_r, nrow, ncol, nz)
  IF (nrow == 4 .AND. ncol == 4 .AND. nz == 7) THEN
    WRITE(*,'(A)') "[PASS] CSC real get_dimensions"
  ELSE
    WRITE(*,'(A)') "[FAIL] CSC real get_dimensions"
    test_passed = .FALSE.
  END IF
  
  ! Test deallocate
  CALL sparse_deallocate(mat_csc_r)
  IF (.NOT. ALLOCATED(mat_csc_r%irow) .AND. &
      .NOT. ALLOCATED(mat_csc_r%pcol) .AND. &
      .NOT. ALLOCATED(mat_csc_r%val) .AND. &
      mat_csc_r%nrow == 0) THEN
    WRITE(*,'(A)') "[PASS] CSC real deallocate"
  ELSE
    WRITE(*,'(A)') "[FAIL] CSC real deallocate"
    test_passed = .FALSE.
  END IF
  
  ! Test 2: Create and initialize CSC complex matrix
  WRITE(*,'(A)') "Test 2: CSC Complex Matrix"
  mat_csc_c%nrow = 3
  mat_csc_c%ncol = 3
  mat_csc_c%nz = 5
  ALLOCATE(mat_csc_c%irow(mat_csc_c%nz))
  ALLOCATE(mat_csc_c%pcol(mat_csc_c%ncol+1))
  ALLOCATE(mat_csc_c%val(mat_csc_c%nz))
  
  mat_csc_c%pcol = (/1, 3, 4, 6/)
  mat_csc_c%irow = (/1, 2, 2, 1, 3/)
  mat_csc_c%val = (/(1.0_dp, 0.0_dp), (0.0_dp, 1.0_dp), &
                    (2.0_dp, -1.0_dp), (3.0_dp, 2.0_dp), (4.0_dp, 0.0_dp)/)
  
  CALL sparse_get_dimensions(mat_csc_c, nrow, ncol, nz)
  IF (nrow == 3 .AND. ncol == 3 .AND. nz == 5) THEN
    WRITE(*,'(A)') "[PASS] CSC complex get_dimensions"
  ELSE
    WRITE(*,'(A)') "[FAIL] CSC complex get_dimensions"
    test_passed = .FALSE.
  END IF
  
  CALL sparse_deallocate(mat_csc_c)
  
  ! Test 3: COO real matrix
  WRITE(*,'(A)') "Test 3: COO Real Matrix"
  mat_coo_r%nrow = 5
  mat_coo_r%ncol = 5
  mat_coo_r%nz = 8
  ALLOCATE(mat_coo_r%irow(mat_coo_r%nz))
  ALLOCATE(mat_coo_r%icol(mat_coo_r%nz))
  ALLOCATE(mat_coo_r%val(mat_coo_r%nz))
  
  CALL sparse_get_dimensions(mat_coo_r, nrow, ncol, nz)
  IF (nrow == 5 .AND. ncol == 5 .AND. nz == 8) THEN
    WRITE(*,'(A)') "[PASS] COO real get_dimensions"
  ELSE
    WRITE(*,'(A)') "[FAIL] COO real get_dimensions"
    test_passed = .FALSE.
  END IF
  
  CALL sparse_deallocate(mat_coo_r)
  
  ! Test 4: CSR real matrix
  WRITE(*,'(A)') "Test 4: CSR Real Matrix"
  mat_csr_r%nrow = 3
  mat_csr_r%ncol = 4
  mat_csr_r%nz = 6
  ALLOCATE(mat_csr_r%prow(mat_csr_r%nrow+1))
  ALLOCATE(mat_csr_r%icol(mat_csr_r%nz))
  ALLOCATE(mat_csr_r%val(mat_csr_r%nz))
  
  CALL sparse_get_dimensions(mat_csr_r, nrow, ncol, nz)
  IF (nrow == 3 .AND. ncol == 4 .AND. nz == 6) THEN
    WRITE(*,'(A)') "[PASS] CSR real get_dimensions"
  ELSE
    WRITE(*,'(A)') "[FAIL] CSR real get_dimensions"
    test_passed = .FALSE.
  END IF
  
  CALL sparse_deallocate(mat_csr_r)
  
  ! Test 5: Kind parameters
  WRITE(*,'(A)') "Test 5: Kind Parameters"
  IF (dp == KIND(1.0d0)) THEN
    WRITE(*,'(A)') "[PASS] dp parameter"
  ELSE
    WRITE(*,'(A)') "[FAIL] dp parameter"
    test_passed = .FALSE.
  END IF
  
  IF (long == 8) THEN
    WRITE(*,'(A)') "[PASS] long parameter"
  ELSE
    WRITE(*,'(A)') "[FAIL] long parameter"
    test_passed = .FALSE.
  END IF
  
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