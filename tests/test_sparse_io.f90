PROGRAM test_sparse_io
  ! Test for sparse_io_mod module
  
  USE sparse_types_mod, ONLY: dp
  USE sparse_io_mod
  IMPLICIT NONE
  
  ! Test variables
  INTEGER :: nrow, ncol, nz
  INTEGER, ALLOCATABLE :: irow(:), pcol(:)
  REAL(kind=dp), ALLOCATABLE :: val(:), A_full(:,:)
  COMPLEX(kind=dp), ALLOCATABLE :: z_val(:)
  CHARACTER(len=256) :: test_file
  INTEGER :: unit, i
  LOGICAL :: test_passed, file_exists
  REAL(kind=dp), PARAMETER :: tol = 1.0e-14_dp
  
  test_passed = .TRUE.
  
  WRITE(*,'(A)') "================================="
  WRITE(*,'(A)') "Sparse I/O Module Test"
  WRITE(*,'(A)') "================================="
  WRITE(*,*)
  
  ! Test 1: Load mini example
  WRITE(*,'(A)') "Test 1: Load mini example"
  CALL load_mini_example(A_full)
  
  IF (SIZE(A_full,1) == 5 .AND. SIZE(A_full,2) == 5) THEN
    ! Check some specific values
    IF (ABS(A_full(1,1) - 1.0_dp) < tol .AND. &
        ABS(A_full(2,4)) < tol .AND. &  ! Should be zero
        ABS(A_full(3,3)) < tol) THEN    ! Should be zero
      WRITE(*,'(A)') "[PASS] Load mini example"
    ELSE
      WRITE(*,'(A)') "[FAIL] Load mini example - incorrect values"
      test_passed = .FALSE.
    END IF
  ELSE
    WRITE(*,'(A)') "[FAIL] Load mini example - incorrect size"
    test_passed = .FALSE.
  END IF
  
  DEALLOCATE(A_full)
  
  ! Test 2: Create and load compressed example
  WRITE(*,'(A)') "Test 2: Load compressed example"
  
  ! Create a test file
  test_file = "test_compressed_matrix.dat"
  unit = 20
  OPEN(unit=unit, file=test_file, status='replace', action='write')
  WRITE(unit,*) 3, 3, 5  ! nrow, ncol, nz
  WRITE(unit,*) 1, 2, 2, 1, 3  ! irow
  WRITE(unit,*) 1, 3, 4, 6     ! pcol
  WRITE(unit,*) 1.0, 2.0, 3.0, 4.0, 5.0  ! val
  CLOSE(unit)
  
  ! Load the file
  CALL load_compressed_example(test_file, nrow, ncol, nz, irow, pcol, val)
  
  IF (nrow == 3 .AND. ncol == 3 .AND. nz == 5) THEN
    IF (ALL(irow == (/1, 2, 2, 1, 3/)) .AND. &
        ALL(pcol == (/1, 3, 4, 6/)) .AND. &
        MAXVAL(ABS(val - (/1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp/))) < tol) THEN
      WRITE(*,'(A)') "[PASS] Load compressed example"
    ELSE
      WRITE(*,'(A)') "[FAIL] Load compressed example - incorrect data"
      test_passed = .FALSE.
    END IF
  ELSE
    WRITE(*,'(A)') "[FAIL] Load compressed example - incorrect dimensions"
    test_passed = .FALSE.
  END IF
  
  ! Clean up test file
  OPEN(unit=unit, file=test_file, status='old')
  CLOSE(unit, status='delete')
  
  IF (ALLOCATED(irow)) DEALLOCATE(irow)
  IF (ALLOCATED(pcol)) DEALLOCATE(pcol)
  IF (ALLOCATED(val)) DEALLOCATE(val)
  
  ! Test 3: Create and load Octave format matrix
  WRITE(*,'(A)') "Test 3: Load Octave matrix"
  
  ! Create a test file in Octave format
  test_file = "test_octave_matrix.dat"
  unit = 20
  OPEN(unit=unit, file=test_file, status='replace', action='write')
  WRITE(unit,*) 3, 3, 4  ! nrow, ncol, nz
  WRITE(unit,*) 1, 1, 1.0  ! row 1, col 1, value
  WRITE(unit,*) 2, 1, 2.0  ! row 2, col 1, value
  WRITE(unit,*) 2, 2, 3.0  ! row 2, col 2, value
  WRITE(unit,*) 3, 3, 4.0  ! row 3, col 3, value
  CLOSE(unit)
  
  ! Load the file
  CALL load_octave_matrices(test_file, nrow, ncol, nz, irow, pcol, val)
  
  IF (nrow == 3 .AND. ncol == 3 .AND. nz == 4) THEN
    ! Octave format is converted to CSC, check column pointers
    IF (pcol(1) == 1 .AND. pcol(2) == 3 .AND. &
        pcol(3) == 4 .AND. pcol(4) == 5) THEN
      WRITE(*,'(A)') "[PASS] Load Octave matrix"
    ELSE
      WRITE(*,'(A)') "[FAIL] Load Octave matrix - incorrect structure"
      test_passed = .FALSE.
    END IF
  ELSE
    WRITE(*,'(A)') "[FAIL] Load Octave matrix - incorrect dimensions"
    test_passed = .FALSE.
  END IF
  
  ! Clean up test file
  OPEN(unit=unit, file=test_file, status='old')
  CLOSE(unit, status='delete')
  
  IF (ALLOCATED(irow)) DEALLOCATE(irow)
  IF (ALLOCATED(pcol)) DEALLOCATE(pcol)
  IF (ALLOCATED(val)) DEALLOCATE(val)
  
  ! Test 4: Create and load complex Octave matrix
  WRITE(*,'(A)') "Test 4: Load complex Octave matrix"
  
  ! Create a test file
  test_file = "test_octave_complex.dat"
  unit = 20
  OPEN(unit=unit, file=test_file, status='replace', action='write')
  WRITE(unit,*) 2, 2, 3  ! nrow, ncol, nz
  WRITE(unit,*) 1, 1, "(1.0,0.0)"
  WRITE(unit,*) 2, 1, "(0.0,1.0)"
  WRITE(unit,*) 1, 2, "(2.0,-1.0)"
  CLOSE(unit)
  
  ! Load the file
  CALL load_octave_matrices(test_file, nrow, ncol, nz, irow, pcol, z_val)
  
  IF (nrow == 2 .AND. ncol == 2 .AND. nz == 3) THEN
    WRITE(*,'(A)') "[PASS] Load complex Octave matrix"
  ELSE
    WRITE(*,'(A)') "[FAIL] Load complex Octave matrix"
    test_passed = .FALSE.
  END IF
  
  ! Clean up test file
  OPEN(unit=unit, file=test_file, status='old')
  CLOSE(unit, status='delete')
  
  IF (ALLOCATED(irow)) DEALLOCATE(irow)
  IF (ALLOCATED(pcol)) DEALLOCATE(pcol)
  IF (ALLOCATED(z_val)) DEALLOCATE(z_val)
  
  ! Test 5: Test find_unit helper
  WRITE(*,'(A)') "Test 5: Find free unit"
  unit = 10
  CALL find_unit(unit)
  IF (unit >= 10) THEN
    WRITE(*,'(A)') "[PASS] Find free unit"
  ELSE
    WRITE(*,'(A)') "[FAIL] Find free unit"
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
  
END PROGRAM test_sparse_io