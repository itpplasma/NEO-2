PROGRAM test_sparse_types
  ! Test for sparse_types_mod module
  
  USE sparse_types_mod
  IMPLICIT NONE
  
  ! Test variables
  LOGICAL :: test_passed
  REAL(kind=dp) :: test_real
  COMPLEX(kind=dp) :: test_complex
  INTEGER(kind=long) :: test_long
  
  test_passed = .TRUE.
  
  WRITE(*,'(A)') "================================="
  WRITE(*,'(A)') "Sparse Types Module Test"
  WRITE(*,'(A)') "================================="
  WRITE(*,*)
  
  ! Test 1: dp parameter
  WRITE(*,'(A)') "Test 1: dp parameter"
  test_real = 1.0_dp
  IF (dp == KIND(1.0d0) .AND. KIND(test_real) == dp) THEN
    WRITE(*,'(A)') "[PASS] dp parameter"
  ELSE
    WRITE(*,'(A)') "[FAIL] dp parameter"
    test_passed = .FALSE.
  END IF
  
  ! Test 2: long parameter  
  WRITE(*,'(A)') "Test 2: long parameter"
  test_long = 1_long
  IF (long == 8 .AND. KIND(test_long) == long) THEN
    WRITE(*,'(A)') "[PASS] long parameter"
  ELSE
    WRITE(*,'(A)') "[FAIL] long parameter"
    test_passed = .FALSE.
  END IF
  
  ! Test 3: Complex with dp
  WRITE(*,'(A)') "Test 3: Complex with dp"
  test_complex = (1.0_dp, 2.0_dp)
  IF (KIND(REAL(test_complex)) == dp .AND. KIND(AIMAG(test_complex)) == dp) THEN
    WRITE(*,'(A)') "[PASS] Complex dp parameter"
  ELSE
    WRITE(*,'(A)') "[FAIL] Complex dp parameter"
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