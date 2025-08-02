PROGRAM test_sparse_solvers
  ! Test for sparse_solvers_mod module
  
  USE sparse_types_mod, ONLY: dp, long
  USE sparse_solvers_mod
  USE sparse_conversion_mod, ONLY: full2sparse
  IMPLICIT NONE
  
  ! Test variables
  INTEGER :: nrow, ncol, nz
  INTEGER, ALLOCATABLE :: irow(:), pcol(:)
  REAL(kind=dp), ALLOCATABLE :: val(:), b(:,:), x(:)
  REAL(kind=dp), ALLOCATABLE :: A_full(:,:)
  COMPLEX(kind=dp), ALLOCATABLE :: z_val(:), z_b(:), z_x(:)
  COMPLEX(kind=dp), ALLOCATABLE :: z_A_full(:,:)
  INTEGER :: iopt, info
  REAL(kind=dp) :: max_abs_err, max_rel_err
  INTEGER :: i, j
  LOGICAL :: test_passed
  REAL(kind=dp), PARAMETER :: tol = 1.0e-12_dp
  
  test_passed = .TRUE.
  
  WRITE(*,'(A)') "================================="
  WRITE(*,'(A)') "Sparse Solvers Module Test"
  WRITE(*,'(A)') "================================="
  WRITE(*,*)
  
  ! Test 1: Direct solver for real system (method 3 - UMFPACK)
  WRITE(*,'(A)') "Test 1: Direct solver for real system (UMFPACK)"
  
  ! Create a simple 3x3 test system
  ALLOCATE(A_full(3,3))
  A_full = 0.0_dp
  A_full(1,1) = 4.0_dp
  A_full(1,2) = 1.0_dp
  A_full(2,1) = 1.0_dp
  A_full(2,2) = 3.0_dp
  A_full(3,3) = 2.0_dp
  
  ! Convert to sparse format
  CALL full2sparse(A_full, irow, pcol, val, nrow, ncol, nz)
  
  ! Create RHS
  ALLOCATE(b(3,1))
  b(:,1) = (/5.0_dp, 4.0_dp, 2.0_dp/)  ! Solution should be x = [1, 1, 1]
  
  ! Set solver method
  sparse_solve_method = 3  ! UMFPACK
  iopt = 0  ! Full solve (factorize + solve + cleanup)
  
  ! Solve the system
  CALL sparse_solve(nrow, ncol, nz, irow, pcol, val, b, iopt)
  
  ! Check solution
  IF (ABS(b(1,1) - 1.0_dp) < tol .AND. &
      ABS(b(2,1) - 1.0_dp) < tol .AND. &
      ABS(b(3,1) - 1.0_dp) < tol) THEN
    WRITE(*,'(A)') "[PASS] Direct solver (UMFPACK)"
  ELSE
    WRITE(*,'(A)') "[FAIL] Direct solver (UMFPACK)"
    WRITE(*,'(A,3F10.6)') "  Solution: ", b(:,1)
    test_passed = .FALSE.
  END IF
  
  DEALLOCATE(A_full, b, irow, pcol, val)
  
  ! Test 2: Direct solver with multiple RHS
  WRITE(*,'(A)') "Test 2: Direct solver with multiple RHS"
  
  ! Create same system
  ALLOCATE(A_full(3,3))
  A_full = 0.0_dp
  A_full(1,1) = 4.0_dp
  A_full(1,2) = 1.0_dp
  A_full(2,1) = 1.0_dp
  A_full(2,2) = 3.0_dp
  A_full(3,3) = 2.0_dp
  
  CALL full2sparse(A_full, irow, pcol, val, nrow, ncol, nz)
  
  ! Create multiple RHS as 2D array
  ALLOCATE(b(3,2))  ! 3 rows, 2 RHS vectors
  b(:,1) = (/5.0_dp, 4.0_dp, 2.0_dp/)    ! First RHS
  b(:,2) = (/8.0_dp, 7.0_dp, 4.0_dp/)    ! Second RHS
  
  iopt = 0
  CALL sparse_solve(nrow, ncol, nz, irow, pcol, val, b, iopt)
  
  IF (ABS(b(1,1) - 1.0_dp) < tol .AND. &
      ABS(b(2,1) - 1.0_dp) < tol .AND. &
      ABS(b(3,1) - 1.0_dp) < tol .AND. &
      ABS(b(1,2) - (17.0_dp/11.0_dp)) < tol .AND. &
      ABS(b(2,2) - (20.0_dp/11.0_dp)) < tol .AND. &
      ABS(b(3,2) - 2.0_dp) < tol) THEN
    WRITE(*,'(A)') "[PASS] Multiple RHS"
  ELSE
    WRITE(*,'(A)') "[FAIL] Multiple RHS"
    WRITE(*,'(A,3F12.8)') "  Solution 1: ", b(:,1)
    WRITE(*,'(A,3F12.8)') "  Solution 2: ", b(:,2)
    WRITE(*,'(A,3F12.8)') "  Expected 1: ", (/1.0_dp, 1.0_dp, 1.0_dp/)
    WRITE(*,'(A,3F12.8)') "  Expected 2: ", (/17.0_dp/11.0_dp, 20.0_dp/11.0_dp, 2.0_dp/)
    test_passed = .FALSE.
  END IF
  
  DEALLOCATE(A_full, b, irow, pcol, val)
  
  ! Test 3: Complex system solver
  WRITE(*,'(A)') "Test 3: Complex system solver"
  
  ALLOCATE(z_A_full(2,2))
  z_A_full(1,1) = (2.0_dp, 0.0_dp)
  z_A_full(1,2) = (0.0_dp, 1.0_dp)
  z_A_full(2,1) = (0.0_dp, -1.0_dp)
  z_A_full(2,2) = (2.0_dp, 0.0_dp)
  
  CALL full2sparse(z_A_full, irow, pcol, z_val, nrow, ncol, nz)
  
  ALLOCATE(z_b(2))
  z_b = (/(2.0_dp, 1.0_dp), (2.0_dp, -1.0_dp)/)  ! Solution: [1, i]
  
  
  iopt = 0
  CALL sparse_solve(nrow, ncol, nz, irow, pcol, z_val, z_b, iopt)
  
  
  IF (ABS(REAL(z_b(1)) - 1.0_dp) < tol .AND. &
      ABS(AIMAG(z_b(1)) - 0.0_dp) < tol .AND. &
      ABS(REAL(z_b(2)) - 1.0_dp) < tol .AND. &
      ABS(AIMAG(z_b(2)) - 0.0_dp) < tol) THEN
    WRITE(*,'(A)') "[PASS] Complex solver"
  ELSE
    WRITE(*,'(A)') "[FAIL] Complex solver"
    test_passed = .FALSE.
  END IF
  
  DEALLOCATE(z_A_full, z_b, irow, pcol, z_val)
  
  ! Test 4: Full matrix interface
  WRITE(*,'(A)') "Test 4: Full matrix interface"
  
  ALLOCATE(A_full(3,3))
  A_full = 0.0_dp
  A_full(1,1) = 2.0_dp
  A_full(2,2) = 3.0_dp
  A_full(3,3) = 4.0_dp
  
  ALLOCATE(b(3,1))
  b(:,1) = (/2.0_dp, 6.0_dp, 12.0_dp/)  ! Solution: [1, 2, 3]
  
  iopt = 0
  CALL sparse_solve(A_full, b, iopt)
  
  IF (ABS(b(1,1) - 1.0_dp) < tol .AND. &
      ABS(b(2,1) - 2.0_dp) < tol .AND. &
      ABS(b(3,1) - 3.0_dp) < tol) THEN
    WRITE(*,'(A)') "[PASS] Full matrix interface"
  ELSE
    WRITE(*,'(A)') "[FAIL] Full matrix interface"
    test_passed = .FALSE.
  END IF
  
  DEALLOCATE(A_full, b)
  
  ! Test 5: Solver with factorization reuse
  WRITE(*,'(A)') "Test 5: Factorization reuse"
  
  ALLOCATE(A_full(3,3))
  A_full = 0.0_dp
  A_full(1,1) = 4.0_dp
  A_full(1,2) = 1.0_dp
  A_full(2,1) = 1.0_dp
  A_full(2,2) = 3.0_dp
  A_full(3,3) = 2.0_dp
  
  CALL full2sparse(A_full, irow, pcol, val, nrow, ncol, nz)
  
  ! First solve with factorization
  ALLOCATE(b(3,1))
  b(:,1) = (/5.0_dp, 4.0_dp, 2.0_dp/)
  iopt = 0  ! Perform factorization
  CALL sparse_solve(nrow, ncol, nz, irow, pcol, val, b, iopt)
  
  ! Second solve reusing factorization
  b(:,1) = (/8.0_dp, 7.0_dp, 4.0_dp/)
  iopt = 1  ! Reuse factorization
  CALL sparse_solve(nrow, ncol, nz, irow, pcol, val, b, iopt)
  
  
  IF (ABS(b(1,1) - (17.0_dp/11.0_dp)) < tol .AND. &
      ABS(b(2,1) - (20.0_dp/11.0_dp)) < tol .AND. &
      ABS(b(3,1) - 2.0_dp) < tol) THEN
    WRITE(*,'(A)') "[PASS] Factorization reuse"
  ELSE
    WRITE(*,'(A)') "[FAIL] Factorization reuse"
    test_passed = .FALSE.
  END IF
  
  DEALLOCATE(A_full, b, irow, pcol, val)
  
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
  
END PROGRAM test_sparse_solvers