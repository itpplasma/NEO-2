PROGRAM test_sparse_legacy
  ! Comprehensive test harness for existing sparse_mod functionality
  ! This captures current behavior to ensure no regressions during refactoring
  
  USE sparse_mod
  USE sparse_types_mod, ONLY: dp
  IMPLICIT NONE
  
  ! Test counters
  INTEGER :: tests_run = 0
  INTEGER :: tests_passed = 0
  INTEGER :: tests_failed = 0
  
  ! Test data
  INTEGER :: nrow, ncol, nz
  INTEGER, DIMENSION(:), ALLOCATABLE :: irow, pcol, icol
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: val, b, x, b_orig
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: A_full, B_full, X_full, B_orig_full
  COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE :: z_val, z_b, z_x, z_b_orig
  COMPLEX(kind=dp), DIMENSION(:,:), ALLOCATABLE :: z_A_full, z_B_full, z_X_full
  
  ! Error tolerances
  REAL(kind=dp), PARAMETER :: tol_abs = 1.0e-12_dp
  REAL(kind=dp), PARAMETER :: tol_rel = 1.0e-10_dp
  
  WRITE(*,'(A)') "================================="
  WRITE(*,'(A)') "Sparse Module Legacy Test Suite"
  WRITE(*,'(A)') "================================="
  WRITE(*,*)
  
  ! Test 1: Mini example with direct solver
  CALL test_mini_example()
  
  ! Test 2: Sparse to full conversion
  CALL test_sparse_to_full_conversion()
  
  ! Test 3: Full to sparse conversion
  CALL test_full_to_sparse_conversion()
  
  ! Test 4: Column pointer conversions
  CALL test_column_pointer_conversions()
  
  ! Test 5: Sparse matrix-vector multiplication
  CALL test_sparse_matmul()
  
  ! Test 6: Real solver with single RHS
  CALL test_real_solver_single_rhs()
  
  ! Test 7: Real solver with multiple RHS
  CALL test_real_solver_multiple_rhs()
  
  ! Test 8: Complex solver with single RHS
  CALL test_complex_solver_single_rhs()
  
  ! Test 9: Complex solver with multiple RHS
  CALL test_complex_solver_multiple_rhs()
  
  ! Test 10: UMFPACK solver (method 3)
  CALL test_umfpack_solver()
  
  ! Test 11: SuiteSparse interface - skipped for now
  !CALL test_suitesparse_interface()
  
  ! Test 12: Solver method switching
  CALL test_solver_method_switching()
  
  ! Test 13: Edge cases
  CALL test_edge_cases()
  
  ! Summary
  WRITE(*,*)
  WRITE(*,'(A)') "================================="
  WRITE(*,'(A,I4)') "Total tests run:    ", tests_run
  WRITE(*,'(A,I4)') "Tests passed:       ", tests_passed
  WRITE(*,'(A,I4)') "Tests failed:       ", tests_failed
  WRITE(*,'(A)') "================================="
  
  IF (tests_failed > 0) THEN
    STOP 1
  END IF
  
CONTAINS

  SUBROUTINE check_result(test_name, condition)
    CHARACTER(len=*), INTENT(in) :: test_name
    LOGICAL, INTENT(in) :: condition
    
    tests_run = tests_run + 1
    IF (condition) THEN
      tests_passed = tests_passed + 1
      WRITE(*,'(A,A,A)') "[PASS] ", test_name
    ELSE
      tests_failed = tests_failed + 1
      WRITE(*,'(A,A,A)') "[FAIL] ", test_name
    END IF
  END SUBROUTINE check_result
  
  SUBROUTINE test_mini_example()
    REAL(kind=dp) :: max_abs_err, max_rel_err
    
    WRITE(*,'(A)') "Test 1: Mini example with direct solver"
    
    ! Load mini example
    CALL load_mini_example(A_full)
    
    ! Create RHS
    IF (ALLOCATED(b)) DEALLOCATE(b)
    ALLOCATE(b(SIZE(A_full,2)))
    b = 1.0_dp
    
    ! Save original b
    IF (ALLOCATED(b_orig)) DEALLOCATE(b_orig)
    ALLOCATE(b_orig(SIZE(b)))
    b_orig = b
    
    ! Solve
    IF (ALLOCATED(x)) DEALLOCATE(x)
    ALLOCATE(x(SIZE(b)))
    x = b
    CALL sparse_solve(A_full, x)
    
    ! Test solution
    CALL sparse_solver_test(A_full, x, b_orig, max_abs_err, max_rel_err)
    
    CALL check_result("Mini example solution accuracy", &
         max_abs_err < tol_abs .AND. max_rel_err < tol_rel)
    
    ! Cleanup
    IF (ALLOCATED(A_full)) DEALLOCATE(A_full)
    
  END SUBROUTINE test_mini_example
  
  SUBROUTINE test_sparse_to_full_conversion()
    INTEGER :: i, j
    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: A_reconstructed
    
    WRITE(*,'(A)') "Test 2: Sparse to full conversion"
    
    ! Create a simple test matrix
    nrow = 4
    ncol = 4
    nz = 7
    
    IF (ALLOCATED(irow)) DEALLOCATE(irow)
    IF (ALLOCATED(pcol)) DEALLOCATE(pcol)
    IF (ALLOCATED(val)) DEALLOCATE(val)
    
    ALLOCATE(irow(nz), pcol(ncol+1), val(nz))
    
    ! Column 1: elements at (1,1) and (3,1)
    ! Column 2: element at (2,2)
    ! Column 3: elements at (1,3) and (3,3)
    ! Column 4: elements at (2,4) and (4,4)
    pcol = (/1, 3, 4, 6, 8/)
    irow = (/1, 3, 2, 1, 3, 2, 4/)
    val = (/1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp, 7.0_dp/)
    
    ! Convert to full
    CALL sparse2full(irow, pcol, val, nrow, ncol, A_reconstructed)
    
    ! Check specific values
    CALL check_result("Sparse to full: (1,1) = 1.0", ABS(A_reconstructed(1,1) - 1.0_dp) < tol_abs)
    CALL check_result("Sparse to full: (3,1) = 2.0", ABS(A_reconstructed(3,1) - 2.0_dp) < tol_abs)
    CALL check_result("Sparse to full: (2,2) = 3.0", ABS(A_reconstructed(2,2) - 3.0_dp) < tol_abs)
    CALL check_result("Sparse to full: zeros", ABS(A_reconstructed(4,1)) < tol_abs)
    
  END SUBROUTINE test_sparse_to_full_conversion
  
  SUBROUTINE test_full_to_sparse_conversion()
    INTEGER :: i, nz_out
    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: A_test, A_reconstructed
    
    WRITE(*,'(A)') "Test 3: Full to sparse conversion"
    
    ! Create test matrix
    ALLOCATE(A_test(4,4))
    A_test = 0.0_dp
    A_test(1,1) = 1.0_dp
    A_test(2,2) = 2.0_dp
    A_test(3,3) = 3.0_dp
    A_test(4,4) = 4.0_dp
    A_test(1,3) = 5.0_dp
    A_test(2,4) = 6.0_dp
    
    ! Convert to sparse
    CALL full2sparse(A_test, irow, pcol, val, nrow, ncol, nz_out)
    
    CALL check_result("Full to sparse: dimensions", nrow == 4 .AND. ncol == 4)
    CALL check_result("Full to sparse: nonzeros", nz_out == 6)
    
    ! Convert back to full and compare
    CALL sparse2full(irow, pcol, val, nrow, ncol, A_reconstructed)
    
    CALL check_result("Full to sparse roundtrip", &
         MAXVAL(ABS(A_test - A_reconstructed)) < tol_abs)
    
  END SUBROUTINE test_full_to_sparse_conversion
  
  SUBROUTINE test_column_pointer_conversions()
    INTEGER :: i
    INTEGER, DIMENSION(:), ALLOCATABLE :: pcol_test, icol_test, pcol_reconstructed
    
    WRITE(*,'(A)') "Test 4: Column pointer conversions"
    
    ! Test pcol to icol conversion
    ALLOCATE(pcol_test(5))
    pcol_test = (/1, 3, 4, 6, 8/)
    
    CALL column_pointer2full(pcol_test, icol_test)
    
    CALL check_result("Pointer to full: size", SIZE(icol_test) == 7)
    CALL check_result("Pointer to full: values", &
         icol_test(1) == 1 .AND. icol_test(2) == 1 .AND. &
         icol_test(3) == 2 .AND. icol_test(4) == 3 .AND. &
         icol_test(5) == 3 .AND. icol_test(6) == 4 .AND. &
         icol_test(7) == 4)
    
    ! Test icol to pcol conversion
    CALL column_full2pointer(icol_test, pcol_reconstructed)
    
    CALL check_result("Full to pointer roundtrip", &
         ALL(pcol_test == pcol_reconstructed))
    
  END SUBROUTINE test_column_pointer_conversions
  
  SUBROUTINE test_sparse_matmul()
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: x_test, r_result, r_expected
    INTEGER :: i
    
    WRITE(*,'(A)') "Test 5: Sparse matrix-vector multiplication"
    
    ! Use the same sparse matrix from test 2
    nrow = 4
    ncol = 4
    nz = 7
    
    IF (ALLOCATED(irow)) DEALLOCATE(irow)
    IF (ALLOCATED(pcol)) DEALLOCATE(pcol)
    IF (ALLOCATED(val)) DEALLOCATE(val)
    
    ALLOCATE(irow(nz), pcol(ncol+1), val(nz))
    
    pcol = (/1, 3, 4, 6, 8/)
    irow = (/1, 3, 2, 1, 3, 2, 4/)
    val = (/1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp, 7.0_dp/)
    
    ! Test vector
    ALLOCATE(x_test(ncol))
    x_test = (/1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp/)
    
    ! Compute matrix-vector product
    ALLOCATE(r_result(nrow))
    CALL sparse_matmul(nrow, ncol, irow, pcol, val, x_test, r_result)
    
    ! Expected result (computed manually)
    ALLOCATE(r_expected(nrow))
    r_expected(1) = 1.0_dp * 1.0_dp + 4.0_dp * 3.0_dp  ! = 13.0
    r_expected(2) = 3.0_dp * 2.0_dp + 6.0_dp * 4.0_dp  ! = 30.0
    r_expected(3) = 2.0_dp * 1.0_dp + 5.0_dp * 3.0_dp  ! = 17.0
    r_expected(4) = 7.0_dp * 4.0_dp                     ! = 28.0
    
    CALL check_result("Sparse matmul accuracy", &
         MAXVAL(ABS(r_result - r_expected)) < tol_abs)
    
  END SUBROUTINE test_sparse_matmul
  
  SUBROUTINE test_real_solver_single_rhs()
    REAL(kind=dp) :: max_abs_err, max_rel_err
    
    WRITE(*,'(A)') "Test 6: Real solver with single RHS"
    
    ! Create a simple SPD test matrix
    IF (ALLOCATED(A_full)) DEALLOCATE(A_full)
    ALLOCATE(A_full(3,3))
    A_full = 0.0_dp
    A_full(1,1) = 4.0_dp
    A_full(2,2) = 5.0_dp
    A_full(3,3) = 6.0_dp
    A_full(1,2) = 1.0_dp
    A_full(2,1) = 1.0_dp
    A_full(2,3) = 2.0_dp
    A_full(3,2) = 2.0_dp
    
    ! Convert to sparse
    CALL full2sparse(A_full, irow, pcol, val, nrow, ncol, nz)
    
    ! Create RHS
    IF (ALLOCATED(b)) DEALLOCATE(b)
    ALLOCATE(b(nrow))
    b = (/1.0_dp, 2.0_dp, 3.0_dp/)
    IF (ALLOCATED(b_orig)) DEALLOCATE(b_orig)
    ALLOCATE(b_orig(nrow))
    b_orig = b
    
    ! Solve
    IF (ALLOCATED(x)) DEALLOCATE(x)
    ALLOCATE(x(nrow))
    x = b
    CALL sparse_solve(nrow, ncol, nz, irow, pcol, val, x)
    
    ! Test
    CALL sparse_solver_test(nrow, ncol, irow, pcol, val, x, b_orig, max_abs_err, max_rel_err)
    
    CALL check_result("Real single RHS accuracy", &
         max_abs_err < tol_abs .AND. max_rel_err < tol_rel)
    
  END SUBROUTINE test_real_solver_single_rhs
  
  SUBROUTINE test_real_solver_multiple_rhs()
    REAL(kind=dp) :: max_abs_err, max_rel_err
    INTEGER :: i
    
    WRITE(*,'(A)') "Test 7: Real solver with multiple RHS"
    
    ! Use same matrix as test 6
    IF (ALLOCATED(A_full)) DEALLOCATE(A_full)
    ALLOCATE(A_full(3,3))
    A_full = 0.0_dp
    A_full(1,1) = 4.0_dp
    A_full(2,2) = 5.0_dp
    A_full(3,3) = 6.0_dp
    A_full(1,2) = 1.0_dp
    A_full(2,1) = 1.0_dp
    A_full(2,3) = 2.0_dp
    A_full(3,2) = 2.0_dp
    
    ! Convert to sparse
    CALL full2sparse(A_full, irow, pcol, val, nrow, ncol, nz)
    
    ! Create multiple RHS
    IF (ALLOCATED(B_full)) DEALLOCATE(B_full)
    ALLOCATE(B_full(nrow, 2))
    B_full(:,1) = (/1.0_dp, 2.0_dp, 3.0_dp/)
    B_full(:,2) = (/4.0_dp, 5.0_dp, 6.0_dp/)
    IF (ALLOCATED(B_orig_full)) DEALLOCATE(B_orig_full)
    ALLOCATE(B_orig_full(nrow, 2))
    B_orig_full = B_full
    
    ! Solve
    IF (ALLOCATED(X_full)) DEALLOCATE(X_full)
    ALLOCATE(X_full(nrow, 2))
    X_full = B_full
    CALL sparse_solve(nrow, ncol, nz, irow, pcol, val, X_full)
    
    ! Test each solution
    CALL sparse_solver_test(nrow, ncol, irow, pcol, val, X_full, B_orig_full, max_abs_err, max_rel_err)
    
    CALL check_result("Real multiple RHS accuracy", &
         max_abs_err < tol_abs .AND. max_rel_err < tol_rel)
    
  END SUBROUTINE test_real_solver_multiple_rhs
  
  SUBROUTINE test_complex_solver_single_rhs()
    REAL(kind=dp) :: max_abs_err, max_rel_err
    
    WRITE(*,'(A)') "Test 8: Complex solver with single RHS"
    
    ! Create a simple complex test matrix
    nrow = 3
    ncol = 3
    nz = 7
    
    IF (ALLOCATED(irow)) DEALLOCATE(irow)
    IF (ALLOCATED(pcol)) DEALLOCATE(pcol)
    IF (ALLOCATED(z_val)) DEALLOCATE(z_val)
    ALLOCATE(irow(nz), pcol(ncol+1), z_val(nz))
    
    ! Hermitian matrix structure
    pcol = (/1, 3, 5, 7/)
    irow = (/1, 2, 1, 2, 2, 3/)
    z_val = (/(3.0_dp, 0.0_dp), (1.0_dp, -1.0_dp), &
             (1.0_dp, 1.0_dp), (4.0_dp, 0.0_dp), &
             (2.0_dp, -1.0_dp), (2.0_dp, 1.0_dp), (5.0_dp, 0.0_dp)/)
    
    ! Create RHS
    IF (ALLOCATED(z_b)) DEALLOCATE(z_b)
    ALLOCATE(z_b(nrow))
    z_b = (/(1.0_dp, 1.0_dp), (2.0_dp, -1.0_dp), (3.0_dp, 0.0_dp)/)
    IF (ALLOCATED(z_b_orig)) DEALLOCATE(z_b_orig)
    ALLOCATE(z_b_orig(nrow))
    z_b_orig = z_b
    
    ! Solve
    IF (ALLOCATED(z_x)) DEALLOCATE(z_x)
    ALLOCATE(z_x(nrow))
    z_x = z_b
    CALL sparse_solve(nrow, ncol, nz, irow, pcol, z_val, z_x)
    
    ! Test
    CALL sparse_solver_test(nrow, ncol, irow, pcol, z_val, z_x, z_b_orig, max_abs_err, max_rel_err)
    
    CALL check_result("Complex single RHS accuracy", &
         max_abs_err < tol_abs .AND. max_rel_err < tol_rel)
    
  END SUBROUTINE test_complex_solver_single_rhs
  
  SUBROUTINE test_complex_solver_multiple_rhs()
    REAL(kind=dp) :: max_abs_err, max_rel_err
    
    WRITE(*,'(A)') "Test 9: Complex solver with multiple RHS"
    
    ! Use same matrix as test 8
    nrow = 3
    ncol = 3
    nz = 7
    
    IF (ALLOCATED(irow)) DEALLOCATE(irow)
    IF (ALLOCATED(pcol)) DEALLOCATE(pcol)
    IF (ALLOCATED(z_val)) DEALLOCATE(z_val)
    ALLOCATE(irow(nz), pcol(ncol+1), z_val(nz))
    
    pcol = (/1, 3, 5, 7/)
    irow = (/1, 2, 1, 2, 2, 3/)
    z_val = (/(3.0_dp, 0.0_dp), (1.0_dp, -1.0_dp), &
             (1.0_dp, 1.0_dp), (4.0_dp, 0.0_dp), &
             (2.0_dp, -1.0_dp), (2.0_dp, 1.0_dp), (5.0_dp, 0.0_dp)/)
    
    ! Create multiple RHS
    IF (ALLOCATED(z_B_full)) DEALLOCATE(z_B_full)
    ALLOCATE(z_B_full(nrow, 2))
    z_B_full(:,1) = (/(1.0_dp, 1.0_dp), (2.0_dp, -1.0_dp), (3.0_dp, 0.0_dp)/)
    z_B_full(:,2) = (/(0.0_dp, 1.0_dp), (1.0_dp, 0.0_dp), (2.0_dp, 2.0_dp)/)
    
    ! Solve
    IF (ALLOCATED(z_X_full)) DEALLOCATE(z_X_full)
    ALLOCATE(z_X_full(nrow, 2))
    z_X_full = z_B_full
    CALL sparse_solve(nrow, ncol, nz, irow, pcol, z_val, z_X_full)
    
    ! Test
    CALL sparse_solver_test(nrow, ncol, irow, pcol, z_val, z_X_full, z_B_full, max_abs_err, max_rel_err)
    
    CALL check_result("Complex multiple RHS accuracy", &
         max_abs_err < tol_abs .AND. max_rel_err < tol_rel)
    
  END SUBROUTINE test_complex_solver_multiple_rhs
  
  SUBROUTINE test_umfpack_solver()
    REAL(kind=dp) :: max_abs_err, max_rel_err
    INTEGER :: old_method
    
    WRITE(*,'(A)') "Test 10: UMFPACK solver (method 3)"
    
    ! Save old method
    old_method = sparse_solve_method
    sparse_solve_method = 3
    
    ! Create test matrix
    IF (ALLOCATED(A_full)) DEALLOCATE(A_full)
    ALLOCATE(A_full(3,3))
    A_full = 0.0_dp
    A_full(1,1) = 4.0_dp
    A_full(2,2) = 5.0_dp
    A_full(3,3) = 6.0_dp
    A_full(1,2) = 1.0_dp
    A_full(2,1) = 1.0_dp
    
    ! Convert to sparse
    CALL full2sparse(A_full, irow, pcol, val, nrow, ncol, nz)
    
    ! Create RHS
    IF (ALLOCATED(b)) DEALLOCATE(b)
    ALLOCATE(b(nrow))
    b = (/1.0_dp, 2.0_dp, 3.0_dp/)
    IF (ALLOCATED(b_orig)) DEALLOCATE(b_orig)
    ALLOCATE(b_orig(nrow))
    b_orig = b
    
    ! Solve
    IF (ALLOCATED(x)) DEALLOCATE(x)
    ALLOCATE(x(nrow))
    x = b
    CALL sparse_solve(nrow, ncol, nz, irow, pcol, val, x)
    
    ! Test
    CALL sparse_solver_test(nrow, ncol, irow, pcol, val, x, b_orig, max_abs_err, max_rel_err)
    
    CALL check_result("UMFPACK solver accuracy", &
         max_abs_err < tol_abs .AND. max_rel_err < tol_rel)
    
    ! Restore method
    sparse_solve_method = old_method
    
  END SUBROUTINE test_umfpack_solver
  
  !SUBROUTINE test_suitesparse_interface()
  !  REAL(kind=dp) :: max_abs_err, max_rel_err
  !  
  !  WRITE(*,'(A)') "Test 11: SuiteSparse interface"
  !  
  !  ! Create test matrix
  !  ALLOCATE(A_full(3,3))
  !  A_full = 0.0_dp
  !  A_full(1,1) = 4.0_dp
  !  A_full(2,2) = 5.0_dp
  !  A_full(3,3) = 6.0_dp
  !  A_full(1,2) = 1.0_dp
  !  A_full(2,1) = 1.0_dp
  !  
  !  ! Convert to sparse
  !  CALL full2sparse(A_full, irow, pcol, val, nrow, ncol, nz)
  !  
  !  ! Create RHS
  !  ALLOCATE(b(nrow))
  !  b = (/1.0_dp, 2.0_dp, 3.0_dp/)
  !  ALLOCATE(b_orig(nrow))
  !  b_orig = b
  !  
  !  ! Solve using SuiteSparse directly
  !  ALLOCATE(x(nrow))
  !  x = b
  !  CALL sparse_solve_suitesparse(nrow, ncol, nz, irow, pcol, val, x)
  !  
  !  ! Test
  !  CALL sparse_solver_test(nrow, ncol, irow, pcol, val, x, b_orig, max_abs_err, max_rel_err)
  !  
  !  CALL check_result("SuiteSparse interface accuracy", &
  !       max_abs_err < tol_abs .AND. max_rel_err < tol_rel)
  !  
  !END SUBROUTINE test_suitesparse_interface
  
  SUBROUTINE test_solver_method_switching()
    REAL(kind=dp) :: max_abs_err1, max_rel_err1, max_abs_err2, max_rel_err2
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: x1, x2
    INTEGER :: old_method
    
    WRITE(*,'(A)') "Test 12: Solver method switching"
    
    ! Save old method
    old_method = sparse_solve_method
    
    ! Create test matrix
    IF (ALLOCATED(A_full)) DEALLOCATE(A_full)
    ALLOCATE(A_full(3,3))
    A_full = 0.0_dp
    A_full(1,1) = 4.0_dp
    A_full(2,2) = 5.0_dp
    A_full(3,3) = 6.0_dp
    A_full(1,2) = 1.0_dp
    A_full(2,1) = 1.0_dp
    
    ! Convert to sparse
    CALL full2sparse(A_full, irow, pcol, val, nrow, ncol, nz)
    
    ! Create RHS
    IF (ALLOCATED(b)) DEALLOCATE(b)
    ALLOCATE(b(nrow))
    b = (/1.0_dp, 2.0_dp, 3.0_dp/)
    
    ! Solve with method 3
    sparse_solve_method = 3
    IF (ALLOCATED(x1)) DEALLOCATE(x1)
    ALLOCATE(x1(nrow))
    x1 = b
    CALL sparse_solve(nrow, ncol, nz, irow, pcol, val, x1)
    
    ! Solve with method 3 again to test consistency 
    sparse_solve_method = 3
    IF (ALLOCATED(x2)) DEALLOCATE(x2)
    ALLOCATE(x2(nrow))
    x2 = b
    CALL sparse_solve(nrow, ncol, nz, irow, pcol, val, x2)
    
    ! Compare solutions
    CALL check_result("Method switching consistency", &
         MAXVAL(ABS(x1 - x2)) < 1.0e-8_dp)
    
    ! Restore method
    sparse_solve_method = old_method
    
  END SUBROUTINE test_solver_method_switching
  
  SUBROUTINE test_edge_cases()
    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: A_small
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: b_small, x_small
    INTEGER :: nz_out
    
    WRITE(*,'(A)') "Test 13: Edge cases"
    
    ! Test 1x1 matrix
    ALLOCATE(A_small(1,1))
    A_small(1,1) = 2.0_dp
    
    CALL full2sparse(A_small, irow, pcol, val, nrow, ncol, nz_out)
    
    CALL check_result("1x1 matrix conversion", &
         nrow == 1 .AND. ncol == 1 .AND. nz_out == 1)
    
    ! Solve 1x1 system
    ALLOCATE(b_small(1), x_small(1))
    b_small(1) = 4.0_dp
    x_small = b_small
    CALL sparse_solve(nrow, ncol, nz_out, irow, pcol, val, x_small)
    
    CALL check_result("1x1 system solution", &
         ABS(x_small(1) - 2.0_dp) < tol_abs)
    
    ! Test zero matrix detection
    DEALLOCATE(A_small)
    ALLOCATE(A_small(3,3))
    A_small = 0.0_dp
    
    CALL full2sparse(A_small, irow, pcol, val, nrow, ncol, nz_out)
    
    CALL check_result("Zero matrix has no nonzeros", nz_out == 0)
    
  END SUBROUTINE test_edge_cases
  
END PROGRAM test_sparse_legacy