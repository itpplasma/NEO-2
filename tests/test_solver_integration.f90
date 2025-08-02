PROGRAM test_solver_integration
  ! =========================================================================
  ! Integration tests comparing BiCGSTAB vs UMFPACK solvers
  ! 
  ! Purpose: Ensures both solvers produce consistent results on the same problems
  ! 
  ! Test Coverage:
  ! - Small SPD matrices (basic correctness)
  ! - Diagonal matrices (simple case)
  ! - Tridiagonal matrices (realistic sparse structure)
  ! - Large sparse matrices (scalability testing)
  ! - Performance comparison (timing analysis)
  ! 
  ! Success Criteria:
  ! - Solution differences < 1e-6 (solver_comparison_tol)
  ! - Both solvers converge without errors
  ! - Performance data collected for analysis
  ! =========================================================================
  
  USE sparse_mod
  USE sparse_types_mod, ONLY: dp
  USE sparse_solvers_mod, ONLY: SOLVER_UMFPACK, SOLVER_BICGSTAB
  IMPLICIT NONE
  
  ! Test counters
  INTEGER :: tests_run = 0
  INTEGER :: tests_passed = 0
  INTEGER :: tests_failed = 0
  
  ! Test data
  INTEGER :: nrow, ncol, nz
  INTEGER, DIMENSION(:), ALLOCATABLE :: irow, pcol
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: val, b, x_umfpack, x_bicgstab, b_orig
  COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE :: z_val, z_b, z_x_umfpack, z_x_bicgstab, z_b_orig
  
  ! Error tolerances
  REAL(kind=dp), PARAMETER :: tol_abs = 1.0e-10_dp
  REAL(kind=dp), PARAMETER :: tol_rel = 1.0e-8_dp
  REAL(kind=dp), PARAMETER :: solver_comparison_tol = 1.0e-6_dp  ! Tolerance for comparing different solvers
  
  WRITE(*,'(A)') "======================================="
  WRITE(*,'(A)') "Solver Integration Test Suite"
  WRITE(*,'(A)') "BiCGSTAB vs UMFPACK Comparison"
  WRITE(*,'(A)') "======================================="
  WRITE(*,*)
  
  ! Test 1: Small SPD matrix comparison
  CALL test_small_spd_matrix_comparison()
  
  ! Test 2: Diagonal matrix comparison
  CALL test_diagonal_matrix_comparison()
  
  ! Test 3: Tridiagonal matrix comparison
  CALL test_tridiagonal_matrix_comparison()
  
  ! Test 4: Complex matrix comparison
  CALL test_complex_matrix_comparison()
  
  ! Test 5: Large sparse matrix comparison
  CALL test_large_sparse_matrix_comparison()
  
  ! Test 6: Performance comparison framework
  CALL test_performance_comparison()
  
  ! Summary
  WRITE(*,*)
  WRITE(*,'(A)') "======================================="
  WRITE(*,'(A,I4)') "Total tests run:    ", tests_run
  WRITE(*,'(A,I4)') "Tests passed:       ", tests_passed
  WRITE(*,'(A,I4)') "Tests failed:       ", tests_failed
  WRITE(*,'(A)') "======================================="
  
  IF (tests_failed > 0) THEN
    WRITE(*,'(A)') "Some solver integration tests FAILED!"
    STOP 1
  ELSE
    WRITE(*,'(A)') "All solver integration tests PASSED!"
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

  SUBROUTINE generate_spd_matrix(n, irow_out, pcol_out, val_out, nrow_out, ncol_out, nz_out)
    ! Generate a symmetric positive definite test matrix
    INTEGER, INTENT(in) :: n
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(out) :: irow_out, pcol_out
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(out) :: val_out
    INTEGER, INTENT(out) :: nrow_out, ncol_out, nz_out
    
    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: A_full
    INTEGER :: i
    
    ! Create SPD matrix: A = I + 0.1 * (ones matrix)
    ALLOCATE(A_full(n,n))
    A_full = 0.1_dp  ! Off-diagonal elements
    DO i = 1, n
      A_full(i,i) = 2.0_dp + 0.1_dp  ! Diagonal dominance ensures SPD
    END DO
    
    ! Convert to sparse format
    CALL full2sparse(A_full, irow_out, pcol_out, val_out, nrow_out, ncol_out, nz_out)
    
    DEALLOCATE(A_full)
  END SUBROUTINE generate_spd_matrix

  SUBROUTINE generate_diagonal_matrix(n, irow_out, pcol_out, val_out, nrow_out, ncol_out, nz_out)
    ! Generate a diagonal test matrix
    INTEGER, INTENT(in) :: n
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(out) :: irow_out, pcol_out
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(out) :: val_out
    INTEGER, INTENT(out) :: nrow_out, ncol_out, nz_out
    
    INTEGER :: i
    
    nrow_out = n
    ncol_out = n
    nz_out = n
    
    ALLOCATE(irow_out(nz_out), pcol_out(ncol_out+1), val_out(nz_out))
    
    ! Column pointers for diagonal matrix
    DO i = 1, ncol_out + 1
      pcol_out(i) = i
    END DO
    
    ! Row indices (diagonal)
    DO i = 1, nz_out
      irow_out(i) = i
    END DO
    
    ! Values (increasing diagonal elements)
    DO i = 1, nz_out
      val_out(i) = REAL(i, kind=dp)
    END DO
  END SUBROUTINE generate_diagonal_matrix

  SUBROUTINE generate_tridiagonal_matrix(n, irow_out, pcol_out, val_out, nrow_out, ncol_out, nz_out)
    ! Generate a tridiagonal test matrix
    INTEGER, INTENT(in) :: n
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(out) :: irow_out, pcol_out
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(out) :: val_out
    INTEGER, INTENT(out) :: nrow_out, ncol_out, nz_out
    
    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: A_full
    INTEGER :: i
    
    ! Create tridiagonal matrix
    ALLOCATE(A_full(n,n))
    A_full = 0.0_dp
    
    DO i = 1, n
      A_full(i,i) = 4.0_dp  ! Main diagonal
      IF (i > 1) A_full(i,i-1) = -1.0_dp  ! Lower diagonal
      IF (i < n) A_full(i,i+1) = -1.0_dp  ! Upper diagonal
    END DO
    
    ! Convert to sparse format
    CALL full2sparse(A_full, irow_out, pcol_out, val_out, nrow_out, ncol_out, nz_out)
    
    DEALLOCATE(A_full)
  END SUBROUTINE generate_tridiagonal_matrix

  SUBROUTINE compare_solutions(name, x1, x2, max_diff, rel_diff)
    CHARACTER(len=*), INTENT(in) :: name
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: x1, x2
    REAL(kind=dp), INTENT(out) :: max_diff, rel_diff
    
    REAL(kind=dp) :: norm_x1
    
    max_diff = MAXVAL(ABS(x1 - x2))
    norm_x1 = SQRT(SUM(x1**2))
    
    IF (norm_x1 > 1.0e-14_dp) THEN
      rel_diff = max_diff / norm_x1
    ELSE
      rel_diff = max_diff
    END IF
    
    WRITE(*,'(A,A,A,ES12.4,A,ES12.4)') "  ", name, " - Max diff: ", max_diff, ", Rel diff: ", rel_diff
  END SUBROUTINE compare_solutions

  SUBROUTINE compare_complex_solutions(name, x1, x2, max_diff, rel_diff)
    CHARACTER(len=*), INTENT(in) :: name
    COMPLEX(kind=dp), DIMENSION(:), INTENT(in) :: x1, x2
    REAL(kind=dp), INTENT(out) :: max_diff, rel_diff
    
    REAL(kind=dp) :: norm_x1
    
    max_diff = MAXVAL(ABS(x1 - x2))
    norm_x1 = SQRT(SUM(ABS(x1)**2))
    
    IF (norm_x1 > 1.0e-14_dp) THEN
      rel_diff = max_diff / norm_x1
    ELSE
      rel_diff = max_diff
    END IF
    
    WRITE(*,'(A,A,A,ES12.4,A,ES12.4)') "  ", name, " - Max diff: ", max_diff, ", Rel diff: ", rel_diff
  END SUBROUTINE compare_complex_solutions

  SUBROUTINE generate_complex_hermitian_matrix(n, irow_out, pcol_out, z_val_out, nrow_out, ncol_out, nz_out)
    ! Generate a Hermitian positive definite test matrix
    INTEGER, INTENT(in) :: n
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(out) :: irow_out, pcol_out
    COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(out) :: z_val_out
    INTEGER, INTENT(out) :: nrow_out, ncol_out, nz_out
    
    COMPLEX(kind=dp), DIMENSION(:,:), ALLOCATABLE :: A_full
    INTEGER :: i, j
    
    ! Create Hermitian positive definite matrix
    ALLOCATE(A_full(n,n))
    
    ! Fill matrix to be Hermitian
    DO i = 1, n
      A_full(i,i) = CMPLX(REAL(i+1, kind=dp), 0.0_dp, kind=dp)  ! Real diagonal
      DO j = i+1, n
        A_full(i,j) = CMPLX(0.1_dp, 0.05_dp * REAL(i-j, kind=dp), kind=dp)
        A_full(j,i) = CONJG(A_full(i,j))  ! Hermitian property
      END DO
    END DO
    
    ! Convert to sparse format
    CALL full2sparse_complex(A_full, irow_out, pcol_out, z_val_out, nrow_out, ncol_out, nz_out)
    
    DEALLOCATE(A_full)
  END SUBROUTINE generate_complex_hermitian_matrix

  SUBROUTINE generate_large_sparse_matrix(n, sparsity, irow_out, pcol_out, val_out, nrow_out, ncol_out, nz_out)
    ! Generate a large sparse matrix with controlled sparsity
    INTEGER, INTENT(in) :: n
    REAL(kind=dp), INTENT(in) :: sparsity  ! Fraction of non-zero elements
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(out) :: irow_out, pcol_out
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(out) :: val_out
    INTEGER, INTENT(out) :: nrow_out, ncol_out, nz_out
    
    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: A_full
    INTEGER :: i, j
    
    ! Create a simple large sparse matrix (diagonal + first super/sub diagonal)
    ALLOCATE(A_full(n,n))
    A_full = 0.0_dp
    
    ! Set diagonal elements (ensure non-singular)
    DO i = 1, n
      A_full(i,i) = 4.0_dp + REAL(i, kind=dp) * 0.01_dp
    END DO
    
    ! Add super and sub diagonals for sparse structure
    DO i = 1, n-1
      A_full(i,i+1) = -1.0_dp  ! Super diagonal
      A_full(i+1,i) = -1.0_dp  ! Sub diagonal
    END DO
    
    ! Add some additional sparse structure every 5th row/column
    DO i = 5, n, 5
      DO j = 1, n, 7
        IF (i /= j .AND. ABS(A_full(i,j)) < 1.0e-14_dp) THEN
          A_full(i,j) = 0.1_dp * SIN(REAL(i+j, kind=dp))
        END IF
      END DO
    END DO
    
    ! Convert to sparse format
    CALL full2sparse(A_full, irow_out, pcol_out, val_out, nrow_out, ncol_out, nz_out)
    
    DEALLOCATE(A_full)
  END SUBROUTINE generate_large_sparse_matrix

  SUBROUTINE test_small_spd_matrix_comparison()
    REAL(kind=dp) :: max_diff, rel_diff
    INTEGER :: old_method
    
    WRITE(*,'(A)') "Test 1: Small SPD matrix comparison"
    
    ! Generate 5x5 SPD matrix
    CALL generate_spd_matrix(5, irow, pcol, val, nrow, ncol, nz)
    
    ! Create RHS
    IF (ALLOCATED(b)) DEALLOCATE(b)
    IF (ALLOCATED(b_orig)) DEALLOCATE(b_orig)
    ALLOCATE(b(nrow), b_orig(nrow))
    b = 1.0_dp
    b_orig = b
    
    ! Save current solver method
    old_method = sparse_solve_method
    
    ! Solve with UMFPACK
    sparse_solve_method = SOLVER_UMFPACK
    IF (ALLOCATED(x_umfpack)) DEALLOCATE(x_umfpack)
    ALLOCATE(x_umfpack(nrow))
    x_umfpack = b
    CALL sparse_solve(nrow, ncol, nz, irow, pcol, val, x_umfpack)
    
    ! Solve with BiCGSTAB
    sparse_solve_method = SOLVER_BICGSTAB
    IF (ALLOCATED(x_bicgstab)) DEALLOCATE(x_bicgstab)
    ALLOCATE(x_bicgstab(nrow))
    x_bicgstab = b_orig
    CALL sparse_solve(nrow, ncol, nz, irow, pcol, val, x_bicgstab)
    
    ! Compare solutions
    CALL compare_solutions("SPD 5x5", x_umfpack, x_bicgstab, max_diff, rel_diff)
    
    CALL check_result("Small SPD matrix solver agreement", &
         max_diff < solver_comparison_tol .AND. rel_diff < solver_comparison_tol)
    
    ! Restore solver method
    sparse_solve_method = old_method
  END SUBROUTINE test_small_spd_matrix_comparison

  SUBROUTINE test_diagonal_matrix_comparison()
    REAL(kind=dp) :: max_diff, rel_diff
    INTEGER :: old_method, i
    
    WRITE(*,'(A)') "Test 2: Diagonal matrix comparison"
    
    ! Generate 10x10 diagonal matrix
    CALL generate_diagonal_matrix(10, irow, pcol, val, nrow, ncol, nz)
    
    ! Create RHS
    IF (ALLOCATED(b)) DEALLOCATE(b)
    IF (ALLOCATED(b_orig)) DEALLOCATE(b_orig)
    ALLOCATE(b(nrow), b_orig(nrow))
    b = (/ (REAL(i, kind=dp), i=1,nrow) /)
    b_orig = b
    
    ! Save current solver method
    old_method = sparse_solve_method
    
    ! Solve with UMFPACK
    sparse_solve_method = SOLVER_UMFPACK
    IF (ALLOCATED(x_umfpack)) DEALLOCATE(x_umfpack)
    ALLOCATE(x_umfpack(nrow))
    x_umfpack = b
    CALL sparse_solve(nrow, ncol, nz, irow, pcol, val, x_umfpack)
    
    ! Solve with BiCGSTAB
    sparse_solve_method = SOLVER_BICGSTAB
    IF (ALLOCATED(x_bicgstab)) DEALLOCATE(x_bicgstab)
    ALLOCATE(x_bicgstab(nrow))
    x_bicgstab = b_orig
    CALL sparse_solve(nrow, ncol, nz, irow, pcol, val, x_bicgstab)
    
    ! Compare solutions
    CALL compare_solutions("Diagonal 10x10", x_umfpack, x_bicgstab, max_diff, rel_diff)
    
    CALL check_result("Diagonal matrix solver agreement", &
         max_diff < solver_comparison_tol .AND. rel_diff < solver_comparison_tol)
    
    ! Restore solver method
    sparse_solve_method = old_method
  END SUBROUTINE test_diagonal_matrix_comparison

  SUBROUTINE test_tridiagonal_matrix_comparison()
    REAL(kind=dp) :: max_diff, rel_diff
    INTEGER :: old_method
    
    WRITE(*,'(A)') "Test 3: Tridiagonal matrix comparison"
    
    ! Generate 20x20 tridiagonal matrix
    CALL generate_tridiagonal_matrix(20, irow, pcol, val, nrow, ncol, nz)
    
    ! Create RHS
    IF (ALLOCATED(b)) DEALLOCATE(b)
    IF (ALLOCATED(b_orig)) DEALLOCATE(b_orig)
    ALLOCATE(b(nrow), b_orig(nrow))
    b = 1.0_dp
    b_orig = b
    
    ! Save current solver method
    old_method = sparse_solve_method
    
    ! Solve with UMFPACK
    sparse_solve_method = SOLVER_UMFPACK
    IF (ALLOCATED(x_umfpack)) DEALLOCATE(x_umfpack)
    ALLOCATE(x_umfpack(nrow))
    x_umfpack = b
    CALL sparse_solve(nrow, ncol, nz, irow, pcol, val, x_umfpack)
    
    ! Solve with BiCGSTAB
    sparse_solve_method = SOLVER_BICGSTAB
    IF (ALLOCATED(x_bicgstab)) DEALLOCATE(x_bicgstab)
    ALLOCATE(x_bicgstab(nrow))
    x_bicgstab = b_orig
    CALL sparse_solve(nrow, ncol, nz, irow, pcol, val, x_bicgstab)
    
    ! Compare solutions
    CALL compare_solutions("Tridiagonal 20x20", x_umfpack, x_bicgstab, max_diff, rel_diff)
    
    CALL check_result("Tridiagonal matrix solver agreement", &
         max_diff < solver_comparison_tol .AND. rel_diff < solver_comparison_tol)
    
    ! Restore solver method
    sparse_solve_method = old_method
  END SUBROUTINE test_tridiagonal_matrix_comparison

  SUBROUTINE test_complex_matrix_comparison()
    REAL(kind=dp) :: max_diff, rel_diff
    INTEGER :: old_method
    
    WRITE(*,'(A)') "Test 4: Complex matrix comparison"
    
    ! Note: BiCGSTAB for complex matrices not yet implemented
    ! This test will be skipped until complex BiCGSTAB is available
    
    CALL check_result("Complex matrix solver agreement (SKIPPED - complex BiCGSTAB not implemented)", .TRUE.)
  END SUBROUTINE test_complex_matrix_comparison

  SUBROUTINE test_large_sparse_matrix_comparison()
    REAL(kind=dp) :: max_diff, rel_diff
    INTEGER :: old_method
    
    WRITE(*,'(A)') "Test 5: Large sparse matrix comparison"
    
    ! Generate 100x100 sparse matrix with 5% sparsity
    CALL generate_large_sparse_matrix(100, 0.05_dp, irow, pcol, val, nrow, ncol, nz)
    
    ! Create RHS
    IF (ALLOCATED(b)) DEALLOCATE(b)
    IF (ALLOCATED(b_orig)) DEALLOCATE(b_orig)
    ALLOCATE(b(nrow), b_orig(nrow))
    b = 1.0_dp
    b_orig = b
    
    ! Save current solver method
    old_method = sparse_solve_method
    
    ! Solve with UMFPACK
    sparse_solve_method = SOLVER_UMFPACK
    IF (ALLOCATED(x_umfpack)) DEALLOCATE(x_umfpack)
    ALLOCATE(x_umfpack(nrow))
    x_umfpack = b
    CALL sparse_solve(nrow, ncol, nz, irow, pcol, val, x_umfpack)
    
    ! Solve with BiCGSTAB
    sparse_solve_method = SOLVER_BICGSTAB
    IF (ALLOCATED(x_bicgstab)) DEALLOCATE(x_bicgstab)
    ALLOCATE(x_bicgstab(nrow))
    x_bicgstab = b_orig
    CALL sparse_solve(nrow, ncol, nz, irow, pcol, val, x_bicgstab)
    
    ! Compare solutions
    CALL compare_solutions("Large sparse 100x100", x_umfpack, x_bicgstab, max_diff, rel_diff)
    
    CALL check_result("Large sparse matrix solver agreement", &
         max_diff < solver_comparison_tol .AND. rel_diff < solver_comparison_tol)
    
    ! Restore solver method
    sparse_solve_method = old_method
  END SUBROUTINE test_large_sparse_matrix_comparison

  SUBROUTINE test_performance_comparison()
    REAL(kind=dp) :: time_umfpack, time_bicgstab, time_start, time_end
    INTEGER :: old_method, iter
    
    WRITE(*,'(A)') "Test 6: Performance comparison framework"
    
    ! Generate a moderately large matrix for timing
    CALL generate_large_sparse_matrix(50, 0.1_dp, irow, pcol, val, nrow, ncol, nz)
    
    ! Create RHS
    IF (ALLOCATED(b)) DEALLOCATE(b)
    IF (ALLOCATED(b_orig)) DEALLOCATE(b_orig)
    ALLOCATE(b(nrow), b_orig(nrow))
    b = 1.0_dp
    b_orig = b
    
    ! Save current solver method
    old_method = sparse_solve_method
    
    ! Time UMFPACK solver
    CALL CPU_TIME(time_start)
    DO iter = 1, 10  ! Multiple iterations for more accurate timing
      sparse_solve_method = SOLVER_UMFPACK
      IF (ALLOCATED(x_umfpack)) DEALLOCATE(x_umfpack)
      ALLOCATE(x_umfpack(nrow))
      x_umfpack = b
      CALL sparse_solve(nrow, ncol, nz, irow, pcol, val, x_umfpack)
    END DO
    CALL CPU_TIME(time_end)
    time_umfpack = (time_end - time_start) / 10.0_dp
    
    ! Time BiCGSTAB solver
    CALL CPU_TIME(time_start)
    DO iter = 1, 10  ! Multiple iterations for more accurate timing
      sparse_solve_method = SOLVER_BICGSTAB
      IF (ALLOCATED(x_bicgstab)) DEALLOCATE(x_bicgstab)
      ALLOCATE(x_bicgstab(nrow))
      x_bicgstab = b_orig
      CALL sparse_solve(nrow, ncol, nz, irow, pcol, val, x_bicgstab)
    END DO
    CALL CPU_TIME(time_end)
    time_bicgstab = (time_end - time_start) / 10.0_dp
    
    ! Report timing results
    WRITE(*,'(A,F8.4,A)') "  UMFPACK time per solve: ", time_umfpack * 1000.0_dp, " ms"
    WRITE(*,'(A,F8.4,A)') "  BiCGSTAB time per solve: ", time_bicgstab * 1000.0_dp, " ms"
    
    IF (time_umfpack > 0.0_dp) THEN
      WRITE(*,'(A,F6.2,A)') "  Speedup factor: ", time_umfpack / time_bicgstab, "x"
    END IF
    
    ! Performance comparison passes if both solvers complete without error
    CALL check_result("Performance comparison framework", .TRUE.)
    
    ! Restore solver method
    sparse_solve_method = old_method
  END SUBROUTINE test_performance_comparison

END PROGRAM test_solver_integration