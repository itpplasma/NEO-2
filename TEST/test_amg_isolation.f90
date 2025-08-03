PROGRAM test_amg_isolation
  !> Test AMG preconditioner in complete isolation to debug issues

  USE nrtype, ONLY: I4B, DP
  USE amg_types_mod, ONLY: amg_hierarchy
  USE amg_precond_mod, ONLY: amg_precond_setup, amg_precond_apply, amg_precond_destroy
  IMPLICIT NONE

  ! Test parameters - use the exact same 20-point problem where IDR(s) works
  INTEGER(I4B), PARAMETER :: n = 20
  INTEGER(I4B) :: i, j, nnz, info
  
  ! Sparse matrix in CSR format (same as IDR(s) working case)
  INTEGER(I4B), ALLOCATABLE :: row_ptr(:), col_idx(:)
  REAL(DP), ALLOCATABLE :: values(:)
  
  ! Vectors for testing
  REAL(DP) :: b(n), x(n), x_amg(n), x_exact(n)
  REAL(DP) :: residual_before, residual_after, reduction_factor
  TYPE(amg_hierarchy) :: amg_hier
  
  WRITE(*,'(A)') '============================================='
  WRITE(*,'(A)') 'Testing AMG Preconditioner in Isolation'
  WRITE(*,'(A)') '============================================='
  WRITE(*,*)
  
  ! Create same simple tridiagonal matrix that IDR(s) solves easily
  ! A = [ 2 -1  0  0 ...]
  !     [-1  2 -1  0 ...]
  !     [ 0 -1  2 -1 ...]
  nnz = 3*n - 2
  ALLOCATE(row_ptr(n+1), col_idx(nnz), values(nnz))
  
  ! Build CSR format
  nnz = 0
  DO i = 1, n
    row_ptr(i) = nnz + 1
    IF (i > 1) THEN
      nnz = nnz + 1
      col_idx(nnz) = i - 1
      values(nnz) = -1.0_DP
    END IF
    nnz = nnz + 1
    col_idx(nnz) = i
    values(nnz) = 2.0_DP
    IF (i < n) THEN
      nnz = nnz + 1
      col_idx(nnz) = i + 1
      values(nnz) = -1.0_DP
    END IF
  END DO
  row_ptr(n+1) = nnz + 1
  
  WRITE(*,'(A,I0)') 'Matrix size: ', n
  WRITE(*,'(A,I0)') 'Nonzeros: ', nnz
  WRITE(*,*)
  
  ! Set up exact solution and RHS
  DO i = 1, n
    x_exact(i) = REAL(i, DP)  ! Simple test solution
  END DO
  
  ! Compute RHS: b = A * x_exact
  b = 0.0_DP
  DO i = 1, n
    DO j = row_ptr(i), row_ptr(i+1) - 1
      b(i) = b(i) + values(j) * x_exact(col_idx(j))
    END DO
  END DO
  
  WRITE(*,'(A)') '1. Testing AMG setup...'
  
  ! Set up AMG hierarchy
  CALL amg_precond_setup(amg_hier, n, nnz, row_ptr, col_idx, values)
  
  ! Check if AMG was set up correctly
  IF (.NOT. ALLOCATED(amg_hier%levels) .OR. amg_hier%n_levels < 1) THEN
    WRITE(*,'(A)') '   [FAILURE] AMG hierarchy not allocated properly'
    STOP 1
  END IF
  
  IF (amg_hier%levels(1)%n /= n) THEN
    WRITE(*,'(A)') '   [FAILURE] AMG hierarchy dimension mismatch'
    WRITE(*,'(A,I0,A,I0)') '   Expected: ', n, ', Got: ', amg_hier%levels(1)%n
    STOP 1
  END IF
  
  WRITE(*,'(A)') '   [SUCCESS] AMG hierarchy built'
  WRITE(*,'(A,I0)') '   Number of levels: ', amg_hier%n_levels
  WRITE(*,*)
  
  ! Test 2: AMG preconditioning effectiveness
  WRITE(*,'(A)') '2. Testing AMG preconditioning effectiveness...'
  
  ! Start with zero initial guess
  x = 0.0_DP
  
  ! Compute initial residual: r = b - A*x = b (since x=0)
  residual_before = SQRT(DOT_PRODUCT(b, b))
  WRITE(*,'(A,ES12.5)') '   Initial residual norm: ', residual_before
  
  ! Apply AMG preconditioner: x_amg = M^(-1) * b
  x_amg = b  ! Copy RHS
  CALL amg_precond_apply(amg_hier, x_amg, info)
  
  IF (info /= 0) THEN
    WRITE(*,'(A,I0)') '   [FAILURE] AMG apply failed with info = ', info
    STOP 1
  END IF
  
  ! Check for NaN/Inf in AMG result
  IF (ANY(x_amg /= x_amg)) THEN
    WRITE(*,'(A)') '   [FAILURE] AMG produced NaN values'
    DO i = 1, n
      IF (x_amg(i) /= x_amg(i)) THEN
        WRITE(*,'(A,I0,A,ES12.5)') '   x_amg(', i, ') = ', x_amg(i)
      END IF
    END DO
    STOP 1
  END IF
  
  ! Compute residual after AMG: r_new = b - A*(x_amg)
  x = 0.0_DP  ! Compute A*x_amg
  DO i = 1, n
    DO j = row_ptr(i), row_ptr(i+1) - 1
      x(i) = x(i) + values(j) * x_amg(col_idx(j))
    END DO
  END DO
  
  x = b - x  ! r_new = b - A*x_amg
  residual_after = SQRT(DOT_PRODUCT(x, x))
  
  reduction_factor = residual_before / residual_after
  
  WRITE(*,'(A,ES12.5)') '   Residual after AMG: ', residual_after
  WRITE(*,'(A,ES12.5)') '   Reduction factor: ', reduction_factor
  
  IF (reduction_factor < 2.0_DP) THEN
    WRITE(*,'(A)') '   [WARNING] AMG not very effective (reduction < 2x)'
  ELSE IF (reduction_factor > 1000.0_DP) THEN
    WRITE(*,'(A)') '   [WARNING] AMG too aggressive (reduction > 1000x)'
    WRITE(*,'(A)') '            This might cause numerical issues in IDR(s)'
  ELSE
    WRITE(*,'(A)') '   [SUCCESS] AMG reduction factor looks reasonable'
  END IF
  
  WRITE(*,*)
  
  ! Test 3: Check AMG solution quality
  WRITE(*,'(A)') '3. Testing AMG solution quality...'
  
  ! Compute error vs exact solution
  x = x_amg - x_exact
  REAL(DP) :: error_norm
  error_norm = SQRT(DOT_PRODUCT(x, x))
  
  WRITE(*,'(A,ES12.5)') '   Solution error norm: ', error_norm
  
  IF (error_norm < 1.0e-2_DP) THEN
    WRITE(*,'(A)') '   [SUCCESS] AMG provides good approximate solution'
  ELSE
    WRITE(*,'(A)') '   [INFO] AMG approximate solution (expected for preconditioner)'
  END IF
  
  WRITE(*,*)
  WRITE(*,'(A)') '4. AMG preconditioner values summary:'
  WRITE(*,'(A,ES12.5)') '   Max AMG output: ', MAXVAL(ABS(x_amg))
  WRITE(*,'(A,ES12.5)') '   Min AMG output: ', MINVAL(ABS(x_amg))
  WRITE(*,'(A,ES12.5)') '   AMG output norm: ', SQRT(DOT_PRODUCT(x_amg, x_amg))
  
  ! Clean up
  CALL amg_precond_destroy(amg_hier)
  DEALLOCATE(row_ptr, col_idx, values)
  
  WRITE(*,*)
  WRITE(*,'(A)') '============================================='
  WRITE(*,'(A)') 'AMG isolation test completed'
  WRITE(*,'(A)') '============================================='

END PROGRAM test_amg_isolation