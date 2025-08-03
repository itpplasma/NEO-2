PROGRAM test_idrs_no_precond
  !> Test IDR(s) without preconditioning on a simple well-conditioned matrix
  
  USE nrtype, ONLY: I4B, DP
  USE idrs_mod, ONLY: idrs_solve_amg_preconditioned
  USE amg_types_mod, ONLY: amg_hierarchy
  IMPLICIT NONE
  
  ! Test parameters
  INTEGER(I4B), PARAMETER :: n = 10
  INTEGER(I4B) :: i, j, nnz
  
  ! Sparse matrix in CSR format
  INTEGER(I4B), ALLOCATABLE :: row_ptr(:), col_idx(:)
  REAL(DP), ALLOCATABLE :: values(:)
  
  ! Vectors
  REAL(DP) :: b(n), x(n), x_exact(n), x_initial(n)
  REAL(DP) :: error, residual_norm
  
  ! IDR(s) parameters
  INTEGER :: max_iter, shadow_dim, iter, info
  REAL(DP) :: tol
  LOGICAL :: converged
  TYPE(amg_hierarchy) :: dummy_amg
  INTEGER :: workspace
  
  WRITE(*,'(A)') '==========================================================='
  WRITE(*,'(A)') 'Testing IDR(s) on Simple Tridiagonal Matrix (No Precond)'
  WRITE(*,'(A)') '==========================================================='
  
  ! Create a simple tridiagonal matrix
  ! A = [ 2 -1  0  0 ...]
  !     [-1  2 -1  0 ...]
  !     [ 0 -1  2 -1 ...]
  !     [... ]
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
  
  ! Set exact solution
  DO i = 1, n
    x_exact(i) = REAL(i, DP)
  END DO
  
  ! Compute RHS: b = A * x_exact
  b = 0.0_DP
  DO i = 1, n
    DO j = row_ptr(i), row_ptr(i+1) - 1
      b(i) = b(i) + values(j) * x_exact(col_idx(j))
    END DO
  END DO
  
  WRITE(*,'(A)') 'Matrix properties:'
  WRITE(*,'(A,I5)') '  Size: ', n
  WRITE(*,'(A,I5)') '  Nonzeros: ', nnz
  WRITE(*,'(A,F8.3)') '  Condition number (approx): ', 4.0_DP
  WRITE(*,*)
  
  ! Test without preconditioning (identity preconditioner)
  ! Initialize dummy AMG hierarchy that acts as identity
  dummy_amg%n_levels = 1
  ALLOCATE(dummy_amg%levels(1))
  dummy_amg%levels(1)%n = n
  dummy_amg%levels(1)%nnz = n
  ALLOCATE(dummy_amg%levels(1)%row_ptr(n+1))
  ALLOCATE(dummy_amg%levels(1)%col_idx(n))
  ALLOCATE(dummy_amg%levels(1)%values(n))
  DO i = 1, n
    dummy_amg%levels(1)%row_ptr(i) = i
    dummy_amg%levels(1)%col_idx(i) = i
    dummy_amg%levels(1)%values(i) = 1.0_DP
  END DO
  dummy_amg%levels(1)%row_ptr(n+1) = n + 1
  
  ! IDR(s) parameters
  max_iter = 50
  tol = 1.0E-10_DP
  shadow_dim = 4
  x_initial = 0.0_DP
  
  WRITE(*,'(A)') 'Testing IDR(s) with identity preconditioner...'
  WRITE(*,'(A,I3)') '  Shadow dimension s = ', shadow_dim
  WRITE(*,'(A,I5)') '  Max iterations = ', max_iter
  WRITE(*,'(A,ES12.5)') '  Tolerance = ', tol
  WRITE(*,*)
  
  ! Call IDR(s)
  CALL idrs_solve_amg_preconditioned(workspace, n, row_ptr, col_idx, values, &
                                    b, x_initial, max_iter, tol, dummy_amg, &
                                    shadow_dim, x, iter, residual_norm, &
                                    converged, info)
  
  ! Check results
  error = 0.0_DP
  DO i = 1, n
    error = error + (x(i) - x_exact(i))**2
  END DO
  error = SQRT(error)
  
  WRITE(*,*)
  WRITE(*,'(A)') 'Results:'
  WRITE(*,'(A,L1)') '  Converged: ', converged
  WRITE(*,'(A,I5)') '  Iterations: ', iter
  WRITE(*,'(A,ES12.5)') '  Final residual: ', residual_norm
  WRITE(*,'(A,ES12.5)') '  Solution error: ', error
  WRITE(*,'(A,I3)') '  Info code: ', info
  
  IF (converged .AND. error < 1.0E-8_DP) THEN
    WRITE(*,*)
    WRITE(*,'(A)') 'SUCCESS: IDR(s) converged to correct solution'
  ELSE
    WRITE(*,*)
    WRITE(*,'(A)') 'FAILURE: IDR(s) did not converge properly'
    WRITE(*,*)
    WRITE(*,'(A)') 'Expected vs Computed solution:'
    DO i = 1, n
      WRITE(*,'(I3,2F12.6)') i, x_exact(i), x(i)
    END DO
  END IF
  
  ! Clean up
  DEALLOCATE(row_ptr, col_idx, values)
  DEALLOCATE(dummy_amg%levels(1)%row_ptr)
  DEALLOCATE(dummy_amg%levels(1)%col_idx)
  DEALLOCATE(dummy_amg%levels(1)%values)
  DEALLOCATE(dummy_amg%levels)
  
END PROGRAM test_idrs_no_precond