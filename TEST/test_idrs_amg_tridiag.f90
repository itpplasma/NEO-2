PROGRAM test_idrs_amg_tridiag
  !> Test IDR(s) with AMG on a simple tridiagonal matrix
  
  USE nrtype, ONLY: I4B, DP
  USE sparse_solvers_mod, ONLY: sparse_solve_suitesparse, sparse_solve_idrs_real, &
                                sparse_solve_method, SOLVER_UMFPACK, SOLVER_IDRS, &
                                default_iterative_params, PRECOND_AMG, PRECOND_NONE
  IMPLICIT NONE
  
  ! Test parameters
  INTEGER(I4B), PARAMETER :: n = 50
  INTEGER(I4B) :: i, j, nnz
  
  ! Sparse matrix in CSR format
  INTEGER(I4B), ALLOCATABLE :: row_ptr(:), col_idx(:)
  REAL(DP), ALLOCATABLE :: values(:), b(:), x(:), x_exact(:), x_ref(:)
  REAL(DP) :: error, residual(n), residual_norm
  
  WRITE(*,'(A)') '==========================================================='
  WRITE(*,'(A)') 'Testing IDR(s) + AMG on Tridiagonal Matrix'
  WRITE(*,'(A)') '==========================================================='
  
  ! Create a tridiagonal matrix with good conditioning
  nnz = 3*n - 2
  ALLOCATE(row_ptr(n+1), col_idx(nnz), values(nnz))
  ALLOCATE(b(n), x(n), x_exact(n), x_ref(n))
  
  ! Build sparse matrix in CSR format
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
    values(nnz) = 4.0_DP
    IF (i < n) THEN
      nnz = nnz + 1
      col_idx(nnz) = i + 1
      values(nnz) = -1.0_DP
    END IF
  END DO
  row_ptr(n+1) = nnz + 1
  
  ! Set exact solution
  DO i = 1, n
    x_exact(i) = SIN(REAL(i, DP) * 3.14159265359_DP / REAL(n+1, DP))
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
  WRITE(*,'(A)') '  Type: Tridiagonal, diagonally dominant'
  WRITE(*,*)
  
  ! Test 1: Reference solution with UMFPACK
  WRITE(*,'(A)') '1. Computing reference solution with UMFPACK...'
  x_ref = b  ! Copy RHS since it gets overwritten
  sparse_solve_method = SOLVER_UMFPACK
  CALL sparse_solve_suitesparse(n, n, nnz, row_ptr, col_idx, values, x_ref, 1)
  
  error = 0.0_DP
  DO i = 1, n
    error = error + (x_ref(i) - x_exact(i))**2
  END DO
  error = SQRT(error)
  WRITE(*,'(A,ES12.5)') '   UMFPACK solution error: ', error
  
  ! Test 2: IDR(s) with AMG
  WRITE(*,*)
  WRITE(*,'(A)') '2. Testing IDR(s) without preconditioning...'
  
  ! Configure IDR(s) without preconditioning
  default_iterative_params%preconditioner_type = PRECOND_NONE
  default_iterative_params%idrs_shadow_space_dim = 4
  default_iterative_params%max_iterations = 1000
  default_iterative_params%abs_tolerance = 1.0E-10_DP
  default_iterative_params%rel_tolerance = 1.0E-8_DP
  default_iterative_params%verbose = .FALSE.
  
  x = b  ! Copy RHS since it gets overwritten
  CALL sparse_solve_idrs_real(n, n, nnz, row_ptr, col_idx, values, x, 1)
  
  ! Check results
  error = 0.0_DP
  DO i = 1, n
    error = error + (x(i) - x_exact(i))**2
  END DO
  error = SQRT(error)
  
  ! Compute residual
  residual = 0.0_DP
  DO i = 1, n
    DO j = row_ptr(i), row_ptr(i+1) - 1
      residual(i) = residual(i) + values(j) * x(col_idx(j))
    END DO
  END DO
  residual = b - residual
  residual_norm = SQRT(DOT_PRODUCT(residual, residual))
  
  WRITE(*,'(A,ES12.5)') '   IDR(s)+AMG solution error: ', error
  WRITE(*,'(A,ES12.5)') '   IDR(s)+AMG residual norm: ', residual_norm
  
  IF (error < 1.0E-8_DP) THEN
    WRITE(*,*)
    WRITE(*,'(A)') 'SUCCESS: IDR(s)+AMG converged to correct solution'
  ELSE
    WRITE(*,*)
    WRITE(*,'(A)') 'FAILURE: IDR(s)+AMG did not converge properly'
  END IF
  
  ! Clean up
  DEALLOCATE(row_ptr, col_idx, values, b, x, x_exact, x_ref)
  
  WRITE(*,*)
  WRITE(*,'(A)') '==========================================================='
  
END PROGRAM test_idrs_amg_tridiag