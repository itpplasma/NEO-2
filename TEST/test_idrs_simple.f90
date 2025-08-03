PROGRAM test_idrs_simple
  !> Test simple IDR(s) implementation
  
  USE nrtype, ONLY: I4B, DP
  USE idrs_simple_mod, ONLY: idrs_simple
  IMPLICIT NONE
  
  ! Test parameters
  INTEGER(I4B), PARAMETER :: n = 10
  INTEGER(I4B) :: i, j, nnz
  
  ! Sparse matrix in CSR format
  INTEGER(I4B), ALLOCATABLE :: row_ptr(:), col_idx(:)
  REAL(DP), ALLOCATABLE :: values(:)
  
  ! Vectors
  REAL(DP) :: b(n), x(n), x_exact(n), x0(n)
  REAL(DP) :: error, residual(n), normr
  
  ! IDR(s) parameters
  INTEGER :: max_iter, s, iter
  REAL(DP) :: tol
  LOGICAL :: converged
  
  WRITE(*,'(A)') '==========================================='
  WRITE(*,'(A)') 'Testing Simple IDR(s) Implementation'
  WRITE(*,'(A)') '==========================================='
  
  ! Create a simple tridiagonal matrix
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
  
  ! Test with different s values
  DO s = 1, 4
    WRITE(*,*)
    WRITE(*,'(A,I2)') 'Testing with s = ', s
    WRITE(*,'(A)') '--------------------------'
    
    ! IDR(s) parameters
    max_iter = 50
    tol = 1.0E-10_DP
    x0 = 0.0_DP
    
    ! Call IDR(s)
    CALL idrs_simple(n, row_ptr, col_idx, values, b, x0, &
                     s, max_iter, tol, x, iter, converged)
    
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
    normr = SQRT(DOT_PRODUCT(residual, residual))
    
    WRITE(*,'(A,L1)') '  Converged: ', converged
    WRITE(*,'(A,I5)') '  Iterations: ', iter
    WRITE(*,'(A,ES12.5)') '  Residual norm: ', normr
    WRITE(*,'(A,ES12.5)') '  Solution error: ', error
    
    IF (converged .AND. error < 1.0E-8_DP) THEN
      WRITE(*,'(A)') '  Status: SUCCESS'
    ELSE
      WRITE(*,'(A)') '  Status: FAILURE'
    END IF
    
  END DO
  
  ! Clean up
  DEALLOCATE(row_ptr, col_idx, values)
  
  WRITE(*,*)
  WRITE(*,'(A)') '==========================================='
  
END PROGRAM test_idrs_simple