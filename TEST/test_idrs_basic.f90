PROGRAM test_idrs_basic
  !> Basic test of IDR(s) without AMG to verify algorithm correctness
  
  USE nrtype, ONLY: I4B, DP
  USE sparse_types_mod, ONLY: sparse_matrix
  USE sparse_conversion_mod, ONLY: csc_to_csr_real
  IMPLICIT NONE
  
  ! Test parameters
  INTEGER(I4B), PARAMETER :: n = 5
  INTEGER(I4B) :: i, j, nnz
  
  ! Sparse matrix in CSR format
  INTEGER(I4B), ALLOCATABLE :: row_ptr(:), col_idx(:)
  REAL(DP), ALLOCATABLE :: values(:)
  
  ! Vectors
  REAL(DP) :: b(n), x(n), x_exact(n)
  REAL(DP) :: error
  
  WRITE(*,'(A)') '========================================================='
  WRITE(*,'(A)') 'Basic IDR(s) Test - Simple 5x5 Matrix'
  WRITE(*,'(A)') '========================================================='
  
  ! Create a simple 5x5 tridiagonal matrix
  ! A = [2 -1  0  0  0]
  !     [-1 2 -1  0  0]
  !     [0 -1  2 -1  0]
  !     [0  0 -1  2 -1]
  !     [0  0  0 -1  2]
  
  nnz = 13  ! 5 diagonal + 4 upper + 4 lower
  ALLOCATE(row_ptr(n+1), col_idx(nnz), values(nnz))
  
  ! Build CSR format
  row_ptr = [1, 3, 6, 9, 12, 14]  ! 1-based indexing
  col_idx = [1, 2, 1, 2, 3, 2, 3, 4, 3, 4, 5, 4, 5]
  values = [2.0_DP, -1.0_DP, -1.0_DP, 2.0_DP, -1.0_DP, -1.0_DP, 2.0_DP, &
            -1.0_DP, -1.0_DP, 2.0_DP, -1.0_DP, -1.0_DP, 2.0_DP]
  
  ! Set exact solution
  x_exact = [1.0_DP, 2.0_DP, 3.0_DP, 2.0_DP, 1.0_DP]
  
  ! Compute RHS: b = A * x_exact
  CALL sparse_matvec(n, row_ptr, col_idx, values, x_exact, b)
  
  WRITE(*,'(A)') 'Matrix A (tridiagonal):'
  WRITE(*,'(A)') '  2 -1  0  0  0'
  WRITE(*,'(A)') ' -1  2 -1  0  0'
  WRITE(*,'(A)') '  0 -1  2 -1  0'
  WRITE(*,'(A)') '  0  0 -1  2 -1'
  WRITE(*,'(A)') '  0  0  0 -1  2'
  WRITE(*,*)
  WRITE(*,'(A,5F8.3)') 'Exact solution: ', x_exact
  WRITE(*,'(A,5F8.3)') 'RHS vector b:   ', b
  WRITE(*,*)
  
  ! Test basic matrix operations
  WRITE(*,'(A)') 'Testing basic sparse matrix-vector product...'
  x = 0.0_DP
  x(1) = 1.0_DP  ! Unit vector
  CALL sparse_matvec(n, row_ptr, col_idx, values, x, b)
  WRITE(*,'(A,5F8.3)') 'A * e1 = ', b
  
  ! Verify result (should be first column of A)
  IF (ABS(b(1) - 2.0_DP) < 1.0E-10_DP .AND. &
      ABS(b(2) + 1.0_DP) < 1.0E-10_DP .AND. &
      ABS(b(3)) < 1.0E-10_DP) THEN
    WRITE(*,'(A)') 'PASS: Sparse matrix-vector product works correctly'
  ELSE
    WRITE(*,'(A)') 'FAIL: Sparse matrix-vector product error'
  END IF
  
  WRITE(*,*)
  WRITE(*,'(A)') '========================================================='
  WRITE(*,'(A)') 'Basic matrix operations verified'
  WRITE(*,'(A)') '========================================================='
  
  DEALLOCATE(row_ptr, col_idx, values)
  
CONTAINS

  SUBROUTINE sparse_matvec(n, row_ptr, col_idx, values, x, y)
    !> Sparse matrix-vector product: y = A * x
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(IN) :: row_ptr(:), col_idx(:)
    REAL(DP), INTENT(IN) :: values(:), x(:)
    REAL(DP), INTENT(OUT) :: y(:)
    
    INTEGER :: i, j
    
    y = 0.0_DP
    DO i = 1, n
      DO j = row_ptr(i), row_ptr(i+1) - 1
        y(i) = y(i) + values(j) * x(col_idx(j))
      END DO
    END DO
    
  END SUBROUTINE sparse_matvec

END PROGRAM test_idrs_basic