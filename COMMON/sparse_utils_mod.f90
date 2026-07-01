MODULE sparse_utils_mod
  ! Sparse matrix utilities for CSR/CSC format conversions and operations
  ! Implements CSC<->CSR conversions, matrix-vector multiplication, diagonal extraction
  
  USE sparse_types_mod, ONLY: dp, long
  IMPLICIT NONE
  
  PRIVATE
  
  ! Public interfaces
  PUBLIC :: csc_to_csr
  PUBLIC :: csr_to_csc
  PUBLIC :: csr_matvec
  PUBLIC :: csr_extract_diagonal
  
  ! Generic interfaces for real and complex versions
  INTERFACE csc_to_csr
    MODULE PROCEDURE csc_to_csr_real
    MODULE PROCEDURE csc_to_csr_complex
  END INTERFACE csc_to_csr
  
  INTERFACE csr_to_csc
    MODULE PROCEDURE csr_to_csc_real
    MODULE PROCEDURE csr_to_csc_complex
  END INTERFACE csr_to_csc
  
  INTERFACE csr_matvec
    MODULE PROCEDURE csr_matvec_real
    MODULE PROCEDURE csr_matvec_complex
  END INTERFACE csr_matvec
  
  INTERFACE csr_extract_diagonal
    MODULE PROCEDURE csr_extract_diagonal_real
    MODULE PROCEDURE csr_extract_diagonal_complex
  END INTERFACE csr_extract_diagonal
  
CONTAINS

  !==============================================================================
  ! CSC to CSR conversion - Real version
  !==============================================================================
  SUBROUTINE csc_to_csr_real(nrow, ncol, nnz, csc_col_ptr, csc_row_idx, csc_val, &
                              csr_row_ptr, csr_col_idx, csr_val)
    INTEGER, INTENT(IN) :: nrow, ncol, nnz
    INTEGER, INTENT(IN) :: csc_col_ptr(ncol+1)
    INTEGER, INTENT(IN) :: csc_row_idx(nnz)
    REAL(KIND=dp), INTENT(IN) :: csc_val(nnz)
    INTEGER, INTENT(OUT) :: csr_row_ptr(nrow+1)
    INTEGER, INTENT(OUT) :: csr_col_idx(nnz)
    REAL(KIND=dp), INTENT(OUT) :: csr_val(nnz)
    
    INTEGER :: i, j, k, row, col
    INTEGER :: row_counts(nrow)
    INTEGER :: row_positions(nrow)
    
    ! Handle empty matrix
    IF (nnz == 0) THEN
      csr_row_ptr = 1
      RETURN
    END IF
    
    ! Count non-zeros per row
    row_counts = 0
    DO col = 1, ncol
      DO k = csc_col_ptr(col), csc_col_ptr(col+1)-1
        row = csc_row_idx(k)
        row_counts(row) = row_counts(row) + 1
      END DO
    END DO
    
    ! Set up row pointers
    csr_row_ptr(1) = 1
    DO i = 1, nrow
      csr_row_ptr(i+1) = csr_row_ptr(i) + row_counts(i)
    END DO
    
    ! Initialize current positions for each row
    DO i = 1, nrow
      row_positions(i) = csr_row_ptr(i)
    END DO
    
    ! Fill CSR arrays
    DO col = 1, ncol
      DO k = csc_col_ptr(col), csc_col_ptr(col+1)-1
        row = csc_row_idx(k)
        j = row_positions(row)
        csr_col_idx(j) = col
        csr_val(j) = csc_val(k)
        row_positions(row) = row_positions(row) + 1
      END DO
    END DO
    
  END SUBROUTINE csc_to_csr_real
  
  !==============================================================================
  ! CSC to CSR conversion - Complex version
  !==============================================================================
  SUBROUTINE csc_to_csr_complex(nrow, ncol, nnz, csc_col_ptr, csc_row_idx, csc_val, &
                                 csr_row_ptr, csr_col_idx, csr_val)
    INTEGER, INTENT(IN) :: nrow, ncol, nnz
    INTEGER, INTENT(IN) :: csc_col_ptr(ncol+1)
    INTEGER, INTENT(IN) :: csc_row_idx(nnz)
    COMPLEX(KIND=dp), INTENT(IN) :: csc_val(nnz)
    INTEGER, INTENT(OUT) :: csr_row_ptr(nrow+1)
    INTEGER, INTENT(OUT) :: csr_col_idx(nnz)
    COMPLEX(KIND=dp), INTENT(OUT) :: csr_val(nnz)
    
    INTEGER :: i, j, k, row, col
    INTEGER :: row_counts(nrow)
    INTEGER :: row_positions(nrow)
    
    ! Handle empty matrix
    IF (nnz == 0) THEN
      csr_row_ptr = 1
      RETURN
    END IF
    
    ! Count non-zeros per row
    row_counts = 0
    DO col = 1, ncol
      DO k = csc_col_ptr(col), csc_col_ptr(col+1)-1
        row = csc_row_idx(k)
        row_counts(row) = row_counts(row) + 1
      END DO
    END DO
    
    ! Set up row pointers
    csr_row_ptr(1) = 1
    DO i = 1, nrow
      csr_row_ptr(i+1) = csr_row_ptr(i) + row_counts(i)
    END DO
    
    ! Initialize current positions for each row
    DO i = 1, nrow
      row_positions(i) = csr_row_ptr(i)
    END DO
    
    ! Fill CSR arrays
    DO col = 1, ncol
      DO k = csc_col_ptr(col), csc_col_ptr(col+1)-1
        row = csc_row_idx(k)
        j = row_positions(row)
        csr_col_idx(j) = col
        csr_val(j) = csc_val(k)
        row_positions(row) = row_positions(row) + 1
      END DO
    END DO
    
  END SUBROUTINE csc_to_csr_complex
  
  !==============================================================================
  ! CSR to CSC conversion - Real version
  !==============================================================================
  SUBROUTINE csr_to_csc_real(nrow, ncol, nnz, csr_row_ptr, csr_col_idx, csr_val, &
                              csc_col_ptr, csc_row_idx, csc_val)
    INTEGER, INTENT(IN) :: nrow, ncol, nnz
    INTEGER, INTENT(IN) :: csr_row_ptr(nrow+1)
    INTEGER, INTENT(IN) :: csr_col_idx(nnz)
    REAL(KIND=dp), INTENT(IN) :: csr_val(nnz)
    INTEGER, INTENT(OUT) :: csc_col_ptr(ncol+1)
    INTEGER, INTENT(OUT) :: csc_row_idx(nnz)
    REAL(KIND=dp), INTENT(OUT) :: csc_val(nnz)
    
    INTEGER :: i, j, k, row, col
    INTEGER :: col_counts(ncol)
    INTEGER :: col_positions(ncol)
    
    ! Handle empty matrix
    IF (nnz == 0) THEN
      csc_col_ptr = 1
      RETURN
    END IF
    
    ! Count non-zeros per column
    col_counts = 0
    DO row = 1, nrow
      DO k = csr_row_ptr(row), csr_row_ptr(row+1)-1
        col = csr_col_idx(k)
        col_counts(col) = col_counts(col) + 1
      END DO
    END DO
    
    ! Set up column pointers
    csc_col_ptr(1) = 1
    DO j = 1, ncol
      csc_col_ptr(j+1) = csc_col_ptr(j) + col_counts(j)
    END DO
    
    ! Initialize current positions for each column
    DO j = 1, ncol
      col_positions(j) = csc_col_ptr(j)
    END DO
    
    ! Fill CSC arrays
    DO row = 1, nrow
      DO k = csr_row_ptr(row), csr_row_ptr(row+1)-1
        col = csr_col_idx(k)
        i = col_positions(col)
        csc_row_idx(i) = row
        csc_val(i) = csr_val(k)
        col_positions(col) = col_positions(col) + 1
      END DO
    END DO
    
  END SUBROUTINE csr_to_csc_real
  
  !==============================================================================
  ! CSR to CSC conversion - Complex version
  !==============================================================================
  SUBROUTINE csr_to_csc_complex(nrow, ncol, nnz, csr_row_ptr, csr_col_idx, csr_val, &
                                 csc_col_ptr, csc_row_idx, csc_val)
    INTEGER, INTENT(IN) :: nrow, ncol, nnz
    INTEGER, INTENT(IN) :: csr_row_ptr(nrow+1)
    INTEGER, INTENT(IN) :: csr_col_idx(nnz)
    COMPLEX(KIND=dp), INTENT(IN) :: csr_val(nnz)
    INTEGER, INTENT(OUT) :: csc_col_ptr(ncol+1)
    INTEGER, INTENT(OUT) :: csc_row_idx(nnz)
    COMPLEX(KIND=dp), INTENT(OUT) :: csc_val(nnz)
    
    INTEGER :: i, j, k, row, col
    INTEGER :: col_counts(ncol)
    INTEGER :: col_positions(ncol)
    
    ! Handle empty matrix
    IF (nnz == 0) THEN
      csc_col_ptr = 1
      RETURN
    END IF
    
    ! Count non-zeros per column
    col_counts = 0
    DO row = 1, nrow
      DO k = csr_row_ptr(row), csr_row_ptr(row+1)-1
        col = csr_col_idx(k)
        col_counts(col) = col_counts(col) + 1
      END DO
    END DO
    
    ! Set up column pointers
    csc_col_ptr(1) = 1
    DO j = 1, ncol
      csc_col_ptr(j+1) = csc_col_ptr(j) + col_counts(j)
    END DO
    
    ! Initialize current positions for each column
    DO j = 1, ncol
      col_positions(j) = csc_col_ptr(j)
    END DO
    
    ! Fill CSC arrays
    DO row = 1, nrow
      DO k = csr_row_ptr(row), csr_row_ptr(row+1)-1
        col = csr_col_idx(k)
        i = col_positions(col)
        csc_row_idx(i) = row
        csc_val(i) = csr_val(k)
        col_positions(col) = col_positions(col) + 1
      END DO
    END DO
    
  END SUBROUTINE csr_to_csc_complex
  
  !==============================================================================
  ! CSR matrix-vector multiplication: y = A*x - Real version
  !==============================================================================
  SUBROUTINE csr_matvec_real(nrow, csr_row_ptr, csr_col_idx, csr_val, x, y)
    INTEGER, INTENT(IN) :: nrow
    INTEGER, INTENT(IN) :: csr_row_ptr(nrow+1)
    INTEGER, INTENT(IN) :: csr_col_idx(:)
    REAL(KIND=dp), INTENT(IN) :: csr_val(:)
    REAL(KIND=dp), INTENT(IN) :: x(:)
    REAL(KIND=dp), INTENT(OUT) :: y(nrow)
    
    INTEGER :: i, k
    REAL(KIND=dp) :: sum
    
    ! Compute y = A*x
    DO i = 1, nrow
      sum = 0.0_dp
      DO k = csr_row_ptr(i), csr_row_ptr(i+1)-1
        sum = sum + csr_val(k) * x(csr_col_idx(k))
      END DO
      y(i) = sum
    END DO
    
  END SUBROUTINE csr_matvec_real
  
  !==============================================================================
  ! CSR matrix-vector multiplication: y = A*x - Complex version
  !==============================================================================
  SUBROUTINE csr_matvec_complex(nrow, csr_row_ptr, csr_col_idx, csr_val, x, y)
    INTEGER, INTENT(IN) :: nrow
    INTEGER, INTENT(IN) :: csr_row_ptr(nrow+1)
    INTEGER, INTENT(IN) :: csr_col_idx(:)
    COMPLEX(KIND=dp), INTENT(IN) :: csr_val(:)
    COMPLEX(KIND=dp), INTENT(IN) :: x(:)
    COMPLEX(KIND=dp), INTENT(OUT) :: y(nrow)
    
    INTEGER :: i, k
    COMPLEX(KIND=dp) :: sum
    
    ! Compute y = A*x
    DO i = 1, nrow
      sum = (0.0_dp, 0.0_dp)
      DO k = csr_row_ptr(i), csr_row_ptr(i+1)-1
        sum = sum + csr_val(k) * x(csr_col_idx(k))
      END DO
      y(i) = sum
    END DO
    
  END SUBROUTINE csr_matvec_complex
  
  !==============================================================================
  ! Extract diagonal from CSR matrix - Real version
  !==============================================================================
  SUBROUTINE csr_extract_diagonal_real(nrow, csr_row_ptr, csr_col_idx, csr_val, diag)
    INTEGER, INTENT(IN) :: nrow
    INTEGER, INTENT(IN) :: csr_row_ptr(nrow+1)
    INTEGER, INTENT(IN) :: csr_col_idx(:)
    REAL(KIND=dp), INTENT(IN) :: csr_val(:)
    REAL(KIND=dp), INTENT(OUT) :: diag(nrow)
    
    INTEGER :: i, k
    LOGICAL :: found
    
    ! Extract diagonal elements
    DO i = 1, nrow
      found = .FALSE.
      DO k = csr_row_ptr(i), csr_row_ptr(i+1)-1
        IF (csr_col_idx(k) == i) THEN
          diag(i) = csr_val(k)
          found = .TRUE.
          EXIT
        END IF
      END DO
      IF (.NOT. found) THEN
        diag(i) = 0.0_dp
      END IF
    END DO
    
  END SUBROUTINE csr_extract_diagonal_real
  
  !==============================================================================
  ! Extract diagonal from CSR matrix - Complex version
  !==============================================================================
  SUBROUTINE csr_extract_diagonal_complex(nrow, csr_row_ptr, csr_col_idx, csr_val, diag)
    INTEGER, INTENT(IN) :: nrow
    INTEGER, INTENT(IN) :: csr_row_ptr(nrow+1)
    INTEGER, INTENT(IN) :: csr_col_idx(:)
    COMPLEX(KIND=dp), INTENT(IN) :: csr_val(:)
    COMPLEX(KIND=dp), INTENT(OUT) :: diag(nrow)
    
    INTEGER :: i, k
    LOGICAL :: found
    
    ! Extract diagonal elements
    DO i = 1, nrow
      found = .FALSE.
      DO k = csr_row_ptr(i), csr_row_ptr(i+1)-1
        IF (csr_col_idx(k) == i) THEN
          diag(i) = csr_val(k)
          found = .TRUE.
          EXIT
        END IF
      END DO
      IF (.NOT. found) THEN
        diag(i) = (0.0_dp, 0.0_dp)
      END IF
    END DO
    
  END SUBROUTINE csr_extract_diagonal_complex
  
END MODULE sparse_utils_mod