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
  PUBLIC :: remap_rc
  
  ! Export specific procedures for sparse_solvers_mod
  PUBLIC :: csc_to_csr_real, csc_to_csr_complex
  
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
  
  INTERFACE remap_rc
    MODULE PROCEDURE remap_rc_real
    MODULE PROCEDURE remap_rc_complex
  END INTERFACE remap_rc
  
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
  
  !==============================================================================
  ! Re-arrange and squeeze sparse matrix elements - Real version
  !==============================================================================
  SUBROUTINE remap_rc_real(nz,nz_sqeezed,irow,icol,amat)
    ! Re-arranges matrix elements which may be unordered and may have
    ! different elements with the same row and column indices is such
    ! a way that column index, icol, forms a non-decreasing sequence
    ! and row index, irow, forms increasing sub-sequences for itervals
    ! with a fixed column index. Sums up elements of the matrix which
    ! have the same row and column indices to one element with these
    ! indices
    !
    ! Arguments:
    ! nz          - (input)  number of elements in irow,icol,amat
    ! nz_sqeezed  - (output) number of elements with different (irow(k),icol(k))
    ! irow        - (inout)  row indices
    ! icol        - (inout)  column indices
    ! amat        - (inout)  matrix values

    INTEGER, INTENT(in)                          :: nz
    INTEGER, INTENT(out)                         :: nz_sqeezed
    INTEGER, DIMENSION(nz), INTENT(inout)        :: irow,icol
    REAL(kind=dp), DIMENSION(nz), INTENT(inout)  :: amat

    INTEGER                            :: ncol,i,j,k,kbeg,kend,ips,iflag,ksq
    INTEGER, DIMENSION(:), ALLOCATABLE :: nrows,icount,ipoi
    INTEGER                            :: ksq_ne0
    INTEGER, DIMENSION(:), ALLOCATABLE :: kne0

    ncol=MAXVAL(icol)
    ALLOCATE(nrows(ncol),icount(ncol),ipoi(nz))
    nrows=0

    ! count number of rows in a given column:
    DO k=1,nz
       j=icol(k)
       nrows(j)=nrows(j)+1
    ENDDO

    ! compute starting index - 1 of rows in a general list for each column:

    icount(1)=0

    DO i=1,ncol-1
       icount(i+1)=icount(i)+nrows(i)
    ENDDO

    ! compute the pointer from the list ordered by columns to a general list

    DO k=1,nz
       j=icol(k)
       icount(j)=icount(j)+1
       ipoi(icount(j))=k
    ENDDO

    ! re-order row indices to non-decreasing sub-sequences

    DO i=1,ncol
       kend=icount(i)
       kbeg=kend-nrows(i)+1
       DO j=1,kend-kbeg
          iflag=0
          DO k=kbeg+1,kend
             IF(irow(ipoi(k)).LT.irow(ipoi(k-1))) THEN
                iflag=1
                ips=ipoi(k)
                ipoi(k)=ipoi(k-1)
                ipoi(k-1)=ips
             ENDIF
          ENDDO
          IF(iflag.EQ.0) EXIT
       ENDDO
    ENDDO

    irow=irow(ipoi)
    icol=icol(ipoi)
    amat=amat(ipoi)

    ! squeese the data - sum up matrix elements with the same indices

    ksq=1
    !
    DO k=2,nz
       IF(irow(k).EQ.irow(k-1).AND.icol(k).EQ.icol(k-1)) THEN
          amat(ksq)=amat(ksq)+amat(k)
       ELSE
          ksq=ksq+1
          irow(ksq)=irow(k)
          icol(ksq)=icol(k)
          amat(ksq)=amat(k)
       ENDIF
    ENDDO

    ! remove zeros from the sparse vector

    ALLOCATE(kne0(ksq))
    ksq_ne0=0
    DO k=1,ksq
       IF(amat(k) .NE. 0.0d0) THEN
          ksq_ne0=ksq_ne0+1
          kne0(ksq_ne0)=k
       ENDIF
    ENDDO
    IF(ksq_ne0 .EQ. 0) THEN
       PRINT *,'sparse_utils_mod/remap_rc: All entries of the sparse vector are zero!'
    ELSE
       irow(1:ksq_ne0)=irow(kne0(1:ksq_ne0))
       icol(1:ksq_ne0)=icol(kne0(1:ksq_ne0))
       amat(1:ksq_ne0)=amat(kne0(1:ksq_ne0))
    ENDIF

    nz_sqeezed=ksq_ne0
    DEALLOCATE(nrows,icount,ipoi,kne0)

  END SUBROUTINE remap_rc_real
  
  !==============================================================================
  ! Re-arrange and squeeze sparse matrix elements - Complex version
  !==============================================================================
  SUBROUTINE remap_rc_complex(nz,nz_sqeezed,irow,icol,amat)
    ! Re-arranges matrix elements which may be unordered and may have
    ! different elements with the same row and column indices is such
    ! a way that column index, icol, forms a non-decreasing sequence
    ! and row index, irow, forms increasing sub-sequences for itervals
    ! with a fixed column index. Sums up elements of the matrix which
    ! have the same row and column indices to one element with these
    ! indices
    !
    ! Arguments:
    ! nz          - (input)  number of elements in irow,icol,amat
    ! nz_sqeezed  - (output) number of elements with different (irow(k),icol(k))
    ! irow        - (inout)  row indices
    ! icol        - (inout)  column indices
    ! amat        - (inout)  matrix values

    INTEGER, INTENT(in)                          :: nz
    INTEGER, INTENT(out)                         :: nz_sqeezed
    INTEGER, DIMENSION(nz), INTENT(inout)        :: irow,icol
    COMPLEX(kind=dp), DIMENSION(nz), INTENT(inout) :: amat

    INTEGER                            :: ncol,i,j,k,kbeg,kend,ips,iflag,ksq
    INTEGER, DIMENSION(:), ALLOCATABLE :: nrows,icount,ipoi
    INTEGER                            :: ksq_ne0
    INTEGER, DIMENSION(:), ALLOCATABLE :: kne0

    ncol=MAXVAL(icol)
    ALLOCATE(nrows(ncol),icount(ncol),ipoi(nz))
    nrows=0

    ! count number of rows in a given column:

    DO k=1,nz
       j=icol(k)
       nrows(j)=nrows(j)+1
    ENDDO

    ! compute starting index - 1 of rows in a general list for each column:

    icount(1)=0
    !
    DO i=1,ncol-1
       icount(i+1)=icount(i)+nrows(i)
    ENDDO

    ! compute the pointer from the list ordered by columns to a general list

    DO k=1,nz
       j=icol(k)
       icount(j)=icount(j)+1
       ipoi(icount(j))=k
    ENDDO

    ! re-order row indices to non-decreasing sub-sequences

    DO i=1,ncol
       kend=icount(i)
       kbeg=kend-nrows(i)+1
       DO j=1,kend-kbeg
          iflag=0
          DO k=kbeg+1,kend
             IF(irow(ipoi(k)).LT.irow(ipoi(k-1))) THEN
                iflag=1
                ips=ipoi(k)
                ipoi(k)=ipoi(k-1)
                ipoi(k-1)=ips
             ENDIF
          ENDDO
          IF(iflag.EQ.0) EXIT
       ENDDO
    ENDDO

    irow=irow(ipoi)
    icol=icol(ipoi)
    amat=amat(ipoi)

    ! squeese the data - sum up matrix elements with the same indices

    ksq=1

    DO k=2,nz
       IF(irow(k).EQ.irow(k-1).AND.icol(k).EQ.icol(k-1)) THEN
          amat(ksq)=amat(ksq)+amat(k)
       ELSE
          ksq=ksq+1
          irow(ksq)=irow(k)
          icol(ksq)=icol(k)
          amat(ksq)=amat(k)
       ENDIF
    ENDDO

    ! remove zeros from the sparse vector

    ALLOCATE(kne0(ksq))
    ksq_ne0=0
    DO k=1,ksq
       IF(amat(k) .NE. (0.d0,0.d0)) THEN
          ksq_ne0=ksq_ne0+1
          kne0(ksq_ne0)=k
       ENDIF
    ENDDO
    IF(ksq_ne0 .EQ. 0) THEN
       PRINT *,'sparse_utils_mod/remap_rc: All entries of the sparse vector are zero!'
    ELSE
       irow(1:ksq_ne0)=irow(kne0(1:ksq_ne0))
       icol(1:ksq_ne0)=icol(kne0(1:ksq_ne0))
       amat(1:ksq_ne0)=amat(kne0(1:ksq_ne0))
    ENDIF

    nz_sqeezed=ksq_ne0
    DEALLOCATE(nrows,icount,ipoi,kne0)

  END SUBROUTINE remap_rc_complex
  
END MODULE sparse_utils_mod