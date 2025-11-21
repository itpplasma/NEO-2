MODULE sparse_conversion_mod
  ! Module containing sparse matrix format conversion routines
  ! Extracted from sparse_mod.f90 for better modularity
  
  USE sparse_types_mod, ONLY: dp
  IMPLICIT NONE
  
  PUBLIC :: column_pointer2full
  PUBLIC :: column_full2pointer
  PUBLIC :: sparse2full
  PUBLIC :: full2sparse
  
  INTERFACE column_pointer2full
    MODULE PROCEDURE col_pointer2full
  END INTERFACE column_pointer2full
  
  INTERFACE column_full2pointer
    MODULE PROCEDURE col_full2pointer
  END INTERFACE column_full2pointer
  
  INTERFACE sparse2full
    MODULE PROCEDURE sp2full, sp2fullComplex
  END INTERFACE sparse2full
  
  INTERFACE full2sparse
    MODULE PROCEDURE full2sp, full2spComplex
  END INTERFACE full2sparse
  
CONTAINS
  
  !-------------------------------------------------------------------------------
  ! Convert column pointer pcol to full column index icol
  !
  ! This converts from Compressed Sparse Column (CSC) format to a full
  ! column index array where each element stores its column number
  !
  ! Input:
  !   pcol - Column pointer array (size: ncol+1)
  !
  ! Output:
  !   icol - Full column index array (size: nz)
  !
  SUBROUTINE col_pointer2full(pcol, icol)
    INTEGER, DIMENSION(:), INTENT(in) :: pcol
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(inout) :: icol
    
    INTEGER :: nz
    INTEGER :: nc_old, c, nc, ncol
    
    ncol = SIZE(pcol,1) - 1
    nz = pcol(ncol+1) - 1
    IF (ALLOCATED(icol)) DEALLOCATE(icol)
    ALLOCATE(icol(nz))
    nc_old = 0
    DO c = 1, ncol
      nc = pcol(c+1) - pcol(c)
      icol(nc_old+1:nc_old+nc) = c
      nc_old = nc_old + nc
    END DO
    
  END SUBROUTINE col_pointer2full
  !-------------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------------
  ! Convert full column index icol to column pointer pcol
  !
  ! This converts from a full column index array to Compressed Sparse Column
  ! (CSC) format pointer array
  !
  ! Input:
  !   icol - Full column index array (size: nz)
  !
  ! Output:
  !   pcol - Column pointer array (size: ncol+1)
  !
  SUBROUTINE col_full2pointer(icol, pcol)
    INTEGER, DIMENSION(:), INTENT(in) :: icol
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(inout) :: pcol
    
    INTEGER :: ncol, nz
    INTEGER :: c_c, c_old, k, c, kc
    
    ncol = MAXVAL(icol)
    nz = SIZE(icol,1)
    
    IF (ALLOCATED(pcol)) DEALLOCATE(pcol)
    ALLOCATE(pcol(ncol+1))
    
    c_c = 1
    pcol(c_c) = 1
    c_old = 0
    DO k = 1, nz
      c = icol(k)
      IF (c .NE. c_old) THEN
        IF (c .GT. c_old + 1) THEN
          DO kc = c_old+1, c
            c_c = c_c + 1
            pcol(c_c) = k
          END DO
        ELSE
          c_c = c_c + 1
          pcol(c_c) = k+1
        END IF
        c_old = c
      ELSE
        pcol(c_c) = k+1
      END IF
    END DO
    IF (c_c .LT. ncol+1) pcol(c_c+1:ncol+1) = pcol(c_c)
    
  END SUBROUTINE col_full2pointer
  !-------------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------------
  ! Convert sparse matrix to full (dense) matrix - Real version
  !
  ! Input:
  !   irow - Row indices (size: nz)
  !   pcol - Column pointers (size: ncol+1) or column indices (size: nz)
  !   val  - Matrix values (size: nz)
  !   nrow - Number of rows
  !   ncol - Number of columns
  !
  ! Output:
  !   A - Full matrix (nrow x ncol)
  !
  SUBROUTINE sp2full(irow, pcol, val, nrow, ncol, A)
    INTEGER, DIMENSION(:), INTENT(in) :: irow, pcol
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: val
    INTEGER, INTENT(in) :: nrow, ncol
    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: A
    
    INTEGER :: nz, n, ir, ic
    INTEGER, DIMENSION(:), ALLOCATABLE :: icol
    
    nz = SIZE(val,1)
    IF (SIZE(pcol,1) .NE. nz) THEN
      IF (SIZE(pcol,1) .EQ. ncol+1) THEN
        CALL column_pointer2full(pcol, icol)
      ELSE
        PRINT *, 'Error in sparse2full: icol is not correct'
        STOP
      END IF
    ELSE
      ALLOCATE(icol(SIZE(pcol)))
      icol = pcol
    END IF
    
    IF (ALLOCATED(A)) DEALLOCATE(A)
    ALLOCATE(A(nrow,ncol))
    A = 0.0_dp
    DO n = 1, nz
      ir = irow(n)
      ic = icol(n)
      A(ir,ic) = A(ir,ic) + val(n)
    END DO
    IF (ALLOCATED(icol)) DEALLOCATE(icol)
    
  END SUBROUTINE sp2full
  !-------------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------------
  ! Convert sparse matrix to full (dense) matrix - Complex version
  !
  ! Input:
  !   irow - Row indices (size: nz)
  !   pcol - Column pointers (size: ncol+1) or column indices (size: nz)
  !   val  - Matrix values (size: nz)
  !   nrow - Number of rows
  !   ncol - Number of columns
  !
  ! Output:
  !   A - Full matrix (nrow x ncol)
  !
  SUBROUTINE sp2fullComplex(irow, pcol, val, nrow, ncol, A)
    INTEGER, DIMENSION(:), INTENT(in) :: irow, pcol
    COMPLEX(kind=dp), DIMENSION(:), INTENT(in) :: val
    INTEGER, INTENT(in) :: nrow, ncol
    COMPLEX(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: A
    
    INTEGER :: nz, n, ir, ic
    INTEGER, DIMENSION(:), ALLOCATABLE :: icol
    
    nz = SIZE(val,1)
    IF (SIZE(pcol,1) .NE. nz) THEN
      IF (SIZE(pcol,1) .EQ. ncol+1) THEN
        CALL column_pointer2full(pcol, icol)
      ELSE
        PRINT *, 'Error in sparse2full: icol is not correct'
        STOP
      END IF
    ELSE
      ALLOCATE(icol(SIZE(pcol)))
      icol = pcol
    END IF
    
    IF (ALLOCATED(A)) DEALLOCATE(A)
    ALLOCATE(A(nrow,ncol))
    A = (0.0_dp, 0.0_dp)
    DO n = 1, nz
      ir = irow(n)
      ic = icol(n)
      A(ir,ic) = A(ir,ic) + val(n)
    END DO
    IF (ALLOCATED(icol)) DEALLOCATE(icol)
    
  END SUBROUTINE sp2fullComplex
  !-------------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------------
  ! Convert full (dense) matrix to sparse format - Real version
  !
  ! Input:
  !   A - Full matrix to convert
  !
  ! Output:
  !   irow   - Row indices of nonzero elements
  !   pcol   - Column pointers (CSC format)
  !   values - Nonzero values
  !   nrow   - Number of rows
  !   ncol   - Number of columns
  !   nz_out - Number of nonzeros (optional)
  !
  SUBROUTINE full2sp(A, irow, pcol, values, nrow, ncol, nz_out)
    REAL(kind=dp), DIMENSION(:,:), INTENT(in) :: A
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(inout) :: irow, pcol
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(inout) :: values
    INTEGER, INTENT(out) :: nrow, ncol
    INTEGER, OPTIONAL, INTENT(out) :: nz_out
    
    INTEGER, DIMENSION(:), ALLOCATABLE :: icol
    INTEGER :: nz, nc, nr, n
    
    nrow = SIZE(A,1)
    ncol = SIZE(A,2)
    
    ! Count nonzeros
    nz = 0
    DO nc = 1, ncol
      DO nr = 1, nrow
        IF (A(nr,nc) .NE. 0.0_dp) nz = nz + 1
      END DO
    END DO
    
    ! Allocate arrays
    IF (ALLOCATED(irow)) DEALLOCATE(irow)
    ALLOCATE(irow(nz))
    IF (ALLOCATED(icol)) DEALLOCATE(icol)
    ALLOCATE(icol(nz))
    IF (ALLOCATED(values)) DEALLOCATE(values)
    ALLOCATE(values(nz))
    
    ! Fill arrays
    n = 0
    DO nc = 1, ncol
      DO nr = 1, nrow
        IF (A(nr,nc) .NE. 0.0_dp) THEN
          n = n + 1
          irow(n) = nr
          icol(n) = nc
          values(n) = A(nr,nc)
        END IF
      END DO
    END DO
    
    ! Convert to column pointer format
    CALL column_full2pointer(icol, pcol)
    
    IF (PRESENT(nz_out)) nz_out = nz
    IF (ALLOCATED(icol)) DEALLOCATE(icol)
    
  END SUBROUTINE full2sp
  !-------------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------------
  ! Convert full (dense) matrix to sparse format - Complex version
  !
  ! Input:
  !   A - Full matrix to convert
  !
  ! Output:
  !   irow   - Row indices of nonzero elements
  !   pcol   - Column pointers (CSC format)
  !   val    - Nonzero values
  !   nrow   - Number of rows
  !   ncol   - Number of columns
  !   nz_out - Number of nonzeros (optional)
  !
  SUBROUTINE full2spComplex(A, irow, pcol, val, nrow, ncol, nz_out)
    COMPLEX(kind=dp), DIMENSION(:,:), INTENT(in) :: A
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(inout) :: irow, pcol
    COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(inout) :: val
    INTEGER, INTENT(out) :: nrow, ncol
    INTEGER, OPTIONAL, INTENT(out) :: nz_out
    
    INTEGER, DIMENSION(:), ALLOCATABLE :: icol
    INTEGER :: nz, nc, nr, n
    
    nrow = SIZE(A,1)
    ncol = SIZE(A,2)
    
    ! Count nonzeros
    nz = 0
    DO nc = 1, ncol
      DO nr = 1, nrow
        IF (A(nr,nc) .NE. (0.0_dp, 0.0_dp)) nz = nz + 1
      END DO
    END DO
    
    ! Allocate arrays
    IF (ALLOCATED(irow)) DEALLOCATE(irow)
    ALLOCATE(irow(nz))
    IF (ALLOCATED(icol)) DEALLOCATE(icol)
    ALLOCATE(icol(nz))
    IF (ALLOCATED(val)) DEALLOCATE(val)
    ALLOCATE(val(nz))
    
    ! Fill arrays
    n = 0
    DO nc = 1, ncol
      DO nr = 1, nrow
        IF (A(nr,nc) .NE. (0.0_dp, 0.0_dp)) THEN
          n = n + 1
          irow(n) = nr
          icol(n) = nc
          val(n) = A(nr,nc)
        END IF
      END DO
    END DO
    
    ! Convert to column pointer format
    CALL column_full2pointer(icol, pcol)
    
    IF (PRESENT(nz_out)) nz_out = nz
    IF (ALLOCATED(icol)) DEALLOCATE(icol)
    
  END SUBROUTINE full2spComplex
  !-------------------------------------------------------------------------------
  
END MODULE sparse_conversion_mod