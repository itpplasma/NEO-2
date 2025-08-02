MODULE sparse_types_mod
  ! Module defining types for sparse matrix representations
  ! Extracted from sparse_mod.f90 for better modularity
  
  IMPLICIT NONE
  
  ! Kind parameters
  INTEGER, PARAMETER :: dp = KIND(1.0d0)
  INTEGER, PARAMETER :: long = 8
  
  ! Sparse matrix storage formats
  INTEGER, PARAMETER :: SPARSE_FORMAT_COO = 1  ! Coordinate format
  INTEGER, PARAMETER :: SPARSE_FORMAT_CSC = 2  ! Compressed Sparse Column
  INTEGER, PARAMETER :: SPARSE_FORMAT_CSR = 3  ! Compressed Sparse Row
  
  ! Type for real sparse matrix in CSC format
  TYPE :: sparse_matrix_csc_real
    INTEGER :: nrow = 0                         ! Number of rows
    INTEGER :: ncol = 0                         ! Number of columns  
    INTEGER :: nz = 0                           ! Number of nonzeros
    INTEGER, ALLOCATABLE :: irow(:)             ! Row indices (size: nz)
    INTEGER, ALLOCATABLE :: pcol(:)             ! Column pointers (size: ncol+1)
    REAL(kind=dp), ALLOCATABLE :: val(:)        ! Values (size: nz)
  END TYPE sparse_matrix_csc_real
  
  ! Type for complex sparse matrix in CSC format
  TYPE :: sparse_matrix_csc_complex
    INTEGER :: nrow = 0                         ! Number of rows
    INTEGER :: ncol = 0                         ! Number of columns
    INTEGER :: nz = 0                           ! Number of nonzeros
    INTEGER, ALLOCATABLE :: irow(:)             ! Row indices (size: nz)
    INTEGER, ALLOCATABLE :: pcol(:)             ! Column pointers (size: ncol+1)
    COMPLEX(kind=dp), ALLOCATABLE :: val(:)     ! Values (size: nz)
  END TYPE sparse_matrix_csc_complex
  
  ! Type for real sparse matrix in CSR format (future use)
  TYPE :: sparse_matrix_csr_real
    INTEGER :: nrow = 0                         ! Number of rows
    INTEGER :: ncol = 0                         ! Number of columns
    INTEGER :: nz = 0                           ! Number of nonzeros  
    INTEGER, ALLOCATABLE :: jcol(:)             ! Column indices (size: nz)
    INTEGER, ALLOCATABLE :: prow(:)             ! Row pointers (size: nrow+1)
    REAL(kind=dp), ALLOCATABLE :: val(:)        ! Values (size: nz)
  END TYPE sparse_matrix_csr_real
  
  ! Type for complex sparse matrix in CSR format (future use)
  TYPE :: sparse_matrix_csr_complex
    INTEGER :: nrow = 0                         ! Number of rows
    INTEGER :: ncol = 0                         ! Number of columns
    INTEGER :: nz = 0                           ! Number of nonzeros
    INTEGER, ALLOCATABLE :: jcol(:)             ! Column indices (size: nz)
    INTEGER, ALLOCATABLE :: prow(:)             ! Row pointers (size: nrow+1)
    COMPLEX(kind=dp), ALLOCATABLE :: val(:)     ! Values (size: nz)
  END TYPE sparse_matrix_csr_complex

  ! Generic interfaces for allocation/deallocation
  INTERFACE allocate_sparse
    MODULE PROCEDURE allocate_csc_real, allocate_csc_complex
    MODULE PROCEDURE allocate_csr_real, allocate_csr_complex
  END INTERFACE allocate_sparse
  
  INTERFACE deallocate_sparse
    MODULE PROCEDURE deallocate_csc_real, deallocate_csc_complex
    MODULE PROCEDURE deallocate_csr_real, deallocate_csr_complex
  END INTERFACE deallocate_sparse

CONTAINS

  ! CSC Real allocation
  SUBROUTINE allocate_csc_real(mat, nrow, ncol, nz)
    TYPE(sparse_matrix_csc_real), INTENT(inout) :: mat
    INTEGER, INTENT(in) :: nrow, ncol, nz
    
    mat%nrow = nrow
    mat%ncol = ncol
    mat%nz = nz
    
    IF (ALLOCATED(mat%irow)) DEALLOCATE(mat%irow)
    IF (ALLOCATED(mat%pcol)) DEALLOCATE(mat%pcol)
    IF (ALLOCATED(mat%val)) DEALLOCATE(mat%val)
    
    ALLOCATE(mat%irow(nz))
    ALLOCATE(mat%pcol(ncol+1))
    ALLOCATE(mat%val(nz))
  END SUBROUTINE allocate_csc_real
  
  ! CSC Complex allocation
  SUBROUTINE allocate_csc_complex(mat, nrow, ncol, nz)
    TYPE(sparse_matrix_csc_complex), INTENT(inout) :: mat
    INTEGER, INTENT(in) :: nrow, ncol, nz
    
    mat%nrow = nrow
    mat%ncol = ncol
    mat%nz = nz
    
    IF (ALLOCATED(mat%irow)) DEALLOCATE(mat%irow)
    IF (ALLOCATED(mat%pcol)) DEALLOCATE(mat%pcol)
    IF (ALLOCATED(mat%val)) DEALLOCATE(mat%val)
    
    ALLOCATE(mat%irow(nz))
    ALLOCATE(mat%pcol(ncol+1))
    ALLOCATE(mat%val(nz))
  END SUBROUTINE allocate_csc_complex
  
  ! CSR Real allocation
  SUBROUTINE allocate_csr_real(mat, nrow, ncol, nz)
    TYPE(sparse_matrix_csr_real), INTENT(inout) :: mat
    INTEGER, INTENT(in) :: nrow, ncol, nz
    
    mat%nrow = nrow
    mat%ncol = ncol
    mat%nz = nz
    
    IF (ALLOCATED(mat%jcol)) DEALLOCATE(mat%jcol)
    IF (ALLOCATED(mat%prow)) DEALLOCATE(mat%prow)
    IF (ALLOCATED(mat%val)) DEALLOCATE(mat%val)
    
    ALLOCATE(mat%jcol(nz))
    ALLOCATE(mat%prow(nrow+1))
    ALLOCATE(mat%val(nz))
  END SUBROUTINE allocate_csr_real
  
  ! CSR Complex allocation
  SUBROUTINE allocate_csr_complex(mat, nrow, ncol, nz)
    TYPE(sparse_matrix_csr_complex), INTENT(inout) :: mat
    INTEGER, INTENT(in) :: nrow, ncol, nz
    
    mat%nrow = nrow
    mat%ncol = ncol
    mat%nz = nz
    
    IF (ALLOCATED(mat%jcol)) DEALLOCATE(mat%jcol)
    IF (ALLOCATED(mat%prow)) DEALLOCATE(mat%prow)
    IF (ALLOCATED(mat%val)) DEALLOCATE(mat%val)
    
    ALLOCATE(mat%jcol(nz))
    ALLOCATE(mat%prow(nrow+1))
    ALLOCATE(mat%val(nz))
  END SUBROUTINE allocate_csr_complex
  
  ! CSC Real deallocation
  SUBROUTINE deallocate_csc_real(mat)
    TYPE(sparse_matrix_csc_real), INTENT(inout) :: mat
    
    IF (ALLOCATED(mat%irow)) DEALLOCATE(mat%irow)
    IF (ALLOCATED(mat%pcol)) DEALLOCATE(mat%pcol)
    IF (ALLOCATED(mat%val)) DEALLOCATE(mat%val)
    
    mat%nrow = 0
    mat%ncol = 0
    mat%nz = 0
  END SUBROUTINE deallocate_csc_real
  
  ! CSC Complex deallocation
  SUBROUTINE deallocate_csc_complex(mat)
    TYPE(sparse_matrix_csc_complex), INTENT(inout) :: mat
    
    IF (ALLOCATED(mat%irow)) DEALLOCATE(mat%irow)
    IF (ALLOCATED(mat%pcol)) DEALLOCATE(mat%pcol)
    IF (ALLOCATED(mat%val)) DEALLOCATE(mat%val)
    
    mat%nrow = 0
    mat%ncol = 0
    mat%nz = 0
  END SUBROUTINE deallocate_csc_complex
  
  ! CSR Real deallocation
  SUBROUTINE deallocate_csr_real(mat)
    TYPE(sparse_matrix_csr_real), INTENT(inout) :: mat
    
    IF (ALLOCATED(mat%jcol)) DEALLOCATE(mat%jcol)
    IF (ALLOCATED(mat%prow)) DEALLOCATE(mat%prow)
    IF (ALLOCATED(mat%val)) DEALLOCATE(mat%val)
    
    mat%nrow = 0
    mat%ncol = 0
    mat%nz = 0
  END SUBROUTINE deallocate_csr_real
  
  ! CSR Complex deallocation
  SUBROUTINE deallocate_csr_complex(mat)
    TYPE(sparse_matrix_csr_complex), INTENT(inout) :: mat
    
    IF (ALLOCATED(mat%jcol)) DEALLOCATE(mat%jcol)
    IF (ALLOCATED(mat%prow)) DEALLOCATE(mat%prow)
    IF (ALLOCATED(mat%val)) DEALLOCATE(mat%val)
    
    mat%nrow = 0
    mat%ncol = 0
    mat%nz = 0
  END SUBROUTINE deallocate_csr_complex
  
END MODULE sparse_types_mod