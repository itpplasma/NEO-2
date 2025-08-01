MODULE sparse_types_mod
<<<<<<< HEAD
  ! Module containing basic type definitions and parameters
=======
  ! Module defining types for sparse matrix representations
>>>>>>> a1e4b81 (refactor: Extract sparse types and conversions into separate modules)
  ! Extracted from sparse_mod.f90 for better modularity
  
  IMPLICIT NONE
  
  ! Kind parameters
  INTEGER, PARAMETER :: dp = KIND(1.0d0)
  INTEGER, PARAMETER :: long = 8
  
<<<<<<< HEAD
=======
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
  
  ! Type for real sparse matrix in COO format
  TYPE :: sparse_matrix_coo_real
    INTEGER :: nrow = 0                         ! Number of rows
    INTEGER :: ncol = 0                         ! Number of columns
    INTEGER :: nz = 0                           ! Number of nonzeros
    INTEGER, ALLOCATABLE :: irow(:)             ! Row indices (size: nz)
    INTEGER, ALLOCATABLE :: icol(:)             ! Column indices (size: nz)
    REAL(kind=dp), ALLOCATABLE :: val(:)        ! Values (size: nz)
  END TYPE sparse_matrix_coo_real
  
  ! Type for complex sparse matrix in COO format
  TYPE :: sparse_matrix_coo_complex
    INTEGER :: nrow = 0                         ! Number of rows
    INTEGER :: ncol = 0                         ! Number of columns
    INTEGER :: nz = 0                           ! Number of nonzeros
    INTEGER, ALLOCATABLE :: irow(:)             ! Row indices (size: nz)
    INTEGER, ALLOCATABLE :: icol(:)             ! Column indices (size: nz)
    COMPLEX(kind=dp), ALLOCATABLE :: val(:)     ! Values (size: nz)
  END TYPE sparse_matrix_coo_complex
  
  ! Type for real sparse matrix in CSR format (for future use)
  TYPE :: sparse_matrix_csr_real
    INTEGER :: nrow = 0                         ! Number of rows
    INTEGER :: ncol = 0                         ! Number of columns
    INTEGER :: nz = 0                           ! Number of nonzeros
    INTEGER, ALLOCATABLE :: prow(:)             ! Row pointers (size: nrow+1)
    INTEGER, ALLOCATABLE :: icol(:)             ! Column indices (size: nz)
    REAL(kind=dp), ALLOCATABLE :: val(:)        ! Values (size: nz)
  END TYPE sparse_matrix_csr_real
  
  ! Type for complex sparse matrix in CSR format (for future use)
  TYPE :: sparse_matrix_csr_complex
    INTEGER :: nrow = 0                         ! Number of rows
    INTEGER :: ncol = 0                         ! Number of columns
    INTEGER :: nz = 0                           ! Number of nonzeros
    INTEGER, ALLOCATABLE :: prow(:)             ! Row pointers (size: nrow+1)
    INTEGER, ALLOCATABLE :: icol(:)             ! Column indices (size: nz)
    COMPLEX(kind=dp), ALLOCATABLE :: val(:)     ! Values (size: nz)
  END TYPE sparse_matrix_csr_complex
  
  ! Generic interfaces for working with different sparse matrix types
  INTERFACE sparse_get_dimensions
    MODULE PROCEDURE get_dimensions_csc_real
    MODULE PROCEDURE get_dimensions_csc_complex
    MODULE PROCEDURE get_dimensions_coo_real
    MODULE PROCEDURE get_dimensions_coo_complex
    MODULE PROCEDURE get_dimensions_csr_real
    MODULE PROCEDURE get_dimensions_csr_complex
  END INTERFACE sparse_get_dimensions
  
  INTERFACE sparse_deallocate
    MODULE PROCEDURE deallocate_csc_real
    MODULE PROCEDURE deallocate_csc_complex
    MODULE PROCEDURE deallocate_coo_real
    MODULE PROCEDURE deallocate_coo_complex
    MODULE PROCEDURE deallocate_csr_real
    MODULE PROCEDURE deallocate_csr_complex
  END INTERFACE sparse_deallocate
  
CONTAINS
  
  ! Get dimensions for CSC real matrix
  SUBROUTINE get_dimensions_csc_real(mat, nrow, ncol, nz)
    TYPE(sparse_matrix_csc_real), INTENT(in) :: mat
    INTEGER, INTENT(out) :: nrow, ncol, nz
    nrow = mat%nrow
    ncol = mat%ncol
    nz = mat%nz
  END SUBROUTINE get_dimensions_csc_real
  
  ! Get dimensions for CSC complex matrix
  SUBROUTINE get_dimensions_csc_complex(mat, nrow, ncol, nz)
    TYPE(sparse_matrix_csc_complex), INTENT(in) :: mat
    INTEGER, INTENT(out) :: nrow, ncol, nz
    nrow = mat%nrow
    ncol = mat%ncol
    nz = mat%nz
  END SUBROUTINE get_dimensions_csc_complex
  
  ! Get dimensions for COO real matrix
  SUBROUTINE get_dimensions_coo_real(mat, nrow, ncol, nz)
    TYPE(sparse_matrix_coo_real), INTENT(in) :: mat
    INTEGER, INTENT(out) :: nrow, ncol, nz
    nrow = mat%nrow
    ncol = mat%ncol
    nz = mat%nz
  END SUBROUTINE get_dimensions_coo_real
  
  ! Get dimensions for COO complex matrix
  SUBROUTINE get_dimensions_coo_complex(mat, nrow, ncol, nz)
    TYPE(sparse_matrix_coo_complex), INTENT(in) :: mat
    INTEGER, INTENT(out) :: nrow, ncol, nz
    nrow = mat%nrow
    ncol = mat%ncol
    nz = mat%nz
  END SUBROUTINE get_dimensions_coo_complex
  
  ! Get dimensions for CSR real matrix
  SUBROUTINE get_dimensions_csr_real(mat, nrow, ncol, nz)
    TYPE(sparse_matrix_csr_real), INTENT(in) :: mat
    INTEGER, INTENT(out) :: nrow, ncol, nz
    nrow = mat%nrow
    ncol = mat%ncol
    nz = mat%nz
  END SUBROUTINE get_dimensions_csr_real
  
  ! Get dimensions for CSR complex matrix
  SUBROUTINE get_dimensions_csr_complex(mat, nrow, ncol, nz)
    TYPE(sparse_matrix_csr_complex), INTENT(in) :: mat
    INTEGER, INTENT(out) :: nrow, ncol, nz
    nrow = mat%nrow
    ncol = mat%ncol
    nz = mat%nz
  END SUBROUTINE get_dimensions_csr_complex
  
  ! Deallocate CSC real matrix
  SUBROUTINE deallocate_csc_real(mat)
    TYPE(sparse_matrix_csc_real), INTENT(inout) :: mat
    IF (ALLOCATED(mat%irow)) DEALLOCATE(mat%irow)
    IF (ALLOCATED(mat%pcol)) DEALLOCATE(mat%pcol)
    IF (ALLOCATED(mat%val)) DEALLOCATE(mat%val)
    mat%nrow = 0
    mat%ncol = 0
    mat%nz = 0
  END SUBROUTINE deallocate_csc_real
  
  ! Deallocate CSC complex matrix
  SUBROUTINE deallocate_csc_complex(mat)
    TYPE(sparse_matrix_csc_complex), INTENT(inout) :: mat
    IF (ALLOCATED(mat%irow)) DEALLOCATE(mat%irow)
    IF (ALLOCATED(mat%pcol)) DEALLOCATE(mat%pcol)
    IF (ALLOCATED(mat%val)) DEALLOCATE(mat%val)
    mat%nrow = 0
    mat%ncol = 0
    mat%nz = 0
  END SUBROUTINE deallocate_csc_complex
  
  ! Deallocate COO real matrix
  SUBROUTINE deallocate_coo_real(mat)
    TYPE(sparse_matrix_coo_real), INTENT(inout) :: mat
    IF (ALLOCATED(mat%irow)) DEALLOCATE(mat%irow)
    IF (ALLOCATED(mat%icol)) DEALLOCATE(mat%icol)
    IF (ALLOCATED(mat%val)) DEALLOCATE(mat%val)
    mat%nrow = 0
    mat%ncol = 0
    mat%nz = 0
  END SUBROUTINE deallocate_coo_real
  
  ! Deallocate COO complex matrix
  SUBROUTINE deallocate_coo_complex(mat)
    TYPE(sparse_matrix_coo_complex), INTENT(inout) :: mat
    IF (ALLOCATED(mat%irow)) DEALLOCATE(mat%irow)
    IF (ALLOCATED(mat%icol)) DEALLOCATE(mat%icol)
    IF (ALLOCATED(mat%val)) DEALLOCATE(mat%val)
    mat%nrow = 0
    mat%ncol = 0
    mat%nz = 0
  END SUBROUTINE deallocate_coo_complex
  
  ! Deallocate CSR real matrix
  SUBROUTINE deallocate_csr_real(mat)
    TYPE(sparse_matrix_csr_real), INTENT(inout) :: mat
    IF (ALLOCATED(mat%prow)) DEALLOCATE(mat%prow)
    IF (ALLOCATED(mat%icol)) DEALLOCATE(mat%icol)
    IF (ALLOCATED(mat%val)) DEALLOCATE(mat%val)
    mat%nrow = 0
    mat%ncol = 0
    mat%nz = 0
  END SUBROUTINE deallocate_csr_real
  
  ! Deallocate CSR complex matrix
  SUBROUTINE deallocate_csr_complex(mat)
    TYPE(sparse_matrix_csr_complex), INTENT(inout) :: mat
    IF (ALLOCATED(mat%prow)) DEALLOCATE(mat%prow)
    IF (ALLOCATED(mat%icol)) DEALLOCATE(mat%icol)
    IF (ALLOCATED(mat%val)) DEALLOCATE(mat%val)
    mat%nrow = 0
    mat%ncol = 0
    mat%nz = 0
  END SUBROUTINE deallocate_csr_complex
  
>>>>>>> a1e4b81 (refactor: Extract sparse types and conversions into separate modules)
END MODULE sparse_types_mod