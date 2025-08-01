MODULE sparse_io_mod
  ! Module containing sparse matrix I/O operations
  ! Extracted from sparse_mod.f90 for better modularity
  
  USE sparse_types_mod, ONLY: dp
  IMPLICIT NONE
  
  PUBLIC :: load_mini_example
  PUBLIC :: load_compressed_example
  PUBLIC :: load_standard_example
  PUBLIC :: load_octave_matrices
  PUBLIC :: find_unit
  
  PRIVATE :: load_mini_ex
  PRIVATE :: load_compressed_ex
  PRIVATE :: load_standard_ex
  PRIVATE :: load_octave_mat
  PRIVATE :: load_octave_matComplex
  
  INTERFACE load_mini_example
    MODULE PROCEDURE load_mini_ex
  END INTERFACE load_mini_example
  
  INTERFACE load_compressed_example
    MODULE PROCEDURE load_compressed_ex
  END INTERFACE load_compressed_example
  
  INTERFACE load_standard_example
    MODULE PROCEDURE load_standard_ex
  END INTERFACE load_standard_example
  
  INTERFACE load_octave_matrices
    MODULE PROCEDURE load_octave_mat, load_octave_matComplex
  END INTERFACE load_octave_matrices
  
CONTAINS
  
  !-------------------------------------------------------------------------------
  ! Find a free unit number for file I/O
  !
  ! Input/Output:
  !   unit - Starting unit number, returns first free unit >= input value
  !
  SUBROUTINE find_unit(unit)
    INTEGER, INTENT(inout) :: unit
    LOGICAL :: opened
    
    DO
      INQUIRE(unit=unit, opened=opened)
      IF (.NOT. opened) EXIT
      unit = unit + 1
    END DO
    
  END SUBROUTINE find_unit
  !-------------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------------
  ! Load a mini example matrix (5x5 test matrix)
  !
  ! Output:
  !   A - Full matrix (5x5) with specific test pattern
  !
  SUBROUTINE load_mini_ex(A)
    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(out) :: A
    
    ALLOCATE(A(5,5))
    A(:,1) = (/1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp/)
    A(:,2) = A(:,1)*5 + 2
    A(:,3) = A(:,2)*7 + 2
    A(:,4) = A(:,3)*2 + 2
    A(:,5) = A(:,4)*9 + 2
    
    A(2,4) = 0.0_dp
    A(3,3) = 0.0_dp
    A(4,2) = 0.0_dp
    
  END SUBROUTINE load_mini_ex
  !-------------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------------
  ! Load sparse matrix from compressed format file
  !
  ! File format:
  !   Line 1: nrow ncol nz
  !   Line 2: irow(1:nz)
  !   Line 3: pcol(1:ncol+1)
  !   Line 4: val(1:nz)
  !
  ! Input:
  !   name - Filename to read
  !
  ! Output:
  !   nrow - Number of rows
  !   ncol - Number of columns
  !   nz   - Number of nonzeros
  !   irow - Row indices (size: nz)
  !   pcol - Column pointers (size: ncol+1)
  !   val  - Values (size: nz)
  !
  SUBROUTINE load_compressed_ex(name, nrow, ncol, nz, irow, pcol, val)
    CHARACTER(LEN=*), INTENT(in) :: name
    INTEGER, INTENT(out) :: nrow, ncol, nz
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(out) :: irow, pcol
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(out) :: val
    
    INTEGER :: unit, i
    
    unit = 10
    CALL find_unit(unit)
    OPEN(unit=unit, file=TRIM(ADJUSTL(name)), action='read')
    
    READ(unit,*) nrow, ncol, nz
    ALLOCATE(irow(nz), pcol(ncol+1), val(nz))
    READ(unit,*) (irow(i), i = 1, nz)
    READ(unit,*) (pcol(i), i = 1, ncol+1)
    READ(unit,*) (val(i), i = 1, nz)
    
    CLOSE(unit=unit)
    
  END SUBROUTINE load_compressed_ex
  !-------------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------------
  ! Load sparse matrix from SuperLU Harwell-Boeing format
  !
  ! Input:
  !   name - Filename to read
  !
  ! Output:
  !   nrow - Number of rows
  !   ncol - Number of columns
  !   nz   - Number of nonzeros
  !   irow - Row indices (size: nz)
  !   pcol - Column pointers (size: ncol+1)
  !   val  - Values (size: nz)
  !
  SUBROUTINE load_standard_ex(name, nrow, ncol, nz, irow, pcol, val)
    CHARACTER(LEN=*), INTENT(in) :: name
    INTEGER, INTENT(out) :: nrow, ncol, nz
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(out) :: irow, pcol
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(out) :: val
    
    INTEGER :: unit, i
    
    CHARACTER(len=72) :: fmt1
    CHARACTER(len=72) :: title
    CHARACTER(len=8)  :: key
    CHARACTER(len=3)  :: mxtype
    CHARACTER(len=16) :: ptrfmt, indfmt
    CHARACTER(len=20) :: valfmt, rhsfmt
    
    INTEGER :: totcrd, ptrcrd, indcrd, valcrd, rhscrd, neltvl
    
    fmt1 = '( A72, A8 / 5I14 / A3, 11X, 4I14 / 2A16, 2A20 )'
    
    unit = 10
    CALL find_unit(unit)
    OPEN(unit=unit, file=TRIM(ADJUSTL(name)), action='read')
    
    READ(unit=unit, fmt=fmt1) &
         title, key, totcrd, ptrcrd, indcrd, valcrd, rhscrd, &
         mxtype, nrow, ncol, nz, neltvl, &
         ptrfmt, indfmt, valfmt, rhsfmt
         
    ALLOCATE(irow(nz), pcol(ncol+1), val(nz))
    READ(unit=unit, fmt=ptrfmt) (pcol(i), i = 1, ncol+1)
    READ(unit=unit, fmt=indfmt) (irow(i), i = 1, nz)
    READ(unit=unit, fmt=valfmt) (val(i), i = 1, nz)
    
    CLOSE(unit=unit)
    
  END SUBROUTINE load_standard_ex
  !-------------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------------
  ! Load real sparse matrix from Octave format
  !
  ! Octave format stores matrix in COO format:
  !   Line 1: nrow ncol nz
  !   Lines 2-nz+1: row col value
  !
  ! This routine converts to CSC format
  !
  ! Input:
  !   name - Filename to read
  !
  ! Output:
  !   nrow - Number of rows
  !   ncol - Number of columns
  !   nz   - Number of nonzeros
  !   irow - Row indices in CSC format (size: nz)
  !   pcol - Column pointers (size: ncol+1)
  !   val  - Values (size: nz)
  !
  SUBROUTINE load_octave_mat(name, nrow, ncol, nz, irow, pcol, val)
    CHARACTER(LEN=*), INTENT(in) :: name
    INTEGER, INTENT(out) :: nrow, ncol, nz
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(out) :: irow, pcol
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(out) :: val
    
    INTEGER :: unit, i, k
    INTEGER, DIMENSION(:), ALLOCATABLE :: octave_pcol
    
    ! Open the input file
    unit = 10
    CALL find_unit(unit)
    OPEN(unit=unit, file=TRIM(ADJUSTL(name)), action='read')
    
    ! Read dimensions
    READ(unit,*) nrow, ncol, nz
    ALLOCATE(irow(nz), pcol(ncol+1), octave_pcol(nz), val(nz))
    
    ! Read the sparse matrix in Octave format (COO)
    DO i = 1, nz
      READ(unit,*) irow(i), octave_pcol(i), val(i)
    END DO
    CLOSE(unit=unit)
    
    ! Convert from COO to CSC format
    ! First step: count entries in each column
    pcol = 0
    pcol(1) = 1
    DO i = 1, nz
      pcol(octave_pcol(i)+1) = pcol(octave_pcol(i)+1) + 1
    END DO
    
    ! Second step: cumulative sum to get column pointers
    DO i = 1, ncol
      pcol(i+1) = pcol(i) + pcol(i+1)
    END DO
    
  END SUBROUTINE load_octave_mat
  !-------------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------------
  ! Load complex sparse matrix from Octave format
  !
  ! Octave format stores matrix in COO format:
  !   Line 1: nrow ncol nz
  !   Lines 2-nz+1: row col value
  !
  ! This routine converts to CSC format
  !
  ! Input:
  !   name - Filename to read
  !
  ! Output:
  !   nrow - Number of rows
  !   ncol - Number of columns
  !   nz   - Number of nonzeros
  !   irow - Row indices in CSC format (size: nz)
  !   pcol - Column pointers (size: ncol+1)
  !   val  - Complex values (size: nz)
  !
  SUBROUTINE load_octave_matComplex(name, nrow, ncol, nz, irow, pcol, val)
    CHARACTER(LEN=*), INTENT(in) :: name
    INTEGER, INTENT(out) :: nrow, ncol, nz
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(out) :: irow, pcol
    COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(out) :: val
    
    INTEGER :: unit, i, k
    INTEGER, DIMENSION(:), ALLOCATABLE :: octave_pcol
    
    ! Open the input file
    unit = 10
    CALL find_unit(unit)
    OPEN(unit=unit, file=TRIM(ADJUSTL(name)), action='read')
    
    ! Read dimensions
    READ(unit,*) nrow, ncol, nz
    ALLOCATE(irow(nz), pcol(ncol+1), octave_pcol(nz), val(nz))
    
    ! Read the sparse matrix in Octave format (COO)
    DO i = 1, nz
      READ(unit,*) irow(i), octave_pcol(i), val(i)
    END DO
    CLOSE(unit=unit)
    
    ! Convert from COO to CSC format
    ! First step: count entries in each column
    pcol = 0
    pcol(1) = 1
    DO i = 1, nz
      pcol(octave_pcol(i)+1) = pcol(octave_pcol(i)+1) + 1
    END DO
    
    ! Second step: cumulative sum to get column pointers
    DO i = 1, ncol
      pcol(i+1) = pcol(i) + pcol(i+1)
    END DO
    
  END SUBROUTINE load_octave_matComplex
  !-------------------------------------------------------------------------------
  
END MODULE sparse_io_mod