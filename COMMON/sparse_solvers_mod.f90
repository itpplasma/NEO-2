MODULE sparse_solvers_mod
  ! Module containing sparse matrix solver operations
  ! Extracted from sparse_mod.f90 for better modularity
  
  USE sparse_types_mod, ONLY: dp, long
  USE sparse_conversion_mod
  IMPLICIT NONE
  
  PUBLIC :: sparse_solve
  PUBLIC :: sparse_solve_suitesparse
  PUBLIC :: sparse_solve_method
  PUBLIC :: factorization_exists
  
  INTEGER :: sparse_solve_method = 3
  LOGICAL :: factorization_exists = .FALSE.
  
  ! SuiteSparse solver data address pointers
  INTEGER(kind=long), PRIVATE :: symbolic, numeric
  INTEGER(kind=long), PRIVATE :: sys = 0
  REAL(kind=dp), PRIVATE :: control(20), info_suitesparse(90)
  
  INTERFACE sparse_solve
    MODULE PROCEDURE sparse_solveReal_b1, sparse_solveReal_b2, sparse_solveReal_A_b1, sparse_solveReal_A_b2, &
                     sparse_solveComplex_b1, sparse_solveComplex_b2, sparse_solveComplex_A_b1, sparse_solveComplex_A_b2
  END INTERFACE sparse_solve
  
  INTERFACE sparse_solve_suitesparse
    MODULE PROCEDURE sparse_solve_suitesparse_b1, sparse_solve_suitesparse_b2_loop, &
                     sparse_solve_suitesparseComplex_b1, sparse_solve_suitesparseComplex_b2_loop
  END INTERFACE sparse_solve_suitesparse
  
  PRIVATE :: sparse_solveReal_b1, sparse_solveReal_b2, sparse_solveReal_A_b1, sparse_solveReal_A_b2
  PRIVATE :: sparse_solveComplex_b1, sparse_solveComplex_b2, sparse_solveComplex_A_b1, sparse_solveComplex_A_b2
  PRIVATE :: sparse_solve_suitesparse_b1, sparse_solve_suitesparse_b2_loop
  PRIVATE :: sparse_solve_suitesparseComplex_b1, sparse_solve_suitesparseComplex_b2_loop
  
CONTAINS
  
  !-------------------------------------------------------------------------------
  ! solves A*x = b for sparse A and 1-D vector b
  ! A is specified through nrow,ncol,nz,irow,pcol,val
  ! results are returned in b
  SUBROUTINE sparse_solveReal_b1(nrow,ncol,nz,irow,pcol,val,b,iopt_in)
    INTEGER, INTENT(in) :: nrow,ncol,nz
    INTEGER, DIMENSION(:), INTENT(in) :: irow,pcol
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: val
    REAL(kind=dp), DIMENSION(:), INTENT(inout) :: b
    INTEGER, INTENT(in), OPTIONAL :: iopt_in
    
    INTEGER :: iopt = 0
    INTEGER :: n
    INTEGER, DIMENSION(:), ALLOCATABLE :: pcoln
    LOGICAL :: pcol_modified
    
    IF (PRESENT(iopt_in)) iopt = iopt_in
    
    pcol_modified = .FALSE.
    ! pcoln to remove "C convention" used in SuiteSparse
    IF (pcol(1) .EQ. 0) THEN
       pcol_modified = .TRUE.
       ALLOCATE(pcoln(ncol+1))
       pcoln = pcol + 1
    END IF
    
    ! check about existing factorization
    IF (factorization_exists .AND. iopt .EQ. 1) THEN ! free memory first
       IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. 3) ) THEN
          IF (pcol_modified) THEN
             CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcoln,val,b,3)
          ELSE
             CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,3)
          END IF
       END IF
    END IF
    IF (.NOT. factorization_exists .AND. iopt .EQ. 2) THEN ! factorize first
       IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. 3) ) THEN
          IF (pcol_modified) THEN
             CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcoln,val,b,1)
          ELSE
             CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,1)
          END IF
       END IF
       factorization_exists = .TRUE.
    END IF
    IF (iopt .EQ. 1) factorization_exists = .TRUE.
    IF (iopt .EQ. 3) factorization_exists = .FALSE.
    
    IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. 3) ) THEN
       IF (pcol_modified) THEN
          CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcoln,val,b,iopt)
       ELSE
          CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,iopt)
       END IF
    ELSE
       PRINT *, 'sparse_solve_method ',sparse_solve_method,'not implemented'
       STOP
    END IF
    
    IF (ALLOCATED(pcoln)) DEALLOCATE(pcoln)
    RETURN
  END SUBROUTINE sparse_solveReal_b1
  !-------------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------------
  ! solves A*x = b for sparse A and 1-D vector b
  ! A is specified through nrow,ncol,nz,irow,pcol,val
  ! results are returned in b
  SUBROUTINE sparse_solveComplex_b1(nrow,ncol,nz,irow,pcol,val,b,iopt_in)
    INTEGER, INTENT(in) :: nrow,ncol,nz
    INTEGER, DIMENSION(:), INTENT(in) :: irow,pcol
    COMPLEX(kind=dp), DIMENSION(:), INTENT(in) :: val
    COMPLEX(kind=dp), DIMENSION(:), INTENT(inout) :: b
    INTEGER, INTENT(in), OPTIONAL :: iopt_in
    
    INTEGER :: iopt = 0
    INTEGER :: n
    INTEGER, DIMENSION(:), ALLOCATABLE :: pcoln
    LOGICAL :: pcol_modified
    
    IF (PRESENT(iopt_in)) iopt = iopt_in
    
    pcol_modified = .FALSE.
    ! pcoln to remove "C convention" used in SuiteSparse
    IF (pcol(1) .EQ. 0) THEN
       pcol_modified = .TRUE.
       ALLOCATE(pcoln(ncol+1))
       pcoln = pcol + 1
    END IF
    
    ! check about existing factorization
    IF (factorization_exists .AND. iopt .EQ. 1) THEN ! free memory first
       IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. 3) ) THEN
          IF (pcol_modified) THEN
             CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcoln,val,b,3)
          ELSE
             CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,3)
          END IF
       END IF
    END IF
    IF (.NOT. factorization_exists .AND. iopt .EQ. 2) THEN ! factorize first
       IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. 3) ) THEN
          IF (pcol_modified) THEN
             CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcoln,val,b,1)
          ELSE
             CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,1)
          END IF
       END IF
       factorization_exists = .TRUE.
    END IF
    IF (iopt .EQ. 1) factorization_exists = .TRUE.
    IF (iopt .EQ. 3) factorization_exists = .FALSE.
    
    IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. 3) ) THEN
       IF (pcol_modified) THEN
          CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcoln,val,b,iopt)
       ELSE
          CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,iopt)
       END IF
    ELSE
       PRINT *, 'sparse_solve_method ',sparse_solve_method,'not implemented'
       STOP
    END IF
    
    IF (ALLOCATED(pcoln)) DEALLOCATE(pcoln)
    
  END SUBROUTINE sparse_solveComplex_b1
  !-------------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------------
  ! solves A*x = b for sparse A and 1-D array b
  ! A is specified through nrow,ncol,nz,irow,pcol,val
  ! results are returned in b
  SUBROUTINE sparse_solveReal_b2(nrow,ncol,nz,irow,pcol,val,b,iopt_in)
    INTEGER, INTENT(in) :: nrow,ncol,nz
    INTEGER, DIMENSION(:), INTENT(in) :: irow,pcol
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: val
    REAL(kind=dp), DIMENSION(:,:), INTENT(inout) :: b
    INTEGER, INTENT(in), OPTIONAL :: iopt_in
    
    INTEGER :: iopt = 0
    INTEGER :: n
    INTEGER, DIMENSION(:), ALLOCATABLE :: pcoln
    LOGICAL :: pcol_modified
    
    IF (PRESENT(iopt_in)) iopt = iopt_in
    
    pcol_modified = .FALSE.
    ! pcoln to remove "C convention" used in SuiteSparse
    IF (pcol(1) .EQ. 0) THEN
       pcol_modified = .TRUE.
       ALLOCATE(pcoln(ncol+1))
       pcoln = pcol + 1
    END IF
    
    ! check about existing factorization
    IF (factorization_exists .AND. iopt .EQ. 1) THEN ! free memory first
       IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. 3) ) THEN
          IF (pcol_modified) THEN
             CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcoln,val,b,3)
          ELSE
             CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,3)
          END IF
       END IF
    END IF
    IF (.NOT. factorization_exists .AND. iopt .EQ. 2) THEN ! factorize first
       IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. 3) ) THEN
          IF (pcol_modified) THEN
             CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcoln,val,b,1)
          ELSE
             CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,1)
          END IF
       END IF
       factorization_exists = .TRUE.
    END IF
    IF (iopt .EQ. 1) factorization_exists = .TRUE.
    IF (iopt .EQ. 3) factorization_exists = .FALSE.
    
    IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. 3) ) THEN
       IF (pcol_modified) THEN
          CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcoln,val,b,iopt)
       ELSE
          CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,iopt)
       END IF
    ELSE
       PRINT *, 'sparse_solve_method ',sparse_solve_method,'not implemented'
       STOP
    END IF
    
    IF (ALLOCATED(pcoln)) DEALLOCATE(pcoln)
    
  END SUBROUTINE sparse_solveReal_b2
  !-------------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------------
  ! solves A*x = b for sparse A and 1-D array b
  ! A is specified through nrow,ncol,nz,irow,pcol,val
  ! results are returned in b
  SUBROUTINE sparse_solveComplex_b2(nrow,ncol,nz,irow,pcol,val,b,iopt_in)
    INTEGER, INTENT(in) :: nrow,ncol,nz
    INTEGER, DIMENSION(:), INTENT(in) :: irow,pcol
    COMPLEX(kind=dp), DIMENSION(:), INTENT(in) :: val
    COMPLEX(kind=dp), DIMENSION(:,:), INTENT(inout) :: b
    INTEGER, INTENT(in), OPTIONAL :: iopt_in
    
    INTEGER :: iopt = 0
    INTEGER :: n
    INTEGER, DIMENSION(:), ALLOCATABLE :: pcoln
    LOGICAL :: pcol_modified
    
    IF (PRESENT(iopt_in)) iopt = iopt_in
    
    pcol_modified = .FALSE.
    ! pcoln to remove "C convention" used in SuiteSparse
    IF (pcol(1) .EQ. 0) THEN
       pcol_modified = .TRUE.
       ALLOCATE(pcoln(ncol+1))
       pcoln = pcol + 1
    END IF
    
    ! check about existing factorization
    IF (factorization_exists .AND. iopt .EQ. 1) THEN ! free memory first
       IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. 3) ) THEN
          IF (pcol_modified) THEN
             CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcoln,val,b,3)
          ELSE
             CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,3)
          END IF
       END IF
    END IF
    IF (.NOT. factorization_exists .AND. iopt .EQ. 2) THEN ! factorize first
       IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. 3) ) THEN
          IF (pcol_modified) THEN
             CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcoln,val,b,1)
          ELSE
             CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,1)
          END IF
       END IF
       factorization_exists = .TRUE.
    END IF
    IF (iopt .EQ. 1) factorization_exists = .TRUE.
    IF (iopt .EQ. 3) factorization_exists = .FALSE.
    
    IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. 3) ) THEN
       IF (pcol_modified) THEN
          CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcoln,val,b,iopt)
       ELSE
          CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,iopt)
       END IF
    ELSE
       PRINT *, 'sparse_solve_method ',sparse_solve_method,'not implemented'
       STOP
    END IF
    
    IF (ALLOCATED(pcoln)) DEALLOCATE(pcoln)
    
  END SUBROUTINE sparse_solveComplex_b2
  !-------------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------------
  ! solves A*x = b for sparse A and 1-D vector b
  ! A is given as a full matrix
  ! results are returned in b
  SUBROUTINE sparse_solveReal_A_b1(A,b,iopt_in)
    REAL(kind=dp), DIMENSION(:,:), INTENT(in) :: A
    REAL(kind=dp), DIMENSION(:), INTENT(inout) :: b
    INTEGER, INTENT(in), OPTIONAL :: iopt_in
    
    INTEGER :: iopt = 0
    INTEGER :: nrow,ncol,nz
    INTEGER, DIMENSION(:), ALLOCATABLE :: irow,pcol
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: val
    
    IF (PRESENT(iopt_in)) iopt = iopt_in
    
    ! check about existing factorization
    IF (factorization_exists .AND. iopt .EQ. 1) THEN ! free memory first
       IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. 3) ) THEN
          CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,3)
       END IF
    END IF
    IF (.NOT. factorization_exists .AND. iopt .EQ. 2) THEN ! factorize first
       IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. 3) ) THEN
          CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,1)
       END IF
       factorization_exists = .TRUE.
    END IF
    IF (iopt .EQ. 1) factorization_exists = .TRUE.
    IF (iopt .EQ. 3) factorization_exists = .FALSE.
    
    IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. 3) ) THEN
       CALL full2sparse(A,irow,pcol,val,nrow,ncol,nz)
       CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,iopt)
    ELSE
       PRINT *, 'sparse_solve_method ',sparse_solve_method,'not implemented'
       STOP
    END IF
    
    IF (ALLOCATED(irow)) DEALLOCATE(irow)
    IF (ALLOCATED(pcol)) DEALLOCATE(pcol)
    IF (ALLOCATED(val))  DEALLOCATE(val)
    
  END SUBROUTINE sparse_solveReal_A_b1
  !-------------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------------
  ! solves A*x = b for sparse A and 1-D vector b
  ! A is given as a full matrix
  ! results are returned in b
  SUBROUTINE sparse_solveComplex_A_b1(A,b,iopt_in)
    COMPLEX(kind=dp), DIMENSION(:,:), INTENT(in) :: A
    COMPLEX(kind=dp), DIMENSION(:), INTENT(inout) :: b
    INTEGER, INTENT(in), OPTIONAL :: iopt_in
    
    INTEGER :: iopt = 0
    INTEGER :: nrow,ncol,nz
    INTEGER, DIMENSION(:), ALLOCATABLE :: irow,pcol
    COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE :: val
    
    IF (PRESENT(iopt_in)) iopt = iopt_in
    
    ! check about existing factorization
    IF (factorization_exists .AND. iopt .EQ. 1) THEN ! free memory first
       IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. 3) ) THEN
          CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,3)
       END IF
    END IF
    IF (.NOT. factorization_exists .AND. iopt .EQ. 2) THEN ! factorize first
       IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. 3) ) THEN
          CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,1)
       END IF
       factorization_exists = .TRUE.
    END IF
    IF (iopt .EQ. 1) factorization_exists = .TRUE.
    IF (iopt .EQ. 3) factorization_exists = .FALSE.
    
    IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. 3) ) THEN
       CALL full2sparse(A,irow,pcol,val,nrow,ncol,nz)
       CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,iopt)
    ELSE
       PRINT *, 'sparse_solve_method ',sparse_solve_method,'not implemented'
       STOP
    END IF
    
    IF (ALLOCATED(irow)) DEALLOCATE(irow)
    IF (ALLOCATED(pcol)) DEALLOCATE(pcol)
    IF (ALLOCATED(val))  DEALLOCATE(val)
    
  END SUBROUTINE sparse_solveComplex_A_b1
  !-------------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------------
  ! solves A*x = b for sparse A and 2-D array b
  ! A is given as a full matrix
  ! results are returned in b
  SUBROUTINE sparse_solveReal_A_b2(A,b,iopt_in)
    REAL(kind=dp), DIMENSION(:,:), INTENT(in) :: A
    REAL(kind=dp), DIMENSION(:,:), INTENT(inout) :: b
    INTEGER, INTENT(in), OPTIONAL :: iopt_in
    
    INTEGER :: iopt = 0
    INTEGER :: nrow,ncol,nz
    INTEGER, DIMENSION(:), ALLOCATABLE :: irow,pcol
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: val
    
    IF (PRESENT(iopt_in)) iopt = iopt_in
    
    ! check about existing factorization
    IF (factorization_exists .AND. iopt .EQ. 1) THEN ! free memory first
       IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. 3) ) THEN
          CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,3)
       END IF
    END IF
    IF (.NOT. factorization_exists .AND. iopt .EQ. 2) THEN ! factorize first
       IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. 3) ) THEN
          CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,1)
       END IF
       factorization_exists = .TRUE.
    END IF
    IF (iopt .EQ. 1) factorization_exists = .TRUE.
    IF (iopt .EQ. 3) factorization_exists = .FALSE.
    
    IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. 3) ) THEN
       CALL full2sparse(A,irow,pcol,val,nrow,ncol,nz)
       CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,iopt)
    ELSE
       PRINT *, 'sparse_solve_method ',sparse_solve_method,'not implemented'
       STOP
    END IF
    
    IF (ALLOCATED(irow)) DEALLOCATE(irow)
    IF (ALLOCATED(pcol)) DEALLOCATE(pcol)
    IF (ALLOCATED(val))  DEALLOCATE(val)
    
  END SUBROUTINE sparse_solveReal_A_b2
  !-------------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------------
  ! solves A*x = b for sparse A and 2-D array b
  ! A is given as a full matrix
  ! results are returned in b
  SUBROUTINE sparse_solveComplex_A_b2(A,b,iopt_in)
    COMPLEX(kind=dp), DIMENSION(:,:), INTENT(in) :: A
    COMPLEX(kind=dp), DIMENSION(:,:), INTENT(inout) :: b
    INTEGER, INTENT(in), OPTIONAL :: iopt_in
    
    INTEGER :: iopt = 0
    INTEGER :: nrow,ncol,nz
    INTEGER, DIMENSION(:), ALLOCATABLE :: irow,pcol
    COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE :: val
    
    IF (PRESENT(iopt_in)) iopt = iopt_in
    
    ! check about existing factorization
    IF (factorization_exists .AND. iopt .EQ. 1) THEN ! free memory first
       IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. 3) ) THEN
          CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,3)
       END IF
    END IF
    IF (.NOT. factorization_exists .AND. iopt .EQ. 2) THEN ! factorize first
       IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. 3) ) THEN
          CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,1)
       END IF
       factorization_exists = .TRUE.
    END IF
    IF (iopt .EQ. 1) factorization_exists = .TRUE.
    IF (iopt .EQ. 3) factorization_exists = .FALSE.
    
    IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. 3) ) THEN
       CALL full2sparse(A,irow,pcol,val,nrow,ncol,nz)
       CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,iopt)
    ELSE
       PRINT *, 'sparse_solve_method ',sparse_solve_method,'not implemented'
       STOP
    END IF
    
    IF (ALLOCATED(irow)) DEALLOCATE(irow)
    IF (ALLOCATED(pcol)) DEALLOCATE(pcol)
    IF (ALLOCATED(val))  DEALLOCATE(val)
    
  END SUBROUTINE sparse_solveComplex_A_b2
  !-------------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------------
  ! solves A*x = b for sparse A and 1-D vector b
  ! A is specified through nrow,ncol,nz,irow,pcol,val
  ! results are returned in b
  ! Routines from SuperLU-Distribution
  SUBROUTINE sparse_solve_suitesparse_b1(nrow,ncol,nz,irow,pcol,val,b,iopt_in)
    INTEGER, INTENT(in) :: nrow,ncol,nz
    INTEGER, DIMENSION(:), INTENT(in) :: irow,pcol
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: val
    REAL(kind=dp), DIMENSION(:), INTENT(inout) :: b
    INTEGER, INTENT(in) :: iopt_in
    
    INTEGER :: n,nc
    INTEGER :: nrhs = 1
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: x
    
    INTEGER :: Ap_len, Ai_len, Ax_len, Az_len
    INTEGER, DIMENSION(:), ALLOCATABLE :: Ap, Ai, p
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: Ax, Az
    
    ! C wrappers of SuiteSparse functions
    INTEGER(kind=long) :: umf4def_
    INTEGER(kind=long) :: umf4sym_
    INTEGER(kind=long) :: umf4num_
    INTEGER(kind=long) :: umf4solr_
    INTEGER(kind=long) :: umf4sol_
    INTEGER(kind=long) :: umf4fnum_
    INTEGER(kind=long) :: umf4fsym_
    
    Ap_len = ncol + 1
    Ai_len = nz
    Ax_len = nz
    Az_len = 0 ! nz
    ALLOCATE(Ap(Ap_len), Ai(Ai_len))
    Ap = pcol
    Ai = irow
    
    IF (iopt_in .EQ. 3) THEN  ! free memory from last solution
       CALL umf4fnum(numeric)
       factorization_exists = .FALSE.
       DEALLOCATE(Ap, Ai)
       RETURN
    END IF
    
    ALLOCATE(x(nrow))
    
    IF (iopt_in .EQ. 0 .OR. iopt_in .EQ. 1) THEN
       IF ( sparse_solve_method .EQ. 3 ) THEN ! SuiteSparse (without (=3) iterative refinement)
          CALL umf4def(control)
          control(1) = 0 ! No output - there are other options, see the manual
       ELSE IF ( sparse_solve_method .EQ. 2 ) THEN ! SuiteSparse (with (=2) iterative refinement)
          CALL umf4def(control)
          control(1) = 0 ! No output - there are other options, see the manual
          control(8) = 10 ! max number of iterative refinement steps
       END IF
       
       ! Pre-order and symbolic analysis
       CALL umf4sym (n, n, Ap, Ai, val, symbolic, control, info_suitesparse)
       IF (info_suitesparse(1) .LT. 0) THEN
          PRINT *, 'Error occurred in umf4sym: ', info_suitesparse(1)
       END IF
       
       CALL umf4num(Ap, Ai, val, symbolic, numeric, control, info_suitesparse)
       IF (info_suitesparse(1) .LT. 0) THEN
          PRINT *, 'Error occurred in umf4num: ', info_suitesparse(1)
       END IF
    END IF
    
    ! Note: in newer versions, the function interfaces has been changed to match the types in an cleaner way
    !       use n instead of nrow
    !       use nc instead of ncol
    n = nrow
    nc = ncol
    
    IF (iopt_in .EQ. 0 .OR. iopt_in .EQ. 2) THEN
       IF ( sparse_solve_method .EQ. 2 ) THEN ! SuiteSparse (with (=2)
          CALL umf4solr (sys, Ap, Ai, val, x, b, numeric, control, info_suitesparse) !iterative refinement
       ELSE IF ( sparse_solve_method .EQ. 3 ) THEN ! SuiteSparse (without (=3) iterative refinement)
          CALL umf4sol (sys, x, b, numeric, control, info_suitesparse)
       END IF
       IF (info_suitesparse(1) .LT. 0) THEN
          PRINT *, 'Error occurred in umf4solr: ', info_suitesparse(1)
       END IF
    END IF
    
    b = x
    
    ! Last, free the storage allocated inside SuiteSparse
    IF (iopt_in .EQ. 0 .OR. iopt_in .EQ. 3) THEN
       CALL umf4fnum (numeric)
       CALL umf4fsym (symbolic)
    END IF
    
    IF (ALLOCATED(Ap))  DEALLOCATE(Ap)
    IF (ALLOCATED(Ai))  DEALLOCATE(Ai)
    IF (ALLOCATED(x))  DEALLOCATE(x)
    
  END SUBROUTINE sparse_solve_suitesparse_b1
  !-------------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------------
  ! solves A*x = b for sparse A and 1-D vector b
  ! A is specified through nrow,ncol,nz,irow,pcol,val
  ! results are returned in b
  ! Routines from SuperLU-Distribution
  SUBROUTINE sparse_solve_suitesparseComplex_b1(nrow,ncol,nz,irow,pcol,val,b,iopt_in)
    INTEGER, INTENT(in) :: nrow,ncol,nz
    INTEGER, DIMENSION(:), INTENT(in) :: irow,pcol
    COMPLEX(kind=dp), DIMENSION(:), INTENT(in) :: val
    COMPLEX(kind=dp), DIMENSION(:), INTENT(inout) :: b
    INTEGER, INTENT(in) :: iopt_in
    
    INTEGER :: n,nc
    INTEGER :: nrhs = 1
    COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE :: x
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: xx, xz, bx, bz
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: valx, valz
    
    INTEGER :: Ap_len, Ai_len, Ax_len, Az_len
    INTEGER, DIMENSION(:), ALLOCATABLE :: Ap, Ai, p
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: Ax, Az
    
    ! C wrappers of SuiteSparse functions
    INTEGER(kind=long) :: umf4zdef_
    INTEGER(kind=long) :: umf4zsym_
    INTEGER(kind=long) :: umf4znum_
    INTEGER(kind=long) :: umf4zsolr_
    INTEGER(kind=long) :: umf4zsol_
    INTEGER(kind=long) :: umf4zfnum_
    INTEGER(kind=long) :: umf4zfsym_
    
    Ap_len = ncol + 1
    Ai_len = nz
    Ax_len = nz
    Az_len = nz
    ALLOCATE(Ap(Ap_len), Ai(Ai_len))
    Ap = pcol
    Ai = irow
    
    ALLOCATE(valx(nz), valz(nz))
    valx = REAL(val)
    valz = AIMAG(val)
    
    ALLOCATE(bx(nrow), bz(nrow))
    bx = REAL(b)
    bz = AIMAG(b)
    
    IF (iopt_in .EQ. 3) THEN  ! free memory from last solution
       CALL umf4zfnum(numeric)
       factorization_exists = .FALSE.
       DEALLOCATE(Ap, Ai, valx, valz, bx, bz)
       RETURN
    END IF
    
    ALLOCATE(xx(nrow), xz(nrow))
    
    IF (iopt_in .EQ. 0 .OR. iopt_in .EQ. 1) THEN
       IF ( sparse_solve_method .EQ. 3 ) THEN ! SuiteSparse (without (=3) iterative refinement)
          CALL umf4zdef(control)
          control(1) = 0 ! No output - there are other options, see the manual
       ELSE IF ( sparse_solve_method .EQ. 2 ) THEN ! SuiteSparse (with (=2) iterative refinement)
          CALL umf4zdef(control)
          control(1) = 0 ! No output - there are other options, see the manual
          control(8) = 10 ! max number of iterative refinement steps
       END IF
       ! Pre-order and symbolic analysis
       CALL umf4zsym (n, n, Ap, Ai, valx, valz, symbolic, control, info_suitesparse)
       IF (info_suitesparse(1) .LT. 0) THEN
          PRINT *, 'Error occurred in umf4zsym: ', info_suitesparse(1)
       END IF
       
       CALL umf4znum(Ap, Ai, valx, valz, symbolic, numeric, control, info_suitesparse)
       IF (info_suitesparse(1) .LT. 0) THEN
          PRINT *, 'Error occurred in umf4znum: ', info_suitesparse(1)
       END IF
    END IF
    
    ! Note: in newer versions, the function interfaces has been changed to match the types in an cleaner way
    !       use n instead of nrow
    !       use nc instead of ncol
    n = nrow
    nc = ncol
    
    IF (iopt_in .EQ. 0 .OR. iopt_in .EQ. 2) THEN
       IF ( sparse_solve_method .EQ. 2 ) THEN ! SuiteSparse (with (=2)
          CALL umf4zsolr (sys, Ap, Ai, valx, valz, xx, xz, bx, bz, numeric, &
               control, info_suitesparse) !iterative refinement
       ELSE IF ( sparse_solve_method .EQ. 3 ) THEN ! SuiteSparse (without (=3) iterative refinement)
          CALL umf4zsol (sys, xx, xz, bx, bz, numeric, control, info_suitesparse)
       END IF
       IF (info_suitesparse(1) .LT. 0) THEN
          PRINT *, 'Error occurred in umf4zsolr: ', info_suitesparse(1)
       END IF
    END IF
    
    b = CMPLX(xx, xz, KIND=dp)
    
    ! Last, free the storage allocated inside SuiteSparse
    IF (iopt_in .EQ. 0 .OR. iopt_in .EQ. 3) THEN
       CALL umf4zfnum (numeric)
       CALL umf4zfsym (symbolic)
    END IF
    
    IF (ALLOCATED(Ap))   DEALLOCATE(Ap)
    IF (ALLOCATED(Ai))   DEALLOCATE(Ai)
    IF (ALLOCATED(xx))   DEALLOCATE(xx)
    IF (ALLOCATED(xz))   DEALLOCATE(xz)
    IF (ALLOCATED(bx))   DEALLOCATE(bx)
    IF (ALLOCATED(bz))   DEALLOCATE(bz)
    IF (ALLOCATED(valx))  DEALLOCATE(valx)
    IF (ALLOCATED(valz))  DEALLOCATE(valz)
    
  END SUBROUTINE sparse_solve_suitesparseComplex_b1
  !-------------------------------------------------------------------------------
  
  
  !-------------------------------------------------------------------------------
  ! solves A*x = b for sparse A and 2-D array b
  ! A is specified through nrow,ncol,nz,irow,pcol,val
  ! results are returned in b
  ! Routines from SuiteSparse-Distribution
  SUBROUTINE sparse_solve_suitesparse_b2_loop(nrow,ncol,nz,irow,pcol,val,b,iopt_in)
    INTEGER, INTENT(in) :: nrow,ncol,nz
    INTEGER, DIMENSION(:), INTENT(in) :: irow,pcol
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: val
    REAL(kind=dp), DIMENSION(:,:), INTENT(inout) :: b
    INTEGER, INTENT(in) :: iopt_in
    
    INTEGER :: n,nc,i,nrhs
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: x,bloc
    
    INTEGER :: Ap_len, Ai_len, Ax_len, Az_len
    INTEGER, DIMENSION(:), ALLOCATABLE :: Ap, Ai, p
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: Ax, Az
    
    ! C wrappers of SuiteSparse functions
    INTEGER(kind=long) :: umf4def_
    INTEGER(kind=long) :: umf4sym_
    INTEGER(kind=long) :: umf4num_
    INTEGER(kind=long) :: umf4solr_
    INTEGER(kind=long) :: umf4sol_
    INTEGER(kind=long) :: umf4fnum_
    INTEGER(kind=long) :: umf4fsym_
    
    Ap_len = ncol + 1
    Ai_len = nz
    Ax_len = nz
    Az_len = 0 ! nz
    ALLOCATE(Ap(Ap_len), Ai(Ai_len))
    Ap = pcol
    Ai = irow
    
    IF (iopt_in .EQ. 3) THEN  ! free memory from last solution
       CALL umf4fnum(numeric)
       factorization_exists = .FALSE.
       DEALLOCATE(Ap, Ai)
       RETURN
    END IF
    
    nrhs = SIZE(b,2)
    ALLOCATE(x(nrow), bloc(nrow))
    
    IF (iopt_in .EQ. 0 .OR. iopt_in .EQ. 1) THEN
       IF ( sparse_solve_method .EQ. 3 ) THEN ! SuiteSparse (without (=3) iterative refinement)
          CALL umf4def(control)
          control(1) = 0 ! No output - there are other options, see the manual
       ELSE IF ( sparse_solve_method .EQ. 2 ) THEN ! SuiteSparse (with (=2) iterative refinement)
          CALL umf4def(control)
          control(1) = 0 ! No output - there are other options, see the manual
          control(8) = 10 ! max number of iterative refinement steps
       END IF
       
       ! Pre-order and symbolic analysis
       CALL umf4sym (n, n, Ap, Ai, val, symbolic, control, info_suitesparse)
       IF (info_suitesparse(1) .LT. 0) THEN
          PRINT *, 'Error occurred in umf4sym: ', info_suitesparse(1)
       END IF
       
       CALL umf4num(Ap, Ai, val, symbolic, numeric, control, info_suitesparse)
       IF (info_suitesparse(1) .LT. 0) THEN
          PRINT *, 'Error occurred in umf4num: ', info_suitesparse(1)
       END IF
    END IF
    
    ! Note: in newer versions, the function interfaces has been changed to match the types in an cleaner way
    !       use n instead of nrow
    !       use nc instead of ncol
    n = nrow
    nc = ncol
    
    IF (iopt_in .EQ. 0 .OR. iopt_in .EQ. 2) THEN
       DO i = 1,nrhs
          bloc = b(:,i)
          IF ( sparse_solve_method .EQ. 2 ) THEN ! SuiteSparse (with (=2)
             CALL umf4solr (sys, Ap, Ai, val, x, bloc, numeric, control, info_suitesparse) !iterative refinement
          ELSE IF ( sparse_solve_method .EQ. 3 ) THEN ! SuiteSparse (without (=3) iterative refinement)
             CALL umf4sol (sys, x, bloc, numeric, control, info_suitesparse)
          END IF
          IF (info_suitesparse(1) .LT. 0) THEN
             PRINT *, 'Error occurred in umf4solr: ', info_suitesparse(1)
          END IF
          b(:,i) = x
       END DO
    END IF
    
    ! Last, free the storage allocated inside SuiteSparse
    IF (iopt_in .EQ. 0 .OR. iopt_in .EQ. 3) THEN
       CALL umf4fnum (numeric)
       CALL umf4fsym (symbolic)
    END IF
    
    IF (ALLOCATED(Ap))  DEALLOCATE(Ap)
    IF (ALLOCATED(Ai))  DEALLOCATE(Ai)
    IF (ALLOCATED(x))  DEALLOCATE(x)
    IF (ALLOCATED(bloc))  DEALLOCATE(bloc)
    
    RETURN
  END SUBROUTINE sparse_solve_suitesparse_b2_loop
  !-------------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------------
  ! solves A*x = b for sparse A and 2-D array b
  ! A is specified through nrow,ncol,nz,irow,pcol,val
  ! results are returned in b
  ! Routines from SuiteSparse-Distribution
  SUBROUTINE sparse_solve_suitesparseComplex_b2_loop(nrow,ncol,nz,irow,pcol,val,b,iopt_in)
    INTEGER, INTENT(in) :: nrow,ncol,nz
    INTEGER, DIMENSION(:), INTENT(in) :: irow,pcol
    COMPLEX(kind=dp), DIMENSION(:), INTENT(in) :: val
    COMPLEX(kind=dp), DIMENSION(:,:), INTENT(inout) :: b
    INTEGER, INTENT(in) :: iopt_in
    
    INTEGER :: n,nc,i,nrhs
    COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE :: x,bloc
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: xx, xz, blocx, blocz
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: valx, valz
    
    INTEGER :: Ap_len, Ai_len, Ax_len, Az_len
    INTEGER, DIMENSION(:), ALLOCATABLE :: Ap, Ai, p
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: Ax, Az
    
    ! C wrappers of SuiteSparse functions
    INTEGER(kind=long) :: umf4zdef_
    INTEGER(kind=long) :: umf4zsym_
    INTEGER(kind=long) :: umf4znum_
    INTEGER(kind=long) :: umf4zsolr_
    INTEGER(kind=long) :: umf4zsol_
    INTEGER(kind=long) :: umf4zfnum_
    INTEGER(kind=long) :: umf4zfsym_
    
    Ap_len = ncol + 1
    Ai_len = nz
    Ax_len = nz
    Az_len = nz
    ALLOCATE(Ap(Ap_len), Ai(Ai_len))
    Ap = pcol
    Ai = irow
    
    ALLOCATE(valx(nz), valz(nz))
    valx = REAL(val)
    valz = AIMAG(val)
    
    IF (iopt_in .EQ. 3) THEN  ! free memory from last solution
       CALL umf4zfnum(numeric)
       factorization_exists = .FALSE.
       DEALLOCATE(Ap, Ai, valx, valz)
       RETURN
    END IF
    
    nrhs = SIZE(b,2)
    ALLOCATE(xx(nrow), xz(nrow), blocx(nrow), blocz(nrow))
    
    IF (iopt_in .EQ. 0 .OR. iopt_in .EQ. 1) THEN
       IF ( sparse_solve_method .EQ. 3 ) THEN ! SuiteSparse (without (=3) iterative refinement)
          CALL umf4zdef(control)
          control(1) = 0 ! No output - there are other options, see the manual
       ELSE IF ( sparse_solve_method .EQ. 2 ) THEN ! SuiteSparse (with (=2) iterative refinement)
          CALL umf4zdef(control)
          control(1) = 0 ! No output - there are other options, see the manual
          control(8) = 10 ! max number of iterative refinement steps
       END IF
       ! Pre-order and symbolic analysis
       CALL umf4zsym (n, n, Ap, Ai, valx, valz, symbolic, control, info_suitesparse)
       IF (info_suitesparse(1) .LT. 0) THEN
          PRINT *, 'Error occurred in umf4zsym: ', info_suitesparse(1)
       END IF
       
       CALL umf4znum(Ap, Ai, valx, valz, symbolic, numeric, control, info_suitesparse)
       IF (info_suitesparse(1) .LT. 0) THEN
          PRINT *, 'Error occurred in umf4znum: ', info_suitesparse(1)
       END IF
    END IF
    
    ! Note: in newer versions, the function interfaces has been changed to match the types in an cleaner way
    !       use n instead of nrow
    !       use nc instead of ncol
    n = nrow
    nc = ncol
    
    IF (iopt_in .EQ. 0 .OR. iopt_in .EQ. 2) THEN
       DO i = 1,nrhs
          blocx = REAL(b(:,i))
          blocz = AIMAG(b(:,i))
          IF ( sparse_solve_method .EQ. 2 ) THEN ! SuiteSparse (with (=2)
             CALL umf4zsolr (sys, Ap, Ai, valx, valz, xx, xz, blocx, blocz, numeric,&
                  control, info_suitesparse) !iterative refinement
          ELSE IF ( sparse_solve_method .EQ. 3 ) THEN ! SuiteSparse (without (=3) iterative refinement)
             CALL umf4zsol (sys, xx, xz, blocx, blocz, numeric, control, info_suitesparse)
          END IF
          IF (info_suitesparse(1) .LT. 0) THEN
             PRINT *, 'Error occurred in umf4zsolr: ', info_suitesparse(1)
          END IF
          b(:,i) = CMPLX(xx, xz, KIND=dp)
       END DO
    END IF
    
    ! Last, free the storage allocated inside SuiteSparse
    IF (iopt_in .EQ. 0 .OR. iopt_in .EQ. 3) THEN
       CALL umf4zfnum (numeric)
       CALL umf4zfsym (symbolic)
    END IF
    
    IF (ALLOCATED(Ap))    DEALLOCATE(Ap)
    IF (ALLOCATED(Ai))    DEALLOCATE(Ai)
    IF (ALLOCATED(xx))    DEALLOCATE(xx)
    IF (ALLOCATED(xz))    DEALLOCATE(xz)
    IF (ALLOCATED(blocx)) DEALLOCATE(blocx)
    IF (ALLOCATED(blocz)) DEALLOCATE(blocz)
    IF (ALLOCATED(valx))  DEALLOCATE(valx)
    IF (ALLOCATED(valz))  DEALLOCATE(valz)
    
  END SUBROUTINE sparse_solve_suitesparseComplex_b2_loop
  !-------------------------------------------------------------------------------
  
END MODULE sparse_solvers_mod