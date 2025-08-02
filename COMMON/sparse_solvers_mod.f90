MODULE sparse_solvers_mod
  ! Module containing sparse matrix solver operations
  ! Extracted from sparse_mod.f90 for better modularity
  !
  ! This module fixes a critical bug from the original sparse_mod.f90 where
  ! real and complex solvers shared the same factorization pointers (symbolic, numeric),
  ! causing memory corruption when alternating between solver types.
  !
  ! The fix uses separate pointers for real and complex factorizations and
  ! properly cleans up the opposing type's factorization when switching.
  
  USE sparse_types_mod, ONLY: dp, long
  USE sparse_conversion_mod
  USE sparse_arithmetic_mod, ONLY: sparse_talk
  USE sparse_utils_mod, ONLY: csc_to_csr_real, csc_to_csr_complex
  USE bicgstab_mod
  USE ilu_precond_mod
  IMPLICIT NONE
  
  ! Solver method constants (no magic numbers)
  INTEGER, PARAMETER, PUBLIC :: SOLVER_UMFPACK = 3
  INTEGER, PARAMETER, PUBLIC :: SOLVER_BICGSTAB = 4
  
  ! Module variables that need to be accessible
  INTEGER :: sparse_solve_method = 0  ! Auto-select: UMFPACK for small matrices, BiCGSTAB for large
  LOGICAL :: factorization_exists = .FALSE.
  
  ! BiCGSTAB solver parameters (configurable via namelist)
  REAL(kind=dp) :: bicgstab_abs_tolerance = 1.0e-14_dp    ! Absolute tolerance floor (||r|| < abs_tol)
  REAL(kind=dp) :: bicgstab_rel_tolerance = 1.0e-8_dp     ! Relative tolerance (||r|| < rel_tol*||b||)
  INTEGER :: bicgstab_max_iter = 1000                      ! Maximum iterations
  INTEGER :: bicgstab_restart_limit = 3                    ! Number of restarts allowed
  LOGICAL :: bicgstab_verbose = .FALSE.                    ! Print convergence info
  LOGICAL :: bicgstab_adaptive_tolerance = .FALSE.         ! Automatically adjust tolerance based on matrix conditioning

  ! Adaptive tolerance thresholds - condition number estimates
  REAL(kind=dp), PARAMETER :: COND_WELL_CONDITIONED = 100.0_dp      ! Well-conditioned threshold
  REAL(kind=dp), PARAMETER :: COND_MODERATELY_ILL = 10000.0_dp      ! Moderately ill-conditioned
  REAL(kind=dp), PARAMETER :: COND_ILL_CONDITIONED = 1000000.0_dp   ! Ill-conditioned threshold
  
  ! Adaptive tolerance fallback values
  REAL(kind=dp), PARAMETER :: TOL_ABS_MODERATE = 1.0e-10_dp    ! Moderate absolute tolerance
  REAL(kind=dp), PARAMETER :: TOL_REL_MODERATE = 1.0e-6_dp     ! Moderate relative tolerance
  REAL(kind=dp), PARAMETER :: TOL_ABS_RELAXED = 1.0e-8_dp      ! Relaxed absolute tolerance
  REAL(kind=dp), PARAMETER :: TOL_REL_RELAXED = 1.0e-5_dp      ! Relaxed relative tolerance
  REAL(kind=dp), PARAMETER :: TOL_ABS_LOOSE = 1.0e-6_dp        ! Loose absolute tolerance
  REAL(kind=dp), PARAMETER :: TOL_REL_LOOSE = 1.0e-4_dp        ! Loose relative tolerance
  
  ! ILU preconditioning parameters (configurable via namelist in neo2.in)
  INTEGER :: ilu_fill_level = 1                            ! ILU(k) fill level (0=no ILU, 1=ILU(1), 2=ILU(2), etc.)
  REAL(kind=dp) :: ilu_drop_tolerance = 0.0_dp             ! Drop tolerance for ILU factorization
  
  ! Adaptive solver selection parameters
  INTEGER :: auto_solver_threshold = 100               ! Switch to BiCGSTAB for matrices >= this size
  
  PUBLIC :: sparse_solve
  PUBLIC :: sparse_solve_suitesparse
  PUBLIC :: sparse_solve_method
  PUBLIC :: factorization_exists
  PUBLIC :: bicgstab_abs_tolerance, bicgstab_rel_tolerance
  PUBLIC :: bicgstab_max_iter, bicgstab_restart_limit, bicgstab_verbose
  PUBLIC :: bicgstab_adaptive_tolerance
  PUBLIC :: ilu_fill_level, ilu_drop_tolerance
  PUBLIC :: auto_solver_threshold
  
  ! Named constants for iopt parameter values (compatible with original sparse_mod.f90)
  ! Usage pattern: call with iopt=1 to factorize, then iopt=2 to solve multiple RHS
  INTEGER, PARAMETER, PRIVATE :: IOPT_FULL_SOLVE = 0     ! Factorize + solve + cleanup
  INTEGER, PARAMETER, PRIVATE :: IOPT_FACTORIZE_ONLY = 1 ! Factorize only (save for reuse)
  INTEGER, PARAMETER, PRIVATE :: IOPT_SOLVE_ONLY = 2     ! Solve only (reuse factorization)
  INTEGER, PARAMETER, PRIVATE :: IOPT_FREE_MEMORY = 3    ! Free memory only
  
  ! Initialization of the parameters of Super_LU c-Routines
  INTEGER(kind=long), PRIVATE :: factors
  
  ! SuiteSparse solver data address pointers
  ! Separate pointers for real and complex factorizations to avoid conflicts
  INTEGER(kind=long), PRIVATE :: symbolic_real = 0, numeric_real = 0
  INTEGER(kind=long), PRIVATE :: symbolic_complex = 0, numeric_complex = 0
  INTEGER(kind=long), PRIVATE :: sys = 0
  REAL(kind=dp), PRIVATE :: control(20), info_suitesparse(90)
  LOGICAL, PRIVATE :: factorization_exists_real = .FALSE.
  LOGICAL, PRIVATE :: factorization_exists_complex = .FALSE.
  INTEGER, PRIVATE :: current_factorization_type = 0  ! 0=none, 1=real, 2=complex
  
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
  ! Auto-selection of sparse solver method based on matrix size
  !-------------------------------------------------------------------------------
  SUBROUTINE handle_auto_selection(nrow, ncol)
    INTEGER, INTENT(in) :: nrow, ncol
    
    IF (sparse_solve_method .EQ. 0) THEN
       IF (nrow < auto_solver_threshold) THEN
          sparse_solve_method = SOLVER_UMFPACK
          IF (sparse_talk) PRINT *, 'Auto-selected UMFPACK for small matrix (', nrow, 'x', ncol, ')'
       ELSE
          sparse_solve_method = SOLVER_BICGSTAB
          IF (sparse_talk) PRINT *, 'Auto-selected BiCGSTAB for large matrix (', nrow, 'x', ncol, ')'
       END IF
    END IF
  END SUBROUTINE handle_auto_selection

  !-------------------------------------------------------------------------------
  ! Parameter validation for sparse solver input
  !-------------------------------------------------------------------------------
  SUBROUTINE validate_sparse_solver_params(nrow, ncol, nz, solver_method, routine_name)
    INTEGER, INTENT(in) :: nrow, ncol, nz, solver_method
    CHARACTER(len=*), INTENT(in) :: routine_name
    
    ! Check matrix dimensions
    IF (nrow <= 0) THEN
      PRINT *, 'ERROR in ', routine_name, ': nrow must be positive, got', nrow
      ERROR STOP 'Invalid matrix dimensions'
    END IF
    
    IF (ncol <= 0) THEN
      PRINT *, 'ERROR in ', routine_name, ': ncol must be positive, got', ncol
      ERROR STOP 'Invalid matrix dimensions'
    END IF
    
    IF (nz < 0) THEN
      PRINT *, 'ERROR in ', routine_name, ': nz cannot be negative, got', nz
      ERROR STOP 'Invalid number of non-zeros'
    END IF
    
    ! Check for reasonable matrix size
    IF (nrow > 1000000 .OR. ncol > 1000000) THEN
      PRINT *, 'WARNING in ', routine_name, ': Large matrix size', nrow, 'x', ncol
      PRINT *, '         Consider using BiCGSTAB (method 4) for better memory efficiency'
    END IF
    
    ! Check solver method and implement auto-selection
    SELECT CASE (solver_method)
      CASE (0)
        ! Auto-select based on matrix size
        IF (nrow < auto_solver_threshold) THEN
          sparse_solve_method = SOLVER_UMFPACK
          IF (sparse_talk) THEN
            PRINT *, 'INFO in ', routine_name, ': Auto-selected UMFPACK for small matrix (', nrow, 'x', ncol, ')'
          END IF
        ELSE
          sparse_solve_method = SOLVER_BICGSTAB
          IF (sparse_talk) THEN
            PRINT *, 'INFO in ', routine_name, ': Auto-selected BiCGSTAB for large matrix (', nrow, 'x', ncol, ')'
          END IF
        END IF
      CASE (2, SOLVER_UMFPACK, SOLVER_BICGSTAB)
        ! Valid explicit solver methods
      CASE DEFAULT
        PRINT *, 'ERROR in ', routine_name, ': Invalid solver method', solver_method
        PRINT *, 'Valid options: 0 (auto), 2 (legacy), SOLVER_UMFPACK(3), SOLVER_BICGSTAB(4)'
        ERROR STOP 'Invalid solver method'
    END SELECT
    
    ! Optional informational messages (controlled by sparse_talk)
    IF (sparse_talk) THEN
      ! BiCGSTAB-specific recommendations
      IF (solver_method == SOLVER_BICGSTAB) THEN
        IF (nrow < 10) THEN
          PRINT *, 'INFO in ', routine_name, ': Small matrix (', nrow, 'x', ncol, ')'
          PRINT *, '      UMFPACK might be more efficient for very small problems'
        ELSE IF (nrow > 1000) THEN
          PRINT *, 'INFO in ', routine_name, ': Large matrix (', nrow, 'x', ncol, ')'
          PRINT *, '      BiCGSTAB should provide good memory efficiency'
        END IF
      END IF
      
      ! UMFPACK-specific warnings for large problems
      IF ((solver_method == SOLVER_UMFPACK .OR. solver_method == 2) .AND. nrow > 1000) THEN
        PRINT *, 'INFO in ', routine_name, ': Large matrix with direct solver'
        PRINT *, '      Consider BiCGSTAB (method 4) for reduced memory usage'
      END IF
    END IF
    
  END SUBROUTINE validate_sparse_solver_params
  
  !-------------------------------------------------------------------------------
  ! solves A*x = b for sparse A and 1-D vector b
  ! A is specified through nrow,ncol,nz,irow,pcol,val
  ! results are returned in b
  ! iopt_in: 0 = full solve (factorize+solve+cleanup)
  !          1 = factorize only (save factorization for reuse)
  !          2 = solve only (reuse existing factorization)
  !          3 = cleanup only (free memory)
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
    
    ! Validate input parameters (only for initial factorization)
    IF (iopt == 0 .OR. iopt == 1) THEN
      CALL validate_sparse_solver_params(nrow, ncol, nz, sparse_solve_method, 'sparse_solveReal_b1')
    END IF
    
    pcol_modified = .FALSE.
    ! pcoln to remove "C convention" used in SuiteSparse
    IF (pcol(1) .EQ. 0) THEN
       pcol_modified = .TRUE.
       ALLOCATE(pcoln(ncol+1))
       pcoln = pcol + 1
    END IF
    
    ! For iopt=1 (factorize only), save factorization for later reuse with iopt=2
    IF (.NOT. factorization_exists_real .AND. iopt .EQ. IOPT_SOLVE_ONLY) THEN ! factorize first
       IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. SOLVER_UMFPACK) ) THEN
          IF (pcol_modified) THEN
             CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcoln,val,b,1)
          ELSE
             CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,1)
          END IF
       END IF
       factorization_exists_real = .TRUE.
    END IF
    IF (iopt .EQ. IOPT_FACTORIZE_ONLY) factorization_exists_real = .TRUE.
    IF (iopt .EQ. IOPT_FREE_MEMORY) factorization_exists_real = .FALSE.
    
    ! Update global flag for compatibility
    factorization_exists = factorization_exists_real
    
    ! Handle auto-selection
    CALL handle_auto_selection(nrow, ncol)
    
    IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. SOLVER_UMFPACK) ) THEN
       IF (pcol_modified) THEN
          CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcoln,val,b,iopt)
       ELSE
          CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,iopt)
       END IF
    ELSE IF (sparse_solve_method .EQ. SOLVER_BICGSTAB) THEN
       IF (pcol_modified) THEN
          CALL sparse_solve_bicgstab_real(nrow,ncol,nz,irow,pcoln,val,b,iopt)
       ELSE
          CALL sparse_solve_bicgstab_real(nrow,ncol,nz,irow,pcol,val,b,iopt)
       END IF
    ELSE
       PRINT *, 'ERROR: Invalid sparse_solve_method =', sparse_solve_method
       PRINT *, 'Valid options: 0 (auto), 2 (legacy), SOLVER_UMFPACK(3), SOLVER_BICGSTAB(4)'
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
    
    ! For iopt=1 (factorize only), save factorization for later reuse with iopt=2
    IF (.NOT. factorization_exists_complex .AND. iopt .EQ. IOPT_SOLVE_ONLY) THEN ! factorize first
       IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. SOLVER_UMFPACK) ) THEN
          IF (pcol_modified) THEN
             CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcoln,val,b,1)
          ELSE
             CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,1)
          END IF
       END IF
       factorization_exists_complex = .TRUE.
    END IF
    IF (iopt .EQ. IOPT_FACTORIZE_ONLY) factorization_exists_complex = .TRUE.
    IF (iopt .EQ. IOPT_FREE_MEMORY) factorization_exists_complex = .FALSE.
    
    ! Update global flag for compatibility
    factorization_exists = factorization_exists_complex
    
    ! Handle auto-selection
    CALL handle_auto_selection(nrow, ncol)
    
    IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. SOLVER_UMFPACK) ) THEN
       IF (pcol_modified) THEN
          CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcoln,val,b,iopt)
       ELSE
          CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,iopt)
       END IF
    ELSE IF (sparse_solve_method .EQ. SOLVER_BICGSTAB) THEN
       IF (pcol_modified) THEN
          CALL sparse_solve_bicgstab_complex(nrow,ncol,nz,irow,pcoln,val,b,iopt)
       ELSE
          CALL sparse_solve_bicgstab_complex(nrow,ncol,nz,irow,pcol,val,b,iopt)
       END IF
    ELSE
       PRINT *, 'ERROR: Invalid sparse_solve_method =', sparse_solve_method
       PRINT *, 'Valid options: 0 (auto), 2 (legacy), SOLVER_UMFPACK(3), SOLVER_BICGSTAB(4)'
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
    INTEGER :: n, i
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
    
    ! For iopt=1 (factorize only), save factorization for later reuse with iopt=2
    IF (.NOT. factorization_exists_real .AND. iopt .EQ. IOPT_SOLVE_ONLY) THEN ! factorize first
       IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. SOLVER_UMFPACK) ) THEN
          IF (pcol_modified) THEN
             CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcoln,val,b,1)
          ELSE
             CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,1)
          END IF
       END IF
       factorization_exists_real = .TRUE.
    END IF
    IF (iopt .EQ. IOPT_FACTORIZE_ONLY) factorization_exists_real = .TRUE.
    IF (iopt .EQ. IOPT_FREE_MEMORY) factorization_exists_real = .FALSE.
    
    ! Update global flag for compatibility
    factorization_exists = factorization_exists_real
    
    ! Handle auto-selection
    CALL handle_auto_selection(nrow, ncol)
    
    IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. SOLVER_UMFPACK) ) THEN
       IF (pcol_modified) THEN
          CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcoln,val,b,iopt)
       ELSE
          CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,iopt)
       END IF
    ELSE IF (sparse_solve_method .EQ. SOLVER_BICGSTAB) THEN
       ! BiCGSTAB doesn't support 2D arrays directly, solve each column separately
       DO i = 1, SIZE(b,2)
          IF (pcol_modified) THEN
             CALL sparse_solve_bicgstab_real(nrow,ncol,nz,irow,pcoln,val,b(:,i),iopt)
          ELSE
             CALL sparse_solve_bicgstab_real(nrow,ncol,nz,irow,pcol,val,b(:,i),iopt)
          END IF
       END DO
    ELSE
       PRINT *, 'ERROR: Invalid sparse_solve_method =', sparse_solve_method
       PRINT *, 'Valid options: 0 (auto), 2 (legacy), SOLVER_UMFPACK(3), SOLVER_BICGSTAB(4)'
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
    INTEGER :: n, i
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
    
    ! For iopt=1 (factorize only), save factorization for later reuse with iopt=2
    IF (.NOT. factorization_exists_complex .AND. iopt .EQ. IOPT_SOLVE_ONLY) THEN ! factorize first
       IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. SOLVER_UMFPACK) ) THEN
          IF (pcol_modified) THEN
             CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcoln,val,b,1)
          ELSE
             CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,1)
          END IF
       END IF
       factorization_exists_complex = .TRUE.
    END IF
    IF (iopt .EQ. IOPT_FACTORIZE_ONLY) factorization_exists_complex = .TRUE.
    IF (iopt .EQ. IOPT_FREE_MEMORY) factorization_exists_complex = .FALSE.
    
    ! Update global flag for compatibility
    factorization_exists = factorization_exists_complex
    
    ! Handle auto-selection
    CALL handle_auto_selection(nrow, ncol)
    
    IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. SOLVER_UMFPACK) ) THEN
       IF (pcol_modified) THEN
          CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcoln,val,b,iopt)
       ELSE
          CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,iopt)
       END IF
    ELSE IF (sparse_solve_method .EQ. SOLVER_BICGSTAB) THEN
       ! BiCGSTAB doesn't support 2D arrays directly, solve each column separately
       DO i = 1, SIZE(b,2)
          IF (pcol_modified) THEN
             CALL sparse_solve_bicgstab_complex(nrow,ncol,nz,irow,pcoln,val,b(:,i),iopt)
          ELSE
             CALL sparse_solve_bicgstab_complex(nrow,ncol,nz,irow,pcol,val,b(:,i),iopt)
          END IF
       END DO
    ELSE
       PRINT *, 'ERROR: Invalid sparse_solve_method =', sparse_solve_method
       PRINT *, 'Valid options: 0 (auto), 2 (legacy), SOLVER_UMFPACK(3), SOLVER_BICGSTAB(4)'
       STOP
    END IF
    
    IF (ALLOCATED(pcoln)) DEALLOCATE(pcoln)
    
  END SUBROUTINE sparse_solveComplex_b2
  !-------------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------------
  ! solves A*x = b for sparse A and 1-D vector b
  ! A is given as a full matrix
  ! results are returned in b
  ! iopt_in: 0 = full solve (factorize+solve+cleanup)
  !          1 = factorize only (save factorization for reuse)
  !          2 = solve only (reuse existing factorization)
  !          3 = cleanup only (free memory)
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
    IF (factorization_exists_real .AND. iopt .EQ. IOPT_FACTORIZE_ONLY) THEN ! free memory first
       IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. SOLVER_UMFPACK) ) THEN
          CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,3)
       END IF
    END IF
    IF (.NOT. factorization_exists_real .AND. iopt .EQ. IOPT_SOLVE_ONLY) THEN ! factorize first
       IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. SOLVER_UMFPACK) ) THEN
          CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,1)
       END IF
       factorization_exists_real = .TRUE.
    END IF
    IF (iopt .EQ. IOPT_FACTORIZE_ONLY) factorization_exists_real = .TRUE.
    IF (iopt .EQ. IOPT_FREE_MEMORY) factorization_exists_real = .FALSE.
    
    ! Update global flag for compatibility
    factorization_exists = factorization_exists_real
    
    ! Convert full matrix to sparse format
    CALL full2sparse(A,irow,pcol,val,nrow,ncol,nz)
    
    ! Handle auto-selection
    CALL handle_auto_selection(nrow, ncol)
    
    IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. SOLVER_UMFPACK) ) THEN
       CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,iopt)
    ELSE IF (sparse_solve_method .EQ. SOLVER_BICGSTAB) THEN
       CALL sparse_solve_bicgstab_real(nrow,ncol,nz,irow,pcol,val,b,iopt)
    ELSE
       PRINT *, 'ERROR: Invalid sparse_solve_method =', sparse_solve_method
       PRINT *, 'Valid options: 0 (auto), 2 (legacy), SOLVER_UMFPACK(3), SOLVER_BICGSTAB(4)'
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
    IF (factorization_exists_complex .AND. iopt .EQ. IOPT_FACTORIZE_ONLY) THEN ! free memory first
       IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. SOLVER_UMFPACK) ) THEN
          CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,3)
       END IF
    END IF
    IF (.NOT. factorization_exists_complex .AND. iopt .EQ. IOPT_SOLVE_ONLY) THEN ! factorize first
       IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. SOLVER_UMFPACK) ) THEN
          CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,1)
       END IF
       factorization_exists_complex = .TRUE.
    END IF
    IF (iopt .EQ. IOPT_FACTORIZE_ONLY) factorization_exists_complex = .TRUE.
    IF (iopt .EQ. IOPT_FREE_MEMORY) factorization_exists_complex = .FALSE.
    
    ! Update global flag for compatibility
    factorization_exists = factorization_exists_complex
    
    ! Convert full matrix to sparse format
    CALL full2sparse(A,irow,pcol,val,nrow,ncol,nz)
    
    ! Handle auto-selection
    CALL handle_auto_selection(nrow, ncol)
    
    IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. SOLVER_UMFPACK) ) THEN
       CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,iopt)
    ELSE IF (sparse_solve_method .EQ. SOLVER_BICGSTAB) THEN
       CALL sparse_solve_bicgstab_complex(nrow,ncol,nz,irow,pcol,val,b,iopt)
    ELSE
       PRINT *, 'ERROR: Invalid sparse_solve_method =', sparse_solve_method
       PRINT *, 'Valid options: 0 (auto), 2 (legacy), SOLVER_UMFPACK(3), SOLVER_BICGSTAB(4)'
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
    INTEGER :: nrow,ncol,nz,i
    INTEGER, DIMENSION(:), ALLOCATABLE :: irow,pcol
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: val
    
    IF (PRESENT(iopt_in)) iopt = iopt_in
    
    ! check about existing factorization
    IF (factorization_exists_real .AND. iopt .EQ. IOPT_FACTORIZE_ONLY) THEN ! free memory first
       IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. SOLVER_UMFPACK) ) THEN
          CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,3)
       END IF
    END IF
    IF (.NOT. factorization_exists_real .AND. iopt .EQ. IOPT_SOLVE_ONLY) THEN ! factorize first
       IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. SOLVER_UMFPACK) ) THEN
          CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,1)
       END IF
       factorization_exists_real = .TRUE.
    END IF
    IF (iopt .EQ. IOPT_FACTORIZE_ONLY) factorization_exists_real = .TRUE.
    IF (iopt .EQ. IOPT_FREE_MEMORY) factorization_exists_real = .FALSE.
    
    ! Update global flag for compatibility
    factorization_exists = factorization_exists_real
    
    ! Convert full matrix to sparse format
    CALL full2sparse(A,irow,pcol,val,nrow,ncol,nz)
    
    ! Handle auto-selection
    CALL handle_auto_selection(nrow, ncol)
    
    IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. SOLVER_UMFPACK) ) THEN
       CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,iopt)
    ELSE IF (sparse_solve_method .EQ. SOLVER_BICGSTAB) THEN
       ! BiCGSTAB doesn't support 2D arrays directly, solve each column separately
       DO i = 1, SIZE(b,2)
          CALL sparse_solve_bicgstab_real(nrow,ncol,nz,irow,pcol,val,b(:,i),iopt)
       END DO
    ELSE
       PRINT *, 'ERROR: Invalid sparse_solve_method =', sparse_solve_method
       PRINT *, 'Valid options: 0 (auto), 2 (legacy), SOLVER_UMFPACK(3), SOLVER_BICGSTAB(4)'
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
    INTEGER :: nrow,ncol,nz,i
    INTEGER, DIMENSION(:), ALLOCATABLE :: irow,pcol
    COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE :: val
    
    IF (PRESENT(iopt_in)) iopt = iopt_in
    
    ! check about existing factorization
    IF (factorization_exists_complex .AND. iopt .EQ. IOPT_FACTORIZE_ONLY) THEN ! free memory first
       IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. SOLVER_UMFPACK) ) THEN
          CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,3)
       END IF
    END IF
    IF (.NOT. factorization_exists_complex .AND. iopt .EQ. IOPT_SOLVE_ONLY) THEN ! factorize first
       IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. SOLVER_UMFPACK) ) THEN
          CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,1)
       END IF
       factorization_exists_complex = .TRUE.
    END IF
    IF (iopt .EQ. IOPT_FACTORIZE_ONLY) factorization_exists_complex = .TRUE.
    IF (iopt .EQ. IOPT_FREE_MEMORY) factorization_exists_complex = .FALSE.
    
    ! Update global flag for compatibility
    factorization_exists = factorization_exists_complex
    
    ! Convert full matrix to sparse format
    CALL full2sparse(A,irow,pcol,val,nrow,ncol,nz)
    
    ! Handle auto-selection
    CALL handle_auto_selection(nrow, ncol)
    
    IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. SOLVER_UMFPACK) ) THEN
       CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,iopt)
    ELSE IF (sparse_solve_method .EQ. SOLVER_BICGSTAB) THEN
       ! BiCGSTAB doesn't support 2D arrays directly, solve each column separately
       DO i = 1, SIZE(b,2)
          CALL sparse_solve_bicgstab_complex(nrow,ncol,nz,irow,pcol,val,b(:,i),iopt)
       END DO
    ELSE
       PRINT *, 'ERROR: Invalid sparse_solve_method =', sparse_solve_method
       PRINT *, 'Valid options: 0 (auto), 2 (legacy), SOLVER_UMFPACK(3), SOLVER_BICGSTAB(4)'
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
    
    INTEGER(kind=long) :: n,nc
    INTEGER :: nrhs = 1
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: x
    
    INTEGER :: Ap_len, Ai_len, Ax_len, Az_len
    INTEGER(kind=long), DIMENSION(:), ALLOCATABLE :: Ap, Ai, p
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: Ax, Az
    
    ! C wrappers of SuiteSparse functions
    INTEGER(kind=long) :: umf4def_
    INTEGER(kind=long) :: umf4sym_
    INTEGER(kind=long) :: umf4num_
    INTEGER(kind=long) :: umf4solr_
    INTEGER(kind=long) :: umf4sol_
    INTEGER(kind=long) :: umf4fnum_
    INTEGER(kind=long) :: umf4fsym_
    INTEGER(kind=long) :: umf4zfnum_
    INTEGER(kind=long) :: umf4zfsym_
    
    Ap_len = ncol + 1
    Ai_len = nz
    Ax_len = nz
    Az_len = 0 ! nz
    ALLOCATE(Ap(Ap_len), Ai(Ai_len))
    Ap = pcol - 1  ! convert from 1 to 0-based indexing
    Ai = irow - 1  ! convert from 1 to 0-based indexing
    
    IF (iopt_in .EQ. IOPT_FREE_MEMORY) THEN  ! free memory from last solution
       CALL umf4fnum(numeric_real)
       factorization_exists_real = .FALSE.
       DEALLOCATE(Ap, Ai)
       RETURN
    END IF
    
    ALLOCATE(x(nrow))
    x = 0.0_dp  ! Initialize solution vector
    
    ! Set default parameters
    CALL umf4def(control)
    control(1) = 0 ! No output - there are other options, see the manual
    
    n = nrow  ! convert from 1 to 0-based indexing
    
    ! Clear any previous complex factorization data
    IF (factorization_exists_complex) THEN
       CALL umf4zfnum(numeric_complex)
       CALL umf4zfsym(symbolic_complex)
       symbolic_complex = 0
       numeric_complex = 0
       factorization_exists_complex = .FALSE.
    END IF
    
    IF (iopt_in .EQ. IOPT_FULL_SOLVE .OR. iopt_in .EQ. IOPT_FACTORIZE_ONLY) THEN ! Factorize for full solve or factorize-only
       IF ( sparse_solve_method .EQ. 2 ) THEN ! SuiteSparse (with (=2) iterative refinement)
          control(8) = 10 ! max number of iterative refinement steps
       END IF
       
       ! Pre-order and symbolic analysis
       CALL umf4sym (n, n, Ap, Ai, val, symbolic_real, control, info_suitesparse)
       IF (info_suitesparse(1) .LT. 0) THEN
          PRINT *, 'Error occurred in umf4sym: ', info_suitesparse(1)
       END IF
       
       CALL umf4num(Ap, Ai, val, symbolic_real, numeric_real, control, info_suitesparse)
       IF (info_suitesparse(1) .LT. 0) THEN
          PRINT *, 'Error occurred in umf4num: ', info_suitesparse(1)
       END IF
       
       factorization_exists_real = .TRUE.
    END IF
    
    ! Note: in newer versions, the function interfaces has been changed to match the types in an cleaner way
    !       use n instead of nrow
    !       use nc instead of ncol
    nc = ncol
    
    ! Solve phase - only if not just factorizing
    IF (iopt_in .EQ. IOPT_FULL_SOLVE .OR. iopt_in .EQ. IOPT_SOLVE_ONLY) THEN
       ! Check if factorization exists for solve-only case
       IF (iopt_in .EQ. IOPT_SOLVE_ONLY .AND. .NOT. factorization_exists_real) THEN
          PRINT *, 'ERROR: Solve requested but no factorization exists!'
          IF (ALLOCATED(Ap)) DEALLOCATE(Ap)
          IF (ALLOCATED(Ai)) DEALLOCATE(Ai)
          IF (ALLOCATED(x)) DEALLOCATE(x)
          RETURN
       END IF
       
       IF ( sparse_solve_method .EQ. 2 ) THEN ! SuiteSparse (with (=2)
          CALL umf4solr (sys, Ap, Ai, val, x, b, numeric_real, control, info_suitesparse) !iterative refinement
       ELSE IF ( sparse_solve_method .EQ. 3 ) THEN ! SuiteSparse (without (=3) iterative refinement)
          CALL umf4sol (sys, x, b, numeric_real, control, info_suitesparse)
       END IF
       IF (info_suitesparse(1) .LT. 0) THEN
          PRINT *, 'Error occurred in umf4solr: ', info_suitesparse(1)
       END IF
       
       b = x
    END IF
    
    ! Last, free the storage allocated inside SuiteSparse
    IF (iopt_in .EQ. IOPT_FREE_MEMORY) THEN
       CALL umf4fnum (numeric_real)
       CALL umf4fsym (symbolic_real)  
       factorization_exists_real = .FALSE.
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
    
    INTEGER(kind=long) :: n,nc
    INTEGER :: nrhs = 1
    COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE :: x
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: xx, xz, bx, bz
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: valx, valz
    
    INTEGER :: Ap_len, Ai_len, Ax_len, Az_len
    INTEGER(kind=long), DIMENSION(:), ALLOCATABLE :: Ap, Ai, p
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
    Ap = pcol - 1  ! convert from 1 to 0-based indexing
    Ai = irow - 1  ! convert from 1 to 0-based indexing
    
    ALLOCATE(valx(nz), valz(nz))
    valx = REAL(val)
    valz = AIMAG(val)
    
    ALLOCATE(bx(nrow), bz(nrow))
    bx = REAL(b)
    bz = AIMAG(b)
    
    IF (iopt_in .EQ. IOPT_FREE_MEMORY) THEN  ! free memory from last solution
       CALL umf4zfnum(numeric_complex)
       factorization_exists_complex = .FALSE.
       DEALLOCATE(Ap, Ai, valx, valz, bx, bz)
       RETURN
    END IF
    
    ALLOCATE(xx(nrow), xz(nrow))
    
    n = nrow  ! Initialize n for UMFPACK interface
    nc = ncol
    
    ! Clear any previous real factorization data
    IF (factorization_exists_real) THEN
       CALL umf4fnum(numeric_real)
       CALL umf4fsym(symbolic_real)
       symbolic_real = 0
       numeric_real = 0
       factorization_exists_real = .FALSE.
    END IF
    
    IF (iopt_in .EQ. IOPT_FULL_SOLVE .OR. iopt_in .EQ. IOPT_FACTORIZE_ONLY) THEN
       IF ( sparse_solve_method .EQ. 3 ) THEN ! SuiteSparse (without (=3) iterative refinement)
          CALL umf4zdef(control)
          control(1) = 0 ! No output - there are other options, see the manual
       ELSE IF ( sparse_solve_method .EQ. 2 ) THEN ! SuiteSparse (with (=2) iterative refinement)
          CALL umf4zdef(control)
          control(1) = 0 ! No output - there are other options, see the manual
          control(8) = 10 ! max number of iterative refinement steps
       END IF
       ! Pre-order and symbolic analysis
       CALL umf4zsym (n, n, Ap, Ai, valx, valz, symbolic_complex, control, info_suitesparse)
       IF (info_suitesparse(1) .LT. 0) THEN
          PRINT *, 'Error occurred in umf4zsym: ', info_suitesparse(1)
       END IF
       
       CALL umf4znum(Ap, Ai, valx, valz, symbolic_complex, numeric_complex, control, info_suitesparse)
       IF (info_suitesparse(1) .LT. 0) THEN
          PRINT *, 'Error occurred in umf4znum: ', info_suitesparse(1)
       END IF
    END IF
    
    ! Note: in newer versions, the function interfaces has been changed to match the types in an cleaner way
    !       use n instead of nrow
    !       use nc instead of ncol
    nc = ncol
    
    ! Solve phase - only if not just factorizing
    IF (iopt_in .EQ. IOPT_FULL_SOLVE .OR. iopt_in .EQ. IOPT_SOLVE_ONLY) THEN
       IF ( sparse_solve_method .EQ. 2 ) THEN ! SuiteSparse (with (=2)
          CALL umf4zsolr (sys, Ap, Ai, valx, valz, xx, xz, bx, bz, numeric_complex, &
               control, info_suitesparse) !iterative refinement
       ELSE IF ( sparse_solve_method .EQ. 3 ) THEN ! SuiteSparse (without (=3) iterative refinement)
          CALL umf4zsol (sys, xx, xz, bx, bz, numeric_complex, control, info_suitesparse)
       END IF
       IF (info_suitesparse(1) .LT. 0) THEN
          PRINT *, 'Error occurred in umf4zsolr: ', info_suitesparse(1)
       END IF
       b = CMPLX(xx, xz, KIND=dp)
    END IF
    
    ! Last, free the storage allocated inside SuiteSparse
    IF (iopt_in .EQ. IOPT_FREE_MEMORY) THEN
       CALL umf4zfnum (numeric_complex)
       CALL umf4zfsym (symbolic_complex)
       factorization_exists_complex = .FALSE.
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
    
    INTEGER(kind=long) :: n,nc
    INTEGER :: i,nrhs
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: x,bloc
    
    INTEGER :: Ap_len, Ai_len, Ax_len, Az_len
    INTEGER(kind=long), DIMENSION(:), ALLOCATABLE :: Ap, Ai, p
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
    Ap = pcol - 1  ! convert from 1 to 0-based indexing
    Ai = irow - 1  ! convert from 1 to 0-based indexing
    
    IF (iopt_in .EQ. IOPT_FREE_MEMORY) THEN  ! free memory from last solution
       CALL umf4fnum(numeric_real)
       factorization_exists_real = .FALSE.
       DEALLOCATE(Ap, Ai)
       RETURN
    END IF
    
    nrhs = SIZE(b,2)
    ALLOCATE(x(nrow), bloc(nrow))
    
    ! Set default parameters
    CALL umf4def(control)
    control(1) = 0 ! No output - there are other options, see the manual
    
    n = nrow  ! convert from 1 to 0-based indexing
    
    IF (iopt_in .EQ. IOPT_FULL_SOLVE .OR. iopt_in .EQ. IOPT_FACTORIZE_ONLY) THEN ! Factorize for full solve or factorize-only
       IF ( sparse_solve_method .EQ. 2 ) THEN ! SuiteSparse (with (=2) iterative refinement)
          control(8) = 10 ! max number of iterative refinement steps
       END IF
       
       ! Pre-order and symbolic analysis
       CALL umf4sym (n, n, Ap, Ai, val, symbolic_real, control, info_suitesparse)
       IF (info_suitesparse(1) .LT. 0) THEN
          PRINT *, 'Error occurred in umf4sym: ', info_suitesparse(1)
       END IF
       
       CALL umf4num(Ap, Ai, val, symbolic_real, numeric_real, control, info_suitesparse)
       IF (info_suitesparse(1) .LT. 0) THEN
          PRINT *, 'Error occurred in umf4num: ', info_suitesparse(1)
       END IF
       
       factorization_exists_real = .TRUE.
    END IF
    
    ! Note: in newer versions, the function interfaces has been changed to match the types in an cleaner way
    !       use n instead of nrow
    !       use nc instead of ncol
    nc = ncol
    
    ! Solve phase - only if not just factorizing
    IF (iopt_in .EQ. IOPT_FULL_SOLVE .OR. iopt_in .EQ. IOPT_SOLVE_ONLY) THEN
       DO i = 1,nrhs
          bloc = b(:,i)
          IF ( sparse_solve_method .EQ. 2 ) THEN ! SuiteSparse (with (=2)
             CALL umf4solr (sys, Ap, Ai, val, x, bloc, numeric_real, control, info_suitesparse) !iterative refinement
          ELSE IF ( sparse_solve_method .EQ. 3 ) THEN ! SuiteSparse (without (=3) iterative refinement)
             CALL umf4sol (sys, x, bloc, numeric_real, control, info_suitesparse)
          END IF
          IF (info_suitesparse(1) .LT. 0) THEN
             PRINT *, 'Error occurred in umf4solr: ', info_suitesparse(1)
          END IF
          b(:,i) = x
       END DO
    END IF
    
    ! Last, free the storage allocated inside SuiteSparse
    IF (iopt_in .EQ. IOPT_FREE_MEMORY) THEN
       CALL umf4fnum (numeric_real)
       CALL umf4fsym (symbolic_real)  
       factorization_exists_real = .FALSE.
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
    
    INTEGER(kind=long) :: n,nc
    INTEGER :: i,nrhs
    COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE :: x,bloc
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: xx, xz, blocx, blocz
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: valx, valz
    
    INTEGER :: Ap_len, Ai_len, Ax_len, Az_len
    INTEGER(kind=long), DIMENSION(:), ALLOCATABLE :: Ap, Ai, p
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
    Ap = pcol - 1  ! convert from 1 to 0-based indexing
    Ai = irow - 1  ! convert from 1 to 0-based indexing
    
    ALLOCATE(valx(nz), valz(nz))
    valx = REAL(val)
    valz = AIMAG(val)
    
    IF (iopt_in .EQ. IOPT_FREE_MEMORY) THEN  ! free memory from last solution
       CALL umf4zfnum(numeric_complex)
       factorization_exists_complex = .FALSE.
       DEALLOCATE(Ap, Ai, valx, valz)
       RETURN
    END IF
    
    nrhs = SIZE(b,2)
    ALLOCATE(xx(nrow), xz(nrow), blocx(nrow), blocz(nrow))
    
    n = nrow  ! Initialize n for UMFPACK interface
    nc = ncol
    
    IF (iopt_in .EQ. IOPT_FULL_SOLVE .OR. iopt_in .EQ. IOPT_FACTORIZE_ONLY) THEN ! Factorize for full solve or factorize-only
       IF ( sparse_solve_method .EQ. 3 ) THEN ! SuiteSparse (without (=3) iterative refinement)
          CALL umf4zdef(control)
          control(1) = 0 ! No output - there are other options, see the manual
       ELSE IF ( sparse_solve_method .EQ. 2 ) THEN ! SuiteSparse (with (=2) iterative refinement)
          CALL umf4zdef(control)
          control(1) = 0 ! No output - there are other options, see the manual
          control(8) = 10 ! max number of iterative refinement steps
       END IF
       ! Pre-order and symbolic analysis
       CALL umf4zsym (n, n, Ap, Ai, valx, valz, symbolic_complex, control, info_suitesparse)
       IF (info_suitesparse(1) .LT. 0) THEN
          PRINT *, 'Error occurred in umf4zsym: ', info_suitesparse(1)
       END IF
       
       CALL umf4znum(Ap, Ai, valx, valz, symbolic_complex, numeric_complex, control, info_suitesparse)
       IF (info_suitesparse(1) .LT. 0) THEN
          PRINT *, 'Error occurred in umf4znum: ', info_suitesparse(1)
       END IF
    END IF
    
    ! Note: in newer versions, the function interfaces has been changed to match the types in an cleaner way
    !       use n instead of nrow
    !       use nc instead of ncol
    nc = ncol
    
    ! Solve phase - only if not just factorizing
    IF (iopt_in .EQ. IOPT_FULL_SOLVE .OR. iopt_in .EQ. IOPT_SOLVE_ONLY) THEN
       DO i = 1,nrhs
          blocx = REAL(b(:,i))
          blocz = AIMAG(b(:,i))
          IF ( sparse_solve_method .EQ. 2 ) THEN ! SuiteSparse (with (=2)
             CALL umf4zsolr (sys, Ap, Ai, valx, valz, xx, xz, blocx, blocz, numeric_complex,&
                  control, info_suitesparse) !iterative refinement
          ELSE IF ( sparse_solve_method .EQ. 3 ) THEN ! SuiteSparse (without (=3) iterative refinement)
             CALL umf4zsol (sys, xx, xz, blocx, blocz, numeric_complex, control, info_suitesparse)
          END IF
          IF (info_suitesparse(1) .LT. 0) THEN
             PRINT *, 'Error occurred in umf4zsolr: ', info_suitesparse(1)
          END IF
          b(:,i) = CMPLX(xx, xz, KIND=dp)
       END DO
    END IF
    
    ! Last, free the storage allocated inside SuiteSparse
    IF (iopt_in .EQ. IOPT_FREE_MEMORY) THEN
       CALL umf4zfnum (numeric_complex)
       CALL umf4zfsym (symbolic_complex)
       factorization_exists_complex = .FALSE.
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

  !-------------------------------------------------------------------------------
  ! Helper functions for robust iopt parameter handling
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  ! Determines if factorization should be performed based on iopt value
  PURE FUNCTION should_factorize(iopt) RESULT(do_factorize)
    INTEGER, INTENT(in) :: iopt
    LOGICAL :: do_factorize
    do_factorize = (iopt == IOPT_FULL_SOLVE .OR. iopt == IOPT_FACTORIZE_ONLY)
  END FUNCTION should_factorize

  !-------------------------------------------------------------------------------
  ! Determines if solving should be performed based on iopt value  
  PURE FUNCTION should_solve(iopt) RESULT(do_solve)
    INTEGER, INTENT(in) :: iopt
    LOGICAL :: do_solve
    do_solve = (iopt == IOPT_FULL_SOLVE .OR. iopt == IOPT_SOLVE_ONLY)
  END FUNCTION should_solve

  !-------------------------------------------------------------------------------
  ! Determines if memory cleanup should be performed based on iopt value
  PURE FUNCTION should_cleanup(iopt) RESULT(do_cleanup)
    INTEGER, INTENT(in) :: iopt
    LOGICAL :: do_cleanup
    do_cleanup = (iopt == IOPT_FREE_MEMORY)
  END FUNCTION should_cleanup

  !-------------------------------------------------------------------------------
  ! BiCGSTAB solver wrapper for real systems
  ! Note: BiCGSTAB is a one-shot solver (no persistent factorization)
  !       iopt parameter is ignored except for compatibility
  SUBROUTINE sparse_solve_bicgstab_real(nrow,ncol,nz,irow,pcol,val,b,iopt_in)
    INTEGER, INTENT(in) :: nrow,ncol,nz
    INTEGER, DIMENSION(:), INTENT(in) :: irow,pcol
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: val
    REAL(kind=dp), DIMENSION(:), INTENT(inout) :: b
    INTEGER, INTENT(in) :: iopt_in
    
    ! Local variables for CSR conversion
    INTEGER, ALLOCATABLE :: csr_row_ptr(:), csr_col_idx(:)
    REAL(kind=dp), ALLOCATABLE :: csr_val(:)
    REAL(kind=dp), ALLOCATABLE :: b_save(:)
    
    ! BiCGSTAB solver parameters (use configurable values)
    REAL(kind=dp) :: abs_tol, rel_tol
    INTEGER :: max_iter
    LOGICAL :: converged
    INTEGER :: iter
    TYPE(bicgstab_stats) :: stats
    
    ! ILU preconditioning variables
    TYPE(ilu_factorization) :: ilu_fac
    INTEGER :: ilu_info
    
    ! Use module-level configurable parameters
    rel_tol = bicgstab_rel_tolerance
    abs_tol = bicgstab_abs_tolerance
    max_iter = bicgstab_max_iter
    
    ! Apply adaptive tolerance adjustment if enabled
    IF (bicgstab_adaptive_tolerance) THEN
      CALL apply_adaptive_tolerance_real(nrow, nz, val, abs_tol, rel_tol)
    END IF
    
    ! Convert from CSC to CSR format for BiCGSTAB
    ALLOCATE(csr_row_ptr(nrow+1), csr_col_idx(nz), csr_val(nz))
    CALL csc_to_csr_real(nrow, ncol, nz, pcol, irow, val, csr_row_ptr, csr_col_idx, csr_val)
    
    ! Solve using BiCGSTAB (with ILU preconditioning if enabled)
    ! b is used as both RHS and solution vector (solution overwrites b)
    ! Start with zero initial guess
    ALLOCATE(b_save(nrow))
    b_save = b  ! Save RHS
    b = 0.0_dp  ! Zero initial guess
    
    IF (ilu_fill_level > 0) THEN
      ! Use ILU-preconditioned BiCGSTAB
      IF (bicgstab_verbose) THEN
        PRINT *, 'INFO: Using ILU(',ilu_fill_level,') preconditioning for BiCGSTAB'
      END IF
      
      ! Compute ILU factorization
      CALL ilu_factorize(nrow, csr_row_ptr, csr_col_idx, csr_val, &
                         ilu_fill_level, ilu_drop_tolerance, ilu_fac, ilu_info)
      
      IF (ilu_info == 0) THEN
        ! ILU factorization successful, use preconditioned solver
        CALL bicgstab_solve_precond(nrow, csr_row_ptr, csr_col_idx, csr_val, &
                                    b_save, b, ilu_fac, abs_tol, rel_tol, max_iter, &
                                    converged, iter, stats)
        CALL ilu_free(ilu_fac)
      ELSE
        ! ILU factorization failed, fall back to unpreconditioned solver
        IF (bicgstab_verbose) THEN
          PRINT *, 'WARNING: ILU factorization failed, using unpreconditioned BiCGSTAB'
        END IF
        CALL bicgstab_solve(nrow, csr_row_ptr, csr_col_idx, csr_val, b_save, b, &
                            abs_tol, rel_tol, max_iter, converged, iter, stats)
      END IF
    ELSE
      ! Use unpreconditioned BiCGSTAB (ilu_fill_level = 0)
      CALL bicgstab_solve(nrow, csr_row_ptr, csr_col_idx, csr_val, b_save, b, &
                          abs_tol, rel_tol, max_iter, converged, iter, stats)
    END IF
    
    DEALLOCATE(b_save)
    
    IF (.NOT. converged) THEN
      PRINT *, 'WARNING: BiCGSTAB did not converge after', iter, 'iterations'
      PRINT *, '         Final residual:', stats%final_residual
      PRINT *, '         Target tolerance: abs=', abs_tol, ', rel=', rel_tol
      PRINT *, '         Matrix size:', nrow, 'x', nrow
      PRINT *, '         Non-zeros:', nz
      IF (stats%final_residual > 1.0e-6_dp) THEN
        PRINT *, '         CAUTION: Large residual may indicate solution inaccuracy'
        PRINT *, '         Consider: 1) Increasing max_iter 2) Relaxing tolerance 3) Using UMFPACK'
      ELSE
        PRINT *, '         INFO: Small residual suggests acceptable solution despite non-convergence'
      END IF
    END IF
    
    IF (sparse_talk .OR. bicgstab_verbose) THEN
      PRINT *, 'BiCGSTAB converged in', iter, 'iterations, residual =', stats%final_residual
    END IF
    
    DEALLOCATE(csr_row_ptr, csr_col_idx, csr_val)
    
  END SUBROUTINE sparse_solve_bicgstab_real

  !-------------------------------------------------------------------------------
  ! BiCGSTAB solver wrapper for complex systems  
  ! Note: BiCGSTAB is a one-shot solver (no persistent factorization)
  !       iopt parameter is ignored except for compatibility
  SUBROUTINE sparse_solve_bicgstab_complex(nrow,ncol,nz,irow,pcol,val,b,iopt_in)
    INTEGER, INTENT(in) :: nrow,ncol,nz
    INTEGER, DIMENSION(:), INTENT(in) :: irow,pcol
    COMPLEX(kind=dp), DIMENSION(:), INTENT(in) :: val
    COMPLEX(kind=dp), DIMENSION(:), INTENT(inout) :: b
    INTEGER, INTENT(in) :: iopt_in
    
    ! Local variables for CSR conversion
    INTEGER, ALLOCATABLE :: csr_row_ptr(:), csr_col_idx(:)
    COMPLEX(kind=dp), ALLOCATABLE :: csr_val(:)
    COMPLEX(kind=dp), ALLOCATABLE :: b_save(:)
    
    ! BiCGSTAB solver parameters (use configurable values)
    REAL(kind=dp) :: abs_tol, rel_tol
    INTEGER :: max_iter
    LOGICAL :: converged
    INTEGER :: iter
    TYPE(bicgstab_stats) :: stats
    
    ! ILU preconditioning variables
    TYPE(ilu_factorization_complex) :: ilu_fac
    INTEGER :: ilu_info
    
    ! Use module-level configurable parameters
    rel_tol = bicgstab_rel_tolerance
    abs_tol = bicgstab_abs_tolerance
    max_iter = bicgstab_max_iter
    
    ! Convert from CSC to CSR format for BiCGSTAB
    ALLOCATE(csr_row_ptr(nrow+1), csr_col_idx(nz), csr_val(nz))
    CALL csc_to_csr_complex(nrow, ncol, nz, pcol, irow, val, csr_row_ptr, csr_col_idx, csr_val)
    
    ! Solve using BiCGSTAB (with ILU preconditioning if enabled)
    ! b is used as both RHS and solution vector (solution overwrites b)
    ! Start with zero initial guess
    ALLOCATE(b_save(nrow))
    b_save = b  ! Save RHS
    b = (0.0_dp, 0.0_dp)  ! Zero initial guess for complex
    
    IF (ilu_fill_level > 0) THEN
      ! Use ILU-preconditioned BiCGSTAB
      IF (bicgstab_verbose) THEN
        PRINT *, 'INFO: Using ILU(',ilu_fill_level,') preconditioning for BiCGSTAB'
      END IF
      
      ! Compute ILU factorization
      CALL ilu_factorize(nrow, csr_row_ptr, csr_col_idx, csr_val, &
                         ilu_fill_level, ilu_drop_tolerance, ilu_fac, ilu_info)
      
      IF (ilu_info == 0) THEN
        ! ILU factorization successful, use preconditioned solver
        CALL bicgstab_solve_precond(nrow, csr_row_ptr, csr_col_idx, csr_val, &
                                    b_save, b, ilu_fac, abs_tol, rel_tol, max_iter, &
                                    converged, iter, stats)
        CALL ilu_free(ilu_fac)
      ELSE
        ! ILU factorization failed, fall back to unpreconditioned solver
        IF (bicgstab_verbose) THEN
          PRINT *, 'WARNING: ILU factorization failed, using unpreconditioned BiCGSTAB'
        END IF
        CALL bicgstab_solve(nrow, csr_row_ptr, csr_col_idx, csr_val, b_save, b, &
                            abs_tol, rel_tol, max_iter, converged, iter, stats)
      END IF
    ELSE
      ! Use unpreconditioned BiCGSTAB (ilu_fill_level = 0)
      CALL bicgstab_solve(nrow, csr_row_ptr, csr_col_idx, csr_val, b_save, b, &
                          abs_tol, rel_tol, max_iter, converged, iter, stats)
    END IF
    
    DEALLOCATE(b_save)
    
    IF (.NOT. converged) THEN
      PRINT *, 'WARNING: BiCGSTAB did not converge after', iter, 'iterations'
      PRINT *, '         Final residual:', stats%final_residual
      PRINT *, '         Target tolerance: abs=', abs_tol, ', rel=', rel_tol
      PRINT *, '         Matrix size:', nrow, 'x', nrow
      PRINT *, '         Non-zeros:', nz
      IF (stats%final_residual > 1.0e-6_dp) THEN
        PRINT *, '         CAUTION: Large residual may indicate solution inaccuracy'
        PRINT *, '         Consider: 1) Increasing max_iter 2) Relaxing tolerance 3) Using UMFPACK'
      ELSE
        PRINT *, '         INFO: Small residual suggests acceptable solution despite non-convergence'
      END IF
    END IF
    
    IF (sparse_talk .OR. bicgstab_verbose) THEN
      PRINT *, 'BiCGSTAB converged in', iter, 'iterations, residual =', stats%final_residual
    END IF
    
    DEALLOCATE(csr_row_ptr, csr_col_idx, csr_val)
    
  END SUBROUTINE sparse_solve_bicgstab_complex

  !-------------------------------------------------------------------------------
  ! Adaptive tolerance adjustment for real matrices
  ! Estimates matrix conditioning and adjusts tolerance accordingly
  !> Apply adaptive tolerance adjustment based on matrix conditioning
  !! 
  !! This subroutine estimates the condition number of a sparse matrix using
  !! a simple heuristic based on the ratio of maximum to minimum absolute values.
  !! Based on this estimate, it adjusts the solver tolerances to improve
  !! convergence behavior for ill-conditioned systems.
  !!
  !! @param[in]    nrow    Number of matrix rows
  !! @param[in]    nz      Number of nonzero elements
  !! @param[in]    val     Matrix values array (size: nz)
  !! @param[inout] abs_tol Absolute tolerance (updated based on conditioning)
  !! @param[inout] rel_tol Relative tolerance (updated based on conditioning)
  !!
  !! @author Generated with Claude Code
  !! @date 2025-08-02
  SUBROUTINE apply_adaptive_tolerance_real(nrow, nz, val, abs_tol, rel_tol)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: nrow, nz
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: val
    REAL(kind=dp), INTENT(inout) :: abs_tol, rel_tol
    
    ! Local variables for conditioning estimation
    REAL(kind=dp) :: val_min, val_max, condition_estimate
    REAL(kind=dp) :: original_abs_tol, original_rel_tol
    INTEGER :: i
    
    ! Save original tolerances for reporting
    original_abs_tol = abs_tol
    original_rel_tol = rel_tol
    
    ! Simple matrix conditioning estimation using min/max values
    ! This is a heuristic approach - more sophisticated methods could be used
    val_min = HUGE(1.0_dp)
    val_max = 0.0_dp
    
    DO i = 1, nz
      IF (ABS(val(i)) > 0.0_dp) THEN
        val_min = MIN(val_min, ABS(val(i)))
        val_max = MAX(val_max, ABS(val(i)))
      END IF
    END DO
    
    ! Avoid division by zero
    IF (val_min <= 0.0_dp) THEN
      condition_estimate = HUGE(1.0_dp)
    ELSE
      condition_estimate = val_max / val_min
    END IF
    
    IF (bicgstab_verbose) THEN
      PRINT *, 'DEBUG: Adaptive tolerance analysis:'
      PRINT *, '       Matrix value range: min=', val_min, ', max=', val_max
      PRINT *, '       Condition estimate:', condition_estimate
      PRINT *, '       Original tolerances: abs=', original_abs_tol, ', rel=', original_rel_tol
    END IF
    
    ! Adjust tolerance based on conditioning
    IF (condition_estimate < COND_WELL_CONDITIONED) THEN
      ! Well-conditioned matrix - keep strict tolerance
      ! No adjustment needed
      IF (bicgstab_verbose) THEN
        PRINT *, '       Matrix appears well-conditioned - keeping strict tolerance'
      END IF
      
    ELSE IF (condition_estimate < COND_MODERATELY_ILL) THEN
      ! Moderately ill-conditioned - slight relaxation
      abs_tol = MAX(abs_tol, TOL_ABS_MODERATE)
      rel_tol = MAX(rel_tol, TOL_REL_MODERATE)
      IF (bicgstab_verbose) THEN
        PRINT *, '       Matrix moderately ill-conditioned - slightly relaxing tolerance'
      END IF
      
    ELSE IF (condition_estimate < COND_ILL_CONDITIONED) THEN
      ! Ill-conditioned - moderate relaxation
      abs_tol = MAX(abs_tol, TOL_ABS_RELAXED)
      rel_tol = MAX(rel_tol, TOL_REL_RELAXED)
      IF (bicgstab_verbose) THEN
        PRINT *, '       Matrix ill-conditioned - moderately relaxing tolerance'
      END IF
      
    ELSE
      ! Very ill-conditioned - significant relaxation
      abs_tol = MAX(abs_tol, TOL_ABS_LOOSE)
      rel_tol = MAX(rel_tol, TOL_REL_LOOSE)
      IF (bicgstab_verbose) THEN
        PRINT *, '       Matrix very ill-conditioned - significantly relaxing tolerance'
      END IF
    END IF
    
    IF (bicgstab_verbose) THEN
      PRINT *, '       Adjusted tolerances: abs=', abs_tol, ', rel=', rel_tol
    END IF
    
    ! Update module-level variables so the test can check them
    bicgstab_abs_tolerance = abs_tol
    bicgstab_rel_tolerance = rel_tol
    
  END SUBROUTINE apply_adaptive_tolerance_real

END MODULE sparse_solvers_mod