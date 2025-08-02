PROGRAM test_ilu_fill_levels
  ! Test different ILU fill levels on spline matrices until BiCGSTAB converges
  ! This explores the effectiveness of ILU(k) preconditioning for k = 0, 1, 2, 3, ...
  
  USE nrtype, ONLY: I4B, DP
  USE sparse_solvers_mod, ONLY: sparse_solve_method, SOLVER_BICGSTAB, &
                                ilu_fill_level, bicgstab_max_iter, &
                                bicgstab_abs_tolerance, bicgstab_rel_tolerance, &
                                bicgstab_verbose
  USE sparse_mod, ONLY: sparse_solve
  IMPLICIT NONE
  
  ! Test parameters
  INTEGER, PARAMETER :: MATRIX_SIZE = 100
  INTEGER, PARAMETER :: MAX_FILL_LEVEL = 5
  REAL(DP), PARAMETER :: TEST_TOLERANCE = 1.0e-10_DP
  
  ! Test variables
  REAL(DP), DIMENSION(MATRIX_SIZE, MATRIX_SIZE) :: A
  REAL(DP), DIMENSION(MATRIX_SIZE) :: b, x_exact, x_bicgstab, residual
  INTEGER :: i, j, fill_k, converged_fill_level
  REAL(DP) :: h, residual_norm, error_norm
  LOGICAL :: converged
  INTEGER :: saved_max_iter, saved_fill_level
  REAL(DP) :: saved_abs_tol, saved_rel_tol
  LOGICAL :: saved_verbose
  
  WRITE(*,'(A)') '============================================='
  WRITE(*,'(A)') 'ILU Fill Level Progression Test'
  WRITE(*,'(A)') 'Testing BiCGSTAB + ILU(k) on Spline Matrices'
  WRITE(*,'(A)') '============================================='
  WRITE(*,*)
  
  ! Save original solver settings
  saved_max_iter = bicgstab_max_iter
  saved_fill_level = ilu_fill_level
  saved_abs_tol = bicgstab_abs_tolerance
  saved_rel_tol = bicgstab_rel_tolerance
  saved_verbose = bicgstab_verbose
  
  ! Configure solver for testing
  sparse_solve_method = SOLVER_BICGSTAB
  bicgstab_max_iter = 1000
  bicgstab_abs_tolerance = TEST_TOLERANCE
  bicgstab_rel_tolerance = TEST_TOLERANCE
  bicgstab_verbose = .FALSE.  ! Keep output clean for progression test
  
  ! Create a challenging spline-like matrix (tridiagonal with additional structure)
  WRITE(*,'(A,I0,A,I0,A)') 'Creating ', MATRIX_SIZE, 'x', MATRIX_SIZE, ' spline-like test matrix'
  h = 1.0_DP / REAL(MATRIX_SIZE + 1, DP)
  A = 0.0_DP
  
  ! Create a discretized second derivative operator with additional coupling
  ! This mimics the structure of spline coefficient matrices
  DO i = 1, MATRIX_SIZE
    ! Main diagonal (stronger than typical)
    A(i, i) = 4.0_DP / (h * h)
    
    ! Off-diagonals
    IF (i > 1) A(i, i-1) = -1.0_DP / (h * h)
    IF (i < MATRIX_SIZE) A(i, i+1) = -1.0_DP / (h * h)
    
    ! Additional coupling for spline-like behavior (makes it more challenging)
    IF (i > 2) A(i, i-2) = 0.1_DP / (h * h)
    IF (i < MATRIX_SIZE-1) A(i, i+2) = 0.1_DP / (h * h)
  END DO
  
  ! Create exact solution and RHS
  DO i = 1, MATRIX_SIZE
    x_exact(i) = SIN(3.14159_DP * REAL(i, DP) * h)
  END DO
  
  ! Compute RHS: b = A * x_exact
  DO i = 1, MATRIX_SIZE
    b(i) = 0.0_DP
    DO j = 1, MATRIX_SIZE
      b(i) = b(i) + A(i, j) * x_exact(j)
    END DO
  END DO
  
  WRITE(*,'(A)') 'Matrix created successfully'
  WRITE(*,'(A,E12.4)') 'Matrix condition estimate (crude): ', MAXVAL(ABS(A)) / MINVAL(ABS(A), MASK=ABS(A)>1e-15)
  WRITE(*,*)
  
  ! Test different ILU fill levels
  converged_fill_level = -1
  
  DO fill_k = 0, MAX_FILL_LEVEL
    WRITE(*,'(A,I0,A)', ADVANCE='NO') 'Testing ILU(', fill_k, ') preconditioning... '
    
    ! Configure ILU fill level
    ilu_fill_level = fill_k
    
    ! Reset solution vector
    x_bicgstab = 0.0_DP
    
    ! Attempt to solve with current fill level
    converged = .TRUE.
    
    ! Use a simple error handling approach
    x_bicgstab = b  ! Copy RHS as initial guess for solver
    
    ! Call sparse solver (this will modify x_bicgstab in place)
    CALL sparse_solve(A, x_bicgstab)
    
    ! Check convergence by computing residual
    residual = 0.0_DP
    DO i = 1, MATRIX_SIZE
      DO j = 1, MATRIX_SIZE
        residual(i) = residual(i) + A(i, j) * x_bicgstab(j)
      END DO
      residual(i) = residual(i) - b(i)
    END DO
    
    ! Compute residual norm
    residual_norm = 0.0_DP
    DO i = 1, MATRIX_SIZE
      residual_norm = residual_norm + residual(i) * residual(i)
    END DO
    residual_norm = SQRT(residual_norm)
    
    ! Compute solution error norm
    error_norm = 0.0_DP
    DO i = 1, MATRIX_SIZE
      error_norm = error_norm + (x_bicgstab(i) - x_exact(i))**2
    END DO
    error_norm = SQRT(error_norm)
    
    ! Check if converged (residual small enough)
    IF (residual_norm < TEST_TOLERANCE * 10.0_DP .AND. error_norm < TEST_TOLERANCE * 100.0_DP) THEN
      WRITE(*,'(A)') '[CONVERGED]'
      WRITE(*,'(A,E12.4)') '    Residual norm: ', residual_norm
      WRITE(*,'(A,E12.4)') '    Solution error: ', error_norm
      IF (converged_fill_level < 0) converged_fill_level = fill_k
    ELSE
      WRITE(*,'(A)') '[FAILED]'
      WRITE(*,'(A,E12.4)') '    Residual norm: ', residual_norm
      WRITE(*,'(A,E12.4)') '    Solution error: ', error_norm
    END IF
    
    WRITE(*,*)
  END DO
  
  ! Summary
  WRITE(*,'(A)') '============================================='
  WRITE(*,'(A)') 'ILU Fill Level Test Results'
  WRITE(*,'(A)') '============================================='
  
  IF (converged_fill_level >= 0) THEN
    WRITE(*,'(A,I0)') 'SUCCESS: BiCGSTAB converged with ILU(', converged_fill_level, ') preconditioning'
    WRITE(*,'(A,I0,A)') 'Minimum fill level required: ', converged_fill_level, ' (higher levels also work)'
    
    IF (converged_fill_level == 0) THEN
      WRITE(*,'(A)') 'NOTE: ILU(0) was sufficient - matrix is reasonably well-conditioned'
    ELSE IF (converged_fill_level == 1) THEN
      WRITE(*,'(A)') 'NOTE: ILU(1) was required - standard preconditioning effective'
    ELSE IF (converged_fill_level == 2) THEN
      WRITE(*,'(A)') 'NOTE: ILU(2) was required - matrix benefits from higher-order preconditioning'
    ELSE
      WRITE(*,'(A,I0,A)') 'NOTE: High fill level ILU(', converged_fill_level, ') was required - challenging matrix'
    END IF
  ELSE
    WRITE(*,'(A,I0,A)') 'FAILURE: BiCGSTAB did not converge even with ILU(', MAX_FILL_LEVEL, ') preconditioning'
    WRITE(*,'(A)') 'This matrix may require:'
    WRITE(*,'(A)') '  - Higher fill levels (k > ', MAX_FILL_LEVEL, ')'
    WRITE(*,'(A)') '  - Different preconditioning strategy'
    WRITE(*,'(A)') '  - Direct solver (UMFPACK)'
    STOP 1
  END IF
  
  WRITE(*,'(A)') '============================================='
  
  ! Restore original solver settings
  bicgstab_max_iter = saved_max_iter
  ilu_fill_level = saved_fill_level
  bicgstab_abs_tolerance = saved_abs_tol
  bicgstab_rel_tolerance = saved_rel_tol
  bicgstab_verbose = saved_verbose
  
END PROGRAM test_ilu_fill_levels