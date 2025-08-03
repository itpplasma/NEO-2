PROGRAM test_spline_ilu_fill_levels
  ! Test different ILU fill levels on the actual strange spline matrix
  ! This uses the real spline coefficient matrix that was causing convergence issues
  
  USE nrtype, ONLY: I4B, DP
  USE sparse_solvers_mod, ONLY: sparse_solve_method, SOLVER_BICGSTAB, SOLVER_UMFPACK, &
                                ilu_fill_level, bicgstab_max_iter, &
                                bicgstab_abs_tolerance, bicgstab_rel_tolerance, &
                                bicgstab_verbose
  IMPLICIT NONE
  
  INTERFACE
    SUBROUTINE splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
         a, b, c, d, m, f)
      USE nrtype, ONLY: I4B, DP
      REAL(DP),                   INTENT(inout) :: c1, cn
      REAL(DP),     DIMENSION(:), INTENT(in)    :: x, y, lambda1
      INTEGER(I4B), DIMENSION(:), INTENT(in)    :: indx
      REAL(DP),     DIMENSION(:), INTENT(out)   :: a, b, c, d
      INTEGER(I4B),               INTENT(in)    :: sw1, sw2
      REAL(DP),                   INTENT(in)    :: m
      INTERFACE
        FUNCTION f(x,m)
          USE nrtype, ONLY : DP
          IMPLICIT NONE
          REAL(DP), INTENT(in) :: x, m
          REAL(DP)             :: f
        END FUNCTION f
      END INTERFACE
    END SUBROUTINE splinecof3_a
  END INTERFACE

  ! Test parameters
  INTEGER(I4B), PARAMETER :: n_intervals = 50
  INTEGER(I4B), PARAMETER :: MAX_FILL_LEVEL = 5
  REAL(DP), PARAMETER :: BASE_TOLERANCE = 1.0e-12_DP
  INTEGER(I4B) :: n_points, i, fill_k, converged_fill_level
  REAL(DP), ALLOCATABLE :: x(:), y(:), lambda1(:)
  INTEGER(I4B), ALLOCATABLE :: indx(:)
  REAL(DP), ALLOCATABLE :: a_ref(:), b_ref(:), c_ref(:), d_ref(:)  ! Reference (UMFPACK)
  REAL(DP), ALLOCATABLE :: a_test(:), b_test(:), c_test(:), d_test(:)  ! Test (BiCGSTAB)
  REAL(DP) :: c1, cn, m
  INTEGER(I4B) :: sw1, sw2
  REAL(DP) :: max_diff_a, max_diff_b, max_diff_c, max_diff_d, total_error
  
  ! Saved settings
  INTEGER :: saved_method, saved_max_iter, saved_fill_level
  REAL(DP) :: saved_abs_tol, saved_rel_tol
  LOGICAL :: saved_verbose
  
  WRITE(*,'(A)') '============================================='
  WRITE(*,'(A)') 'ILU Fill Level Test on Real Spline Matrix'
  WRITE(*,'(A)') 'Testing BiCGSTAB + ILU(k) on Actual Spline Problem'
  WRITE(*,'(A)') '============================================='
  WRITE(*,*)
  
  ! Save original solver settings
  saved_method = sparse_solve_method
  saved_max_iter = bicgstab_max_iter
  saved_fill_level = ilu_fill_level
  saved_abs_tol = bicgstab_abs_tolerance
  saved_rel_tol = bicgstab_rel_tolerance
  saved_verbose = bicgstab_verbose
  
  ! Setup the EXACT same problematic spline from test_spline_bicgstab_accuracy.f90
  n_points = n_intervals * 5  ! This gives 250 points, but the failing matrix is 348x348
  ! Let's use the exact setup that creates the 348x348 failing matrix
  n_points = 58  ! This should create the problematic 348x348 matrix
  
  ALLOCATE(x(n_points), y(n_points), lambda1(n_points), indx(n_points))
  ALLOCATE(a_ref(n_points), b_ref(n_points), c_ref(n_points), d_ref(n_points))
  ALLOCATE(a_test(n_points), b_test(n_points), c_test(n_points), d_test(n_points))
  
  ! Create the EXACT same problematic spline data that fails
  DO i = 1, n_points
    x(i) = REAL(i-1, DP) / REAL(n_points-1, DP) * 10.0_DP  ! Scale to match failing case
    y(i) = SIN(2.0_DP * 3.14159_DP * x(i)) + 0.5_DP * EXP(-x(i))  ! Complex function
    lambda1(i) = 1.0e-6_DP  ! Very small smoothing - creates ill-conditioning!
    indx(i) = i
  END DO
  
  c1 = 0.0_DP
  cn = 0.0_DP
  sw1 = 2  ! Boundary condition type (from working test)
  sw2 = 4  ! Boundary condition type (from working test)
  m = 0.0_DP
  
  WRITE(*,'(A,I0)') 'Problem size: ', n_points, ' points'
  WRITE(*,'(A)') 'Spline function: sin(2πx) + 0.1*sin(25x) with small smoothing'
  WRITE(*,*)
  
  ! Get reference solution with UMFPACK
  WRITE(*,'(A)') 'Computing reference solution with UMFPACK...'
  sparse_solve_method = SOLVER_UMFPACK
  CALL splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                    a_ref, b_ref, c_ref, d_ref, m, test_function)
  WRITE(*,'(A)') 'Reference solution computed successfully'
  WRITE(*,*)
  
  ! Configure BiCGSTAB settings
  sparse_solve_method = SOLVER_BICGSTAB
  bicgstab_max_iter = 2000
  bicgstab_abs_tolerance = BASE_TOLERANCE
  bicgstab_rel_tolerance = BASE_TOLERANCE
  bicgstab_verbose = .FALSE.
  
  ! Test different ILU fill levels
  converged_fill_level = -1
  
  DO fill_k = 0, MAX_FILL_LEVEL
    WRITE(*,'(A,I0,A)', ADVANCE='NO') 'Testing ILU(', fill_k, ') preconditioning... '
    
    ! Configure ILU fill level
    ilu_fill_level = fill_k
    
    ! Reset boundary conditions (splinecof3_a modifies these)
    c1 = 0.0_DP
    cn = 0.0_DP
    
    ! Attempt to solve with current fill level
    CALL splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                      a_test, b_test, c_test, d_test, m, test_function)
    
    ! Compute differences from reference solution
    max_diff_a = 0.0_DP
    max_diff_b = 0.0_DP
    max_diff_c = 0.0_DP
    max_diff_d = 0.0_DP
    
    DO i = 1, n_points
      max_diff_a = MAX(max_diff_a, ABS(a_test(i) - a_ref(i)))
      max_diff_b = MAX(max_diff_b, ABS(b_test(i) - b_ref(i)))
      max_diff_c = MAX(max_diff_c, ABS(c_test(i) - c_ref(i)))
      max_diff_d = MAX(max_diff_d, ABS(d_test(i) - d_ref(i)))
    END DO
    
    total_error = max_diff_a + max_diff_b + max_diff_c + max_diff_d
    
    ! Check convergence (should match UMFPACK to reasonable precision)
    IF (total_error < BASE_TOLERANCE * 1000.0_DP) THEN
      WRITE(*,'(A)') '[CONVERGED]'
      WRITE(*,'(A,E12.4)') '    Total coefficient error: ', total_error
      WRITE(*,'(A,E12.4,A,E12.4,A,E12.4,A,E12.4)') '    Max errors: a=', max_diff_a, &
                                                     ', b=', max_diff_b, ', c=', max_diff_c, ', d=', max_diff_d
      IF (converged_fill_level < 0) converged_fill_level = fill_k
    ELSE
      WRITE(*,'(A)') '[FAILED]'
      WRITE(*,'(A,E12.4)') '    Total coefficient error: ', total_error
      WRITE(*,'(A,E12.4,A,E12.4,A,E12.4,A,E12.4)') '    Max errors: a=', max_diff_a, &
                                                     ', b=', max_diff_b, ', c=', max_diff_c, ', d=', max_diff_d
    END IF
    
    WRITE(*,*)
  END DO
  
  ! Summary
  WRITE(*,'(A)') '============================================='
  WRITE(*,'(A)') 'Real Spline Matrix ILU Test Results'
  WRITE(*,'(A)') '============================================='
  
  IF (converged_fill_level >= 0) THEN
    WRITE(*,'(A,I0)') 'SUCCESS: BiCGSTAB converged with ILU(', converged_fill_level, ') preconditioning'
    WRITE(*,'(A,I0,A)') 'Minimum fill level required: ', converged_fill_level, ' for the strange spline matrix'
    
    IF (converged_fill_level == 0) THEN
      WRITE(*,'(A)') 'SURPRISING: ILU(0) was sufficient - this spline matrix is better conditioned than expected'
    ELSE IF (converged_fill_level == 1) THEN
      WRITE(*,'(A)') 'EXPECTED: ILU(1) required - confirms standard preconditioning choice'
    ELSE IF (converged_fill_level == 2) THEN
      WRITE(*,'(A)') 'INTERESTING: ILU(2) required - this spline matrix benefits from higher-order preconditioning'
    ELSE
      WRITE(*,'(A,I0,A)') 'CHALLENGING: High fill level ILU(', converged_fill_level, ') required - very difficult spline matrix'
    END IF
  ELSE
    WRITE(*,'(A,I0,A)') 'FAILURE: BiCGSTAB did not converge even with ILU(', MAX_FILL_LEVEL, ') preconditioning'
    WRITE(*,'(A)') 'This real spline matrix is extremely challenging and may require:'
    WRITE(*,'(A)') '  - Even higher fill levels (k > ', MAX_FILL_LEVEL, ')'
    WRITE(*,'(A)') '  - Alternative preconditioning (AMG, etc.)'
    WRITE(*,'(A)') '  - Direct solver (UMFPACK) for reliability'
    WRITE(*,'(A)') '  - Different discretization or smoothing parameters'
    STOP 1
  END IF
  
  WRITE(*,'(A)') '============================================='
  
  ! Restore original solver settings
  sparse_solve_method = saved_method
  bicgstab_max_iter = saved_max_iter
  ilu_fill_level = saved_fill_level
  bicgstab_abs_tolerance = saved_abs_tol
  bicgstab_rel_tolerance = saved_rel_tol
  bicgstab_verbose = saved_verbose
  
  DEALLOCATE(x, y, lambda1, indx)
  DEALLOCATE(a_ref, b_ref, c_ref, d_ref, a_test, b_test, c_test, d_test)

CONTAINS

  FUNCTION test_function(x, m) RESULT(result)
    REAL(DP), INTENT(in) :: x, m
    REAL(DP) :: result
    result = 1.0_DP  ! Simple test function
  END FUNCTION test_function
  
END PROGRAM test_spline_ilu_fill_levels