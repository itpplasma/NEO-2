PROGRAM test_bicgstab_l_spline
  ! Test BiCGSTAB(ℓ) with different stabilization parameters and ILU fill levels
  ! on the exact failing spline matrix to find optimal configuration
  
  USE nrtype, ONLY: I4B, DP
  USE sparse_solvers_mod, ONLY: sparse_solve_method, SOLVER_BICGSTAB, SOLVER_UMFPACK, &
                                ilu_fill_level, bicgstab_max_iter, &
                                bicgstab_abs_tolerance, bicgstab_rel_tolerance, &
                                bicgstab_stabilization_param, bicgstab_verbose
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

  ! Test parameters for systematic exploration
  INTEGER(I4B), PARAMETER :: n_intervals = 50
  INTEGER(I4B), PARAMETER :: MAX_L_PARAM = 4      ! Test BiCGSTAB(1) through BiCGSTAB(4)
  INTEGER(I4B), PARAMETER :: MAX_ILU_LEVEL = 5    ! Test ILU(0) through ILU(5)
  REAL(DP), PARAMETER :: BASE_TOLERANCE = 1.0e-10_DP
  
  ! Test variables
  INTEGER(I4B) :: n_points, i, l_val, ilu_k
  REAL(DP), ALLOCATABLE :: x(:), y(:), lambda1(:)
  INTEGER(I4B), ALLOCATABLE :: indx(:)
  REAL(DP), ALLOCATABLE :: a_ref(:), b_ref(:), c_ref(:), d_ref(:)  ! Reference (UMFPACK)
  REAL(DP), ALLOCATABLE :: a_test(:), b_test(:), c_test(:), d_test(:)  ! Test (BiCGSTAB(ℓ))
  REAL(DP) :: c1, cn, m
  INTEGER(I4B) :: sw1, sw2
  REAL(DP) :: max_diff_a, max_diff_b, max_diff_c, max_diff_d, total_error
  LOGICAL :: converged
  INTEGER(I4B) :: optimal_l, optimal_k
  REAL(DP) :: best_error
  
  ! Result tracking
  LOGICAL :: success_matrix(1:MAX_L_PARAM, 0:MAX_ILU_LEVEL)
  REAL(DP) :: error_matrix(1:MAX_L_PARAM, 0:MAX_ILU_LEVEL)
  
  ! Analysis variables
  INTEGER :: success_count, successes_for_l, successes_for_k, min_k_for_l
  REAL(DP) :: min_error_for_l
  
  ! Saved settings
  INTEGER :: saved_method, saved_max_iter, saved_fill_level, saved_l_param
  REAL(DP) :: saved_abs_tol, saved_rel_tol
  LOGICAL :: saved_verbose
  
  WRITE(*,'(A)') '============================================='
  WRITE(*,'(A)') 'BiCGSTAB(ℓ) + ILU(k) Systematic Test'
  WRITE(*,'(A)') 'Testing on Pathological Spline Matrix'
  WRITE(*,'(A)') '============================================='
  WRITE(*,*)
  
  ! Save original solver settings
  saved_method = sparse_solve_method
  saved_max_iter = bicgstab_max_iter
  saved_fill_level = ilu_fill_level
  saved_l_param = bicgstab_stabilization_param
  saved_abs_tol = bicgstab_abs_tolerance
  saved_rel_tol = bicgstab_rel_tolerance
  saved_verbose = bicgstab_verbose
  
  ! Setup the exact same problematic spline from test_spline_bicgstab_accuracy.f90
  n_points = n_intervals * 5  ! This creates the challenging problem
  
  ALLOCATE(x(n_points), y(n_points), lambda1(n_points), indx(n_points))
  ALLOCATE(a_ref(n_points), b_ref(n_points), c_ref(n_points), d_ref(n_points))
  ALLOCATE(a_test(n_points), b_test(n_points), c_test(n_points), d_test(n_points))
  
  ! Create the exact same challenging spline data
  DO i = 1, n_points
    x(i) = REAL(i-1, DP) / REAL(n_points-1, DP) * 10.0_DP
    y(i) = SIN(2.0_DP * 3.14159_DP * x(i)) + 0.5_DP * EXP(-x(i))
    lambda1(i) = 1.0e-6_DP  ! Very small smoothing - creates the pathological matrix!
    indx(i) = i
  END DO
  
  c1 = 0.0_DP
  cn = 0.0_DP
  sw1 = 2  ! From working test
  sw2 = 4  ! From working test
  m = 0.0_DP
  
  WRITE(*,'(A,I0,A)') 'Problem size: ', n_points, ' points'
  WRITE(*,'(A)') 'Spline function: sin(2πx) + 0.5*exp(-x) with λ=1e-6 (pathological!)'
  WRITE(*,*)
  
  ! Get reference solution with UMFPACK
  WRITE(*,'(A)') 'Computing reference solution with UMFPACK...'
  sparse_solve_method = SOLVER_UMFPACK
  CALL splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                    a_ref, b_ref, c_ref, d_ref, m, test_function)
  WRITE(*,'(A)') 'Reference solution computed successfully'
  WRITE(*,*)
  
  ! Configure BiCGSTAB settings for systematic testing
  sparse_solve_method = SOLVER_BICGSTAB
  bicgstab_max_iter = 3000  ! Allow more iterations for difficult cases
  bicgstab_abs_tolerance = BASE_TOLERANCE
  bicgstab_rel_tolerance = BASE_TOLERANCE
  bicgstab_verbose = .TRUE.   ! Enable verbose to debug algorithm path
  
  ! Initialize tracking
  success_matrix = .FALSE.
  error_matrix = HUGE(1.0_DP)
  best_error = HUGE(1.0_DP)
  optimal_l = -1
  optimal_k = -1
  
  WRITE(*,'(A)') 'Systematic BiCGSTAB(ℓ) + ILU(k) Testing:'
  WRITE(*,'(A)') 'ℓ\\k     ILU(0)    ILU(1)    ILU(2)    ILU(3)    ILU(4)    ILU(5)'
  WRITE(*,'(A)') '--- ---------- ---------- ---------- ---------- ---------- ----------'
  
  ! Systematic test: all combinations of ℓ and k
  DO l_val = 1, MAX_L_PARAM
    WRITE(*,'(I0,A)', ADVANCE='NO') l_val, '   '
    
    DO ilu_k = 0, MAX_ILU_LEVEL
      ! Configure solver parameters
      bicgstab_stabilization_param = l_val
      ilu_fill_level = ilu_k
      
      ! Reset boundary conditions (splinecof3_a modifies these)
      c1 = 0.0_DP
      cn = 0.0_DP
      
      ! Attempt to solve with current configuration
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
      error_matrix(l_val, ilu_k) = total_error
      
      ! Check convergence (should match UMFPACK to reasonable precision)
      IF (total_error < BASE_TOLERANCE * 1000.0_DP) THEN
        success_matrix(l_val, ilu_k) = .TRUE.
        WRITE(*,'(A)', ADVANCE='NO') '    PASS  '
        
        ! Track best configuration
        IF (total_error < best_error) THEN
          best_error = total_error
          optimal_l = l_val
          optimal_k = ilu_k
        END IF
      ELSE
        WRITE(*,'(A)', ADVANCE='NO') '    FAIL  '
      END IF
      
    END DO
    
    WRITE(*,*)  ! End line for this ℓ value
  END DO
  
  WRITE(*,*)
  
  ! Analysis and summary
  WRITE(*,'(A)') '============================================='
  WRITE(*,'(A)') 'BiCGSTAB(ℓ) + ILU(k) Test Results Analysis'
  WRITE(*,'(A)') '============================================='
  
  ! Count successes
  success_count = COUNT(success_matrix)
  
  IF (success_count > 0) THEN
    WRITE(*,'(A,I0,A,I0,A)') 'SUCCESS: ', success_count, ' out of ', &
                             MAX_L_PARAM * (MAX_ILU_LEVEL + 1), ' configurations worked'
    WRITE(*,'(A,I0,A,I0,A)') 'OPTIMAL: BiCGSTAB(', optimal_l, ') + ILU(', optimal_k, ')'
    WRITE(*,'(A,E12.4)') 'Best error achieved:', best_error
    WRITE(*,*)
    
    ! Analysis by ℓ parameter
    WRITE(*,'(A)') 'Analysis by stabilization parameter ℓ:'
    DO l_val = 1, MAX_L_PARAM
      successes_for_l = COUNT(success_matrix(l_val, :))
      WRITE(*,'(A,I0,A,I0,A,I0,A)') '  BiCGSTAB(', l_val, '): ', successes_for_l, ' / ', &
                                     MAX_ILU_LEVEL + 1, ' ILU levels work'
      
      IF (successes_for_l > 0) THEN
        ! Find minimum k that works for this ℓ
        min_k_for_l = -1
        DO ilu_k = 0, MAX_ILU_LEVEL
          IF (success_matrix(l_val, ilu_k)) THEN
            min_k_for_l = ilu_k
            EXIT
          END IF
        END DO
        WRITE(*,'(A,I0)') '    Minimum ILU level needed: k = ', min_k_for_l
      END IF
    END DO
    
    WRITE(*,*)
    
    ! Analysis by ILU level
    WRITE(*,'(A)') 'Analysis by ILU fill level k:'
    DO ilu_k = 0, MAX_ILU_LEVEL
      successes_for_k = COUNT(success_matrix(:, ilu_k))
      WRITE(*,'(A,I0,A,I0,A,I0,A)') '  ILU(', ilu_k, '): ', successes_for_k, ' / ', &
                                     MAX_L_PARAM, ' BiCGSTAB(ℓ) variants work'
    END DO
    
    WRITE(*,*)
    
    ! Recommendations
    WRITE(*,'(A)') 'RECOMMENDATIONS:'
    IF (optimal_l == 1) THEN
      WRITE(*,'(A)') '  - Standard BiCGSTAB(1) sufficient with proper ILU preconditioning'
    ELSE IF (optimal_l == 2) THEN
      WRITE(*,'(A)') '  - BiCGSTAB(2) provides better stability for this matrix'
    ELSE
      WRITE(*,'(A,I0,A)') '  - Higher-order BiCGSTAB(', optimal_l, ') required for this pathological case'
    END IF
    
    IF (optimal_k == 0) THEN
      WRITE(*,'(A)') '  - No ILU preconditioning needed (matrix better conditioned than expected)'
    ELSE IF (optimal_k == 1) THEN
      WRITE(*,'(A)') '  - ILU(1) preconditioning confirms standard choice'
    ELSE
      WRITE(*,'(A,I0,A)') '  - Higher fill level ILU(', optimal_k, ') required for this ill-conditioned matrix'
    END IF
    
  ELSE
    WRITE(*,'(A)') 'COMPLETE FAILURE: No BiCGSTAB(ℓ) + ILU(k) combination worked!'
    WRITE(*,'(A)') 'This spline matrix is extremely pathological.'
    WRITE(*,'(A)') ''
    WRITE(*,'(A)') 'Recommendations:'
    WRITE(*,'(A)') '  - Use UMFPACK (direct solver) for reliability'
    WRITE(*,'(A)') '  - Consider GMRES with ILU preconditioning'
    WRITE(*,'(A)') '  - Modify spline parameters (larger λ) to reduce ill-conditioning'
    WRITE(*,*)
    
    ! Show worst errors for debugging
    WRITE(*,'(A)') 'Error analysis (smallest errors achieved):'
    DO l_val = 1, MAX_L_PARAM
      min_error_for_l = MINVAL(error_matrix(l_val, :))
      WRITE(*,'(A,I0,A,E12.4)') '  BiCGSTAB(', l_val, ') best error: ', min_error_for_l
    END DO
    
    STOP 1
  END IF
  
  WRITE(*,'(A)') '============================================='
  
  ! Restore original solver settings
  sparse_solve_method = saved_method
  bicgstab_max_iter = saved_max_iter
  ilu_fill_level = saved_fill_level
  bicgstab_stabilization_param = saved_l_param
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
  
END PROGRAM test_bicgstab_l_spline