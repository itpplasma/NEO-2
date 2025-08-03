PROGRAM test_idrs_404_spline
  !> Test IDR(s) solver without preconditioning on 404x404 spline matrix
  
  USE nrtype, ONLY: I4B, DP
  USE sparse_solvers_mod, ONLY: sparse_solve_idrs_real, default_iterative_params, &
                                PRECOND_NONE, PRECOND_AMG, sparse_solve_method, SOLVER_UMFPACK, SOLVER_IDRS
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

  ! Test parameters - let's try smaller, better-conditioned problems
  INTEGER(I4B), PARAMETER :: n_points = 67  ! Back to 404x404 problem
  REAL(DP), PARAMETER :: lambda_small = 1.0e-6_DP  ! Original challenging conditioning
  
  ! Variables
  INTEGER(I4B) :: i, info, actual_matrix_size
  REAL(DP), ALLOCATABLE :: x(:), y(:), lambda1(:)
  INTEGER(I4B), ALLOCATABLE :: indx(:)
  REAL(DP), ALLOCATABLE :: a_idrs(:), b_idrs(:), c_idrs(:), d_idrs(:)
  REAL(DP), ALLOCATABLE :: a_umf(:), b_umf(:), c_umf(:), d_umf(:)  ! UMFPACK reference
  REAL(DP) :: c1, cn, m, total_error
  INTEGER(I4B) :: sw1, sw2
  
  WRITE(*,'(A)') '========================================================='
  WRITE(*,'(A)') 'Testing IDR(s) on 404x404 Challenging Spline Matrix'
  WRITE(*,'(A)') '========================================================='
  
  ! Allocate arrays - using the exact same structure as test_spline_amg
  ALLOCATE(x(n_points), y(n_points), lambda1(n_points), indx(n_points))
  
  ! Create the exact same problematic spline data as test_spline_amg
  DO i = 1, n_points
    x(i) = REAL(i-1, DP) / REAL(n_points-1, DP) * 10.0_DP
    ! Complex oscillatory function that creates ill-conditioning
    y(i) = SIN(2.0_DP * 3.14159_DP * x(i)) + 0.1_DP * SIN(25.0_DP * x(i))
    lambda1(i) = lambda_small  ! Very small smoothing - critical for ill-conditioning!
    indx(i) = i
  END DO
  
  ! Boundary conditions - same as test_spline_amg
  c1 = 0.0_DP
  cn = 0.0_DP
  sw1 = 2  ! Natural spline
  sw2 = 4  ! Natural spline  
  m = 0.0_DP
  
  WRITE(*,'(A,I3)') 'Data points: ', n_points
  WRITE(*,'(A,ES10.3)') 'Smoothing parameter (lambda): ', lambda_small
  WRITE(*,*)
  
  ! Allocate result arrays
  ALLOCATE(a_idrs(n_points), b_idrs(n_points), c_idrs(n_points), d_idrs(n_points))
  ALLOCATE(a_umf(n_points), b_umf(n_points), c_umf(n_points), d_umf(n_points))
  
  ! Test 1: Reference solution with UMFPACK
  WRITE(*,'(A)') '1. Computing reference solution with UMFPACK...'
  sparse_solve_method = SOLVER_UMFPACK
  CALL splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                    a_umf, b_umf, c_umf, d_umf, m, dummy_function)
  
  ! Check UMFPACK results
  IF (ANY(a_umf /= a_umf) .OR. ANY(b_umf /= b_umf) .OR. &
      ANY(c_umf /= c_umf) .OR. ANY(d_umf /= d_umf)) THEN
    WRITE(*,'(A)') '   [FAILURE] UMFPACK produced NaN values'
  ELSE
    WRITE(*,'(A)') '   [SUCCESS] UMFPACK computed solution'
    WRITE(*,'(A,ES12.5)') '   Max |a| coefficient: ', MAXVAL(ABS(a_umf))
    WRITE(*,'(A,ES12.5)') '   Max |b| coefficient: ', MAXVAL(ABS(b_umf))
    WRITE(*,'(A,ES12.5)') '   Max |c| coefficient: ', MAXVAL(ABS(c_umf))
    WRITE(*,'(A,ES12.5)') '   Max |d| coefficient: ', MAXVAL(ABS(d_umf))
  END IF
  WRITE(*,*)
  
  ! Test 2: IDR(s) with AMG preconditioning (test the fix!)
  WRITE(*,'(A)') '2. Testing IDR(s) WITH AMG preconditioning (FIXED VERSION)...'
  
  ! Configure IDR(s) with AMG preconditioning
  default_iterative_params%preconditioner_type = PRECOND_AMG
  default_iterative_params%idrs_shadow_space_dim = 4
  default_iterative_params%max_iterations = 1000
  default_iterative_params%abs_tolerance = 1.0E-10_DP
  default_iterative_params%rel_tolerance = 1.0E-8_DP
  default_iterative_params%verbose = .TRUE.  ! Enable verbose to see convergence
  
  ! Call spline with IDR(s)+AMG
  sparse_solve_method = SOLVER_IDRS
  CALL splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                    a_idrs, b_idrs, c_idrs, d_idrs, m, dummy_function)
  
  ! Check if spline coefficients are reasonable
  IF (ANY(a_idrs /= a_idrs) .OR. ANY(b_idrs /= b_idrs) .OR. &
      ANY(c_idrs /= c_idrs) .OR. ANY(d_idrs /= d_idrs)) THEN
    WRITE(*,'(A)') '   [FAILURE] IDR(s)+AMG produced NaN values'
  ELSE
    WRITE(*,'(A)') '   [SUCCESS] IDR(s)+AMG computed coefficients without NaN!'
    
    ! Check magnitude of coefficients
    WRITE(*,'(A,ES12.5)') '   Max |a| coefficient: ', MAXVAL(ABS(a_idrs))
    WRITE(*,'(A,ES12.5)') '   Max |b| coefficient: ', MAXVAL(ABS(b_idrs))
    WRITE(*,'(A,ES12.5)') '   Max |c| coefficient: ', MAXVAL(ABS(c_idrs))
    WRITE(*,'(A,ES12.5)') '   Max |d| coefficient: ', MAXVAL(ABS(d_idrs))
    
    ! Compare with UMFPACK reference
    WRITE(*,*)
    WRITE(*,'(A)') '3. Comparison with UMFPACK reference:'
    WRITE(*,'(A,ES12.5)') '   Max difference in a: ', MAXVAL(ABS(a_idrs - a_umf))
    WRITE(*,'(A,ES12.5)') '   Max difference in b: ', MAXVAL(ABS(b_idrs - b_umf))
    WRITE(*,'(A,ES12.5)') '   Max difference in c: ', MAXVAL(ABS(c_idrs - c_umf))
    WRITE(*,'(A,ES12.5)') '   Max difference in d: ', MAXVAL(ABS(d_idrs - d_umf))
    
    ! Check if solutions are close
    total_error = MAXVAL(ABS(a_idrs - a_umf)) + MAXVAL(ABS(b_idrs - b_umf)) + &
                  MAXVAL(ABS(c_idrs - c_umf)) + MAXVAL(ABS(d_idrs - d_umf))
    
    IF (total_error < 1.0e-3_DP) THEN
      WRITE(*,'(A)') '   [SUCCESS] IDR(s)+AMG results match UMFPACK within tolerance'
    ELSE
      WRITE(*,'(A,ES12.5)') '   [FAILURE] IDR(s)+AMG differs from UMFPACK, total error = ', total_error
    END IF
  END IF
  
  ! Clean up
  DEALLOCATE(x, y, lambda1, indx)
  DEALLOCATE(a_idrs, b_idrs, c_idrs, d_idrs)
  
  WRITE(*,'(A)') '========================================================='

CONTAINS

  FUNCTION dummy_function(x, m) RESULT(result)
    REAL(DP), INTENT(in) :: x, m
    REAL(DP) :: result
    result = 1.0_DP  ! Same as test_spline_amg.f90
  END FUNCTION dummy_function

END PROGRAM test_idrs_404_spline