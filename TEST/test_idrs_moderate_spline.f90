PROGRAM test_idrs_moderate_spline
  !> Test IDR(s) on moderately-conditioned spline problems

  USE nrtype, ONLY: I4B, DP
  USE sparse_solvers_mod, ONLY: sparse_solve_idrs_real, default_iterative_params, &
                                PRECOND_NONE, sparse_solve_method, SOLVER_UMFPACK, SOLVER_IDRS
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

  ! Test different problem sizes and conditioning
  INTEGER(I4B), PARAMETER :: n_tests = 4
  INTEGER(I4B), PARAMETER :: test_sizes(n_tests) = [10, 20, 30, 50]
  REAL(DP), PARAMETER :: test_lambdas(n_tests) = [1.0e-2_DP, 1.0e-3_DP, 1.0e-4_DP, 1.0e-5_DP]
  CHARACTER(20), PARAMETER :: test_names(n_tests) = ['Well-conditioned  ', &
                                                     'Moderate condition', &
                                                     'Poor condition    ', &
                                                     'Very poor condition']
  
  ! Variables
  INTEGER(I4B) :: test_id, i, n_points
  REAL(DP) :: lambda_val, total_error
  REAL(DP), ALLOCATABLE :: x(:), y(:), lambda1(:)
  INTEGER(I4B), ALLOCATABLE :: indx(:)
  REAL(DP), ALLOCATABLE :: a_idrs(:), b_idrs(:), c_idrs(:), d_idrs(:)
  REAL(DP), ALLOCATABLE :: a_umf(:), b_umf(:), c_umf(:), d_umf(:)
  REAL(DP) :: c1, cn, m
  INTEGER(I4B) :: sw1, sw2
  LOGICAL :: idrs_success
  
  WRITE(*,'(A)') '======================================================='
  WRITE(*,'(A)') 'Testing IDR(s) on Various Spline Problems'
  WRITE(*,'(A)') '======================================================='
  WRITE(*,*)
  
  ! Configure IDR(s) once
  default_iterative_params%preconditioner_type = PRECOND_NONE
  default_iterative_params%idrs_shadow_space_dim = 4
  default_iterative_params%max_iterations = 1000
  default_iterative_params%abs_tolerance = 1.0E-10_DP
  default_iterative_params%rel_tolerance = 1.0E-8_DP
  default_iterative_params%verbose = .FALSE.
  
  ! Test loop
  DO test_id = 1, n_tests
    n_points = test_sizes(test_id)
    lambda_val = test_lambdas(test_id)
    
    WRITE(*,'(A,I0,A,A,A)') 'Test ', test_id, ': ', TRIM(test_names(test_id)), ' spline'
    WRITE(*,'(A,I0,A,ES9.2)') '  Size: ', n_points, ' points, lambda = ', lambda_val
    
    ! Allocate arrays
    ALLOCATE(x(n_points), y(n_points), lambda1(n_points), indx(n_points))
    ALLOCATE(a_idrs(n_points), b_idrs(n_points), c_idrs(n_points), d_idrs(n_points))
    ALLOCATE(a_umf(n_points), b_umf(n_points), c_umf(n_points), d_umf(n_points))
    
    ! Create simple sine wave data
    DO i = 1, n_points
      x(i) = REAL(i-1, DP) / REAL(n_points-1, DP) * 6.28_DP  ! 0 to 2*pi
      y(i) = SIN(x(i))
      lambda1(i) = lambda_val
      indx(i) = i
    END DO
    
    ! Boundary conditions
    c1 = 0.0_DP
    cn = 0.0_DP
    sw1 = 2
    sw2 = 4
    m = 0.0_DP
    
    ! Get UMFPACK reference
    sparse_solve_method = SOLVER_UMFPACK
    CALL splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                      a_umf, b_umf, c_umf, d_umf, m, simple_function)
    
    ! Test IDR(s)
    sparse_solve_method = SOLVER_IDRS
    CALL splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                      a_idrs, b_idrs, c_idrs, d_idrs, m, simple_function)
    
    ! Check results
    IF (ANY(a_idrs /= a_idrs) .OR. ANY(b_idrs /= b_idrs) .OR. &
        ANY(c_idrs /= c_idrs) .OR. ANY(d_idrs /= d_idrs)) THEN
      WRITE(*,'(A)') '  [FAILURE] IDR(s) produced NaN'
      idrs_success = .FALSE.
    ELSE
      ! Calculate error vs UMFPACK
      total_error = MAXVAL(ABS(a_idrs - a_umf)) + MAXVAL(ABS(b_idrs - b_umf)) + &
                    MAXVAL(ABS(c_idrs - c_umf)) + MAXVAL(ABS(d_idrs - d_umf))
      
      IF (total_error < 1.0e-6_DP) THEN
        WRITE(*,'(A,ES9.2)') '  [SUCCESS] IDR(s) error vs UMFPACK: ', total_error
        idrs_success = .TRUE.
      ELSE
        WRITE(*,'(A,ES9.2)') '  [FAILURE] IDR(s) large error vs UMFPACK: ', total_error
        idrs_success = .FALSE.
      END IF
    END IF
    
    ! Clean up for next test
    DEALLOCATE(x, y, lambda1, indx)
    DEALLOCATE(a_idrs, b_idrs, c_idrs, d_idrs, a_umf, b_umf, c_umf, d_umf)
    
    WRITE(*,*)
  END DO
  
  WRITE(*,'(A)') '======================================================='

CONTAINS

  FUNCTION simple_function(x, m) RESULT(result)
    REAL(DP), INTENT(in) :: x, m
    REAL(DP) :: result
    result = 1.0_DP
  END FUNCTION simple_function

END PROGRAM test_idrs_moderate_spline