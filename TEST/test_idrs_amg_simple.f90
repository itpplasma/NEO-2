PROGRAM test_idrs_amg_simple
  !> Test native sparse IDR(s)+AMG on a well-conditioned spline problem
  
  USE nrtype, ONLY: I4B, DP
  USE sparse_solvers_mod, ONLY: sparse_solve_method, &
                                SOLVER_UMFPACK, SOLVER_IDRS, &
                                PRECOND_AMG, &
                                default_iterative_params
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

  ! Test parameters - small, well-conditioned spline 
  INTEGER(I4B), PARAMETER :: n_points = 10
  REAL(DP), PARAMETER :: lambda_good = 1.0e-3_DP  ! Well-conditioned
  
  ! Variables
  INTEGER(I4B) :: i
  REAL(DP), ALLOCATABLE :: x(:), y(:), lambda1(:)
  INTEGER(I4B), ALLOCATABLE :: indx(:)
  REAL(DP), ALLOCATABLE :: a_ref(:), b_ref(:), c_ref(:), d_ref(:)  ! Reference (UMFPACK)
  REAL(DP), ALLOCATABLE :: a_idrs(:), b_idrs(:), c_idrs(:), d_idrs(:)  ! IDR(s)+AMG
  REAL(DP) :: c1, cn, m
  INTEGER(I4B) :: sw1, sw2
  REAL(DP) :: error_idrs
  LOGICAL :: idrs_success
  
  ! Saved settings
  INTEGER :: saved_method, saved_precond
  LOGICAL :: saved_verbose
  
  WRITE(*,'(A)') '========================================================='
  WRITE(*,'(A)') 'IDR(s)+AMG Simple Test: Well-Conditioned Small Spline' 
  WRITE(*,'(A)') '========================================================='
  WRITE(*,*)
  
  ! Save original solver settings
  saved_method = sparse_solve_method
  saved_precond = default_iterative_params%preconditioner_type
  saved_verbose = default_iterative_params%verbose
  
  ! Allocate arrays
  ALLOCATE(x(n_points), y(n_points), lambda1(n_points), indx(n_points))
  ALLOCATE(a_ref(n_points), b_ref(n_points), c_ref(n_points), d_ref(n_points))
  ALLOCATE(a_idrs(n_points), b_idrs(n_points), c_idrs(n_points), d_idrs(n_points))
  
  ! Create well-conditioned spline data
  DO i = 1, n_points
    x(i) = REAL(i-1, DP) / REAL(n_points-1, DP) * 5.0_DP
    ! Simple smooth function
    y(i) = SIN(x(i)) + 0.1_DP * x(i)
    lambda1(i) = lambda_good  ! Well-conditioned smoothing
    indx(i) = i
  END DO
  
  ! Boundary conditions
  c1 = 0.0_DP
  cn = 0.0_DP
  sw1 = 2  ! Natural spline
  sw2 = 4  ! Natural spline  
  m = 0.0_DP
  
  WRITE(*,'(A,I0,A)') 'Problem size: ', n_points, ' data points'
  WRITE(*,'(A,ES10.3)') 'Smoothing parameter (lambda): ', lambda_good
  WRITE(*,'(A)') 'This creates a small, well-conditioned system matrix'
  WRITE(*,*)
  
  ! Test 1: Reference solution with UMFPACK
  WRITE(*,'(A)') '1. Computing reference solution with UMFPACK...'
  sparse_solve_method = SOLVER_UMFPACK
  CALL splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                    a_ref, b_ref, c_ref, d_ref, m, test_function)
  WRITE(*,'(A)') '   [SUCCESS] UMFPACK computed reference solution'
  WRITE(*,*)
  
  ! Test 2: IDR(s) + AMG preconditioning
  WRITE(*,'(A)') '2. Testing IDR(s) + AMG preconditioning...'
  sparse_solve_method = SOLVER_IDRS
  default_iterative_params%preconditioner_type = PRECOND_AMG
  default_iterative_params%verbose = .TRUE.
  default_iterative_params%max_iterations = 200
  default_iterative_params%abs_tolerance = 1.0e-10_DP
  default_iterative_params%rel_tolerance = 1.0e-8_DP
  default_iterative_params%idrs_shadow_space_dim = 4
  
  CALL splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                    a_idrs, b_idrs, c_idrs, d_idrs, m, test_function)
  
  ! Check IDR(s)+AMG accuracy
  error_idrs = MAXVAL(ABS(a_idrs - a_ref)) + &
               MAXVAL(ABS(b_idrs - b_ref)) + &
               MAXVAL(ABS(c_idrs - c_ref)) + &
               MAXVAL(ABS(d_idrs - d_ref))
  idrs_success = error_idrs < 1.0e-6_DP
  
  IF (idrs_success) THEN
    WRITE(*,'(A,ES10.3)') '   [SUCCESS] IDR(s)+AMG error: ', error_idrs
  ELSE
    WRITE(*,'(A,ES10.3)') '   [FAILURE] IDR(s)+AMG error: ', error_idrs
  END IF
  WRITE(*,*)
  
  ! Restore original settings
  sparse_solve_method = saved_method
  default_iterative_params%preconditioner_type = saved_precond
  default_iterative_params%verbose = saved_verbose
  
  ! Final results
  WRITE(*,'(A)') '========================================================='
  WRITE(*,'(A)') 'IDR(s)+AMG Simple Test Results'
  WRITE(*,'(A)') '========================================================='
  IF (idrs_success) THEN
    WRITE(*,'(A)') 'SUCCESS: Native sparse IDR(s)+AMG works on well-conditioned matrix!'
    STOP 0
  ELSE
    WRITE(*,'(A)') 'FAILURE: IDR(s)+AMG failed on well-conditioned matrix'
    STOP 1
  END IF

CONTAINS
  FUNCTION test_function(x, m) RESULT(result)
    REAL(DP), INTENT(in) :: x, m
    REAL(DP) :: result
    result = 1.0_DP  ! Simple test function
  END FUNCTION test_function
  
END PROGRAM test_idrs_amg_simple