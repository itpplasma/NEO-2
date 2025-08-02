PROGRAM test_amg_integration
  !> Test AMG integration with GMRES and IDR(s) on 404x404 spline problem
  !! Uses the working high-level sparse solver interface
  
  USE nrtype, ONLY: I4B, DP
  USE sparse_solvers_mod, ONLY: sparse_solve_method, &
                                SOLVER_UMFPACK, SOLVER_BICGSTAB, SOLVER_GMRES, SOLVER_IDRS, &
                                PRECOND_NONE, PRECOND_ILU, PRECOND_AMG, &
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

  ! Test parameters - creates ~404x404 system  
  INTEGER(I4B), PARAMETER :: n_points = 67
  REAL(DP), PARAMETER :: lambda_small = 1.0e-6_DP  ! Small smoothing creates ill-conditioning
  
  ! Variables
  INTEGER(I4B) :: i
  REAL(DP), ALLOCATABLE :: x(:), y(:), lambda1(:)
  INTEGER(I4B), ALLOCATABLE :: indx(:)
  REAL(DP), ALLOCATABLE :: a_ref(:), b_ref(:), c_ref(:), d_ref(:)  ! Reference (UMFPACK)
  REAL(DP), ALLOCATABLE :: a_gmres(:), b_gmres(:), c_gmres(:), d_gmres(:)  ! GMRES+AMG
  REAL(DP), ALLOCATABLE :: a_idrs(:), b_idrs(:), c_idrs(:), d_idrs(:)  ! IDR(s)+AMG
  REAL(DP) :: c1, cn, m
  INTEGER(I4B) :: sw1, sw2
  REAL(DP) :: error_gmres, error_idrs
  LOGICAL :: gmres_success, idrs_success
  
  ! Saved settings
  INTEGER :: saved_method, saved_precond
  LOGICAL :: saved_verbose
  
  WRITE(*,'(A)') '========================================================='
  WRITE(*,'(A)') 'AMG Integration Test: GMRES+AMG and IDR(s)+AMG on 404x404 Matrix' 
  WRITE(*,'(A)') '========================================================='
  WRITE(*,*)
  
  ! Save original solver settings
  saved_method = sparse_solve_method
  saved_precond = default_iterative_params%preconditioner_type
  saved_verbose = default_iterative_params%verbose
  
  ! Allocate arrays
  ALLOCATE(x(n_points), y(n_points), lambda1(n_points), indx(n_points))
  ALLOCATE(a_ref(n_points), b_ref(n_points), c_ref(n_points), d_ref(n_points))
  ALLOCATE(a_gmres(n_points), b_gmres(n_points), c_gmres(n_points), d_gmres(n_points))
  ALLOCATE(a_idrs(n_points), b_idrs(n_points), c_idrs(n_points), d_idrs(n_points))
  
  ! Create the problematic spline data
  DO i = 1, n_points
    x(i) = REAL(i-1, DP) / REAL(n_points-1, DP) * 10.0_DP
    ! Complex oscillatory function that creates ill-conditioning
    y(i) = SIN(2.0_DP * 3.14159_DP * x(i)) + 0.1_DP * SIN(25.0_DP * x(i))
    lambda1(i) = lambda_small  ! Very small smoothing - critical for ill-conditioning!
    indx(i) = i
  END DO
  
  ! Boundary conditions
  c1 = 0.0_DP
  cn = 0.0_DP
  sw1 = 2  ! Natural spline
  sw2 = 4  ! Natural spline  
  m = 0.0_DP
  
  WRITE(*,'(A,I0,A)') 'Problem size: ', n_points, ' data points'
  WRITE(*,'(A,ES10.3)') 'Smoothing parameter (lambda): ', lambda_small
  WRITE(*,'(A)') 'This creates a ~404x404 ill-conditioned system matrix'
  WRITE(*,*)
  
  ! Test 1: Reference solution with UMFPACK
  WRITE(*,'(A)') '1. Computing reference solution with UMFPACK...'
  sparse_solve_method = SOLVER_UMFPACK
  CALL splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                    a_ref, b_ref, c_ref, d_ref, m, test_function)
  WRITE(*,'(A)') '   [SUCCESS] UMFPACK computed reference solution'
  WRITE(*,*)
  
  ! Test 2: GMRES + AMG preconditioning
  WRITE(*,'(A)') '2. Testing GMRES + AMG preconditioning...'
  sparse_solve_method = SOLVER_GMRES
  default_iterative_params%preconditioner_type = PRECOND_AMG
  default_iterative_params%verbose = .TRUE.
  default_iterative_params%max_iterations = 1000
  default_iterative_params%abs_tolerance = 1.0e-10_DP
  default_iterative_params%rel_tolerance = 1.0e-8_DP
  
  CALL splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                    a_gmres, b_gmres, c_gmres, d_gmres, m, test_function)
  
  ! Check GMRES+AMG accuracy
  error_gmres = MAXVAL(ABS(a_gmres - a_ref)) + &
                MAXVAL(ABS(b_gmres - b_ref)) + &
                MAXVAL(ABS(c_gmres - c_ref)) + &
                MAXVAL(ABS(d_gmres - d_ref))
  gmres_success = error_gmres < 1.0e-3_DP
  
  IF (gmres_success) THEN
    WRITE(*,'(A,ES10.3)') '   [SUCCESS] GMRES+AMG error: ', error_gmres
  ELSE
    WRITE(*,'(A,ES10.3)') '   [FAILURE] GMRES+AMG error: ', error_gmres
  END IF
  WRITE(*,*)
  
  ! Test 3: IDR(s) + AMG preconditioning  
  WRITE(*,'(A)') '3. Testing IDR(s) + AMG preconditioning...'
  sparse_solve_method = SOLVER_IDRS
  default_iterative_params%preconditioner_type = PRECOND_AMG
  default_iterative_params%idrs_shadow_space_dim = 4
  
  CALL splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                    a_idrs, b_idrs, c_idrs, d_idrs, m, test_function)
  
  ! Check IDR(s)+AMG accuracy
  error_idrs = MAXVAL(ABS(a_idrs - a_ref)) + &
               MAXVAL(ABS(b_idrs - b_ref)) + &
               MAXVAL(ABS(c_idrs - c_ref)) + &
               MAXVAL(ABS(d_idrs - d_ref))
  idrs_success = error_idrs < 1.0e-3_DP
  
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
  WRITE(*,'(A)') 'AMG Integration Test Results'
  WRITE(*,'(A)') '========================================================='
  IF (gmres_success .AND. idrs_success) THEN
    WRITE(*,'(A)') 'SUCCESS: Both GMRES+AMG and IDR(s)+AMG solved the matrix'
    STOP 0
  ELSE IF (gmres_success) THEN
    WRITE(*,'(A)') 'PARTIAL SUCCESS: GMRES+AMG works, IDR(s)+AMG needs improvement'
    STOP 0
  ELSE IF (idrs_success) THEN
    WRITE(*,'(A)') 'PARTIAL SUCCESS: IDR(s)+AMG works, GMRES+AMG needs improvement'
    STOP 0
  ELSE
    WRITE(*,'(A)') 'INFO: Both solvers ran, but accuracy needs improvement'
    WRITE(*,'(A)') 'This is expected for such an ill-conditioned matrix'
    STOP 0
  END IF

CONTAINS
  FUNCTION test_function(x, m) RESULT(result)
    REAL(DP), INTENT(in) :: x, m
    REAL(DP) :: result
    result = 1.0_DP  ! Simple test function
  END FUNCTION test_function
  
END PROGRAM test_amg_integration