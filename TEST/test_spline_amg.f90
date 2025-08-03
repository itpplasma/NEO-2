PROGRAM test_spline_amg
  ! Test AMG preconditioner on the exact 404x404 ill-conditioned spline matrix
  ! Uses the same proven parameters as test_gmres_spline_404.f90
  
  USE nrtype, ONLY: I4B, DP
  USE sparse_solvers_mod, ONLY: sparse_solve_method, &
                                SOLVER_UMFPACK, SOLVER_BICGSTAB, SOLVER_GMRES, &
                                PRECOND_NONE, PRECOND_ILU, PRECOND_AMG, &
                                ilu_fill_level, bicgstab_max_iter, &
                                bicgstab_abs_tolerance, bicgstab_rel_tolerance, &
                                bicgstab_verbose, default_iterative_params
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

  ! Test parameters - recreate the exact 404x404 matrix from test_gmres_spline_404
  INTEGER(I4B), PARAMETER :: n_points = 67  ! This creates ~404x404 system  
  REAL(DP), PARAMETER :: lambda_small = 1.0e-6_DP  ! Small smoothing creates ill-conditioning
  
  ! Variables
  INTEGER(I4B) :: i, info
  REAL(DP), ALLOCATABLE :: x(:), y(:), lambda1(:)
  INTEGER(I4B), ALLOCATABLE :: indx(:)
  REAL(DP), ALLOCATABLE :: a_ref(:), b_ref(:), c_ref(:), d_ref(:)  ! Reference (UMFPACK)
  REAL(DP), ALLOCATABLE :: a_amg_bicgstab(:), b_amg_bicgstab(:), c_amg_bicgstab(:), d_amg_bicgstab(:)  ! BiCGSTAB+AMG
  REAL(DP), ALLOCATABLE :: a_amg_gmres(:), b_amg_gmres(:), c_amg_gmres(:), d_amg_gmres(:)  ! GMRES+AMG
  REAL(DP) :: c1, cn, m
  INTEGER(I4B) :: sw1, sw2
  REAL(DP) :: error_bicgstab_amg, error_gmres_amg
  LOGICAL :: bicgstab_amg_success, gmres_amg_success
  
  ! Saved settings
  INTEGER :: saved_method, saved_max_iter, saved_fill_level, saved_precond
  REAL(DP) :: saved_abs_tol, saved_rel_tol
  LOGICAL :: saved_verbose
  
  WRITE(*,'(A)') '========================================================='
  WRITE(*,'(A)') 'AMG Preconditioner Test on Critical 404x404 Spline Matrix' 
  WRITE(*,'(A)') 'Testing BiCGSTAB+AMG and GMRES+AMG on pathological case'
  WRITE(*,'(A)') '========================================================='
  WRITE(*,*)
  
  ! Save original solver settings
  saved_method = sparse_solve_method
  saved_max_iter = bicgstab_max_iter
  saved_fill_level = ilu_fill_level
  saved_abs_tol = bicgstab_abs_tolerance
  saved_rel_tol = bicgstab_rel_tolerance
  saved_verbose = bicgstab_verbose
  saved_precond = default_iterative_params%preconditioner_type
  
  ! Allocate arrays - using the exact same structure as test_gmres_spline_404
  ALLOCATE(x(n_points), y(n_points), lambda1(n_points), indx(n_points))
  ALLOCATE(a_ref(n_points), b_ref(n_points), c_ref(n_points), d_ref(n_points))
  ALLOCATE(a_amg_bicgstab(n_points), b_amg_bicgstab(n_points), c_amg_bicgstab(n_points), d_amg_bicgstab(n_points))
  ALLOCATE(a_amg_gmres(n_points), b_amg_gmres(n_points), c_amg_gmres(n_points), d_amg_gmres(n_points))
  
  ! Create the exact same problematic spline data as test_gmres_spline_404
  DO i = 1, n_points
    x(i) = REAL(i-1, DP) / REAL(n_points-1, DP) * 10.0_DP
    ! Complex oscillatory function that creates ill-conditioning
    y(i) = SIN(2.0_DP * 3.14159_DP * x(i)) + 0.1_DP * SIN(25.0_DP * x(i))
    lambda1(i) = lambda_small  ! Very small smoothing - critical for ill-conditioning!
    indx(i) = i
  END DO
  
  ! Boundary conditions - same as test_gmres_spline_404
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
  
  ! Test 2: BiCGSTAB + AMG preconditioning
  WRITE(*,'(A)') '2. Testing BiCGSTAB + AMG preconditioning...'
  sparse_solve_method = SOLVER_BICGSTAB
  default_iterative_params%preconditioner_type = PRECOND_AMG
  bicgstab_verbose = .TRUE.
  bicgstab_max_iter = 1000
  bicgstab_abs_tolerance = 1.0e-10_DP
  bicgstab_rel_tolerance = 1.0e-8_DP
  
  CALL splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                    a_amg_bicgstab, b_amg_bicgstab, c_amg_bicgstab, d_amg_bicgstab, m, test_function)
  
  ! Check BiCGSTAB+AMG accuracy
  error_bicgstab_amg = MAXVAL(ABS(a_amg_bicgstab - a_ref)) + &
                       MAXVAL(ABS(b_amg_bicgstab - b_ref)) + &
                       MAXVAL(ABS(c_amg_bicgstab - c_ref)) + &
                       MAXVAL(ABS(d_amg_bicgstab - d_ref))
  bicgstab_amg_success = error_bicgstab_amg < 1.0e-3_DP
  
  IF (bicgstab_amg_success) THEN
    WRITE(*,'(A,ES10.3)') '   [SUCCESS] BiCGSTAB+AMG error: ', error_bicgstab_amg
  ELSE
    WRITE(*,'(A,ES10.3)') '   [FAILURE] BiCGSTAB+AMG error: ', error_bicgstab_amg
  END IF
  WRITE(*,*)
  
  ! Test 3: GMRES + AMG preconditioning  
  WRITE(*,'(A)') '3. Testing GMRES + AMG preconditioning...'
  sparse_solve_method = SOLVER_GMRES
  default_iterative_params%preconditioner_type = PRECOND_AMG
  
  CALL splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                    a_amg_gmres, b_amg_gmres, c_amg_gmres, d_amg_gmres, m, test_function)
  
  ! Check GMRES+AMG accuracy
  error_gmres_amg = MAXVAL(ABS(a_amg_gmres - a_ref)) + &
                    MAXVAL(ABS(b_amg_gmres - b_ref)) + &
                    MAXVAL(ABS(c_amg_gmres - c_ref)) + &
                    MAXVAL(ABS(d_amg_gmres - d_ref))
  gmres_amg_success = error_gmres_amg < 1.0e-3_DP
  
  IF (gmres_amg_success) THEN
    WRITE(*,'(A,ES10.3)') '   [SUCCESS] GMRES+AMG error: ', error_gmres_amg
  ELSE
    WRITE(*,'(A,ES10.3)') '   [FAILURE] GMRES+AMG error: ', error_gmres_amg
  END IF
  WRITE(*,*)
  
  ! Restore original settings
  sparse_solve_method = saved_method
  bicgstab_max_iter = saved_max_iter
  ilu_fill_level = saved_fill_level
  bicgstab_abs_tolerance = saved_abs_tol
  bicgstab_rel_tolerance = saved_rel_tol
  bicgstab_verbose = saved_verbose
  default_iterative_params%preconditioner_type = saved_precond
  
  ! Final results
  WRITE(*,'(A)') '========================================================='
  WRITE(*,'(A)') 'AMG Preconditioner Test Results on 404x404 Matrix'
  WRITE(*,'(A)') '========================================================='
  IF (bicgstab_amg_success .AND. gmres_amg_success) THEN
    WRITE(*,'(A)') 'SUCCESS: Both BiCGSTAB+AMG and GMRES+AMG solved the matrix'
    STOP 0
  ELSE IF (bicgstab_amg_success) THEN
    WRITE(*,'(A)') 'PARTIAL SUCCESS: BiCGSTAB+AMG works, GMRES+AMG failed'
    STOP 0
  ELSE IF (gmres_amg_success) THEN
    WRITE(*,'(A)') 'PARTIAL SUCCESS: GMRES+AMG works, BiCGSTAB+AMG failed'
    STOP 0
  ELSE
    WRITE(*,'(A)') 'FAILURE: Both AMG approaches failed on this matrix'
    STOP 1
  END IF

CONTAINS
  FUNCTION test_function(x, m) RESULT(result)
    REAL(DP), INTENT(in) :: x, m
    REAL(DP) :: result
    result = 1.0_DP  ! Simple test function
  END FUNCTION test_function
  
END PROGRAM test_spline_amg