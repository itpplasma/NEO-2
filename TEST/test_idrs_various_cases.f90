PROGRAM test_idrs_various_cases
  !> Test IDR(s) on various problem sizes and conditioning

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

  ! Test parameters
  INTEGER(I4B) :: test_id, i, n_points
  REAL(DP) :: lambda_val, total_error
  REAL(DP), ALLOCATABLE :: x(:), y(:), lambda1(:)
  INTEGER(I4B), ALLOCATABLE :: indx(:)
  REAL(DP), ALLOCATABLE :: a_idrs(:), b_idrs(:), c_idrs(:), d_idrs(:)
  REAL(DP), ALLOCATABLE :: a_umf(:), b_umf(:), c_umf(:), d_umf(:)
  REAL(DP) :: c1, cn, m
  INTEGER(I4B) :: sw1, sw2
  
  WRITE(*,'(A)') '============================================='
  WRITE(*,'(A)') 'Testing IDR(s) on Various Problem Cases'
  WRITE(*,'(A)') '============================================='
  WRITE(*,*)
  
  ! Configure IDR(s) once
  default_iterative_params%preconditioner_type = PRECOND_NONE
  default_iterative_params%idrs_shadow_space_dim = 4
  default_iterative_params%max_iterations = 500
  default_iterative_params%abs_tolerance = 1.0E-10_DP
  default_iterative_params%rel_tolerance = 1.0E-8_DP
  default_iterative_params%verbose = .FALSE.
  
  ! Test 1: Small well-conditioned problem
  test_id = 1
  n_points = 10
  lambda_val = 1.0e-2_DP
  
  WRITE(*,'(A,I0)') 'Test ', test_id
  WRITE(*,'(A)') 'Small well-conditioned spline (10 points, lambda=1e-2)'
  
  ALLOCATE(x(n_points), y(n_points), lambda1(n_points), indx(n_points))
  ALLOCATE(a_idrs(n_points), b_idrs(n_points), c_idrs(n_points), d_idrs(n_points))
  ALLOCATE(a_umf(n_points), b_umf(n_points), c_umf(n_points), d_umf(n_points))
  
  DO i = 1, n_points
    x(i) = REAL(i-1, DP)
    y(i) = SIN(x(i))
    lambda1(i) = lambda_val
    indx(i) = i
  END DO
  
  c1 = 0.0_DP; cn = 0.0_DP; sw1 = 2; sw2 = 4; m = 0.0_DP
  
  ! UMFPACK reference
  sparse_solve_method = SOLVER_UMFPACK
  CALL splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                    a_umf, b_umf, c_umf, d_umf, m, simple_function)
  
  ! IDR(s) test
  sparse_solve_method = SOLVER_IDRS
  CALL splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                    a_idrs, b_idrs, c_idrs, d_idrs, m, simple_function)
  
  total_error = MAXVAL(ABS(a_idrs - a_umf)) + MAXVAL(ABS(b_idrs - b_umf)) + &
                MAXVAL(ABS(c_idrs - c_umf)) + MAXVAL(ABS(d_idrs - d_umf))
  
  IF (total_error < 1.0e-8_DP) THEN
    WRITE(*,'(A,ES9.2)') '  [SUCCESS] Error vs UMFPACK: ', total_error
  ELSE
    WRITE(*,'(A,ES9.2)') '  [FAILURE] Error vs UMFPACK: ', total_error
  END IF
  
  DEALLOCATE(x, y, lambda1, indx, a_idrs, b_idrs, c_idrs, d_idrs, a_umf, b_umf, c_umf, d_umf)
  WRITE(*,*)
  
  ! Test 2: Medium moderately-conditioned problem
  test_id = 2
  n_points = 25
  lambda_val = 1.0e-4_DP
  
  WRITE(*,'(A,I0)') 'Test ', test_id
  WRITE(*,'(A)') 'Medium moderately-conditioned spline (25 points, lambda=1e-4)'
  
  ALLOCATE(x(n_points), y(n_points), lambda1(n_points), indx(n_points))
  ALLOCATE(a_idrs(n_points), b_idrs(n_points), c_idrs(n_points), d_idrs(n_points))
  ALLOCATE(a_umf(n_points), b_umf(n_points), c_umf(n_points), d_umf(n_points))
  
  DO i = 1, n_points
    x(i) = REAL(i-1, DP) / REAL(n_points-1, DP) * 3.14_DP
    y(i) = SIN(2.0_DP * x(i))
    lambda1(i) = lambda_val
    indx(i) = i
  END DO
  
  c1 = 0.0_DP; cn = 0.0_DP; sw1 = 2; sw2 = 4; m = 0.0_DP
  
  ! UMFPACK reference
  sparse_solve_method = SOLVER_UMFPACK
  CALL splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                    a_umf, b_umf, c_umf, d_umf, m, simple_function)
  
  ! IDR(s) test
  sparse_solve_method = SOLVER_IDRS
  CALL splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                    a_idrs, b_idrs, c_idrs, d_idrs, m, simple_function)
  
  total_error = MAXVAL(ABS(a_idrs - a_umf)) + MAXVAL(ABS(b_idrs - b_umf)) + &
                MAXVAL(ABS(c_idrs - c_umf)) + MAXVAL(ABS(d_idrs - d_umf))
  
  IF (total_error < 1.0e-6_DP) THEN
    WRITE(*,'(A,ES9.2)') '  [SUCCESS] Error vs UMFPACK: ', total_error
  ELSE
    WRITE(*,'(A,ES9.2)') '  [FAILURE] Error vs UMFPACK: ', total_error
  END IF
  
  DEALLOCATE(x, y, lambda1, indx, a_idrs, b_idrs, c_idrs, d_idrs, a_umf, b_umf, c_umf, d_umf)
  WRITE(*,*)
  
  ! Test 3: Larger poorly-conditioned problem
  test_id = 3
  n_points = 40
  lambda_val = 1.0e-5_DP
  
  WRITE(*,'(A,I0)') 'Test ', test_id
  WRITE(*,'(A)') 'Larger poorly-conditioned spline (40 points, lambda=1e-5)'
  
  ALLOCATE(x(n_points), y(n_points), lambda1(n_points), indx(n_points))
  ALLOCATE(a_idrs(n_points), b_idrs(n_points), c_idrs(n_points), d_idrs(n_points))
  ALLOCATE(a_umf(n_points), b_umf(n_points), c_umf(n_points), d_umf(n_points))
  
  DO i = 1, n_points
    x(i) = REAL(i-1, DP) / REAL(n_points-1, DP) * 6.28_DP
    y(i) = SIN(x(i)) + 0.1_DP * SIN(10.0_DP * x(i))  ! More oscillatory
    lambda1(i) = lambda_val
    indx(i) = i
  END DO
  
  c1 = 0.0_DP; cn = 0.0_DP; sw1 = 2; sw2 = 4; m = 0.0_DP
  
  ! UMFPACK reference
  sparse_solve_method = SOLVER_UMFPACK
  CALL splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                    a_umf, b_umf, c_umf, d_umf, m, simple_function)
  
  ! IDR(s) test
  sparse_solve_method = SOLVER_IDRS
  CALL splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                    a_idrs, b_idrs, c_idrs, d_idrs, m, simple_function)
  
  total_error = MAXVAL(ABS(a_idrs - a_umf)) + MAXVAL(ABS(b_idrs - b_umf)) + &
                MAXVAL(ABS(c_idrs - c_umf)) + MAXVAL(ABS(d_idrs - d_umf))
  
  IF (total_error < 1.0e-4_DP) THEN
    WRITE(*,'(A,ES9.2)') '  [SUCCESS] Error vs UMFPACK: ', total_error
  ELSE
    WRITE(*,'(A,ES9.2)') '  [FAILURE] Error vs UMFPACK: ', total_error
  END IF
  
  DEALLOCATE(x, y, lambda1, indx, a_idrs, b_idrs, c_idrs, d_idrs, a_umf, b_umf, c_umf, d_umf)
  WRITE(*,*)
  
  WRITE(*,'(A)') '============================================='

CONTAINS

  FUNCTION simple_function(x, m) RESULT(result)
    REAL(DP), INTENT(in) :: x, m
    REAL(DP) :: result
    result = 1.0_DP
  END FUNCTION simple_function

END PROGRAM test_idrs_various_cases