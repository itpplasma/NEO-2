PROGRAM test_amg_realistic
  ! Test AMG on realistic, well-conditioned matrices first
  ! Start simple and work up to more complex cases
  
  USE nrtype, ONLY: I4B, DP
  USE sparse_solvers_mod, ONLY: sparse_solve_method, SOLVER_BICGSTAB, SOLVER_GMRES, &
                                SOLVER_UMFPACK, PRECOND_NONE, PRECOND_ILU, PRECOND_AMG, &
                                default_iterative_params, bicgstab_max_iter, &
                                bicgstab_abs_tolerance, bicgstab_rel_tolerance, &
                                bicgstab_verbose, ilu_fill_level
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

  LOGICAL :: all_tests_passed
  INTEGER :: test_count, passed_count
  
  ! Saved settings
  INTEGER :: saved_method, saved_max_iter, saved_fill_level, saved_precond
  REAL(DP) :: saved_abs_tol, saved_rel_tol
  LOGICAL :: saved_verbose
  
  WRITE(*,'(A)') '========================================================='
  WRITE(*,'(A)') 'AMG Realistic Tests - Progressive Difficulty'
  WRITE(*,'(A)') 'Testing AMG on well-conditioned problems first'
  WRITE(*,'(A)') '========================================================='
  WRITE(*,*)
  
  ! Save current settings
  saved_method = sparse_solve_method
  saved_max_iter = bicgstab_max_iter
  saved_fill_level = ilu_fill_level
  saved_abs_tol = bicgstab_abs_tolerance
  saved_rel_tol = bicgstab_rel_tolerance
  saved_verbose = bicgstab_verbose
  saved_precond = default_iterative_params%preconditioner_type
  
  all_tests_passed = .TRUE.
  test_count = 0
  passed_count = 0
  
  ! Test 1: Small well-conditioned spline (easy case)
  WRITE(*,'(A)') '=== Test 1: Small Well-Conditioned Spline ==='
  test_count = test_count + 1
  IF (test_small_spline()) THEN
    WRITE(*,'(A)') 'PASS: Small spline test succeeded'
    passed_count = passed_count + 1
  ELSE
    WRITE(*,'(A)') 'FAIL: Small spline test failed'
    all_tests_passed = .FALSE.
  END IF
  WRITE(*,*)
  
  ! Test 2: Medium spline with moderate conditioning
  WRITE(*,'(A)') '=== Test 2: Medium Spline (Moderate Lambda) ==='
  test_count = test_count + 1
  IF (test_medium_spline()) THEN
    WRITE(*,'(A)') 'PASS: Medium spline test succeeded'
    passed_count = passed_count + 1
  ELSE
    WRITE(*,'(A)') 'FAIL: Medium spline test failed'
    all_tests_passed = .FALSE.
  END IF
  WRITE(*,*)
  
  ! Test 3: Larger spline with smaller lambda (more challenging)
  WRITE(*,'(A)') '=== Test 3: Larger Spline (Smaller Lambda) ==='
  test_count = test_count + 1
  IF (test_large_spline()) THEN
    WRITE(*,'(A)') 'PASS: Large spline test succeeded'
    passed_count = passed_count + 1
  ELSE
    WRITE(*,'(A)') 'FAIL: Large spline test failed'
    all_tests_passed = .FALSE.
  END IF
  WRITE(*,*)
  
  ! Test 4: Compare AMG vs ILU performance
  WRITE(*,'(A)') '=== Test 4: AMG vs ILU Performance Comparison ==='
  test_count = test_count + 1
  IF (test_amg_vs_ilu()) THEN
    WRITE(*,'(A)') 'PASS: AMG vs ILU comparison succeeded'
    passed_count = passed_count + 1
  ELSE
    WRITE(*,'(A)') 'FAIL: AMG vs ILU comparison failed'
    all_tests_passed = .FALSE.
  END IF
  WRITE(*,*)
  
  ! Restore settings
  sparse_solve_method = saved_method
  bicgstab_max_iter = saved_max_iter
  ilu_fill_level = saved_fill_level
  bicgstab_abs_tolerance = saved_abs_tol
  bicgstab_rel_tolerance = saved_rel_tol
  bicgstab_verbose = saved_verbose
  default_iterative_params%preconditioner_type = saved_precond
  
  ! Final results
  WRITE(*,'(A)') '========================================================='
  WRITE(*,'(A,I0,A,I0,A)') 'Results: ', passed_count, ' of ', test_count, ' tests passed'
  IF (all_tests_passed) THEN
    WRITE(*,'(A)') 'SUCCESS: All AMG realistic tests passed'
    STOP 0
  ELSE
    WRITE(*,'(A)') 'FAILURE: Some AMG realistic tests failed'
    STOP 1
  END IF

CONTAINS

  FUNCTION test_small_spline() RESULT(success)
    LOGICAL :: success
    INTEGER, PARAMETER :: n_points = 10  ! Small test case
    INTEGER, PARAMETER :: n_intervals = 8
    REAL(DP), PARAMETER :: lambda_moderate = 1.0e-3_DP  ! Well-conditioned
    
    REAL(DP), ALLOCATABLE :: x(:), y(:), lambda1(:)
    INTEGER(I4B), ALLOCATABLE :: indx(:)
    REAL(DP), ALLOCATABLE :: a_ref(:), b_ref(:), c_ref(:), d_ref(:)
    REAL(DP), ALLOCATABLE :: a_amg(:), b_amg(:), c_amg(:), d_amg(:)
    REAL(DP) :: c1, cn, m
    INTEGER(I4B) :: sw1, sw2
    REAL(DP) :: error_amg
    INTEGER :: i
    
    WRITE(*,'(A,I0,A)') '  Problem size: ', n_points, ' data points'
    WRITE(*,'(A,ES10.3)') '  Lambda (smoothing): ', lambda_moderate
    
    ! Allocate arrays
    ALLOCATE(x(n_points), y(n_points), lambda1(n_intervals), indx(n_intervals))
    ALLOCATE(a_ref(n_intervals), b_ref(n_intervals), c_ref(n_intervals), d_ref(n_intervals))
    ALLOCATE(a_amg(n_intervals), b_amg(n_intervals), c_amg(n_intervals), d_amg(n_intervals))
    
    ! Create simple, smooth test data
    DO i = 1, n_points
      x(i) = REAL(i-1, DP) / REAL(n_points-1, DP) * 2.0_DP
      y(i) = SIN(x(i))  ! Simple smooth function
    END DO
    
    ! Create intervals and lambda
    DO i = 1, n_intervals
      indx(i) = 1 + (i-1) * (n_points-1) / (n_intervals-1)
      lambda1(i) = lambda_moderate
    END DO
    indx(n_intervals) = n_points
    
    ! Boundary conditions
    c1 = 0.0_DP
    cn = 0.0_DP
    sw1 = 2  ! Natural spline
    sw2 = 4  ! Natural spline
    m = 0.0_DP
    
    ! Get reference solution with UMFPACK
    sparse_solve_method = SOLVER_UMFPACK
    CALL splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                      a_ref, b_ref, c_ref, d_ref, m, test_function)
    
    ! Test BiCGSTAB + AMG
    sparse_solve_method = SOLVER_BICGSTAB
    default_iterative_params%preconditioner_type = PRECOND_AMG
    bicgstab_verbose = .TRUE.   ! Enable debug output
    bicgstab_max_iter = 200
    bicgstab_abs_tolerance = 1.0e-10_DP
    bicgstab_rel_tolerance = 1.0e-8_DP
    
    CALL splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                      a_amg, b_amg, c_amg, d_amg, m, test_function)
    
    ! Check accuracy
    error_amg = MAXVAL(ABS(a_amg - a_ref)) + &
                MAXVAL(ABS(b_amg - b_ref)) + &
                MAXVAL(ABS(c_amg - c_ref)) + &
                MAXVAL(ABS(d_amg - d_ref))
    
    success = error_amg < 1.0e-6_DP
    
    ! Check for NaN/Inf
    success = success .AND. ALL(a_amg == a_amg) .AND. ALL(b_amg == b_amg) .AND. &
              ALL(c_amg == c_amg) .AND. ALL(d_amg == d_amg)
    success = success .AND. ALL(ABS(a_amg) < HUGE(1.0_DP)) .AND. ALL(ABS(b_amg) < HUGE(1.0_DP)) .AND. &
              ALL(ABS(c_amg) < HUGE(1.0_DP)) .AND. ALL(ABS(d_amg) < HUGE(1.0_DP))
    
    WRITE(*,'(A,ES10.3)') '  AMG Error: ', error_amg
    
    DEALLOCATE(x, y, lambda1, indx)
    DEALLOCATE(a_ref, b_ref, c_ref, d_ref)
    DEALLOCATE(a_amg, b_amg, c_amg, d_amg)
    
  END FUNCTION test_small_spline
  
  FUNCTION test_medium_spline() RESULT(success)
    LOGICAL :: success
    INTEGER, PARAMETER :: n_points = 25  ! Medium test case
    INTEGER, PARAMETER :: n_intervals = 20
    REAL(DP), PARAMETER :: lambda_moderate = 1.0e-4_DP  ! Moderately conditioned
    
    REAL(DP), ALLOCATABLE :: x(:), y(:), lambda1(:)
    INTEGER(I4B), ALLOCATABLE :: indx(:)
    REAL(DP), ALLOCATABLE :: a_ref(:), b_ref(:), c_ref(:), d_ref(:)
    REAL(DP), ALLOCATABLE :: a_amg(:), b_amg(:), c_amg(:), d_amg(:)
    REAL(DP) :: c1, cn, m
    INTEGER(I4B) :: sw1, sw2
    REAL(DP) :: error_amg
    INTEGER :: i
    
    WRITE(*,'(A,I0,A)') '  Problem size: ', n_points, ' data points'
    WRITE(*,'(A,ES10.3)') '  Lambda (smoothing): ', lambda_moderate
    
    ! Allocate arrays
    ALLOCATE(x(n_points), y(n_points), lambda1(n_intervals), indx(n_intervals))
    ALLOCATE(a_ref(n_intervals), b_ref(n_intervals), c_ref(n_intervals), d_ref(n_intervals))
    ALLOCATE(a_amg(n_intervals), b_amg(n_intervals), c_amg(n_intervals), d_amg(n_intervals))
    
    ! Create test data with some complexity
    DO i = 1, n_points
      x(i) = REAL(i-1, DP) / REAL(n_points-1, DP) * 4.0_DP
      y(i) = SIN(x(i)) + 0.3_DP * SIN(3.0_DP * x(i))  ! Slightly more complex
    END DO
    
    ! Create intervals and lambda
    DO i = 1, n_intervals
      indx(i) = 1 + (i-1) * (n_points-1) / (n_intervals-1)
      lambda1(i) = lambda_moderate
    END DO
    indx(n_intervals) = n_points
    
    ! Boundary conditions
    c1 = 0.0_DP
    cn = 0.0_DP
    sw1 = 2
    sw2 = 4
    m = 0.0_DP
    
    ! Get reference solution
    sparse_solve_method = SOLVER_UMFPACK
    CALL splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                      a_ref, b_ref, c_ref, d_ref, m, test_function)
    
    ! Test BiCGSTAB + AMG
    sparse_solve_method = SOLVER_BICGSTAB
    default_iterative_params%preconditioner_type = PRECOND_AMG
    bicgstab_verbose = .FALSE.
    bicgstab_max_iter = 300
    bicgstab_abs_tolerance = 1.0e-10_DP
    bicgstab_rel_tolerance = 1.0e-8_DP
    
    CALL splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                      a_amg, b_amg, c_amg, d_amg, m, test_function)
    
    ! Check accuracy
    error_amg = MAXVAL(ABS(a_amg - a_ref)) + &
                MAXVAL(ABS(b_amg - b_ref)) + &
                MAXVAL(ABS(c_amg - c_ref)) + &
                MAXVAL(ABS(d_amg - d_ref))
    
    success = error_amg < 1.0e-5_DP  ! Slightly relaxed tolerance
    
    ! Check for NaN/Inf
    success = success .AND. ALL(a_amg == a_amg) .AND. ALL(b_amg == b_amg) .AND. &
              ALL(c_amg == c_amg) .AND. ALL(d_amg == d_amg)
    success = success .AND. ALL(ABS(a_amg) < HUGE(1.0_DP)) .AND. ALL(ABS(b_amg) < HUGE(1.0_DP)) .AND. &
              ALL(ABS(c_amg) < HUGE(1.0_DP)) .AND. ALL(ABS(d_amg) < HUGE(1.0_DP))
    
    WRITE(*,'(A,ES10.3)') '  AMG Error: ', error_amg
    
    DEALLOCATE(x, y, lambda1, indx)
    DEALLOCATE(a_ref, b_ref, c_ref, d_ref)
    DEALLOCATE(a_amg, b_amg, c_amg, d_amg)
    
  END FUNCTION test_medium_spline
  
  FUNCTION test_large_spline() RESULT(success)
    LOGICAL :: success
    INTEGER, PARAMETER :: n_points = 50  ! Larger test case
    INTEGER, PARAMETER :: n_intervals = 40
    REAL(DP), PARAMETER :: lambda_small = 1.0e-5_DP  ! Smaller lambda, more challenging
    
    REAL(DP), ALLOCATABLE :: x(:), y(:), lambda1(:)
    INTEGER(I4B), ALLOCATABLE :: indx(:)
    REAL(DP), ALLOCATABLE :: a_ref(:), b_ref(:), c_ref(:), d_ref(:)
    REAL(DP), ALLOCATABLE :: a_amg(:), b_amg(:), c_amg(:), d_amg(:)
    REAL(DP) :: c1, cn, m
    INTEGER(I4B) :: sw1, sw2
    REAL(DP) :: error_amg
    INTEGER :: i
    
    WRITE(*,'(A,I0,A)') '  Problem size: ', n_points, ' data points'
    WRITE(*,'(A,ES10.3)') '  Lambda (smoothing): ', lambda_small
    
    ! Allocate arrays
    ALLOCATE(x(n_points), y(n_points), lambda1(n_intervals), indx(n_intervals))
    ALLOCATE(a_ref(n_intervals), b_ref(n_intervals), c_ref(n_intervals), d_ref(n_intervals))
    ALLOCATE(a_amg(n_intervals), b_amg(n_intervals), c_amg(n_intervals), d_amg(n_intervals))
    
    ! Create test data
    DO i = 1, n_points
      x(i) = REAL(i-1, DP) / REAL(n_points-1, DP) * 6.0_DP
      y(i) = EXP(-x(i)) * SIN(2.0_DP * x(i))  ! Exponentially decaying oscillation
    END DO
    
    ! Create intervals and lambda
    DO i = 1, n_intervals
      indx(i) = 1 + (i-1) * (n_points-1) / (n_intervals-1)
      lambda1(i) = lambda_small
    END DO
    indx(n_intervals) = n_points
    
    ! Boundary conditions
    c1 = 0.0_DP
    cn = 0.0_DP
    sw1 = 2
    sw2 = 4
    m = 0.0_DP
    
    ! Get reference solution
    sparse_solve_method = SOLVER_UMFPACK
    CALL splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                      a_ref, b_ref, c_ref, d_ref, m, test_function)
    
    ! Test BiCGSTAB + AMG
    sparse_solve_method = SOLVER_BICGSTAB
    default_iterative_params%preconditioner_type = PRECOND_AMG
    bicgstab_verbose = .FALSE.
    bicgstab_max_iter = 500
    bicgstab_abs_tolerance = 1.0e-9_DP
    bicgstab_rel_tolerance = 1.0e-7_DP
    
    CALL splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                      a_amg, b_amg, c_amg, d_amg, m, test_function)
    
    ! Check accuracy (more relaxed for challenging case)
    error_amg = MAXVAL(ABS(a_amg - a_ref)) + &
                MAXVAL(ABS(b_amg - b_ref)) + &
                MAXVAL(ABS(c_amg - c_ref)) + &
                MAXVAL(ABS(d_amg - d_ref))
    
    success = error_amg < 1.0e-4_DP  ! Relaxed tolerance for challenging case
    
    ! Check for NaN/Inf
    success = success .AND. ALL(a_amg == a_amg) .AND. ALL(b_amg == b_amg) .AND. &
              ALL(c_amg == c_amg) .AND. ALL(d_amg == d_amg)
    success = success .AND. ALL(ABS(a_amg) < HUGE(1.0_DP)) .AND. ALL(ABS(b_amg) < HUGE(1.0_DP)) .AND. &
              ALL(ABS(c_amg) < HUGE(1.0_DP)) .AND. ALL(ABS(d_amg) < HUGE(1.0_DP))
    
    WRITE(*,'(A,ES10.3)') '  AMG Error: ', error_amg
    
    DEALLOCATE(x, y, lambda1, indx)
    DEALLOCATE(a_ref, b_ref, c_ref, d_ref)
    DEALLOCATE(a_amg, b_amg, c_amg, d_amg)
    
  END FUNCTION test_large_spline
  
  FUNCTION test_amg_vs_ilu() RESULT(success)
    LOGICAL :: success
    INTEGER, PARAMETER :: n_points = 30
    INTEGER, PARAMETER :: n_intervals = 25
    REAL(DP), PARAMETER :: lambda_test = 1.0e-4_DP
    
    REAL(DP), ALLOCATABLE :: x(:), y(:), lambda1(:)
    INTEGER(I4B), ALLOCATABLE :: indx(:)
    REAL(DP), ALLOCATABLE :: a_ref(:), b_ref(:), c_ref(:), d_ref(:)
    REAL(DP), ALLOCATABLE :: a_amg(:), b_amg(:), c_amg(:), d_amg(:)
    REAL(DP), ALLOCATABLE :: a_ilu(:), b_ilu(:), c_ilu(:), d_ilu(:)
    REAL(DP) :: c1, cn, m
    INTEGER(I4B) :: sw1, sw2
    REAL(DP) :: error_amg, error_ilu
    INTEGER :: i
    
    WRITE(*,'(A,I0,A)') '  Problem size: ', n_points, ' data points'
    WRITE(*,'(A,ES10.3)') '  Lambda (smoothing): ', lambda_test
    
    ! Allocate arrays
    ALLOCATE(x(n_points), y(n_points), lambda1(n_intervals), indx(n_intervals))
    ALLOCATE(a_ref(n_intervals), b_ref(n_intervals), c_ref(n_intervals), d_ref(n_intervals))
    ALLOCATE(a_amg(n_intervals), b_amg(n_intervals), c_amg(n_intervals), d_amg(n_intervals))
    ALLOCATE(a_ilu(n_intervals), b_ilu(n_intervals), c_ilu(n_intervals), d_ilu(n_intervals))
    
    ! Create test data
    DO i = 1, n_points
      x(i) = REAL(i-1, DP) / REAL(n_points-1, DP) * 3.0_DP
      y(i) = SIN(x(i)) * EXP(-0.2_DP * x(i))
    END DO
    
    ! Create intervals and lambda
    DO i = 1, n_intervals
      indx(i) = 1 + (i-1) * (n_points-1) / (n_intervals-1)
      lambda1(i) = lambda_test
    END DO
    indx(n_intervals) = n_points
    
    ! Boundary conditions
    c1 = 0.0_DP
    cn = 0.0_DP
    sw1 = 2
    sw2 = 4
    m = 0.0_DP
    
    ! Get reference solution
    sparse_solve_method = SOLVER_UMFPACK
    CALL splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                      a_ref, b_ref, c_ref, d_ref, m, test_function)
    
    ! Test BiCGSTAB + AMG
    sparse_solve_method = SOLVER_BICGSTAB
    default_iterative_params%preconditioner_type = PRECOND_AMG
    bicgstab_verbose = .FALSE.
    bicgstab_max_iter = 300
    bicgstab_abs_tolerance = 1.0e-10_DP
    bicgstab_rel_tolerance = 1.0e-8_DP
    
    CALL splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                      a_amg, b_amg, c_amg, d_amg, m, test_function)
    
    ! Test BiCGSTAB + ILU for comparison
    default_iterative_params%preconditioner_type = PRECOND_ILU
    ilu_fill_level = 1
    
    CALL splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                      a_ilu, b_ilu, c_ilu, d_ilu, m, test_function)
    
    ! Calculate errors
    error_amg = MAXVAL(ABS(a_amg - a_ref)) + &
                MAXVAL(ABS(b_amg - b_ref)) + &
                MAXVAL(ABS(c_amg - c_ref)) + &
                MAXVAL(ABS(d_amg - d_ref))
    
    error_ilu = MAXVAL(ABS(a_ilu - a_ref)) + &
                MAXVAL(ABS(b_ilu - b_ref)) + &
                MAXVAL(ABS(c_ilu - c_ref)) + &
                MAXVAL(ABS(d_ilu - d_ref))
    
    WRITE(*,'(A,ES10.3)') '  AMG Error: ', error_amg
    WRITE(*,'(A,ES10.3)') '  ILU Error: ', error_ilu
    
    ! Both should be accurate, AMG ideally better or comparable
    success = (error_amg < 1.0e-5_DP) .AND. (error_ilu < 1.0e-5_DP)
    
    ! Check for NaN/Inf in AMG solution
    success = success .AND. ALL(a_amg == a_amg) .AND. ALL(b_amg == b_amg) .AND. &
              ALL(c_amg == c_amg) .AND. ALL(d_amg == d_amg)
    
    IF (error_amg < error_ilu) THEN
      WRITE(*,'(A)') '  AMG outperformed ILU!'
    ELSE IF (error_amg < 2.0_DP * error_ilu) THEN
      WRITE(*,'(A)') '  AMG performance comparable to ILU'
    ELSE
      WRITE(*,'(A)') '  ILU performed better than AMG'
    END IF
    
    DEALLOCATE(x, y, lambda1, indx)
    DEALLOCATE(a_ref, b_ref, c_ref, d_ref)
    DEALLOCATE(a_amg, b_amg, c_amg, d_amg)
    DEALLOCATE(a_ilu, b_ilu, c_ilu, d_ilu)
    
  END FUNCTION test_amg_vs_ilu

  FUNCTION test_function(x, m) RESULT(result)
    REAL(DP), INTENT(in) :: x, m
    REAL(DP) :: result
    result = 1.0_DP  ! Simple test function
  END FUNCTION test_function
  
END PROGRAM test_amg_realistic