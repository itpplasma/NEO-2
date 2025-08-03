PROGRAM test_idrs_amg_debug
  !> Debug IDR(s)+AMG by direct comparison with GMRES+AMG on same problem

  USE nrtype, ONLY: I4B, DP
  USE sparse_solvers_mod, ONLY: sparse_solve_idrs_real, sparse_solve_gmres_real, &
                                default_iterative_params, PRECOND_AMG, PRECOND_NONE, &
                                sparse_solve_method, SOLVER_UMFPACK, SOLVER_IDRS, SOLVER_GMRES
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

  ! Test on a simple problem where both should work
  INTEGER(I4B), PARAMETER :: n_points = 15  ! Small enough for AMG to work
  REAL(DP), PARAMETER :: lambda_val = 1.0e-4_DP  ! Well-conditioned
  
  ! Variables
  INTEGER(I4B) :: i
  REAL(DP), ALLOCATABLE :: x(:), y(:), lambda1(:)
  INTEGER(I4B), ALLOCATABLE :: indx(:)
  REAL(DP), ALLOCATABLE :: a_umf(:), b_umf(:), c_umf(:), d_umf(:)      ! UMFPACK reference
  REAL(DP), ALLOCATABLE :: a_gmres(:), b_gmres(:), c_gmres(:), d_gmres(:)  ! GMRES+AMG
  REAL(DP), ALLOCATABLE :: a_idrs(:), b_idrs(:), c_idrs(:), d_idrs(:)      ! IDR(s)+AMG
  REAL(DP) :: c1, cn, m
  INTEGER(I4B) :: sw1, sw2
  REAL(DP) :: error_gmres, error_idrs
  
  WRITE(*,'(A)') '========================================='
  WRITE(*,'(A)') 'Debug IDR(s)+AMG vs GMRES+AMG'
  WRITE(*,'(A)') 'Same problem, different iterative solver'
  WRITE(*,'(A)') '========================================='
  WRITE(*,*)
  
  ! Allocate arrays
  ALLOCATE(x(n_points), y(n_points), lambda1(n_points), indx(n_points))
  ALLOCATE(a_umf(n_points), b_umf(n_points), c_umf(n_points), d_umf(n_points))
  ALLOCATE(a_gmres(n_points), b_gmres(n_points), c_gmres(n_points), d_gmres(n_points))
  ALLOCATE(a_idrs(n_points), b_idrs(n_points), c_idrs(n_points), d_idrs(n_points))
  
  ! Create simple well-conditioned data
  DO i = 1, n_points
    x(i) = REAL(i-1, DP) / REAL(n_points-1, DP) * 3.14_DP
    y(i) = SIN(x(i))
    lambda1(i) = lambda_val
    indx(i) = i
  END DO
  
  c1 = 0.0_DP; cn = 0.0_DP; sw1 = 2; sw2 = 4; m = 0.0_DP
  
  WRITE(*,'(A,I0)') 'Problem size: ', n_points, ' points'
  WRITE(*,'(A,ES10.3)') 'Lambda (smoothing): ', lambda_val
  WRITE(*,*)
  
  ! Step 1: Get UMFPACK reference
  WRITE(*,'(A)') '1. Computing UMFPACK reference...'
  sparse_solve_method = SOLVER_UMFPACK
  CALL splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                    a_umf, b_umf, c_umf, d_umf, m, test_function)
  WRITE(*,'(A)') '   [SUCCESS] UMFPACK computed reference'
  WRITE(*,*)
  
  ! Step 2: Test GMRES+AMG
  WRITE(*,'(A)') '2. Testing GMRES+AMG...'
  
  ! Configure GMRES+AMG
  default_iterative_params%preconditioner_type = PRECOND_AMG
  default_iterative_params%max_iterations = 500
  default_iterative_params%abs_tolerance = 1.0E-10_DP
  default_iterative_params%rel_tolerance = 1.0E-8_DP
  default_iterative_params%verbose = .TRUE.
  
  sparse_solve_method = SOLVER_GMRES
  CALL splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                    a_gmres, b_gmres, c_gmres, d_gmres, m, test_function)
  
  ! Check GMRES+AMG results
  IF (ANY(a_gmres /= a_gmres) .OR. ANY(b_gmres /= b_gmres) .OR. &
      ANY(c_gmres /= c_gmres) .OR. ANY(d_gmres /= d_gmres)) THEN
    WRITE(*,'(A)') '   [FAILURE] GMRES+AMG produced NaN'
    error_gmres = 1.0e10_DP
  ELSE
    error_gmres = MAXVAL(ABS(a_gmres - a_umf)) + MAXVAL(ABS(b_gmres - b_umf)) + &
                  MAXVAL(ABS(c_gmres - c_umf)) + MAXVAL(ABS(d_gmres - d_umf))
    WRITE(*,'(A,ES10.3)') '   [SUCCESS] GMRES+AMG error: ', error_gmres
  END IF
  WRITE(*,*)
  
  ! Step 3: Test IDR(s)+AMG with same settings
  WRITE(*,'(A)') '3. Testing IDR(s)+AMG (same AMG settings)...'
  
  ! Keep same AMG settings, just switch to IDR(s)
  default_iterative_params%idrs_shadow_space_dim = 4
  
  sparse_solve_method = SOLVER_IDRS
  CALL splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                    a_idrs, b_idrs, c_idrs, d_idrs, m, test_function)
  
  ! Check IDR(s)+AMG results
  IF (ANY(a_idrs /= a_idrs) .OR. ANY(b_idrs /= b_idrs) .OR. &
      ANY(c_idrs /= c_idrs) .OR. ANY(d_idrs /= d_idrs)) THEN
    WRITE(*,'(A)') '   [FAILURE] IDR(s)+AMG produced NaN - THIS IS THE BUG!'
    error_idrs = 1.0e10_DP
    
    ! Print some diagnostic info
    WRITE(*,'(A)') '   Diagnostic info:'
    WRITE(*,'(A,ES12.5)') '   Max |a_idrs|: ', MAXVAL(ABS(a_idrs), MASK=(a_idrs==a_idrs))
    WRITE(*,'(A,ES12.5)') '   Max |b_idrs|: ', MAXVAL(ABS(b_idrs), MASK=(b_idrs==b_idrs))
    WRITE(*,'(A,ES12.5)') '   Max |c_idrs|: ', MAXVAL(ABS(c_idrs), MASK=(c_idrs==c_idrs))
    WRITE(*,'(A,ES12.5)') '   Max |d_idrs|: ', MAXVAL(ABS(d_idrs), MASK=(d_idrs==d_idrs))
  ELSE
    error_idrs = MAXVAL(ABS(a_idrs - a_umf)) + MAXVAL(ABS(b_idrs - b_umf)) + &
                 MAXVAL(ABS(c_idrs - c_umf)) + MAXVAL(ABS(d_idrs - d_umf))
    WRITE(*,'(A,ES10.3)') '   [SUCCESS] IDR(s)+AMG error: ', error_idrs
  END IF
  WRITE(*,*)
  
  ! Step 4: Test IDR(s) without AMG for comparison
  WRITE(*,'(A)') '4. Testing IDR(s) without preconditioning (for comparison)...'
  
  default_iterative_params%preconditioner_type = PRECOND_NONE
  default_iterative_params%verbose = .FALSE.  ! Reduce output
  
  sparse_solve_method = SOLVER_IDRS
  CALL splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                    a_idrs, b_idrs, c_idrs, d_idrs, m, test_function)
  
  IF (ANY(a_idrs /= a_idrs) .OR. ANY(b_idrs /= b_idrs) .OR. &
      ANY(c_idrs /= c_idrs) .OR. ANY(d_idrs /= d_idrs)) THEN
    WRITE(*,'(A)') '   [FAILURE] IDR(s) without precond produced NaN'
  ELSE
    error_idrs = MAXVAL(ABS(a_idrs - a_umf)) + MAXVAL(ABS(b_idrs - b_umf)) + &
                 MAXVAL(ABS(c_idrs - c_umf)) + MAXVAL(ABS(d_idrs - d_umf))
    WRITE(*,'(A,ES10.3)') '   [SUCCESS] IDR(s) no precond error: ', error_idrs
  END IF
  WRITE(*,*)
  
  ! Summary
  WRITE(*,'(A)') '========================================='
  WRITE(*,'(A)') 'Summary: AMG Compatibility Test'
  WRITE(*,'(A)') '========================================='
  IF (error_gmres < 1.0e-6_DP) THEN
    WRITE(*,'(A)') 'GMRES+AMG: WORKS'
  ELSE
    WRITE(*,'(A)') 'GMRES+AMG: FAILS'
  END IF
  
  IF (error_idrs < 1.0e-6_DP) THEN
    WRITE(*,'(A)') 'IDR(s)+AMG: WORKS'
  ELSE
    WRITE(*,'(A)') 'IDR(s)+AMG: FAILS - NEED TO DEBUG THIS!'
  END IF
  
  ! Clean up
  DEALLOCATE(x, y, lambda1, indx)
  DEALLOCATE(a_umf, b_umf, c_umf, d_umf, a_gmres, b_gmres, c_gmres, d_gmres)
  DEALLOCATE(a_idrs, b_idrs, c_idrs, d_idrs)

CONTAINS

  FUNCTION test_function(x, m) RESULT(result)
    REAL(DP), INTENT(in) :: x, m
    REAL(DP) :: result
    result = 1.0_DP
  END FUNCTION test_function

END PROGRAM test_idrs_amg_debug