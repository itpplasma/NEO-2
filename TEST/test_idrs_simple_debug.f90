PROGRAM test_idrs_simple_debug
  !> Debug IDR(s) with a simple 5x5 matrix

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

  ! Simple small problem first
  INTEGER(I4B), PARAMETER :: n_points = 5  ! This creates small system
  REAL(DP), PARAMETER :: lambda_large = 1.0e-2_DP  ! Large smoothing for well-conditioning
  
  ! Variables
  INTEGER(I4B) :: i
  REAL(DP), ALLOCATABLE :: x(:), y(:), lambda1(:)
  INTEGER(I4B), ALLOCATABLE :: indx(:)
  REAL(DP), ALLOCATABLE :: a_idrs(:), b_idrs(:), c_idrs(:), d_idrs(:)
  REAL(DP), ALLOCATABLE :: a_umf(:), b_umf(:), c_umf(:), d_umf(:)
  REAL(DP) :: c1, cn, m
  INTEGER(I4B) :: sw1, sw2
  
  WRITE(*,'(A)') '========================================='
  WRITE(*,'(A)') 'Debug IDR(s) with Simple 5-Point Problem'
  WRITE(*,'(A)') '========================================='
  
  ! Allocate arrays
  ALLOCATE(x(n_points), y(n_points), lambda1(n_points), indx(n_points))
  ALLOCATE(a_idrs(n_points), b_idrs(n_points), c_idrs(n_points), d_idrs(n_points))
  ALLOCATE(a_umf(n_points), b_umf(n_points), c_umf(n_points), d_umf(n_points))
  
  ! Create simple linear data
  DO i = 1, n_points
    x(i) = REAL(i-1, DP) 
    y(i) = x(i) * 2.0_DP + 1.0_DP  ! Simple linear function y = 2x + 1
    lambda1(i) = lambda_large  ! Well-conditioned
    indx(i) = i
  END DO
  
  ! Boundary conditions
  c1 = 0.0_DP
  cn = 0.0_DP
  sw1 = 2
  sw2 = 4
  m = 0.0_DP
  
  WRITE(*,'(A,I3)') 'Data points: ', n_points
  WRITE(*,'(A,ES10.3)') 'Smoothing parameter (lambda): ', lambda_large
  WRITE(*,*) 'Simple data: y = 2x + 1'
  WRITE(*,*)
  
  ! Test 1: Reference solution with UMFPACK
  WRITE(*,'(A)') '1. Computing reference solution with UMFPACK...'
  sparse_solve_method = SOLVER_UMFPACK
  CALL splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                    a_umf, b_umf, c_umf, d_umf, m, simple_function)
  
  WRITE(*,'(A)') '   [SUCCESS] UMFPACK computed solution'
  WRITE(*,'(A,ES12.5)') '   Max |a| coefficient: ', MAXVAL(ABS(a_umf))
  WRITE(*,'(A,ES12.5)') '   Max |b| coefficient: ', MAXVAL(ABS(b_umf))
  WRITE(*,'(A,ES12.5)') '   Max |c| coefficient: ', MAXVAL(ABS(c_umf))
  WRITE(*,'(A,ES12.5)') '   Max |d| coefficient: ', MAXVAL(ABS(d_umf))
  WRITE(*,*)
  
  ! Test 2: IDR(s) solution
  WRITE(*,'(A)') '2. Testing IDR(s)...'
  
  ! Configure IDR(s) with verbose output to see what happens
  default_iterative_params%preconditioner_type = PRECOND_NONE
  default_iterative_params%idrs_shadow_space_dim = 4
  default_iterative_params%max_iterations = 100
  default_iterative_params%abs_tolerance = 1.0E-10_DP
  default_iterative_params%rel_tolerance = 1.0E-8_DP
  default_iterative_params%verbose = .TRUE.  ! Enable verbose output
  
  sparse_solve_method = SOLVER_IDRS
  CALL splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                    a_idrs, b_idrs, c_idrs, d_idrs, m, simple_function)
  
  WRITE(*,'(A)') '   IDR(s) Results:'
  WRITE(*,'(A,ES12.5)') '   Max |a| coefficient: ', MAXVAL(ABS(a_idrs))
  WRITE(*,'(A,ES12.5)') '   Max |b| coefficient: ', MAXVAL(ABS(b_idrs))
  WRITE(*,'(A,ES12.5)') '   Max |c| coefficient: ', MAXVAL(ABS(c_idrs))
  WRITE(*,'(A,ES12.5)') '   Max |d| coefficient: ', MAXVAL(ABS(d_idrs))
  
  ! Compare solutions
  WRITE(*,*)
  WRITE(*,'(A)') '3. Comparison:'
  WRITE(*,'(A,ES12.5)') '   Max difference in a: ', MAXVAL(ABS(a_idrs - a_umf))
  WRITE(*,'(A,ES12.5)') '   Max difference in b: ', MAXVAL(ABS(b_idrs - b_umf))
  WRITE(*,'(A,ES12.5)') '   Max difference in c: ', MAXVAL(ABS(c_idrs - c_umf))
  WRITE(*,'(A,ES12.5)') '   Max difference in d: ', MAXVAL(ABS(d_idrs - d_umf))
  
  ! Clean up
  DEALLOCATE(x, y, lambda1, indx)
  DEALLOCATE(a_idrs, b_idrs, c_idrs, d_idrs, a_umf, b_umf, c_umf, d_umf)
  
  WRITE(*,'(A)') '========================================='

CONTAINS

  FUNCTION simple_function(x, m) RESULT(result)
    REAL(DP), INTENT(in) :: x, m
    REAL(DP) :: result
    result = 1.0_DP
  END FUNCTION simple_function

END PROGRAM test_idrs_simple_debug