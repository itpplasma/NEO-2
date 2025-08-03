PROGRAM test_bicgstab_l_validation
  ! Quick validation test for BiCGSTAB(ℓ) on well-conditioned matrices
  ! to verify the implementation works correctly before testing pathological cases
  
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

  ! Test parameters for well-conditioned matrix
  INTEGER(I4B), PARAMETER :: n_intervals = 10  ! Small, well-conditioned
  
  ! Test variables
  INTEGER(I4B) :: n_points, i, l_val
  REAL(DP), ALLOCATABLE :: x(:), y(:), lambda1(:)
  INTEGER(I4B), ALLOCATABLE :: indx(:)
  REAL(DP), ALLOCATABLE :: a_ref(:), b_ref(:), c_ref(:), d_ref(:)
  REAL(DP), ALLOCATABLE :: a_test(:), b_test(:), c_test(:), d_test(:)
  REAL(DP) :: c1, cn, m
  INTEGER(I4B) :: sw1, sw2
  REAL(DP) :: max_diff_a, max_diff_b, max_diff_c, max_diff_d, total_error
  
  ! Saved settings
  INTEGER :: saved_method, saved_max_iter, saved_fill_level, saved_l_param
  REAL(DP) :: saved_abs_tol, saved_rel_tol
  LOGICAL :: saved_verbose
  
  WRITE(*,'(A)') '==========================================='
  WRITE(*,'(A)') 'BiCGSTAB(ℓ) Validation Test'
  WRITE(*,'(A)') 'Testing on Well-Conditioned Spline Matrix'
  WRITE(*,'(A)') '==========================================='
  WRITE(*,*)
  
  ! Save original solver settings
  saved_method = sparse_solve_method
  saved_max_iter = bicgstab_max_iter
  saved_fill_level = ilu_fill_level
  saved_l_param = bicgstab_stabilization_param
  saved_abs_tol = bicgstab_abs_tolerance
  saved_rel_tol = bicgstab_rel_tolerance
  saved_verbose = bicgstab_verbose
  
  ! Setup well-conditioned spline (larger λ for good conditioning)
  n_points = n_intervals * 5
  
  ALLOCATE(x(n_points), y(n_points), lambda1(n_points), indx(n_points))
  ALLOCATE(a_ref(n_points), b_ref(n_points), c_ref(n_points), d_ref(n_points))
  ALLOCATE(a_test(n_points), b_test(n_points), c_test(n_points), d_test(n_points))
  
  ! Create well-conditioned spline data
  DO i = 1, n_points
    x(i) = REAL(i-1, DP) / REAL(n_points-1, DP) * 5.0_DP
    y(i) = SIN(x(i)) + 0.1_DP * COS(3.0_DP * x(i))
    lambda1(i) = 1.0e-3_DP  ! Reasonable smoothing - well-conditioned
    indx(i) = i
  END DO
  
  c1 = 0.0_DP
  cn = 0.0_DP
  sw1 = 2
  sw2 = 4
  m = 0.0_DP
  
  WRITE(*,'(A,I0,A)') 'Problem size: ', n_points, ' points'
  WRITE(*,'(A)') 'Spline function: sin(x) + 0.1*cos(3x) with λ=1e-3 (well-conditioned)'
  WRITE(*,*)
  
  ! Get reference solution with UMFPACK
  WRITE(*,'(A)') 'Computing reference solution with UMFPACK...'
  sparse_solve_method = SOLVER_UMFPACK
  CALL splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                    a_ref, b_ref, c_ref, d_ref, m, test_function)
  WRITE(*,'(A)') 'Reference solution computed successfully'
  WRITE(*,*)
  
  ! Configure BiCGSTAB settings
  sparse_solve_method = SOLVER_BICGSTAB
  bicgstab_max_iter = 1000
  bicgstab_abs_tolerance = 1.0e-12_DP
  bicgstab_rel_tolerance = 1.0e-10_DP
  bicgstab_verbose = .FALSE.
  ilu_fill_level = 1
  
  WRITE(*,'(A)') 'Testing BiCGSTAB(ℓ) variants:'
  WRITE(*,*)
  
  ! Test BiCGSTAB(ℓ) for ℓ = 1, 2, 3, 4
  DO l_val = 1, 4
    bicgstab_stabilization_param = l_val
    
    ! Reset boundary conditions
    c1 = 0.0_DP
    cn = 0.0_DP
    
    WRITE(*,'(A,I0,A)', ADVANCE='NO') 'BiCGSTAB(', l_val, '): '
    
    ! Solve with current ℓ value
    CALL splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                      a_test, b_test, c_test, d_test, m, test_function)
    
    ! Compute differences from reference solution
    max_diff_a = MAXVAL(ABS(a_test - a_ref))
    max_diff_b = MAXVAL(ABS(b_test - b_ref))
    max_diff_c = MAXVAL(ABS(c_test - c_ref))
    max_diff_d = MAXVAL(ABS(d_test - d_ref))
    
    total_error = max_diff_a + max_diff_b + max_diff_c + max_diff_d
    
    IF (total_error < 1.0e-8_DP) THEN
      WRITE(*,'(A,E10.3)') 'PASS (error = ', total_error, ')'
    ELSE
      WRITE(*,'(A,E10.3)') 'FAIL (error = ', total_error, ')'
    END IF
  END DO
  
  WRITE(*,*)
  WRITE(*,'(A)') 'BiCGSTAB(ℓ) validation completed'
  
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
    result = 1.0_DP
  END FUNCTION test_function
  
END PROGRAM test_bicgstab_l_validation