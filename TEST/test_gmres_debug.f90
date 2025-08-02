PROGRAM test_gmres_debug
  ! Debug dense GMRES implementation
  
  USE nrtype, ONLY: I4B, DP
  USE gmres_mod, ONLY: gmres_workspace, create_gmres_workspace, destroy_gmres_workspace, &
                       gmres_solve_structured
  IMPLICIT NONE
  
  TYPE(gmres_workspace) :: workspace
  INTEGER(I4B), PARAMETER :: n = 10, restart = 5
  REAL(DP), ALLOCATABLE :: A(:,:), b(:), x(:), x_exact(:), x_initial(:)
  REAL(DP) :: residual_norm, error_norm, tol
  INTEGER(I4B) :: iterations, info, i, j
  LOGICAL :: converged
  
  WRITE(*,*) 'Debug GMRES on simple diagonal system'
  
  ! Allocate arrays
  ALLOCATE(A(n,n), b(n), x(n), x_exact(n), x_initial(n))
  
  ! Create diagonal matrix
  A = 0.0_DP
  DO i = 1, n
    A(i,i) = 2.0_DP
  END DO
  
  ! Setup exact solution
  x_exact = 1.0_DP
  
  ! Compute RHS: b = A * x_exact
  b = 0.0_DP
  DO i = 1, n
    DO j = 1, n
      b(i) = b(i) + A(i,j) * x_exact(j)
    END DO
  END DO
  
  WRITE(*,*) 'Matrix diagonal:', (A(i,i), i=1,MIN(5,n))
  WRITE(*,*) 'RHS vector:', (b(i), i=1,MIN(5,n))
  WRITE(*,*) 'Norm of b:', SQRT(DOT_PRODUCT(b,b))
  
  ! Create workspace
  CALL create_gmres_workspace(workspace, n, restart)
  
  ! Set initial guess
  x_initial = 0.0_DP
  
  ! Solve
  tol = 1.0e-10_DP
  CALL gmres_solve_structured(workspace, A, b, x_initial, 50, tol, &
                              x, iterations, residual_norm, converged, info)
  
  WRITE(*,*) 'Converged:', converged
  WRITE(*,*) 'Iterations:', iterations
  WRITE(*,*) 'Residual norm:', residual_norm
  WRITE(*,*) 'Solution:', (x(i), i=1,MIN(5,n))
  WRITE(*,*) 'Error norm:', SQRT(SUM((x - x_exact)**2))
  
  ! Cleanup
  CALL destroy_gmres_workspace(workspace)
  DEALLOCATE(A, b, x, x_exact, x_initial)
  
END PROGRAM test_gmres_debug