MODULE gmres_mod
  !> Minimal GMRES stubs - implementation not yet available
  !! This module provides minimal type definitions and stub procedures
  !! to satisfy dependencies until full GMRES implementation is complete.

  USE nrtype, ONLY: I4B, DP
  IMPLICIT NONE
  PRIVATE

  ! Minimal type definitions to satisfy sparse_solvers_mod dependencies
  TYPE, PUBLIC :: gmres_workspace
    INTEGER :: placeholder = 0
  END TYPE gmres_workspace

  PUBLIC :: gmres_solve_real, gmres_solve_complex, gmres_solve_structured_preconditioned, &
            gmres_solve_structured, create_gmres_workspace, destroy_gmres_workspace

CONTAINS

  SUBROUTINE gmres_solve_real(A, b, x, stats)
    REAL(DP), INTENT(IN) :: A(:,:)
    REAL(DP), INTENT(IN) :: b(:)
    REAL(DP), INTENT(INOUT) :: x(:)
    INTEGER, INTENT(OUT), OPTIONAL :: stats
    
    ! Stub implementation - just return error
    STOP 'ERROR: GMRES solver not implemented yet'
  END SUBROUTINE gmres_solve_real

  SUBROUTINE gmres_solve_complex(A, b, x, stats)
    COMPLEX(DP), INTENT(IN) :: A(:,:)
    COMPLEX(DP), INTENT(IN) :: b(:)
    COMPLEX(DP), INTENT(INOUT) :: x(:)
    INTEGER, INTENT(OUT), OPTIONAL :: stats
    
    ! Stub implementation - just return error
    STOP 'ERROR: GMRES solver not implemented yet'
  END SUBROUTINE gmres_solve_complex

  SUBROUTINE gmres_solve_structured_preconditioned(workspace, n, row_ptr, col_idx, values, &
                                                   b, x, max_iter, tol, precond, &
                                                   result, iter, residual, converged, info, &
                                                   use_preconditioner)
    USE ilu_precond_mod, ONLY: ilu_factorization
    TYPE(gmres_workspace), INTENT(INOUT) :: workspace
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(IN) :: row_ptr(:), col_idx(:)
    REAL(DP), INTENT(IN) :: values(:), b(:)
    REAL(DP), INTENT(INOUT) :: x(:)
    INTEGER, INTENT(IN) :: max_iter
    REAL(DP), INTENT(IN) :: tol
    TYPE(ilu_factorization), INTENT(IN) :: precond  ! ILU preconditioner
    REAL(DP), INTENT(OUT) :: result(:)
    INTEGER, INTENT(OUT) :: iter
    REAL(DP), INTENT(OUT) :: residual
    LOGICAL, INTENT(OUT) :: converged
    INTEGER, INTENT(OUT) :: info
    LOGICAL, INTENT(IN), OPTIONAL :: use_preconditioner
    
    ! Stub implementation - just return error
    STOP 'ERROR: GMRES solver not implemented yet'
  END SUBROUTINE gmres_solve_structured_preconditioned

  SUBROUTINE create_gmres_workspace(workspace, n, restart_dim)
    TYPE(gmres_workspace), INTENT(OUT) :: workspace
    INTEGER, INTENT(IN) :: n, restart_dim
    
    ! Stub implementation - do nothing
    workspace%placeholder = 0
  END SUBROUTINE create_gmres_workspace

  SUBROUTINE destroy_gmres_workspace(workspace)
    TYPE(gmres_workspace), INTENT(INOUT) :: workspace
    
    ! Stub implementation - do nothing
    workspace%placeholder = 0
  END SUBROUTINE destroy_gmres_workspace

  SUBROUTINE gmres_solve_structured(workspace, A, b, x_initial, max_iter, tol, &
                                   x, iter, residual_norm, converged, info)
    TYPE(gmres_workspace), INTENT(INOUT) :: workspace
    REAL(DP), INTENT(IN) :: A(:,:), b(:), x_initial(:)
    INTEGER, INTENT(IN) :: max_iter
    REAL(DP), INTENT(IN) :: tol
    REAL(DP), INTENT(INOUT) :: x(:)
    INTEGER, INTENT(OUT) :: iter
    REAL(DP), INTENT(OUT) :: residual_norm
    LOGICAL, INTENT(OUT) :: converged
    INTEGER, INTENT(OUT) :: info
    
    ! Stub implementation - just return error
    STOP 'ERROR: GMRES solver not implemented yet'
  END SUBROUTINE gmres_solve_structured

END MODULE gmres_mod