MODULE bicgstab_mod
  ! BiCGSTAB (Bi-Conjugate Gradient Stabilized) iterative solver module
  ! Implements BiCGSTAB algorithm with optional preconditioning
  ! Provides robust iterative solution for sparse linear systems
  
  USE sparse_types_mod, ONLY: dp, long
  USE sparse_utils_mod
  USE ilu_precond_mod
  IMPLICIT NONE
  
  PRIVATE
  
  ! Public types
  PUBLIC :: bicgstab_stats
  
  ! Public interfaces
  PUBLIC :: bicgstab_solve
  PUBLIC :: bicgstab_solve_precond
  
  ! Type for storing solver statistics
  TYPE :: bicgstab_stats
    REAL(kind=dp) :: initial_residual    ! Initial residual norm
    REAL(kind=dp) :: final_residual      ! Final residual norm
    INTEGER :: iterations                 ! Number of iterations performed
    LOGICAL :: converged                  ! Convergence flag
    INTEGER :: restart_count              ! Number of restarts
    REAL(kind=dp) :: solve_time           ! Total solve time (seconds)
  END TYPE bicgstab_stats
  
  ! Generic interfaces
  INTERFACE bicgstab_solve
    MODULE PROCEDURE bicgstab_solve_real
    MODULE PROCEDURE bicgstab_solve_complex
  END INTERFACE bicgstab_solve
  
  INTERFACE bicgstab_solve_precond
    MODULE PROCEDURE bicgstab_solve_precond_real
    MODULE PROCEDURE bicgstab_solve_precond_complex
  END INTERFACE bicgstab_solve_precond
  
CONTAINS

  !==============================================================================
  ! BiCGSTAB solver without preconditioning - Real version
  !==============================================================================
  SUBROUTINE bicgstab_solve_real(n, csr_row_ptr, csr_col_idx, csr_val, b, x, &
                                  tol, max_iter, converged, iter, stats)
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(IN) :: csr_row_ptr(n+1)
    INTEGER, INTENT(IN) :: csr_col_idx(:)
    REAL(kind=dp), INTENT(IN) :: csr_val(:)
    REAL(kind=dp), INTENT(IN) :: b(n)
    REAL(kind=dp), INTENT(INOUT) :: x(n)
    REAL(kind=dp), INTENT(IN) :: tol
    INTEGER, INTENT(IN) :: max_iter
    LOGICAL, INTENT(OUT) :: converged
    INTEGER, INTENT(OUT) :: iter
    TYPE(bicgstab_stats), INTENT(OUT) :: stats
    
    ! Local variables
    REAL(kind=dp), ALLOCATABLE :: r(:), r0(:), p(:), v(:), s(:), t(:)
    REAL(kind=dp) :: rho, rho_old, alpha, omega, beta
    REAL(kind=dp) :: rnorm, bnorm, relative_tol
    INTEGER :: i
    REAL(kind=dp) :: r0_dot_r, r0_dot_v, t_dot_s, t_dot_t
    REAL(kind=dp) :: start_time, end_time
    
    ! Basic input validation
    IF (n <= 0) ERROR STOP "BiCGSTAB: n must be positive"
    IF (tol <= 0.0_dp) ERROR STOP "BiCGSTAB: tolerance must be positive"
    IF (max_iter <= 0) ERROR STOP "BiCGSTAB: max_iter must be positive"
    
    ! Initialize statistics
    CALL CPU_TIME(start_time)
    stats%iterations = 0
    stats%restart_count = 0
    converged = .FALSE.
    iter = 0
    
    ! Allocate work arrays
    ALLOCATE(r(n), r0(n), p(n), v(n), s(n), t(n))
    
    ! Compute initial residual r = b - A*x
    CALL csr_matvec(n, csr_row_ptr, csr_col_idx, csr_val, x, r)
    r = b - r
    
    ! Compute norms
    rnorm = SQRT(DOT_PRODUCT(r, r))
    bnorm = SQRT(DOT_PRODUCT(b, b))
    stats%initial_residual = rnorm
    
    ! Check for zero RHS
    IF (bnorm < 1.0e-15_dp) THEN
      x = 0.0_dp
      converged = .TRUE.
      iter = 0
      stats%final_residual = 0.0_dp
      stats%converged = .TRUE.
      CALL CPU_TIME(end_time)
      stats%solve_time = end_time - start_time
      DEALLOCATE(r, r0, p, v, s, t)
      RETURN
    END IF
    
    ! Set relative tolerance
    relative_tol = tol * bnorm
    
    ! Check initial convergence
    IF (rnorm < relative_tol) THEN
      converged = .TRUE.
      iter = 0
      stats%final_residual = rnorm
      stats%converged = .TRUE.
      CALL CPU_TIME(end_time)
      stats%solve_time = end_time - start_time
      DEALLOCATE(r, r0, p, v, s, t)
      RETURN
    END IF
    
    ! Choose r0 (arbitrary, but must have r0.r != 0)
    r0 = r
    
    ! Initialize BiCGSTAB iteration
    rho = DOT_PRODUCT(r0, r)
    p = r
    
    ! Main BiCGSTAB loop
    DO i = 1, max_iter
      ! Compute v = A*p
      CALL csr_matvec(n, csr_row_ptr, csr_col_idx, csr_val, p, v)
      
      ! Compute alpha = rho / (r0.v)
      r0_dot_v = DOT_PRODUCT(r0, v)
      IF (ABS(r0_dot_v) < 1.0e-15_dp) THEN
        ! Breakdown - restart with current x
        stats%restart_count = stats%restart_count + 1
        r0 = r
        rho = DOT_PRODUCT(r0, r)
        p = r
        CYCLE
      END IF
      alpha = rho / r0_dot_v
      
      ! Compute s = r - alpha*v
      s = r - alpha * v
      
      ! Check convergence of s
      rnorm = SQRT(DOT_PRODUCT(s, s))
      IF (rnorm < relative_tol) THEN
        x = x + alpha * p
        converged = .TRUE.
        iter = i
        stats%final_residual = rnorm
        EXIT
      END IF
      
      ! Compute t = A*s
      CALL csr_matvec(n, csr_row_ptr, csr_col_idx, csr_val, s, t)
      
      ! Compute omega = (t.s) / (t.t)
      t_dot_s = DOT_PRODUCT(t, s)
      t_dot_t = DOT_PRODUCT(t, t)
      IF (t_dot_t < 1.0e-15_dp) THEN
        ! Breakdown
        omega = 0.0_dp
      ELSE
        omega = t_dot_s / t_dot_t
      END IF
      
      ! Update x
      x = x + alpha * p + omega * s
      
      ! Update r
      r = s - omega * t
      
      ! Check convergence
      rnorm = SQRT(DOT_PRODUCT(r, r))
      IF (rnorm < relative_tol) THEN
        converged = .TRUE.
        iter = i
        stats%final_residual = rnorm
        EXIT
      END IF
      
      ! Check for breakdown
      IF (ABS(omega) < 1.0e-15_dp) THEN
        ! Restart
        stats%restart_count = stats%restart_count + 1
        r0 = r
        rho = DOT_PRODUCT(r0, r)
        p = r
        CYCLE
      END IF
      
      ! Prepare for next iteration
      rho_old = rho
      rho = DOT_PRODUCT(r0, r)
      
      IF (ABS(rho) < 1.0e-15_dp) THEN
        ! Breakdown - restart
        stats%restart_count = stats%restart_count + 1
        r0 = r
        rho = DOT_PRODUCT(r0, r)
        p = r
        CYCLE
      END IF
      
      beta = (rho / rho_old) * (alpha / omega)
      p = r + beta * (p - omega * v)
      
      iter = i
    END DO
    
    ! Final statistics
    stats%iterations = iter
    stats%converged = converged
    IF (.NOT. converged) THEN
      stats%final_residual = rnorm
    END IF
    
    CALL CPU_TIME(end_time)
    stats%solve_time = end_time - start_time
    
    ! Cleanup
    DEALLOCATE(r, r0, p, v, s, t)
    
  END SUBROUTINE bicgstab_solve_real
  
  !==============================================================================
  ! BiCGSTAB solver without preconditioning - Complex version
  !==============================================================================
  SUBROUTINE bicgstab_solve_complex(n, csr_row_ptr, csr_col_idx, csr_val, b, x, &
                                     tol, max_iter, converged, iter, stats)
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(IN) :: csr_row_ptr(n+1)
    INTEGER, INTENT(IN) :: csr_col_idx(:)
    COMPLEX(kind=dp), INTENT(IN) :: csr_val(:)
    COMPLEX(kind=dp), INTENT(IN) :: b(n)
    COMPLEX(kind=dp), INTENT(INOUT) :: x(n)
    REAL(kind=dp), INTENT(IN) :: tol
    INTEGER, INTENT(IN) :: max_iter
    LOGICAL, INTENT(OUT) :: converged
    INTEGER, INTENT(OUT) :: iter
    TYPE(bicgstab_stats), INTENT(OUT) :: stats
    
    ! Local variables
    COMPLEX(kind=dp), ALLOCATABLE :: r(:), r0(:), p(:), v(:), s(:), t(:)
    COMPLEX(kind=dp) :: rho, rho_old, alpha, omega, beta
    REAL(kind=dp) :: rnorm, bnorm, relative_tol
    INTEGER :: i
    COMPLEX(kind=dp) :: r0_dot_r, r0_dot_v, t_dot_s, t_dot_t
    REAL(kind=dp) :: start_time, end_time
    
    ! Basic input validation
    IF (n <= 0) ERROR STOP "BiCGSTAB: n must be positive"
    IF (tol <= 0.0_dp) ERROR STOP "BiCGSTAB: tolerance must be positive"
    IF (max_iter <= 0) ERROR STOP "BiCGSTAB: max_iter must be positive"
    
    ! Initialize statistics
    CALL CPU_TIME(start_time)
    stats%iterations = 0
    stats%restart_count = 0
    converged = .FALSE.
    iter = 0
    
    ! Allocate work arrays
    ALLOCATE(r(n), r0(n), p(n), v(n), s(n), t(n))
    
    ! Compute initial residual r = b - A*x
    CALL csr_matvec(n, csr_row_ptr, csr_col_idx, csr_val, x, r)
    r = b - r
    
    ! Compute norms
    rnorm = SQRT(SUM(ABS(r)**2))
    bnorm = SQRT(SUM(ABS(b)**2))
    stats%initial_residual = rnorm
    
    ! Check for zero RHS
    IF (bnorm < 1.0e-15_dp) THEN
      x = (0.0_dp, 0.0_dp)
      converged = .TRUE.
      iter = 0
      stats%final_residual = 0.0_dp
      stats%converged = .TRUE.
      CALL CPU_TIME(end_time)
      stats%solve_time = end_time - start_time
      DEALLOCATE(r, r0, p, v, s, t)
      RETURN
    END IF
    
    ! Set relative tolerance
    relative_tol = tol * bnorm
    
    ! Check initial convergence
    IF (rnorm < relative_tol) THEN
      converged = .TRUE.
      iter = 0
      stats%final_residual = rnorm
      stats%converged = .TRUE.
      CALL CPU_TIME(end_time)
      stats%solve_time = end_time - start_time
      DEALLOCATE(r, r0, p, v, s, t)
      RETURN
    END IF
    
    ! Choose r0 (arbitrary, but must have r0.r != 0)
    r0 = r
    
    ! Initialize BiCGSTAB iteration
    rho = SUM(CONJG(r0) * r)
    p = r
    
    ! Main BiCGSTAB loop
    DO i = 1, max_iter
      ! Compute v = A*p
      CALL csr_matvec(n, csr_row_ptr, csr_col_idx, csr_val, p, v)
      
      ! Compute alpha = rho / (r0.v)
      r0_dot_v = SUM(CONJG(r0) * v)
      IF (ABS(r0_dot_v) < 1.0e-15_dp) THEN
        ! Breakdown - restart with current x
        stats%restart_count = stats%restart_count + 1
        r0 = r
        rho = SUM(CONJG(r0) * r)
        p = r
        CYCLE
      END IF
      alpha = rho / r0_dot_v
      
      ! Compute s = r - alpha*v
      s = r - alpha * v
      
      ! Check convergence of s
      rnorm = SQRT(SUM(ABS(s)**2))
      IF (rnorm < relative_tol) THEN
        x = x + alpha * p
        converged = .TRUE.
        iter = i
        stats%final_residual = rnorm
        EXIT
      END IF
      
      ! Compute t = A*s
      CALL csr_matvec(n, csr_row_ptr, csr_col_idx, csr_val, s, t)
      
      ! Compute omega = (t.s) / (t.t)
      t_dot_s = SUM(CONJG(t) * s)
      t_dot_t = SUM(CONJG(t) * t)
      IF (REAL(t_dot_t) < 1.0e-15_dp) THEN
        ! Breakdown
        omega = (0.0_dp, 0.0_dp)
      ELSE
        omega = t_dot_s / t_dot_t
      END IF
      
      ! Update x
      x = x + alpha * p + omega * s
      
      ! Update r
      r = s - omega * t
      
      ! Check convergence
      rnorm = SQRT(SUM(ABS(r)**2))
      IF (rnorm < relative_tol) THEN
        converged = .TRUE.
        iter = i
        stats%final_residual = rnorm
        EXIT
      END IF
      
      ! Check for breakdown
      IF (ABS(omega) < 1.0e-15_dp) THEN
        ! Restart
        stats%restart_count = stats%restart_count + 1
        r0 = r
        rho = SUM(CONJG(r0) * r)
        p = r
        CYCLE
      END IF
      
      ! Prepare for next iteration
      rho_old = rho
      rho = SUM(CONJG(r0) * r)
      
      IF (ABS(rho) < 1.0e-15_dp) THEN
        ! Breakdown - restart
        stats%restart_count = stats%restart_count + 1
        r0 = r
        rho = SUM(CONJG(r0) * r)
        p = r
        CYCLE
      END IF
      
      beta = (rho / rho_old) * (alpha / omega)
      p = r + beta * (p - omega * v)
      
      iter = i
    END DO
    
    ! Final statistics
    stats%iterations = iter
    stats%converged = converged
    IF (.NOT. converged) THEN
      stats%final_residual = rnorm
    END IF
    
    CALL CPU_TIME(end_time)
    stats%solve_time = end_time - start_time
    
    ! Cleanup
    DEALLOCATE(r, r0, p, v, s, t)
    
  END SUBROUTINE bicgstab_solve_complex
  
  !==============================================================================
  ! BiCGSTAB solver with preconditioning - Real version
  !==============================================================================
  SUBROUTINE bicgstab_solve_precond_real(n, csr_row_ptr, csr_col_idx, csr_val, b, x, &
                                          ilu_fac, tol, max_iter, converged, iter, stats)
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(IN) :: csr_row_ptr(n+1)
    INTEGER, INTENT(IN) :: csr_col_idx(:)
    REAL(kind=dp), INTENT(IN) :: csr_val(:)
    REAL(kind=dp), INTENT(IN) :: b(n)
    REAL(kind=dp), INTENT(INOUT) :: x(n)
    TYPE(ilu_factorization), INTENT(IN) :: ilu_fac
    REAL(kind=dp), INTENT(IN) :: tol
    INTEGER, INTENT(IN) :: max_iter
    LOGICAL, INTENT(OUT) :: converged
    INTEGER, INTENT(OUT) :: iter
    TYPE(bicgstab_stats), INTENT(OUT) :: stats
    
    ! Local variables
    REAL(kind=dp), ALLOCATABLE :: r(:), r0(:), p(:), v(:), s(:), t(:)
    REAL(kind=dp), ALLOCATABLE :: z(:), y(:)  ! Preconditioned vectors
    REAL(kind=dp) :: rho, rho_old, alpha, omega, beta
    REAL(kind=dp) :: rnorm, bnorm, relative_tol
    INTEGER :: i
    REAL(kind=dp) :: r0_dot_z, r0_dot_v, t_dot_s, t_dot_t
    REAL(kind=dp) :: start_time, end_time
    
    ! Basic input validation
    IF (n <= 0) ERROR STOP "BiCGSTAB: n must be positive"
    IF (tol <= 0.0_dp) ERROR STOP "BiCGSTAB: tolerance must be positive"
    IF (max_iter <= 0) ERROR STOP "BiCGSTAB: max_iter must be positive"
    
    ! Initialize statistics
    CALL CPU_TIME(start_time)
    stats%iterations = 0
    stats%restart_count = 0
    converged = .FALSE.
    iter = 0
    
    ! Allocate work arrays
    ALLOCATE(r(n), r0(n), p(n), v(n), s(n), t(n), z(n), y(n))
    
    ! Compute initial residual r = b - A*x
    CALL csr_matvec(n, csr_row_ptr, csr_col_idx, csr_val, x, r)
    r = b - r
    
    ! Compute norms
    rnorm = SQRT(DOT_PRODUCT(r, r))
    bnorm = SQRT(DOT_PRODUCT(b, b))
    stats%initial_residual = rnorm
    
    ! Check for zero RHS
    IF (bnorm < 1.0e-15_dp) THEN
      x = 0.0_dp
      converged = .TRUE.
      iter = 0
      stats%final_residual = 0.0_dp
      stats%converged = .TRUE.
      CALL CPU_TIME(end_time)
      stats%solve_time = end_time - start_time
      DEALLOCATE(r, r0, p, v, s, t, z, y)
      RETURN
    END IF
    
    ! Set relative tolerance
    relative_tol = tol * bnorm
    
    ! Check initial convergence
    IF (rnorm < relative_tol) THEN
      converged = .TRUE.
      iter = 0
      stats%final_residual = rnorm
      stats%converged = .TRUE.
      CALL CPU_TIME(end_time)
      stats%solve_time = end_time - start_time
      DEALLOCATE(r, r0, p, v, s, t, z, y)
      RETURN
    END IF
    
    ! Choose r0 (arbitrary, but must have r0.M^{-1}r != 0)
    r0 = r
    
    ! Apply preconditioner: z = M^{-1}*r
    CALL ilu_solve(ilu_fac, r, z)
    
    ! Initialize BiCGSTAB iteration
    rho = DOT_PRODUCT(r0, z)
    p = z
    
    ! Main BiCGSTAB loop
    DO i = 1, max_iter
      ! Compute v = A*p
      CALL csr_matvec(n, csr_row_ptr, csr_col_idx, csr_val, p, v)
      
      ! Compute alpha = rho / (r0.v)
      r0_dot_v = DOT_PRODUCT(r0, v)
      IF (ABS(r0_dot_v) < 1.0e-15_dp) THEN
        ! Breakdown - restart with current x
        stats%restart_count = stats%restart_count + 1
        r0 = r
        CALL ilu_solve(ilu_fac, r, z)
        rho = DOT_PRODUCT(r0, z)
        p = z
        CYCLE
      END IF
      alpha = rho / r0_dot_v
      
      ! Compute s = r - alpha*v
      s = r - alpha * v
      
      ! Check convergence of s
      rnorm = SQRT(DOT_PRODUCT(s, s))
      IF (rnorm < relative_tol) THEN
        x = x + alpha * p
        converged = .TRUE.
        iter = i
        stats%final_residual = rnorm
        EXIT
      END IF
      
      ! Apply preconditioner: y = M^{-1}*s
      CALL ilu_solve(ilu_fac, s, y)
      
      ! Compute t = A*y
      CALL csr_matvec(n, csr_row_ptr, csr_col_idx, csr_val, y, t)
      
      ! Compute omega = (t.s) / (t.t)
      t_dot_s = DOT_PRODUCT(t, s)
      t_dot_t = DOT_PRODUCT(t, t)
      IF (t_dot_t < 1.0e-15_dp) THEN
        ! Breakdown
        omega = 0.0_dp
      ELSE
        omega = t_dot_s / t_dot_t
      END IF
      
      ! Update x
      x = x + alpha * p + omega * y
      
      ! Update r
      r = s - omega * t
      
      ! Check convergence
      rnorm = SQRT(DOT_PRODUCT(r, r))
      IF (rnorm < relative_tol) THEN
        converged = .TRUE.
        iter = i
        stats%final_residual = rnorm
        EXIT
      END IF
      
      ! Check for breakdown
      IF (ABS(omega) < 1.0e-15_dp) THEN
        ! Restart
        stats%restart_count = stats%restart_count + 1
        r0 = r
        CALL ilu_solve(ilu_fac, r, z)
        rho = DOT_PRODUCT(r0, z)
        p = z
        CYCLE
      END IF
      
      ! Apply preconditioner: z = M^{-1}*r
      CALL ilu_solve(ilu_fac, r, z)
      
      ! Prepare for next iteration
      rho_old = rho
      rho = DOT_PRODUCT(r0, z)
      
      IF (ABS(rho) < 1.0e-15_dp) THEN
        ! Breakdown - restart
        stats%restart_count = stats%restart_count + 1
        r0 = r
        CALL ilu_solve(ilu_fac, r, z)
        rho = DOT_PRODUCT(r0, z)
        p = z
        CYCLE
      END IF
      
      beta = (rho / rho_old) * (alpha / omega)
      p = z + beta * (p - omega * v)
      
      iter = i
    END DO
    
    ! Final statistics
    stats%iterations = iter
    stats%converged = converged
    IF (.NOT. converged) THEN
      stats%final_residual = rnorm
    END IF
    
    CALL CPU_TIME(end_time)
    stats%solve_time = end_time - start_time
    
    ! Cleanup
    DEALLOCATE(r, r0, p, v, s, t, z, y)
    
  END SUBROUTINE bicgstab_solve_precond_real
  
  !==============================================================================
  ! BiCGSTAB solver with preconditioning - Complex version
  !==============================================================================
  SUBROUTINE bicgstab_solve_precond_complex(n, csr_row_ptr, csr_col_idx, csr_val, b, x, &
                                             ilu_fac, tol, max_iter, converged, iter, stats)
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(IN) :: csr_row_ptr(n+1)
    INTEGER, INTENT(IN) :: csr_col_idx(:)
    COMPLEX(kind=dp), INTENT(IN) :: csr_val(:)
    COMPLEX(kind=dp), INTENT(IN) :: b(n)
    COMPLEX(kind=dp), INTENT(INOUT) :: x(n)
    TYPE(ilu_factorization_complex), INTENT(IN) :: ilu_fac
    REAL(kind=dp), INTENT(IN) :: tol
    INTEGER, INTENT(IN) :: max_iter
    LOGICAL, INTENT(OUT) :: converged
    INTEGER, INTENT(OUT) :: iter
    TYPE(bicgstab_stats), INTENT(OUT) :: stats
    
    ! Local variables
    COMPLEX(kind=dp), ALLOCATABLE :: r(:), r0(:), p(:), v(:), s(:), t(:)
    COMPLEX(kind=dp), ALLOCATABLE :: z(:), y(:)  ! Preconditioned vectors
    COMPLEX(kind=dp) :: rho, rho_old, alpha, omega, beta
    REAL(kind=dp) :: rnorm, bnorm, relative_tol
    INTEGER :: i
    COMPLEX(kind=dp) :: r0_dot_z, r0_dot_v, t_dot_s, t_dot_t
    REAL(kind=dp) :: start_time, end_time
    
    ! Basic input validation
    IF (n <= 0) ERROR STOP "BiCGSTAB: n must be positive"
    IF (tol <= 0.0_dp) ERROR STOP "BiCGSTAB: tolerance must be positive"
    IF (max_iter <= 0) ERROR STOP "BiCGSTAB: max_iter must be positive"
    
    ! Initialize statistics
    CALL CPU_TIME(start_time)
    stats%iterations = 0
    stats%restart_count = 0
    converged = .FALSE.
    iter = 0
    
    ! Allocate work arrays
    ALLOCATE(r(n), r0(n), p(n), v(n), s(n), t(n), z(n), y(n))
    
    ! Compute initial residual r = b - A*x
    CALL csr_matvec(n, csr_row_ptr, csr_col_idx, csr_val, x, r)
    r = b - r
    
    ! Compute norms
    rnorm = SQRT(SUM(ABS(r)**2))
    bnorm = SQRT(SUM(ABS(b)**2))
    stats%initial_residual = rnorm
    
    ! Check for zero RHS
    IF (bnorm < 1.0e-15_dp) THEN
      x = (0.0_dp, 0.0_dp)
      converged = .TRUE.
      iter = 0
      stats%final_residual = 0.0_dp
      stats%converged = .TRUE.
      CALL CPU_TIME(end_time)
      stats%solve_time = end_time - start_time
      DEALLOCATE(r, r0, p, v, s, t, z, y)
      RETURN
    END IF
    
    ! Set relative tolerance
    relative_tol = tol * bnorm
    
    ! Check initial convergence
    IF (rnorm < relative_tol) THEN
      converged = .TRUE.
      iter = 0
      stats%final_residual = rnorm
      stats%converged = .TRUE.
      CALL CPU_TIME(end_time)
      stats%solve_time = end_time - start_time
      DEALLOCATE(r, r0, p, v, s, t, z, y)
      RETURN
    END IF
    
    ! Choose r0 (arbitrary, but must have r0.M^{-1}r != 0)
    r0 = r
    
    ! Apply preconditioner: z = M^{-1}*r
    CALL ilu_solve(ilu_fac, r, z)
    
    ! Initialize BiCGSTAB iteration
    rho = SUM(CONJG(r0) * z)
    p = z
    
    ! Main BiCGSTAB loop
    DO i = 1, max_iter
      ! Compute v = A*p
      CALL csr_matvec(n, csr_row_ptr, csr_col_idx, csr_val, p, v)
      
      ! Compute alpha = rho / (r0.v)
      r0_dot_v = SUM(CONJG(r0) * v)
      IF (ABS(r0_dot_v) < 1.0e-15_dp) THEN
        ! Breakdown - restart with current x
        stats%restart_count = stats%restart_count + 1
        r0 = r
        CALL ilu_solve(ilu_fac, r, z)
        rho = SUM(CONJG(r0) * z)
        p = z
        CYCLE
      END IF
      alpha = rho / r0_dot_v
      
      ! Compute s = r - alpha*v
      s = r - alpha * v
      
      ! Check convergence of s
      rnorm = SQRT(SUM(ABS(s)**2))
      IF (rnorm < relative_tol) THEN
        x = x + alpha * p
        converged = .TRUE.
        iter = i
        stats%final_residual = rnorm
        EXIT
      END IF
      
      ! Apply preconditioner: y = M^{-1}*s
      CALL ilu_solve(ilu_fac, s, y)
      
      ! Compute t = A*y
      CALL csr_matvec(n, csr_row_ptr, csr_col_idx, csr_val, y, t)
      
      ! Compute omega = (t.s) / (t.t)
      t_dot_s = SUM(CONJG(t) * s)
      t_dot_t = SUM(CONJG(t) * t)
      IF (REAL(t_dot_t) < 1.0e-15_dp) THEN
        ! Breakdown
        omega = (0.0_dp, 0.0_dp)
      ELSE
        omega = t_dot_s / t_dot_t
      END IF
      
      ! Update x
      x = x + alpha * p + omega * y
      
      ! Update r
      r = s - omega * t
      
      ! Check convergence
      rnorm = SQRT(SUM(ABS(r)**2))
      IF (rnorm < relative_tol) THEN
        converged = .TRUE.
        iter = i
        stats%final_residual = rnorm
        EXIT
      END IF
      
      ! Check for breakdown
      IF (ABS(omega) < 1.0e-15_dp) THEN
        ! Restart
        stats%restart_count = stats%restart_count + 1
        r0 = r
        CALL ilu_solve(ilu_fac, r, z)
        rho = SUM(CONJG(r0) * z)
        p = z
        CYCLE
      END IF
      
      ! Apply preconditioner: z = M^{-1}*r
      CALL ilu_solve(ilu_fac, r, z)
      
      ! Prepare for next iteration
      rho_old = rho
      rho = SUM(CONJG(r0) * z)
      
      IF (ABS(rho) < 1.0e-15_dp) THEN
        ! Breakdown - restart
        stats%restart_count = stats%restart_count + 1
        r0 = r
        CALL ilu_solve(ilu_fac, r, z)
        rho = SUM(CONJG(r0) * z)
        p = z
        CYCLE
      END IF
      
      beta = (rho / rho_old) * (alpha / omega)
      p = z + beta * (p - omega * v)
      
      iter = i
    END DO
    
    ! Final statistics
    stats%iterations = iter
    stats%converged = converged
    IF (.NOT. converged) THEN
      stats%final_residual = rnorm
    END IF
    
    CALL CPU_TIME(end_time)
    stats%solve_time = end_time - start_time
    
    ! Cleanup
    DEALLOCATE(r, r0, p, v, s, t, z, y)
    
  END SUBROUTINE bicgstab_solve_precond_complex

END MODULE bicgstab_mod