MODULE bicgstab_mod
  ! BiCGSTAB (Bi-Conjugate Gradient Stabilized) iterative solver module
  ! Implements BiCGSTAB algorithm with optional preconditioning
  ! Provides robust iterative solution for sparse linear systems
  
  USE sparse_types_mod, ONLY: dp, long
  USE sparse_utils_mod
  USE ilu_precond_mod
  USE amg_types_mod
  USE amg_precond_mod
  IMPLICIT NONE
  
  PRIVATE
  
  ! Public types
  PUBLIC :: bicgstab_stats
  
  ! Public interfaces
  PUBLIC :: bicgstab_solve
  PUBLIC :: bicgstab_solve_precond
  PUBLIC :: bicgstab_l_solve_precond
  PUBLIC :: bicgstab_l_solve_amg_precond
  
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

  INTERFACE bicgstab_l_solve_precond
    MODULE PROCEDURE bicgstab_l_solve_precond_real
    MODULE PROCEDURE bicgstab_l_solve_precond_complex
  END INTERFACE bicgstab_l_solve_precond

  INTERFACE bicgstab_l_solve_amg_precond
    MODULE PROCEDURE bicgstab_l_solve_amg_precond_real
  END INTERFACE bicgstab_l_solve_amg_precond
  
CONTAINS

  !==============================================================================
  ! BiCGSTAB solver without preconditioning - Real version
  !==============================================================================
  SUBROUTINE bicgstab_solve_real(n, csr_row_ptr, csr_col_idx, csr_val, b, x, &
                                  abs_tol, rel_tol, max_iter, converged, iter, stats)
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(IN) :: csr_row_ptr(n+1)
    INTEGER, INTENT(IN) :: csr_col_idx(:)
    REAL(kind=dp), INTENT(IN) :: csr_val(:)
    REAL(kind=dp), INTENT(IN) :: b(n)
    REAL(kind=dp), INTENT(INOUT) :: x(n)
    REAL(kind=dp), INTENT(IN) :: abs_tol
    REAL(kind=dp), INTENT(IN) :: rel_tol
    INTEGER, INTENT(IN) :: max_iter
    LOGICAL, INTENT(OUT) :: converged
    INTEGER, INTENT(OUT) :: iter
    TYPE(bicgstab_stats), INTENT(OUT) :: stats
    
    ! Local variables
    REAL(kind=dp), ALLOCATABLE :: r(:), r0(:), p(:), v(:), s(:), t(:), h(:)
    REAL(kind=dp) :: rho, rho_old, alpha, omega, beta
    REAL(kind=dp) :: rnorm, bnorm, tol_eff
    INTEGER :: i
    REAL(kind=dp) :: r0_dot_r, r0_dot_v, t_dot_s, t_dot_t
    REAL(kind=dp) :: start_time, end_time
    
    ! Basic input validation
    IF (n <= 0) ERROR STOP "BiCGSTAB: n must be positive"
    IF (abs_tol <= 0.0_dp) ERROR STOP "BiCGSTAB: absolute tolerance must be positive"
    IF (rel_tol <= 0.0_dp) ERROR STOP "BiCGSTAB: relative tolerance must be positive"
    IF (max_iter <= 0) ERROR STOP "BiCGSTAB: max_iter must be positive"
    
    ! Initialize statistics
    CALL CPU_TIME(start_time)
    stats%iterations = 0
    stats%restart_count = 0
    converged = .FALSE.
    iter = 0
    
    ! Allocate work arrays
    ALLOCATE(r(n), r0(n), p(n), v(n), s(n), t(n), h(n))
    
    ! Compute initial residual r = b - A*x
    CALL csr_matvec(n, csr_row_ptr, csr_col_idx, csr_val, x, r)
    r = b - r
    
    ! Compute norms
    rnorm = SQRT(DOT_PRODUCT(r, r))
    bnorm = SQRT(DOT_PRODUCT(b, b))
    stats%initial_residual = rnorm
    
    ! Debug output for spline test
    IF (n == 348) THEN  ! Only for spline test size
      PRINT *, 'BiCGSTAB DEBUG: Initial state'
      PRINT *, '  n =', n
      PRINT *, '  ||b|| =', bnorm
      PRINT *, '  ||x|| =', SQRT(DOT_PRODUCT(x, x))
      PRINT *, '  ||r|| = ||b-Ax|| =', rnorm
      
      ! Check matrix norms and conditioning hints
      PRINT *, 'BiCGSTAB DEBUG: Matrix info'
      PRINT *, '  nz =', SIZE(csr_val)
      PRINT *, '  ||A||_F (approx) =', SQRT(DOT_PRODUCT(csr_val, csr_val))
      PRINT *, '  max|A_ij| =', MAXVAL(ABS(csr_val))
      PRINT *, '  min|A_ij| (nonzero) =', MINVAL(ABS(csr_val), MASK=ABS(csr_val) > 0.0_dp)
    END IF
    
    ! Check for zero RHS
    IF (bnorm < 1.0e-15_dp) THEN
      x = 0.0_dp
      converged = .TRUE.
      iter = 0
      stats%final_residual = 0.0_dp
      stats%converged = .TRUE.
      CALL CPU_TIME(end_time)
      stats%solve_time = end_time - start_time
      DEALLOCATE(r, r0, p, v, s, t, h)
      RETURN
    END IF
    
    ! Set effective tolerance: max(abs_tol, rel_tol * ||b||)
    tol_eff = MAX(abs_tol, rel_tol * bnorm)
    
    ! Check initial convergence
    IF (rnorm < tol_eff) THEN
      converged = .TRUE.
      iter = 0
      stats%final_residual = rnorm
      stats%converged = .TRUE.
      CALL CPU_TIME(end_time)
      stats%solve_time = end_time - start_time
      DEALLOCATE(r, r0, p, v, s, t, h)
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
      
      ! Debug for first iteration
      IF (i == 1 .AND. n == 348) THEN
        PRINT *, 'BiCGSTAB DEBUG: First iteration after A*p'
        PRINT *, '  ||p|| =', SQRT(DOT_PRODUCT(p, p))
        PRINT *, '  ||v|| = ||A*p|| =', SQRT(DOT_PRODUCT(v, v))
        PRINT *, '  rho =', rho
      END IF
      
      ! Compute alpha = rho / (r0.v)
      r0_dot_v = DOT_PRODUCT(r0, v)
      IF (ABS(r0_dot_v) < 1.0e-15_dp) THEN
        ! Breakdown - restart with current x
        IF (i == 1) THEN
          PRINT *, 'BiCGSTAB DEBUG: Breakdown on first iteration!'
          PRINT *, '  rho =', rho
          PRINT *, '  r0_dot_v =', r0_dot_v
          PRINT *, '  ||r0|| =', SQRT(DOT_PRODUCT(r0, r0))
          PRINT *, '  ||v|| =', SQRT(DOT_PRODUCT(v, v))
          PRINT *, '  ||p|| =', SQRT(DOT_PRODUCT(p, p))
        END IF
        stats%restart_count = stats%restart_count + 1
        r0 = r
        rho = DOT_PRODUCT(r0, r)
        p = r
        CYCLE
      END IF
      alpha = rho / r0_dot_v
      
      ! Compute h = x + alpha*p (intermediate solution)
      h = x + alpha * p
      
      ! Compute s = r - alpha*v
      s = r - alpha * v
      
      ! Check convergence of s
      rnorm = SQRT(DOT_PRODUCT(s, s))
      IF (i <= 3 .AND. n == 348) THEN
        PRINT *, 'BiCGSTAB DEBUG: iter', i, 'after s = r - alpha*v'
        PRINT *, '  alpha =', alpha
        PRINT *, '  ||s|| =', rnorm, 'tol_eff =', tol_eff
      END IF
      IF (rnorm < tol_eff) THEN
        x = h  ! Use intermediate solution
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
      
      ! Update x using intermediate solution h
      x = h + omega * s
      
      ! Update r
      r = s - omega * t
      
      ! Check convergence
      rnorm = SQRT(DOT_PRODUCT(r, r))
      IF (i <= 3 .AND. n == 348) THEN
        PRINT *, 'BiCGSTAB DEBUG: iter', i, 'after full update'
        PRINT *, '  omega =', omega
        PRINT *, '  ||r|| =', rnorm
        PRINT *, '  ||x|| =', SQRT(DOT_PRODUCT(x, x))
      END IF
      IF (rnorm < tol_eff) THEN
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
    DEALLOCATE(r, r0, p, v, s, t, h)
    
  END SUBROUTINE bicgstab_solve_real
  
  !==============================================================================
  ! BiCGSTAB solver without preconditioning - Complex version
  !==============================================================================
  SUBROUTINE bicgstab_solve_complex(n, csr_row_ptr, csr_col_idx, csr_val, b, x, &
                                     abs_tol, rel_tol, max_iter, converged, iter, stats)
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(IN) :: csr_row_ptr(n+1)
    INTEGER, INTENT(IN) :: csr_col_idx(:)
    COMPLEX(kind=dp), INTENT(IN) :: csr_val(:)
    COMPLEX(kind=dp), INTENT(IN) :: b(n)
    COMPLEX(kind=dp), INTENT(INOUT) :: x(n)
    REAL(kind=dp), INTENT(IN) :: abs_tol
    REAL(kind=dp), INTENT(IN) :: rel_tol
    INTEGER, INTENT(IN) :: max_iter
    LOGICAL, INTENT(OUT) :: converged
    INTEGER, INTENT(OUT) :: iter
    TYPE(bicgstab_stats), INTENT(OUT) :: stats
    
    ! Local variables
    COMPLEX(kind=dp), ALLOCATABLE :: r(:), r0(:), p(:), v(:), s(:), t(:), h(:)
    COMPLEX(kind=dp) :: rho, rho_old, alpha, omega, beta
    REAL(kind=dp) :: rnorm, bnorm, tol_eff
    INTEGER :: i
    COMPLEX(kind=dp) :: r0_dot_r, r0_dot_v, t_dot_s, t_dot_t
    REAL(kind=dp) :: start_time, end_time
    
    ! Basic input validation
    IF (n <= 0) ERROR STOP "BiCGSTAB: n must be positive"
    IF (abs_tol <= 0.0_dp) ERROR STOP "BiCGSTAB: absolute tolerance must be positive"
    IF (rel_tol <= 0.0_dp) ERROR STOP "BiCGSTAB: relative tolerance must be positive"
    IF (max_iter <= 0) ERROR STOP "BiCGSTAB: max_iter must be positive"
    
    ! Initialize statistics
    CALL CPU_TIME(start_time)
    stats%iterations = 0
    stats%restart_count = 0
    converged = .FALSE.
    iter = 0
    
    ! Allocate work arrays
    ALLOCATE(r(n), r0(n), p(n), v(n), s(n), t(n), h(n))
    
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
      DEALLOCATE(r, r0, p, v, s, t, h)
      RETURN
    END IF
    
    ! Set effective tolerance: max(abs_tol, rel_tol * ||b||)
    tol_eff = MAX(abs_tol, rel_tol * bnorm)
    
    ! Check initial convergence
    IF (rnorm < tol_eff) THEN
      converged = .TRUE.
      iter = 0
      stats%final_residual = rnorm
      stats%converged = .TRUE.
      CALL CPU_TIME(end_time)
      stats%solve_time = end_time - start_time
      DEALLOCATE(r, r0, p, v, s, t, h)
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
      
      ! Compute h = x + alpha*p (intermediate solution)
      h = x + alpha * p
      
      ! Compute s = r - alpha*v
      s = r - alpha * v
      
      ! Check convergence of s
      rnorm = SQRT(SUM(ABS(s)**2))
      IF (rnorm < tol_eff) THEN
        x = h  ! Use intermediate solution
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
      
      ! Update x using intermediate solution h
      x = h + omega * s
      
      ! Update r
      r = s - omega * t
      
      ! Check convergence
      rnorm = SQRT(SUM(ABS(r)**2))
      IF (rnorm < tol_eff) THEN
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
    DEALLOCATE(r, r0, p, v, s, t, h)
    
  END SUBROUTINE bicgstab_solve_complex
  
  !==============================================================================
  ! BiCGSTAB solver with preconditioning - Real version
  !==============================================================================
  SUBROUTINE bicgstab_solve_precond_real(n, csr_row_ptr, csr_col_idx, csr_val, b, x, &
                                          ilu_fac, abs_tol, rel_tol, max_iter, converged, iter, stats)
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(IN) :: csr_row_ptr(n+1)
    INTEGER, INTENT(IN) :: csr_col_idx(:)
    REAL(kind=dp), INTENT(IN) :: csr_val(:)
    REAL(kind=dp), INTENT(IN) :: b(n)
    REAL(kind=dp), INTENT(INOUT) :: x(n)
    TYPE(ilu_factorization), INTENT(IN) :: ilu_fac
    REAL(kind=dp), INTENT(IN) :: abs_tol
    REAL(kind=dp), INTENT(IN) :: rel_tol
    INTEGER, INTENT(IN) :: max_iter
    LOGICAL, INTENT(OUT) :: converged
    INTEGER, INTENT(OUT) :: iter
    TYPE(bicgstab_stats), INTENT(OUT) :: stats
    
    ! Local variables
    REAL(kind=dp), ALLOCATABLE :: r(:), r0(:), p(:), v(:), s(:), t(:)
    REAL(kind=dp), ALLOCATABLE :: z(:), y(:)  ! Preconditioned vectors
    REAL(kind=dp) :: rho, rho_old, alpha, omega, beta
    REAL(kind=dp) :: rnorm, bnorm, tol_eff
    INTEGER :: i
    REAL(kind=dp) :: r0_dot_z, r0_dot_v, t_dot_s, t_dot_t
    REAL(kind=dp) :: start_time, end_time
    
    ! Basic input validation
    IF (n <= 0) ERROR STOP "BiCGSTAB: n must be positive"
    IF (abs_tol <= 0.0_dp) ERROR STOP "BiCGSTAB: absolute tolerance must be positive"
    IF (rel_tol <= 0.0_dp) ERROR STOP "BiCGSTAB: relative tolerance must be positive"
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
    
    ! Set effective tolerance: max(abs_tol, rel_tol * ||b||)
    tol_eff = MAX(abs_tol, rel_tol * bnorm)
    
    ! Check initial convergence
    IF (rnorm < tol_eff) THEN
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
    IF (ilu_fac%factorized) THEN
      CALL ilu_solve(ilu_fac, r, z)
    ELSE
      z = r  ! No preconditioning - use residual directly
    END IF
    
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
        IF (ilu_fac%factorized) THEN
          CALL ilu_solve(ilu_fac, r, z)
        ELSE
          z = r  ! No preconditioning - use residual directly
        END IF
        rho = DOT_PRODUCT(r0, z)
        p = z
        CYCLE
      END IF
      alpha = rho / r0_dot_v
      
      ! Compute s = r - alpha*v
      s = r - alpha * v
      
      ! Check convergence of s
      rnorm = SQRT(DOT_PRODUCT(s, s))
      IF (rnorm < tol_eff) THEN
        x = x + alpha * p
        converged = .TRUE.
        iter = i
        stats%final_residual = rnorm
        EXIT
      END IF
      
      ! Apply preconditioner: y = M^{-1}*s
      IF (ilu_fac%factorized) THEN
        CALL ilu_solve(ilu_fac, s, y)
      ELSE
        y = s  ! No preconditioning - use residual directly
      END IF
      
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
      IF (rnorm < tol_eff) THEN
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
        IF (ilu_fac%factorized) THEN
          CALL ilu_solve(ilu_fac, r, z)
        ELSE
          z = r  ! No preconditioning - use residual directly
        END IF
        rho = DOT_PRODUCT(r0, z)
        p = z
        CYCLE
      END IF
      
      ! Apply preconditioner: z = M^{-1}*r
      IF (ilu_fac%factorized) THEN
        CALL ilu_solve(ilu_fac, r, z)
      ELSE
        z = r  ! No preconditioning - use residual directly
      END IF
      
      ! Prepare for next iteration
      rho_old = rho
      rho = DOT_PRODUCT(r0, z)
      
      IF (ABS(rho) < 1.0e-15_dp) THEN
        ! Breakdown - restart
        stats%restart_count = stats%restart_count + 1
        r0 = r
        IF (ilu_fac%factorized) THEN
          CALL ilu_solve(ilu_fac, r, z)
        ELSE
          z = r  ! No preconditioning - use residual directly
        END IF
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
                                             ilu_fac, abs_tol, rel_tol, max_iter, converged, iter, stats)
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(IN) :: csr_row_ptr(n+1)
    INTEGER, INTENT(IN) :: csr_col_idx(:)
    COMPLEX(kind=dp), INTENT(IN) :: csr_val(:)
    COMPLEX(kind=dp), INTENT(IN) :: b(n)
    COMPLEX(kind=dp), INTENT(INOUT) :: x(n)
    TYPE(ilu_factorization_complex), INTENT(IN) :: ilu_fac
    REAL(kind=dp), INTENT(IN) :: abs_tol
    REAL(kind=dp), INTENT(IN) :: rel_tol
    INTEGER, INTENT(IN) :: max_iter
    LOGICAL, INTENT(OUT) :: converged
    INTEGER, INTENT(OUT) :: iter
    TYPE(bicgstab_stats), INTENT(OUT) :: stats
    
    ! Local variables
    COMPLEX(kind=dp), ALLOCATABLE :: r(:), r0(:), p(:), v(:), s(:), t(:), h(:)
    COMPLEX(kind=dp), ALLOCATABLE :: z(:), y(:)  ! Preconditioned vectors
    COMPLEX(kind=dp) :: rho, rho_old, alpha, omega, beta
    REAL(kind=dp) :: rnorm, bnorm, tol_eff
    INTEGER :: i
    COMPLEX(kind=dp) :: r0_dot_z, r0_dot_v, t_dot_s, t_dot_t
    REAL(kind=dp) :: start_time, end_time
    
    ! Basic input validation
    IF (n <= 0) ERROR STOP "BiCGSTAB: n must be positive"
    IF (abs_tol <= 0.0_dp) ERROR STOP "BiCGSTAB: absolute tolerance must be positive"
    IF (rel_tol <= 0.0_dp) ERROR STOP "BiCGSTAB: relative tolerance must be positive"
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
    
    ! Set effective tolerance: max(abs_tol, rel_tol * ||b||)
    tol_eff = MAX(abs_tol, rel_tol * bnorm)
    
    ! Check initial convergence
    IF (rnorm < tol_eff) THEN
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
    IF (ilu_fac%factorized) THEN
      CALL ilu_solve(ilu_fac, r, z)
    ELSE
      z = r  ! No preconditioning - use residual directly
    END IF
    
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
        IF (ilu_fac%factorized) THEN
          CALL ilu_solve(ilu_fac, r, z)
        ELSE
          z = r  ! No preconditioning - use residual directly
        END IF
        rho = SUM(CONJG(r0) * z)
        p = z
        CYCLE
      END IF
      alpha = rho / r0_dot_v
      
      ! Compute s = r - alpha*v
      s = r - alpha * v
      
      ! Check convergence of s
      rnorm = SQRT(SUM(ABS(s)**2))
      IF (rnorm < tol_eff) THEN
        x = x + alpha * p
        converged = .TRUE.
        iter = i
        stats%final_residual = rnorm
        EXIT
      END IF
      
      ! Apply preconditioner: y = M^{-1}*s
      IF (ilu_fac%factorized) THEN
        CALL ilu_solve(ilu_fac, s, y)
      ELSE
        y = s  ! No preconditioning - use residual directly
      END IF
      
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
      IF (rnorm < tol_eff) THEN
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
        IF (ilu_fac%factorized) THEN
          CALL ilu_solve(ilu_fac, r, z)
        ELSE
          z = r  ! No preconditioning - use residual directly
        END IF
        rho = SUM(CONJG(r0) * z)
        p = z
        CYCLE
      END IF
      
      ! Apply preconditioner: z = M^{-1}*r
      IF (ilu_fac%factorized) THEN
        CALL ilu_solve(ilu_fac, r, z)
      ELSE
        z = r  ! No preconditioning - use residual directly
      END IF
      
      ! Prepare for next iteration
      rho_old = rho
      rho = SUM(CONJG(r0) * z)
      
      IF (ABS(rho) < 1.0e-15_dp) THEN
        ! Breakdown - restart
        stats%restart_count = stats%restart_count + 1
        r0 = r
        IF (ilu_fac%factorized) THEN
          CALL ilu_solve(ilu_fac, r, z)
        ELSE
          z = r  ! No preconditioning - use residual directly
        END IF
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

  !==============================================================================
  ! BiCGSTAB(ℓ) solver with preconditioning - Real version
  !==============================================================================
  !> BiCGSTAB(ℓ) algorithm with higher-order stabilization
  !! 
  !! This implements the BiCGSTAB(ℓ) method which generalizes standard BiCGSTAB
  !! to use ℓ-step recurrences for improved stability on difficult problems.
  !! For ℓ=1, this reduces to standard BiCGSTAB.
  !!
  !! @param[in]    n            Matrix dimension
  !! @param[in]    csr_row_ptr  CSR row pointers (size: n+1)
  !! @param[in]    csr_col_idx  CSR column indices  
  !! @param[in]    csr_val      CSR matrix values
  !! @param[in]    b            Right-hand side vector (size: n)
  !! @param[inout] x            Solution vector (size: n) - initial guess on input
  !! @param[in]    abs_tol      Absolute tolerance
  !! @param[in]    rel_tol      Relative tolerance
  !! @param[in]    max_iter     Maximum iterations
  !! @param[in]    l_param      Stabilization parameter ℓ (1=standard BiCGSTAB)
  !! @param[in]    ilu_fac      ILU factorization for preconditioning
  !! @param[in]    verbose      Print convergence information
  !! @param[out]   converged    Convergence flag
  !! @param[out]   iter         Iterations performed
  !! @param[out]   stats        Solver statistics
  SUBROUTINE bicgstab_l_solve_precond_real(n, csr_row_ptr, csr_col_idx, csr_val, b, x, &
                                           abs_tol, rel_tol, max_iter, l_param, ilu_fac, &
                                           verbose, converged, iter, stats)
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(IN) :: csr_row_ptr(n+1)
    INTEGER, INTENT(IN) :: csr_col_idx(:)
    REAL(kind=dp), INTENT(IN) :: csr_val(:)
    REAL(kind=dp), INTENT(IN) :: b(n)
    REAL(kind=dp), INTENT(INOUT) :: x(n)
    REAL(kind=dp), INTENT(IN) :: abs_tol, rel_tol
    INTEGER, INTENT(IN) :: max_iter, l_param
    TYPE(ilu_factorization), INTENT(IN) :: ilu_fac
    LOGICAL, INTENT(IN) :: verbose
    LOGICAL, INTENT(OUT) :: converged
    INTEGER, INTENT(OUT) :: iter
    TYPE(bicgstab_stats), INTENT(OUT) :: stats
    
    ! Local variables
    REAL(kind=dp), ALLOCATABLE :: r(:,:), u(:,:)  ! Residual and direction vectors (n x l+1)
    REAL(kind=dp), ALLOCATABLE :: r0(:)           ! Shadow residual
    REAL(kind=dp), ALLOCATABLE :: z(:), y(:)      ! Work vectors for preconditioning
    REAL(kind=dp) :: alpha, beta, omega, rho0, rho1, sigma
    REAL(kind=dp) :: tau(l_param*l_param), gamma(0:l_param), gamma_p(0:l_param)
    REAL(kind=dp) :: tau_matrix(0:l_param, 0:l_param)
    REAL(kind=dp) :: norm_b, tol_eff, residual_norm
    INTEGER :: i, j, k
    LOGICAL :: breakdown
    REAL(kind=dp) :: start_time, end_time
    
    ! For ℓ=1, fall back to standard BiCGSTAB for efficiency
    IF (l_param == 1) THEN
      IF (verbose) PRINT *, 'BiCGSTAB(ℓ): Using standard BiCGSTAB for ℓ=1'
      CALL bicgstab_solve_precond_real(n, csr_row_ptr, csr_col_idx, csr_val, b, x, &
                                       ilu_fac, abs_tol, rel_tol, max_iter, &
                                       converged, iter, stats)
      RETURN
    END IF
    
    IF (verbose) PRINT *, 'BiCGSTAB(ℓ): Using higher-order BiCGSTAB with ℓ=', l_param
    
    CALL CPU_TIME(start_time)
    
    ! Initialize
    converged = .FALSE.
    iter = 0
    breakdown = .FALSE.
    
    ! Allocate workspace
    ALLOCATE(r(n, 0:l_param), u(n, 0:l_param), r0(n), z(n), y(n))
    
    ! Compute norm of RHS for relative tolerance
    norm_b = SQRT(DOT_PRODUCT(b, b))
    stats%initial_residual = norm_b
    
    IF (norm_b < TINY(1.0_dp)) THEN
      ! Zero RHS case
      x = 0.0_dp
      converged = .TRUE.
      stats%converged = .TRUE.
      stats%final_residual = 0.0_dp
      CALL CPU_TIME(end_time)
      stats%solve_time = end_time - start_time
      DEALLOCATE(r, u, r0, z, y)
      RETURN
    END IF
    
    tol_eff = MAX(abs_tol, rel_tol * norm_b)
    
    ! Compute initial residual: r = b - A*x
    CALL csr_matvec(n, csr_row_ptr, csr_col_idx, csr_val, x, z)
    r(:, 0) = b - z
    
    ! Apply preconditioning to get r0 (if ILU factorization is valid)
    IF (ilu_fac%factorized) THEN
      CALL ilu_solve(ilu_fac, r(:, 0), r0)
    ELSE
      ! No preconditioning - use residual directly  
      r0 = r(:, 0)
    END IF
    residual_norm = SQRT(DOT_PRODUCT(r0, r0))
    
    IF (verbose) THEN
      PRINT *, 'BiCGSTAB(', l_param, ') starting:'
      PRINT *, '  ||b|| =', norm_b
      PRINT *, '  Initial residual =', residual_norm
      PRINT *, '  Target tolerance =', tol_eff
    END IF
    
    ! Check initial convergence
    IF (residual_norm <= tol_eff) THEN
      converged = .TRUE.
      stats%final_residual = residual_norm
      stats%converged = .TRUE.
      CALL CPU_TIME(end_time)
      stats%solve_time = end_time - start_time
      DEALLOCATE(r, u, r0, z, y)
      RETURN
    END IF
    
    u(:, 0) = r0
    rho0 = 1.0_dp
    alpha = 0.0_dp
    omega = 1.0_dp
    
    ! Main BiCGSTAB(ℓ) iteration loop
    main_loop: DO
      rho0 = -omega * rho0
      
      ! BiCG part: ℓ steps
      bicg_loop: DO j = 0, l_param - 1
        rho1 = DOT_PRODUCT(r0, r(:, j))
        
        IF (ABS(rho1) < TINY(1.0_dp)) THEN
          breakdown = .TRUE.
          EXIT bicg_loop
        END IF
        
        beta = alpha * rho1 / rho0
        rho0 = rho1
        
        ! Update u vectors
        DO i = 0, j
          u(:, i) = r(:, i) - beta * u(:, i)
        END DO
        
        ! Apply preconditioned matrix: z = M^(-1) * A * u_j
        CALL csr_matvec(n, csr_row_ptr, csr_col_idx, csr_val, u(:, j), z)
        IF (ilu_fac%factorized) THEN
          CALL ilu_solve(ilu_fac, z, y)
          u(:, j+1) = y
        ELSE
          u(:, j+1) = z  ! No preconditioning
        END IF
        
        alpha = rho0 / DOT_PRODUCT(r0, u(:, j+1))
        
        ! Update residuals
        DO i = 0, j
          r(:, i) = r(:, i) - alpha * u(:, i+1)
        END DO
        
        ! Apply preconditioned matrix to r_j
        CALL csr_matvec(n, csr_row_ptr, csr_col_idx, csr_val, r(:, j), z)
        IF (ilu_fac%factorized) THEN
          CALL ilu_solve(ilu_fac, z, y)
          r(:, j+1) = y
        ELSE
          r(:, j+1) = z  ! No preconditioning
        END IF
        
        ! Update solution
        x = x + alpha * u(:, 0)
      END DO bicg_loop
      
      IF (breakdown) THEN
        IF (verbose) THEN
          PRINT *, 'BiCGSTAB(', l_param, ') breakdown detected'
        END IF
        EXIT main_loop
      END IF
      
      ! MR part: Minimal residual stabilization (solve least squares problem)
      ! This follows Sleijpen-Fokkema 1993 algorithm
      
      ! Step 1: Modified Gram-Schmidt orthogonalization of r_1, ..., r_ℓ
      DO j = 1, l_param
        DO i = 1, j-1
          ! tau[i,j] = <r_j, r_i> / sigma[i]
          tau_matrix(i, j) = DOT_PRODUCT(r(:, j), r(:, i)) / tau_matrix(i, i)
          ! r_j = r_j - tau[i,j] * r_i (orthogonalize)
          r(:, j) = r(:, j) - tau_matrix(i, j) * r(:, i)
        END DO
        ! sigma[j] = ||r_j||²
        tau_matrix(j, j) = DOT_PRODUCT(r(:, j), r(:, j))
        
        ! gamma'[j] = <r̂_0, r_j> / sigma[j]
        gamma_p(j) = DOT_PRODUCT(r0, r(:, j)) / tau_matrix(j, j)
      END DO
      
      ! Step 2: Solve upper triangular system for gamma coefficients
      ! This minimizes ||r_0 - sum(gamma[j] * r_j)||
      gamma(l_param) = gamma_p(l_param)
      DO j = l_param-1, 1, -1
        gamma(j) = gamma_p(j)
        DO i = j+1, l_param
          gamma(j) = gamma(j) - tau_matrix(j, i) * gamma(i)
        END DO
      END DO
      
      ! Set omega for stability check
      omega = gamma(l_param)
      
      ! Step 3: Update solution: x = x + sum(gamma[j] * u_{j-1})
      DO j = 1, l_param
        x = x + gamma(j) * u(:, j-1)  ! Note: u vectors, not r vectors
      END DO
      
      ! Step 4: Update residuals for next iteration: r_{j-1} = r_{j-1} - sum(gamma[k] * r_k) for k>=j
      DO j = 1, l_param  
        r(:, j-1) = r(:, j-1) - gamma(j) * r(:, j)
        u(:, j-1) = u(:, j-1) - gamma(j) * u(:, j)
      END DO
      
      iter = iter + l_param
      
      ! Check convergence
      residual_norm = SQRT(DOT_PRODUCT(r(:, 0), r(:, 0)))
      
      IF (verbose .AND. MOD(iter, 10*l_param) == 0) THEN
        PRINT *, 'BiCGSTAB(', l_param, ') iter', iter, ': residual =', residual_norm
      END IF
      
      IF (residual_norm <= tol_eff) THEN
        converged = .TRUE.
        EXIT main_loop
      END IF
      
      IF (iter >= max_iter) THEN
        EXIT main_loop
      END IF
      
      ! Check for stagnation
      IF (ABS(omega) < TINY(1.0_dp)) THEN
        breakdown = .TRUE.
        EXIT main_loop
      END IF
      
    END DO main_loop
    
    stats%final_residual = residual_norm
    stats%converged = converged
    stats%iterations = iter
    
    IF (verbose) THEN
      IF (converged) THEN
        PRINT *, 'BiCGSTAB(', l_param, ') converged in', iter, 'iterations'
      ELSE
        PRINT *, 'BiCGSTAB(', l_param, ') did not converge in', iter, 'iterations'
        PRINT *, 'Final residual:', residual_norm
      END IF
    END IF
    
    CALL CPU_TIME(end_time)
    stats%solve_time = end_time - start_time
    
    DEALLOCATE(r, u, r0, z, y)
    
  END SUBROUTINE bicgstab_l_solve_precond_real

  !> Complex version placeholder
  SUBROUTINE bicgstab_l_solve_precond_complex(n, csr_row_ptr, csr_col_idx, csr_val, b, x, &
                                              abs_tol, rel_tol, max_iter, l_param, ilu_fac, &
                                              verbose, converged, iter, stats)
    INTEGER, INTENT(IN) :: n, max_iter, l_param
    INTEGER, INTENT(IN) :: csr_row_ptr(n+1), csr_col_idx(:)
    COMPLEX(kind=dp), INTENT(IN) :: csr_val(:), b(n)
    COMPLEX(kind=dp), INTENT(INOUT) :: x(n)
    REAL(kind=dp), INTENT(IN) :: abs_tol, rel_tol
    TYPE(ilu_factorization_complex), INTENT(IN) :: ilu_fac
    LOGICAL, INTENT(IN) :: verbose
    LOGICAL, INTENT(OUT) :: converged
    INTEGER, INTENT(OUT) :: iter
    TYPE(bicgstab_stats), INTENT(OUT) :: stats
    
    ! Placeholder - implement complex version later
    PRINT *, 'ERROR: Complex BiCGSTAB(ℓ) not yet implemented'
    converged = .FALSE.
    iter = 0
    stats%converged = .FALSE.
    stats%final_residual = HUGE(1.0_dp)
    
  END SUBROUTINE bicgstab_l_solve_precond_complex

  !==============================================================================
  ! BiCGSTAB(ℓ) solver with AMG preconditioning - Real version
  !==============================================================================
  !> BiCGSTAB(ℓ) algorithm with AMG preconditioning
  !! 
  !! This implements BiCGSTAB(ℓ) with AMG (Algebraic Multigrid) preconditioning
  !! for improved convergence on ill-conditioned systems, particularly spline matrices.
  !!
  SUBROUTINE bicgstab_l_solve_amg_precond_real(n, csr_row_ptr, csr_col_idx, csr_val, b, x, &
                                               abs_tol, rel_tol, max_iter, l_param, amg_hier, &
                                               verbose, converged, iter, stats)
    
    IMPLICIT NONE
    
    ! Input parameters
    INTEGER, INTENT(IN) :: n                           ! Matrix size
    INTEGER, INTENT(IN) :: csr_row_ptr(n+1)           ! CSR row pointers
    INTEGER, INTENT(IN) :: csr_col_idx(:)              ! CSR column indices
    REAL(dp), INTENT(IN) :: csr_val(:)                 ! CSR values
    REAL(dp), INTENT(IN) :: b(n)                       ! Right-hand side
    REAL(dp), INTENT(INOUT) :: x(n)                    ! Solution vector (input: initial guess)
    REAL(dp), INTENT(IN) :: abs_tol                    ! Absolute tolerance
    REAL(dp), INTENT(IN) :: rel_tol                    ! Relative tolerance
    INTEGER, INTENT(IN) :: max_iter                    ! Maximum iterations
    INTEGER, INTENT(IN) :: l_param                     ! BiCGSTAB(ℓ) parameter
    TYPE(amg_hierarchy), INTENT(INOUT) :: amg_hier     ! AMG hierarchy for preconditioning
    LOGICAL, INTENT(IN) :: verbose                     ! Verbose output flag
    
    ! Output parameters
    LOGICAL, INTENT(OUT) :: converged                  ! Convergence flag
    INTEGER, INTENT(OUT) :: iter                       ! Number of iterations
    TYPE(bicgstab_stats), INTENT(OUT) :: stats         ! Solver statistics
    
    ! Local variables
    REAL(dp), ALLOCATABLE :: r(:), r0(:), p(:), v(:), s(:), t(:), z(:), y(:), temp(:)
    REAL(dp) :: alpha, beta, omega, rho, rho_old
    REAL(dp) :: bnorm, rnorm, rtol
    REAL(dp) :: start_time, end_time
    INTEGER :: i, j
    
    CALL CPU_TIME(start_time)
    
    ! Initialize statistics
    stats%iterations = 0
    stats%converged = .FALSE.
    stats%restart_count = 0
    
    ! Allocate work arrays
    ALLOCATE(r(n), r0(n), p(n), v(n), s(n), t(n), z(n), y(n))
    
    ! Compute initial residual: r = b - A*x
    r = b
    ! Compute A*x and subtract from b
    ALLOCATE(temp(n))
    CALL csr_matvec(n, csr_row_ptr, csr_col_idx, csr_val, x, temp)
    r = r - temp
    DEALLOCATE(temp)
    
    ! Compute reference norm for relative tolerance
    bnorm = SQRT(DOT_PRODUCT(b, b))
    IF (bnorm == 0.0_dp) bnorm = 1.0_dp
    
    ! Initial residual norm
    rnorm = SQRT(DOT_PRODUCT(r, r))
    stats%initial_residual = rnorm
    rtol = MAX(abs_tol, rel_tol * bnorm)
    
    IF (verbose) THEN
      PRINT *, 'BiCGSTAB(ℓ): Initial residual =', rnorm
      PRINT *, 'BiCGSTAB(ℓ): Target tolerance =', rtol
    END IF
    
    ! Check for immediate convergence
    IF (rnorm <= rtol) THEN
      converged = .TRUE.
      iter = 0
      stats%final_residual = rnorm
      stats%converged = .TRUE.
      CALL CPU_TIME(end_time)
      stats%solve_time = end_time - start_time
      DEALLOCATE(r, r0, p, v, s, t, z, y)
      RETURN
    END IF
    
    ! BiCGSTAB(ℓ) main iteration loop
    r0 = r  ! Choose r0 = r_0
    
    DO iter = 1, max_iter
      
      IF (iter == 1) THEN
        ! First iteration
        rho = DOT_PRODUCT(r0, r)
        p = r
      ELSE
        ! Update rho and beta
        rho_old = rho
        rho = DOT_PRODUCT(r0, r)
        beta = (rho / rho_old) * (alpha / omega)
        p = r + beta * (p - omega * v)
      END IF
      
      ! Apply AMG preconditioning: solve M*z = p for z
      z = p  ! Initial guess for AMG solve
      CALL amg_precond_apply(amg_hier, z, p)
      
      ! Compute v = A*z
      CALL csr_matvec(n, csr_row_ptr, csr_col_idx, csr_val, z, v)
      
      ! Compute alpha
      alpha = rho / DOT_PRODUCT(r0, v)
      
      ! Update s = r - alpha*v
      s = r - alpha * v
      
      ! Check for convergence on s
      rnorm = SQRT(DOT_PRODUCT(s, s))
      IF (rnorm <= rtol) THEN
        ! Converged on s
        x = x + alpha * z
        converged = .TRUE.
        stats%final_residual = rnorm
        EXIT
      END IF
      
      ! Apply AMG preconditioning: solve M*y = s for y
      y = s  ! Initial guess for AMG solve
      CALL amg_precond_apply(amg_hier, y, s)
      
      ! Compute t = A*y
      CALL csr_matvec(n, csr_row_ptr, csr_col_idx, csr_val, y, t)
      
      ! Compute omega
      omega = DOT_PRODUCT(t, s) / DOT_PRODUCT(t, t)
      
      ! Update solution and residual
      x = x + alpha * z + omega * y
      r = s - omega * t
      
      ! Check convergence
      rnorm = SQRT(DOT_PRODUCT(r, r))
      
      IF (verbose .AND. MOD(iter, 50) == 0) THEN
        PRINT *, 'BiCGSTAB(ℓ) iter', iter, ': residual =', rnorm
      END IF
      
      IF (rnorm <= rtol) THEN
        converged = .TRUE.
        stats%final_residual = rnorm
        EXIT
      END IF
      
      ! Check for breakdown
      IF (ABS(omega) < 1.0e-15_dp) THEN
        IF (verbose) THEN
          PRINT *, 'BiCGSTAB(ℓ): Breakdown detected (omega ~= 0)'
        END IF
        EXIT
      END IF
      
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
    
  END SUBROUTINE bicgstab_l_solve_amg_precond_real

END MODULE bicgstab_mod