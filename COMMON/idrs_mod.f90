MODULE idrs_mod
  !> Exact port of Julia IterativeSolvers.jl IDR(s) implementation
  !! Based on the algorithm from:
  !! - Sonneveld & van Gijzen (2008) "IDR(s): a family of simple and fast algorithms
  !!   for solving large nonsymmetric linear systems" SIAM J. Sci. Comput.
  !! - Van Gijzen & Sonneveld (2011) "Algorithm 913: An Elegant IDR(s) Variant that
  !!   Efficiently Exploits Bi-orthogonality Properties" ACM Trans. Math. Software
  !! 
  !! This implementation uses left preconditioning: solves P^(-1) A x = P^(-1) b
  
  USE nrtype, ONLY: I4B, DP
  USE amg_types_mod, ONLY: amg_hierarchy
  USE amg_precond_mod, ONLY: amg_precond_apply
  IMPLICIT NONE
  PRIVATE
  
  ! Parameters
  REAL(DP), PARAMETER :: EPS_DP = EPSILON(1.0_DP)
  
  ! Public interface
  PUBLIC :: idrs_solve_amg_preconditioned

CONTAINS

  SUBROUTINE idrs_solve_amg_preconditioned(workspace, n, row_ptr, col_idx, values, &
                                          b, x_initial, max_iter, tol, amg_hier, &
                                          shadow_dim, x, iter, residual_norm, &
                                          converged, info)
    !> IDR(s) with left AMG preconditioning
    !! Solves: A x = b  using left-preconditioned IDR(s)
    !! Actually solves: P^(-1) A x = P^(-1) b where P is the AMG preconditioner
    TYPE(*), INTENT(INOUT) :: workspace  ! For interface compatibility
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(IN) :: row_ptr(:), col_idx(:)
    REAL(DP), INTENT(IN) :: values(:), b(:), x_initial(:)
    INTEGER, INTENT(IN) :: max_iter, shadow_dim
    REAL(DP), INTENT(IN) :: tol
    TYPE(amg_hierarchy), INTENT(INOUT) :: amg_hier
    REAL(DP), INTENT(OUT) :: x(:)
    INTEGER, INTENT(OUT) :: iter
    REAL(DP), INTENT(OUT) :: residual_norm
    LOGICAL, INTENT(OUT) :: converged
    INTEGER, INTENT(OUT) :: info
    
    ! Local variables
    REAL(DP), ALLOCATABLE :: r(:), t(:), v(:), w(:)  ! Working vectors
    REAL(DP), ALLOCATABLE :: P(:,:)                  ! Shadow space (n x s)
    REAL(DP), ALLOCATABLE :: G(:,:), U(:,:)         ! Direction matrices (n x s)
    REAL(DP), ALLOCATABLE :: M(:,:)                  ! Small matrix (s x s)
    REAL(DP), ALLOCATABLE :: f(:), c(:)              ! Small vectors (s)
    
    ! Scalar variables
    REAL(DP) :: normr, normr0, tolb, omega
    REAL(DP) :: alpha, beta, rho
    INTEGER :: i, j, k
    LOGICAL :: verbose
    
    ! Enable verbose output for debugging
    verbose = .FALSE.
    
    ! Initialize
    info = 0
    converged = .FALSE.
    iter = 0
    
    ! Allocate working arrays
    ALLOCATE(r(n), t(n), v(n), w(n))
    ALLOCATE(P(n, shadow_dim))
    ALLOCATE(G(n, shadow_dim), U(n, shadow_dim))
    ALLOCATE(M(shadow_dim, shadow_dim))
    ALLOCATE(f(shadow_dim), c(shadow_dim))
    
    ! Initialize solution
    x = x_initial
    
    ! Compute initial residual: r = b - A*x
    CALL sparse_matvec(n, row_ptr, col_idx, values, x, v)
    r = b - v
    
    ! Initial residual norm (using unpreconditioned residual for true error)
    normr = SQRT(DOT_PRODUCT(r, r))
    normr0 = normr
    residual_norm = normr
    
    ! Set tolerance
    tolb = tol * normr0
    
    IF (verbose) THEN
      WRITE(*,'(A,ES12.5)') 'IDR(s): Initial residual = ', normr
    END IF
    
    ! Check for immediate convergence
    IF (normr <= tolb) THEN
      converged = .TRUE.
      RETURN
    END IF
    
    ! Generate random orthonormal shadow space P
    CALL initialize_P_random(P, n, shadow_dim)
    
    ! Initialize
    U = 0.0_DP
    G = 0.0_DP
    M = 0.0_DP
    omega = 1.0_DP
    
    ! Main IDR(s) iteration
    DO WHILE (iter < max_iter .AND. .NOT. converged)
      
      ! Apply preconditioner: v = P^(-1) * r
      v = 0.0_DP
      CALL amg_precond_apply(amg_hier, v, r)
      
      ! Compute f = P^T * v (where P is shadow space, not preconditioner!)
      DO i = 1, shadow_dim
        f(i) = DOT_PRODUCT(P(:, i), v)
      END DO
      
      ! s steps of IDR
      DO k = 1, shadow_dim
        
        ! Solve small system for c
        IF (k > 1) THEN
          CALL solve_small_system(M(1:k-1, 1:k-1), f(1:k-1), c(1:k-1), k-1)
          
          ! Update v = v - sum(G(:,j) * c(j))
          DO j = 1, k-1
            v = v - c(j) * G(:, j)
          END DO
          
          ! Update U(:,k) = omega * v + sum(U(:,j) * c(j))
          U(:, k) = omega * v
          DO j = 1, k-1
            U(:, k) = U(:, k) + c(j) * U(:, j)
          END DO
        ELSE
          ! First iteration: no system to solve
          U(:, k) = omega * v
        END IF
        
        ! Compute G(:,k) = P^(-1) * A * U(:,k)
        ! First: w = A * U(:,k)
        CALL sparse_matvec(n, row_ptr, col_idx, values, U(:, k), w)
        ! Then: G(:,k) = P^(-1) * w
        G(:, k) = 0.0_DP
        CALL amg_precond_apply(amg_hier, G(:, k), w)
        
        ! Bi-orthogonalize: M(i,k) = P(:,i)^T * G(:,k)
        DO i = 1, shadow_dim
          M(i, k) = DOT_PRODUCT(P(:, i), G(:, k))
        END DO
        
        ! Check for breakdown
        IF (ABS(M(k, k)) < EPS_DP * normr0) THEN
          IF (verbose) THEN
            WRITE(*,'(A,I3,A,ES12.5)') 'IDR(s) breakdown: M(', k, ',', k, ') = ', M(k, k)
          END IF
          info = 2
          EXIT
        END IF
        
        ! Compute alpha and update solution/residual
        alpha = f(k) / M(k, k)
        x = x + alpha * U(:, k)
        r = r - alpha * w  ! Note: update with unpreconditioned w = A*U(:,k)
        
        ! Update f for next iterations
        IF (k < shadow_dim) THEN
          DO i = k+1, shadow_dim
            f(i) = f(i) - alpha * M(i, k)
          END DO
        END IF
        
        ! Update iteration count
        iter = iter + 1
        
        ! Check convergence
        normr = SQRT(DOT_PRODUCT(r, r))
        residual_norm = normr
        
        IF (normr <= tolb) THEN
          converged = .TRUE.
          EXIT
        END IF
        
        IF (verbose .AND. MOD(iter, 10) == 0) THEN
          WRITE(*,'(A,I5,A,ES12.5)') 'IDR(s) iter ', iter, ': residual = ', normr
        END IF
        
      END DO  ! k loop
      
      IF (converged .OR. info /= 0) EXIT
      
      ! Additional (s+1)-th step
      IF (iter < max_iter) THEN
        
        ! Solve M * c = f
        CALL solve_small_system(M, f, c, shadow_dim)
        
        ! Update solution and residual
        ! v = sum(U(:,j) * c(j))
        v = 0.0_DP
        DO j = 1, shadow_dim
          v = v + c(j) * U(:, j)
        END DO
        x = x + v
        
        ! w = sum(G(:,j) * c(j))  (this is P^(-1)*A*v)
        w = 0.0_DP
        DO j = 1, shadow_dim
          w = w + c(j) * G(:, j)
        END DO
        
        ! Need unpreconditioned A*v for residual update
        CALL sparse_matvec(n, row_ptr, col_idx, values, v, t)
        r = r - t
        
        ! Apply preconditioner to new residual
        v = 0.0_DP
        CALL amg_precond_apply(amg_hier, v, r)
        
        ! Compute omega using Julia's approach
        ! w = A * v
        CALL sparse_matvec(n, row_ptr, col_idx, values, v, w)
        ! t = P^(-1) * A * v
        t = 0.0_DP
        CALL amg_precond_apply(amg_hier, t, w)
        
        ! Julia's omega function
        CALL compute_omega_julia(t, r, omega)
        
        ! Update solution and residual
        x = x + omega * v
        r = r - omega * w  ! w = A*v (unpreconditioned)
        
        iter = iter + 1
        
        ! Check convergence
        normr = SQRT(DOT_PRODUCT(r, r))
        residual_norm = normr
        
        IF (normr <= tolb) THEN
          converged = .TRUE.
        END IF
        
        IF (verbose .AND. MOD(iter, 10) == 0) THEN
          WRITE(*,'(A,I5,A,ES12.5,A,F8.4)') 'IDR(s) iter ', iter, &
                 ': residual = ', normr, ', omega = ', omega
        END IF
        
      END IF
      
    END DO  ! Main loop
    
    ! Final status
    IF (verbose) THEN
      IF (converged) THEN
        WRITE(*,'(A,I5,A,ES12.5)') 'IDR(s) converged in ', iter, &
               ' iterations, final residual = ', residual_norm
      ELSE
        WRITE(*,'(A,I5,A,ES12.5)') 'IDR(s) did not converge after ', iter, &
               ' iterations, final residual = ', residual_norm
      END IF
    END IF
    
    ! Set return code
    IF (.NOT. converged .AND. iter >= max_iter) info = 1
    
    ! Clean up
    DEALLOCATE(r, t, v, w)
    DEALLOCATE(P, G, U, M, f, c)
    
  END SUBROUTINE idrs_solve_amg_preconditioned
  
  SUBROUTINE initialize_P_random(P, n, s)
    !> Initialize shadow space P with random orthonormal vectors
    REAL(DP), INTENT(OUT) :: P(:,:)
    INTEGER, INTENT(IN) :: n, s
    
    INTEGER :: i, j, k
    INTEGER, ALLOCATABLE :: seed(:)
    REAL(DP) :: norm_val, dot_val
    
    ! Initialize random seed
    CALL RANDOM_SEED(SIZE=k)
    ALLOCATE(seed(k))
    seed = 12345
    CALL RANDOM_SEED(PUT=seed)
    DEALLOCATE(seed)
    
    ! Fill with random values
    CALL RANDOM_NUMBER(P)
    
    ! Shift to [-1, 1]
    P = 2.0_DP * P - 1.0_DP
    
    ! Orthonormalize using modified Gram-Schmidt
    DO j = 1, s
      ! Orthogonalize against previous vectors
      DO k = 1, j-1
        dot_val = DOT_PRODUCT(P(:, j), P(:, k))
        P(:, j) = P(:, j) - dot_val * P(:, k)
      END DO
      
      ! Normalize
      norm_val = SQRT(DOT_PRODUCT(P(:, j), P(:, j)))
      IF (norm_val > EPS_DP) THEN
        P(:, j) = P(:, j) / norm_val
      ELSE
        ! If dependent, set to unit vector
        P(:, j) = 0.0_DP
        IF (j <= n) P(j, j) = 1.0_DP
      END IF
    END DO
    
  END SUBROUTINE initialize_P_random
  
  SUBROUTINE solve_small_system(A, b, x, n)
    !> Solve small dense system Ax = b using Gaussian elimination with pivoting
    REAL(DP), INTENT(IN) :: A(:,:), b(:)
    REAL(DP), INTENT(OUT) :: x(:)
    INTEGER, INTENT(IN) :: n
    
    REAL(DP), ALLOCATABLE :: LU(:,:), work(:)
    INTEGER :: i, j, k, pivot
    REAL(DP) :: factor, temp
    
    IF (n == 0) RETURN
    
    ! Allocate working arrays
    ALLOCATE(LU(n, n), work(n))
    
    ! Copy to working arrays
    LU = A(1:n, 1:n)
    work = b(1:n)
    
    ! Forward elimination with partial pivoting
    DO k = 1, n-1
      ! Find pivot
      pivot = k
      DO i = k+1, n
        IF (ABS(LU(i, k)) > ABS(LU(pivot, k))) pivot = i
      END DO
      
      ! Swap rows if needed
      IF (pivot /= k) THEN
        DO j = k, n
          temp = LU(k, j)
          LU(k, j) = LU(pivot, j)
          LU(pivot, j) = temp
        END DO
        temp = work(k)
        work(k) = work(pivot)
        work(pivot) = temp
      END IF
      
      ! Eliminate column
      IF (ABS(LU(k, k)) > EPS_DP) THEN
        DO i = k+1, n
          factor = LU(i, k) / LU(k, k)
          DO j = k+1, n
            LU(i, j) = LU(i, j) - factor * LU(k, j)
          END DO
          work(i) = work(i) - factor * work(k)
        END DO
      END IF
    END DO
    
    ! Back substitution
    x = 0.0_DP
    DO i = n, 1, -1
      IF (ABS(LU(i, i)) > EPS_DP) THEN
        x(i) = work(i)
        DO j = i+1, n
          x(i) = x(i) - LU(i, j) * x(j)
        END DO
        x(i) = x(i) / LU(i, i)
      END IF
    END DO
    
    DEALLOCATE(LU, work)
    
  END SUBROUTINE solve_small_system
  
  SUBROUTINE sparse_matvec(n, row_ptr, col_idx, values, x, y)
    !> Sparse matrix-vector product: y = A * x
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(IN) :: row_ptr(:), col_idx(:)
    REAL(DP), INTENT(IN) :: values(:), x(:)
    REAL(DP), INTENT(OUT) :: y(:)
    
    INTEGER :: i, j
    
    y = 0.0_DP
    DO i = 1, n
      DO j = row_ptr(i), row_ptr(i+1) - 1
        y(i) = y(i) + values(j) * x(col_idx(j))
      END DO
    END DO
    
  END SUBROUTINE sparse_matvec
  
  SUBROUTINE compute_omega_julia(t, s, omega)
    !> Compute omega exactly as in Julia IterativeSolvers.jl
    REAL(DP), INTENT(IN) :: t(:), s(:)
    REAL(DP), INTENT(OUT) :: omega
    
    REAL(DP), PARAMETER :: angle = 0.7071067811865476_DP  ! sqrt(2)/2
    REAL(DP) :: ns, nt, ts, rho
    
    ns = SQRT(DOT_PRODUCT(s, s))
    nt = SQRT(DOT_PRODUCT(t, t))
    ts = DOT_PRODUCT(t, s)
    
    IF (nt > EPS_DP .AND. ns > EPS_DP) THEN
      rho = ABS(ts) / (nt * ns)
      omega = ts / (nt * nt)
      
      IF (rho < angle) THEN
        omega = omega * angle / rho
      END IF
    ELSE
      omega = 0.0_DP
    END IF
    
  END SUBROUTINE compute_omega_julia

END MODULE idrs_mod