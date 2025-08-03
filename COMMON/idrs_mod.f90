MODULE idrs_mod
  !> Exact line-by-line port of Julia IterativeSolvers.jl IDR(s) implementation
  !! Based on the algorithm from:
  !! - Sonneveld & van Gijzen (2008) "IDR(s): a family of simple and fast algorithms
  !!   for solving large nonsymmetric linear systems" SIAM J. Sci. Comput.
  !! - Van Gijzen & Sonneveld (2011) "Algorithm 913: An Elegant IDR(s) Variant that
  !!   Efficiently Exploits Bi-orthogonality Properties" ACM Trans. Math. Software
  !! 
  !! This is an exact port of the Julia implementation from IterativeSolvers.jl
  
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
    !> IDR(s) with AMG preconditioning - exact port of Julia implementation
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
    
    ! Local variables matching Julia implementation
    REAL(DP), ALLOCATABLE :: R(:), V(:), Q(:), Z(:)   ! Working vectors
    REAL(DP), ALLOCATABLE :: P(:,:)                   ! Shadow space (n x s)
    REAL(DP), ALLOCATABLE :: U(:,:), G(:,:)          ! Search directions (n x s)
    REAL(DP), ALLOCATABLE :: M(:,:)                   ! Small matrix (s x s)
    REAL(DP), ALLOCATABLE :: f(:), c(:)               ! Small vectors (s)
    
    ! Scalar variables
    REAL(DP) :: normR, normR0, tolb
    REAL(DP) :: omega_val, alpha, beta, gamma
    REAL(DP) :: rho, ns, nt, ts, angle
    INTEGER :: i, j, k, step, s
    LOGICAL :: verbose
    
    ! Initialize
    s = shadow_dim
    info = 0
    converged = .FALSE.
    iter = 0
    verbose = .FALSE.
    
    ! Allocate working arrays
    ALLOCATE(R(n), V(n), Q(n), Z(n))
    ALLOCATE(P(n, s))
    ALLOCATE(U(n, s), G(n, s))
    ALLOCATE(M(s, s))
    ALLOCATE(f(s), c(s))
    
    ! Initialize solution
    x = x_initial
    
    ! Line 115: R = C - A*X (compute initial residual)
    CALL sparse_matvec(n, row_ptr, col_idx, values, x, R)
    R = b - R
    
    ! Line 116: normR = norm(R)
    normR = SQRT(DOT_PRODUCT(R, R))
    normR0 = normR
    
    ! Line 117: tol = max(reltol * normR, abstol)
    tolb = tol * normR0
    
    IF (verbose) THEN
      WRITE(*,'(A)') '=== idrs ==='
      WRITE(*,'(A4,A,A4,A,A7)') 'iter', '    ', 'step', '    ', 'resnorm'
    END IF
    
    ! Line 130: Z = zero(C)
    Z = 0.0_DP
    
    ! Line 132: P = T[rand!(copy(C)) for k in 1:s]
    CALL initialize_P_random(P, n, s)
    
    ! Line 133-134: U = T[copy(Z) for k in 1:s], G = T[copy(Z) for k in 1:s]
    U = 0.0_DP
    G = 0.0_DP
    
    ! Line 135-136: Q = copy(Z), V = copy(Z)
    Q = 0.0_DP
    V = 0.0_DP
    
    ! Line 138: M = Matrix{eltype(C)}(I,s,s)
    M = 0.0_DP
    DO i = 1, s
      M(i, i) = 1.0_DP
    END DO
    
    ! Line 139-140: f = zeros(eltype(C),s), c = zeros(eltype(C),s)
    f = 0.0_DP
    c = 0.0_DP
    
    ! Line 142: omega::eltype(C) = 1
    omega_val = 1.0_DP
    
    ! Check initial residual for NaN
    IF (normR /= normR) THEN
      IF (verbose) THEN
        WRITE(*,'(A)') 'ERROR: Initial residual is NaN!'
      END IF
      info = -97
      RETURN
    END IF
    
    ! Main iteration loop
    step = 1
    DO WHILE (iter < max_iter)
      
      IF (verbose) THEN
        WRITE(*,'(A,I4,A,I2,A,ES12.5)') 'Iter ', iter+1, ', step ', step, ', normR=', normR
      END IF
      
      ! Check for NaN in residual
      IF (normR /= normR) THEN
        IF (verbose) THEN
          WRITE(*,'(A)') 'ERROR: Residual became NaN during iteration!'
        END IF
        info = -96
        EXIT
      END IF
      
      ! Line 167: if it.normR < it.tol || iter > it.maxiter
      IF (normR < tolb) THEN
        converged = .TRUE.
        EXIT
      END IF
      
      ! Line 176: if step in 1:s
      IF (step >= 1 .AND. step <= s) THEN
        
        ! Line 177-181: if step == 1, compute f[i] = dot(P[i], R)
        IF (step == 1) THEN
          DO i = 1, s
            f(i) = DOT_PRODUCT(P(:, i), R)
          END DO
        END IF
        
        k = step
        
        ! Line 186: c = LowerTriangular(M[k:s,k:s])\f[k:s]
        ! Solve lower triangular system
        DO i = k, s
          c(i-k+1) = f(i)
          DO j = k, i-1
            c(i-k+1) = c(i-k+1) - M(i, j) * c(j-k+1)
          END DO
          c(i-k+1) = c(i-k+1) / M(i, i)
        END DO
        
        ! Line 187-188: V .= c[1] .* G[k], Q .= c[1] .* U[k]
        IF (k <= s) THEN
          V = c(1) * G(:, k)
          Q = c(1) * U(:, k)
        END IF
        
        ! Line 190-193: for i = k+1:s
        DO i = k+1, s
          V = V + c(i-k+1) * G(:, i)
          Q = Q + c(i-k+1) * U(:, i)
        END DO
        
        ! Line 196: V .= R .- V
        V = R - V
        
        ! Line 199: ldiv!(Pl, V) - Apply preconditioner
        IF (verbose) THEN
          WRITE(*,'(A,I3,A,ES12.5)') '  Before AMG: step=', step, ', norm(V)=', SQRT(DOT_PRODUCT(V, V))
        END IF
        
        ! Apply AMG preconditioning: solve M*Z = V for Z (correct way)
        Z = 0.0_DP  ! Zero initial guess like GMRES
        CALL amg_precond_apply(amg_hier, Z, V)
        V = Z  ! Copy preconditioned result back to V
        
        IF (verbose) THEN
          WRITE(*,'(A,I3,A,ES12.5)') '  After AMG: step=', step, ', norm(V)=', SQRT(DOT_PRODUCT(V, V))
          ! Check for NaN in V
          IF (ANY(V /= V)) THEN
            WRITE(*,'(A)') '  ERROR: NaN detected in V after AMG preconditioning!'
            info = -99
            RETURN
          END IF
        END IF
        
        ! Line 201: U[k] .= Q .+ it.omega .* V
        U(:, k) = Q + omega_val * V
        
        ! Line 202: mul!(G[k], A, U[k])
        CALL sparse_matvec(n, row_ptr, col_idx, values, U(:, k), G(:, k))
        
        ! Line 206-210: Bi-orthogonalise the new basis vectors
        DO i = 1, k-1
          alpha = DOT_PRODUCT(P(:, i), G(:, k)) / M(i, i)
          G(:, k) = G(:, k) - alpha * G(:, i)
          U(:, k) = U(:, k) - alpha * U(:, i)
        END DO
        
        ! Line 214-216: New column of M = P'*G
        DO i = k, s
          M(i, k) = DOT_PRODUCT(P(:, i), G(:, k))
        END DO
        
        ! Line 220-222: Make r orthogonal to q_i
        IF (verbose) THEN
          WRITE(*,'(A,I3,A,ES12.5,A,ES12.5)') '  Division: k=', k, ', f(k)=', f(k), ', M(k,k)=', M(k, k)
        END IF
        
        IF (ABS(M(k, k)) < EPS_DP) THEN
          IF (verbose) THEN
            WRITE(*,'(A)') '  ERROR: M(k,k) is near zero - potential division by zero!'
          END IF
          info = -98
          RETURN
        END IF
        
        beta = f(k) / M(k, k)
        
        IF (verbose) THEN
          WRITE(*,'(A,ES12.5)') '  beta=', beta
        END IF
        
        R = R - beta * G(:, k)
        x = x + beta * U(:, k)
        
        ! Line 224: it.normR = norm(R)
        normR = SQRT(DOT_PRODUCT(R, R))
        residual_norm = normR
        
        ! Line 235-237: if k < s, update f
        IF (k < s) THEN
          DO i = k+1, s
            f(i) = f(i) - beta * M(i, k)
          END DO
        END IF
        
        step = step + 1
        
      ELSE IF (step == s + 1) THEN
        ! Line 239: elseif step == s + 1
        
        ! Line 243: copyto!(V, R)
        V = R
        
        ! Line 246: ldiv!(Pl, V) - Apply preconditioner
        IF (verbose) THEN
          WRITE(*,'(A,ES12.5)') '  Before 2nd AMG: norm(V)=', SQRT(DOT_PRODUCT(V, V))
        END IF
        
        ! Apply AMG preconditioning: solve M*Z = V for Z (correct way)
        Z = 0.0_DP  ! Zero initial guess like GMRES
        CALL amg_precond_apply(amg_hier, Z, V)
        V = Z  ! Copy preconditioned result back to V
        
        IF (verbose) THEN
          WRITE(*,'(A,ES12.5)') '  After 2nd AMG: norm(V)=', SQRT(DOT_PRODUCT(V, V))
          ! Check for NaN in V
          IF (ANY(V /= V)) THEN
            WRITE(*,'(A)') '  ERROR: NaN detected in V after 2nd AMG preconditioning!'
            info = -99
            RETURN
          END IF
        END IF
        
        ! Line 248: mul!(Q, A, V)
        CALL sparse_matvec(n, row_ptr, col_idx, values, V, Q)
        
        ! Line 249: it.omega = omega(Q, R)
        ! Julia omega function implementation
        angle = SQRT(2.0_DP) / 2.0_DP
        ns = SQRT(DOT_PRODUCT(R, R))
        nt = SQRT(DOT_PRODUCT(Q, Q))
        ts = DOT_PRODUCT(Q, R)
        
        IF (verbose) THEN
          WRITE(*,'(A,ES12.5,A,ES12.5,A,ES12.5)') '  Omega calc: ns=', ns, ', nt=', nt, ', ts=', ts
        END IF
        
        IF (nt > EPS_DP .AND. ns > EPS_DP) THEN
          rho = ABS(ts) / (nt * ns)
          omega_val = ts / (nt * nt)
          IF (verbose) THEN
            WRITE(*,'(A,ES12.5,A,ES12.5)') '  rho=', rho, ', omega_val(initial)=', omega_val
          END IF
          IF (rho < angle) THEN
            omega_val = omega_val * angle / rho
            IF (verbose) THEN
              WRITE(*,'(A,ES12.5)') '  omega_val(adjusted)=', omega_val
            END IF
          END IF
        ELSE
          IF (verbose) THEN
            WRITE(*,'(A)') '  WARNING: nt or ns too small, keeping omega=1'
          END IF
        END IF
        
        ! Line 250-251: R .-= it.omega .* Q, X .+= it.omega .* V
        R = R - omega_val * Q
        x = x + omega_val * V
        
        ! Line 253: it.normR = norm(R)
        normR = SQRT(DOT_PRODUCT(R, R))
        residual_norm = normR
        
        ! Line 264: nextstep = 1
        step = 1
      END IF
      
      iter = iter + 1
      
      IF (verbose) THEN
        WRITE(*,'(I3,A,I3,A,ES12.5)') iter, '    ', step-1, '    ', normR
      END IF
      
    END DO  ! Main loop
    
    ! Final status
    IF (verbose) THEN
      WRITE(*,*)
    END IF
    
    ! Set return code
    IF (.NOT. converged .AND. iter >= max_iter) info = 1
    
    ! Clean up
    DEALLOCATE(R, V, Q, Z)
    DEALLOCATE(P, U, G, M, f, c)
    
  END SUBROUTINE idrs_solve_amg_preconditioned
  
  SUBROUTINE initialize_P_random(P, n, s)
    !> Initialize shadow space P with random values (matching Julia's rand!)
    REAL(DP), INTENT(OUT) :: P(:,:)
    INTEGER, INTENT(IN) :: n, s
    
    INTEGER :: i, j, k
    INTEGER, ALLOCATABLE :: seed(:)
    
    ! Initialize random seed for reproducibility
    CALL RANDOM_SEED(SIZE=k)
    ALLOCATE(seed(k))
    seed = 12345
    CALL RANDOM_SEED(PUT=seed)
    DEALLOCATE(seed)
    
    ! Julia's rand! generates uniform random numbers in [0,1]
    CALL RANDOM_NUMBER(P)
    
  END SUBROUTINE initialize_P_random
  
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

END MODULE idrs_mod