MODULE idrs_mod
  !> Exact replication of Julia IterativeSolvers.jl IDR(s) implementation
  !! Replicates idrs.jl algorithms for Induced Dimension Reduction
  !! MIT License compliance from IterativeSolvers.jl
  
  USE nrtype, ONLY: I4B, DP
  IMPLICIT NONE
  PRIVATE
  
  ! Exact replication of Julia IDRSIterable structure
  TYPE, PUBLIC :: idrs_workspace
    REAL(DP), ALLOCATABLE :: X(:)        ! Solution vector
    REAL(DP), ALLOCATABLE :: R(:)        ! Residual vector
    REAL(DP), ALLOCATABLE :: P(:,:)      ! Shadow space vectors s x n
    REAL(DP), ALLOCATABLE :: U(:,:)      ! Search directions s x n
    REAL(DP), ALLOCATABLE :: G(:,:)      ! A*U vectors s x n
    REAL(DP), ALLOCATABLE :: Q(:)        ! Workspace vector
    REAL(DP), ALLOCATABLE :: V(:)        ! Workspace vector
    REAL(DP), ALLOCATABLE :: Z(:)        ! Zero vector
    REAL(DP), ALLOCATABLE :: M(:,:)      ! s x s matrix
    REAL(DP), ALLOCATABLE :: f(:)        ! s vector
    REAL(DP), ALLOCATABLE :: c(:)        ! s vector
    
    ! Smoothing vectors (optional)
    REAL(DP), ALLOCATABLE :: X_s(:)      ! Smoothed solution
    REAL(DP), ALLOCATABLE :: R_s(:)      ! Smoothed residual
    REAL(DP), ALLOCATABLE :: T_s(:)      ! Smoothing workspace
    
    INTEGER :: n                         ! Problem size
    INTEGER :: s                         ! Shadow space dimension
    INTEGER :: maxiter                   ! Maximum iterations
    INTEGER :: iter                      ! Current iteration
    INTEGER :: step                      ! Current step (1 to s+1)
    REAL(DP) :: omega                    ! Relaxation parameter
    REAL(DP) :: normR                    ! Current residual norm
    REAL(DP) :: tol                      ! Tolerance
    REAL(DP) :: abstol, reltol           ! Absolute and relative tolerance
    LOGICAL :: converged                 ! Convergence flag
    LOGICAL :: smoothing                 ! Enable smoothing
  END TYPE idrs_workspace
  
  PUBLIC :: idrs_solve_real, idrs_solve_structured, create_idrs_workspace, &
            destroy_idrs_workspace
  
CONTAINS

  SUBROUTINE idrs_solve_real(A, b, x, s, stats)
    REAL(DP), INTENT(IN) :: A(:,:)
    REAL(DP), INTENT(IN) :: b(:)
    REAL(DP), INTENT(INOUT) :: x(:)
    INTEGER, INTENT(IN), OPTIONAL :: s
    INTEGER, INTENT(OUT), OPTIONAL :: stats
    
    TYPE(idrs_workspace) :: workspace
    INTEGER :: shadow_dim, iter, info
    REAL(DP) :: residual_norm
    LOGICAL :: converged
    
    ! Set default shadow space dimension
    shadow_dim = 8
    IF (PRESENT(s)) shadow_dim = s
    
    ! Create workspace and solve
    CALL create_idrs_workspace(workspace, SIZE(b), shadow_dim)
    CALL idrs_solve_structured(workspace, A, b, x, 1000, 1.0E-8_DP, &
                               iter, residual_norm, converged, info)
    CALL destroy_idrs_workspace(workspace)
    
    IF (PRESENT(stats)) stats = info
  END SUBROUTINE idrs_solve_real
  
  SUBROUTINE create_idrs_workspace(workspace, n, s)
    TYPE(idrs_workspace), INTENT(OUT) :: workspace
    INTEGER, INTENT(IN) :: n, s
    
    INTEGER :: i, j
    
    ! Set dimensions
    workspace%n = n
    workspace%s = s
    
    ! Allocate main vectors
    ALLOCATE(workspace%X(n))
    ALLOCATE(workspace%R(n))
    ALLOCATE(workspace%Q(n))
    ALLOCATE(workspace%V(n))
    ALLOCATE(workspace%Z(n))
    
    ! Allocate matrix arrays
    ALLOCATE(workspace%P(n, s))
    ALLOCATE(workspace%U(n, s))
    ALLOCATE(workspace%G(n, s))
    
    ! Allocate small matrices/vectors
    ALLOCATE(workspace%M(s, s))
    ALLOCATE(workspace%f(s))
    ALLOCATE(workspace%c(s))
    
    ! Initialize shadow space P with random vectors (Julia: rand!(copy(C)))
    ! Using simple pseudo-random initialization
    DO j = 1, s
      DO i = 1, n
        workspace%P(i, j) = SIN(REAL(i*j, DP)) + COS(REAL(i+j, DP))
      END DO
      ! Normalize
      workspace%P(:, j) = workspace%P(:, j) / NORM2(workspace%P(:, j))
    END DO
    
    ! Initialize M as identity matrix (Julia: Matrix{eltype(C)}(I,s,s))
    workspace%M = 0.0_DP
    DO i = 1, s
      workspace%M(i, i) = 1.0_DP
    END DO
    
    ! Initialize other arrays
    workspace%U = 0.0_DP
    workspace%G = 0.0_DP
    workspace%Z = 0.0_DP
    workspace%f = 0.0_DP
    workspace%c = 0.0_DP
    
    ! Initialize state
    workspace%omega = 1.0_DP
    workspace%iter = 1
    workspace%step = 1
    workspace%converged = .FALSE.
    workspace%smoothing = .FALSE.
    
  END SUBROUTINE create_idrs_workspace
  
  SUBROUTINE destroy_idrs_workspace(workspace)
    TYPE(idrs_workspace), INTENT(INOUT) :: workspace
    
    ! Deallocate all arrays
    IF (ALLOCATED(workspace%X)) DEALLOCATE(workspace%X)
    IF (ALLOCATED(workspace%R)) DEALLOCATE(workspace%R)
    IF (ALLOCATED(workspace%P)) DEALLOCATE(workspace%P)
    IF (ALLOCATED(workspace%U)) DEALLOCATE(workspace%U)
    IF (ALLOCATED(workspace%G)) DEALLOCATE(workspace%G)
    IF (ALLOCATED(workspace%Q)) DEALLOCATE(workspace%Q)
    IF (ALLOCATED(workspace%V)) DEALLOCATE(workspace%V)
    IF (ALLOCATED(workspace%Z)) DEALLOCATE(workspace%Z)
    IF (ALLOCATED(workspace%M)) DEALLOCATE(workspace%M)
    IF (ALLOCATED(workspace%f)) DEALLOCATE(workspace%f)
    IF (ALLOCATED(workspace%c)) DEALLOCATE(workspace%c)
    IF (ALLOCATED(workspace%X_s)) DEALLOCATE(workspace%X_s)
    IF (ALLOCATED(workspace%R_s)) DEALLOCATE(workspace%R_s)
    IF (ALLOCATED(workspace%T_s)) DEALLOCATE(workspace%T_s)
  END SUBROUTINE destroy_idrs_workspace
  
  SUBROUTINE idrs_solve_structured(workspace, A, b, x, max_iter, tol, &
                                  iter, residual_norm, converged, info)
    TYPE(idrs_workspace), INTENT(INOUT) :: workspace
    REAL(DP), INTENT(IN) :: A(:,:), b(:)
    REAL(DP), INTENT(INOUT) :: x(:)
    INTEGER, INTENT(IN) :: max_iter
    REAL(DP), INTENT(IN) :: tol
    INTEGER, INTENT(OUT) :: iter
    REAL(DP), INTENT(OUT) :: residual_norm
    LOGICAL, INTENT(OUT) :: converged
    INTEGER, INTENT(OUT) :: info
    
    ! Exact replication of Julia idrs_method! function
    INTEGER :: i
    
    ! Initialize workspace
    workspace%X = x
    workspace%maxiter = max_iter
    workspace%abstol = 0.0_DP
    workspace%reltol = tol
    
    ! Julia: R = C - A*X (line 115)
    CALL matrix_vector_product(A, workspace%X, workspace%R)
    workspace%R = b - workspace%R
    
    ! Julia: normR = norm(R) (line 116)
    workspace%normR = NORM2(workspace%R)
    
    ! Julia: tol = max(reltol * normR, abstol) (line 117)
    workspace%tol = MAX(workspace%reltol * workspace%normR, workspace%abstol)
    
    ! Check initial convergence
    workspace%converged = (workspace%normR < workspace%tol)
    
    info = 0
    iter = 0
    
    ! Main IDR(s) iteration loop
    DO WHILE (iter < max_iter .AND. .NOT. workspace%converged)
      
      ! Julia iterate function (lines 163-272)
      CALL idrs_iterate(workspace, A)
      iter = iter + 1
      
    END DO
    
    ! Output results  
    x = workspace%X
    residual_norm = workspace%normR
    converged = workspace%converged
    
    IF (iter >= max_iter .AND. .NOT. converged) info = 1
    
  END SUBROUTINE idrs_solve_structured
  
  SUBROUTINE idrs_iterate(workspace, A)
    ! Exact replication of Julia iterate function (lines 163-272)
    TYPE(idrs_workspace), INTENT(INOUT) :: workspace
    REAL(DP), INTENT(IN) :: A(:,:)
    
    INTEGER :: i, j, k
    REAL(DP) :: alpha, beta, gamma
    
    ! Check convergence (Julia lines 167-174)
    IF (workspace%normR < workspace%tol) THEN
      workspace%converged = .TRUE.
      RETURN
    END IF
    
    ! Julia: if step in 1:s (line 176)
    IF (workspace%step >= 1 .AND. workspace%step <= workspace%s) THEN
      
      ! Julia: if step == 1 (line 177)
      IF (workspace%step == 1) THEN
        ! Julia: f[i] = dot(P[i], R) (line 179)
        DO i = 1, workspace%s
          workspace%f(i) = DOT_PRODUCT(workspace%P(:, i), workspace%R)
        END DO
      END IF
      
      k = workspace%step
      
      ! Julia: Solve small system (lines 186-193)
      ! c = LowerTriangular(M[k:s,k:s])\f[k:s]
      CALL solve_lower_triangular(workspace%M(k:workspace%s, k:workspace%s), &
                                  workspace%f(k:workspace%s), &
                                  workspace%c(1:workspace%s-k+1))
      
      ! Julia: V .= c[1] .* G[k] (line 187)
      workspace%V = workspace%c(1) * workspace%G(:, k)
      workspace%Q = workspace%c(1) * workspace%U(:, k)
      
      ! Julia: for i = k+1:s (lines 190-193)
      DO i = k+1, workspace%s
        workspace%V = workspace%V + workspace%c(i-k+1) * workspace%G(:, i)
        workspace%Q = workspace%Q + workspace%c(i-k+1) * workspace%U(:, i)
      END DO
      
      ! Julia: V .= R .- V (line 196)
      workspace%V = workspace%R - workspace%V
      
      ! No preconditioning for now (Julia line 199: ldiv!(Pl, V))
      
      ! Julia: U[k] .= Q .+ omega .* V (line 201)
      workspace%U(:, k) = workspace%Q + workspace%omega * workspace%V
      
      ! Julia: mul!(G[k], A, U[k]) (line 202)
      CALL matrix_vector_product(A, workspace%U(:, k), workspace%G(:, k))
      
      ! Julia: Bi-orthogonalise (lines 206-211)
      DO i = 1, k-1
        alpha = DOT_PRODUCT(workspace%P(:, i), workspace%G(:, k)) / workspace%M(i, i)
        workspace%G(:, k) = workspace%G(:, k) - alpha * workspace%G(:, i)
        workspace%U(:, k) = workspace%U(:, k) - alpha * workspace%U(:, i)
      END DO
      
      ! Julia: New column of M (lines 214-216)
      DO i = k, workspace%s
        workspace%M(i, k) = DOT_PRODUCT(workspace%P(:, i), workspace%G(:, k))
      END DO
      
      ! Julia: Make r orthogonal to q_i (lines 220-223)
      beta = workspace%f(k) / workspace%M(k, k)
      workspace%R = workspace%R - beta * workspace%G(:, k)
      workspace%X = workspace%X + beta * workspace%U(:, k)
      
      ! Julia: normR = norm(R) (line 224)
      workspace%normR = NORM2(workspace%R)
      
      ! Julia: Update f for next steps (lines 235-237)
      IF (k < workspace%s) THEN
        DO i = k+1, workspace%s
          workspace%f(i) = workspace%f(i) - beta * workspace%M(i, k)
        END DO
      END IF
      
      workspace%step = workspace%step + 1
      
    ! Julia: elseif step == s + 1 (line 239)
    ELSE IF (workspace%step == workspace%s + 1) THEN
      
      ! Julia: copyto!(V, R) (line 243)
      workspace%V = workspace%R
      
      ! No preconditioning (Julia line 246: ldiv!(Pl, V))
      
      ! Julia: mul!(Q, A, V) (line 248)
      CALL matrix_vector_product(A, workspace%V, workspace%Q)
      
      ! Julia: omega = omega(Q, R) (line 249)
      workspace%omega = compute_omega(workspace%Q, workspace%R)
      
      ! Julia: R .-= omega .* Q (line 250)
      workspace%R = workspace%R - workspace%omega * workspace%Q
      workspace%X = workspace%X + workspace%omega * workspace%V
      
      ! Julia: normR = norm(R) (line 253)
      workspace%normR = NORM2(workspace%R)
      
      workspace%step = 1
      
    END IF
    
    ! Check convergence
    workspace%converged = (workspace%normR <= workspace%tol)
    
  END SUBROUTINE idrs_iterate
  
  SUBROUTINE solve_lower_triangular(L, b, x)
    ! Solve Lx = b where L is lower triangular
    REAL(DP), INTENT(IN) :: L(:,:), b(:)
    REAL(DP), INTENT(OUT) :: x(:)
    
    INTEGER :: i, j, n
    
    n = SIZE(b)
    
    ! Forward substitution
    DO i = 1, n
      x(i) = b(i)
      DO j = 1, i-1
        x(i) = x(i) - L(i, j) * x(j)
      END DO
      x(i) = x(i) / L(i, i)
    END DO
    
  END SUBROUTINE solve_lower_triangular
  
  FUNCTION compute_omega(t, s) RESULT(omega_val)
    ! Exact replication of Julia omega function (lines 70-81)
    REAL(DP), INTENT(IN) :: t(:), s(:)
    REAL(DP) :: omega_val
    
    REAL(DP), PARAMETER :: angle = SQRT(2.0_DP) / 2.0_DP
    REAL(DP) :: ns, nt, ts, rho
    
    ! Julia: ns = norm(s), nt = norm(t), ts = dot(t,s)
    ns = NORM2(s)
    nt = NORM2(t)
    ts = DOT_PRODUCT(t, s)
    
    ! Julia: rho = abs(ts/(nt*ns))
    rho = ABS(ts / (nt * ns))
    
    ! Julia: omega = ts/(nt*nt)
    omega_val = ts / (nt * nt)
    
    ! Julia: if rho < angle
    IF (rho < angle) THEN
      omega_val = omega_val * angle / rho
    END IF
    
  END FUNCTION compute_omega
  
  SUBROUTINE matrix_vector_product(A, x, Ax)
    ! Standard matrix-vector product
    REAL(DP), INTENT(IN) :: A(:,:), x(:)
    REAL(DP), INTENT(OUT) :: Ax(:)
    
    INTEGER :: i, j
    
    DO i = 1, SIZE(A, 1)
      Ax(i) = 0.0_DP
      DO j = 1, SIZE(A, 2)
        Ax(i) = Ax(i) + A(i, j) * x(j)
      END DO
    END DO
    
  END SUBROUTINE matrix_vector_product

END MODULE idrs_mod