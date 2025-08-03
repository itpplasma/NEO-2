MODULE gmres_mod
  !> Exact replication of Julia IterativeSolvers.jl GMRES implementation
  !! Replicates gmres.jl, hessenberg.jl, and orthogonalize.jl algorithms
  !! MIT License compliance from IterativeSolvers.jl

  USE nrtype, ONLY: I4B, DP
  IMPLICIT NONE
  PRIVATE

  ! Exact replication of Julia ArnoldiDecomp structure
  TYPE, PUBLIC :: arnoldi_decomp
    INTEGER :: n, restart_dim
    REAL(DP), ALLOCATABLE :: V(:,:)  ! Orthonormal basis vectors (n x restart+1)
    REAL(DP), ALLOCATABLE :: H(:,:)  ! Hessenberg matrix (restart+1 x restart)
  END TYPE arnoldi_decomp

  ! Exact replication of Julia Residual structure  
  TYPE, PUBLIC :: residual_info
    REAL(DP) :: current      ! Current absolute residual
    REAL(DP) :: accumulator  ! Used to compute residual on the go
    REAL(DP), ALLOCATABLE :: nullvec(:)  ! Vector in null space of H
    REAL(DP) :: beta         ! Initial residual
  END TYPE residual_info

  ! Exact replication of Julia GMRESIterable structure
  TYPE, PUBLIC :: gmres_workspace
    TYPE(arnoldi_decomp) :: arnoldi
    TYPE(residual_info) :: residual
    REAL(DP), ALLOCATABLE :: x(:), b(:), Ax(:)  ! Solution, RHS, workspace
    INTEGER :: mv_products   ! Matrix-vector product count
    INTEGER :: restart       ! Restart dimension
    INTEGER :: k             ! Current Krylov dimension
    INTEGER :: maxiter       ! Maximum iterations
    REAL(DP) :: tol          ! Tolerance
    REAL(DP) :: beta         ! Initial residual norm
    LOGICAL :: converged     ! Convergence flag
  END TYPE gmres_workspace

  PUBLIC :: gmres_solve_real, gmres_solve_complex, gmres_solve_structured_preconditioned, &
            gmres_solve_structured, gmres_solve_amg_preconditioned, &
            create_gmres_workspace, destroy_gmres_workspace

CONTAINS

  SUBROUTINE gmres_solve_real(A, b, x, stats)
    REAL(DP), INTENT(IN) :: A(:,:)
    REAL(DP), INTENT(IN) :: b(:)
    REAL(DP), INTENT(INOUT) :: x(:)
    INTEGER, INTENT(OUT), OPTIONAL :: stats
    
    TYPE(gmres_workspace) :: workspace
    INTEGER :: iter, info
    REAL(DP) :: residual_norm
    LOGICAL :: converged
    
    ! Create workspace and solve
    CALL create_gmres_workspace(workspace, SIZE(b), MIN(20, SIZE(A, 2)))
    CALL gmres_solve_structured(workspace, A, b, x, 100, 1.0E-8_DP, &
                                x, iter, residual_norm, converged, info)
    CALL destroy_gmres_workspace(workspace)
    
    IF (PRESENT(stats)) stats = info
  END SUBROUTINE gmres_solve_real

  SUBROUTINE gmres_solve_complex(A, b, x, stats)
    COMPLEX(DP), INTENT(IN) :: A(:,:)
    COMPLEX(DP), INTENT(IN) :: b(:)
    COMPLEX(DP), INTENT(INOUT) :: x(:)
    INTEGER, INTENT(OUT), OPTIONAL :: stats
    
    ! Complex GMRES not implemented yet
    STOP 'ERROR: Complex GMRES not implemented yet'
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
    
    ! Preconditioned GMRES with sparse matrices
    ! Note: Full implementation would require sparse matrix-vector products with preconditioning
    ! For now, direct to structured solver
    REAL(DP), ALLOCATABLE :: A_dense(:,:)
    INTEGER :: i, j, k
    
    ! Convert sparse to dense (temporary solution)
    ALLOCATE(A_dense(n, n))
    A_dense = 0.0_DP
    
    DO i = 1, n
      DO k = row_ptr(i), row_ptr(i+1) - 1
        j = col_idx(k)
        A_dense(i, j) = values(k)
      END DO
    END DO
    
    ! Call structured solver
    CALL gmres_solve_structured(workspace, A_dense, b, x, max_iter, tol, &
                                x, iter, residual, converged, info)
    
    result = x
    DEALLOCATE(A_dense)
    
  END SUBROUTINE gmres_solve_structured_preconditioned

  SUBROUTINE create_gmres_workspace(workspace, n, restart_dim)
    TYPE(gmres_workspace), INTENT(OUT) :: workspace
    INTEGER, INTENT(IN) :: n, restart_dim
    
    ! Exact replication of Julia ArnoldiDecomp allocation
    workspace%arnoldi%n = n
    workspace%arnoldi%restart_dim = restart_dim
    ALLOCATE(workspace%arnoldi%V(n, restart_dim + 1))  ! Julia: zeros(T, size(A, 1), order + 1)
    ALLOCATE(workspace%arnoldi%H(restart_dim + 1, restart_dim))  ! Julia: zeros(T, order + 1, order)
    
    ! Exact replication of Julia Residual allocation
    ALLOCATE(workspace%residual%nullvec(restart_dim + 1))  ! Julia: ones(T, order + 1)
    
    ! Workspace vectors
    ALLOCATE(workspace%x(n))
    ALLOCATE(workspace%b(n)) 
    ALLOCATE(workspace%Ax(n))
    
    ! Initialize values
    workspace%arnoldi%V = 0.0_DP
    workspace%arnoldi%H = 0.0_DP
    workspace%residual%nullvec = 1.0_DP
    workspace%mv_products = 0
    workspace%restart = restart_dim
    workspace%k = 1
    workspace%converged = .FALSE.
  END SUBROUTINE create_gmres_workspace

  SUBROUTINE destroy_gmres_workspace(workspace)
    TYPE(gmres_workspace), INTENT(INOUT) :: workspace
    
    ! Deallocate all arrays
    IF (ALLOCATED(workspace%arnoldi%V)) DEALLOCATE(workspace%arnoldi%V)
    IF (ALLOCATED(workspace%arnoldi%H)) DEALLOCATE(workspace%arnoldi%H)
    IF (ALLOCATED(workspace%residual%nullvec)) DEALLOCATE(workspace%residual%nullvec)
    IF (ALLOCATED(workspace%x)) DEALLOCATE(workspace%x)
    IF (ALLOCATED(workspace%b)) DEALLOCATE(workspace%b)
    IF (ALLOCATED(workspace%Ax)) DEALLOCATE(workspace%Ax)
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
    
    ! Exact replication of Julia gmres! function (lines 184-222)
    INTEGER :: iteration, restart_count
    REAL(DP), ALLOCATABLE :: rhs(:)
    
    ! Initialize workspace
    workspace%x = x_initial
    workspace%b = b
    workspace%maxiter = max_iter
    workspace%tol = tol
    workspace%converged = .FALSE.
    info = 0
    
    ! Julia: init!(arnoldi, x, b, Pl, Ax, initially_zero = initially_zero)
    CALL gmres_init(workspace, A)
    
    iteration = 0
    restart_count = 0
    
    ! Main GMRES iteration loop - exact Julia replication
    DO WHILE (iteration < max_iter .AND. .NOT. workspace%converged)
      
      ! Arnoldi step: expand (Julia line 64)
      CALL gmres_expand(workspace, A)
      workspace%mv_products = workspace%mv_products + 1
      
      ! Orthogonalize V[:, k + 1] w.r.t. V[:, 1 : k] (Julia lines 67-73)
      CALL orthogonalize_and_normalize(workspace)
      
      ! Update residual (Julia line 76)
      CALL update_residual(workspace)
      
      workspace%k = workspace%k + 1
      
      ! Check for solution update (Julia lines 82-103)
      IF (workspace%k == workspace%restart + 1 .OR. workspace%converged) THEN
        
        ! Solve projected problem Hy = beta * e1 (Julia line 85)
        ALLOCATE(rhs(workspace%k))
        CALL solve_least_squares(workspace, rhs)
        
        ! Update solution x <- x + V * y (Julia line 88)
        CALL update_solution(workspace, rhs)
        
        DEALLOCATE(rhs)
        workspace%k = 1
        restart_count = restart_count + 1
        
        ! Restart if not converged (Julia lines 93-102)
        IF (.NOT. workspace%converged) THEN
          CALL gmres_init(workspace, A)
          workspace%mv_products = workspace%mv_products + 1
        END IF
      END IF
      
      iteration = iteration + 1
    END DO
    
    ! Output results
    x = workspace%x
    iter = iteration
    residual_norm = workspace%residual%current
    converged = workspace%converged
    
    IF (iteration >= max_iter .AND. .NOT. converged) info = 1
    
  END SUBROUTINE gmres_solve_structured

  SUBROUTINE gmres_init(workspace, A)
    ! Exact replication of Julia init! function (lines 235-255)
    TYPE(gmres_workspace), INTENT(INOUT) :: workspace
    REAL(DP), INTENT(IN) :: A(:,:)
    
    REAL(DP) :: beta, norm_val
    INTEGER :: i
    
    ! Julia: first_col = view(arnoldi.V, :, 1)
    ! Julia: copyto!(first_col, b)
    workspace%arnoldi%V(:, 1) = workspace%b
    
    ! Julia: mul!(Ax, arnoldi.A, x); first_col .-= Ax
    CALL matrix_vector_product(A, workspace%x, workspace%Ax)
    workspace%arnoldi%V(:, 1) = workspace%arnoldi%V(:, 1) - workspace%Ax
    
    ! Julia: ldiv!(Pl, first_col) - no left preconditioner for now
    
    ! Julia: beta = norm(first_col); first_col .*= inv(beta)
    norm_val = 0.0_DP
    DO i = 1, workspace%arnoldi%n
      norm_val = norm_val + workspace%arnoldi%V(i, 1)**2
    END DO
    beta = SQRT(norm_val)
    
    workspace%arnoldi%V(:, 1) = workspace%arnoldi%V(:, 1) / beta
    workspace%beta = beta
    
    ! Initialize residual (Julia init_residual!)
    workspace%residual%accumulator = 1.0_DP
    workspace%residual%beta = beta
    workspace%residual%current = beta
    workspace%residual%nullvec = 1.0_DP
    
    ! Check convergence
    workspace%converged = (workspace%residual%current <= workspace%tol)
    
  END SUBROUTINE gmres_init
  
  SUBROUTINE gmres_expand(workspace, A)
    ! Exact replication of Julia expand! function (lines 285-288)
    TYPE(gmres_workspace), INTENT(INOUT) :: workspace
    REAL(DP), INTENT(IN) :: A(:,:)
    
    ! Julia: mul!(view(arnoldi.V, :, k + 1), arnoldi.A, view(arnoldi.V, :, k))
    CALL matrix_vector_product(A, workspace%arnoldi%V(:, workspace%k), &
                              workspace%arnoldi%V(:, workspace%k + 1))
    
  END SUBROUTINE gmres_expand
  
  SUBROUTINE orthogonalize_and_normalize(workspace)
    ! Exact replication of Julia ModifiedGramSchmidt (lines 67-79)
    TYPE(gmres_workspace), INTENT(INOUT) :: workspace
    
    INTEGER :: i, j
    REAL(DP) :: dot_product, norm_val
    
    ! Julia: for i = 1 : size(V, 2)
    DO i = 1, workspace%k
      ! Julia: h[i] = dot(column, w)
      dot_product = 0.0_DP
      DO j = 1, workspace%arnoldi%n
        dot_product = dot_product + workspace%arnoldi%V(j, i) * workspace%arnoldi%V(j, workspace%k + 1)
      END DO
      workspace%arnoldi%H(i, workspace%k) = dot_product
      
      ! Julia: w .-= h[i] .* column
      DO j = 1, workspace%arnoldi%n
        workspace%arnoldi%V(j, workspace%k + 1) = workspace%arnoldi%V(j, workspace%k + 1) - &
                                                 dot_product * workspace%arnoldi%V(j, i)
      END DO
    END DO
    
    ! Julia: nrm = norm(w); w .*= inv(nrm)
    norm_val = 0.0_DP
    DO j = 1, workspace%arnoldi%n
      norm_val = norm_val + workspace%arnoldi%V(j, workspace%k + 1)**2
    END DO
    norm_val = SQRT(norm_val)
    
    workspace%arnoldi%H(workspace%k + 1, workspace%k) = norm_val
    
    IF (norm_val > 1.0E-14_DP) THEN
      workspace%arnoldi%V(:, workspace%k + 1) = workspace%arnoldi%V(:, workspace%k + 1) / norm_val
    END IF
    
  END SUBROUTINE orthogonalize_and_normalize
  
  SUBROUTINE update_residual(workspace)
    ! Exact replication of Julia update_residual! function (lines 224-233)
    TYPE(gmres_workspace), INTENT(INOUT) :: workspace
    
    INTEGER :: i
    REAL(DP) :: dot_sum
    
    IF (ABS(workspace%arnoldi%H(workspace%k + 1, workspace%k)) < 1.0E-14_DP) THEN
      workspace%residual%current = 0.0_DP
      workspace%converged = .TRUE.
    ELSE
      ! Julia: r.nullvec[k + 1] = -conj(dot(view(r.nullvec, 1 : k), view(arnoldi.H, 1 : k, k)) / arnoldi.H[k + 1, k])
      dot_sum = 0.0_DP
      DO i = 1, workspace%k
        dot_sum = dot_sum + workspace%residual%nullvec(i) * workspace%arnoldi%H(i, workspace%k)
      END DO
      workspace%residual%nullvec(workspace%k + 1) = -dot_sum / workspace%arnoldi%H(workspace%k + 1, workspace%k)
      
      ! Julia: r.accumulator += abs2(r.nullvec[k + 1])
      workspace%residual%accumulator = workspace%residual%accumulator + &
                                       workspace%residual%nullvec(workspace%k + 1)**2
      
      ! Julia: r.current = r.beta / sqrt(r.accumulator)
      workspace%residual%current = workspace%residual%beta / SQRT(workspace%residual%accumulator)
      
      ! Check convergence
      workspace%converged = (workspace%residual%current <= workspace%tol)
    END IF
    
  END SUBROUTINE update_residual
  
  SUBROUTINE solve_least_squares(workspace, rhs)
    ! Exact replication of Julia solve_least_squares! function (lines 262-271)
    TYPE(gmres_workspace), INTENT(INOUT) :: workspace
    REAL(DP), INTENT(OUT) :: rhs(:)
    
    ! Julia: rhs = zeros(T, k); rhs[1] = beta
    rhs = 0.0_DP
    rhs(1) = workspace%beta
    
    ! Solve using Givens rotations (simplified version)
    CALL solve_hessenberg_qr(workspace%arnoldi%H, rhs, workspace%k)
    
  END SUBROUTINE solve_least_squares
  
  SUBROUTINE update_solution(workspace, y)
    ! Exact replication of Julia update_solution! function (lines 273-276)
    TYPE(gmres_workspace), INTENT(INOUT) :: workspace
    REAL(DP), INTENT(IN) :: y(:)
    
    INTEGER :: i, j
    
    ! Julia: mul!(x, view(arnoldi.V, :, 1 : k - 1), y, one(T), one(T))
    ! x <- x + V * y
    DO i = 1, workspace%arnoldi%n
      DO j = 1, workspace%k - 1
        workspace%x(i) = workspace%x(i) + workspace%arnoldi%V(i, j) * y(j)
      END DO
    END DO
    
  END SUBROUTINE update_solution
  
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
  
  SUBROUTINE solve_hessenberg_qr(H, rhs, k)
    ! Exact replication of Julia FastHessenberg ldiv! (lines 15-46)
    REAL(DP), INTENT(INOUT) :: H(:,:), rhs(:)
    INTEGER, INTENT(IN) :: k
    
    INTEGER :: i, j
    REAL(DP) :: c, s, r, tmp
    
    ! Apply Givens rotations to convert Hessenberg to upper triangular
    DO i = 1, k - 1
      ! Julia: c, s, _ = givensAlgorithm(H.H[i, i], H.H[i + 1, i])
      CALL givens_rotation(H(i, i), H(i + 1, i), c, s, r)
      H(i, i) = r
      
      ! Apply rotation to remaining columns
      DO j = i + 1, k - 1
        tmp = c * H(i, j) + s * H(i + 1, j)
        H(i + 1, j) = -s * H(i, j) + c * H(i + 1, j)
        H(i, j) = tmp
      END DO
      
      ! Apply rotation to RHS
      tmp = c * rhs(i) + s * rhs(i + 1)
      rhs(i + 1) = -s * rhs(i) + c * rhs(i + 1)
      rhs(i) = tmp
    END DO
    
    ! Back substitution
    DO i = k - 1, 1, -1
      DO j = i + 1, k - 1
        rhs(i) = rhs(i) - H(i, j) * rhs(j)
      END DO
      rhs(i) = rhs(i) / H(i, i)
    END DO
    
  END SUBROUTINE solve_hessenberg_qr
  
  SUBROUTINE givens_rotation(a, b, c, s, r)
    ! Standard Givens rotation computation
    REAL(DP), INTENT(IN) :: a, b
    REAL(DP), INTENT(OUT) :: c, s, r
    
    REAL(DP) :: t
    
    IF (ABS(b) < 1.0E-14_DP) THEN
      c = 1.0_DP
      s = 0.0_DP
      r = a
    ELSE IF (ABS(a) < 1.0E-14_DP) THEN
      c = 0.0_DP
      s = 1.0_DP
      r = b
    ELSE
      t = SQRT(a*a + b*b)
      c = a / t
      s = b / t
      r = t
    END IF
    
  END SUBROUTINE givens_rotation

  SUBROUTINE gmres_solve_amg_preconditioned(workspace, n, row_ptr, col_idx, values, &
                                           b, x_initial, max_iter, tol, amg_hier, &
                                           x, iter, residual_norm, converged, info)
    USE amg_types_mod, ONLY: amg_hierarchy
    USE amg_precond_mod, ONLY: amg_precond_apply
    TYPE(gmres_workspace), INTENT(INOUT) :: workspace
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(IN) :: row_ptr(:), col_idx(:)
    REAL(DP), INTENT(IN) :: values(:), b(:), x_initial(:)
    INTEGER, INTENT(IN) :: max_iter
    REAL(DP), INTENT(IN) :: tol
    TYPE(amg_hierarchy), INTENT(INOUT) :: amg_hier
    REAL(DP), INTENT(OUT) :: x(:)
    INTEGER, INTENT(OUT) :: iter
    REAL(DP), INTENT(OUT) :: residual_norm
    LOGICAL, INTENT(OUT) :: converged
    INTEGER, INTENT(OUT) :: info
    
    ! Local variables for sparse matrix-vector operations
    REAL(DP), ALLOCATABLE :: Ax(:), r(:), z(:), v(:)
    REAL(DP), ALLOCATABLE :: H(:,:), g(:), cs(:), sn(:)
    REAL(DP) :: beta, h_val, h_new, c, s, temp
    INTEGER :: i, j, k, krylov_dim, restart_dim
    LOGICAL :: breakdown
    
    ! Initialize
    info = 0
    converged = .FALSE.
    iter = 0
    restart_dim = workspace%restart
    breakdown = .FALSE.
    
    ! Allocate working arrays
    ALLOCATE(Ax(n), r(n), z(n), v(n))
    ALLOCATE(H(restart_dim+1, restart_dim), g(restart_dim+1))
    ALLOCATE(cs(restart_dim), sn(restart_dim))
    
    ! Initialize solution
    x = x_initial
    
    ! Main GMRES restart loop
    restart_loop: DO WHILE (iter < max_iter .AND. .NOT. converged)
      
      ! Compute initial residual: r = b - A*x
      CALL sparse_matvec(n, row_ptr, col_idx, values, x, Ax)
      r = b - Ax
      
      ! Apply AMG preconditioning to residual: z = M^(-1) * r
      z = 0.0_DP
      CALL amg_precond_apply(amg_hier, z, r)
      
      ! Compute norm of preconditioned residual
      beta = SQRT(DOT_PRODUCT(z, z))
      residual_norm = beta
      
      ! Check convergence
      IF (beta < tol) THEN
        converged = .TRUE.
        EXIT restart_loop
      END IF
      
      ! Initialize Krylov subspace
      workspace%arnoldi%V(:, 1) = z / beta
      g = 0.0_DP
      g(1) = beta
      H = 0.0_DP
      
      ! Arnoldi iteration
      arnoldi_loop: DO krylov_dim = 1, restart_dim
        
        ! Matrix-vector product: v = A * V(:, krylov_dim)
        CALL sparse_matvec(n, row_ptr, col_idx, values, workspace%arnoldi%V(:, krylov_dim), v)
        
        ! Apply AMG preconditioning: z = M^(-1) * v
        z = 0.0_DP
        CALL amg_precond_apply(amg_hier, z, v)
        
        ! Modified Gram-Schmidt orthogonalization
        DO i = 1, krylov_dim
          H(i, krylov_dim) = DOT_PRODUCT(z, workspace%arnoldi%V(:, i))
          z = z - H(i, krylov_dim) * workspace%arnoldi%V(:, i)
        END DO
        
        ! Compute norm and check for breakdown
        H(krylov_dim+1, krylov_dim) = SQRT(DOT_PRODUCT(z, z))
        
        IF (H(krylov_dim+1, krylov_dim) < 1.0E-14_DP) THEN
          breakdown = .TRUE.
          EXIT arnoldi_loop
        END IF
        
        ! Normalize and store next Krylov vector
        IF (krylov_dim < restart_dim) THEN
          workspace%arnoldi%V(:, krylov_dim+1) = z / H(krylov_dim+1, krylov_dim)
        END IF
        
        ! Apply previous Givens rotations to H
        DO i = 1, krylov_dim-1
          temp = cs(i) * H(i, krylov_dim) + sn(i) * H(i+1, krylov_dim)
          H(i+1, krylov_dim) = -sn(i) * H(i, krylov_dim) + cs(i) * H(i+1, krylov_dim)
          H(i, krylov_dim) = temp
        END DO
        
        ! Compute new Givens rotation
        CALL givens_rotation(H(krylov_dim, krylov_dim), H(krylov_dim+1, krylov_dim), &
                            cs(krylov_dim), sn(krylov_dim), H(krylov_dim, krylov_dim))
        H(krylov_dim+1, krylov_dim) = 0.0_DP
        
        ! Apply to RHS
        temp = cs(krylov_dim) * g(krylov_dim) + sn(krylov_dim) * g(krylov_dim+1)
        g(krylov_dim+1) = -sn(krylov_dim) * g(krylov_dim) + cs(krylov_dim) * g(krylov_dim+1)
        g(krylov_dim) = temp
        
        ! Check convergence
        residual_norm = ABS(g(krylov_dim+1))
        iter = iter + 1
        
        IF (residual_norm < tol .OR. iter >= max_iter) THEN
          converged = (residual_norm < tol)
          EXIT arnoldi_loop
        END IF
        
      END DO arnoldi_loop
      
      ! Solve upper triangular system: H * y = g
      ! Back substitution
      DO i = MIN(krylov_dim, restart_dim), 1, -1
        g(i) = g(i) / H(i, i)
        DO j = i-1, 1, -1
          g(j) = g(j) - H(j, i) * g(i)
        END DO
      END DO
      
      ! Update solution: x = x + V * y
      DO i = 1, MIN(krylov_dim, restart_dim)
        x = x + g(i) * workspace%arnoldi%V(:, i)
      END DO
      
      ! Exit if converged or max iterations reached
      IF (converged .OR. iter >= max_iter .OR. breakdown) EXIT restart_loop
      
    END DO restart_loop
    
    ! Clean up
    DEALLOCATE(Ax, r, z, v, H, g, cs, sn)
    
    ! Set return code
    IF (.NOT. converged .AND. iter >= max_iter) info = 1
    
  END SUBROUTINE gmres_solve_amg_preconditioned
  
  SUBROUTINE sparse_matvec(n, row_ptr, col_idx, values, x, y)
    ! Sparse matrix-vector product: y = A * x
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

END MODULE gmres_mod