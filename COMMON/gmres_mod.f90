MODULE gmres_mod
  !> GMRES iterative solver with ILU preconditioning
  !! 
  !! This module implements the Generalized Minimal Residual (GMRES) method
  !! for solving sparse linear systems Ax = b. GMRES is particularly effective
  !! for non-symmetric and ill-conditioned matrices where BiCGSTAB may fail.
  !!
  !! Key features:
  !! - Restarted GMRES(m) with configurable restart parameter
  !! - Support for ILU(k) preconditioning
  !! - Arnoldi orthogonalization process
  !! - QR decomposition via Givens rotations
  !! - Comprehensive convergence monitoring
  !!
  !! @author Generated with Claude Code
  !! @date 2025-08-02

  USE nrtype, ONLY: I4B, DP
  USE ilu_precond_mod
  IMPLICIT NONE
  PRIVATE

  ! Arnoldi decomposition type (based on IterativeSolvers.jl ArnoldiDecomp)
  TYPE, PUBLIC :: arnoldi_decomp
    REAL(DP), ALLOCATABLE :: V(:,:)     ! Orthonormal basis vectors (n x restart+1)
    REAL(DP), ALLOCATABLE :: H(:,:)     ! Upper Hessenberg matrix (restart+1 x restart)  
    INTEGER(I4B) :: order               ! Restart dimension
  END TYPE arnoldi_decomp
  
  ! GMRES workspace type (based on IterativeSolvers.jl GMRESIterable)
  TYPE, PUBLIC :: gmres_workspace
    TYPE(arnoldi_decomp) :: arnoldi                    ! Embedded Arnoldi decomposition
    REAL(DP), ALLOCATABLE :: givens_c(:), givens_s(:) ! Givens rotation coefficients
    REAL(DP), ALLOCATABLE :: rhs_qr(:)                ! QR decomposition RHS vector
    REAL(DP) :: residual_norm                          ! Current residual norm
    INTEGER(I4B) :: k                                  ! Current Krylov subspace dimension
  END TYPE gmres_workspace
  
  ! GMRES statistics type
  TYPE, PUBLIC :: gmres_stats
    INTEGER :: iterations = 0              ! Total iterations performed
    INTEGER :: restarts = 0               ! Number of restarts
    REAL(kind=dp) :: initial_residual = 0.0_dp    ! ||r0||
    REAL(kind=dp) :: final_residual = 0.0_dp      ! ||r_final||
    REAL(kind=dp) :: solve_time = 0.0_dp           ! Wall clock time
    LOGICAL :: converged = .FALSE.                 ! Convergence flag
    CHARACTER(len=64) :: message = ''              ! Status message
  END TYPE gmres_stats

  PUBLIC :: gmres_solve_real, gmres_solve_complex
  PUBLIC :: create_arnoldi_decomp, destroy_arnoldi_decomp
  PUBLIC :: create_gmres_workspace, destroy_gmres_workspace
  PUBLIC :: modified_gram_schmidt, arnoldi_expand, arnoldi_step
  PUBLIC :: check_orthogonality, arnoldi_expand_csr
  PUBLIC :: compute_givens_rotation, apply_givens_rotation
  PUBLIC :: qr_decompose_hessenberg, solve_least_squares_qr, update_qr_factorization
  PUBLIC :: gmres_solve_structured, gmres_iterate
  PUBLIC :: initialize_gmres, finalize_gmres_iteration
  PUBLIC :: gmres_solve_structured_preconditioned

CONTAINS

  !> Solve real sparse linear system using GMRES with optional preconditioning
  !!
  !! @param[in]    n         Matrix dimension
  !! @param[in]    row_ptr   CSR row pointers (size: n+1)
  !! @param[in]    col_idx   CSR column indices (size: nnz)
  !! @param[in]    val       CSR matrix values (size: nnz)
  !! @param[in]    b         Right-hand side vector (size: n)
  !! @param[inout] x         Solution vector (size: n) - initial guess on input
  !! @param[in]    restart   GMRES restart parameter (m)
  !! @param[in]    max_iter  Maximum total iterations
  !! @param[in]    abs_tol   Absolute tolerance
  !! @param[in]    rel_tol   Relative tolerance
  !! @param[in]    ilu_fac   ILU factorization (optional)
  !! @param[in]    use_ilu   Use ILU preconditioning flag
  !! @param[in]    verbose   Print convergence information
  !! @param[out]   stats     GMRES statistics
  !! @param[out]   info      Return code (0=success, >0=failure)
  SUBROUTINE gmres_solve_real(n, row_ptr, col_idx, val, b, x, restart, max_iter, &
                              abs_tol, rel_tol, ilu_fac, use_ilu, verbose, stats, info)
    IMPLICIT NONE
    INTEGER(I4B), INTENT(in) :: n
    INTEGER(I4B), DIMENSION(:), INTENT(in) :: row_ptr, col_idx
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: val, b
    REAL(kind=dp), DIMENSION(:), INTENT(inout) :: x
    INTEGER(I4B), INTENT(in) :: restart, max_iter
    REAL(kind=dp), INTENT(in) :: abs_tol, rel_tol
    TYPE(ilu_factorization), INTENT(in), OPTIONAL :: ilu_fac
    LOGICAL, INTENT(in) :: use_ilu, verbose
    TYPE(gmres_stats), INTENT(out) :: stats
    INTEGER(I4B), INTENT(out) :: info

    ! Local variables for GMRES algorithm
    REAL(kind=dp), ALLOCATABLE :: V(:,:)           ! Krylov basis vectors (n x restart+1)
    REAL(kind=dp), ALLOCATABLE :: H(:,:)           ! Upper Hessenberg matrix
    REAL(kind=dp), ALLOCATABLE :: g(:)             ! RHS of least squares problem
    REAL(kind=dp), ALLOCATABLE :: y(:)             ! Solution of least squares problem
    REAL(kind=dp), ALLOCATABLE :: c(:), s(:)       ! Givens rotation coefficients
    REAL(kind=dp), ALLOCATABLE :: r(:), w(:), z(:) ! Work vectors
    
    ! Local variables
    INTEGER(I4B) :: iter, restart_iter, i, j, k
    REAL(kind=dp) :: beta, norm_b, tol_eff, residual_norm
    REAL(kind=dp) :: temp, h_ik, h_kk
    REAL(kind=dp) :: start_time, end_time
    LOGICAL :: converged
    
    ! Initialize
    CALL CPU_TIME(start_time)
    info = 0
    converged = .FALSE.
    stats%iterations = 0
    stats%restarts = 0
    
    ! Allocate workspace
    ALLOCATE(V(n, restart+1), H(restart+1, restart), g(restart+1), y(restart))
    ALLOCATE(c(restart), s(restart), r(n), w(n), z(n))
    
    ! Compute ||b|| for relative tolerance
    norm_b = SQRT(DOT_PRODUCT(b, b))
    stats%initial_residual = norm_b
    
    IF (norm_b < TINY(1.0_dp)) THEN
      ! Zero RHS - solution is zero
      x = 0.0_dp
      stats%converged = .TRUE.
      stats%message = 'Zero RHS - trivial solution'
      CALL CPU_TIME(end_time)
      stats%solve_time = end_time - start_time
      DEALLOCATE(V, H, g, y, c, s, r, w, z)
      RETURN
    END IF
    
    tol_eff = MAX(abs_tol, rel_tol * norm_b)
    
    IF (verbose) THEN
      PRINT *, 'GMRES: Starting solve'
      PRINT *, '  Matrix size:', n, 'x', n
      PRINT *, '  Restart parameter:', restart
      PRINT *, '  ||b|| =', norm_b
      PRINT *, '  Effective tolerance:', tol_eff
    END IF

    ! Main GMRES restart loop
    restart_loop: DO 
      stats%restarts = stats%restarts + 1
      
      ! Compute initial residual r = b - A*x
      CALL matvec_csr_real(n, row_ptr, col_idx, val, x, w)
      r = b - w
      
      ! Apply preconditioning: solve M*z = r
      IF (use_ilu .AND. PRESENT(ilu_fac)) THEN
        CALL ilu_solve(ilu_fac, r, z)
      ELSE
        z = r
      END IF
      
      beta = SQRT(DOT_PRODUCT(z, z))
      residual_norm = beta
      
      IF (verbose .AND. stats%restarts == 1) THEN
        PRINT *, '  Initial residual norm:', residual_norm
      END IF
      
      ! Check convergence
      IF (residual_norm <= tol_eff) THEN
        converged = .TRUE.
        EXIT restart_loop
      END IF
      
      ! Initialize first Krylov vector
      V(:, 1) = z / beta
      g = 0.0_dp
      g(1) = beta
      
      ! Arnoldi process
      arnoldi_loop: DO k = 1, restart
        iter = k
        stats%iterations = stats%iterations + 1
        
        ! Apply matrix: w = A * V(:,k)
        CALL matvec_csr_real(n, row_ptr, col_idx, val, V(:, k), w)
        
        ! Apply preconditioning: solve M*z = w  
        IF (use_ilu .AND. PRESENT(ilu_fac)) THEN
          CALL ilu_solve(ilu_fac, w, z)
        ELSE
          z = w
        END IF
        
        ! Modified Gram-Schmidt orthogonalization
        DO j = 1, k
          H(j, k) = DOT_PRODUCT(z, V(:, j))
          z = z - H(j, k) * V(:, j)
        END DO
        
        H(k+1, k) = SQRT(DOT_PRODUCT(z, z))
        
        IF (H(k+1, k) > TINY(1.0_dp)) THEN
          V(:, k+1) = z / H(k+1, k)
        END IF
        
        ! Apply previous Givens rotations to new column of H
        DO j = 1, k-1
          temp = c(j) * H(j, k) + s(j) * H(j+1, k)
          H(j+1, k) = -s(j) * H(j, k) + c(j) * H(j+1, k)
          H(j, k) = temp
        END DO
        
        ! Compute new Givens rotation
        IF (ABS(H(k+1, k)) < TINY(1.0_dp)) THEN
          c(k) = 1.0_dp
          s(k) = 0.0_dp
        ELSE
          temp = SQRT(H(k, k)**2 + H(k+1, k)**2)
          c(k) = H(k, k) / temp
          s(k) = H(k+1, k) / temp
        END IF
        
        ! Apply new Givens rotation
        H(k, k) = c(k) * H(k, k) + s(k) * H(k+1, k)
        H(k+1, k) = 0.0_dp
        
        ! Update RHS of least squares problem
        temp = c(k) * g(k) + s(k) * g(k+1)
        g(k+1) = -s(k) * g(k) + c(k) * g(k+1)
        g(k) = temp
        
        residual_norm = ABS(g(k+1))
        
        IF (verbose) THEN
          PRINT *, '    Restart', stats%restarts, 'iter', k, ': residual =', residual_norm
        END IF
        
        ! Check convergence
        IF (residual_norm <= tol_eff) THEN
          converged = .TRUE.
          EXIT arnoldi_loop
        END IF
        
        ! Check maximum iterations
        IF (stats%iterations >= max_iter) THEN
          EXIT arnoldi_loop
        END IF
      END DO arnoldi_loop
      
      ! Solve upper triangular system H*y = g by back substitution
      DO i = iter, 1, -1
        y(i) = g(i)
        DO j = i+1, iter
          y(i) = y(i) - H(i, j) * y(j)
        END DO
        y(i) = y(i) / H(i, i)
      END DO
      
      ! Update solution: x = x + V*y
      DO i = 1, iter
        x = x + y(i) * V(:, i)
      END DO
      
      IF (converged .OR. stats%iterations >= max_iter) THEN
        EXIT restart_loop
      END IF
      
    END DO restart_loop
    
    stats%final_residual = residual_norm
    stats%converged = converged
    
    IF (converged) THEN
      stats%message = 'Converged successfully'
      IF (verbose) THEN
        PRINT *, 'GMRES: Converged in', stats%iterations, 'iterations'
        PRINT *, '  Final residual:', residual_norm
      END IF
    ELSE
      stats%message = 'Maximum iterations reached'
      info = 1
      IF (verbose) THEN
        PRINT *, 'GMRES: Did not converge in', stats%iterations, 'iterations'
        PRINT *, '  Final residual:', residual_norm
      END IF
    END IF
    
    CALL CPU_TIME(end_time)
    stats%solve_time = end_time - start_time
    
    DEALLOCATE(V, H, g, y, c, s, r, w, z)
    
  END SUBROUTINE gmres_solve_real

  !> Complex version - placeholder for now
  SUBROUTINE gmres_solve_complex(n, row_ptr, col_idx, val, b, x, restart, max_iter, &
                                 abs_tol, rel_tol, ilu_fac, use_ilu, verbose, stats, info)
    INTEGER(I4B), INTENT(in) :: n
    INTEGER(I4B), DIMENSION(:), INTENT(in) :: row_ptr, col_idx
    COMPLEX(kind=dp), DIMENSION(:), INTENT(in) :: val, b
    COMPLEX(kind=dp), DIMENSION(:), INTENT(inout) :: x
    INTEGER(I4B), INTENT(in) :: restart, max_iter
    REAL(kind=dp), INTENT(in) :: abs_tol, rel_tol
    TYPE(ilu_factorization_complex), INTENT(in), OPTIONAL :: ilu_fac
    LOGICAL, INTENT(in) :: use_ilu, verbose
    TYPE(gmres_stats), INTENT(out) :: stats
    INTEGER(I4B), INTENT(out) :: info
    
    ! Placeholder - implement complex version later
    PRINT *, 'ERROR: Complex GMRES not yet implemented'
    info = -1
    stats%converged = .FALSE.
    stats%message = 'Complex GMRES not implemented'
    
  END SUBROUTINE gmres_solve_complex

  !> Matrix-vector multiplication for CSR format
  SUBROUTINE matvec_csr_real(n, row_ptr, col_idx, val, x, y)
    INTEGER(I4B), INTENT(in) :: n
    INTEGER(I4B), DIMENSION(:), INTENT(in) :: row_ptr, col_idx
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: val, x
    REAL(kind=dp), DIMENSION(:), INTENT(out) :: y
    
    INTEGER(I4B) :: i, j
    
    DO i = 1, n
      y(i) = 0.0_dp
      DO j = row_ptr(i), row_ptr(i+1) - 1
        y(i) = y(i) + val(j) * x(col_idx(j))
      END DO
    END DO
    
  END SUBROUTINE matvec_csr_real

  ! ============================================================================
  ! Structured data management routines (IterativeSolvers.jl template approach)
  ! ============================================================================

  SUBROUTINE create_arnoldi_decomp(arnoldi, n, restart)
    ! Create and initialize Arnoldi decomposition workspace
    TYPE(arnoldi_decomp), INTENT(out) :: arnoldi
    INTEGER(I4B), INTENT(in) :: n        ! Matrix size
    INTEGER(I4B), INTENT(in) :: restart  ! Restart dimension
    
    ! Store restart dimension
    arnoldi%order = restart
    
    ! Allocate V matrix: n x (restart+1)
    ! Extra column for the (restart+1)-th Arnoldi vector
    ALLOCATE(arnoldi%V(n, restart + 1))
    arnoldi%V = 0.0_DP
    
    ! Allocate H matrix: (restart+1) x restart
    ! Upper Hessenberg matrix from Arnoldi process
    ALLOCATE(arnoldi%H(restart + 1, restart))
    arnoldi%H = 0.0_DP
    
  END SUBROUTINE create_arnoldi_decomp
  
  SUBROUTINE destroy_arnoldi_decomp(arnoldi)
    ! Clean up Arnoldi decomposition workspace
    TYPE(arnoldi_decomp), INTENT(inout) :: arnoldi
    
    IF (ALLOCATED(arnoldi%V)) DEALLOCATE(arnoldi%V)
    IF (ALLOCATED(arnoldi%H)) DEALLOCATE(arnoldi%H)
    arnoldi%order = 0
    
  END SUBROUTINE destroy_arnoldi_decomp
  
  SUBROUTINE create_gmres_workspace(workspace, n, restart)
    ! Create and initialize GMRES workspace
    TYPE(gmres_workspace), INTENT(out) :: workspace
    INTEGER(I4B), INTENT(in) :: n        ! Matrix size  
    INTEGER(I4B), INTENT(in) :: restart  ! Restart dimension
    
    ! Initialize embedded Arnoldi decomposition
    CALL create_arnoldi_decomp(workspace%arnoldi, n, restart)
    
    ! Allocate Givens rotation coefficients (restart elements)
    ALLOCATE(workspace%givens_c(restart))
    ALLOCATE(workspace%givens_s(restart))
    workspace%givens_c = 0.0_DP
    workspace%givens_s = 0.0_DP
    
    ! Allocate QR RHS vector (restart+1 elements)
    ALLOCATE(workspace%rhs_qr(restart + 1))
    workspace%rhs_qr = 0.0_DP
    
    ! Initialize scalar values
    workspace%residual_norm = 0.0_DP
    workspace%k = 0
    
  END SUBROUTINE create_gmres_workspace
  
  SUBROUTINE destroy_gmres_workspace(workspace)
    ! Clean up GMRES workspace
    TYPE(gmres_workspace), INTENT(inout) :: workspace
    
    ! Clean up embedded Arnoldi decomposition
    CALL destroy_arnoldi_decomp(workspace%arnoldi)
    
    ! Deallocate Givens rotation arrays
    IF (ALLOCATED(workspace%givens_c)) DEALLOCATE(workspace%givens_c)
    IF (ALLOCATED(workspace%givens_s)) DEALLOCATE(workspace%givens_s)
    
    ! Deallocate QR RHS vector
    IF (ALLOCATED(workspace%rhs_qr)) DEALLOCATE(workspace%rhs_qr)
    
    ! Reset scalar values
    workspace%residual_norm = 0.0_DP
    workspace%k = 0
    
  END SUBROUTINE destroy_gmres_workspace

  ! ============================================================================
  ! Arnoldi orthogonalization routines (IterativeSolvers.jl template approach)
  ! ============================================================================

  SUBROUTINE modified_gram_schmidt(Q, v, h_coeffs)
    ! Modified Gram-Schmidt orthogonalization of vector v against columns of Q
    ! Based on IterativeSolvers.jl orthogonalize_and_normalize! function
    REAL(DP), INTENT(in) :: Q(:,:)              ! Existing orthonormal basis
    REAL(DP), INTENT(inout) :: v(:)             ! Vector to orthogonalize (modified in place)
    REAL(DP), INTENT(out) :: h_coeffs(:)        ! Orthogonalization coefficients
    
    INTEGER(I4B) :: k, j
    
    k = SIZE(Q, 2)  ! Number of existing basis vectors
    
    ! Modified Gram-Schmidt: orthogonalize against each basis vector
    DO j = 1, k
      h_coeffs(j) = DOT_PRODUCT(v, Q(:, j))     ! Projection coefficient
      v = v - h_coeffs(j) * Q(:, j)             ! Remove component
    END DO
    
  END SUBROUTINE modified_gram_schmidt
  
  SUBROUTINE arnoldi_expand(A, v_in, v_out)
    ! Matrix-vector multiplication: v_out = A * v_in
    ! Dense matrix version for testing
    REAL(DP), INTENT(in) :: A(:,:)              ! Dense matrix
    REAL(DP), INTENT(in) :: v_in(:)             ! Input vector
    REAL(DP), INTENT(out) :: v_out(:)           ! Output vector
    
    v_out = MATMUL(A, v_in)
    
  END SUBROUTINE arnoldi_expand
  
  SUBROUTINE arnoldi_expand_csr(n, row_ptr, col_idx, val, v_in, v_out)
    ! Matrix-vector multiplication for CSR sparse format
    ! v_out = A * v_in where A is in CSR format
    INTEGER(I4B), INTENT(in) :: n                  ! Matrix dimension
    INTEGER(I4B), INTENT(in) :: row_ptr(:)         ! CSR row pointers
    INTEGER(I4B), INTENT(in) :: col_idx(:)         ! CSR column indices
    REAL(DP), INTENT(in) :: val(:)                 ! CSR values
    REAL(DP), INTENT(in) :: v_in(:)                ! Input vector
    REAL(DP), INTENT(out) :: v_out(:)              ! Output vector
    
    INTEGER(I4B) :: i, j
    
    v_out = 0.0_DP
    DO i = 1, n
      DO j = row_ptr(i), row_ptr(i+1) - 1
        v_out(i) = v_out(i) + val(j) * v_in(col_idx(j))
      END DO
    END DO
    
  END SUBROUTINE arnoldi_expand_csr
  
  SUBROUTINE arnoldi_step(A, arnoldi, k)
    ! Complete Arnoldi step: expand, orthogonalize, normalize
    ! Based on IterativeSolvers.jl expand! and orthogonalize_and_normalize!
    REAL(DP), INTENT(in) :: A(:,:)                 ! Matrix (dense for testing)
    TYPE(arnoldi_decomp), INTENT(inout) :: arnoldi ! Arnoldi workspace
    INTEGER(I4B), INTENT(in) :: k                  ! Current step
    
    REAL(DP), ALLOCATABLE :: w(:)                  ! Work vector
    REAL(DP) :: norm_w
    INTEGER(I4B) :: n
    
    n = SIZE(arnoldi%V, 1)
    ALLOCATE(w(n))
    
    ! Arnoldi expansion: w = A * V(:, k)
    CALL arnoldi_expand(A, arnoldi%V(:, k), w)
    
    ! Modified Gram-Schmidt orthogonalization against V(:, 1:k)
    CALL modified_gram_schmidt(arnoldi%V(:, 1:k), w, arnoldi%H(1:k, k))
    
    ! Compute norm and store in H(k+1, k)
    norm_w = SQRT(DOT_PRODUCT(w, w))
    arnoldi%H(k+1, k) = norm_w
    
    ! Normalize and store as next Arnoldi vector
    IF (norm_w > TINY(1.0_DP)) THEN
      arnoldi%V(:, k+1) = w / norm_w
    ELSE
      ! Breakdown - use random vector (shouldn't happen in tests)
      arnoldi%V(:, k+1) = 0.0_DP
      IF (k+1 <= SIZE(arnoldi%V, 1)) arnoldi%V(k+1, k+1) = 1.0_DP
    END IF
    
    DEALLOCATE(w)
    
  END SUBROUTINE arnoldi_step
  
  SUBROUTINE check_orthogonality(Q, max_error)
    ! Check orthogonality of matrix Q: measure || Q^T * Q - I ||
    REAL(DP), INTENT(in) :: Q(:,:)                 ! Matrix to check
    REAL(DP), INTENT(out) :: max_error             ! Maximum orthogonality error
    
    REAL(DP), ALLOCATABLE :: QtQ(:,:)              ! Q^T * Q
    INTEGER(I4B) :: k, i, j
    
    k = SIZE(Q, 2)
    ALLOCATE(QtQ(k, k))
    
    ! Compute Q^T * Q
    DO i = 1, k
      DO j = 1, k
        QtQ(i, j) = DOT_PRODUCT(Q(:, i), Q(:, j))
      END DO
    END DO
    
    ! Subtract identity matrix and find maximum error
    max_error = 0.0_DP
    DO i = 1, k
      DO j = 1, k
        IF (i == j) THEN
          max_error = MAX(max_error, ABS(QtQ(i, j) - 1.0_DP))
        ELSE
          max_error = MAX(max_error, ABS(QtQ(i, j)))
        END IF
      END DO
    END DO
    
    DEALLOCATE(QtQ)
    
  END SUBROUTINE check_orthogonality

  ! ============================================================================
  ! QR decomposition routines via Givens rotations (IterativeSolvers.jl approach)
  ! ============================================================================

  SUBROUTINE compute_givens_rotation(a, b, c, s)
    ! Compute Givens rotation to eliminate b: [c s; -s c] * [a; b] = [r; 0]
    ! Based on IterativeSolvers.jl Givens rotation generation
    REAL(DP), INTENT(in) :: a, b               ! Input values
    REAL(DP), INTENT(out) :: c, s              ! Givens rotation coefficients
    
    REAL(DP) :: r, t
    
    IF (ABS(b) < TINY(1.0_DP)) THEN
      ! b is effectively zero
      c = 1.0_DP
      s = 0.0_DP
    ELSE IF (ABS(a) < TINY(1.0_DP)) THEN
      ! a is effectively zero
      c = 0.0_DP
      s = 1.0_DP
    ELSE
      ! General case - use stable Givens rotation computation
      IF (ABS(b) > ABS(a)) THEN
        t = a / b
        s = 1.0_DP / SQRT(1.0_DP + t*t)
        c = s * t
      ELSE
        t = b / a
        c = 1.0_DP / SQRT(1.0_DP + t*t)
        s = c * t
      END IF
      
      ! Ensure c >= 0 for numerical stability
      IF (c < 0.0_DP) THEN
        c = -c
        s = -s
      END IF
    END IF
    
  END SUBROUTINE compute_givens_rotation
  
  SUBROUTINE apply_givens_rotation(H, k, c, s)
    ! Apply Givens rotation to rows k and k+1 of matrix H
    ! Eliminates H(k+1, k) by rotating rows k and k+1
    REAL(DP), INTENT(inout) :: H(:,:)          ! Matrix to transform
    INTEGER(I4B), INTENT(in) :: k              ! Row index
    REAL(DP), INTENT(in) :: c, s               ! Givens rotation coefficients
    
    REAL(DP) :: temp
    INTEGER(I4B) :: j, n_cols
    
    n_cols = SIZE(H, 2)
    
    ! Apply rotation to all columns: [c s; -s c] * [H(k,:); H(k+1,:)]
    DO j = 1, n_cols
      temp = c * H(k, j) + s * H(k+1, j)
      H(k+1, j) = -s * H(k, j) + c * H(k+1, j)
      H(k, j) = temp
    END DO
    
  END SUBROUTINE apply_givens_rotation
  
  SUBROUTINE qr_decompose_hessenberg(H, g, c, s)
    ! QR decomposition of upper Hessenberg matrix via Givens rotations
    ! Also applies same rotations to RHS vector g
    REAL(DP), INTENT(inout) :: H(:,:)          ! Hessenberg matrix (modified to R)
    REAL(DP), INTENT(inout) :: g(:)            ! RHS vector (transformed)
    REAL(DP), INTENT(out) :: c(:), s(:)        ! Givens rotation coefficients
    
    INTEGER(I4B) :: k, m
    REAL(DP) :: temp
    
    m = SIZE(H, 2)  ! Number of columns
    
    ! Apply Givens rotations to eliminate subdiagonal
    DO k = 1, m
      IF (k+1 <= SIZE(H, 1)) THEN
        ! Compute rotation to eliminate H(k+1, k)
        CALL compute_givens_rotation(H(k, k), H(k+1, k), c(k), s(k))
        
        ! Apply rotation to matrix H
        CALL apply_givens_rotation(H, k, c(k), s(k))
        
        ! Apply same rotation to RHS vector g
        temp = c(k) * g(k) + s(k) * g(k+1)
        g(k+1) = -s(k) * g(k) + c(k) * g(k+1)
        g(k) = temp
      END IF
    END DO
    
  END SUBROUTINE qr_decompose_hessenberg
  
  SUBROUTINE solve_least_squares_qr(R, g, y)
    ! Solve upper triangular system R * y = g via back substitution
    ! Used after QR decomposition to solve GMRES least-squares subproblem
    REAL(DP), INTENT(in) :: R(:,:)             ! Upper triangular matrix R
    REAL(DP), INTENT(in) :: g(:)               ! Transformed RHS vector
    REAL(DP), INTENT(out) :: y(:)              ! Solution vector
    
    INTEGER(I4B) :: i, j, n
    
    n = SIZE(R, 1)
    
    ! Back substitution: solve R * y = g
    DO i = n, 1, -1
      y(i) = g(i)
      
      ! Subtract off-diagonal contributions
      DO j = i+1, n
        y(i) = y(i) - R(i, j) * y(j)
      END DO
      
      ! Divide by diagonal element
      IF (ABS(R(i, i)) > TINY(1.0_DP)) THEN
        y(i) = y(i) / R(i, i)
      ELSE
        ! Near-singular case - could indicate breakdown
        y(i) = 0.0_DP
      END IF
    END DO
    
  END SUBROUTINE solve_least_squares_qr
  
  SUBROUTINE update_qr_factorization(g, k, c, s)
    ! Update QR factorization incrementally - apply rotation to RHS vector
    ! Used in incremental GMRES QR construction
    REAL(DP), INTENT(inout) :: g(:)            ! RHS vector to update
    INTEGER(I4B), INTENT(in) :: k              ! Position of rotation
    REAL(DP), INTENT(in) :: c, s               ! Givens rotation coefficients
    
    REAL(DP) :: temp
    
    ! Apply Givens rotation to positions k and k+1 of vector g
    temp = c * g(k) + s * g(k+1)
    g(k+1) = -s * g(k) + c * g(k+1)
    g(k) = temp
    
  END SUBROUTINE update_qr_factorization

  ! ============================================================================
  ! Complete GMRES algorithm routines (IterativeSolvers.jl template approach)
  ! ============================================================================

  SUBROUTINE gmres_solve_structured(workspace, A, b, x_initial, max_iter, tol, &
                                    x, iterations, residual_norm, converged, info)
    ! High-level structured GMRES solver interface
    ! Based on IterativeSolvers.jl gmres! function
    TYPE(gmres_workspace), INTENT(inout) :: workspace
    REAL(DP), INTENT(in) :: A(:,:)             ! Matrix (dense for testing)
    REAL(DP), INTENT(in) :: b(:)               ! RHS vector
    REAL(DP), INTENT(in) :: x_initial(:)       ! Initial guess
    INTEGER(I4B), INTENT(in) :: max_iter       ! Maximum iterations
    REAL(DP), INTENT(in) :: tol                ! Convergence tolerance
    REAL(DP), INTENT(out) :: x(:)              ! Solution vector
    INTEGER(I4B), INTENT(out) :: iterations    ! Total iterations performed
    REAL(DP), INTENT(out) :: residual_norm     ! Final residual norm
    LOGICAL, INTENT(out) :: converged          ! Convergence flag
    INTEGER(I4B), INTENT(out) :: info          ! Return code
    
    REAL(DP) :: tol_eff, norm_b
    INTEGER(I4B) :: restart_count
    
    ! Initialize
    info = 0
    converged = .FALSE.
    iterations = 0
    restart_count = 0
    x = x_initial
    
    ! Compute effective tolerance
    norm_b = SQRT(DOT_PRODUCT(b, b))
    tol_eff = tol * norm_b
    
    ! Check for zero RHS
    IF (norm_b < TINY(1.0_DP)) THEN
      x = 0.0_DP
      converged = .TRUE.
      residual_norm = 0.0_DP
      RETURN
    END IF
    
    ! Main GMRES restart loop
    DO WHILE (.NOT. converged .AND. iterations < max_iter)
      restart_count = restart_count + 1
      
      ! Initialize for this restart cycle
      CALL initialize_gmres(workspace, A, b, x, residual_norm)
      
      ! Check initial convergence
      IF (residual_norm <= tol_eff) THEN
        converged = .TRUE.
        EXIT
      END IF
      
      ! GMRES iteration within restart cycle
      DO WHILE (workspace%k < workspace%arnoldi%order .AND. &
               iterations < max_iter .AND. .NOT. converged)
        
        iterations = iterations + 1
        CALL gmres_iterate(workspace, A, converged, residual_norm, info)
        
        IF (residual_norm <= tol_eff) THEN
          converged = .TRUE.
        END IF
      END DO
      
      ! Update solution at end of cycle or convergence
      CALL finalize_gmres_iteration(workspace, x)
    END DO
    
  END SUBROUTINE gmres_solve_structured
  
  SUBROUTINE initialize_gmres(workspace, A, b, x, initial_residual_norm)
    ! Initialize GMRES workspace for new solve or restart
    ! Based on IterativeSolvers.jl init! function
    TYPE(gmres_workspace), INTENT(inout) :: workspace
    REAL(DP), INTENT(in) :: A(:,:)             ! Matrix
    REAL(DP), INTENT(in) :: b(:)               ! RHS vector
    REAL(DP), INTENT(in) :: x(:)               ! Current solution estimate
    REAL(DP), INTENT(out) :: initial_residual_norm
    
    REAL(DP), ALLOCATABLE :: r(:)              ! Residual vector
    INTEGER(I4B) :: n
    
    n = SIZE(b)
    ALLOCATE(r(n))
    
    ! Compute initial residual r = b - A*x
    r = b - MATMUL(A, x)
    
    ! Compute residual norm
    initial_residual_norm = SQRT(DOT_PRODUCT(r, r))
    workspace%residual_norm = initial_residual_norm
    
    ! Initialize first Arnoldi vector V(:,1) = r / ||r||
    workspace%arnoldi%V(:, 1) = r / initial_residual_norm
    
    ! Initialize RHS vector for QR problem: g = ||r|| * e_1
    workspace%rhs_qr = 0.0_DP
    workspace%rhs_qr(1) = initial_residual_norm
    
    ! Reset iteration counter
    workspace%k = 0
    
    ! Reset Givens rotations
    workspace%givens_c = 0.0_DP
    workspace%givens_s = 0.0_DP
    
    DEALLOCATE(r)
    
  END SUBROUTINE initialize_gmres
  
  SUBROUTINE gmres_iterate(workspace, A, converged, residual_norm, info)
    ! Perform one GMRES iteration (expand Arnoldi basis by one vector)
    ! Based on IterativeSolvers.jl iterate function
    TYPE(gmres_workspace), INTENT(inout) :: workspace
    REAL(DP), INTENT(in) :: A(:,:)             ! Matrix
    LOGICAL, INTENT(out) :: converged          ! Convergence flag
    REAL(DP), INTENT(out) :: residual_norm     ! Current residual norm
    INTEGER(I4B), INTENT(out) :: info          ! Return code
    
    REAL(DP), ALLOCATABLE :: w(:)              ! Work vector
    REAL(DP) :: temp
    INTEGER(I4B) :: n, j
    
    info = 0
    converged = .FALSE.
    n = SIZE(workspace%arnoldi%V, 1)
    workspace%k = workspace%k + 1
    
    ALLOCATE(w(n))
    
    ! Arnoldi expansion: w = A * V(:,k)
    w = MATMUL(A, workspace%arnoldi%V(:, workspace%k))
    
    ! Modified Gram-Schmidt orthogonalization
    DO j = 1, workspace%k
      workspace%arnoldi%H(j, workspace%k) = DOT_PRODUCT(w, workspace%arnoldi%V(:, j))
      w = w - workspace%arnoldi%H(j, workspace%k) * workspace%arnoldi%V(:, j)
    END DO
    
    ! Compute new H entry and normalize
    workspace%arnoldi%H(workspace%k + 1, workspace%k) = SQRT(DOT_PRODUCT(w, w))
    
    IF (workspace%arnoldi%H(workspace%k + 1, workspace%k) > TINY(1.0_DP)) THEN
      workspace%arnoldi%V(:, workspace%k + 1) = &
        w / workspace%arnoldi%H(workspace%k + 1, workspace%k)
    ELSE
      ! Breakdown - exact solution found
      converged = .TRUE.
    END IF
    
    ! Apply previous Givens rotations to new column of H
    DO j = 1, workspace%k - 1
      temp = workspace%givens_c(j) * workspace%arnoldi%H(j, workspace%k) + &
             workspace%givens_s(j) * workspace%arnoldi%H(j + 1, workspace%k)
      workspace%arnoldi%H(j + 1, workspace%k) = &
        -workspace%givens_s(j) * workspace%arnoldi%H(j, workspace%k) + &
         workspace%givens_c(j) * workspace%arnoldi%H(j + 1, workspace%k)
      workspace%arnoldi%H(j, workspace%k) = temp
    END DO
    
    ! Compute new Givens rotation
    CALL compute_givens_rotation( &
      workspace%arnoldi%H(workspace%k, workspace%k), &
      workspace%arnoldi%H(workspace%k + 1, workspace%k), &
      workspace%givens_c(workspace%k), &
      workspace%givens_s(workspace%k))
    
    ! Apply new Givens rotation
    workspace%arnoldi%H(workspace%k, workspace%k) = &
      workspace%givens_c(workspace%k) * workspace%arnoldi%H(workspace%k, workspace%k) + &
      workspace%givens_s(workspace%k) * workspace%arnoldi%H(workspace%k + 1, workspace%k)
    workspace%arnoldi%H(workspace%k + 1, workspace%k) = 0.0_DP
    
    ! Update RHS of least squares problem
    temp = workspace%givens_c(workspace%k) * workspace%rhs_qr(workspace%k) + &
           workspace%givens_s(workspace%k) * workspace%rhs_qr(workspace%k + 1)
    workspace%rhs_qr(workspace%k + 1) = &
      -workspace%givens_s(workspace%k) * workspace%rhs_qr(workspace%k) + &
       workspace%givens_c(workspace%k) * workspace%rhs_qr(workspace%k + 1)
    workspace%rhs_qr(workspace%k) = temp
    
    ! Compute residual norm (cheap - just |g_{k+1}|)
    residual_norm = ABS(workspace%rhs_qr(workspace%k + 1))
    workspace%residual_norm = residual_norm
    
    DEALLOCATE(w)
    
  END SUBROUTINE gmres_iterate
  
  SUBROUTINE finalize_gmres_iteration(workspace, x)
    ! Update solution vector after GMRES iterations
    ! Solve least squares problem and update x = x + V*y
    TYPE(gmres_workspace), INTENT(inout) :: workspace
    REAL(DP), INTENT(inout) :: x(:)            ! Solution vector to update
    
    REAL(DP), ALLOCATABLE :: y(:)              ! Least squares solution
    INTEGER(I4B) :: i, j, k
    
    k = workspace%k
    IF (k == 0) RETURN  ! Nothing to do
    
    ALLOCATE(y(k))
    
    ! Solve upper triangular system H*y = g by back substitution
    DO i = k, 1, -1
      y(i) = workspace%rhs_qr(i)
      DO j = i + 1, k
        y(i) = y(i) - workspace%arnoldi%H(i, j) * y(j)
      END DO
      y(i) = y(i) / workspace%arnoldi%H(i, i)
    END DO
    
    ! Update solution: x = x + V*y
    DO i = 1, k
      x = x + y(i) * workspace%arnoldi%V(:, i)
    END DO
    
    DEALLOCATE(y)
    
  END SUBROUTINE finalize_gmres_iteration

  ! ============================================================================
  ! Preconditioned GMRES for sparse matrices (CSR format)
  ! ============================================================================

  SUBROUTINE gmres_solve_structured_preconditioned(workspace, n, row_ptr, col_idx, val, &
                                                   b, x_initial, max_iter, tol, ilu_fac, &
                                                   x, iterations, residual_norm, converged, info, &
                                                   use_preconditioner, precond_side)
    ! Preconditioned GMRES solver for sparse matrices in CSR format
    ! Integrates with existing ILU preconditioning module
    TYPE(gmres_workspace), INTENT(inout) :: workspace
    INTEGER(I4B), INTENT(in) :: n               ! Matrix dimension
    INTEGER(I4B), INTENT(in) :: row_ptr(n+1)    ! CSR row pointers
    INTEGER(I4B), INTENT(in) :: col_idx(:)      ! CSR column indices
    REAL(DP), INTENT(in) :: val(:)              ! CSR values
    REAL(DP), INTENT(in) :: b(:)                ! RHS vector
    REAL(DP), INTENT(in) :: x_initial(:)        ! Initial guess
    INTEGER(I4B), INTENT(in) :: max_iter        ! Maximum iterations
    REAL(DP), INTENT(in) :: tol                 ! Convergence tolerance
    TYPE(ilu_factorization), INTENT(in) :: ilu_fac  ! ILU factorization
    REAL(DP), INTENT(out) :: x(:)               ! Solution vector
    INTEGER(I4B), INTENT(out) :: iterations     ! Total iterations performed
    REAL(DP), INTENT(out) :: residual_norm      ! Final residual norm
    LOGICAL, INTENT(out) :: converged           ! Convergence flag
    INTEGER(I4B), INTENT(out) :: info           ! Return code
    LOGICAL, INTENT(in), OPTIONAL :: use_preconditioner  ! Use ILU preconditioning
    CHARACTER(len=*), INTENT(in), OPTIONAL :: precond_side  ! 'left' or 'right'
    
    ! Local variables
    LOGICAL :: use_precond
    CHARACTER(len=10) :: precond_location
    REAL(DP) :: tol_eff, norm_b
    INTEGER(I4B) :: restart_count
    
    ! Initialize
    info = 0
    converged = .FALSE.
    iterations = 0
    restart_count = 0
    x = x_initial
    
    ! Handle optional arguments
    use_precond = .TRUE.
    IF (PRESENT(use_preconditioner)) use_precond = use_preconditioner
    
    precond_location = 'left'
    IF (PRESENT(precond_side)) precond_location = precond_side
    
    ! Compute effective tolerance
    norm_b = SQRT(DOT_PRODUCT(b, b))
    tol_eff = tol * norm_b
    
    ! Check for zero RHS
    IF (norm_b < TINY(1.0_DP)) THEN
      x = 0.0_DP
      converged = .TRUE.
      residual_norm = 0.0_DP
      RETURN
    END IF
    
    ! Main GMRES restart loop with preconditioning
    DO WHILE (.NOT. converged .AND. iterations < max_iter)
      restart_count = restart_count + 1
      
      ! Initialize for this restart cycle
      CALL initialize_gmres_preconditioned(workspace, n, row_ptr, col_idx, val, &
                                          b, x, use_precond, ilu_fac, residual_norm)
      
      ! Update tolerance based on initial residual (Julia style)
      IF (restart_count == 1) THEN
        tol_eff = MAX(tol * residual_norm, tol)  ! max(reltol * residual.current, abstol)
      END IF
      
      ! Check initial convergence
      IF (residual_norm <= tol_eff) THEN
        converged = .TRUE.
        EXIT
      END IF
      
      ! GMRES iteration within restart cycle
      DO WHILE (workspace%k < workspace%arnoldi%order .AND. &
               iterations < max_iter .AND. .NOT. converged)
        
        iterations = iterations + 1
        CALL gmres_iterate_preconditioned(workspace, n, row_ptr, col_idx, val, &
                                         use_precond, ilu_fac, converged, residual_norm, info)
        
        IF (residual_norm <= tol_eff) THEN
          converged = .TRUE.
        END IF
      END DO
      
      ! Update solution at end of cycle or convergence
      CALL finalize_gmres_iteration(workspace, x)
    END DO
    
  END SUBROUTINE gmres_solve_structured_preconditioned
  
  SUBROUTINE initialize_gmres_preconditioned(workspace, n, row_ptr, col_idx, val, &
                                            b, x, use_precond, ilu_fac, initial_residual_norm)
    ! Initialize GMRES workspace with preconditioning support
    TYPE(gmres_workspace), INTENT(inout) :: workspace
    INTEGER(I4B), INTENT(in) :: n               ! Matrix dimension
    INTEGER(I4B), INTENT(in) :: row_ptr(n+1)    ! CSR row pointers
    INTEGER(I4B), INTENT(in) :: col_idx(:)      ! CSR column indices
    REAL(DP), INTENT(in) :: val(:)              ! CSR values
    REAL(DP), INTENT(in) :: b(:)                ! RHS vector
    REAL(DP), INTENT(in) :: x(:)                ! Current solution estimate
    LOGICAL, INTENT(in) :: use_precond          ! Use preconditioning
    TYPE(ilu_factorization), INTENT(in) :: ilu_fac  ! ILU factorization
    REAL(DP), INTENT(out) :: initial_residual_norm
    
    REAL(DP), ALLOCATABLE :: r(:), z(:)         ! Residual and preconditioned residual
    
    ALLOCATE(r(n), z(n))
    
    ! Compute initial residual r = b - A*x
    CALL matvec_csr_real(n, row_ptr, col_idx, val, x, r)
    r = b - r
    
    ! Apply preconditioning if requested: solve M*z = r
    IF (use_precond .AND. ilu_fac%factorized) THEN
      CALL ilu_solve(ilu_fac, r, z)
    ELSE
      z = r
    END IF
    
    ! Compute residual norm (of preconditioned residual)
    initial_residual_norm = SQRT(DOT_PRODUCT(z, z))
    workspace%residual_norm = initial_residual_norm
    
    
    ! Initialize first Arnoldi vector V(:,1) = z / ||z||
    workspace%arnoldi%V(:, 1) = z / initial_residual_norm
    
    ! Initialize RHS vector for QR problem: g = ||z|| * e_1
    workspace%rhs_qr = 0.0_DP
    workspace%rhs_qr(1) = initial_residual_norm
    
    ! Reset iteration counter
    workspace%k = 0
    
    ! Reset Givens rotations
    workspace%givens_c = 0.0_DP
    workspace%givens_s = 0.0_DP
    
    DEALLOCATE(r, z)
    
  END SUBROUTINE initialize_gmres_preconditioned
  
  SUBROUTINE gmres_iterate_preconditioned(workspace, n, row_ptr, col_idx, val, &
                                          use_precond, ilu_fac, converged, residual_norm, info)
    ! Perform one preconditioned GMRES iteration
    TYPE(gmres_workspace), INTENT(inout) :: workspace
    INTEGER(I4B), INTENT(in) :: n               ! Matrix dimension
    INTEGER(I4B), INTENT(in) :: row_ptr(n+1)    ! CSR row pointers
    INTEGER(I4B), INTENT(in) :: col_idx(:)      ! CSR column indices
    REAL(DP), INTENT(in) :: val(:)              ! CSR values
    LOGICAL, INTENT(in) :: use_precond          ! Use preconditioning
    TYPE(ilu_factorization), INTENT(in) :: ilu_fac  ! ILU factorization
    LOGICAL, INTENT(out) :: converged           ! Convergence flag
    REAL(DP), INTENT(out) :: residual_norm      ! Current residual norm
    INTEGER(I4B), INTENT(out) :: info           ! Return code
    
    REAL(DP), ALLOCATABLE :: w(:), z(:)         ! Work vectors
    REAL(DP) :: temp
    INTEGER(I4B) :: j
    
    info = 0
    converged = .FALSE.
    workspace%k = workspace%k + 1
    
    ALLOCATE(w(n), z(n))
    
    ! Arnoldi expansion: w = A * V(:,k)
    CALL matvec_csr_real(n, row_ptr, col_idx, val, workspace%arnoldi%V(:, workspace%k), w)
    
    ! Apply preconditioning: solve M*z = w
    IF (use_precond .AND. ilu_fac%factorized) THEN
      CALL ilu_solve(ilu_fac, w, z)
    ELSE
      z = w
    END IF
    
    ! Modified Gram-Schmidt orthogonalization
    DO j = 1, workspace%k
      workspace%arnoldi%H(j, workspace%k) = DOT_PRODUCT(z, workspace%arnoldi%V(:, j))
      z = z - workspace%arnoldi%H(j, workspace%k) * workspace%arnoldi%V(:, j)
    END DO
    
    ! Compute new H entry and normalize
    workspace%arnoldi%H(workspace%k + 1, workspace%k) = SQRT(DOT_PRODUCT(z, z))
    
    IF (workspace%arnoldi%H(workspace%k + 1, workspace%k) > TINY(1.0_DP)) THEN
      workspace%arnoldi%V(:, workspace%k + 1) = &
        z / workspace%arnoldi%H(workspace%k + 1, workspace%k)
    ELSE
      ! Breakdown - exact solution found
      converged = .TRUE.
    END IF
    
    ! Apply previous Givens rotations to new column of H
    DO j = 1, workspace%k - 1
      temp = workspace%givens_c(j) * workspace%arnoldi%H(j, workspace%k) + &
             workspace%givens_s(j) * workspace%arnoldi%H(j + 1, workspace%k)
      workspace%arnoldi%H(j + 1, workspace%k) = &
        -workspace%givens_s(j) * workspace%arnoldi%H(j, workspace%k) + &
         workspace%givens_c(j) * workspace%arnoldi%H(j + 1, workspace%k)
      workspace%arnoldi%H(j, workspace%k) = temp
    END DO
    
    ! Compute new Givens rotation
    CALL compute_givens_rotation( &
      workspace%arnoldi%H(workspace%k, workspace%k), &
      workspace%arnoldi%H(workspace%k + 1, workspace%k), &
      workspace%givens_c(workspace%k), &
      workspace%givens_s(workspace%k))
    
    ! Apply new Givens rotation
    workspace%arnoldi%H(workspace%k, workspace%k) = &
      workspace%givens_c(workspace%k) * workspace%arnoldi%H(workspace%k, workspace%k) + &
      workspace%givens_s(workspace%k) * workspace%arnoldi%H(workspace%k + 1, workspace%k)
    workspace%arnoldi%H(workspace%k + 1, workspace%k) = 0.0_DP
    
    ! Update RHS of least squares problem
    temp = workspace%givens_c(workspace%k) * workspace%rhs_qr(workspace%k) + &
           workspace%givens_s(workspace%k) * workspace%rhs_qr(workspace%k + 1)
    workspace%rhs_qr(workspace%k + 1) = &
      -workspace%givens_s(workspace%k) * workspace%rhs_qr(workspace%k) + &
       workspace%givens_c(workspace%k) * workspace%rhs_qr(workspace%k + 1)
    workspace%rhs_qr(workspace%k) = temp
    
    ! Compute residual norm (cheap - just |g_{k+1}|)
    residual_norm = ABS(workspace%rhs_qr(workspace%k + 1))
    workspace%residual_norm = residual_norm
    
    DEALLOCATE(w, z)
    
  END SUBROUTINE gmres_iterate_preconditioned

END MODULE gmres_mod