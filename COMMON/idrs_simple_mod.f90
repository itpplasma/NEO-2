MODULE idrs_simple_mod
  !> Simple unpreconditioned IDR(s) implementation for testing
  !! Based on Algorithm 1 from Sonneveld & van Gijzen (2008)
  
  USE nrtype, ONLY: I4B, DP
  IMPLICIT NONE
  PRIVATE
  
  ! Parameters
  REAL(DP), PARAMETER :: EPS_DP = EPSILON(1.0_DP)
  
  ! Public interface
  PUBLIC :: idrs_simple

CONTAINS

  SUBROUTINE idrs_simple(n, row_ptr, col_idx, values, b, x0, &
                        s, max_iter, tol, x, iter, converged)
    !> Simple IDR(s) without preconditioning
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(IN) :: row_ptr(:), col_idx(:)
    REAL(DP), INTENT(IN) :: values(:), b(:), x0(:)
    INTEGER, INTENT(IN) :: s, max_iter
    REAL(DP), INTENT(IN) :: tol
    REAL(DP), INTENT(OUT) :: x(:)
    INTEGER, INTENT(OUT) :: iter
    LOGICAL, INTENT(OUT) :: converged
    
    ! Local variables
    REAL(DP), ALLOCATABLE :: r(:), v(:), t(:)
    REAL(DP), ALLOCATABLE :: P(:,:), G(:,:), U(:,:), M(:,:)
    REAL(DP), ALLOCATABLE :: f(:), c(:), om(:)
    REAL(DP) :: normb, normr, tolb, omega, rho
    INTEGER :: i, j, k, m
    
    ! Allocate arrays
    ALLOCATE(r(n), v(n), t(n))
    ALLOCATE(P(n,s), G(n,s), U(n,s), M(s,s))
    ALLOCATE(f(s), c(s), om(s+1))
    
    ! Initialize
    x = x0
    iter = 0
    converged = .FALSE.
    
    ! Initial residual r = b - A*x
    CALL matvec(n, row_ptr, col_idx, values, x, r)
    r = b - r
    
    normb = NORM2(b)
    normr = NORM2(r)
    tolb = tol * normb
    
    IF (normr <= tolb) THEN
      converged = .TRUE.
      RETURN
    END IF
    
    ! Initialize shadow space P (random orthonormal)
    CALL init_P(P, n, s)
    
    ! Initialize
    om = 1.0_DP
    U = 0.0_DP
    G = 0.0_DP
    
    ! Main iteration
    DO WHILE (iter < max_iter .AND. .NOT. converged)
      
      ! Compute f = P^T * r
      DO i = 1, s
        f(i) = DOT_PRODUCT(P(:,i), r)
      END DO
      
      ! IDR cycle
      DO k = 1, s
        
        ! Solve for c
        c = 0.0_DP
        DO i = 1, k-1
          c(i) = f(i)
          DO j = 1, i-1
            c(i) = c(i) - M(i,j) * c(j)
          END DO
          c(i) = c(i) / M(i,i)
        END DO
        
        ! Update v
        v = r
        DO j = 1, k-1
          v = v - c(j) * G(:,j)
        END DO
        v = v / om(k)
        
        ! Update U(:,k)
        U(:,k) = om(k) * v
        DO j = 1, k-1
          U(:,k) = U(:,k) + c(j) * U(:,j)
        END DO
        
        ! Compute G(:,k) = A * U(:,k)
        CALL matvec(n, row_ptr, col_idx, values, U(:,k), G(:,k))
        
        ! Update M
        DO i = 1, s
          M(i,k) = DOT_PRODUCT(P(:,i), G(:,k))
        END DO
        
        ! Update solution and residual
        DO m = 1, k
          f(m) = f(m) - M(m,k) * om(k+1)
        END DO
        r = r - om(k+1) * G(:,k)
        x = x + om(k+1) * U(:,k)
        
        iter = iter + 1
        normr = NORM2(r)
        
        IF (normr <= tolb) THEN
          converged = .TRUE.
          EXIT
        END IF
        
        IF (k < s .AND. ABS(f(k+1)) < EPS_DP * normr) THEN
          ! Early termination
          EXIT
        END IF
        
      END DO ! k loop
      
      IF (converged) EXIT
      
      ! Compute new omega
      IF (iter < max_iter) THEN
        
        ! Solve M(1:s,1:s) * c = f
        c = f
        DO i = 1, s
          DO j = 1, i-1
            c(i) = c(i) - M(i,j) * c(j)
          END DO
          c(i) = c(i) / M(i,i)
        END DO
        
        ! Update solution
        v = 0.0_DP
        DO j = 1, s
          v = v + c(j) * U(:,j)
        END DO
        x = x + v
        
        ! Update residual
        t = 0.0_DP
        DO j = 1, s
          t = t + c(j) * G(:,j)
        END DO
        r = r - t
        
        ! Compute omega
        CALL matvec(n, row_ptr, col_idx, values, r, v)
        rho = DOT_PRODUCT(r, v)
        omega = DOT_PRODUCT(r, r) / rho
        
        IF (ABS(omega) < 0.01_DP) omega = 0.01_DP
        
        ! Update x and r
        x = x + omega * r
        r = r - omega * v
        
        om(2:s+1) = omega
        om(1) = 1.0_DP
        
        iter = iter + 1
        normr = NORM2(r)
        
        IF (normr <= tolb) THEN
          converged = .TRUE.
        END IF
        
      END IF
      
    END DO ! Main loop
    
    ! Cleanup
    DEALLOCATE(r, v, t, P, G, U, M, f, c, om)
    
  END SUBROUTINE idrs_simple
  
  SUBROUTINE init_P(P, n, s)
    !> Initialize random orthonormal P
    REAL(DP), INTENT(OUT) :: P(:,:)
    INTEGER, INTENT(IN) :: n, s
    INTEGER :: i, j, k
    REAL(DP) :: norm
    
    ! Random values
    CALL RANDOM_NUMBER(P)
    P = 2.0_DP * P - 1.0_DP
    
    ! Orthonormalize
    DO j = 1, s
      DO k = 1, j-1
        P(:,j) = P(:,j) - DOT_PRODUCT(P(:,j), P(:,k)) * P(:,k)
      END DO
      norm = NORM2(P(:,j))
      IF (norm > EPS_DP) THEN
        P(:,j) = P(:,j) / norm
      END IF
    END DO
    
  END SUBROUTINE init_P
  
  SUBROUTINE matvec(n, row_ptr, col_idx, values, x, y)
    !> Sparse matrix-vector product
    INTEGER, INTENT(IN) :: n, row_ptr(:), col_idx(:)
    REAL(DP), INTENT(IN) :: values(:), x(:)
    REAL(DP), INTENT(OUT) :: y(:)
    INTEGER :: i, j
    
    y = 0.0_DP
    DO i = 1, n
      DO j = row_ptr(i), row_ptr(i+1) - 1
        y(i) = y(i) + values(j) * x(col_idx(j))
      END DO
    END DO
    
  END SUBROUTINE matvec

END MODULE idrs_simple_mod