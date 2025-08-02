MODULE ilu_precond_mod
  ! ILU (Incomplete LU) preconditioner module
  ! Implements ILU(k) factorization with drop tolerance
  ! Provides forward/backward substitution solvers
  
  USE sparse_types_mod, ONLY: dp, long
  USE sparse_utils_mod
  IMPLICIT NONE
  
  PRIVATE
  
  ! Public types
  PUBLIC :: ilu_factorization
  PUBLIC :: ilu_factorization_complex
  
  ! Public interfaces
  PUBLIC :: ilu_factorize
  PUBLIC :: ilu_solve
  PUBLIC :: ilu_free
  
  ! Type for storing ILU factorization (real)
  TYPE :: ilu_factorization
    INTEGER :: n                          ! Matrix dimension
    INTEGER :: L_nnz, U_nnz              ! Number of non-zeros in L and U
    INTEGER, ALLOCATABLE :: L_row_ptr(:) ! CSR row pointers for L
    INTEGER, ALLOCATABLE :: L_col_idx(:) ! CSR column indices for L
    REAL(kind=dp), ALLOCATABLE :: L_val(:) ! CSR values for L
    INTEGER, ALLOCATABLE :: U_row_ptr(:) ! CSR row pointers for U
    INTEGER, ALLOCATABLE :: U_col_idx(:) ! CSR column indices for U
    REAL(kind=dp), ALLOCATABLE :: U_val(:) ! CSR values for U
    LOGICAL :: factorized = .FALSE.      ! Factorization status
  END TYPE ilu_factorization
  
  ! Type for storing ILU factorization (complex)
  TYPE :: ilu_factorization_complex
    INTEGER :: n                          ! Matrix dimension
    INTEGER :: L_nnz, U_nnz              ! Number of non-zeros in L and U
    INTEGER, ALLOCATABLE :: L_row_ptr(:) ! CSR row pointers for L
    INTEGER, ALLOCATABLE :: L_col_idx(:) ! CSR column indices for L
    COMPLEX(kind=dp), ALLOCATABLE :: L_val(:) ! CSR values for L
    INTEGER, ALLOCATABLE :: U_row_ptr(:) ! CSR row pointers for U
    INTEGER, ALLOCATABLE :: U_col_idx(:) ! CSR column indices for U
    COMPLEX(kind=dp), ALLOCATABLE :: U_val(:) ! CSR values for U
    LOGICAL :: factorized = .FALSE.      ! Factorization status
  END TYPE ilu_factorization_complex
  
  ! Generic interfaces
  INTERFACE ilu_factorize
    MODULE PROCEDURE ilu_factorize_real
    MODULE PROCEDURE ilu_factorize_complex
  END INTERFACE ilu_factorize
  
  INTERFACE ilu_solve
    MODULE PROCEDURE ilu_solve_real
    MODULE PROCEDURE ilu_solve_complex
  END INTERFACE ilu_solve
  
  INTERFACE ilu_free
    MODULE PROCEDURE ilu_free_real
    MODULE PROCEDURE ilu_free_complex
  END INTERFACE ilu_free
  
CONTAINS

  !==============================================================================
  ! ILU(k) factorization - Real version
  !==============================================================================
  SUBROUTINE ilu_factorize_real(n, csr_row_ptr, csr_col_idx, csr_val, &
                                 fill_level, drop_tol, ilu_fac, info)
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(IN) :: csr_row_ptr(n+1)
    INTEGER, INTENT(IN) :: csr_col_idx(:)
    REAL(kind=dp), INTENT(IN) :: csr_val(:)
    INTEGER, INTENT(IN) :: fill_level
    REAL(kind=dp), INTENT(IN) :: drop_tol
    TYPE(ilu_factorization), INTENT(OUT) :: ilu_fac
    INTEGER, INTENT(OUT) :: info
    
    ! Local variables
    INTEGER :: i, j, k, kk, jj, jpos
    INTEGER :: nnz_in, nnz_ilu
    REAL(kind=dp) :: pivot, multiplier, sum
    INTEGER, ALLOCATABLE :: ilu_row_ptr(:), ilu_col_idx(:)
    REAL(kind=dp), ALLOCATABLE :: ilu_val(:), w(:)
    INTEGER, ALLOCATABLE :: jw(:)
    INTEGER :: max_nnz
    
    ! Initialize
    info = 0
    ilu_fac%n = n
    ilu_fac%factorized = .FALSE.
    
    ! Count input non-zeros
    nnz_in = csr_row_ptr(n+1) - 1
    
    ! For ILU(k), estimate maximum possible non-zeros
    IF (fill_level == 0) THEN
      max_nnz = nnz_in  ! ILU(0) has same pattern as A
    ELSE
      max_nnz = MIN(n*n, nnz_in * (fill_level + 1) * 2)  ! Rough upper bound
    END IF
    
    ! Allocate ILU storage
    ALLOCATE(ilu_row_ptr(n+1), ilu_col_idx(max_nnz), ilu_val(max_nnz))
    ALLOCATE(w(n), jw(n))
    
    ! Initialize
    ilu_row_ptr(1) = 1
    nnz_ilu = 0
    
    ! Main ILU factorization loop
    DO i = 1, n
      ! Initialize work arrays
      w = 0.0_dp
      jw = 0
      
      ! Copy row i of A into work array
      DO k = csr_row_ptr(i), csr_row_ptr(i+1)-1
        j = csr_col_idx(k)
        w(j) = csr_val(k)
        jw(j) = j  ! Mark position as non-zero
      END DO
      
      ! Eliminate using previous rows
      DO k = 1, i-1
        ! Check if element (i,k) exists
        IF (w(k) /= 0.0_dp) THEN
          ! Find pivot U(k,k) in row k of ILU
          pivot = 0.0_dp
          DO kk = ilu_row_ptr(k), ilu_row_ptr(k+1)-1
            IF (ilu_col_idx(kk) == k) THEN
              pivot = ilu_val(kk)
              EXIT
            END IF
          END DO
          
          IF (ABS(pivot) > 1.0e-15_dp) THEN
            ! Compute multiplier L(i,k) = A(i,k) / U(k,k)
            multiplier = w(k) / pivot
            w(k) = multiplier  ! Store L(i,k)
            
            ! Update row i: A(i,j) = A(i,j) - L(i,k) * U(k,j)
            DO kk = ilu_row_ptr(k), ilu_row_ptr(k+1)-1
              j = ilu_col_idx(kk)
              IF (j > k) THEN
                ! For ILU(0), only update if (i,j) exists in original pattern
                IF (fill_level == 0) THEN
                  IF (jw(j) /= 0) THEN
                    w(j) = w(j) - multiplier * ilu_val(kk)
                  END IF
                ELSE
                  ! For ILU(k), allow fill-in
                  w(j) = w(j) - multiplier * ilu_val(kk)
                  IF (jw(j) == 0) jw(j) = j  ! Mark new fill-in
                END IF
              END IF
            END DO
          ELSE
            ! Zero pivot - set L(i,k) = 0
            w(k) = 0.0_dp
          END IF
        END IF
      END DO
      
      ! Check diagonal element
      IF (ABS(w(i)) < 1.0e-15_dp) THEN
        IF (jw(i) == 0) THEN
          info = -3  ! Structural zero on diagonal
          DEALLOCATE(ilu_row_ptr, ilu_col_idx, ilu_val, w, jw)
          RETURN
        END IF
        info = 1  ! Near-zero pivot warning
        w(i) = 1.0e-12_dp  ! Replace with small value
      END IF
      
      ! Store row i in ILU structure, applying drop tolerance
      DO j = 1, n
        IF (jw(j) /= 0 .AND. ABS(w(j)) > drop_tol) THEN
          nnz_ilu = nnz_ilu + 1
          IF (nnz_ilu > max_nnz) THEN
            info = -2  ! Insufficient storage
            DEALLOCATE(ilu_row_ptr, ilu_col_idx, ilu_val, w, jw)
            RETURN
          END IF
          ilu_col_idx(nnz_ilu) = j
          ilu_val(nnz_ilu) = w(j)
        END IF
      END DO
      
      ilu_row_ptr(i+1) = nnz_ilu + 1
    END DO
    
    ! Split into L and U
    CALL split_lu_from_ilu(n, nnz_ilu, ilu_row_ptr, ilu_col_idx, ilu_val, ilu_fac)
    
    ! Cleanup
    DEALLOCATE(ilu_row_ptr, ilu_col_idx, ilu_val, w, jw)
    
    ilu_fac%factorized = .TRUE.
    
  END SUBROUTINE ilu_factorize_real
  
  !==============================================================================
  ! ILU(k) factorization - Complex version
  !==============================================================================
  SUBROUTINE ilu_factorize_complex(n, csr_row_ptr, csr_col_idx, csr_val, &
                                    fill_level, drop_tol, ilu_fac, info)
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(IN) :: csr_row_ptr(n+1)
    INTEGER, INTENT(IN) :: csr_col_idx(:)
    COMPLEX(kind=dp), INTENT(IN) :: csr_val(:)
    INTEGER, INTENT(IN) :: fill_level
    REAL(kind=dp), INTENT(IN) :: drop_tol
    TYPE(ilu_factorization_complex), INTENT(OUT) :: ilu_fac
    INTEGER, INTENT(OUT) :: info
    
    ! Local variables
    INTEGER :: i, j, k, kk, jj, jpos
    INTEGER :: nnz_in, nnz_ilu
    COMPLEX(kind=dp) :: pivot, multiplier
    INTEGER, ALLOCATABLE :: ilu_row_ptr(:), ilu_col_idx(:)
    COMPLEX(kind=dp), ALLOCATABLE :: ilu_val(:), w(:)
    INTEGER, ALLOCATABLE :: jw(:)
    INTEGER :: max_nnz
    
    ! Initialize
    info = 0
    ilu_fac%n = n
    ilu_fac%factorized = .FALSE.
    
    ! Count input non-zeros
    nnz_in = csr_row_ptr(n+1) - 1
    
    ! For ILU(k), estimate maximum possible non-zeros
    IF (fill_level == 0) THEN
      max_nnz = nnz_in  ! ILU(0) has same pattern as A
    ELSE
      max_nnz = MIN(n*n, nnz_in * (fill_level + 1) * 2)  ! Rough upper bound
    END IF
    
    ! Allocate ILU storage
    ALLOCATE(ilu_row_ptr(n+1), ilu_col_idx(max_nnz), ilu_val(max_nnz))
    ALLOCATE(w(n), jw(n))
    
    ! Initialize
    ilu_row_ptr(1) = 1
    nnz_ilu = 0
    
    ! Main ILU factorization loop
    DO i = 1, n
      ! Initialize work arrays
      w = (0.0_dp, 0.0_dp)
      jw = 0
      
      ! Copy row i of A into work array
      DO k = csr_row_ptr(i), csr_row_ptr(i+1)-1
        j = csr_col_idx(k)
        w(j) = csr_val(k)
        jw(j) = j  ! Mark position as non-zero
      END DO
      
      ! Eliminate using previous rows
      DO k = 1, i-1
        ! Check if element (i,k) exists
        IF (ABS(w(k)) > 0.0_dp) THEN
          ! Find pivot U(k,k) in row k of ILU
          pivot = (0.0_dp, 0.0_dp)
          DO kk = ilu_row_ptr(k), ilu_row_ptr(k+1)-1
            IF (ilu_col_idx(kk) == k) THEN
              pivot = ilu_val(kk)
              EXIT
            END IF
          END DO
          
          IF (ABS(pivot) > 1.0e-15_dp) THEN
            ! Compute multiplier L(i,k) = A(i,k) / U(k,k)
            multiplier = w(k) / pivot
            w(k) = multiplier  ! Store L(i,k)
            
            ! Update row i: A(i,j) = A(i,j) - L(i,k) * U(k,j)
            DO kk = ilu_row_ptr(k), ilu_row_ptr(k+1)-1
              j = ilu_col_idx(kk)
              IF (j > k) THEN
                ! For ILU(0), only update if (i,j) exists in original pattern
                IF (fill_level == 0) THEN
                  IF (jw(j) /= 0) THEN
                    w(j) = w(j) - multiplier * ilu_val(kk)
                  END IF
                ELSE
                  ! For ILU(k), allow fill-in
                  w(j) = w(j) - multiplier * ilu_val(kk)
                  IF (jw(j) == 0) jw(j) = j  ! Mark new fill-in
                END IF
              END IF
            END DO
          ELSE
            ! Zero pivot - set L(i,k) = 0
            w(k) = (0.0_dp, 0.0_dp)
          END IF
        END IF
      END DO
      
      ! Check diagonal element
      IF (ABS(w(i)) < 1.0e-15_dp) THEN
        IF (jw(i) == 0) THEN
          info = -3  ! Structural zero on diagonal
          DEALLOCATE(ilu_row_ptr, ilu_col_idx, ilu_val, w, jw)
          RETURN
        END IF
        info = 1  ! Near-zero pivot warning
        w(i) = (1.0e-12_dp, 0.0_dp)  ! Replace with small value
      END IF
      
      ! Store row i in ILU structure, applying drop tolerance
      DO j = 1, n
        IF (jw(j) /= 0 .AND. ABS(w(j)) > drop_tol) THEN
          nnz_ilu = nnz_ilu + 1
          IF (nnz_ilu > max_nnz) THEN
            info = -2  ! Insufficient storage
            DEALLOCATE(ilu_row_ptr, ilu_col_idx, ilu_val, w, jw)
            RETURN
          END IF
          ilu_col_idx(nnz_ilu) = j
          ilu_val(nnz_ilu) = w(j)
        END IF
      END DO
      
      ilu_row_ptr(i+1) = nnz_ilu + 1
    END DO
    
    ! Split into L and U
    CALL split_lu_from_ilu_complex(n, nnz_ilu, ilu_row_ptr, ilu_col_idx, ilu_val, ilu_fac)
    
    ! Cleanup
    DEALLOCATE(ilu_row_ptr, ilu_col_idx, ilu_val, w, jw)
    
    ilu_fac%factorized = .TRUE.
    
  END SUBROUTINE ilu_factorize_complex
  
  !==============================================================================
  ! Forward/backward substitution - Real version
  !==============================================================================
  SUBROUTINE ilu_solve_real(ilu_fac, b, x)
    TYPE(ilu_factorization), INTENT(IN) :: ilu_fac
    REAL(kind=dp), INTENT(IN) :: b(:)
    REAL(kind=dp), INTENT(OUT) :: x(:)
    
    INTEGER :: i, k
    REAL(kind=dp), ALLOCATABLE :: y(:)
    INTEGER :: n
    
    n = ilu_fac%n
    ALLOCATE(y(n))
    
    ! Forward substitution: L*y = b
    DO i = 1, n
      y(i) = b(i)
      DO k = ilu_fac%L_row_ptr(i), ilu_fac%L_row_ptr(i+1)-1
        IF (ilu_fac%L_col_idx(k) < i) THEN
          y(i) = y(i) - ilu_fac%L_val(k) * y(ilu_fac%L_col_idx(k))
        END IF
      END DO
    END DO
    
    ! Backward substitution: U*x = y
    DO i = n, 1, -1
      x(i) = y(i)
      DO k = ilu_fac%U_row_ptr(i), ilu_fac%U_row_ptr(i+1)-1
        IF (ilu_fac%U_col_idx(k) > i) THEN
          x(i) = x(i) - ilu_fac%U_val(k) * x(ilu_fac%U_col_idx(k))
        ELSE IF (ilu_fac%U_col_idx(k) == i) THEN
          x(i) = x(i) / ilu_fac%U_val(k)
        END IF
      END DO
    END DO
    
    DEALLOCATE(y)
    
  END SUBROUTINE ilu_solve_real
  
  !==============================================================================
  ! Forward/backward substitution - Complex version
  !==============================================================================
  SUBROUTINE ilu_solve_complex(ilu_fac, b, x)
    TYPE(ilu_factorization_complex), INTENT(IN) :: ilu_fac
    COMPLEX(kind=dp), INTENT(IN) :: b(:)
    COMPLEX(kind=dp), INTENT(OUT) :: x(:)
    
    INTEGER :: i, k
    COMPLEX(kind=dp), ALLOCATABLE :: y(:)
    INTEGER :: n
    
    n = ilu_fac%n
    ALLOCATE(y(n))
    
    ! Forward substitution: L*y = b
    DO i = 1, n
      y(i) = b(i)
      DO k = ilu_fac%L_row_ptr(i), ilu_fac%L_row_ptr(i+1)-1
        IF (ilu_fac%L_col_idx(k) < i) THEN
          y(i) = y(i) - ilu_fac%L_val(k) * y(ilu_fac%L_col_idx(k))
        END IF
      END DO
    END DO
    
    ! Backward substitution: U*x = y
    DO i = n, 1, -1
      x(i) = y(i)
      DO k = ilu_fac%U_row_ptr(i), ilu_fac%U_row_ptr(i+1)-1
        IF (ilu_fac%U_col_idx(k) > i) THEN
          x(i) = x(i) - ilu_fac%U_val(k) * x(ilu_fac%U_col_idx(k))
        ELSE IF (ilu_fac%U_col_idx(k) == i) THEN
          x(i) = x(i) / ilu_fac%U_val(k)
        END IF
      END DO
    END DO
    
    DEALLOCATE(y)
    
  END SUBROUTINE ilu_solve_complex
  
  !==============================================================================
  ! Free ILU factorization - Real version
  !==============================================================================
  SUBROUTINE ilu_free_real(ilu_fac)
    TYPE(ilu_factorization), INTENT(INOUT) :: ilu_fac
    
    IF (ALLOCATED(ilu_fac%L_row_ptr)) DEALLOCATE(ilu_fac%L_row_ptr)
    IF (ALLOCATED(ilu_fac%L_col_idx)) DEALLOCATE(ilu_fac%L_col_idx)
    IF (ALLOCATED(ilu_fac%L_val)) DEALLOCATE(ilu_fac%L_val)
    IF (ALLOCATED(ilu_fac%U_row_ptr)) DEALLOCATE(ilu_fac%U_row_ptr)
    IF (ALLOCATED(ilu_fac%U_col_idx)) DEALLOCATE(ilu_fac%U_col_idx)
    IF (ALLOCATED(ilu_fac%U_val)) DEALLOCATE(ilu_fac%U_val)
    
    ilu_fac%factorized = .FALSE.
    ilu_fac%L_nnz = 0
    ilu_fac%U_nnz = 0
    
  END SUBROUTINE ilu_free_real
  
  !==============================================================================
  ! Free ILU factorization - Complex version
  !==============================================================================
  SUBROUTINE ilu_free_complex(ilu_fac)
    TYPE(ilu_factorization_complex), INTENT(INOUT) :: ilu_fac
    
    IF (ALLOCATED(ilu_fac%L_row_ptr)) DEALLOCATE(ilu_fac%L_row_ptr)
    IF (ALLOCATED(ilu_fac%L_col_idx)) DEALLOCATE(ilu_fac%L_col_idx)
    IF (ALLOCATED(ilu_fac%L_val)) DEALLOCATE(ilu_fac%L_val)
    IF (ALLOCATED(ilu_fac%U_row_ptr)) DEALLOCATE(ilu_fac%U_row_ptr)
    IF (ALLOCATED(ilu_fac%U_col_idx)) DEALLOCATE(ilu_fac%U_col_idx)
    IF (ALLOCATED(ilu_fac%U_val)) DEALLOCATE(ilu_fac%U_val)
    
    ilu_fac%factorized = .FALSE.
    ilu_fac%L_nnz = 0
    ilu_fac%U_nnz = 0
    
  END SUBROUTINE ilu_free_complex
  
  !==============================================================================
  ! Helper function to determine level of fill
  !==============================================================================
  LOGICAL FUNCTION level_of_fill(i, j, k, max_level)
    INTEGER, INTENT(IN) :: i, j, k, max_level
    
    ! Simple level-of-fill criterion
    ! For ILU(k), allow fill if path length <= k
    ! This is a simplified version - actual implementation could be more sophisticated
    level_of_fill = .TRUE.  ! For now, allow fill up to max_level
    
  END FUNCTION level_of_fill
  
  !==============================================================================
  ! Split ILU factorization into separate L and U - Real version
  !==============================================================================
  SUBROUTINE split_lu_from_ilu(n, nnz, ilu_row_ptr, ilu_col_idx, ilu_val, ilu_fac)
    INTEGER, INTENT(IN) :: n, nnz
    INTEGER, INTENT(IN) :: ilu_row_ptr(n+1), ilu_col_idx(nnz)
    REAL(kind=dp), INTENT(IN) :: ilu_val(nnz)
    TYPE(ilu_factorization), INTENT(INOUT) :: ilu_fac
    
    INTEGER :: i, j, k, l_count, u_count
    
    ! Count entries in L and U
    l_count = 0
    u_count = 0
    
    DO i = 1, n
      DO k = ilu_row_ptr(i), ilu_row_ptr(i+1)-1
        IF (ilu_col_idx(k) < i) THEN
          l_count = l_count + 1  ! Strict lower
        ELSE
          u_count = u_count + 1  ! Diagonal and upper
        END IF
      END DO
      l_count = l_count + 1  ! Unit diagonal for L
    END DO
    
    ilu_fac%L_nnz = l_count
    ilu_fac%U_nnz = u_count
    
    ! Allocate storage
    ALLOCATE(ilu_fac%L_row_ptr(n+1), ilu_fac%L_col_idx(l_count), ilu_fac%L_val(l_count))
    ALLOCATE(ilu_fac%U_row_ptr(n+1), ilu_fac%U_col_idx(u_count), ilu_fac%U_val(u_count))
    
    ! Fill L and U
    l_count = 0
    u_count = 0
    
    DO i = 1, n
      ilu_fac%L_row_ptr(i) = l_count + 1
      ilu_fac%U_row_ptr(i) = u_count + 1
      
      ! Process row i
      DO k = ilu_row_ptr(i), ilu_row_ptr(i+1)-1
        j = ilu_col_idx(k)
        IF (j < i) THEN
          ! Strict lower - goes to L
          l_count = l_count + 1
          ilu_fac%L_col_idx(l_count) = j
          ilu_fac%L_val(l_count) = ilu_val(k)
        ELSE
          ! Diagonal and upper - goes to U
          u_count = u_count + 1
          ilu_fac%U_col_idx(u_count) = j
          ilu_fac%U_val(u_count) = ilu_val(k)
        END IF
      END DO
      
      ! Add unit diagonal to L
      l_count = l_count + 1
      ilu_fac%L_col_idx(l_count) = i
      ilu_fac%L_val(l_count) = 1.0_dp
    END DO
    
    ilu_fac%L_row_ptr(n+1) = l_count + 1
    ilu_fac%U_row_ptr(n+1) = u_count + 1
    
  END SUBROUTINE split_lu_from_ilu
  
  !==============================================================================
  ! Split combined LU into separate L and U - Real version (legacy)
  !==============================================================================
  SUBROUTINE split_lu(n, nnz, ju, jlu, alu, ilu_fac)
    INTEGER, INTENT(IN) :: n, nnz
    INTEGER, INTENT(IN) :: ju(n+1), jlu(nnz)
    REAL(kind=dp), INTENT(IN) :: alu(nnz)
    TYPE(ilu_factorization), INTENT(INOUT) :: ilu_fac
    
    INTEGER :: i, j, k, l_count, u_count
    
    ! Count entries in L and U
    l_count = n  ! Diagonal entries (unit diagonal for L)
    u_count = 0
    
    DO i = 1, n
      DO k = ju(i), ju(i+1)-1
        IF (jlu(k) < i) THEN
          l_count = l_count + 1
        ELSE
          u_count = u_count + 1
        END IF
      END DO
    END DO
    
    ilu_fac%L_nnz = l_count
    ilu_fac%U_nnz = u_count
    
    ! Allocate storage
    ALLOCATE(ilu_fac%L_row_ptr(n+1), ilu_fac%L_col_idx(l_count), ilu_fac%L_val(l_count))
    ALLOCATE(ilu_fac%U_row_ptr(n+1), ilu_fac%U_col_idx(u_count), ilu_fac%U_val(u_count))
    
    ! Fill L and U
    l_count = 0
    u_count = 0
    
    DO i = 1, n
      ilu_fac%L_row_ptr(i) = l_count + 1
      ilu_fac%U_row_ptr(i) = u_count + 1
      
      ! L part (including unit diagonal)
      DO k = ju(i), ju(i+1)-1
        IF (jlu(k) < i) THEN
          l_count = l_count + 1
          ilu_fac%L_col_idx(l_count) = jlu(k)
          ilu_fac%L_val(l_count) = alu(k)
        END IF
      END DO
      
      ! Unit diagonal for L
      l_count = l_count + 1
      ilu_fac%L_col_idx(l_count) = i
      ilu_fac%L_val(l_count) = 1.0_dp
      
      ! U part
      DO k = ju(i), ju(i+1)-1
        IF (jlu(k) >= i) THEN
          u_count = u_count + 1
          ilu_fac%U_col_idx(u_count) = jlu(k)
          ilu_fac%U_val(u_count) = alu(k)
        END IF
      END DO
    END DO
    
    ilu_fac%L_row_ptr(n+1) = l_count + 1
    ilu_fac%U_row_ptr(n+1) = u_count + 1
    
  END SUBROUTINE split_lu
  
  !==============================================================================
  ! Split ILU factorization into separate L and U - Complex version
  !==============================================================================
  SUBROUTINE split_lu_from_ilu_complex(n, nnz, ilu_row_ptr, ilu_col_idx, ilu_val, ilu_fac)
    INTEGER, INTENT(IN) :: n, nnz
    INTEGER, INTENT(IN) :: ilu_row_ptr(n+1), ilu_col_idx(nnz)
    COMPLEX(kind=dp), INTENT(IN) :: ilu_val(nnz)
    TYPE(ilu_factorization_complex), INTENT(INOUT) :: ilu_fac
    
    INTEGER :: i, j, k, l_count, u_count
    
    ! Count entries in L and U
    l_count = 0
    u_count = 0
    
    DO i = 1, n
      DO k = ilu_row_ptr(i), ilu_row_ptr(i+1)-1
        IF (ilu_col_idx(k) < i) THEN
          l_count = l_count + 1  ! Strict lower
        ELSE
          u_count = u_count + 1  ! Diagonal and upper
        END IF
      END DO
      l_count = l_count + 1  ! Unit diagonal for L
    END DO
    
    ilu_fac%L_nnz = l_count
    ilu_fac%U_nnz = u_count
    
    ! Allocate storage
    ALLOCATE(ilu_fac%L_row_ptr(n+1), ilu_fac%L_col_idx(l_count), ilu_fac%L_val(l_count))
    ALLOCATE(ilu_fac%U_row_ptr(n+1), ilu_fac%U_col_idx(u_count), ilu_fac%U_val(u_count))
    
    ! Fill L and U
    l_count = 0
    u_count = 0
    
    DO i = 1, n
      ilu_fac%L_row_ptr(i) = l_count + 1
      ilu_fac%U_row_ptr(i) = u_count + 1
      
      ! Process row i
      DO k = ilu_row_ptr(i), ilu_row_ptr(i+1)-1
        j = ilu_col_idx(k)
        IF (j < i) THEN
          ! Strict lower - goes to L
          l_count = l_count + 1
          ilu_fac%L_col_idx(l_count) = j
          ilu_fac%L_val(l_count) = ilu_val(k)
        ELSE
          ! Diagonal and upper - goes to U
          u_count = u_count + 1
          ilu_fac%U_col_idx(u_count) = j
          ilu_fac%U_val(u_count) = ilu_val(k)
        END IF
      END DO
      
      ! Add unit diagonal to L
      l_count = l_count + 1
      ilu_fac%L_col_idx(l_count) = i
      ilu_fac%L_val(l_count) = (1.0_dp, 0.0_dp)
    END DO
    
    ilu_fac%L_row_ptr(n+1) = l_count + 1
    ilu_fac%U_row_ptr(n+1) = u_count + 1
    
  END SUBROUTINE split_lu_from_ilu_complex
  
  !==============================================================================
  ! Split combined LU into separate L and U - Complex version
  !==============================================================================
  SUBROUTINE split_lu_complex(n, nnz, ju, jlu, alu, ilu_fac)
    INTEGER, INTENT(IN) :: n, nnz
    INTEGER, INTENT(IN) :: ju(n+1), jlu(nnz)
    COMPLEX(kind=dp), INTENT(IN) :: alu(nnz)
    TYPE(ilu_factorization_complex), INTENT(INOUT) :: ilu_fac
    
    INTEGER :: i, j, k, l_count, u_count
    
    ! Count entries in L and U
    l_count = n  ! Diagonal entries
    u_count = 0
    
    DO i = 1, n
      DO k = ju(i), ju(i+1)-1
        IF (jlu(k) < i) THEN
          l_count = l_count + 1
        ELSE
          u_count = u_count + 1
        END IF
      END DO
    END DO
    
    ilu_fac%L_nnz = l_count
    ilu_fac%U_nnz = u_count
    
    ! Allocate storage
    ALLOCATE(ilu_fac%L_row_ptr(n+1), ilu_fac%L_col_idx(l_count), ilu_fac%L_val(l_count))
    ALLOCATE(ilu_fac%U_row_ptr(n+1), ilu_fac%U_col_idx(u_count), ilu_fac%U_val(u_count))
    
    ! Fill L and U
    l_count = 0
    u_count = 0
    
    DO i = 1, n
      ilu_fac%L_row_ptr(i) = l_count + 1
      ilu_fac%U_row_ptr(i) = u_count + 1
      
      ! L part
      DO k = ju(i), ju(i+1)-1
        IF (jlu(k) < i) THEN
          l_count = l_count + 1
          ilu_fac%L_col_idx(l_count) = jlu(k)
          ilu_fac%L_val(l_count) = alu(k)
        END IF
      END DO
      
      ! Unit diagonal for L
      l_count = l_count + 1
      ilu_fac%L_col_idx(l_count) = i
      ilu_fac%L_val(l_count) = (1.0_dp, 0.0_dp)
      
      ! U part
      DO k = ju(i), ju(i+1)-1
        IF (jlu(k) >= i) THEN
          u_count = u_count + 1
          ilu_fac%U_col_idx(u_count) = jlu(k)
          ilu_fac%U_val(u_count) = alu(k)
        END IF
      END DO
    END DO
    
    ilu_fac%L_row_ptr(n+1) = l_count + 1
    ilu_fac%U_row_ptr(n+1) = u_count + 1
    
  END SUBROUTINE split_lu_complex
  
END MODULE ilu_precond_mod