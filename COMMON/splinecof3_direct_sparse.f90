!> Direct sparse implementation that builds matrix in COO format and converts to CSC
module splinecof3_direct_sparse_mod
  use nrtype, only : I4B, DP
  use sparse_mod, only: sparse_solve
  use inter_interfaces, only: calc_opt_lambda3
  implicit none
  
  private
  public :: splinecof3_direct_sparse, splinecof3_direct_sparse_get_coo
  
  ! Module variables to store COO matrix for inspection
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE, SAVE :: last_irow_coo, last_icol_coo
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE :: last_val_coo, last_rhs_coo
  INTEGER(I4B), SAVE :: last_nnz = 0, last_n = 0
  
contains

  !> Direct sparse implementation matching splinecof3_a algorithm
  SUBROUTINE splinecof3_direct_sparse(x, y, c1, cn, lambda1, indx, sw1, sw2, &
       a, b, c, d, m, f)
    REAL(DP),                   INTENT(INOUT) :: c1, cn
    REAL(DP),     DIMENSION(:), INTENT(IN)    :: x, y, lambda1
    INTEGER(I4B), DIMENSION(:), INTENT(IN)    :: indx
    REAL(DP),     DIMENSION(:), INTENT(OUT)   :: a, b, c, d
    INTEGER(I4B),               INTENT(IN)    :: sw1, sw2
    REAL(DP),                   INTENT(IN)    :: m
    INTERFACE
       FUNCTION f(x,m)
         use nrtype, only : DP
         IMPLICIT NONE
         REAL(DP), INTENT(IN) :: x, m
         REAL(DP)             :: f
       END FUNCTION f
    END INTERFACE

    ! Local variables
    INTEGER(I4B) :: len_indx, VAR, size_dimension
    INTEGER(I4B) :: i, j, k, l, ii, ie, nnz, idx, max_nnz
    INTEGER(I4B) :: i_alloc, mu1, mu2, nu1, nu2, sig1, sig2, rho1, rho2
    INTEGER(I4B) :: nrow, ncol, pos, len_x
    REAL(DP) :: h, h_j, x_h
    REAL(DP) :: help_a, help_b, help_c, help_d, help_i, help_inh
    REAL(DP), DIMENSION(:), ALLOCATABLE :: lambda, omega, inh
    ! COO format arrays
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: irow_coo, icol_coo
    REAL(DP), DIMENSION(:), ALLOCATABLE :: val_coo
    ! CSC format arrays
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: irow_csc, pcol_csc
    REAL(DP), DIMENSION(:), ALLOCATABLE :: val_csc
    ! Helper arrays
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: col_count
    character(200) :: error_message

    ! Initialize variables
    VAR = 7
    len_x = SIZE(x)
    len_indx = SIZE(indx)
    size_dimension = VAR * len_indx - 2
    nrow = size_dimension
    ncol = size_dimension

    ! Validation checks
    if ( .NOT. ( size(x) == size(y) ) ) then
      write (*,*) 'splinecof3_direct_sparse: assertion 1 failed'
      stop 'program terminated'
    end if
    if ( .NOT. ( size(a) == size(b) .AND. size(a) == size(c) &
         .AND.   size(a) == size(d) .AND. size(a) == size(indx) &
         .AND.   size(a) == size(lambda1) ) ) then
      write (*,*) 'splinecof3_direct_sparse: assertion 2 failed'
      stop 'program terminated'
    end if

    do i = 1, len_x-1
      if (x(i) >= x(i+1)) then
        print *, 'SPLINECOF3_DIRECT_SPARSE: error i, x(i), x(i+1)', &
             i, x(i), x(i+1)
        stop 'SPLINECOF3_DIRECT_SPARSE: error  wrong order of x(i)'
      end if
    end do
    do i = 1, len_indx-1
      if (indx(i) < 1) then
        print *, 'SPLINECOF3_DIRECT_SPARSE: error i, indx(i)', i, indx(i)
        stop 'SPLINECOF3_DIRECT_SPARSE: error  indx(i) < 1'
      end if
      if (indx(i) >= indx(i+1)) then
        print *, 'SPLINECOF3_DIRECT_SPARSE: error i, indx(i), indx(i+1)', &
              i, indx(i), indx(i+1)
        stop 'SPLINECOF3_DIRECT_SPARSE: error  wrong order of indx(i)'
      end if
      if (indx(i) > len_x) then
        print *, 'SPLINECOF3_DIRECT_SPARSE: error i, indx(i), indx(i+1)', &
              i, indx(i), indx(i+1)
        stop 'SPLINECOF3_DIRECT_SPARSE: error  indx(i) > len_x'
      end if
    end do
    if (indx(len_indx) < 1) then
      print *, 'SPLINECOF3_DIRECT_SPARSE: error len_indx, indx(len_indx)', &
            len_indx, indx(len_indx)
      stop 'SPLINECOF3_DIRECT_SPARSE: error  indx(max) < 1'
    end if
    if (indx(len_indx) > len_x) then
      print *, 'SPLINECOF3_DIRECT_SPARSE: error len_indx, indx(len_indx)', &
            len_indx, indx(len_indx)
      stop 'SPLINECOF3_DIRECT_SPARSE: error  indx(max) > len_x'
    end if

    if (sw1 == sw2) then
      stop 'SPLINECOF3_DIRECT_SPARSE: error  two identical boundary conditions'
    end if

    ! Allocate work arrays
    ALLOCATE(lambda(len_indx), omega(len_indx), inh(size_dimension), &
             stat = i_alloc, errmsg=error_message)
    if(i_alloc /= 0) then
      write(*,*) 'splinecof3_direct_sparse: Allocation failed:', trim(error_message)
      stop
    end if

    ! Process boundary conditions
    IF (DABS(c1) > 1.0E30) THEN
      c1 = 0.0D0
    END IF
    IF (DABS(cn) > 1.0E30) THEN
      cn = 0.0D0
    END IF

    ! Calculate optimal weights for smoothing (lambda)
    IF ( MAXVAL(lambda1) < 0.0D0 ) THEN
      CALL calc_opt_lambda3(x, y, omega)
    ELSE
      omega  = lambda1
    END IF
    lambda = 1.0D0 - omega
    
    ! Initialize RHS vector
    inh = 0.0D0

    ! Set boundary condition switches
    mu1  = 0; mu2  = 0
    nu1  = 0; nu2  = 0
    sig1 = 0; sig2 = 0
    rho1 = 0; rho2 = 0

    SELECT CASE(sw1)
    CASE(1); mu1 = 1
    CASE(2); nu1 = 1
    CASE(3); sig1 = 1
    CASE(4); rho1 = 1
    END SELECT

    SELECT CASE(sw2)
    CASE(1); mu2 = 1
    CASE(2); nu2 = 1
    CASE(3); sig2 = 1
    CASE(4); rho2 = 1
    END SELECT

    ! Estimate maximum non-zeros (extremely conservative to prevent overflow)
    ! Each equation can have up to ~15 non-zeros, size_dimension equations total
    max_nnz = 20 * size_dimension
    
    ! Allocate COO format arrays
    ALLOCATE(irow_coo(max_nnz), icol_coo(max_nnz), val_coo(max_nnz), &
             stat = i_alloc)
    if(i_alloc /= 0) stop 'Allocation for COO arrays failed!'
    
    ! Build the sparse matrix in COO format
    idx = 0
    i = 0

    ! Boundary condition 1
    i = i + 1
    ! For sparse matrices, only add non-zero entries
    IF (mu1 /= 0) THEN
       idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = 2; val_coo(idx) = DBLE(mu1)
    END IF
    IF (nu1 /= 0) THEN
       idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = 3; val_coo(idx) = DBLE(nu1)
    END IF
    IF (sig1 /= 0) THEN
       idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = (len_indx-1)*VAR + 2; val_coo(idx) = DBLE(sig1)
    END IF
    IF (rho1 /= 0) THEN
       idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = (len_indx-1)*VAR + 3; val_coo(idx) = DBLE(rho1)
    END IF
    inh(i) = c1

    ! Main loop over intervals
    DO j = 1, VAR*(len_indx-1)-1, VAR
       ii = indx((j-1)/VAR+1)
       ie = indx((j-1)/VAR+2) - 1
       h  = x(ie+1) - x(ii)

       ! Continuity conditions - A_i, B_i, C_i
       ! A_i continuity
       i = i + 1
       idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j; val_coo(idx) = 1.0D0
       idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j+1; val_coo(idx) = h
       idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j+2; val_coo(idx) = h*h
       idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j+3; val_coo(idx) = h*h*h
       idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j+VAR; val_coo(idx) = -1.0D0

       ! B_i continuity
       i = i + 1
       idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j+1; val_coo(idx) = 1.0D0
       idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j+2; val_coo(idx) = 2.0D0*h
       idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j+3; val_coo(idx) = 3.0D0*h*h
       idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j+VAR+1; val_coo(idx) = -1.0D0

       ! C_i continuity
       i = i + 1
       idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j+2; val_coo(idx) = 1.0D0
       idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j+3; val_coo(idx) = 3.0D0*h
       idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j+VAR+2; val_coo(idx) = -1.0D0

       ! Fitting conditions - compute coefficients
       help_a = 0.0D0; help_b = 0.0D0; help_c = 0.0D0; help_d = 0.0D0
       help_i = 0.0D0
       
       DO l = ii, ie
          h_j = x(l) - x(ii)
          x_h = f(x(l),m) * f(x(l),m)
          help_a = help_a + x_h
          help_b = help_b + h_j * x_h
          help_c = help_c + h_j * h_j * x_h
          help_d = help_d + h_j * h_j * h_j * x_h
          help_i = help_i + f(x(l),m) * y(l)
       END DO

       ! delta a_i
       i = i + 1
       IF (ABS(help_a) > 1D-15) THEN
          idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j; val_coo(idx) = omega((j-1)/VAR+1) * help_a
       END IF
       IF (ABS(help_b) > 1D-15) THEN
          idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j+1; val_coo(idx) = omega((j-1)/VAR+1) * help_b
       END IF
       IF (ABS(help_c) > 1D-15) THEN
          idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j+2; val_coo(idx) = omega((j-1)/VAR+1) * help_c
       END IF
       IF (ABS(help_d) > 1D-15) THEN
          idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j+3; val_coo(idx) = omega((j-1)/VAR+1) * help_d
       END IF
       idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j+4; val_coo(idx) = 1.0D0
       IF (j > 1) THEN
          idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j-VAR+4; val_coo(idx) = -1.0D0
       END IF
       inh(i) = omega((j-1)/VAR+1) * help_i

       ! delta b_i
       i = i + 1
       help_a = 0.0D0
       help_b = 0.0D0
       help_c = 0.0D0
       help_d = 0.0D0
       help_i = 0.0D0
       DO l = ii, ie
          h_j = x(l) - x(ii)
          x_h = f(x(l),m) * f(x(l),m)
          help_a = help_a + h_j * x_h
          help_b = help_b + h_j * h_j * x_h
          help_c = help_c + h_j * h_j * h_j * x_h
          help_d = help_d + h_j * h_j * h_j * h_j * x_h
          help_i = help_i + h_j * f(x(l),m) * y(l)
       END DO
       IF (ABS(help_a) > 1D-15) THEN
          idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j; val_coo(idx) = omega((j-1)/VAR+1) * help_a
       END IF
       IF (ABS(help_b) > 1D-15) THEN
          idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j+1; val_coo(idx) = omega((j-1)/VAR+1) * help_b
       END IF
       IF (ABS(help_c) > 1D-15) THEN
          idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j+2; val_coo(idx) = omega((j-1)/VAR+1) * help_c
       END IF
       IF (ABS(help_d) > 1D-15) THEN
          idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j+3; val_coo(idx) = omega((j-1)/VAR+1) * help_d
       END IF
       idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j+4; val_coo(idx) = h
       idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j+5; val_coo(idx) = 1.0D0
       IF (j == 1) THEN
          IF (mu1 == 1) THEN
             idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = (len_indx-1)*VAR+4; val_coo(idx) = DBLE(mu1)
          END IF
          IF (mu2 == 1) THEN
             idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = (len_indx-1)*VAR+5; val_coo(idx) = DBLE(mu2)
          END IF
       ELSE
          idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j-VAR+5; val_coo(idx) = -1.0D0
       END IF
       inh(i) = omega((j-1)/VAR+1) * help_i

       ! delta c_i
       i = i + 1
       help_a = 0.0D0
       help_b = 0.0D0
       help_c = 0.0D0
       help_d = 0.0D0
       help_i = 0.0D0
       DO l = ii, ie
          h_j = x(l) - x(ii)
          x_h = f(x(l),m) * f(x(l),m)
          help_a = help_a + h_j * h_j * x_h
          help_b = help_b + h_j * h_j * h_j * x_h
          help_c = help_c + h_j * h_j * h_j * h_j * x_h
          help_d = help_d + h_j * h_j * h_j * h_j * h_j * x_h
          help_i = help_i + h_j * h_j * f(x(l),m) * y(l)
       END DO
       IF (ABS(help_a) > 1D-15) THEN
          idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j; val_coo(idx) = omega((j-1)/VAR+1) * help_a
       END IF
       IF (ABS(help_b) > 1D-15) THEN
          idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j+1; val_coo(idx) = omega((j-1)/VAR+1) * help_b
       END IF
       IF (ABS(help_c) > 1D-15) THEN
          idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j+2; val_coo(idx) = omega((j-1)/VAR+1) * help_c
       END IF
       IF (ABS(help_d) > 1D-15) THEN
          idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j+3; val_coo(idx) = omega((j-1)/VAR+1) * help_d
       END IF
       idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j+4; val_coo(idx) = h * h
       idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j+5; val_coo(idx) = 2.0D0 * h
       idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j+6; val_coo(idx) = 1.0D0
       IF (j == 1) THEN
          IF (nu1 == 1) THEN
             idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = (len_indx-1)*VAR+4; val_coo(idx) = DBLE(nu1)
          END IF
          IF (nu2 == 1) THEN
             idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = (len_indx-1)*VAR+5; val_coo(idx) = DBLE(nu2)
          END IF
       ELSE
          idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j-VAR+6; val_coo(idx) = -1.0D0
       END IF
       inh(i) = omega((j-1)/VAR+1) * help_i

       ! delta DELTA d_i
       i = i + 1
       help_a = 0.0D0
       help_b = 0.0D0
       help_c = 0.0D0
       help_d = 0.0D0
       help_i = 0.0D0
       DO l = ii, ie
          h_j = x(l) - x(ii)
          x_h = f(x(l),m) * f(x(l),m)
          help_a = help_a + h_j * h_j * h_j * x_h
          help_b = help_b + h_j * h_j * h_j * h_j * x_h
          help_c = help_c + h_j * h_j * h_j * h_j * h_j * x_h
          help_d = help_d + h_j * h_j * h_j * h_j * h_j * h_j * x_h
          help_i = help_i + h_j * h_j * h_j * f(x(l),m) * y(l)
       END DO
       IF (ABS(help_a) > 1D-15) THEN
          idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j; val_coo(idx) = omega((j-1)/VAR+1) * help_a
       END IF
       IF (ABS(help_b) > 1D-15) THEN
          idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j+1; val_coo(idx) = omega((j-1)/VAR+1) * help_b
       END IF
       IF (ABS(help_c) > 1D-15) THEN
          idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j+2; val_coo(idx) = omega((j-1)/VAR+1) * help_c
       END IF
       IF (ABS(help_d) > 1D-15 .OR. ABS(lambda((j-1)/VAR+1)) > 1D-15) THEN
          idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j+3; val_coo(idx) = omega((j-1)/VAR+1) * help_d + lambda((j-1)/VAR+1)
       END IF
       idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j+4; val_coo(idx) = h * h * h
       idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j+5; val_coo(idx) = 3.0D0 * h * h
       idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j+6; val_coo(idx) = 3.0D0 * h
       inh(i) = omega((j-1)/VAR+1) * help_i
    END DO

    ! Last segment special conditions
    j = VAR*(len_indx-1)+1
    ii = indx(len_indx)
    ie = SIZE(x)
    
    ! delta a_{N-1}
    i = i + 1
    help_a = 0.0D0
    help_inh = 0.0D0
    DO l = ii, ie
       help_a = help_a + f(x(l),m) * f(x(l),m)
       help_inh = help_inh + f(x(l),m) * y(l)
    END DO
    IF (ABS(help_a) > 1D-15) THEN
       idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = (len_indx-1)*VAR+1; val_coo(idx) = omega(len_indx) * help_a
    END IF
    idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = (len_indx-2)*VAR+5; val_coo(idx) = -1.0D0
    inh(i) = omega(len_indx) * help_inh
    
    ! delta b_{N-1}
    i = i + 1
    idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = (len_indx-2)*VAR+6; val_coo(idx) = -1.0D0
    IF (sig1 == 1) THEN
       idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = (len_indx-1)*VAR+4; val_coo(idx) = DBLE(sig1)
    END IF
    IF (sig2 == 1) THEN
       idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = (len_indx-1)*VAR+5; val_coo(idx) = DBLE(sig2)
    END IF
    
    ! delta c_{N-1}
    i = i + 1
    idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = (len_indx-2)*VAR+7; val_coo(idx) = -1.0D0
    IF (rho1 == 1) THEN
       idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = (len_indx-1)*VAR+4; val_coo(idx) = DBLE(rho1)
    END IF
    IF (rho2 == 1) THEN
       idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = (len_indx-1)*VAR+5; val_coo(idx) = DBLE(rho2)
    END IF
    
    ! Boundary condition 2
    i = i + 1
    ! For sparse matrices, only add non-zero entries
    IF (mu2 /= 0) THEN
       idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = 2; val_coo(idx) = DBLE(mu2)
    END IF
    IF (nu2 /= 0) THEN
       idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = 3; val_coo(idx) = DBLE(nu2)
    END IF
    IF (sig2 /= 0) THEN
       idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = (len_indx-1)*VAR + 2; val_coo(idx) = DBLE(sig2)
    END IF
    IF (rho2 /= 0) THEN
       idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = (len_indx-1)*VAR + 3; val_coo(idx) = DBLE(rho2)
    END IF
    inh(i) = cn

    ! Total non-zeros
    nnz = idx
    
    IF (nnz == 0) THEN
       WRITE(0,*) 'ERROR: No non-zero entries in matrix!'
       STOP
    END IF
    
    IF (nnz > max_nnz) THEN
       WRITE(0,*) 'CRITICAL ERROR: Buffer overflow detected!'
       WRITE(0,*) 'Actual non-zeros:', nnz, ' > estimated max:', max_nnz
       WRITE(0,*) 'This indicates memory corruption has occurred.'
       WRITE(0,*) 'Increase max_nnz estimate in splinecof3_direct_sparse.f90'
       STOP 'Memory safety violation detected'
    END IF

    ! Store COO matrix for inspection
    IF (ALLOCATED(last_irow_coo)) DEALLOCATE(last_irow_coo)
    IF (ALLOCATED(last_icol_coo)) DEALLOCATE(last_icol_coo)
    IF (ALLOCATED(last_val_coo)) DEALLOCATE(last_val_coo)
    IF (ALLOCATED(last_rhs_coo)) DEALLOCATE(last_rhs_coo)
    ALLOCATE(last_irow_coo(nnz), last_icol_coo(nnz), last_val_coo(nnz), &
             last_rhs_coo(size_dimension))
    last_irow_coo(1:nnz) = irow_coo(1:nnz)
    last_icol_coo(1:nnz) = icol_coo(1:nnz)
    last_val_coo(1:nnz) = val_coo(1:nnz)
    last_rhs_coo = inh
    last_nnz = nnz
    last_n = size_dimension

    ! Now convert from COO to CSC format
    ! First count entries per column
    ALLOCATE(col_count(ncol), pcol_csc(ncol+1), stat = i_alloc)
    if(i_alloc /= 0) stop 'Allocation for column counts failed!'
    
    col_count = 0
    DO k = 1, nnz
       IF (icol_coo(k) < 1 .OR. icol_coo(k) > ncol) THEN
          WRITE(*,*) 'ERROR: Invalid column index', icol_coo(k), 'at entry', k
          WRITE(*,*) '  Valid range: 1 to', ncol
          STOP
       END IF
       col_count(icol_coo(k)) = col_count(icol_coo(k)) + 1
    END DO
    
    ! Build column pointer
    pcol_csc(1) = 1
    DO j = 1, ncol
       pcol_csc(j+1) = pcol_csc(j) + col_count(j)
    END DO
    
    ! Allocate CSC arrays
    ALLOCATE(irow_csc(nnz), val_csc(nnz), stat = i_alloc)
    if(i_alloc /= 0) stop 'Allocation for CSC arrays failed!'
    
    ! Reset column count for second pass
    col_count = 0
    
    ! Fill CSC arrays (this sorts by column)
    DO k = 1, nnz
       j = icol_coo(k)
       pos = pcol_csc(j) + col_count(j)
       irow_csc(pos) = irow_coo(k)
       val_csc(pos) = val_coo(k)
       col_count(j) = col_count(j) + 1
    END DO
    
    ! Call sparse_solve with CSC format
    CALL sparse_solve(nrow, ncol, nnz, irow_csc, pcol_csc, val_csc, inh)
    
    ! Extract solution
    DO i = 1, len_indx
       a(i) = inh((i-1)*VAR+1)
       b(i) = inh((i-1)*VAR+2)
       c(i) = inh((i-1)*VAR+3)
       d(i) = inh((i-1)*VAR+4)
    END DO
    
    ! Clean up
    DEALLOCATE(irow_coo, icol_coo, val_coo, irow_csc, pcol_csc, val_csc, &
               col_count, lambda, omega, inh)

  END SUBROUTINE splinecof3_direct_sparse

  !> Get the last computed COO matrix for inspection
  SUBROUTINE splinecof3_direct_sparse_get_coo(irow, icol, val, rhs, nnz, n)
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: irow, icol
    REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: val, rhs
    INTEGER(I4B), INTENT(OUT) :: nnz, n
    
    nnz = last_nnz
    n = last_n
    IF (nnz > 0 .AND. ALLOCATED(last_irow_coo)) THEN
       ALLOCATE(irow(nnz), icol(nnz), val(nnz), rhs(n))
       irow = last_irow_coo
       icol = last_icol_coo
       val = last_val_coo
       rhs = last_rhs_coo
    END IF
  END SUBROUTINE splinecof3_direct_sparse_get_coo

end module splinecof3_direct_sparse_mod

! Wrapper subroutine to match interface expectations
SUBROUTINE splinecof3_direct_sparse_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
     a, b, c, d, m, f)
  use splinecof3_direct_sparse_mod, only: splinecof3_direct_sparse
  use nrtype, only : I4B, DP
  REAL(DP),                   INTENT(INOUT) :: c1, cn
  REAL(DP),     DIMENSION(:), INTENT(IN)    :: x, y, lambda1
  INTEGER(I4B), DIMENSION(:), INTENT(IN)    :: indx
  REAL(DP),     DIMENSION(:), INTENT(OUT)   :: a, b, c, d
  INTEGER(I4B),               INTENT(IN)    :: sw1, sw2
  REAL(DP),                   INTENT(IN)    :: m
  INTERFACE
     FUNCTION f(x,m)
       use nrtype, only : DP
       IMPLICIT NONE
       REAL(DP), INTENT(IN) :: x, m
       REAL(DP)             :: f
     END FUNCTION f
  END INTERFACE

  CALL splinecof3_direct_sparse(x, y, c1, cn, lambda1, indx, sw1, sw2, &
       a, b, c, d, m, f)
END SUBROUTINE splinecof3_direct_sparse_a