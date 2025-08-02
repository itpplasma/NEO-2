!> Direct sparse implementation that builds matrix in COO format and converts to CSC
module splinecof3_direct_sparse_mod
  use nrtype, only : I4B, DP
  use sparse_mod, only: sparse_solve
  use inter_interfaces, only: calc_opt_lambda3
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
  implicit none
  
  private
  public :: splinecof3_direct_sparse, splinecof3_assemble_matrix
  
contains

  !> Add a matrix entry if non-zero (counting mode just increments counter)
  SUBROUTINE add_entry(counting, idx, i, j, val, irow, icol, vals)
    LOGICAL, INTENT(IN) :: counting
    INTEGER(I4B), INTENT(INOUT) :: idx
    INTEGER(I4B), INTENT(IN) :: i, j
    REAL(DP), INTENT(IN) :: val
    INTEGER(I4B), DIMENSION(:), INTENT(INOUT), OPTIONAL :: irow, icol
    REAL(DP), DIMENSION(:), INTENT(INOUT), OPTIONAL :: vals
    
    ! Add entry following original dense implementation behavior
    IF (should_include_element(val)) THEN
       idx = idx + 1
       IF (.NOT. counting) THEN
          irow(idx) = i
          icol(idx) = j
          vals(idx) = val
       END IF
    END IF
  END SUBROUTINE add_entry
  
  !> Check if matrix element should be included (matches original dense implementation)
  LOGICAL FUNCTION should_include_element(val)
    REAL(DP), INTENT(IN) :: val
    ! Match the behavior of full2sparse which excludes exact zeros
    ! This ensures exact numerical compatibility with dense->sparse conversion
    should_include_element = (val /= 0.0_DP)
  END FUNCTION should_include_element
  
  !> Add boundary condition entries
  SUBROUTINE add_boundary_condition_1(counting, idx, i, mu1, nu1, sig1, rho1, &
                                       len_indx, VAR, irow, icol, vals)
    LOGICAL, INTENT(IN) :: counting
    INTEGER(I4B), INTENT(INOUT) :: idx, i
    INTEGER(I4B), INTENT(IN) :: mu1, nu1, sig1, rho1, len_indx, VAR
    INTEGER(I4B), DIMENSION(:), INTENT(INOUT), OPTIONAL :: irow, icol
    REAL(DP), DIMENSION(:), INTENT(INOUT), OPTIONAL :: vals
    
    i = i + 1
    ! Add ALL boundary parameters unconditionally to match original dense implementation
    CALL add_entry(counting, idx, i, 2, DBLE(mu1), irow, icol, vals)
    CALL add_entry(counting, idx, i, 3, DBLE(nu1), irow, icol, vals)
    CALL add_entry(counting, idx, i, (len_indx-1)*VAR + 2, DBLE(sig1), irow, icol, vals)
    CALL add_entry(counting, idx, i, (len_indx-1)*VAR + 3, DBLE(rho1), irow, icol, vals)
  END SUBROUTINE add_boundary_condition_1
  
  !> Add continuity conditions
  SUBROUTINE add_continuity_conditions(counting, idx, i, j, h, VAR, irow, icol, vals)
    LOGICAL, INTENT(IN) :: counting
    INTEGER(I4B), INTENT(INOUT) :: idx, i
    INTEGER(I4B), INTENT(IN) :: j, VAR
    REAL(DP), INTENT(IN) :: h
    INTEGER(I4B), DIMENSION(:), INTENT(INOUT), OPTIONAL :: irow, icol
    REAL(DP), DIMENSION(:), INTENT(INOUT), OPTIONAL :: vals
    
    ! A_i continuity
    i = i + 1
    idx = idx + 1
    IF (.NOT. counting) THEN
       irow(idx) = i; icol(idx) = j; vals(idx) = 1.0D0
    END IF
    idx = idx + 1
    IF (.NOT. counting) THEN
       irow(idx) = i; icol(idx) = j+1; vals(idx) = h
    END IF
    idx = idx + 1
    IF (.NOT. counting) THEN
       irow(idx) = i; icol(idx) = j+2; vals(idx) = h*h
    END IF
    idx = idx + 1
    IF (.NOT. counting) THEN
       irow(idx) = i; icol(idx) = j+3; vals(idx) = h*h*h
    END IF
    idx = idx + 1
    IF (.NOT. counting) THEN
       irow(idx) = i; icol(idx) = j+VAR; vals(idx) = -1.0D0
    END IF
    
    ! B_i continuity
    i = i + 1
    idx = idx + 1
    IF (.NOT. counting) THEN
       irow(idx) = i; icol(idx) = j+1; vals(idx) = 1.0D0
    END IF
    idx = idx + 1
    IF (.NOT. counting) THEN
       irow(idx) = i; icol(idx) = j+2; vals(idx) = 2.0D0*h
    END IF
    idx = idx + 1
    IF (.NOT. counting) THEN
       irow(idx) = i; icol(idx) = j+3; vals(idx) = 3.0D0*h*h
    END IF
    idx = idx + 1
    IF (.NOT. counting) THEN
       irow(idx) = i; icol(idx) = j+VAR+1; vals(idx) = -1.0D0
    END IF
    
    ! C_i continuity
    i = i + 1
    idx = idx + 1
    IF (.NOT. counting) THEN
       irow(idx) = i; icol(idx) = j+2; vals(idx) = 1.0D0
    END IF
    idx = idx + 1
    IF (.NOT. counting) THEN
       irow(idx) = i; icol(idx) = j+3; vals(idx) = 3.0D0*h
    END IF
    idx = idx + 1
    IF (.NOT. counting) THEN
       irow(idx) = i; icol(idx) = j+VAR+2; vals(idx) = -1.0D0
    END IF
  END SUBROUTINE add_continuity_conditions
  
  !> Compute fitting coefficients for an interval
  SUBROUTINE compute_fitting_coeffs(ii, ie, x, y, m, f, help_a, help_b, help_c, help_d, help_i)
    INTEGER(I4B), INTENT(IN) :: ii, ie
    REAL(DP), DIMENSION(:), INTENT(IN) :: x, y
    REAL(DP), INTENT(IN) :: m
    INTERFACE
       FUNCTION f(x,m)
         use nrtype, only : DP
         IMPLICIT NONE
         REAL(DP), INTENT(IN) :: x, m
         REAL(DP)             :: f
       END FUNCTION f
    END INTERFACE
    REAL(DP), INTENT(OUT) :: help_a, help_b, help_c, help_d, help_i
    
    INTEGER(I4B) :: l
    REAL(DP) :: h_j, x_h
    
    help_a = 0.0D0; help_b = 0.0D0; help_c = 0.0D0; help_d = 0.0D0; help_i = 0.0D0
    
    DO l = ii, ie
       h_j = x(l) - x(ii)
       x_h = f(x(l),m) * f(x(l),m)
       help_a = help_a + x_h
       help_b = help_b + h_j * x_h
       help_c = help_c + h_j * h_j * x_h
       help_d = help_d + h_j * h_j * h_j * x_h
       help_i = help_i + f(x(l),m) * y(l)
    END DO
  END SUBROUTINE compute_fitting_coeffs
  
  !> Process one interval's matrix entries
  SUBROUTINE process_interval(counting, idx, i, j, ii, ie, x, y, m, f, omega, lambda, &
                               mu1, mu2, VAR, indx, irow, icol, vals, inh)
    LOGICAL, INTENT(IN) :: counting
    INTEGER(I4B), INTENT(INOUT) :: idx, i
    INTEGER(I4B), INTENT(IN) :: j, ii, ie, mu1, mu2, VAR
    REAL(DP), DIMENSION(:), INTENT(IN) :: x, y, omega, lambda
    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
    REAL(DP), INTENT(IN) :: m
    INTERFACE
       FUNCTION f(x,m)
         use nrtype, only : DP
         IMPLICIT NONE
         REAL(DP), INTENT(IN) :: x, m
         REAL(DP)             :: f
       END FUNCTION f
    END INTERFACE
    INTEGER(I4B), DIMENSION(:), INTENT(INOUT), OPTIONAL :: irow, icol
    REAL(DP), DIMENSION(:), INTENT(INOUT), OPTIONAL :: vals, inh
    
    REAL(DP) :: help_a, help_b, help_c, help_d, help_i, h_j, x_h, h
    INTEGER(I4B) :: interval_idx, l, len_indx
    
    interval_idx = (j-1)/VAR + 1
    len_indx = SIZE(indx)
    h = x(indx((j-1)/VAR+2)) - x(ii)
    
    ! Delta a_i
    CALL compute_fitting_coeffs(ii, ie, x, y, m, f, help_a, help_b, help_c, help_d, help_i)
    i = i + 1
    CALL add_entry(counting, idx, i, j, omega(interval_idx) * help_a, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+1, omega(interval_idx) * help_b, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+2, omega(interval_idx) * help_c, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+3, omega(interval_idx) * help_d, irow, icol, vals)
    idx = idx + 1
    IF (.NOT. counting) THEN
       irow(idx) = i; icol(idx) = j+4; vals(idx) = 1.0D0
    END IF
    IF (j > 1) THEN
       idx = idx + 1
       IF (.NOT. counting) THEN
          irow(idx) = i; icol(idx) = j-VAR+4; vals(idx) = -1.0D0
       END IF
    END IF
    IF (.NOT. counting) inh(i) = omega(interval_idx) * help_i
    
    ! Delta b_i
    i = i + 1
    help_a = 0.0D0; help_b = 0.0D0; help_c = 0.0D0; help_d = 0.0D0; help_i = 0.0D0
    DO l = ii, ie
       h_j = x(l) - x(ii)
       x_h = f(x(l),m) * f(x(l),m)
       help_a = help_a + h_j * x_h
       help_b = help_b + h_j * h_j * x_h
       help_c = help_c + h_j * h_j * h_j * x_h
       help_d = help_d + h_j * h_j * h_j * h_j * x_h
       help_i = help_i + h_j * f(x(l),m) * y(l)
    END DO
    CALL add_entry(counting, idx, i, j, omega(interval_idx) * help_a, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+1, omega(interval_idx) * help_b, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+2, omega(interval_idx) * help_c, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+3, omega(interval_idx) * help_d, irow, icol, vals)
    idx = idx + 1
    IF (.NOT. counting) THEN
       irow(idx) = i; icol(idx) = j+4; vals(idx) = h
    END IF
    idx = idx + 1
    IF (.NOT. counting) THEN
       irow(idx) = i; icol(idx) = j+5; vals(idx) = 1.0D0
    END IF
    IF (j == 1) THEN
       CALL add_entry(counting, idx, i, (len_indx-1)*VAR+4, DBLE(mu1), irow, icol, vals)
       CALL add_entry(counting, idx, i, (len_indx-1)*VAR+5, DBLE(mu2), irow, icol, vals)
    ELSE
       idx = idx + 1
       IF (.NOT. counting) THEN
          irow(idx) = i; icol(idx) = j-VAR+5; vals(idx) = -1.0D0
       END IF
    END IF
    IF (.NOT. counting) inh(i) = omega(interval_idx) * help_i
    
    ! Delta c_i
    i = i + 1
    help_a = 0.0D0; help_b = 0.0D0; help_c = 0.0D0; help_d = 0.0D0; help_i = 0.0D0
    DO l = ii, ie
       h_j = x(l) - x(ii)
       x_h = f(x(l),m) * f(x(l),m)
       help_a = help_a + h_j * h_j * h_j * x_h
       help_b = help_b + h_j * h_j * h_j * h_j * x_h
       help_c = help_c + h_j * h_j * h_j * h_j * h_j * x_h
       help_d = help_d + h_j * h_j * h_j * h_j * h_j * h_j * x_h
       help_i = help_i + h_j * h_j * h_j * f(x(l),m) * y(l)
    END DO
    CALL add_entry(counting, idx, i, j, omega(interval_idx) * help_a, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+1, omega(interval_idx) * help_b, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+2, omega(interval_idx) * help_c, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+3, omega(interval_idx) * help_d + lambda(interval_idx), irow, icol, vals)
    idx = idx + 1
    IF (.NOT. counting) THEN
       irow(idx) = i; icol(idx) = j+4; vals(idx) = h * h * h
    END IF
    idx = idx + 1
    IF (.NOT. counting) THEN
       irow(idx) = i; icol(idx) = j+5; vals(idx) = 3.0D0 * h * h
    END IF
    idx = idx + 1
    IF (.NOT. counting) THEN
       irow(idx) = i; icol(idx) = j+6; vals(idx) = 3.0D0 * h
    END IF
    IF (.NOT. counting) inh(i) = omega(interval_idx) * help_i
    
  END SUBROUTINE process_interval

  !> Process first interval fitting conditions exactly as in dense reference
  SUBROUTINE process_first_interval(counting, idx, i, j, ii, ie, x, y, m, f, omega, lambda, &
                                   mu1, mu2, nu1, nu2, VAR, len_indx, indx, irow, icol, vals, inh)
    LOGICAL, INTENT(IN) :: counting
    INTEGER(I4B), INTENT(INOUT) :: idx, i
    INTEGER(I4B), INTENT(IN) :: j, ii, ie, mu1, mu2, nu1, nu2, VAR, len_indx
    REAL(DP), DIMENSION(:), INTENT(IN) :: x, y, omega, lambda
    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
    REAL(DP), INTENT(IN) :: m
    INTERFACE
       FUNCTION f(x,m)
         use nrtype, only : DP
         IMPLICIT NONE
         REAL(DP), INTENT(IN) :: x, m
         REAL(DP)             :: f
       END FUNCTION f
    END INTERFACE
    INTEGER(I4B), DIMENSION(:), INTENT(INOUT), OPTIONAL :: irow, icol
    REAL(DP), DIMENSION(:), INTENT(INOUT), OPTIONAL :: vals, inh
    
    REAL(DP) :: help_a, help_b, help_c, help_d, help_i, h_j, x_h, h
    INTEGER(I4B) :: interval_idx, l
    
    interval_idx = (j-1)/VAR + 1
    h = x(indx((j-1)/VAR+2)) - x(ii)
    
    ! delta a_i
    i = i + 1
    help_a = 0.0D0; help_b = 0.0D0; help_c = 0.0D0; help_d = 0.0D0; help_i = 0.0D0
    DO l = ii, ie
       h_j = x(l) - x(ii)
       x_h = f(x(l),m) * f(x(l),m)
       help_a = help_a + x_h
       help_b = help_b + h_j * x_h
       help_c = help_c + h_j * h_j * x_h
       help_d = help_d + h_j * h_j * h_j * x_h
       help_i = help_i + f(x(l),m) * y(l)
    END DO
    ! Always add fitting coefficients (even if small) to match dense structure
    idx = idx + 1
    IF (.NOT. counting) THEN
       irow(idx) = i; icol(idx) = j+0; vals(idx) = omega(interval_idx) * help_a
    END IF
    idx = idx + 1
    IF (.NOT. counting) THEN
       irow(idx) = i; icol(idx) = j+1; vals(idx) = omega(interval_idx) * help_b
    END IF
    idx = idx + 1
    IF (.NOT. counting) THEN
       irow(idx) = i; icol(idx) = j+2; vals(idx) = omega(interval_idx) * help_c
    END IF
    idx = idx + 1
    IF (.NOT. counting) THEN
       irow(idx) = i; icol(idx) = j+3; vals(idx) = omega(interval_idx) * help_d
    END IF
    CALL add_entry(counting, idx, i, j+4, 1.0D0, irow, icol, vals)
    IF (.NOT. counting) inh(i) = omega(interval_idx) * help_i
    
    ! delta b_i
    i = i + 1
    help_a = 0.0D0; help_b = 0.0D0; help_c = 0.0D0; help_d = 0.0D0; help_i = 0.0D0
    DO l = ii, ie
       h_j = x(l) - x(ii)
       x_h = f(x(l),m) * f(x(l),m)
       help_a = help_a + h_j * x_h
       help_b = help_b + h_j * h_j * x_h
       help_c = help_c + h_j * h_j * h_j * x_h
       help_d = help_d + h_j * h_j * h_j * h_j * x_h
       help_i = help_i + h_j * f(x(l),m) * y(l)
    END DO
    CALL add_entry(counting, idx, i, j+0, omega(interval_idx) * help_a, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+1, omega(interval_idx) * help_b, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+2, omega(interval_idx) * help_c, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+3, omega(interval_idx) * help_d, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+4, h, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+5, 1.0D0, irow, icol, vals)
    CALL add_entry(counting, idx, i, (len_indx-1)*VAR+4, DBLE(mu1), irow, icol, vals)
    CALL add_entry(counting, idx, i, (len_indx-1)*VAR+5, DBLE(mu2), irow, icol, vals)
    IF (.NOT. counting) inh(i) = omega(interval_idx) * help_i
    
    ! delta c_i
    i = i + 1
    help_a = 0.0D0; help_b = 0.0D0; help_c = 0.0D0; help_d = 0.0D0; help_i = 0.0D0
    DO l = ii, ie
       h_j = x(l) - x(ii)
       x_h = f(x(l),m) * f(x(l),m)
       help_a = help_a + h_j * h_j * x_h
       help_b = help_b + h_j * h_j * h_j * x_h
       help_c = help_c + h_j * h_j * h_j * h_j * x_h
       help_d = help_d + h_j * h_j * h_j * h_j * h_j * x_h
       help_i = help_i + h_j * h_j * f(x(l),m) * y(l)
    END DO
    CALL add_entry(counting, idx, i, j+0, omega(interval_idx) * help_a, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+1, omega(interval_idx) * help_b, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+2, omega(interval_idx) * help_c, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+3, omega(interval_idx) * help_d, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+4, h*h, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+5, 2.0D0*h, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+6, 1.0D0, irow, icol, vals)
    CALL add_entry(counting, idx, i, (len_indx-1)*VAR+4, DBLE(nu1), irow, icol, vals)
    CALL add_entry(counting, idx, i, (len_indx-1)*VAR+5, DBLE(nu2), irow, icol, vals)
    IF (.NOT. counting) inh(i) = omega(interval_idx) * help_i
    
    ! delta DELTA d_i
    i = i + 1
    help_a = 0.0D0; help_b = 0.0D0; help_c = 0.0D0; help_d = 0.0D0; help_i = 0.0D0
    DO l = ii, ie
       h_j = x(l) - x(ii)
       x_h = f(x(l),m) * f(x(l),m)
       help_a = help_a + h_j * h_j * h_j * x_h
       help_b = help_b + h_j * h_j * h_j * h_j * x_h
       help_c = help_c + h_j * h_j * h_j * h_j * h_j * x_h
       help_d = help_d + h_j * h_j * h_j * h_j * h_j * h_j * x_h
       help_i = help_i + h_j * h_j * h_j * f(x(l),m) * y(l)
    END DO
    CALL add_entry(counting, idx, i, j+0, omega(interval_idx) * help_a, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+1, omega(interval_idx) * help_b, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+2, omega(interval_idx) * help_c, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+3, omega(interval_idx) * help_d + lambda(interval_idx), irow, icol, vals)
    CALL add_entry(counting, idx, i, j+4, h*h*h, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+5, 3.0D0*h*h, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+6, 3.0D0*h, irow, icol, vals)
    IF (.NOT. counting) inh(i) = omega(interval_idx) * help_i
    
  END SUBROUTINE process_first_interval

  !> Process middle interval fitting conditions exactly as in dense reference
  SUBROUTINE process_middle_interval(counting, idx, i, j, ii, ie, x, y, m, f, omega, lambda, &
                                    VAR, indx, irow, icol, vals, inh)
    LOGICAL, INTENT(IN) :: counting
    INTEGER(I4B), INTENT(INOUT) :: idx, i
    INTEGER(I4B), INTENT(IN) :: j, ii, ie, VAR
    REAL(DP), DIMENSION(:), INTENT(IN) :: x, y, omega, lambda
    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
    REAL(DP), INTENT(IN) :: m
    INTERFACE
       FUNCTION f(x,m)
         use nrtype, only : DP
         IMPLICIT NONE
         REAL(DP), INTENT(IN) :: x, m
         REAL(DP)             :: f
       END FUNCTION f
    END INTERFACE
    INTEGER(I4B), DIMENSION(:), INTENT(INOUT), OPTIONAL :: irow, icol
    REAL(DP), DIMENSION(:), INTENT(INOUT), OPTIONAL :: vals, inh
    
    REAL(DP) :: help_a, help_b, help_c, help_d, help_i, h_j, x_h, h
    INTEGER(I4B) :: interval_idx, l
    
    interval_idx = (j-1)/VAR + 1
    h = x(indx((j-1)/VAR+2)) - x(ii)
    
    ! delta a_i
    i = i + 1
    help_a = 0.0D0; help_b = 0.0D0; help_c = 0.0D0; help_d = 0.0D0; help_i = 0.0D0
    DO l = ii, ie
       h_j = x(l) - x(ii)
       x_h = f(x(l),m) * f(x(l),m)
       help_a = help_a + x_h
       help_b = help_b + h_j * x_h
       help_c = help_c + h_j * h_j * x_h
       help_d = help_d + h_j * h_j * h_j * x_h
       help_i = help_i + f(x(l),m) * y(l)
    END DO
    CALL add_entry(counting, idx, i, j+0, omega(interval_idx) * help_a, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+1, omega(interval_idx) * help_b, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+2, omega(interval_idx) * help_c, irow, icol, vals)  
    CALL add_entry(counting, idx, i, j+3, omega(interval_idx) * help_d, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+4, 1.0D0, irow, icol, vals)
    CALL add_entry(counting, idx, i, j-VAR+4, -1.0D0, irow, icol, vals)
    IF (.NOT. counting) inh(i) = omega(interval_idx) * help_i
    
    ! delta b_i
    i = i + 1
    help_a = 0.0D0; help_b = 0.0D0; help_c = 0.0D0; help_d = 0.0D0; help_i = 0.0D0
    DO l = ii, ie
       h_j = x(l) - x(ii)
       x_h = f(x(l),m) * f(x(l),m)
       help_a = help_a + h_j * x_h
       help_b = help_b + h_j * h_j * x_h
       help_c = help_c + h_j * h_j * h_j * x_h
       help_d = help_d + h_j * h_j * h_j * h_j * x_h
       help_i = help_i + h_j * f(x(l),m) * y(l)
    END DO
    CALL add_entry(counting, idx, i, j+0, omega(interval_idx) * help_a, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+1, omega(interval_idx) * help_b, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+2, omega(interval_idx) * help_c, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+3, omega(interval_idx) * help_d, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+4, h, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+5, 1.0D0, irow, icol, vals)
    CALL add_entry(counting, idx, i, j-VAR+5, -1.0D0, irow, icol, vals)
    IF (.NOT. counting) inh(i) = omega(interval_idx) * help_i
    
    ! delta c_i
    i = i + 1
    help_a = 0.0D0; help_b = 0.0D0; help_c = 0.0D0; help_d = 0.0D0; help_i = 0.0D0
    DO l = ii, ie
       h_j = x(l) - x(ii)
       x_h = f(x(l),m) * f(x(l),m)
       help_a = help_a + h_j * h_j * x_h
       help_b = help_b + h_j * h_j * h_j * x_h
       help_c = help_c + h_j * h_j * h_j * h_j * x_h
       help_d = help_d + h_j * h_j * h_j * h_j * h_j * x_h
       help_i = help_i + h_j * h_j * f(x(l),m) * y(l)
    END DO
    CALL add_entry(counting, idx, i, j+0, omega(interval_idx) * help_a, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+1, omega(interval_idx) * help_b, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+2, omega(interval_idx) * help_c, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+3, omega(interval_idx) * help_d, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+4, h*h, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+5, 2.0D0*h, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+6, 1.0D0, irow, icol, vals)
    CALL add_entry(counting, idx, i, j-VAR+6, -1.0D0, irow, icol, vals)
    IF (.NOT. counting) inh(i) = omega(interval_idx) * help_i
    
    ! delta DELTA d_i
    i = i + 1
    help_a = 0.0D0; help_b = 0.0D0; help_c = 0.0D0; help_d = 0.0D0; help_i = 0.0D0
    DO l = ii, ie
       h_j = x(l) - x(ii)
       x_h = f(x(l),m) * f(x(l),m)
       help_a = help_a + h_j * h_j * h_j * x_h
       help_b = help_b + h_j * h_j * h_j * h_j * x_h
       help_c = help_c + h_j * h_j * h_j * h_j * h_j * x_h
       help_d = help_d + h_j * h_j * h_j * h_j * h_j * h_j * x_h
       help_i = help_i + h_j * h_j * h_j * f(x(l),m) * y(l)
    END DO
    CALL add_entry(counting, idx, i, j+0, omega(interval_idx) * help_a, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+1, omega(interval_idx) * help_b, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+2, omega(interval_idx) * help_c, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+3, omega(interval_idx) * help_d + lambda(interval_idx), irow, icol, vals)
    CALL add_entry(counting, idx, i, j+4, h*h*h, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+5, 3.0D0*h*h, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+6, 3.0D0*h, irow, icol, vals)
    IF (.NOT. counting) inh(i) = omega(interval_idx) * help_i
    
  END SUBROUTINE process_middle_interval

  !> Build matrix in two passes: count non-zeros, then fill
  SUBROUTINE build_matrix_two_pass(counting, idx, i, x, y, m, f, lambda, omega, &
                                   indx, mu1, mu2, nu1, nu2, sig1, sig2, rho1, rho2, &
                                   c1, cn, VAR, len_indx, irow, icol, vals, inh)
    LOGICAL, INTENT(IN) :: counting
    INTEGER(I4B), INTENT(INOUT) :: idx, i
    REAL(DP), DIMENSION(:), INTENT(IN) :: x, y, lambda, omega
    REAL(DP), INTENT(IN) :: m, c1, cn
    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
    INTEGER(I4B), INTENT(IN) :: mu1, mu2, nu1, nu2, sig1, sig2, rho1, rho2
    INTEGER(I4B), INTENT(IN) :: VAR, len_indx
    INTERFACE
       FUNCTION f(x,m)
         use nrtype, only : DP
         IMPLICIT NONE
         REAL(DP), INTENT(IN) :: x, m
         REAL(DP)             :: f
       END FUNCTION f
    END INTERFACE
    INTEGER(I4B), DIMENSION(:), INTENT(INOUT), OPTIONAL :: irow, icol
    REAL(DP), DIMENSION(:), INTENT(INOUT), OPTIONAL :: vals, inh
    
    INTEGER(I4B) :: j, ii, ie, l
    REAL(DP) :: h, h_j, x_h, help_a, help_b, help_c, help_d, help_i, help_inh
    
    ! Initialize
    idx = 0
    i = 0
    
    ! Boundary condition 1 - Always add these entries (even if zero) to match dense structure
    i = i + 1
    idx = idx + 1
    IF (.NOT. counting) THEN
       irow(idx) = i; icol(idx) = 2; vals(idx) = DBLE(mu1)
    END IF
    idx = idx + 1
    IF (.NOT. counting) THEN
       irow(idx) = i; icol(idx) = 3; vals(idx) = DBLE(nu1)
    END IF
    idx = idx + 1
    IF (.NOT. counting) THEN
       irow(idx) = i; icol(idx) = (len_indx-1)*VAR + 2; vals(idx) = DBLE(sig1)
    END IF
    idx = idx + 1
    IF (.NOT. counting) THEN
       irow(idx) = i; icol(idx) = (len_indx-1)*VAR + 3; vals(idx) = DBLE(rho1)
    END IF
    IF (.NOT. counting) inh(i) = c1
    
    ! Coefs for first point
    j = 1
    ii = indx((j-1)/VAR+1)
    ie = indx((j-1)/VAR+2) - 1
    h = x(indx((j-1)/VAR+2)) - x(ii)
    
    ! A_i
    i = i + 1
    CALL add_entry(counting, idx, i, j+0, 1.0D0, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+1, h, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+2, h*h, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+3, h*h*h, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+VAR+0, -1.0D0, irow, icol, vals)
    
    ! B_i
    i = i + 1
    CALL add_entry(counting, idx, i, j+1, 1.0D0, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+2, 2.0D0*h, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+3, 3.0D0*h*h, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+VAR+1, -1.0D0, irow, icol, vals)
    
    ! C_i
    i = i + 1
    CALL add_entry(counting, idx, i, j+2, 1.0D0, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+3, 3.0D0*h, irow, icol, vals)
    CALL add_entry(counting, idx, i, j+VAR+2, -1.0D0, irow, icol, vals)
    
    ! delta a_i, b_i, c_i, d_i for first interval - exactly as in dense reference
    CALL process_first_interval(counting, idx, i, j, ii, ie, x, y, m, f, omega, lambda, &
                               mu1, mu2, nu1, nu2, VAR, len_indx, indx, irow, icol, vals, inh)
    
    ! Coefs for points 2 to len_indx-1 - exactly matching dense loop structure
    DO j = VAR+1, VAR*(len_indx-1)-1, VAR
       ii = indx((j-1)/VAR+1)
       ie = indx((j-1)/VAR+2) - 1
       h = x(indx((j-1)/VAR+2)) - x(ii)
       
       ! A_i
       i = i + 1
       CALL add_entry(counting, idx, i, j+0, 1.0D0, irow, icol, vals)
       CALL add_entry(counting, idx, i, j+1, h, irow, icol, vals)
       CALL add_entry(counting, idx, i, j+2, h*h, irow, icol, vals)
       CALL add_entry(counting, idx, i, j+3, h*h*h, irow, icol, vals)
       CALL add_entry(counting, idx, i, j+VAR+0, -1.0D0, irow, icol, vals)
       
       ! B_i
       i = i + 1
       CALL add_entry(counting, idx, i, j+1, 1.0D0, irow, icol, vals)
       CALL add_entry(counting, idx, i, j+2, 2.0D0*h, irow, icol, vals)
       CALL add_entry(counting, idx, i, j+3, 3.0D0*h*h, irow, icol, vals)
       CALL add_entry(counting, idx, i, j+VAR+1, -1.0D0, irow, icol, vals)
       
       ! C_i  
       i = i + 1
       CALL add_entry(counting, idx, i, j+2, 1.0D0, irow, icol, vals)
       CALL add_entry(counting, idx, i, j+3, 3.0D0*h, irow, icol, vals)
       CALL add_entry(counting, idx, i, j+VAR+2, -1.0D0, irow, icol, vals)
       
       ! delta a_i, b_i, c_i, d_i for middle intervals
       CALL process_middle_interval(counting, idx, i, j, ii, ie, x, y, m, f, omega, lambda, &
                                   VAR, indx, irow, icol, vals, inh)
    END DO
    
    ! Last point - exactly as in dense reference
    ii = indx(len_indx)
    ie = ii
    help_a = 0.0D0
    help_inh = 0.0D0
    l = ii
    help_a = help_a + f(x(l),m) * f(x(l),m)
    help_inh = help_inh + f(x(l),m) * y(l)
    
    ! delta a_i
    i = i + 1
    CALL add_entry(counting, idx, i, (len_indx-1)*VAR+1, omega(len_indx) * help_a, irow, icol, vals)
    CALL add_entry(counting, idx, i, (len_indx-2)*VAR+5, omega(len_indx) * (-1.0D0), irow, icol, vals)
    IF (.NOT. counting) inh(i) = omega(len_indx) * help_inh
    
    ! delta b_i
    i = i + 1
    CALL add_entry(counting, idx, i, (len_indx-2)*VAR+6, -1.0D0, irow, icol, vals)
    CALL add_entry(counting, idx, i, (len_indx-1)*VAR+4, DBLE(sig1), irow, icol, vals)
    CALL add_entry(counting, idx, i, (len_indx-1)*VAR+5, DBLE(sig2), irow, icol, vals)
    
    ! delta c_i
    i = i + 1
    CALL add_entry(counting, idx, i, (len_indx-2)*VAR+7, -1.0D0, irow, icol, vals)
    CALL add_entry(counting, idx, i, (len_indx-1)*VAR+4, DBLE(rho1), irow, icol, vals)
    CALL add_entry(counting, idx, i, (len_indx-1)*VAR+5, DBLE(rho2), irow, icol, vals)
    
    ! Boundary condition 2 - use add_entry to handle zero exclusion consistently
    i = i + 1
    CALL add_entry(counting, idx, i, 2, DBLE(mu2), irow, icol, vals)
    CALL add_entry(counting, idx, i, 3, DBLE(nu2), irow, icol, vals)
    CALL add_entry(counting, idx, i, (len_indx-1)*VAR + 2, DBLE(sig2), irow, icol, vals)
    CALL add_entry(counting, idx, i, (len_indx-1)*VAR + 3, DBLE(rho2), irow, icol, vals)
    IF (.NOT. counting) inh(i) = cn
    
  END SUBROUTINE build_matrix_two_pass

  !> Build matrix using original proven approach (single pass)
  SUBROUTINE build_matrix_original(idx, i, x, y, m, f, lambda, omega, &
                                   indx, mu1, mu2, nu1, nu2, sig1, sig2, rho1, rho2, &
                                   c1, cn, VAR, len_indx, irow_coo, icol_coo, val_coo, inh)
    INTEGER(I4B), INTENT(INOUT) :: idx, i
    REAL(DP), DIMENSION(:), INTENT(IN) :: x, y, lambda, omega
    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
    REAL(DP), INTENT(IN) :: m, c1, cn
    INTEGER(I4B), INTENT(IN) :: mu1, mu2, nu1, nu2, sig1, sig2, rho1, rho2, VAR, len_indx
    INTERFACE
       FUNCTION f(x,m)
         use nrtype, only : DP
         IMPLICIT NONE
         REAL(DP), INTENT(IN) :: x, m
         REAL(DP)             :: f
       END FUNCTION f
    END INTERFACE
    INTEGER(I4B), DIMENSION(:), INTENT(INOUT) :: irow_coo, icol_coo
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: val_coo, inh

    INTEGER(I4B) :: j, ii, ie, l
    REAL(DP) :: help_a, help_b, help_c, help_d, help_i, help_inh, h, h_j, x_h

    ! Boundary condition 1
    i = i + 1
    ! Add ALL boundary parameters unconditionally to match original dense implementation
    ! The should_include_element check will handle zero exclusion
    IF (should_include_element(DBLE(mu1))) THEN
       idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = 2; val_coo(idx) = DBLE(mu1)
    END IF
    IF (should_include_element(DBLE(nu1))) THEN
       idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = 3; val_coo(idx) = DBLE(nu1)
    END IF
    IF (should_include_element(DBLE(sig1))) THEN
       idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = (len_indx-1)*VAR + 2; val_coo(idx) = DBLE(sig1)
    END IF
    IF (should_include_element(DBLE(rho1))) THEN
       idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = (len_indx-1)*VAR + 3; val_coo(idx) = DBLE(rho1)
    END IF
    inh(i) = c1

    ! Main loop for each interval
    DO j = 1, VAR*(len_indx-1), VAR
       ii = indx((j-1)/VAR+1)
       ie = indx((j-1)/VAR+2) - 1
       h = x(indx((j-1)/VAR+2)) - x(ii)

       ! delta a_i
       i = i + 1
       help_a = 0.0D0
       help_b = 0.0D0
       help_c = 0.0D0
       help_d = 0.0D0
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
       ! Add fitting coefficients - matches original dense implementation behavior
       IF (should_include_element(omega((j-1)/VAR+1) * help_a)) THEN
          idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j; val_coo(idx) = omega((j-1)/VAR+1) * help_a
       END IF
       IF (should_include_element(omega((j-1)/VAR+1) * help_b)) THEN
          idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j+1; val_coo(idx) = omega((j-1)/VAR+1) * help_b
       END IF
       IF (should_include_element(omega((j-1)/VAR+1) * help_c)) THEN
          idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j+2; val_coo(idx) = omega((j-1)/VAR+1) * help_c
       END IF
       IF (should_include_element(omega((j-1)/VAR+1) * help_d)) THEN
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
       ! Add fitting coefficients - matches original dense implementation behavior
       IF (should_include_element(omega((j-1)/VAR+1) * help_a)) THEN
          idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j; val_coo(idx) = omega((j-1)/VAR+1) * help_a
       END IF
       IF (should_include_element(omega((j-1)/VAR+1) * help_b)) THEN
          idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j+1; val_coo(idx) = omega((j-1)/VAR+1) * help_b
       END IF
       IF (should_include_element(omega((j-1)/VAR+1) * help_c)) THEN
          idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j+2; val_coo(idx) = omega((j-1)/VAR+1) * help_c
       END IF
       IF (should_include_element(omega((j-1)/VAR+1) * help_d)) THEN
          idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j+3; val_coo(idx) = omega((j-1)/VAR+1) * help_d
       END IF
       idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j+4; val_coo(idx) = h
       idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j+5; val_coo(idx) = 1.0D0
       IF (j == 1) THEN
          IF (should_include_element(DBLE(nu1))) THEN
             idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = (len_indx-1)*VAR+4; val_coo(idx) = DBLE(nu1)
          END IF
          IF (should_include_element(DBLE(nu2))) THEN
             idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = (len_indx-1)*VAR+5; val_coo(idx) = DBLE(nu2)
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
          help_a = help_a + h_j * h_j * h_j * x_h
          help_b = help_b + h_j * h_j * h_j * h_j * x_h
          help_c = help_c + h_j * h_j * h_j * h_j * h_j * x_h
          help_d = help_d + h_j * h_j * h_j * h_j * h_j * h_j * x_h
          help_i = help_i + h_j * h_j * h_j * f(x(l),m) * y(l)
       END DO
       ! Add fitting coefficients - use small threshold to avoid numerical issues
       IF (ABS(omega((j-1)/VAR+1) * help_a) > 1.0D-15) THEN
          idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j; val_coo(idx) = omega((j-1)/VAR+1) * help_a
       END IF
       IF (ABS(omega((j-1)/VAR+1) * help_b) > 1.0D-15) THEN
          idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j+1; val_coo(idx) = omega((j-1)/VAR+1) * help_b
       END IF
       IF (ABS(omega((j-1)/VAR+1) * help_c) > 1.0D-15) THEN
          idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = j+2; val_coo(idx) = omega((j-1)/VAR+1) * help_c
       END IF
       IF (should_include_element(omega((j-1)/VAR+1) * help_d + lambda((j-1)/VAR+1))) THEN
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
    ie = ii  ! Last point only, matching original algorithm
    
    ! delta a_{N-1}
    i = i + 1
    help_a = 0.0D0
    help_inh = 0.0D0
    l = ii
    help_a = help_a + f(x(l),m) * f(x(l),m)
    help_inh = help_inh + f(x(l),m) * y(l)
    IF (should_include_element(omega(len_indx) * help_a)) THEN
       idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = (len_indx-1)*VAR+1; val_coo(idx) = omega(len_indx) * help_a
    END IF
    idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = (len_indx-2)*VAR+5; val_coo(idx) = -1.0D0
    inh(i) = omega(len_indx) * help_inh
    
    ! delta b_{N-1}
    i = i + 1
    idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = (len_indx-2)*VAR+6; val_coo(idx) = -1.0D0
    IF (should_include_element(DBLE(sig1))) THEN
       idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = (len_indx-1)*VAR+4; val_coo(idx) = DBLE(sig1)
    END IF
    IF (should_include_element(DBLE(sig2))) THEN
       idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = (len_indx-1)*VAR+5; val_coo(idx) = DBLE(sig2)
    END IF
    
    ! delta c_{N-1}
    i = i + 1
    idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = (len_indx-2)*VAR+7; val_coo(idx) = -1.0D0
    IF (should_include_element(DBLE(rho1))) THEN
       idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = (len_indx-1)*VAR+4; val_coo(idx) = DBLE(rho1)
    END IF
    IF (should_include_element(DBLE(rho2))) THEN
       idx = idx + 1; irow_coo(idx) = i; icol_coo(idx) = (len_indx-1)*VAR+5; val_coo(idx) = DBLE(rho2)
    END IF
    
    ! Boundary condition 2
    i = i + 1
    ! Only add non-zero boundary condition entries
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

  END SUBROUTINE build_matrix_original

  !> Direct sparse implementation matching splinecof3_a algorithm
  !>
  !> IMPORTANT NOTE ON BOUNDARY CONDITIONS:
  !> For clamped end conditions (sw2=3), this implementation has a known limitation:
  !> - The constraint should enforce S'(x_n) = cn (derivative at last data point)
  !> - Instead, it sets b(n-1) = cn, where b(n-1) represents S'(x_{n-1})
  !> - This is mathematically incorrect but maintains compatibility with all other
  !>   implementations in NEO-2 (original dense, fast path)
  !> - The sparse matrix construction naturally produces this behavior, matching
  !>   the original dense implementation exactly
  !> - The spline will NOT have the correct derivative at x_n, but this appears
  !>   sufficient for NEO-2's practical applications
  !>
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
    INTEGER(I4B) :: i, j, k, nnz, idx, nnz_max, ii, ie, l, neq
    INTEGER(I4B) :: i_alloc, mu1, mu2, nu1, nu2, sig1, sig2, rho1, rho2
    INTEGER(I4B) :: nrow, ncol, pos, len_x
    REAL(DP) :: help_a, help_b, help_c, help_d, help_i, h, h_j, x_h
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
    logical :: consecutive_indices

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

    ! Calculate system size exactly as in dense reference implementation
    size_dimension = VAR * len_indx - 2
    neq = size_dimension
    ncol = size_dimension
    nrow = size_dimension
    
    ! Use two-pass approach: first count exact non-zeros, then allocate and fill
    idx = 0
    i = 0
    CALL build_matrix_two_pass(.TRUE., idx, i, x, y, m, f, lambda, omega, &
                              indx, mu1, mu2, nu1, nu2, sig1, sig2, rho1, rho2, &
                              c1, cn, VAR, len_indx)
    nnz = idx
    
    ! Allocate with exact count (no waste)
    ALLOCATE(irow_coo(nnz), icol_coo(nnz), val_coo(nnz), stat=i_alloc, errmsg=error_message)
    if(i_alloc /= 0) then
      write(*,'(A,I0)') 'SPLINECOF3_DIRECT_SPARSE: COO allocation failed (error code: ', i_alloc, ')'
      write(*,'(A)') 'Error message: ' // trim(error_message)
      write(*,'(A,I0)') 'Attempted to allocate arrays of size nnz=', nnz
      write(*,'(A,F0.2,A)') 'Estimated memory required: ', real(nnz*3)*8.0/1024.0/1024.0, ' MB'
      error stop 'SPLINECOF3_DIRECT_SPARSE: Memory allocation failure for COO arrays'
    end if

    ! Second pass: fill the arrays
    idx = 0
    i = 0
    CALL build_matrix_two_pass(.FALSE., idx, i, x, y, m, f, lambda, omega, &
                              indx, mu1, mu2, nu1, nu2, sig1, sig2, rho1, rho2, &
                              c1, cn, VAR, len_indx, irow_coo, icol_coo, val_coo, inh)
    nnz = idx

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
    ! NOTE: Spline matrices are often ill-conditioned. UMFPACK (method 3) is recommended
    ! over BiCGSTAB (method 4) for better stability and accuracy.
    CALL sparse_solve(nrow, ncol, nnz, irow_csc, pcol_csc, val_csc, inh)
    
    ! Extract solution and check for NaN/Inf
    DO i = 1, len_indx
       a(i) = inh((i-1)*VAR+1)
       b(i) = inh((i-1)*VAR+2)
       c(i) = inh((i-1)*VAR+3)
       d(i) = inh((i-1)*VAR+4)
       
       ! Check for NaN or Inf in solution using IEEE intrinsics
       IF (.NOT. ieee_is_finite(a(i)) .OR. .NOT. ieee_is_finite(b(i)) .OR. &
           .NOT. ieee_is_finite(c(i)) .OR. .NOT. ieee_is_finite(d(i))) THEN
          WRITE(*,'(A,I0)') 'ERROR: Non-finite values in spline coefficients at interval ', i
          WRITE(*,'(A,4ES15.6)') '  Spline coefficients [a,b,c,d]: ', a(i), b(i), c(i), d(i)
          IF (ieee_is_nan(a(i)) .OR. ieee_is_nan(b(i)) .OR. ieee_is_nan(c(i)) .OR. ieee_is_nan(d(i))) THEN
            WRITE(*,*) '  NaN detected - likely caused by:'
            WRITE(*,*) '    - Singular or ill-conditioned matrix'
            WRITE(*,*) '    - Invalid boundary conditions or lambda weights'
            WRITE(*,*) '    - Duplicate or improperly ordered x values'
          ELSE
            WRITE(*,*) '  Infinite values detected - likely caused by:'
            WRITE(*,*) '    - Numerical overflow in matrix construction'
            WRITE(*,*) '    - Extreme values in input data or boundary conditions'
          END IF
          WRITE(*,'(A,I0,A,I0)') '  Problem size: len_x=', len_x, ', len_indx=', len_indx
          WRITE(*,'(A,2ES15.6)') '  Boundary conditions c1, cn: ', c1, cn
          WRITE(*,'(A,I0,A,I0)') '  Boundary condition types sw1, sw2: ', sw1, ', sw2: ', sw2
          ERROR STOP 'SPLINECOF3_DIRECT_SPARSE: Non-finite spline coefficients'
       END IF
    END DO
    
    
    ! Follow spline_cof convention: set n-th element to zero
    ! Note: arrays are size len_indx, not len_x when indx is a subset
    IF (len_indx == len_x) THEN
       a(len_x) = 0.0_DP
       b(len_x) = 0.0_DP
       c(len_x) = 0.0_DP
       d(len_x) = 0.0_DP
    END IF
    
    ! Clean up
    DEALLOCATE(irow_coo, icol_coo, val_coo, irow_csc, pcol_csc, val_csc, &
               col_count, lambda, omega, inh)

  END SUBROUTINE splinecof3_direct_sparse

  !> Extract matrix assembly logic from splinecof3_direct_sparse
  !> Returns the assembled COO matrix and RHS vector without solving
  SUBROUTINE splinecof3_assemble_matrix(x, y, c1, cn, lambda1, indx, sw1, sw2, &
       m, f, nrow, ncol, nnz, irow_coo, icol_coo, val_coo, rhs)
    REAL(DP),                   INTENT(INOUT) :: c1, cn
    REAL(DP),     DIMENSION(:), INTENT(IN)    :: x, y, lambda1
    INTEGER(I4B), DIMENSION(:), INTENT(IN)    :: indx
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
    
    ! Output: COO matrix and RHS
    INTEGER(I4B), INTENT(OUT) :: nrow, ncol, nnz
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: irow_coo, icol_coo
    REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: val_coo, rhs
    
    ! Local variables (copied from splinecof3_direct_sparse)
    INTEGER(I4B) :: len_indx, VAR, size_dimension
    INTEGER(I4B) :: i, idx, i_alloc
    INTEGER(I4B) :: mu1, mu2, nu1, nu2, sig1, sig2, rho1, rho2
    INTEGER(I4B) :: len_x
    REAL(DP), DIMENSION(:), ALLOCATABLE :: lambda, omega, inh
    character(200) :: error_message

    ! Initialize variables (copied from splinecof3_direct_sparse)
    VAR = 7
    len_x = SIZE(x)
    len_indx = SIZE(indx)
    size_dimension = VAR * len_indx - 2
    nrow = size_dimension
    ncol = size_dimension

    ! Validation checks (copied from splinecof3_direct_sparse)
    if ( .NOT. ( size(x) == size(y) ) ) then
      write (*,*) 'splinecof3_assemble_matrix: assertion 1 failed'
      stop 'program terminated'
    end if
    if ( .NOT. ( size(indx) == size(lambda1) ) ) then
      write (*,*) 'splinecof3_assemble_matrix: assertion 2 failed'
      stop 'program terminated'
    end if

    do i = 1, len_x-1
      if (x(i) >= x(i+1)) then
        print *, 'SPLINECOF3_ASSEMBLE_MATRIX: error i, x(i), x(i+1)', &
             i, x(i), x(i+1)
        stop 'SPLINECOF3_ASSEMBLE_MATRIX: error  wrong order of x(i)'
      end if
    end do
    do i = 1, len_indx-1
      if (indx(i) < 1) then
        print *, 'SPLINECOF3_ASSEMBLE_MATRIX: error i, indx(i)', i, indx(i)
        stop 'SPLINECOF3_ASSEMBLE_MATRIX: error  indx(i) < 1'
      end if
      if (indx(i) >= indx(i+1)) then
        print *, 'SPLINECOF3_ASSEMBLE_MATRIX: error i, indx(i), indx(i+1)', &
              i, indx(i), indx(i+1)
        stop 'SPLINECOF3_ASSEMBLE_MATRIX: error  wrong order of indx(i)'
      end if
      if (indx(i) > len_x) then
        print *, 'SPLINECOF3_ASSEMBLE_MATRIX: error i, indx(i), indx(i+1)', &
              i, indx(i), indx(i+1)
        stop 'SPLINECOF3_ASSEMBLE_MATRIX: error  indx(i) > len_x'
      end if
    end do
    if (indx(len_indx) < 1) then
      print *, 'SPLINECOF3_ASSEMBLE_MATRIX: error len_indx, indx(len_indx)', &
            len_indx, indx(len_indx)
      stop 'SPLINECOF3_ASSEMBLE_MATRIX: error  indx(max) < 1'
    end if
    if (indx(len_indx) > len_x) then
      print *, 'SPLINECOF3_ASSEMBLE_MATRIX: error len_indx, indx(len_indx)', &
            len_indx, indx(len_indx)
      stop 'SPLINECOF3_ASSEMBLE_MATRIX: error  indx(max) > len_x'
    end if

    if (sw1 == sw2) then
      stop 'SPLINECOF3_ASSEMBLE_MATRIX: error  two identical boundary conditions'
    end if

    ! Allocate work arrays (copied from splinecof3_direct_sparse)
    ALLOCATE(lambda(len_indx), omega(len_indx), inh(size_dimension), &
             stat = i_alloc, errmsg=error_message)
    if(i_alloc /= 0) then
      write(*,*) 'splinecof3_assemble_matrix: Allocation failed:', trim(error_message)
      stop
    end if

    ! Process boundary conditions (copied from splinecof3_direct_sparse)
    IF (DABS(c1) > 1.0E30) THEN
      c1 = 0.0D0
    END IF
    IF (DABS(cn) > 1.0E30) THEN
      cn = 0.0D0
    END IF

    ! Calculate optimal weights for smoothing (copied from splinecof3_direct_sparse)
    IF ( MAXVAL(lambda1) < 0.0D0 ) THEN
      CALL calc_opt_lambda3(x, y, omega)
    ELSE
      omega  = lambda1
    END IF
    lambda = 1.0D0 - omega
    
    ! Initialize RHS vector
    inh = 0.0D0

    ! Set boundary condition switches (copied from splinecof3_direct_sparse)
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
    
    ! Use two-pass approach to count exact non-zeros, then allocate and fill
    idx = 0
    i = 0
    CALL build_matrix_two_pass(.TRUE., idx, i, x, y, m, f, lambda, omega, &
                              indx, mu1, mu2, nu1, nu2, sig1, sig2, rho1, rho2, &
                              c1, cn, VAR, len_indx)
    nnz = idx
    
    ! Allocate COO arrays with exact count
    ALLOCATE(irow_coo(nnz), icol_coo(nnz), val_coo(nnz), rhs(size_dimension), &
             stat=i_alloc, errmsg=error_message)
    if(i_alloc /= 0) then
      write(*,'(A,I0)') 'SPLINECOF3_ASSEMBLE_MATRIX: COO allocation failed (error code: ', i_alloc, ')'
      write(*,'(A)') 'Error message: ' // trim(error_message)
      write(*,'(A,I0)') 'Attempted to allocate arrays of size nnz=', nnz
      error stop 'SPLINECOF3_ASSEMBLE_MATRIX: Memory allocation failure for COO arrays'
    end if

    ! Second pass: fill the arrays  
    idx = 0
    i = 0
    CALL build_matrix_two_pass(.FALSE., idx, i, x, y, m, f, lambda, omega, &
                              indx, mu1, mu2, nu1, nu2, sig1, sig2, rho1, rho2, &
                              c1, cn, VAR, len_indx, irow_coo, icol_coo, val_coo, inh)
    nnz = idx
    
    ! Copy RHS to output
    rhs = inh
    
    ! Clean up work arrays
    DEALLOCATE(lambda, omega, inh)
    
  END SUBROUTINE splinecof3_assemble_matrix

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