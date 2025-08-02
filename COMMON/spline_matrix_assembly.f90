!> Module for spline matrix assembly routines
!> Separates matrix construction from solving for testability
module spline_matrix_assembly_mod
  use nrtype, only : I4B, DP
  use inter_interfaces, only: calc_opt_lambda3
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  
  private
  public :: assemble_spline_matrix_sparse_coo
  public :: assemble_spline_matrix_fast_tridiagonal
  public :: compare_sparse_matrices
  public :: extract_tridiagonal_from_sparse
  
contains

  !> Assemble sparse matrix for spline system (returns COO format)
  !> This extracts the matrix assembly logic from splinecof3_direct_sparse
  subroutine assemble_spline_matrix_sparse_coo(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                                               m, f, nrow, ncol, nnz, irow_coo, icol_coo, &
                                               val_coo, rhs)
    real(DP), dimension(:), intent(in) :: x, y, lambda1
    real(DP), intent(inout) :: c1, cn
    real(DP), intent(in) :: m
    integer(I4B), dimension(:), intent(in) :: indx
    integer(I4B), intent(in) :: sw1, sw2
    interface
       function f(x,m)
         use nrtype, only : DP
         implicit none
         real(DP), intent(in) :: x, m
         real(DP) :: f
       end function f
    end interface
    integer(I4B), intent(out) :: nrow, ncol, nnz
    integer(I4B), allocatable, dimension(:), intent(out) :: irow_coo, icol_coo
    real(DP), allocatable, dimension(:), intent(out) :: val_coo, rhs
    
    ! Local variables matching original implementation
    integer(I4B), parameter :: VAR = 7
    integer(I4B) :: len_indx, len_x, i_alloc, idx, i
    real(DP), allocatable :: omega(:), lambda(:)
    integer(I4B) :: mu1, mu2, nu1, nu2, sig1, sig2, rho1, rho2
    character(200) :: error_message
    
    ! Get dimensions
    len_x = size(x)
    len_indx = size(indx)
    nrow = VAR * len_indx - 2
    ncol = nrow
    
    ! Process boundary conditions (matching original)
    if (dabs(c1) > 1.0E30) then
      c1 = 0.0D0
    end if
    if (dabs(cn) > 1.0E30) then
      cn = 0.0D0
    end if
    
    ! Allocate workspace
    allocate(omega(len_indx), lambda(len_indx), rhs(nrow), stat=i_alloc, errmsg=error_message)
    if (i_alloc /= 0) then
      write(*,*) 'assemble_spline_matrix_sparse_coo: Allocation failed:', trim(error_message)
      stop
    end if
    
    ! Calculate optimal weights for smoothing
    if (maxval(lambda1) < 0.0_DP) then
      call calc_opt_lambda3(x, y, omega)
    else
      omega = lambda1
    end if
    lambda = 1.0_DP - omega
    
    ! Initialize RHS
    rhs = 0.0_DP
    
    ! Set boundary condition parameters
    call set_boundary_params(sw1, sw2, mu1, mu2, nu1, nu2, sig1, sig2, rho1, rho2)
    
    ! First pass: count actual non-zeros using build_matrix_two_pass
    idx = 0
    i = 0
    call build_matrix_two_pass(.TRUE., idx, i, x, y, m, f, lambda, omega, &
                              indx, mu1, mu2, nu1, nu2, sig1, sig2, rho1, rho2, &
                              c1, cn, VAR, len_indx)
    nnz = idx
    
    ! Allocate COO arrays
    allocate(irow_coo(nnz), icol_coo(nnz), val_coo(nnz))
    
    ! Second pass: fill arrays
    idx = 0
    i = 0
    call build_matrix_two_pass(.FALSE., idx, i, x, y, m, f, lambda, omega, &
                              indx, mu1, mu2, nu1, nu2, sig1, sig2, rho1, rho2, &
                              c1, cn, VAR, len_indx, irow_coo, icol_coo, val_coo, rhs)
    
    deallocate(omega, lambda)
    
  end subroutine assemble_spline_matrix_sparse_coo
  
  !> Add continuity conditions between intervals
  subroutine add_continuity_conditions(counting, idx, i, j, h, VAR, irow, icol, vals)
    logical, intent(in) :: counting
    integer(I4B), intent(inout) :: idx, i
    integer(I4B), intent(in) :: j, VAR
    real(DP), intent(in) :: h
    integer(I4B), dimension(:), intent(inout), optional :: irow, icol
    real(DP), dimension(:), intent(inout), optional :: vals
    
    ! A_i continuity
    i = i + 1
    idx = idx + 1
    if (.not. counting) then
       irow(idx) = i; icol(idx) = j; vals(idx) = 1.0D0
    end if
    idx = idx + 1
    if (.not. counting) then
       irow(idx) = i; icol(idx) = j+1; vals(idx) = h
    end if
    idx = idx + 1
    if (.not. counting) then
       irow(idx) = i; icol(idx) = j+2; vals(idx) = h*h
    end if
    idx = idx + 1
    if (.not. counting) then
       irow(idx) = i; icol(idx) = j+3; vals(idx) = h*h*h
    end if
    idx = idx + 1
    if (.not. counting) then
       irow(idx) = i; icol(idx) = j+VAR; vals(idx) = -1.0D0
    end if
    
    ! B_i continuity
    i = i + 1
    idx = idx + 1
    if (.not. counting) then
       irow(idx) = i; icol(idx) = j+1; vals(idx) = 1.0D0
    end if
    idx = idx + 1
    if (.not. counting) then
       irow(idx) = i; icol(idx) = j+2; vals(idx) = 2.0D0*h
    end if
    idx = idx + 1
    if (.not. counting) then
       irow(idx) = i; icol(idx) = j+3; vals(idx) = 3.0D0*h*h
    end if
    idx = idx + 1
    if (.not. counting) then
       irow(idx) = i; icol(idx) = j+VAR+1; vals(idx) = -1.0D0
    end if
    
    ! C_i continuity
    i = i + 1
    idx = idx + 1
    if (.not. counting) then
       irow(idx) = i; icol(idx) = j+2; vals(idx) = 1.0D0
    end if
    idx = idx + 1
    if (.not. counting) then
       irow(idx) = i; icol(idx) = j+3; vals(idx) = 3.0D0*h
    end if
    idx = idx + 1
    if (.not. counting) then
       irow(idx) = i; icol(idx) = j+VAR+2; vals(idx) = -1.0D0
    end if
  end subroutine add_continuity_conditions
  
  !> Add boundary condition entries
  subroutine add_boundary_condition_1(counting, idx, i, mu1, nu1, sig1, rho1, &
                                       len_indx, VAR, irow, icol, vals)
    logical, intent(in) :: counting
    integer(I4B), intent(inout) :: idx, i
    integer(I4B), intent(in) :: mu1, nu1, sig1, rho1, len_indx, VAR
    integer(I4B), dimension(:), intent(inout), optional :: irow, icol
    real(DP), dimension(:), intent(inout), optional :: vals
    
    i = i + 1
    if (mu1 /= 0) call add_entry(counting, idx, i, 2, dble(mu1), irow, icol, vals)
    if (nu1 /= 0) call add_entry(counting, idx, i, 3, dble(nu1), irow, icol, vals)
    if (sig1 /= 0) call add_entry(counting, idx, i, (len_indx-1)*VAR + 2, dble(sig1), irow, icol, vals)
    if (rho1 /= 0) call add_entry(counting, idx, i, (len_indx-1)*VAR + 3, dble(rho1), irow, icol, vals)
  end subroutine add_boundary_condition_1
  
  !> Add second boundary condition
  subroutine add_boundary_condition_2(counting, idx, i, mu2, nu2, sig1, sig2, rho1, rho2, &
                                      len_indx, VAR, cn, irow, icol, vals, inh)
    logical, intent(in) :: counting
    integer(I4B), intent(inout) :: idx, i
    integer(I4B), intent(in) :: mu2, nu2, sig1, sig2, rho1, rho2, len_indx, VAR
    real(DP), intent(in) :: cn
    integer(I4B), dimension(:), intent(inout), optional :: irow, icol
    real(DP), dimension(:), intent(inout), optional :: vals, inh
    
    ! delta b_i
    i = i + 1
    call add_entry(counting, idx, i, (len_indx-2)*VAR+6, -1.0D0, irow, icol, vals)
    call add_entry(counting, idx, i, (len_indx-1)*VAR+4, dble(sig1), irow, icol, vals)
    call add_entry(counting, idx, i, (len_indx-1)*VAR+5, dble(sig2), irow, icol, vals)
    
    ! delta c_i
    i = i + 1
    call add_entry(counting, idx, i, (len_indx-2)*VAR+7, -1.0D0, irow, icol, vals)
    call add_entry(counting, idx, i, (len_indx-1)*VAR+4, dble(rho1), irow, icol, vals)
    call add_entry(counting, idx, i, (len_indx-1)*VAR+5, dble(rho2), irow, icol, vals)
    
    ! Boundary condition 2
    i = i + 1
    idx = idx + 1
    if (.not. counting) then
       irow(idx) = i; icol(idx) = 2; vals(idx) = dble(mu2)
    end if
    idx = idx + 1
    if (.not. counting) then
       irow(idx) = i; icol(idx) = 3; vals(idx) = dble(nu2)
    end if
    idx = idx + 1
    if (.not. counting) then
       irow(idx) = i; icol(idx) = (len_indx-1)*VAR + 2; vals(idx) = dble(sig2)
    end if
    idx = idx + 1
    if (.not. counting) then
       irow(idx) = i; icol(idx) = (len_indx-1)*VAR + 3; vals(idx) = dble(rho2)
    end if
    if (.not. counting) inh(i) = cn
  end subroutine add_boundary_condition_2
  
  !> Helper: Set boundary condition parameters
  subroutine set_boundary_params(sw1, sw2, mu1, mu2, nu1, nu2, sig1, sig2, rho1, rho2)
    integer(I4B), intent(in) :: sw1, sw2
    integer(I4B), intent(out) :: mu1, mu2, nu1, nu2, sig1, sig2, rho1, rho2
    
    ! First boundary condition
    select case(sw1)
    case(1)
      mu1 = 1; nu1 = 0; sig1 = 0; rho1 = 0
    case(2)
      mu1 = 0; nu1 = 1; sig1 = 0; rho1 = 0
    case(3)
      mu1 = 0; nu1 = 0; sig1 = 1; rho1 = 0
    case(4)
      mu1 = 0; nu1 = 0; sig1 = 0; rho1 = 1
    end select
    
    ! Second boundary condition
    select case(sw2)
    case(1)
      mu2 = 1; nu2 = 0; sig2 = 0; rho2 = 0
    case(2)
      mu2 = 0; nu2 = 1; sig2 = 0; rho2 = 0
    case(3)
      mu2 = 0; nu2 = 0; sig2 = 1; rho2 = 0
    case(4)
      mu2 = 0; nu2 = 0; sig2 = 0; rho2 = 1
    end select
    
  end subroutine set_boundary_params
  
  !> Build matrix using two-pass approach (from splinecof3_direct_sparse)
  !> First pass counts entries, second pass fills arrays
  subroutine build_matrix_two_pass(counting, idx, i, x, y, m, f, lambda, omega, &
                                   indx, mu1, mu2, nu1, nu2, sig1, sig2, rho1, rho2, &
                                   c1, cn, VAR, len_indx, irow, icol, vals, inh)
    logical, intent(in) :: counting
    integer(I4B), intent(inout) :: idx, i
    real(DP), dimension(:), intent(in) :: x, y, lambda, omega
    real(DP), intent(in) :: m, c1, cn
    integer(I4B), dimension(:), intent(in) :: indx
    integer(I4B), intent(in) :: mu1, mu2, nu1, nu2, sig1, sig2, rho1, rho2
    integer(I4B), intent(in) :: VAR, len_indx
    interface
       function f(x,m)
         use nrtype, only : DP
         implicit none
         real(DP), intent(in) :: x, m
         real(DP) :: f
       end function f
    end interface
    integer(I4B), dimension(:), intent(inout), optional :: irow, icol
    real(DP), dimension(:), intent(inout), optional :: vals, inh
    
    integer(I4B) :: j, ii, ie, l
    real(DP) :: help_a, help_inh, h
    
    ! Include all the matrix assembly logic from splinecof3_direct_sparse
    ! This is a simplified version - the full implementation would include
    ! all boundary conditions, continuity conditions, and fitting conditions
    
    ! Boundary condition 1
    call add_boundary_condition_1(counting, idx, i, mu1, nu1, sig1, rho1, &
                                  len_indx, VAR, irow, icol, vals)
    if (.not. counting) inh(i) = c1
    
    ! Process each interval
    do j = 1, VAR*(len_indx-1), VAR
       ii = indx((j-1)/VAR+1)
       ie = indx((j-1)/VAR+2) - 1
       h = x(indx((j-1)/VAR+2)) - x(ii)
       
       if (j == 1) then
          ! First interval - special handling
          call process_first_interval(counting, idx, i, j, ii, ie, x, y, m, f, omega, lambda, &
                                     mu1, mu2, nu1, nu2, VAR, len_indx, indx, irow, icol, vals, inh)
       else
          ! Continuity conditions
          call add_continuity_conditions(counting, idx, i, j, h, VAR, irow, icol, vals)
          
          ! Middle interval fitting conditions
          call process_middle_interval(counting, idx, i, j, ii, ie, x, y, m, f, omega, lambda, &
                                      VAR, indx, irow, icol, vals, inh)
       end if
    end do
    
    ! Process last point
    ii = indx(len_indx)
    ie = ii
    help_a = 0.0D0
    help_inh = 0.0D0
    l = ii
    help_a = help_a + f(x(l),m) * f(x(l),m)
    help_inh = help_inh + f(x(l),m) * y(l)
    
    ! Last point conditions
    i = i + 1
    call add_entry(counting, idx, i, (len_indx-1)*VAR+1, omega(len_indx) * help_a, irow, icol, vals)
    call add_entry(counting, idx, i, (len_indx-2)*VAR+5, omega(len_indx) * (-1.0D0), irow, icol, vals)
    if (.not. counting) inh(i) = omega(len_indx) * help_inh
    
    ! Boundary condition 2
    call add_boundary_condition_2(counting, idx, i, mu2, nu2, sig1, sig2, rho1, rho2, &
                                  len_indx, VAR, cn, irow, icol, vals, inh)
    
  end subroutine build_matrix_two_pass
  
  !> Check if matrix element should be included (matches original dense implementation)
  logical function should_include_element(val)
    real(DP), intent(in) :: val
    ! Original dense implementation adds ALL elements unconditionally
    should_include_element = .TRUE.
  end function should_include_element
  
  !> Add a matrix entry if non-zero (counting mode just increments counter)
  subroutine add_entry(counting, idx, i, j, val, irow, icol, vals)
    logical, intent(in) :: counting
    integer(I4B), intent(inout) :: idx
    integer(I4B), intent(in) :: i, j
    real(DP), intent(in) :: val
    integer(I4B), dimension(:), intent(inout), optional :: irow, icol
    real(DP), dimension(:), intent(inout), optional :: vals
    
    ! Add entry following original dense implementation behavior
    if (should_include_element(val)) then
       idx = idx + 1
       if (.not. counting) then
          irow(idx) = i
          icol(idx) = j
          vals(idx) = val
       end if
    end if
  end subroutine add_entry

  !> Process first interval (simplified - would include full fitting logic)
  subroutine process_first_interval(counting, idx, i, j, ii, ie, x, y, m, f, omega, lambda, &
                                   mu1, mu2, nu1, nu2, VAR, len_indx, indx, irow, icol, vals, inh)
    logical, intent(in) :: counting
    integer(I4B), intent(inout) :: idx, i
    integer(I4B), intent(in) :: j, ii, ie, mu1, mu2, nu1, nu2, VAR, len_indx
    real(DP), dimension(:), intent(in) :: x, y, omega, lambda
    integer(I4B), dimension(:), intent(in) :: indx
    real(DP), intent(in) :: m
    interface
       function f(x,m)
         use nrtype, only : DP
         implicit none
         real(DP), intent(in) :: x, m
         real(DP) :: f
       end function f
    end interface
    integer(I4B), dimension(:), intent(inout), optional :: irow, icol
    real(DP), dimension(:), intent(inout), optional :: vals, inh
    
    ! This is a simplified placeholder - the full implementation would include
    ! all the fitting coefficient calculations from the original
    i = i + 4  ! Skip 4 fitting equations for this interval
    
  end subroutine process_first_interval
  
  !> Process middle interval (simplified - would include full fitting logic)
  subroutine process_middle_interval(counting, idx, i, j, ii, ie, x, y, m, f, omega, lambda, &
                                    VAR, indx, irow, icol, vals, inh)
    logical, intent(in) :: counting
    integer(I4B), intent(inout) :: idx, i
    integer(I4B), intent(in) :: j, ii, ie, VAR
    real(DP), dimension(:), intent(in) :: x, y, omega, lambda
    integer(I4B), dimension(:), intent(in) :: indx
    real(DP), intent(in) :: m
    interface
       function f(x,m)
         use nrtype, only : DP
         implicit none
         real(DP), intent(in) :: x, m
         real(DP) :: f
       end function f
    end interface
    integer(I4B), dimension(:), intent(inout), optional :: irow, icol
    real(DP), dimension(:), intent(inout), optional :: vals, inh
    
    ! This is a simplified placeholder - the full implementation would include
    ! all the fitting coefficient calculations from the original
    i = i + 4  ! Skip 4 fitting equations for this interval
    
  end subroutine process_middle_interval
  
  !> Assemble tridiagonal matrix for fast spline path
  subroutine assemble_spline_matrix_fast_tridiagonal(x, y, c1, cn, sw1, sw2, n, &
                                                     diag, super_diag, sub_diag, rhs)
    real(DP), dimension(:), intent(in) :: x, y
    real(DP), intent(in) :: c1, cn
    integer(I4B), intent(in) :: sw1, sw2, n
    real(DP), dimension(:), allocatable, intent(out) :: diag, super_diag, sub_diag, rhs
    
    integer(I4B) :: i
    real(DP), allocatable :: h(:), alpha(:)
    logical :: natural_start, natural_end, clamped_start, clamped_end
    
    ! Determine boundary condition types
    natural_start = (sw1 == 2)
    natural_end = (sw2 == 4)
    clamped_start = (sw1 == 1)
    clamped_end = (sw2 == 3)
    
    ! Allocate arrays
    allocate(h(n-1), alpha(n), diag(n), super_diag(n-1), sub_diag(n-1), rhs(n))
    
    ! Step 1: Compute h_i = x_{i+1} - x_i
    do i = 1, n-1
      h(i) = x(i+1) - x(i)
    end do
    
    ! Step 2: Compute alpha values
    alpha(1) = 0.0_DP
    do i = 2, n-1
      alpha(i) = 3.0_DP/h(i)*(y(i+1) - y(i)) - 3.0_DP/h(i-1)*(y(i) - y(i-1))
    end do
    alpha(n) = 0.0_DP
    
    ! Step 3: Set up tridiagonal system based on boundary conditions
    if (clamped_start) then
      alpha(1) = 3.0_DP*(y(2) - y(1))/h(1) - 3.0_DP*c1
      diag(1) = 2.0_DP*h(1)
      super_diag(1) = h(1)
      rhs(1) = alpha(1)
    else  ! natural_start
      diag(1) = 1.0_DP
      super_diag(1) = 0.0_DP
      rhs(1) = 0.0_DP
    end if
    
    ! Middle rows
    do i = 2, n-1
      sub_diag(i-1) = h(i-1)
      diag(i) = 2.0_DP*(h(i-1) + h(i))
      if (i < n-1) super_diag(i) = h(i)
      rhs(i) = alpha(i)
    end do
    
    ! Last row
    if (clamped_end) then
      alpha(n) = 3.0_DP*cn - 3.0_DP*(y(n) - y(n-1))/h(n-1)
      sub_diag(n-1) = h(n-1)
      diag(n) = 2.0_DP*h(n-1)
      rhs(n) = alpha(n)
    else  ! natural_end
      sub_diag(n-1) = 0.0_DP
      diag(n) = 1.0_DP
      rhs(n) = 0.0_DP
    end if
    
    deallocate(h, alpha)
    
  end subroutine assemble_spline_matrix_fast_tridiagonal
  
  !> Extract tridiagonal entries from sparse matrix for comparison
  subroutine extract_tridiagonal_from_sparse(n, nnz, irow, icol, vals, &
                                             diag, super_diag, sub_diag)
    integer(I4B), intent(in) :: n, nnz
    integer(I4B), dimension(:), intent(in) :: irow, icol
    real(DP), dimension(:), intent(in) :: vals
    real(DP), dimension(:), allocatable, intent(out) :: diag, super_diag, sub_diag
    
    integer(I4B) :: k
    
    allocate(diag(n), super_diag(n-1), sub_diag(n-1))
    diag = 0.0_DP
    super_diag = 0.0_DP
    sub_diag = 0.0_DP
    
    ! Extract diagonal and off-diagonal elements
    do k = 1, nnz
      if (irow(k) == icol(k)) then
        diag(irow(k)) = vals(k)
      else if (irow(k) == icol(k) - 1) then
        super_diag(irow(k)) = vals(k)
      else if (irow(k) == icol(k) + 1) then
        sub_diag(icol(k)) = vals(k)
      end if
    end do
    
  end subroutine extract_tridiagonal_from_sparse
  
  !> Compare two sparse matrices in COO format
  function compare_sparse_matrices(n1, irow1, icol1, val1, &
                                  n2, irow2, icol2, val2, tol) result(matches)
    integer(I4B), intent(in) :: n1, n2
    integer(I4B), dimension(:), intent(in) :: irow1, icol1, irow2, icol2
    real(DP), dimension(:), intent(in) :: val1, val2
    real(DP), intent(in) :: tol
    logical :: matches
    
    integer(I4B) :: i, j, found_idx
    real(DP) :: max_diff
    
    matches = .true.
    max_diff = 0.0_DP
    
    ! First check if same number of non-zeros
    if (n1 /= n2) then
      write(*,'(A,I0,A,I0)') 'Different number of non-zeros: ', n1, ' vs ', n2
      matches = .false.
      return
    end if
    
    ! Check each entry in matrix 1 exists in matrix 2 with same value
    do i = 1, n1
      found_idx = 0
      do j = 1, n2
        if (irow1(i) == irow2(j) .and. icol1(i) == icol2(j)) then
          found_idx = j
          exit
        end if
      end do
      
      if (found_idx == 0) then
        write(*,'(A,I0,A,I0,A)') 'Entry (', irow1(i), ',', icol1(i), ') not found in second matrix'
        matches = .false.
      else
        if (abs(val1(i) - val2(found_idx)) > tol) then
          max_diff = max(max_diff, abs(val1(i) - val2(found_idx)))
          if (matches) then  ! First difference
            write(*,'(A,I0,A,I0,A,E15.6)') 'First difference at (', irow1(i), ',', icol1(i), '): ', &
                  abs(val1(i) - val2(found_idx))
            write(*,'(A,E15.6,A,E15.6)') '  Matrix 1: ', val1(i), ', Matrix 2: ', val2(found_idx)
          end if
          matches = .false.
        end if
      end if
    end do
    
    if (.not. matches) then
      write(*,'(A,E15.6)') 'Maximum element difference: ', max_diff
    end if
    
  end function compare_sparse_matrices
  
end module spline_matrix_assembly_mod