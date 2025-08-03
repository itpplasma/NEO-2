module amg_precond_mod
  use amg_types_mod
  use amg_smoothed_aggregation_mod
  use amg_cycles_mod
  use sparse_types_mod, only: dp
  implicit none
  
  private
  public :: amg_precond_setup, amg_precond_apply, amg_precond_destroy
  public :: amg_precond_info, adapt_amg_parameters
  
contains

  subroutine amg_precond_setup(hierarchy, n, nnz, row_ptr, col_idx, values, params)
    type(amg_hierarchy), intent(out) :: hierarchy
    integer, intent(in) :: n, nnz
    integer, intent(in) :: row_ptr(n+1), col_idx(nnz)
    real(dp), intent(in) :: values(nnz)
    type(amg_params), intent(in), optional :: params
    
    type(amg_level) :: dummy_level
    
    ! Set parameters
    if (present(params)) then
      hierarchy%params = params
    else
      ! Use defaults from amg_types_mod
      hierarchy%params = amg_params()
    end if
    
    ! Enable verbose for debugging
    hierarchy%params%verbose = .true.
    
    ! Adaptive parameter selection based on matrix properties
    call adapt_amg_parameters(hierarchy%params, n, nnz, row_ptr, col_idx, values)
    
    ! Build AMG hierarchy
    select case(hierarchy%params%method)
    case(AMG_SMOOTHED_AGGREGATION)
      call sa_amg_setup(hierarchy, dummy_level, n, nnz, row_ptr, col_idx, values)
      
    case(AMG_CLASSICAL)
      ! Classical AMG not yet implemented
      error stop "Classical AMG not yet implemented"
      
    case default
      error stop "Unknown AMG method in amg_precond_setup"
    end select
    
    if (hierarchy%params%verbose) then
      call amg_precond_info(hierarchy)
    end if
    
  end subroutine amg_precond_setup
  
  subroutine amg_precond_apply(hierarchy, x, b)
    type(amg_hierarchy), intent(inout) :: hierarchy
    real(dp), intent(inout) :: x(:)
    real(dp), intent(in) :: b(:)
    
    integer :: info
    
    ! Apply AMG as a preconditioner (one or more cycles)
    ! For preconditioning, typically use just one V-cycle
    select case(hierarchy%params%cycle_type)
    case(AMG_VCYCLE)
      call amg_vcycle_apply(hierarchy, 1, x, b)
      
    case(AMG_WCYCLE)
      call amg_wcycle_apply(hierarchy, 1, x, b)
      
    case default
      error stop "Unknown cycle type in amg_precond_apply"
    end select
    
  end subroutine amg_precond_apply
  
  subroutine amg_precond_destroy(hierarchy)
    type(amg_hierarchy), intent(inout) :: hierarchy
    
    integer :: i
    
    if (allocated(hierarchy%levels)) then
      do i = 1, hierarchy%n_levels
        if (allocated(hierarchy%levels(i)%row_ptr)) deallocate(hierarchy%levels(i)%row_ptr)
        if (allocated(hierarchy%levels(i)%col_idx)) deallocate(hierarchy%levels(i)%col_idx)
        if (allocated(hierarchy%levels(i)%values)) deallocate(hierarchy%levels(i)%values)
        if (allocated(hierarchy%levels(i)%P_row_ptr)) deallocate(hierarchy%levels(i)%P_row_ptr)
        if (allocated(hierarchy%levels(i)%P_col_idx)) deallocate(hierarchy%levels(i)%P_col_idx)
        if (allocated(hierarchy%levels(i)%P_values)) deallocate(hierarchy%levels(i)%P_values)
        if (allocated(hierarchy%levels(i)%R_row_ptr)) deallocate(hierarchy%levels(i)%R_row_ptr)
        if (allocated(hierarchy%levels(i)%R_col_idx)) deallocate(hierarchy%levels(i)%R_col_idx)
        if (allocated(hierarchy%levels(i)%R_values)) deallocate(hierarchy%levels(i)%R_values)
        if (allocated(hierarchy%levels(i)%S_row_ptr)) deallocate(hierarchy%levels(i)%S_row_ptr)
        if (allocated(hierarchy%levels(i)%S_col_idx)) deallocate(hierarchy%levels(i)%S_col_idx)
        if (allocated(hierarchy%levels(i)%S_values)) deallocate(hierarchy%levels(i)%S_values)
        if (allocated(hierarchy%levels(i)%diagonal)) deallocate(hierarchy%levels(i)%diagonal)
        if (allocated(hierarchy%levels(i)%res_work)) deallocate(hierarchy%levels(i)%res_work)
        if (allocated(hierarchy%levels(i)%sol_work)) deallocate(hierarchy%levels(i)%sol_work)
        if (allocated(hierarchy%levels(i)%aggregates)) deallocate(hierarchy%levels(i)%aggregates)
      end do
      deallocate(hierarchy%levels)
    end if
    
    hierarchy%n_levels = 0
    hierarchy%setup_done = .false.
    
  end subroutine amg_precond_destroy
  
  subroutine amg_precond_info(hierarchy)
    type(amg_hierarchy), intent(in) :: hierarchy
    
    integer :: i
    character(len=20) :: method_name
    
    select case(hierarchy%params%method)
    case(AMG_SMOOTHED_AGGREGATION)
      method_name = "Smoothed Aggregation"
    case(AMG_CLASSICAL)
      method_name = "Classical Ruge-Stuben"
    case default
      method_name = "Unknown"
    end select
    
    print '(A)', "=========================================="
    print '(A)', "AMG Preconditioner Information:"
    print '(A)', "=========================================="
    print '(A,A)', "Method: ", trim(method_name)
    print '(A,I0)', "Number of levels: ", hierarchy%n_levels
    print '(A,F6.3)', "Grid complexity: ", hierarchy%grid_complexity
    print '(A,F6.3)', "Operator complexity: ", hierarchy%operator_complexity
    print '(A)', ""
    print '(A)', "Level information:"
    print '(A)', "Level    Size    Nonzeros    Aggregates"
    print '(A)', "----------------------------------------"
    
    do i = 1, hierarchy%n_levels
      if (i < hierarchy%n_levels) then
        print '(I5,I8,I12,I12)', i, hierarchy%levels(i)%n, &
              hierarchy%levels(i)%nnz, hierarchy%levels(i)%n_aggregates
      else
        print '(I5,I8,I12,A12)', i, hierarchy%levels(i)%n, &
              hierarchy%levels(i)%nnz, "    -"
      end if
    end do
    
    print '(A)', "=========================================="
    
  end subroutine amg_precond_info
  
  subroutine adapt_amg_parameters(params, n, nnz, row_ptr, col_idx, values)
    type(amg_params), intent(inout) :: params
    integer, intent(in) :: n, nnz
    integer, intent(in) :: row_ptr(n+1), col_idx(nnz)
    real(dp), intent(in) :: values(nnz)
    
    real(dp) :: sparsity_ratio, max_off_diag, diag_dominance
    real(dp) :: cond_estimate, bandwidth_ratio
    integer :: i, j, k, bandwidth, max_bandwidth
    logical :: is_symmetric, is_spd
    
    ! Analyze matrix properties
    sparsity_ratio = real(nnz, dp) / real(n*n, dp)
    max_off_diag = 0.0_dp
    diag_dominance = 0.0_dp
    max_bandwidth = 0
    is_symmetric = .true.
    is_spd = .true.
    
    ! Check matrix properties
    do i = 1, n
      do k = row_ptr(i), row_ptr(i+1)-1
        j = col_idx(k)
        
        if (i == j) then
          ! Diagonal entry
          if (values(k) <= 0.0_dp) is_spd = .false.
        else
          ! Off-diagonal entry
          max_off_diag = max(max_off_diag, abs(values(k)))
          max_bandwidth = max(max_bandwidth, abs(i - j))
          
          ! Simple symmetry check (approximate)
          if (abs(values(k)) > 1.0e-12_dp) then
            ! This is a simplified check - full symmetry would require
            ! checking if A(j,i) exists and equals A(i,j)
          end if
        end if
      end do
    end do
    
    bandwidth_ratio = real(max_bandwidth, dp) / real(n, dp)
    
    ! Adaptive parameter selection based on matrix characteristics
    
    ! For sparse, banded matrices (typical of finite differences)
    if (sparsity_ratio < 0.01_dp .and. bandwidth_ratio < 0.1_dp) then
      params%strength_threshold = 0.25_dp
      params%max_coarse_ratio = 0.5_dp
      params%prolongation_damping = 1.33_dp
      
    ! For dense-ish or ill-conditioned matrices (splines, etc.)
    else if (sparsity_ratio > 0.05_dp .or. bandwidth_ratio > 0.3_dp) then
      params%strength_threshold = 0.1_dp  ! Lower threshold for more connections
      params%max_coarse_ratio = 0.7_dp    ! Allow more aggressive coarsening
      params%prolongation_damping = 0.67_dp  ! More conservative damping
      params%max_levels = 10  ! Limit levels for ill-conditioned problems
      
    ! For very large, sparse matrices
    else if (n > 10000 .and. sparsity_ratio < 0.001_dp) then
      params%strength_threshold = 0.5_dp  ! Aggressive thresholding
      params%max_coarse_ratio = 0.3_dp    
      params%prolongation_damping = 1.5_dp
      
    ! Default case - use moderate settings
    else
      params%strength_threshold = 0.25_dp
      params%max_coarse_ratio = 0.6_dp
      params%prolongation_damping = 1.0_dp
    end if
    
    ! Adjust for suspected spline matrices (dense rows, smooth operators)
    if (bandwidth_ratio > 0.5_dp .and. sparsity_ratio > 0.1_dp) then
      params%strength_threshold = 0.05_dp  ! Very permissive for splines
      params%coarsest_size = max(20, n/100)  ! Smaller coarsest level
      params%presmoother_steps = 1  ! Fewer smoothing steps
      params%postsmoother_steps = 1
    end if
    
  end subroutine adapt_amg_parameters
  
end module amg_precond_mod