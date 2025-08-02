module amg_cycles_mod
  use amg_types_mod
  use amg_smoothers_mod
  use sparse_types_mod, only: dp
  implicit none
  
  private
  public :: amg_solve, amg_vcycle_apply, amg_wcycle_apply
  
contains

  subroutine amg_solve(hierarchy, x, b, info)
    type(amg_hierarchy), intent(inout) :: hierarchy
    real(dp), intent(inout) :: x(:)
    real(dp), intent(in) :: b(:)
    integer, intent(out) :: info
    
    integer :: iter
    real(dp) :: rnorm, bnorm, tol
    real(dp), allocatable :: r(:)
    
    info = 0
    
    if (.not. hierarchy%setup_done) then
      info = -1
      return
    end if
    
    allocate(r(hierarchy%levels(1)%n))
    
    ! Compute initial residual
    call compute_residual(hierarchy%levels(1), x, b, r)
    rnorm = norm2(r)
    bnorm = norm2(b)
    
    if (bnorm > 0.0_dp) then
      tol = hierarchy%params%tolerance * bnorm
    else
      tol = hierarchy%params%tolerance
    end if
    
    if (hierarchy%params%verbose) then
      print '(A,I4,A,ES12.5)', "AMG: Initial residual = ", 0, ", ||r|| = ", rnorm
    end if
    
    ! Main iteration loop
    do iter = 1, hierarchy%params%max_iter
      ! Apply one multigrid cycle
      select case(hierarchy%params%cycle_type)
      case(AMG_VCYCLE)
        call amg_vcycle_apply(hierarchy, 1, x, b)
        
      case(AMG_WCYCLE)
        call amg_wcycle_apply(hierarchy, 1, x, b)
        
      case default
        error stop "Unknown cycle type in amg_solve"
      end select
      
      ! Check convergence
      call compute_residual(hierarchy%levels(1), x, b, r)
      rnorm = norm2(r)
      
      if (hierarchy%params%verbose) then
        print '(A,I4,A,ES12.5)', "AMG: Iteration ", iter, ", ||r|| = ", rnorm
      end if
      
      if (rnorm < tol) then
        if (hierarchy%params%verbose) then
          print '(A,I4,A)', "AMG: Converged in ", iter, " iterations"
        end if
        exit
      end if
    end do
    
    if (iter > hierarchy%params%max_iter) then
      info = 1  ! Did not converge
      if (hierarchy%params%verbose) then
        print '(A)', "AMG: Warning - maximum iterations reached"
      end if
    end if
    
    deallocate(r)
    
  end subroutine amg_solve
  
  recursive subroutine amg_vcycle_apply(hierarchy, level_idx, x, b)
    type(amg_hierarchy), intent(inout) :: hierarchy
    integer, intent(in) :: level_idx
    real(dp), intent(inout) :: x(:)
    real(dp), intent(in) :: b(:)
    
    real(dp), allocatable :: r(:), e(:), b_coarse(:)
    
    if (level_idx == hierarchy%n_levels) then
      ! Coarsest level - solve directly
      call coarse_solve(hierarchy%levels(level_idx), x, b)
    else
      ! Pre-smoothing
      call amg_smooth(hierarchy%levels(level_idx), x, b, hierarchy%params, .true.)
      
      ! Compute residual r = b - A*x
      allocate(r(hierarchy%levels(level_idx)%n))
      call compute_residual(hierarchy%levels(level_idx), x, b, r)
      
      ! Restrict residual to coarse grid
      allocate(b_coarse(hierarchy%levels(level_idx+1)%n))
      call restrict_vector(hierarchy%levels(level_idx), r, b_coarse)
      
      ! Solve coarse grid equation A_c * e_c = r_c
      allocate(e(hierarchy%levels(level_idx+1)%n))
      e = 0.0_dp
      call amg_vcycle_apply(hierarchy, level_idx+1, e, b_coarse)
      
      ! Prolongate error to fine grid and correct
      call prolongate_and_correct(hierarchy%levels(level_idx), e, x)
      
      ! Post-smoothing
      call amg_smooth(hierarchy%levels(level_idx), x, b, hierarchy%params, .false.)
      
      deallocate(r, e, b_coarse)
    end if
    
  end subroutine amg_vcycle_apply
  
  recursive subroutine amg_wcycle_apply(hierarchy, level_idx, x, b)
    type(amg_hierarchy), intent(inout) :: hierarchy
    integer, intent(in) :: level_idx
    real(dp), intent(inout) :: x(:)
    real(dp), intent(in) :: b(:)
    
    real(dp), allocatable :: r(:), e(:), b_coarse(:)
    integer :: i
    
    if (level_idx == hierarchy%n_levels) then
      ! Coarsest level - solve directly
      call coarse_solve(hierarchy%levels(level_idx), x, b)
    else
      ! Pre-smoothing
      call amg_smooth(hierarchy%levels(level_idx), x, b, hierarchy%params, .true.)
      
      ! Compute residual r = b - A*x
      allocate(r(hierarchy%levels(level_idx)%n))
      call compute_residual(hierarchy%levels(level_idx), x, b, r)
      
      ! Restrict residual to coarse grid
      allocate(b_coarse(hierarchy%levels(level_idx+1)%n))
      call restrict_vector(hierarchy%levels(level_idx), r, b_coarse)
      
      ! W-cycle: Apply cycle twice on coarse grid
      allocate(e(hierarchy%levels(level_idx+1)%n))
      e = 0.0_dp
      do i = 1, 2
        call amg_wcycle_apply(hierarchy, level_idx+1, e, b_coarse)
      end do
      
      ! Prolongate error to fine grid and correct
      call prolongate_and_correct(hierarchy%levels(level_idx), e, x)
      
      ! Post-smoothing
      call amg_smooth(hierarchy%levels(level_idx), x, b, hierarchy%params, .false.)
      
      deallocate(r, e, b_coarse)
    end if
    
  end subroutine amg_wcycle_apply
  
  subroutine compute_residual(level, x, b, r)
    type(amg_level), intent(in) :: level
    real(dp), intent(in) :: x(:), b(:)
    real(dp), intent(out) :: r(:)
    
    integer :: i, j, k
    real(dp) :: sum
    
    ! r = b - A*x
    do i = 1, level%n
      sum = 0.0_dp
      do k = level%row_ptr(i), level%row_ptr(i+1)-1
        j = level%col_idx(k)
        sum = sum + level%values(k) * x(j)
      end do
      r(i) = b(i) - sum
    end do
    
  end subroutine compute_residual
  
  subroutine restrict_vector(level, v_fine, v_coarse)
    type(amg_level), intent(in) :: level
    real(dp), intent(in) :: v_fine(:)
    real(dp), intent(out) :: v_coarse(:)
    
    integer :: i, j, k
    real(dp) :: sum
    
    ! v_coarse = R * v_fine
    do i = 1, level%n_coarse
      sum = 0.0_dp
      do k = level%R_row_ptr(i), level%R_row_ptr(i+1)-1
        j = level%R_col_idx(k)
        sum = sum + level%R_values(k) * v_fine(j)
      end do
      v_coarse(i) = sum
    end do
    
  end subroutine restrict_vector
  
  subroutine prolongate_and_correct(level, v_coarse, x_fine)
    type(amg_level), intent(in) :: level
    real(dp), intent(in) :: v_coarse(:)
    real(dp), intent(inout) :: x_fine(:)
    
    integer :: i, j, k
    real(dp) :: sum
    
    ! x_fine = x_fine + P * v_coarse
    do i = 1, level%n_fine
      sum = 0.0_dp
      do k = level%P_row_ptr(i), level%P_row_ptr(i+1)-1
        j = level%P_col_idx(k)
        sum = sum + level%P_values(k) * v_coarse(j)
      end do
      x_fine(i) = x_fine(i) + sum
    end do
    
  end subroutine prolongate_and_correct
  
  subroutine coarse_solve(level, x, b)
    type(amg_level), intent(in) :: level
    real(dp), intent(inout) :: x(:)
    real(dp), intent(in) :: b(:)
    
    ! For small coarse problems, use many iterations of smoother
    ! In production, could use direct solver for very small problems
    integer :: i
    
    do i = 1, 20
      call gauss_seidel_forward(level, x, b, 1)
    end do
    
  end subroutine coarse_solve
  
end module amg_cycles_mod