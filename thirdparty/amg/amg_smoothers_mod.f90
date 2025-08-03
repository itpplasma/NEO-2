module amg_smoothers_mod
  use amg_types_mod
  use sparse_types_mod, only: dp
  implicit none
  
  private
  public :: amg_smooth, gauss_seidel_forward, gauss_seidel_backward
  public :: symmetric_gauss_seidel, jacobi_smooth
  
contains

  subroutine amg_smooth(level, x, b, params, presweep)
    type(amg_level), intent(in) :: level
    real(dp), intent(inout) :: x(:)
    real(dp), intent(in) :: b(:)
    type(amg_params), intent(in) :: params
    logical, intent(in) :: presweep
    
    integer :: n_sweeps
    
    if (presweep) then
      n_sweeps = params%presmoother_steps
    else
      n_sweeps = params%postsmoother_steps
    end if
    
    select case(params%smoother)
    case(AMG_JACOBI)
      call jacobi_smooth(level, x, b, n_sweeps, params%jacobi_weight)
      
    case(AMG_GAUSS_SEIDEL)
      call gauss_seidel_forward(level, x, b, n_sweeps)
      
    case(AMG_SYM_GAUSS_SEIDEL)
      call symmetric_gauss_seidel(level, x, b, n_sweeps)
      
    case default
      error stop "Unknown smoother type in amg_smooth"
    end select
    
  end subroutine amg_smooth
  
  subroutine gauss_seidel_forward(level, x, b, n_sweeps)
    type(amg_level), intent(in) :: level
    real(dp), intent(inout) :: x(:)
    real(dp), intent(in) :: b(:)
    integer, intent(in) :: n_sweeps
    
    integer :: sweep, i, j, k
    real(dp) :: rsum, diag
    
    do sweep = 1, n_sweeps
      do i = 1, level%n
        rsum = b(i)
        diag = 0.0_dp
        
        ! Compute residual: r = b - A*x
        do k = level%row_ptr(i), level%row_ptr(i+1)-1
          j = level%col_idx(k)
          if (i == j) then
            diag = level%values(k)
          else
            rsum = rsum - level%values(k) * x(j)
          end if
        end do
        
        ! Update x(i) = x(i) + r(i)/A(i,i)
        ! Handle zero diagonal case (keep old value)
        if (diag /= 0.0_dp) then
          x(i) = rsum / diag
        end if
      end do
    end do
    
  end subroutine gauss_seidel_forward
  
  subroutine gauss_seidel_backward(level, x, b, n_sweeps)
    type(amg_level), intent(in) :: level
    real(dp), intent(inout) :: x(:)
    real(dp), intent(in) :: b(:)
    integer, intent(in) :: n_sweeps
    
    integer :: sweep, i, j, k
    real(dp) :: rsum, diag
    
    do sweep = 1, n_sweeps
      do i = level%n, 1, -1
        rsum = b(i)
        diag = 0.0_dp
        
        ! Compute residual: r = b - A*x
        do k = level%row_ptr(i), level%row_ptr(i+1)-1
          j = level%col_idx(k)
          if (i == j) then
            diag = level%values(k)
          else
            rsum = rsum - level%values(k) * x(j)
          end if
        end do
        
        ! Update x(i) = x(i) + r(i)/A(i,i)
        ! Handle zero diagonal case (keep old value)
        if (diag /= 0.0_dp) then
          x(i) = rsum / diag
        end if
      end do
    end do
    
  end subroutine gauss_seidel_backward
  
  subroutine symmetric_gauss_seidel(level, x, b, n_sweeps)
    type(amg_level), intent(in) :: level
    real(dp), intent(inout) :: x(:)
    real(dp), intent(in) :: b(:)
    integer, intent(in) :: n_sweeps
    
    integer :: sweep
    
    do sweep = 1, n_sweeps
      ! Forward sweep
      call gauss_seidel_forward(level, x, b, 1)
      ! Backward sweep
      call gauss_seidel_backward(level, x, b, 1)
    end do
    
  end subroutine symmetric_gauss_seidel
  
  subroutine jacobi_smooth(level, x, b, n_sweeps, omega)
    type(amg_level), intent(in) :: level
    real(dp), intent(inout) :: x(:)
    real(dp), intent(in) :: b(:)
    integer, intent(in) :: n_sweeps
    real(dp), intent(in) :: omega
    
    integer :: sweep, i, j, k
    real(dp) :: rsum, diag
    real(dp), allocatable :: x_new(:)
    
    allocate(x_new(level%n))
    
    do sweep = 1, n_sweeps
      x_new = x
      
      do i = 1, level%n
        rsum = b(i)
        diag = 0.0_dp
        
        ! Compute residual: r = b - A*x
        do k = level%row_ptr(i), level%row_ptr(i+1)-1
          j = level%col_idx(k)
          if (i == j) then
            diag = level%values(k)
          else
            rsum = rsum - level%values(k) * x(j)
          end if
        end do
        
        ! Update with relaxation: x_new = x + omega * r/diag
        if (diag /= 0.0_dp) then
          x_new(i) = x(i) + omega * rsum / diag
        end if
      end do
      
      x = x_new
    end do
    
    deallocate(x_new)
    
  end subroutine jacobi_smooth
  
end module amg_smoothers_mod