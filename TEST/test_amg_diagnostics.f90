program test_amg_diagnostics
  use sparse_types_mod
  use amg_types_mod  
  use amg_precond_mod
  use spline_mod
  implicit none

  integer, parameter :: n_data = 25
  real(dp), parameter :: lambda = 1.0e-4_dp
  
  integer :: n, nnz, info
  integer, allocatable :: row_ptr(:), col_idx(:)
  real(dp), allocatable :: values(:), x_data(:), y_data(:)
  type(amg_hierarchy) :: hierarchy
  
  print *, "=== AMG Coarse Operator Diagnostics ==="
  print *, "Building spline matrix with n_data =", n_data
  print *, "Lambda =", lambda
  
  ! Create spline data
  allocate(x_data(n_data), y_data(n_data))
  call create_test_spline_data(n_data, x_data, y_data)
  
  ! Build spline matrix  
  call create_spline_matrix(x_data, y_data, lambda, n, nnz, row_ptr, col_idx, values, info)
  if (info /= 0) then
    print *, "ERROR: Failed to create spline matrix"
    stop 1
  end if
  
  print *, "Spline matrix: n =", n, "nnz =", nnz
  
  ! Setup AMG with verbose output
  call amg_precond_setup(hierarchy, n, nnz, row_ptr, col_idx, values)
  
  print *, "=== Analysis Complete ==="
  
  ! Cleanup
  call amg_precond_destroy(hierarchy)
  deallocate(row_ptr, col_idx, values, x_data, y_data)

contains

  subroutine create_test_spline_data(n_points, x, y)
    integer, intent(in) :: n_points
    real(dp), intent(out) :: x(:), y(:)
    integer :: i
    real(dp) :: t
    
    do i = 1, n_points
      t = real(i-1, dp) / real(n_points-1, dp)
      x(i) = t
      y(i) = sin(4.0_dp * 3.14159_dp * t) + 0.1_dp * cos(20.0_dp * 3.14159_dp * t)
    end do
  end subroutine create_test_spline_data

end program test_amg_diagnostics