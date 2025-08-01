!> Fast cubic spline implementation using LAPACK tridiagonal solver
!> Based on standard textbook formulations for natural, clamped, and mixed boundary conditions
module splinecof3_fast_mod
  use nrtype, only : I4B, DP
  implicit none
  
  private
  public :: splinecof3_general_fast
  
contains

  !> General fast cubic spline using tridiagonal solver
  !> 
  !> Supports: (2,4) natural, (1,3) clamped, (1,4) mixed, (2,3) mixed
  SUBROUTINE splinecof3_general_fast(x, y, c1, cn, sw1, sw2, a, b, c, d)
    real(DP), dimension(:), intent(in) :: x, y
    real(DP), intent(in) :: c1, cn
    integer(I4B), intent(in) :: sw1, sw2
    real(DP), dimension(:), intent(out) :: a, b, c, d

    integer(I4B) :: info, n, i
    real(DP), allocatable :: h(:), alpha(:), l(:), mu(:), z(:), c_work(:)
    logical :: natural_start, natural_end, clamped_start, clamped_end

    n = size(x)
    
    ! Determine boundary condition types
    natural_start = (sw1 == 2)   ! S''(x1) = 0
    natural_end = (sw2 == 4)     ! S''(xn) = 0  
    clamped_start = (sw1 == 1)   ! S'(x1) = c1
    clamped_end = (sw2 == 3)     ! S'(xn) = cn
    
    ! Validate supported combinations
    if (.not. ((sw1 == 2 .and. sw2 == 4) .or. &  ! Natural
               (sw1 == 1 .and. sw2 == 3) .or. &  ! Clamped
               (sw1 == 1 .and. sw2 == 4) .or. &  ! Mixed: clamped start, natural end
               (sw1 == 2 .and. sw2 == 3))) then  ! Mixed: natural start, clamped end
      write(*,'(A,2I0)') 'splinecof3_general_fast: ERROR - Unsupported boundary combination sw1=', sw1, ', sw2=', sw2
      error stop 'splinecof3_general_fast: Invalid boundary conditions'
    end if
    
    ! Validate inputs
    if (size(y) /= n .or. size(a) /= n .or. size(b) /= n .or. &
        size(c) /= n .or. size(d) /= n .or. n < 3) then
      error stop 'splinecof3_general_fast: Array size mismatch or insufficient points'
    end if
    
    do i = 1, n-1
      if (x(i) >= x(i+1)) then
        error stop 'splinecof3_general_fast: Non-monotonic x values'
      end if
    end do

    ! Allocate work arrays
    allocate(h(n-1), alpha(n), l(n), mu(n), z(n), c_work(n))
    
    ! Step 1: Compute h_i = x_{i+1} - x_i
    do i = 1, n-1
      h(i) = x(i+1) - x(i)
    end do
    
    ! Step 2: Compute alpha values based on boundary conditions
    alpha(1) = 0.0_DP  ! Will be set based on boundary condition
    do i = 2, n-1
      alpha(i) = 3.0_DP/h(i)*(y(i+1) - y(i)) - 3.0_DP/h(i-1)*(y(i) - y(i-1))
    end do
    alpha(n) = 0.0_DP  ! Will be set based on boundary condition
    
    ! Step 3: Set up tridiagonal system based on boundary conditions
    if (clamped_start) then
      alpha(1) = 3.0_DP*(y(2) - y(1))/h(1) - 3.0_DP*c1
      l(1) = 2.0_DP*h(1)
      mu(1) = 0.5_DP
      z(1) = alpha(1)/l(1)
    else  ! natural_start
      l(1) = 1.0_DP
      mu(1) = 0.0_DP
      z(1) = 0.0_DP
    end if
    
    ! Step 4: Forward elimination
    do i = 2, n-1
      if (clamped_start .or. i > 2) then
        l(i) = 2.0_DP*(x(i+1) - x(i-1)) - h(i-1)*mu(i-1)
        mu(i) = h(i)/l(i)
        z(i) = (alpha(i) - h(i-1)*z(i-1))/l(i)
      else  ! i = 2 and natural_start
        l(i) = 2.0_DP*(x(i+1) - x(i-1))
        mu(i) = h(i)/l(i)
        z(i) = alpha(i)/l(i)
      end if
    end do
    
    ! Step 5: Set final values based on end boundary condition
    if (clamped_end) then
      alpha(n) = 3.0_DP*cn - 3.0_DP*(y(n) - y(n-1))/h(n-1)
      l(n) = h(n-1)*(2.0_DP - mu(n-1))
      z(n) = (alpha(n) - h(n-1)*z(n-1))/l(n)
      c_work(n) = z(n)
    else  ! natural_end
      l(n) = 1.0_DP
      z(n) = 0.0_DP
      c_work(n) = 0.0_DP
    end if
    
    ! Step 6: Back substitution
    if (natural_end) then
      c_work(n-1) = z(n-1)
    else  ! clamped_end
      c_work(n-1) = z(n-1) - mu(n-1)*c_work(n)
    end if
    
    do i = n-2, 1, -1
      if (natural_start .and. i == 1) then
        c_work(i) = 0.0_DP
      else
        c_work(i) = z(i) - mu(i)*c_work(i+1)
      end if
    end do
    
    ! Step 7: Compute spline coefficients
    ! a_i = y_i
    a(1:n-1) = y(1:n-1)
    
    ! c_i = c_work_i (second derivatives)
    c(1:n-1) = c_work(1:n-1)
    
    ! b_i and d_i
    do i = 1, n-1
      d(i) = (c_work(i+1) - c_work(i))/(3.0_DP*h(i))
      b(i) = (y(i+1) - y(i))/h(i) - h(i)*(c_work(i+1) + 2.0_DP*c_work(i))/3.0_DP
    end do
    
    ! Override b values for clamped boundaries
    if (clamped_start) then
      b(1) = c1
    end if
    if (clamped_end) then
      b(n-1) = cn
    end if
    
    ! Follow spline_cof convention: set n-th element to zero
    a(n) = 0.0_DP
    b(n) = 0.0_DP
    c(n) = 0.0_DP
    d(n) = 0.0_DP

    ! Clean up
    deallocate(h, alpha, l, mu, z, c_work)
    
  END SUBROUTINE splinecof3_general_fast

end module splinecof3_fast_mod