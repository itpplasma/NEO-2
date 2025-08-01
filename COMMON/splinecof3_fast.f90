!> Fast natural cubic spline implementation using LAPACK tridiagonal solver
!> This module provides a high-performance implementation for the special case
!> of natural cubic splines (zero second derivatives at endpoints) on uniformly
!> spaced or arbitrary data points.
module splinecof3_fast_mod
  use nrtype, only : I4B, DP
  implicit none
  
  private
  public :: splinecof3_general_fast
  
contains


  !> General fast cubic spline using the proven natural spline as base
  !> 
  !> This reuses the working natural spline implementation with boundary condition modifications
  !> Supports: (2,4) natural, (1,3) clamped, (1,4) mixed, (2,3) mixed
  SUBROUTINE splinecof3_general_fast(x, y, c1, cn, sw1, sw2, a, b, c, d)
    real(DP), dimension(:), intent(in) :: x, y
    real(DP), intent(in) :: c1, cn
    integer(I4B), intent(in) :: sw1, sw2
    real(DP), dimension(:), intent(out) :: a, b, c, d

    integer(I4B) :: info, n, i
    real(DP), allocatable :: h(:), r(:), dl(:), ds(:), cs(:)
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
    
    ! Follow spline_cof convention: coefficient arrays have size n, but only use n-1 elements
    if (size(y) /= n .or. size(a) /= n .or. size(b) /= n .or. &
        size(c) /= n .or. size(d) /= n .or. n < 3) then
      error stop 'splinecof3_general_fast: Array size mismatch or insufficient points'
    end if
    
    do i = 1, n-1
      if (x(i) >= x(i+1)) then
        error stop 'splinecof3_general_fast: Non-monotonic x values'
      end if
    end do

    ! Allocate work arrays using reference implementation sizing
    allocate(h(n-1), r(n-1), dl(n-3), ds(n-2), cs(n-2))
    
    ! Base setup from working natural spline reference
    h = x(2:) - x(1:n-1)
    r = y(2:) - y(1:n-1)
    
    dl = h(2:n-2)
    ds = 2.0_DP*(h(1:n-2) + h(2:))
    
    ! RHS from working natural spline reference
    cs = 3.0_DP*(r(2:)/h(2:) - r(1:n-2)/h(1:n-2))
    
    ! Simple boundary condition modifications (minimal changes from natural case)
    if (clamped_start) then
      ! Modify first equation RHS for clamped start boundary condition
      cs(1) = cs(1) - 3.0_DP * ((y(2)-y(1))/h(1) - c1)
    end if
    
    if (clamped_end) then
      ! Modify last equation RHS for clamped end boundary condition  
      cs(n-2) = cs(n-2) - 3.0_DP * (cn - (y(n)-y(n-1))/h(n-1))
    end if

    ! Reuse exact same solver call as working natural spline
    call dptsv(n-2, 1, ds, dl, cs, n-2, info)

    ! Reuse exact same error handling as working natural spline
    if (info /= 0) then
      if (info < 0) then
        write(*,'(A,I0,A)') 'splinecof3_general_fast: LAPACK dptsv error - illegal value in argument ', -info, '.'
        error stop 'splinecof3_general_fast: Invalid argument to dptsv'
      else
        write(*,'(A,I0,A)') 'splinecof3_general_fast: LAPACK dptsv error - diagonal element ', info, ' is zero.'
        write(*,*) 'The tridiagonal system is singular and cannot be solved.'
        error stop 'splinecof3_general_fast: Singular tridiagonal system in dptsv'
      end if
    end if

    ! Extract coefficients using working natural spline reference formulas
    a(1:n-1) = y(1:n-1)
    
    ! b coefficients with boundary condition handling
    if (clamped_start) then
      b(1) = c1  ! Specified first derivative
    else
      b(1) = r(1)/h(1) - h(1)/3.0_DP*cs(1)  ! Natural start
    end if
    
    b(2:n-2) = r(2:n-2)/h(2:n-2) - h(2:n-2)/3.0_DP*(cs(2:n-2) + 2.0_DP*cs(1:n-3))
    
    if (clamped_end) then
      b(n-1) = cn  ! Specified first derivative
    else
      b(n-1) = r(n-1)/h(n-1) - h(n-1)/3.0_DP*(2.0_DP*cs(n-2))  ! Natural end
    end if
    
    ! c coefficients
    if (natural_start) then
      c(1) = 0.0_DP  ! Natural boundary
    else
      ! For clamped start, compute boundary second derivative
      c(1) = 3.0_DP/h(1)*(r(1)/h(1) - c1) - cs(1)/2.0_DP
    end if
    
    c(2:n-1) = cs  ! Interior second derivatives from tridiagonal solve
    
    ! d coefficients
    d(1) = 1.0_DP/(3.0_DP*h(1))*cs(1)
    d(2:n-2) = 1.0_DP/(3.0_DP*h(2:n-2))*(cs(2:n-2) - cs(1:n-3))
    d(n-1) = 1.0_DP/(3.0_DP*h(n-1))*(-cs(n-2))

    ! Follow spline_cof convention: set n-th element to zero
    a(n) = 0.0_DP
    b(n) = 0.0_DP
    c(n) = 0.0_DP
    d(n) = 0.0_DP

    ! Clean up
    deallocate(h, r, dl, ds, cs)
    
  END SUBROUTINE splinecof3_general_fast

end module splinecof3_fast_mod