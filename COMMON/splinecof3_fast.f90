!> Fast natural cubic spline implementation using LAPACK tridiagonal solver
!> This module provides a high-performance implementation for the special case
!> of natural cubic splines (zero second derivatives at endpoints) on uniformly
!> spaced or arbitrary data points.
module splinecof3_fast_mod
  use nrtype, only : I4B, DP
  implicit none
  
  private
  public :: splinecof3_fast
  
contains

  !> Fast natural cubic spline coefficient calculation
  !> 
  !> This routine implements natural cubic spline interpolation using LAPACK's
  !> optimized tridiagonal solver (dptsv). It's significantly faster than the
  !> general sparse matrix approach for the natural boundary condition case.
  !>
  !> INPUT:
  !>   x(:) - knot positions (must be strictly increasing)
  !>   y(:) - function values at knots
  !>
  !> OUTPUT:
  !>   a(:) - spline coefficients (size n-1)
  !>   b(:) - spline coefficients (size n-1) 
  !>   c(:) - spline coefficients (size n-1)
  !>   d(:) - spline coefficients (size n-1)
  !>
  !> The spline in interval i is: 
  !>   S(t) = a(i) + b(i)*(t-x(i)) + c(i)*(t-x(i))^2 + d(i)*(t-x(i))^3
  !>   for t in [x(i), x(i+1)]
  SUBROUTINE splinecof3_fast(x, y, a, b, c, d)
    real(DP), dimension(:), intent(in) :: x, y
    real(DP), dimension(:), intent(out) :: a, b, c, d

    integer(I4B) :: info, n, i
    real(DP), allocatable :: r(:), h(:), dl(:), ds(:), cs(:)
    character(100) :: error_msg

    n = size(x)
    
    ! Validate input arrays
    if (size(y) /= n) then
      write(*,'(A)') 'splinecof3_fast: ERROR - x and y arrays must have same size'
      error stop 'splinecof3_fast: Array size mismatch'
    end if
    
    if (size(a) /= n-1 .or. size(b) /= n-1 .or. size(c) /= n-1 .or. size(d) /= n-1) then
      write(*,'(A)') 'splinecof3_fast: ERROR - output arrays must have size n-1'
      write(*,'(A,I0,A,I0)') 'Expected size: ', n-1, ', got sizes: a=', size(a)
      error stop 'splinecof3_fast: Output array size mismatch'
    end if
    
    if (n < 3) then
      write(*,'(A,I0)') 'splinecof3_fast: ERROR - need at least 3 points, got ', n
      error stop 'splinecof3_fast: Insufficient data points'
    end if
    
    ! Check that x values are strictly increasing
    do i = 1, n-1
      if (x(i) >= x(i+1)) then
        write(*,'(A,I0,A,ES15.6,A,ES15.6)') 'splinecof3_fast: ERROR - x values not increasing at i=', &
             i, ': x(', x(i), ') >= x(', x(i+1)
        error stop 'splinecof3_fast: Non-monotonic x values'
      end if
    end do

    ! Allocate work arrays
    allocate(h(n-1), r(n-1), dl(n-3), ds(n-2), cs(n-2))

    ! Calculate intervals and differences
    h = x(2:n) - x(1:n-1)
    r = y(2:n) - y(1:n-1)

    ! Set up tridiagonal system for natural cubic spline
    ! The system solves for second derivatives at interior points
    if (n > 2) then
      dl = h(2:n-2)                                    ! sub-diagonal
      ds = 2.0_DP * (h(1:n-2) + h(2:n-1))             ! main diagonal  
      cs = 3.0_DP * (r(2:n-1)/h(2:n-1) - r(1:n-2)/h(1:n-2))  ! RHS

      ! Solve tridiagonal system using LAPACK
      call dptsv(n-2, 1, ds, dl, cs, n-2, info)

      if (info /= 0) then
        if (info < 0) then
          write(*,'(A,I0,A)') 'splinecof3_fast: LAPACK dptsv error - illegal value in argument ', -info
          error stop 'splinecof3_fast: Invalid argument to dptsv'
        else
          write(*,'(A,I0,A)') 'splinecof3_fast: LAPACK dptsv error - diagonal element ', info, ' is zero'
          write(*,'(A)') 'The tridiagonal system is singular and cannot be solved.'
          write(*,'(A)') 'This may indicate duplicate x values or other data issues.'
          error stop 'splinecof3_fast: Singular tridiagonal system in dptsv'
        end if
      end if
    end if

    ! Calculate spline coefficients
    ! a(i) = y(i) for each interval
    a(1:n-1) = y(1:n-1)
    
    ! c(i) values: 0 at endpoints (natural boundary), cs from solver for interior
    c(1) = 0.0_DP
    if (n > 2) then
      c(2:n-1) = cs(1:n-2)
    end if
    
    ! b(i) and d(i) coefficients from natural cubic spline formulas
    b(1) = r(1)/h(1) - h(1)/3.0_DP * c(2)
    if (n > 2) then
      do i = 2, n-2
        b(i) = r(i)/h(i) - h(i)/3.0_DP * (c(i+1) + 2.0_DP*c(i))
      end do
    end if
    b(n-1) = r(n-1)/h(n-1) - h(n-1)/3.0_DP * (2.0_DP*c(n-1))
    
    d(1) = c(2)/(3.0_DP*h(1))
    if (n > 2) then
      do i = 2, n-2
        d(i) = (c(i+1) - c(i))/(3.0_DP*h(i))
      end do
    end if
    d(n-1) = -c(n-1)/(3.0_DP*h(n-1))

    ! Clean up
    deallocate(h, r, dl, ds, cs)
    
  END SUBROUTINE splinecof3_fast

end module splinecof3_fast_mod