module radial_lagrange_cof
  ! Optional local cubic Lagrange replacement for the radial Boozer splines.
  ! NEO-2 evaluates them at one surface, so the banded solve over all ns
  ! surfaces per mode is unnecessary. Interval k gets the cubic through nodes
  ! k-1..k+2 (shifted inwards at the edges), stored as the a, b, c, d
  ! coefficients of p(x - x(k)) that splint_horner3 evaluates.
  use nrtype, only : I4B, DP
  use neo_spline_data, only : lsw_lagrange_boozer

  implicit none

  private
  public :: lagrange_cof3, radial_cof3, radial_cof3_hi_driv

  abstract interface
    function lead_function(x, m)
      import :: DP
      real(DP), intent(in) :: x
      real(DP), intent(in) :: m
      real(DP) :: lead_function
    end function lead_function
  end interface

contains

  subroutine lagrange_cof3(x, y, m, a, b, c, d, indx, f)
    ! Coefficients of y/f(x, m) at x(indx). As in splinecof3_hi_driv, the
    ! first node is skipped for m /= 0, where y/f is undefined at the axis.
    real(DP), dimension(:), intent(in) :: x, y
    real(DP), intent(in) :: m
    real(DP), dimension(:), intent(out) :: a, b, c, d
    integer(I4B), dimension(:), intent(in) :: indx
    procedure(lead_function) :: f

    real(DP) :: xn(size(indx)), gn(size(indx)), t(4), dd(4), cf(4)
    integer(I4B) :: n, k, j, i, kmin, j0, np

    n = size(indx)
    kmin = 1
    if (m /= 0.0_DP) kmin = 2
    np = min(4, n - kmin + 1)

    xn = x(indx)
    gn(1) = 0.0_DP
    do j = kmin, n
      gn(j) = y(indx(j)) / f(xn(j), m)
    end do

    do k = 1, n
      j0 = min(max(k - 1, kmin), n - np + 1)
      t(1:np) = xn(j0:j0+np-1) - xn(k)
      ! Newton divided differences, then expand into powers of x - xn(k)
      dd(1:np) = gn(j0:j0+np-1)
      do j = 2, np
        do i = np, j, -1
          dd(i) = (dd(i) - dd(i-1)) / (t(i) - t(i-j+1))
        end do
      end do
      cf = 0.0_DP
      cf(1) = dd(np)
      do j = np - 1, 1, -1
        do i = np, 2, -1
          cf(i) = cf(i-1) - t(j) * cf(i)
        end do
        cf(1) = dd(j) - t(j) * cf(1)
      end do
      a(k) = cf(1)
      b(k) = cf(2)
      c(k) = cf(3)
      d(k) = cf(4)
    end do
  end subroutine lagrange_cof3

  subroutine radial_cof3(x, y, c1, cn, lambda, indx, sw1, sw2, a, b, c, d, m, f)
    ! splinecof3, or lagrange_cof3 if lsw_lagrange_boozer.
    use inter_interfaces, only : splinecof3
    real(DP), dimension(:), intent(in) :: x, y, lambda
    real(DP), intent(inout) :: c1, cn
    integer(I4B), dimension(:), intent(in) :: indx
    integer(I4B), intent(in) :: sw1, sw2
    real(DP), dimension(:), intent(out) :: a, b, c, d
    real(DP), intent(in) :: m
    procedure(lead_function) :: f

    if (lsw_lagrange_boozer) then
      call lagrange_cof3(x, y, m, a, b, c, d, indx, f)
    else
      call splinecof3(x, y, c1, cn, lambda, indx, sw1, sw2, a, b, c, d, m, f)
    end if
  end subroutine radial_cof3

  subroutine radial_cof3_hi_driv(x, y, m, a, b, c, d, indx, f)
    ! splinecof3_hi_driv, or column-wise lagrange_cof3 if lsw_lagrange_boozer.
    use inter_interfaces, only : splinecof3_hi_driv
    real(DP), dimension(:), intent(in) :: x, m
    real(DP), dimension(:,:), intent(in) :: y
    real(DP), dimension(:,:), intent(out) :: a, b, c, d
    integer(I4B), dimension(:), intent(in) :: indx
    procedure(lead_function) :: f

    integer(I4B) :: i

    if (lsw_lagrange_boozer) then
      do i = 1, size(y, 2)
        call lagrange_cof3(x, y(:,i), m(i), a(:,i), b(:,i), c(:,i), d(:,i), &
          indx, f)
      end do
    else
      call splinecof3_hi_driv(x, y, m, a, b, c, d, indx, f)
    end if
  end subroutine radial_cof3_hi_driv

end module radial_lagrange_cof
