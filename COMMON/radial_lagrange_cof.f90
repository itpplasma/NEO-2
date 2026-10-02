module radial_lagrange_cof
  ! Local Lagrange alternative to the radial cubic splines of the Boozer
  ! data (neo_init_spline, neo_prep_b00, neo_init_spline_pert).
  !
  ! The splines cost a banded solve over all ns surfaces per Fourier mode,
  ! although NEO-2 evaluates them at one flux surface. Here each interval
  ! [x(k), x(k+1)] gets the cubic through the four nodes k-1..k+2 (shifted
  ! inwards at the boundaries). The result is stored in the same a, b, c, d
  ! arrays that splint_horner3 evaluates, y = f(x,m) * p(x - x(k)), so all
  ! evaluation code stays unchanged. Values and first radial derivatives
  ! converge as h^4 and h^3, like the cubic spline; the first derivative is
  ! discontinuous across nodes, which does not matter at a fixed surface.
  use nrtype, only : I4B, DP, splinecof_compatibility

  implicit none

  private
  public :: lsw_lagrange_boozer, lagrange_cof3, lagrange_cof3_hi_driv
  public :: radial_cof3, radial_cof3_hi_driv

  ! Use local Lagrange polynomials instead of radial cubic splines
  ! (namelist switch, default off).
  logical, save :: lsw_lagrange_boozer = .false.

  integer(I4B), parameter :: npoly = 4
    ! Opt-in transport convergence has been checked for radial spacings
    ! 0.001 and 0.002 in normalized toroidal flux (k differences < 1e-9).
    ! The wider 0.004 and 0.008 studies lose accuracy; retain splines there.
    real(DP), parameter :: max_lagrange_spacing = 0.002_DP
    logical, save :: fallback_reported = .false.

  abstract interface
    function test_function(x, m)
      import :: DP
      real(DP), intent(in) :: x
      real(DP), intent(in) :: m
      real(DP) :: test_function
    end function test_function
  end interface

contains

    logical function use_local_lagrange(x, indx) result(admitted)
        real(DP), intent(in) :: x(:)
        integer(I4B), intent(in) :: indx(:)
        integer(I4B) :: j
        real(DP) :: spacing

        admitted = .false.
        if (.not. lsw_lagrange_boozer) return
        if (any(indx < 1) .or. any(indx > size(x))) return
        ! Weighted modes exclude the first node, so five nodes are needed
        ! to retain a cubic stencil. Check every stencil conservatively.
        if (size(indx) >= npoly + 1) then
            admitted = .true.
            do j = 1, size(indx) - 1
                spacing = x(indx(j + 1)) - x(indx(j))
                if (spacing <= 0.0_DP .or. spacing > max_lagrange_spacing &
                    + 64.0_DP * epsilon(1.0_DP)) then
                    admitted = .false.
                    exit
                end if
            end do
        end if
        if (.not. admitted) then
            if (.not. fallback_reported) then
                print *, 'LSW_LAGRANGE_BOOZER: grid outside validated spacing; ', &
                    'using splines (maximum supported normalized ds = 0.002)'
                fallback_reported = .true.
            end if
        end if
    end function use_local_lagrange

  subroutine lagrange_cof3(x, y, m, a, b, c, d, indx, f)
    ! Piecewise cubic Lagrange coefficients of y(x) / f(x, m) at the nodes
    ! x(indx). As in splinecof3_hi_driv, the first node is skipped for
    ! m /= 0 because y/f is undefined at the axis.
    real(DP), dimension(:), intent(in) :: x, y
    real(DP), intent(in) :: m
    real(DP), dimension(:), intent(out) :: a, b, c, d
    integer(I4B), dimension(:), intent(in) :: indx
    procedure(test_function) :: f

    real(DP), dimension(size(indx)) :: xn, gn
    real(DP), dimension(npoly) :: t, dd, cf
    integer(I4B) :: n, k, j, i, kmin, j0, np

    if (splinecof_compatibility) &
      error stop 'lagrange_cof3: splinecof_compatibility is not supported'

    n = size(indx)
    if (size(a) /= n .or. size(b) /= n .or. size(c) /= n .or. size(d) /= n) &
      error stop 'lagrange_cof3: size mismatch'

    kmin = 1
    if (m /= 0.0_DP) kmin = 2
    np = min(npoly, n - kmin + 1)
    if (np < 2) error stop 'lagrange_cof3: too few nodes'

    xn = x(indx)
    do j = kmin, n
      gn(j) = y(indx(j)) / f(xn(j), m)
    end do
    if (kmin == 2) gn(1) = 0.0_DP

    do k = 1, n
      j0 = min(max(k - 1, kmin), n - np + 1)
      t(1:np) = xn(j0:j0+np-1) - xn(k)
      ! Newton divided differences
      dd(1:np) = gn(j0:j0+np-1)
      do j = 2, np
        do i = np, j, -1
          dd(i) = (dd(i) - dd(i-1)) / (t(i) - t(i-j+1))
        end do
      end do
      ! Expand the Newton form into monomials in h = x - xn(k)
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

  subroutine lagrange_cof3_hi_driv(x, y, m, a, b, c, d, indx, f)
    ! Column-wise lagrange_cof3, same interface as splinecof3_hi_driv.
    real(DP), dimension(:), intent(in) :: x
    real(DP), dimension(:,:), intent(in) :: y
    real(DP), dimension(:), intent(in) :: m
    real(DP), dimension(:,:), intent(out) :: a, b, c, d
    integer(I4B), dimension(:), intent(in) :: indx
    procedure(test_function) :: f

    integer(I4B) :: i

    do i = 1, size(y, 2)
      call lagrange_cof3(x, y(:,i), m(i), a(:,i), b(:,i), c(:,i), d(:,i), &
        indx, f)
    end do
  end subroutine lagrange_cof3_hi_driv

  subroutine radial_cof3(x, y, a, b, c, d, indx, f)
    ! Radial profile without smoothing and natural boundary conditions
    ! (as in neo_init_spline), or local Lagrange if lsw_lagrange_boozer.
    use inter_interfaces, only : splinecof3
    real(DP), dimension(:), intent(in) :: x, y
    real(DP), dimension(:), intent(out) :: a, b, c, d
    integer(I4B), dimension(:), intent(in) :: indx
    procedure(test_function) :: f

    real(DP), dimension(size(indx)) :: lambda
    real(DP), parameter :: m0 = 0.0_DP
    real(DP) :: c1, cn
    integer(I4B), parameter :: sw1 = 2, sw2 = 4

    if (use_local_lagrange(x, indx)) then
      call lagrange_cof3(x, y, m0, a, b, c, d, indx, f)
    else
      lambda = 1.0_DP
      c1 = 0.0_DP
      cn = 0.0_DP
      call splinecof3(x, y, c1, cn, lambda, indx, sw1, sw2, a, b, c, d, m0, f)
    end if
  end subroutine radial_cof3

  subroutine radial_cof3_hi_driv(x, y, m, a, b, c, d, indx, f)
    ! splinecof3_hi_driv, or local Lagrange if lsw_lagrange_boozer.
    use inter_interfaces, only : splinecof3_hi_driv
    real(DP), dimension(:), intent(in) :: x
    real(DP), dimension(:,:), intent(in) :: y
    real(DP), dimension(:), intent(in) :: m
    real(DP), dimension(:,:), intent(out) :: a, b, c, d
    integer(I4B), dimension(:), intent(in) :: indx
    procedure(test_function) :: f

    if (use_local_lagrange(x, indx)) then
      call lagrange_cof3_hi_driv(x, y, m, a, b, c, d, indx, f)
    else
      call splinecof3_hi_driv(x, y, m, a, b, c, d, indx, f)
    end if
  end subroutine radial_cof3_hi_driv

end module radial_lagrange_cof
