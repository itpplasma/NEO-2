!> Find the poloidal angle of a given (R,Z) point on a closed curve
!> (R(theta), Z(theta)), theta in [0, 2 pi).
!>
!> Used by neo_magfie::calc_thetaB_RZloc to locate a cylindrical point on a
!> flux surface. Solving R(theta) = R0 and Z(theta) = Z0 separately is not
!> well posed: on a closed curve each equation has two roots, so separate
!> iterations can converge to different angles. Instead the squared distance
!> D(theta) = (R - R0)^2 + (Z - Z0)^2 is minimised globally: the zeros of
!> g(theta) = dD/dtheta / 2 = (R - R0) R' + (Z - Z0) Z' at which g changes
!> sign from - to + are bracketed by a uniform scan and refined by a
!> safeguarded secant/bisection iteration; the closest one is kept, and the
!> scan is refined until the result no longer changes. Only R, Z and their
!> theta derivatives are needed.
module theta_rz_solver_mod
  use nrtype, only: dp, twopi
  implicit none
  private

  public :: find_theta_of_rz

  abstract interface
     subroutine rz_of_theta(theta, R, R_tb, Z, Z_tb)
       import :: dp
       real(dp), intent(in)  :: theta
       real(dp), intent(out) :: R, R_tb, Z, Z_tb
     end subroutine rz_of_theta
  end interface
  public :: rz_of_theta

contains

  !> theta    : angle in [0, 2 pi) of the curve point closest to (R0, Z0)
  !> dist     : distance between that point and (R0, Z0)
  !> extent   : (max R - min R) + (max Z - min Z) of the curve from the scan,
  !>            a length scale for judging dist
  !> nscan    : optional number of initial scan points (default 360)
  !>
  !> Every scan interval in which g changes sign from - to + (a local minimum
  !> of D) is refined, and the refined candidate with the smallest distance
  !> is kept. The scan resolution is doubled (at most 8 times) and the
  !> closest candidate over all resolutions is returned; the doubling stops
  !> early only once a point on the curve (distance <= 1e-10 extent) is found
  !> at two successive resolutions, so minima narrower than the scan step are
  !> not missed.
  subroutine find_theta_of_rz(rz, R0, Z0, theta, dist, extent, nscan)
    procedure(rz_of_theta)        :: rz
    real(dp), intent(in)          :: R0, Z0
    real(dp), intent(out)         :: theta, dist, extent
    integer, intent(in), optional :: nscan

    integer, parameter :: max_doublings = 8
    real(dp), parameter :: tol_same = 1.0e-9_dp, tol_on_curve = 1.0e-10_dp
    integer :: n, idouble
    real(dp) :: theta_n, dist_n, theta_prev

    n = 360
    if (present(nscan)) n = max(nscan, 8)

    call solve_at_resolution(rz, R0, Z0, n, theta, dist, extent)
    theta_prev = theta
    do idouble = 1, max_doublings
       n = 2 * n
       call solve_at_resolution(rz, R0, Z0, n, theta_n, dist_n, extent)
       if (dist_n < dist) then
          theta = theta_n
          dist = dist_n
       end if
       ! early exit only for a point on the curve found at two resolutions
       if (dist <= tol_on_curve * extent .and. &
            angle_distance(theta_n, theta_prev) < tol_same) exit
       theta_prev = theta_n
    end do
  end subroutine find_theta_of_rz

  subroutine solve_at_resolution(rz, R0, Z0, n, theta, dist, extent)
    procedure(rz_of_theta) :: rz
    real(dp), intent(in)   :: R0, Z0
    integer, intent(in)    :: n
    real(dp), intent(out)  :: theta, dist, extent

    integer :: i, ip
    real(dp) :: dtheta, R, R_tb, Z, Z_tb, d2min, tc, dc
    real(dp) :: Rmin, Rmax, Zmin, Zmax
    real(dp), allocatable :: g(:), d2(:)

    allocate(g(0:n-1), d2(0:n-1))
    dtheta = twopi / real(n, dp)
    Rmin = huge(1.0_dp); Rmax = -huge(1.0_dp)
    Zmin = huge(1.0_dp); Zmax = -huge(1.0_dp)
    do i = 0, n - 1
       call rz(real(i, dp) * dtheta, R, R_tb, Z, Z_tb)
       g(i) = (R - R0) * R_tb + (Z - Z0) * Z_tb
       d2(i) = (R - R0)**2 + (Z - Z0)**2
       Rmin = min(Rmin, R); Rmax = max(Rmax, R)
       Zmin = min(Zmin, Z); Zmax = max(Zmax, Z)
    end do
    extent = (Rmax - Rmin) + (Zmax - Zmin)

    ! fallback if no sign change is resolved: best sample
    i = minloc(d2, dim=1) - 1
    theta = real(i, dp) * dtheta
    d2min = d2(i)
    do i = 0, n - 1
       ip = modulo(i + 1, n)
       if (g(i) < 0.0_dp .and. g(ip) >= 0.0_dp) then
          call refine(rz, R0, Z0, real(i, dp) * dtheta, &
               real(i + 1, dp) * dtheta, g(i), g(ip), tc, dc)
          if (dc < d2min) then
             d2min = dc
             theta = tc
          end if
       end if
    end do

    theta = modulo(theta, twopi)
    if (theta >= twopi) theta = 0.0_dp
    call rz(theta, R, R_tb, Z, Z_tb)
    dist = sqrt((R - R0)**2 + (Z - Z0)**2)
  end subroutine solve_at_resolution

  !> Root of g in [a, b] with g(a) < 0 <= g(b): secant steps, replaced by
  !> bisection when they leave the bracket, plus a forced bisection every
  !> third step so that the bracket shrinks geometrically. Returns the root
  !> t and the squared distance d2 there.
  subroutine refine(rz, R0, Z0, a_in, b_in, ga_in, gb_in, t, d2)
    procedure(rz_of_theta) :: rz
    real(dp), intent(in)   :: R0, Z0, a_in, b_in, ga_in, gb_in
    real(dp), intent(out)  :: t, d2

    integer, parameter :: kmax = 200
    real(dp), parameter :: tol_theta = 1.0e-13_dp
    integer :: k
    real(dp) :: a, b, ga, gb, c, gc, R, R_tb, Z, Z_tb

    a = a_in; b = b_in; ga = ga_in; gb = gb_in
    t = b
    if (gb /= 0.0_dp) then
       do k = 1, kmax
          if (mod(k, 3) == 0) then
             c = 0.5_dp * (a + b)
          else
             c = b - gb * (b - a) / (gb - ga)
             if (.not. (c > a .and. c < b)) c = 0.5_dp * (a + b)
          end if
          call rz(modulo(c, twopi), R, R_tb, Z, Z_tb)
          gc = (R - R0) * R_tb + (Z - Z0) * Z_tb
          t = c
          if (gc == 0.0_dp) exit
          if (gc < 0.0_dp) then
             a = c; ga = gc
          else
             b = c; gb = gc
          end if
          if (b - a < tol_theta) exit
       end do
    end if
    call rz(modulo(t, twopi), R, R_tb, Z, Z_tb)
    d2 = (R - R0)**2 + (Z - Z0)**2
  end subroutine refine

  real(dp) function angle_distance(a, b)
    real(dp), intent(in) :: a, b
    angle_distance = abs(modulo(a - b + 0.5_dp * twopi, twopi) - 0.5_dp * twopi)
  end function angle_distance

end module theta_rz_solver_mod
