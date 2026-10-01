!> Find the poloidal angle of a given (R,Z) point on a closed curve
!> (R(theta), Z(theta)), theta in [0, 2 pi).
!>
!> Used by neo_magfie::calc_thetaB_RZloc to locate a cylindrical point on a
!> flux surface. Solving R(theta) = R0 and Z(theta) = Z0 separately is not
!> well posed: on a closed curve each equation has two roots, so separate
!> iterations can converge to different angles. Instead the squared distance
!> D(theta) = (R - R0)^2 + (Z - Z0)^2 is minimised: a uniform scan selects the
!> global minimum, then g(theta) = dD/dtheta / 2 = (R - R0) R' + (Z - Z0) Z'
!> is driven to zero by a safeguarded secant/bisection iteration inside the
!> bracket around the scan minimum. Only R, Z and their theta derivatives are
!> needed.
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
  !> nscan    : optional number of scan points (default 360)
  subroutine find_theta_of_rz(rz, R0, Z0, theta, dist, extent, nscan)
    procedure(rz_of_theta)        :: rz
    real(dp), intent(in)          :: R0, Z0
    real(dp), intent(out)         :: theta, dist, extent
    integer, intent(in), optional :: nscan

    integer, parameter :: kmax = 200
    real(dp), parameter :: tol_theta = 1.0e-13_dp
    integer :: n, i, imin, k
    real(dp) :: dtheta, th, R, R_tb, Z, Z_tb, d2, d2min
    real(dp) :: Rmin, Rmax, Zmin, Zmax
    real(dp) :: a, b, ga, gb, c, gc

    n = 360
    if (present(nscan)) n = max(nscan, 8)
    dtheta = twopi / real(n, dp)

    imin = 0
    d2min = huge(1.0_dp)
    Rmin = huge(1.0_dp); Rmax = -huge(1.0_dp)
    Zmin = huge(1.0_dp); Zmax = -huge(1.0_dp)
    do i = 0, n - 1
       th = real(i, dp) * dtheta
       call rz(th, R, R_tb, Z, Z_tb)
       d2 = (R - R0)**2 + (Z - Z0)**2
       if (d2 < d2min) then
          d2min = d2
          imin = i
       end if
       Rmin = min(Rmin, R); Rmax = max(Rmax, R)
       Zmin = min(Zmin, Z); Zmax = max(Zmax, Z)
    end do
    extent = (Rmax - Rmin) + (Zmax - Zmin)

    ! Bracket [a, b] around the scan minimum; g < 0 left of a minimum of D
    ! and g > 0 right of it. Angles are kept unwrapped inside the bracket.
    a = real(imin - 1, dp) * dtheta
    b = real(imin + 1, dp) * dtheta
    ga = g_of(a)
    gb = g_of(b)
    theta = real(imin, dp) * dtheta
    if (ga <= 0.0_dp .and. gb >= 0.0_dp .and. ga < gb) then
       do k = 1, kmax
          ! secant step, replaced by bisection if it leaves the bracket
          c = b - gb * (b - a) / (gb - ga)
          if (.not. (c > a .and. c < b)) c = 0.5_dp * (a + b)
          gc = g_of(c)
          if (gc == 0.0_dp) then
             a = c; b = c
          else if (gc < 0.0_dp) then
             a = c; ga = gc
          else
             b = c; gb = gc
          end if
          theta = c
          if (b - a < tol_theta) exit
          ! guarantee bracket shrinkage when the secant stalls on one side
          if (mod(k, 3) == 0) then
             c = 0.5_dp * (a + b)
             gc = g_of(c)
             if (gc < 0.0_dp) then
                a = c; ga = gc
             else
                b = c; gb = gc
             end if
             theta = c
          end if
       end do
    end if

    theta = modulo(theta, twopi)
    call rz(theta, R, R_tb, Z, Z_tb)
    dist = sqrt((R - R0)**2 + (Z - Z0)**2)

  contains

    real(dp) function g_of(t)
      real(dp), intent(in) :: t
      real(dp) :: Rt, Rt_tb, Zt, Zt_tb
      call rz(modulo(t, twopi), Rt, Rt_tb, Zt, Zt_tb)
      g_of = (Rt - R0) * Rt_tb + (Zt - Z0) * Zt_tb
    end function g_of

  end subroutine find_theta_of_rz

end module theta_rz_solver_mod
