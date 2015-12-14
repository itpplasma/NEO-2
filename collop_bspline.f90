module collop_bspline
  use nrtype, only : dp, pi
  use collisionality_mod, only : phi_x_max, collop_bspline_dist, collop_bspline_order
  !use fgsl
  use gsl_bspline_routines_mod

  implicit none
  
  real(fgsl_double), dimension(:), allocatable :: x, b
  real(fgsl_double), dimension(:,:), allocatable :: db
  real(fgsl_double) :: xd
  integer :: i, j, knots
  integer :: nder
  logical, private :: bspline_initialized = .false.
  
contains

  !**********************************************************
  ! Splines
  !**********************************************************
  subroutine init_phi_bspline(lagmax, legmax)
    integer :: lagmax, legmax
    integer :: k
    real(kind=dp) :: xp, gam_all, x_del

    if (.not. bspline_initialized) then
       bspline_initialized = .true.
       
       knots = lagmax - collop_bspline_order + 3
       nder  = collop_bspline_order  ! Maximum derivative

       if (.not. allocated(x))  allocate(x(knots))
       if (.not. allocated(b))  allocate(b(knots+collop_bspline_order-2))
       if (.not. allocated(db)) allocate(db(nder+1, knots+collop_bspline_order-2))

       gam_all = 0d0
       do k = 1, knots-1
          gam_all = gam_all + collop_bspline_dist**k
       end do

       x_del = phi_x_max / gam_all
       do k = 1, knots-1
          !x_dat(k) = k*x_del
          x(k+1) = x(k) + x_del * collop_bspline_dist**k
       end do
       x(1)     = 0d0
       x(knots) = phi_x_max

       write (*,*) "Knots: ", x
       call init_bspline(collop_bspline_order, knots)
       call set_bspline_knots(x)
       
    end if
  end subroutine init_phi_bspline

  function phi_bspline(m, x) result(phi)
    integer       :: m
    real(kind=dp) :: x, phi

    call bspline_eval(x, b)
    phi = b(m+1)
    
  end function phi_bspline

  function d_phi_bspline(m, x) result(d_phi)
    integer       :: m
    real(kind=dp) :: x, d_phi

    call bspline_deriv_eval(x, nder, db)
    d_phi = db(2, m+1)

  end function d_phi_bspline

  function dd_phi_bspline(m, x) result(dd_phi)
    integer       :: m
    real(kind=dp) :: x, dd_phi
    
    call bspline_deriv_eval(x, nder, db)
    dd_phi = db(3, m+1)
    
  end function dd_phi_bspline

end module collop_bspline
