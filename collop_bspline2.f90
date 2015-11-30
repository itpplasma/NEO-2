module collop_bspline2
  use nrtype, only : dp, pi
  use collisionality_mod, only : phi_x_max
  !use fgsl
  use gsl_bspline_routines_mod

  implicit none
  
  real(fgsl_double), dimension(:), allocatable :: x, b
  real(fgsl_double), dimension(:,:), allocatable :: db
  real(fgsl_double) :: xd
  integer :: i, j, order, knots
  integer :: nder
  logical, private :: bspline_initialized = .false.
  
contains

  !**********************************************************
  ! Splines
  !**********************************************************
  subroutine init_phi_bspline2(lagmax, legmax)
    integer :: lagmax, legmax
    real(kind=dp) :: xp
    !integer :: order, knots

    if (.not. bspline_initialized) then
       bspline_initialized = .true.
       order = 3
       knots = lagmax - order + 3
       nder  = order  ! Maximum derivative

       if (.not. allocated(x))  allocate(x(knots))
       if (.not. allocated(b))  allocate(b(knots+order-2))
       if (.not. allocated(db)) allocate(db(nder+1, knots+order-2))

       do i = 1, knots
          x(i) = phi_x_max**2/(knots-1) * (i-1)
       end do
       x(knots) = phi_x_max**2
       write (*,*) "Knots: ", x
       call init_bspline(order, knots)
       call set_bspline_knots(x)
       
    end if
  end subroutine init_phi_bspline2

  function phi_bspline2(m, x) result(phi)
    integer       :: m
    real(kind=dp) :: x, phi

    b = 0d0
    if (x**2 .le. phi_x_max**2) then
       call bspline_eval(x**2, b)
       phi = b(m+1)
    else
       ! Taylor expansion
       call bspline_eval(phi_x_max**2, b)
       call bspline_deriv_eval(phi_x_max**2, nder, db)
       phi = b(m+1) + db(2, m+1)*(x**2-phi_x_max**2)+db(3,m+1)/2*(x**2-phi_x_max**2)**2! + &
             !db(4,m+1)/6*(x-phi_x_max)**3 + db(5,m+1)/24*(x-phi_x_max)**4
    end if
    !write (*,*) x, phi
    !call bspline_deriv_eval(xd, nder, db)
  end function phi_bspline2

  function d_phi_bspline2(m, x) result(d_phi)
    integer       :: m
    real(kind=dp) :: x, d_phi

    db = 0d0
    if (x**2 .le. phi_x_max**2) then
       call bspline_deriv_eval(x**2, nder, db)
       d_phi = db(2, m+1)
    else
       call bspline_deriv_eval(phi_x_max**2, nder, db)
       d_phi = db(2, m+1) + db(3,m+1)*(x**2-phi_x_max**2)! + &
             !db(4,m+1)/2*(x-phi_x_max)**2 + db(5,m+1)/6*(x-phi_x_max)**3
    end if
  end function d_phi_bspline2

  function dd_phi_bspline2(m, x) result(dd_phi)
    integer       :: m
    real(kind=dp) :: x, dd_phi

    db = 0d0
    if (x**2 .le. phi_x_max**2) then
       call bspline_deriv_eval(x**2, nder, db)
       dd_phi = db(3, m+1)
    else
       call bspline_deriv_eval(phi_x_max**2, nder, db)
       dd_phi = db(3,m+1)! + db(4,m+1)*(x-phi_x_max) + db(5,m+1)/2*(x-phi_x_max)**2 
    end if
    !write (*,*) dd_phi
  end function dd_phi_bspline2

end module collop_bspline2
