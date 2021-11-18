module collop_bspline
  !use fgsl
  use gsl_bspline_routines_mod

  implicit none
  
  REAL(fgsl_double), DIMENSION(:), ALLOCATABLE, target :: xknots
  real(fgsl_double), dimension(:), allocatable :: b, db
  real(fgsl_double) :: xd
  integer :: i, j, knots
  integer :: nder
  logical, private :: bspline_initialized = .false.
  
contains
  subroutine check()
    use collisionality_mod, only : collop_bspline_order
    !use collop, only : collop_base_exp, collop_base_prj

    ! This check leads to circular dependencies and is thus not used until
    ! it was found out how to solve the problem.
    !if ((collop_base_prj .eq. 11) .or. (collop_base_exp .eq. 11)) then
      if (collop_bspline_order < 2) then
        stop 'ERROR: collop_bspline_order should be >= 2'
      end if
    !end if
  end subroutine check

  !**********************************************************
  ! Splines
  !**********************************************************
  subroutine init_phi_bspline(lagmax, legmax)
    use nrtype, only : dp
    use collisionality_mod, only : phi_x_max, collop_bspline_dist, collop_bspline_order, &
      &  collop_bspline_taylor

    integer :: lagmax, legmax
    integer :: k
    real(kind=dp) :: xp, gam_all, x_del

    taylorExpansion = collop_bspline_taylor

    if (.not. bspline_initialized) then
       bspline_initialized = .true.

       knots = lagmax - collop_bspline_order + 3
       nder  = collop_bspline_order  ! Maximum derivative

       if (.not. allocated(xknots))  allocate(xknots(knots))
       if (.not. allocated(b))  allocate(b(knots+collop_bspline_order-2))
       if (.not. allocated(db)) allocate(db(knots+collop_bspline_order-2))

       gam_all = 0d0
       do k = 1, knots-1
          gam_all = gam_all + collop_bspline_dist**k
       end do
       
       xknots(1) = 0d0
       x_del = phi_x_max / gam_all
       do k = 1, knots-1
          !x_dat(k) = k*x_del
          xknots(k+1) = xknots(k) + x_del * collop_bspline_dist**k
       end do
      
       xknots(knots) = phi_x_max

       write (*,*) "Knots: ", xknots
       call init_bspline(collop_bspline_order, knots)
       call set_bspline_knots(xknots)

       ! Testing

       !write (*,*) phi_bspline(2, 4d0)
       !write (*,*) d_phi_bspline(2, 4d0)
       !write (*,*) dd_phi_bspline(2, 4d0)

       !stop
      call check()
    end if
  end subroutine init_phi_bspline

  function phi_bspline(m, x) result(phi)
    use nrtype, only : dp

    integer       :: m
    real(kind=dp) :: x, phi

    call bspline_eval(x, b)
    !call bspline_deriv_eval(x, 0, b)
    phi = b(m+1)
    
  end function phi_bspline

  function d_phi_bspline(m, x) result(d_phi)
    use nrtype, only : dp 

    integer       :: m
    real(kind=dp) :: x, d_phi

    call bspline_deriv_eval(x, 1, db)
    d_phi = db(m+1)

  end function d_phi_bspline

  function dd_phi_bspline(m, x) result(dd_phi)
    use nrtype, only : dp

    integer       :: m
    real(kind=dp) :: x, dd_phi
    
    call bspline_deriv_eval(x, 2, db)
    dd_phi = db(m+1)
    
  end function dd_phi_bspline

end module collop_bspline
