!-----------------------------------------------------------------------------------------!
! module: gsl_bspline_routines_mod                                                        !
! authors: Gernot Kapper, Andreas F. Martitsch                                            !
! date: 06.10.2016                                                                        !
! version: 0.2.1                                                                          !
!-----------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------!
! History                                                                                 !
! 0.1 - Initial version
! 0.2 - Redesign
!-----------------------------------------------------------------------------------------!

module gsl_bspline_routines_mod
  use fgsl

  implicit none

  !**********************************************************
  ! FGSL Workspace
  !**********************************************************
  type(fgsl_bspline_workspace), private :: sw
  type(fgsl_bspline_deriv_workspace), private :: dw

  !**********************************************************
  ! FGSL variables
  !**********************************************************
  integer(fgsl_size_t), private :: fgsl_order, fgsl_nbreak
  integer(fgsl_size_t), private :: fgsl_ncbf
  type(fgsl_vector), private :: fgsl_xd
  type(fgsl_vector), private :: fgsl_by, fgsl_by_end
  type(fgsl_matrix), private :: fgsl_dby, fgsl_dby_end
  integer(fgsl_size_t)       :: fgsl_max_nder

  !**********************************************************
  ! Fortran variables
  !**********************************************************
  real(kind=kind(1d0))                                  :: xd_beg, xd_end
  real(fgsl_double), dimension(:,:), allocatable        :: dby
  real(fgsl_double), dimension(:), allocatable, private :: by
  real(kind=kind(1d0)), dimension(:), allocatable       :: by_end
  real(kind=kind(1d0)), dimension(:,:), allocatable     :: dby_end

  logical :: taylorExpansion = .true.
  !logical :: taylorExpansion = .false.

contains

  subroutine fgsl_check(status)
    integer(fgsl_int) :: status

    if (status .ne. 0) then
       write (*,*) "FGSL Error: ", status, fgsl_strerror(status)
       stop
    end if
    
  end subroutine fgsl_check
  
  subroutine init_bspline(order, nbreak)
    integer :: order, nbreak

    fgsl_order    = order
    fgsl_nbreak   = nbreak
    fgsl_ncbf     = nbreak + order - 2
    fgsl_max_nder = order
    
    sw = fgsl_bspline_alloc(fgsl_order, fgsl_nbreak)
    dw = fgsl_bspline_deriv_alloc(fgsl_order)

    if (allocated(dby)) deallocate(dby)
    allocate(dby(1:fgsl_max_nder+1, 1:nbreak+order-2))
    fgsl_dby = fgsl_matrix_init(1.0_fgsl_double)
    call fgsl_check(fgsl_matrix_align(dby, fgsl_max_nder+1, fgsl_max_nder+1, fgsl_ncbf, fgsl_dby))
    
  end subroutine init_bspline

  subroutine set_bspline_knots(xd)
    real(fgsl_double), dimension(:) :: xd

    !**********************************************************
    ! Initialize knots
    !**********************************************************
    xd_beg = xd(1)
    xd_end = xd(fgsl_nbreak)
    
    fgsl_xd = fgsl_vector_init(1.0_fgsl_double)
    call fgsl_check(fgsl_vector_align(xd, fgsl_nbreak, fgsl_xd, fgsl_nbreak, 0_fgsl_size_t, 1_fgsl_size_t))
    call fgsl_check(fgsl_bspline_knots(fgsl_xd, sw))

    fgsl_by  = fgsl_vector_init(1.0_fgsl_double)

    if (allocated(by)) deallocate(by)
    allocate(by(fgsl_ncbf))
    call fgsl_check(fgsl_vector_align(by, fgsl_ncbf, fgsl_by, fgsl_ncbf, 0_fgsl_size_t, 1_fgsl_size_t))

    !write (*,*) "ncoefs:"
    !write (*,*) fgsl_bspline_ncoeffs(sw)

    !**********************************************************
    ! Initialize Taylor expansion coefficients
    !**********************************************************
    if (allocated(by_end)) deallocate(by_end)
    allocate(by_end(fgsl_ncbf))
    fgsl_by_end = fgsl_vector_init(1.0_fgsl_double)
    call fgsl_check(fgsl_vector_align(by_end, fgsl_ncbf, fgsl_by_end, fgsl_ncbf, 0_fgsl_size_t, 1_fgsl_size_t))
    call fgsl_check(fgsl_bspline_eval(xd_end, fgsl_by_end, sw))
        
    if (allocated(dby_end)) deallocate(dby_end)
    allocate(dby_end(0:fgsl_max_nder, 1:fgsl_ncbf))
    fgsl_dby_end = fgsl_matrix_init(1.0_fgsl_double)
    call fgsl_check(fgsl_matrix_align(dby_end, fgsl_max_nder+1, fgsl_max_nder+1, fgsl_ncbf, fgsl_dby_end))
    call fgsl_check(fgsl_bspline_deriv_eval(xd_end, fgsl_max_nder, fgsl_dby_end, sw, dw))
        
  end subroutine set_bspline_knots

  subroutine set_bspline_knots_uniform(a, b)
    real(fgsl_double) :: a, b
    call fgsl_check(fgsl_bspline_knots_uniform(a, b, sw))
    
  end subroutine set_bspline_knots_uniform

  recursive subroutine bspline_eval(x, by_loc) 
    real(fgsl_double) :: x
    real(fgsl_double), dimension(:) :: by_loc
    integer(fgsl_size_t) :: k, m
    real(fgsl_double)    :: gam
    
    if (x .le. xd_end) then
       call fgsl_check(fgsl_bspline_eval(x, fgsl_by, sw))
       by_loc = fgsl_by
    elseif (taylorExpansion) then
       !write (*,*) "Taylor Expansion of BSplines", x
       by_loc = by_end
       do m = 1, fgsl_ncbf
          gam = 1d0
          do k = 1, fgsl_order-1
             gam = k * gam
             by_loc(m) = by_loc(m) + dby_end(k, m)/gam * (x - xd_end)**k
          end do
       end do
    else
       by_loc = 0.0_fgsl_double
       !write (*,*) "Stop, BSpline not defined outside interval."
       !stop
    end if
    
  end subroutine bspline_eval

  recursive subroutine bspline_deriv_eval(x, nder, dby_loc)
    real(fgsl_double) :: x
    integer(fgsl_int) :: nder
    integer(fgsl_size_t) :: fgsl_nder
    real(fgsl_double), dimension(:) :: dby_loc
    integer(fgsl_size_t) :: m, k
    real(fgsl_double) :: gam

    if (x .le. xd_end) then
       fgsl_nder = nder
       call fgsl_check(fgsl_bspline_deriv_eval(x, fgsl_nder, fgsl_dby, sw, dw))
       dby_loc = dby(nder+1, :)
    elseif (taylorExpansion) then
       !write (*,*) "Taylor Expansion of BSplines (deriv)", nder, x
       dby_loc = dby_end(nder, :)
       do m = 1, fgsl_ncbf
          gam = 1d0
          do k = 1, fgsl_order - nder - 1
             gam = k * gam
             dby_loc(m) = dby_loc(m) + dby_end(k+nder, m)/gam * (x - xd_end)**k
          end do
       end do
    else
       dby_loc = 0.0_fgsl_double
       !write (*,*) "Stop, BSpline not defined outside interval."
       !stop            
    end if
    
  end subroutine bspline_deriv_eval
  
end module gsl_bspline_routines_mod
