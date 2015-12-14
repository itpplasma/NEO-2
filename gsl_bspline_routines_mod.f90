!-----------------------------------------------------------------------------------------!
! module: gsl_bspline_routines_mod                                                        !
! authors: Gernot Kapper, Andreas F. Martitsch                                            !
! date: 02.11.2015                                                                        !
! version: 0.1                                                                            !
!-----------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------!
! History                                                                                 !
! 0.1 - Initial version                                                                   !
!-----------------------------------------------------------------------------------------!

module gsl_bspline_routines_mod
  use fgsl

  implicit none

  type(fgsl_bspline_workspace), private :: sw
  type(fgsl_bspline_deriv_workspace), private :: dw
  integer(fgsl_size_t), private :: fgsl_order, fgsl_nbreak
  integer(fgsl_size_t), private :: fgsl_ncbf
  type(fgsl_vector), private :: xd_fgsl
  integer(fgsl_int), private :: status
  type(fgsl_vector), private :: b
  type(fgsl_matrix), private :: db
  real(fgsl_double), dimension(:), allocatable, private :: bv
  real(kind=kind(1d0)) :: xd_beg, xd_end
  integer :: max_nder
  real(kind=kind(1d0)), dimension(:,:), allocatable :: rdb
  
  logical :: taylorExpansion = .true.

contains

  subroutine init_bspline(order, nbreak)
    integer :: order, nbreak

    fgsl_order  = order
    fgsl_nbreak = nbreak
    fgsl_ncbf   = nbreak + order - 2
    max_nder    = order
    
    sw = fgsl_bspline_alloc(fgsl_order, fgsl_nbreak)
    dw = fgsl_bspline_deriv_alloc(fgsl_order)

    if (allocated(rdb)) deallocate(rdb)
    allocate(rdb(0:max_nder, 1:nbreak+order-2))
    
  end subroutine init_bspline

  subroutine set_bspline_knots(xd)
    real(fgsl_double), dimension(:) :: xd

    xd_beg = xd(1)
    xd_end = xd(fgsl_nbreak)
    
    xd_fgsl = fgsl_vector_init(1.0_fgsl_double)
    status = fgsl_vector_align(xd, fgsl_nbreak, xd_fgsl, fgsl_nbreak, 0_fgsl_size_t, 1_fgsl_size_t)
    status = fgsl_bspline_knots(xd_fgsl, sw)

    b  = fgsl_vector_init(1.0_fgsl_double)

    if (allocated(bv)) deallocate(bv)
    allocate(bv(fgsl_ncbf))
    status = fgsl_vector_align(bv, fgsl_ncbf, b, fgsl_ncbf, 0_fgsl_size_t, 1_fgsl_size_t)
    !status = fgsl_vector_align(bp, b)

    write (*,*) "ncoefs:"
    write (*,*) fgsl_bspline_ncoeffs(sw)
    
  end subroutine set_bspline_knots

  subroutine set_bspline_knots_uniform(a, b)
    real(fgsl_double) :: a, b
    status = fgsl_bspline_knots_uniform(a, b, sw)
    
  end subroutine set_bspline_knots_uniform

  recursive subroutine bspline_eval(x, by) 
    real(fgsl_double) :: x
    real(fgsl_double), dimension(:) :: by
    integer :: k, m
    
    if (x .le. xd_end) then
       status = fgsl_bspline_eval(x, b, sw)
       by = b
    elseif (taylorExpansion) then
       !status = fgsl_bspline_eval(x, b, sw)

       call bspline_eval(xd_end, by)
       call bspline_deriv_eval(xd_end, max_nder, rdb)
       
       do m = 1, fgsl_ncbf
          do k = 1, fgsl_order
             by(m) = by(m) + rdb(k, m)/gamma(k+1d0) * (x - xd_end)**k
          end do
       end do
    else
       write (*,*) "Stop, BSpline not defined outside interval."
       stop
    end if
    
  end subroutine bspline_eval

  recursive subroutine bspline_deriv_eval(x, nder, dby)
    real(fgsl_double) :: x
    integer(fgsl_int) :: nder
    integer(fgsl_size_t) :: fgsl_nder
    real(fgsl_double), dimension(:,:) :: dby
    integer :: n, m, k

    if (x .le. xd_end) then
       fgsl_nder = nder
       db = fgsl_matrix_init(1.0_fgsl_double)
       status = fgsl_matrix_align(dby, fgsl_nder+1, fgsl_nder+1, fgsl_ncbf, db)
       status = fgsl_bspline_deriv_eval(x, fgsl_nder, db, sw, dw)
       call fgsl_matrix_free(db)
    elseif (taylorExpansion) then
       call bspline_deriv_eval(xd_end, max_nder, rdb)
       dby = rdb
       do n = 1, nder
          do m = 1, fgsl_ncbf
             do k = 1, fgsl_order - (n-1)
                dby(n,m) = dby(n,m) + dby(k+(n), m)/gamma(k+1d0) * (x - xd_end)**k
             end do
          end do
       end do
    else
       write (*,*) "Stop, BSpline not defined outside interval."
       stop            
    end if
    
  end subroutine bspline_deriv_eval
  
end module gsl_bspline_routines_mod
