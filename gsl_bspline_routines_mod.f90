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

contains

  subroutine init_bspline(order, nbreak)
    integer :: order, nbreak

    fgsl_order  = order
    fgsl_nbreak = nbreak
    fgsl_ncbf   = nbreak + order - 2
    
    sw = fgsl_bspline_alloc(fgsl_order, fgsl_nbreak)
    dw = fgsl_bspline_deriv_alloc(fgsl_order)
    
  end subroutine init_bspline

  subroutine set_bspline_knots(xd)
    real(fgsl_double), dimension(:) :: xd
    
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

  subroutine bspline_eval(x, by)
    real(fgsl_double) :: x
    real(fgsl_double), dimension(:) :: by 

    status = fgsl_bspline_eval(x, b, sw)
    by = b
    
  end subroutine bspline_eval

  subroutine bspline_deriv_eval(x, nder, dby)
    real(fgsl_double) :: x
    integer(fgsl_int) :: nder
    integer(fgsl_size_t) :: fgsl_nder
    real(fgsl_double), dimension(:,:) :: dby

    fgsl_nder = nder
    db = fgsl_matrix_init(1.0_fgsl_double)
    status = fgsl_matrix_align(dby, fgsl_nder+1, fgsl_nder+1, fgsl_ncbf, db)
    
    status = fgsl_bspline_deriv_eval(x, fgsl_nder, db, sw, dw)
    
    call fgsl_matrix_free(db)
    
  end subroutine bspline_deriv_eval
  
end module gsl_bspline_routines_mod
