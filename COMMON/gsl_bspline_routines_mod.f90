!-----------------------------------------------------------------------------------------!
! module: gsl_bspline_routines_mod                                                        !
! authors: TU Graz ITPcp Plasma, Andreas F. Martitsch                                     !
! date: 06.10.2016                                                                        !
! version: 0.3.0                                                                          !
!-----------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------!
! History                                                                                 !
! 0.1 - Initial version
! 0.2 - Redesign
! 0.3 - Backend replaced by fortnum_bspline; public interface unchanged
!-----------------------------------------------------------------------------------------!

module gsl_bspline_routines_mod
  use, intrinsic :: iso_fortran_env, only: fgsl_double => real64
  use fortnum_bspline, only: bspline_workspace_t, bspline_init, bspline_set_knots, &
    &  bspline_eval_basis, bspline_eval_deriv
  use fortnum_status, only: fortnum_status_t, FORTNUM_OK

  implicit none

  !**********************************************************
  ! B-spline workspace (caller-owned, replaces FGSL handle)
  !**********************************************************
  type(bspline_workspace_t), private :: sw

  !**********************************************************
  ! Spline geometry
  !**********************************************************
  integer, private :: bspline_order, bspline_nbreak
  integer, private :: bspline_ncbf
  integer, private :: bspline_max_nder

  !**********************************************************
  ! Fortran variables
  !**********************************************************
  real(fgsl_double)                                  :: xd_beg, xd_end
  real(fgsl_double), dimension(:,:), allocatable, target :: dby
  real(fgsl_double), dimension(:), allocatable, target, private :: by
  real(fgsl_double), dimension(:), allocatable, target :: by_end
  real(fgsl_double), dimension(:,:), allocatable , target :: dby_end

  logical :: taylorExpansion = .true.

contains

  subroutine fortnum_check(status)
    type(fortnum_status_t) :: status

    if (status%code .ne. FORTNUM_OK) then
       write (*,*) "fortnum B-spline error: ", status%code, trim(status%msg)
       stop
    end if

  end subroutine fortnum_check

  subroutine init_bspline(order, nbreak)
    integer :: order, nbreak
    type(fortnum_status_t) :: status

    bspline_order    = order
    bspline_nbreak   = nbreak
    bspline_ncbf     = nbreak + order - 2
    bspline_max_nder = order

    call bspline_init(sw, order, nbreak, status)
    call fortnum_check(status)

    if (allocated(dby)) deallocate(dby)
    allocate(dby(0:bspline_max_nder, 1:bspline_ncbf))

  end subroutine init_bspline

  subroutine set_bspline_knots(xd)
    real(fgsl_double), dimension(:), target :: xd
    type(fortnum_status_t) :: status

    !**********************************************************
    ! Initialize knots
    !**********************************************************
    xd_beg = xd(1)
    xd_end = xd(bspline_nbreak)

    call bspline_set_knots(sw, xd(1:bspline_nbreak), status)
    call fortnum_check(status)

    if (allocated(by)) deallocate(by)
    allocate(by(bspline_ncbf))

    !**********************************************************
    ! Initialize Taylor expansion coefficients
    !**********************************************************
    if (allocated(by_end)) deallocate(by_end)
    allocate(by_end(bspline_ncbf))
    call bspline_eval_basis(sw, xd_end, by_end, status)
    call fortnum_check(status)

    if (allocated(dby_end)) deallocate(dby_end)
    allocate(dby_end(0:bspline_max_nder, 1:bspline_ncbf))
    call bspline_eval_deriv(sw, xd_end, bspline_max_nder, dby_end, status)
    call fortnum_check(status)

  end subroutine set_bspline_knots

  recursive subroutine bspline_eval(x, by_loc)
    real(fgsl_double) :: x
    real(fgsl_double), dimension(:) :: by_loc
    type(fortnum_status_t) :: status
    integer :: k, m
    real(fgsl_double) :: gam

    if (x .le. xd_end) then
       call bspline_eval_basis(sw, x, by, status)
       call fortnum_check(status)
       by_loc = by
    elseif (taylorExpansion) then
       by_loc = by_end
       do m = 1, bspline_ncbf
          gam = 1d0
          do k = 1, bspline_order-1
             gam = k * gam
             by_loc(m) = by_loc(m) + dby_end(k, m)/gam * (x - xd_end)**k
          end do
       end do
    else
       by_loc = 0.0_fgsl_double
    end if

  end subroutine bspline_eval

  recursive subroutine bspline_deriv_eval(x, nder, dby_loc)
    real(fgsl_double) :: x
    integer :: nder
    real(fgsl_double), dimension(:) :: dby_loc
    type(fortnum_status_t) :: status
    integer :: m, k
    real(fgsl_double) :: gam

    if (x .le. xd_end) then
       call bspline_eval_deriv(sw, x, bspline_max_nder, dby, status)
       call fortnum_check(status)
       dby_loc = dby(nder, :)
    elseif (taylorExpansion) then
       dby_loc = dby_end(nder, :)
       do m = 1, bspline_ncbf
          gam = 1d0
          do k = 1, bspline_order - nder - 1
             gam = k * gam
             dby_loc(m) = dby_loc(m) + dby_end(k+nder, m)/gam * (x - xd_end)**k
          end do
       end do
    else
       dby_loc = 0.0_fgsl_double
    end if

  end subroutine bspline_deriv_eval

end module gsl_bspline_routines_mod
