module collop_spline
  use nrtype, only : dp, pi
  use collisionality_mod, only : v_max_resolution
  use inter_interfaces, only : splinecof3, splint_horner3, tf, tfzero, tfone
  implicit none

  real(kind=dp), dimension(:),   allocatable, private :: x_dat
  real(kind=dp), dimension(:,:), allocatable, private :: y_dat

  real(kind=dp), dimension(:,:), allocatable, private :: a_spl, b_spl, c_spl, d_spl
  
contains

  !**********************************************************
  ! Splines
  !**********************************************************
  subroutine init_phi_spline(lagmax, legmax)
    integer :: lagmax, legmax
    integer :: k, m
    real(kind=dp) :: x_del
    integer :: swd, sw1, sw2
    real(kind=dp) :: m0, c1, cn
    real(kind=dp), dimension(:), allocatable :: lambda
    integer, dimension(:), allocatable :: sp_index

    if (allocated(x_dat)) deallocate(x_dat)
    if (allocated(y_dat)) deallocate(y_dat)
    
    allocate(x_dat(0:lagmax))
    allocate(y_dat(0:lagmax, 0:lagmax))
    x_dat = 0d0
    y_dat = 0d0
    x_del = v_max_resolution / lagmax
    do k = 0, lagmax
       x_dat(k) = k*x_del
       y_dat(k,k) = 1.0d0
    end do

    sw1 = 2
    sw2 = 4
    c1 = 0d0
    cn = 0d0
    m0 = 0d0
    allocate(lambda(lagmax+1))
    lambda = 1d0

    allocate(sp_index(lagmax+1))
    sp_index = (/ (k, k=1, lagmax+1) /)

    if (allocated(a_spl)) deallocate(a_spl)
    if (allocated(b_spl)) deallocate(b_spl)
    if (allocated(c_spl)) deallocate(c_spl)
    if (allocated(d_spl)) deallocate(d_spl)
    allocate(a_spl(0:lagmax, 0:lagmax))
    allocate(b_spl(0:lagmax, 0:lagmax))
    allocate(c_spl(0:lagmax, 0:lagmax))
    allocate(d_spl(0:lagmax, 0:lagmax))

    do m = 0, lagmax

       call splinecof3(x_dat, y_dat(m,:), c1, cn, lambda, sp_index, sw1, sw2, &
            a_spl(m,:), b_spl(m,:), c_spl(m,:), d_spl(m,:), m0, tf)

       !write (*,*) x_dat, y_dat(m,:), c1, cn, lambda, sp_index, sw1, sw2, a_spl(m,:), &
       !     b_spl(m,:), c_spl(m,:), d_spl(m,:), m0
       
    end do

    deallocate(sp_index)
    deallocate(lambda)
    
  end subroutine init_phi_spline

  function phi_spline(m, x) result(phi)
    integer       :: m
    real(kind=dp) :: x, phi
    integer       :: swd
    real(kind=dp) :: m0
    real(kind=dp) :: y, yp, ypp, yppp
    
    swd = 1
    m0 = 0d0
    call splint_horner3(x_dat, a_spl(m,:), b_spl(m,:), c_spl(m,:), d_spl(m,:), swd, m0, &
         x, tfone, tfzero, tfzero, tfzero, y, yp, ypp, yppp)
    phi = y
  end function phi_spline

  function d_phi_spline(m, x) result(d_phi)
    integer       :: m
    real(kind=dp) :: x, d_phi
    integer       :: swd
    real(kind=dp) :: m0
    real(kind=dp) :: y, yp, ypp, yppp
    
    swd = 1
    m0 = 0d0

    call splint_horner3(x_dat, a_spl(m,:), b_spl(m,:), c_spl(m,:), d_spl(m,:), swd, m0, &
         x, tfone, tfzero, tfzero, tfzero, y, yp, ypp, yppp)
    d_phi = yp
  end function d_phi_spline

  function dd_phi_spline(m, x) result(dd_phi)
    integer       :: m
    real(kind=dp) :: x, dd_phi
    integer       :: swd
    real(kind=dp) :: m0
    real(kind=dp) :: y, yp, ypp, yppp
    
    swd = 1
    m0 = 0d0

    call splint_horner3(x_dat, a_spl(m,:), b_spl(m,:), c_spl(m,:), d_spl(m,:), swd, m0, &
         x, tfone, tfzero, tfzero, tfzero, y, yp, ypp, yppp)
    dd_phi = ypp
  end function dd_phi_spline
  
end module collop_spline
