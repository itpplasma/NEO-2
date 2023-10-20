module collop_definitions
  use nrtype, only : dp, pi
  use gsl_integration_routines_mod
  use collop_laguerre
  use collop_polynomial
  implicit none

  !**********************************************************
  ! Weighting
  !**********************************************************
  real(kind=dp)    :: alpha  = 0
  real(kind=dp)    :: beta   = 0

  !**********************************************************
  ! Species
  !**********************************************************
  real(kind=dp)    :: m_a
  real(kind=dp)    :: m_b
  character(len=3) :: tag_a
  character(len=3) :: tag_b

  !**********************************************************
  ! Test function
  !**********************************************************
  real(kind=dp), dimension(:,:),   allocatable :: coefleg
   
  !**********************************************************
  ! Matrix size
  !**********************************************************
  integer :: lagmax
  integer :: legmax

  !**********************************************************
  ! Integration settings
  !**********************************************************
  real(kind=dp) :: epsabs = 1d-10
  real(kind=dp) :: epsrel = 1d-10
  integer       :: sw_qag_rule = 2

  abstract interface
     subroutine init_phi_interface(lagmax, legmax)
       integer :: lagmax, legmax
     end subroutine init_phi_interface

     function phi_interface(m,x)
       use nrtype, only : dp
       integer :: m
       real(kind=dp) :: x, phi_interface
     end function phi_interface

     function d_phi_interface(m,x)
       use nrtype, only : dp
       integer :: m
       real(kind=dp) :: x, d_phi_interface
     end function d_phi_interface

     function dd_phi_interface(m,x)
       use nrtype, only : dp
       integer :: m
       real(kind=dp) :: x, dd_phi_interface
     end function dd_phi_interface
  end interface

  !**********************************************************
  ! Function pointers for different base functions
  !**********************************************************
  procedure(init_phi_interface), pointer :: init_phi => null()
  procedure(phi_interface),      pointer :: phi      => null()
  procedure(d_phi_interface),    pointer :: d_phi    => null()
  procedure(dd_phi_interface),   pointer :: dd_phi   => null()
  
contains
  
  subroutine init_legendre(n)
    ! Computes coefficients of Legendre polynomials of orders from 0 to n
    !
    ! Input parameters:
    !           Formal: n            - maximum order of Legendre polynomials
    ! Output parameters:
    !           Formal: coefleg(i,j) - j-th coefficient of Legendre polynomial
    !                                  of the order i

    integer :: n,i,j

    double precision :: frontfac,rearfac

    if(allocated(coefleg)) return
    write (*,*) "Initializing Legendre coefficients..."
    allocate(coefleg(0:n,0:n))

    coefleg=0.d0
    coefleg(0,0)=1.d0
    coefleg(1,1)=1.d0
    frontfac=1.d0

    do i=2,n
       frontfac=frontfac*(2.d0-1.d0/dble(i))
       rearfac=frontfac
       coefleg(i,i)=rearfac
       do j=i-2,0,-2
          rearfac=-rearfac*dble(j+1)*dble(j+2)/dble(i-j)/dble(i+j+1)
          coefleg(i,j)=rearfac
       enddo
    enddo
  end subroutine init_legendre

  function G(x) result(y)
    real(kind=dp) :: x
    real(kind=dp) :: y

    y = (erf(x) - x*d_erf(x)) / (2*x**2)

  end function G

  function d_erf(x) result(y)
    real(kind=dp) :: x
    real(kind=dp) :: y

    y = 2d0/sqrt(pi) * exp(-x**2)
  end function d_erf

  function dd_erf(x) result(y)
    real(kind=dp) :: x
    real(kind=dp) :: y

    y = -4d0*x/sqrt(pi) * exp(-x**2)
  end function dd_erf

  function d_G(x) result(y)
    real(kind=dp) :: x
    real(kind=dp) :: y

    y = (4*x**2*d_erf(x) - 4*x*erf(x) - 2*x**3*dd_erf(x))
    y = y / (2*x**2)**2
  end function d_G

end module collop_definitions
