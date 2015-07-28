module collop_definitions
  use nrtype, only : dp, pi
  use gsl_integration_routines_mod
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
  character(len=128) :: phi_m_desc = 'Associated Laguerre polynomials (Default)'
  integer            :: phi_m_tag  = 0
  real(kind=dp), dimension(:),     allocatable :: phi_hm
  real(kind=dp), dimension(:,:),   allocatable :: coefleg
  real(kind=dp), dimension(:,:,:), allocatable :: gencoeflag, phitrans
  
  !**********************************************************
  ! Matrix size
  !**********************************************************
  integer, parameter :: lagmax = 5
  integer, parameter :: legmax = 5

  !**********************************************************
  ! Integration settings
  !**********************************************************
  real(kind=dp) :: epsabs = 1d-10
  real(kind=dp) :: epsrel = 1d-10
  integer       :: sw_qag_rule = 2
  !real(kind=dp) :: x_max    = 20
contains
  
  subroutine init_legendre(n)
    !
    ! Computes coefficients of Legendre polynomials of orders from 0 to n

    !
    ! Input parameters:
    !           Formal: n            - maximum order of Legendre polynomials
    ! Output parameters:
    !           Formal: coefleg(i,j) - j-th coefficient of Legendre polynomial
    !                                  of the order i
    !
    integer :: n,i,j
    !
    double precision :: frontfac,rearfac
    !
    if(allocated(coefleg)) return
    write (*,*) "Initializing Legendre coefficients..."
    allocate(coefleg(0:n,0:n))
    !
    coefleg=0.d0
    coefleg(0,0)=1.d0
    coefleg(1,1)=1.d0
    frontfac=1.d0
    !
    do i=2,n
       frontfac=frontfac*(2.d0-1.d0/dble(i))
       rearfac=frontfac
       coefleg(i,i)=rearfac
       do j=i-2,0,-2
          rearfac=-rearfac*dble(j+1)*dble(j+2)/dble(i-j)/dble(i+j+1)
          coefleg(i,j)=rearfac
       enddo
    enddo
    write (*,*) "Done."
  end subroutine init_legendre
  
  subroutine init_laguerre(lagmax, legmax)
    integer :: lagmax, legmax
    integer :: i, j, k, l
    real(kind=dp) :: bincoef

    if (allocated(gencoeflag)) return
    write (*,*) "Initializing Laguerre coefficients..."
    allocate(gencoeflag(0:lagmax, 0:lagmax, 0:legmax))
    gencoeflag=0.d0
    gencoeflag(0,0,:)=1.d0

    do l=1, legmax
       do i=1, lagmax
          do k=0, i
             bincoef=1.d0
             do j=1,i-k
                bincoef=bincoef*(1.d0+(k+l+0.5d0)/j)
             enddo
             do j=1,k
                bincoef=-bincoef/j
             enddo
             gencoeflag(i,k,l)=bincoef
          enddo
       enddo
    enddo
    write (*,*) "Done."
  end subroutine init_laguerre

  subroutine init_phi(n)
    integer :: n
    integer :: m
    if (allocated(phi_hm)) deallocate(phi_hm)
    allocate(phi_hm(0:n))

    do m = 0, n
       phi_hm(m) = 1d0/sqrt(gamma(m + 2.5d0) / (2d0 * gamma(m+1d0)))
    end do

  end subroutine init_phi

  function phi(m, x)
    integer       :: m, k
    real(kind=dp) :: phi
    real(kind=dp) :: x, plag, xpow, add
    real(kind=dp) :: x2

    x2 = x**2
    plag=gencoeflag(m,0,1)
    xpow=1.d0

    do k=1,m
       add=gencoeflag(m,k,1)*xpow
       plag=plag+add*x2
       xpow=xpow*x2
    enddo

    phi = plag * pi**(3d0/4d0) * phi_hm(m)
  end function phi

  function d_phi(m, x)
    integer       :: m, k
    real(kind=dp) :: d_phi
    real(kind=dp) :: x, dplag, xpow, add
    real(kind=dp) :: x2

    x2 = x**2
    dplag=0.d0
    xpow=1.d0

    do k=1,m
       add=gencoeflag(m,k,1)*xpow
       dplag=dplag+k*add
       xpow=xpow*x2
    enddo
    
    d_phi = 2d0 * x * dplag * pi**(3d0/4d0) * phi_hm(m)
  end function d_phi
  
  function dd_phi(m, x)
    integer :: m,k
    real(kind=dp) :: dd_phi
    real(kind=dp) :: x, dplag, ddplag, xpow,add
    real(kind=dp) :: x2

    x2 = x**2
    dplag=0.d0
    ddplag=0d0
    xpow=1.d0

    do k=1,m
       add=gencoeflag(m,k,1)*xpow
       dplag=dplag+k*add
       ddplag = ddplag + k*(k-1)*add/x2
       xpow=xpow*x2
    enddo

    dd_phi = (2d0 * dplag + 4d0*x**2*ddplag) * pi**(3d0/4d0) * phi_hm(m)
  end function dd_phi

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
