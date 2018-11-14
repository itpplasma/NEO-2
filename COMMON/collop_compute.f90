module collop_compute

  use hdf5_tools
  !use nrtype, only : dp, pi
  use gsl_integration_routines_mod
  use gsl_specialfunctions_mod
  use collop_laguerre
  use collop_polynomial
  use collop_spline
  use collop_bspline
  use collop_nbi
  !use collop_bspline2
  !use collop_spline4
  use ieee_arithmetic
  
  implicit none

  !**********************************************************
  ! Species definitions (cgs)
  !**********************************************************
  real(kind=dp) :: m_ele = 9.109382150d-28
  real(kind=dp) :: m_pro = 1.672621637d-24
  real(kind=dp) :: m_alp = 6.644656200d-24
  real(kind=dp) :: m_d   = 3.343583719d-24
  real(kind=dp) :: m_C   = 19.94406876d-24
  real(kind=dp) :: m_W   = 305.2734971d-24
  real(kind=dp) :: c     = 2.9979d10
  real(kind=dp) :: ev    = 1.6022d-12

  !**********************************************************
  ! Thermal velocity ratio
  !**********************************************************
  real(kind=dp) :: gamma_ab

  !**********************************************************
  ! Collision operator matrices
  !**********************************************************
  integer :: ispec
  real(kind=dp), dimension(:,:),   allocatable :: asource_s, anumm_s, denmm_s
  real(kind=dp), dimension(:,:,:), allocatable :: I1_mmp_s, I2_mmp_s, I3_mmp_s, I4_mmp_s, ailmm_s
  real(kind=dp), dimension(:,:),   allocatable :: M_transform, M_transform_inv
  real(kind=dp), dimension(:),     allocatable :: C_m
  real(kind=dp), dimension(:),     allocatable :: weightenerg_offset

  !**********************************************************
  ! Profile
  !**********************************************************
  real(kind=dp) :: T_a
  real(kind=dp) :: T_b
  real(kind=dp) :: v_ta
  real(kind=dp) :: v_tb

  !**********************************************************
  ! Weighting
  !**********************************************************
  real(kind=dp)    :: alpha  = 0d0
  real(kind=dp)    :: beta   = 0d0

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
  ! Relativistic parameters
  !**********************************************************
  real(kind=dp) :: rmu, norm_maxwell
  integer       :: isw_relativistic
  real(kind=dp) :: a_00_offset = 1.00d0
  real(kind=dp) :: a_02_offset = 1.50d0
  real(kind=dp) :: a_22_offset = 3.75d0
  
  !**********************************************************
  ! Matrix size
  !**********************************************************
  integer :: lagmax
  integer :: legmax

  !**********************************************************
  ! Integration settings
  !**********************************************************
  real(kind=dp) :: epsabs = 1d-13
  real(kind=dp) :: epsrel = 1d-13
  integer       :: sw_qag_rule = 5
  logical       :: integral_cutoff = .true.
  real(kind=dp) :: x_cutoff = 3.0d3
  logical       :: lsw_interval_sep = .true.
  logical       :: lsw_split_interval = .true. ! split intervals manually into further sub-intervals
  integer       :: num_sub_intervals = 5 ! number of sub-intervals between nodes (active if lsw_split_interval)
  integer       :: num_sub_intervals_cutoff = 40 ! number of sub-intervals between last node and cutoff (active if lsw_split_interval)
  ! Taylor expansion of (non-relativistic) matrix elements around x=0 (numerical stability):
  logical       :: lsw_expand_kernel =.false. ! Taylor expansion of integral-part subintegrands
  logical       :: lsw_expand_G =.true. ! Taylor expansion of Chandrasekhar function
  ! corresponding domain of Taylor expansion (0 <= x <= xmax_expand_kernel)
  real(kind=dp) :: xmax_expand_kernel = 1.0d-5
  ! numerical stabilization of integral-part computation (I2 + I3 + I4)
  logical       :: lsw_stabilize_Immp=.true.
  
  !**********************************************************
  ! Pre-computed matrix elements
  !**********************************************************
  !character(len=100) :: matelem_name='MatElem_aa_hatfun_lag10_xmax6.h5'
  !character(len=100) :: matelem_name='MatElem_aa_hatfun_lag5_xmax6.h5'
  !character(len=100) :: matelem_name='MatElem_aa_hatfun_lag14_xmax5.h5'
  !character(len=100) :: matelem_name='MatElem_aa_hatfun_lag30_xmax4.h5'
  character(len=100) :: matelem_name='MatElem_aa_hatfun_lag28_xmax4.h5'
  !character(len=100) :: matelem_name='MatElem_aa_hatfun_lag24_xmax4.h5'
  !character(len=100) :: matelem_name='MatElem_aa_hatfun_lag22_xmax4.h5'
  !character(len=100) :: matelem_name='MatElem_aa_hatfun_lag20_xmax4.h5'
  !character(len=100) :: matelem_name='MatElem_aa_hatfun_lag19_xmax4.h5'
  !character(len=100) :: matelem_name='MatElem_aa_hatfun_lag17_xmax4.h5'
  !character(len=100) :: matelem_name='MatElem_aa_hatfun_lag15_xmax4.h5'
  !character(len=100) :: matelem_name='MatElem_aa_hatfun_lag14_xmax4.h5'
  !character(len=100) :: matelem_name='MatElem_aa_hatfun_lag12_xmax4.h5'
  !character(len=100) :: matelem_name='MatElem_aa_hatfun_lag10_xmax4.h5'
  !character(len=100) :: matelem_name='MatElem_aa_hatfun_lag10_xmax4.h5'
  !character(len=100) :: matelem_name='MatElem_aa_hatfun_lag7_xmax4.h5'
  !character(len=100) :: matelem_name='MatElem_aa_hatfun_lag5_xmax4.h5'
  !character(len=100) :: matelem_name='MatElem_aa_hatfun_lag10_xmax3.h5'
  !character(len=100) :: matelem_name='MatElem_aa_hatfun_lag7_xmax3.h5'
  !character(len=100) :: matelem_name='MatElem_aa_hatfun_lag5_xmax3.h5'
  !character(len=100) :: matelem_name='MatElem_aa_hatfun_lag10_xmax2.h5'
  !character(len=100) :: matelem_name='MatElem_aa_hatfun_lag7_xmax2.h5'
  !character(len=100) :: matelem_name='MatElem_aa_hatfun_lag5_xmax2.h5'
  integer(HID_T) :: h5id_matelem
  logical :: precomp=.false.
  logical :: make_ortho=.true.

  interface chop
     module procedure chop_0
     module procedure chop_1
     module procedure chop_2
     module procedure chop_3
     module procedure chop_4
     module procedure chop_5
  end interface chop

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
  procedure(init_phi_interface), pointer :: init_phi_prj => null()
  procedure(phi_interface),      pointer :: phi_prj      => null()
  procedure(d_phi_interface),    pointer :: d_phi_prj    => null()
  procedure(dd_phi_interface),   pointer :: dd_phi_prj   => null()

  procedure(init_phi_interface), pointer :: init_phi_exp => null()
  procedure(phi_interface),      pointer :: phi_exp      => null()
  procedure(d_phi_interface),    pointer :: d_phi_exp    => null()
  procedure(dd_phi_interface),   pointer :: dd_phi_exp   => null()

  !**********************************************************
  ! jump discontinuities
  !**********************************************************
  logical       :: d_phi_discont = .false.
  real(kind=dp), dimension(:), allocatable :: t_vec

  interface integrate
     module procedure integrate_a_to_infinity, integrate_a_to_b
  end interface integrate

  interface integrate_param
     module procedure integrate_a_to_infinity_param, integrate_a_to_b_param
  end interface integrate_param

contains

!     ..................................................................
!
  subroutine DBESK(X,N,BK,IER)
!
      double precision X, BK
      double precision T, B
      double precision G0, G1, GJ
      integer N,IER,L,J
!
      dimension T(12)
      BK=.0
      if(N)10,11,11
   10 IER=1
      return
   11 if(X)12,12,20
   12 IER=2
      return
   20 continue
   22 IER=0
      if(X-1.)36,36,25
   25 B=1./X
      T(1)=B
      do 26 L=2,12
   26 T(L)=T(L-1)*B
      if(N-1)27,29,27
!
!     COMPUTE KO USING POLYNOMIAL APPROXIMATION
!
   27 G0=(1.2533141-.1566642*T(1)+.08811128*T(2)-.09139095*T(3)   &
      +.1344596*T(4)-.2299850*T(5)+.3792410*T(6)-.5247277*T(7)    &
      +.5575368*T(8)-.4262633*T(9)+.2184518*T(10)-.06680977*T(11) &
      +.009189383*T(12))
      if(N)20,28,29
   28 BK=G0
      return
!
!     COMPUTE K1 USING POLYNOMIAL APPROXIMATION
!
   29 G1=(1.2533141+.4699927*T(1)-.1468583*T(2)+.1280427*T(3)     &
      -.1736432*T(4)+.2847618*T(5)-.4594342*T(6)+.6283381*T(7)    &
      -.6632295*T(8)+.5050239*T(9)-.2581304*T(10)+.07880001*T(11) &
      -.01082418*T(12))
      if(N-1)20,30,31
   30 BK=G1
      return
!
!     FROM KO,K1 COMPUTE KN USING RECURRENCE RELATION
!
   31 do 35 J=2,N
      GJ=2.*(FLOAT(J)-1.)*G1/X+G0
      if(GJ-1.0D300)33,33,32
   32 IER=4
      GO TO 34
   33 G0=G1
   35 G1=GJ
   34 BK=GJ
      return
!   
   36 IER=3
      BK=0.
      return
  end
!
  
  subroutine init_collop(collop_base_prj, collop_base_exp, scalprod_alpha, scalprod_beta)
    use rkstep_mod, only : lag, leg
    use collop_bspline, only : xknots, nder
    integer :: collop_base_prj, collop_base_exp
    real(kind=dp) :: scalprod_alpha
    real(kind=dp) :: scalprod_beta

    alpha = scalprod_alpha
    beta  = scalprod_beta
    lagmax = lag
    legmax = leg

    select case (collop_base_prj)
    case (0)
       write (*,*) "Using Laguerre polynomials as collision operator projection base."
       init_phi_prj => init_phi_laguerre
       phi_prj      => phi_laguerre
       d_phi_prj    => d_phi_laguerre
       dd_phi_prj   => dd_phi_laguerre
    case (1)
       write (*,*) "Using standard polynomials as collision operator projection base."
       init_phi_prj => init_phi_polynomial
       phi_prj      => phi_polynomial
       d_phi_prj    => d_phi_polynomial
       dd_phi_prj   => dd_phi_polynomial     
    case (2)
       write (*,*) "Using squared polynomials as collision operator projection base."
       init_phi_prj => init_phi_polynomial_2
       phi_prj      => phi_polynomial_2
       d_phi_prj    => d_phi_polynomial_2
       dd_phi_prj   => dd_phi_polynomial_2
    case (3)
       write (*,*) "Using squared polynomials without zeroth order as collision operator projection base."
       init_phi_prj => init_phi_polynomial_3
       phi_prj      => phi_polynomial_3
       d_phi_prj    => d_phi_polynomial_3
       dd_phi_prj   => dd_phi_polynomial_3
    case (10)
       write (*,*) "Using splines as collision operator projection base."
       init_phi_prj => init_phi_spline
       phi_prj      => phi_spline
       d_phi_prj    => d_phi_spline
       dd_phi_prj   => dd_phi_spline
    case (11)
       write (*,*) "Using B-Splines as collision operator projection base."
       init_phi_prj => init_phi_bspline
       phi_prj      => phi_bspline
       d_phi_prj    => d_phi_bspline
       dd_phi_prj   => dd_phi_bspline
    !case (12)
    !   write (*,*) "Using 4th-order Splines as collision operator projection base."
    !   init_phi_prj => init_phi_spline4
    !   phi_prj      => phi_spline4
    !   d_phi_prj    => d_phi_spline4
    !   dd_phi_prj   => dd_phi_spline4
    !case (13)
    !   write (*,*) "Using B-Splines with x**2 as collision operator projection base."
    !   init_phi_prj => init_phi_bspline2
    !   phi_prj      => phi_bspline2
    !   d_phi_prj    => d_phi_bspline2
    !   dd_phi_prj   => dd_phi_bspline2
    case (100)
       write (*,*) "Using hat functions as collision operator projection base."
       precomp=.true.
    case default
       write (*,*) "Undefined collision operator projection base ", collop_base_prj
       stop
    end select

    select case (collop_base_exp)
    case (0)
       write (*,*) "Using Laguerre polynomials as collision operator expansion base."
       init_phi_exp => init_phi_laguerre
       phi_exp      => phi_laguerre
       d_phi_exp    => d_phi_laguerre
       dd_phi_exp   => dd_phi_laguerre
    case (1)
       write (*,*) "Using standard polynomials as collision operator expansion base."
       init_phi_exp => init_phi_polynomial
       phi_exp      => phi_polynomial
       d_phi_exp    => d_phi_polynomial
       dd_phi_exp   => dd_phi_polynomial     
    case (2)
       write (*,*) "Using squared polynomials as collision operator expansion base."
       init_phi_exp => init_phi_polynomial_2
       phi_exp      => phi_polynomial_2
       d_phi_exp    => d_phi_polynomial_2
       dd_phi_exp   => dd_phi_polynomial_2
    case (3)
       write (*,*) "Using squared polynomials without zeroth order as collision operator expansion base."
       init_phi_exp => init_phi_polynomial_3
       phi_exp      => phi_polynomial_3
       d_phi_exp    => d_phi_polynomial_3
       dd_phi_exp   => dd_phi_polynomial_3
    case (10)
       write (*,*) "Using splines as collision operator expansion base."
       init_phi_exp => init_phi_spline
       phi_exp      => phi_spline
       d_phi_exp    => d_phi_spline
       dd_phi_exp   => dd_phi_spline
    case (11)
       write (*,*) "Using B-Splines as collision operator expansion base."
       init_phi_exp => init_phi_bspline
       phi_exp      => phi_bspline
       d_phi_exp    => d_phi_bspline
       dd_phi_exp   => dd_phi_bspline
    !case (12)
    !   write (*,*) "Using 4th-order Splines as collision operator expansion base."
    !   init_phi_exp => init_phi_spline4
    !   phi_exp      => phi_spline4
    !   d_phi_exp    => d_phi_spline4
    !   dd_phi_exp   => dd_phi_spline4
    !case (13)
    !   write (*,*) "Using B-Splines with x**2 as collision operator expansion base."
    !   init_phi_exp => init_phi_bspline2
    !   phi_exp      => phi_bspline2
    !   d_phi_exp    => d_phi_bspline2
    !   dd_phi_exp   => dd_phi_bspline2
    case (100)
       write (*,*) "Using hat functions as collision operator expansion base."
    case default
       write (*,*) "Undefined collision operator expansion base ", collop_base_exp
       stop
    end select
  
    call init_legendre(legmax)
    call init_phi_prj(lagmax, legmax)
    call init_phi_exp(lagmax, legmax)

    !**********************************************************
    ! Special for B-Splines
    ! Inform integration routines about knots and continuity of
    ! bases functions
    !**********************************************************
    if(lsw_interval_sep) then
       if ((collop_base_prj .eq. 11) .or. (collop_base_exp .eq. 11)) then
          if (allocated(t_vec)) deallocate(t_vec)
          allocate(t_vec(lbound(xknots,1):ubound(xknots,1)))
          t_vec = xknots
          !write (*,*) t_vec

          !**********************************************************
          ! Detect if hat functions
          !**********************************************************
          if (nder .eq. 2) then  

             !**********************************************************
             ! Cutoff integration after last hat
             !**********************************************************
             integral_cutoff = .true.
             x_cutoff = phi_x_max

             !**********************************************************
             ! Avoid second derivative in integral part (I4)
             !**********************************************************
             d_phi_discont = .true.
          end if

       end if
    end if
 
  end subroutine init_collop

  subroutine chop_0(x)
    real(kind=dp) :: x, chop

    !if (abs(x) .lt. 1d-10) then
    !   x = 0d0
       !elseif (abs(x - 1d0) .lt. 1d-10) then
       !   x = 1d0
    !end if
    x = x
    
  end subroutine chop_0

  subroutine chop_1(x)
    real(kind=dp), dimension(:) :: x
    integer :: k

    do k = lbound(x,1), ubound(x,1)
       call chop(x(k))
    end do

  end subroutine chop_1

  subroutine chop_2(x)
    real(kind=dp), dimension(:,:) :: x
    integer :: k

    do k = lbound(x,1), ubound(x,1)
       call chop(x(k,:))
    end do

  end subroutine chop_2

  subroutine chop_3(x)
    real(kind=dp), dimension(:,:,:) :: x
    integer :: k

    do k = lbound(x,1), ubound(x,1)
       call chop(x(k,:,:))
    end do

  end subroutine chop_3

  subroutine chop_4(x)
    real(kind=dp), dimension(:,:,:,:) :: x
    integer :: k

    do k = lbound(x,1), ubound(x,1)
       call chop(x(k,:,:,:))
    end do

  end subroutine chop_4

  subroutine chop_5(x)
    real(kind=dp), dimension(:,:,:,:,:) :: x
    integer :: k

    do k = lbound(x,1), ubound(x,1)
       call chop(x(k,:,:,:,:))
    end do

  end subroutine chop_5

  function dmuk2ovk2(mu)
    real(kind=dp) :: mu, dmuk2ovk2
    
    if (mu .ge. 5d2) then
       dmuk2ovk2 = -1.5d0 - 15d0/(8d0*mu) + 15d0/(8d0*mu**2)                    &
            - 135d0/(128d0 * mu**3) - 45d0/(32d0*mu**4) + 7425d0/(1024d0*mu**5) &
            - 675d0/(32d0*mu**6) + 1905525d0/(32768*mu**7)
    else
       dmuk2ovk2 = -3d0 + mu - (mu*besselk(1,mu))/besselk(2,mu)
    end if
  end function dmuk2ovk2
  
  function ddmuk2ovk2(mu)
    real(kind=dp) :: mu, ddmuk2ovk2

    if (mu .ge. 5d2) then
       ddmuk2ovk2 = 15d0/4d0 + 75d0/(8d0 * mu) - 495d0/(64d0 * mu**2) &
            + 45d0/(128d0 * mu**3) + 9585d0 / (512d0 * mu**4)         &
            - 65475d0/(1024d0 * mu**5) + 2942325d0/(16384d0 * mu**6)  &
            - 17380575d0/(32768d0*mu**7)
    else
       ddmuk2ovk2 = (2*mu*(6 + (-3 + mu)*mu)*besselk(0,mu) &
            + (24 + mu*(-12 + (7 - 2*mu)*mu))*besselk(1,mu))/(mu*besselk(2,mu))
    end if
  end function ddmuk2ovk2
 
  subroutine karney(z,gamma,j0_1,j0_11,j0_2,j0_02,j0_022,               &
       j1_0,j1_1,j1_2,j1_02,j1_11,j1_022)
    !
    implicit none
    !
    double precision :: gamma,j0_1, j0_11, j0_2,j0_02,j0_022,           &
         j1_0,j1_1,j1_2,j1_02,j1_11,j1_022
    double precision :: z,sigma,z2
    !
    z2=z**2
    !
    if(z.lt.1d-2) then
       gamma=1.d0+z2*(0.5d0+z2*(-0.125d0+z2/16.d0))
       j0_1=1d0
       j0_11=z2/6d0 - z2**2/15d0 + 4*z**6/105d0 + 8*z**8/315d0 + 64*z**10/3465d0
       j0_2=gamma
       j0_02=z2*(1.d0/6.d0+z2*(-0.05d0+z2*(3.d0/112.d0-5.d0*z2/288.d0)))
       j0_022=z2**2*(1.d0/120.d0+z2*(-3.d0/560.d0+z2*(5.d0/1344.d0 &
            -35.d0*z2/12672.d0)))
       j1_0=z*(1.d0/3.d0+z2*(-2.d0/15.d0+z2*(8.d0/105.d0-16.d0*z2/315.d0)))
       j1_1=z*(1.d0/3.d0+z2*(-0.1d0+z2*(3.d0/56.d0-5.d0*z2/144.d0)))
       j1_2=z/3.d0
       j1_02=z*z2*(1.d0/30.d0+z2*(-2.d0/105.d0+z2*(4.d0/315.d0 &
            -32.d0*z2/3465.d0)))
       j1_11=z*z2*(1.d0/30.d0+z2*(-3.d0/140.d0+z2*(5.d0/336.d0 &
            -35.d0*z2/3168.d0)))
       j1_022=z*z2**2*(1.d0/840.d0+z2*(-1.d0/945.d0+z2*(1.d0/1155.d0 &
            -32.d0*z2/45045.d0)))
       return
    endif
    !
    gamma=dsqrt(z2+1.d0)
    sigma=dlog(z+gamma)
    !
    j0_1=1d0
    j0_11=(sigma*gamma - z)/(2d0*z)
    !
    j0_2=gamma
    j0_02=(z*gamma-sigma)/(4.d0*z)
    j0_022=((3.d0+2.d0*z2)*sigma-3.d0*z*gamma)/(32.d0*z)
    !
    j1_0=(gamma*sigma-z)/z2
    j1_1=(z*gamma-sigma)/(2.d0*z2)
    j1_2=z/3.d0
    j1_02=(3.d0*z+z**3-3.d0*gamma*sigma)/(12.d0*z2)
    j1_11=((3.d0+2.d0*z2)*sigma-3.d0*z*gamma)/(8.d0*z2)
    j1_022=((15.d0+6.d0*z2)*gamma*sigma-15.d0*z-11.d0*z**3)/288.d0/z2
    !
    !
    return
  end subroutine karney

  function D_thth(rmu, x)
    real(kind=dp) :: rmu, x, D_thth
    real(kind=dp) :: z, gam
    real(kind=dp) :: j0_1,j0_11,j0_2,j0_02,j0_022,j1_0,j1_1,j1_2,j1_02,j1_11,j1_022
    real(kind=dp) :: I1, I2
    
    z   = sqrt(2/rmu) * x
    gam = sqrt(1+z**2)
    call karney(z, gam,j0_1,j0_11,j0_2,j0_02,j0_022,j1_0,j1_1,j1_2,j1_02,j1_11,j1_022)

    I1 = integrate(D_thth1, 0d0, x)
    I2 = integrate(D_thth2, x)

    D_thth = 4/sqrt(pi) * norm_maxwell * (I1 + I2)
  contains

    function D_thth1(xp)
      real(kind=dp) :: xp, D_thth1
      real(kind=dp) :: zp, gamp, exp_maxwellp
      real(kind=dp) :: j0_1p,j0_11p,j0_2p,j0_02p,j0_022p,j1_0p,j1_1p,j1_2p,j1_02p,j1_11p,j1_022p
      
      zp  = sqrt(2/rmu) * xp
      gamp = sqrt(1+zp**2)
      exp_maxwellp = 2d0/(sqrt(1+2*xp**2/rmu) + 1)
      
      call karney(zp, gamp, j0_1p,j0_11p,j0_2p,j0_02p,j0_022p,j1_0p,j1_1p,j1_2p,j1_02p,j1_11p,j1_022p)
      D_thth1 = ((0.5*j0_2p - (z**(-2) + 1d0/gam**2d0)*j0_02p)  &
           + 4d0/(gam**2)*z**(-2)*j0_022p)*gam/gamp*xp**2/x &
           * exp(-exp_maxwellp * xp**2)

    end function D_thth1

    function D_thth2(xp)
      real(kind=dp) :: xp, D_thth2
      real(kind=dp) :: zp, gamp, exp_maxwellp
      
      zp  = sqrt(2/rmu) * xp
      gamp = sqrt(1+zp**2)
      exp_maxwellp = 2d0/(sqrt(1+2*xp**2/rmu) + 1)
      
      D_thth2 = (0.5*gamp**2/gam**2*j0_2 - xp**2/x**2*(zp**(-2) + 1d0/gam**2d0)*j0_02 &
           + 4d0/(gam**2)*z**(-2)*j0_022)*gam/gamp*xp &
           * exp(-exp_maxwellp * xp**2)

    end function D_thth2
    
  end function D_thth

  function D_uu(rmu, x)
    real(kind=dp) :: rmu, x, D_uu
    real(kind=dp) :: z, gam
    real(kind=dp) :: j0_1,j0_11,j0_2,j0_02,j0_022,j1_0,j1_1,j1_2,j1_02,j1_11,j1_022
    real(kind=dp) :: I1, I2
    
    z   = sqrt(2/rmu) * x
    gam = sqrt(1+z**2)
    call karney(z, gam, j0_1,j0_11,j0_2,j0_02,j0_022,j1_0,j1_1,j1_2,j1_02,j1_11,j1_022)

    I1 = integrate(D_uu1, 0d0, x)
    I2 = integrate(D_uu2, x)
    
    D_uu = 4d0/sqrt(pi) * norm_maxwell * (I1 + I2)
    !write (400,*) I, I1, I2, norm_maxwell, 4d0/sqrt(pi)* norm_maxwell * (I1 + I2), D_uu
  contains

    function D_uu1(xp)
      real(kind=dp) :: xp, D_uu1
      real(kind=dp) :: zp, gamp, exp_maxwellp
      real(kind=dp) :: j0_1p,j0_11p,j0_2p,j0_02p,j0_022p,j1_0p,j1_1p,j1_2p,j1_02p,j1_11p,j1_022p
      
      zp  = sqrt(2/rmu) * xp
      gamp = sqrt(1+zp**2)
      exp_maxwellp = 2d0/(sqrt(1+2*xp**2/rmu) + 1)
      
      call karney(zp, gamp,j0_1p,j0_11p,j0_2p,j0_02p,j0_022p,j1_0p,j1_1p,j1_2p,j1_02p,j1_11p,j1_022p)
      D_uu1 = ((2*gam**2*j0_02p - 8*j0_022p)*gam*z**(-2) * (xp**2/(x*gamp))) * exp(-exp_maxwellp * xp**2)
    end function D_uu1

    function D_uu2(xp)
      real(kind=dp) :: xp, D_uu2
      real(kind=dp) :: zp, gamp, exp_maxwellp
      
      zp  = sqrt(2/rmu) * xp
      gamp = sqrt(1+zp**2)
      exp_maxwellp = 2d0/(sqrt(1+2*xp**2/rmu) + 1)
      
      D_uu2 = (2*gamp**2*j0_02 - 8*j0_022) * z**(-2) * gam/gamp * xp * exp(-exp_maxwellp * xp**2)
    end function D_uu2
    
  end function D_uu

  function D_uu_b(rmu, x, mp)
    real(kind=dp) :: rmu, x, D_uu_b
    integer       :: mp
    real(kind=dp) :: z, gam
    real(kind=dp) :: j0_1,j0_11,j0_2,j0_02,j0_022,j1_0,j1_1,j1_2,j1_02,j1_11,j1_022
    real(kind=dp) :: I1, I2
    
    z   = sqrt(2/rmu) * x
    gam = sqrt(1+z**2)
    call karney(z, gam,j0_1,j0_11,j0_2,j0_02,j0_022,j1_0,j1_1,j1_2,j1_02,j1_11,j1_022)

    I1 = integrate(D_uu_b1, 0d0, x)
    I2 = integrate(D_uu_b2, x)
    
    D_uu_b = 4d0/sqrt(pi) * norm_maxwell * (I1 + I2)
    !write (400,*) I, I1, I2, norm_maxwell, 4d0/sqrt(pi)* norm_maxwell * (I1 + I2), D_uu
  contains

    function D_uu_b1(xp)
      real(kind=dp) :: xp, D_uu_b1
      real(kind=dp) :: zp, gamp, exp_maxwellp
      real(kind=dp) :: j0_1p,j0_11p,j0_2p,j0_02p,j0_022p,j1_0p,j1_1p,j1_2p,j1_02p,j1_11p,j1_022p
      
      zp  = sqrt(2/rmu) * xp
      gamp = sqrt(1+zp**2)
      exp_maxwellp = 2d0/(sqrt(1+2*xp**2/rmu) + 1)
      
      call karney(zp, gamp,j0_1p,j0_11p,j0_2p,j0_02p,j0_022p,j1_0p,j1_1p,j1_2p,j1_02p,j1_11p,j1_022p)
      D_uu_b1 = ((2*gam**2*j0_02p - 8*j0_022p)*gam*z**(-2) * (xp**2/(x*gamp))) * exp(-exp_maxwellp * xp**2) * phi_exp(mp,xp)
    end function D_uu_b1

    function D_uu_b2(xp)
      real(kind=dp) :: xp, D_uu_b2
      real(kind=dp) :: zp, gamp, exp_maxwellp
      
      zp  = sqrt(2/rmu) * xp
      gamp = sqrt(1+zp**2)
      exp_maxwellp = 2d0/(sqrt(1+2*xp**2/rmu) + 1)
      
      D_uu_b2 = (2*gamp**2*j0_02 - 8*j0_022) * z**(-2) * gam/gamp * xp * exp(-exp_maxwellp * xp**2) * phi_exp(mp,xp)
    end function D_uu_b2
    
  end function D_uu_b  

  function F_u_b(rmu, x, mp)
    real(kind=dp) :: rmu, x, F_u_b
    integer       :: mp
    real(kind=dp) :: z, gam
    real(kind=dp) :: j0_1,j0_11,j0_2,j0_02,j0_022,j1_0,j1_1,j1_2,j1_02,j1_11,j1_022
    real(kind=dp) :: I1, I2
    z   = sqrt(2/rmu) * x
    gam = sqrt(1+z**2)
    call karney(z, gam,j0_1,j0_11,j0_2,j0_02,j0_022,j1_0,j1_1,j1_2,j1_02,j1_11,j1_022)
    
    I1 = integrate(F_u_b1, 0d0, x)
    I2 = integrate(F_u_b2, x)
    
    F_u_b = -4d0/sqrt(pi) * norm_maxwell * (I1 + I2)
    !write (400,*) I, I1, I2, norm_maxwell, 4d0/sqrt(pi)* norm_maxwell * (I1 + I2), D_uu
  contains

    function F_u_b1(xp)
      real(kind=dp) :: xp, F_u_b1
      real(kind=dp) :: zp, gamp, exp_maxwellp
      real(kind=dp) :: j0_1p,j0_11p,j0_2p,j0_02p,j0_022p,j1_0p,j1_1p,j1_2p,j1_02p,j1_11p,j1_022p
      
      zp  = sqrt(2/rmu) * xp
      gamp = sqrt(1+zp**2)
      exp_maxwellp = 2d0/(sqrt(1+2*xp**2/rmu) + 1)
  
      call karney(zp, gamp,j0_1p,j0_11p,j0_2p,j0_02p,j0_022p,j1_0p,j1_1p,j1_2p,j1_02p,j1_11p,j1_022p)
      
      F_u_b1 = ((gam**2 * j0_1p - 2*j0_11p) * xp**2/(x**2*gamp)) * exp(-exp_maxwellp * xp**2) * phi_exp(mp,xp)
    end function F_u_b1

    function F_u_b2(xp)
      real(kind=dp) :: xp, F_u_b2
      real(kind=dp) :: zp, gamp, exp_maxwellp
      
      zp  = sqrt(2/rmu) * xp
      gamp = sqrt(1+zp**2)
      exp_maxwellp = 2d0/(sqrt(1+2*xp**2/rmu) + 1)
      
      F_u_b2 = (4d0*xp/x * j0_02) * exp(-exp_maxwellp * xp**2) * phi_exp(mp,xp)
    end function F_u_b2
    
  end function F_u_b

  function C_I_rel1_mmp(rmu, l, m, mp)
    real(kind=dp) :: rmu, C_I_rel1_mmp
    integer       :: l, m, mp

    if (l .eq. 0) then
       C_I_rel1_mmp = C_I_l0_rel1(rmu, m, mp)
    elseif (l .eq. 1) then 
       C_I_rel1_mmp = 1.5d0 * 3.d0/(pi**1.5d0) * integrate(k, 0d0)
    end if
  contains

    function k(x)
      real(kind=dp) :: x, k

      k = C_I_l1_rel1(rmu, x, m, mp)
      
    end function k

  end function C_I_rel1_mmp

  function C_I_l0_rel1(rmu, m, mp)
    real(kind=dp) :: rmu, C_I_l0_rel1
    integer       :: m, mp

    !write (*,*) F_u_b(rmu, 10d0, 0), D_uu_b(rmu, 10d0, 0)
    !stop
    C_I_l0_rel1 = 3.d0/(8d0*pi) * norm_maxwell * integrate(k, 0d0)

  contains

    function k(x)
      real(kind=dp) :: k, x
      real(kind=dp) :: z, gam, exp_maxwell
      
      z   = sqrt(2/rmu) * x
      gam = sqrt(1+z**2)
      exp_maxwell = 2d0/(sqrt(1+2*x**2/rmu) + 1)
      k = exp(-exp_maxwell * x**2) * x**2 * (phi_prj(m,x) + x*d_phi_prj(m,x)) &
          * ((2d0*x/gam) * D_uu_b(rmu, x, mp) + F_u_b(rmu,x, mp))
           
    end function k
    
  end function C_I_l0_rel1
  
  function C_I_l1_rel1(rmu, x, m, mp)
    real(kind=dp) :: rmu, C_I_l1_rel1
    integer       :: m, mp
    real(kind=dp) :: gam, k_l1, x
    real(kind=dp) :: I1, I2, exp_maxwell
    real(kind=dp) :: j0_1,j0_11,j0_2,j0_02,j0_022,j1_0,j1_1,j1_2,j1_02,j1_11,j1_022
    real(kind=dp) :: z

    z   = sqrt(2d0/rmu) * x
    gam = sqrt(1+z**2)
    call karney(z,gam,j0_1,j0_11,j0_2,j0_02,j0_022,j1_0,j1_1,j1_2,j1_02,j1_11,j1_022)
    
    exp_maxwell = 2d0/(sqrt(1+2*x**2/rmu) + 1)
    I1 = integrate(k1, 0d0, x)
    I2 = integrate(k2, x)
    k_l1 = norm_maxwell * (1d0/gam * exp(-exp_maxwell * x**2) * phi_exp(mp,x) + I1 + I2)

    C_I_l1_rel1 = norm_maxwell * x**3 * exp(-exp_maxwell * x**2) * phi_prj(m,x) * k_l1
    
  contains

    function k1(xp)
      real(kind=dp) :: xp, k1
      real(kind=dp) :: zp, gamp, exp_maxwellp
      real(kind=dp) :: j0_1p,j0_11p,j0_2p,j0_02p,j0_022p,j1_0p,j1_1p,j1_2p,j1_02p,j1_11p,j1_022p
      
      zp  = sqrt(2d0/rmu) * xp
      gamp = sqrt(1+zp**2)
      exp_maxwellp = 2d0/(sqrt(1+2*xp**2/rmu) + 1)
      
      call karney(zp,gamp,j0_1p,j0_11p,j0_2p,j0_02p,j0_022p,j1_0p,j1_1p,j1_2p,j1_02p,j1_11p,j1_022p)

      k1 = sqrt(2d0) * xp**2/(gam*gamp) * (&
           1.0d0/x**2 * (2.0d0/sqrt(rmu)*j1_1p &
           + j1_2p*sqrt(rmu) &
           - 10d0*j1_02p*sqrt(rmu)) &
           + gam/x**2*(-2*j1_1p*sqrt(rmu) &
           + 4d0*j1_11p*sqrt(rmu) &
           + 6d0*sqrt(rmu)**3*j1_02p &
           - 24d0*sqrt(rmu)**3*j1_022p)  &
           + (2d0*j1_0p/sqrt(rmu)) &
           + gam*(4d0*j1_02p*sqrt(rmu)) &
           ) &
           * exp(-exp_maxwellp*xp**2) * phi_exp(mp,xp)
      
    end function k1

    function k2(xp)
      real(kind=dp) :: xp,k2
      real(kind=dp) :: zp, gamp, exp_maxwellp

      zp  = sqrt(2d0/rmu) * xp
      gamp = sqrt(1+zp**2)
      exp_maxwellp = 2d0/(sqrt(1+2*xp**2/rmu) + 1)

      k2 = sqrt(2d0) * xp**2/(gam*gamp) * (&
           1.0d0/xp**2*(2.0d0/sqrt(rmu)*j1_1 &
           + j1_2*sqrt(rmu) &
           - 10d0*j1_02*sqrt(rmu)) &
           + gamp/xp**2*(-2*j1_1*sqrt(rmu) &
           + 4d0*j1_11*sqrt(rmu) &
           + 6d0*sqrt(rmu)**3 * j1_02 & 
           - 24d0*sqrt(rmu)**3 * j1_022)  &
           + (2d0*j1_0/sqrt(rmu)) &
           + gamp*(4d0*j1_02*sqrt(rmu)) &
           ) &
           * exp(-exp_maxwellp*xp**2) * phi_exp(mp,xp)
      
    end function k2
    
  end function C_I_l1_rel1

  function C_I_rel2_mmp(rmu, l, m, mp)
    real(kind=dp) :: rmu, C_I_rel2_mmp
    integer       :: l, m, mp

    C_I_rel2_mmp = (2d0*l+1)/2d0 * 3.d0/(pi**1.5d0) * integrate(kernel_rel2, 0d0)

    !x = 0d0
    !do k = 0, 20000
    !   x = k * 10d0/20000
    !   write (600,*) x, kernel_rel2(x)
    !end do
    !write (*,*) C_I_rel2_mmp
    !stop
  contains

    function kernel_rel2(x)
      real(kind=dp) :: x, kernel_rel2

      kernel_rel2 = C_I_rel2_mmp_kernel(x, rmu, l, m, mp)
      !write (600, *) x, kernel_rel2
      
    end function kernel_rel2
    
  end function C_I_rel2_mmp
  
  function C_I_rel2_mmp_kernel(x, rmu, l, m, mp)
    use rel_kernels_mod
    real(kind=dp) :: x, rmu, C_I_rel2_mmp_kernel
    integer       :: l, m, mp
    real(kind=dp) :: z, I10, I11, I20, I21
    real(kind=dp) :: I1, I2, exp_maxwell

    z   = sqrt(2d0/rmu) * x
    gam = sqrt(1+z**2)
    exp_maxwell = 2d0/(sqrt(1+2*x**2/rmu) + 1)
    
    I10 = integrate(k1, 0d0, x)
    I11 = integrate(k1, x)
    I1  = (phi_prj(m,x) + x*d_phi_prj(m,x)) * x**2 * norm_maxwell**2 * exp(-exp_maxwell * x**2) * (I10 + I11)
    I20 = integrate(k2, 0d0, x)
    I21 = integrate(k2, x)
    I2  = x**3 * phi_prj(m,x) * norm_maxwell**2 * exp(-exp_maxwell * x**2) * (I20 + I21)    
    
    C_I_rel2_mmp_kernel = I1 + I2
    
    contains
    
    function k1(xp)
      real(kind=dp) :: xp, k1
      real(kind=dp) :: R_11, R_10, R_01, R_00, exp
      real(kind=dp) :: zp, gamp, exp_maxwellp

      zp  = sqrt(2d0/rmu) * xp
      gamp = sqrt(1+zp**2)
      exp_maxwellp = 2d0/(sqrt(1+2*xp**2/rmu) + 1)

      call rel_kernels(l,z,zp,R_11,R_10,R_01,R_00)

      k1 = xp**2 * exp(-exp_maxwellp * xp**2) * (sqrt(2d0/rmu) * R_11 * d_phi_exp(mp, xp) + 2d0/rmu * R_10 * phi_exp(mp, xp))
    end function k1

    function k2(xp)
      real(kind=dp) :: xp, k2
      real(kind=dp) :: R_11, R_10, R_01, R_00
      real(kind=dp) :: zp, gamp, exp_maxwellp

      zp  = sqrt(2d0/rmu) * xp
      gamp = sqrt(1+zp**2)
      exp_maxwellp = 2d0/(sqrt(1+2*xp**2/rmu) + 1)

      call rel_kernels(l,z,zp,R_11,R_10,R_01,R_00)

      k2 = xp**2 * exp(-exp_maxwellp * xp**2) * (2d0/rmu * R_01 * d_phi_exp(mp, xp) + (2d0/rmu)**1.5d0 * R_00 * phi_exp(mp, xp))
    end function k2
    
  end function C_I_rel2_mmp_kernel
    
  recursive function integrate_a_to_b(func1d, a, b) result(y)
    real(kind=dp) :: a, b
    real(kind=dp) :: y
    real(kind=dp), dimension(2) :: res_int
    integer :: k, kmax
    real(kind=dp), dimension(128) :: x_inter ! For performance issues is this preallocated with a sufficient number of knots
    logical :: binknots
    real(kind=dp) :: x_sub_low, x_sub_up, x_sub_del ! Split integration domain into sub-intervals
    integer :: k_sub, n_sub ! Split integration domain into sub-intervals
    
    interface  
       function func1d(x)
         use nrtype, only : dp
         real(kind=dp) :: func1d
         real(kind=dp) :: x
       end function func1d
    end interface

!!$    ! For debug
!!$    res_int = fint1d_qag(func1d, a, b, epsabs, epsrel, sw_qag_rule)
!!$    y = res_int(1)
!!$    return
    
    if (allocated(t_vec)) then
       y = 0d0
       kmax = 0
       binknots = .false.
       
       if (a .lt. t_vec(1)) then
          kmax = kmax + 1
          x_inter(kmax) = a
       end if
       
       do k = lbound(t_vec, 1), ubound(t_vec, 1)
          if (a .ge. t_vec(k)) then
             cycle
          elseif (kmax .eq. 0) then
             kmax = kmax + 1
             x_inter(kmax) = a
          end if
          
          if (t_vec(k) .ge. b) then
             kmax = kmax + 1
             x_inter(kmax) = b
             binknots = .true.
             exit
          end if

          kmax = kmax + 1
          x_inter(kmax) = t_vec(k)

       end do

       if (.not. binknots) then
          kmax = kmax + 1
          x_inter(kmax) = b
       end if

       !write (*,*) "***"
       !write (*,*) t_vec
       !write (*,*) a, b
       !write (*,*) x_inter(1:kmax)

       if (lsw_split_interval) then
          do k = 1, kmax-1
             n_sub = num_sub_intervals
             if ((.not. binknots) .and. (k .eq. kmax-1)) n_sub = num_sub_intervals_cutoff
             x_sub_del = (x_inter(k+1) - x_inter(k))/dble(n_sub)
             do k_sub = 1,n_sub
                x_sub_low = x_inter(k) + (k_sub-1) * x_sub_del
                x_sub_up = x_inter(k) + k_sub * x_sub_del
                !
                !res_int = fint1d_qag(func1d, x_sub_low , x_sub_up, epsabs, epsrel, sw_qag_rule)
                !res_int = fint1d_qags(func1d, x_sub_low , x_sub_up, epsabs, epsrel)
                res_int = fint1d_cquad(func1d, x_sub_low , x_sub_up, epsabs, epsrel)
                y = y + res_int(1)
             end do
          end do
       else
          do k = 1, kmax-1
             !res_int = fint1d_qag(func1d, x_inter(k), x_inter(k+1), epsabs, epsrel, sw_qag_rule)
             !res_int = fint1d_qags(func1d, x_inter(k), x_inter(k+1), epsabs, epsrel)
             res_int = fint1d_cquad(func1d, x_inter(k), x_inter(k+1), epsabs, epsrel)
             y = y + res_int(1)
          end do
       end if
    else
       !res_int = fint1d_qag(func1d, a, b, epsabs, epsrel, sw_qag_rule)
       !res_int = fint1d_qags(func1d, a, b, epsabs, epsrel)
       res_int = fint1d_cquad(func1d, a, b, epsabs, epsrel)
       y = res_int(1)
    end if
       
  end function integrate_a_to_b

  recursive function integrate_a_to_b_param(func1d, param, a, b) result(y)
    real(kind=dp) :: a, b, param
    real(kind=dp) :: y
    real(kind=dp), dimension(2) :: res_int
    integer :: k, kmax
    real(kind=dp), dimension(128) :: x_inter ! For performance issues is this preallocated with a sufficient number of knots
    logical :: binknots
    real(kind=dp) :: x_sub_low, x_sub_up, x_sub_del ! Split integration domain into sub-intervals
    integer :: k_sub, n_sub ! Split integration domain into sub-intervals
    
    interface  
       function func1d(x, param1)
         use nrtype, only : dp
         real(kind=dp) :: func1d
         real(kind=dp) :: x, param1
       end function func1d
    end interface

!!$    ! For debug
!!$    res_int = fint1d_qag(func1d, param, a, b, epsabs, epsrel, sw_qag_rule)
!!$    y = res_int(1)
!!$    return
    
    if (allocated(t_vec)) then
       y = 0d0
       kmax = 0
       binknots = .false.
       
       if (a .lt. t_vec(1)) then
          kmax = kmax + 1
          x_inter(kmax) = a
       end if
       
       do k = lbound(t_vec, 1), ubound(t_vec, 1)
          if (a .ge. t_vec(k)) then
             cycle
          elseif (kmax .eq. 0) then
             kmax = kmax + 1
             x_inter(kmax) = a
          end if
          
          if (t_vec(k) .ge. b) then
             kmax = kmax + 1
             x_inter(kmax) = b
             binknots = .true.
             exit
          end if

          kmax = kmax + 1
          x_inter(kmax) = t_vec(k)

       end do

       if (.not. binknots) then
          kmax = kmax + 1
          x_inter(kmax) = b
       end if

       !write (*,*) "***"
       !write (*,*) t_vec
       !write (*,*) a, b
       !write (*,*) x_inter(1:kmax)

       if (lsw_split_interval) then
          do k = 1, kmax-1
             n_sub = num_sub_intervals
             if ((.not. binknots) .and. (k .eq. kmax-1)) n_sub = num_sub_intervals_cutoff
             x_sub_del = (x_inter(k+1) - x_inter(k))/dble(n_sub)
             do k_sub = 1,n_sub
                x_sub_low = x_inter(k) + (k_sub-1) * x_sub_del
                x_sub_up = x_inter(k) + k_sub * x_sub_del
                !
                !res_int = fint1d_qag(func1d, param, x_sub_low , x_sub_up, epsabs, epsrel, sw_qag_rule)
                !res_int = fint1d_qags(func1d, param, x_sub_low , x_sub_up, epsabs, epsrel)
                res_int = fint1d_cquad(func1d, param, x_sub_low , x_sub_up, epsabs, epsrel)

                y = y + res_int(1)
             end do
          end do
       else
          do k = 1, kmax-1
             !res_int = fint1d_qag(func1d, param, x_inter(k), x_inter(k+1), epsabs, epsrel, sw_qag_rule)
             !res_int = fint1d_qags(func1d, param, x_inter(k), x_inter(k+1), epsabs, epsrel)
             res_int = fint1d_cquad(func1d, param, x_inter(k), x_inter(k+1), epsabs, epsrel)
             y = y + res_int(1)
          end do
       end if
    else
       !res_int = fint1d_qag(func1d, param, a, b, epsabs, epsrel, sw_qag_rule)
       !res_int = fint1d_qags(func1d, param, a, b, epsabs, epsrel)
       res_int = fint1d_cquad(func1d, param, a, b, epsabs, epsrel)
       y = res_int(1)
    end if
       
  end function integrate_a_to_b_param

  recursive function integrate_a_to_infinity(func1d, a) result(y)
    real(kind=dp) :: a
    real(kind=dp) :: y
    real(kind=dp), dimension(2) :: res_int
    integer :: k
    logical :: in_interval
    real(kind=dp) :: x_sub_low, x_sub_up, x_sub_del ! Split integration domain into sub-intervals
    integer :: k_sub ! Split integration domain into sub-intervals
    real(kind=dp) :: x_cutoff_local
    
    interface  
       function func1d(x)
         use nrtype, only : dp
         real(kind=dp) :: func1d
         real(kind=dp) :: x
       end function func1d
    end interface

!!$    ! For debug
!!$    if (integral_cutoff) then
!!$       res_int = fint1d_qag(func1d, a, x_cutoff, epsabs, epsrel, sw_qag_rule)
!!$    else
!!$       res_int = fint1d_qagiu(func1d, a, epsabs, epsrel)
!!$    end if
!!$    y = res_int(1)
!!$    return
    if (integral_cutoff) then
      x_cutoff_local = x_cutoff
      do while (abs(func1d(x_cutoff_local-1.0)) < 1.0d-250 .and. (x_cutoff_local - a) > 10.0)
        x_cutoff_local = x_cutoff_local-1.0
      end do
!~       write(*,*) 'Using x_cutoff_local=', x_cutoff_local, 'with |f(x)|=', abs(func1d(x_cutoff_local))
    end if

    if (allocated(t_vec)) then
       ! Assuming that t_vec is ordered, this should be the same as a .le. maxval(t_vec)
       if (a .le. t_vec(ubound(t_vec,1))) then
          y = 0d0
          in_interval = .false.
          !write (*,*) "Integration", a, t_vec

          ! Integrate the part up to maximum of the parallel velocity grid.
          do k = lbound(t_vec, 1), ubound(t_vec, 1)
             if (.not. in_interval) then
                if (a .ge. t_vec(k)) then
                   !write (*,*) "Cycle"
                   cycle
                else
                   if (lsw_split_interval) then
                      x_sub_del = (t_vec(k)-a)/dble(num_sub_intervals)
                      do k_sub = 1,num_sub_intervals
                         x_sub_low = a + (k_sub-1) * x_sub_del
                         x_sub_up = a + k_sub * x_sub_del
                         !
                         !res_int = fint1d_qag(func1d, x_sub_low , x_sub_up, epsabs, epsrel, sw_qag_rule)
                         !res_int = fint1d_qags(func1d, x_sub_low , x_sub_up, epsabs, epsrel)
                         res_int = fint1d_cquad(func1d, x_sub_low , x_sub_up, epsabs, epsrel)

                         y = y + res_int(1)
                      end do
                   else
                      !res_int = fint1d_qag(func1d, a, t_vec(k), epsabs, epsrel, sw_qag_rule)
                      !res_int = fint1d_qags(func1d, a, t_vec(k), epsabs, epsrel)
                      res_int = fint1d_cquad(func1d, a, t_vec(k), epsabs, epsrel)

                      y = y + res_int(1)
                   end if
                   !write (*,*) "1. Int", a, t_vec(k), res_int
                   in_interval = .true.
                end if
             end if
             if (in_interval) then
                if (k .lt. ubound(t_vec,1)) then
                   if (lsw_split_interval) then
                      x_sub_del = (t_vec(k+1) - t_vec(k))/dble(num_sub_intervals)
                      do k_sub = 1,num_sub_intervals
                         x_sub_low = t_vec(k) + (k_sub-1) * x_sub_del
                         x_sub_up = t_vec(k) + k_sub * x_sub_del
                         !
                         !res_int = fint1d_qag(func1d, x_sub_low , x_sub_up, epsabs, epsrel, sw_qag_rule)
                         !res_int = fint1d_qags(func1d, x_sub_low , x_sub_up, epsabs, epsrel)
                         res_int = fint1d_cquad(func1d, x_sub_low , x_sub_up, epsabs, epsrel)

                         y = y + res_int(1)
                      end do
                   else
                      !res_int = fint1d_qag(func1d, t_vec(k), t_vec(k+1), epsabs, epsrel, sw_qag_rule)
                      !res_int = fint1d_qags(func1d, t_vec(k), t_vec(k+1), epsabs, epsrel)
                      res_int = fint1d_cquad(func1d, t_vec(k), t_vec(k+1), epsabs, epsrel)

                      y = y + res_int(1)
                   end if
                   !write (*,*) "Rest int", t_vec(k), t_vec(k+1), y
                end if
             end if

          end do

          ! Integrate from the maximum of the parallel velocity grid to
          ! infinity (or at least some high enough x value).
          if (integral_cutoff) then
             if (lsw_split_interval) then
                x_sub_del = (x_cutoff_local - t_vec(ubound(t_vec,1)))/dble(num_sub_intervals_cutoff)
                do k_sub = 1,num_sub_intervals_cutoff
                   x_sub_low = t_vec(ubound(t_vec,1)) + (k_sub-1) * x_sub_del
                   x_sub_up = t_vec(ubound(t_vec,1)) + k_sub * x_sub_del
                   !
                   !res_int = fint1d_qag(func1d, x_sub_low , x_sub_up, epsabs, epsrel, sw_qag_rule)
                   !res_int = fint1d_qags(func1d, x_sub_low , x_sub_up, epsabs, epsrel)
                   res_int = fint1d_cquad(func1d, x_sub_low , x_sub_up, epsabs, epsrel)

                   y = y + res_int(1)
                end do
             else
                !res_int = fint1d_qag(func1d, t_vec(ubound(t_vec,1)), x_cutoff_local, epsabs, epsrel, sw_qag_rule)
                !res_int = fint1d_qags(func1d, t_vec(ubound(t_vec,1)), x_cutoff_local, epsabs, epsrel)
                res_int = fint1d_cquad(func1d, t_vec(ubound(t_vec,1)), x_cutoff_local, epsabs, epsrel)

                y = y + res_int(1)
             end if
          else
             res_int = fint1d_qagiu(func1d, t_vec(ubound(t_vec,1)), epsabs, epsrel)
             !res_int = fint1d_qags(func1d, t_vec(ubound(t_vec,1)), x_cutoff_local, epsabs, epsrel)
             !res_int = fint1d_cquad(func1d, t_vec(ubound(t_vec,1)), x_cutoff_local, epsabs, epsrel)

             y = y + res_int(1)
          end if
          !write (*,*) "Final int", t_vec(ubound(t_vec,1)), y
       else
          !write (*,*) "int a outside domain", a, integral_cutoff, x_cutoff_local
          if (integral_cutoff) then
             if (lsw_split_interval) then
                y = 0d0
                x_sub_del = (x_cutoff_local - a)/dble(num_sub_intervals_cutoff)
                do k_sub = 1,num_sub_intervals_cutoff
                   x_sub_low = a + (k_sub-1) * x_sub_del
                   x_sub_up = a + k_sub * x_sub_del
                   !
                   !res_int = fint1d_qag(func1d, x_sub_low , x_sub_up, epsabs, epsrel, sw_qag_rule)
                   !res_int = fint1d_qags(func1d, x_sub_low , x_sub_up, epsabs, epsrel)
                   res_int = fint1d_cquad(func1d, x_sub_low , x_sub_up, epsabs, epsrel)

                   y = y + res_int(1)
                end do
             else
                !res_int = fint1d_qag(func1d, a, x_cutoff_local, epsabs, epsrel, sw_qag_rule)
                !res_int = fint1d_qags(func1d, a, x_cutoff_local, epsabs, epsrel)
                res_int = fint1d_cquad(func1d, a, x_cutoff_local, epsabs, epsrel)

                y = res_int(1)
             end if
          else
             res_int = fint1d_qagiu(func1d, a, epsabs, epsrel)
             !res_int = fint1d_qags(func1d, a, x_cutoff_local, epsabs, epsrel)
             !res_int = fint1d_cquad(func1d, a, x_cutoff_local, epsabs, epsrel)

             y = res_int(1)
          end if
       end if
       
    else
       
       if (integral_cutoff) then
          !res_int = fint1d_qag(func1d, a, x_cutoff_local, epsabs, epsrel, sw_qag_rule)
          !res_int = fint1d_qags(func1d, a, x_cutoff_local, epsabs, epsrel)
          res_int = fint1d_cquad(func1d, a, x_cutoff_local, epsabs, epsrel)

       else
          res_int = fint1d_qagiu(func1d, a, epsabs, epsrel)
          !res_int = fint1d_qags(func1d, a, x_cutoff_local, epsabs, epsrel)
          !res_int = fint1d_cquad(func1d, a, x_cutoff_local, epsabs, epsrel)
       end if
       y = res_int(1)
       
    end if
    
  end function integrate_a_to_infinity

  recursive function integrate_a_to_infinity_param(func1d, param, a) result(y)
    real(kind=dp) :: a, param
    real(kind=dp) :: y
    real(kind=dp), dimension(2) :: res_int
    integer :: k
    logical :: in_interval
    real(kind=dp) :: x_sub_low, x_sub_up, x_sub_del ! Split integration domain into sub-intervals
    integer :: k_sub ! Split integration domain into sub-intervals
    real(kind=dp) :: x_cutoff_local
    
    interface  
       function func1d(x, param1)
         use nrtype, only : dp
         real(kind=dp) :: func1d
         real(kind=dp) :: x, param1
       end function func1d
    end interface

!!$    ! For debug
!!$    if (integral_cutoff) then
!!$       res_int = fint1d_qag(func1d, param, a, x_cutoff, epsabs, epsrel, sw_qag_rule)
!!$    else
!!$       res_int = fint1d_qagiu(func1d, param, a, epsabs, epsrel)
!!$    end if
!!$    y = res_int(1)
!!$    return
    if (integral_cutoff) then
      x_cutoff_local = x_cutoff
      do while (abs(func1d(x_cutoff_local-1.0, param)) < 1.0d-250 .and. (x_cutoff_local - a) > 10.0)
        x_cutoff_local = x_cutoff_local-1.0
      end do
!~       write(*,*) 'Using x_cutoff_local=', x_cutoff_local, 'with |f(x)|=', abs(func1d(x_cutoff_local))
    end if

    if (allocated(t_vec)) then
       
       if (a .le. t_vec(ubound(t_vec,1))) then
          y = 0d0
          in_interval = .false.
          !write (*,*) "Integration", a, t_vec
          do k = lbound(t_vec, 1), ubound(t_vec, 1)
             if (.not. in_interval) then
                if (a .ge. t_vec(k)) then
                   !write (*,*) "Cycle"
                   cycle
                else
                   if (lsw_split_interval) then
                      x_sub_del = (t_vec(k)-a)/dble(num_sub_intervals)
                      do k_sub = 1,num_sub_intervals
                         x_sub_low = a + (k_sub-1) * x_sub_del
                         x_sub_up = a + k_sub * x_sub_del
                         !
                         !res_int = fint1d_qag(func1d, param, x_sub_low , x_sub_up, epsabs, epsrel, sw_qag_rule)
                         !res_int = fint1d_qags(func1d, param, x_sub_low , x_sub_up, epsabs, epsrel)
                         res_int = fint1d_cquad(func1d, param, x_sub_low , x_sub_up, epsabs, epsrel)

                         y = y + res_int(1)
                      end do
                   else
                      !res_int = fint1d_qag(func1d, param, a, t_vec(k), epsabs, epsrel, sw_qag_rule)
                      !res_int = fint1d_qags(func1d, param, a, t_vec(k), epsabs, epsrel)
                      res_int = fint1d_cquad(func1d, param, a, t_vec(k), epsabs, epsrel)

                      y = y + res_int(1)
                   end if
                   !write (*,*) "1. Int", a, t_vec(k), res_int
                   in_interval = .true.
                end if
             end if
             if (in_interval) then
                if (k .lt. ubound(t_vec,1)) then
                   if (lsw_split_interval) then
                      x_sub_del = (t_vec(k+1) - t_vec(k))/dble(num_sub_intervals)
                      do k_sub = 1,num_sub_intervals
                         x_sub_low = t_vec(k) + (k_sub-1) * x_sub_del
                         x_sub_up = t_vec(k) + k_sub * x_sub_del
                         !
                         !res_int = fint1d_qag(func1d, param, x_sub_low , x_sub_up, epsabs, epsrel, sw_qag_rule)
                         !res_int = fint1d_qags(func1d, param, x_sub_low , x_sub_up, epsabs, epsrel)
                         res_int = fint1d_cquad(func1d, param, x_sub_low , x_sub_up, epsabs, epsrel)

                         y = y + res_int(1)
                      end do
                   else
                      !res_int = fint1d_qag(func1d, param, t_vec(k), t_vec(k+1), epsabs, epsrel, sw_qag_rule)
                      !res_int = fint1d_qags(func1d, param, t_vec(k), t_vec(k+1), epsabs, epsrel)
                      res_int = fint1d_cquad(func1d, param, t_vec(k), t_vec(k+1), epsabs, epsrel)

                      y = y + res_int(1)
                   end if
                   !write (*,*) "Rest int", t_vec(k), t_vec(k+1), y
                end if
             end if

          end do
          if (integral_cutoff) then
             if (lsw_split_interval) then
                x_sub_del = (x_cutoff_local - t_vec(ubound(t_vec,1)))/dble(num_sub_intervals_cutoff)
                do k_sub = 1,num_sub_intervals_cutoff
                   x_sub_low = t_vec(ubound(t_vec,1)) + (k_sub-1) * x_sub_del
                   x_sub_up = t_vec(ubound(t_vec,1)) + k_sub * x_sub_del
                   !
                   !res_int = fint1d_qag(func1d, param, x_sub_low , x_sub_up, epsabs, epsrel, sw_qag_rule)
                   !res_int = fint1d_qags(func1d, param, x_sub_low , x_sub_up, epsabs, epsrel)
                   res_int = fint1d_cquad(func1d, param, x_sub_low , x_sub_up, epsabs, epsrel)

                   y = y + res_int(1)
                end do
             else
                !res_int = fint1d_qag(func1d, param, t_vec(ubound(t_vec,1)), x_cutoff_local, epsabs, epsrel, sw_qag_rule)
                !res_int = fint1d_qags(func1d, param, t_vec(ubound(t_vec,1)), x_cutoff_local, epsabs, epsrel)
                res_int = fint1d_cquad(func1d, param, t_vec(ubound(t_vec,1)), x_cutoff_local, epsabs, epsrel)

                y = y + res_int(1)
             end if
          else
             res_int = fint1d_qagiu(func1d, param, t_vec(ubound(t_vec,1)), epsabs, epsrel)
             !res_int = fint1d_qags(func1d, param, t_vec(ubound(t_vec,1)), x_cutoff_local, epsabs, epsrel)
             !res_int = fint1d_cquad(func1d, param, t_vec(ubound(t_vec,1)), x_cutoff_local, epsabs, epsrel)

             y = y + res_int(1)
          end if
          !write (*,*) "Final int", t_vec(ubound(t_vec,1)), y
       else
          !write (*,*) "int a outside domain", a, integral_cutoff, x_cutoff_local
          if (integral_cutoff) then
             if (lsw_split_interval) then
                y = 0d0
                x_sub_del = (x_cutoff_local - a)/dble(num_sub_intervals_cutoff)
                do k_sub = 1,num_sub_intervals_cutoff
                   x_sub_low = a + (k_sub-1) * x_sub_del
                   x_sub_up = a + k_sub * x_sub_del
                   !
                   !res_int = fint1d_qag(func1d, param, x_sub_low , x_sub_up, epsabs, epsrel, sw_qag_rule)
                   !res_int = fint1d_qags(func1d, param, x_sub_low , x_sub_up, epsabs, epsrel)
                   res_int = fint1d_cquad(func1d, param, x_sub_low , x_sub_up, epsabs, epsrel)

                   y = y + res_int(1)
                end do
             else
                !res_int = fint1d_qag(func1d, param, a, x_cutoff_local, epsabs, epsrel, sw_qag_rule)
                !res_int = fint1d_qags(func1d, param, a, x_cutoff_local, epsabs, epsrel)
                res_int = fint1d_cquad(func1d, param, a, x_cutoff_local, epsabs, epsrel)

                y = res_int(1)
             end if
          else
             res_int = fint1d_qagiu(func1d, param, a, epsabs, epsrel)
             !res_int = fint1d_qags(func1d, param, a, x_cutoff_local, epsabs, epsrel)
             !res_int = fint1d_cquad(func1d, param, a, x_cutoff_local, epsabs, epsrel)

             y = res_int(1)
          end if
       end if
       
    else
       
       if (integral_cutoff) then
          !res_int = fint1d_qag(func1d, param, a, x_cutoff_local, epsabs, epsrel, sw_qag_rule)
          !res_int = fint1d_qags(func1d, param, a, x_cutoff_local, epsabs, epsrel)
          res_int = fint1d_cquad(func1d, param, a, x_cutoff_local, epsabs, epsrel)

       else
          res_int = fint1d_qagiu(func1d, param, a, epsabs, epsrel)
          !res_int = fint1d_qags(func1d, param, a, x_cutoff_local, epsabs, epsrel)
          !res_int = fint1d_cquad(func1d, param, a, x_cutoff_local, epsabs, epsrel)

       end if
       y = res_int(1)
       
    end if
    
  end function integrate_a_to_infinity_param

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
  end subroutine init_legendre

  function G(x) result(y)
    real(kind=dp) :: x
    real(kind=dp) :: y

    if (lsw_expand_G .and. (x .le. xmax_expand_kernel)) then
       ! use 5th order Taylor expansion around x=0 (numerical stability)
       y = 2d0*x/(3d0*sqrt(pi))-2d0*(x**3)/(5d0*sqrt(pi))+(x**5)/(7d0*sqrt(pi))
    else
      ! use default expression

      !> According to the taylor expansion, the function behaves smoothly
      !> near zero, thus return the limit for x = 0.
      if (x .eq. 0.0d0) then
        y = 0.0d0
      else
        y = (erf(x) - x*d_erf(x)) / (2*x**2)
      end if
    end if

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

    !> According to the taylor expansion from G(x), the derivative
    !> goes smoothly to the value \[ 2/(3\sqrt(\pi) \], despite being
    !> not defined. Thus if x is zero, simply return this limit.
    if (x .eq. 0.0d0) then
      y = 2d0/(3d0*sqrt(pi))
    else
      y = (4*x**2*d_erf(x) - 4*x*erf(x) - 2*x**3*dd_erf(x))
      y = y / (2*x**2)**2
    end if
  end function d_G

  function nu_D_hat(x) result(nu)
    real(kind=dp) :: x, y, nu

    ! At the momement only for self-collisions !!!!!!!
    y = x * 1d0

    nu = 3*sqrt(pi)/4d0 * (erf(y) - G(y))/(x**3)
  end function nu_D_hat

  subroutine inv(A)
    real(dp), dimension(:,:), intent(inout) :: A
    integer,  dimension(size(A,1)) :: ipiv
    real(dp), dimension(size(A,1)) :: work  
    integer :: n, info

    n = size(A,1)

    !**********************************************************
    ! LU Factorization (LAPACK)
    !**********************************************************
    call DGETRF(n, n, A, n, ipiv, info)

    !**********************************************************
    ! Check state
    !**********************************************************
    if (info /= 0) then
       stop 'Error: Matrix is numerically singular!'
    end if

    !**********************************************************
    ! Compute inverse matrix (LAPACK)
    !**********************************************************
    call DGETRI(n, A, n, ipiv, work, n, info)

    !**********************************************************
    ! Check state
    !**********************************************************
    if (info /= 0) then
       stop 'Error: Matrix inversion failed!'
    end if
  end subroutine inv

  subroutine compute_Minv(M)
    real(kind=dp), dimension(:,:), intent(out) :: M
    integer :: mm, mp

    write (*,*) "Computing phi transformation matrix..."
    if (precomp) then ! load pre-computed M_transform
       call h5_open(trim(adjustl(matelem_name)), h5id_matelem)
       call h5_get(h5id_matelem,'Ammp',M)
       call h5_close(h5id_matelem)
       !DO mm=1,SIZE(Minv,1)
       !   PRINT *,(Minv(mm,mp),mp=1,SIZE(Minv,2))
       !END DO
       !STOP
    else
       if (isw_relativistic .eq. 0) then
          do mm = 0, lagmax
             do mp = 0, lagmax
                M(mm+1,mp+1) = integrate(phim_phimp, 0d0)
             end do
          end do
       elseif (isw_relativistic .ge. 1) then
          do mm = 0, lagmax
             do mp = 0, lagmax
                M(mm+1,mp+1) = norm_maxwell * integrate(phim_phimp_rel1, 0d0)
                !M(mm+1,mp+1) = integrate(phim_phimp, 0d0)
             end do
          end do          
       end if
    end if

    if (make_ortho) then ! make DKE orthogonal w.r.t. to derivative along field line
      M_transform = 0.0d0
      do mm = 0, lagmax
        M_transform(mm, mm) = 1.0d0
      end do
    else
       M_transform = M
    end if

    
    ! DEBUG - MAKE MATRIX ALWAYS UNIT MATRIX
!!$    write (*,*) "WARNING, INVERSE MATRIX NOT COMPUTED!!!!!!"
!!$    do mm = 0, lagmax
!!$       do mp = 0, lagmax
!!$          if (mm .eq. mp) then
!!$             M(mm+1,mp+1) = 1d0
!!$          else
!!$             M(mm+1,mp+1) = 0d0
!!$          end if
!!$       end do
!!$    end do

    
    !**********************************************************
    ! M -> inv(M)
    !**********************************************************
    call inv(M)
    call chop(M)

  contains

    function phim_phimp(x)
      real(kind=dp) :: x, phim_phimp

      phim_phimp =  pi**(-3d0/2d0) * x**(4+alpha) * exp(-(1+beta)*x**2) * phi_prj(mm,x) * phi_exp(mp,x)

    end function phim_phimp

    function phim_phimp_rel1(x)
      real(kind=dp) :: x, phim_phimp_rel1
      real(kind=dp) :: exp_maxwell, gam
      
      exp_maxwell = 2d0/(sqrt(1+2*x**2/rmu) + 1)
      gam = sqrt(1 + 2*x**2/rmu)

      !phim_phimp_rel1 =  pi**(-3d0/2d0) * x**(4+alpha) * exp(-(1+beta)*x**2) * phi_prj(mm,x) * phi_exp(mp,x)
      phim_phimp_rel1 = 1d0/gam * pi**(-3d0/2d0) * x**(4+alpha) * exp(-(exp_maxwell+beta)*x**2) * phi_prj(mm,x) * phi_exp(mp,x)

    end function phim_phimp_rel1
  end subroutine compute_Minv

  subroutine compute_sources(asource_s, weightlag_s, weightden_s, weightparflow_s, weightenerg_s)
    real(kind=dp), dimension(:,:), intent(out) :: asource_s, weightlag_s
    real(kind=dp), dimension(:), intent(out)   :: weightden_s, weightparflow_s, weightenerg_s
    real(kind=dp) :: res_int
    integer :: m, k, j

    write (*,*) "Computing sources..."

    !if (allocated(asource_s)) deallocate(asource_s)
    !allocate(asource_s(0:lagmax, 1:3))
    !write (*,*) lbound(asource_s), ubound(asource_s)
    if (precomp) then ! load pre-computed sources
       call h5_open(trim(adjustl(matelem_name)), h5id_matelem)
       call h5_get(h5id_matelem,'a1m',asource_s(:,1))
       call h5_get(h5id_matelem,'a2m',asource_s(:,2))
       call h5_get(h5id_matelem,'a3m',asource_s(:,3))
       call h5_close(h5id_matelem)
       !do m=1,size(asource_s,1)
       !   print *,asource_s(m,1),asource_s(m,2),asource_s(m,3)
       !end do
       !stop
    else
       if (isw_relativistic .eq. 0) then 
          do k = 1, 3
             do m = 0, lagmax
                asource_s(m+1, k) = integrate(am, 0d0)
             end do
          end do
       elseif (isw_relativistic .ge.1) then
          do k = 1, 3
             do m = 0, lagmax
                asource_s(m+1, k) = norm_maxwell * integrate(am_rel1, 0d0)
             end do
          end do          
       end if
    end if

    if (make_ortho) then ! make DKE orthogonal w.r.t. to derivative along field line
       do k = 1, 3
          asource_s(:,k) = matmul(M_transform_inv, asource_s(:,k))
       end do
    end if

    call chop(asource_s)
    write (*,*) "Computing weighting coefficients..."

    if (allocated(C_m)) deallocate(C_m)
    allocate(C_m(0:lagmax))

    if (precomp) then ! load pre-computed weighting coefficients
       call h5_open(trim(adjustl(matelem_name)), h5id_matelem)
       call h5_get(h5id_matelem,'b1m',weightlag_s(1,:))
       call h5_get(h5id_matelem,'b2m',weightlag_s(2,:))
       call h5_get(h5id_matelem,'b3m',weightlag_s(3,:))
       call h5_close(h5id_matelem)
       !do m=1,size(weightlag_s,2)
       !   PRINT *,weightlag_s(1,m),weightlag_s(2,m),weightlag_s(3,m)
       !end do
       !stop
    else
       if (isw_relativistic .eq. 0 ) then
          do j = 1, 3
             do m = 0, lagmax
                res_int = integrate(bm, 0d0)
                weightlag_s(j,m+1) = 1d0/sqrt(pi) * res_int
                if (j .eq. 3) C_m(m) = res_int
             end do
          end do
       elseif (isw_relativistic .ge. 1) then
          do j = 1, 3
             do m = 0, lagmax
                res_int = norm_maxwell * integrate(bm_rel1, 0d0)
                weightlag_s(j,m+1) = 1d0/sqrt(pi) * res_int
                if (j .eq. 3) C_m(m) = res_int
                !write (*,*) "weightlag: ", weightlag_s(j,m+1)
             end do
          end do
       end if
    end if

    ! remove null-space of axisymmetric solution (tokamak mode, energy conservation)
    do m = 0, lagmax
       res_int = integrate(weightenerg_kernel, 0d0)
       weightenerg_s(m+1) = res_int
    end do

    if (make_ortho) then
       weightenerg_s = matmul(M_transform_inv, weightenerg_s)
    end if
    
    ! Was used for pre-conditioned iterations (not necessary/depricated)
    ! weightparflow for computation of bvec_parflow (must be orthogonal?)
    if (make_ortho) then
       weightparflow_s = asource_s(:,1)
    else
       weightparflow_s = matmul(M_transform_inv, asource_s(:,1))
    end if
    
    call chop(weightlag_s)

    if (isw_relativistic .eq. 0) then
       do m = 0, lagmax
          weightden_s(m+1) = 1d0/sqrt(pi) * integrate(weightden_kernel, 0d0)
          weightenerg_offset(m) = 1d0/sqrt(pi) * integrate(weightenerg_offset_kernel, 0d0)
       end do
    elseif (isw_relativistic .ge. 1) then
       do m = 0, lagmax
          weightden_s(m+1) = norm_maxwell * 1d0/sqrt(pi) * integrate(weightden_kernel_rel1, 0d0)
          weightenerg_offset(m) = norm_maxwell * 1d0/sqrt(pi) * &
               integrate(weightenerg_offset_kernel_rel1, 0d0)       
       end do       
    end if

  contains

    function am(x)
      real(kind=dp) :: x, am

      am = pi**(-3d0/2d0) * x**(4+alpha) * exp(-(1+beta)*x**2) * phi_prj(m, x) * x**(2*k - 1 - 5*kdelta(3,k))
    end function am

    function am_rel1(x)
      real(kind=dp) :: x, am_rel1
      real(kind=dp) :: exp_maxwell, gam, gam_fac

      exp_maxwell = 2d0/(sqrt(1+2*x**2/rmu) + 1)
      gam = sqrt(1+2*x**2/rmu)
     
      if (k .ne. 2) then
         gam_fac = 1.0d0/gam
      else
         gam_fac = 1.d0/gam * 2.0d0/(gam+1d0)
      end if
      
      am_rel1 = gam_fac * pi**(-3d0/2d0) * x**(4+alpha) * exp(-(exp_maxwell+beta)*x**2) &
           * phi_prj(m, x) * x**(2*k - 1 - 5*kdelta(3,k))
    end function am_rel1

    function bm(x)
      real(kind=dp) :: x, bm

      bm = exp(-x**2) * x**(2*(j+1)-5*kdelta(3,j)) * phi_exp(m,x)
    end function bm

    function bm_rel1(x)
      real(kind=dp) :: x, bm_rel1
      real(kind=dp) :: exp_maxwell, gam
      real(kind=dp) :: gam_fac
      
      exp_maxwell = 2d0/(sqrt(1+2*x**2/rmu) + 1)
      gam = sqrt(1+2*x**2/rmu)
      
      if (j .ne. 2) then
         gam_fac = 1.0d0/gam
      else
         gam_fac = 1.0d0/gam * 2.0d0/(gam+1d0)
      end if
      
      bm_rel1 = gam_fac * exp(-exp_maxwell*x**2) * x**(2*(j+1)-5*kdelta(3,j)) * phi_exp(m,x)
    end function bm_rel1

    function weightenerg_kernel(x)
      real(kind=dp) :: x, weightenerg_kernel

      weightenerg_kernel = &
           pi**(-3d0/2d0) * x**(4+alpha) * exp(-(1+beta)*x**2) * phi_prj(m, x) * (x**2 - 1.5d0)
    end function weightenerg_kernel
    
    function weightden_kernel(x)
      real(kind=dp) :: x, weightden_kernel
      weightden_kernel = exp(-x**2) * x**2 * phi_exp(m,x)
    end function weightden_kernel

    function weightden_kernel_rel1(x)
      real(kind=dp) :: x, weightden_kernel_rel1
      real(kind=dp) :: exp_maxwell, gam

      exp_maxwell = 2d0/(sqrt(1+2*x**2/rmu) + 1)
      gam = sqrt(1+2*x**2/rmu)
      weightden_kernel_rel1 =   exp(-exp_maxwell*x**2) * x**2 * phi_exp(m,x)
    end function weightden_kernel_rel1

    function weightenerg_offset_kernel(x)
      real(kind=dp) :: x, weightenerg_offset_kernel

      weightenerg_offset_kernel = exp(-x**2) * x**4 * phi_exp(m,x)  
    end function weightenerg_offset_kernel
    
    function weightenerg_offset_kernel_rel1(x)
      real(kind=dp) :: x, weightenerg_offset_kernel_rel1
      real(kind=dp) :: exp_maxwell, gam

      exp_maxwell = 2d0/(sqrt(1+2*x**2/rmu) + 1)
      gam = sqrt(1+2*x**2/rmu)
      weightenerg_offset_kernel_rel1 =   2.0d0/(gam+1d0) * exp(-exp_maxwell*x**2) * x**4 * phi_exp(m,x)  
    end function weightenerg_offset_kernel_rel1
    
    function kdelta(a,b)
      integer :: a, b
      integer :: kdelta

      kdelta = 0
      if (a .eq. b) kdelta = 1
    end function kdelta

  end subroutine compute_sources

  subroutine compute_xmmp(x1mm_s,x2mm_s)
    real(kind=dp), dimension(:,:) :: x1mm_s, x2mm_s
    real(kind=dp) :: cnorm
    integer :: m, mp

    write(*,*) "Computing x1mm and x2mm ..."
    if (precomp) then ! load pre-computed x1mm / x2mm
       call h5_open(trim(adjustl(matelem_name)), h5id_matelem)
       call h5_get(h5id_matelem,'x1mmp',x1mm_s(:,:))
       call h5_get(h5id_matelem,'x2mmp',x2mm_s(:,:))
       call h5_close(h5id_matelem)
       !do m=1,size(x1mm_s,1)
       !   print *,(x1mm_s(m,mp),mp=1,size(x1mm_s,2))
       !end do
       !print *,' '
       !do m=1,size(x2mm_s,1)
       !   print *,(x2mm_s(m,mp),mp=1,size(x2mm_s,2))
       !end do
       !stop
    else
       cnorm=1.d0/(pi*sqrt(pi))

       do m = 0, lagmax
          do mp = 0, lagmax
             x1mm_s(m+1, mp+1) = cnorm * integrate(integrand_x1mmp, 0d0)
             x2mm_s(m+1, mp+1) = cnorm * integrate(integrand_x2mmp, 0d0)
          end do
       end do
    end if

    if (make_ortho) then ! make DKE orthogonal w.r.t. to derivative along field line
       x1mm_s = matmul(M_transform_inv, x1mm_s)
       x2mm_s = matmul(M_transform_inv, x2mm_s)
    end if

    call chop(x1mm_s)
    call chop(x2mm_s)

  contains

    function integrand_x1mmp(x)
      real(kind=dp) :: x
      real(kind=dp) :: integrand_x1mmp

      integrand_x1mmp= x**(alpha) * exp(-(1+beta)*x**2) * &
           x**(3) * phi_prj(m, x) * phi_exp(mp, x)

    end function integrand_x1mmp

    function integrand_x2mmp(x)
      real(kind=dp) :: x
      real(kind=dp) :: integrand_x2mmp

      integrand_x2mmp= x**(alpha) * exp(-(1+beta)*x**2) * &
           x**(5) * phi_prj(m, x) * phi_exp(mp, x)

    end function integrand_x2mmp

  end subroutine compute_xmmp
  
  subroutine compute_lorentz(anumm_s)
    real(kind=dp), dimension(:,:) :: anumm_s
    integer :: m, mp

    !if (allocated(anumm_s)) deallocate(anumm_s)
    !allocate(anumm_s(0:lagmax, 0:lagmax))

    write (*,*) "Computing Lorentz part..."

    if (precomp) then ! load pre-computed Lorentz part
       call h5_open(trim(adjustl(matelem_name)), h5id_matelem)
       call h5_get(h5id_matelem,'NuabHat',anumm_s(:,:))
       call h5_close(h5id_matelem)
       !do m=1,size(anumm_s,1)
       !   print *,(anumm_s(m,mp),mp=1,size(anumm_s,2))
       !end do
       !stop
    else
       if (isw_relativistic .eq. 0) then
          do m = 0, lagmax
             do mp = 0, lagmax
                anumm_s(m+1, mp+1) = 3d0/(4d0 * pi) * integrate(integrand, 0d0)
             end do
          end do
       elseif (isw_relativistic .ge. 1) then
          do m = 0, lagmax
             do mp = 0, lagmax
                anumm_s(m+1, mp+1) = norm_maxwell * 6d0/(4d0 * pi) * integrate(integrand_rel1, 0d0)
             end do
          end do          
       end if
    end if
  
    if (make_ortho) then ! make DKE orthogonal w.r.t. to derivative along field line
       anumm_s = matmul(M_transform_inv, anumm_s)
    end if

    call chop(anumm_s)
    !write (*,*) "Done."
    
  contains

    function integrand(x)
      real(kind=dp) :: x, y
      real(kind=dp) :: integrand

      y = x * gamma_ab
      integrand = x**(alpha) * exp(-(1+beta)*x**2) * phi_prj(m, x) * (erf(y) - G(y)) * phi_exp(mp, x)

    end function integrand

    function integrand_rel1(x)
      real(kind=dp) :: x
      real(kind=dp) :: integrand_rel1  
      real(kind=dp) :: exp_maxwell, z, gam

      z   = sqrt(2/rmu) * x
      gam = sqrt(1+z**2)      
      exp_maxwell = 2d0/(sqrt(1+2*x**2/rmu) + 1)
      integrand_rel1 = exp(-exp_maxwell * x**2) * x * D_thth(rmu, x) * phi_prj(m,x) * phi_exp(mp, x)

    end function integrand_rel1
    
  end subroutine compute_lorentz

  subroutine compute_lorentz_inf(anumm_s)
    real(kind=dp), dimension(:,:) :: anumm_s
    integer :: m, mp

    !if (allocated(anumm_s)) deallocate(anumm_s)
    !allocate(anumm_s(0:lagmax, 0:lagmax))

    write (*,*) "Computing Lorentz part (gamma -> inf)..."

    if (isw_relativistic .eq. 0) then
       do m = 0, lagmax
          do mp = 0, lagmax
             anumm_s(m+1, mp+1) = 3d0/(4d0 * pi) * integrate(integrand_inf, 0d0)
          end do
       end do
    elseif (isw_relativistic .ge. 1) then
       do m = 0, lagmax
          do mp = 0, lagmax
             anumm_s(m+1, mp+1) = 3d0/(4d0 * pi) * norm_maxwell * integrate(integrand_inf_rel1, 0d0)
          end do
       end do
    end if
    
    if (make_ortho) then ! make DKE orthogonal w.r.t. to derivative along field line
       anumm_s = matmul(M_transform_inv, anumm_s)
    end if

    call chop(anumm_s)

    !write (*,*) "Done."

  contains

    function integrand_inf(x)
      real(kind=dp) :: x, y
      real(kind=dp) :: integrand_inf

      y = x * gamma_ab
      integrand_inf = x**(alpha) * exp(-(1+beta)*x**2) * phi_prj(m, x) * phi_exp(mp, x)

    end function integrand_inf
    
    function integrand_inf_rel1(x)
      real(kind=dp) :: x
      real(kind=dp) :: integrand_inf_rel1
      real(kind=dp) :: exp_maxwell, z, gam

      z   = sqrt(2/rmu) * x
      gam = sqrt(1+z**2)
      exp_maxwell = 2d0/(sqrt(1+2*x**2/rmu) + 1)
      integrand_inf_rel1 = gam * exp(-exp_maxwell*x**2) * phi_prj(m, x) * phi_exp(mp, x)
      
    end function integrand_inf_rel1
    
  end subroutine compute_lorentz_inf

  subroutine compute_energyscattering(denmm_s)
    real(kind=dp), dimension(:,:) :: denmm_s
    integer :: m, mp

    write (*,*) "Computing energy scattering part..."

    if (precomp) then ! load pre-computed energy scattering part
       call h5_open(trim(adjustl(matelem_name)), h5id_matelem)
       call h5_get(h5id_matelem,'DabHat',denmm_s(:,:))
       call h5_close(h5id_matelem)
       !do m=1,size(denmm_s,1)
       !   print *,(denmm_s(m,mp),mp=1,size(denmm_s,2))
       !end do
       !stop
    else
       if (isw_relativistic .eq. 0) then
          do m = 0, lagmax
             do mp = 0, lagmax
                denmm_s(m+1, mp+1) = 3d0/(4d0 * pi) * integrate(integrand, 0d0)
                !write (*,*) denmm_s(m+1, mp+1)
                !stop
                !write (*,*) "denmm_s", m, mp, integrand(2d0), denmm_s(m, mp)
                !write (*,*) G(2d0), d_G(2d0), gamma_ab, phi(m, 2d0), d_phi(m, 2d0), dd_phi(m,2d0)
             end do
          end do
       elseif (isw_relativistic .ge. 1) then
          do m = 0, lagmax
             do mp = 0, lagmax
                denmm_s(m+1, mp+1) = -3d0/(4d0 * pi) * norm_maxwell * integrate(integrand_rel1, 0d0)
                !write (*,*) denmm_s(m+1, mp+1)
                !stop
                !write (*,*) "denmm_s", m, mp, integrand(2d0), denmm_s(m, mp)
                !write (*,*) G(2d0), d_G(2d0), gamma_ab, phi(m, 2d0), d_phi(m, 2d0), dd_phi(m,2d0)
             end do
          end do          
       end if
    end if
    
    if (make_ortho) then ! make DKE orthogonal w.r.t. to derivative along field line
       denmm_s = matmul(M_transform_inv, denmm_s)
    end if

    call chop(denmm_s)

  contains

    function integrand(x)
      real(kind=dp) :: x, y
      real(kind=dp) :: integrand
      real(kind=dp) :: D_1, D_2

      y = x * gamma_ab

      ! First approach - Use second derivates
      !D_1 = d_G(y) * gamma_ab * x * d_phi_exp(mp, x) &
      !     + G(y) * (-2*x * x * d_phi_exp(mp, x) + x * dd_phi_exp(mp, x) &
      !     + d_phi_exp(mp, x))

      !D_2 = 2 * (1 - T_a/T_b) * (d_G(y) * gamma_ab * x**2 *phi_exp(mp, x) + &
      !     G(y) * (2*x*phi_exp(mp, x) - 2*x**3*phi_exp(mp, x) + x**2 * d_phi_exp(mp, x)))
      !integrand = x**(1+alpha)*exp(-(1+beta)*x**2) * phi_prj(m, x) * (D_1 - D_2)

      ! Second approach - Use only first derivative
      D_1 = x**(1+alpha) * exp(-beta*x**2) * phi_prj(m, x) * (T_a/T_b - 1) * (exp(-x**2) * x &
           * (phi_exp(mp, x) * (-2*(-1+x**2)*G(y) + y*d_G(y)) + x*G(y)*d_phi_exp(mp, x)))
      ! is the same as ( alpha = beta = 0)
      !D_1 = (1-T_a/T_b)*G(y)*(x**2)*exp(-x**2)*phi_exp(mp, x)*(phi_prj(m, x) + x*d_phi_prj(m, x))
      
      D_2 = (exp(-beta*x**2) * x**alpha * ((1 + alpha - 2*beta*x**2) * phi_prj(m, x) + x*d_phi_prj(m, x))) &
           * (x * G(y) * exp(-x**2) * d_phi_exp(mp, x))
      integrand = 2*D_1 - D_2
    end function integrand

    function integrand_rel1(x)
      real(kind=dp) :: x
      real(kind=dp) :: integrand_rel1
      real(kind=dp) :: exp_maxwell, z, gam
      z   = sqrt(2/rmu) * x
      gam = sqrt(1+z**2)
      
      exp_maxwell = 2d0/(sqrt(1+2*x**2/rmu) + 1)
      integrand_rel1 = exp(-exp_maxwell*x**2) * x**2 * D_uu(rmu, x) &
           * d_phi_prj(mp,x) * (phi_prj(m,x) + x*d_phi_prj(m,x))

      !write (*,*) D_uu(rmu,2d0), integrand_rel1
    end function integrand_rel1
  end subroutine compute_energyscattering

  subroutine compute_integralpart(ailmm_s)
    real(kind=dp), dimension(:,:,:) :: ailmm_s
    integer                         :: l, m, mp

    write (*,*) "Computing momentum conservation part..."

    if (precomp) then ! load pre-computed momentum conservation part
       call h5_open(trim(adjustl(matelem_name)), h5id_matelem)
       if (allocated(I1_mmp_s)) deallocate(I1_mmp_s) ! I1_mmp
       allocate(I1_mmp_s(0:lagmax, 0:lagmax, 0:legmax))
       call h5_get(h5id_matelem,'Immp1',I1_mmp_s(:,:,:))
       if (allocated(I2_mmp_s)) deallocate(I2_mmp_s) ! I2_mmp
       allocate(I2_mmp_s(0:lagmax, 0:lagmax, 0:legmax))
       call h5_get(h5id_matelem,'Immp2',I2_mmp_s(:,:,:))
       if (allocated(I3_mmp_s)) deallocate(I3_mmp_s) ! I3_mmp
       allocate(I3_mmp_s(0:lagmax, 0:lagmax, 0:legmax))
       call h5_get(h5id_matelem,'Immp3',I3_mmp_s(:,:,:))
       if (allocated(I4_mmp_s)) deallocate(I4_mmp_s) ! I4_mmp
       allocate(I4_mmp_s(0:lagmax, 0:lagmax, 0:legmax))
       call h5_get(h5id_matelem,'Immp4',I4_mmp_s(:,:,:))
       call h5_close(h5id_matelem)
       !do m=1,size(I1_mmp_s,1)
       !   print *,(I1_mmp_s(m-1,mp-1,5),mp=1,size(I1_mmp_s,2))
       !end do
       !print *,' '
       !do m=1,size(I2_mmp_s,1)
       !   print *,(I2_mmp_s(m-1,mp-1,5),mp=1,size(I2_mmp_s,2))
       !end do
       !print *,' '
       !do m=1,size(I3_mmp_s,1)
       !   print *,(I3_mmp_s(m-1,mp-1,5),mp=1,size(I3_mmp_s,2))
       !end do
       !print *,' '
       !do m=1,size(I4_mmp_s,1)
       !   print *,(I4_mmp_s(m-1,mp-1,5),mp=1,size(I4_mmp_s,2))
       !end do
       !stop
    else
       if (isw_relativistic .eq. 0) then
          !print *,'*********'
          !print *,'I1_mmp:'
          call compute_I1_mmp_s()
          !call disp_gsl_integration_error()
          !print *,'I2_mmp:'
          call compute_I2_mmp_s()
          !call disp_gsl_integration_error()
          !print *,'I3_mmp:'
          call compute_I3_mmp_s()
          !call disp_gsl_integration_error()
          !print *,'I4_mmp:'
          call compute_I4_mmp_s()
          !call disp_gsl_integration_error()
          !print *,'*********'
          
          !if (allocated(ailmm_s)) deallocate(ailmm_s)
          !allocate(ailmm_s(0:lagmax, 0:lagmax, 0:legmax))
          if (make_ortho) then ! make DKE orthogonal w.r.t. to derivative along field line
             do l = 0, legmax
                !ailmm_s(:,:,l+1) = I1_mmp_s(:,:,l) + I2_mmp_s(:,:,l) + I3_mmp_s(:,:,l) + I4_mmp_s(:,:,l)
                ailmm_s(:,:,l+1) = matmul(M_transform_inv, I1_mmp_s(:,:,l) + I2_mmp_s(:,:,l) + I3_mmp_s(:,:,l) + I4_mmp_s(:,:,l))
             end do
          else
             do l = 0, legmax
                ailmm_s(:,:,l+1) = I1_mmp_s(:,:,l) + I2_mmp_s(:,:,l) + I3_mmp_s(:,:,l) + I4_mmp_s(:,:,l)
             end do
          end if

       elseif (isw_relativistic .eq. 1) then
          ailmm_s = 0d0
          call compute_I1_mmp_s()
          call compute_I2_mmp_s()
          call compute_I3_mmp_s()
          call compute_I4_mmp_s()
          do l = 0, legmax
             do m = 0, lagmax
                do mp = 0, lagmax
                   if (l .le. 1) then
                      ailmm_s(m+1,mp+1,l+1) = C_I_rel1_mmp(rmu, l, m, mp)
                   else
                      ailmm_s(:,:,l+1) = I1_mmp_s(:,:,l) + I2_mmp_s(:,:,l) + I3_mmp_s(:,:,l) + I4_mmp_s(:,:,l)
                   end if
                end do
             end do
             if (make_ortho) then
                ailmm_s(:,:,l+1) = matmul(M_transform_inv, ailmm_s(:,:,l+1))
             end if
          end do
       elseif (isw_relativistic .eq. 2) then

          do l = 0, legmax
             do m = 0, lagmax
                do mp = 0, lagmax
                   write (*,'(A, I2, A, I2, A, I2,A)') " Computing matrix element (leg = ", l, ", m = ", m, ", mp = ", mp, ")"
                   ailmm_s(m+1,mp+1,l+1) = C_I_rel2_mmp(rmu, l, m, mp)

                   if (ieee_is_nan(ailmm_s(m+1, mp+1, l+1))) then
                      write (*,*) "*** Matrix element is NaN. Check integration settings. ***"
                      call disp_gsl_integration_error()
                      stop
                   end if
                end do
             end do
             if (make_ortho) then
                ailmm_s(:,:,l+1) = matmul(M_transform_inv, ailmm_s(:,:,l+1))
             end if

          end do
       end if
    end if


    !call chop(ailmm_s)

  end subroutine compute_integralpart
   
  subroutine compute_nbisource(tag_a, tag_b, m_a0, m_b0, T_a0, T_b0, Inbi_lmmp_s)
    use  collisionality_mod, only : lsw_nbi

    character(len=*) :: tag_a, tag_b
    real(kind=dp)    :: m_a0, m_b0
    real(kind=dp)    :: T_a0, T_b0
    real(kind=dp), dimension(:,:) :: Inbi_lmmp_s

    m_a = m_a0
    m_b = m_b0
    T_a = T_a0
    T_b = T_b0

    v_ta = sqrt(2*T_a / m_a)
    v_tb = sqrt(2*T_b / m_b)
    gamma_ab = v_ta/v_tb
    write (*,'(A,A,A,A,A,1E13.6)') " Computing NBI source for ", tag_a, "-", tag_b, " with gamma_ab =", gamma_ab

    ! TODO for NBI (most likely not completed!!)
    ! 1. Read data file
    ! 2. Write linear spline routine for read data including first and second derivative
    ! 3. Change the function pointers of phi_exp, d_phi_exp, and dd_phi_exp to the splined data
    ! 4. Look into compute_integralpart() how the four terms of the integral part are computed and summed up
    ! 5. Do the same in this function, but with modified legmax
    ! 6. Sum up the matrix elements into Inbi_lmmp_s
    ! 7. Set back to the original pointers of phi_exp, d_phi_exp, and dd_phi_exp
    
    
  end subroutine compute_nbisource
  
  subroutine compute_I1_mmp_s()
    integer :: l, m, mp

    if (allocated(I1_mmp_s)) deallocate(I1_mmp_s)
    allocate(I1_mmp_s(0:lagmax, 0:lagmax, 0:legmax))

    do l = 0, legmax
       do m = 0, lagmax
          do mp = 0, lagmax
             I1_mmp_s(m, mp, l) = 3d0/pi**(3d0/2d0) * (2d0*l+1)/2d0 * gamma_ab**3 * m_a/m_b * integrate(int_I1_mmp_s, 0d0)
          end do
       end do
    end do

  contains

    function int_I1_mmp_s(x) result(integrand)
      real(kind=dp) :: x
      real(kind=dp) :: integrand

!!$      integrand = x**(3+alpha) * exp(-(beta + 1)*x**2) * exp(-(gamma_ab*x)**2) * phi_prj(m, x) * phi_exp(mp, x)
      integrand = x**(3+alpha) * exp(-(beta + 1)*x**2) * exp(-(gamma_ab*x)**2) * phi_prj(m, x) * phi_exp(mp, gamma_ab*x)
    end function int_I1_mmp_s

  end subroutine compute_I1_mmp_s

  subroutine compute_I2_mmp_s()
    integer :: l, m, mp

    if (allocated(I2_mmp_s)) deallocate(I2_mmp_s)
    allocate(I2_mmp_s(0:lagmax, 0:lagmax, 0:legmax))

    do l = 0, legmax
       do m = 0, lagmax
          do mp = 0, lagmax
             I2_mmp_s(m, mp, l) = -3d0/pi**(1.5d0) * gamma_ab**3 * integrate(I_phi, 0d0)
             !print *,'l,m,mp: ',l,m,mp
             !call disp_gsl_integration_error()
          end do
       end do
    end do

  contains

    function I_phi(x)
      real(kind=dp) :: x, I1, I2, I_phi

      if (lsw_expand_kernel .and. (x .le. xmax_expand_kernel)) then
         ! use 2nd order Taylor expansion around x=0 (numerical stability)
         I1 = ((l+1d0)*(l+2d0)/(l+3d0)-l)*(x**2/2d0)*phi_exp(mp, 0d0)
         I2 = integrate(I_phi_2, x)

         I_phi = x**(3+alpha) * exp(-(beta+1)*x**2) * phi_prj(m,x) * (I1 + x**l * I2)
      else
         if (lsw_stabilize_Immp) then
          !> Integrating kernel 1 from zero to zero will lead to problems,
          !> as this will result in a 0/0, as then also the parameter is zero.
          !> Thus handle this case explicitly.
          !> Also kernel 2 makes problems for x=0, (and param=0), thus
          !> this is also handled explicitly, and for param=0 the
          !> kernel is zero, except at x=0, which is only a single point
          !> and thus should make no difference. Also this fits the not
          !> stabilized case below, where I2 would be multiplied with 0.
          if (x .eq. 0.0d0) then
            I1 = 0.0d0
            I2 = 0.0d0
          else
            I1 = integrate_param(I_phi_1_param, x, 0d0, x)
            I2 = integrate_param(I_phi_2_param, x, x)
          end if

            I_phi = x**(3+alpha) * exp(-(beta+1)*x**2) * phi_prj(m,x) * (I1 + I2)
         else
            ! use default expression
            I1 = integrate(I_phi_1, 0d0, x)
            I2 = integrate(I_phi_2, x)

            I_phi = x**(3+alpha) * exp(-(beta+1)*x**2) * phi_prj(m,x) * (x**(-l-1) * I1 + x**l * I2)
         end if
      end if

    end function I_phi

    function I_phi_1(xp)
      real(kind=dp) :: xp, yp, I_phi_1

      yp = gamma_ab * xp
!!$      I_phi_1 = exp(-yp**2) * phi_exp(mp, xp) * xp**(l+2)
      I_phi_1 = exp(-yp**2) * phi_exp(mp, yp) * xp**(l+2)
    end function I_phi_1

    function I_phi_2(xp)
      real(kind=dp) :: xp, yp, I_phi_2

      yp = gamma_ab * xp
!!$      I_phi_2 = exp(-yp**2) * phi_exp(mp, xp) * xp**(-l+1)
      I_phi_2 = exp(-yp**2) * phi_exp(mp, yp) * xp**(-l+1)
    end function I_phi_2

    function I_phi_1_param(xp,xval)
      real(kind=dp) :: xp, yp, I_phi_1_param, xval

      yp = gamma_ab * xp
      I_phi_1_param = exp(-yp**2) * phi_exp(mp, yp) * xp * ((xp/xval)**(l+1))
    end function I_phi_1_param

    function I_phi_2_param(xp,xval)
      real(kind=dp) :: xp, yp, I_phi_2_param, xval

      yp = gamma_ab * xp
      I_phi_2_param = exp(-yp**2) * phi_exp(mp, yp) * xp * ((xval/xp)**l)
    end function I_phi_2_param
    
  end subroutine compute_I2_mmp_s

  subroutine compute_I3_mmp_s()
    integer :: l, m, mp

    if (allocated(I3_mmp_s)) deallocate(I3_mmp_s)
    allocate(I3_mmp_s(0:lagmax, 0:lagmax, 0:legmax))

    do l = 0, legmax
       do m = 0, lagmax
          do mp = 0, lagmax
             I3_mmp_s(m, mp, l) = 3d0/pi**(1.5d0) * (1-m_a/m_b) * gamma_ab**3 * integrate(I_phi, 0d0)
          end do
       end do
    end do

  contains

    function I_phi(x)
      real(kind=dp) :: x, I1, I2, I_phi


      if (lsw_expand_kernel .and. (x .le. xmax_expand_kernel)) then
         ! use 2nd order Taylor expansion around x=0 (numerical stability)
         I1 = ((l+1d0)*(l+2d0)/(l+3d0)-l)*(x**2/2d0)*phi_exp(mp, 0d0)
         I2 = integrate(I_phi_2, x)

         I_phi = ((x**(3+alpha)) * exp(-(beta+1)*(x**2)) * &
              (alpha-2*(beta+1)*(x**2)+4) * phi_prj(m,x) + &
              (x**(4+alpha))  * exp(-(beta+1)*(x**2)) * d_phi_prj(m,x)) * &
              (I1 + (x**l) * I2)
      else
         if (lsw_stabilize_Immp) then
            I1 = integrate_param(I_phi_1_param, x, 0d0, x)
            I2 = integrate_param(I_phi_2_param, x, x)

            I_phi = ((x**(3+alpha)) * exp(-(beta+1)*(x**2)) * &
                 (alpha-2*(beta+1)*(x**2)+4) * phi_prj(m,x) + &
                 (x**(4+alpha))  * exp(-(beta+1)*(x**2)) * d_phi_prj(m,x)) * &
                 (I1 + I2)
         else
            ! use default expression
            I1 = integrate(I_phi_1, 0d0, x)
            I2 = integrate(I_phi_2, x)

            I_phi = ((x**(3+alpha)) * exp(-(beta+1)*(x**2)) * &
                 (alpha-2*(beta+1)*(x**2)+4) * phi_prj(m,x) + &
                 (x**(4+alpha))  * exp(-(beta+1)*(x**2)) * d_phi_prj(m,x)) * &
                 ((x**(-l-1)) * I1 + (x**l) * I2)
         end if
      end if

    end function I_phi

    function I_phi_1(xp)
      real(kind=dp) :: xp, yp, I_phi_1

      yp = gamma_ab * xp
!!$      I_phi_1 = exp(-yp**2) * phi_exp(mp, xp) * xp**(l+2)
      I_phi_1 = exp(-yp**2) * phi_exp(mp, yp) * xp**(l+2)
    end function I_phi_1

    function I_phi_2(xp)
      real(kind=dp) :: xp, yp, I_phi_2

      yp = gamma_ab * xp
!!$      I_phi_2 = exp(-yp**2) * phi_exp(mp, xp) * xp**(-l+1)
      I_phi_2 = exp(-yp**2) * phi_exp(mp, yp) * xp**(-l+1)
    end function I_phi_2

    function I_phi_1_param(xp,xval)
      real(kind=dp) :: xp, yp, I_phi_1_param, xval

      yp = gamma_ab * xp
      I_phi_1_param = exp(-yp**2) * phi_exp(mp, yp) * xp * ((xp/xval)**(l+1))
    end function I_phi_1_param

    function I_phi_2_param(xp,xval)
      real(kind=dp) :: xp, yp, I_phi_2_param, xval

      yp = gamma_ab * xp
      I_phi_2_param = exp(-yp**2) * phi_exp(mp, yp) * xp * ((xval/xp)**l)
    end function I_phi_2_param
    
  end subroutine compute_I3_mmp_s

  subroutine compute_I4_mmp_s()
    integer :: l, m, mp
    real(kind=dp) :: ti, tip1, I4_1

    if (allocated(I4_mmp_s)) deallocate(I4_mmp_s)
    allocate(I4_mmp_s(0:lagmax, 0:lagmax, 0:legmax))
       
    if (.not. d_phi_discont) then
       do l = 0, legmax
          do m = 0, lagmax
             do mp = 0, lagmax
                I4_mmp_s(m, mp, l) = 3d0/pi**(1.5d0) * gamma_ab**3 * integrate(I_psi, 0d0)
             end do
          end do
       end do
    else

        do l = 0, legmax
          do m = 0, lagmax
             do mp = 0, lagmax
!!$                if (integral_cutoff) then
!!$                   res_int = fint1d_qag(I_psi, 0d0, x_cutoff, epsabs, epsrel, sw_qag_rule)
!!$                else
!!$                   res_int = fint1d_qagiu(I_psi, 0d0, epsabs, epsrel)
!!$                end if
!!$                I4_mmp_s(m, mp, l) = 3d0/pi**(1.5d0) * gamma_ab**3 * res_int(1)

                I4_mmp_s(m, mp, l) = 0d0
                do i = lbound(t_vec, 1), ubound(t_vec, 1)-1
                   ti = t_vec(i)
                   tip1 = t_vec(i+1)
                   
                   I4_1 = I4_part(tip1 - 1d-10) - I4_part(ti + 1d-10)

                   !write (*,*) ti, tip1, I4_part(ti), I4_part(tip1)
                   
                   !res_int = fint1d_qag(I_psi, ti, tip1, epsabs, epsrel, sw_qag_rule)
                   I4_mmp_s(m, mp, l) = I4_mmp_s(m, mp, l) - 3d0/pi**(1.5d0) * gamma_ab**3 * (I4_1 - integrate(I_psi, ti, tip1))

                   !write (*,*) d_phi_prj(1, 1.999999d0),   d_phi_prj(1, 2d0), d_phi_prj(1, 2.0001d0)
                   !stop
                   
                end do
                
             end do
          end do
       end do      
       
    end if
    
  contains

    function I4_part(x)
      real(kind=dp) :: I4_part, x

      I4_part = ((x**(4 + alpha)*((5 + alpha - 2*(1 + beta)*x**2) * phi_prj(m,x) &
           + x*d_phi_prj(m, x))) *  exp(-((1 + beta)*x**2))) * K(x)

    end function I4_part
    
    function K(x)
      real(kind=dp) :: K, x
      real(kind=dp) :: I1, I2, I3, I4

      if (lsw_expand_kernel .and. (x .le. xmax_expand_kernel)) then
         ! use 2nd order Taylor expansion around x=0 (numerical stability)
         
         !res_int = fint1d_qag(I_psi_1, 0d0, x, epsabs, epsrel, sw_qag_rule)
         I1 = ( 4d0 - 2d0*l + (l+1d0)*(l+2d0)*(-2d0-l+(l+3d0)*(l+4d0)/(l+4d0)) ) * &
              (x**4/24d0) * phi_exp(mp, 0d0) / (2d0*l+3d0)

         !res_int = fint1d_qag(I_psi_2, 0d0, x, epsabs, epsrel, sw_qag_rule)
         I2 = ( 8d0 - 2d0*l + (l-1d0)*l*(-l+(l+1d0)*(l+2d0)/(l+3d0)) ) * &
              (x**4/24d0) * phi_exp(mp, 0d0) / (2d0*l-1d0)

         !if (integral_cutoff) then
         !   res_int = fint1d_qag(I_psi_3, x, x_cutoff, epsabs, epsrel, sw_qag_rule)
         !else
         !   res_int = fint1d_qagiu(I_psi_3, x, epsabs, epsrel)
         !end if
         I3 = integrate(I_psi_3, x)

         !if (integral_cutoff) then
         !   res_int = fint1d_qag(I_psi_4, x, x_cutoff, epsabs, epsrel, sw_qag_rule)
         !else
         !   res_int = fint1d_qagiu(I_psi_4, x, epsabs, epsrel)
         !end if
         I4 = integrate(I_psi_4, x)

         K = x**(-l-1)/(2*l+3)*I1 -  x**(-l+1)/(2*l-1)*I2 +  x**(l+2)/(2*l+3)*I3 - x**l/(2*l-1)*I4
      else
         if (lsw_stabilize_Immp) then
            I1 = integrate_param(I_psi_1_param, x, 0d0, x)
            I2 = integrate_param(I_psi_2_param, x, 0d0, x)
            I3 = integrate_param(I_psi_3_param, x, x)
            I4 = integrate_param(I_psi_4_param, x, x)

            K = I1/(2*l+3) -  I2/(2*l-1) +  I3/(2*l+3) - I4/(2*l-1)
         else
            ! use default expression

            !res_int = fint1d_qag(I_psi_1, 0d0, x, epsabs, epsrel, sw_qag_rule)
            I1 = integrate(I_psi_1, 0d0, x)

            !res_int = fint1d_qag(I_psi_2, 0d0, x, epsabs, epsrel, sw_qag_rule)
            I2 = integrate(I_psi_2, 0d0, x)

            !if (integral_cutoff) then
            !   res_int = fint1d_qag(I_psi_3, x, x_cutoff, epsabs, epsrel, sw_qag_rule)
            !else
            !   res_int = fint1d_qagiu(I_psi_3, x, epsabs, epsrel)
            !end if
            I3 = integrate(I_psi_3, x)

            !if (integral_cutoff) then
            !   res_int = fint1d_qag(I_psi_4, x, x_cutoff, epsabs, epsrel, sw_qag_rule)
            !else
            !   res_int = fint1d_qagiu(I_psi_4, x, epsabs, epsrel)
            !end if
            I4 = integrate(I_psi_4, x)

            K = x**(-l-1)/(2*l+3)*I1 -  x**(-l+1)/(2*l-1)*I2 +  x**(l+2)/(2*l+3)*I3 - x**l/(2*l-1)*I4
         end if
      end if

    end function K
    
    function I_psi(x)
      real(kind=dp) :: x, I_psi

      ! First approach - use second derivates
      !  c = 1 + beta
!!$      I_psi = exp(-c*x**2)*x**(3+alpha) * ((20+9*alpha + &
!!$           alpha**2 - 2*(11+2*alpha) * c*x**2 + 4*c**2*x**4) * phi_prj(m,x) + &
!!$           x*(2*(5+alpha-2*c*x**2) * d_phi_prj(m, x) + x*dd_phi_prj(m,x))) * &
!!$           (x**(-l-1)/(2*l+3)*I1 - x**(-l+1)/(2*l-1)*I2 + x**(l+2)/(2*l+3)*I3 - x**l/(2*l-1)*I4)

      ! Second approach - avoid second derivates
      !K_temp = K(x)
      !I4_1 = ((x**(4 + alpha)*((5 + alpha - 2*(1 + beta)*x**2) * phi_prj(m,x) &
      !       + x*d_phi_prj(m, x))) *  exp(-((1 + beta)*x**2))) * K_temp

      I_psi = ((x**(3 + alpha) * ((20 + 9*alpha + alpha**2 - 2*(11 + 2*alpha)*(1 + beta)*x**2 + 4*(1 + beta)**2*x**4)*phi_prj(m,x) &
           + x*(2*(5 + alpha - 2*(1 + beta)*x**2)*d_phi_prj(m,x) + x*dd_phi_prj(m,x)))) * exp(-(1 + beta)*x**2)) * K(x)
      
    end function I_psi

    function I_psi_1(xp)
      real(kind=dp) :: xp, yp, I_psi_1

      yp = gamma_ab * xp
!!$      I_psi_1 = xp**(l+4)*exp(-yp**2)*phi_exp(mp,xp)
      I_psi_1 = xp**(l+4)*exp(-yp**2)*phi_exp(mp,yp)
    end function I_psi_1

    function I_psi_2(xp)
      real(kind=dp) :: xp, yp, I_psi_2

      yp = gamma_ab * xp
!!$      I_psi_2 = xp**(l+2)*exp(-yp**2)*phi_exp(mp,xp)
      I_psi_2 = xp**(l+2)*exp(-yp**2)*phi_exp(mp,yp)
    end function I_psi_2

    function I_psi_3(xp)
      real(kind=dp) :: xp, yp, I_psi_3

      yp = gamma_ab * xp
!!$      I_psi_3 = xp**(-l+1)*exp(-yp**2)*phi_exp(mp,xp)
      I_psi_3 = xp**(-l+1)*exp(-yp**2)*phi_exp(mp,yp)
    end function I_psi_3

    function I_psi_4(xp)
      real(kind=dp) :: xp, yp, I_psi_4

      yp = gamma_ab * xp
!!$      I_psi_4 = xp**(-l+3)*exp(-yp**2)*phi_exp(mp,xp)
      I_psi_4 = xp**(-l+3)*exp(-yp**2)*phi_exp(mp,yp)
    end function I_psi_4

    function I_psi_1_param(xp,xval)
      real(kind=dp) :: xp, yp, I_psi_1_param, xval

      yp = gamma_ab * xp
      I_psi_1_param = (xp**3) * (xp/xval)**(l+1) * exp(-yp**2) * phi_exp(mp,yp)
    end function I_psi_1_param

    function I_psi_2_param(xp,xval)
      real(kind=dp) :: xp, yp, I_psi_2_param, xval

      yp = gamma_ab * xp
      I_psi_2_param = (xp**3) * ((xp/xval)**(l-1)) * exp(-yp**2) * phi_exp(mp,yp)
    end function I_psi_2_param

    function I_psi_3_param(xp,xval)
      real(kind=dp) :: xp, yp, I_psi_3_param, xval

      yp = gamma_ab * xp
      I_psi_3_param = (xp**3) * ((xval/xp)**(l+2)) * exp(-yp**2) * phi_exp(mp,yp)
    end function I_psi_3_param

    function I_psi_4_param(xp,xval)
      real(kind=dp) :: xp, yp, I_psi_4_param, xval

      yp = gamma_ab * xp
      I_psi_4_param = (xp**3) * ((xval/xp)**l) * exp(-yp**2) * phi_exp(mp,yp)
    end function I_psi_4_param
    
  end subroutine compute_I4_mmp_s

  subroutine compute_source(asource_s, weightlag_s, bzero_s, weightparflow_s, weightenerg_s, Amm_s)
    real(kind=dp), dimension(:,:), intent(out) :: Amm_s, asource_s, weightlag_s
    real(kind=dp), dimension(:), intent(out)   :: bzero_s, weightparflow_s, weightenerg_s

    if (allocated(M_transform)) deallocate(M_transform)
    allocate(M_transform(0:lagmax, 0:lagmax))
    
    if (allocated(M_transform_inv)) deallocate(M_transform_inv)
    allocate(M_transform_inv(0:lagmax, 0:lagmax))

    if (allocated(weightenerg_offset)) deallocate(weightenerg_offset)
    allocate(weightenerg_offset(0:lagmax))
    
    call compute_Minv(M_transform_inv)
    call disp_gsl_integration_error()
    call compute_sources(asource_s, weightlag_s, bzero_s, weightparflow_s, weightenerg_s)
    call disp_gsl_integration_error()
    Amm_s=M_transform

  end subroutine compute_source

  subroutine compute_collop_rel(isw_rel, T_e, asource_s, weightlag_s, bzero_s, weightparflow_s, &
       weightenerg_s, Amm_s, anumm_s, anumm_inf_s, denmm_s, ailmm_s)
    real(kind=dp), dimension(:,:) :: asource_s, weightlag_s, Amm_s
    real(kind=dp), dimension(:)   :: bzero_s, weightparflow_s, weightenerg_s
    real(kind=dp), dimension(:,:)   :: anumm_s, anumm_inf_s, denmm_s
    real(kind=dp), dimension(:,:,:) :: ailmm_s
    real(kind=dp) :: T_e!, rmu_beg, rmu_end
    integer :: ierr, n, isw_rel!, rmu_steps

    !**********************************************************
    ! Preparations for relativistic formulas
    !**********************************************************
    isw_relativistic = isw_rel
    rmu = (c**2 * m_ele)/(eV*T_e)
    
    n=2
    call DBESK(rmu,n,norm_maxwell,ierr)
    norm_maxwell = sqrt(pi/2d0) * 1/norm_maxwell
    
    write (*,*) "Inverse relativistic velocity mu = ", rmu
    write (*,*) "Relativistic Maxwell normalization norm_maxwell = ", norm_maxwell

    !**********************************************************
    ! Define constants for offset correction
    !**********************************************************
    a_00_offset = 1.00d0
    a_02_offset = (-1d0) * dmuk2ovk2(rmu)
    a_22_offset = ddmuk2ovk2(rmu)
    
!!$    rmu_beg = 0d0
!!$    rmu_end = 1d4
!!$    rmu_steps = 100000
!!$    do i = 1, rmu_steps
!!$
!!$       rmu = (rmu_end - rmu_beg)/rmu_steps * i
!!$       write (110,*) rmu, dk2ovk2(rmu), ddk2ovk2(rmu)
!!$       
!!$    end do
!!$       
!!$    stop
    
    if (allocated(M_transform)) deallocate(M_transform)
    allocate(M_transform(0:lagmax, 0:lagmax))
    
    if (allocated(M_transform_inv)) deallocate(M_transform_inv)
    allocate(M_transform_inv(0:lagmax, 0:lagmax))

    if (allocated(weightenerg_offset)) deallocate(weightenerg_offset)
    allocate(weightenerg_offset(0:lagmax))

    call compute_Minv(M_transform_inv)
    call disp_gsl_integration_error()
    call compute_sources(asource_s, weightlag_s, bzero_s, weightparflow_s, weightenerg_s)
    call disp_gsl_integration_error()
    Amm_s=M_transform
    
    gamma_ab = 1.0d0
    m_a = 1d0
    m_b = 1d0
    T_a = 1d0
    T_b = 1d0
    write (*,'(A,A,A,A,A,ES13.6)') " Computing relativistic collision operator"
    
    call compute_lorentz(anumm_s)
    call disp_gsl_integration_error()
    call compute_lorentz_inf(anumm_inf_s)
    call disp_gsl_integration_error()
    call compute_energyscattering(denmm_s)
    call disp_gsl_integration_error()
    call compute_integralpart(ailmm_s)
    call disp_gsl_integration_error()

    !write (*,*) 2d0/3d0 * ailmm_s(:,2,2) / (anumm_s(:,2) - denmm_s(:,2))
 
  end subroutine compute_collop_rel
  
  subroutine compute_collop(tag_a, tag_b, m_a0, m_b0, T_a0, T_b0, anumm_s, denmm_s, ailmm_s)
    character(len=*) :: tag_a, tag_b
    real(kind=dp)    :: m_a0, m_b0
    real(kind=dp)    :: T_a0, T_b0
    real(kind=dp), dimension(:,:)   :: anumm_s, denmm_s
    real(kind=dp), dimension(:,:,:) :: ailmm_s

    m_a = m_a0
    m_b = m_b0
    T_a = T_a0
    T_b = T_b0

    v_ta = sqrt(2*T_a / m_a)
    v_tb = sqrt(2*T_b / m_b)
    gamma_ab = v_ta/v_tb
    write (*,'(A,A,A,A,A,1E13.6)') " Computing collision operator for ", tag_a, "-", tag_b, " with gamma_ab =", gamma_ab

    !**********************************************************
    ! Define constants for offset correction
    !**********************************************************
    a_00_offset = 1.00d0
    a_02_offset = 1.50d0
    a_22_offset = 3.75d0
    
    call compute_lorentz(anumm_s)
    call disp_gsl_integration_error()
    call compute_energyscattering(denmm_s)
    call disp_gsl_integration_error()    
    call compute_integralpart(ailmm_s)
    call disp_gsl_integration_error()
    
  end subroutine compute_collop
  
  subroutine compute_collop_inf(tag_a, tag_b, m_a0, m_b0, T_a0, T_b0, anumm_s, anumm_inf_s, denmm_s, ailmm_s)
    character(len=*) :: tag_a, tag_b
    real(kind=dp)    :: m_a0, m_b0
    real(kind=dp)    :: T_a0, T_b0
    real(kind=dp), dimension(:,:)   :: anumm_s, anumm_inf_s, denmm_s
    real(kind=dp), dimension(:,:,:) :: ailmm_s

    m_a = m_a0
    m_b = m_b0
    T_a = T_a0
    T_b = T_b0

    v_ta = sqrt(2*T_a / m_a)
    v_tb = sqrt(2*T_b / m_b)
    gamma_ab = v_ta/v_tb
    write (*,'(A,A,A,A,A,ES13.6)') " Computing collision operator for ", tag_a, "-", tag_b, " with gamma_ab =", gamma_ab

    !**********************************************************
    ! Define constants for offset correction
    !**********************************************************
    a_00_offset = 1.00d0
    a_02_offset = 1.50d0
    a_22_offset = 3.75d0
    
    call compute_lorentz(anumm_s)
    call disp_gsl_integration_error()
    call compute_lorentz_inf(anumm_inf_s)
    call disp_gsl_integration_error()
    call compute_energyscattering(denmm_s)
    call disp_gsl_integration_error()
    call compute_integralpart(ailmm_s)
    call disp_gsl_integration_error()

  end subroutine compute_collop_inf

  subroutine compute_collop_lorentz(tag_a, tag_b, m_a0, m_b0, T_a0, T_b0, anumm_s)
    character(len=*) :: tag_a, tag_b
    real(kind=dp)    :: m_a0, m_b0
    real(kind=dp)    :: T_a0, T_b0
    real(kind=dp), dimension(:,:)   :: anumm_s

    m_a = m_a0
    m_b = m_b0
    T_a = T_a0
    T_b = T_b0

    v_ta = sqrt(2*T_a / m_a)
    v_tb = sqrt(2*T_b / m_b)
    gamma_ab = v_ta/v_tb
    write (*,'(A,A,A,A,A,ES13.6)') " Computing collision operator for ", tag_a, "-", tag_b, " with gamma_ab =", gamma_ab

    !**********************************************************
    ! Define constants for offset correction
    !**********************************************************
    a_00_offset = 1.00d0
    a_02_offset = 1.50d0
    a_22_offset = 3.75d0
    
    call compute_lorentz(anumm_s)
    call disp_gsl_integration_error()

  end subroutine compute_collop_lorentz

end module collop_compute
