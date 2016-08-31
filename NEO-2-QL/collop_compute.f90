module collop_compute

  use hdf5_tools
  !use nrtype, only : dp, pi
  use gsl_integration_routines_mod
  use collop_laguerre
  use collop_polynomial
  
  implicit none
  
  !**********************************************************
  ! Species definitions (cgs)
  !**********************************************************
  real(kind=dp) :: m_ele = 9.109382150d-28
  real(kind=dp) :: m_pro = 1.672621637d-24
  real(kind=dp) :: m_alp = 6.644656200d-24
  real(kind=dp) :: m_d   = 3.343583719d-24
  real(kind=dp) :: m_C   = 19.94406876d-24
  
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
  !real(kind=dp) :: x_max    = 20

  !**********************************************************
  ! Pre-computed matrix elements
  !**********************************************************
  !character(len=100) :: matelem_name='MatElem_aa_hatfun_lag10_xmax6.h5'
  !character(len=100) :: matelem_name='MatElem_aa_hatfun_lag5_xmax6.h5'
  !character(len=100) :: matelem_name='MatElem_aa_hatfun_lag30_xmax4.h5'
  !character(len=100) :: matelem_name='MatElem_aa_hatfun_lag20_xmax4.h5'
  !character(len=100) :: matelem_name='MatElem_aa_hatfun_lag15_xmax4.h5'
  character(len=100) :: matelem_name='MatElem_aa_hatfun_lag10_xmax4.h5'
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
  
contains

  subroutine init_collop(collop_base_prj, collop_base_exp, scalprod_alpha, scalprod_beta)
    use rkstep_mod, only : lag, leg
    integer :: collop_base_prj, collop_base_exp
    real(kind=dp) :: scalprod_alpha
    real(kind=dp) :: scalprod_beta

    alpha = scalprod_alpha
    beta  = scalprod_beta
    lagmax = lag
    legmax = leg
    
    if (collop_base_prj .eq. 0) then
       write (*,*) "Using Laguerre polynomials as collision operator projection base."
       init_phi_prj => init_phi_laguerre
       phi_prj      => phi_laguerre
       d_phi_prj    => d_phi_laguerre
       dd_phi_prj   => dd_phi_laguerre
    elseif (collop_base_prj .eq. 1) then
        write (*,*) "Using standard polynomials as collision operator projection base."
       init_phi_prj => init_phi_polynomial
       phi_prj      => phi_polynomial
       d_phi_prj    => d_phi_polynomial
       dd_phi_prj   => dd_phi_polynomial     
    elseif (collop_base_prj .eq. 2) then
       write (*,*) "Using squared polynomials as collision operator projection base."
       init_phi_prj => init_phi_polynomial_2
       phi_prj      => phi_polynomial_2
       d_phi_prj    => d_phi_polynomial_2
       dd_phi_prj   => dd_phi_polynomial_2
    elseif (collop_base_prj .eq. 3) then
       write (*,*) "Using squared polynomials without zeroth order as collision operator projection base."
       init_phi_prj => init_phi_polynomial_3
       phi_prj      => phi_polynomial_3
       d_phi_prj    => d_phi_polynomial_3
       dd_phi_prj   => dd_phi_polynomial_3
    elseif (collop_base_prj .eq. 10) then
       write (*,*) "Using hat functions as collision operator projection base."
       precomp=.true.
       make_ortho=.false.
    else
       write (*,*) "Undefined collision operator projection base ", collop_base_prj
       stop
    end if

    if (collop_base_exp .eq. 0) then
       write (*,*) "Using Laguerre polynomials as collision operator expansion base."
       init_phi_exp => init_phi_laguerre
       phi_exp      => phi_laguerre
       d_phi_exp    => d_phi_laguerre
       dd_phi_exp   => dd_phi_laguerre
    elseif (collop_base_exp .eq. 1) then
        write (*,*) "Using standard polynomials as collision operator expansion base."
       init_phi_exp => init_phi_polynomial
       phi_exp      => phi_polynomial
       d_phi_exp    => d_phi_polynomial
       dd_phi_exp   => dd_phi_polynomial     
    elseif (collop_base_exp .eq. 2) then
       write (*,*) "Using squared polynomials as collision operator expansion base."
       init_phi_exp => init_phi_polynomial_2
       phi_exp      => phi_polynomial_2
       d_phi_exp    => d_phi_polynomial_2
       dd_phi_exp   => dd_phi_polynomial_2
    elseif (collop_base_exp .eq. 3) then
       write (*,*) "Using squared polynomials without zeroth order as collision operator expansion base."
       init_phi_exp => init_phi_polynomial_3
       phi_exp      => phi_polynomial_3
       d_phi_exp    => d_phi_polynomial_3
       dd_phi_exp   => dd_phi_polynomial_3
    elseif (collop_base_exp .eq. 10) then
       write (*,*) "Using hat functions as collision operator expansion base."
    else
       write (*,*) "Undefined collision operator expansion base ", collop_base_exp
       stop
    end if
    
    call init_legendre(legmax)
    call init_phi_prj(lagmax, legmax)
    call init_phi_exp(lagmax, legmax)
    
  end subroutine init_collop

  subroutine chop_0(x)
    real(kind=dp) :: x, chop

    if (abs(x) .lt. epsabs) then
       x = 0d0
    elseif (abs(x - 1d0) .lt. epsabs) then
       x = 1d0
    end if

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

  function nu_D_hat(x) result(nu)
    real(kind=dp) :: x, y, nu

    ! At the momement only for self-collisions
    y = x * 1d0

    nu = 3*sqrt(pi)/4d0 * (erf(y) - G(y))/(x**3)
  end function nu_D_hat

  subroutine inv(A)
    real(dp), dimension(:,:) :: A
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
  
  subroutine compute_Minv(Minv)
    real(kind=dp), dimension(:,:) :: Minv
    real(kind=dp), dimension(2)   :: res_int
    integer :: mm, mp

    write (*,*) "Computing phi transformation matrix..."

    if (precomp) then ! load pre-computed M_transform
       call h5_open(trim(adjustl(matelem_name)), h5id_matelem)
       call h5_get(h5id_matelem,'Ammp',Minv)
       call h5_close(h5id_matelem)
       !DO mm=1,SIZE(Minv,1)
       !   PRINT *,(Minv(mm,mp),mp=1,SIZE(Minv,2))
       !END DO
       !STOP
    else
       do mm = 0, lagmax
          do mp = 0, lagmax
             res_int = fint1d_qagiu(phim_phimp, 0d0, epsabs, epsrel)
             Minv(mm+1,mp+1) = res_int(1)
          end do
       end do
    end if

    if (make_ortho) then ! make DKE orthogonal w.r.t. to derivative along field line
       do mm = 0, lagmax
          do mp = 0, lagmax
             if (mm .eq. mp) then
                M_transform(mm,mp)=1.0d0
             else
                M_transform(mm,mp)=0.0d0
             end if
          end do
       end do
    else
       M_transform = Minv
    end if

    call inv(Minv)
    call chop(Minv)

  contains

    function phim_phimp(x)
      real(kind=dp) :: x, phim_phimp

      phim_phimp =  pi**(-3d0/2d0) * x**(4+alpha) * exp(-(1+beta)*x**2) * phi_prj(mm,x) * phi_exp(mp,x)
      
    end function phim_phimp
  end subroutine compute_Minv

  subroutine compute_sources(asource_s, weightlag_s)
    real(kind=dp), dimension(:,:) :: asource_s, weightlag_s
    real(kind=dp), dimension(2) :: res_int
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
       do k = 1, 3
          do m = 0, lagmax
             res_int = fint1d_qagiu(am, 0d0, epsabs, epsrel)
             asource_s(m+1, k) = res_int(1)
          end do
       end do
    end if 

    if (make_ortho) then ! make DKE orthogonal w.r.t. to derivative along field line
       do k = 1, 3
          asource_s(:,k) = matmul(M_transform_inv, asource_s(:,k))
       end do
    end if

    call chop(asource_s)
    !write (*,*) "Done."
    
    write (*,*) "Computing weighting coefficients..."

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
       do j = 1, 3
          do m = 0, lagmax
             !write (*,*) j, m, lbound(weightlag_s), ubound(weightlag_s)
             res_int = fint1d_qagiu(bm, 0d0, epsabs, epsrel)
             weightlag_s(j,m+1) = 1d0/sqrt(pi) * res_int(1)
          end do
       end do
    end if

    ! weightlag for computation of bvec_parflow
    if (make_ortho) then ! make DKE orthogonal w.r.t. to derivative along field line
       weightlag_s(4,:) = asource_s(:,1)
    else
       weightlag_s(4,:) = matmul(M_transform_inv, asource_s(:,1))
    end if

    call chop(weightlag_s)
    !write (*,*) "Done."

    contains

      function am(x)
        real(kind=dp) :: x, am

        am = pi**(-3d0/2d0) * x**(4+alpha) * exp(-(1+beta)*x**2) * phi_prj(m, x) * x**(2*k - 1 - 5*kdelta(3,k))
      end function am

      function bm(x)
        real(kind=dp) :: x, bm
        bm = exp(-x**2) * x**(2*(j+1)-5*kdelta(3,j)) * phi_prj(m,x)
      end function bm

      function kdelta(a,b)
        integer :: a, b
        integer :: kdelta

        kdelta = 0
        if (a .eq. b) kdelta = 1
      end function kdelta

    end subroutine compute_sources

    subroutine compute_xmmp(x1mm_s,x2mm_s)
      real(kind=dp), dimension(:,:) :: x1mm_s, x2mm_s
      real(kind=dp), dimension(2) :: res_int
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
               res_int = fint1d_qagiu(integrand_x1mmp, 0d0, epsabs, epsrel)
               x1mm_s(m+1, mp+1) = cnorm * res_int(1)
               res_int = fint1d_qagiu(integrand_x2mmp, 0d0, epsabs, epsrel)
               x2mm_s(m+1, mp+1) = cnorm * res_int(1)
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
      real(kind=dp), dimension(2) :: res_int
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
         do m = 0, lagmax
            do mp = 0, lagmax
               res_int = fint1d_qagiu(integrand, 0d0, epsabs, epsrel)
               anumm_s(m+1, mp+1) = 3d0/(4d0 * pi) * res_int(1)
            end do
         end do
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
    
  end subroutine compute_lorentz

  subroutine compute_lorentz_inf(anumm_s)
    real(kind=dp), dimension(:,:) :: anumm_s
    real(kind=dp), dimension(2) :: res_int
    integer :: m, mp

    !if (allocated(anumm_s)) deallocate(anumm_s)
    !allocate(anumm_s(0:lagmax, 0:lagmax))

    write (*,*) "Computing Lorentz part (gamma -> inf)..."

    do m = 0, lagmax
       do mp = 0, lagmax
          res_int = fint1d_qagiu(integrand_inf, 0d0, epsabs, epsrel)
          anumm_s(m+1, mp+1) = 3d0/(4d0 * pi) * res_int(1)
       end do
    end do

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
      integrand_inf = x**(alpha) * exp(-(1+beta)*x**2) * phi_prj(m, x) * (1 - 0) * phi_exp(mp, x)

    end function integrand_inf

  end subroutine compute_lorentz_inf
   
  subroutine compute_energyscattering(denmm_s)
    real(kind=dp), dimension(:,:) :: denmm_s
    integer :: m, mp
    real(kind=dp), dimension(2) :: res_int

    !if (allocated(denmm_s)) deallocate(denmm_s)
    !allocate(denmm_s(0:lagmax, 0:lagmax))

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
       do m = 0, lagmax
          do mp = 0, lagmax
             res_int = fint1d_qagiu(integrand, 0d0, epsabs, epsrel)
             !res_int = fint1d_qag(integrand, 0d0, 100d0, epsabs, epsrel, 2)       
             denmm_s(m+1, mp+1) = 3d0/(4d0 * pi) * res_int(1)

             !write (*,*) "denmm_s", m, mp, integrand(2d0), denmm_s(m, mp)
             !write (*,*) G(2d0), d_G(2d0), gamma_ab, phi(m, 2d0), d_phi(m, 2d0), dd_phi(m,2d0)
          end do
       end do
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

      D_1 = d_G(y) * gamma_ab * x * d_phi_exp(mp, x) &
           + G(y) * (-2*x * x * d_phi_exp(mp, x) + x * dd_phi_exp(mp, x) &
           + d_phi_exp(mp, x))

      D_2 = 2 * (1 - T_a/T_b) * (d_G(y) * gamma_ab * x**2 *phi_exp(mp, x) + &
           G(y) * (2*x*phi_exp(mp, x) - 2*x**3*phi_exp(mp, x) + x**2 * d_phi_exp(mp, x)))
      
      integrand = x**(1+alpha)*exp(-(1+beta)*x**2) * phi_prj(m, x) * (D_1 - D_2)
    end function integrand
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
       call compute_I1_mmp_s()
       call compute_I2_mmp_s()
       call compute_I3_mmp_s()
       call compute_I4_mmp_s()
    end if

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
    
    
    call chop(ailmm_s)
    !write (*,*) "Done."
    
  end subroutine compute_integralpart

  subroutine compute_I1_mmp_s()
    integer :: l, m, mp
    real(kind=dp), dimension(2) :: res_int

    if (allocated(I1_mmp_s)) deallocate(I1_mmp_s)
    allocate(I1_mmp_s(0:lagmax, 0:lagmax, 0:legmax))
    
    do l = 0, legmax
       do m = 0, lagmax
          do mp = 0, lagmax  
             !res_int = fint1d_qag(int_I1_mmp_s, 0d0, 10d0, epsabs, epsrel, sw_qag_rule)
             res_int = fint1d_qagiu(int_I1_mmp_s, 0d0, epsabs, epsrel)

             I1_mmp_s(m, mp, l) = 3d0/pi**(3d0/2d0) * (2d0*l+1)/2d0 * gamma_ab**3 * m_a/m_b * res_int(1)
             
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
    real(kind=dp), dimension(2) :: res_int

    if (allocated(I2_mmp_s)) deallocate(I2_mmp_s)
    allocate(I2_mmp_s(0:lagmax, 0:lagmax, 0:legmax))
    
    do l = 0, legmax
       do m = 0, lagmax
          do mp = 0, lagmax  
             !res_int = fint1d_qag(I_phi, 0d0, 10d0, epsabs, epsrel, sw_qag_rule)
             res_int = fint1d_qagiu(I_phi, 0d0, epsabs, epsrel)

             I2_mmp_s(m, mp, l) = -3d0/pi**(1.5d0) * gamma_ab**3 * res_int(1)
          end do
       end do
    end do
    
  contains

    function I_phi(x)
      real(kind=dp) :: x, I1, I2, I_phi
      real(kind=dp), dimension(2) :: res_int
    
      res_int = fint1d_qag(I_phi_1, 0d0, x, epsabs, epsrel, sw_qag_rule)
      I1 = res_int(1)

      !res_int = fint1d_qag(I_phi_2, x, 10d0, epsabs, epsrel, sw_qag_rule)
      res_int = fint1d_qagiu(I_phi_2, x, epsabs, epsrel)

      I2 = res_int(1)
      
      I_phi = x**(3+alpha) * exp(-(beta+1)*x**2) * phi_prj(m,x) * (x**(-l-1) * I1 + x**l * I2)
      
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
    
  end subroutine compute_I2_mmp_s

  subroutine compute_I3_mmp_s
    integer :: l, m, mp
    real(kind=dp), dimension(2) :: res_int

    if (allocated(I3_mmp_s)) deallocate(I3_mmp_s)
    allocate(I3_mmp_s(0:lagmax, 0:lagmax, 0:legmax))

    do l = 0, legmax
       do m = 0, lagmax
          do mp = 0, lagmax  
             !res_int = fint1d_qag(I_phi, 0d0, 10d0, epsabs, epsrel, sw_qag_rule)
             res_int = fint1d_qagiu(I_phi, 0d0, epsabs, epsrel)

             I3_mmp_s(m, mp, l) = 3d0/pi**(1.5d0) * (1-m_a/m_b) * gamma_ab**3 * res_int(1)
          end do
       end do
    end do

  contains

    function I_phi(x)
      real(kind=dp) :: x, I1, I2, I_phi
      real(kind=dp), dimension(2) :: res_int
    
      res_int = fint1d_qag(I_phi_1, 0d0, x, epsabs, epsrel, sw_qag_rule)
      I1 = res_int(1)

      !res_int = fint1d_qag(I_phi_2, x, 10d0, epsabs, epsrel, sw_qag_rule)
      res_int = fint1d_qagiu(I_phi_2, x, epsabs, epsrel)

      I2 = res_int(1)
      
      I_phi = ((x**(3+alpha)) * exp(-(beta+1)*(x**2)) * &
           (alpha-2*(beta+1)*(x**2)+4) * phi_prj(m,x) + &
           (x**(4+alpha))  * exp(-(beta+1)*(x**2)) * d_phi_prj(m,x)) * &
           ((x**(-l-1)) * I1 + (x**l) * I2)      
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
    
  end subroutine compute_I3_mmp_s

  subroutine compute_I4_mmp_s()
    integer :: l, m, mp
    real(kind=dp), dimension(2) :: res_int

    if (allocated(I4_mmp_s)) deallocate(I4_mmp_s)
    allocate(I4_mmp_s(0:lagmax, 0:lagmax, 0:legmax))
    
    do l = 0, legmax
       do m = 0, lagmax
          do mp = 0, lagmax  
             !res_int = fint1d_qag(I_psi, 0d0, 10d0, epsabs, epsrel, sw_qag_rule)
             res_int = fint1d_qagiu(I_psi, 0d0, epsabs, epsrel)
             I4_mmp_s(m, mp, l) = 3d0/pi**(1.5d0) * gamma_ab**3 * res_int(1)
          end do
       end do
    end do

  contains

    function I_psi(x)
      real(kind=dp) :: x, I1, I2, I3, I4, I_psi
      real(kind=dp), dimension(2) :: res_int
      real(kind=dp) :: c

      c = 1 + beta
      
      res_int = fint1d_qag(I_psi_1, 0d0, x, epsabs, epsrel, sw_qag_rule)
      I1 = res_int(1)

      res_int = fint1d_qag(I_psi_2, 0d0, x, epsabs, epsrel, sw_qag_rule)
      I2 = res_int(1)

      !res_int = fint1d_qag(I_psi_3, x, 20d0, epsabs, epsrel, sw_qag_rule)
      res_int = fint1d_qagiu(I_psi_3, x, epsabs, epsrel)
      I3 = res_int(1)
      
      !res_int = fint1d_qag(I_psi_4, x, 20d0, epsabs, epsrel, sw_qag_rule)
      res_int = fint1d_qagiu(I_psi_4, x, epsabs, epsrel)
      I4 = res_int(1)
      
      I_psi = exp(-c*x**2)*x**(3+alpha) * ((20+9*alpha + &
           alpha**2 - 2*(11+2*alpha) * c*x**2 + 4*c**2*x**4) * phi_prj(m,x) + &
           x*(2*(5+alpha-2*c*x**2) * d_phi_prj(m, x) + x*dd_phi_prj(m,x))) * &
           (x**(-l-1)/(2*l+3)*I1 - x**(-l+1)/(2*l-1)*I2 + x**(l+2)/(2*l+3)*I3 - x**l/(2*l-1)*I4)
      
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
    
  end subroutine compute_I4_mmp_s
  
  subroutine compute_source(asource_s, weightlag_s, Amm_s)
    real(kind=dp), dimension(:,:) :: asource_s, weightlag_s, Amm_s

    if (allocated(M_transform_inv)) deallocate(M_transform_inv)
    allocate(M_transform_inv(0:lagmax, 0:lagmax))
    if (allocated(M_transform)) deallocate(M_transform)
    allocate(M_transform(0:lagmax, 0:lagmax))
    
    call compute_Minv(M_transform_inv)
    call compute_sources(asource_s, weightlag_s)
    Amm_s=M_transform

  end subroutine compute_source
  
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

    call compute_lorentz(anumm_s)
    call compute_energyscattering(denmm_s)
    call compute_integralpart(ailmm_s)

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
    write (*,'(A,A,A,A,A,1E13.6)') " Computing collision operator for ", tag_a, "-", tag_b, " with gamma_ab =", gamma_ab

    call compute_lorentz(anumm_s)
    call compute_lorentz_inf(anumm_inf_s)
    call compute_energyscattering(denmm_s)
    call compute_integralpart(ailmm_s)

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
    write (*,'(A,A,A,A,A,1E13.6)') " Computing collision operator for ", tag_a, "-", tag_b, " with gamma_ab =", gamma_ab

    call compute_lorentz(anumm_s)

  end subroutine compute_collop_lorentz
  
end module collop_compute
