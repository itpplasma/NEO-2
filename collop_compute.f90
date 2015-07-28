module collop_compute

  use collop_definitions
  use hdf5_tools
  
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

  !**********************************************************
  ! Profile
  !**********************************************************
  real(kind=dp) :: n_a
  real(kind=dp) :: n_b
  real(kind=dp) :: T_a
  real(kind=dp) :: T_b
  real(kind=dp) :: v_ta
  real(kind=dp) :: v_tb

contains

  subroutine init_collop
    call init_laguerre(lagmax, legmax)
    call init_legendre(legmax)
    call init_phi(lagmax)

    v_ta = sqrt(2*T_a / m_a)
    v_tb = sqrt(2*T_b / m_b)
    gamma_ab = v_ta/v_tb

    write (*,'(A,E13.6)') " Computing collision operator for gamma_ab =", gamma_ab
    
  end subroutine init_collop
  
  subroutine compute_sources
    real(kind=dp), dimension(2) :: res_int
    integer :: m, k

    write (*,*) "Computing sources..."
    
    if (allocated(asource_s)) deallocate(asource_s)
    allocate(asource_s(0:lagmax, 1:3))
    
    do k = 1, 3
       do m = 0, lagmax
          res_int = fint1d_qagiu(integrand, 0d0, epsabs, epsrel)
          asource_s(m, k) = res_int(1)
       end do
    end do

    write (*,*) "Done."
    
    contains

      function integrand(x) result(y)
        real(kind=dp) :: x, y

        y = pi**(-3d0/2d0) * x**(4+alpha) * exp(-(1+beta)*x**2) * phi(m, x) * x**(2*k - 1 - 5*kdelta(3,k))
      end function integrand
     
    function kdelta(a,b)
      integer :: a, b
      integer :: kdelta

      kdelta = 0
      if (a .eq. b) kdelta = 1
    end function kdelta

  end subroutine compute_sources

  subroutine compute_lorentz
    real(kind=dp), dimension(2) :: res_int
    integer :: m, mp

    if (allocated(anumm_s)) deallocate(anumm_s)
    allocate(anumm_s(0:lagmax, 0:lagmax))

    write (*,*) "Computing Lorentz part..."

    do m = 0, lagmax
       do mp = 0, lagmax
          res_int = fint1d_qagiu(integrand, 0d0, epsabs, epsrel)
          anumm_s(m, mp) = 3d0/(4d0 * pi) * res_int(1) 
       end do
    end do
    
    write (*,*) "Done."

  contains

    function integrand(x)
      real(kind=dp) :: x, y
      real(kind=dp) :: integrand

      y = x * gamma_ab
      integrand = x**(alpha) * exp(-(1+beta)*x**2) * phi(m, x) * (erf(y) - G(y)) * phi(mp, x)

    end function integrand
    
  end subroutine compute_lorentz
   
  subroutine compute_energyscattering
    integer :: m, mp
    real(kind=dp), dimension(2) :: res_int

    if (allocated(denmm_s)) deallocate(denmm_s)
    allocate(denmm_s(0:lagmax, 0:lagmax))

    do m = 0, lagmax
       do mp = 0, lagmax
          res_int = fint1d_qagiu(integrand, 0d0, epsabs, epsrel)
          !res_int = fint1d_qag(integrand, 0d0, 100d0, epsabs, epsrel, 2)       
          denmm_s(m, mp) = 3d0/(4d0 * pi) * res_int(1)

          !write (*,*) "denmm_s", m, mp, integrand(2d0), denmm_s(m, mp)
          !write (*,*) G(2d0), d_G(2d0), gamma_ab, phi(m, 2d0), d_phi(m, 2d0), dd_phi(m,2d0)
       end do
    end do
    
  contains
    
    function integrand(x)
      real(kind=dp) :: x, y
      real(kind=dp) :: integrand
      real(kind=dp) :: D_1, D_2
      
      y = x * gamma_ab

      D_1 = d_G(y) * gamma_ab * exp(-x**2) * x * d_phi(mp, x) &
           + G(y) * (-2*x * exp(-x**2) * x * d_phi(mp, x) + exp(-x**2) * x * dd_phi(mp, x) &
           + exp(-x**2) * d_phi(mp, x))

      D_2 = 2 * (1 - T_a/T_b) * (d_G(y) * gamma_ab * x**2 * exp(-x**2) * phi(mp, x) + &
           G(y) * (2*x*exp(-x**2)*phi(mp, x) - 2*x**3*exp(-x**2)*phi(mp, x) + x**2 * exp(-x**2) * d_phi(mp, x)))
      
      integrand = x**(1+alpha)*exp(-beta*x**2) * phi(m, x) * (D_1 - D_2)
    end function integrand
  end subroutine compute_energyscattering

  subroutine compute_integralpart

    write (*,*) "Computing Integral part..."

    call compute_I1_mmp_s()
    call compute_I2_mmp_s()
    call compute_I3_mmp_s()
    call compute_I4_mmp_s()
    
    if (allocated(ailmm_s)) deallocate(ailmm_s)
    allocate(ailmm_s(0:lagmax, 0:lagmax, 0:legmax))

    ailmm_s = I1_mmp_s + I2_mmp_s + I3_mmp_s + I4_mmp_s

    write (*,*) "Done."
    
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

      integrand = x**(3+alpha) * exp(-(beta + 1)*x**2) * exp(-(gamma_ab*x)**2) * phi(m, x) * phi(mp, x)
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
      
      I_phi = x**(3+alpha) * exp(-(beta+1)*x**2) * phi(m,x) * (x**(-l-1) * I1 + x**l * I2)
      
    end function I_phi
    
    function I_phi_1(xp)
      real(kind=dp) :: xp, yp, I_phi_1

      yp = gamma_ab * xp
      I_phi_1 = exp(-yp**2) * phi(mp, xp) * xp**(l+2)
    end function I_phi_1
    
    function I_phi_2(xp)
      real(kind=dp) :: xp, yp, I_phi_2

      yp = gamma_ab * xp
      I_phi_2 = exp(-yp**2) * phi(mp, xp) * xp**(-l+1)
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
      
      I_phi = ((x**(3+alpha)) * exp(-(beta+1)*(x**2)) * (alpha-2*(beta+1)*(x**2)+4) * phi(m,x) + &
               (x**(4+alpha))  * exp(-(beta+1)*(x**2)) * d_phi(m,x)) * &
               ((x**(-l-1)) * I1 + (x**l) * I2)

      
    end function I_phi

    function I_phi_1(xp)
      real(kind=dp) :: xp, yp, I_phi_1

      yp = gamma_ab * xp
      I_phi_1 = exp(-yp**2) * phi(mp, xp) * xp**(l+2)
    end function I_phi_1
    
    function I_phi_2(xp)
      real(kind=dp) :: xp, yp, I_phi_2

      yp = gamma_ab * xp
      I_phi_2 = exp(-yp**2) * phi(mp, xp) * xp**(-l+1)
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
      
      I_psi = exp(-c*x**2)*x**(3+alpha) * ((20+9*alpha + alpha**2 - 2*(11+2*alpha) * c*x**2 + 4*c**2*x**4) * phi(m,x) + &
              x*(2*(5+alpha-2*c*x**2) * d_phi(m, x) + x*dd_phi(m,x))) * &
              (x**(-l-1)/(2*l+3)*I1 - x**(-l+1)/(2*l-1)*I2 + x**(l+2)/(2*l+3)*I3 - x**l/(2*l-1)*I4)
      
    end function I_psi

    function I_psi_1(xp)
      real(kind=dp) :: xp, yp, I_psi_1
      
      yp = gamma_ab * xp
      I_psi_1 = xp**(l+4)*exp(-yp**2)*phi(mp,xp)
    end function I_psi_1
    
    function I_psi_2(xp)
      real(kind=dp) :: xp, yp, I_psi_2

      yp = gamma_ab * xp
      I_psi_2 = xp**(l+2)*exp(-yp**2)*phi(mp,xp)
    end function I_psi_2

    function I_psi_3(xp)
      real(kind=dp) :: xp, yp, I_psi_3

      yp = gamma_ab * xp
      I_psi_3 = xp**(-l+1)*exp(-yp**2)*phi(mp,xp)
    end function I_psi_3
    
    function I_psi_4(xp)
      real(kind=dp) :: xp, yp, I_psi_4

      yp = gamma_ab * xp
      I_psi_4 = xp**(-l+3)*exp(-yp**2)*phi(mp,xp)
    end function I_psi_4
    
  end subroutine compute_I4_mmp_s

  subroutine write_collop(h5filename)
    character(len=*) :: h5filename
    integer(HID_T) :: h5id_collop, h5id_meta, h5id_species

    call h5_open_rw(h5filename, h5id_collop)
    call h5_define_group(h5id_collop, trim(tag_a) //'-'// trim(tag_b), h5id_species)
    call h5_define_group(h5id_species, 'meta', h5id_meta)

    call h5_add(h5id_meta, 'lagmax', lagmax)
    call h5_add(h5id_meta, 'legmax', legmax)
    call h5_add(h5id_meta, 'alpha', alpha)
    call h5_add(h5id_meta, 'beta', beta)
    call h5_add(h5id_meta, 'm_a', m_a)
    call h5_add(h5id_meta, 'm_b', m_b)
    call h5_add(h5id_meta, 'T_a', T_a)
    call h5_add(h5id_meta, 'T_b', T_b)
    call h5_add(h5id_meta, 'gamma_ab', gamma_ab)
    call h5_add(h5id_meta, 'phi_m_desc', phi_m_desc)
    call h5_add(h5id_meta, 'phi_m_tag', phi_m_tag)
    call h5_add(h5id_meta, 'tag_a', tag_a)
    call h5_add(h5id_meta, 'tag_b', tag_b)
    call h5_close_group(h5id_meta)

    call h5_add(h5id_species, 'a_m', asource_s, lbound(asource_s), ubound(asource_s))
    call h5_add(h5id_species, 'nu_hat_mmp', anumm_s, lbound(anumm_s), ubound(anumm_s))
    call h5_add(h5id_species, 'D_hat_mmp', denmm_s, lbound(denmm_s), ubound(denmm_s))
    call h5_add(h5id_species, 'I_lmmp', ailmm_s, lbound(ailmm_s), ubound(ailmm_s))
    call h5_add(h5id_species, 'I1_mmp', I1_mmp_s, lbound(I1_mmp_s), ubound(I1_mmp_s))
    call h5_add(h5id_species, 'I2_mmp', I2_mmp_s, lbound(I2_mmp_s), ubound(I2_mmp_s))
    call h5_add(h5id_species, 'I3_mmp', I3_mmp_s, lbound(I3_mmp_s), ubound(I3_mmp_s))
    call h5_add(h5id_species, 'I4_mmp', I4_mmp_s, lbound(I4_mmp_s), ubound(I4_mmp_s))
    call h5_close_group(h5id_species)

    call h5_close(h5id_collop)

  end subroutine write_collop
  
  subroutine compute_collop()
    
    call init_collop()

    call compute_sources()
    call compute_lorentz()
    call compute_energyscattering()
    call compute_integralpart()
    
  end subroutine compute_collop
end module collop_compute
