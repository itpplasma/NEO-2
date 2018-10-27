module collop_polynomial
  use nrtype, only : dp

  implicit none
  
contains

  !**********************************************************
  ! Standard polynomials
  !**********************************************************
  subroutine init_phi_polynomial(lagmax, legmax)
    integer :: lagmax, legmax

  end subroutine init_phi_polynomial

  function phi_polynomial(m, x)
    integer       :: m
    real(kind=dp) :: x, phi_polynomial

    !phi_polynomial = sqrt(pi**1.5d0) * sqrt(2d0/gamma(5d0/2d0 + m)) * x**(m)
    phi_polynomial = x**m
  end function phi_polynomial

  function d_phi_polynomial(m, x)
    integer       :: m
    real(kind=dp) :: x, d_phi_polynomial

    !d_phi_polynomial =  sqrt(pi**1.5d0) * sqrt(2d0/gamma(5d0/2d0 + m)) * m * x**(m-1)
    d_phi_polynomial = m * x**(m-1)
  end function d_phi_polynomial

  function dd_phi_polynomial(m, x)
    integer :: m
    real(kind=dp) :: x, dd_phi_polynomial

    !dd_phi_polynomial =  sqrt(pi**1.5d0) * sqrt(2d0/gamma(5d0/2d0 + m)) * m*(m-1) * x**(m-2)
    dd_phi_polynomial = m*(m-1) * x**(m-2)
  end function dd_phi_polynomial

  !**********************************************************
  ! Square polynomial
  !**********************************************************
  subroutine init_phi_polynomial_2(lagmax, legmax)
    integer :: lagmax, legmax

  end subroutine init_phi_polynomial_2

  function phi_polynomial_2(m, x)
    integer       :: m
    real(kind=dp) :: x, phi_polynomial_2

    !phi_polynomial = sqrt(pi**1.5d0) * sqrt(2d0/gamma(5d0/2d0 + m)) * x**(m)
    phi_polynomial_2 = x**(2*m)
  end function phi_polynomial_2

  function d_phi_polynomial_2(m, x)
    integer       :: m
    real(kind=dp) :: x, d_phi_polynomial_2

    !d_phi_polynomial =  sqrt(pi**1.5d0) * sqrt(2d0/gamma(5d0/2d0 + m)) * m * x**(m-1)
    d_phi_polynomial_2 = (2*m) * x**(2*m-1)
  end function d_phi_polynomial_2

  function dd_phi_polynomial_2(m, x)
    integer :: m
    real(kind=dp) :: x, dd_phi_polynomial_2

    !dd_phi_polynomial =  sqrt(pi**1.5d0) * sqrt(2d0/gamma(5d0/2d0 + m)) * m*(m-1) * x**(m-2)
    dd_phi_polynomial_2 = 2*m*(2*m-1) * x**(2*m-2)
  end function dd_phi_polynomial_2

  !**********************************************************
  ! Square polynomial without zeroth order
  !**********************************************************
  subroutine init_phi_polynomial_3(lagmax, legmax)
    integer :: lagmax, legmax

  end subroutine init_phi_polynomial_3

  function phi_polynomial_3(m, x)
    integer       :: m
    real(kind=dp) :: x, phi_polynomial_3

    !phi_polynomial = sqrt(pi**1.5d0) * sqrt(2d0/gamma(5d0/2d0 + m)) * x**(m)
    phi_polynomial_3 = x**(2*(m+1))
  end function phi_polynomial_3

  function d_phi_polynomial_3(m, x)
    integer       :: m
    real(kind=dp) :: x, d_phi_polynomial_3

    !d_phi_polynomial =  sqrt(pi**1.5d0) * sqrt(2d0/gamma(5d0/2d0 + m)) * m * x**(m-1)
    d_phi_polynomial_3 = 2*(m+1) * x**(2*m+1)
  end function d_phi_polynomial_3

  function dd_phi_polynomial_3(m, x)
    integer :: m
    real(kind=dp) :: x, dd_phi_polynomial_3

    !dd_phi_polynomial =  sqrt(pi**1.5d0) * sqrt(2d0/gamma(5d0/2d0 + m)) * m*(m-1) * x**(m-2)
    dd_phi_polynomial_3 = 2*(m+1)*(2*m+1) * x**(2*m)
  end function dd_phi_polynomial_3
  
end module collop_polynomial
