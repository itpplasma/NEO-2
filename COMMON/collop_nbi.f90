module collop_nbi
  use nrtype, only : dp, pi

  implicit none
  
contains

  subroutine init_phi_nbi(lagmax, legmax)
    integer :: lagmax, legmax

    ! Load data
    
  end subroutine init_phi_nbi

  function phi_nbi(m, x) result(phi)
    integer       :: m
    real(kind=dp) :: x, phi

    ! Return interpolated NBI source data
    
  end function phi_nbi

  function d_phi_nbi(m, x) result(d_phi)
    integer       :: m
    real(kind=dp) :: x, d_phi

    ! Return first derivative of phi_nbi
    
  end function d_phi_nbi

  function dd_phi_nbi(m, x) result(dd_phi)
    integer       :: m
    real(kind=dp) :: x, dd_phi
    
  end function dd_phi_nbi

end module collop_nbi
