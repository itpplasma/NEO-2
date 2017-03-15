!-----------------------------------------------------------------------------------------!
! module: gsl_specialfunctions_mod                                                        !
! authors: Gernot Kapper, Andreas F. Martitsch                                            !
! date: 14.03.2017                                                                        !
! version: 0.1                                                                            !
!-----------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------!
! History                                                                                 !
! 0.1 - Initial version for modified Bessel function of second kind                       !
!-----------------------------------------------------------------------------------------!

module gsl_specialfunctions_mod
  use fgsl

  implicit none

  contains
  
  function besselk(n, z)
    integer :: n
    real(fgsl_double) :: z, besselk
    
    besselk = fgsl_sf_bessel_kcn(n, z)
    
  end function besselk
  
end module gsl_specialfunctions_mod
