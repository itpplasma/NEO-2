!-----------------------------------------------------------------------------------------!
! module: gsl_specialfunctions_mod                                                        !
! authors: TU Graz ITPcp Plasma, Andreas F. Martitsch                                            !
! date: 14.03.2017                                                                        !
! version: 0.2                                                                            !
!-----------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------!
! History                                                                                 !
! 0.1 - Initial version for modified Bessel function of second kind                       !
! 0.2 - Backend moved from FGSL to fortnum (bessel_kn)                                     !
!-----------------------------------------------------------------------------------------!

module gsl_specialfunctions_mod
  use, intrinsic :: iso_fortran_env, only: wp => real64
  use fortnum_special, only: bessel_kn

  implicit none

  contains

  function besselk(n, z)
    integer, intent(in) :: n
    real(wp), intent(in) :: z
    real(wp) :: besselk

    besselk = bessel_kn(n, z)

  end function besselk

end module gsl_specialfunctions_mod
