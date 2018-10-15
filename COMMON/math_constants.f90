!> Despite the name this module contains mathematical and physical
!> constants, as it is sometimes hard to distinguish, e.g. for
!> conversion factors.
module math_constants
  use nrtype, only : DP, dpc

  real(DP), parameter :: PI=3.141592653589793238462643383279502884197_dp
  real(DP), parameter :: PIO2=1.57079632679489661923132169163975144209858_dp
  real(DP), parameter :: TWOPI=6.283185307179586476925286766559005768394_dp
  real(DP), parameter :: SQRT2=1.41421356237309504880168872420969807856967_dp
  !> Defined as \f$ \gamma = \int\limits_{1}^{\infty} 1/floor(x) - 1/x dx \f$
  !> or alternatively \f$ \sum\limits_{k=1}^{\infty} ( 1/k - \ln ((k+1)/k)\f$
  !> or \f$ \Gamma'(1) = - \gamma \f$.
  real(DP), parameter :: EULER=0.5772156649015328606065120900824024310422_dp

  real(DP), parameter :: PI_D=3.141592653589793238462643383279502884197_dp
  real(DP), parameter :: PIO2_D=1.57079632679489661923132169163975144209858_dp
  real(DP), parameter :: TWOPI_D=6.283185307179586476925286766559005768394_dp

  complex(DPC), parameter :: IMUN = (0.d0,1.d0)


  ! Define physical constants (cgs-units)
  REAL(kind=dp), PARAMETER, PUBLIC :: c=2.9979e10_dp       ! speed of light
  REAL(kind=dp), PARAMETER, PUBLIC :: e=4.8032e-10_dp      ! elementary charge
  REAL(kind=dp), PARAMETER, PUBLIC :: u=1.660539040e-24_dp ! atomic mass unit

  REAL(kind=dp), PARAMETER                 :: mc_o_e        = 5.6856793d-8 ! cgs
  REAL(kind=dp), PARAMETER                 :: mc_o_2e       = 2.8428397d-8 ! cgs
  REAL(kind=dp), PARAMETER                 :: e_o_mc        = 1.7588048d7  ! cgs
  REAL(kind=dp), PARAMETER                 :: b_convfac     = 1.0d4  ! SI to cgs
  REAL(kind=dp), PARAMETER                 :: i_convfac     = 1.0d6  ! SI to cgs
  REAL(kind=dp), PARAMETER                 :: sqg11_convfac = 1.0d6  ! SI to cgs
  REAL(kind=dp), PARAMETER                 :: len_convfac   = 1.0d2  ! SI to cgs
  REAL(kind=dp), PARAMETER                 :: te_to_vte     = 4.19d7 ! ev to cm/s

end module math_constants
