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

end module math_constants
