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

  complex(DPC), parameter :: IMUN = (0.d0,1.d0)

  ! Define physical constants (cgs-units)
  real(kind=dp), parameter, public :: c=2.9979e10_dp       !> speed of light, [cm/s]
  real(kind=dp), parameter, public :: e=4.8032e-10_dp      !> elementary charge
  real(kind=dp), parameter, public :: u=1.660539040e-24_dp !> atomic mass unit

  ! Masses of some particles/atoms. Atoms get designated with m\_ and
  ! symbol (e.g. H for hydrogen).
  real(kind=dp), parameter, public :: m_ele = 9.109382150d-28 !> Electron mass
  real(kind=dp), parameter, public :: m_pro = 1.672621637d-24 !> Proton mass
  real(kind=dp), parameter, public :: m_alp = 6.644656200d-24 !> alpha particle mass
  real(kind=dp), parameter, public :: m_D   = 3.343583719d-24 !> Deuterium mass
  real(kind=dp), parameter, public :: m_C   = 19.94406876d-24 !> Carbon mass
  real(kind=dp), parameter, public :: m_W   = 305.2734971d-24 !> Tungsten mass


  real(kind=dp), parameter, public :: mc_o_e        = 5.6856793d-8 !> cgs
  real(kind=dp), parameter, public :: mc_o_2e       = 2.8428397d-8 !> cgs
  real(kind=dp), parameter, public :: e_o_mc        = 1.7588048d7  !> cgs
  real(kind=dp), parameter, public :: b_convfac     = 1.0d4  !> SI to cgs for magnetic fields, i.e. convert from T to Gauss.
  real(kind=dp), parameter, public :: i_convfac     = 1.0d6  !> SI to cgs
  real(kind=dp), parameter, public :: sqg11_convfac = 1.0d6  !> SI to cgs
  real(kind=dp), parameter, public :: len_convfac   = 1.0d2  !> SI to cgs for length, i.e. convert from m to cm.
  real(kind=dp), parameter, public :: te_to_vte     = 4.19d7 !> Thermal energy to thermal velocity, i.e. eV to cm/s.
  real(kind=dp), parameter, public :: ev_to_cgs     = 1.6022d-12 !> eV to gcm^2/s^2, the cgs equivalent of joule.

end module math_constants
