module nrtype
  public I4B, I2B, I1B
  public longint
  public SP, DP, SPC, DPC
  public LGT
  public PI, PIO2, TWOPI, SQRT2, EULER, PI_D, PIO2_D, TWOPI_D
  public sprs2_sp, sprs2_dp

  ! Definition of types taken from Numerical Recipes
  integer, parameter :: I4B = SELECTED_INT_KIND(9)
  integer, parameter :: I2B = SELECTED_INT_KIND(4)
  integer, parameter :: I1B = SELECTED_INT_KIND(2)
  ! lahey, g95, intel
  integer, parameter    :: longint = 8
  ! nag
  !integer, parameter    :: longint = 4
  integer, parameter :: SP = kind(1.0)
  integer, parameter :: DP = kind(1.0D0)
  integer, parameter :: SPC = kind((1.0,1.0))
  integer, parameter :: DPC = kind((1.0D0,1.0D0))
  integer, parameter :: LGT = kind(.TRUE.)

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

  type sprs2_sp
     integer(I4B) :: n,len
     real(SP), dimension(:), pointer :: val
     integer(I4B), dimension(:), pointer :: irow
     integer(I4B), dimension(:), pointer :: jcol
  end type sprs2_sp

  type sprs2_dp
     integer(I4B) :: n,len
     real(DP), dimension(:), pointer :: val
     integer(I4B), dimension(:), pointer :: irow
     integer(I4B), dimension(:), pointer :: jcol
  end type sprs2_dp
end module nrtype
