MODULE nrtype
! Definition of types taken from Numerical Recipes
  INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
  INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
  INTEGER, PARAMETER :: I1B = SELECTED_INT_KIND(2)
  integer, parameter :: longint = 8 !< \todo Replace with one of the above.
  INTEGER, PARAMETER :: SP = KIND(1.0)
  INTEGER, PARAMETER :: DP = KIND(1.0D0)
  INTEGER, PARAMETER :: SPC = KIND((1.0,1.0))
  INTEGER, PARAMETER :: DPC = KIND((1.0D0,1.0D0))
  INTEGER, PARAMETER :: LGT = KIND(.TRUE.)
  REAL(DP), PARAMETER :: PI=3.141592653589793238462643383279502884197_dp
  REAL(DP), PARAMETER :: PIO2=1.57079632679489661923132169163975144209858_dp
  REAL(DP), PARAMETER :: TWOPI=6.283185307179586476925286766559005768394_dp
  REAL(DP), PARAMETER :: SQRT2=1.41421356237309504880168872420969807856967_dp
  REAL(DP), PARAMETER :: EULER=0.5772156649015328606065120900824024310422_dp
  REAL(DP), PARAMETER :: PI_D=3.141592653589793238462643383279502884197_dp
  REAL(DP), PARAMETER :: PIO2_D=1.57079632679489661923132169163975144209858_dp
  REAL(DP), PARAMETER :: TWOPI_D=6.283185307179586476925286766559005768394_dp
  !> Type for sparse quadratic matrix with single precision entries.
  !> Storing row/column index, as well as matrix size (n from
  !> 'n x n Matrix') and number of entries.
  TYPE sprs2_sp
     INTEGER(I4B) :: n,len
     REAL(SP), DIMENSION(:), POINTER :: val
     INTEGER(I4B), DIMENSION(:), POINTER :: irow
     INTEGER(I4B), DIMENSION(:), POINTER :: jcol
  END TYPE sprs2_sp
  !> As sprs2_sp, but with double precision matrix values.
  TYPE sprs2_dp
     INTEGER(I4B) :: n,len
     REAL(DP), DIMENSION(:), POINTER :: val
     INTEGER(I4B), DIMENSION(:), POINTER :: irow
     INTEGER(I4B), DIMENSION(:), POINTER :: jcol
  END TYPE sprs2_dp

  !> Variable to be able to use linear interpolation (=true) for spline
  !> coefficients. Value true is used by nfp from tools/create_surfaces.
  logical :: splinecof_compatibility = .false.

END MODULE nrtype
