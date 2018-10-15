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
