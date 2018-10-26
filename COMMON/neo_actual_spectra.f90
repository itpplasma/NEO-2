!> Actual spectra on flux surface
module neo_actual_spectra
  use nrtype, only : dp
  real(kind=dp),    dimension(:),     allocatable :: s_rmnc, s_zmnc, s_lmnc
  real(kind=dp),    dimension(:),     allocatable :: s_bmnc
  !! Modifications by Andreas F. Martitsch (06.08.2014)
  ! Additional data from Boozer files without Stellarator symmetry
  real(kind=dp),    dimension(:),     allocatable :: s_rmns, s_zmns, s_lmns
  real(kind=dp),    dimension(:),     allocatable :: s_bmns
  !! End Modifications by Andreas F. Martitsch (06.08.2014)
end module neo_actual_spectra
