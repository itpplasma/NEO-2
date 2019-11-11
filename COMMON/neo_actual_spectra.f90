MODULE neo_actual_spectra
! Actual spectra on flux surface
  USE neo_precision
  REAL(kind=dp),    DIMENSION(:),     ALLOCATABLE :: s_rmnc, s_zmnc, s_lmnc
  REAL(kind=dp),    DIMENSION(:),     ALLOCATABLE :: s_bmnc
  !! Modifications by Andreas F. Martitsch (06.08.2014)
  ! Additional data from Boozer files without Stellarator symmetry
  REAL(kind=dp),    DIMENSION(:),     ALLOCATABLE :: s_rmns, s_zmns, s_lmns
  REAL(kind=dp),    DIMENSION(:),     ALLOCATABLE :: s_bmns
  !! End Modifications by Andreas F. Martitsch (06.08.2014)
END MODULE neo_actual_spectra
