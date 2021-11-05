MODULE neo_actual_spectra
! Actual spectra on flux surface
  use nrtype, only : dp
  REAL(kind=dp),    DIMENSION(:),     ALLOCATABLE :: s_rmnc, s_zmnc, s_lmnc
  REAL(kind=dp),    DIMENSION(:),     ALLOCATABLE :: s_bmnc
END MODULE neo_actual_spectra
