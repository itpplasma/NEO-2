MODULE neo_actual_fluxs
! Actual data on flux surface
  use nrtype, only : dp
  REAL(kind=dp)                                   :: s_es
  REAL(kind=dp)                                   :: s_iota
  REAL(kind=dp)                                   :: s_pprime
  REAL(kind=dp)                                   :: s_sqrtg00
  REAL(kind=dp)                                   :: s_curr_tor, s_curr_pol
  REAL(kind=dp)                                   :: s_b00, s_b00_s
  INTEGER                                         :: through_fourier
END MODULE neo_actual_fluxs
