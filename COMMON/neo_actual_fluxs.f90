module neo_actual_fluxs
! Actual data on flux surface
  use nrtype, only : dp
  real(kind=dp)  :: s_es
  real(kind=dp)  :: s_iota
  real(kind=dp)  :: s_pprime
  real(kind=dp)  :: s_sqrtg00
  real(kind=dp)  :: s_curr_tor, s_curr_pol
  real(kind=dp)  :: s_b00, s_b00_s
  integer        :: through_fourier
end module neo_actual_fluxs
