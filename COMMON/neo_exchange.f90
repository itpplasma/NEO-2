module neo_exchange
  use nrtype, only : dp

  real(kind=dp)              :: b_min, b_max
  integer                    :: nper
  real(kind=dp)              :: rt0, rt0_g
  real(kind=dp)              :: bmref, bmref_g, bmref_a
  integer                    :: nstep_per, nstep_min, nstep_max
  integer                    :: write_integrate
  integer                    :: write_diagnostic
  integer                    :: write_cur_inte
  integer                    :: write_pla_inte
  real(kind=dp)              :: acc_req
  integer                    :: no_bins
  integer                    :: psi_ind
  integer                    :: calc_nstep_max
  real(kind=dp)              :: theta_bmin, phi_bmin
  real(kind=dp)              :: theta_bmax, phi_bmax
  real(kind=dp),    dimension(:),     allocatable :: iota
  real(kind=dp),    dimension(:),     allocatable :: curr_pol
  real(kind=dp),    dimension(:),     allocatable :: curr_tor
  real(kind=dp)              :: fac
  integer                    :: calc_cur, calc_eps, calc_pla, calc_van
  integer                    :: hit_rat, nfp_rat, nfl_rat
  real(kind=dp)              :: delta_theta_rat
  real(kind=dp)              :: delta_cur_fac
  integer                    :: cutoff_cur_int
  integer                    :: write_cur_disp
end module neo_exchange
