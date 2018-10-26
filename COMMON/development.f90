module development
  use nrtype, only : dp

  integer :: iprintflag,proptag_old=-100000000,nstfp
  real(kind=dp) :: deleta
  real(kind=dp), dimension(:),   allocatable :: alam_l,alam_r
  real(kind=dp), dimension(:),   allocatable :: delta_eta_l,delta_eta_r
  real(kind=dp), dimension(:),   allocatable :: phi_stfp
  real(kind=dp), dimension(:,:), allocatable :: stfp_ov_hstep
  integer :: solver_talk
  integer :: switch_off_asymp
  ! accuracy factor in ripple solver
  real(kind=dp) :: ripple_solver_accurfac
  ! asymptotical stuff
  integer :: asymp_margin_zero
  integer :: asymp_margin_npass
  real(kind=dp) :: asymp_pardeleta
end module
