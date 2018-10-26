! some Winny stuff about flint
module flint_mod
  use nrtype, only : dp

  integer :: plot_gauss
  integer :: plot_prop
  real(kind=dp), dimension(:), allocatable :: phiarr
  integer, dimension(:), allocatable :: phi_divide
  real(kind=dp) :: hphi_mult
  integer :: phi_split_mode,phi_place_mode,phi_split_min
  integer :: max_solver_try
  ! WINNY
  real(kind=dp) :: bsfunc_local_err_max_mult
  real(kind=dp) :: bsfunc_max_mult_reach
  real(kind=dp) :: boundary_dist_limit_factor

  integer :: bsfunc_modelfunc_num
  integer :: bsfunc_divide
  integer :: bsfunc_ignore_trap_levels

  real(kind=dp) :: bsfunc_local_shield_factor
  logical :: bsfunc_shield
  logical :: bsfunc_lambda_loc_res
  real(kind=dp) :: eta_savemem_dist1, eta_savemem_dist2, eta_savemem_sigma_mult

end module flint_mod
