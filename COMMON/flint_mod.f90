! some Winny stuff about flint
MODULE flint_mod
  INTEGER :: plot_gauss
  INTEGER :: plot_prop
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: phiarr
  INTEGER, DIMENSION(:), ALLOCATABLE :: phi_divide
  DOUBLE PRECISION :: hphi_mult
  INTEGER :: phi_split_mode,phi_place_mode,phi_split_min
  INTEGER :: max_solver_try
  ! WINNY
  DOUBLE PRECISION :: bsfunc_local_err_max_mult
  DOUBLE PRECISION :: bsfunc_max_mult_reach
  DOUBLE PRECISION :: boundary_dist_limit_factor

  INTEGER :: bsfunc_modelfunc_num
  INTEGER :: bsfunc_divide
  INTEGER :: bsfunc_ignore_trap_levels

  DOUBLE PRECISION :: bsfunc_local_shield_factor
  LOGICAL :: bsfunc_shield
  LOGICAL :: bsfunc_lambda_loc_res
  DOUBLE PRECISION :: eta_savemem_dist1, eta_savemem_dist2, eta_savemem_sigma_mult
END MODULE flint_mod
