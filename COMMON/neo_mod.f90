!!$  MODULE ripple_solver_mod
!!$    INTEGER :: npart_loc
!!$    INTEGER :: npass_l,npass_r,mhill
!!$    INTEGER :: npart_halfband
!!$    DOUBLE PRECISION :: qflux,qcurr,qflux_new,qcurr_new
!!$    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE   :: flux_p,flux_m
!!$    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE   :: curr_p,curr_m
!!$    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE   :: source_p,source_m
!!$    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: amat_plus_minus
!!$    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: amat_minus_plus
!!$    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: amat_plus_plus
!!$    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: amat_minus_minus
!!$  END MODULE ripple_solver_mod

! Winny module for ripple added
MODULE ripple_mod
  INTEGER :: ripple_counter
  INTEGER,          DIMENSION(:), ALLOCATABLE :: imax1_ripple
  INTEGER,          DIMENSION(:), ALLOCATABLE :: imax2_ripple
  INTEGER,          DIMENSION(:), ALLOCATABLE :: imin_ripple
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: col_ripple
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: eta_x0
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: eta_s
  DOUBLE PRECISION                            :: bmax_abs,bmin_abs
  INTEGER :: propagator_counter
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: phi1_prop
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: phi2_prop  
  INTEGER, DIMENSION(:), ALLOCATABLE :: ibeg_prop,iend_prop
  INTEGER, DIMENSION(:), ALLOCATABLE :: ripplenumber_prop
  INTEGER :: ibeg
  INTEGER :: iend
  INTEGER :: ibeg_prop_act
  INTEGER :: iend_prop_act
  INTEGER :: imax_eta
END MODULE ripple_mod

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

!!$MODULE join_ripples_simple_mod
!!$  INTEGER :: o_l,o_r,n_l,n_r,n_laguerra
!!$  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: o_amat_p_p,n_amat_p_p
!!$  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: o_amat_p_m,n_amat_p_m
!!$  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: o_amat_m_m,n_amat_m_m
!!$  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: o_amat_m_p,n_amat_m_p
!!$  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: o_source_p,n_source_p
!!$  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: o_source_m,n_source_m
!!$  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: o_flux_p,n_flux_p
!!$  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: o_flux_m,n_flux_m
!!$  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: o_qflux,n_qflux
!!$END MODULE join_ripples_simple_mod
