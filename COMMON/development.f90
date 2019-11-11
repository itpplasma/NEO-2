MODULE development
  INTEGER :: iprintflag,proptag_old=-100000000,nstfp
  DOUBLE PRECISION :: deleta
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: alam_l,alam_r
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: delta_eta_l,delta_eta_r
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: phi_stfp
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: stfp_ov_hstep
  INTEGER :: solver_talk
  INTEGER :: switch_off_asymp
  ! accuracy factor in ripple solver
  DOUBLE PRECISION :: ripple_solver_accurfac
  ! asymptotical stuff
  INTEGER :: asymp_margin_zero
  INTEGER :: asymp_margin_npass
  DOUBLE PRECISION :: asymp_pardeleta
END MODULE development
