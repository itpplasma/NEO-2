MODULE neo_van
  use nrtype, only : dp
  REAL(kind=dp)                            :: v_phi0, v_theta0
  REAL(kind=dp)                            :: bmin_tol
  INTEGER                                  :: v_nper, v_steps
  INTEGER                                  :: v_num_mm
  INTEGER                                  :: no_minima
  INTEGER,  DIMENSION(:), ALLOCATABLE      :: li_minima
  INTEGER                                  :: no_gamma
  INTEGER                                  :: tau_num
  INTEGER                                  :: tau_max_iter
  REAL(kind=dp)                            :: lambda_fac
  REAL(kind=dp)                            :: temp_e
  REAL(kind=dp)                            :: gamma_eps
  REAL(kind=dp)                            :: phi_eps
END MODULE neo_van
