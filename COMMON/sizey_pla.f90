MODULE sizey_pla
  USE neo_precision
  ! Definition for rk4d_pla
  INTEGER            ::  npart_pla
  INTEGER            ::  ndim_pla
  INTEGER, PARAMETER ::  npq_pla = 3
  REAL(kind=dp)      ::  lamup_pla
  REAL(kind=dp)      ::  lambda_alpha
  REAL(kind=dp)      ::  nufac_pla
END MODULE sizey_pla
