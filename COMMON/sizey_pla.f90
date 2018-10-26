module sizey_pla
  use nrtype, only : dp

  ! Definition for rk4d_pla
  integer            ::  npart_pla
  integer            ::  ndim_pla
  integer, parameter ::  npq_pla = 3
  real(kind=dp)      ::  lamup_pla
  real(kind=dp)      ::  lambda_alpha
  real(kind=dp)      ::  nufac_pla
end module sizey_pla
