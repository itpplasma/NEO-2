module neo_output
  use nrtype, only : dp
  real(kind=dp), dimension(:), allocatable ::  epspar
  real(kind=dp)                            ::  epstot,ctrone,ctrtot
  real(kind=dp)                            ::  epstothat
  real(kind=dp)                            ::  bareph,barept,drdpsi
  real(kind=dp)                            ::  yps
  integer                                  ::  nintfp
  integer                                  ::  ierr
  real(kind=dp)                            ::  lambda_b
  real(kind=dp)                            ::  lambda_b1, lambda_b2
  real(kind=dp)                            ::  lambda_ps1, lambda_ps2
  real(kind=dp)                            ::  lambda_del
  real(kind=dp)                            ::  avnabpsi,rfint
  real(kind=dp)                            ::  avb2, f_c, f_p
  real(kind=dp)                            ::  lambda_pla
  real(kind=dp)                            ::  delta_cur_max, typ_cur_len
end module neo_output
