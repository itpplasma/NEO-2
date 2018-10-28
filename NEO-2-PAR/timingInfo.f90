module timinginfo
  use nrtype, only : dp

  real(kind=dp), public, save :: dgesvTime = 0
  integer, save :: dgesvCalls = 0

end module timinginfo
