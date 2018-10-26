module rk4_kin_mod
  use nrtype, only : dp

  real(kind=dp), dimension(:), allocatable :: y,dydx,yt,dyt,dym
end module
