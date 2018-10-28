
SUBROUTINE rk4_kin(x,h)

  use nrtype, only : dp
  USE rk4_kin_mod

  IMPLICIT NONE

  real(kind=dp), parameter :: eps = 1.d-10

  INTEGER          :: n
  real(kind=dp) :: x,h,xh

  external :: rhs_kin

  n=SIZE(y)
  xh=x+h

  call odeint_allroutines(y,n,x,xh,eps,rhs_kin)

  x=xh

  RETURN
END SUBROUTINE rk4_kin
