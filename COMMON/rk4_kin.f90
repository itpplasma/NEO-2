SUBROUTINE rk4_kin(x,h)

  USE rk4_kin_mod

  IMPLICIT NONE

  double precision, parameter :: eps = 1.d-10

  INTEGER          :: n
  DOUBLE PRECISION :: x,h,xh

  external :: rhs_kin

  n=SIZE(y)
  xh=x+h

  call odeint_allroutines(y,n,x,xh,eps,rhs_kin)

  x=xh

END SUBROUTINE rk4_kin
