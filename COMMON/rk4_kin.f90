SUBROUTINE rk4_kin(x,h)

  USE rk4_kin_mod
  USE odeint_allroutines_sub, only: odeint_allroutines, compute_derivative
  USE rhs_kin_sub, only: rhs_kin

  IMPLICIT NONE

  double precision, parameter :: eps = 1.d-10

  INTEGER          :: n
  DOUBLE PRECISION :: x,h,xh

  n=SIZE(y)
  xh=x+h

  call odeint_allroutines(y,n,x,xh,eps,rhs_kin)

  x=xh

END SUBROUTINE rk4_kin
