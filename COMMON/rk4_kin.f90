SUBROUTINE rk4_kin(x,h)

  USE rk4_kin_mod
  USE odeint_allroutines_sub, only: odeint_allroutines, odeint_has_failed
  USE rhs_kin_sub, only: rhs_kin

  IMPLICIT NONE

  double precision, parameter :: eps = 1.d-10

  INTEGER          :: n
  DOUBLE PRECISION :: x,h,xh
  DOUBLE PRECISION :: scale_floor(SIZE(y))

  n=SIZE(y)
  xh=x+h

  ! The running integrals y(11:14) start at zero with zero derivative and
  ! need an absolute error scale, otherwise odeint can fail at roundoff.
  scale_floor = 0.d0
  scale_floor(11:14) = 1.d0

  call odeint_allroutines(y,n,x,xh,eps,rhs_kin,scale_floor=scale_floor)
  if (odeint_has_failed()) then
    print *, 'rk4_kin: field line integration failed from phi = ', x, &
         ' to phi = ', xh
    error stop 'rk4_kin: field line integration failed'
  end if

  x=xh

END SUBROUTINE rk4_kin
