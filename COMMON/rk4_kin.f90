SUBROUTINE rk4_kin(x,h)

  USE rk4_kin_mod
  USE rhs_kin_sub, only: rhs_kin

  IMPLICIT NONE

  double precision, parameter :: eps = 1.d-10

  DOUBLE PRECISION :: x,h,xh
  LOGICAL          :: failed

  xh=x+h

  call odeint_kin_retry(y,x,xh,eps,rhs_kin,failed)
  if (failed) then
    print *, 'rk4_kin: field line integration failed from phi = ', x, &
         ' to phi = ', xh
    error stop 'rk4_kin: field line integration failed'
  end if

  x=xh

END SUBROUTINE rk4_kin
