SUBROUTINE rk4_kin(x,h)

  USE rk4_kin_mod
  USE odeint_robust_mod, only: odeint_robust
  USE rhs_kin_sub, only: rhs_kin

  IMPLICIT NONE

  double precision, parameter :: eps = 1.d-10

  INTEGER          :: ierr
  DOUBLE PRECISION :: x,h,xh

  xh=x+h

  call odeint_robust(y,x,xh,eps,rhs_kin,ierr)
  if (ierr .ne. 0) then
    print *, 'rk4_kin: field line integration failed from phi = ', x, &
         ' to phi = ', xh
    error stop 'rk4_kin: field line integration failed'
  end if

  x=xh

END SUBROUTINE rk4_kin
