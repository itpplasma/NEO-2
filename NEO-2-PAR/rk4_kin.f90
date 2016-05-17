!
  SUBROUTINE rk4_kin(x,h)
!
  USE rk4_kin_mod
!
  IMPLICIT NONE
!
  INTEGER          :: n,i
  DOUBLE PRECISION :: x,h,hh,h6,xh
!
  n=SIZE(y)
!
  hh=h*0.5d0
  h6=h/6.d0
  xh=x+hh
!
  CALL rhs_kin(x,y,dydx)
!
  yt=y+hh*dydx
!
  CALL rhs_kin(xh,yt,dyt)
!
  yt=y+hh*dyt
!
  CALL rhs_kin(xh,yt,dym)
!
  yt=y+h*dym
  dym=dyt+dym
!
  CALL rhs_kin(x+h,yt,dyt)
!
  y=y+h6*(dydx+dyt+2.d0*dym)
!
  x=x+h
!
  RETURN
  END
