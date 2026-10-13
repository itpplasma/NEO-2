MODULE rk4_kin_mod
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: y,dydx,yt,dyt,dym

CONTAINS

  !> odeint_allroutines from x1 to x2; on failure retry once with the zero
  !> components of y(11:14) shifted by one. Returns failed if the retry fails.
  SUBROUTINE odeint_kin_retry(yv,x1,x2,eps,derivs,failed)
    USE odeint_allroutines_sub, only: odeint_allroutines, odeint_has_failed

    DOUBLE PRECISION, INTENT(INOUT) :: yv(:)
    DOUBLE PRECISION, INTENT(IN)    :: x1,x2,eps
    LOGICAL,          INTENT(OUT)   :: failed
    INTERFACE
      SUBROUTINE derivs(x,y,dydx)
        DOUBLE PRECISION, INTENT(IN)  :: x
        DOUBLE PRECISION, INTENT(IN)  :: y(:)
        DOUBLE PRECISION, INTENT(OUT) :: dydx(:)
      END SUBROUTINE derivs
    END INTERFACE

    DOUBLE PRECISION :: y_start(SIZE(yv)), shift(11:14)

    y_start=yv
    call odeint_allroutines(yv,SIZE(yv),x1,x2,eps,derivs)
    failed = odeint_has_failed()
    if (failed) then
      ! The running integrals y(11:14) start at zero with zero derivative, so the
      ! relative error test fails down to roundoff and odeint gives up. They do
      ! not enter rhs_kin: retry with the zero ones shifted to an absolute scale.
      shift = merge(1.d0, 0.d0, y_start(11:14) .eq. 0.d0)
      yv = y_start
      yv(11:14) = yv(11:14) + shift
      call odeint_allroutines(yv,SIZE(yv),x1,x2,eps,derivs)
      yv(11:14) = yv(11:14) - shift
      failed = odeint_has_failed()
    end if
  END SUBROUTINE odeint_kin_retry

END MODULE rk4_kin_mod
