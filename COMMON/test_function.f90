
! testfunction for spline tf(x, m)
!           1. derivative tfp(x, m)
!           2. derivative tfpp(x, m)
!           3. derivative tfppp(x, m)

FUNCTION tf(x,m)
!
! calculate testfunction for spline tf(x, m)
!
!-----------------------------------------------------------------------
! Modules
!-----------------------------------------------------------------------

  USE inter_precision

!-----------------------------------------------------------------------

  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: x
  REAL(DP), INTENT(IN) :: m
  REAL(DP)             :: tf

!-----------------------------------------------------------------------

  IF (m .NE. 0.0D0) THEN
     IF (x == 0.0D0) THEN
        tf = 0.0D0
     ELSE
        tf = x**m
     END IF
  ELSE
     tf = 1.0D0
  END IF

END FUNCTION tf

FUNCTION tfp(x,m)
!
! calculate 1. derivative of testfunction for spline tfp(x, m)
!
!-----------------------------------------------------------------------
! Modules
!-----------------------------------------------------------------------

  USE inter_precision

!-----------------------------------------------------------------------

  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: x
  REAL(DP), INTENT(IN) :: m
  REAL(DP)             :: tfp

!-----------------------------------------------------------------------

  IF ((m - 1.0D0) .NE. 0.0D0) THEN
     IF (x == 0.0D0) THEN
        tfp = 0.0D0
     ELSE
        tfp = m * x**(m-1.0D0)
     END IF
  ELSE
     tfp = 1.0D0
  END IF

END FUNCTION tfp


FUNCTION tfpp(x,m)
!
! calculate 2. derivative of testfunction for spline tfpp(x, m)
!
!-----------------------------------------------------------------------
! Modules
!-----------------------------------------------------------------------

  USE inter_precision

!-----------------------------------------------------------------------

  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: x
  REAL(DP), INTENT(IN) :: m
  REAL(DP)             :: tfpp

!-----------------------------------------------------------------------

  IF ((m - 2.0D0) .NE. 0.0D0) THEN
     IF (x == 0.0D0) THEN
        tfpp = 0.0D0
     ELSE
        tfpp = m * (m - 1.0D0) * x**(m-2.0D0)
     END IF
  ELSE
     tfpp = 1.0D0
  END IF

END FUNCTION tfpp


FUNCTION tfppp(x,m)
!
! calculate 3. derivative of testfunction for spline tfppp(x, m)
!
!-----------------------------------------------------------------------
! Modules
!-----------------------------------------------------------------------

  USE inter_precision

!-----------------------------------------------------------------------

  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: x
  REAL(DP), INTENT(IN) :: m
  REAL(DP)             :: tfppp

!-----------------------------------------------------------------------

  IF ((m - 3.0D0) .NE. 0.0D0) THEN
     IF (x == 0.0D0) THEN
        tfppp = 0.0D0
     ELSE
        tfppp = m * (m - 1.0D0) * (m - 2.0D0) * x**(m-3.0D0)
     END IF
  ELSE
     tfppp = 1.0D0
  END IF

END FUNCTION tfppp

FUNCTION tfone(x,m)
! calculate testfunction for spline tf(x, m); here just 1
!-----------------------------------------------------------------------
! Modules
!-----------------------------------------------------------------------
  USE inter_precision
!-----------------------------------------------------------------------
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: x
  REAL(DP), INTENT(IN) :: m
  REAL(DP)             :: tfone
!-----------------------------------------------------------------------
  tfone = 1.0D0
END FUNCTION tfone
FUNCTION tfzero(x,m)
! calculate testfunction for spline tf(x, m); here just 1
!-----------------------------------------------------------------------
! Modules
!-----------------------------------------------------------------------
  USE inter_precision
!-----------------------------------------------------------------------
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: x
  REAL(DP), INTENT(IN) :: m
  REAL(DP)             :: tfzero
!-----------------------------------------------------------------------
  tfzero = 0.0D0
END FUNCTION tfzero
