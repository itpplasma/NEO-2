!
MODULE polleg_mod

  PUBLIC  polleg
  PRIVATE polleg_1
  INTERFACE polleg
     MODULE PROCEDURE polleg_1
  END INTERFACE

  PUBLIC  binomial
  PRIVATE binomial_1
  INTERFACE binomial
     MODULE PROCEDURE binomial_1
  END INTERFACE

CONTAINS
!
  SUBROUTINE polleg_1(n,coefleg)
!
! Computes coefficients of Legendre polynomials of orders from 0 to n

!
! Input parameters:
!           Formal: n            - maximum order of Legendre polynomials
! Output parameters:
!           Formal: coefleg(i,j) - j-th coefficient of Legendre polynomial
!                                  of the order i
  IMPLICIT NONE
!
  INTEGER :: n,i,j
!
  DOUBLE PRECISION :: frontfac,rearfac
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: coefleg
!
  IF(ALLOCATED(coefleg)) DEALLOCATE(coefleg)
  ALLOCATE(coefleg(0:n,0:n))
!
  coefleg=0.d0
  coefleg(0,0)=1.d0
  coefleg(1,1)=1.d0
  frontfac=1.d0
!
  DO i=2,n
    frontfac=frontfac*(2.d0-1.d0/DBLE(i))
    rearfac=frontfac
    coefleg(i,i)=rearfac
    DO j=i-2,0,-2
      rearfac=-rearfac*DBLE(j+1)*DBLE(j+2)/DBLE(i-j)/DBLE(i+j+1)
      coefleg(i,j)=rearfac
    ENDDO
  ENDDO
!
  RETURN
END SUBROUTINE polleg_1
!
  SUBROUTINE binomial_1(n,coefbin)
!
! Computes binomial coefficients of orders from 0 to n
!
! Input parameters:
!           Formal: n            - maximum power of the binom
! Output parameters:
!           Formal: coefbin(i,j) - j-th coefficient of the binom
  IMPLICIT NONE
!
  INTEGER :: n,i,j
!
  DOUBLE PRECISION :: frontfac,factforw,factbackw
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: coefbin
!
  IF(ALLOCATED(coefbin)) DEALLOCATE(coefbin)
  ALLOCATE(coefbin(0:n,0:n))
!
  coefbin=0.d0
  coefbin(0,0)=1.d0
  frontfac=1.d0
!
  DO i=1,n
    frontfac=frontfac*DBLE(i)
    factforw=1.d0
    factbackw=frontfac*DBLE(i+1)
    DO j=0,i
      IF(j.GT.0) factforw=factforw*DBLE(j)
      factbackw=factbackw/DBLE(i-j+1)
      coefbin(i,j)=frontfac/factforw/factbackw
    ENDDO
  ENDDO
!
  RETURN
END SUBROUTINE binomial_1
END MODULE
