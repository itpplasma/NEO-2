!> Fast path implementation for natural cubic splines (m=0, sw1=2, sw2=4)
module splinecof3_fast_mod
  use nrtype, only : I4B, DP
  implicit none
  
contains

  !> Fast natural cubic spline implementation
  SUBROUTINE splinecof3_fast(x, y, c1, cn, lambda1, indx, sw1, sw2, &
       a, b, c, d, m, f)
    REAL(DP),                   INTENT(INOUT) :: c1, cn
    REAL(DP),     DIMENSION(:), INTENT(IN)    :: x, y, lambda1
    INTEGER(I4B), DIMENSION(:), INTENT(IN)    :: indx
    REAL(DP),     DIMENSION(:), INTENT(OUT)   :: a, b, c, d
    INTEGER(I4B),               INTENT(IN)    :: sw1, sw2
    REAL(DP),                   INTENT(IN)    :: m
    INTERFACE
       FUNCTION f(x,m)
         use nrtype, only : DP
         IMPLICIT NONE
         REAL(DP), INTENT(IN) :: x, m
         REAL(DP)             :: f
       END FUNCTION f
    END INTERFACE

    ! Local variables
    INTEGER(I4B) :: len_x, len_indx, i, j, n
    REAL(DP) :: h_i, h_i1
    REAL(DP), DIMENSION(:), ALLOCATABLE :: h, alpha, mu, z, l

    ! Check if fast path conditions are met
    IF (m /= 0 .OR. sw1 /= 2 .OR. sw2 /= 4) THEN
       WRITE(*,*) 'splinecof3_fast: Invalid conditions for fast path'
       WRITE(*,*) 'm=', m, ' sw1=', sw1, ' sw2=', sw2
       STOP 'Fast path requires m=0, sw1=2, sw2=4'
    END IF

    len_x = SIZE(x)
    len_indx = SIZE(indx)
    n = len_indx

    ! Allocate work arrays
    ALLOCATE(h(n), alpha(n), mu(n), z(n), l(n))

    ! Compute intervals h
    DO i = 1, n-1
       h(i) = x(indx(i+1)) - x(indx(i))
    END DO

    ! For natural splines, c1 = cn = 0
    c1 = 0.0D0
    cn = 0.0D0

    ! Compute alpha values
    alpha(1) = 0.0D0
    DO i = 2, n-1
       alpha(i) = 3.0D0/h(i) * (y(indx(i+1)) - y(indx(i))) - &
                  3.0D0/h(i-1) * (y(indx(i)) - y(indx(i-1)))
    END DO
    alpha(n) = 0.0D0

    ! Forward elimination
    l(1) = 1.0D0
    mu(1) = 0.0D0
    z(1) = 0.0D0

    DO i = 2, n-1
       l(i) = 2.0D0 * (x(indx(i+1)) - x(indx(i-1))) - h(i-1) * mu(i-1)
       mu(i) = h(i) / l(i)
       z(i) = (alpha(i) - h(i-1) * z(i-1)) / l(i)
    END DO

    l(n) = 1.0D0
    z(n) = 0.0D0

    ! Back substitution
    c(n) = 0.0D0
    DO j = n-1, 1, -1
       c(j) = z(j) - mu(j) * c(j+1)
    END DO

    ! Compute remaining coefficients
    DO i = 1, n-1
       h_i = h(i)
       a(i) = y(indx(i))
       b(i) = (y(indx(i+1)) - y(indx(i)))/h_i - h_i*(c(i+1) + 2.0D0*c(i))/3.0D0
       d(i) = (c(i+1) - c(i))/(3.0D0*h_i)
    END DO

    ! Last segment
    a(n) = y(indx(n))
    b(n) = 0.0D0
    c(n) = 0.0D0
    d(n) = 0.0D0

    ! Clean up
    DEALLOCATE(h, alpha, mu, z, l)

  END SUBROUTINE splinecof3_fast

end module splinecof3_fast_mod