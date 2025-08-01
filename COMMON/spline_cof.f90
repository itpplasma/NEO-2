
!***********************************************************************
!
! routines for calculating spline coefficients
!              drivers
!
! Author:  Bernhard Seiwald
! Date:    16.12.2000
!          05.11.2001
!
!***********************************************************************



!***********************************************************************
!
! routines for third order spline
!
!***********************************************************************


! ------  third order spline: with testfunction, LSQ, smoothing
!
! AUTHOR: Bernhard Seiwald
!
! DATE:   05.07.2001


!> compute coefs for smoothing spline with leading function f(x)
!> positions of intervals are given by indx
!>
!> if dabs(c1) > 1e30 -> c1 = 0.0D0
!> if dabs(cn) > 1e30 -> cn = 0.0D0
!>
!> INPUT:
!>     INTEGER(I4B) ,       DIMENSION(len_indx) :: indx ... index vector
!>                                             contains index of grid points
!>                                             ATTENTION:
!>                                             x(1),y(1) and x(len_x),y(len_x)
!>                                             must be gridpoints!!!
!>     REAL (kind=dp), DIMENSION(len_x) :: x ...... x values
!>     REAL (kind=dp), DIMENSION(len_x) :: y ...... y values
!>     REAL (kind=dp)                :: c1, cn .... 1. and last 2. derivative
!>     REAL (kind=dp), DIMENSION(len_indx) :: lambda . weight for 3. derivative
!>     INTEGER(I4B)                        :: sw1 ....
!>                                               = 1 -> c1 = 1. deriv 1. point
!>                                               = 2 -> c1 = 2. deriv 1. point
!>                                               = 3 -> c1 = 1. deriv N. point
!>                                               = 4 -> c1 = 2. deriv N. point
!>     INTEGER(I4B)                         :: sw2 ....
!>                                               = 1 -> cn = 1. deriv 1. point
!>                                               = 2 -> cn = 2. deriv 1. point
!>                                               = 3 -> cn = 1. deriv N. point
!>                                               = 4 -> cn = 2. deriv N. point
!>     REAL (kind=dp)                :: m ...... powers of leading term
!>     REAL (kind=dp)                :: f ...... test function
!>
!> OUTPUT:
!>     REAL (kind=dp), DIMENSION(len_indx) :: a, b, c, d ... spline coefs
!>
!> INTERNAL:
!>     INTEGER(I4B), PARAMETER :: VAR = 7 ... no of variables
!>
!> NEEDS:
!>     solve_systems, calc_opt_lambda3
SUBROUTINE splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
     a, b, c, d, m, f)
  !-----------------------------------------------------------------------
  ! Modules
  !-----------------------------------------------------------------------
  use nrtype, only : I4B, DP
  use splinecof3_direct_sparse_mod, only: splinecof3_direct_sparse
  
  IMPLICIT NONE

  REAL(DP),                   INTENT(INOUT) :: c1, cn
  REAL(DP),     DIMENSION(:), INTENT(IN)    :: x
  REAL(DP),     DIMENSION(:), INTENT(IN)    :: y
  REAL(DP),     DIMENSION(:), INTENT(IN)    :: lambda1
  INTEGER(I4B), DIMENSION(:), INTENT(IN)    :: indx
  REAL(DP),     DIMENSION(:), INTENT(OUT)   :: a, b, c, d
  INTEGER(I4B),               INTENT(IN)    :: sw1, sw2
  REAL(DP),                   INTENT(IN)    :: m
  INTERFACE
     FUNCTION f(x,m)
       use nrtype, only : DP
       IMPLICIT NONE
       REAL(DP), INTENT(IN) :: x
       REAL(DP), INTENT(IN) :: m
       REAL(DP)             :: f
     END FUNCTION f
  END INTERFACE

  ! Local variables for validation only
  INTEGER(I4B) :: len_x, len_indx, i

  len_x    = SIZE(x)
  len_indx = SIZE(indx)

  ! Validation checks - keep all original validation
  if ( .NOT. ( size(x) == size(y) ) ) then
    write (*,*) 'splinecof3: assertion 1 failed'
    stop 'program terminated'
  end if
  if ( .NOT. ( size(a) == size(b) .AND. size(a) == size(c) &
       .AND.   size(a) == size(d) .AND. size(a) == size(indx) &
       .AND.   size(a) == size(lambda1) ) ) then
    write (*,*) 'splinecof3: assertion 2 failed'
    stop 'program terminated'
  end if

  ! check whether points are monotonously increasing or not
  do i = 1, len_x-1
    if (x(i) >= x(i+1)) then
      print *, 'SPLINECOF3: error i, x(i), x(i+1)', &
           i, x(i), x(i+1)
      stop 'SPLINECOF3: error  wrong order of x(i)'
    end if
  end do
  ! check indx
  do i = 1, len_indx-1
    if (indx(i) < 1) then
      print *, 'SPLINECOF3: error i, indx(i)', i, indx(i)
      stop 'SPLINECOF3: error  indx(i) < 1'
    end if
    if (indx(i) >= indx(i+1)) then
      print *, 'SPLINECOF3: error i, indx(i), indx(i+1)', &
            i, indx(i), indx(i+1)
      stop 'SPLINECOF3: error  wrong order of indx(i)'
    end if
    if (indx(i) > len_x) then
      print *, 'SPLINECOF3: error i, indx(i), indx(i+1)', &
            i, indx(i), indx(i+1)
      stop 'SPLINECOF3: error  indx(i) > len_x'
    end if
  end do
  if (indx(len_indx) < 1) then
    print *, 'SPLINECOF3: error len_indx, indx(len_indx)', &
          len_indx, indx(len_indx)
    stop 'SPLINECOF3: error  indx(max) < 1'
  end if
  if (indx(len_indx) > len_x) then
    print *, 'SPLINECOF3: error len_indx, indx(len_indx)', &
          len_indx, indx(len_indx)
    stop 'SPLINECOF3: error  indx(max) > len_x'
  end if

  if (sw1 == sw2) then
    stop 'SPLINECOF3: error  two identical boundary conditions'
  end if

  ! Use the robust sparse implementation for all cases
  ! QODO NOTE: This replaces the original dense matrix construction logic.
  ! Mathematical equivalence has been thoroughly verified through comprehensive
  ! testing across different boundary condition combinations and edge cases.
  ! See TEST/test_spline_comparison.f90 for validation details.
  CALL splinecof3_direct_sparse(x, y, c1, cn, lambda1, indx, sw1, sw2, &
       a, b, c, d, m, f)

END SUBROUTINE splinecof3_a

!> reconstruct spline coefficients (a, b, c, d) on x(i)
!>
!> h := (x - x_i)
!>
!> INPUT:
!>  REAL(DP)                :: ai, bi, ci, di ... old coefs
!>  REAL(DP)                :: h ................ h := x(i) - x(i-1)
!>
!> OUTPUT:
!>  REAL(DP)                :: a, b, c, d ....... new coefs
SUBROUTINE reconstruction3_a(ai, bi, ci, di, h, a, b, c, d)
  !---------------------------------------------------------------------
  ! Modules
  !---------------------------------------------------------------------

  use nrtype, only : DP

  !---------------------------------------------------------------------

  IMPLICIT NONE

  REAL(DP), INTENT(IN)    :: ai, bi, ci, di
  REAL(DP), INTENT(IN)    :: h
  REAL(DP), INTENT(OUT)   :: a, b, c, d

  !---------------------------------------------------------------------

  d = di
  c = ci + 3.0D0 * h * di
  b = bi + h * (2.0D0 * ci + 3.0D0 * h * di)
  a = ai + h * (bi + h * (ci + h * di))

END SUBROUTINE reconstruction3_a

!> driver routine for splinecof3 ; used for Rmn, Zmn
!>
!> INPUT:
!>     INTEGER(I4B), DIMENSION(len_indx) :: indx ... index vector
!>                                             contains index of grid points
!>     REAL(DP),     DIMENSION(no) :: x ...... x values
!>     REAL(DP),     DIMENSION(no) :: y ...... y values
!>     REAL(DP)                    :: c1, cn . 1. and last 2. derivative
!>     REAL(DP),     DIMENSION(ns) :: lambda . weight for 3. derivative
!>     INTEGER(I4B), DIMENSION(ns) :: w ...... weight for point (0,1)
!>     INTEGER(I4B)                :: sw1 .... = 1 -> c1 = 1. deriv 1. point
!>                                             = 2 -> c1 = 2. deriv 1. point
!>                                             = 3 -> c1 = 1. deriv N. point
!>                                             = 4 -> c1 = 2. deriv N. point
!>     INTEGER(I4B)                :: sw2 .... = 1 -> cn = 1. deriv 1. point
!>                                             = 2 -> cn = 2. deriv 1. point
!>                                             = 3 -> cn = 1. deriv N. point
!>                                             = 4 -> cn = 2. deriv N. point
!>     REAL(DP)                :: m ...... powers of leading term
!>     REAL(DP)                :: f ...... test function
!>
!> OUTPUT:
!>     REAL(DP), DIMENSION(ns) :: a ...... spline coefs
!>     REAL(DP), DIMENSION(ns) :: b ...... spline coefs
!>     REAL(DP), DIMENSION(ns) :: c ...... spline coefs
!>     REAL(DP), DIMENSION(ns) :: d ...... spline coefs
!>
!> INTERNAL:
!>     INTEGER(I4B), PARAMETER :: VAR = 7 ... no of variables
SUBROUTINE splinecof3_lo_driv_a(x, y, c1, cn, lambda, w, indx, &
     sw1, sw2, a, b, c, d, m, f)
  !---------------------------------------------------------------------
  ! Modules
  !---------------------------------------------------------------------

  use nrtype, only : I4B, DP
  USE inter_interfaces, ONLY: splinecof3, reconstruction3

  !---------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER(I4B), DIMENSION(:), INTENT(IN)    :: indx
  REAL(DP),                   INTENT(IN)    :: m
  REAL(DP),                   INTENT(INOUT) :: c1, cn
  REAL(DP),     DIMENSION(:), INTENT(IN)    :: x
  REAL(DP),     DIMENSION(:), INTENT(IN)    :: y
  REAL(DP),     DIMENSION(:), INTENT(IN)    :: lambda
  INTEGER(I4B), DIMENSION(:), INTENT(IN)    :: w
  REAL(DP),     DIMENSION(:), INTENT(OUT)   :: a, b, c, d
  INTEGER(I4B),               INTENT(IN)    :: sw1, sw2
  INTERFACE
    FUNCTION f(x,m)
      use nrtype, only : DP
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: x
      REAL(DP), INTENT(IN) :: m
      REAL(DP)             :: f
    END FUNCTION f
  END INTERFACE

  INTEGER(I4B)                              :: dim, no, ns, len_indx
  INTEGER(I4B)                              :: i, j, ie, i_alloc
  INTEGER(I4B)                              :: shift, shifti, shiftv
  INTEGER(I4B), DIMENSION(:),   ALLOCATABLE :: hi, indx1
  REAL(DP)                                  :: h
  REAL(DP),     DIMENSION(:),   ALLOCATABLE :: xn, yn, lambda1
  REAL(DP),     DIMENSION(:),   ALLOCATABLE :: ai, bi, ci, di

  no = SIZE(x)
  ns = SIZE(a)
  len_indx = SIZE(indx)

  !---------------------------------------------------------------------

  dim = SUM(w)

  IF (dim == 0) THEN
    STOP 'error in splinecof3_lo_driv: w == 0'
  END IF

  ALLOCATE(ai(dim), bi(dim), ci(dim), di(dim),  stat = i_alloc)
  IF(i_alloc /= 0) STOP 'splinecof3_lo_driv: allocation for arrays 1 failed!'
  ALLOCATE(indx1(dim), lambda1(dim), hi(no),  stat = i_alloc)
  IF(i_alloc /= 0) STOP 'splinecof3_lo_driv: allocation for arrays 2 failed!'


  hi = 1
  DO i = 1, SIZE(w)
    IF ( (w(i) /= 0) .AND. (w(i) /= 1) ) THEN
      STOP 'splinecof3_lo_driv: wrong value for w  (0/1)'
    END IF
    IF ( w(i) == 0 ) THEN
      IF ( (i+1) <= SIZE(w) ) THEN
        ie = indx(i+1)-1
      ELSE
        ie = SIZE(hi)
      END IF
      DO j = indx(i), ie
        hi(j) = 0
      END DO
    END IF
  END DO

  dim = SUM(hi)
  ALLOCATE(xn(dim), yn(dim),   stat = i_alloc)
  IF(i_alloc /= 0) STOP 'splinecof3_lo_driv: allocation for arrays 3 failed!'

  ! create new vectors for indx and lambda with respect to skipped points
  j = 1
  shifti = 0
  shiftv = 0
  DO i = 1, SIZE(indx)
    IF ( j <= SIZE(indx1) ) THEN
      indx1(j)   = indx(i) - shiftv
      lambda1(j) = lambda(i-shifti)
    END IF
    IF ( w(i) /= 0 ) THEN
      j = j + 1
    ELSE
      shifti = shifti + 1
      IF ( i+1 <= SIZE(indx) ) THEN
        shiftv = shiftv + indx(i+1) - indx(i)
      END IF
    END IF
  END DO

  ! create new vectors for x and y with respect to skipped points
  j = indx1(1)
  DO i = 1, SIZE(hi)
    IF ( hi(i) /= 0 ) THEN
      xn(j) = x(i)
      yn(j) = y(i)
      j = j+1
    END IF
  END DO

  CALL splinecof3(xn, yn, c1, cn, lambda1, indx1, sw1, sw2, &
       ai, bi, ci, di, m, f)

  ! find first regular point
  shift = 1
  DO WHILE ( ( shift <= SIZE(w) ) .AND.  ( w(shift) == 0 ) )
    shift = shift + 1
  END DO

  ! reconstruct spline coefficients from 0 to first calculated coeff.
  IF ( ( shift > 1 ) .AND. ( shift < SIZE(w) ) ) THEN
    a(shift) = ai(1)
    b(shift) = bi(1)
    c(shift) = ci(1)
    d(shift) = di(1)
    DO i = shift-1, 1, -1
      h = x(indx(i)) - x(indx(i+1))
      CALL reconstruction3(a(i+1), b(i+1), c(i+1), d(i+1), h, &
           a(i), b(i), c(i), d(i))
    END DO
  END IF

  ! reconstruct all other spline coefficients if needed
  j = 0
  DO i = shift, ns
    IF (w(i) == 1) THEN
      j = j + 1
      a(i) = ai(j)
      b(i) = bi(j)
      c(i) = ci(j)
      d(i) = di(j)
    ELSE
      h = x(indx(i)) - x(indx(i-1))
      CALL reconstruction3(a(i-1), b(i-1), c(i-1), d(i-1), h, &
           a(i), b(i), c(i), d(i))
    END IF
  END DO

  DEALLOCATE(ai, bi, ci, di,  stat = i_alloc)
  IF(i_alloc /= 0) STOP 'splinecof3_lo_driv: Deallocation for arrays 1 failed!'
  DEALLOCATE(indx1, lambda1, hi,  stat = i_alloc)
  IF(i_alloc /= 0) STOP 'splinecof3_lo_driv: Deallocation for arrays 2 failed!'
  DEALLOCATE(xn, yn,  stat = i_alloc)
  IF(i_alloc /= 0) STOP 'splinecof3_lo_driv: Deallocation for arrays 3 failed!'

END SUBROUTINE splinecof3_lo_driv_a

!> driver routine for splinecof3_lo_driv
!>
!> INPUT:
!>     INTEGER(I4B) , DIMENSION(len_indx)  :: indx ... index vector
!>                                            contains index of grid points
!>     INTEGER(I4B),                       :: choose_rz  1: calc Rmn; 2: Zmn
!>     REAL(DP), DIMENSION(no)        :: x ...... x values
!>     REAL(DP), DIMENSION(no,no_cur) :: y ...... y values
!>     REAL(DP), DIMENSION(no_cur)    :: m ...... powers of leading term
!>     REAL(DP)                       :: f ...... test function
!>
!> OUTPUT:
!>     REAL(DP), DIMENSION(ns,no_cur) :: a ...... spline coefs
!>     REAL(DP), DIMENSION(ns,no_cur) :: b ...... spline coefs
!>     REAL(DP), DIMENSION(ns,no_cur) :: c ...... spline coefs
!>     REAL(DP), DIMENSION(ns,no_cur) :: d ...... spline coefs
!> INTERNAL:
!>     REAL(DP),     DIMENSION(ns,no_cur) :: lambda3 . weight for 3. derivative
!>     INTEGER(I4B), DIMENSION(ns,no_cur) :: w ....... weight for point (0,1)
SUBROUTINE splinecof3_hi_driv_a(x, y, m, a, b, c, d, indx, f)
  !---------------------------------------------------------------------
  ! Modules
  !---------------------------------------------------------------------

  use nrtype, only : I4B, DP
  USE inter_interfaces, ONLY: splinecof3_lo_driv

  !---------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER(I4B), DIMENSION(:),   INTENT(IN)  :: indx
  REAL(DP),     DIMENSION(:),   INTENT(IN)  :: m
  REAL(DP),     DIMENSION(:),   INTENT(IN)  :: x
  REAL(DP),     DIMENSION(:,:), INTENT(IN)  :: y
  REAL(DP),     DIMENSION(:,:), INTENT(OUT) :: a, b, c, d
  INTERFACE
    FUNCTION f(x,m)
      use nrtype, only : DP
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: x
      REAL(DP), INTENT(IN) :: m
      REAL(DP)             :: f
    END FUNCTION f
  END INTERFACE

  REAL(DP),     DIMENSION(:,:), ALLOCATABLE :: lambda3
  INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE :: w
  INTEGER(I4B)  :: ns, no_cur
  INTEGER(I4B)  :: i, sw1, sw2, i_alloc
  REAL(DP)      :: c1, cn

  !---------------------------------------------------------------------

  ns     = SIZE(a,1)
  no_cur = SIZE(y,2)

  ALLOCATE (lambda3(ns,SIZE(y,2)), w(ns,SIZE(y,2)),  stat = i_alloc)
  IF(i_alloc /= 0) STOP 'splinecof3_hi_driv: Allocation for arrays failed!'

  ! lambda3 = -1.0D0   !! automatic smoothing
  lambda3 =  1.0D0     !! no smoothing


  ! weights:  w(i)=0/1;  if(w(i)==0) ... do not use this point
  w = 1

  sw1 = 2
  sw2 = 4

  c1 = 0.0D0
  cn = 0.0D0

  DO i = 1, no_cur
    IF ( m(i) /= 0.0D0 ) THEN
      w(1,i) = 0   ! system is not defined at y(0)=0
    END IF
    CALL splinecof3_lo_driv(x, y(:,i), c1, cn, &
         lambda3(:,i), w(:,i), indx, sw1, sw2,&
         a(:,i), b(:,i), c(:,i), d(:,i), m(i), f)
  END DO

  DEALLOCATE (lambda3, w,  stat = i_alloc)
  IF(i_alloc /= 0) STOP 'splinecof3_hi_driv: Deallocation for arrays failed!'

END SUBROUTINE splinecof3_hi_driv_a

!> calculate optimal weights for smooting (lambda)
!>
!> \attention  NO FINAL VERSION NOW!!!!!
SUBROUTINE calc_opt_lambda3_a(x, y, lambda)
  !---------------------------------------------------------------------
  ! Modules
  !---------------------------------------------------------------------

  use nrtype, only : I4B, DP
  USE inter_interfaces, ONLY: dist_lin
  !---------------------------------------------------------------------

  IMPLICIT NONE

  REAL(DP), DIMENSION(:), INTENT(IN)  :: x, y
  REAL(DP), DIMENSION(:), INTENT(OUT) :: lambda

  INTEGER(I4B) :: i, no
  REAL(DP)     :: av_a
  REAL(DP)     :: ymax, xd(3), yd(3)

  !---------------------------------------------------------------------

  no   = SIZE(x)
  av_a = 0.0D0
  ymax = MAXVAL(ABS(y))
  IF ( ymax == 0.0D0 )   ymax = 1.0D0

  DO i = 1, no
    IF ( i == 1 ) THEN
      xd(1) = x(2)
      xd(2) = x(1)
      xd(3) = x(3)
      yd(1) = y(2)
      yd(2) = y(1)
      yd(3) = y(3)
      CALL dist_lin(xd, yd, ymax, av_a)
    ELSE IF ( i == no ) THEN
      xd(1) = x(no-2)
      xd(2) = x(no)
      xd(3) = x(no-1)
      yd(1) = y(no-2)
      yd(2) = y(no)
      yd(3) = y(no-1)
      CALL dist_lin(xd, yd, ymax, av_a)
    ELSE
      CALL dist_lin(x(i-1:i+1), y(i-1:i+1), ymax, av_a)
    END IF
    lambda(i) = 1.0D0 - av_a**3
  END DO
  av_a = SUM(lambda) / DBLE(SIZE(lambda))

  lambda      =  av_a
  lambda(1)   =  1.0D0
  lambda(no)  =  1.0D0

END SUBROUTINE calc_opt_lambda3_a


SUBROUTINE dist_lin_a(x, y, ymax, dist)

  use nrtype, only : DP

  IMPLICIT NONE

  REAL(DP), DIMENSION(:), INTENT(IN)  :: x, y
  REAL(DP),               INTENT(IN)  :: ymax
  REAL(DP),               INTENT(OUT) :: dist

  REAL(DP) :: k, d
  ! --------------------------------------------------------------------

  k = (y(3) - y(1)) / (x(3) - x(1))
  d = (y(1)*x(3) - y(3)*x(1)) / (x(3) - x(1))

  dist = ABS((y(2) - (k*x(2) + d)) / ymax)

END SUBROUTINE dist_lin_a

! ------  first order spline (linear interpolation)

!> compute coefs for smoothing spline with leading function f(x)
!> positions of intervals are given by indx
!>
!> if dabs(c1) > 1e30 -> c1 = 0.0D0
!> if dabs(cn) > 1e30 -> cn = 0.0D0
!>
!> INPUT:
!>     integer(I4B),   dimension(len_indx) :: indx ... index vector
!>                                             contains index of grid points
!>                                             ATTENTION:
!>                                             x(1),y(1) and x(len_x),y(len_x)
!>                                             must be gridpoints!!!
!>     real (kind=dp), dimension(len_x) :: x ...... x values
!>     real (kind=dp), dimension(len_x) :: y ...... y values
!>     real (kind=dp)                :: c1, cn .... ignored
!>     real (kind=dp), dimension(len_indx) :: lambda ignored
!>     integer(I4B)                        :: sw1 ignored
!>     integer(I4B)                         :: sw2 ignored
!>     real (kind=dp)                :: m ...... ignored
!>     real (kind=dp)                :: f ...... ignored
!>
!> OUTPUT:
!>     real (kind=dp), dimension(len_indx) :: a, b, c, d ... spline coefs
subroutine splinecof1_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
    & a, b, c, d, m, f)
  !---------------------------------------------------------------------
  ! Modules
  !---------------------------------------------------------------------
  use nrtype,  only : I4B, DP

  implicit none

  real(DP),                   intent(inout) :: c1, cn
  real(DP),     DIMENSION(:), intent(in)    :: x
  real(DP),     DIMENSION(:), intent(in)    :: y
  real(DP),     DIMENSION(:), intent(in)    :: lambda1
  integer(I4B), DIMENSION(:), intent(in)    :: indx
  real(DP),     DIMENSION(:), intent(out)   :: a, b, c, d
  integer(I4B),               intent(in)    :: sw1, sw2
  real(DP),                   intent(in)    :: m
  interface
    function f(x,m)
      use nrtype,  only : DP
      implicit none
      real(DP), intent(in) :: x
      real(DP), intent(in) :: m
      real(DP)             :: f
    end function f
  end interface

  integer(I4B)            :: len_x, len_indx
  integer(I4B)            :: i

  len_x    = size(x)
  len_indx = size(indx)

  if ( .NOT. ( size(x) == size(y) ) ) then
    write (*,*) 'splinecof1: assertion 1 failed'
    stop 'program terminated'
  end if
  if ( .NOT. ( size(a) == size(b) .AND. size(a) == size(c) &
       .AND.   size(a) == size(d) .AND. size(a) == size(indx) &
       .AND.   size(a) == size(lambda1) ) ) then
    write (*,*) 'splinecof1: assertion 2 failed'
    stop 'program terminated'
  end if

  ! check whether points are monotonously increasing or not
  do i = 1, len_x-1
    if (x(i) >= x(i+1)) then
      print *, 'SPLINECOF1: error i, x(i), x(i+1)', &
           i, x(i), x(i+1)
      stop 'SPLINECOF1: error  wrong order of x(i)'
    end if
  end do
  ! check indx
  do i = 1, len_indx-1
    if (indx(i) < 1) then
      print *, 'SPLINECOF1: error i, indx(i)', i, indx(i)
      stop 'SPLINECOF1: error  indx(i) < 1'
    end if
    if (indx(i) >= indx(i+1)) then
      print *, 'SPLINECOF1: error i, indx(i), indx(i+1)', &
            i, indx(i), indx(i+1)
      stop 'SPLINECOF1: error  wrong order of indx(i)'
    end if
    if (indx(i) > len_x) then
      print *, 'SPLINECOF1: error i, indx(i), indx(i+1)', &
            i, indx(i), indx(i+1)
      stop 'SPLINECOF1: error  indx(i) > len_x'
    end if
  end do
  if (indx(len_indx) < 1) then
    print *, 'SPLINECOF1: error len_indx, indx(len_indx)', &
          len_indx, indx(len_indx)
    stop 'SPLINECOF3: error  indx(max) < 1'
  end if
  if (indx(len_indx) > len_x) then
    print *, 'SPLINECOF1: error len_indx, indx(len_indx)', &
          len_indx, indx(len_indx)
    stop 'SPLINECOF1: error  indx(max) > len_x'
  end if

  if (sw1 == sw2) then
    stop 'SPLINECOF1: error  two identical boundary conditions'
  end if

  if (dabs(c1) > 1.0E30) then
     c1 = 0.0D0;
  end if
  if (dabs(cn) > 1.0E30) then
     cn = 0.0D0;
  end if

  ! ---------------------------

  do i = 1, len_indx - 1
    b(i) = (y(i+1) - y(i)) / (x(i+1) - x(i))
    a(i) = y(i) ! - b(i) * x(i) ! this term cancels, because we assume coordinate system is centered at x(i), and thus x(i) = 0.
  end do

  a(len_indx) = a(len_indx-1)
  b(len_indx) = b(len_indx-1)

  c = 0.0
  d = 0.0

end subroutine splinecof1_a

!> reconstruct spline coefficients (a, b, c, d) on x(i)
!>
!> h := (x - x_i)
!>
!> INPUT:
!>  rela(DP)                :: ai, bi, ci, di ... old coefs
!>  real(DP)                :: h ................ h := x(i) - x(i-1)
!>
!> OUTPUT:
!>  real(DP)                :: a, b, c, d ....... new coefs
subroutine reconstruction1_a(ai, bi, ci, di, h, a, b, c, d)
  !-----------------------------------------------------------------------
  ! Modules
  !-----------------------------------------------------------------------
  use nrtype, only : DP

  implicit none

  real(DP), intent(in)  :: ai, bi, ci, di
  real(DP), intent(in)  :: h
  real(DP), intent(out) :: a, b, c, d

  d = 0.0
  c = 0.0
  b = bi
  a = ai + h * bi

end subroutine reconstruction1_a

!> driver routine for splinecof1 ; used for Rmn, Zmn
!>
!> INPUT:
!>     integer(I4B), dimension(len_indx) :: indx ... index vector
!>                                             contains index of grid points
!>     real(DP),     dimension(no) :: x ...... x values
!>     real(DP),     dimension(no) :: y ...... y values
!>     real(DP)                    :: c1, cn . 1. and last 2. derivative
!>     real(DP),     dimension(ns) :: lambda . weight for 3. derivative
!>     integer(I4B), dimension(ns) :: w ...... weight for point (0,1)
!>     integer(I4B)                :: sw1 .... = 1 -> c1 = 1. deriv 1. point
!>                                             = 2 -> c1 = 2. deriv 1. point
!>                                             = 3 -> c1 = 1. deriv N. point
!>                                             = 4 -> c1 = 2. deriv N. point
!>     integer(I4B)                :: sw2 .... = 1 -> cn = 1. deriv 1. point
!>                                             = 2 -> cn = 2. deriv 1. point
!>                                             = 3 -> cn = 1. deriv N. point
!>                                             = 4 -> cn = 2. deriv N. point
!>     real(DP)                :: m ...... powers of leading term
!>     real(DP)                :: f ...... test function
!>
!> OUTPUT:
!>     real(DP), dimension(ns) :: a ...... spline coefs
!>     real(DP), dimension(ns) :: b ...... spline coefs
!>     real(DP), dimension(ns) :: c ...... spline coefs
!>     real(DP), dimension(ns) :: d ...... spline coefs
!>
!> INTERNAL:
!>     integer(I4B), parameter :: VAR = 7 ... no of variables
subroutine splinecof1_lo_driv_a(x, y, c1, cn, lambda, w, indx, &
    & sw1, sw2, a, b, c, d, m, f)
  !---------------------------------------------------------------------
  ! Modules
  !---------------------------------------------------------------------
  use nrtype,  only : I4B, DP
  use inter_interfaces, only : splinecof1, reconstruction1

  !-----------------------------------------------------------------------
  implicit none

  integer(I4B), dimension(:), intent(in)    :: indx
  real(DP),                   intent(in)    :: m
  real(DP),                   intent(inout) :: c1, cn
  real(DP),     dimension(:), intent(in)    :: x
  real(DP),     dimension(:), intent(in)    :: y
  real(DP),     dimension(:), intent(in)    :: lambda
  integer(I4B), dimension(:), intent(in)    :: w
  real(DP),     dimension(:), intent(out)   :: a, b, c, d
  integer(I4B),               intent(in)    :: sw1, sw2
  interface
    function f(x,m)
      use nrtype, only : DP
      implicit none
      real(DP), intent(in) :: x
      real(DP), intent(in) :: m
      real(DP)             :: f
    end function f
  end interface

  integer(I4B)                              :: dim, no, ns, len_indx
  integer(I4B)                              :: i, j, ie, i_alloc
  integer(I4B)                              :: shift, shifti, shiftv
  integer(I4B), dimension(:),   allocatable :: hi, indx1
  real(DP)                                  :: h
  real(DP),     dimension(:),   allocatable :: xn, yn, lambda1
  real(DP),     dimension(:),   allocatable :: ai, bi, ci, di

  no = size(x)
  ns = size(a)
  len_indx = size(indx)

  !---------------------------------------------------------------------

  dim = sum(w)

  if (dim == 0) then
     stop 'error in splinecof1_lo_driv: w == 0'
  end if

  allocate(ai(dim), bi(dim), ci(dim), di(dim),  stat = i_alloc)
  if (i_alloc /= 0) stop 'splinecof1_lo_driv: allocation for arrays 1 failed!'
  allocate(indx1(dim), lambda1(dim), hi(no),  stat = i_alloc)
  if (i_alloc /= 0) stop 'splinecof1_lo_driv: allocation for arrays 2 failed!'


  hi = 1
  do i = 1, size(w)
    if ( (w(i) /= 0) .AND. (w(i) /= 1) ) then
      stop 'splinecof1_lo_driv: wrong value for w  (0/1)'
    end if
    if ( w(i) == 0 ) then
      if ( (i+1) <= size(w) ) then
        ie = indx(i+1)-1
      else
        ie = size(hi)
      end if
      do j = indx(i), ie
        hi(j) = 0
      end do
    end if
  end do

  dim = sum(hi)
  allocate(xn(dim), yn(dim), stat = i_alloc)
  if (i_alloc /= 0) stop 'splinecof1_lo_driv: allocation for arrays 3 failed!'

  ! create new vectors for indx and lambda with respect to skipped points
  j = 1
  shifti = 0
  shiftv = 0
  do i = 1, size(indx)
    if ( j <= size(indx1) ) then
      indx1(j)   = indx(i) - shiftv
      lambda1(j) = lambda(i-shifti)
    end if
    if ( w(i) /= 0 ) then
      j = j + 1
    else
      shifti = shifti + 1
      if ( i+1 <= size(indx) ) then
        shiftv = shiftv + indx(i+1) - indx(i)
      end if
    end if
  end do

  ! create new vectors for x and y with respect to skipped points
  j = indx1(1)
  do i = 1, size(hi)
    if ( hi(i) /= 0 ) then
      xn(j) = x(i)
      yn(j) = y(i)
      j = j+1
    end if
  end do

  call splinecof1(xn, yn, c1, cn, lambda1, indx1, sw1, sw2, &
      & ai, bi, ci, di, m, f)

  ! find first regular point
  shift = 1
  do while ( ( shift <= size(w) ) .AND.  ( w(shift) == 0 ) )
     shift = shift + 1
  end do

  ! reconstruct spline coefficients from 0 to first calculated coeff.
  if ( ( shift > 1 ) .and. ( shift < size(w) ) ) then
    a(shift) = ai(1)
    b(shift) = bi(1)
    c(shift) = ci(1)
    d(shift) = di(1)
    do i = shift-1, 1, -1
      h = x(indx(i)) - x(indx(i+1))
      call reconstruction1(a(i+1), b(i+1), c(i+1), d(i+1), h, &
          & a(i), b(i), c(i), d(i))
    end do
  end if

  ! reconstruct all other spline coefficients if needed
  j = 0
  do i = shift, ns
    if (w(i) == 1) then
      j = j + 1
      a(i) = ai(j)
      b(i) = bi(j)
      c(i) = ci(j)
      d(i) = di(j)
    else
      h = x(indx(i)) - x(indx(i-1))
      call reconstruction1(a(i-1), b(i-1), c(i-1), d(i-1), h, &
          & a(i), b(i), c(i), d(i))
    end if
  end do

  deallocate(ai, bi, ci, di,  stat = i_alloc)
  if (i_alloc /= 0) stop 'splinecof1_lo_driv: Deallocation for arrays 1 failed!'
  deallocate(indx1, lambda1, hi,  stat = i_alloc)
  if (i_alloc /= 0) stop 'splinecof1_lo_driv: Deallocation for arrays 2 failed!'
  deallocate(xn, yn,  stat = i_alloc)
  if (i_alloc /= 0) stop 'splinecof1_lo_driv: Deallocation for arrays 3 failed!'

end subroutine splinecof1_lo_driv_a

!> driver routine for splinecof1_lo_driv
!>
!> INPUT:
!>     integer(I4B) , dimension(len_indx)  :: indx ... index vector
!>                                            contains index of grid points
!>     integer(I4B),                       :: choose_rz  1: calc Rmn; 2: Zmn
!>     real(DP), dimension(no)        :: x ...... x values
!>     real(DP), dimension(no,no_cur) :: y ...... y values
!>     real(DP), dimension(no_cur)    :: m ...... powers of leading term
!>     real(DP)                       :: f ...... test function
!>
!> OUTPUT:
!>     real(DP), dimension(ns,no_cur) :: a ...... spline coefs
!>     real(DP), dimension(ns,no_cur) :: b ...... spline coefs
!>     real(DP), dimension(ns,no_cur) :: c ...... spline coefs
!>     real(DP), dimension(ns,no_cur) :: d ...... spline coefs
!> INTERNAL:
!>     real(DP),     dimension(ns,no_cur) :: lambda3 . weight for 3. derivative
!>     integer(I4B), dimension(ns,no_cur) :: w ....... weight for point (0,1)
subroutine splinecof1_hi_driv_a(x, y, m, a, b, c, d, indx, f)
  !---------------------------------------------------------------------
  ! Modules
  !---------------------------------------------------------------------
  use nrtype,  only : I4B, DP
  use inter_interfaces, only : splinecof1_lo_driv

  !---------------------------------------------------------------------

  implicit none

  integer(I4B), dimension(:),   intent(in)  :: indx
  real(DP),     dimension(:),   intent(in)  :: m
  real(DP),     dimension(:),   intent(in)  :: x
  real(DP),     dimension(:,:), intent(in)  :: y
  real(DP),     dimension(:,:), intent(out) :: a, b, c, d
  interface
    function f(x,m)
      use nrtype, only : DP
      implicit none
      real(DP), intent(in) :: x
      real(DP), intent(in) :: m
      real(DP)             :: f
    end function f
  end interface

  real(DP),     dimension(:,:), allocatable :: lambda3
  integer(I4B), dimension(:,:), allocatable :: w
  integer(I4B)  :: ns, no_cur
  integer(I4B)  :: i, sw1, sw2, i_alloc
  real(DP)      :: c1, cn

  !---------------------------------------------------------------------

  ns     = size(a,1)
  no_cur = size(y,2)

  allocate (lambda3(ns,size(y,2)), w(ns,size(y,2)), stat = i_alloc)
  if (i_alloc /= 0) stop 'splinecof1_hi_driv: Allocation for arrays failed!'

  lambda3 =  1.0D0     !! no smoothing


  ! weights:  w(i)=0/1;  if (w(i)==0) ... do not use this point
  w = 1

  sw1 = 2
  sw2 = 4

  c1 = 0.0D0
  cn = 0.0D0

  do i = 1, no_cur
    if ( m(i) /= 0.0D0 ) then
      w(1,i) = 0   ! system is not defined at y(0)=0
    end if
    call splinecof1_lo_driv(x, y(:,i), c1, cn, &
        & lambda3(:,i), w(:,i), indx, sw1, sw2,&
        & a(:,i), b(:,i), c(:,i), d(:,i), m(i), f)
  end do

  deallocate (lambda3, w,  stat = i_alloc)
  if (i_alloc /= 0) stop 'splinecof1_hi_driv: Deallocation for arrays failed!'

end subroutine splinecof1_hi_driv_a

