
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


SUBROUTINE splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
     a, b, c, d, m, f)

!-----------------------------------------------------------------------
! 
! compute coefs for smoothing spline with leading function f(x)
! positions of intervals are given by indx 
!
! if dabs(c1) > 1e30 -> c1 = 0.0D0
! if dabs(cn) > 1e30 -> cn = 0.0D0
!
! INPUT:
!     INTEGER(I4B) ,       DIMENSION(len_indx) :: indx ... index vector
!                                             contains index of grid points
!                                             ATTENTION: 
!                                             x(1),y(1) and x(len_x),y(len_x)
!                                             must be gridpoints!!!
!     REAL (kind=dp), DIMENSION(len_x) :: x ...... x values
!     REAL (kind=dp), DIMENSION(len_x) :: y ...... y values
!     REAL (kind=dp)                :: c1, cn .... 1. and last 2. derivative
!     REAL (kind=dp), DIMENSION(len_indx) :: lambda . weight for 3. derivative
!     INTEGER(I4B)                        :: sw1 .... 
!                                               = 1 -> c1 = 1. deriv 1. point
!                                               = 2 -> c1 = 2. deriv 1. point
!                                               = 3 -> c1 = 1. deriv N. point
!                                               = 4 -> c1 = 2. deriv N. point
!     INTEGER(I4B)                         :: sw2 .... 
!                                               = 1 -> cn = 1. deriv 1. point
!                                               = 2 -> cn = 2. deriv 1. point
!                                               = 3 -> cn = 1. deriv N. point
!                                               = 4 -> cn = 2. deriv N. point
!     REAL (kind=dp)                :: m ...... powers of leading term
!     REAL (kind=dp)                :: f ...... test function
!
! OUTPUT:
!     REAL (kind=dp), DIMENSION(len_indx) :: a, b, c, d ... spline coefs
!
! INTERNAL:
!     INTEGER(I4B), PARAMETER :: VAR = 7 ... no of variables
!
! NEEDS:
!     solve_systems, calc_opt_lambda3
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
! Modules
!-----------------------------------------------------------------------

  USE inter_precision,  ONLY: I4B, DP
  USE inter_interfaces, ONLY: calc_opt_lambda3
  !! Modifications by Andreas F. Martitsch (06.08.2014)
  !Replace standard solver from Lapack with sparse solver
  !(Bad performance for more than 1000 flux surfaces ~ (3*nsurf)^2)
  USE solve_systems
  USE sparse_mod, ONLY : sparse_solve
  !! End Modifications by Andreas F. Martitsch (06.08.2014)

!-----------------------------------------------------------------------

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
       USE inter_precision, ONLY: DP
       IMPLICIT NONE
       REAL(DP), INTENT(IN) :: x
       REAL(DP), INTENT(IN) :: m
       REAL(DP)             :: f
     END FUNCTION f
  END INTERFACE

  INTEGER(I4B), PARAMETER :: VAR = 7
  INTEGER(I4B)            :: dim
  INTEGER(I4B)            :: i_alloc, info
  INTEGER(I4B)            :: len_x, len_indx
  INTEGER(I4B)            :: i, j, l, ii, ie
  INTEGER(I4B)            :: mu1, mu2, nu1, nu2
  INTEGER(I4B)            :: sig1, sig2, rho1, rho2
  INTEGER(I4B), DIMENSION(:),   ALLOCATABLE :: indx_lu
  REAL(DP)                :: h, h_j, x_h, help_i, help_inh
  REAL(DP)                :: help_a, help_b, help_c, help_d
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: MA
  REAL(DP), DIMENSION(:),   ALLOCATABLE :: inh, simqa, lambda, omega

  len_x    = SIZE(x)
  len_indx = SIZE(indx)
  dim = VAR * len_indx - 2

  ALLOCATE(MA(dim, dim),  stat = i_alloc)
  IF(i_alloc /= 0) STOP 'splinecof3: Allocation for arrays 1 failed!'
  ALLOCATE(inh(dim), indx_lu(dim),  stat = i_alloc)
  IF(i_alloc /= 0) STOP 'splinecof3: Allocation for arrays 2 failed!'
  ALLOCATE(simqa(dim*dim),  stat = i_alloc)
  IF(i_alloc /= 0) STOP 'splinecof3: Allocation for arrays 3 failed!'
  ALLOCATE(lambda(SIZE(lambda1)),  stat = i_alloc)
  IF(i_alloc /= 0) STOP 'splinecof3: Allocation for lambda failed!'
  ALLOCATE(omega(SIZE(lambda1)),  stat = i_alloc)
  IF(i_alloc /= 0) STOP 'splinecof3: Allocation for omega failed!'
!-----------------------------------------------------------------------

  IF ( .NOT. ( SIZE(x) == SIZE(y) ) ) THEN
     WRITE (*,*) 'splinecof3: assertion 1 failed'
     STOP 'program terminated'
  END IF
  IF ( .NOT. ( SIZE(a) == SIZE(b) .AND. SIZE(a) == SIZE(c) &
       .AND.   SIZE(a) == SIZE(d) .AND. SIZE(a) == SIZE(indx) &
       .AND.   SIZE(a) == SIZE(lambda1) ) ) THEN
     WRITE (*,*) 'splinecof3: assertion 2 failed'
     STOP 'program terminated'
  END IF


  IF (DABS(c1) > 1.0E30) THEN
     c1 = 0.0D0;
  END IF
  IF (DABS(cn) > 1.0E30) THEN
     cn = 0.0D0;
  END IF

! setting all to zero
  MA(:,:) = 0.0D0
  inh(:)  = 0.0D0

! check whether points are monotonously increasing or not
  DO i = 1, len_x-1
     IF (x(i) >= x(i+1)) THEN
        PRINT *, 'SPLINECOF3: error i, x(i), x(i+1)', &
             i, x(i), x(i+1)
        STOP 'SPLINECOF3: error  wrong order of x(i)'
     END IF
  END DO
! check indx
  DO i = 1, len_indx-1
     IF (indx(i) < 1) THEN
        PRINT *, 'SPLINECOF3: error i, indx(i)', i, indx(i)
        STOP 'SPLINECOF3: error  indx(i) < 1'
     END IF
     IF (indx(i) >= indx(i+1)) THEN
        PRINT *, 'SPLINECOF3: error i, indx(i), indx(i+1)', &
             i, indx(i), indx(i+1)
        STOP 'SPLINECOF3: error  wrong order of indx(i)'
     END IF
     IF (indx(i) > len_x) THEN
        PRINT *, 'SPLINECOF3: error i, indx(i), indx(i+1)', &
             i, indx(i), indx(i+1)
        STOP 'SPLINECOF3: error  indx(i) > len_x'
     END IF
  END DO
  IF (indx(len_indx) < 1) THEN
     PRINT *, 'SPLINECOF3: error len_indx, indx(len_indx)', &
          len_indx, indx(len_indx)
     STOP 'SPLINECOF3: error  indx(max) < 1'
  END IF
  IF (indx(len_indx) > len_x) THEN
     PRINT *, 'SPLINECOF3: error len_indx, indx(len_indx)', &
          len_indx, indx(len_indx)
     STOP 'SPLINECOF3: error  indx(max) > len_x'
  END IF

! calculate optimal weights for smooting (lambda)
  IF ( MAXVAL(lambda1) < 0.0D0 ) THEN
     CALL calc_opt_lambda3(x, y, omega)
  ELSE
     omega  = lambda1
  END IF
  lambda = 1.0D0 - omega

  IF (sw1 == sw2) THEN
     STOP 'SPLINECOF3: error  two identical boundary conditions'
  END IF

  IF (sw1 == 1) THEN
     mu1  = 1
     nu1  = 0
     sig1 = 0
     rho1 = 0
  ELSE IF (sw1 == 2) THEN
     mu1  = 0
     nu1  = 1
     sig1 = 0
     rho1 = 0
  ELSE IF (sw1 == 3) THEN
     mu1  = 0
     nu1  = 0
     sig1 = 1
     rho1 = 0
  ELSE IF (sw1 == 4) THEN
     mu1  = 0
     nu1  = 0
     sig1 = 0
     rho1 = 1
  ELSE
     STOP 'SPLINECOF3: error  in using boundary condition 1'
  END IF

  IF (sw2 == 1) THEN
     mu2  = 1
     nu2  = 0
     sig2 = 0
     rho2 = 0
  ELSE IF (sw2 == 2) THEN
     mu2  = 0
     nu2  = 1
     sig2 = 0
     rho2 = 0
  ELSE IF (sw2 == 3) THEN
     mu2  = 0
     nu2  = 0
     sig2 = 1
     rho2 = 0
  ELSE IF (sw2 == 4) THEN
     mu2  = 0
     nu2  = 0
     sig2 = 0
     rho2 = 1
  ELSE
     STOP 'SPLINECOF3: error  in using boundary condition 2'
  END IF


! coefs for first point
  i  = 0
  j  = 1
  ii = indx((j-1)/VAR+1)
  ie = indx((j-1)/VAR+2) - 1
  h  = x(indx((j-1)/VAR+2)) - x(ii)

! boundary condition 1
  i = i + 1
  MA(i, 2) = DBLE(mu1)
  MA(i, 3) = DBLE(nu1)
  MA(i, (len_indx-1)*VAR + 2) = DBLE(sig1)
  MA(i, (len_indx-1)*VAR + 3) = DBLE(rho1)
  inh(i) = c1

! A_i
  i = i + 1
  MA(i, j+0  +0) =  1.0D0
  MA(i, j+0  +1) =  h
  MA(i, j+0  +2) =  h * h
  MA(i, j+0  +3) =  h * h * h
  MA(i, j+VAR+0) = -1.0D0
! B_i
  i = i + 1
  MA(i, j+0  +1) =  1.0D0
  MA(i, j+0  +2) =  2.0D0 * h
  MA(i, j+0  +3) =  3.0D0 * h * h
  MA(i, j+VAR+1) = -1.0D0
! C_i
  i = i + 1
  MA(i, j+0  +2) =  1.0D0
  MA(i, j+0  +3) =  3.0D0 * h
  MA(i, j+VAR+2) = -1.0D0
! delta a_i
  i = i + 1
  help_a = 0.0D0
  help_b = 0.0D0
  help_c = 0.0D0
  help_d = 0.0D0
  help_i = 0.0D0
  DO l = ii, ie
     h_j = x(l) - x(ii)
     x_h    = f(x(l),m) * f(x(l),m)
     help_a = help_a + x_h
     help_b = help_b + h_j * x_h
     help_c = help_c + h_j * h_j * x_h
     help_d = help_d + h_j * h_j * h_j * x_h
     help_i = help_i + f(x(l),m) * y(l)
  END DO  ! DO l = ii, ie
  MA(i, j+0  +0) =  omega((j-1)/VAR+1) * help_a
  MA(i, j+0  +1) =  omega((j-1)/VAR+1) * help_b
  MA(i, j+0  +2) =  omega((j-1)/VAR+1) * help_c
  MA(i, j+0  +3) =  omega((j-1)/VAR+1) * help_d
  MA(i, j+0  +4) =  1.0D0
  inh(i)         =  omega((j-1)/VAR+1) * help_i
! delta b_i
  i = i + 1
  help_a = 0.0D0
  help_b = 0.0D0
  help_c = 0.0D0
  help_d = 0.0D0
  help_i = 0.0D0
  DO l = ii, ie
     h_j = x(l) - x(ii)
     x_h    = f(x(l),m) * f(x(l),m)
     help_a = help_a + h_j * x_h
     help_b = help_b + h_j * h_j * x_h
     help_c = help_c + h_j * h_j * h_j * x_h
     help_d = help_d + h_j * h_j * h_j * h_j * x_h
     help_i = help_i + h_j * f(x(l),m) * y(l)
  END DO  ! DO l = ii, ie
  MA(i, j+0  +0) =  omega((j-1)/VAR+1) * help_a
  MA(i, j+0  +1) =  omega((j-1)/VAR+1) * help_b
  MA(i, j+0  +2) =  omega((j-1)/VAR+1) * help_c
  MA(i, j+0  +3) =  omega((j-1)/VAR+1) * help_d
  MA(i, j+0  +4) =  h
  MA(i, j+0  +5) =  1.0D0
  MA(i, (len_indx-1)*VAR+4) = DBLE(mu1)
  MA(i, (len_indx-1)*VAR+5) = DBLE(mu2)
  inh(i)         =  omega((j-1)/VAR+1) * help_i
! delta c_i
  i = i + 1
  help_a = 0.0D0
  help_b = 0.0D0
  help_c = 0.0D0
  help_d = 0.0D0
  help_i = 0.0D0
  DO l = ii, ie
     h_j = x(l) - x(ii)
     x_h    = f(x(l),m) * f(x(l),m)
     help_a = help_a + h_j * h_j * x_h
     help_b = help_b + h_j * h_j * h_j * x_h
     help_c = help_c + h_j * h_j * h_j * h_j * x_h
     help_d = help_d + h_j * h_j * h_j * h_j * h_j * x_h
     help_i = help_i + h_j * h_j * f(x(l),m) * y(l)
  END DO  ! DO l = ii, ie
  MA(i, j+0  +0) =  omega((j-1)/VAR+1) * help_a
  MA(i, j+0  +1) =  omega((j-1)/VAR+1) * help_b
  MA(i, j+0  +2) =  omega((j-1)/VAR+1) * help_c
  MA(i, j+0  +3) =  omega((j-1)/VAR+1) * help_d
  MA(i, j+0  +4) =  h * h
  MA(i, j+0  +5) =  2.0D0 * h
  MA(i, j+0  +6) =  1.0D0
  MA(i, (len_indx-1)*VAR+4) = DBLE(nu1)
  MA(i, (len_indx-1)*VAR+5) = DBLE(nu2)
  inh(i)         =  omega((j-1)/VAR+1) * help_i
! delta DELTA d_i
  i = i + 1
  help_a = 0.0D0
  help_b = 0.0D0
  help_c = 0.0D0
  help_d = 0.0D0
  help_i = 0.0D0
  DO l = ii, ie
     h_j = x(l) - x(ii)
     x_h    = f(x(l),m) * f(x(l),m)
     help_a = help_a + h_j * h_j * h_j * x_h
     help_b = help_b + h_j * h_j * h_j * h_j * x_h
     help_c = help_c + h_j * h_j * h_j * h_j * h_j * x_h
     help_d = help_d + h_j * h_j * h_j * h_j * h_j * h_j * x_h
     help_i = help_i + h_j * h_j * h_j * f(x(l),m) * y(l)
  END DO  ! DO l = ii, ie
  MA(i, j+0  +0) =  omega((j-1)/VAR+1) * help_a
  MA(i, j+0  +1) =  omega((j-1)/VAR+1) * help_b
  MA(i, j+0  +2) =  omega((j-1)/VAR+1) * help_c
  MA(i, j+0  +3) =  omega((j-1)/VAR+1) * help_d + lambda((j-1)/VAR+1)
  MA(i, j+0  +4) =  h * h * h
  MA(i, j+0  +5) =  3.0D0 * h * h
  MA(i, j+0  +6) =  3.0D0 * h
  inh(i)         =  omega((j-1)/VAR+1) * help_i

! coefs for point 2 to len_x_points-1
  DO j = VAR+1, VAR*(len_indx-1)-1, VAR
     ii = indx((j-1)/VAR+1)
     ie = indx((j-1)/VAR+2) - 1
     h  = x(indx((j-1)/VAR+2)) - x(ii)
! A_i
     i = i + 1
     MA(i, j+0  +0) =  1.0D0
     MA(i, j+0  +1) =  h
     MA(i, j+0  +2) =  h * h
     MA(i, j+0  +3) =  h * h * h
     MA(i, j+VAR+0) = -1.0D0
! B_i
     i = i + 1
     MA(i, j+0  +1) =  1.0D0
     MA(i, j+0  +2) =  2.0D0 * h
     MA(i, j+0  +3) =  3.0D0 * h * h
     MA(i, j+VAR+1) = -1.0D0
! C_i
     i = i + 1
     MA(i, j+0  +2) =  1.0D0
     MA(i, j+0  +3) =  3.0D0 * h
     MA(i, j+VAR+2) = -1.0D0
! delta a_i
     i = i + 1
     help_a = 0.0D0
     help_b = 0.0D0
     help_c = 0.0D0
     help_d = 0.0D0
     help_i = 0.0D0
     DO l = ii, ie
        h_j = x(l) - x(ii)
        x_h    = f(x(l),m) * f(x(l),m)
        help_a = help_a + x_h
        help_b = help_b + h_j * x_h
        help_c = help_c + h_j * h_j * x_h
        help_d = help_d + h_j * h_j * h_j * x_h
        help_i = help_i + f(x(l),m) * y(l)
     END DO   ! DO l = ii, ie
     MA(i, j+0  +0) =  omega((j-1)/VAR+1) * help_a
     MA(i, j+0  +1) =  omega((j-1)/VAR+1) * help_b
     MA(i, j+0  +2) =  omega((j-1)/VAR+1) * help_c
     MA(i, j+0  +3) =  omega((j-1)/VAR+1) * help_d
     MA(i, j+0  +4) =  1.0D0
     MA(i, j-VAR+4) = -1.0D0
     inh(i)         =  omega((j-1)/VAR+1) * help_i
! delta b_i
     i = i + 1
     help_a = 0.0D0
     help_b = 0.0D0
     help_c = 0.0D0
     help_d = 0.0D0
     help_i = 0.0D0
     DO l = ii, ie
        h_j = x(l) - x(ii)
        x_h    = f(x(l),m) * f(x(l),m)
        help_a = help_a + h_j * x_h
        help_b = help_b + h_j * h_j * x_h
        help_c = help_c + h_j * h_j * h_j * x_h
        help_d = help_d + h_j * h_j * h_j * h_j * x_h
        help_i = help_i + h_j * f(x(l),m) * y(l)
     END DO  ! DO l = ii, ie
     MA(i, j+0  +0) =  omega((j-1)/VAR+1) * help_a
     MA(i, j+0  +1) =  omega((j-1)/VAR+1) * help_b
     MA(i, j+0  +2) =  omega((j-1)/VAR+1) * help_c
     MA(i, j+0  +3) =  omega((j-1)/VAR+1) * help_d
     MA(i, j+0  +4) =  h
     MA(i, j+0  +5) =  1.0D0
     MA(i, j-VAR+5) = -1.0D0
     inh(i)         =  omega((j-1)/VAR+1) * help_i
! delta c_i
     i = i + 1
     help_a = 0.0D0
     help_b = 0.0D0
     help_c = 0.0D0
     help_d = 0.0D0
     help_i = 0.0D0
     DO l = ii, ie
        h_j = x(l) - x(ii)
        x_h    = f(x(l),m) * f(x(l),m)
        help_a = help_a + h_j * h_j * x_h
        help_b = help_b + h_j * h_j * h_j * x_h
        help_c = help_c + h_j * h_j * h_j * h_j * x_h
        help_d = help_d + h_j * h_j * h_j * h_j * h_j * x_h
        help_i = help_i + h_j * h_j * f(x(l),m) * y(l)
     END DO  ! DO l = ii, ie
     MA(i, j+0  +0) =  omega((j-1)/VAR+1) * help_a
     MA(i, j+0  +1) =  omega((j-1)/VAR+1) * help_b
     MA(i, j+0  +2) =  omega((j-1)/VAR+1) * help_c
     MA(i, j+0  +3) =  omega((j-1)/VAR+1) * help_d
     MA(i, j+0  +4) =  h * h
     MA(i, j+0  +5) =  2.0D0 * h
     MA(i, j+0  +6) =  1.0D0
     MA(i, j-VAR+6) = -1.0D0
     inh(i)         =  omega((j-1)/VAR+1) * help_i
! delta DELTA d_i
     i = i + 1
     help_a = 0.0D0
     help_b = 0.0D0
     help_c = 0.0D0
     help_d = 0.0D0
     help_i = 0.0D0
     DO l = ii, ie
        h_j = x(l) - x(ii)
        x_h    = f(x(l),m) * f(x(l),m)
        help_a = help_a + h_j * h_j * h_j * x_h
        help_b = help_b + h_j * h_j * h_j * h_j * x_h
        help_c = help_c + h_j * h_j * h_j * h_j * h_j * x_h
        help_d = help_d + h_j * h_j * h_j * h_j * h_j * h_j * x_h
        help_i = help_i + h_j * h_j * h_j * f(x(l),m) * y(l)
     END DO  ! DO l = ii, ie
     MA(i, j+0  +0) =  omega((j-1)/VAR+1) * help_a
     MA(i, j+0  +1) =  omega((j-1)/VAR+1) * help_b
     MA(i, j+0  +2) =  omega((j-1)/VAR+1) * help_c
     MA(i, j+0  +3) =  omega((j-1)/VAR+1) * help_d + lambda((j-1)/VAR+1)
     MA(i, j+0  +4) =  h * h * h
     MA(i, j+0  +5) =  3.0D0 * h * h
     MA(i, j+0  +6) =  3.0D0 * h
     inh(i)         =  omega((j-1)/VAR+1) * help_i
  END DO  ! DO j = VAR+1, VAR*(len_indx-1)-1, VAR

! last point
! delta a_i
  i = i + 1
  ii = indx((j-1)/VAR+1)
  ie = ii
  help_a   = 0.0D0
  help_inh = 0.0D0
  l = ii
  help_a   = help_a   + f(x(l),m) * f(x(l),m)
  help_inh = help_inh + f(x(l),m) * y(l)
!!$  MA(i, (len_indx-1)*VAR+1) = help_a     ! allowed, if omega != 0
!!$  MA(i, (len_indx-2)*VAR+5) = -1.0D0
!!$  inh(i)                    = help_inh
  MA(i, (len_indx-1)*VAR+1) = omega((j-1)/VAR+1) * help_a
  MA(i, (len_indx-2)*VAR+5) = omega((j-1)/VAR+1) * (-1.0D0)
  inh(i)                    = omega((j-1)/VAR+1) * help_inh
! delta b_i
  i = i + 1
  MA(i, (len_indx-2)*VAR+6) = -1.0D0
  MA(i, (len_indx-1)*VAR+4) =  DBLE(sig1)
  MA(i, (len_indx-1)*VAR+5) =  DBLE(sig2)
! delta c_i
  i = i + 1
  MA(i, (len_indx-2)*VAR+7) = -1.0D0
  MA(i, (len_indx-1)*VAR+4) =  DBLE(rho1)
  MA(i, (len_indx-1)*VAR+5) =  DBLE(rho2)

!!$! boundary condition 1
!!$  i = i + 1
!!$  MA(i, 2) = DBLE(mu1)
!!$  MA(i, 3) = DBLE(nu1)
!!$  MA(i, (len_indx-1)*VAR + 2) = DBLE(sig1)
!!$  MA(i, (len_indx-1)*VAR + 3) = DBLE(rho1)
!!$  inh(i) = c1

! boundary condition 2
  i = i + 1
  MA(i, 2) = DBLE(mu2)
  MA(i, 3) = DBLE(nu2)
  MA(i, (len_indx-1)*VAR + 2) = DBLE(sig2)
  MA(i, (len_indx-1)*VAR + 3) = DBLE(rho2)
  inh(i) = cn

! make column vector simqa out of matrix MA
!!$  k = 1
!!$  DO i=1, dim
!!$     DO j=1,dim
!!$        simqa(k) = MA(j,i)
!!$        k = k + 1
!!$     END DO
!!$  END DO

! --- for checking matrix ---
!!$!  PRINT *, dim
!!$  DO i=1, dim
!!$!     PRINT *, 'i = ',i
!!$     DO j=1,dim
!!$        PRINT *, i, j, MA(i,j)
!!$!        PRINT *, MA(i,:)
!!$     END DO
!!$  END DO
!!$!  PRINT *,'inhom'
!!$!  DO i=1, dim
!!$!    PRINT *, i, inh(i)
!!$!    PRINT *, inh(i)
!!$!  END DO
! ---------------------------

! solve system
  !! Modifications by Andreas F. Martitsch (06.08.2014)
  !Replace standard solver from Lapack with sparse solver
  !(Bad performance for more than 1000 flux surfaces ~ (3*nsurf)^2)
  !PRINT *,"Before - solve system for spline coefficients"
  !CALL solve_eqsys(MA, inh, info)
  !IF (info /= 0) STOP 'splinecof3_a: Singular matrix in solve_eqsys()!'
  CALL sparse_solve(MA, inh)
  !PRINT *,"After - solve system for spline coefficients"
  !! Modifications by Andreas F. Martitsch (06.08.2014)
  
!  PRINT *, ' A B C D'
! take a(), b(), c(), d()
  DO i = 1, len_indx
     a(i) = inh((i-1)*VAR+1) 
     b(i) = inh((i-1)*VAR+2) 
     c(i) = inh((i-1)*VAR+3)
     d(i) = inh((i-1)*VAR+4) 
!     PRINT *, a(i), b(i), c(i), d(i)
  END DO
  

  DEALLOCATE(MA,  stat = i_alloc)
  IF(i_alloc /= 0) STOP 'splinecof3: Deallocation for arrays 1 failed!'
  DEALLOCATE(inh, indx_lu,  stat = i_alloc)
  IF(i_alloc /= 0) STOP 'splinecof3: Deallocation for arrays 2 failed!'
  DEALLOCATE(simqa,  stat = i_alloc)
  IF(i_alloc /= 0) STOP 'splinecof3: Deallocation for arrays 3 failed!'
  DEALLOCATE(lambda,  stat = i_alloc)
  IF(i_alloc /= 0) STOP 'splinecof3: Deallocation for lambda failed!'
  DEALLOCATE(omega,  stat = i_alloc)
  IF(i_alloc /= 0) STOP 'splinecof3: Deallocation for omega failed!'

END SUBROUTINE splinecof3_a




SUBROUTINE reconstruction3_a(ai, bi, ci, di, h, a, b, c, d)

!-----------------------------------------------------------------------
! 
! reconstruct spline coefficients (a, b, c, d) on x(i)
! h := (x - x_i)
!
! INPUT:
!  REAL(DP)                :: ai, bi, ci, di ... old coefs
!  REAL(DP)                :: h ................ h := x(i) - x(i-1)
!
! OUTPUT:
!  REAL(DP)                :: a, b, c, d ....... new coefs
!
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
! Modules
!-----------------------------------------------------------------------

  USE inter_precision, ONLY: DP

!-----------------------------------------------------------------------

  IMPLICIT NONE

  REAL(DP), INTENT(IN)    :: ai, bi, ci, di
  REAL(DP), INTENT(IN)    :: h
  REAL(DP), INTENT(OUT)   :: a, b, c, d 

!-----------------------------------------------------------------------

  d = di
  c = ci + 3.0D0 * h * di
  b = bi + h * (2.0D0 * ci + 3.0D0 * h * di)
  a = ai + h * (bi + h * (ci + h * di))

END SUBROUTINE reconstruction3_a




SUBROUTINE splinecof3_lo_driv_a(x, y, c1, cn, lambda, w, indx, &
     sw1, sw2, a, b, c, d, m, f)

!-----------------------------------------------------------------------
! 
! driver routine for splinecof3 ; used for Rmn, Zmn
!
! INPUT:
!     INTEGER(I4B) ,       DIMENSION(len_indx)  :: indx ... index vector
!                                             contains index of grid points
!     REAL(DP),     DIMENSION(no) :: x ...... x values
!     REAL(DP),     DIMENSION(no) :: y ...... y values
!     REAL(DP)                    :: c1, cn . 1. and last 2. derivative
!     REAL(DP),     DIMENSION(ns) :: lambda . weight for 3. derivative
!     INTEGER(I4B), DIMENSION(ns) :: w ...... weight for point (0,1)
!     INTEGER(I4B)                :: sw1 .... = 1 -> c1 = 1. deriv 1. point
!                                             = 2 -> c1 = 2. deriv 1. point
!                                             = 3 -> c1 = 1. deriv N. point
!                                             = 4 -> c1 = 2. deriv N. point
!     INTEGER(I4B)                :: sw2 .... = 1 -> cn = 1. deriv 1. point
!                                             = 2 -> cn = 2. deriv 1. point
!                                             = 3 -> cn = 1. deriv N. point
!                                             = 4 -> cn = 2. deriv N. point
!     REAL(DP)                :: m ...... powers of leading term
!     REAL(DP)                :: f ...... test function
!
! OUTPUT:
!     REAL(DP), DIMENSION(ns) :: a ...... spline coefs
!     REAL(DP), DIMENSION(ns) :: b ...... spline coefs
!     REAL(DP), DIMENSION(ns) :: c ...... spline coefs
!     REAL(DP), DIMENSION(ns) :: d ...... spline coefs
!
! INTERNAL:
!     INTEGER(I4B), PARAMETER :: VAR = 7 ... no of variables
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
! Modules
!-----------------------------------------------------------------------

  USE inter_precision,  ONLY: I4B, DP
  USE inter_interfaces, ONLY: splinecof3, reconstruction3

!-----------------------------------------------------------------------

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
       USE inter_precision, ONLY: DP
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

!-----------------------------------------------------------------------

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




SUBROUTINE splinecof3_hi_driv_a(x, y, m, a, b, c, d, indx, f)

!-----------------------------------------------------------------------
! 
! driver routine for splinecof3_lo_driv
!
! INPUT:
!     INTEGER(I4B) , DIMENSION(len_indx)  :: indx ... index vector
!                                            contains index of grid points
!     INTEGER(I4B),                       :: choose_rz  1: calc Rmn; 2: Zmn
!     REAL(DP), DIMENSION(no)        :: x ...... x values
!     REAL(DP), DIMENSION(no,no_cur) :: y ...... y values
!     REAL(DP), DIMENSION(no_cur)    :: m ...... powers of leading term
!     REAL(DP)                       :: f ...... test function
!
! OUTPUT:
!     REAL(DP), DIMENSION(ns,no_cur) :: a ...... spline coefs
!     REAL(DP), DIMENSION(ns,no_cur) :: b ...... spline coefs
!     REAL(DP), DIMENSION(ns,no_cur) :: c ...... spline coefs
!     REAL(DP), DIMENSION(ns,no_cur) :: d ...... spline coefs
! INTERNAL:
!     REAL(DP),     DIMENSION(ns,no_cur) :: lambda3 . weight for 3. derivative
!     INTEGER(I4B), DIMENSION(ns,no_cur) :: w ....... weight for point (0,1)
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
! Modules
!-----------------------------------------------------------------------

  USE inter_precision,  ONLY: I4B, DP
  USE inter_interfaces, ONLY: splinecof3_lo_driv

!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER(I4B), DIMENSION(:),   INTENT(IN)  :: indx
  REAL(DP),     DIMENSION(:),   INTENT(IN)  :: m
  REAL(DP),     DIMENSION(:),   INTENT(IN)  :: x
  REAL(DP),     DIMENSION(:,:), INTENT(IN)  :: y
  REAL(DP),     DIMENSION(:,:), INTENT(OUT) :: a, b, c, d
  INTERFACE
     FUNCTION f(x,m)
       USE inter_precision, ONLY: DP
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

!-----------------------------------------------------------------------

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



SUBROUTINE calc_opt_lambda3_a(x, y, lambda)
!   NO FINAL VERSION NOW!!!!!
! calculate optimal weights for smooting (lambda)
!
!-----------------------------------------------------------------------
! Modules
!-----------------------------------------------------------------------

  USE inter_precision,  ONLY: I4B, DP
  USE inter_interfaces, ONLY: dist_lin
!-----------------------------------------------------------------------

  IMPLICIT NONE

  REAL(DP), DIMENSION(:), INTENT(IN)  :: x, y
  REAL(DP), DIMENSION(:), INTENT(OUT) :: lambda

  INTEGER(I4B) :: i, no
  REAL(DP)     :: av_a
  REAL(DP)     :: ymax, xd(3), yd(3)

!-----------------------------------------------------------------------

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
!!$     IF (x(i) < 0.2) THEN
!!$     lambda(i) = 1.0D0 - av_a**(1.5D0)
     lambda(i) = 1.0D0 - av_a**3
!!$        lambda(i) = 1.0D0 - av_a**(2.5)
!!$     ELSE
!!$        lambda(i) = 1.0D0 - av_a**(5.5)
!!$     END IF
  END DO
  av_a = SUM(lambda) / DBLE(SIZE(lambda))

  lambda      =  av_a
  lambda(1)   =  1.0D0
  lambda(no)  =  1.0D0

END SUBROUTINE calc_opt_lambda3_a


SUBROUTINE dist_lin_a(x, y, ymax, dist)

  USE inter_precision, ONLY: DP

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
