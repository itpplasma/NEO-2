
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
SUBROUTINE splinecof3_original_dense(x, y, c1, cn, lambda1, indx, sw1, sw2, &
     a, b, c, d, m, f)
  !-----------------------------------------------------------------------
  ! Modules
  !-----------------------------------------------------------------------

  use nrtype, only : I4B, DP
  USE inter_interfaces, ONLY: calc_opt_lambda3
  !! Modifications by Andreas F. Martitsch (06.08.2014)
  !Replace standard solver from Lapack with sparse solver
  !(Bad performance for more than 1000 flux surfaces ~ (3*nsurf)^2)
  USE sparse_mod, ONLY : sparse_solve
  !! End Modifications by Andreas F. Martitsch (06.08.2014)

  !---------------------------------------------------------------------

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

  INTEGER(I4B), PARAMETER :: VAR = 7
  INTEGER(I4B)            :: size_dimension
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
  character(200) :: error_message

  len_x    = SIZE(x)
  len_indx = SIZE(indx)
  size_dimension = VAR * len_indx - 2

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

  ALLOCATE(MA(size_dimension, size_dimension),  stat = i_alloc, errmsg=error_message)
  if(i_alloc /= 0) then
    write(*,*) 'splinecof3: Allocation for array ma failed with error message:'
    write(*,*) trim(error_message)
    write(*,*) 'size should be ', size_dimension, ' x ', size_dimension
    stop
  end if
  ALLOCATE(inh(size_dimension), indx_lu(size_dimension),  stat = i_alloc, errmsg=error_message)
  if(i_alloc /= 0) then
    write(*,*) 'splinecof3: Allocation for arrays inh and indx_lu failed with error message:'
    write(*,*) trim(error_message)
    write(*,*) 'size should be ', size_dimension
    stop
  end if
  ALLOCATE(simqa(size_dimension*size_dimension),  stat = i_alloc, errmsg=error_message)
  if(i_alloc /= 0) then
    write(*,*) 'splinecof3: Allocation for array simqa failed with error message:'
    write(*,*) trim(error_message)
    write(*,*) 'size should be ', size_dimension*size_dimension
    stop
  end if
  ALLOCATE(lambda(SIZE(lambda1)),  stat = i_alloc, errmsg=error_message)
  if(i_alloc /= 0) then
    write(*,*) 'splinecof3: Allocation for array lambda failed with error message:'
    write(*,*) trim(error_message)
    write(*,*) 'size should be ', size(lambda1)
    stop
  end if
  ALLOCATE(omega(SIZE(lambda1)),  stat = i_alloc, errmsg=error_message)
  if(i_alloc /= 0) then
    write(*,*) 'splinecof3: Allocation for array omega failed with message:'
    write(*,*) trim(error_message)
    write(*,*) 'size should be ', size(lambda1)
    stop
  end if
  !---------------------------------------------------------------------


  IF (DABS(c1) > 1.0E30) THEN
    c1 = 0.0D0;
  END IF
  IF (DABS(cn) > 1.0E30) THEN
    cn = 0.0D0;
  END IF

  ! setting all to zero
  MA(:,:) = 0.0D0
  inh(:)  = 0.0D0

  ! calculate optimal weights for smooting (lambda)
  IF ( MAXVAL(lambda1) < 0.0D0 ) THEN
    CALL calc_opt_lambda3(x, y, omega)
  ELSE
    omega  = lambda1
  END IF
  lambda = 1.0D0 - omega

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

  ! boundary condition 2
  i = i + 1
  MA(i, 2) = DBLE(mu2)
  MA(i, 3) = DBLE(nu2)
  MA(i, (len_indx-1)*VAR + 2) = DBLE(sig2)
  MA(i, (len_indx-1)*VAR + 3) = DBLE(rho2)
  inh(i) = cn

! ---------------------------

  ! solve system
  CALL sparse_solve(MA, inh)

  ! take a(), b(), c(), d()
  DO i = 1, len_indx
     a(i) = inh((i-1)*VAR+1)
     b(i) = inh((i-1)*VAR+2)
     c(i) = inh((i-1)*VAR+3)
     d(i) = inh((i-1)*VAR+4)
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

END SUBROUTINE splinecof3_original_dense
