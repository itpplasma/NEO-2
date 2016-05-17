MODULE spline_mod

CONTAINS

SUBROUTINE spl2d(nx,ny,hx,hy,mx,my,f,spl)

! Makes a 2-dimensional cubic spline of function f(x,y)
!
! Input:  nx, ny              number of values in x and y  
!         hx, hy              step size in x and y (aequidistant)
!         mx, my              spline mode (0: standard, 1: periodic)
!         f(nx,ny)            f(x,y)-values
! Output: spl                 Array with spline parameters

  USE neo_precision

  IMPLICIT NONE

  INTEGER,                             INTENT(in)  :: nx, ny, mx, my
  REAL(kind=dp),                       INTENT(in)  :: hx, hy
  REAL(kind=dp), DIMENSION(nx,ny)    , INTENT(in)  :: f
  REAL(kind=dp), DIMENSION(4,4,nx,ny), INTENT(out) :: spl
  
  REAL(kind=dp), DIMENSION(:),         ALLOCATABLE :: bi, ci, di, s
  INTEGER                                          :: i, j, k

  ALLOCATE ( bi(nx), ci(nx), di(nx), s(nx) )
  DO j = 1,ny
     DO i = 1,nx
        s(i) = f(i,j)
     END DO
     IF (mx .EQ. 0) THEN
        CALL splreg(nx,hx,s,bi,ci,di)
     ELSE
        CALL splper(nx,hx,s,bi,ci,di)
     ENDIF
     DO i = 1,nx
        spl(1,1,i,j) = s(i)
        spl(2,1,i,j) = bi(i)
        spl(3,1,i,j) = ci(i)
        spl(4,1,i,j) = di(i)
     END DO
  END DO
  DEALLOCATE ( bi, ci, di, s )

  ALLOCATE ( bi(ny), ci(ny), di(ny), s(ny) )
  DO k = 1,4
     DO i = 1,nx
        DO j = 1,ny
           s(j) = spl(k,1,i,j)
        END DO
        IF (my .EQ. 0) THEN
           CALL splreg(ny,hy,s,bi,ci,di)
        ELSE
           CALL splper(ny,hy,s,bi,ci,di)
        ENDIF
        DO j=1,ny
           spl(k,2,i,j)=bi(j)
           spl(k,3,i,j)=ci(j)
           spl(k,4,i,j)=di(j)
        END DO
     END DO
  END DO
  DEALLOCATE ( bi, ci, di, s )

  RETURN
END SUBROUTINE spl2d
!=====================================================
SUBROUTINE eva2d(nx,ny,ix,iy,dx,dy,spl,spval)

! Evaluates a 2-dimensional cubic spline of function f(x,y)
!
! Input:  nx, ny              number of values in x and y  
!         ix, iy              pointer into the spline array spl
!         dx, dy              distance from x(ix) and y(iy)
!         spl                 array with spline data
! Output: spval               evaluated function value

  USE neo_precision
  
  IMPLICIT NONE
  
  INTEGER,                             INTENT(in)  :: nx, ny, ix, iy
  REAL(kind=dp),                       INTENT(in)  :: dx, dy
  REAL(kind=dp), DIMENSION(4,4,nx,ny), INTENT(in)  :: spl
  REAL(kind=dp),                       INTENT(out) :: spval

  REAL(kind=dp), DIMENSION(4)                      :: a
  INTEGER                                          :: l

  DO l=1,4
     a(l) = spl(1,l,ix,iy) + dx*(spl(2,l,ix,iy) +              & 
          dx*(spl(3,l,ix,iy) + dx* spl(4,l,ix,iy)))
  END DO
  spval = a(1)+dy*(a(2)+dy*(a(3)+dy*a(4)))

  RETURN
END SUBROUTINE eva2d
!=====================================================
SUBROUTINE eva2d_fd(nx,ny,ix,iy,dx,dy,spl,spval)

! Evaluates the first derivatives of 2-dimensional cubic spline of function f(x,y)
!
! Input:  nx, ny              number of values in x and y  
!         ix, iy              pointer into the spline array spl
!         dx, dy              distance from x(ix) and y(iy)
!         spl                 array with spline data
! Output: spval(2)            evaluated function value
!                             spval(1) = df/dx
!                             spval(2) = df/dy

  USE neo_precision
  
  IMPLICIT NONE
  
  INTEGER,                             INTENT(in)  :: nx, ny, ix, iy
  REAL(kind=dp),                       INTENT(in)  :: dx, dy
  REAL(kind=dp), DIMENSION(4,4,nx,ny), INTENT(in)  :: spl
  REAL(kind=dp), DIMENSION(2),         INTENT(out) :: spval

  INTEGER                                          :: i,j
  REAL(kind=dp)                                    :: muli, mulj

  spval = 0.0_dp

! df/dx
  DO i=2,4
     IF (i == 2) THEN
        muli = 1.0_dp
     ELSE
        muli = dx**(i-2)
     END IF
     muli = muli * (i-1)
     DO j=1,4
       IF (j == 1) THEN
           mulj = 1.0_dp
        ELSE
           mulj = dy**(j-1)
        END IF
        spval(1) = spval(1) + spl(i,j,ix,iy) * muli * mulj
     END DO
  END DO

! df/dy
  DO i=1,4
     IF (i == 1) THEN
        muli = 1.0_dp
     ELSE
        muli = dx**(i-1)
     END IF
     DO j=2,4
        IF (j == 2) THEN
           mulj = 1.0_dp
        ELSE
           mulj = dy**(j-2)
        END IF
        mulj = mulj * (j-1)
        spval(2) = spval(2) + spl(i,j,ix,iy) * muli * mulj 
     END DO
  END DO

  RETURN
END SUBROUTINE eva2d_fd
!=====================================================
SUBROUTINE eva2d_sd(nx,ny,ix,iy,dx,dy,spl,spval)

! Evaluates the second derivatives of 2-dimensional cubic spline of function f(x,y)
!
! Input:  nx, ny              number of values in x and y  
!         ix, iy              pointer into the spline array spl
!         dx, dy              distance from x(ix) and y(iy)
!         spl                 array with spline data
! Output: spval(3)            evaluated function values
!                             spval(1) = d^2f/dx^2
!                             spval(2) = d^2f/(dxdy)
!                             spval(3) = d^2f/dy^2

  USE neo_precision
  
  IMPLICIT NONE
  
  INTEGER,                             INTENT(in)  :: nx, ny, ix, iy
  REAL(kind=dp),                       INTENT(in)  :: dx, dy
  REAL(kind=dp), DIMENSION(4,4,nx,ny), INTENT(in)  :: spl
  REAL(kind=dp), DIMENSION(3),         INTENT(out) :: spval

  INTEGER                                          :: i,j
  REAL(kind=dp)                                    :: muli, mulj

  spval = 0.0_dp

! d^2f/dx^2
  DO i=3,4
     IF (i == 3) THEN
        muli = 1.0_dp
     ELSE
        muli = dx**(i-3)
     END IF
     muli = muli * (i-1) * (i-2)
     DO j=1,4
        IF (j == 1) THEN
           mulj = 1.0_dp
        ELSE
           mulj = dy**(j-1)
        END IF
        spval(1) = spval(1) + spl(i,j,ix,iy) * muli * mulj
     END DO
  END DO

! d^2f/(dxdy)
  DO i=2,4
     IF (i == 2) THEN
        muli = 1.0_dp
     ELSE
        muli = dx**(i-2)
     END IF
     muli = muli * (i-1)
     DO j=2,4
        IF (j == 2) THEN
           mulj = 1.0_dp
        ELSE
           mulj = dy**(j-2)
        END IF
        mulj = mulj * (j-1)
        spval(2) = spval(2) + spl(i,j,ix,iy) * muli * mulj
     END DO
  END DO

! d^2f/dy^2
  DO i=1,4
     IF (i == 1) THEN
        muli = 1.0_dp
     ELSE
        muli = dx**(i-1)
     END IF
     DO j=3,4
        IF (j == 3) THEN
           mulj = 1.0_dp
        ELSE
           mulj = dy**(j-3)
        END IF
        mulj = mulj * (j-1) * (j-2)
        spval(3) = spval(3) + spl(i,j,ix,iy) * muli * mulj
     END DO
  END DO

  RETURN
END SUBROUTINE eva2d_sd
!=====================================================
SUBROUTINE splreg(n,h,y,bi,ci,di)

! Makes a cubic spline of function y(x)
!
! Input:  n                   number of values in y  
!         h                   step size in x (aequidistant)
!         y(n)                y-values
! Output: bi(n),ci(n),di(n)   Spline parameters

  USE neo_precision

  IMPLICIT NONE

  INTEGER,                     INTENT(in)  :: n
  REAL(kind=dp),               INTENT(in)  :: h
  REAL(kind=dp), DIMENSION(n), INTENT(in)  :: y
  REAL(kind=dp), DIMENSION(n), INTENT(out) :: bi, ci, di

  REAL(kind=dp)                            :: ak1, ak2, am1, am2, c, e, c1
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: al, bt
  INTEGER                                  :: k, n2, i, i5

  ALLOCATE ( al(n), bt(n) )

  ak1 = 0.d0
  ak2 = 0.d0
  am1 = 0.d0
  am2 = 0.d0
  k = n-1
  al(1) = ak1
  bt(1) = am1
  n2 = n-2
  c = -4.d0*h
  DO i = 1,n2
     e = -3.d0*((y(i+2)-y(i+1))-(y(i+1)-y(i)))/h
     c1 = c-al(i)*h
     al(i+1) = h/c1
     bt(i+1) = (h*bt(i)+e)/c1
  END DO
  ci(n) = (am2+ak2*bt(k))/(1.d0-al(k)*ak2)
  DO i = 1,k
     i5 = n-i
     ci(i5) = al(i5)*ci(i5+1)+bt(i5)
  END DO
  n2 = n-1
  DO i = 1,n2
     bi(i) = (y(i+1)-y(i))/h-h*(ci(i+1)+2.d0*ci(i))/3.d0
     di(i) = (ci(i+1)-ci(i))/h/3.d0
  END DO
  DEALLOCATE ( al, bt )

  RETURN
END SUBROUTINE splreg
!=====================================================
SUBROUTINE splper(n,h,y,bi,ci,di)

! Makes a cubic spline of periodic function y(x)
!
! Input:  n                   number of values in y
!         h                   step size in x (aequidistant)
!         y(n)                y-values
! Output: bi(n),ci(n),di(n)   Spline parameters

  USE neo_precision

  IMPLICIT NONE

  INTEGER,                     INTENT(in)  :: n
  REAL(kind=dp),               INTENT(in)  :: h
  REAL(kind=dp), DIMENSION(n), INTENT(in)  :: y
  REAL(kind=dp), DIMENSION(n), INTENT(out) :: bi, ci, di

  REAL(kind=dp)                            :: psi, ss
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: bmx, yl
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: amx1, amx2, amx3
  INTEGER                                  :: nmx, n1, n2, i, i1

  ALLOCATE ( bmx(n), yl(n), amx1(n), amx2(n), amx3(n) )
  
  bmx(1) = 1.d30

  nmx=n-1
  n1=nmx-1
  n2=nmx-2
  psi=3.d0/h/h

  CALL spfper(n,amx1,amx2,amx3)

  bmx(nmx) = (y(nmx+1)-2.d0*y(nmx)+y(nmx-1))*psi
  bmx(1)   =(y(2)-y(1)-y(nmx+1)+y(nmx))*psi
  DO i = 3,nmx
     bmx(i-1) = (y(i)-2.d0*y(i-1)+y(i-2))*psi
  END DO
  yl(1) = bmx(1)/amx1(1)
  DO i = 2,n1
     i1 = i-1
     yl(i) = (bmx(i)-yl(i1)*amx2(i1))/amx1(i)
  END DO
  ss = 0.d0
  DO i = 1,n1
     ss = ss+yl(i)*amx3(i)
  END DO
  yl(nmx) = (bmx(nmx)-ss)/amx1(nmx)
  bmx(nmx) = yl(nmx)/amx1(nmx)
  bmx(n1) = (yl(n1)-amx2(n1)*bmx(nmx))/amx1(n1)
  DO i = n2,1,-1
     bmx(i) = (yl(i)-amx3(i)*bmx(nmx)-amx2(i)*bmx(i+1))/amx1(i)
  END DO
  DO i = 1,nmx
     ci(i) = bmx(i)
  END DO

  DO i = 1,n1
     bi(i) = (y(i+1)-y(i))/h-h*(ci(i+1)+2.d0*ci(i))/3.d0
     di(i) = (ci(i+1)-ci(i))/h/3.d0
  END DO
  bi(nmx) = (y(n)-y(n-1))/h-h*(ci(1)+2.d0*ci(nmx))/3.d0
  di(nmx) = (ci(1)-ci(nmx))/h/3.d0
!
! Fix of problems at upper periodicity boundary
!
  bi(n) = bi(1)
  ci(n) = ci(1)
  di(n) = di(1)

  DEALLOCATE ( bmx, yl, amx1, amx2, amx3 )

  RETURN
END SUBROUTINE splper
!=====================================================
SUBROUTINE spfper(np1,amx1,amx2,amx3)

! Helper routine for splfi

  USE neo_precision

  IMPLICIT NONE

  INTEGER,                       INTENT(in)  :: np1
  REAL(kind=dp), DIMENSION(np1), INTENT(out) :: amx1, amx2, amx3
  REAL(kind=dp)                              :: beta, ss
  INTEGER                                    :: n, n1, i, i1

  n = np1-1

  n1 = n-1
  amx1(1) = 2.d0
  amx2(1) = 0.5d0
  amx3(1) = 0.5d0
  amx1(2) = SQRT(15.d0)/2.d0
  amx2(2) = 1.d0/amx1(2)
  amx3(2) = -.25d0/amx1(2)
  beta = 3.75d0
  DO i = 3,n1
     i1 = i-1
     beta = 4.d0-1.d0/beta
     amx1(i) = SQRT(beta)
     amx2(i) = 1.d0/amx1(i)
     amx3(i) = -amx3(i1)/amx1(i)/amx1(i1)
  END DO
  amx3(n1) = amx3(n1)+1.d0/amx1(n1)
  amx2(n1) = amx3(n1)
  ss = 0.0d0
  DO i = 1,n1
     ss = ss+amx3(i)*amx3(i)
  END DO
  amx1(n) = SQRT(4.d0-ss)

  RETURN
END SUBROUTINE spfper
!=====================================================
SUBROUTINE poi2d(hx,hy,mx,my,                         &
     xmin,xmax,ymin,ymax,                             &
     x,y,ix,iy,dx,dy,ierr)
! Creates Pointers for eva2d
!
! Input:  hx, hy              increment in x and y  
!         mx, my              standard (0) or periodic (1) spline
!         xmin, xmax          Minimum and maximum x
!         ymin, ymax          Minimum and maximum y
!         x, y                x and y values for spline avaluation
! Output: spval               evaluated function value
!         ix, iy              pointer into the spline array spl
!         dx, dy              distance from x(ix) and y(iy)
!         ierr                error (> 0)

  USE neo_precision

  IMPLICIT NONE

  REAL(kind=dp),                       INTENT(in)  :: hx, hy
  INTEGER,                             INTENT(in)  :: mx, my
  REAL(kind=dp),                       INTENT(in)  :: xmin, xmax, ymin, ymax
  REAL(kind=dp),                       INTENT(in)  :: x, y

  INTEGER,                             INTENT(out) :: ix, iy
  REAL(kind=dp),                       INTENT(out) :: dx, dy
  INTEGER,                             INTENT(out) :: ierr

  REAL(kind=dp)                                    :: dxx, x1, dyy, y1
  REAL(kind=dp)                                    :: dxmax, dymax

  ierr = 0

  dxx = x-xmin
  IF (mx .EQ. 0) THEN
     IF (dxx .LT. 0.d0) THEN
        ierr = 1
        RETURN
     END IF
     IF (x .GT. xmax) THEN
        ierr = 2
        RETURN
     END IF
  ELSE
     dxmax = xmax - xmin
     IF(dxx .LT. 0.d0) THEN
        dxx = dxx+DBLE(float(1+INT(ABS(dxx/dxmax))))*dxmax
     ELSE IF(dxx .GT. dxmax) THEN
        dxx = dxx-DBLE(float(INT(ABS(dxx/dxmax))))*dxmax
     END IF
  END IF
  x1 = dxx/hx
  ix = INT(x1)
  dx = hx*(x1-DBLE(float(ix)))
  ix = ix+1

  dyy = y-ymin
  IF (my .EQ. 0) THEN
     IF (dyy .LT. 0.d0) THEN
        ierr = 3
        RETURN
     END IF
     IF (y .GT. ymax) THEN
        ierr = 4
        RETURN
     END IF
  ELSE
     dymax = ymax - ymin
     IF(dyy .LT. 0.d0) THEN
        dyy = dyy+DBLE(float(1+INT(ABS(dyy/dymax))))*dymax
     ELSE IF(dyy .GT. dymax) THEN
        dyy = dyy-DBLE(float(INT(ABS(dyy/dymax))))*dymax
     END IF
  END IF
  y1 = dyy/hy
  iy = INT(y1)
  dy = hy*(y1-DBLE(float(iy)))
  iy = iy+1

  RETURN
END SUBROUTINE poi2d

END MODULE spline_mod





