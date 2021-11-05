!***********************************************************************
!
! routines for spline interpolation
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

!> \todo own module for these routines and those in spline_cof? or into pspline?
!> \todo abstract interface for the functions, to shorten parameter definitions

subroutine splint_horner3_a(xa,a,b,c,d,swd,m,x_in,f,fp,fpp,fppp,&
                          y,yp,ypp,yppp)
  ! Computes value y(x) dy/dx(x) from cubic spline
  !
  ! Attention - fastest routine; no check at all
  !
  ! Input:
  !         xa(n)         x-values
  !         a(n),b(n),c(n),d(n)  coefs from spline
  !         swd           Switch for derivatives (0: no / 1: yes)
  !         m             powers of leading term
  !         x_in          x-value for y(x_in) and yp(x_in)
  !         f             'leading function' for spline
  !         fp            'leading function' for spline, 1. derivative
  !         fpp           'leading function' for spline, 2. derivative
  !         fppp          'leading function' for spline, 3. derivative
  ! Output:
  !         y             y-value at x_in
  !         yp            dy/dx-value at x_in
  !         ypp           d2y/dx2-value at x_in
  !         yppp          d3y/dx3-value at x_in

  !---------------------------------------------------------------------
  ! Modules
  !---------------------------------------------------------------------
  use nrtype, only : I4B, DP, splinecof_compatibility

  !---------------------------------------------------------------------
  implicit none

  integer(I4B),               intent(in)  :: swd
  real(DP),                   intent(in)  :: m
  real(DP),     dimension(:), intent(in)  :: xa, a, b, c, d
  real(DP),                   intent(in)  :: x_in
  real(DP),                   intent(out) :: y, yp, ypp, yppp
  interface
    function f(x,m)
      use nrtype, only : DP
      implicit none
      real(DP), intent(in) :: x
      real(DP), intent(in) :: m
      real(DP)             :: f
    end function f
  end interface
  interface
    function fp(x,m)
      use nrtype, only : DP
      implicit none
      real(DP), intent(in) :: x
      real(DP), intent(in) :: m
      real(DP)             :: fp
    end function fp
  end interface
  interface
    function fpp(x,m)
      use nrtype, only : DP
      implicit none
      real(DP), intent(in) :: x
      real(DP), intent(in) :: m
      real(DP)             :: fpp
    end function fpp
  end interface
  interface
    function fppp(x,m)
      use nrtype, only : DP
      implicit none
      real(DP), intent(in) :: x
      real(DP), intent(in) :: m
      real(DP)             :: fppp
    end function fppp
  end interface

  integer(I4B) :: klo, khi, n
  integer(I4B) :: k
  real(DP)     :: h, p, p1, p2, p3
  real(DP)     :: x
!-----------------------------------------------------------------------
  real(DP)     :: delta
!-----------------------------------------------------------------------

  n = size(a)

  if (.NOT. ( n == size(b) .AND. n == size(c) &
      .AND. n == size(d) .AND. n == size(xa) ) ) then
    write (*,*) 'splint_horner3: assertion 1 failed'
    write (*,*) size(a), size(b), size(c), size(d), size(xa)
    stop 'program terminated'
  end if

  x  = x_in

  klo=1
  khi=n

  ! Bisection to find k value for which xa(klo) < x < xa(klo+1)
  do while ( (khi-klo) .GT. 1 )
    k=(khi+klo)/2
    if (xa(k).GT.x) then
      khi=k
    else
      klo=k
    end if
  end do

  ! Checks to see if bisection was sucessfull.
  if ((klo < 0) .OR. (klo > n)) then
    print *, 'splint_horner3: n, klo: ', n, klo
    stop
  end if
  if ((khi < 0) .OR. (khi > n)) then
    print *, 'splint_horner3: n, khi: ', n, khi
    stop
  end if

  h = x - xa(klo)

  !! for equidistant points only!
  ! h = xa(2) - xa(1)
  ! klo = int((x - xa(1))/h) + 1
  ! if (klo .le. 0) klo = 1
  ! if (klo .ge. n) klo = n - 1
  ! khi = klo + 1

  if (splinecof_compatibility) then
    !  Linear interpolation
    delta = xa(khi) - xa(klo)
    p = a(klo) + h * (b(klo) + delta * (c(klo) + delta * d(klo)))
  else
    ! Nonlinear interpolation
    p = a(klo) + h * (b(klo) + h * (c(klo) + h * d(klo)))
  end if

  y = f(x,m) * p

  if ( swd .NE. 0 ) then
    p1   = b(klo) + h * (2.0D0 * c(klo) + 3.0D0 * d(klo) * h)
    p2   = 2.0D0 * c(klo) + 6.0D0 * d(klo) * h
    p3   = 6.0D0 * d(klo)
    yp   = fp(x_in,m) * p + f(x_in,m) * p1
    ypp  = fpp(x_in,m) * p + 2.0D0 * fp(x_in,m) * p1 + f(x_in,m) * p2
    yppp = fppp(x_in,m) * p + 3.0D0 * fpp(x_in,m) * p1 &
          + 3.0D0 * fp(x_in,m) * p2 + f(x_in,m) * p3
  else
    yp   = 0.0D0
    ypp  = 0.0D0
    yppp = 0.0D0
  end if

end subroutine splint_horner3_a


subroutine splint_horner3_driv_s_a(svec,a,b,c,d,swd,ixm,ixn,s,theta,phi,&
                                   f,fp,fpp,fppp,y,ys,yt,yp)
  ! driver routine for splint_horner3
  ! y  =  y_mn  *   sin(m*theta - n*phi)
  ! ys =  y_mn' *   sin(m*theta - n*phi)
  ! yt =  y_mn  * m*cos(m*theta - n*phi)
  ! yp = -y_mn  * n*cos(m*theta - n*phi)
  !
  ! Input:
  !         svec()        s-values, dimension (ns)
  !         a(),b(),c(),d()  coefs from spline, dimension (ns,no_cur)
  !         swd           Switch for derivatives (0: no / 1: yes)
  !         ixm           powers of leading term, mode  sin(m*theta-n*phi)
  !         ixn           mode  sin(m*theta-n*phi)
  !         s             s-value for y(s) and ys(s) and...
  !         theta         angle
  !         phi           angle
  !         f             'leading function' for spline
  !         fp            'leading function' for spline, 1. derivative
  !         fpp           'leading function' for spline, 2. derivative
  !         fppp          'leading function' for spline, 3. derivative
  ! Output:
  !         y             y
  !         ys            dy/ds
  !         yt            d2y/dteta
  !         yp            d3y/dphi

  !---------------------------------------------------------------------
  ! Modules
  !---------------------------------------------------------------------
  use nrtype,  only : I4B, DP
  use inter_interfaces, only : splint_horner3

  !---------------------------------------------------------------------
  implicit none

  real(DP),     dimension(:),   intent(in)  :: svec
  real(DP),     dimension(:,:), intent(in)  :: a, b, c, d
  integer(I4B),                 intent(in)  :: swd
  real(DP),     dimension(:),   intent(in)  :: ixm
  integer(I4B), dimension(:),   intent(in)  :: ixn
  real(DP),                     intent(in)  :: s, theta, phi
  real(DP),                     intent(out) :: y, ys, yt, yp
  interface
    function f(x,m)
      use nrtype, only : DP
      implicit none
      real(DP), intent(in) :: x
      real(DP), intent(in) :: m
      real(DP)             :: f
    end function f
  end interface
  interface
    function fp(x,m)
      use nrtype, only : DP
      implicit none
      real(DP), intent(in) :: x
      real(DP), intent(in) :: m
      real(DP)             :: fp
    end function fp
  end interface
  interface
    function fpp(x,m)
      use nrtype, only : DP
      implicit none
      real(DP), intent(in) :: x
      real(DP), intent(in) :: m
      real(DP)             :: fpp
    end function fpp
  end interface
  interface
    function fppp(x,m)
      use nrtype, only : DP
      implicit none
      real(DP), intent(in) :: x
      real(DP), intent(in) :: m
      real(DP)             :: fppp
    end function fppp
  end interface

  integer(I4B) :: no_cur, i
  real(DP)     :: arg, si, co
  real(DP)     :: ay, ays, ayss, aysss

!-----------------------------------------------------------------------

  no_cur = size(a,2)

  y  = 0.0D0
  ys = 0.0D0
  yt = 0.0D0
  yp = 0.0D0

  do i = 1, no_cur
    call splint_horner3(svec, a(:,i), b(:,i), c(:,i), d(:,i), swd, &
                     & ixm(i), s, f, fp, fpp, fppp, ay, ays, ayss, aysss)

    arg = ixm(i) * theta  -  ixn(i) * phi
    si  = sin(arg)
    co  = cos(arg)
    y  =  y  +          ay  * si
    ys =  ys +          ays * si
    yt =  yt + ixm(i) * ay  * co
    yp =  yp - ixn(i) * ay  * co
  end do

end subroutine splint_horner3_driv_s_a


subroutine splint_horner3_driv_c_a(svec,a,b,c,d,swd,ixm,ixn,s,theta,phi,&
                                   f,fp,fpp,fppp,y,ys,yt,yp)
  ! driver routine for splint_horner3
  ! y  =  y_mn  *   cos(m*theta - n*phi)
  ! ys =  y_mn' *   cos(m*theta - n*phi)
  ! yt = -y_mn  * m*sin(m*theta - n*phi)
  ! yp =  y_mn  * n*sin(m*theta - n*phi)
  !
  ! Input:
  !         svec()        s-values, dimension (ns)
  !         a(),b(),c(),d()  coefs from spline, dimension (ns,no_cur)
  !         swd           Switch for derivatives (0: no / 1: yes)
  !         ixm           powers of leading term, mode  sin(m*theta-n*phi)
  !         ixn           mode  sin(m*theta-n*phi)
  !         s             s-value for y(s) and ys(s) and...
  !         theta         angle
  !         phi           angle
  !         f             'leading function' for spline
  !         fp            'leading function' for spline, 1. derivative
  !         fpp           'leading function' for spline, 2. derivative
  !         fppp          'leading function' for spline, 3. derivative
  ! Output:
  !         y             y
  !         ys            dy/ds
  !         yt            d2y/dteta
  !         yp            d3y/dphi

  !---------------------------------------------------------------------
  ! Modules
  !---------------------------------------------------------------------
  use nrtype,  only : I4B, DP
  use inter_interfaces, only : splint_horner3

  !---------------------------------------------------------------------
  implicit none

  real(DP),     dimension(:),   intent(in)  :: svec
  real(DP),     dimension(:,:), intent(in)  :: a, b, c, d
  integer(I4B),                 intent(in)  :: swd
  real(DP),     dimension(:),   intent(in)  :: ixm
  integer(I4B), dimension(:),   intent(in)  :: ixn
  real(DP),                     intent(in)  :: s, theta, phi
  real(DP),                     intent(out) :: y, ys, yt, yp
  interface
    function f(x,m)
      use nrtype, only : DP
      implicit none
      real(DP), intent(in) :: x
      real(DP), intent(in) :: m
      real(DP)             :: f
    end function f
  end interface
  interface
    function fp(x,m)
      use nrtype, only : DP
      implicit none
      real(DP), intent(in) :: x
      real(DP), intent(in) :: m
      real(DP)             :: fp
    end function fp
  end interface
  interface
    function fpp(x,m)
      use nrtype, only : DP
      implicit none
      real(DP), intent(in) :: x
      real(DP), intent(in) :: m
      real(DP)             :: fpp
    end function fpp
  end interface
  interface
    function fppp(x,m)
      use nrtype, only : DP
      implicit none
      real(DP), intent(in) :: x
      real(DP), intent(in) :: m
      real(DP)             :: fppp
    end function fppp
  end interface

  integer(I4B) :: no_cur, i
  real(DP)     :: arg, si, co
  real(DP)     :: ay, ays, ayss, aysss

  !---------------------------------------------------------------------

  no_cur = size(a,2)

  y  = 0.0D0
  ys = 0.0D0
  yt = 0.0D0
  yp = 0.0D0

  do i = 1, no_cur
    call splint_horner3(svec, a(:,i), b(:,i), c(:,i), d(:,i), swd, &
                     & ixm(i), s, f, fp, fpp, fppp, ay, ays, ayss, aysss)

    arg = ixm(i) * theta  -  ixn(i) * phi
    si  = sin(arg)
    co  = cos(arg)
    y  =  y  +          ay  * co
    ys =  ys +          ays * co
    yt =  yt - ixm(i) * ay  * si
    yp =  yp + ixn(i) * ay  * si
  end do

end subroutine splint_horner3_driv_c_a


!> Computes value y(x) dy/dx(x) from linear spline (linear interpolation)
!>
!> Attention - fastest routine; no check at all
!>
!> Input:
!>         xa(n)         x-values
!>         a(n),b(n),c(n),d(n)  coefs from spline
!>         swd           Switch for derivatives (0: no / 1: yes)
!>         m             powers of leading term
!>         x_in          x-value for y(x_in) and yp(x_in)
!>         f             'leading function' for spline
!>         fp            'leading function' for spline, 1. derivative
!>         fpp           'leading function' for spline, 2. derivative
!>         fppp          'leading function' for spline, 3. derivative
!> Output:
!>         y             y-value at x_in
!>         yp            dy/dx-value at x_in
!>         ypp           d2y/dx2-value at x_in (= 0)
!>         yppp          d3y/dx3-value at x_in (= 0)
subroutine splint_horner1_a(xa, a, b, c, d, swd, m, x_in, f, fp, fpp, fppp, &
                        & y, yp, ypp, yppp)
  !---------------------------------------------------------------------
  ! Modules
  !---------------------------------------------------------------------
  use nrtype, only : I4B, DP

  implicit none

  integer(I4B),               intent(in)  :: swd
  real(DP),                   intent(in)  :: m
  real(DP),     dimension(:), intent(in)  :: xa, a, b, c, d
  real(DP),                   intent(in)  :: x_in
  real(DP),                   intent(out) :: y, yp, ypp, yppp
  interface
    function f(x,m)
      use nrtype, only : DP
      implicit none
      real(DP), intent(in) :: x
      real(DP), intent(in) :: m
      real(DP)             :: f
    end function f
  end interface
  interface
    function fp(x,m)
      use nrtype, only : DP
      implicit none
      real(DP), intent(in) :: x
      real(DP), intent(in) :: m
      real(DP)             :: fp
    end function fp
  end interface
  interface
    function fpp(x,m)
      use nrtype, only : DP
      implicit none
      real(DP), intent(in) :: x
      real(DP), intent(in) :: m
      real(DP)             :: fpp
    end function fpp
  end interface
  interface
    function fppp(x,m)
      use nrtype, only : DP
      implicit none
      real(DP), intent(in) :: x
      real(DP), intent(in) :: m
      real(DP)             :: fppp
    end function fppp
  end interface

  integer(I4B) :: klo, khi, n
  integer(I4B) :: k
  real(DP)     :: h
  real(DP)     :: x
  !---------------------------------------------------------------------

  n = size(a)

  if (.not. ( n == size(b) .and. n == size(c) &
      & .and. n == size(d) .and. n == size(xa) ) ) then
    write (*,*) 'splint_horner3: assertion 1 failed'
    write (*,*) size(a), size(b), size(c), size(d), size(xa)
    stop 'program terminated'
  end if

  x  = x_in

  klo = 1
  khi = n

  ! Bisection to find k value for which xa(klo) < x < xa(klo+1)
  do while ( (khi-klo) > 1 )
    k = (khi+klo)/2
    if (xa(k) > x) then
      khi=k
    else
      klo=k
    end if
  end do

  ! Checks to see if bisection was sucessfull.
  if ((klo < 0) .or. (klo > n)) then
    write(*,*) 'splint_horner3: n, klo: ', n, klo
    stop
  end if
  if ((khi < 0) .or. (khi > n)) then
    write(*,*) 'splint_horner3: n, khi: ', n, khi
    stop
  end if

  h = x - xa(klo)

  y = a(klo) + h * b(klo)

  if ( swd .NE. 0 ) then
    yp   = b(klo)
  else
    yp   = 0.0D0
  end if

  ypp  = 0.0D0
  yppp = 0.0D0

end subroutine splint_horner1_a
