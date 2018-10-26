module extremum_mod

  use magnetics_mod
  use plagrange_mod

  implicit none

  integer, parameter :: nlagrange = 5

  public find_extremum
  private find_ext
  interface find_extremum
     module procedure find_ext
  end interface

contains

  subroutine find_ext(fieldperiod,x1i,x2i,dxi,x,y)
    use nrtype, only : dp

    type(fieldperiod_struct), pointer :: fieldperiod
    real(kind=dp), intent(in)   :: x1i,x2i,dxi
    real(kind=dp), intent(out)  :: x,y

    real(kind=dp) :: A,d,x_hlp,d2,dummy,dx,x1_p,x2_p
    real(kind=dp) :: x1_in,x2_in,x1,x2,x3,x4
    real(kind=dp) :: fx1,fx2,fxm,fx3,fx4
    integer       :: n,k
    logical       :: L

    x1 = x1i
    x2 = x2i

    A = 2.0d0 / (1.0d0 + SQRT(5.0d0))
    n = 0
    dx = ABS(dxi)

    if (x2 .LT. x1) then
      x_hlp = x2
      x2 = x1
      x1 = x_hlp
    end if

    x1_in = x1
    x2_in = x2

    ! second derivative
    call plagrange_interp(fieldperiod,x1,nlagrange,fx1,dummy)
    call plagrange_interp(fieldperiod,x2,nlagrange,fx2,dummy)
    call plagrange_interp(fieldperiod,(x1+x2)/2.0d0,nlagrange,fxm,dummy)
    d2 = (fx1+fx2-2.0d0*fxm) / ((x2-x1)/2.0d0)**2

    do k = 1,2

      d = x2 - x1
      x3 = x1 + A * d
      x4 = x2 - A * d

      call plagrange_interp(fieldperiod,x3,nlagrange,fx3,dummy)
      call plagrange_interp(fieldperiod,x4,nlagrange,fx4,dummy)

      do while (ABS(d)/ABS(MAX(x1_in,x2_in)) .GT. dx)
        n = n + 1
        if (d2 < 0) then
          L = fx4 < fx3
        else
          L = fx3 < fx4
        end if
        x1_p = x1
        x2_p = x2
        if (L) then
          x1 = x4
          x4 = x3
          fx4 = fx3
          d = x2 - x1
          x3 = x1 + A * d
          call plagrange_interp(fieldperiod,x3,nlagrange,fx3,dummy)
        else
          x2 = x3
          x3 = x4
          fx3 = fx4
          d = x2 - x1
          x4 = x2 - A * d
          call plagrange_interp(fieldperiod,x4,nlagrange,fx4,dummy)
        end if
        !if (x1 .lt. 0.1d0) print *, x1,x2
        if (x1_p .EQ. x1 .AND. x2_p .EQ. x2) then
          exit
        end if
      end do
      !print *, 'n ',n,d2,k

      call plagrange_interp(fieldperiod,x1,nlagrange,fx1,dummy)
      call plagrange_interp(fieldperiod,x2,nlagrange,fx2,dummy)

      if (d2 .LT. 0) then
        L = fx2 < fx1
      else
        L = fx1 < fx2
      end if

      if (L) then
        x = x1
        y = fx1
      else
        x = x2
        y = fx2
      end if

      if (x .GT. x1_in .AND. x .LT. x2_in) then
        exit
      elseif (x .LE. x1_in) then
        x2 = x1_in
        x1 = x1_in - (x2_in - x1_in)
      else
        x1 = x2_in
        x2 = x2_in + (x2_in - x1_in)
      end if
    end do

  end subroutine find_ext

end module extremum_mod
