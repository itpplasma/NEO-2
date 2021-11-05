
! testfunction for spline tf(x, m)
!           1. derivative tfp(x, m)
!           2. derivative tfpp(x, m)
!           3. derivative tfppp(x, m)

function tf(x,m)
  ! calculate testfunction for spline tf(x, m)

  !---------------------------------------------------------------------
  ! Modules
  !---------------------------------------------------------------------
  use nrtype, only : DP

  !---------------------------------------------------------------------
  implicit none

  real(DP), intent(in) :: x
  real(DP), intent(in) :: m
  real(DP)             :: tf

  !---------------------------------------------------------------------
  if (m .NE. 0.0D0) then
     if (x == 0.0D0) then
        tf = 0.0D0
     else
        tf = x**m
     end if
  else
     tf = 1.0D0
  end if

end function tf


function tfp(x,m)
  ! calculate 1. derivative of testfunction for spline tfp(x, m)

  !---------------------------------------------------------------------
  ! Modules
  !---------------------------------------------------------------------
  use nrtype, only : DP

  !---------------------------------------------------------------------
  implicit none

  real(DP), intent(in) :: x
  real(DP), intent(in) :: m
  real(DP)             :: tfp

  !---------------------------------------------------------------------
  if ((m - 1.0D0) .NE. 0.0D0) then
     if (x == 0.0D0) then
        tfp = 0.0D0
     else
        tfp = m * x**(m-1.0D0)
     end if
  else
     tfp = 1.0D0
  end if

end function tfp


function tfpp(x,m)
  ! calculate 2. derivative of testfunction for spline tfpp(x, m)

  !---------------------------------------------------------------------
  ! Modules
  !---------------------------------------------------------------------
  use nrtype, only : DP

  !---------------------------------------------------------------------
  implicit none

  real(DP), intent(in) :: x
  real(DP), intent(in) :: m
  real(DP)             :: tfpp

  !---------------------------------------------------------------------
  if ((m - 2.0D0) .NE. 0.0D0) then
     if (x == 0.0D0) then
        tfpp = 0.0D0
     else
        tfpp = m * (m - 1.0D0) * x**(m-2.0D0)
     end if
  else
     tfpp = 1.0D0
  end if

end function tfpp


function tfppp(x,m)
  ! calculate 3. derivative of testfunction for spline tfppp(x, m)

  !---------------------------------------------------------------------
  ! Modules
  !---------------------------------------------------------------------
  use nrtype, only : DP

  !---------------------------------------------------------------------
  implicit none

  real(DP), intent(in) :: x
  real(DP), intent(in) :: m
  real(DP)             :: tfppp

  !---------------------------------------------------------------------
  if ((m - 3.0D0) .NE. 0.0D0) then
     if (x == 0.0D0) then
        tfppp = 0.0D0
     else
        tfppp = m * (m - 1.0D0) * (m - 2.0D0) * x**(m-3.0D0)
     end if
  else
     tfppp = 1.0D0
  end if

end function tfppp


function tfone(x,m)
  ! calculate testfunction for spline tf(x, m); here just 1

  !---------------------------------------------------------------------
  ! Modules
  !---------------------------------------------------------------------
  use nrtype, only : DP

  !---------------------------------------------------------------------
  implicit none

  real(DP), intent(in) :: x
  real(DP), intent(in) :: m
  real(DP)             :: tfone

  !---------------------------------------------------------------------
  tfone = 1.0D0
end function tfone


function tfzero(x,m)
  ! calculate testfunction for spline tf(x, m); here just 1

  !---------------------------------------------------------------------
  ! Modules
  !---------------------------------------------------------------------
  use nrtype, only : DP

  !---------------------------------------------------------------------
  implicit none

  real(DP), intent(in) :: x
  real(DP), intent(in) :: m
  real(DP)             :: tfzero

  !---------------------------------------------------------------------
  tfzero = 0.0D0
end function tfzero
