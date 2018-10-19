
subroutine kin_allocate(lalloc)

  use size_mod, only : ndim
  use rk4_kin_mod, only : y, dydx, yt, dyt, dym

  implicit none

  logical, intent(in) :: lalloc

  if (lalloc) then
    allocate(y(ndim), dydx(ndim), yt(ndim), dyt(ndim), dym(ndim))
  else
    deallocate(y, dydx, yt, dyt, dym)
  end if

end subroutine kin_allocate
