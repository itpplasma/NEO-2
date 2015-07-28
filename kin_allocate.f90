!
  subroutine kin_allocate(ialloc)
!
  use size_mod
  use rk4_kin_mod
!
  implicit none
!
  integer :: ialloc
!
  if(ialloc.eq.1) then
    allocate(y(ndim),dydx(ndim),yt(ndim),dyt(ndim),dym(ndim))
  else
    deallocate(y,dydx,yt,dyt,dym)
  endif
!
  return
  end
