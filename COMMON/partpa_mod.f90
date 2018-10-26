module partpa_mod
  use nrtype, only : dp

  ! Winny ipmin,iminr added
  integer                                  :: ipmax,ipmin,istart,iminr
  integer                                  :: npassing_start,npass_max
  real(kind=dp)                            :: pardeb0,bmod0,dirint,hxeta
  real(kind=dp), dimension(:), allocatable :: eta
  ! Winny for talking between
  integer :: nhalfstep
  ! Winny end
end module
