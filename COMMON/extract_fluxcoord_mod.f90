module extract_fluxcoord_mod
  use nrtype, only : dp

  integer, save :: load_extract_fluxcoord=1
  integer :: nphinorm
  real(kind=dp) :: psif_extract,theta_extract,psifmin,hpsif
  real(kind=dp) :: psifmax,phifmax,sigcos
  real(kind=dp), dimension(:), allocatable :: phinorm_arr
end module extract_fluxcoord_mod
