module bdivfree_mod
  use nrtype, only : dp

  integer :: nr,nz,ntor,icp
  integer, dimension(:,:), allocatable :: ipoint
  real(kind=dp) :: rmin,zmin,hr,hz,pmin,pfac
  real(kind=dp), dimension(:),       allocatable :: rpoi,zpoi
  real(kind=dp), dimension(:,:,:),   allocatable :: apav,rbpav_coef
  real(kind=dp), dimension(:,:,:,:), allocatable :: aznre,aznim,arnre,arnim
end module bdivfree_mod
