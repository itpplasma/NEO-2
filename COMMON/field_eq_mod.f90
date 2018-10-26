module field_eq_mod
  use nrtype, only : dp

  integer :: icall_eq=0
  integer :: nrad,nzet,icp,nwindow_r,nwindow_z

  real(kind=dp) :: psib,btf,rtf,hrad,hzet
  real(kind=dp), dimension(:,:), allocatable    :: psi, psi0
  real(kind=dp), dimension(:,:,:), allocatable  :: splpsi
  real(kind=dp), dimension(:), allocatable      :: rad, zet, xi,f
  integer, dimension(:), allocatable           :: imi,ima,jmi,jma
  integer, dimension(:,:), allocatable         :: ipoint
  real(kind=dp) :: psif,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2
end module field_eq_mod
