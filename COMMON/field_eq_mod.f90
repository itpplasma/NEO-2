module field_eq_mod
  integer :: icall_eq=0
  integer :: nrad,nzet,icp,nwindow_r,nwindow_z
  real(kind=8), parameter                      :: pi=3.14159265358979d0
  real(kind=8) :: psib,btf,rtf,hrad,hzet
  real(kind=8), dimension(:,:), allocatable    :: psi, psi0
  real(kind=8), dimension(:,:,:), allocatable  :: splpsi
  real(kind=8), dimension(:), allocatable      :: rad, zet, xi,f
  integer, dimension(:), allocatable           :: imi,ima,jmi,jma
  integer, dimension(:,:), allocatable         :: ipoint
  double precision :: psif,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2
end module field_eq_mod
