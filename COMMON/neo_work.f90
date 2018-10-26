!> Working parameters
module neo_work
  use nrtype, only : dp

  real(kind=dp) ::   theta_start
  real(kind=dp) ::   theta_end
  real(kind=dp) ::   theta_int

  real(kind=dp) ::   phi_start
  real(kind=dp) ::   phi_end
  real(kind=dp) ::   phi_int

! real(kind=dp),    dimension(:,:,:),   allocatable :: cosval,sinval
  real(kind=dp),    dimension(:,:),     allocatable :: cosmth,sinmth
  real(kind=dp),    dimension(:,:),     allocatable :: cosnph,sinnph
  real(kind=dp),    dimension(:),       allocatable :: theta_arr,phi_arr
  real(kind=dp),    dimension(:,:),     allocatable :: r,z,l,b
  real(kind=dp),    dimension(:,:),     allocatable :: r_tb,z_tb,p_tb,b_tb
  real(kind=dp),    dimension(:,:),     allocatable :: r_pb,z_pb,p_pb,b_pb
  real(kind=dp),    dimension(:,:),     allocatable :: gtbtb,gpbpb,gtbpb
  real(kind=dp),    dimension(:,:),     allocatable :: isqrg,sqrg11,kg,pard
  real(kind=dp),    dimension(:,:),     allocatable :: bqtphi
  real(kind=dp),    dimension(:,:),     allocatable :: r_nabpsi
  real(kind=dp),    dimension(:,:),     allocatable :: psi_r, psi_z
end module neo_work
