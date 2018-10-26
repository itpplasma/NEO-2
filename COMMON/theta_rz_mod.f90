module theta_rz_mod
  use nrtype, only : dp

  integer :: icall=0
  integer :: nsqp,nlab,nthe,icp_pt
  integer, dimension(:,:), allocatable :: ipoint_pt
  real(kind=dp) :: hsqpsi,hlabel,htheqt,psiaxis,sigma_qt,raxis,zaxis
  real(kind=dp), dimension(:,:),   allocatable :: spllabel
  real(kind=dp), dimension(:,:,:), allocatable :: splthet
  real(kind=dp), dimension(:),     allocatable :: sqpsi,flab,theqt
end module theta_rz_mod
