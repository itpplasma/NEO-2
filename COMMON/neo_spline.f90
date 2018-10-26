!> Spline arrays
module neo_spline
  use nrtype, only : dp

  integer, parameter ::   mt = 1
  integer, parameter ::   mp = 1
  integer            ::   theta_ind, phi_ind
  integer            ::   ierr

  real(kind=dp) ::  theta_d, phi_d

  !> Spline array for modb
  real(kind=dp),    dimension(:,:,:,:), allocatable :: b_spl

  !> Spline array for geodesic curviture
  real(kind=dp),    dimension(:,:,:,:), allocatable :: k_spl

  !> Spline array for sqrg11
  real(kind=dp),    dimension(:,:,:,:), allocatable :: g_spl

  !> Spline array for parallel derivative
  real(kind=dp),    dimension(:,:,:,:), allocatable :: p_spl

  !> Spline array for quasi-toroidal phi component of b
  real(kind=dp),    dimension(:,:,:,:), allocatable :: q_spl

  !> Spline array for r_nabpsi
  real(kind=dp),    dimension(:,:,:,:), allocatable :: r_spl
end module neo_spline
