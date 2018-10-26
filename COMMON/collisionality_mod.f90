module collisionality_mod
  use nrtype, only : dp

  integer          :: isw_lorentz
  integer          :: isw_integral
  integer          :: isw_energy
  integer          :: isw_axisymm
  integer          :: isw_relativistic
  real(kind=dp) :: T_e

  integer          :: isw_momentum
  integer          :: vel_num
  integer          :: vel_distri_swi
  integer          :: nvel
  real(kind=dp) :: vel_max

  real(kind=dp) :: collpar
  real(kind=dp) :: conl_over_mfp
  real(kind=dp) :: coeps
  real(kind=dp), allocatable, dimension(:) :: y_axi_averages
  real(kind=dp), allocatable, dimension(:) :: vel_array

  !**********************************************************
  ! Changes required for multispecies support
  !**********************************************************
  logical          :: lsw_multispecies
  ! isw_coul_log = 0: Coulomb logarithm set as species independent (overrides values for n_spec)
  ! isw_coul_log = 1: Coulomb logarithm computed for each species using n_spec, T_spec
  !                   (overrides values for collisionality parameters)
  integer          :: isw_coul_log
  integer          :: num_spec ! number of species
  integer, dimension(:), allocatable :: species_tag
  real(kind=dp), dimension(:), allocatable :: conl_over_mfp_spec
  real(kind=dp), dimension(:), allocatable :: collpar_spec
  real(kind=dp), dimension(:), allocatable :: z_spec ! species charge number
  real(kind=dp), dimension(:), allocatable :: m_spec ! species mass [g]
  real(kind=dp), dimension(:), allocatable :: n_spec ! species density [cm^-3]
  real(kind=dp), dimension(:), allocatable :: T_spec ! species temperature [erg]

  real(kind=dp) :: collpar_min,collpar_max,v_min_resolution,v_max_resolution
  real(kind=dp) :: phi_x_max
  real(kind=dp) :: collop_bspline_dist
  integer          :: collop_bspline_order
  logical          :: collop_bspline_taylor

  real(kind=dp) :: m_nbi, T_nbi
  logical          :: lsw_nbi

end module collisionality_mod
