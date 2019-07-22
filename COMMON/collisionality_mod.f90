MODULE collisionality_mod
  INTEGER          :: isw_lorentz
  INTEGER          :: isw_integral
  INTEGER          :: isw_energy
  INTEGER          :: isw_axisymm
  INTEGER          :: isw_relativistic
  DOUBLE PRECISION :: T_e

  INTEGER          :: isw_momentum
  INTEGER          :: vel_num
  INTEGER          :: vel_distri_swi
  INTEGER          :: nvel
  DOUBLE PRECISION :: vel_max

  DOUBLE PRECISION :: collpar
  DOUBLE PRECISION :: conl_over_mfp
  DOUBLE PRECISION :: coeps
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: y_axi_averages
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: vel_array

  !**********************************************************
  ! Changes required for multispecies support
  !**********************************************************
  LOGICAL          :: lsw_multispecies
  ! isw_coul_log = 0: Coulomb logarithm set as species independent (overrides values for n_spec)
  ! isw_coul_log = 1: Coulomb logarithm computed for each species using n_spec, T_spec
  !                   (overrides values for collisionality parameters)
  INTEGER          :: isw_coul_log
  INTEGER          :: num_spec ! number of species
  INTEGER, DIMENSION(:), ALLOCATABLE :: species_tag
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: conl_over_mfp_spec
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: collpar_spec
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: z_spec ! species charge number
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: m_spec ! species mass [g]
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: n_spec ! species density [cm^-3]
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: T_spec ! species temperature [erg]

  DOUBLE PRECISION :: collpar_min,collpar_max,v_min_resolution,v_max_resolution
  DOUBLE PRECISION :: phi_x_max
  DOUBLE PRECISION :: collop_bspline_dist
  INTEGER          :: collop_bspline_order
  LOGICAL          :: collop_bspline_taylor

  double precision :: m_nbi, T_nbi
  LOGICAL          :: lsw_nbi

END MODULE collisionality_mod
