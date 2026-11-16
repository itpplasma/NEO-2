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

CONTAINS

  ! Coulomb logarithm of the pair (a,b), isw_coul_log = 1: e-e 14.9, e-i 15.2 (electron n, T),
  ! i-i 17.3 with 1.5 ln T (test-species n, T); n in 1e20 m^-3, T in keV (Wesson, Tokamaks).
  DOUBLE PRECISION FUNCTION coulomb_log(a, b)
    INTEGER, INTENT(in) :: a, b
    DOUBLE PRECISION, PARAMETER :: erg_per_keV = 1.602176634d-9
    INTEGER :: k
    k = a
    IF (z_spec(b) .LT. 0d0) k = b
    IF (z_spec(a) .LT. 0d0 .AND. z_spec(b) .LT. 0d0) THEN
      coulomb_log = 14.9d0 - 0.5d0*LOG(n_spec(k)*1d-14) + LOG(T_spec(k)/erg_per_keV)
    ELSE IF (z_spec(a) .LT. 0d0 .OR. z_spec(b) .LT. 0d0) THEN
      coulomb_log = 15.2d0 - 0.5d0*LOG(n_spec(k)*1d-14) + LOG(T_spec(k)/erg_per_keV)
    ELSE
      coulomb_log = 17.3d0 - 0.5d0*LOG(n_spec(k)*1d-14) + 1.5d0*LOG(T_spec(k)/erg_per_keV)
    END IF
    coulomb_log = coulomb_log - LOG(ABS(z_spec(a)*z_spec(b)))
  END FUNCTION coulomb_log

END MODULE collisionality_mod
