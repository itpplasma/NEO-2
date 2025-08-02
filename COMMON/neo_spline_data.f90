MODULE neo_spline_data
  ! Splines along s
  use nrtype, only : dp, I4B
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: a_rmnc,b_rmnc,c_rmnc,d_rmnc
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: a_zmnc,b_zmnc,c_zmnc,d_zmnc
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: a_lmnc,b_lmnc,c_lmnc,d_lmnc
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: a_bmnc,b_bmnc,c_bmnc,d_bmnc
  !! Modifications by Andreas F. Martitsch (06.08.2014)
  ! Additional data from Boozer files without Stellarator symmetry
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: a_rmns,b_rmns,c_rmns,d_rmns
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: a_zmns,b_zmns,c_zmns,d_zmns
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: a_lmns,b_lmns,c_lmns,d_lmns
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: a_bmns,b_bmns,c_bmns,d_bmns
  !! End Modifications by Andreas F. Martitsch (06.08.2014)

  REAL(kind=dp), DIMENSION(:),   ALLOCATABLE :: a_iota,b_iota
  REAL(kind=dp), DIMENSION(:),   ALLOCATABLE :: c_iota,d_iota
  REAL(kind=dp), DIMENSION(:),   ALLOCATABLE :: a_pprime,b_pprime
  REAL(kind=dp), DIMENSION(:),   ALLOCATABLE :: c_pprime,d_pprime
  REAL(kind=dp), DIMENSION(:),   ALLOCATABLE :: a_sqrtg00,b_sqrtg00
  REAL(kind=dp), DIMENSION(:),   ALLOCATABLE :: c_sqrtg00,d_sqrtg00
  REAL(kind=dp), DIMENSION(:),   ALLOCATABLE :: a_curr_tor,b_curr_tor
  REAL(kind=dp), DIMENSION(:),   ALLOCATABLE :: c_curr_tor,d_curr_tor
  REAL(kind=dp), DIMENSION(:),   ALLOCATABLE :: a_curr_pol,b_curr_pol
  REAL(kind=dp), DIMENSION(:),   ALLOCATABLE :: c_curr_pol,d_curr_pol

  REAL(kind=dp), DIMENSION(:),   ALLOCATABLE :: r_m, r_mhalf
  INTEGER(I4B),  DIMENSION(:),   ALLOCATABLE :: sp_index

  logical, save :: lsw_linear_boozer

END MODULE neo_spline_data
