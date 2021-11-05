MODULE neo_spline_data
  ! Splines along s
  use nrtype, only : dp, i4b
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: a_rmnc,b_rmnc,c_rmnc,d_rmnc
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: a_zmnc,b_zmnc,c_zmnc,d_zmnc
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: a_lmnc,b_lmnc,c_lmnc,d_lmnc
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: a_bmnc,b_bmnc,c_bmnc,d_bmnc

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

END MODULE neo_spline_data
