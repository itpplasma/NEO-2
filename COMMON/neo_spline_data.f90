!> Splines along s
module neo_spline_data
  use nrtype, only : dp, i4b

  real(kind=dp), dimension(:,:), allocatable :: a_rmnc,b_rmnc,c_rmnc,d_rmnc
  real(kind=dp), dimension(:,:), allocatable :: a_zmnc,b_zmnc,c_zmnc,d_zmnc
  real(kind=dp), dimension(:,:), allocatable :: a_lmnc,b_lmnc,c_lmnc,d_lmnc
  real(kind=dp), dimension(:,:), allocatable :: a_bmnc,b_bmnc,c_bmnc,d_bmnc
  !! Modifications by Andreas F. Martitsch (06.08.2014)
  ! Additional data from Boozer files without Stellarator symmetry
  real(kind=dp), dimension(:,:), allocatable :: a_rmns,b_rmns,c_rmns,d_rmns
  real(kind=dp), dimension(:,:), allocatable :: a_zmns,b_zmns,c_zmns,d_zmns
  real(kind=dp), dimension(:,:), allocatable :: a_lmns,b_lmns,c_lmns,d_lmns
  real(kind=dp), dimension(:,:), allocatable :: a_bmns,b_bmns,c_bmns,d_bmns
  !! End Modifications by Andreas F. Martitsch (06.08.2014)

  real(kind=dp), dimension(:),   allocatable :: a_iota,b_iota
  real(kind=dp), dimension(:),   allocatable :: c_iota,d_iota
  real(kind=dp), dimension(:),   allocatable :: a_pprime,b_pprime
  real(kind=dp), dimension(:),   allocatable :: c_pprime,d_pprime
  real(kind=dp), dimension(:),   allocatable :: a_sqrtg00,b_sqrtg00
  real(kind=dp), dimension(:),   allocatable :: c_sqrtg00,d_sqrtg00
  real(kind=dp), dimension(:),   allocatable :: a_curr_tor,b_curr_tor
  real(kind=dp), dimension(:),   allocatable :: c_curr_tor,d_curr_tor
  real(kind=dp), dimension(:),   allocatable :: a_curr_pol,b_curr_pol
  real(kind=dp), dimension(:),   allocatable :: c_curr_pol,d_curr_pol

  real(kind=dp), dimension(:),   allocatable :: r_m, r_mhalf
  integer(I4B),  dimension(:),   allocatable :: sp_index

end module neo_spline_data
