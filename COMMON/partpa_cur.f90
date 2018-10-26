module partpa_cur
  use nrtype, only : dp
  ! Exchange between flint_cur and rhs_cur
  use sizey_cur

  real(kind=dp)                            :: bmod0
  real(kind=dp)                            :: gamma_cur
  real(kind=dp), dimension(:), allocatable :: y_part, yfac, sqyfac
  real(kind=dp), dimension(:), allocatable :: yfac_xin, yfac_xid
  real(kind=dp), dimension(:), allocatable :: delta_cur
  real(kind=dp), dimension(:), allocatable :: k_fac1, k_fac2
  real(kind=dp), dimension(:), allocatable :: contrif
  real(kind=dp), dimension(:), allocatable :: fcontrif
  integer                                  :: write_curint
end module partpa_cur
