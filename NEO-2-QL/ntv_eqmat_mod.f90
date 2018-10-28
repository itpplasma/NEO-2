module ntv_eqmat_mod
  use nrtype, only : dp, dpc

  integer                                      :: nz_symm,nz_asymm,nz_regper
  integer                                      :: nz_per_pos,nz_per_neg
  integer,          dimension(:),  allocatable :: irow_symm,icol_symm
  integer,          dimension(:),  allocatable :: irow_regper,icol_regper
  integer,          dimension(:),  allocatable :: irow_asymm,icol_asymm
  integer,          dimension(:),  allocatable :: irow_per_pos,icol_per_pos
  integer,          dimension(:),  allocatable :: irow_per_neg,icol_per_neg
  real(kind=dp), dimension(:),     allocatable :: amat_symm
  real(kind=dp), dimension(:),     allocatable :: amat_regper
  complex(kind=dpc), dimension(:), allocatable :: amat_asymm
  real(kind=dp), dimension(:,:),   allocatable :: f0_coll,f0_ttmp
  real(kind=dp), dimension(:,:,:), allocatable :: f0_coll_all,f0_ttmp_all
end module ntv_eqmat_mod
