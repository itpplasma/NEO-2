MODULE neo_input
! Input from data files (Boozer)
  use nrtype, only : dp
  INTEGER, DIMENSION(:),     ALLOCATABLE :: ixm, ixn
  INTEGER, DIMENSION(:),     ALLOCATABLE :: pixm, pixn
  INTEGER, DIMENSION(:),     ALLOCATABLE :: i_m, i_n

  REAL(kind=dp),    DIMENSION(:),     ALLOCATABLE :: es
  REAL(kind=dp),    DIMENSION(:),     ALLOCATABLE :: pprime
  REAL(kind=dp),    DIMENSION(:),     ALLOCATABLE :: sqrtg00

  REAL(kind=dp),    DIMENSION(:,:),   ALLOCATABLE :: rmnc,  zmnc,  lmnc
  REAL(kind=dp),    DIMENSION(:,:),   ALLOCATABLE :: bmnc
  !! Modifications by Andreas F. Martitsch (06.08.2014)
  ! Additional data from Boozer files without Stellarator symmetry
  REAL(kind=dp),    DIMENSION(:,:),   ALLOCATABLE :: rmns,  zmns,  lmns
  REAL(kind=dp),    DIMENSION(:,:),   ALLOCATABLE :: bmns
  !! End Modifications by Andreas F. Martitsch (06.08.2014)
  REAL(kind=dp),    DIMENSION(:),     ALLOCATABLE :: b00

  REAL(kind=dp) :: flux, psi_pr

  INTEGER  :: m0b, n0b
  INTEGER  :: ns, mnmax, nfp
  INTEGER  :: m_max, n_max
END MODULE neo_input
