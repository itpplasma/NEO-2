module rkstep_mod
  use nrtype, only : dp

  integer :: legmax,leg,lag
  !**********************************************************
  ! Changes required for multispecies support
  !**********************************************************
  ! Old Version
  ! real(kind=dp), dimension(:,:),     allocatable :: anumm,denmm,asource
  ! real(kind=dp), dimension(:,:,:),   allocatable :: ailmm
  ! real(kind=dp), dimension(:,:),     allocatable :: weightlag

  ! New Version - up to now only differential part
  CHARACTER(len=3), dimension(:), allocatable       :: species_tags
  real(kind=dp), dimension(:,:),     allocatable :: Amm
  real(kind=dp), dimension(:,:),     allocatable :: asource
  real(kind=dp), dimension(:,:,:),   allocatable, TARGET :: anumm_a
  real(kind=dp), dimension(:,:,:),   allocatable, TARGET :: denmm_a
  real(kind=dp), dimension(:,:,:,:), allocatable, TARGET :: anumm_aa
  real(kind=dp), dimension(:,:,:,:), allocatable, TARGET :: denmm_aa
  real(kind=dp), dimension(:,:,:,:,:),allocatable, TARGET:: ailmm_aa
  real(kind=dp), dimension(:,:),     POINTER     :: anumm
  real(kind=dp), dimension(:,:),     allocatable :: anumm_lag
  real(kind=dp), dimension(:,:),     POINTER     :: denmm
  real(kind=dp), dimension(:,:,:),   POINTER     :: ailmm
  real(kind=dp), dimension(:,:),     allocatable :: weightlag
  real(kind=dp), dimension(:),       allocatable :: weightden
  real(kind=dp), dimension(:),       allocatable :: weightparflow
  real(kind=dp), dimension(:),       allocatable :: weightenerg
  real(kind=dp), dimension(:,:,:),   allocatable :: Inbi_lmmp_a
  !**********************************************************
  ! WINNY - for flint
  ! real(kind=dp), dimension(:),       allocatable :: collision_sigma_multiplier
  ! WINNY - for flint
  !**********************************************************
  real(kind=dp)               :: epserr_sink       ! Regularization
  DOUBLE COMPLEX                 :: epserr_sink_cmplx ! Regularization
  real(kind=dp)               :: epserr_iter
  integer                        :: niter
  real(kind=dp), dimension(3) :: fluxes
end module
