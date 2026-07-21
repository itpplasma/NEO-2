MODULE rkstep_mod
  INTEGER :: legmax,leg,lag

  ! New Version - up to now only differential part
  CHARACTER(len=3), DIMENSION(:), ALLOCATABLE       :: species_tags
  DOUBLE PRECISION, DIMENSION(:,:),     ALLOCATABLE :: Amm
  DOUBLE PRECISION, DIMENSION(:,:),     ALLOCATABLE :: asource
  ! Velocity moments of the single-helicity misalignment drive
  ! (column 1: weight x^-1, ExB piece / A1 channel;
  !  column 2: weight x^2, vpar corrugation piece / A2 channel)
  DOUBLE PRECISION, DIMENSION(:,:),     ALLOCATABLE :: asource_hel
  DOUBLE PRECISION, DIMENSION(:,:,:),   ALLOCATABLE, TARGET :: anumm_a
  DOUBLE PRECISION, DIMENSION(:,:,:),   ALLOCATABLE, TARGET :: denmm_a
  DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE, TARGET :: anumm_aa
  DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE, TARGET :: denmm_aa
  DOUBLE PRECISION, DIMENSION(:,:,:,:,:),ALLOCATABLE, TARGET:: ailmm_aa
  DOUBLE PRECISION, DIMENSION(:,:),     POINTER     :: anumm
  DOUBLE PRECISION, DIMENSION(:,:),     ALLOCATABLE :: anumm_lag
  DOUBLE PRECISION, DIMENSION(:,:),     POINTER     :: denmm
  DOUBLE PRECISION, DIMENSION(:,:,:),   POINTER     :: ailmm
  DOUBLE PRECISION, DIMENSION(:,:),     ALLOCATABLE :: weightlag
  DOUBLE PRECISION, DIMENSION(:),       ALLOCATABLE :: weightden
  DOUBLE PRECISION, DIMENSION(:),       ALLOCATABLE :: weightparflow
  DOUBLE PRECISION, DIMENSION(:),       ALLOCATABLE :: weightenerg
  double precision, dimension(:,:,:),   allocatable :: Inbi_lmmp_a

  DOUBLE PRECISION               :: epserr_sink       ! Regularization
  COMPLEX(kind=kind(1d0))        :: epserr_sink_cmplx ! Regularization
  DOUBLE PRECISION               :: epserr_iter       ! Relative error of integral part iterations
  INTEGER                        :: niter
  DOUBLE PRECISION, DIMENSION(3) :: fluxes
END MODULE rkstep_mod
