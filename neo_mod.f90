!
  MODULE size_mod
    INTEGER :: ndim0,npart,ndim,nmat
  END MODULE
!
  MODULE partpa_mod
    ! Winny ipmin,iminr added
    INTEGER                                     :: ipmax,ipmin,istart,iminr
    INTEGER                                     :: npassing_start,npass_max
    DOUBLE PRECISION                            :: pardeb0,bmod0,dirint,hxeta
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: eta 
    ! Winny for talking between 
    INTEGER :: nhalfstep
    ! Winny end
  END MODULE
!
  MODULE rk4_kin_mod
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: y,dydx,yt,dyt,dym
  END MODULE
!

!!$  MODULE ripple_solver_mod
!!$    INTEGER :: npart_loc
!!$    INTEGER :: npass_l,npass_r,mhill
!!$    INTEGER :: npart_halfband
!!$    DOUBLE PRECISION :: qflux,qcurr,qflux_new,qcurr_new
!!$    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE   :: flux_p,flux_m
!!$    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE   :: curr_p,curr_m
!!$    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE   :: source_p,source_m
!!$    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: amat_plus_minus
!!$    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: amat_minus_plus
!!$    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: amat_plus_plus
!!$    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: amat_minus_minus
!!$  END MODULE ripple_solver_mod
!
  MODULE collisionality_mod
    INTEGER          :: isw_lorentz
    INTEGER          :: isw_integral
    INTEGER          :: isw_energy
    INTEGER          :: isw_axisymm

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
  END MODULE
!
  MODULE rkstep_mod
    INTEGER :: legmax,leg,lag
    !**********************************************************
    ! Change for multispecies support
    !**********************************************************
    ! Old Version
    ! DOUBLE PRECISION, DIMENSION(:,:),     ALLOCATABLE :: anumm,denmm,asource
    ! DOUBLE PRECISION, DIMENSION(:,:,:),   ALLOCATABLE :: ailmm
    ! DOUBLE PRECISION, DIMENSION(:,:),     ALLOCATABLE :: weightlag

    ! New Version - up to now only differential part
    CHARACTER(len=3), DIMENSION(:), ALLOCATABLE       :: species_tags
    DOUBLE PRECISION, DIMENSION(:,:),     ALLOCATABLE :: Amm
    DOUBLE PRECISION, DIMENSION(:,:),     ALLOCATABLE :: asource
    DOUBLE PRECISION, DIMENSION(:,:,:),   ALLOCATABLE, TARGET :: anumm_a
    DOUBLE PRECISION, DIMENSION(:,:,:),   ALLOCATABLE, TARGET :: denmm_a
    DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE, TARGET :: anumm_aa
    DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE, TARGET :: denmm_aa
    DOUBLE PRECISION, DIMENSION(:,:,:,:,:),ALLOCATABLE, TARGET:: ailmm_aa
    DOUBLE PRECISION, DIMENSION(:,:),     ALLOCATABLE :: anumm_lag
    DOUBLE PRECISION, DIMENSION(:,:),     POINTER     :: anumm
    DOUBLE PRECISION, DIMENSION(:,:),     POINTER     :: denmm
    DOUBLE PRECISION, DIMENSION(:,:,:),   POINTER     :: ailmm
    DOUBLE PRECISION, DIMENSION(:,:),     ALLOCATABLE :: weightlag
    !**********************************************************
    
    DOUBLE PRECISION, DIMENSION(3) :: fluxes
  END MODULE
!
MODULE development
INTEGER :: iprintflag,proptag_old=-100000000,nstfp
DOUBLE PRECISION :: deleta
DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: alam_l,alam_r
DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: delta_eta_l,delta_eta_r
DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: phi_stfp
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: stfp_ov_hstep
INTEGER :: solver_talk
INTEGER :: switch_off_asymp
! accuracy factor in ripple solver
DOUBLE PRECISION :: ripple_solver_accurfac
! asymptotical stuff
INTEGER :: asymp_margin_zero
INTEGER :: asymp_margin_npass
DOUBLE PRECISION :: asymp_pardeleta
END MODULE

! Winny module for ripple added
MODULE ripple_mod
  INTEGER :: ripple_counter
  INTEGER,          DIMENSION(:), ALLOCATABLE :: imax1_ripple
  INTEGER,          DIMENSION(:), ALLOCATABLE :: imax2_ripple
  INTEGER,          DIMENSION(:), ALLOCATABLE :: imin_ripple
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: col_ripple
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: eta_x0
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: eta_s
  DOUBLE PRECISION                            :: bmax_abs,bmin_abs
  INTEGER :: propagator_counter
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: phi1_prop
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: phi2_prop  
  INTEGER, DIMENSION(:), ALLOCATABLE :: ibeg_prop,iend_prop
  INTEGER, DIMENSION(:), ALLOCATABLE :: ripplenumber_prop
  INTEGER :: ibeg
  INTEGER :: iend
  INTEGER :: ibeg_prop_act
  INTEGER :: iend_prop_act
  INTEGER :: imax_eta
END MODULE ripple_mod

! some Winny stuff about flint
MODULE flint_mod
  INTEGER :: plot_gauss
  INTEGER :: plot_prop
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: phiarr
  INTEGER, DIMENSION(:), ALLOCATABLE :: phi_divide
  DOUBLE PRECISION :: hphi_mult
  INTEGER :: phi_split_mode,phi_place_mode,phi_split_min
  INTEGER :: max_solver_try
  ! WINNY
  DOUBLE PRECISION :: bsfunc_local_err_max_mult
  DOUBLE PRECISION :: bsfunc_max_mult_reach
  DOUBLE PRECISION :: boundary_dist_limit_factor

  INTEGER :: bsfunc_modelfunc_num
  INTEGER :: bsfunc_divide
  INTEGER :: bsfunc_ignore_trap_levels

  DOUBLE PRECISION :: bsfunc_local_shield_factor
  LOGICAL :: bsfunc_shield
END MODULE flint_mod

!!$MODULE join_ripples_simple_mod
!!$  INTEGER :: o_l,o_r,n_l,n_r,n_laguerra
!!$  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: o_amat_p_p,n_amat_p_p
!!$  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: o_amat_p_m,n_amat_p_m
!!$  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: o_amat_m_m,n_amat_m_m
!!$  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: o_amat_m_p,n_amat_m_p
!!$  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: o_source_p,n_source_p
!!$  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: o_source_m,n_source_m
!!$  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: o_flux_p,n_flux_p
!!$  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: o_flux_m,n_flux_m
!!$  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: o_qflux,n_qflux
!!$END MODULE join_ripples_simple_mod
