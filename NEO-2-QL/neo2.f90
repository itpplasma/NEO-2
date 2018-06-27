PROGRAM neo2

  !**********************************************************
  ! MPI Support
  !**********************************************************
  USE mpiprovider_module
  USE hdf5_tools

  USE size_mod
  !USE partpa_mod, ONLY : hxeta
  USE flint_mod, ONLY : plot_gauss,plot_prop,phi_split_mode,        &
       phi_place_mode,phi_split_min,hphi_mult,max_solver_try,       &
       bsfunc_local_err_max_mult,bsfunc_max_mult_reach,             &
       bsfunc_modelfunc_num,bsfunc_divide,                          &
       bsfunc_ignore_trap_levels,boundary_dist_limit_factor,        &
       bsfunc_local_shield_factor,bsfunc_shield
  USE device_mod
  USE collisionality_mod, ONLY : conl_over_mfp,isw_lorentz,         &
       isw_integral,isw_energy,isw_axisymm,                         &
       isw_momentum,vel_distri_swi,vel_num,vel_max,                 &
       nvel,vel_array,v_max_resolution,v_min_resolution,            &
       phi_x_max, collop_bspline_order, collop_bspline_dist,        &
       isw_relativistic, T_e, lsw_multispecies, isw_coul_log,       &
       num_spec, species_tag, conl_over_mfp_spec, z_spec, m_spec,   &
       T_spec, n_spec, collop_bspline_taylor, lsw_nbi, m_nbi, T_nbi
  USE propagator_mod, ONLY : reconstruct_prop_dist,                 &
       prop_diagphys,prop_overwrite,                                &
       prop_diagnostic,prop_binary,                                 &
       prop_timing,prop_join_ends,prop_fluxsplitmode,               &
       prop_write,prop_reconstruct,prop_ripple_plot,                &
       prop_reconstruct_levels
  USE magnetics_mod, ONLY : mag_talk,mag_infotalk
  USE mag_interface_mod, ONLY : mag_local_sigma, hphi_lim,          &
       mag_magfield,mag_nperiod_min,mag_save_memory,                &
       magnetic_device,mag_cycle_ripples,mag_start_special,         &
       aiota_tokamak,mag_close_fieldline,mag_ripple_contribution,   &
       mag_coordinates,boozer_s,boozer_theta_beg,boozer_phi_beg,    &
       mag_dbhat_min,mag_dphi_inf_min,mag_inflection_mult,          &
       mag_symmetric,mag_symmetric_shorten,                         &
       sigma_shield_factor,split_inflection_points
  USE binarysplit_mod, ONLY : bsfunc_message,bsfunc_modelfunc,      &
       bsfunc_total_err, bsfunc_local_err, bsfunc_min_distance,     &
       bsfunc_max_index, bsfunc_max_splitlevel,                     &
       bsfunc_sigma_mult, bsfunc_sigma_min, bsfunc_local_solver   
  USE binarysplit_int, ONLY : linspace
  !! Modifications by Andreas F. Martitsch (15.07.2014)
  ! Path for the collision operator matrices is now specified via neo2.in
  ! (necessary for computations with Condor)  
  USE collop, ONLY : collop_construct, collop_deconstruct,          &
       collop_load, collop_unload, z_eff, collop_path,              &
  !! End Modifications by Andreas F. Martitsch (15.07.2014)
       collop_base_prj, collop_base_exp, scalprod_alpha,            &
       scalprod_beta, lsw_read_precom, lsw_write_precom !! Added lsw_read_precom
       !! and lsw_write_precom by Michael Draxler (25.08.2017)
  USE rkstep_mod, ONLY : lag,leg,legmax                            
      
  USE development, ONLY : solver_talk,switch_off_asymp, &
       asymp_margin_zero, asymp_margin_npass, asymp_pardeleta,      &
       ripple_solver_accurfac
  USE sparse_mod, ONLY : sparse_talk,sparse_solve_method,sparse_example
  !! Modification by Andreas F. Martitsch (14.07.2015)
  ! Extra input for NTV computations
  USE ntv_mod, ONLY : isw_ntv_mode, isw_qflux_NA, in_file_pert,     &
       MtOvR, B_rho_L_loc, xstart_cyl, isw_ripple_solver,           &
       isw_calc_Er, isw_calc_MagDrift, species_tag_Vphi,            &
       isw_Vphi_loc, Vphi, R_Vphi, Z_Vphi, boozer_theta_Vphi,       &
       dn_spec_ov_ds, dT_spec_ov_ds
  !! End Modification by Andreas F. Martitsch (14.07.2015)
  !! Modifications by Andreas F. Martitsch (17.03.2016)
  ! derivative of iota for non-local NTV computations
  ! (with magnetic shear)
  USE neo_magfie_mod, ONLY : isw_mag_shear
  !! End Modifications by Andreas F. Martitsch (17.03.2016)
  USE neo_sub_mod, ONLY : neo_read_control ! only used for preparation of multi-spec input
  USE neo_control, ONLY: in_file, inp_swi, lab_swi

  !************************************
  ! HDF5
  !************************************
  USE hdf5_tools
  USE hdf5_tools_f2003
  !
  IMPLICIT NONE

  INTEGER, PARAMETER :: dp = KIND(1.0d0)
  INTEGER, PARAMETER :: MAXDIM = 1000

  LOGICAL :: opened

  !**********************************************************
  ! Include version information
  !**********************************************************
  INCLUDE "version.f90"
  !**********************************************************

  
  
  REAL(kind=dp), PARAMETER :: pi=3.14159265358979_dp

  !**********************************************************
  ! Include version information
  !**********************************************************
  !INCLUDE "cmake_version.f90"
  !INCLUDE "version.f90"

  !************************************************
  ! HDF5
  !************************************************
  INTEGER(HID_T) :: h5_config_id
  INTEGER(HID_T) :: h5_config_group
  CHARACTER(8)  :: date
  CHARACTER(10) :: time
  CHARACTER(50) :: datetimestring
  REAL(kind=dp) :: rand_num
  CHARACTER(6)  :: rand_hash
  CHARACTER(32) :: fieldname

  INTEGER(HID_T)  :: h5id_taginfo, h5id_propfile, h5id_final, h5id_surf, h5id_neo2
  INTEGER(HID_T)  :: h5id_propagators, h5id_prop
  CHARACTER(512)  :: surfname
  CHARACTER(512)  :: h5_filename
  CHARACTER(1024) :: cwd
  INTEGER         :: k,l
  INTEGER         :: tag_first, tag_last
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: cg0_1_num_prop, cg2_1_num_prop
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: cg0_2_num_prop, cg2_2_num_prop
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: cg0_3_num_prop, cg2_3_num_prop, denom_mflint_prop
  REAL(kind=dp)   :: cg0_1_avg, cg2_1_avg
  REAL(kind=dp)   :: cg0_2_avg, cg2_2_avg
  REAL(kind=dp)   :: cg0_3_avg, cg2_3_avg
  !**********************************************************
  
  REAL(kind=dp) :: rbeg,zbeg
  REAL(kind=dp) :: phimi
  REAL(kind=dp) :: xetama,xetami

  REAL(kind=dp) :: eta_s_lim
  ! REAL(kind=dp) :: z_eff

  INTEGER :: proptag_first,proptag_last,proptag_start,proptag_end
  INTEGER :: proptag_begin,proptag_final
  INTEGER :: uw
  INTEGER :: ialloc
  INTEGER :: nstep,nperiod
  INTEGER :: eta_part,lambda_equi
  INTEGER :: bin_split_mode
  INTEGER :: eta_part_global,eta_part_trapped
  REAL(kind=dp) :: eta_part_globalfac,eta_part_globalfac_p,eta_part_globalfac_t
  REAL(kind=dp) :: eta_alpha_p,eta_alpha_t
  ! ---------------------------------------------------------------------------
  !
  ! the input is read from 2 files
  !   neo2.def   with default values
  !   neo2.in    overwrites defaultvalues if present
  !
  !  they are both in namelist-format (rather free and convenient format)
  !
  ! settings for namelist
  INTEGER :: u1=10
  INTEGER :: ios
  INTEGER :: jf
  CHARACTER(len=20), DIMENSION(2) :: fnames
  !! Modification by Andreas F. Martitsch (23.08.2015)
  ! multi-species part:
  ! -> prepare multi-species computation
  INTEGER :: ind_spec, ind_boozer_s, ind_char, ctr_spec
  INTEGER :: isw_multispecies_init, OMP_NUM_THREADS
  INTEGER(HID_T) :: h5id_multispec_in
  CHARACTER(len=40) :: fname_multispec, fname_multispec_in
  CHARACTER(len=40) :: fname_exec, fname_exec_precom
  CHARACTER(len=200) :: dir_name, cmd_line
  INTEGER :: num_radial_pts, num_species_all
  INTEGER, DIMENSION(:), ALLOCATABLE :: rel_stages_prof, species_tag_prof
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: boozer_s_prof, Vphi_prof
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: R_Vphi_prof, Z_Vphi_prof
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: boozer_theta_Vphi_prof
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: T_prof, n_prof
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: dT_ov_ds_prof, dn_ov_ds_prof
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: kappa_prof
  REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: species_def_prof
  ! -> read species-relevant info into a large array (dynamic allocation not supported)
  INTEGER, DIMENSION(:), ALLOCATABLE :: species_tag_vec   ! vector with species tags
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: conl_over_mfp_vec ! collisionality parameter
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: z_vec ! species charge number 
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: m_vec ! species mass
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: T_vec ! species temperature
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: n_vec ! species density (used only for isw_coul_log > 0)
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: dT_vec_ov_ds
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: dn_vec_ov_ds
  !! End Modification by Andreas F. Martitsch (23.08.2015)
  ! groups for namelist
  !! Modification by Andreas F. Martitsch (21.02.2017)
  ! multi-species part:
  NAMELIST /multi_spec/                                                       &
       lsw_multispecies, isw_multispecies_init, fname_multispec_in,           &
       isw_coul_log, num_spec, species_tag_vec, conl_over_mfp_vec, z_vec,     &
       m_vec, T_vec, n_vec, isw_calc_Er, isw_calc_MagDrift,                   &
       species_tag_Vphi, isw_Vphi_loc, Vphi, R_Vphi, Z_Vphi,                  &
       boozer_theta_Vphi, dn_vec_ov_ds, dT_vec_ov_ds
  !! End Modification by Andreas F. Martitsch (21.02.2017)
  NAMELIST /settings/                                                         &
       phimi,nstep,nperiod,xetami,xetama,ndim0,zbeg,rbeg,                     &
       proptag_begin,proptag_final,mag_start_special,                         &
       mag_magfield,mag_nperiod_min,                                          &
       mag_save_memory,magnetic_device,mag_cycle_ripples,                     &
       aiota_tokamak,mag_close_fieldline,eta_part_global,eta_part_globalfac,  &
       eta_part_globalfac_p,eta_part_globalfac_t,                             &
       eta_alpha_p,eta_alpha_t,eta_part_trapped,                              &
       mag_coordinates,boozer_s,boozer_theta_beg,boozer_phi_beg,              &
       mag_dbhat_min,mag_dphi_inf_min,mag_inflection_mult,                    & 
       solver_talk,switch_off_asymp,                                          &
       asymp_margin_zero,asymp_margin_npass,asymp_pardeleta,                  &
       ripple_solver_accurfac,                                                &
       sparse_talk,sparse_solve_method, OMP_NUM_THREADS,                      &
       mag_symmetric,mag_symmetric_shorten
  NAMELIST /collision/                                                        &
       conl_over_mfp,lag,leg,legmax,z_eff,isw_lorentz,                        &
       isw_integral,isw_energy,isw_axisymm,                                   &
       isw_momentum,vel_distri_swi,vel_num,vel_max,collop_path,               &      
       collop_base_prj, collop_base_exp, scalprod_alpha, scalprod_beta,       &
       phi_x_max, collop_bspline_order, collop_bspline_dist,                  &
       v_min_resolution, v_max_resolution, isw_relativistic, T_e,             &
       lsw_read_precom, lsw_write_precom, collop_bspline_taylor, lsw_nbi,     &
       m_nbi, T_nbi
  NAMELIST /binsplit/                                                         &
       eta_s_lim,eta_part,lambda_equi,phi_split_mode,phi_place_mode,          &
       phi_split_min,max_solver_try,                                          &
       hphi_mult,bin_split_mode,bsfunc_message,bsfunc_modelfunc,              &
       bsfunc_total_err,bsfunc_local_err,bsfunc_min_distance,                 &
       bsfunc_max_index,bsfunc_max_splitlevel,                                &
       bsfunc_sigma_mult, bsfunc_sigma_min, bsfunc_local_solver,              &
       mag_local_sigma,mag_ripple_contribution,                               &
       bsfunc_local_err_max_mult,bsfunc_max_mult_reach,                       &
       bsfunc_modelfunc_num,bsfunc_divide,                                    &
       bsfunc_ignore_trap_levels,boundary_dist_limit_factor,                  &
       bsfunc_local_shield_factor,bsfunc_shield,sigma_shield_factor,          &
       split_inflection_points
  NAMELIST /propagator/                                                       &
       prop_diagphys,prop_overwrite,                                          &
       prop_diagnostic,prop_binary,prop_timing,prop_join_ends,                &
       prop_fluxsplitmode,                                                    &
       mag_talk,mag_infotalk,                                                 &
       hphi_lim,                                                              &
       prop_write,prop_reconstruct,prop_ripple_plot,                          &
       prop_reconstruct_levels
  NAMELIST /plotting/                                                         &
       plot_gauss,plot_prop
  !! Modification by Andreas F. Martitsch (14.07.2015)
  ! Extra input for NTV computation
  NAMELIST /ntv_input/                                                        &
       isw_ntv_mode, isw_qflux_NA, in_file_pert, MtOvR, B_rho_L_loc,          &
       isw_ripple_solver, isw_mag_shear
  !! End Modification by Andreas F. Martitsch (14.07.2015)
  
  ! ---------------------------------------------------------------------------
  ! filenames (default file and specific input file) for namelist
  fnames = (/'neo2.def','neo2.in '/)
  ! file-names of multi-species input and startup-script for NEO-2
  fname_multispec = 'neo2.in'
  fname_exec = 'run_neo2.sh'
  fname_exec_precom = 'run_neo2_precom.sh'
  ! ---------------------------------------------------------------------------
  ! defaults
  !
  ! settings
    !! Modification by Andreas F. Martitsch (21.02.2017)
  ! multi-species part:
  lsw_multispecies = .FALSE.
  isw_multispecies_init = 0
  fname_multispec_in = ''
  ! isw_coul_log = 0: Coulomb logarithm set as species independent (overrides values for n_spec)
  ! isw_coul_log = 1: Coulomb logarithm computed for each species using n_spec, T_spec
  !                   (overrides values for collisionality parameters)
  isw_coul_log = 0
  num_spec = 1
  IF(ALLOCATED(species_tag_vec)) DEALLOCATE(species_tag_vec)
  ALLOCATE(species_tag_vec(MAXDIM))
  species_tag_vec = (/ (ind_spec,ind_spec=0,MAXDIM-1) /)
  IF(ALLOCATED(conl_over_mfp_vec)) DEALLOCATE(conl_over_mfp_vec)
  ALLOCATE(conl_over_mfp_vec(MAXDIM))
  conl_over_mfp_vec = 1.0d-3
  IF(ALLOCATED(z_vec)) DEALLOCATE(z_vec)
  ALLOCATE(z_vec(MAXDIM))
  z_vec = 1.0d0
  IF(ALLOCATED(m_vec)) DEALLOCATE(m_vec)
  ALLOCATE(m_vec(MAXDIM))
  m_vec = 1.672621637d-24 ! proton mass [g]
  IF(ALLOCATED(T_vec)) DEALLOCATE(T_vec)
  ALLOCATE(T_vec(MAXDIM))
  T_vec = 1.6d-9 ! temperature [erg] (=1keV)
  IF(ALLOCATED(n_vec)) DEALLOCATE(n_vec)
  ALLOCATE(n_vec(MAXDIM))
  n_vec = 1.0d13 ! density  [cm^-3]
  isw_calc_Er = 0
  isw_calc_MagDrift = 0
  species_tag_Vphi = 0 ! species index should be the same for all flux surfaces
  isw_Vphi_loc = 0 ! 0: <V_\varphi>, 1: V_{\varphi}(R,Z)
  Vphi = 0.0d0
  R_Vphi = 0.0d0 ! only used for isw_Vphi_loc=1
  Z_Vphi = 0.0d0 ! only used for isw_Vphi_loc=1
  boozer_theta_Vphi = 0.0d0 ! only used for isw_Vphi_loc=2
  IF(ALLOCATED(dT_vec_ov_ds)) DEALLOCATE(dT_vec_ov_ds)
  ALLOCATE(dT_vec_ov_ds(MAXDIM))
  dT_vec_ov_ds = 0.0d0 ! temperature [erg] (=1keV)
  IF(ALLOCATED(dn_vec_ov_ds)) DEALLOCATE(dn_vec_ov_ds)
  ALLOCATE(dn_vec_ov_ds(MAXDIM))
  dn_vec_ov_ds = 0.0d0 ! density  [cm^-3]
  !! End Modification by Andreas F. Martitsch (21.02.2017)
  mag_magfield = 1
  magnetic_device = 1
  mag_nperiod_min = 300
  mag_save_memory = 1
  mag_cycle_ripples = 0
  mag_close_fieldline = 1
  mag_ripple_contribution = 1
  mag_dbhat_min = 0.1d0
  mag_dphi_inf_min = 0.05d0
  mag_inflection_mult = 3.0d0
  solver_talk = 0
  switch_off_asymp = 0
  asymp_margin_zero = 10
  asymp_margin_npass = 4
  asymp_pardeleta = 10.0d0
  ripple_solver_accurfac = 3.0d0
  phimi=0.d0 
  nstep=480
  nperiod=500
  xetami=0.0d0
  xetama=1.300001d0
  eta_part_global = 0
  eta_part_trapped = 10
  eta_part_globalfac = 3.0_dp
  eta_part_globalfac_p = 3.0_dp
  eta_part_globalfac_t = 3.0_dp
  eta_alpha_p = 4.0_dp
  eta_alpha_t = 1.0_dp
  ndim0=14
  zbeg=0.d0
  rbeg=181.d0
  proptag_begin=0
  proptag_final=0
  mag_start_special=0
  aiota_tokamak = 1.0d0/3.0d0
  mag_coordinates = 0
  boozer_s = 0.5_dp
  boozer_theta_beg = 0.0_dp
  boozer_phi_beg = 0.0_dp
  sparse_talk = .FALSE.
  !  sparse_solve_method = 0
  OMP_NUM_THREADS = 1
  ! collision 
  conl_over_mfp = 1.0d-3
  lag=10
  leg=20
  legmax=20
  z_eff=1.d0
  isw_lorentz = 1
  isw_integral = 0
  isw_energy = 0
  isw_axisymm = 0
  isw_momentum = 0
  isw_relativistic = 0  ! 0 is non-relativistic, 1 is Braams and Karney, 2 is high order Legendre
  T_e = 1 ! Default 1eV
  vel_distri_swi = 0
  vel_num = 10
  vel_max = 5.0d0
  !! Modifications by Andreas F. Martitsch (15.07.2014)
  ! Default path for the collision operator matrices
  collop_path = '/afs/itp.tugraz.at/proj/plasma/DOCUMENTS/Neo2/data-MatrixElements/'
  !! End Modifications by Andreas F. Martitsch (15.07.2014)
  collop_base_prj = 0
  collop_base_exp = 0
  scalprod_alpha = 0d0
  scalprod_beta  = 0d0
  v_min_resolution = 0.1d0
  v_max_resolution = 5.0d0
  phi_x_max        = 5.0d0
  collop_bspline_order = 4
  collop_bspline_dist  = 1d0
  collop_bspline_taylor = .true.
  lsw_read_precom = .FALSE. !! Added lsw_read_precom and lsw_write_precom 
  lsw_write_precom = .FALSE.   !! by Michael Draxler (25.08.2017)
  ! binsplit
  eta_s_lim = 1.2d0
  eta_part = 100
  lambda_equi = 0
  phi_split_mode = 2
  phi_place_mode = 2
  phi_split_min = 1
  max_solver_try = 1
  hphi_mult = 1.0d0 
  bin_split_mode = 1
  bsfunc_message = 0
  bsfunc_modelfunc = 1
  bsfunc_modelfunc_num = 1
  bsfunc_ignore_trap_levels = 0
  boundary_dist_limit_factor = 1.e-2
  bsfunc_local_shield_factor = 1.0d0
  bsfunc_shield = .FALSE.
  bsfunc_divide = 0  
  bsfunc_total_err = 1.0d-1
  bsfunc_local_err = 1.0d-2
  bsfunc_local_err_max_mult = 1.0d0
  bsfunc_max_mult_reach = 3.0d0
  bsfunc_min_distance = 0.0d0 
  bsfunc_max_index = 20*eta_part
  bsfunc_max_splitlevel = 32    
  bsfunc_sigma_mult = 1.0_dp
  bsfunc_sigma_min = 0.0_dp
  bsfunc_local_solver = 0
  sigma_shield_factor = 3.0d0
  split_inflection_points = .TRUE.
  mag_local_sigma = 0
  mag_symmetric = .FALSE.
  mag_symmetric_shorten = .FALSE.
  ! propagator
  prop_diagphys = 1
  prop_overwrite   = 1
  prop_diagnostic = 1
  prop_binary = 0
  prop_timing = 1 
  prop_join_ends = 0
  prop_fluxsplitmode = 1
  prop_write = 0
  prop_reconstruct = 0
  prop_ripple_plot = 0
  prop_reconstruct_levels = 0
  mag_talk = .TRUE. 
  mag_infotalk = .TRUE.
  hphi_lim = 1.0d-6
  ! plotting
  plot_gauss = 0 
  plot_prop  = 0
  !! Modification by Andreas F. Martitsch (14.07.2015)
  ! ntv_input
  isw_ntv_mode = 0
  isw_qflux_NA = 0
  MtOvR = 0.0d0
  B_rho_L_loc = 0.0d0
  isw_ripple_solver = 1
  isw_mag_shear = 0
  !! End Modification by Andreas F. Martitsch (14.07.2015)

  lsw_nbi = .false.
  T_nbi  = 70d3
  m_nbi  = 3.343583719d-24
  
  CALL h5_init()
  
  ! reading
  DO jf = 1,SIZE(fnames)
     IF(jf .EQ. 1) CYCLE ! skip neo2.def (Andreas F. Martitsch - 21.10.2015)
     OPEN(unit=u1,file=fnames(jf),status='old',iostat=ios)
     if (ios .ne. 0) then
        PRINT *, 'WARNING: File ',fnames(jf),' cannot be OPENED!'
        PRINT *, ''
        STOP
     ELSE
        ! Read variables from group settings
        !! Modification by Andreas F. Martitsch (21.02.2017)
        ! multi-species part:
        READ(u1,nml=multi_spec,iostat=ios)
        IF (ios .NE. 0) THEN
           PRINT *, 'WARNING: group multi_spec in ',fnames(jf),' cannot be READ!'
           PRINT *, ''
           !STOP
        END IF
        !! End Modification by Andreas F. Martitsch (21.02.2017)
        REWIND(u1) ! start reading file from beginning (Andreas F. Martitsch - 23.02.2017)
        READ(u1,nml=settings,iostat=ios)
        IF (ios .NE. 0) THEN
           PRINT *, 'WARNING: group settings in ',fnames(jf),' cannot be READ!'
           PRINT *, ''
           STOP
        END IF
        REWIND(u1) ! start reading file from beginning (Andreas F. Martitsch - 23.02.2017)
        READ(u1,nml=collision,iostat=ios)
        IF (ios .NE. 0) THEN
           PRINT *, 'WARNING: group collision in ',fnames(jf),' cannot be READ!'
           PRINT *, ''
           STOP
        END IF
        REWIND(u1) ! start reading file from beginning (Andreas F. Martitsch - 23.02.2017)
        READ(u1,nml=binsplit,iostat=ios)
        IF (ios .NE. 0) THEN
           PRINT *, 'WARNING: group binsplit in ',fnames(jf),' cannot be READ!'
           PRINT *, ''
           STOP
        END IF
        REWIND(u1) ! start reading file from beginning (Andreas F. Martitsch - 23.02.2017)
        READ(u1,nml=propagator,iostat=ios)
        IF (ios .NE. 0) THEN
           PRINT *, 'WARNING: group propagator in ',fnames(jf),' cannot be READ!'
           PRINT *, ''
           STOP
        END IF
        REWIND(u1) ! start reading file from beginning (Andreas F. Martitsch - 23.02.2017)
        READ(u1,nml=plotting,iostat=ios)
        IF (ios .NE. 0) THEN
           PRINT *, 'WARNING: group plotting in ',fnames(jf),' cannot be READ!'
           PRINT *, ''
           STOP
        END IF
        !! Modification by Andreas F. Martitsch (17.07.2014)
        ! ntv_input
        REWIND(u1) ! start reading file from beginning (Andreas F. Martitsch - 23.02.2017)
        READ(u1,nml=ntv_input,iostat=ios)
        IF (ios .NE. 0) THEN
           PRINT *, 'WARNING: group ntv_input in ',fnames(jf),' cannot be READ!'
           PRINT *, ''
           STOP
        END IF
        !! End Modification by Andreas F. Martitsch (17.07.2014)
     END IF
     CLOSE(unit=u1)
  END DO
  ! PAUSE
  
  !! Modification by Andreas F. Martitsch (20.02.2017)
  ! Prepare  multi-species computations for a given profile
  ! -> prepare input files, directories
  IF(lsw_multispecies .AND. (isw_multispecies_init .GT. 0)) THEN
     PRINT *,'Prepare  multi-species computations for a given profile'
     isw_multispecies_init = 0
     !
     ! get file-name of axisymmetric equilibrium (neo.in)
     ! -> in_file
     CALL neo_read_control()
     !
     ! read multi-species input (HDF5 file)
     !
     ! open file
     CALL h5_open(TRIM(ADJUSTL(fname_multispec_in)), h5id_multispec_in)
     ! get size of arrays
     CALL h5_get(h5id_multispec_in,'num_radial_pts',num_radial_pts)
     CALL h5_get(h5id_multispec_in,'num_species',num_species_all)
     ! get boozer_s-profile
     IF(ALLOCATED(boozer_s_prof)) DEALLOCATE(boozer_s_prof)
     ALLOCATE(boozer_s_prof(num_radial_pts))
     CALL h5_get(h5id_multispec_in,'boozer_s',boozer_s_prof)
     !PRINT *,boozer_s_prof
     ! get species definition (species charge number and mass)
     IF(ALLOCATED(species_def_prof)) DEALLOCATE(species_def_prof)
     ALLOCATE(species_def_prof(num_radial_pts,num_species_all,2))
     CALL h5_get(h5id_multispec_in,'species_def',species_def_prof)
     !PRINT *,species_def_prof(:,1)
     !PRINT *,species_def_prof(:,2)
     ! get species tag and number of relevant ionization stages
     IF(ALLOCATED(species_tag_prof)) DEALLOCATE(species_tag_prof)
     ALLOCATE(species_tag_prof(num_species_all))
     CALL h5_get(h5id_multispec_in,'species_tag',species_tag_prof)
     !PRINT *,species_tag_prof
     IF(ALLOCATED(rel_stages_prof)) DEALLOCATE(rel_stages_prof)
     ALLOCATE(rel_stages_prof(num_radial_pts))
     CALL h5_get(h5id_multispec_in,'rel_stages',rel_stages_prof)
     !PRINT *,rel_stages_prof
     ! get density and temperature profiles
     IF(ALLOCATED(n_prof)) DEALLOCATE(n_prof)
     ALLOCATE(n_prof(num_radial_pts,num_species_all))
     CALL h5_get(h5id_multispec_in,'n_prof',n_prof)
     !DO ind_boozer_s=1,num_radial_pts
     !   WRITE(*,*) (n_prof(ind_boozer_s,ind_spec),ind_spec=1,num_species_all)
     !END DO
     IF(ALLOCATED(T_prof)) DEALLOCATE(T_prof)
     ALLOCATE(T_prof(num_radial_pts,num_species_all))
     CALL h5_get(h5id_multispec_in,'T_prof',T_prof)
     !DO ind_boozer_s=1,num_radial_pts
     !   WRITE(*,*) (T_prof(ind_boozer_s,ind_spec),ind_spec=1,num_species_all)
     !END DO
     ! get gradients of density and temperature profiles
     IF(ALLOCATED(dn_ov_ds_prof)) DEALLOCATE(dn_ov_ds_prof)
     ALLOCATE(dn_ov_ds_prof(num_radial_pts,num_species_all))
     CALL h5_get(h5id_multispec_in,'dn_ov_ds_prof',dn_ov_ds_prof)
     !DO ind_boozer_s=1,num_radial_pts
     !   WRITE(*,*) (dn_ov_ds_prof(ind_boozer_s,ind_spec),ind_spec=1,num_species_all)
     !END DO
     IF(ALLOCATED(dT_ov_ds_prof)) DEALLOCATE(dT_ov_ds_prof)
     ALLOCATE(dT_ov_ds_prof(num_radial_pts,num_species_all))
     CALL h5_get(h5id_multispec_in,'dT_ov_ds_prof',dT_ov_ds_prof)
     !DO ind_boozer_s=1,num_radial_pts
     !   WRITE(*,*) (dT_ov_ds_prof(ind_boozer_s,ind_spec),ind_spec=1,num_species_all)
     !END DO
     ! get collisionality profile
     IF(ALLOCATED(kappa_prof)) DEALLOCATE(kappa_prof)
     ALLOCATE(kappa_prof(num_radial_pts,num_species_all))
     CALL h5_get(h5id_multispec_in,'kappa_prof',kappa_prof)
     !DO ind_boozer_s=1,num_radial_pts
     !   WRITE(*,*) (kappa_prof(ind_boozer_s,ind_spec),ind_spec=1,num_species_all)
     !END DO
     ! get measured toroidal rotation profile and its species-tag
     IF(ALLOCATED(Vphi_prof)) DEALLOCATE(Vphi_prof)
     ALLOCATE(Vphi_prof(num_radial_pts))
     CALL h5_get(h5id_multispec_in,'Vphi',Vphi_prof)
     CALL h5_get(h5id_multispec_in,'species_tag_Vphi',species_tag_Vphi)
     CALL h5_get(h5id_multispec_in,'isw_Vphi_loc',isw_Vphi_loc)
     IF (isw_Vphi_loc .EQ. 1) THEN
        PRINT *,"neo2.f90: Warning switch isw_Vphi_loc=1 is not tested!"
        STOP
        !
        IF(ALLOCATED(R_Vphi_prof)) DEALLOCATE(R_Vphi_prof)
        ALLOCATE(R_Vphi_prof(num_radial_pts))
        CALL h5_get(h5id_multispec_in,'R_Vphi',R_Vphi_prof)
        IF(ALLOCATED(Z_Vphi_prof)) DEALLOCATE(Z_Vphi_prof)
        ALLOCATE(Z_Vphi_prof(num_radial_pts))
        CALL h5_get(h5id_multispec_in,'Z_Vphi',Z_Vphi_prof)
     ELSE IF(isw_Vphi_loc .EQ. 2) THEN
        IF(ALLOCATED(boozer_theta_Vphi_prof)) DEALLOCATE(boozer_theta_Vphi_prof)
        ALLOCATE(boozer_theta_Vphi_prof(num_radial_pts))
        CALL h5_get(h5id_multispec_in,'boozer_theta_Vphi',boozer_theta_Vphi_prof)
     ELSE IF (isw_Vphi_loc.LT.0 .OR. isw_Vphi_loc.GT.2) THEN
        PRINT *,"neo2.f90: Undefined state of switch isw_Vphi_loc (= 0 / 1 / 2)!"
        STOP
     END IF
     ! close file
     CALL h5_close(h5id_multispec_in)
     !
     ! prepare directories for NEO-2 runs
     !
     DO ind_boozer_s = 1,num_radial_pts
        !
        ! get number of relevant species for each radial point
        num_spec = rel_stages_prof(ind_boozer_s)
        IF (num_spec .EQ. 0) CYCLE ! try next radial point
        !
        ! directory name
        WRITE(dir_name,fmt='(F10.5)') boozer_s_prof(ind_boozer_s)
        ind_char=INDEX(dir_name,'.')
        dir_name=dir_name(:ind_char-1) // 'p' // dir_name(ind_char+1:)
        dir_name = 'es_' // TRIM(ADJUSTL(dir_name))
        !PRINT *,dir_name
        !
        ! shell command: create directories
        cmd_line = &
             'if [ ! -d ' // TRIM(ADJUSTL(dir_name)) // ' ]; then mkdir ' // &
             TRIM(ADJUSTL(dir_name)) // '; fi'
        CALL execute_command_LINE(cmd_line)
        !
        ! go to directory
        CALL chdir(TRIM(ADJUSTL(dir_name)))
        !
        ! shell command: link input-files to current directory
        cmd_line = &
             'if [ ! -e ' // TRIM(ADJUSTL(in_file)) // ' ]; then ln -s ../' // &
             TRIM(ADJUSTL(in_file)) // ' . ; fi'
        CALL execute_command_LINE(cmd_line)
        cmd_line = &
             'if [ ! -e ' // TRIM(ADJUSTL(in_file_pert)) // ' ]; then ln -s ../' // &
             TRIM(ADJUSTL(in_file_pert)) // ' . ; fi'
        CALL execute_command_LINE(cmd_line)
        cmd_line = &
             'if [ ! -e neo.in ]; then ln -s ../neo.in . ; fi'
        CALL execute_command_LINE(cmd_line)
        cmd_line = &
             'if [ ! -e neo_2.x ]; then ln -s ../neo_2.x . ; fi'
        CALL execute_command_LINE(cmd_line)
        !
        ! write start-up script for NEO-2 run
        OPEN(unit=u1,file=fname_exec,action='write',iostat=ios)
        IF (ios .NE. 0) THEN
           PRINT *, 'WARNING: File ',fname_exec,' cannot be OPENED!'
           PRINT *, ''
           STOP
        ELSE
           WRITE(u1,fmt='(A)',iostat=ios) '#! /bin/bash'
           IF (ios .NE. 0) THEN
              PRINT *, 'WARNING: ',fname_exec,' cannot be WRITTEN!'
              PRINT *, ''
              STOP
           END IF
           WRITE(cmd_line,fmt='(I2)') OMP_NUM_THREADS
           cmd_line = 'OMP_NUM_THREADS=' // TRIM(ADJUSTL(cmd_line)) // &
                ' mpiexec -mca orte_tmpdir_base "/tmp/" -x OMP_NUM_THREADS -hostfile hosts -np '
           WRITE(cmd_line,fmt='(A,I3,A10)') TRIM(ADJUSTL(cmd_line)),num_spec,' ./neo_2.x'
           WRITE(u1,fmt='(A)',iostat=ios) TRIM(ADJUSTL(cmd_line))
           IF (ios .NE. 0) THEN
              PRINT *, 'WARNING: ',fname_exec,' cannot be WRITTEN!'
              PRINT *, ''
              STOP
           END IF
        END IF
        CLOSE(unit=u1)
        cmd_line = 'chmod u+x ' // TRIM(ADJUSTL(fname_exec))
        CALL execute_command_LINE(cmd_line)
        !
        ! write start-up script for NEO-2 pre-run (pre-computation of matrix elements)
        OPEN(unit=u1,file=fname_exec_precom,action='write',iostat=ios)
        IF (ios .NE. 0) THEN
           PRINT *, 'WARNING: File ',fname_exec,' cannot be OPENED!'
           PRINT *, ''
           STOP
        ELSE
           WRITE(u1,fmt='(A)',iostat=ios) '#! /bin/bash'
           IF (ios .NE. 0) THEN
              PRINT *, 'WARNING: ',fname_exec,' cannot be WRITTEN!'
              PRINT *, ''
              STOP
           END IF
           WRITE(cmd_line,fmt='(I2)') OMP_NUM_THREADS
           cmd_line = 'OMP_NUM_THREADS=' // TRIM(ADJUSTL(cmd_line)) // &
                ' mpiexec -mca orte_tmpdir_base "/tmp/" -x OMP_NUM_THREADS -np '
           WRITE(cmd_line,fmt='(A,I3,A10)') TRIM(ADJUSTL(cmd_line)),num_spec,' ./neo_2.x'
           WRITE(u1,fmt='(A)',iostat=ios) TRIM(ADJUSTL(cmd_line))
           IF (ios .NE. 0) THEN
              PRINT *, 'WARNING: ',fname_exec,' cannot be WRITTEN!'
              PRINT *, ''
              STOP
           END IF
        END IF
        CLOSE(unit=u1)
        cmd_line = 'chmod u+x ' // TRIM(ADJUSTL(fname_exec_precom))
        CALL execute_command_LINE(cmd_line)
        !
        ! prepare multi-species input
        !
        ! specify radial point
        boozer_s = boozer_s_prof(ind_boozer_s)
        Vphi = Vphi_prof(ind_boozer_s)
        IF (isw_Vphi_loc .EQ. 1) THEN
           PRINT *,"neo2.f90: Warning switch isw_Vphi_loc=1 is not tested!"
           STOP
           !
           R_Vphi = R_Vphi_prof(ind_boozer_s)
           Z_Vphi = Z_Vphi_prof(ind_boozer_s)
        ELSE IF(isw_Vphi_loc .EQ. 2) THEN
           boozer_theta_Vphi = boozer_theta_Vphi_prof(ind_boozer_s)
        ELSE IF (isw_Vphi_loc.LT.0 .OR. isw_Vphi_loc.GT.2) THEN
           PRINT *,"neo2.f90: Undefined state of switch isw_Vphi_loc (= 0 / 1 / 2)!"
           STOP
        END IF
        !
        ! allocate species-tags
        IF(ALLOCATED(species_tag_vec)) DEALLOCATE(species_tag_vec)
        ALLOCATE(species_tag_vec(num_spec))
        ! allocate conl_over_mfp_vec
        IF(ALLOCATED(conl_over_mfp_vec)) DEALLOCATE(conl_over_mfp_vec)
        ALLOCATE(conl_over_mfp_vec(num_spec))
        ! allocate z_vec
        IF(ALLOCATED(z_vec)) DEALLOCATE(z_vec)
        ALLOCATE(z_vec(num_spec))
        ! allocate m_vec
        IF(ALLOCATED(m_vec)) DEALLOCATE(m_vec)
        ALLOCATE(m_vec(num_spec))
        ! allocate T_vec
        IF(ALLOCATED(T_vec)) DEALLOCATE(T_vec)
        ALLOCATE(T_vec(num_spec))
        ! allocate n_vec
        IF(ALLOCATED(n_vec)) DEALLOCATE(n_vec)
        ALLOCATE(n_vec(num_spec))
        ! allocate dT_vec_ov_ds
        IF(ALLOCATED(dT_vec_ov_ds)) DEALLOCATE(dT_vec_ov_ds)
        ALLOCATE(dT_vec_ov_ds(num_spec))
        ! allocate dn_vec_ov_ds
        IF(ALLOCATED(dn_vec_ov_ds)) DEALLOCATE(dn_vec_ov_ds)
        ALLOCATE(dn_vec_ov_ds(num_spec))
        !
        ! fill the arrays
        ctr_spec = 0
        DO ind_spec = 1,num_species_all
           IF (n_prof(ind_boozer_s,ind_spec) .LE. 0.0_dp) CYCLE
           ctr_spec = ctr_spec + 1
           IF (ctr_spec .GT. num_spec) THEN
              PRINT *,"neo2.f90: Error during preparation of &
                   &multi-species computations!"
              PRINT *,"Number of density-values inconsistent &
                   &with number of relevant species!"
              STOP
           END IF
           species_tag_vec(ctr_spec) = species_tag_prof(ind_spec)
           conl_over_mfp_vec(ctr_spec) = -kappa_prof(ind_boozer_s,ind_spec)
           z_vec(ctr_spec) = species_def_prof(ind_boozer_s,ind_spec,1)
           m_vec(ctr_spec) = species_def_prof(ind_boozer_s,ind_spec,2)
           T_vec(ctr_spec) = T_prof(ind_boozer_s,ind_spec)
           n_vec(ctr_spec) = n_prof(ind_boozer_s,ind_spec)
           dT_vec_ov_ds(ctr_spec) = dT_ov_ds_prof(ind_boozer_s,ind_spec)
           dn_vec_ov_ds(ctr_spec) = dn_ov_ds_prof(ind_boozer_s,ind_spec)
        END DO
        !
        ! write multi-species input file (namelist)
        !
        OPEN(unit=u1,file=fname_multispec,action='write',iostat=ios)
        IF (ios .NE. 0) THEN
           PRINT *, 'WARNING: File ',fname_multispec,' cannot be OPENED!'
           PRINT *, ''
           STOP
        ELSE
           ! write variables into groups
           ! multi-species part:
           WRITE(u1,nml=multi_spec,iostat=ios)
           IF (ios .NE. 0) THEN
              PRINT *, 'WARNING: group multi_spec in ',fname_multispec,' cannot be WRITTEN!'
              PRINT *, ''
              STOP
           END IF
           WRITE(u1,nml=settings,iostat=ios)
           IF (ios .NE. 0) THEN
              PRINT *, 'WARNING: group settings in ',fname_multispec,' cannot be WRITTEN!'
              PRINT *, ''
              STOP
           END IF
           WRITE(u1,nml=collision,iostat=ios)
           IF (ios .NE. 0) THEN
              PRINT *, 'WARNING: group collision in ',fname_multispec,' cannot be WRITTEN!'
              PRINT *, ''
              STOP
           END IF
           WRITE(u1,nml=binsplit,iostat=ios)
           IF (ios .NE. 0) THEN
              PRINT *, 'WARNING: group binsplit in ',fname_multispec,' cannot be WRITTEN!'
              PRINT *, ''
              STOP
           END IF
           WRITE(u1,nml=propagator,iostat=ios)
           IF (ios .NE. 0) THEN
              PRINT *, 'WARNING: group propagator in ',fname_multispec,' cannot be WRITTEN!'
              PRINT *, ''
              STOP
           END IF
           WRITE(u1,nml=plotting,iostat=ios)
           IF (ios .NE. 0) THEN
              PRINT *, 'WARNING: group plotting in ',fname_multispec,' cannot be WRITTEN!'
              PRINT *, ''
              STOP
           END IF
           !! Modification by Andreas F. Martitsch (17.07.2014)
           ! ntv_input
           WRITE(u1,nml=ntv_input,iostat=ios)
           IF (ios .NE. 0) THEN
              PRINT *, 'WARNING: group ntv_input in ',fname_multispec,' cannot be WRITTEN!'
              PRINT *, ''
              STOP
           END IF
           !! End Modification by Andreas F. Martitsch (17.07.2014)
        END IF
        CLOSE(unit=u1)
        !
        ! go back to initial directory
        CALL chdir('..')
        !
     END DO
     !   
     STOP
  END IF
  !! End Modification by Andreas F. Martitsch (20.02.2017)

  !! Modification by Andreas F. Martitsch (23.08.2015)
  ! multi-species part:
  ! -> read species-relevant info into a large array (allocatable not supported)
  !
  IF(ALLOCATED(species_tag)) DEALLOCATE(species_tag)
  ALLOCATE(species_tag(0:num_spec-1))
  species_tag(0:num_spec-1)=species_tag_vec(1:num_spec)
  !
  IF(ALLOCATED(conl_over_mfp_spec)) DEALLOCATE(conl_over_mfp_spec)
  ALLOCATE(conl_over_mfp_spec(0:num_spec-1))
  conl_over_mfp_spec(0:num_spec-1)=conl_over_mfp_vec(1:num_spec)
  !
  IF(ALLOCATED(z_spec)) DEALLOCATE(z_spec)
  ALLOCATE(z_spec(0:num_spec-1))
  z_spec(0:num_spec-1)=z_vec(1:num_spec)
  !
  IF(ALLOCATED(m_spec)) DEALLOCATE(m_spec)
  ALLOCATE(m_spec(0:num_spec-1))
  m_spec(0:num_spec-1)=m_vec(1:num_spec)
  !
  IF(ALLOCATED(T_spec)) DEALLOCATE(T_spec)
  ALLOCATE(T_spec(0:num_spec-1))
  T_spec(0:num_spec-1)=T_vec(1:num_spec)
  !
  IF(ALLOCATED(n_spec)) DEALLOCATE(n_spec)
  ALLOCATE(n_spec(0:num_spec-1))
  n_spec(0:num_spec-1)=n_vec(1:num_spec)
  !
  IF(ALLOCATED(dT_spec_ov_ds)) DEALLOCATE(dT_spec_ov_ds)
  ALLOCATE(dT_spec_ov_ds(0:num_spec-1))
  dT_spec_ov_ds(0:num_spec-1)=dT_vec_ov_ds(1:num_spec)
  !
  IF(ALLOCATED(dn_spec_ov_ds)) DEALLOCATE(dn_spec_ov_ds)
  ALLOCATE(dn_spec_ov_ds(0:num_spec-1))
  dn_spec_ov_ds(0:num_spec-1)=dn_vec_ov_ds(1:num_spec) 
  !
  ! print multi-species input
  IF(lsw_multispecies) THEN
     PRINT *,'isw_coul_log       : ',isw_coul_log
     PRINT *,'num_spec           : ',num_spec
     PRINT *,'species_tag        : ',species_tag
     PRINT *,'conl_over_mfp_spec : ',conl_over_mfp_spec
     PRINT *,'z_spec             : ',z_spec
     PRINT *,'m_spec             : ',m_spec
     PRINT *,'T_spec             : ',T_spec
     PRINT *,'n_spec             : ',n_spec
     !STOP
  END IF
  !! End Modification by Andreas F. Martitsch (23.08.2015) 

  IF (mag_magfield .EQ. 0) THEN ! homogeneous case
     PRINT *, 'WARNING: some input quantities modified - homogeneous case!'
     phi_split_mode = 1
     phi_place_mode = 1
     bin_split_mode = 0
     mag_coordinates = 0 ! cylindrical
     mag_symmetric = .FALSE.
  END IF
  ! ---------------------------------------------------------------------------
  ! end of reading
  ! ---------------------------------------------------------------------------

  !**********************************************************
  ! Initialize MPI module
  !**********************************************************
  CALL mpro%init()
  
  
  !****************************************************
  !  Git version check
  !*****************************************************
  IF (mpro%isMaster()) THEN
    CALL write_version_info()
  END IF

 
 
  !! Modification by Andreas F. Martitsch (31.07.2014)
  ! Save here starting point of the field line for cylindircal
  ! coordinates (used for normalizations for final NTV output)
  xstart_cyl = (/rbeg,phimi,zbeg/)
  !! End Modification by Andreas F. Martitsch (31.07.2014)


!!$  ! ---------------------------------------------------------------------------
!!$  ! test sparse solver
!!$  sparse_talk = .TRUE.
!!$  sparse_solve_method = 1
!!$  CALL sparse_example(2)
!!$  STOP
!!$  ! ---------------------------------------------------------------------------

  IF (prop_reconstruct .EQ. 1) THEN
     PRINT *, 'Reconstruction run!'
     CALL reconstruct_prop_dist
     PRINT *, 'No further calculations!'
     STOP
  END IF

  ! ---------------------------------------------------------------------------
  ! some settings
  ! nmat=npart*npart
  ndim=ndim0
  ! allocation of some arrays (should be moved)
  ! this part was not touched
  ialloc=1
  CALL kin_allocate(ialloc)
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  ! prepare the whole configuration
  CALL flint_prepare(phimi,rbeg,zbeg,nstep,nperiod,bin_split_mode,eta_s_lim)

  ! ---------------------------------------------------------------------------
  ! matrix elements
  ! ---------------------------------------------------------------------------
  !IF (isw_integral .EQ. 0 .AND. isw_energy .EQ. 0) THEN
  !   isw_lorentz = 1
  !END IF
  IF (isw_momentum .EQ. 0) THEN ! Laguerre
     CALL collop_construct
     CALL collop_load
     nvel = lag
  ELSEIF (isw_momentum .EQ. 1) THEN ! Grid
     nvel = vel_num
     IF (vel_distri_swi .EQ. 0) THEN
        CALL linspace(0.0_dp,vel_max,vel_num+1,vel_array)
        !print *, 'vel_array ',lbound(vel_array,1),ubound(vel_array,1)
        !print *, vel_array
     ELSE
        PRINT *, 'vel_distri_swi = ',vel_distri_swi,' not implemented!'
        STOP
     END IF
  ELSE
     PRINT *, 'isw_momentum = ',isw_momentum,' not implemented!'
     STOP
  END IF
  ! ---------------------------------------------------------------------------
  ! erase arrays
  ! ---------------------------------------------------------------------------
  IF (mpro%isMaster() .AND. prop_reconstruct .EQ. 0) THEN
     ! find free unit
     uw = 100
     DO
        INQUIRE(unit=uw,opened=opened)
        IF(.NOT. opened) EXIT
        uw = uw + 100
     END DO
     OPEN(uw,file='evolve.dat',status='replace')
     CLOSE(uw)
  END IF
  ! ---------------------------------------------------------------------------



     ! ---------------------------------------------------------------------------
!!$  ! THIS PART WAS MOVED BEFORE COLLOP
!!$  ! ---------------------------------------------------------------------------
!!$  ! some settings
!!$  ! nmat=npart*npart
!!$  ndim=ndim0
!!$  ! allocation of some arrays (should be moved)
!!$  ! this part was not touched
!!$  ialloc=1
!!$  CALL kin_allocate(ialloc)
!!$  ! ---------------------------------------------------------------------------
!!$
!!$  ! ---------------------------------------------------------------------------
!!$  ! prepare the whole configuration
!!$  CALL flint_prepare(phimi,rbeg,zbeg,nstep,nperiod,bin_split_mode,eta_s_lim)
     CALL flint_prepare_2(bin_split_mode,eta_s_lim)

     !*********************************************************
     ! Write information about the run to a HDF5 file
     !*********************************************************
     IF (mpro%isMaster()) THEN
        CALL write_run_info()
     END IF



  ! ---------------------------------------------------------------------------
  ! this is just for christian, sergie please switch it off
  ! CALL sort_theta
  ! nr,nz,nphi
  !CALL write_volume_data(40,40,100,'w7as_vol.dat')
  !CALL write_surface_data('w7as_sur_181.dat')



  ! ---------------------------------------------------------------------------
  ! these are the tags of the first and last fieldpropagator
  ! of the actual fieldline (we have only one at the moment)
  !  fieldline has a first and a last child : fieldperiod
  !  each fieldperiod has a first and a last child : fieldpropagator
  !  for each structure there is a tag which numbers it in
  !   ascending order
  fieldperiod => fieldline%ch_fir
  DO WHILE (fieldperiod%extra .EQ. 1) 
     fieldperiod => fieldperiod%next
  END DO
  proptag_first = fieldperiod%ch_fir%tag

  fieldperiod => fieldline%ch_las
  DO WHILE (fieldperiod%extra .EQ. 1) 
     fieldperiod => fieldperiod%prev
  END DO
  proptag_last = fieldperiod%ch_las%tag
  ! here one can pick whatever one likes between proptag_first and
  !  proptag_last (now from input file)
  IF (proptag_begin .GT. proptag_last) proptag_begin = proptag_last
  IF (proptag_begin .LT. proptag_first) THEN
     proptag_start = proptag_first
     IF (proptag_final .LE. 0 .OR. proptag_final .GT. proptag_last) THEN
        proptag_end = proptag_last
     ELSE
        proptag_end = proptag_final
     END IF
  ELSE
     proptag_start = proptag_begin
     IF (proptag_final .LE. 0 .OR. proptag_final .GT. proptag_last) THEN
        proptag_end = proptag_start - 1
        IF (proptag_end .LT. proptag_first) proptag_end = proptag_last
     ELSE
        proptag_end = proptag_final
     END IF
  END IF
  !
  !IF (proptag_start .LE. proptag_end) THEN
  ! ------------------------------------------------------------------------
  ! real computation
  CALL flint(eta_part_globalfac,eta_part_globalfac_p,eta_part_globalfac_t, &
       eta_alpha_p,eta_alpha_t,                                            &
       xetami,xetama,eta_part,lambda_equi,                                 &
       eta_part_global,eta_part_trapped,                                   &
       bin_split_mode,                                                     &
       proptag_start,proptag_end)
  ! ------------------------------------------------------------------------
  !ELSE
  !   PRINT *, 'NOTHING TO COMPUTE'
  !END IF


  ! ---------------------------------------------------------------------------
  ! final deallocation of device and all its children
  !PRINT *, 'Before destruct_magnetics'
  !CALL destruct_magnetics(device)
  ! ---------------------------------------------------------------------------
  ! final deallocation of device
  !PRINT *, 'Beforecollop_unload'
  IF (isw_momentum .EQ. 0) THEN
     CALL collop_unload
     !PRINT *, 'Beforecollop_deconstruct'
     CALL collop_deconstruct
  END IF

  !*******************************************
  ! Deinitialize MPI module
  !*******************************************
  CALL mpro%deinit(.FALSE.)
  CALL h5_deinit()


  !! Modification by Andreas F. Martitsch (17.07.2014)
  ! Uncomment "STOP" to see IEEE floating-point exceptions (underflow is present).
  ! This FORTRAN 2008 feature is implemented in gfortran-4.9.2.
  ! See change-log entry (https://gcc.gnu.org/gcc-4.9/changes.html):
  ! "When STOP or ERROR STOP are used to terminate the execution and any exception
  ! (but inexact) is signaling, a warning is printed to ERROR_UNIT, indicating which
  ! exceptions are signaling. The -ffpe-summary= command-line option can be used to
  ! fine-tune for which exceptions the warning should be shown.
  STOP
  !! End Modification by Andreas F. Martitsch (17.07.2014)
  

CONTAINS

  SUBROUTINE write_run_info()
    IF (prop_reconstruct .EQ. 0 ) THEN
       CALL h5_create('neo2_config.h5', h5_config_id, 1)

       CALL h5_define_group(h5_config_id, 'metadata', h5_config_group)
       !CALL h5_add(h5_config_group, 'NEO-2 Version', Neo2_Version)
       !CALL h5_add(h5_config_group, 'MPILib Version', MyMPILib_Version)
       !CALL h5_add(h5_config_group, 'CMake_Compiler', CMake_Compiler)
       !CALL h5_add(h5_config_group, 'CMake_Compiler_Version', CMake_Compiler_Version)
       !CALL h5_add(h5_config_group, 'CMake_Build_Type', CMake_Build_Type)
       !CALL h5_add(h5_config_group, 'CMake_Flags', CMake_Flags)
       !CALL h5_add(h5_config_group, 'CMake_Flags_Release', CMake_Flags_Release)
       !CALL h5_add(h5_config_group, 'CMake_Flags_Debug', CMake_Flags_Debug)
       !CALL h5_add(h5_config_group, 'CMake_System', CMake_System)
       !CALL h5_add(h5_config_group, 'CMake_SuiteSparse_Dir', CMake_SuiteSparse_Dir)
       !CALL h5_add(h5_config_group, 'CMake_Blas_Lib', CMake_Blas_Lib)
       CALL h5_close_group(h5_config_group)

       CALL h5_define_group(h5_config_id, 'neo', h5_config_group)
       CALL h5_add(h5_config_group, 'in_file', in_file, 'Boozer file')
       CALL h5_add(h5_config_group, 'lab_swi', lab_swi)
       CALL h5_add(h5_config_group, 'inp_swi', inp_swi)
       CALL h5_close_group(h5_config_group)

       CALL h5_define_group(h5_config_id, 'settings', h5_config_group)
       CALL h5_add(h5_config_group, 'phimi', phimi, 'Beginning of period', 'Rad')
       CALL h5_add(h5_config_group, 'nstep', nstep, 'Number of integration steps per period', 'Rad')
       CALL h5_add(h5_config_group, 'nperiod', nperiod, 'Number of periods')
       CALL h5_add(h5_config_group, 'mag_nperiod_min', mag_nperiod_min)
       CALL h5_add(h5_config_group, 'mag_magfield', mag_magfield)
       CALL h5_add(h5_config_group, 'boozer_theta_beg', boozer_theta_beg)
       CALL h5_add(h5_config_group, 'boozer_phi_beg', boozer_phi_beg)
       CALL h5_add(h5_config_group, 'mag_start_special', mag_start_special)
       CALL h5_add(h5_config_group, 'mag_cycle_ripples', mag_cycle_ripples)
       CALL h5_add(h5_config_group, 'mag_close_fieldline', mag_close_fieldline)
       CALL h5_add(h5_config_group, 'aiota_tokamak', aiota_tokamak)
       CALL h5_add(h5_config_group, 'eta_part_global', eta_part_global)
       CALL h5_add(h5_config_group, 'eta_part_globalfac', eta_part_globalfac)
       CALL h5_add(h5_config_group, 'eta_part_globalfac_p',eta_part_globalfac_p )
       CALL h5_add(h5_config_group, 'eta_part_globalfac_t', eta_part_globalfac_t)
       CALL h5_add(h5_config_group, 'eta_part_trapped',eta_part_trapped )
       CALL h5_add(h5_config_group, 'eta_alpha_p', eta_alpha_p)
       CALL h5_add(h5_config_group, 'eta_alpha_t', eta_alpha_t)
       CALL h5_add(h5_config_group, 'sparse_solve_method', sparse_solve_method)
       CALL h5_add(h5_config_group, 'magnetic_device', magnetic_device, 'Magnetic device (0: Tokamak, 1: W7-AS)')
       CALL h5_add(h5_config_group, 'mag_coordinates', mag_coordinates, '0: Cylindrical, 1: Boozer')
       CALL h5_add(h5_config_group, 'boozer_s', boozer_s, 'Flux surface')
       CALL h5_add(h5_config_group, 'xetami', xetami)
       CALL h5_add(h5_config_group, 'xetama', xetama)
       CALL h5_add(h5_config_group, 'ndim0', ndim0)
       CALL h5_add(h5_config_group, 'zbeg', zbeg)
       CALL h5_add(h5_config_group, 'rbeg', rbeg)
       CALL h5_add(h5_config_group, 'proptag_begin', proptag_begin)
       CALL h5_add(h5_config_group, 'proptag_final', proptag_final)
       CALL h5_add(h5_config_group, 'mag_save_memory', mag_save_memory)
       CALL h5_add(h5_config_group, 'mag_dbhat_min', mag_dbhat_min)
       CALL h5_add(h5_config_group, 'mag_dphi_inf_min', mag_dphi_inf_min)
       CALL h5_add(h5_config_group, 'mag_inflection_mult', mag_inflection_mult)
       CALL h5_add(h5_config_group, 'solver_talk', solver_talk)
       CALL h5_add(h5_config_group, 'switch_off_asymp', switch_off_asymp) 
       CALL h5_add(h5_config_group, 'asymp_margin_zero', asymp_margin_zero)
       CALL h5_add(h5_config_group, 'asymp_margin_npass', asymp_margin_npass)
       CALL h5_add(h5_config_group, 'asymp_pardeleta', asymp_pardeleta)
       CALL h5_add(h5_config_group, 'ripple_solver_accurfac', ripple_solver_accurfac)
       CALL h5_add(h5_config_group, 'sparse_talk', sparse_talk)
       CALL h5_add(h5_config_group, 'mag_symmetric', mag_symmetric)
       CALL h5_add(h5_config_group, 'mag_symmetric_shorten', mag_symmetric_shorten)
       CALL h5_close_group(h5_config_group)

       CALL h5_define_group(h5_config_id, 'collision', h5_config_group)
       CALL h5_add(h5_config_group, 'conl_over_mfp', conl_over_mfp, 'Collisionality parameter')
       CALL h5_add(h5_config_group, 'lag', lag, 'Number of Laguerre polynomials')
       CALL h5_add(h5_config_group, 'leg', leg, 'Number of Legendre polynomials')
       CALL h5_add(h5_config_group, 'legmax', legmax, 'Maximum number of Legendre polynomials')
       CALL h5_add(h5_config_group, 'z_eff', z_eff, 'Effective charge')
       CALL h5_add(h5_config_group, 'isw_lorentz', isw_lorentz, '')
       CALL h5_add(h5_config_group, 'isw_integral', isw_integral, '')
       CALL h5_add(h5_config_group, 'isw_energy', isw_energy, '')
       CALL h5_add(h5_config_group, 'isw_axisymm', isw_axisymm, '')
       CALL h5_add(h5_config_group, 'collop_path', collop_path, 'Path to collision operator matrix')
       CALL h5_add(h5_config_group, 'collop_base_prj', collop_base_prj, 'Projection base of collision operator')
       CALL h5_add(h5_config_group, 'collop_base_exp', collop_base_exp, 'Expansion base of collision operator')
       CALL h5_add(h5_config_group, 'v_min_resolution', v_min_resolution, 'Minimum velocity for level placement')
       CALL h5_add(h5_config_group, 'v_max_resolution', v_max_resolution, 'Maximum velocity for level placement')
       CALL h5_add(h5_config_group, 'phi_x_max', phi_x_max, 'Maximum velocity for base function')
       CALL h5_add(h5_config_group, 'collop_bspline_order', collop_bspline_order, 'BSpline order')
       CALL h5_add(h5_config_group, 'collop_bspline_dist', collop_bspline_dist, 'BSpline knots distribution factor')
       CALL h5_add(h5_config_group, 'isw_relativistic', isw_relativistic)
       CALL h5_add(h5_config_group, 'T_e', T_e)

       CALL h5_close_group(h5_config_group)

       CALL h5_define_group(h5_config_id, 'binsplit', h5_config_group)
       CALL h5_add(h5_config_group, 'eta_s_lim', eta_s_lim)
       CALL h5_add(h5_config_group, 'eta_part',eta_part )
       CALL h5_add(h5_config_group, 'lambda_equi', lambda_equi)
       CALL h5_add(h5_config_group, 'phi_split_mode', phi_split_mode)
       CALL h5_add(h5_config_group, 'phi_place_mode', phi_place_mode)
       CALL h5_add(h5_config_group, 'phi_split_min', phi_split_min)
       CALL h5_add(h5_config_group, 'max_solver_try', max_solver_try)
       CALL h5_add(h5_config_group, 'hphi_mult', hphi_mult)
       CALL h5_add(h5_config_group, 'bin_split_mode', bin_split_mode)
       CALL h5_add(h5_config_group, 'bsfunc_message', bsfunc_message)
       CALL h5_add(h5_config_group, 'bsfunc_modelfunc', bsfunc_modelfunc)
       CALL h5_add(h5_config_group, 'bsfunc_local_err', bsfunc_local_err)
       CALL h5_add(h5_config_group, 'bsfunc_min_distance',bsfunc_min_distance )
       CALL h5_add(h5_config_group, 'bsfunc_max_index',bsfunc_max_index )
       CALL h5_add(h5_config_group, 'bsfunc_max_splitlevel', bsfunc_max_splitlevel)
       CALL h5_add(h5_config_group, 'bsfunc_sigma_mult',bsfunc_sigma_mult)
       CALL h5_add(h5_config_group, 'bsfunc_sigma_min', bsfunc_sigma_min)
       CALL h5_add(h5_config_group, 'bsfunc_local_solver', bsfunc_local_solver)
       CALL h5_add(h5_config_group, 'mag_local_sigma',mag_local_sigma )
       CALL h5_add(h5_config_group, 'bsfunc_divide', bsfunc_divide)
       CALL h5_add(h5_config_group, 'mag_ripple_contribution', mag_ripple_contribution)
       CALL h5_close_group(h5_config_group)

       CALL h5_define_group(h5_config_id, 'propagator', h5_config_group)
       CALL h5_add(h5_config_group, 'prop_diagphys', prop_diagphys)
       CALL h5_add(h5_config_group, 'prop_overwrite', prop_overwrite)
       CALL h5_add(h5_config_group, 'prop_diagnostic', prop_diagnostic)
       CALL h5_add(h5_config_group, 'prop_binary', prop_binary)
       CALL h5_add(h5_config_group, 'prop_timing', prop_timing)
       CALL h5_add(h5_config_group, 'prop_join_ends',prop_join_ends )
       CALL h5_add(h5_config_group, 'prop_fluxsplitmode', prop_fluxsplitmode)
       CALL h5_add(h5_config_group, 'mag_talk', mag_talk)
       CALL h5_add(h5_config_group, 'mag_infotalk',mag_infotalk )
       CALL h5_add(h5_config_group, 'hphi_lim', hphi_lim)
       CALL h5_add(h5_config_group, 'prop_write', prop_write)
       CALL h5_add(h5_config_group, 'prop_reconstruct',prop_reconstruct )
       !CALL h5_add(h5_config_group, 'prop_fileformat',prop_fileformat )
       CALL h5_add(h5_config_group, 'prop_ripple_plot',prop_ripple_plot )
       CALL h5_close_group(h5_config_group)

       CALL h5_define_group(h5_config_id, 'plotting', h5_config_group)
       CALL h5_add(h5_config_group, 'plot_gauss', plot_gauss)
       CALL h5_add(h5_config_group, 'plot_prop', plot_prop)
       CALL h5_close_group(h5_config_group)        

       CALL h5_close(h5_config_id)

    END IF

    !**********************************************************
    ! Save timestamps
    !**********************************************************
    CALL h5_open_rw('neo2_config.h5', h5_config_id)
    CALL h5_open_group(h5_config_id, 'metadata', h5_config_group)

    CALL random_NUMBER(rand_num)
    WRITE (rand_hash,'(Z6.6)') CEILING(16**6 * rand_num)

    WRITE (fieldname, '(A,I1)') "Run_", prop_reconstruct
    CALL h5_delete(h5_config_group, fieldname)
    CALL h5_add(h5_config_group, fieldname, rand_hash)

    CALL date_and_TIME(date,time)
    WRITE (datetimestring, '(A, A)') date, time

    WRITE (fieldname, '(A,I1)') "Timestamp_start_", prop_reconstruct
    CALL h5_delete(h5_config_group, fieldname)
    CALL h5_add(h5_config_group, fieldname, datetimestring)
    WRITE (fieldname, '(A,I1)') "Nodes_", prop_reconstruct
    CALL h5_delete(h5_config_group, fieldname)
    CALL h5_add(h5_config_group, fieldname, mpro%getNumProcs())

    CALL h5_close_group(h5_config_group)
    CALL h5_close(h5_config_id)

  END SUBROUTINE write_run_info

  SUBROUTINE write_version_info()

    WRITE (*,*) ''
    WRITE (*,*) "---------- NEO-2 Git Revision ----------"
    WRITE (*,*) Neo2_Version
    WRITE (*,*) Neo2_Version_Date
    WRITE (*,*) "----------------------------------------"
    WRITE (*,*) ''
    IF (len_TRIM(Neo2_Version_Additional) /= 0) THEN
       WRITE (*,*) "#################################### NEO-2 Git Additional Information ####################################"
       WRITE (*,*) Neo2_Version_Additional
       WRITE (*,*) "##########################################################################################################"
       WRITE (*,*) ''
    END IF
    
  END SUBROUTINE write_version_info  

END PROGRAM neo2
