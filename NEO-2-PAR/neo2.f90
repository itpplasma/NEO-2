PROGRAM neo2

  !**********************************************************
  ! MPI SUPPORT
  !**********************************************************
  USE mpiprovider_module

  USE size_mod
  !USE partpa_mod, ONLY : hxeta
  USE flint_mod, ONLY : plot_gauss,plot_prop,phi_split_mode,        &
       phi_place_mode,phi_split_min,hphi_mult,max_solver_try,       &
       bsfunc_local_err_max_mult,bsfunc_max_mult_reach,             &
       bsfunc_modelfunc_num,bsfunc_divide,                          &
       bsfunc_ignore_trap_levels,boundary_dist_limit_factor,        &
       bsfunc_local_shield_factor,bsfunc_shield,                    &
       bsfunc_lambda_loc_res, eta_savemem_dist1, eta_savemem_dist2, &
       eta_savemem_sigma_mult
  USE device_mod
  USE collisionality_mod, ONLY : conl_over_mfp,isw_lorentz,         &
       isw_integral,isw_energy,isw_axisymm,                         &
       isw_momentum,vel_distri_swi,vel_num,vel_max,                 &
       nvel,vel_array,v_max_resolution,v_min_resolution,            &
       phi_x_max, collop_bspline_order, collop_bspline_dist,        &
       isw_relativistic, T_e, lsw_multispecies, num_spec,           &
       conl_over_mfp_spec, z_spec, collop_bspline_taylor, lsw_nbi,  &
       T_nbi, m_nbi
  USE propagator_mod, ONLY : reconstruct_prop_dist,   &
       prop_diagphys,prop_overwrite,                                &
       prop_diagnostic,prop_binary,                                 &
       prop_timing,prop_join_ends,prop_fluxsplitmode,               &
       prop_write,prop_reconstruct,prop_ripple_plot,                &
       prop_reconstruct_levels, prop_fileformat,                    &
       prop_finaljoin_mode, lsw_save_dentf, lsw_save_enetf,         &
       lsw_save_spitf
  USE magnetics_mod, ONLY : mag_talk,mag_infotalk,mag_write_hdf5,   &
       h5_magnetics_file_name
  USE mag_interface_mod, ONLY : mag_local_sigma, hphi_lim,          &
       mag_magfield,mag_nperiod_min,mag_save_memory,                &
       magnetic_device,mag_cycle_ripples,mag_start_special,         &
       aiota_tokamak,mag_close_fieldline,mag_ripple_contribution,   &
       mag_coordinates,boozer_s,boozer_theta_beg,boozer_phi_beg,    &
       mag_dbhat_min,mag_dphi_inf_min,mag_inflection_mult,          &
       mag_symmetric,mag_symmetric_shorten,                         &
       sigma_shield_factor,split_inflection_points,                 &
       split_at_period_boundary
  USE binarysplit_mod, ONLY : bsfunc_message,bsfunc_modelfunc,      &
       bsfunc_total_err, bsfunc_local_err, bsfunc_min_distance,     &
       bsfunc_max_index, bsfunc_max_splitlevel,                     &
       bsfunc_sigma_mult, bsfunc_sigma_min, bsfunc_local_solver
  USE binarysplit_int, ONLY : linspace
  USE collop, ONLY : collop_construct, collop_deconstruct,          &
       collop_load, collop_unload, z_eff, collop_path,              &
       collop_base_prj, collop_base_exp, scalprod_alpha, scalprod_beta, &
       collop_only_precompute, lsw_read_precom, lsw_write_precom !! Added lsw_read_precom
       !! and lsw_write_precom by Michael Draxler (25.08.2017)
  USE rkstep_mod, ONLY : lag,leg,legmax, epserr_sink, epserr_iter, &
       niter
  USE development, ONLY : solver_talk,switch_off_asymp, &
       asymp_margin_zero, asymp_margin_npass, asymp_pardeleta,      &
       ripple_solver_accurfac
  USE sparse_mod, ONLY : sparse_talk,sparse_solve_method,sparse_example
  USE neo_control, ONLY: in_file, inp_swi, lab_swi

  !************************************
  ! HDF5
  !************************************
  USE hdf5_tools
  USE hdf5_tools_f2003

  IMPLICIT NONE

  INTEGER, PARAMETER :: dp = KIND(1.0d0)

  LOGICAL :: opened

  !**********************************************************
  ! Include version information
  !**********************************************************
  INCLUDE "cmake_version.f90"
  INCLUDE "version.f90"

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

  REAL(kind=dp), PARAMETER :: pi=3.14159265358979_dp

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
  ! -> read species-relevant info into a large array (allocatable not supported)
  REAL(kind=dp), DIMENSION(1000) :: conl_over_mfp_vec
  REAL(kind=dp), DIMENSION(1000) :: z_vec
  !! End Modification by Andreas F. Martitsch (23.08.2015)  
  ! groups for namelist
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
       sparse_talk,sparse_solve_method,mag_symmetric,mag_symmetric_shorten,   &
       epserr_sink, epserr_iter, niter
  NAMELIST /collision/                                                        &
       conl_over_mfp,lag,leg,legmax,z_eff,isw_lorentz,                        &
       isw_integral,isw_energy,isw_axisymm,                                   &
       isw_momentum,vel_distri_swi,vel_num,vel_max,                           &
       collop_path, collop_base_prj, collop_base_exp,                         &
       scalprod_alpha, scalprod_beta, v_min_resolution, v_max_resolution,     &
       phi_x_max, collop_bspline_order, collop_bspline_dist,                  &
       num_spec, conl_over_mfp_vec, z_vec, collop_only_precompute,            &
       isw_relativistic, T_e, lsw_multispecies, collop_bspline_taylor,        &
       lsw_read_precom, lsw_write_precom, lsw_nbi, T_nbi, m_nbi
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
       split_inflection_points,split_at_period_boundary,                      &
       bsfunc_lambda_loc_res,mag_dbhat_min, eta_savemem_dist1,                &
       eta_savemem_dist2, eta_savemem_sigma_mult
  NAMELIST /propagator/                                                       &
       prop_diagphys,prop_overwrite,                                          &
       prop_diagnostic,prop_binary,prop_timing,prop_join_ends,                &
       prop_fluxsplitmode,                                                    &
       mag_talk,mag_infotalk,mag_write_hdf5,                                  &
       hphi_lim,                                                              &
       prop_write,prop_reconstruct,prop_ripple_plot,                          &
       prop_reconstruct_levels,                                               &
       prop_fileformat, prop_finaljoin_mode,                                  &
       lsw_save_dentf, lsw_save_enetf, lsw_save_spitf
  NAMELIST /plotting/                                                         &
       plot_gauss,plot_prop
  ! ---------------------------------------------------------------------------
  ! filenames (default file and specific input file) for namelist
  fnames = (/'neo2.def','neo2.in '/) ! removed because it is not used
  !fnames = (/'neo2.in '/)

  ! ---------------------------------------------------------------------------
  ! defaults
  !
  ! call omp_set_num_threads(4)
  !
  ! settings
  mag_magfield = 1
  magnetic_device = 1
  mag_nperiod_min = 300
  mag_save_memory = 1
  mag_cycle_ripples = 0
  mag_close_fieldline = 1
  mag_ripple_contribution = 2 ! 1 is not in use any more
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
  epserr_sink = 0d0          ! Regularization
  epserr_iter = 1d-3
  niter = 10
  !  sparse_solve_method = 0
  ! collision
  collop_path = '/afs/itp.tugraz.at/proj/plasma/DOCUMENTS/Neo2/data-MatrixElements/'
  collop_base_prj = 0
  collop_base_exp = 0
  scalprod_alpha = 0d0
  scalprod_beta  = 0d0
  conl_over_mfp = 1.0d-3
  v_min_resolution = 0.1d0
  v_max_resolution = 5.0d0
  phi_x_max        = 5.0d0
  collop_bspline_order = 4
  collop_bspline_dist  = 1d0
  collop_bspline_taylor = .true.
  !! Modification by Andreas F. Martitsch (25.08.2015)
  !  multi-species part
  num_spec = 1
  conl_over_mfp_vec = 0.0d0
  z_vec = 1.0d0
  !! End Modification by Andreas F. Martitsch (25.08.2015)
  lsw_multispecies = .FALSE.
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
  bsfunc_lambda_loc_res = .FALSE.
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
  split_inflection_points = .true.
  split_at_period_boundary = .false.
  eta_savemem_dist1 = 0.1d0
  eta_savemem_dist2 = 0.1d0
  eta_savemem_sigma_mult  = 1d0
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
  prop_fileformat = 0      ! 0... ACSII, 1... HDF5
  prop_reconstruct = 0
  prop_ripple_plot = 0
  prop_reconstruct_levels = 0
  prop_finaljoin_mode = 0
  mag_talk = .TRUE.
  mag_infotalk = .TRUE.
  mag_write_hdf5 = .FALSE.
  hphi_lim = 1.0d-6
  ! plotting
  plot_gauss = 0
  plot_prop  = 0
  lsw_save_dentf = .TRUE.
  lsw_save_enetf = .TRUE.
  lsw_save_spitf = .TRUE.

  lsw_nbi = .false.
  T_nbi   = 70d3
  m_nbi   = 3.343583719d-24
  
  !**********************************
  ! Init HDF5 Fortran interface
  !**********************************
  CALL h5_init()

  !**********************************************************
  ! Init MPI
  !**********************************************************
  CALL mpro%init()
  
  !**********************************************************
  ! Read config files
  !**********************************************************
  DO jf = 1,SIZE(fnames)
     IF(jf .EQ. 1) CYCLE ! skip neo2.def (Andreas F. Martitsch - 21.10.2015)
     OPEN(unit=u1,file=fnames(jf),status='old',iostat=ios)
     IF (ios .NE. 0) THEN
        PRINT *, 'WARNING: File ',fnames(jf),' cannot be OPENED!'
        PRINT *, ''
        STOP
     ELSE
        ! Read variables from group settings
        READ(u1,nml=settings,iostat=ios)
        IF (ios .NE. 0) THEN
           PRINT *, 'WARNING: group settings in ',fnames(jf),' cannot be READ!'
           PRINT *, ''
           STOP
        END IF
        READ(u1,nml=collision,iostat=ios)
        IF (ios .NE. 0) THEN
           PRINT *, 'WARNING: group collision in ',fnames(jf),' cannot be READ!'
           PRINT *, ''
           STOP
        END IF
        READ(u1,nml=binsplit,iostat=ios)
        IF (ios .NE. 0) THEN
           PRINT *, 'WARNING: group binsplit in ',fnames(jf),' cannot be READ!'
           PRINT *, ''
           STOP
        END IF
        READ(u1,nml=propagator,iostat=ios)
        IF (ios .NE. 0) THEN
           PRINT *, 'WARNING: group propagator in ',fnames(jf),' cannot be READ!'
           PRINT *, ''
           STOP
        END IF
        READ(u1,nml=plotting,iostat=ios)
        IF (ios .NE. 0) THEN
           PRINT *, 'WARNING: group plotting in ',fnames(jf),' cannot be READ!'
           PRINT *, ''
           STOP
        END IF
     END IF
     CLOSE(unit=u1)
  END DO

  !! Modification by Andreas F. Martitsch (23.08.2015)
  ! multi-species part:
  ! -> read species-relevant info into a large array (allocatable not supported)
  IF(ALLOCATED(conl_over_mfp_spec)) DEALLOCATE(conl_over_mfp_spec)
  ALLOCATE(conl_over_mfp_spec(0:num_spec-1))
  conl_over_mfp_spec(0:num_spec-1)=conl_over_mfp_vec(1:num_spec)
  !
  IF(ALLOCATED(z_spec)) DEALLOCATE(z_spec)
  ALLOCATE(z_spec(0:num_spec-1))
  z_spec(0:num_spec-1)=z_vec(1:num_spec)
  !
  !PRINT *,conl_over_mfp_spec
  !PRINT *,z_spec
  !STOP
  !! End Modification by Andreas F. Martitsch (23.08.2015) 

  IF (mag_magfield .EQ. 0) THEN ! homogeneous case
     PRINT *, 'WARNING: some input quantities modified - homogeneous case!'
     phi_split_mode = 1
     phi_place_mode = 1
     bin_split_mode = 0
     mag_coordinates = 0 ! cylindrical
     mag_symmetric = .FALSE.
  END IF
  !**********************************************************
  ! End of reading config files
  !**********************************************************

  !**********************************************************
  ! Only precomputation of collision operator
  !**********************************************************
  IF (collop_only_precompute) THEN
     WRITE (*,*) "Precomputation of collision operator..."
     
     CALL collop_construct
     CALL collop_load
     
     WRITE (*,*) "Done."     
     STOP
  END IF

  
  !**********************************************************
  ! Write run information to HDF5
  !**********************************************************
  IF (mpro%isMaster()) THEN
     CALL write_version_info()
  END IF

!!$  ! ---------------------------------------------------------------------------
!!$  ! test sparse solver
!!$  sparse_talk = .TRUE.
!!$  sparse_solve_method = 1
!!$  CALL sparse_example(2)
!!$  STOP
!!$  ! ---------------------------------------------------------------------------  

  !**********************************************************
  ! Check reconstruction switch
  !**********************************************************
  IF (prop_reconstruct .EQ. 0 .OR. prop_reconstruct .EQ. 2) THEN
     IF (prop_reconstruct .EQ. 0) THEN
        IF ( mpro%isMaster() ) THEN
           ! close HDF5-File for magnetics
           OPEN(unit=1234, iostat=ios, file=h5_magnetics_file_name, status='old')
           IF (ios .EQ. 0) CLOSE(unit=1234, status='delete')
        ELSE
           mag_write_hdf5 = .FALSE.
        END IF
     ELSE
        mag_write_hdf5 = .FALSE.
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
     IF (prop_reconstruct .EQ. 0) THEN
        ! find free unit
        uw = 100
        DO
           INQUIRE(unit=uw,opened=opened)
           IF(.NOT. opened) EXIT
           uw = uw + 100
        END DO

     END IF

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
     IF (mpro%isMaster() .AND. prop_fileformat .EQ. 1) THEN
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
     !**********************************************************
     ! Save runtime to HDF5 neo2_config file
     !**********************************************************
     IF (mpro%isMaster()) THEN

        IF (prop_fileformat .EQ. 1) THEN

           !**********************************************************
           ! Open taginfo
           !**********************************************************
           CALL h5_open_rw('taginfo.h5', h5id_taginfo)
           CALL h5_get(h5id_taginfo, 'tag_first', tag_first)
           CALL h5_get(h5id_taginfo, 'tag_last',  tag_last)
           CALL h5_close(h5id_taginfo)     

           !**********************************************************
           ! Merge evolve-files
           !**********************************************************
           CALL h5_create("evolve.h5", h5id_prop)

           IF (isw_axisymm .EQ. 1) THEN
              ! Tokamak mode
              tag_first = 3
              tag_last  = 3
           END IF

!!$           DO k = tag_first, tag_last
!!$              DO l = tag_first, tag_last
!!$                 WRITE (h5_filename, '(I0,A,I0)') k, "_", l
!!$
!!$                 OPEN(unit=1234, iostat=ios, file="evolve_" // TRIM(h5_filename) // ".h5", status='old')
!!$                 CLOSE(unit=1234)
!!$
!!$                 !**********************************************************
!!$                 ! Check if file exists
!!$                 !**********************************************************
!!$                 IF (ios .EQ. 0) THEN
!!$
!!$                    CALL h5_open("evolve_" // TRIM(h5_filename) // ".h5", h5id_propfile)
!!$                    CALL h5_copy(h5id_propfile, '/', h5id_prop, "/" // TRIM(h5_filename))
!!$                    CALL h5_close(h5id_propfile)
!!$
!!$                    ! Delete file
!!$                    OPEN(unit=1234, iostat=ios, file="evolve_" // TRIM(h5_filename) // ".h5", status='old')
!!$                    CLOSE(unit=1234, status='delete')
!!$                 END IF
!!$              END DO
!!$           END DO

           CALL h5_close(h5id_prop)

           CALL h5_open_rw('neo2_config.h5', h5_config_id)
           CALL h5_open_group(h5_config_id, 'metadata', h5_config_group)

           CALL date_and_TIME(date,time)
           WRITE (datetimestring, '(A, A)') date, time
           WRITE (fieldname, '(A,I1)') "Timestamp_stop_", prop_reconstruct
           CALL h5_delete(h5_config_group, fieldname)
           CALL h5_add(h5_config_group, fieldname, datetimestring)

           CALL h5_close_group(h5_config_group)
           CALL h5_close(h5_config_id)
        END IF

     END IF
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
     
  ELSEIF (prop_reconstruct .EQ. 1) THEN
     IF (isw_axisymm .EQ. 1) THEN
        WRITE (*,*) "Skipping reconstruction for axisymmetric mode"
        !stop
     END IF

     IF (mpro%isMaster()) THEN
        PRINT *, 'Reconstruction run!'

        !**********************************************************
        ! If HDF5 is enabled, save information about run
        !**********************************************************
        IF (mpro%isMaster() .AND. prop_fileformat .EQ. 1) THEN
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
        END IF

        CALL reconstruct_prop_dist

        !**********************************************************
        ! If HDF5 is enabled, save runtime in file
        !**********************************************************
        IF (mpro%isMaster() .AND. prop_fileformat .EQ. 1) THEN
           CALL h5_open_rw('neo2_config.h5', h5_config_id)
           CALL h5_open_group(h5_config_id, 'metadata', h5_config_group)

           CALL date_and_TIME(date,time)
           WRITE (datetimestring, '(A, A)') date, time
           WRITE (fieldname, '(A,I1)') "Timestamp_stop_", prop_reconstruct
           CALL h5_delete(h5_config_group, fieldname)
           CALL h5_add(h5_config_group, fieldname, datetimestring)

           CALL h5_close_group(h5_config_group)
           CALL h5_close(h5_config_id)
        END IF


        PRINT *, 'No further calculations!'
     END IF
     CALL mpro%barrier()

  ELSEIF (prop_reconstruct .EQ. 3) THEN
     !**********************************************************
     ! Reconstruction 3: Collect HDF5 files
     !**********************************************************
     CALL prop_reconstruct_3()
  END IF

  !**********************************************************
  ! Deinitialize HDF5
  !**********************************************************
  CALL h5_deinit()

  !**********************************************************
  ! Deinitialize MPI
  !**********************************************************
  CALL mpro%deinit()

  !**********************************************************
  ! Quit NEO-2 with STOP to show more compiler warnings
  !**********************************************************
  STOP

CONTAINS

  SUBROUTINE prop_reconstruct_3()

    IF (mpro%isMaster()) THEN
       WRITE (*,*) "Merging HDF5 files..."

       !**********************************************************
       ! Open taginfo
       !**********************************************************
       CALL h5_open('taginfo.h5', h5id_taginfo)
       CALL h5_get(h5id_taginfo, 'tag_first', tag_first)
       CALL h5_get(h5id_taginfo, 'tag_last',  tag_last)
       CALL h5_close(h5id_taginfo)

       IF (isw_axisymm .EQ. 1) THEN
          ! Tokamak mode
          tag_first = 3
          tag_last  = 3
       END IF

       !**********************************************************
       ! Get current working directory
       ! This might only work with gfortran
       !**********************************************************
       CALL getcwd(cwd)
       WRITE (surfname,'(A)') TRIM(ADJUSTL(cwd(INDEX(cwd, '/', .TRUE.)+1:)))
       WRITE (*,*) "Using " // TRIM(ADJUSTL(surfname)) // " as surfname."

       !**********************************************************
       ! Create result file
       !**********************************************************
       CALL h5_create('final.h5', h5id_final, 2)
       CALL h5_define_group(h5id_final, surfname, h5id_surf)
       CALL h5_define_group(h5id_surf, 'NEO-2', h5id_neo2)
       CALL h5_define_group(h5id_neo2, 'propagators', h5id_propagators)

       IF (ALLOCATED(cg0_1_num_prop)) DEALLOCATE(cg0_1_num_prop)
       IF (ALLOCATED(cg2_1_num_prop)) DEALLOCATE(cg2_1_num_prop)
       IF (ALLOCATED(cg0_2_num_prop)) DEALLOCATE(cg0_2_num_prop)
       IF (ALLOCATED(cg2_2_num_prop)) DEALLOCATE(cg2_2_num_prop)
       IF (ALLOCATED(cg0_3_num_prop)) DEALLOCATE(cg0_3_num_prop)
       IF (ALLOCATED(cg2_3_num_prop)) DEALLOCATE(cg2_3_num_prop)
       IF (ALLOCATED(denom_mflint_prop)) DEALLOCATE(denom_mflint_prop)

       ALLOCATE(cg0_1_num_prop(tag_first:tag_last))
       ALLOCATE(cg2_1_num_prop(tag_first:tag_last))
       ALLOCATE(cg0_2_num_prop(tag_first:tag_last))
       ALLOCATE(cg2_2_num_prop(tag_first:tag_last))
       ALLOCATE(cg0_3_num_prop(tag_first:tag_last))
       ALLOCATE(cg2_3_num_prop(tag_first:tag_last))

       ALLOCATE(denom_mflint_prop(tag_first:tag_last))

       !**********************************************************
       ! Iterate over all propagators and merge files
       !**********************************************************
       DO k = tag_first, tag_last
          WRITE (h5_filename, '(I0)') k
          CALL h5_define_group(h5id_propagators, h5_filename, h5id_prop)
          
          IF (lsw_save_spitf) THEN
             CALL h5_open("spitf_" // TRIM(h5_filename) // ".h5", h5id_propfile)
             CALL h5_copy(h5id_propfile, '/', h5id_prop, "spitf")
             CALL h5_close(h5id_propfile)
          END IF

          IF (lsw_save_dentf) THEN
             CALL h5_open("dentf_" // TRIM(h5_filename) // ".h5", h5id_propfile)
             CALL h5_copy(h5id_propfile, '/', h5id_prop, "dentf")
             CALL h5_close(h5id_propfile)
          END IF

          IF (lsw_save_enetf) THEN
             CALL h5_open("enetf_" // TRIM(h5_filename) // ".h5", h5id_propfile)
             CALL h5_copy(h5id_propfile, '/', h5id_prop, "enetf")
             CALL h5_close(h5id_propfile)
          END IF
          
          CALL h5_open("phi_mesh_" // TRIM(h5_filename) // ".h5", h5id_propfile)
          CALL h5_get(h5id_propfile, 'cg0_1_num', cg0_1_num_prop(k))
          CALL h5_get(h5id_propfile, 'cg2_1_num', cg2_1_num_prop(k))
          CALL h5_get(h5id_propfile, 'cg0_2_num', cg0_2_num_prop(k))
          CALL h5_get(h5id_propfile, 'cg2_2_num', cg2_2_num_prop(k))
          CALL h5_get(h5id_propfile, 'cg0_3_num', cg0_3_num_prop(k))
          CALL h5_get(h5id_propfile, 'cg2_3_num', cg2_3_num_prop(k))
          CALL h5_get(h5id_propfile, 'denom_mflint', denom_mflint_prop(k))
          CALL h5_copy(h5id_propfile, '/', h5id_prop, "phi_mesh")

          CALL h5_close(h5id_propfile)

          CALL h5_open("sizeplot_etalev_" // TRIM(h5_filename) // ".h5", h5id_propfile)
          CALL h5_copy(h5id_propfile, '/', h5id_prop, "sizeplot_etalev")
          CALL h5_close(h5id_propfile)

          CALL h5_close_group(h5id_prop)
       END DO

       cg0_1_avg = SUM(cg0_1_num_prop(tag_first:tag_last)) / SUM(denom_mflint_prop(tag_first:tag_last))
       cg2_1_avg = SUM(cg2_1_num_prop(tag_first:tag_last)) / SUM(denom_mflint_prop(tag_first:tag_last))
       cg0_2_avg = SUM(cg0_2_num_prop(tag_first:tag_last)) / SUM(denom_mflint_prop(tag_first:tag_last))
       cg2_2_avg = SUM(cg2_2_num_prop(tag_first:tag_last)) / SUM(denom_mflint_prop(tag_first:tag_last))
       cg0_3_avg = SUM(cg0_3_num_prop(tag_first:tag_last)) / SUM(denom_mflint_prop(tag_first:tag_last))
       cg2_3_avg = SUM(cg2_3_num_prop(tag_first:tag_last)) / SUM(denom_mflint_prop(tag_first:tag_last))

       WRITE (*,*) "cg0_1 = ", cg0_1_avg
       WRITE (*,*) "cg2_1 = ", cg2_1_avg
       WRITE (*,*) "cg0_2 = ", cg0_2_avg
       WRITE (*,*) "cg2_2 = ", cg2_2_avg
       WRITE (*,*) "cg0_3 = ", cg0_3_avg
       WRITE (*,*) "cg2_3 = ", cg2_3_avg

       !**********************************************************
       ! Merge additional files
       !**********************************************************
       CALL h5_open("efinal.h5", h5id_propfile)
       CALL h5_copy(h5id_propfile, '/', h5id_neo2, "efinal")
       CALL h5_close(h5id_propfile)

       CALL h5_open("fulltransp.h5", h5id_propfile)
       CALL h5_copy(h5id_propfile, '/', h5id_neo2, "fulltransp")
       CALL h5_close(h5id_propfile)

       CALL h5_open("neo2_config.h5", h5id_propfile)
       CALL h5_copy(h5id_propfile, '/', h5id_neo2, "neo2_config")
       CALL h5_close(h5id_propfile)

       CALL h5_open_rw("taginfo.h5", h5id_propfile)
       CALL h5_delete(h5id_propfile, 'cg0_1_avg')
       CALL h5_add(h5id_propfile, 'cg0_1_avg', cg0_1_avg)
       CALL h5_delete(h5id_propfile, 'cg2_1_avg')
       CALL h5_add(h5id_propfile, 'cg2_1_avg', cg2_1_avg)

       CALL h5_delete(h5id_propfile, 'cg0_2_avg')
       CALL h5_add(h5id_propfile, 'cg0_2_avg', cg0_2_avg)
       CALL h5_delete(h5id_propfile, 'cg2_2_avg')
       CALL h5_add(h5id_propfile, 'cg2_2_avg', cg2_2_avg)

       CALL h5_delete(h5id_propfile, 'cg0_3_avg')
       CALL h5_add(h5id_propfile, 'cg0_3_avg', cg0_3_avg)
       CALL h5_delete(h5id_propfile, 'cg2_3_avg')
       CALL h5_add(h5id_propfile, 'cg2_3_avg', cg2_3_avg)

       CALL h5_copy(h5id_propfile, '/', h5id_neo2, "taginfo")
       CALL h5_close(h5id_propfile)

       CALL h5_close_group(h5id_propagators)
       CALL h5_close_group(h5id_neo2)
       CALL h5_close_group(h5id_surf)
       CALL h5_close(h5id_final)

       WRITE (*,*) "Result file created. Deleting old files...."

       !**********************************************************
       ! Delete single HDF5 files
       !**********************************************************
       OPEN(unit=1234, iostat=ios, file="propagator_0_0.h5", status='old')
       IF (ios .EQ. 0) CLOSE(unit=1234, status='delete')

       !OPEN(unit=1234, iostat=ios, file="efinal.h5", status='old')
       !IF (ios .EQ. 0) CLOSE(unit=1234, status='delete')

       !OPEN(unit=1234, iostat=ios, file="fulltransp.h5", status='old')
       !IF (ios .EQ. 0) CLOSE(unit=1234, status='delete')

       !open(unit=1234, iostat=ios, file="neo2_config.h5", status='old')
       !if (ios .eq. 0) close(unit=1234, status='delete')

       OPEN(unit=1234, iostat=ios, file="taginfo.h5", status='old')
       IF (ios .EQ. 0) CLOSE(unit=1234, status='delete')    

       DO k = tag_first, tag_last
          DO l = tag_first, tag_last
             WRITE (h5_filename, '(I0,A,I0)') k, "_", l

             OPEN(unit=1234, iostat=ios, file="propagator_" // TRIM(h5_filename) // ".h5", status='old')
             IF (ios .EQ. 0) CLOSE(unit=1234, status='delete')

             OPEN(unit=1234, iostat=ios, file="propagator_boundary_" // TRIM(h5_filename) // ".h5", status='old')
             IF (ios .EQ. 0) CLOSE(unit=1234, status='delete')

             OPEN(unit=1234, iostat=ios, file="reconstruct_" // TRIM(h5_filename) // ".h5", status='old')
             IF (ios .EQ. 0) CLOSE(unit=1234, status='delete')

             OPEN(unit=1234, iostat=ios, file="binarysplit_" // TRIM(h5_filename) // ".h5", status='old')
             IF (ios .EQ. 0) CLOSE(unit=1234, status='delete')
          END DO

       END DO

       DO k = tag_first, tag_last
          WRITE (h5_filename, '(I0)') k

          IF (lsw_save_spitf) THEN
             OPEN(unit=1234, iostat=ios, file="spitf_" // TRIM(h5_filename) // ".h5", status='old')
             IF (ios .EQ. 0) CLOSE(unit=1234, status='delete')
          END IF

          IF (lsw_save_enetf) THEN
             OPEN(unit=1234, iostat=ios, file="enetf_" // TRIM(h5_filename) // ".h5", status='old')
             IF (ios .EQ. 0) CLOSE(unit=1234, status='delete')
          END IF

          IF (lsw_save_dentf) THEN
             OPEN(unit=1234, iostat=ios, file="dentf_" // TRIM(h5_filename) // ".h5", status='old')
             IF (ios .EQ. 0) CLOSE(unit=1234, status='delete')
          END IF
          
          OPEN(unit=1234, iostat=ios, file="phi_mesh_" // TRIM(h5_filename) // ".h5", status='old')
          IF (ios .EQ. 0) CLOSE(unit=1234, status='delete')

          OPEN(unit=1234, iostat=ios, file="sizeplot_etalev_" // TRIM(h5_filename) // ".h5", status='old')
          IF (ios .EQ. 0) CLOSE(unit=1234, status='delete')

       END DO
       WRITE (*,*) "Done."

    END IF

    CALL mpro%barrier()

  END SUBROUTINE prop_reconstruct_3

  SUBROUTINE write_version_info()

    WRITE (*,*) ''
    WRITE (*,*) "---------- NEO-2 Git Revision ----------"
    WRITE (*,*) Neo2_Version
    WRITE (*,*) "----------------------------------------"
    WRITE (*,*) ''
    IF (len_TRIM(Neo2_Version_Additional) /= 0) THEN
       WRITE (*,*) "#################################### NEO-2 Git Additional Information ####################################"
       WRITE (*,*) Neo2_Version_Additional
       WRITE (*,*) "##########################################################################################################"
       WRITE (*,*) ''
    END IF

  END SUBROUTINE write_version_info

  SUBROUTINE write_run_info()
    IF (prop_reconstruct .EQ. 0 ) THEN
       CALL h5_create('neo2_config.h5', h5_config_id, 1)

       CALL h5_define_group(h5_config_id, 'metadata', h5_config_group)
       CALL h5_add(h5_config_group, 'NEO-2 Version', Neo2_Version)
       CALL h5_add(h5_config_group, 'MPILib Version', MyMPILib_Version)
       CALL h5_add(h5_config_group, 'CMake_Compiler', CMake_Compiler)
       CALL h5_add(h5_config_group, 'CMake_Compiler_Version', CMake_Compiler_Version)
       CALL h5_add(h5_config_group, 'CMake_Build_Type', CMake_Build_Type)
       CALL h5_add(h5_config_group, 'CMake_Flags', CMake_Flags)
       CALL h5_add(h5_config_group, 'CMake_Flags_Release', CMake_Flags_Release)
       CALL h5_add(h5_config_group, 'CMake_Flags_Debug', CMake_Flags_Debug)
       CALL h5_add(h5_config_group, 'CMake_System', CMake_System)
       CALL h5_add(h5_config_group, 'CMake_SuiteSparse_Dir', CMake_SuiteSparse_Dir)
       CALL h5_add(h5_config_group, 'CMake_Blas_Lib', CMake_Blas_Lib)
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
       CALL h5_add(h5_config_group, 'epserr_iter', epserr_iter)
       CALL h5_add(h5_config_group, 'epserr_sink', epserr_sink)
       CALL h5_add(h5_config_group, 'niter', niter)

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
       CALL h5_add(h5_config_group, 'collop_bspline_taylor', collop_bspline_taylor, 'BSpline taylor expansion')
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
       CALL h5_add(h5_config_group, 'prop_reconstruct',prop_reconstruct)
       CALL h5_add(h5_config_group, 'prop_fileformat',prop_fileformat)
       CALL h5_add(h5_config_group, 'prop_ripple_plot',prop_ripple_plot)
       CALL h5_add(h5_config_group, 'prop_finaljoin_mode',prop_finaljoin_mode)

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

END PROGRAM neo2
