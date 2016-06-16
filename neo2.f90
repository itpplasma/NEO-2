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
       nvel,vel_array
  USE propagator_mod, ONLY : reconstruct_prop_dist,   &
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
       scalprod_beta, num_spec, conl_over_mfp_spec, z_spec
  USE rkstep_mod, ONLY : lag,leg,legmax                            
      
  USE development, ONLY : solver_talk,switch_off_asymp, &
       asymp_margin_zero, asymp_margin_npass, asymp_pardeleta,      &
       ripple_solver_accurfac
  USE sparse_mod, ONLY : sparse_talk,sparse_solve_method,sparse_example
  !! Modification by Andreas F. Martitsch (14.07.2015)
  ! Extra input for NTV computations
  USE ntv_mod, ONLY : isw_ntv_mode, isw_qflux_NA, in_file_pert,&
       MtOvR, B_rho_L_loc, xstart_cyl, isw_ripple_solver
  !! End Modification by Andreas F. Martitsch (14.07.2015)
  !
  IMPLICIT NONE

  INTEGER, PARAMETER :: dp = KIND(1.0d0)

  LOGICAL :: opened

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
       sparse_talk,sparse_solve_method,mag_symmetric,mag_symmetric_shorten
  !! Modifications by Andreas F. Martitsch (15.07.2014)
  ! Path for the collision operator matrices is now specified via neo2.in
  ! (necessary for computations with Condor)
  NAMELIST /collision/                                                        &
       conl_over_mfp,lag,leg,legmax,z_eff,isw_lorentz,                        &
       isw_integral,isw_energy,isw_axisymm,                                   &
       isw_momentum,vel_distri_swi,vel_num,vel_max,collop_path,               &      
       !! End Modifications by Andreas F. Martitsch (15.07.2014)
       collop_base_prj, collop_base_exp, scalprod_alpha, scalprod_beta,       &
       num_spec, conl_over_mfp_vec, z_vec
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
       isw_ripple_solver
  !! End Modification by Andreas F. Martitsch (14.07.2015)  
  ! ---------------------------------------------------------------------------
  ! filenames (default file and specific input file) for namelist
  fnames = (/'neo2.def','neo2.in '/)
  ! ---------------------------------------------------------------------------
  ! defaults
  !
  ! settings
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
  !! Modification by Andreas F. Martitsch (25.08.2015)
  !  multi-species part
  num_spec = 1
  conl_over_mfp_vec = 0.0d0
  z_vec = 1.0d0
  !! End Modification by Andreas F. Martitsch (25.08.2015)
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
  !! End Modification by Andreas F. Martitsch (14.07.2015)

  CALL h5_init()
  
  ! reading
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
        !! Modification by Andreas F. Martitsch (17.07.2014)
        ! ntv_input
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
  !! Modification by Andreas F. Martitsch (23.08.2015)
  ! multi-species part:
  ! -> read species-relevant info into a large array (allocatable not supported)
  IF(ALLOCATED(conl_over_mfp_spec)) DEALLOCATE(conl_over_mfp_spec)
  ALLOCATE(conl_over_mfp_spec(0:num_spec-1))
  conl_over_mfp_spec(0:num_spec-1)=conl_over_mfp_vec(1:num_spec)
  IF(num_spec .EQ. 1) conl_over_mfp_spec(0)=conl_over_mfp
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
  ! ---------------------------------------------------------------------------
  ! end of reading
  ! ---------------------------------------------------------------------------

  !**********************************************************
  ! Initialize MPI module
  !**********************************************************
  CALL mpro%init()
 
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
     OPEN(uw,file='evolve.dat',status='replace')
     CLOSE(uw)
  END IF
  ! ---------------------------------------------------------------------------



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
  
END PROGRAM neo2
