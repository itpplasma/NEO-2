PROGRAM neo2

! Includes the classes for MPI support, which need a compiler supporting Fortran 2003 standards
! To compile the -cpp flag has to be added for c preprocessor
#if defined(MPI_SUPPORT)
  USE mpiprovider_module
  USE parallelstorage_module, ONLY : globalstorage
#endif

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
  USE collop, ONLY : collop_construct, collop_deconstruct,          &
       collop_load, collop_unload, z_eff
  USE rkstep_mod, ONLY : lag,leg,legmax
  USE development, ONLY : solver_talk,switch_off_asymp, &
       asymp_margin_zero, asymp_margin_npass, asymp_pardeleta,      &
       ripple_solver_accurfac
  USE sparse_mod, ONLY : sparse_talk,sparse_solve_method,sparse_example

  IMPLICIT NONE

  INTEGER, PARAMETER :: dp = KIND(1.0d0)

  LOGICAL :: opened

  ! Used for MPI support
  ! This string is used to give every client an own evolve.dat file
  character(len=32) :: strEvolveFilename
  include "version.f90"

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
  NAMELIST /collision/                                                        &
       conl_over_mfp,lag,leg,legmax,z_eff,isw_lorentz,                        &
       isw_integral,isw_energy,isw_axisymm,                                   &
       isw_momentum,vel_distri_swi,vel_num,vel_max
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
  bsfunc_shield = .false.
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
  mag_local_sigma = 0
  mag_symmetric = .false.
  mag_symmetric_shorten = .false.
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

  ! reading
  DO jf = 1,SIZE(fnames)
     OPEN(unit=u1,file=fnames(jf),status='old',iostat=ios)
     IF (ios .NE. 0) THEN
        PRINT *, 'WARNING: File ',fnames(jf),' cannot be OPENED!'
        PRINT *, ''
     ELSE
        ! Read variables from group settings
        READ(u1,nml=settings,iostat=ios)
        IF (ios .NE. 0) THEN
           PRINT *, 'WARNING: group settings in ',fnames(jf),' cannot be READ!'
           PRINT *, ''
        END IF
        READ(u1,nml=collision,iostat=ios)
        IF (ios .NE. 0) THEN
           PRINT *, 'WARNING: group collision in ',fnames(jf),' cannot be READ!'
           PRINT *, ''
        END IF
        READ(u1,nml=binsplit,iostat=ios)
        IF (ios .NE. 0) THEN
           PRINT *, 'WARNING: group binsplit in ',fnames(jf),' cannot be READ!'
           PRINT *, ''
        END IF
        READ(u1,nml=propagator,iostat=ios)
        IF (ios .NE. 0) THEN
           PRINT *, 'WARNING: group propagator in ',fnames(jf),' cannot be READ!'
           PRINT *, ''
        END IF
        READ(u1,nml=plotting,iostat=ios)
        IF (ios .NE. 0) THEN
           PRINT *, 'WARNING: group plotting in ',fnames(jf),' cannot be READ!'
           PRINT *, ''
        END IF
     END IF
     CLOSE(unit=u1)
  END DO
  ! PAUSE

  IF (mag_magfield .EQ. 0) THEN ! homogeneous case
     PRINT *, 'WARNING: some input quantities modified - homogeneous case!'
     phi_split_mode = 1
     phi_place_mode = 1
     bin_split_mode = 0
     mag_coordinates = 0 ! cylindrical
     mag_symmetric = .false.
  END IF
  ! ---------------------------------------------------------------------------
  ! end of reading
  ! ---------------------------------------------------------------------------

  ! MPI support
#if defined(MPI_SUPPORT)
  call mpro%init()

  if (mpro%getRank() == 1) then
     write (*,*) "NEO-2 Version ", Neo2_Version
  end if

#endif


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
  if (isw_momentum .eq. 0) then ! Laguerre
     CALL collop_construct
     CALL collop_load
     nvel = lag
  elseif (isw_momentum .eq. 1) then ! Grid
     nvel = vel_num
     if (vel_distri_swi .eq. 0) then
        CALL linspace(0.0_dp,vel_max,vel_num+1,vel_array)
        !print *, 'vel_array ',lbound(vel_array,1),ubound(vel_array,1)
        !print *, vel_array
     else
        print *, 'vel_distri_swi = ',vel_distri_swi,' not implemented!'
        stop
     end if
  else
     print *, 'isw_momentum = ',isw_momentum,' not implemented!'
     stop
  end if
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

#if defined(MPI_SUPPORT)
        ! Give every client an own evolve.dat file
        write (globalstorage%evolveFilename, "(A, I3.3, A)"), 'evolve', mpro%getRank(), '.dat'
        OPEN(uw,file=globalstorage%evolveFilename, status='replace')
        CLOSE(uw)
#else
        ! Sequential behaviour
        OPEN(uw,file='evolve.dat',status='replace')
        CLOSE(uw)
#endif

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

  ! MPI support
#if defined(MPI_SUPPORT)
  call mpro%deinit()
#endif
  ! ---------------------------------------------------------------------------
  ! final deallocation of device and all its children
  !PRINT *, 'Before destruct_magnetics'
  !CALL destruct_magnetics(device)
  ! ---------------------------------------------------------------------------
     ! final deallocation of device
     !PRINT *, 'Beforecollop_unload'
     if (isw_momentum .eq. 0) then
        CALL collop_unload
        !PRINT *, 'Beforecollop_deconstruct'
        CALL collop_deconstruct
     end if

  STOP

END PROGRAM neo2
