SUBROUTINE flint_prepare(phimi,rbeg,zbeg,nstep,nperiod,bin_split_mode,eta_s_lim)
  ! this routine is called before flint and does the setup for magnetics
  ! it creates the 
  !  device (stevvo, ....)
  !  surface
  !  fieldline and its children, grandchildren, grandgrandchildren
  !    fieldperiod
  !    fieldpropagator
  !    fieldripple
  
  ! input
  USE size_mod, ONLY : ndim0
  ! input/output
  USE rk4_kin_mod, ONLY : y
  ! output
  USE partpa_mod,  ONLY : bmod0
  ! new module which will handle all magnetics
  USE magnetics_mod
  USE device_mod
  USE mag_interface_mod
  USE collisionality_mod, ONLY : collpar,conl_over_mfp,collpar_min,collpar_max,&
       lsw_multispecies, num_spec, conl_over_mfp_spec, collpar_spec
  USE compute_aiota_mod, ONLY : compute_aiota
  !! Modifications by Andreas F. Martitsch (09.03.2014)
  ! Collection of subroutines (mag.f90) converted to a module.
  ! This allows one to make use of generic interfaces 
  ! for module procedures (e.g., for a different number of
  ! arguments).
  ! Note: This requires changes in "flint_prepare" and 
  ! "write_volume_data" (both flint.f90), and "rhs_kin"
  ! and "magdata_for_particles". 
  USE mag_sub, ONLY: mag
  !! End Modifications by Andreas F. Martitsch (09.03.2014)
  USE mpiprovider_module
  use math_constants, only : pi
  use nrtype, only : dp

  IMPLICIT NONE


  ! parameter list
  REAL(kind=dp),    INTENT(in) :: phimi
  REAL(kind=dp),    INTENT(in) :: rbeg
  REAL(kind=dp),    INTENT(inout) :: zbeg ! changed when EFIT is used
  INTEGER,          INTENT(in) :: nstep
  INTEGER,          INTENT(in) :: nperiod
  INTEGER,          INTENT(in) :: bin_split_mode
  REAL(kind=dp),    INTENT(in) :: eta_s_lim
  ! end of parameter list

  ! Winny's new stuff
  REAL(kind=dp)                         :: xstart(3)
  INTEGER                                  :: i

  ! only for testing
  REAL(kind=dp) :: ehlp
  REAL(kind=dp), ALLOCATABLE :: aeta_x0(:), aeta_s(:)
  TYPE(dnumber_struct), POINTER :: eta_x0 => NULL()
  TYPE(dnumber_struct), POINTER :: eta_s => NULL()
  INTEGER :: test_unit = 200
  LOGICAL :: opened

  ! only to call mag magnetics (should be avoided)
  REAL(kind=dp)                 :: bmod,sqrtg
  REAL(kind=dp), DIMENSION(3)   :: x,bder,hcovar,hctrvr
  REAL(kind=dp), DIMENSION(3,3) :: hcoder,hctder

  ! efit
  INTEGER :: ierr
  !

  ! species index
  INTEGER :: ispec, ispecp
  
  IF (mag_coordinates .EQ. 0) THEN
     IF (mag_magfield .EQ. 3) THEN ! EFIT
        ! only to provide efit_raxis,efit_zaxis,aiota_tokamak
        x(1) = rbeg
        x(2) = 0.0d0
        x(3) = 0.0d0
        CALL compute_aiota(rbeg,efit_raxis,efit_zaxis,aiota_tokamak,ierr)
        !print *, 'aiota_tokamak = ', aiota_tokamak
        IF (ierr .EQ. 1) THEN
           PRINT *, 'Wrong radius for EFIT - I better stop'
           STOP
        END IF
        zbeg = efit_zaxis
        x(3) = zbeg
        CALL mag(x,bmod,sqrtg,bder,hcovar,hctrvr,hcoder,hctder)
     END IF
     ! cylindrical coordinates
     ! this is only here because of this call to mag
     ! for computation of bmod0 (we should get rid of this)
     !  bmod0 should come from somewhere else
     y(1)=rbeg       ! R
     y(2)=zbeg       ! Z
     y(3)=1.d0
     y(4:ndim0)=0.d0
     x(1)=y(1)   ! R 
     x(2)=phimi  ! phi
     x(3)=y(2)   ! Z
     ! this creates the device (look in magnetics for device_struct)
     !print *, 'flint: before make_magnetics'
     CALL make_magnetics('W7_AS')
     ! this is here only for computation of bmod0 (avoid)
     !print *, 'flint: before mag'
     CALL mag(x,bmod,sqrtg,bder,hcovar,hctrvr,hcoder,hctder)
     bmod0=bmod ! seems to be the first value 
     ! this now creates surface and fieldline
     !print *, 'flint: before make_magnetics(bmod0,nperiod,nstep,ndim0)'
     CALL make_magnetics(bmod0,nperiod,nstep,ndim0)
     xstart = (/rbeg,phimi,zbeg/) 
     !print *, 'flint: before make_magnetics(xstart)'
     CALL make_magnetics(xstart)
     !print *, 'flint: after make_magnetics'
  ELSE
     ! boozer
     y(1)=boozer_theta_beg
     y(2)=1.d0
     y(3)=-1.d0
     y(4:ndim0)=0.d0
     x(1)=boozer_s
     x(2)=boozer_phi_beg
     x(3)=boozer_theta_beg
     ! this creates the device (look in magnetics for device_struct)
     CALL make_magnetics('Boozer')
     ! this is here only for computation of bmod0 (avoid)
     CALL mag(x,bmod,sqrtg,bder,hcovar,hctrvr,hcoder,hctder)
     !bmod0=bmod ! seems to be the first value
     ! this now creates surface and fieldline
     bmod0 = boozer_bmod0
     CALL make_magnetics(bmod0,nperiod,nstep,ndim0)
     xstart = (/boozer_s,boozer_phi_beg,boozer_theta_beg/) 
     CALL make_magnetics(xstart)
  END IF

  IF (.NOT. lsw_multispecies) THEN
     ! negative input for conl_over_mfp should provide collpar directly
     IF (conl_over_mfp .GT. 0.0d0) THEN
        collpar=4.d0/(2.d0*pi*device%r0)*conl_over_mfp
     ELSE
        collpar=-conl_over_mfp
     END IF
     collpar_max = collpar
     collpar_min = collpar
  ELSE
     IF(ALLOCATED(collpar_spec)) DEALLOCATE(collpar_spec)
     ALLOCATE(collpar_spec(0:num_spec-1))
     ! negative input for conl_over_mfp_spec should provide collpar_spec directly
     IF (ALL(conl_over_mfp_spec .GT. 0.0d0)) THEN
        collpar_spec=4.d0/(2.d0*pi*device%r0)*conl_over_mfp_spec
        !PRINT *,'1'
        !STOP
     ELSE IF (ALL(conl_over_mfp_spec .LE. 0.0d0)) THEN
        collpar_spec=-conl_over_mfp_spec
        !PRINT *,'2'
        !STOP
     ELSE
        PRINT *,"flint.f90: Values of collisionality &
             &parameter conl_over_spec not consistent!"
        PRINT *,"flint.f90: All entries of conl_over_spec &
             &must be >0 or <=0"
        STOP
     END IF
     ! species index
     ispec = mpro%getRank()
     ! set collisionality parameter for species 'ispec'
     collpar = collpar_spec(ispec)
     conl_over_mfp = conl_over_mfp_spec(ispec)
     ! set default parameters for eta-grid refinement
     ! analogue to single-species version:
     ! collpar_max and collpar_min are the sane for all species
     ! (ToDo: species-dependent refinement - re-discretization)
     collpar_max = collpar_spec(0) 
     collpar_min = collpar_spec(0)
  END IF

END SUBROUTINE flint_prepare

subroutine flint_prepare_2(bin_split_mode,eta_s_lim)
  USE collisionality_mod, ONLY : collpar_min
  USE  mag_interface_mod, ONLY : ripple_eta_magnetics
  use nrtype, only : dp

  IMPLICIT NONE

  ! parameter
  INTEGER,          INTENT(in) :: bin_split_mode
  REAL(kind=dp),    INTENT(in) :: eta_s_lim
  ! find the eta values for splitting
  !  for this the collision parameter is needed
  IF (bin_split_mode .EQ. 1) THEN
     CALL ripple_eta_magnetics(collpar_min,eta_s_lim)
  END IF
  ! print *, 'After ripple_eta_magnetics'
END SUBROUTINE flint_prepare_2


SUBROUTINE flint(eta_part_globalfac,eta_part_globalfac_p,eta_part_globalfac_t, &
     eta_alpha_p,eta_alpha_t,                                                  &
     xetami,xetama,eta_part,lambda_equi,                                       &
     eta_part_global,eta_part_trapped,                                         &
     bin_split_mode,proptag_start,proptag_end)
  ! this routine now walks through the specified fieldpropagators
  !  (see main routine about how to set proptag_start and proptag_end)
  USE size_mod, ONLY : npart
  USE flint_mod, ONLY : plot_gauss,plot_prop,phi_split_mode,        &
       phi_place_mode,phi_split_min,hphi_mult,max_solver_try,       &
       bsfunc_local_err_max_mult,bsfunc_max_mult_reach,             &
       bsfunc_modelfunc_num,bsfunc_divide,                          &
       bsfunc_ignore_trap_levels,boundary_dist_limit_factor,        &
       bsfunc_local_shield_factor,bsfunc_shield,                    &
       bsfunc_lambda_loc_res, eta_savemem_dist1, eta_savemem_dist2, &
       eta_savemem_sigma_mult
  USE magnetics_mod
  USE device_mod
  USE mag_interface_mod, ONLY : mag_save_memory,mag_start_special,magnetic_device, &
       mag_magfield,ripple_prop_joiner,boozer_phi_beg,boozer_theta_beg
  USE propagator_mod, ONLY : prop_ibegperiod,prop_count_call,propagator_solver, &
       prop_reconstruct, prop_reconstruct_levels, prop_write, prop_ctaginfo
  USE collisionality_mod, ONLY : isw_lorentz,isw_axisymm,y_axi_averages, &
       isw_momentum
  !vel_distri_swi,vel_num,vel_max,nvel,vel_array

  USE rkstep_mod, ONLY : lag ! anumm => anumm_lag
  USE collisionality_mod, ONLY : collpar,collpar_min,collpar_max
  use math_constants, only : twopi
  use nrtype, only : dp
  !
  ! types and routines for splitting
  USE binarysplit_mod
  !
  !************************************
  ! HDF5
  !************************************
  USE hdf5_tools
  USE hdf5_tools_f2003

  !**************************************
  ! MPI Support
  ! Source has to be compiled with -cpp flag for C preprocessor directives
  !**************************************
#if defined(MPI_SUPPORT)
  USE neo2scheduler_module
  USE parallelStorage_module
#endif
  ! ---

  IMPLICIT NONE
  !
  ! parameter list
  REAL(kind=dp), INTENT(in)       :: eta_part_globalfac,eta_part_globalfac_p,eta_part_globalfac_t
  REAL(kind=dp), INTENT(in)       :: eta_alpha_p,eta_alpha_t
  REAL(kind=dp), INTENT(in)  :: xetami
  REAL(kind=dp), INTENT(inout)  :: xetama
  INTEGER, INTENT(in) :: eta_part,lambda_equi
  INTEGER, INTENT(in) :: bin_split_mode
  INTEGER, INTENT(inout) :: proptag_start,proptag_end
  INTEGER, INTENT(in)       :: eta_part_global,eta_part_trapped

  ! locals
  INTEGER :: rippletag_old,rippletag,proptag
  INTEGER :: iend,iendperiod,ierr_solv,ierr_join
  INTEGER :: i,j,count_solv
  INTEGER :: x2_ub,i_min_sav,clear_old_ripple
  INTEGER :: proptag_first,proptag_last
  INTEGER :: start_at_begin

  INTEGER :: i_construct,i_construct_s,i_construct_e,i_construct_d
  INTEGER :: eta_min_loc(1)

  INTEGER :: lag_sigma,ilag
  real(kind=dp), DIMENSION(:), ALLOCATABLE :: collision_sigma_multiplier
  real(kind=dp) :: mult_sigma, mult_sigma_mod
  TYPE(binarysplit) :: eta_bs,eta_bs_loc,   eta_bs_store
  REAL(kind=dp), ALLOCATABLE :: eta_x0(:),eta_s(:)
  REAL(kind=dp), ALLOCATABLE :: eta_x0_loc(:),eta_s_loc(:),eta_x0_hlp(:)
  INTEGER,       ALLOCATABLE :: eta_shield_loc(:)
  REAL(kind=dp), ALLOCATABLE :: eta_cl(:)
  REAL(kind=dp), ALLOCATABLE :: eta_ori(:),eta_split(:),eta_split_loc(:)
  REAL(kind=dp), ALLOCATABLE :: lam_ori(:)
  REAL(kind=dp), ALLOCATABLE :: eta_x0_val(:),eta_s_val(:)
  REAL(kind=dp), ALLOCATABLE :: eta_x0_sqrt(:),eta_s_sqrt(:)

  REAL(kind=dp), ALLOCATABLE :: x2_sav(:),bhat_sav(:)
  REAL(kind=dp) :: eta_min_relevant, save_bsfunc_err                       !<-in
  REAL(kind=dp) :: eta_min_global
  REAL(kind=dp) :: eta_b_abs_max
  REAL(kind=dp) :: eta_b_abs_min
  REAL(kind=dp) :: eta_s_loc_min
  REAL(kind=dp) :: t_start,d_p,d_t

  REAL(kind=dp) :: eta_s_relevant
  REAL(kind=dp) :: xe1,xe2,xe3,xe4,de1,de2
  INTEGER       :: eta_part_1,eta_part_2
  REAL(kind=dp), ALLOCATABLE :: xhlp_1(:),xhlp_2(:),xhlp_3(:)

  CHARACTER(len=20)  :: c_ripple_tag,c_period_tag,c_propagator_tag
  CHARACTER(len=100) :: c_filename
  CHARACTER(len=5)   :: c_extprop      = '.prop'
  CHARACTER(len=16)  :: c_fieldprop    = 'fieldpropagator_'
  CHARACTER(len=20)  :: c_fieldpropmod = 'fieldpropagator_mod_'
  CHARACTER(len=8)   :: c_etapropori   = 'eta_ori_'
  CHARACTER(len=10)  :: c_etapropsplit = 'eta_split_'

  TYPE(fieldpropagator_struct), POINTER :: plotpropagator

  INTEGER :: ibmf,idiv

  TYPE(fieldripple_struct), POINTER :: fieldripple_left
  TYPE(fieldripple_struct), POINTER :: fieldripple_right
  REAL(kind=dp) :: previous_dist,boundary_dist,boundary_dist_limit
  REAL(kind=dp) :: next_dist
  REAL(kind=dp) :: old_eta,new_eta,move_eta,opp_limit,opp_limit_safe
  INTEGER :: eb_ind
  INTEGER :: eb_ind_left,eb_ind_right
  INTEGER :: boundary_counter_fixed = 0
  INTEGER :: boundary_counter_partly_fixed = 0
  INTEGER :: boundary_counter_not_fixed = 0

  INTEGER :: boundary_fix_mode
  INTEGER :: b_loop_1,b_loop_2
  INTEGER :: boundary_fix_counter
  INTEGER :: boundary_fix_counter_max = 3
  INTEGER :: boundary_eta_fix_counter
  LOGICAL :: first_ripple,boundary_has_to_be_fixed
  TYPE(dnumber_struct), POINTER :: pointer_eta_close_to_boundary => NULL()
  TYPE(dnumber_struct), POINTER :: pointer_all_bmax => NULL()
  TYPE(dnumber_struct), POINTER :: pointer_ripple_close_to_boundary => NULL()
  TYPE(dnumber_struct), POINTER :: pointer_ripple_all_bmax => NULL()
  TYPE(dnumber_struct), POINTER :: pointer_index_eta_close_to_boundary => NULL()
  TYPE(dnumber_struct), POINTER :: pointer_boundary_dist => NULL()
  TYPE(dnumber_struct), POINTER :: pointer_eta_move => NULL()
  TYPE(dnumber_struct), POINTER :: pointer_next_dist => NULL()
  TYPE(dnumber_struct), POINTER :: pointer_boundary_dist_limit => NULL()
  REAL(kind=dp), ALLOCATABLE :: array_eta_close_to_boundary(:)
  REAL(kind=dp), ALLOCATABLE :: array_all_bmax(:)
  REAL(kind=dp), ALLOCATABLE :: array_ripple_close_to_boundary(:)
  REAL(kind=dp), ALLOCATABLE :: array_ripple_all_bmax(:)
  REAL(kind=dp), ALLOCATABLE :: array_index_eta_close_to_boundary(:)
  REAL(kind=dp), ALLOCATABLE :: array_boundary_dist(:)
  REAL(kind=dp), ALLOCATABLE :: array_eta_move(:)
  REAL(kind=dp), ALLOCATABLE :: array_next_dist(:)
  REAL(kind=dp), ALLOCATABLE :: array_boundary_dist_limit(:)
  LOGICAL,       ALLOCATABLE :: array_eta_used(:)

  REAL(kind=dp), ALLOCATABLE :: eta_shield_hlp(:)
  INTEGER,       ALLOCATABLE :: eta_shield(:)
  REAL(kind=dp), ALLOCATABLE :: eta_type_hlp(:)
  INTEGER,       ALLOCATABLE :: eta_type(:)
  REAL(kind=dp) :: eta_trapped_passing,sigma_trapped_passing 
  REAL(kind=dp) :: eta_highest_local_max 
  REAL(kind=dp) :: save_bsfunc_local_err
  integer :: save_binarysplit_fsplitdepth
  integer :: bsfunc_modelfunc_store

  REAL(kind=dp) :: fieldline_phi_l,fieldline_phi_r,y_conv_factor

  REAL(kind=dp) :: phi_l,phi_r,phi_per,theta_l,theta_r

  !*********************************************
  ! MPI Support
  !*********************************************
#if defined(MPI_SUPPORT)
  ! Define the scheduler
  TYPE(neo2scheduler) :: sched
#endif

  !**********************************************
  ! HDF5 Support
  !**********************************************
  INTEGER(HID_T), DIMENSION(1:10) :: h5ids
  INTEGER                         :: k

  !************************************************
  ! HDF5
  !************************************************
  INTEGER(HID_T) :: h5_config_id
  INTEGER(HID_T) :: h5_config_group

  !**********************************************************
  ! Level decay in trapped domain
  !**********************************************************
  real(kind=dp)  :: eta_savemem_sigma, eta_savemem_sigma_fact
  real(kind=dp)  :: bsfunc_local_err_orig, bsfunc_local_err_decay_sigma, bsfunc_local_err_decay_min_mult

  !print *, 'flint: begin of program'
  ! this is not very sophisticated at the moment
  !  mainly puts the fieldpropagator pointer to the first one
  surface => device%ch_fir
  fieldline => surface%ch_fir
  fieldperiod => fieldline%ch_fir
  fieldpropagator => fieldperiod%ch_fir
  fieldripple => fieldpropagator%ch_act
  proptag_first = fieldline%ch_fir%ch_fir%tag
  proptag_last  = fieldline%ch_las%ch_las%tag

  fieldperiod => fieldline%ch_fir
  IF (mag_magfield .EQ. 0) THEN ! homogeneous
     fieldpropagator => fieldperiod%ch_fir
  ELSE
     fieldpropagator => fieldperiod%ch_ext
  END IF
  !print *, 'flint: before fieldpropagator do'
  DO
     IF (fieldpropagator%tag .EQ. fieldline%abs_max_ptag) EXIT
     IF (ASSOCIATED(fieldpropagator%next)) THEN
        fieldpropagator => fieldpropagator%next
     ELSE
        EXIT
     END IF
  END DO
  !print *, 'flint:  after fieldpropagator do'
  fieldperiod => fieldpropagator%parent
  fieldripple => fieldpropagator%ch_act
!!$  PRINT *, 'abs_max ',fieldline%abs_max_ptag,fieldline%b_abs_max
!!$  PRINT *, 'abs_min ',fieldline%abs_min_ptag,fieldline%b_abs_min
!!$  PRINT *, 'ripple  ',fieldripple%b_max_l,fieldripple%b_max_r
!!$  PRINT *, 'ripple  ',fieldripple%prev%b_max_l,fieldripple%prev%b_max_r
!!$  PRINT *, 'prop-coords    ',fieldripple%pa_fir%coords%x1(0), & 
!!$       fieldripple%pa_fir%coords%x2(0), &
!!$       fieldripple%pa_fir%coords%x3(0)
!!$  PRINT *, 'prop-mdata    ',fieldripple%pa_fir%mdata%bhat(0) 
!!$  PAUSE

  eta_s_relevant=0.d0
  eta_min_global = 1.0d0/fieldline%b_abs_max
  xetama = 1.0_dp / fieldline%b_abs_min / 0.999_dp

  IF (bin_split_mode .NE. 0) THEN
     ! Moved outside the following if-construct
     !print *, 'flint: before extract_array'
     CALL extract_array(fieldripple%eta_x0,eta_x0,1)
     CALL extract_array(fieldripple%eta_x0,eta_x0_hlp,1)
     CALL extract_array(fieldripple%eta_s, eta_s, 1)
     !print *, 'flint:  after extract_array'
     !call disp(fieldripple%eta_x0)

     ALLOCATE(eta_x0_val(1))
     ALLOCATE(eta_s_val(1))
     DO i_construct = 1,2
        eta_min_loc = MINLOC(eta_x0_hlp)
        eta_x0_hlp(eta_min_loc(1)) = 1000.0_dp
        eta_x0_val(1) = eta_x0(eta_min_loc(1))
        eta_s_val(1)  = eta_s(eta_min_loc(1))
        eta_s_relevant=MAX(eta_s_relevant,eta_s_val(1))
     END DO
     ! do it also for the previous ripple
     fieldripple => fieldripple%prev
     !PRINT *, 'EXTRACT 2'
     CALL extract_array(fieldripple%eta_x0,eta_x0,1)
     CALL extract_array(fieldripple%eta_x0,eta_x0_hlp,1)
     CALL extract_array(fieldripple%eta_s, eta_s, 1)
     !PRINT *, 'END EXTRACT 2'
     DO i_construct = 1,2
        eta_min_loc = MINLOC(eta_x0_hlp)
        eta_x0_hlp(eta_min_loc(1)) = 1000.0_dp
        eta_x0_val(1) = eta_x0(eta_min_loc(1))
        eta_s_val(1)  = eta_s(eta_min_loc(1))
        eta_s_relevant=MAX(eta_s_relevant,eta_s_val(1))
     END DO
     DEALLOCATE(eta_x0_val,eta_s_val)
     DEALLOCATE(eta_x0,eta_x0_hlp,eta_s)
     ! Moved outside the following if-construct - End
     ! eta_s_relevant = 0.0_dp
  END IF

  IF (eta_part_global .LT. 0 ) THEN
     ! create linspace of eta_ori values (equidistanz in eta)
     CALL linspace(xetami,xetama,eta_part,eta_ori)
     ! equidistant in lambda
     IF (lambda_equi .EQ. 1) THEN
        CALL linspace(1.0_dp,0.0_dp,eta_part,lam_ori)
        ! 0.999 is put for safety
        eta_ori = (1.0_dp - lam_ori**2) / (fieldline%b_abs_min * 0.999_dp)
        DEALLOCATE(lam_ori)
     END IF
  ELSE
     ! new stuff with base levels around the absolute maximum of B
     !PRINT *, 'EXTRACT 1'

     d_p = MIN( 0.5_dp*(eta_min_global-xetami), eta_part_globalfac_p * eta_s_relevant)
     d_t = MIN( 0.5_dp*(xetama-eta_min_global), eta_part_globalfac_t * eta_s_relevant)

     IF (eta_part_globalfac_t .GT. 0.0_dp) THEN
        d_p = MIN( d_p, eta_part_globalfac_p/eta_part_globalfac_t*d_t)
     END IF
     IF (eta_part_globalfac_p .GT. 0.0_dp) THEN
        d_t = MIN( d_p, eta_part_globalfac_t/eta_part_globalfac_p*d_p)
     END IF

     IF (eta_part_global .EQ. 0) THEN
        d_p = 0.0_dp
        d_t = 0.0_dp
     END IF

     xe2 = eta_min_global - d_p
     xe3 = eta_min_global + d_t

     IF (lambda_equi .EQ. 1) THEN
        xe1 = 1.0_dp
        xe2 = SQRT(1.0_dp - xe2 * fieldline%b_abs_min * 0.999_dp)
        xe3 = SQRT(1.0_dp - xe3 * fieldline%b_abs_min * 0.999_dp)
        xe4 = 0.0_dp
     ELSE
        xe1 = xetami
        xe4 = xetama
     END IF

     !PRINT *, 'xe1-x4 ',xe1,xe2,xe3,xe4

     de1 = ABS(xe2 - xe1)
     de2 = ABS(xe4 - xe3)
     IF (eta_part_trapped .EQ. 0) THEN
        eta_part_1 = MAX(FLOOR(de1/(de1+de2)*DBLE(eta_part)),3)
        eta_part_2 = MAX(eta_part - eta_part_1,3)
     ELSE
        eta_part_1 = MAX(3,eta_part-eta_part_trapped)
        eta_part_2 = MAX(3,eta_part_trapped)
     END IF

     ! old passing
     !CALL linspace(xe1,xe2,eta_part_1,xhlp_1) ! passing
     ! new passing
     ! xe2: is the border of the boundary layer to the passing region
     ! eta_min_global: is the eta at the global maximum of B
     ! eta_alpha_p = 1 (is the old linear behaviour)
     t_start = (1.0_dp - xe2/eta_min_global)**(1.0_dp/eta_alpha_p)
     CALL linspace(1.0_dp,t_start,eta_part_1,xhlp_1)
     xhlp_1(0) = 0.0_dp
     xhlp_1(1:UBOUND(xhlp_1,1)) = eta_min_global * (1.0_dp - xhlp_1(1:UBOUND(xhlp_1,1))**eta_alpha_p)
     ! end of passing

     CALL linspace(xe2,xe3,eta_part_global,xhlp_2) ! boundary - linear

     ! old trapped
     ! CALL linspace(xe3,xe4,eta_part_2,xhlp_3) ! trapped
     ! new trapped
     ! xe3: is the border of the boundary layer to the trapping region
     ! xe4: is the maximum eta 
     ! eta_min_global: is the eta at the global maximum of B
     ! eta_alpha_t = 1 (is the old linear behaviour)
     t_start = ((xe3 - eta_min_global) / (xe4 - eta_min_global))**(1.0_dp/eta_alpha_t)
     CALL linspace(t_start,1.0_dp,eta_part_2,xhlp_3)
     xhlp_3 = eta_min_global + xhlp_3**eta_alpha_t * (xe4 - eta_min_global)
     ! end of trapped

     !PRINT *, 'xhlp_1',xhlp_1
     !PRINT *, 'xhlp_2',xhlp_2
     !PRINT *, 'xhlp_3',xhlp_3

     IF (eta_part_global .EQ. 0) THEN
        ALLOCATE( eta_ori(0:eta_part_1+eta_part_2-2) )
        eta_ori(0:eta_part_1-1) = xhlp_1
        eta_ori(eta_part_1:eta_part_1+eta_part_2-2) = xhlp_3(1:eta_part_2-1)
     ELSE
        ALLOCATE( eta_ori(0:eta_part_1+eta_part_global+eta_part_2-3) )
        eta_ori(0:eta_part_1-1) = xhlp_1
        eta_ori(eta_part_1:eta_part_1+eta_part_global-2) = xhlp_2(1:eta_part_global-1)
        eta_ori(eta_part_1+eta_part_global-1:eta_part_1+eta_part_global+eta_part_2-3) = xhlp_3(1:eta_part_2-1)
     END IF
     IF (lambda_equi .EQ. 1) THEN
        eta_ori = (1.0_dp - eta_ori**2) / (fieldline%b_abs_min * 0.999_dp) ! make eta from lambda
     END IF
     DEALLOCATE (xhlp_1,xhlp_2,xhlp_3)
  END IF

  !PRINT *, 'eta_ori',eta_ori
  !#if !defined(MPI_SUPPORT)
  call h5_magnetics('eta_ori',eta_ori)

  OPEN(9999,file='eta_ori.dat')
  DO i = LBOUND(eta_ori,1),UBOUND(eta_ori,1)
     WRITE(9999,*) i, eta_ori(i)
  END DO
  CLOSE(9999)
  !#endif
  !PAUSE



  ! modify proptag_start, proptag_end
  IF (mag_start_special .GT. 0) THEN
     IF (mag_start_special .EQ. 1) proptag_start = fieldline%abs_max_ptag
     IF (mag_start_special .EQ. 2) proptag_start = fieldline%abs_min_ptag
     IF (mag_start_special .EQ. 3) proptag_start = proptag_start

     IF (proptag_start .LT. proptag_first) proptag_start = proptag_first
     IF (proptag_start .GT. proptag_last)  proptag_start = proptag_last
     proptag_end = proptag_start - 1
     IF (proptag_end .LT. proptag_first) proptag_end = proptag_last
  END IF

  PRINT *, 'proptag_first ',proptag_first
  PRINT *, 'proptag_start ',proptag_start
  PRINT *, 'proptag_end   ',proptag_end
  PRINT *, 'proptag_last  ',proptag_last
  !PAUSE

  ! go to the first propagator which is wanted
  fieldperiod => fieldline%ch_fir 
  fieldpropagator => fieldperiod%ch_fir
  DO WHILE (fieldpropagator%tag .LT. proptag_first)
     IF (.NOT. ASSOCIATED(fieldpropagator%next)) EXIT
     fieldpropagator => fieldpropagator%next
  END DO
  ! internally used for plotting
  plotpropagator => fieldpropagator
  fieldripple => fieldpropagator%ch_act
  eta_b_abs_max = 1.0_dp / fieldline%b_abs_max
  eta_b_abs_min = 1.0_dp / fieldline%b_abs_min
  
  ! for modification of sigma values for binarysplit
  IF (isw_lorentz .EQ. 1 .OR. isw_momentum .EQ. 1) THEN
     lag_sigma = 0
     ! This is for the grid-version because anumm then is not allocated
     !IF (.NOT. ALLOCATED(anumm)) ALLOCATE(anumm(0:0,0:0))
     !anumm(0,0) = 1.0_dp
     
     IF (ALLOCATED(collision_sigma_multiplier)) deallocate(collision_sigma_multiplier)
     ALLOCATE(collision_sigma_multiplier(0:0))
     collision_sigma_multiplier(0) = 1.0_dp 
  ELSE
     mult_sigma = 1.618033988749895d0
     lag_sigma = ceiling( log10( sqrt(collpar_max/collpar_min) ) / log10(mult_sigma) )
     IF (ALLOCATED(collision_sigma_multiplier)) deallocate(collision_sigma_multiplier)
     ALLOCATE(collision_sigma_multiplier(0:lag_sigma))
     mult_sigma_mod = 10**( log10( sqrt(collpar_max/collpar_min) ) / dble(lag_sigma) )
     !write (400,*) "mult_sigma_mod", mult_sigma_mod
     !write (400,*) "lag_sigma", lag_sigma
     collision_sigma_multiplier(0) = 1.0_dp
     do ilag = 1, lag_sigma
        collision_sigma_multiplier(ilag) = collision_sigma_multiplier(ilag-1) * mult_sigma_mod
        !write (400,*) "ilag", ilag, "collision_sigma_multiplier", collision_sigma_multiplier(ilag)
     end do
     !print *, collpar_min,collpar_max,sqrt(collpar_max/collpar_min)
     !print *, collision_sigma_multiplier
     !print *, ' '
  END IF

  rippletag_old = 0
  ! now go the last propagator and do the job
  iend = 0
  iendperiod = 0
  prop_ibegperiod = 1
  prop_count_call = 0
  print *, 'Setting up propagators'
  bsfunc_local_err_orig = bsfunc_local_err
  allprops: DO      ! WHILE (fieldpropagator%tag .LE. proptag_end)
     ! information about propagator
     CALL info_magnetics(fieldpropagator)
     CALL info_magnetics(fieldpropagator%parent)
     ! PAUSE
     fieldperiod => fieldpropagator%parent
     fieldripple => fieldpropagator%ch_act
     rippletag = fieldripple%tag
     proptag = fieldpropagator%tag
     WRITE(c_propagator_tag,*) proptag
     WRITE(c_ripple_tag,*) rippletag
     WRITE(c_period_tag,*) fieldperiod%tag

     write (*,*)   "Level placement for propagator", proptag

     !if (proptag .eq. 5) stop
     
     !PRINT *, 'eta_ori ',LBOUND(eta_ori)

     clear_old_ripple = 0
     ! now process ripple if it is not already processed
     ! PRINT *, 'prop_tag, ripple_tag, period_tag',fieldpropagator%tag,fieldripple%tag,fieldperiod%tag
     newripple: IF (rippletag .NE. rippletag_old) THEN
        clear_old_ripple = 1
        rippletag_old = rippletag
        fieldripple%bin_split_mode = bin_split_mode
        bin_split: IF (bin_split_mode .EQ. 1) THEN ! do the binary split           
           ! extract eta_x0 and eta_s from ripple
           CALL extract_array(fieldripple%eta_x0,eta_x0,1)
           CALL extract_array(fieldripple%eta_x0,eta_x0_hlp,1)
           CALL extract_array(fieldripple%eta_s, eta_s, 1)
           CALL extract_array(fieldripple%eta_shield, eta_shield_hlp, 1)
           IF (ALLOCATED(eta_shield)) DEALLOCATE(eta_shield)
           ALLOCATE(eta_shield(LBOUND(eta_shield_hlp,1):UBOUND(eta_shield_hlp,1)))
           eta_shield = INT(eta_shield_hlp)
           DEALLOCATE(eta_shield_hlp)
           CALL extract_array(fieldripple%eta_type, eta_type_hlp, 1)
           IF (ALLOCATED(eta_type)) DEALLOCATE(eta_type)
           ALLOCATE(eta_type(LBOUND(eta_type_hlp,1):UBOUND(eta_type_hlp,1)))
           eta_type = INT(eta_type_hlp)
           DEALLOCATE(eta_type_hlp)
           ! modify the sigma-values eta_s
           ! eta_s = eta_s * bsfunc_sigma_mult
           fix_sigma: DO i = LBOUND(eta_s,1),UBOUND(eta_s,1)
              eta_s(i) = MAX(eta_s(i),bsfunc_sigma_min)
           END DO fix_sigma

           ! get the local stuff
           IF (bsfunc_local_solver .GE. 1) THEN
              CALL extract_array(fieldripple%eta_cl,eta_cl,1)
              IF (ALLOCATED(eta_x0_loc)) DEALLOCATE(eta_x0_loc)
              IF (ALLOCATED(eta_s_loc))  DEALLOCATE(eta_s_loc)           
              IF (ALLOCATED(eta_shield_loc))  DEALLOCATE(eta_shield_loc)           
              ALLOCATE( eta_x0_loc(UBOUND(eta_cl,1)) )
              ALLOCATE( eta_s_loc(UBOUND(eta_cl,1)) )
              ALLOCATE( eta_shield_loc(UBOUND(eta_cl,1)) )
              DO i = LBOUND(eta_s,1),UBOUND(eta_cl,1)
                 eta_x0_loc(i) = eta_x0(INT(eta_cl(i)))
                 eta_x0_hlp(INT(eta_cl(i))) = 1000.0_dp
                 eta_s_loc(i)  = eta_s(INT(eta_cl(i)))
                 eta_shield_loc(i)  = eta_shield(INT(eta_cl(i)))
              END DO
              ! modify the local sigma-values eta_s_loc
              ! eta_s_loc = eta_s_loc * bsfunc_sigma_mult
              fix_sigma_local: DO i = LBOUND(eta_s_loc,1),UBOUND(eta_s_loc,1)
                 eta_s_loc(i) = MAX(eta_s_loc(i),bsfunc_sigma_min)
              END DO fix_sigma_local
              !PRINT *, 'eta_x0_loc ',eta_x0_loc
              !PRINT *, 'eta_s_loc  ',eta_s_loc              
           END IF

           ALLOCATE(eta_x0_val(1))
           ALLOCATE(eta_s_val(1))
           IF (bsfunc_local_solver .EQ. 0) THEN
              ! GLOBAL ONLY
              CALL construct_bsfunc(eta_x0,eta_s) 
              CALL construct_binarysplit(eta_ori,eta_bs)
              CALL find_binarysplit(eta_bs,eta_x0)
              CALL get_binarysplit(eta_bs,eta_split,'x')
              CALL get_binarysplit(eta_bs,fieldripple%eta,'x')

              ! store it together with the ripple
              fieldripple%eta_bs = eta_bs

              prop_reconstruct_levels = 0

           ELSEIF (bsfunc_local_solver .EQ. 1) THEN
              ! ONLY LOCAL
              prop_reconstruct_levels = 0

              ! all splitting for one level only
              CALL deconstruct_binarysplit(eta_bs)
              CALL deconstruct_binarysplit(eta_bs_loc)
              CALL construct_binarysplit(eta_ori,eta_bs_loc)

              ! Local 
              DO i_construct = 1,UBOUND(eta_x0_loc,1)
                 eta_x0_val(1) = eta_x0_loc(i_construct)
                 DO ilag = 0,lag_sigma
                    eta_s_val(1)  = eta_s_loc(i_construct) * collision_sigma_multiplier(ilag)
                    bsfunc_modelfunc = 1                                     !<-in
                    CALL construct_bsfunc(eta_x0_val,eta_s_val)               
                    CALL find_binarysplit(eta_bs_loc,eta_x0_val)
                    IF (bsfunc_sigma_mult .NE. 1.0_dp) THEN
                       eta_s_val(1)  =  eta_s_val(1) * bsfunc_sigma_mult
                       CALL construct_bsfunc(eta_x0_val,eta_s_val)               
                       CALL find_binarysplit(eta_bs_loc,eta_x0_val)
                    END IF
                 END DO
              END DO
              CALL get_binarysplit(eta_bs_loc,eta_split_loc,'x')
              CALL get_binarysplit(eta_bs_loc,fieldripple%eta_loc,'x')
              fieldripple%eta_bs_loc = eta_bs_loc
              eta_bs = eta_bs_loc
              CALL get_binarysplit(eta_bs,eta_split,'x')
              CALL get_binarysplit(eta_bs,fieldripple%eta,'x')
              fieldripple%eta_bs = eta_bs

           ELSEIF (bsfunc_local_solver .EQ. 2) THEN

              ! First LOCAL, then the ABSOLUTE MAXIMUM
              ! all splitting for one level only
              CALL deconstruct_binarysplit(eta_bs)
              CALL deconstruct_binarysplit(eta_bs_loc)
              CALL construct_binarysplit(eta_ori,eta_bs_loc)

              ! Local 
              DO i_construct = 1,UBOUND(eta_x0_loc,1)
                 eta_x0_val(1) = eta_x0_loc(i_construct)
                 DO ilag = 0,lag_sigma
                    eta_s_val(1)  = eta_s_loc(i_construct) * collision_sigma_multiplier(ilag)
                    bsfunc_modelfunc = 1                                     !<-in
                    CALL construct_bsfunc(eta_x0_val,eta_s_val)               
                    CALL find_binarysplit(eta_bs_loc,eta_x0_val)
                    IF (bsfunc_sigma_mult .NE. 1.0_dp) THEN
                       eta_s_val(1)  =  eta_s_val(1) * bsfunc_sigma_mult
                       CALL construct_bsfunc(eta_x0_val,eta_s_val)               
                       CALL find_binarysplit(eta_bs_loc,eta_x0_val)
                    END IF
                 END DO
              END DO
              CALL get_binarysplit(eta_bs_loc,eta_split_loc,'x')
              CALL get_binarysplit(eta_bs_loc,fieldripple%eta_loc,'x')
              fieldripple%eta_bs_loc = eta_bs_loc

              ! Then absolute Maximum
              eta_bs = eta_bs_loc              
              eta_min_relevant = eta_min_global + eta_part_globalfac * eta_s_relevant
              DO i_construct = 1,2
                 eta_min_loc = MINLOC(eta_x0_hlp)
                 eta_x0_hlp(eta_min_loc(1)) = 1000.0_dp
                 eta_x0_val(1) = eta_x0(eta_min_loc(1))
                 DO ilag = 0,lag_sigma
                    eta_s_val(1)  = eta_s(eta_min_loc(1)) * collision_sigma_multiplier(ilag)
                    bsfunc_modelfunc = 1
                    CALL construct_bsfunc(eta_x0_val,eta_s_val) 
                    CALL find_binarysplit(eta_bs,eta_x0_val)
                    IF (bsfunc_sigma_mult .NE. 1.0_dp) THEN
                       eta_s_val(1)  =  eta_s_val(1) * bsfunc_sigma_mult
                       CALL construct_bsfunc(eta_x0_val,eta_s_val)               
                       CALL find_binarysplit(eta_bs,eta_x0_val)
                    END IF
                 END DO
              END DO
              ! store the stuff
              CALL get_binarysplit(eta_bs,eta_split,'x')
              CALL get_binarysplit(eta_bs,fieldripple%eta,'x')
              fieldripple%eta_bs = eta_bs

           ELSEIF (bsfunc_local_solver .EQ. 3) THEN
              prop_reconstruct_levels = 0

              ! print *, '(A) bsfunc_local_solver = ',bsfunc_local_solver
              !print *, 'ripple loc', fieldripple%tag
              !PRINT *, 'eta_ori ', eta_ori
              !PRINT *, 'eta_x0_loc ', eta_x0_loc
              !PRINT *, 'eta_s_loc ', eta_s_loc
              !PRINT *, 'eta_x0 ', eta_x0
              !PRINT *, 'eta_s ', eta_s

              ! First LOCAL, then the ABSOLUTE MAXIMUM, then the REST
              ! all splitting for one level only
              CALL deconstruct_binarysplit(eta_bs)
              CALL deconstruct_binarysplit(eta_bs_loc)
              CALL construct_binarysplit(eta_ori,eta_bs_loc)

              !PRINT *, 'eta_ori ', eta_ori
              ! PRINT *, 'collision_sigma_multiplier ', collision_sigma_multiplier
              save_bsfunc_err=bsfunc_local_err
              
              ! Local 
              loc_construct: DO i_construct = 1,UBOUND(eta_x0_loc,1)
                 eta_x0_val(1) = eta_x0_loc(i_construct)
                 ! minimum local eta_s 
                 IF (i_construct .EQ. 1) THEN ! left side of ripple
                    eta_s_loc_min = ABS( &
                         fieldripple%b_max_l / fieldripple%d2bp_max_l / &
                         fieldpropagator%mdata%h_phi(LBOUND(fieldpropagator%mdata%h_phi,1))**2 &
                         )
                 ELSE ! right side of ripple
                    eta_s_loc_min = ABS( &
                         fieldripple%b_max_r / fieldripple%d2bp_max_r / &
                         fieldpropagator%mdata%h_phi(UBOUND(fieldpropagator%mdata%h_phi,1))**2 &
                         )
                 END IF
                 eta_s_loc_min = eta_x0_val(1) * SQRT(eta_s_loc_min) * collpar
                 ! minimum local eta_s - end
                 loc_laguerre: DO ilag = 0,lag_sigma
                    eta_s_val(1)  = eta_s_loc(i_construct) * collision_sigma_multiplier(ilag)

                    !**********************************************************
                    ! Additional lines for increasing eta_s_val in trapped domain
                    !**********************************************************
                    save_bsfunc_err  = bsfunc_local_err_orig
                    bsfunc_local_err = bsfunc_local_err_orig
                    eta_savemem_sigma_fact = 1d0
                    
                    if (eta_x0_val(1) .gt. (eta_b_abs_max + eta_savemem_dist1)) then
                       if (eta_x0_val(1) .lt. (eta_b_abs_max + eta_savemem_dist1 + eta_savemem_dist2)) then
                          !eta_savemem_sigma = sqrt((eta_savemem_dist2)**2 &
                          !     / log(eta_savemem_sigma_mult))
                          !eta_s_val(1) = eta_s_val(1) * &
                          !     exp(((eta_x0_val(1) - eta_b_abs_max + eta_savemem_dist1)/eta_savemem_sigma)**2)

                          !write (*,*) "Modifying sigma to ", eta_s_val(1), &
                          !     "by factor", exp(((eta_x0_val(1) - eta_b_abs_max + eta_savemem_dist1) / &
                          !     eta_savemem_sigma)**2), eta_x0_val(1) - eta_b_abs_max, eta_savemem_sigma

                          eta_savemem_sigma_fact = (eta_savemem_sigma_mult-1d0) / eta_savemem_dist2 * &
                               (eta_x0_val(1) - eta_b_abs_max - eta_savemem_dist1) + 1
                       else
                          eta_savemem_sigma_fact = eta_savemem_sigma_mult
                       end if
                    end if
                    !write (*,*) "Modifying eta_s_val by factor", eta_savemem_sigma_fact
                    !eta_s_val(1) = eta_s_val(1) * eta_savemem_sigma_fact

                    ! print *, ilag,eta_s_val(1), eta_x0_val(1)
                    loc_divide: DO idiv = 0,bsfunc_divide
                       loc_modelfunc: do ibmf = 1,bsfunc_modelfunc_num
                          bsfunc_modelfunc = ibmf
                          IF ( (eta_x0_val(1) .GE. eta_b_abs_max - bsfunc_max_mult_reach * eta_s_val(1)) .AND. &
                               (eta_x0_val(1) .LE. eta_b_abs_max + bsfunc_max_mult_reach * eta_s_val(1)) ) &
                               bsfunc_local_err = save_bsfunc_err * bsfunc_local_err_max_mult
                          ! print *, 'I am in local: ',fieldripple%tag,eta_x0_val,eta_s_val
                          call construct_bsfunc(eta_x0_val,eta_s_val)
                          call find_binarysplit(eta_bs_loc,eta_x0_val)
                          bsfunc_local_err = save_bsfunc_err
                          IF (bsfunc_sigma_mult .NE. 1.0_dp) THEN
                             eta_s_val(1)  =  eta_s_val(1) * bsfunc_sigma_mult
                             !   IF ( (eta_x0_val(1) .GE. eta_b_abs_max - bsfunc_max_mult_reach * eta_s_val(1)) .AND. &
                             !        (eta_x0_val(1) .LE. eta_b_abs_max + bsfunc_max_mult_reach * eta_s_val(1)) ) &
                             !        bsfunc_local_err = save_bsfunc_err * bsfunc_local_err_max_mult
                             CALL construct_bsfunc(eta_x0_val,eta_s_val)               
                             CALL find_binarysplit(eta_bs_loc,eta_x0_val)
                             !   bsfunc_local_err = save_bsfunc_err
                             eta_s_val(1)  =  eta_s_val(1) / bsfunc_sigma_mult
                          END IF
                       end do loc_modelfunc
                       eta_s_val(1) = eta_s_val(1) / 1.618033988749895d0
                       !IF (eta_s_val(1) .LT. eta_s_loc_min) EXIT loc_divide
                    END DO loc_divide
                    if ( eta_s_val(1) .gt. 2.0d0 * eta_x0_val(1) ) exit loc_laguerre
                 END DO loc_laguerre
                 !print *, ' '
              END DO loc_construct

              save_bsfunc_err  = bsfunc_local_err_orig
              bsfunc_local_err = bsfunc_local_err_orig
              
              CALL get_binarysplit(eta_bs_loc,eta_split_loc,'x')
              !PRINT *, 'eta_split_loc ', eta_split_loc
              !!! PRINT *, 'eta_split_loc before sqrt: ', ubound(eta_split_loc),' at ',fieldripple%tag
              

              ! WINNNY Sqrt for local lambda dependence around minimum
              if (bsfunc_lambda_loc_res) then
                 bsfunc_modelfunc_store =  bsfunc_modelfunc
                 bsfunc_modelfunc = 3
                 if (allocated(eta_x0_sqrt)) deallocate(eta_x0_sqrt)
                 allocate(eta_x0_sqrt(1))
                 if (allocated(eta_s_sqrt)) deallocate(eta_s_sqrt)
                 allocate(eta_s_sqrt(1))

                 eta_x0_sqrt(1) = 1.0d0 / fieldripple%b_min 
                 sqrt_modelfunc: DO ibmf = 1,bsfunc_modelfunc_num
                    bsfunc_modelfunc = ibmf
                    eta_s_sqrt(1) = 1/fieldripple%b_min - min(1/fieldripple%b_max_r,1/fieldripple%b_max_l)
                    if ( eta_s_sqrt(1) .ge. 0.02d0 ) then
                       CALL construct_bsfunc(eta_x0_sqrt,eta_s_sqrt)               
                       CALL find_binarysplit(eta_bs_loc,eta_x0_sqrt)
                       lambda_idiv: DO 
                          eta_s_sqrt(1) = eta_s_sqrt(1) / 1.618033988749895d0
                          if ( eta_s_sqrt(1) .lt. 0.02d0) exit lambda_idiv
                          CALL construct_bsfunc(eta_x0_sqrt,eta_s_sqrt)               
                          CALL find_binarysplit(eta_bs_loc,eta_x0_sqrt)
                       end DO lambda_idiv
                    end if
                 end DO sqrt_modelfunc
                 deallocate(eta_x0_sqrt)
                 deallocate(eta_s_sqrt)
                 bsfunc_modelfunc = bsfunc_modelfunc_store
              end if
              ! END WINNY Sqrt

              CALL get_binarysplit(eta_bs_loc,eta_split_loc,'x')
              CALL get_binarysplit(eta_bs_loc,fieldripple%eta_loc,'x')
              fieldripple%eta_bs_loc = eta_bs_loc

              !PRINT *, 'eta_split_loc ', eta_split_loc
              !PRINT *, ' ' 
              !PRINT *, 'eta_split_loc after:  ', ubound(eta_split_loc)
              !!! PRINT *, 'eta_split_loc before max:  ', ubound(eta_split_loc),' at ',fieldripple%tag
              
              ! Then absolute Maximum
              eta_bs = eta_bs_loc              
              eta_min_relevant = eta_min_global + eta_part_globalfac * eta_s_relevant
              DO i_construct = 1,2
                 eta_min_loc = MINLOC(eta_x0_hlp)
                 eta_x0_hlp(eta_min_loc(1)) = 1000.0_dp
                 eta_x0_val(1) = eta_x0(eta_min_loc(1))
                 max_laguerre : DO ilag = 0,lag_sigma
                    eta_s_val(1)  = eta_s(eta_min_loc(1)) * collision_sigma_multiplier(ilag)
                    DO ibmf = 1,bsfunc_modelfunc_num
                       bsfunc_modelfunc = ibmf
                       IF ( (eta_x0_val(1) .GE. eta_b_abs_max - bsfunc_max_mult_reach * eta_s_val(1)) .AND. &
                            (eta_x0_val(1) .LE. eta_b_abs_max + bsfunc_max_mult_reach * eta_s_val(1)) ) &
                            bsfunc_local_err = save_bsfunc_err * bsfunc_local_err_max_mult
                       ! PRINT *, 'I am in absolute: ',fieldripple%tag,eta_x0_val,eta_s_val
                       call construct_bsfunc(eta_x0_val,eta_s_val)
                       CALL find_binarysplit(eta_bs,eta_x0_val)
                       bsfunc_local_err = save_bsfunc_err
                       IF (bsfunc_sigma_mult .NE. 1.0_dp) THEN
                          eta_s_val(1)  =  eta_s_val(1) * bsfunc_sigma_mult
                          !   IF ( (eta_x0_val(1) .GE. eta_b_abs_max - bsfunc_max_mult_reach * eta_s_val(1)) .AND. &
                          !        (eta_x0_val(1) .LE. eta_b_abs_max + bsfunc_max_mult_reach * eta_s_val(1)) ) &
                          !        bsfunc_local_err = save_bsfunc_err * bsfunc_local_err_max_mult
                          CALL construct_bsfunc(eta_x0_val,eta_s_val)               
                          CALL find_binarysplit(eta_bs,eta_x0_val)
                          !   bsfunc_local_err = save_bsfunc_err
                          eta_s_val(1)  =  eta_s_val(1) / bsfunc_sigma_mult
                       END IF
                    END DO
                    if ( eta_s_val(1) .gt. 2.0d0 * eta_x0_val(1) ) exit max_laguerre
                 END DO max_laguerre
              END DO

              CALL get_binarysplit(eta_bs,eta_split,'x')
              !!! PRINT *, 'eta_split    before  rest: ', ubound(eta_split),' at ',fieldripple%tag

              ! then the rest
              if (.true.) then
                 DO i_construct = LBOUND(eta_x0,1),UBOUND(eta_x0,1)
                    IF (eta_x0_hlp(i_construct) .LT. 1000._dp) THEN
                       eta_x0_val(1) = eta_x0(i_construct)
                       rest_laguerre : DO ilag = 0,lag_sigma
                          eta_s_val(1)  = eta_s(i_construct) * collision_sigma_multiplier(ilag)
                          
                          !**********************************************************
                          ! Additional lines for increasing eta_s_val in trapped domain
                          !**********************************************************
                          save_bsfunc_err  = bsfunc_local_err_orig
                          bsfunc_local_err = bsfunc_local_err_orig
                          !**********************************************************
                          ! Additional lines for increasing eta_s_val in trapped domain
                          !**********************************************************
                          save_bsfunc_err  = bsfunc_local_err_orig
                          bsfunc_local_err = bsfunc_local_err_orig
                          eta_savemem_sigma_fact = 1d2

                          if (eta_x0_val(1) .gt. (eta_b_abs_max + eta_savemem_dist1)) then
                             if (eta_x0_val(1) .lt. (eta_b_abs_max + eta_savemem_dist1 + eta_savemem_dist2)) then
                                !eta_savemem_sigma = sqrt((eta_savemem_dist2)**2 &
                                !     / log(eta_savemem_sigma_mult))
                                !eta_s_val(1) = eta_s_val(1) * &
                                !     exp(((eta_x0_val(1) - eta_b_abs_max + eta_savemem_dist1)/eta_savemem_sigma)**2)

                                !write (*,*) "Modifying sigma to ", eta_s_val(1), &
                                !     "by factor", exp(((eta_x0_val(1) - eta_b_abs_max + eta_savemem_dist1) / &
                                !     eta_savemem_sigma)**2), eta_x0_val(1) - eta_b_abs_max, eta_savemem_sigma

                                eta_savemem_sigma_fact = (eta_savemem_sigma_mult-1d0) / eta_savemem_dist2 * &
                                     (eta_x0_val(1) - eta_b_abs_max - eta_savemem_dist1) + 1
                             else
                                eta_savemem_sigma_fact = eta_savemem_sigma_mult
                             end if
                          end if
                          !write (*,*) "Modifying eta_s_val by factor", eta_savemem_sigma_fact
                          !eta_s_val(1) = eta_s_val(1) * eta_savemem_sigma_fact/1d2

                          IF(bsfunc_ignore_trap_levels .EQ. 1 .AND. eta_x0_val(1) .GT. eta_min_relevant) CYCLE
                          DO ibmf = 1,bsfunc_modelfunc_num
                             bsfunc_modelfunc = ibmf
                             IF ( (eta_x0_val(1) .GE. eta_b_abs_max - bsfunc_max_mult_reach * eta_s_val(1)) .AND. &
                                  (eta_x0_val(1) .LE. eta_b_abs_max + bsfunc_max_mult_reach * eta_s_val(1)) ) &
                                  bsfunc_local_err = save_bsfunc_err * bsfunc_local_err_max_mult
                             ! PRINT *, 'I am in rest: ',fieldripple%tag,eta_x0_val,eta_s_val
                             CALL construct_bsfunc(eta_x0_val,eta_s_val) 
                             CALL find_binarysplit(eta_bs,eta_x0_val)
                             bsfunc_local_err = save_bsfunc_err
                             IF (bsfunc_sigma_mult .NE. 1.0_dp) THEN
                                eta_s_val(1)  =  eta_s_val(1) * bsfunc_sigma_mult
                                !   IF ( (eta_x0_val(1) .GE. eta_b_abs_max - bsfunc_max_mult_reach * eta_s_val(1)) .AND. &
                                !        (eta_x0_val(1) .LE. eta_b_abs_max + bsfunc_max_mult_reach * eta_s_val(1)) ) &
                                !        bsfunc_local_err = save_bsfunc_err * bsfunc_local_err_max_mult
                                CALL construct_bsfunc(eta_x0_val,eta_s_val)               
                                CALL find_binarysplit(eta_bs,eta_x0_val)
                                !   bsfunc_local_err = save_bsfunc_err
                                eta_s_val(1)  =  eta_s_val(1) / bsfunc_sigma_mult
                             END IF
                          END DO
                          if ( eta_s_val(1) .gt. 2.0d0 * eta_x0_val(1) ) exit rest_laguerre
                       END DO rest_laguerre
                    END IF
                 END DO
              end if

              bsfunc_local_err=save_bsfunc_err
              ! store the stuff
              CALL get_binarysplit(eta_bs,eta_split,'x')
              CALL get_binarysplit(eta_bs,fieldripple%eta,'x')
              fieldripple%eta_bs = eta_bs
              !!! PRINT *, 'eta_split     after  rest: ', ubound(eta_split),' at ',fieldripple%tag
              !!! PRINT *, ' '

              !PRINT *, 'eta_split_loc ', eta_split_loc
              !PRINT *, 'eta_split ', eta_split
              !PRINT *, 'ubound loc ',UBOUND(eta_split_loc,1)
              !PRINT *, 'ubound     ',UBOUND(eta_split,1)
              !PRINT *, ' '
              !print *, 'End - bsfunc_local_solver = ',bsfunc_local_solver
              !PAUSE



           ELSEIF (bsfunc_local_solver .EQ. 4) THEN
              ! This is the mode with the new shielded exit (prepared for in mag_interface)
              ! Actually it is a copy of bsfunc_local_solver = 3 
              ! with shielding for bootstrap current
              prop_reconstruct_levels = 0
              bsfunc_modelfunc = 1
              bsfunc_base_distance = MAXVAL(eta_ori(1:UBOUND(eta_ori,1))-eta_ori(0:UBOUND(eta_ori,1)-1))
              !print *, 'bsfunc_base_distance ',bsfunc_base_distance
              ! First LOCAL, then the ABSOLUTE MAXIMUM, then the REST
              ! all splitting for one level only
              CALL deconstruct_binarysplit(eta_bs)
              CALL deconstruct_binarysplit(eta_bs_loc)
              CALL construct_binarysplit(eta_ori,eta_bs_loc)

              !print *, 'tag    ',fieldripple%tag
              save_binarysplit_fsplitdepth = binarysplit_fsplitdepth
              binarysplit_fsplitdepth = 3
              !PRINT *, 'eta_ori'
              !PRINT *, eta_ori
              save_bsfunc_err = bsfunc_local_err
              save_bsfunc_local_err = bsfunc_local_err
              ! Local 
              i_construct_s = 1
              i_construct_e = UBOUND(eta_x0_loc,1)
              i_construct_d = 1
              loc_construct4: DO i_construct = i_construct_s,i_construct_e,i_construct_d
                 eta_x0_val(1) = eta_x0_loc(i_construct)
                 IF (eta_shield_loc(i_construct) .EQ. 2) THEN
                    bsfunc_local_err = save_bsfunc_local_err * bsfunc_local_shield_factor
                 ELSE
                    bsfunc_local_err = save_bsfunc_local_err
                 END IF
                 loc_laguerre4: DO ilag = 0,lag_sigma
                    eta_s_val(1)  = eta_s_loc(i_construct) * collision_sigma_multiplier(ilag)
                    !print *, 'loca   ',fieldripple%tag,eta_x0_val,eta_s_val,bsfunc_local_err
                    CALL multiple_binarysplit(eta_bs_loc,eta_x0_val,eta_s_val)
                 END DO loc_laguerre4
                 bsfunc_local_err = save_bsfunc_local_err
              END DO loc_construct4
              !CALL printsummary_binarysplit(eta_bs_loc)

              CALL get_binarysplit(eta_bs_loc,eta_split_loc,'x')
              CALL get_binarysplit(eta_bs_loc,fieldripple%eta_loc,'x')
              fieldripple%eta_bs_loc = eta_bs_loc
              eta_bs = eta_bs_loc              


              !print *, 'bsfunc_local_err ',bsfunc_local_err

              !PRINT *, 'eta_split_loc'
              !PRINT *, eta_split_loc

!!$              ! Then absolute Maximum
!!$              eta_min_relevant = eta_min_global + eta_part_globalfac * eta_s_relevant
!!$              DO i_construct = 1,UBOUND(eta_x0,1)
!!$                 eta_min_loc = MINLOC(eta_x0_hlp)
!!$                 eta_x0_val(1) = eta_x0(eta_min_loc(1))
!!$                 if (i_construct .eq. 1) eta_trapped_passing=eta_x0_val(1) !<=NEW
!!$                 if (eta_x0_val(1) .gt. eta_trapped_passing + 1.d-11) exit
!!$                 eta_x0_hlp(eta_min_loc(1)) = 1000.0_dp
!!$                 !sigma_trapped_passing=eta_s(eta_min_loc(1)) !<=NEW
!!$                 DO ilag = 0,lag_sigma
!!$                    eta_s_val(1)  = eta_s(eta_min_loc(1)) * collision_sigma_multiplier(ilag)
!!$                    !print *, 'amax   ',eta_x0_val,eta_s_val,bsfunc_local_err
!!$                    DO ibmf = bsfunc_modelfunc_num,1,-1
!!$                       bsfunc_modelfunc = ibmf
!!$                       CALL construct_bsfunc(eta_x0_val,eta_s_val) 
!!$                       CALL find_binarysplit(eta_bs,eta_x0_val)
!!$                       !bsfunc_local_err = save_bsfunc_err
!!$                       !IF (bsfunc_sigma_mult .NE. 1.0_dp) THEN
!!$                       !   eta_s_val(1)  =  eta_s_val(1) * bsfunc_sigma_mult
!!$                       !   CALL construct_bsfunc(eta_x0_val,eta_s_val)               
!!$                       !   CALL find_binarysplit(eta_bs,eta_x0_val)
!!$                       !   bsfunc_local_err = save_bsfunc_err
!!$                       !   eta_s_val(1)  =  eta_s_val(1) / bsfunc_sigma_mult
!!$                       END IF
!!$                    END DO
!!$                 END DO
!!$              END DO

              ! then the rest
              DO i_construct = LBOUND(eta_x0,1),UBOUND(eta_x0,1)
                 IF (eta_type(i_construct) .EQ. 4) CYCLE ! inflection level
                 IF (eta_type(i_construct) .EQ. 5) CYCLE ! local inflection level
                 IF (eta_type(i_construct) .EQ. 6) CYCLE ! intersection with boundary
                 IF (eta_x0_hlp(i_construct) .GE. 1000._dp) CYCLE
                 eta_x0_val(1) = eta_x0(i_construct)
                 DO ilag = 0,lag_sigma
                    eta_s_val(1)  = eta_s(i_construct) * collision_sigma_multiplier(ilag)
                    IF (eta_shield(i_construct) .EQ. 1 .AND. bsfunc_shield) CYCLE !<= NEW Winny
                    !print *, 'rest   ',eta_x0_val,eta_s_val,bsfunc_local_err
                    ! IF(bsfunc_ignore_trap_levels .EQ. 1 .AND. eta_x0_val(1) .GT. eta_min_relevant) CYCLE
                    CALL multiple_binarysplit(eta_bs,eta_x0_val,eta_s_val)
!!$                    DO ibmf = bsfunc_modelfunc_num,1,-1
!!$                       bsfunc_modelfunc = ibmf
!!$                       CALL construct_bsfunc(eta_x0_val,eta_s_val) 
!!$                       CALL find_binarysplit(eta_bs,eta_x0_val)
!!$                       !bsfunc_local_err = save_bsfunc_err
!!$                       !IF (bsfunc_sigma_mult .NE. 1.0_dp) THEN
!!$                       !   eta_s_val(1)  =  eta_s_val(1) * bsfunc_sigma_mult
!!$                       !   CALL construct_bsfunc(eta_x0_val,eta_s_val)               
!!$                       !   CALL find_binarysplit(eta_bs,eta_x0_val)
!!$                       !   !   bsfunc_local_err = save_bsfunc_err
!!$                       !   eta_s_val(1)  =  eta_s_val(1) / bsfunc_sigma_mult
!!$                       !END IF
!!$                    END DO
                 END DO
              END DO
              bsfunc_local_err = save_bsfunc_local_err

              ! inflection levels
              bsfunc_local_err = save_bsfunc_local_err * bsfunc_local_shield_factor
              DO i_construct = LBOUND(eta_x0,1),UBOUND(eta_x0,1)
                 IF (eta_type(i_construct) .NE. 4) CYCLE ! no inflection level
                 eta_x0_val(1) = eta_x0(i_construct)
                 DO ilag = 0,lag_sigma
                    eta_s_val(1)  = eta_s(i_construct) * collision_sigma_multiplier(ilag)
                    IF (eta_shield(i_construct) .EQ. 1 .AND. bsfunc_shield) CYCLE !<= NEW Winny
                    !print *, 'inf    ',eta_x0_val,eta_s_val,bsfunc_local_err
                    CALL multiple_binarysplit(eta_bs,eta_x0_val,eta_s_val)
                    !DO ibmf = bsfunc_modelfunc_num,1,-1
                    !   bsfunc_modelfunc = ibmf
                    !   CALL construct_bsfunc(eta_x0_val,eta_s_val) 
                    !   CALL find_binarysplit(eta_bs,eta_x0_val)
                    !END DO
                 END DO
              END DO
              bsfunc_local_err = save_bsfunc_local_err

              ! local inflection levels
              bsfunc_local_err = save_bsfunc_local_err * bsfunc_local_shield_factor
              DO i_construct = LBOUND(eta_x0,1),UBOUND(eta_x0,1)
                 IF (eta_type(i_construct) .NE. 5) CYCLE ! no local inflection level
                 eta_x0_val(1) = eta_x0(i_construct)
                 DO ilag = 0,lag_sigma
                    eta_s_val(1)  = eta_s(i_construct) * collision_sigma_multiplier(ilag)
                    !print *, 'infl   ',eta_x0_val,eta_s_val,bsfunc_local_err
                    !print *, 'local inflection ',eta_s_val
                    IF (eta_shield(i_construct) .EQ. 1 .AND. bsfunc_shield) CYCLE !<= NEW Winny
                    CALL multiple_binarysplit(eta_bs,eta_x0_val,eta_s_val)
                    !DO ibmf = bsfunc_modelfunc_num,1,-1
                    !   bsfunc_modelfunc = ibmf
                    !   CALL construct_bsfunc(eta_x0_val,eta_s_val) 
                    !   CALL find_binarysplit(eta_bs,eta_x0_val)
                    !END DO
                 END DO
              END DO
              bsfunc_local_err = save_bsfunc_local_err

              ! local intersection levels
              bsfunc_local_err = save_bsfunc_local_err * bsfunc_local_shield_factor
              DO i_construct = LBOUND(eta_x0,1),UBOUND(eta_x0,1)
                 IF (eta_type(i_construct) .NE. 6) CYCLE ! no local inflection level
                 ! print *, 'flint inter ',fieldripple%tag
                 eta_x0_val(1) = eta_x0(i_construct)
                 DO ilag = 0,lag_sigma
                    eta_s_val(1)  = eta_s(i_construct) * collision_sigma_multiplier(ilag)
                    !print *, 'intl   ',eta_x0_val,eta_s_val,bsfunc_local_err
                    !print *, eta_s_val(1)
                    CALL multiple_binarysplit(eta_bs,eta_x0_val,eta_s_val)
                    !bsfunc_modelfunc = 2
                    !CALL construct_bsfunc(eta_x0_val,eta_s_val) 
                    !CALL find_binarysplit(eta_bs,eta_x0_val)
                 END DO
              END DO
              bsfunc_local_err = save_bsfunc_local_err


              bsfunc_local_err=save_bsfunc_err
              ! store the stuff
              CALL get_binarysplit(eta_bs,eta_split,'x')
              CALL get_binarysplit(eta_bs,fieldripple%eta,'x')
              fieldripple%eta_bs = eta_bs

              !pause

           ELSEIF (bsfunc_local_solver .EQ. 5) THEN
              ! This is the mode with the shielded exit by Sergei
              ! Actually it is a copy of bsfunc_local_solver = 3 
              ! with shielding for bootstrap current
              prop_reconstruct_levels = 0

              ! First LOCAL, then the ABSOLUTE MAXIMUM, then the REST
              ! all splitting for one level only
              CALL deconstruct_binarysplit(eta_bs)
              CALL deconstruct_binarysplit(eta_bs_loc)
              CALL construct_binarysplit(eta_ori,eta_bs_loc)

              !PRINT *, 'eta_ori'
              !PRINT *, eta_ori

              save_bsfunc_err=bsfunc_local_err 
              ! Local 
              loc_construct5: DO i_construct = 1,UBOUND(eta_x0_loc,1)
                 eta_x0_val(1) = eta_x0_loc(i_construct)
                 ! minimum local eta_s 
                 IF (i_construct .EQ. 1) THEN ! left side of ripple
                    eta_s_loc_min = ABS( &
                         fieldripple%b_max_l / fieldripple%d2bp_max_l / &
                         fieldpropagator%mdata%h_phi(LBOUND(fieldpropagator%mdata%h_phi,1))**2 &
                         )
                    eta_highest_local_max=eta_x0_val(1) !<=NEW
                 ELSE ! right side of ripple
                    eta_s_loc_min = ABS( &
                         fieldripple%b_max_r / fieldripple%d2bp_max_r / &
                         fieldpropagator%mdata%h_phi(UBOUND(fieldpropagator%mdata%h_phi,1))**2 &
                         )
                    eta_highest_local_max=MIN(eta_highest_local_max,eta_x0_val(1)) !<=NEW
                 END IF
                 eta_s_loc_min = eta_x0_val(1) * SQRT(eta_s_loc_min) * collpar
                 ! minimum local eta_s - end
                 loc_laguerre5: DO ilag = 0,lag_sigma
                    eta_s_val(1)  = eta_s_loc(i_construct) * collision_sigma_multiplier(ilag)
                    loc_divide5: DO idiv = 0,bsfunc_divide
                       loc_modelfunc5: DO ibmf = 1,bsfunc_modelfunc_num
                          bsfunc_modelfunc = ibmf
                          IF ( (eta_x0_val(1) .GE. eta_b_abs_max - bsfunc_max_mult_reach * eta_s_val(1)) .AND. &
                               (eta_x0_val(1) .LE. eta_b_abs_max + bsfunc_max_mult_reach * eta_s_val(1)) ) &
                               bsfunc_local_err = save_bsfunc_err * bsfunc_local_err_max_mult
                          CALL construct_bsfunc(eta_x0_val,eta_s_val)               
                          CALL find_binarysplit(eta_bs_loc,eta_x0_val)
                          bsfunc_local_err = save_bsfunc_err
                          IF (bsfunc_sigma_mult .NE. 1.0_dp) THEN
                             eta_s_val(1)  =  eta_s_val(1) * bsfunc_sigma_mult
                             !   IF ( (eta_x0_val(1) .GE. eta_b_abs_max - bsfunc_max_mult_reach * eta_s_val(1)) .AND. &
                             !        (eta_x0_val(1) .LE. eta_b_abs_max + bsfunc_max_mult_reach * eta_s_val(1)) ) &
                             !        bsfunc_local_err = save_bsfunc_err * bsfunc_local_err_max_mult
                             CALL construct_bsfunc(eta_x0_val,eta_s_val)               
                             CALL find_binarysplit(eta_bs_loc,eta_x0_val)
                             !   bsfunc_local_err = save_bsfunc_err
                             eta_s_val(1)  =  eta_s_val(1) / bsfunc_sigma_mult
                          END IF
                       END DO loc_modelfunc5
                       eta_s_val(1) = eta_s_val(1) / 1.618033988749895d0
                       !                       IF (eta_s_val(1) .LT. eta_s_loc_min) EXIT loc_divide
                    END DO loc_divide5
                 END DO loc_laguerre5
              END DO loc_construct5
              CALL get_binarysplit(eta_bs_loc,eta_split_loc,'x')
              CALL get_binarysplit(eta_bs_loc,fieldripple%eta_loc,'x')
              fieldripple%eta_bs_loc = eta_bs_loc

              !PRINT *, 'eta_split_loc'
              !PRINT *, eta_split_loc

              ! Then absolute Maximum
              eta_bs = eta_bs_loc              
              eta_min_relevant = eta_min_global + eta_part_globalfac * eta_s_relevant
              DO i_construct = 1,2
                 eta_min_loc = MINLOC(eta_x0_hlp)
                 eta_x0_hlp(eta_min_loc(1)) = 1000.0_dp
                 eta_x0_val(1) = eta_x0(eta_min_loc(1))
                 eta_trapped_passing=eta_x0_val(1) !<=NEW
                 sigma_trapped_passing=eta_s(eta_min_loc(1)) !<=NEW
                 DO ilag = 0,lag_sigma
                    eta_s_val(1)  = eta_s(eta_min_loc(1)) * collision_sigma_multiplier(ilag)
                    DO ibmf = 1,bsfunc_modelfunc_num
                       bsfunc_modelfunc = ibmf
                       IF ( (eta_x0_val(1) .GE. eta_b_abs_max - bsfunc_max_mult_reach * eta_s_val(1)) .AND. &
                            (eta_x0_val(1) .LE. eta_b_abs_max + bsfunc_max_mult_reach * eta_s_val(1)) ) &
                            bsfunc_local_err = save_bsfunc_err * bsfunc_local_err_max_mult
                       !PRINT *, 'I am in absolute'
                       CALL construct_bsfunc(eta_x0_val,eta_s_val) 
                       CALL find_binarysplit(eta_bs,eta_x0_val)
                       bsfunc_local_err = save_bsfunc_err
                       IF (bsfunc_sigma_mult .NE. 1.0_dp) THEN
                          eta_s_val(1)  =  eta_s_val(1) * bsfunc_sigma_mult
                          !   IF ( (eta_x0_val(1) .GE. eta_b_abs_max - bsfunc_max_mult_reach * eta_s_val(1)) .AND. &
                          !        (eta_x0_val(1) .LE. eta_b_abs_max + bsfunc_max_mult_reach * eta_s_val(1)) ) &
                          !        bsfunc_local_err = save_bsfunc_err * bsfunc_local_err_max_mult
                          CALL construct_bsfunc(eta_x0_val,eta_s_val)               
                          CALL find_binarysplit(eta_bs,eta_x0_val)
                          !   bsfunc_local_err = save_bsfunc_err
                          eta_s_val(1)  =  eta_s_val(1) / bsfunc_sigma_mult
                       END IF
                    END DO
                 END DO
              END DO

              ! then the rest
              DO i_construct = LBOUND(eta_x0,1),UBOUND(eta_x0,1)
                 IF (eta_x0_hlp(i_construct) .LT. 1000._dp) THEN
                    eta_x0_val(1) = eta_x0(i_construct)
                    DO ilag = 0,lag_sigma
                       eta_s_val(1)  = eta_s(i_construct) * collision_sigma_multiplier(ilag)
                       ! new exit by Sergei
                       IF(ABS(eta_x0_val(1)-eta_trapped_passing).GT. & !<=NEW
                            3.d0*(sigma_trapped_passing+eta_s_val(1))) THEN !<=NEW
                          CYCLE !<=NEW
                       ELSEIF(eta_x0_val(1)-3.d0*eta_s_val(1).LT.eta_highest_local_max) THEN !<=NEW
                          CYCLE !<=NEW
                       ENDIF !<=NEW
                       ! end - new exit by Sergei                       
                       IF(bsfunc_ignore_trap_levels .EQ. 1 .AND. eta_x0_val(1) .GT. eta_min_relevant) CYCLE
                       DO ibmf = 1,bsfunc_modelfunc_num
                          bsfunc_modelfunc = ibmf
                          IF ( (eta_x0_val(1) .GE. eta_b_abs_max - bsfunc_max_mult_reach * eta_s_val(1)) .AND. &
                               (eta_x0_val(1) .LE. eta_b_abs_max + bsfunc_max_mult_reach * eta_s_val(1)) ) &
                               bsfunc_local_err = save_bsfunc_err * bsfunc_local_err_max_mult
                          !PRINT *, 'I am in rest'
                          CALL construct_bsfunc(eta_x0_val,eta_s_val) 
                          CALL find_binarysplit(eta_bs,eta_x0_val)
                          bsfunc_local_err = save_bsfunc_err
                          IF (bsfunc_sigma_mult .NE. 1.0_dp) THEN
                             eta_s_val(1)  =  eta_s_val(1) * bsfunc_sigma_mult
                             !   IF ( (eta_x0_val(1) .GE. eta_b_abs_max - bsfunc_max_mult_reach * eta_s_val(1)) .AND. &
                             !        (eta_x0_val(1) .LE. eta_b_abs_max + bsfunc_max_mult_reach * eta_s_val(1)) ) &
                             !        bsfunc_local_err = save_bsfunc_err * bsfunc_local_err_max_mult
                             CALL construct_bsfunc(eta_x0_val,eta_s_val)               
                             CALL find_binarysplit(eta_bs,eta_x0_val)
                             !   bsfunc_local_err = save_bsfunc_err
                             eta_s_val(1)  =  eta_s_val(1) / bsfunc_sigma_mult
                          END IF
                       END DO
                    END DO
                 END IF
              END DO

              bsfunc_local_err=save_bsfunc_err
              ! store the stuff
              CALL get_binarysplit(eta_bs,eta_split,'x')
              CALL get_binarysplit(eta_bs,fieldripple%eta,'x')
              fieldripple%eta_bs = eta_bs

           ELSEIF (bsfunc_local_solver .EQ. 11) THEN
              prop_reconstruct_levels = 0
              ! only historical was 1
              ! FIRST LOCAL
              ! create info for gauss
              CALL construct_bsfunc(eta_x0_loc,eta_s_loc) 

              CALL construct_binarysplit(eta_ori,eta_bs_loc)
              ! this does the split with knowledge about maxima
              !  without eta_x0 it would do the split without this knowledge

              CALL find_binarysplit(eta_bs_loc,eta_x0_loc) 
              ! this gets the new eta-split into an array and into fieldripple%eta

              CALL get_binarysplit(eta_bs_loc,eta_split_loc,'x')
              CALL get_binarysplit(eta_bs_loc,fieldripple%eta_loc,'x')

              ! THEN GLOBAL
              CALL construct_bsfunc(eta_x0,eta_s) 
              eta_bs = eta_bs_loc
              CALL find_binarysplit(eta_bs,eta_x0)
              CALL get_binarysplit(eta_bs,eta_split,'x')
              CALL get_binarysplit(eta_bs,fieldripple%eta,'x')

              ! store it together with the ripple
              fieldripple%eta_bs = eta_bs
              fieldripple%eta_bs_loc = eta_bs_loc
           ELSEIF (bsfunc_local_solver .EQ. 12) THEN
              prop_reconstruct_levels = 0
              ! only historical - was 2
              ! ! First LOCAL, then GLOBAL (LOCAL IGNORED AFTERWARDS)
              ! CALL construct_bsfunc(eta_x0_loc,eta_s_loc)               
              ! CALL construct_binarysplit(eta_ori,eta_bs_loc)
              ! CALL find_binarysplit(eta_bs_loc,eta_x0_loc) 
              ! eta_bs = eta_bs_loc
              ! CALL construct_bsfunc(eta_x0,eta_s) 
              ! CALL find_binarysplit(eta_bs,eta_x0)
              ! CALL get_binarysplit(eta_bs,eta_split,'x')
              ! CALL get_binarysplit(eta_bs,fieldripple%eta,'x')
              ! fieldripple%eta_bs = eta_bs



              ! First LOCAL, then the absolute maximum, then the rest
              ! all splitting for one level only

              CALL construct_binarysplit(eta_ori,eta_bs_loc)

              DO i_construct = 1,UBOUND(eta_x0_loc,1)
                 eta_x0_val(1) = eta_x0_loc(i_construct)
                 eta_s_val(1)  = eta_s_loc(i_construct)
                 CALL construct_bsfunc(eta_x0_val,eta_s_val)               
                 CALL find_binarysplit(eta_bs_loc,eta_x0_val) 
              END DO

              CALL get_binarysplit(eta_bs_loc,eta_split_loc,'x')
              CALL get_binarysplit(eta_bs_loc,fieldripple%eta_loc,'x')
              fieldripple%eta_bs_loc = eta_bs_loc

              ! Then absolute Maximum
              eta_bs = eta_bs_loc

              DO i_construct = 1,2
                 eta_min_loc = MINLOC(eta_x0_hlp)
                 eta_x0_hlp(eta_min_loc(1)) = 1000.0_dp
                 eta_x0_val(1) = eta_x0(eta_min_loc(1))
                 eta_s_val(1)  = eta_s(eta_min_loc(1))
                 CALL construct_bsfunc(eta_x0_val,eta_s_val) 
                 CALL find_binarysplit(eta_bs,eta_x0_val)
              END DO

              ! then the rest
              DO i_construct = LBOUND(eta_x0,1),UBOUND(eta_x0,1)
                 IF (eta_x0_hlp(i_construct) .LT. 1000._dp) THEN
                    eta_x0_val(1) = eta_x0(i_construct)
                    eta_s_val(1)  = eta_s(i_construct)
                    CALL construct_bsfunc(eta_x0_val,eta_s_val) 
                    CALL find_binarysplit(eta_bs,eta_x0_val)
                 END IF
              END DO
              ! store the stuff
              CALL get_binarysplit(eta_bs,eta_split,'x')
              CALL get_binarysplit(eta_bs,fieldripple%eta,'x')
              fieldripple%eta_bs = eta_bs
           ELSEIF (bsfunc_local_solver .EQ. 13) THEN
              prop_reconstruct_levels = 0
              ! was 3
              ! First LOCAL, then the absolute maximum, then the rest
              ! all splitting for one level only

              CALL deconstruct_binarysplit(eta_bs)
              CALL deconstruct_binarysplit(eta_bs_loc)
              CALL construct_binarysplit(eta_ori,eta_bs_loc)

              ! local 
              DO i_construct = 1,UBOUND(eta_x0_loc,1)
                 eta_x0_val(1) = eta_x0_loc(i_construct)
                 DO ilag = 0,lag_sigma
                    eta_s_val(1)  = eta_s_loc(i_construct) * collision_sigma_multiplier(ilag)
                    bsfunc_modelfunc = 1                                     !<-in
                    CALL construct_bsfunc(eta_x0_val,eta_s_val)               
                    CALL find_binarysplit(eta_bs_loc,eta_x0_val)
                    eta_s_val(1)  =  eta_s_val(1) * 2.0_dp
                    CALL construct_bsfunc(eta_x0_val,eta_s_val)               
                    CALL find_binarysplit(eta_bs_loc,eta_x0_val)
!!$                    ! WINNY - I switch this off for a test
!!$                    bsfunc_modelfunc = 2                                     !<-in
!!$                    CALL construct_bsfunc(eta_x0_val,eta_s_val)              !<-in 
!!$                    CALL find_binarysplit(eta_bs_loc,eta_x0_val)             !<-in
!!$                    bsfunc_modelfunc = 3                                     !<-in
!!$                    CALL construct_bsfunc(eta_x0_val,eta_s_val)              !<-in 
!!$                    CALL find_binarysplit(eta_bs_loc,eta_x0_val)             !<-in
!!$                    ! WINNY - end

                 END DO
              END DO
              CALL get_binarysplit(eta_bs_loc,eta_split_loc,'x')
              CALL get_binarysplit(eta_bs_loc,fieldripple%eta_loc,'x')
              fieldripple%eta_bs_loc = eta_bs_loc

              ! Then absolute Maximum
              eta_bs = eta_bs_loc

              eta_min_relevant = eta_min_global + 3.d0 * eta_s_relevant
              ! could be
              !eta_min_relevant = eta_min_global + eta_part_globalfac * eta_s_relevant

              ! eta_min_relevant=1000.d0                                    !<-in
              DO i_construct = 1,2
                 eta_min_loc = MINLOC(eta_x0_hlp)
                 eta_x0_hlp(eta_min_loc(1)) = 1000.0_dp
                 eta_x0_val(1) = eta_x0(eta_min_loc(1))
                 DO ilag = 0,lag_sigma
                    eta_s_val(1)  = eta_s(eta_min_loc(1)) * collision_sigma_multiplier(ilag)
                    !eta_min_relevant=MIN(eta_min_relevant,                &  !<-in
                    !                                      eta_x0_val(1)+2.d0*eta_s_val(1))    !<-in
                    !                     eta_x0_val(1)+3.d0*eta_s_val(1))    !<-in
                    bsfunc_modelfunc = 1                                     !<-in
                    !CALL construct_bsfunc(eta_x0_val,eta_s_val) 
                    !CALL find_binarysplit(eta_bs,eta_x0_val)
!!$                    ! WINNY - I switch this off for a test                    
!!$                    bsfunc_modelfunc = 2                                     !<-in
!!$                    CALL construct_bsfunc(eta_x0_val,eta_s_val)              !<-in 
!!$                    CALL find_binarysplit(eta_bs,eta_x0_val)             !<-in ERROR
!!$                    bsfunc_modelfunc = 3                                     !<-in
!!$                    CALL construct_bsfunc(eta_x0_val,eta_s_val)              !<-in 
!!$                    CALL find_binarysplit(eta_bs,eta_x0_val)             !<-in ERROR
!!$                    ! WINNY - end
                 END DO
              END DO

              ! then the rest
              save_bsfunc_err=bsfunc_local_err                            !<-in
              DO i_construct = LBOUND(eta_x0,1),UBOUND(eta_x0,1)
                 IF (eta_x0_hlp(i_construct) .LT. 1000._dp) THEN
                    eta_x0_val(1) = eta_x0(i_construct)
                    DO ilag = 0,lag_sigma
                       eta_s_val(1)  = eta_s(i_construct) * collision_sigma_multiplier(ilag)
                       IF(eta_x0_val(1) .GT.                               &  !<-in
                            eta_min_relevant) CYCLE                            !<-in
                       bsfunc_modelfunc = 1                                  !<-in
                       !CALL construct_bsfunc(eta_x0_val,eta_s_val) 
                       !CALL find_binarysplit(eta_bs,eta_x0_val)
!!$                       ! WINNY - I switch this off for a test
!!$                       bsfunc_modelfunc = 2                                  !<-in
!!$                       CALL construct_bsfunc(eta_x0_val,eta_s_val)           !<-in 
!!$                       CALL find_binarysplit(eta_bs,eta_x0_val)          !<-in ERROR
!!$                       bsfunc_modelfunc = 3                                  !<-in
!!$                       CALL construct_bsfunc(eta_x0_val,eta_s_val)           !<-in 
!!$                       CALL find_binarysplit(eta_bs,eta_x0_val)          !<-in ERROR
!!$                       ! WINNY - end
                    END DO
                 END IF
              END DO
              bsfunc_local_err=save_bsfunc_err                            !<-in
              ! store the stuff
              CALL get_binarysplit(eta_bs,eta_split,'x')
              CALL get_binarysplit(eta_bs,fieldripple%eta,'x')
              fieldripple%eta_bs = eta_bs
           ELSE
              PRINT *, 'bsfunc_local_solver = ',bsfunc_local_solver,' not implemented'
              STOP
           END IF
           DEALLOCATE(eta_x0_val,eta_s_val)

!!$           ! EXPERIMENT
!!$           IF (proptag .EQ. proptag_start) THEN
!!$              eta_bs_store = eta_bs

!!$           END IF
!!$           IF (proptag .EQ. proptag_end) THEN
!!$              eta_bs = eta_bs_store
!!$              CALL get_binarysplit(eta_bs,eta_split,'x')
!!$              CALL get_binarysplit(eta_bs,fieldripple%eta,'x')
!!$              fieldripple%eta_bs = eta_bs
!!$           END IF
!!$           ! EXPERIMENT



           ! if one is interested in plots of eta distributions
           ! 
           !    plot "gauss_func.dat" u 1:2 w l, "gauss_split.dat" u 1:2 w p
           !
           !  program stops with pause
           ! information about ripple

           !CALL printsummary_binarysplit(eta_bs)
           !PRINT *, '----------------------------------------'

           IF (plot_gauss .EQ. 1) THEN
              CALL info_magnetics(fieldripple)
              !CALL printsummary_binarysplit(eta_bs)
              PRINT *, '---------------------------'
              WRITE(c_filename,'(100A)') 'etaori.dat'
              OPEN(1000,file=c_filename)
              DO i = LBOUND(eta_ori,1),UBOUND(eta_ori,1)
                 WRITE(1000,*) eta_ori(i)
              END DO
              CLOSE(1000)
              WRITE(c_filename,'(100A)') c_etapropori,TRIM(ADJUSTL(c_propagator_tag)),c_extprop
              OPEN(1000,file=c_filename)
              DO i = LBOUND(eta_ori,1),UBOUND(eta_ori,1)
                 WRITE(1000,*) eta_ori(i)
              END DO
              CLOSE(1000)

              OPEN(1000,file='etadis.dat')
              DO i = LBOUND(eta_x0,1),UBOUND(eta_x0,1)
                 !PRINT *, eta_x0(i),eta_s(i)
                 WRITE(1000,*) eta_x0(i),eta_s(i)
              END DO
              CLOSE(1000)

              WRITE(c_filename,'(100A)') 'etasplit.dat'
              OPEN(1001,file=c_filename)
              DO i = LBOUND(eta_split,1),UBOUND(eta_split,1)
                 WRITE(1001,*) eta_split(i)
              END DO
              CLOSE(1001)
              WRITE(c_filename,'(100A)') c_etapropsplit,TRIM(ADJUSTL(c_propagator_tag)),c_extprop
              OPEN(1001,file=c_filename)
              DO i = LBOUND(eta_split,1),UBOUND(eta_split,1)
                 WRITE(1001,*) eta_split(i)
              END DO
              CLOSE(1001)

              IF (bsfunc_local_solver .EQ. 1 .OR. bsfunc_local_solver .EQ. 2) THEN
                 OPEN(1001,file='etasplitloc.dat')
                 DO i = LBOUND(eta_split_loc,1),UBOUND(eta_split_loc,1)
                    WRITE(1001,*) eta_split_loc(i)
                 END DO
                 CLOSE(1001)
              ELSE
                 OPEN(1001,file='etasplitloc.dat')
                 CLOSE(1001)                 
              END IF
              PRINT *, '---------------------------'
              CALL plotfile_binarysplit(eta_bs,'gauss_func.dat',9,1000)
              CALL plotfile_binarysplit(eta_bs,'gauss_split.dat',9)
              CALL plotfile_binarysplit(eta_bs,'gauss_inter.dat',9,'c')
              PRINT *, 'plotfile_binarysplit finished!'
              !PAUSE
           END IF
        ELSE
           ! put eta_ori in fieldripple
           IF (ALLOCATED(fieldripple%eta)) DEALLOCATE(fieldripple%eta)
           ALLOCATE(fieldripple%eta(0:UBOUND(eta_ori,1)))
           fieldripple%eta = eta_ori
           prop_reconstruct_levels = 0
        END IF bin_split
        ! deallocation of binarysplit
        CALL deconstruct_binarysplit(eta_bs)
        CALL deconstruct_binarysplit(eta_bs_loc)


!!$        PRINT *, 'eta_ori', eta_ori
!!$        OPEN(1000,file='eta_ori.dat')
!!$        WRITE(1000,*) eta_ori
!!$        CLOSE(1000)
!!$        OPEN(1000,file='eta_split.dat')
!!$        WRITE(1000,*) eta_split
!!$        CLOSE(1000)
!!$        OPEN(1000,file='eta_x0.dat')
!!$        WRITE(1000,*) eta_x0
!!$        CLOSE(1000)



     END IF newripple


     IF (fieldpropagator%tag .EQ. proptag_last) THEN 
        iend = 1
        iendperiod = 1
     ELSE
        IF (ASSOCIATED(fieldpropagator%next)) THEN
           IF (fieldpropagator%parent%tag .NE. fieldpropagator%next%parent%tag) THEN 
              iendperiod = 1
           ELSE
              iendperiod = 0
           END IF
        ELSE
           iendperiod = 1
        END IF
     END IF

     ! PRINT *, 'End of Tag'
     ! go to the next propagator or exit
     IF (fieldpropagator%tag .EQ. proptag_last) EXIT allprops
     IF (.NOT.(ASSOCIATED(fieldpropagator%next))) THEN
        start_at_begin = 1
        fieldpropagator => fieldline%ch_fir%ch_fir
     ELSE
        IF (fieldpropagator%next%tag .LE. fieldline%ch_las%ch_las%tag) THEN
           start_at_begin = 0
           fieldpropagator => fieldpropagator%next
        ELSE
           start_at_begin = 1
           fieldpropagator => fieldline%ch_fir%ch_fir
        END IF
     END IF
  END DO allprops
  PRINT *, 'Setting up propagators - End'
  ! PAUSE

!!$  IF (mag_save_memory .EQ. 1) THEN
!!$     CALL deconstruct_binarysplit(fieldripple%eta_bs)
!!$     CALL deconstruct_binarysplit(fieldripple%eta_bs_loc)
!!$  END IFg

  IF (bin_split_mode .NE. 0) THEN
     ! check for the boundary problem
     boundary_fix_mode = 3
     !boundary_fix_counter_max = 50
     IF (boundary_fix_mode .EQ. 2) THEN ! this was the global solution which failed
        PRINT *, '--------------------------------------------'
        boundaryfixcounter: DO boundary_fix_counter = 1,boundary_fix_counter_max+1
           PRINT *, 'Setting up propagators for boundary check:',' boundary_fix_counter ',boundary_fix_counter
           fieldperiod => fieldline%ch_fir
           fieldpropagator => fieldperiod%ch_fir
           fieldripple => fieldpropagator%ch_act
           first_ripple = .TRUE.
           allripplesetadetect: DO      
              rippletag = fieldripple%tag
              ! left
              IF (first_ripple .AND. boundary_fix_counter .EQ. 1) THEN
                 CALL set_new(pointer_all_bmax,fieldripple%b_max_l)
                 CALL set_new(pointer_ripple_all_bmax,DBLE(rippletag))
                 first_ripple = .FALSE.
              END IF
              DO i = LBOUND(fieldripple%eta,1)+1,UBOUND(fieldripple%eta,1)-1
                 boundary_dist = 1.0d0 / fieldripple%b_max_l - fieldripple%eta(i)
                 previous_dist = fieldripple%eta(i) - fieldripple%eta(i-1)
                 next_dist = fieldripple%eta(i+1) - fieldripple%eta(i)
                 boundary_dist_limit = next_dist * boundary_dist_limit_factor
                 IF (boundary_dist .GE. 0.0d0 .AND. boundary_dist .LE. boundary_dist_limit) THEN
                    !if (abs(boundary_dist) .le. boundary_dist_limit) then
                    PRINT *, ' left:   ','rippletag: ',rippletag, &
                         ' proptag ',fieldripple%pa_fir%tag,' - ' ,fieldripple%pa_las%tag ,' level' ,i
                    PRINT *, ' 1.0d0/fieldripple%b_max_l ',1.0d0/fieldripple%b_max_l
                    PRINT *, ' fieldripple%eta(i)        ',fieldripple%eta(i)
                    PRINT *, ' boundary_dist             ',boundary_dist
                    PRINT *, ' boundary_dist_limit       ',boundary_dist_limit
                    CALL set_new(pointer_index_eta_close_to_boundary,DBLE(i))
                    CALL set_new(pointer_boundary_dist,boundary_dist)
                    CALL set_new(pointer_eta_move,2.0d0 * boundary_dist_limit)
                    CALL set_new(pointer_next_dist,next_dist)
                    CALL set_new(pointer_boundary_dist_limit,boundary_dist_limit)
                    CALL set_new(pointer_eta_close_to_boundary,fieldripple%eta(i))
                    CALL set_new(pointer_ripple_close_to_boundary,DBLE(rippletag))
                 END IF
              END DO
              ! right
              first_ripple = .FALSE.
              IF (boundary_fix_counter .EQ. 1) THEN
                 CALL set_new(pointer_all_bmax,fieldripple%b_max_r)
                 CALL set_new(pointer_ripple_all_bmax,DBLE(rippletag))
              END IF
              DO i = LBOUND(fieldripple%eta,1)+1,UBOUND(fieldripple%eta,1)-1
                 boundary_dist = 1.0d0 / fieldripple%b_max_r - fieldripple%eta(i)
                 previous_dist = fieldripple%eta(i) - fieldripple%eta(i-1)
                 next_dist = fieldripple%eta(i+1) - fieldripple%eta(i)
                 boundary_dist_limit = next_dist * boundary_dist_limit_factor
                 IF (boundary_dist .GE. 0.0d0 .AND. boundary_dist .LE. boundary_dist_limit) THEN
                    !if (abs(boundary_dist) .le. boundary_dist_limit) then
                    PRINT *, ' right:  ','rippletag: ',rippletag, &
                         ' proptag ',fieldripple%pa_fir%tag,' - ',fieldripple%pa_las%tag ,' level' ,i
                    PRINT *, ' 1.0d0/fieldripple%b_max_r ',1.0d0/fieldripple%b_max_r
                    PRINT *, ' fieldripple%eta(i)        ',fieldripple%eta(i)
                    PRINT *, ' boundary_dist             ',boundary_dist
                    PRINT *, ' boundary_dist_limit       ',boundary_dist_limit
                    CALL set_new(pointer_index_eta_close_to_boundary,DBLE(i))
                    CALL set_new(pointer_boundary_dist,boundary_dist)
                    CALL set_new(pointer_eta_move,2.0d0 * boundary_dist_limit)
                    CALL set_new(pointer_next_dist,next_dist)
                    CALL set_new(pointer_boundary_dist_limit,boundary_dist_limit)
                    CALL set_new(pointer_eta_close_to_boundary,fieldripple%eta(i))
                    CALL set_new(pointer_ripple_close_to_boundary,DBLE(rippletag))
                 END IF
              END DO

              IF ( ASSOCIATED(fieldripple%next) ) THEN
                 fieldripple => fieldripple%next
              ELSE
                 EXIT allripplesetadetect
              END IF
           END DO allripplesetadetect

           CALL extract_array(pointer_eta_close_to_boundary,array_eta_close_to_boundary)
           CALL extract_array(pointer_all_bmax,array_all_bmax)
           CALL extract_array(pointer_ripple_close_to_boundary,array_ripple_close_to_boundary)
           CALL extract_array(pointer_ripple_all_bmax,array_ripple_all_bmax)
           CALL extract_array(pointer_index_eta_close_to_boundary,array_index_eta_close_to_boundary)
           CALL extract_array(pointer_boundary_dist,array_boundary_dist)
           CALL extract_array(pointer_boundary_dist_limit,array_boundary_dist_limit)
           CALL extract_array(pointer_eta_move,array_eta_move)
           CALL extract_array(pointer_next_dist,array_next_dist)

           IF (ALLOCATED(array_eta_close_to_boundary)) THEN
              boundary_has_to_be_fixed = .TRUE.
           ELSE
              boundary_has_to_be_fixed = .FALSE.
           END IF

           IF (boundary_fix_counter .EQ. 1) THEN
              PRINT *, array_all_bmax
           END IF

           CALL delete_all(pointer_eta_close_to_boundary)
           CALL delete_all(pointer_all_bmax)
           CALL delete_all(pointer_ripple_close_to_boundary)
           CALL delete_all(pointer_ripple_all_bmax)
           CALL delete_all(pointer_index_eta_close_to_boundary)
           CALL delete_all(pointer_boundary_dist)
           CALL delete_all(pointer_boundary_dist_limit)
           CALL delete_all(pointer_eta_move)
           CALL delete_all(pointer_next_dist)
           PRINT *, 'Setting up propagators for boundary check - End'
           ! check for the boundary problem - end

           IF (boundary_has_to_be_fixed .AND. boundary_fix_counter .GT. boundary_fix_counter_max) THEN
              PRINT *, ''
              PRINT *, 'eta conflicts with boundary could not be solved!'
              STOP
           END IF

           IF (.NOT. boundary_has_to_be_fixed) THEN
              PRINT *, ' all conflicts are fixed'
              EXIT boundaryfixcounter
           END IF
!!$        print *, 'array_all_bmax'
!!$        print *, array_all_bmax
!!$        print *, 'array_ripple_all_bmax'
!!$        print *, int(array_ripple_all_bmax)
!!$        print *, 'array_eta_close_to_boundary'
!!$        if (boundary_has_to_be_fixed) then
!!$           print *, array_eta_close_to_boundary
!!$           print *, 'array_ripple_close_to_boundary'
!!$           print *, int(array_ripple_close_to_boundary)
!!$           print *, 'array_index_eta_close_to_boundary'
!!$           print *, int(array_index_eta_close_to_boundary)
!!$           print *, 'array_boundary_dist'
!!$           print *, array_boundary_dist
!!$           print *, 'array_boundary_dist_limit'
!!$           print *, array_boundary_dist_limit
!!$           print *, 'array_eta_move'
!!$           print *, array_eta_move
!!$        end if

           ! fix the boundary problem
           ! go to the first ripple
           PRINT *, 'Fixing the boundary problem'
           fieldperiod => fieldline%ch_fir 
           fieldpropagator => fieldperiod%ch_fir
           fieldripple => fieldpropagator%ch_act
           ! which etas have to be processed
           ALLOCATE(array_eta_used(LBOUND(array_eta_close_to_boundary,1):UBOUND(array_eta_close_to_boundary,1)))
           array_eta_used = .TRUE.
           DO b_loop_1 = LBOUND(array_eta_close_to_boundary,1),UBOUND(array_eta_close_to_boundary,1)
              IF (array_eta_used(b_loop_1)) THEN
                 DO b_loop_2 = LBOUND(array_eta_close_to_boundary,1),UBOUND(array_eta_close_to_boundary,1)
                    IF (array_eta_used(b_loop_2) .AND. b_loop_1 .NE. b_loop_2) THEN
                       IF (array_eta_close_to_boundary(b_loop_1) .EQ. array_eta_close_to_boundary(b_loop_2)) THEN
                          !print *, 'same eta found'
                          array_eta_used(b_loop_2) = .FALSE.
                          array_eta_move(b_loop_1) = MAX(array_eta_move(b_loop_1),array_eta_move(b_loop_2))
                          ! for rare cases here should be a check whether this for sure causes no conflict
                          ! not necessary because there is now a check in the loop
                       END IF
                    END IF
                 END DO
              END IF
           END DO
           !print *, array_eta_used
           !print *, array_eta_close_to_boundary
           !print *, array_eta_move
           !pause
           ! here now all etas which are in conflict with boundaries should be known
           ! they have to be moved now in all ripples
           boundary_eta_fix_counter = 0
           allripplesfixetabound: DO
              rippletag = fieldripple%tag
              DO b_loop_1 = LBOUND(array_eta_close_to_boundary,1),UBOUND(array_eta_close_to_boundary,1)
                 IF (array_eta_used(b_loop_1)) THEN
                    DO b_loop_2 = LBOUND(fieldripple%eta,1),UBOUND(fieldripple%eta,1)
                       IF (array_eta_close_to_boundary(b_loop_1) .EQ. fieldripple%eta(b_loop_2)) THEN
                          !print *, 'eta ',array_eta_close_to_boundary(b_loop_1),' found in ',rippletag
                          boundary_eta_fix_counter = boundary_eta_fix_counter + 1
                          fieldripple%eta(b_loop_2) = fieldripple%eta(b_loop_2) + array_eta_move(b_loop_1)
                       END IF
                    END DO
                 END IF
              END DO
              IF ( ASSOCIATED(fieldripple%next) ) THEN
                 fieldripple => fieldripple%next
              ELSE
                 EXIT  allripplesfixetabound
              END IF
           END DO allripplesfixetabound
           ! deallocate arrays
           IF (ALLOCATED(array_eta_close_to_boundary)) DEALLOCATE(array_eta_close_to_boundary)
           IF (ALLOCATED(array_ripple_close_to_boundary)) DEALLOCATE(array_ripple_close_to_boundary)
           IF (ALLOCATED(array_index_eta_close_to_boundary)) DEALLOCATE(array_index_eta_close_to_boundary)
           IF (ALLOCATED(array_boundary_dist)) DEALLOCATE(array_boundary_dist)
           IF (ALLOCATED(array_eta_move)) DEALLOCATE(array_eta_move)
           IF (ALLOCATED(array_next_dist)) DEALLOCATE(array_next_dist)
           IF (ALLOCATED(array_boundary_dist_limit)) DEALLOCATE(array_boundary_dist_limit)
           IF (ALLOCATED(array_eta_used)) DEALLOCATE(array_eta_used)
           PRINT *, ' boundary_eta_fix_counter = ',boundary_eta_fix_counter
           PRINT *, 'Fixing the boundary problem - End'
           PRINT *, '--------------------------------------------'        
        END DO boundaryfixcounter

        IF (ALLOCATED(array_all_bmax)) DEALLOCATE(array_all_bmax)
        IF (ALLOCATED(array_ripple_all_bmax)) DEALLOCATE(array_ripple_all_bmax)

        PRINT *, '--------------------------------------------'        

     END IF

     IF (boundary_fix_mode .EQ. 3) THEN ! new stuff - this is the movement in minus direction
        ! check for the boundary problem
        ! go to the first propagator which is wanted
        fieldperiod => fieldline%ch_fir 
        fieldpropagator => fieldperiod%ch_fir
        DO WHILE (fieldpropagator%tag .LT. proptag_first)
           IF (.NOT. ASSOCIATED(fieldpropagator%next)) EXIT
           fieldpropagator => fieldpropagator%next
        END DO
        PRINT *, '--------------------------------------------'
        PRINT *, 'Setting up propagators for boundary check'
        allpropsbound3: DO      
           fieldperiod => fieldpropagator%parent
           fieldripple => fieldpropagator%ch_act
           rippletag = fieldripple%tag
           proptag = fieldpropagator%tag
           ! leftmost propagator in ripple
           IF (proptag .EQ. fieldripple%pa_fir%tag) THEN
              DO i = LBOUND(fieldripple%eta,1)+1,UBOUND(fieldripple%eta,1)-1
                 boundary_dist = 1.0d0 / fieldripple%b_max_l - fieldripple%eta(i)
                 previous_dist = fieldripple%eta(i) - fieldripple%eta(i-1)
                 next_dist = fieldripple%eta(i+1) - fieldripple%eta(i)
                 boundary_dist_limit = previous_dist * boundary_dist_limit_factor
                 IF (boundary_dist .GE. 0.0d0 .AND. boundary_dist .LE. boundary_dist_limit) THEN
                    !print *, ' left:   ','rippletag: ',rippletag, &
                    !     ' proptag ',fieldripple%pa_fir%tag,' - ',fieldripple%pa_las%tag,' level' ,i
                    !print *, '         ',fieldripple%eta(i),1.0d0/fieldripple%b_max_l,boundary_dist
                    !print *, '         ',boundary_dist,boundary_dist_limit
                    fieldripple%eta_boundary_left = fieldripple%eta(i)
                    fieldripple%eta_boundary_modification_left = boundary_dist_limit
                    fieldripple%eta_boundary_index_left = i
                 END IF
              END DO
           END IF

           ! rightmost propagator in ripple
           IF (proptag .EQ. fieldripple%pa_las%tag) THEN
              DO i = LBOUND(fieldripple%eta,1)+1,UBOUND(fieldripple%eta,1)-1
                 boundary_dist = 1.0d0 / fieldripple%b_max_r - fieldripple%eta(i)
                 previous_dist = fieldripple%eta(i) - fieldripple%eta(i-1)
                 next_dist = fieldripple%eta(i+1) - fieldripple%eta(i)
                 boundary_dist_limit = previous_dist * boundary_dist_limit_factor
                 IF (boundary_dist .GE. 0.0d0 .AND. boundary_dist .LE. boundary_dist_limit) THEN
                    !print *, ' right:  ','rippletag: ',rippletag, &
                    !     ' proptag ',fieldripple%pa_fir%tag,' - ',fieldripple%pa_las%tag,' level' ,i
                    !print *, '          ',fieldripple%eta(i),1.0d0/fieldripple%b_max_r,boundary_dist
                    !print *, '          ',boundary_dist,boundary_dist_limit
                    fieldripple%eta_boundary_right = fieldripple%eta(i)
                    fieldripple%eta_boundary_modification_right = boundary_dist_limit
                    fieldripple%eta_boundary_index_right = i
                 END IF
              END DO
           END IF

           IF (fieldpropagator%tag .EQ. proptag_last) THEN 
              iend = 1
              iendperiod = 1
           ELSE
              IF (ASSOCIATED(fieldpropagator%next)) THEN
                 IF (fieldpropagator%parent%tag .NE. fieldpropagator%next%parent%tag) THEN 
                    iendperiod = 1
                 ELSE
                    iendperiod = 0
                 END IF
              ELSE
                 iendperiod = 1
              END IF
           END IF

           IF (fieldpropagator%tag .EQ. proptag_last) EXIT allpropsbound3
           IF (.NOT.(ASSOCIATED(fieldpropagator%next))) THEN
              start_at_begin = 1
              fieldpropagator => fieldline%ch_fir%ch_fir
           ELSE
              IF (fieldpropagator%next%tag .LE. fieldline%ch_las%ch_las%tag) THEN
                 start_at_begin = 0
                 fieldpropagator => fieldpropagator%next
              ELSE
                 start_at_begin = 1
                 fieldpropagator => fieldline%ch_fir%ch_fir
              END IF
           END IF
        END DO allpropsbound3

        PRINT *, 'Setting up propagators for boundary check - End'
        ! check for the boundary problem - end

        ! fix the boundary problem
        ! go to the first ripple
        fieldperiod => fieldline%ch_fir 
        fieldpropagator => fieldperiod%ch_fir
        fieldripple => fieldpropagator%ch_act
        PRINT *, 'Fixing the boundary problem'

        allripplesbound3: DO      
           rippletag = fieldripple%tag

           opp_limit_safe = 1.e-4
           IF ( fieldripple%eta_boundary_index_left .NE. 0 .AND. fieldripple%eta_boundary_index_right .EQ. 0) THEN
              ! only left side affected
              eb_ind_left = fieldripple%eta_boundary_index_left
              move_eta = fieldripple%eta_boundary_modification_left
              old_eta = fieldripple%eta(eb_ind_left)
              new_eta = old_eta - move_eta
              opp_limit = 1.0d0 / fieldripple%b_max_r
              IF (old_eta .GE. opp_limit .AND. new_eta .LT. opp_limit) THEN ! problem on other side
                 new_eta = opp_limit + opp_limit_safe
                 IF (new_eta .LT. old_eta) THEN
                    fieldripple%eta(eb_ind_left) = new_eta
                    boundary_counter_partly_fixed = boundary_counter_partly_fixed + 1
                    PRINT *, ' ripple ',fieldripple%tag,' propagator ',fieldripple%pa_fir%tag,' left  - partly fixed'
                 ELSE
                    boundary_counter_not_fixed = boundary_counter_not_fixed + 1
                    PRINT *, ' ripple ',fieldripple%tag,' propagator ',fieldripple%pa_fir%tag,' left  - not fixed'
                 END IF
              ELSE ! no problem on other side
                 fieldripple%eta(eb_ind_left) = new_eta
                 boundary_counter_fixed = boundary_counter_fixed + 1
                 PRINT *, ' ripple ',fieldripple%tag,' propagator ',fieldripple%pa_fir%tag,' left  - fixed'
              END IF
           ELSEIF ( fieldripple%eta_boundary_index_left .EQ. 0 .AND. fieldripple%eta_boundary_index_right .NE. 0) THEN
              ! only right side affected
              eb_ind_right = fieldripple%eta_boundary_index_right
              move_eta = fieldripple%eta_boundary_modification_right
              old_eta = fieldripple%eta(eb_ind_right)
              new_eta = old_eta - move_eta
              opp_limit = 1.0d0 / fieldripple%b_max_l
              IF (old_eta .GE. opp_limit .AND. new_eta .LT. opp_limit) THEN ! problem on other side
                 new_eta = opp_limit + opp_limit_safe
                 IF (new_eta .LT. old_eta) THEN
                    fieldripple%eta(eb_ind_right) = new_eta
                    boundary_counter_partly_fixed = boundary_counter_partly_fixed + 1
                    PRINT *, ' ripple ',fieldripple%tag,' propagator ',fieldripple%pa_las%tag,' right - partly fixed'
                 ELSE
                    boundary_counter_not_fixed = boundary_counter_not_fixed + 1
                    PRINT *, ' ripple ',fieldripple%tag,' propagator ',fieldripple%pa_las%tag,' right - not fixed'
                 END IF
              ELSE ! no problem on other side
                 fieldripple%eta(eb_ind_right) = new_eta
                 boundary_counter_fixed = boundary_counter_fixed + 1
                 PRINT *, ' ripple ',fieldripple%tag,' propagator ',fieldripple%pa_las%tag,' right - fixed'
              END IF
           ELSEIF ( fieldripple%eta_boundary_index_left .NE. 0 .AND. fieldripple%eta_boundary_index_right .NE. 0) THEN
              ! both sides affected
              eb_ind_left = fieldripple%eta_boundary_index_left
              eb_ind_right = fieldripple%eta_boundary_index_right
              IF (eb_ind_left .EQ. eb_ind_left) THEN ! same index on both sides
                 move_eta = MAX(fieldripple%eta_boundary_modification_right,fieldripple%eta_boundary_modification_left)
                 new_eta = fieldripple%eta(eb_ind_right) - move_eta
                 fieldripple%eta(eb_ind_right) = new_eta
                 boundary_counter_fixed = boundary_counter_fixed + 1
                 PRINT *, ' ripple ',fieldripple%tag,' propagator ', &
                      fieldripple%pa_fir%tag,'-',fieldripple%pa_las%tag,' left,right - fixed'
              ELSE ! different indices
                 ! right
                 move_eta = fieldripple%eta_boundary_modification_right
                 old_eta = fieldripple%eta(eb_ind_right)
                 new_eta = old_eta - move_eta
                 opp_limit = 1.0d0 / fieldripple%b_max_l
                 IF (old_eta .GE. opp_limit .AND. new_eta .LT. opp_limit) THEN ! problem on other side
                    new_eta = opp_limit + opp_limit_safe
                    IF (new_eta .LT. old_eta) THEN
                       fieldripple%eta(eb_ind_right) = new_eta
                       boundary_counter_partly_fixed = boundary_counter_partly_fixed + 1
                       PRINT *, ' ripple ',fieldripple%tag,' propagator ',fieldripple%pa_las%tag,' right - partly fixed'
                    ELSE
                       boundary_counter_not_fixed = boundary_counter_not_fixed + 1
                       PRINT *, ' ripple ',fieldripple%tag,' propagator ',fieldripple%pa_las%tag,' right - not fixed'
                    END IF
                 ELSE ! no problem on other side
                    fieldripple%eta(eb_ind_right) = new_eta
                    boundary_counter_fixed = boundary_counter_fixed + 1
                    PRINT *, ' ripple ',fieldripple%tag,' propagator ',fieldripple%pa_las%tag,' right - fixed'
                 END IF
                 ! left
                 move_eta = fieldripple%eta_boundary_modification_left
                 old_eta = fieldripple%eta(eb_ind_left)
                 new_eta = old_eta - move_eta
                 opp_limit = 1.0d0 / fieldripple%b_max_r
                 IF (old_eta .GE. opp_limit .AND. new_eta .LT. opp_limit) THEN ! problem on other side
                    new_eta = opp_limit + opp_limit_safe
                    IF (new_eta .LT. old_eta) THEN
                       fieldripple%eta(eb_ind_left) = new_eta
                       boundary_counter_partly_fixed = boundary_counter_partly_fixed + 1
                       PRINT *, ' ripple ',fieldripple%tag,' propagator ',fieldripple%pa_fir%tag,' left  - partly fixed'
                    ELSE
                       boundary_counter_not_fixed = boundary_counter_not_fixed + 1
                       PRINT *, ' ripple ',fieldripple%tag,' propagator ',fieldripple%pa_fir%tag,' left  - not fixed'
                    END IF
                 ELSE ! no problem on other side
                    fieldripple%eta(eb_ind_left) = new_eta
                    boundary_counter_fixed = boundary_counter_fixed + 1
                    PRINT *, ' ripple ',fieldripple%tag,' propagator ',fieldripple%pa_fir%tag,' left  - fixed'
                 END IF
              END IF
           END IF

           IF ( ASSOCIATED(fieldripple%next) ) THEN
              fieldripple => fieldripple%next
           ELSE
              EXIT allripplesbound3
           END IF

        END DO allripplesbound3
        PRINT *, ' boundary_counter_fixed:        ',boundary_counter_fixed
        PRINT *, ' boundary_counter_partly_fixed: ',boundary_counter_partly_fixed
        PRINT *, ' boundary_counter_not_fixed:    ',boundary_counter_not_fixed
        PRINT *, 'Fixing the boundary problem - End'
        PRINT *, '--------------------------------------------'
     END IF
     ! End of new code
     ! fix the boundary problem - end

     IF (boundary_fix_mode .EQ. 1) THEN ! old stuff - this is the movement in plus direction
        ! check for the boundary problem
        ! go to the first propagator which is wanted
        fieldperiod => fieldline%ch_fir 
        fieldpropagator => fieldperiod%ch_fir
        DO WHILE (fieldpropagator%tag .LT. proptag_first)
           IF (.NOT. ASSOCIATED(fieldpropagator%next)) EXIT
           fieldpropagator => fieldpropagator%next
        END DO
        PRINT *, ' '
        PRINT *, '--------------------------------------------'
        PRINT *, 'Setting up propagators for boundary check'
        allpropsbound: DO      
           fieldperiod => fieldpropagator%parent
           fieldripple => fieldpropagator%ch_act
           rippletag = fieldripple%tag
           proptag = fieldpropagator%tag
           ! leftmost propagator in ripple
           IF (proptag .EQ. fieldripple%pa_fir%tag) THEN
              DO i = LBOUND(fieldripple%eta,1)+1,UBOUND(fieldripple%eta,1)-1
                 boundary_dist = 1.0d0 / fieldripple%b_max_l - fieldripple%eta(i)
                 previous_dist = fieldripple%eta(i) - fieldripple%eta(i-1)
                 next_dist = fieldripple%eta(i+1) - fieldripple%eta(i)
                 boundary_dist_limit = next_dist * boundary_dist_limit_factor
                 IF (boundary_dist .GE. 0.0d0 .AND. boundary_dist .LE. boundary_dist_limit) THEN
                    !if (abs(boundary_dist) .le. boundary_dist_limit) then
                    PRINT *, 'left:   ','rippletag: ',rippletag, &
                         ' proptag ',fieldripple%pa_fir%tag,' - ',fieldripple%pa_las%tag,' level' ,i
                    PRINT *, '        ',fieldripple%eta(i),1.0d0/fieldripple%b_max_l,boundary_dist
!!$              open(9890,file='eta_left.dat')
!!$              write(9890,*)  fieldripple%eta
!!$              close(9890)
!!$              open(9891,file='eta_prev.dat')
!!$              write(9891,*)  fieldripple%prev%eta
!!$              close(9891)
                    fieldripple%eta_boundary_left = fieldripple%eta(i)
                    fieldripple%eta_boundary_modification_left = 2.0d0 * boundary_dist_limit
                    fieldripple%eta_boundary_index_left = i
                    !fieldripple%eta(i) = fieldripple%eta(i) + 2.0d0 * boundary_dist_limit
                    !boundary_dist = 1.0d0 / fieldripple%b_max_l - fieldripple%eta(i)
                    !print *, 'leftc:  ',rippletag,proptag,i,fieldripple%eta(i),1.0d0/fieldripple%b_max_l,boundary_dist
!!$              pause
                 END IF
              END DO
           END IF

           ! rightmost propagator in ripple
           IF (proptag .EQ. fieldripple%pa_las%tag) THEN
              DO i = LBOUND(fieldripple%eta,1)+1,UBOUND(fieldripple%eta,1)-1
                 boundary_dist = 1.0d0 / fieldripple%b_max_r - fieldripple%eta(i)
                 previous_dist = fieldripple%eta(i) - fieldripple%eta(i-1)
                 next_dist = fieldripple%eta(i+1) - fieldripple%eta(i)
                 boundary_dist_limit = next_dist * boundary_dist_limit_factor
                 IF (boundary_dist .GE. 0.0d0 .AND. boundary_dist .LE. boundary_dist_limit) THEN
                    !if (abs(boundary_dist) .le. boundary_dist_limit) then
                    PRINT *, 'right:  ','rippletag: ',rippletag, &
                         ' proptag ',fieldripple%pa_fir%tag,' - ',fieldripple%pa_las%tag,' level' ,i
                    PRINT *, '        ',fieldripple%eta(i),1.0d0/fieldripple%b_max_l,boundary_dist
!!$              open(9990,file='eta_right.dat')
!!$              write(9990,*)  fieldripple%eta
!!$              close(9990)
!!$              open(9991,file='eta_next.dat')
!!$              write(9991,*)  fieldripple%next%eta
!!$              close(9991)
                    fieldripple%eta_boundary_right = fieldripple%eta(i)
                    fieldripple%eta_boundary_modification_right = 2.0d0 * boundary_dist_limit
                    fieldripple%eta_boundary_index_right = i
                    !fieldripple%eta(i) = fieldripple%eta(i) + 2.0d0 * boundary_dist_limit
                    !boundary_dist = 1.0d0 / fieldripple%b_max_r - fieldripple%eta(i)
                    !print *, 'rightc: ',rippletag,proptag,i,fieldripple%eta(i),1.0d0/fieldripple%b_max_r,boundary_dist
!!$              pause
                 END IF
              END DO
           END IF

           IF (fieldpropagator%tag .EQ. proptag_last) THEN 
              iend = 1
              iendperiod = 1
           ELSE
              IF (ASSOCIATED(fieldpropagator%next)) THEN
                 IF (fieldpropagator%parent%tag .NE. fieldpropagator%next%parent%tag) THEN 
                    iendperiod = 1
                 ELSE
                    iendperiod = 0
                 END IF
              ELSE
                 iendperiod = 1
              END IF
           END IF

           IF (fieldpropagator%tag .EQ. proptag_last) EXIT allpropsbound
           IF (.NOT.(ASSOCIATED(fieldpropagator%next))) THEN
              start_at_begin = 1
              fieldpropagator => fieldline%ch_fir%ch_fir
           ELSE
              IF (fieldpropagator%next%tag .LE. fieldline%ch_las%ch_las%tag) THEN
                 start_at_begin = 0
                 fieldpropagator => fieldpropagator%next
              ELSE
                 start_at_begin = 1
                 fieldpropagator => fieldline%ch_fir%ch_fir
              END IF
           END IF
        END DO allpropsbound

        PRINT *, '--------------------------------------------'
        PRINT *, 'Setting up propagators for boundary check - End'
!!$  pause
        ! check for the boundary problem - end

        ! fix the boundary problem
        ! go to the first ripple
        fieldperiod => fieldline%ch_fir 
        fieldpropagator => fieldperiod%ch_fir
        fieldripple => fieldpropagator%ch_act
        PRINT *, ' '
        PRINT *, '--------------------------------------------'
        PRINT *, 'Fixing the boundary problem'

        allripplesbound: DO      
           rippletag = fieldripple%tag

           ! left and right ripples
           IF ( ASSOCIATED(fieldripple%prev) ) THEN
              fieldripple_left => fieldripple%prev
           ELSE
              fieldripple_left => fieldline%ch_las%ch_las%ch_act ! last
           END IF
           IF ( ASSOCIATED(fieldripple%next) ) THEN
              fieldripple_right => fieldripple%next
           ELSE
              fieldripple_right => fieldline%ch_fir%ch_fir%ch_act ! first
           END IF

           ! left
           IF ( fieldripple%eta_boundary_index_left .NE. 0 ) THEN
              ! left side affected
              PRINT *, 'fieldripple ',fieldripple%tag,' fieldpropagator ',fieldripple%pa_fir%tag
              IF ( fieldripple%eta_boundary_index_left .EQ. fieldripple%eta_boundary_index_right ) THEN
                 ! both sides equal
                 IF ( fieldripple%eta_boundary_left .EQ. fieldripple_left%eta_boundary_right &
                      .AND. &
                      fieldripple%eta_boundary_right .EQ. fieldripple_right%eta_boundary_left & 
                      ) THEN
                    ! can be moved
                    eb_ind = fieldripple%eta_boundary_index_left
                    fieldripple%eta(eb_ind) = fieldripple%eta(eb_ind) + fieldripple%eta_boundary_modification_left
                    boundary_counter_fixed = boundary_counter_fixed + 1
                    PRINT *, 'fieldripple ',fieldripple%tag,' left,right - fixed'
                 ELSE
                    ! can not be moved
                    boundary_counter_not_fixed = boundary_counter_not_fixed + 1
                    PRINT *, 'fieldripple ',fieldripple%tag,' left,right - not fixed'
                 END IF
              ELSE
                 ! only left side
                 IF ( fieldripple%eta_boundary_left .EQ. fieldripple_left%eta_boundary_right ) THEN
                    ! can be moved
                    eb_ind = fieldripple%eta_boundary_index_left
                    fieldripple%eta(eb_ind) = fieldripple%eta(eb_ind) + fieldripple%eta_boundary_modification_left
                    boundary_counter_fixed = boundary_counter_fixed + 1
                    PRINT *, 'fieldripple ',fieldripple%tag,' left - fixed'
                 ELSE
                    ! can not be moved
                    boundary_counter_not_fixed = boundary_counter_not_fixed + 1
                    PRINT *, 'fieldripple ',fieldripple%tag,' left - not fixed'
                 END IF
              END IF
           END IF

           ! right
           IF ( fieldripple%eta_boundary_index_right .NE. 0 &
                .AND. &
                fieldripple%eta_boundary_index_right .NE. fieldripple%eta_boundary_index_left &
                ) THEN
              ! right side affected
              PRINT *, 'fieldripple ',fieldripple%tag,' fieldpropagator ',fieldripple%pa_las%tag
              IF ( fieldripple%eta_boundary_right .EQ. fieldripple_right%eta_boundary_left ) THEN
                 ! can be moved
                 eb_ind = fieldripple%eta_boundary_index_right
                 fieldripple%eta(eb_ind) = fieldripple%eta(eb_ind) + fieldripple%eta_boundary_modification_right
                 boundary_counter_fixed = boundary_counter_fixed + 1
                 PRINT *, 'fieldripple ',fieldripple%tag,' right - fixed'
              ELSE
                 ! can not be moved
                 boundary_counter_not_fixed = boundary_counter_not_fixed + 1
                 PRINT *, 'fieldripple ',fieldripple%tag,' right - not fixed'
              END IF
           END IF

           IF ( ASSOCIATED(fieldripple%next) ) THEN
              fieldripple => fieldripple%next
           ELSE
              EXIT allripplesbound
           END IF

        END DO allripplesbound
        PRINT *, ' '
        PRINT *, 'boundary_counter_fixed:     ',boundary_counter_fixed
        PRINT *, 'boundary_counter_not_fixed: ',boundary_counter_not_fixed
        PRINT *, ' '
        PRINT *, 'Fixing the boundary problem - End'
        PRINT *, '--------------------------------------------'
     END IF
     ! End of now the best code
     ! fix the boundary problem - end
  END IF ! bin_split_mode .ne. 0
  ! Now the magentic field in the propagators is fixed

!!$    ! Testing of Ripple and print
!!$    print *, 'write all eta levels'
!!$    OPEN(unit=5001,file='ripple_level.dat')
!!$    OPEN(unit=5002,file='ripple_level_num.dat')
!!$    fieldperiod => fieldline%ch_fir
!!$    fieldpropagator => fieldperiod%ch_fir
!!$    fieldripple => fieldpropagator%ch_act
!!$    DO
!!$       !print *, 'ripple ',fieldripple%tag
!!$       write(5001,*) fieldripple%tag
!!$       do i = lbound(fieldripple%eta,1),ubound(fieldripple%eta,1)
!!$          write(5001,*) fieldripple%eta(i)
!!$       end do
!!$       write(5002,*) fieldripple%tag,ubound(fieldripple%eta,1)
!!$       IF(ASSOCIATED(fieldripple%next)) THEN
!!$          fieldripple => fieldripple%next
!!$       ELSE
!!$          EXIT
!!$       END IF
!!$    end DO
!!$    close(unit=5001)
!!$    close(unit=5002)
!!$    print *, 'all eta levels written'
!!$
!!$    ! go to the first propagator which is wanted
!!$    fieldperiod => fieldline%ch_fir 
!!$    fieldpropagator => fieldperiod%ch_fir
!!$    OPEN(unit=5001,file='propagator_tags.dat')
!!$    do
!!$       write (5001,*) fieldpropagator%parent%tag,fieldpropagator%tag,fieldpropagator%ch_act%tag, &
!!$            fieldpropagator%phi_l,fieldpropagator%phi_r
!!$       IF (.NOT. ASSOCIATED(fieldpropagator%next)) EXIT
!!$       fieldpropagator => fieldpropagator%next
!!$    end do
!!$    close(unit=5001)
!!$    !pause

  ! go to the first propagator which is wanted
  fieldperiod => fieldline%ch_fir 
  fieldpropagator => fieldperiod%ch_fir
  DO WHILE (fieldpropagator%tag .LT. proptag_first)
     IF (.NOT. ASSOCIATED(fieldpropagator%next)) EXIT
     fieldpropagator => fieldpropagator%next
  END DO
  ! internally used for plotting
  plotpropagator => fieldpropagator
  fieldripple => fieldpropagator%ch_act



  rippletag_old = 0
  ! now go the last propagator and do the job
  iend = 0
  iendperiod = 0
  prop_ibegperiod = 1
  prop_count_call = 0
  PRINT *, 'Fixing Magnetics in all Propagators'
  allpropsfixmag: DO      ! WHILE (fieldpropagator%tag .LE. proptag_end)
     ! information about propagator
     CALL info_magnetics(fieldpropagator)
     CALL info_magnetics(fieldpropagator%parent)
     ! PAUSE
     fieldperiod => fieldpropagator%parent
     fieldripple => fieldpropagator%ch_act
     rippletag = fieldripple%tag
     proptag = fieldpropagator%tag
     WRITE(c_propagator_tag,*) proptag
     WRITE(c_ripple_tag,*) rippletag
     WRITE(c_period_tag,*) fieldperiod%tag

     ! now modify propagator if you want
     !  files prop.dat and propm.dat contain columns
     !   R,phi,Z,bhat,geodcu,h_phi,dlogbdphi
     !  before (prop.dat) and after (propm.dat) modification
     !
     !   plot "prop.dat" u 2:4 w l, "propm.dat" u 2:4 w p
     !
     !  program stops with pause
     IF (plot_prop .EQ. 1) THEN
        CALL info_magnetics(fieldpropagator)
        CALL info_magnetics(fieldripple)
        WRITE(c_filename,'(100A)') 'prop.dat'
        CALL plot_magnetics(plotpropagator,proptag,proptag,c_filename)
        WRITE(c_filename,'(100A)') c_fieldprop,TRIM(ADJUSTL(c_propagator_tag)),c_extprop
        CALL plot_magnetics(plotpropagator,proptag,proptag,c_filename)
        PRINT *, 'plot_magnetics finished!'
     END IF

!!$     IF (mag_save_memory .EQ. 1) THEN
!!$        i_min_sav = fieldpropagator%i_min
!!$        x2_ub = UBOUND(fieldpropagator%coords%x2,1)
!!$        ALLOCATE(x2_sav(0:x2_ub))
!!$        x2_sav   = fieldpropagator%coords%x2
!!$        ALLOCATE(bhat_sav(0:x2_ub))
!!$        bhat_sav = fieldpropagator%mdata%bhat
!!$     END IF

     IF (plot_prop .EQ. 1) PRINT *, 'Before modify_propagator'

     count_solv = 0 

     !PRINT *, 'Before modify_propagator, Tag ',fieldpropagator%tag
     !PRINT *, 'phi_split_mode,phi_split_min ', phi_split_mode,phi_split_min
     !PRINT *, 'hphi_mult,count_solv ',hphi_mult,count_solv
     !PRINT *, 'UBOUND(eta_split,1) ',UBOUND(eta_split,1)
     !PRINT *, 'eta_split ',eta_split
     !PRINT *, 'UBOUND(eta_split_loc,1) ',UBOUND(eta_split_loc,1)
     !PRINT *, 'eta_split_loc ',eta_split_loc
     !PRINT *, 'fieldripple%eta ',fieldripple%eta

     IF (bin_split_mode .EQ. 1) THEN
        IF (prop_reconstruct_levels .EQ. 1) THEN        
           CALL modify_propagator(phi_split_mode,phi_place_mode,phi_split_min, &
                UBOUND(fieldripple%eta_loc,1),fieldripple%eta_loc,hphi_mult,count_solv)
           !PRINT *, '1'
        ELSE
           !PRINT *, '2'
           !PRINT *, 'ubound     ',UBOUND(eta_split,1)
           !PRINT *, eta_split
           CALL modify_propagator(phi_split_mode,phi_place_mode,phi_split_min, &
                UBOUND(fieldripple%eta,1),fieldripple%eta,hphi_mult,count_solv)
           !PAUSE
        END IF
     ELSE
        CALL modify_propagator(phi_split_mode,phi_place_mode,phi_split_min, &
             UBOUND(eta_ori,1),eta_ori,hphi_mult,count_solv)
        !PRINT *, '3'
     END IF

     IF (plot_prop .EQ. 1) PRINT *, 'After modify_propagator'

     !PRINT *, 'ubound ',UBOUND(fieldpropagator%coords%x2,1)
     !PRINT *, 'After modify_propagator'
     !PAUSE


     IF (plot_prop .EQ. 1) THEN
        !CALL info_magnetics(fieldpropagator)
        !CALL info_magnetics(fieldripple)
        plotpropagator => fieldpropagator
        WRITE(c_filename,'(100A)') 'propm.dat'
        CALL plot_magnetics(plotpropagator,proptag,proptag,c_filename)
        WRITE(c_filename,'(100A)') c_fieldprop,TRIM(ADJUSTL(c_propagator_tag)),c_extprop
        CALL plot_magnetics(plotpropagator,proptag,proptag,c_filename)
     END IF

!!$     if (fieldpropagator%tag .ge. 2 .and. fieldpropagator%tag .le. 4) then
!!$        do i = lbound(fieldpropagator%coords%x1,1),ubound(fieldpropagator%coords%x1,1)
!!$           write(7000,*) &
!!$                fieldpropagator%coords%x1(i),fieldpropagator%coords%x2(i),fieldpropagator%coords%x3(i)
!!$           write(7001,*) &
!!$                fieldpropagator%mdata%bhat(i),fieldpropagator%mdata%geodcu(i),fieldpropagator%mdata%h_phi(i), &
!!$                fieldpropagator%mdata%dlogbdphi(i)
!!$        end do
!!$     end if

     IF (fieldpropagator%tag .EQ. proptag_last) THEN 
        iend = 1
        iendperiod = 1
     ELSE
        IF (ASSOCIATED(fieldpropagator%next)) THEN
           IF (fieldpropagator%parent%tag .NE. fieldpropagator%next%parent%tag) THEN 
              iendperiod = 1
           ELSE
              iendperiod = 0
           END IF
        ELSE
           iendperiod = 1
        END IF
     END IF


!!$     CALL propagator_solver(                                  &
!!$          iend,iendperiod,bin_split_mode,eta_ori,             &
!!$          ierr_solv,ierr_join                                 &
!!$          )

     !IF (plot_prop .EQ. 1) pause

!!$     ! do some memory saving and put old information back
!!$     IF (mag_save_memory .EQ. 1) THEN
!!$        CALL set_magnetics_data(fieldpropagator,'sav')
!!$        CALL set_magnetics_data(fieldpropagator%coords%x2,x2_sav)
!!$        DEALLOCATE(x2_sav)
!!$        DEALLOCATE(bhat_sav)
!!$        CALL set_magnetics_data(fieldpropagator%mdata%bhat,bhat_sav)
!!$        fieldpropagator%i_min = i_min_sav
!!$
!!$        IF (clear_old_ripple .EQ. 1 .AND. ASSOCIATED(fieldripple%prev)) THEN
!!$           CALL deconstruct_binarysplit(fieldripple%prev%eta_bs)
!!$        END IF
!!$     END IF

     ! PRINT *, 'End of Tag'
     ! go to the next propagator or exit
     IF (fieldpropagator%tag .EQ. proptag_last) EXIT allpropsfixmag
     IF (.NOT.(ASSOCIATED(fieldpropagator%next))) THEN
        start_at_begin = 1
        fieldpropagator => fieldline%ch_fir%ch_fir
     ELSE
        IF (fieldpropagator%next%tag .LE. fieldline%ch_las%ch_las%tag) THEN
           start_at_begin = 0
           fieldpropagator => fieldpropagator%next
        ELSE
           start_at_begin = 1
           fieldpropagator => fieldline%ch_fir%ch_fir
        END IF
     END IF
  END DO allpropsfixmag
  PRINT *, 'Fixing Magnetics in all Propagators - End'
  PRINT *, '--------------------------------------------'
  !pause
  ! End of fixing magnetics in all propagators 

  ! Write taginfo.prop
  !IF ( (prop_write .EQ. 1 .OR. prop_write .EQ. 2) .AND. prop_reconstruct .EQ. 0) THEN

  !**********************************************************
  ! 15.01.2016 - taginfo is written also if prop_write .eq. 0
  !**********************************************************
  IF (prop_reconstruct .EQ. 0) THEN

     ! --- MPI SUPPORT ---

     ! Only master (or in sequential mode) writes the taginfo.prop file
     IF (mpro%isMaster()) THEN

        !*******************************
        ! HDF5
        !*******************************
        IF (prop_fileformat .EQ. 1) THEN

           call h5_create('taginfo.h5', h5id)

           call h5_add(h5id, 'prop_write', prop_write)
           call h5_add(h5id, 'tag_first', fieldline%ch_fir%ch_fir%tag)
           call h5_add(h5id, 'tag_last',  fieldline%ch_las%ch_las%tag)

           if (mpro%isParallel()) then
              call h5_add(h5id, 'parallel_storage', 1)
           else
              call h5_add(h5id, 'parallel_storage', 0)
           end if

           phi_per = twopi / device%nfp

           call h5_add(h5id, 'aiota',      surface%aiota )
           call h5_add(h5id, 'bmod0',      surface%bmod0)
           call h5_add(h5id, 'b_abs_min',  surface%b_abs_min)
           call h5_add(h5id, 'b_abs_max',  surface%b_abs_max )
           call h5_add(h5id, 'phi_per',    phi_per)
           call h5_add(h5id, 'nfp',        device%nfp)
           call h5_add(h5id, 'boozer_phi_beg',    boozer_phi_beg)
           call h5_add(h5id, 'boozer_theta_beg',  boozer_theta_beg )

           call h5_define_unlimited_array(h5id, 'tag', H5T_NATIVE_INTEGER, h5ids(1))
           call h5_define_unlimited_array(h5id, 'parent_tag', H5T_NATIVE_INTEGER, h5ids(2))
           call h5_define_unlimited_array(h5id, 'fieldperiod_phi_l', H5T_NATIVE_DOUBLE, h5ids(3))
           call h5_define_unlimited_array(h5id, 'phi_l',   H5T_NATIVE_DOUBLE, h5ids(4))
           call h5_define_unlimited_array(h5id, 'phi_r',   H5T_NATIVE_DOUBLE, h5ids(5))
           call h5_define_unlimited_array(h5id, 'theta_l', H5T_NATIVE_DOUBLE, h5ids(6))
           call h5_define_unlimited_array(h5id, 'theta_r', H5T_NATIVE_DOUBLE, h5ids(7))

           fieldperiod => fieldline%ch_fir 
           fieldpropagator => fieldperiod%ch_fir
           k = 0
           allprops_taginfo_nc: do while (fieldpropagator%tag .le. fieldline%ch_las%ch_las%tag)
              k = k + 1

              fieldperiod => fieldpropagator%parent

              phi_l = fieldpropagator%phi_l
              phi_r = fieldpropagator%phi_r
              theta_l = boozer_theta_beg + surface%aiota*(phi_l-boozer_phi_beg)
              theta_r = boozer_theta_beg + surface%aiota*(phi_r-boozer_phi_beg)
              if (fieldpropagator%tag .eq. fieldpropagator%parent%ch_fir%tag) then
                 phi_l = 0.0_dp
              else
                 phi_l = modulo(phi_l-boozer_phi_beg,phi_per)
              end if
              if (fieldpropagator%tag .eq. fieldpropagator%parent%ch_las%tag) then
                 phi_r = phi_per
              else
                 phi_r = modulo(phi_r-boozer_phi_beg,phi_per)
              end if

              call h5_append(h5ids(1), fieldpropagator%tag, k)
              call h5_append(h5ids(2), fieldpropagator%parent%tag, k)
              call h5_append(h5ids(3), fieldperiod%phi_l, k)
              call h5_append(h5ids(4), phi_l, k)
              call h5_append(h5ids(5), phi_r, k)
              call h5_append(h5ids(6), theta_l, k)
              call h5_append(h5ids(7), theta_r, k)

              IF (ASSOCIATED(fieldpropagator%next)) THEN
                 fieldpropagator => fieldpropagator%next
              else
                 exit allprops_taginfo_nc
              end IF
           end DO allprops_taginfo_nc

           call h5_close(h5id)

        else !  if (prop_fileformat .eq. 1) 
           ! ASCII

           CALL unit_propagator
           OPEN(unit=prop_unit,file=prop_ctaginfo,status='replace', &
                form=prop_format,action='write')
           WRITE(prop_unit,*) prop_write
!!$     WRITE(prop_unit,*) prop_first_tag   ! UNSOLVED PROBLEM
!!$     WRITE(prop_unit,*) prop_last_tag    ! UNSOLVED PROBLEM
           WRITE(prop_unit,*) fieldline%ch_fir%ch_fir%tag
           WRITE(prop_unit,*) fieldline%ch_las%ch_las%tag
           ! This was the end of the old write
           ! --- MPI SUPPORT ---
           ! Indicate if NEO-2 ran in parallelized mode
#if defined(MPI_SUPPORT)
           IF (mpro%isParallel()) THEN
              WRITE(prop_unit,*) .TRUE.
           ELSE
              WRITE(prop_unit,*) .FALSE.
           END IF
#else
           WRITE(prop_unit,*) .FALSE.
#endif
           ! ---

           phi_per = twopi / device%nfp

           WRITE(prop_unit,*) surface%aiota
           WRITE(prop_unit,*) surface%bmod0
           WRITE(prop_unit,*) surface%b_abs_min
           WRITE(prop_unit,*) surface%b_abs_max
           WRITE(prop_unit,*) phi_per
           WRITE(prop_unit,*) device%nfp
           WRITE(prop_unit,*) boozer_phi_beg
           WRITE(prop_unit,*) boozer_theta_beg

           fieldperiod => fieldline%ch_fir 
           fieldpropagator => fieldperiod%ch_fir
           allprops_taginfo: DO WHILE (fieldpropagator%tag .LE. fieldline%ch_las%ch_las%tag)

              fieldperiod => fieldpropagator%parent

              phi_l = fieldpropagator%phi_l
              phi_r = fieldpropagator%phi_r
              theta_l = boozer_theta_beg + surface%aiota*(phi_l-boozer_phi_beg)
              theta_r = boozer_theta_beg + surface%aiota*(phi_r-boozer_phi_beg)
              IF (fieldpropagator%tag .EQ. fieldpropagator%parent%ch_fir%tag) THEN
                 phi_l = 0.0_dp
              ELSE
                 phi_l = MODULO(phi_l-boozer_phi_beg,phi_per)
              END IF
              IF (fieldpropagator%tag .EQ. fieldpropagator%parent%ch_las%tag) THEN
                 phi_r = phi_per
              ELSE
                 phi_r = MODULO(phi_r-boozer_phi_beg,phi_per)
              END IF

              WRITE(prop_unit,*) fieldpropagator%tag, fieldpropagator%parent%tag,&
                   fieldperiod%phi_l, &
                   phi_l, phi_r, theta_l, theta_r
              IF (ASSOCIATED(fieldpropagator%next)) THEN
                 fieldpropagator => fieldpropagator%next
              ELSE
                 EXIT allprops_taginfo
              END IF
           END DO allprops_taginfo

           CLOSE(unit=prop_unit)

        END IF
     END IF
     !fieldpropagator%tag

     ! stop
     ! One can activate the stop here, if one just wants a new
     ! taginfo.prop written. Then no harm is done to all
     ! propagator output files


  END IF

  ! write hdf5 information
  call h5_magnetics(device)


  ! Now do the real computation

  ! go to the first propagator which is wanted
  fieldperiod => fieldline%ch_fir 
  fieldpropagator => fieldperiod%ch_fir
  DO WHILE (fieldpropagator%tag .LT. proptag_start)
     IF (.NOT. ASSOCIATED(fieldpropagator%next)) EXIT
     fieldpropagator => fieldpropagator%next
  END DO
  ! internally used for plotting
  plotpropagator => fieldpropagator
  fieldripple => fieldpropagator%ch_act
  !rippletag_old = 0
  ! now go the last propagator and do the job
  iend = 0
  iendperiod = 0
  prop_ibegperiod = 1
  prop_count_call = 0
!!$  IF (magnetic_device .EQ. 0) THEN
!!$     PRINT *, 'Make one propagator for tokamak'
!!$     CALL ripple_prop_joiner(fieldpropagator%parent%parent)
!!$     stop
!!$  end IF
  PRINT *, 'Do the real computation - ripple_solver'
  IF ( (magnetic_device .EQ. 0 .AND. isw_axisymm .EQ. 1) .OR. mag_magfield .EQ. 0 ) THEN
     IF (mag_magfield .EQ. 0) THEN
        PRINT *, 'Use only one propagator for homogeneous field'
     ELSE
        PRINT *, 'Use only one propagator for tokamak'
     END IF
     DO WHILE(ASSOCIATED(fieldripple%prev)) 
        fieldripple => fieldripple%prev
     END DO

     fieldripple => fieldripple%next
     PRINT *, 'fieldripple%tag            ',fieldripple%tag
     PRINT *, 'fieldripple%pa_fir%tag     ',fieldripple%pa_fir%tag
     PRINT *, 'fieldripple%pa_las%tag     ',fieldripple%pa_las%tag
     fieldpropagator => fieldripple%pa_fir


     ! fix the y-vector
     ! begin of fieldline integration
     fieldperiod => fieldpropagator%parent%parent%ch_fir
     fieldline_phi_l = fieldperiod%phi_l
     ! end of fieldline integration
     fieldperiod => fieldpropagator%parent%parent%ch_las
     fieldline_phi_r = fieldperiod%phi_r
     ! length of relevant propagator (ripple) / length of fieldline 
     y_conv_factor = (fieldpropagator%phi_r - fieldpropagator%phi_l) / (fieldline_phi_r - fieldline_phi_l)
     ALLOCATE(y_axi_averages(SIZE(fieldperiod%mdata%yend)))
     ! this is used in propagator.f90 for final output to handle the fact
     ! that only one propagator is computed and not all props are computed
     ! and joined - in some sense it is a hack
     y_axi_averages = fieldperiod%mdata%yend
     y_axi_averages(6:10)  = y_axi_averages(6:10) * y_conv_factor
     y_axi_averages(11:14) = y_axi_averages(11:14) * y_conv_factor**2

     !print *, 'y_conv_factor ',y_conv_factor

     fieldperiod => fieldpropagator%parent
     fieldripple => fieldpropagator%ch_act
     rippletag = fieldripple%tag
     proptag = fieldpropagator%tag
     WRITE(c_propagator_tag,*) proptag
     WRITE(c_ripple_tag,*) rippletag
     WRITE(c_period_tag,*) fieldperiod%tag
     count_solv = 0
     iend = 1
     iendperiod = 1
     !print *, 'flint: before propagator_solver'
     CALL propagator_solver(                                  &
          iend,iendperiod,bin_split_mode,eta_ori,             &
          ierr_solv,ierr_join                                 &
          )
     !print *, 'flint:  after propagator_solver'
  ELSE

     ! --- MPI SUPPORT
#if defined(MPI_SUPPORT)
     ! This is for the case that the program is built with MPI support, but is started in a sequential way
     if (.not. mpro%isParallel()) then
#endif
        ! Sequential program mode
        call propagator_solver(proptag_start,proptag_end,bin_split_mode,eta_ori, parallelMode = .false.)
#if defined(MPI_SUPPORT)
     else
        ! Run program in parallel mode

        ! Set some variables the scheduler needs to access
        ! This should be replaced in the future
        globalstorage%bin_split_mode = bin_split_mode
        globalstorage%eta_ori = eta_ori
        globalstorage%fieldperiod => fieldperiod
        globalstorage%fieldline => fieldline

        ! Initialize the scheduler
        sched%configFilename = 'neo2.in'
        call sched%init()

        ! Prepare runs the initial workunit, if defined
        call sched%prepare()

        ! Schedule will run initMaster on the master process to create the workunits and
        ! will set the clients into ready-mode to receive commands
        ! call omp_set_num_threads(1)
        call sched%schedule()

        ! Deallocate memory and stop the clients
        call sched%deinit()

     end if
#endif
     ! ---



!!$     allprops_comp: DO      ! WHILE (fieldpropagator%tag .LE. proptag_end)
!!$        ! information about propagator
!!$        CALL info_magnetics(fieldpropagator)
!!$        CALL info_magnetics(fieldpropagator%parent)
!!$        ! PAUSE
!!$        fieldperiod => fieldpropagator%parent
!!$        fieldripple => fieldpropagator%ch_act
!!$        rippletag = fieldripple%tag
!!$        proptag = fieldpropagator%tag
!!$        WRITE(c_propagator_tag,*) proptag
!!$        WRITE(c_ripple_tag,*) rippletag
!!$        WRITE(c_period_tag,*) fieldperiod%tag
!!$
!!$        !PRINT *, 'eta_ori ',LBOUND(eta_ori)
!!$        clear_old_ripple = 0
!!$        newripple_comp: IF (rippletag .NE. rippletag_old) THEN
!!$           clear_old_ripple = 1
!!$           rippletag_old = rippletag
!!$           fieldripple%bin_split_mode = bin_split_mode
!!$        END IF newripple_comp
!!$
!!$        IF (fieldpropagator%tag .EQ. proptag_end) THEN 
!!$           iend = 1
!!$           iendperiod = 1
!!$        ELSE
!!$           IF (ASSOCIATED(fieldpropagator%next)) THEN
!!$              IF (fieldpropagator%parent%tag .NE. fieldpropagator%next%parent%tag) THEN 
!!$                 iendperiod = 1
!!$              ELSE
!!$                 iendperiod = 0
!!$              END IF
!!$           ELSE
!!$              iendperiod = 1
!!$           END IF
!!$        END IF
!!$
!!$        count_solv = 0
!!$        CALL propagator_solver(                                  &
!!$             iend,iendperiod,bin_split_mode,eta_ori,             &
!!$             ierr_solv,ierr_join                                 &
!!$             )
!!$        ! go to the next propagator or exit
!!$        IF (fieldpropagator%tag .EQ. proptag_end) EXIT allprops_comp
!!$        IF (.NOT.(ASSOCIATED(fieldpropagator%next))) THEN
!!$           start_at_begin = 1
!!$           fieldpropagator => fieldline%ch_fir%ch_fir
!!$        ELSE
!!$           IF (fieldpropagator%next%tag .LE. fieldline%ch_las%ch_las%tag) THEN
!!$              start_at_begin = 0
!!$              fieldpropagator => fieldpropagator%next
!!$           ELSE
!!$              start_at_begin = 1
!!$              fieldpropagator => fieldline%ch_fir%ch_fir
!!$           END IF
!!$        END IF
!!$     END DO allprops_comp
  END IF
  !PRINT *, 'Final End i Flint before deallocation!'



  ! deallocation of arrays which are used locally
  IF (ALLOCATED(eta_s))      DEALLOCATE(eta_s)
  IF (ALLOCATED(eta_x0))     DEALLOCATE(eta_x0)
  IF (ALLOCATED(eta_s_loc))  DEALLOCATE(eta_s_loc)
  IF (ALLOCATED(eta_shield_loc))  DEALLOCATE(eta_shield_loc)
  IF (ALLOCATED(eta_x0_loc)) DEALLOCATE(eta_x0_loc)
  IF (ALLOCATED(eta_ori))    DEALLOCATE(eta_ori)
  IF (ALLOCATED(eta_split))  DEALLOCATE(eta_split)
  IF (ALLOCATED(eta_type))  DEALLOCATE(eta_type)

  ! deallocation of binarysplit
  !CALL deconstruct_binarysplit(eta_bs)
  !CALL deconstruct_binarysplit(eta_bs_loc)
  RETURN
END SUBROUTINE flint

SUBROUTINE modify_propagator(phi_split_mode,phi_place_mode,phi_split_min, &
     u_eta,eta_split,hphi_mult,count_solv)
  ! this is a helper routine which computes the new phi values
  !  according to phi_split_mode
  !   phi_split_mode = 1 : halfstep (as Sergei did before)
  !   phi_split_mode = 2 : places phi's according to eta-values
  !              so that renormalization in ripplesolver can
  !              work better
  !   phi_split_mode = 3 : divide phi-intervalls according to phi_divide
  ! 
  !   phi_place_mode = 1 : puts only one point between automatically placed phi's
  !   phi_place_mode = 2 : puts points according to hphi * hphi_mult
  !                        always an odd numer of points for Sergei
  !
  ! returns new "magnetics" to the fieldpropagator structure
  ! input/output
  USE rk4_kin_mod, ONLY : y
  ! input/output
  USE partpa_mod, ONLY : ipmax,ipmin
  !
  USE flint_mod, ONLY : phiarr,plot_prop
  USE magnetics_mod
  USE device_mod
  USE binarysplit_mod
  USE plagrange_mod
  use nrtype, only : dp


  IMPLICIT NONE

  !! Modifications by Andreas F. Martitsch (11.06.2014)
  ! Optional output (necessary for modeling the magnetic rotation)
  ! Note: This requires changes in "modify_propagator"
  ! (flint.f90; block using the routine commented out) and
  ! "mag_interface_mod" (mag_interface.f90).
  INTERFACE
    SUBROUTINE magdata_for_particles(phi,bhat,geodcu,h_phi,dlogbdphi,&
          bcovar_s_hat,dlogbds,dbcovar_s_hat_dphi)
      use nrtype, only : dp
      real(kind=dp), INTENT(in)            :: phi
      real(kind=dp), INTENT(out)           :: geodcu,bhat,h_phi,dlogbdphi
      real(kind=dp), OPTIONAL, INTENT(out) :: bcovar_s_hat, &
           & dlogbds, dbcovar_s_hat_dphi
    END SUBROUTINE magdata_for_particles
  END INTERFACE
  !! End Modifications by Andreas F. Martitsch (11.06.2014)

  ! parameter list
  INTEGER, INTENT(in) :: phi_split_mode
  INTEGER, INTENT(in) :: phi_place_mode,phi_split_min
  INTEGER, INTENT(in) :: u_eta
  INTEGER, INTENT(in) :: count_solv
  REAL(kind=dp), INTENT(in) :: eta_split(0:u_eta)
  REAL(kind=dp), INTENT(in)  :: hphi_mult

  ! internal
  INTEGER :: ub,i,imin
!!$  INTEGER :: s_ybeg
  INTEGER :: proptag
  INTEGER,       ALLOCATABLE :: phi_eta_ind(:,:)
!!$  REAL(kind=dp), ALLOCATABLE :: ybeg(:)
  REAL(kind=dp), ALLOCATABLE :: o_eta_split(:)
  REAL(kind=dp) :: phibeg,phiend,hphi,phi
  REAL(kind=dp) :: bhat1,b1
!!$  REAL(kind=dp) :: geodcu1,h_phi1,dlogbdphi1
!!$  REAL(kind=dp) :: b2
!!$  ! using these dnumber_struct's for things whose final size is not
!!$  !  known, see (magnetics_mod)
!!$  !   set_new
!!$  !   extract_array
!!$  !   delete_all
  TYPE(fieldpropagator_struct), POINTER :: plotpropagator
!!$  TYPE(dnumber_struct), POINTER :: x1,x2,x3
!!$  TYPE(dnumber_struct), POINTER :: bhat,geodcu,h_phi,dlogbdphi

  INTEGER :: phi_placer_status

!!$  INTEGER :: k,k_min,k_max,ubprev,nprev
!!$  INTEGER :: loc_size
!!$  REAL(kind=dp), ALLOCATABLE :: lag_fac(:)
!!$  REAL(kind=dp), ALLOCATABLE :: phi_loc(:),dphi_loc(:)
!!$  REAL(kind=dp), ALLOCATABLE :: bhat_loc(:)
!!$  REAL(kind=dp), ALLOCATABLE :: x1_loc(:)
!!$  REAL(kind=dp), ALLOCATABLE :: x3_loc(:)
!!$  REAL(kind=dp), ALLOCATABLE :: geodcu_loc(:)
!!$  REAL(kind=dp), ALLOCATABLE :: h_phi_loc(:)
!!$  REAL(kind=dp), ALLOCATABLE :: dlogbdphi_loc(:)

  INTEGER :: ubn
  INTEGER :: nlagrange = 5
  REAL(kind=dp) dummy
  
!!$  NULLIFY( x1 )
!!$  NULLIFY( x2 )
!!$  NULLIFY( x3 )
!!$  NULLIFY( bhat )
!!$  NULLIFY( geodcu )
!!$  NULLIFY( h_phi )
!!$  NULLIFY( dlogbdphi )
  
  fieldperiod => fieldpropagator%parent


  ! get starting information from propagator
!!$  s_ybeg = SIZE(fieldpropagator%mdata%ybeg,1)
  ub = UBOUND(fieldpropagator%coords%x2,1)

!!$  ALLOCATE(ybeg(s_ybeg))
!!$  ybeg = fieldpropagator%mdata%ybeg
  phibeg = fieldpropagator%coords%x2(0)
  phiend = fieldpropagator%coords%x2(ub)

  ! allocate phi_eta_ind (storage of phi values belonging to eta)
  ALLOCATE(phi_eta_ind(0:u_eta,2))
  phi_eta_ind(:,1) = 0
  phi_eta_ind(:,2) = -1

  ! constuct new phi-values
  IF (phi_split_mode .EQ. 1) THEN
     ! split steps into half
     !ub = 2*ub
     CALL linspace(phibeg,phiend,2*ub+1,phiarr)
  ELSE IF (phi_split_mode .EQ. 2) THEN
     ! split according to eta_split
     ALLOCATE (o_eta_split(1:UBOUND(eta_split,1)))
     o_eta_split = 1.0_dp/eta_split(1:)
     CALL phi_placer(phi_place_mode,phi_split_min,u_eta,o_eta_split(:), &
          phi_eta_ind,hphi_mult,phi_placer_status)
     DEALLOCATE( o_eta_split )
     IF (phi_placer_status .EQ. 1) THEN
        ! there is no eta-value available for this phi_range (intervall too small)
        ! split steps into half as in phi_split_mode .eq. 1
        CALL linspace(phibeg,phiend,2*ub+1,phiarr)
     END IF
  ELSE IF (phi_split_mode .EQ. 3) THEN
     ! new mode according to phi_divide
     CALL phi_divider(u_eta,phi_eta_ind)
  ELSE
     PRINT *, 'not implemented'
     STOP
  END IF

  ! now make the new RK-steps for pre-computed phi-values
  ubn = UBOUND(phiarr,1)
  ! ub  = UBOUND(fieldpropagator%coords%x1,1)
  ipmin =  0
  ipmax =  0
  imin  =  0
  !y     = ybeg
  !phi   = phibeg
  b1 = fieldpropagator%mdata%bhat(0)

  IF (ALLOCATED(fieldpropagator%coords%x1)) DEALLOCATE(fieldpropagator%coords%x1)
  ALLOCATE(fieldpropagator%coords%x1(0:ubn))
  IF (ALLOCATED(fieldpropagator%coords%x2)) DEALLOCATE(fieldpropagator%coords%x2)
  ALLOCATE(fieldpropagator%coords%x2(0:ubn))
  IF (ALLOCATED(fieldpropagator%coords%x3)) DEALLOCATE(fieldpropagator%coords%x3)
  ALLOCATE(fieldpropagator%coords%x3(0:ubn))
  IF (ALLOCATED(fieldpropagator%mdata%bhat)) DEALLOCATE(fieldpropagator%mdata%bhat)
  ALLOCATE(fieldpropagator%mdata%bhat(0:ubn))
  IF (ALLOCATED(fieldpropagator%mdata%geodcu)) DEALLOCATE(fieldpropagator%mdata%geodcu)
  ALLOCATE(fieldpropagator%mdata%geodcu(0:ubn))
  IF (ALLOCATED(fieldpropagator%mdata%h_phi)) DEALLOCATE(fieldpropagator%mdata%h_phi)
  ALLOCATE(fieldpropagator%mdata%h_phi(0:ubn))
  IF (ALLOCATED(fieldpropagator%mdata%dlogbdphi)) DEALLOCATE(fieldpropagator%mdata%dlogbdphi)
  ALLOCATE(fieldpropagator%mdata%dlogbdphi(0:ubn))


  allphi: DO i = 0,ubn
!!$     !CALL magdata_for_particles(phi,y,bhat1,geodcu1,h_phi1,dlogbdphi1)
!!$     phi = phiarr(i)
!!$     allorigphi: DO k = 0,ub
!!$        IF ( phi .LE. fieldpropagator%coords%x2(k) ) EXIT
!!$     END DO allorigphi
!!$  
!!$     loc_size = 6
!!$     ALLOCATE(phi_loc(loc_size))
!!$     IF (ub .GE. 5) THEN
!!$        k_min = MAX(k-3,0)
!!$        k_max = MIN(k+2,ub)
!!$        IF (k_min .EQ. 0)  k_max = k_min + 5
!!$        IF (k_max .EQ. ub) k_min = k_max - 5
!!$        phi_loc       = fieldpropagator%coords%x2(k_min:k_max)
!!$     ELSE
!!$        ubprev = UBOUND(fieldpropagator%prev%coords%x2,1)
!!$        nprev  = 6-ub-1
!!$        phi_loc(1:nprev)       = fieldpropagator%prev%coords%x2(ubprev-nprev:ubprev-1)
!!$        phi_loc(6-ub:6)        = fieldpropagator%coords%x2(0:ub)        
!!$     END IF
!!$     ALLOCATE(dphi_loc(loc_size-1))
!!$     dphi_loc = phi_loc(2:loc_size) - phi_loc(1:loc_size-1)
!!$
!!$     IF (MINVAL(dphi_loc) .LT. 1.d-7) THEN
!!$        loc_size = 2
!!$        DEALLOCATE(phi_loc)
!!$        ALLOCATE(phi_loc(loc_size))
!!$        IF (k .LT. ub) THEN
!!$           k_min = k
!!$           k_max = k + 1
!!$        ELSE
!!$           k_min = k - 1
!!$           k_max = k
!!$        END IF
!!$        phi_loc       = fieldpropagator%coords%x2(k_min:k_max)
!!$     END IF
!!$     DEALLOCATE(dphi_loc)
!!$
!!$     ALLOCATE( bhat_loc(loc_size))
!!$     ALLOCATE( x1_loc(loc_size))
!!$     ALLOCATE( x3_loc(loc_size))
!!$     ALLOCATE( geodcu_loc(loc_size))
!!$     ALLOCATE( h_phi_loc(loc_size))
!!$     ALLOCATE( dlogbdphi_loc(loc_size))
!!$     ALLOCATE( lag_fac(loc_size))
!!$     
!!$     IF (loc_size .EQ. 2 .OR. (loc_size .EQ. 6 .AND. ub .GE. 5)) THEN
!!$        bhat_loc      = fieldpropagator%mdata%bhat(k_min:k_max)
!!$        x1_loc        = fieldpropagator%coords%x1(k_min:k_max)
!!$        x3_loc        = fieldpropagator%coords%x3(k_min:k_max)
!!$        geodcu_loc    = fieldpropagator%mdata%geodcu(k_min:k_max)
!!$        h_phi_loc     = fieldpropagator%mdata%h_phi(k_min:k_max)
!!$        dlogbdphi_loc = fieldpropagator%mdata%dlogbdphi(k_min:k_max)
!!$     ELSE
!!$        bhat_loc(1:nprev)      = fieldpropagator%prev%mdata%bhat(ubprev-nprev:ubprev-1)
!!$        bhat_loc(6-ub:6)       = fieldpropagator%mdata%bhat(0:ub)
!!$        x1_loc(1:nprev)        = fieldpropagator%prev%coords%x1(ubprev-nprev:ubprev-1)
!!$        x1_loc(6-ub:6)         = fieldpropagator%coords%x1(0:ub)
!!$        x3_loc(1:nprev)        = fieldpropagator%prev%coords%x3(ubprev-nprev:ubprev-1)
!!$        x3_loc(6-ub:6)         = fieldpropagator%coords%x3(0:ub)
!!$        geodcu_loc(1:nprev)    = fieldpropagator%prev%mdata%geodcu(ubprev-nprev:ubprev-1)
!!$        geodcu_loc(6-ub:6)     = fieldpropagator%mdata%geodcu(0:ub)
!!$        h_phi_loc(1:nprev)     = fieldpropagator%prev%mdata%h_phi(ubprev-nprev:ubprev-1)
!!$        h_phi_loc(6-ub:6)      = fieldpropagator%mdata%h_phi(0:ub)
!!$        dlogbdphi_loc(1:nprev) = fieldpropagator%prev%mdata%dlogbdphi(ubprev-nprev:ubprev-1)
!!$        dlogbdphi_loc(6-ub:6)  = fieldpropagator%mdata%dlogbdphi(0:ub)
!!$     END IF
!!$
!!$     !PRINT *, 'phi_loc ',phi_loc
!!$     !PRINT *, 'bhat_loc ',bhat_loc
!!$     !PRINT *, 'phi ',phi
!!$     IF (loc_size .EQ. 6) THEN
!!$        CALL lagrange_coefs5(phi,phi_loc,lag_fac)
!!$     ELSE ! linear
!!$        lag_fac(2) = (phi - phi_loc(1)) / (phi_loc(2) - phi_loc(1))
!!$        lag_fac(1) = 1.0_dp - lag_fac(2)
!!$     END IF
!!$     bhat1 = SUM(lag_fac*bhat_loc)

     fieldpropagator%coords%x2(i) = phiarr(i)
     CALL plagrange_interp(fieldperiod,        &
          fieldpropagator%coords%x2(i),        &
          nlagrange,                           &
          fieldpropagator%coords%x1(i),        &
          fieldpropagator%coords%x3(i),        &
          fieldpropagator%mdata%bhat(i),       &
          fieldpropagator%mdata%geodcu(i),     &
          fieldpropagator%mdata%h_phi(i),      &
          fieldpropagator%mdata%dlogbdphi(i)   &
          )

     bhat1 = fieldpropagator%mdata%bhat(i)
     IF (bhat1 .LT. b1) THEN
        imin = i
        b1 = bhat1
     END IF

!!$     CALL set_new(x1,SUM(lag_fac*x1_loc))
!!$     CALL set_new(x2,phi)
!!$     CALL set_new(x3,SUM(lag_fac*x3_loc))
!!$     CALL set_new(bhat,bhat1)
!!$     CALL set_new(geodcu,SUM(lag_fac*geodcu_loc))
!!$     CALL set_new(h_phi,SUM(lag_fac*h_phi_loc))
!!$     CALL set_new(dlogbdphi,SUM(lag_fac*dlogbdphi_loc))
!!$     DEALLOCATE(lag_fac,phi_loc,bhat_loc,x1_loc,x3_loc,geodcu_loc,h_phi_loc,dlogbdphi_loc)

  END DO allphi
  
!!$  ! put it into the fieldpropagator
!!$  CALL extract_array(x1,fieldpropagator%coords%x1,0)
!!$  CALL extract_array(x2,fieldpropagator%coords%x2,0)
!!$  CALL extract_array(x3,fieldpropagator%coords%x3,0)
!!$  CALL extract_array(bhat,fieldpropagator%mdata%bhat,0)
!!$  CALL extract_array(geodcu,fieldpropagator%mdata%geodcu,0)
!!$  CALL extract_array(h_phi,fieldpropagator%mdata%h_phi,0)
!!$  CALL extract_array(dlogbdphi,fieldpropagator%mdata%dlogbdphi,0)

  fieldpropagator%i_min = imin
!!$  fieldpropagator%b_min = fieldpropagator%mdata%bhat(imin)
!!$  fieldpropagator%b_l   = fieldpropagator%mdata%bhat(0)
!!$  fieldpropagator%b_r   = fieldpropagator%mdata%bhat(ubn)

  ! put the phi-values which belong to eta into fieldpropagator
  IF (ALLOCATED(fieldpropagator%phi_eta_ind)) &
       DEALLOCATE(fieldpropagator%phi_eta_ind)
  ALLOCATE(fieldpropagator%phi_eta_ind(0:u_eta,2))
  fieldpropagator%phi_eta_ind = phi_eta_ind

!!$  !check placement of phi
!!$  IF (fieldpropagator%tag .EQ. 82) THEN
!!$     PRINT *, ' '
!!$     PRINT *, 'tag: ',fieldpropagator%tag
!!$     PRINT *, 'CHECK - START: ',0,ub
!!$     DO i = 0, u_eta
!!$        IF (phi_eta_ind(i,1) .NE. 0) THEN
!!$           PRINT *, 'CHECK 1        ',i,phi_eta_ind(i,1), &
!!$                1.0_dp - eta_split(i)*fieldpropagator%mdata%bhat(phi_eta_ind(i,1))
!!$        END IF
!!$        IF (phi_eta_ind(i,2) .NE. 0 .AND. phi_eta_ind(i,2) .NE. ub) THEN
!!$           PRINT *, 'CHECK 2        ',i,phi_eta_ind(i,2), &
!!$                1.0_dp - eta_split(i)*fieldpropagator%mdata%bhat(phi_eta_ind(i,2))
!!$        END IF
!!$     END DO
!!$
!!$     OPEN(123,file='bhat_mfl_mod.dat')
!!$     DO i=0,ubn
!!$        WRITE(123,*) fieldpropagator%coords%x2(i),fieldpropagator%mdata%bhat(i)
!!$     ENDDO
!!$     CLOSE(123)
!!$
!!$     OPEN(123,file='eta_mod.dat')
!!$     DO i=0,u_eta
!!$        WRITE(123,*) fieldpropagator%coords%x2(0),eta_split(i)
!!$        WRITE(123,*) fieldpropagator%coords%x2(ubn),eta_split(i)
!!$        WRITE(123,*) ' '
!!$     ENDDO
!!$     CLOSE(123)
!!$
!!$     PAUSE
!!$  END IF
  
!!$  IF (phi_split_mode .EQ. 3) THEN
!!$     
!!$     OPEN(unit=9999,file='propdiv.dat')
!!$     DO i = 0,ubn
!!$        WRITE (9999,*) fieldpropagator%coords%x2(i),fieldpropagator%mdata%bhat(i)
!!$     END DO
!!$     CLOSE(unit=9999)
!!$     
!!$  END IF

  IF (plot_prop .EQ. 1 .AND. count_solv .GE. 0) THEN
     !CALL info_magnetics(fieldpropagator)
     !CALL info_magnetics(fieldripple)
     plotpropagator => fieldpropagator
     proptag = plotpropagator%tag
     CALL plot_magnetics(plotpropagator,proptag,proptag,'propmm.dat')
     !PAUSE
  END IF

!!$  ! cleaning of dnumber_struct (only used locally)
!!$  CALL delete_all(x1)
!!$  CALL delete_all(x2)
!!$  CALL delete_all(x3)
!!$  CALL delete_all(bhat)
!!$  CALL delete_all(geodcu)
!!$  CALL delete_all(h_phi)
!!$  CALL delete_all(dlogbdphi)
  
  ! final cleaning of locally used arrays
!!$  IF (ALLOCATED(ybeg))        DEALLOCATE(ybeg)
  IF (ALLOCATED(phiarr))      DEALLOCATE(phiarr)
  IF (ALLOCATED(phi_eta_ind)) DEALLOCATE(phi_eta_ind)


END SUBROUTINE modify_propagator


SUBROUTINE phi_placer(phi_place_mode,phi_split_min,u_eta,eta_m1, &
     phi_eta_ind,hphi_mult,phi_placer_status)
  ! helping routine which places phi-values according to eta_values
  !  eta_m1 stands for 1/eta
  USE device_mod
  USE magnetics_mod, ONLY : extract_array,set_new,delete_all,dnumber_struct
  USE flint_mod, ONLY : phiarr
  USE plagrange_mod
  use math_constants, only : twopi
  use nrtype, only : dp
 
  IMPLICIT NONE

  ! parameter list
  INTEGER, INTENT(in) :: phi_place_mode,phi_split_min
  INTEGER, INTENT(in) :: u_eta
  REAL(kind=dp), INTENT(in)    :: eta_m1(1:u_eta)
  INTEGER,       INTENT(inout) :: phi_eta_ind(0:u_eta,2)
  REAL(kind=dp), INTENT(in)    :: hphi_mult
  
  ! pointer used for "unknown length"
  TYPE(dnumber_struct), POINTER :: phi_p => NULL()

  ! local stuff
  INTEGER :: imin,ub,ip,ips,ipe,iphi,i_eta,k,ie_dir,i,m,n
  INTEGER :: phi_count,arr_phi_count
  INTEGER :: nfp,nstep
  ! REAL(kind=dp)              :: x(3),y(3),t(3) ! quadratic
  REAL(kind=dp), ALLOCATABLE :: phi(:), barr(:)
  REAL(kind=dp), ALLOCATABLE :: arr_phi_sl(:), arr_phi_el(:), arr_b_d(:)
  INTEGER,       ALLOCATABLE :: arr_i_eta(:),arr_ip(:)
  REAL(kind=dp) :: phibeg,phiend,delta_phi,philast
  REAL(kind=dp) :: p_k,p_kp1,p_km1,p_u,p_l,p_m
  REAL(kind=dp) :: b_k,b_kp1,b_km1,b_u,b_l,b_m
  REAL(kind=dp) :: b_d
  REAL(kind=dp) :: a,b,c,d_p,d,aa,bb
  REAL(kind=dp) :: b1,b2,e_k
  REAL(kind=dp) :: bp_k,bp_kp1
  REAL(kind=dp) :: phi_el,phi_sl,hphi,dpl
  REAL(kind=dp) :: delta_p = 1.d-13 ! 13
  REAL(kind=dp) :: delta_b = 1.d-4  ! 8

  !REAL(kind=dp) :: phi_loc(6), barr_loc(6)
  !INTEGER :: ubprev,nprev
  !INTEGER :: k_min,k_max
  !REAL(kind=dp) :: lag_fac(6)

  INTEGER :: phi_placer_status
  INTEGER :: bin_counter

  INTEGER :: nlagrange = 5
  REAL(kind=dp) :: dummy

  !PRINT *, 'ubound placer ', UBOUND(eta_m1)

  ! first setup
  phi_placer_status = 0
  ub   = UBOUND(fieldpropagator%coords%x2,1)
  fieldperiod => fieldpropagator%parent

  imin = fieldpropagator%i_min
  ALLOCATE(phi(0:ub))
  ALLOCATE(barr(0:ub))
  phi    = fieldpropagator%coords%x2
  barr   = fieldpropagator%mdata%bhat
  IF (ALLOCATED(phiarr)) DEALLOCATE(phiarr)
  phibeg = fieldpropagator%coords%x2(0)
  phiend = fieldpropagator%coords%x2(ub)
  !hphi = fieldpropagator%coords%x2(1) - phibeg
  nfp   = fieldpropagator%parent%parent%parent%parent%nfp ! device
  nstep = fieldpropagator%parent%parent%parent%nstep      ! surface
  hphi  = (twopi / nfp) / nstep
  hphi  = hphi * hphi_mult

  ALLOCATE( arr_phi_sl(2*u_eta) )
  ALLOCATE( arr_phi_el(2*u_eta) )
  ALLOCATE( arr_b_d(2*u_eta) )
  ALLOCATE( arr_i_eta(2*u_eta) )
  ALLOCATE( arr_ip(2*u_eta) )
  
  ! first phi at begin of propagator
  phi_sl = phibeg
  !CALL set_new(phi_p,phibeg)
  
  ! eta_m1 stands for 1/eta
  !
  ! prop can exist of two parts
  !  towards minimum (1) and from minimum to end (2)
  phi_count = 0
  arr_phi_count = 0
  !IF (fieldpropagator%tag .EQ. 82) OPEN(123,file='phi_placer.dat')
  parts: DO ip = 1,2
     IF (ip .EQ. 1) THEN ! to the minimum
        ips = 0 ! 1
        ipe = imin-1
        b1  = barr(0)
        b2  = barr(imin)
        ie_dir = 1
        ! find first relevant eta_m1 on way to min
        i_eta = -1
        DO i = u_eta, 1, -1
           ! eta outside of interesting region
           IF(eta_m1(i) .LE. b2) phi_eta_ind(i,2) = 0
           IF(eta_m1(i) .LT. b1 .AND. eta_m1(i) .GT. b2) THEN
              i_eta = i
           END IF
           ! i_eta should now be placed at the first value
           ! which is inside the trapped region, or -1 if 
           ! nothing was found
        END DO
     ELSE ! from the minimum
        ips = imin
        ipe = ub - 1
        b1  = barr(ub)
        b2  = barr(imin)
        ie_dir = -1
        i_eta = -1
        ! find first relevant eta_m1 on way from min to end
        DO i = 1, u_eta
           IF(eta_m1(i) .LT. b1 .AND. eta_m1(i) .GT. b2) THEN
              i_eta = i
           END IF
        END DO
        ! now one can continue with eta-values up the hill in B
     END IF
     ! walk through all phi's and look for intersection with 1/eta
     k = ips
     philoop: DO ! k = ips,ipe
        IF (ipe .EQ. ips) THEN
           EXIT philoop
        END IF
        IF (i_eta .GE. 0) THEN ! otherwise no eta
           e_k   = eta_m1(i_eta)

           ! here i find something
           IF ( (ip.EQ.1 .AND. barr(k).GT.e_k .AND. barr(k+1).LE.e_k) .OR. &
                (ip.EQ.2 .AND. barr(k).LT.e_k .AND. barr(k+1).GE.e_k) ) THEN

!!$              ! prepare for 5-th order Lagrange polynomial
!!$              ! d_p    = phi(k+1) - phi(k)
!!$              IF (ub .GE. 5) THEN
!!$                 k_min = MAX(k-2,0)
!!$                 k_max = MIN(k+3,ub)
!!$                 IF (k_min .EQ. 0)  k_max = k_min + 5
!!$                 IF (k_max .EQ. ub) k_min = k_max - 5
!!$                 phi_loc  = phi(k_min:k_max)
!!$                 barr_loc = barr(k_min:k_max)
!!$              ELSE
!!$                 ubprev = UBOUND(fieldpropagator%prev%coords%x2,1)
!!$                 nprev  = 6-ub-1
!!$                 phi_loc(1:nprev)  = fieldpropagator%prev%coords%x2(ubprev-nprev:ubprev-1)
!!$                 barr_loc(1:nprev) = fieldpropagator%prev%mdata%bhat(ubprev-nprev:ubprev-1)
!!$                 phi_loc(6-ub:6)   = phi(0:ub)
!!$                 barr_loc(6-ub:6)  = barr(0:ub)
!!$              END IF
              ! do a binary search
              p_l = phi(k)
              p_u = phi(k+1)
              b_l = barr(k)
              b_u = barr(k+1)
              bin_counter = 0
              binsearch: DO WHILE ( (p_u-p_l) > delta_p )
                 bin_counter = bin_counter + 1
                 p_m = (p_l + p_u) / 2.0_dp
                 !CALL lagrange_coefs5(p_m,phi_loc,lag_fac)
                 !b_m = SUM(lag_fac*barr_loc)
                 CALL plagrange_interp(fieldperiod,p_m,nlagrange,b_m,dummy)
                 IF ( (ip.EQ.1 .AND. b_m.LE.e_k) .OR. &
                      (ip.EQ.2 .AND. b_m.GE.e_k) ) THEN 
                    p_u = p_m
                    b_u = b_m
                 ELSE
                    p_l = p_m
                    b_l = b_m
                 END IF
                 IF (bin_counter .EQ. 40) EXIT
              END DO binsearch

              IF (ip.EQ.1) THEN 
                 phi_el = p_u ! towards minimum
                 b_d    = b_u 
              ELSE
                 phi_el = p_l
                 b_d    = b_l
              END IF

              !IF (fieldpropagator%tag .EQ. 82) WRITE(123,*) e_k,1.0_dp/e_k,phi_el,b_d

              IF (b_d > e_k) THEN
                 PRINT *, 'WARNING', fieldpropagator%tag
                 PRINT *, 'ip,phi_el,b_d,e_k ',ip,phi_el,b_d,e_k
                 PRINT *, 'dummy             ',dummy
                 phi_el = phi_el + (e_k - 1.d-13 - b_d) / dummy
                 CALL plagrange_interp(fieldperiod,phi_el,nlagrange,b_d,dummy)
                 PRINT *, 'ip,phi_el,b_d,e_k ',ip,phi_el,b_d,e_k
                 !PAUSE
              END IF
              ! PRINT *, 'phi_el ',phi_el

              ! now one has to move phi by a small value to ensure
              ! $(1-\eta \Bhat)>0$
              ! derivative
              !d_p = phi_el - phi(k)
              !b_d = 3.0_dp*a*d_p**2 + 2.0_dp*b*d_p + c
              
              ! store the information
              arr_phi_count = arr_phi_count + 1
              arr_phi_sl(arr_phi_count) = phi_sl
              arr_phi_el(arr_phi_count) = phi_el
              arr_b_d(arr_phi_count) = b_d
              arr_i_eta(arr_phi_count) = i_eta
              arr_ip(arr_phi_count) = ip

              ! new starting phi
              phi_sl = phi_el
              i_eta = i_eta + ie_dir ! new index for eta
              e_k   = eta_m1(i_eta)
              ! here i find something
              IF ( (ip.EQ.1 .AND. barr(k).GT.e_k .AND. barr(k+1).LE.e_k) .OR. &
                   (ip.EQ.2 .AND. barr(k).LT.e_k .AND. barr(k+1).GE.e_k) ) THEN
                 ! we have to stay in the same phi intervall to do it
                 ! for another eta
                 k = k
              ELSE
                 ! we can move to the next intervall
                 k = k + 1
              END IF

           ELSE ! end of finding a new intersection
              k = k + 1
           END IF ! end of finding a new intersection
        ELSE ! no relevant eta exists
           EXIT philoop
        END IF 
        IF (k .GT. ipe) EXIT philoop
     END DO philoop
  END DO parts
  !IF (fieldpropagator%tag .EQ. 82) CLOSE(123)

  IF (arr_phi_count .EQ. 0) THEN
     phi_placer_status = 1
     RETURN
  END IF

  ! store the end
  arr_phi_count = arr_phi_count + 1
  arr_phi_sl(arr_phi_count) = phi_sl
  arr_phi_el(arr_phi_count) = phiend
  arr_b_d(arr_phi_count) = 0.0_dp
  arr_ip(arr_phi_count) = arr_ip(arr_phi_count-1)
  IF (ip .EQ. 1) THEN
     arr_i_eta(arr_phi_count) = arr_i_eta(arr_phi_count-1) + 1
  ELSE
     arr_i_eta(arr_phi_count) = arr_i_eta(arr_phi_count-1) - 1
  END IF
  ! PRINT *, 'Size arr_i_eta: ',SIZE(arr_i_eta)
  ! PRINT *, arr_phi_count

  DO k = 1, arr_phi_count
     phi_sl = arr_phi_sl(k)
     phi_el = arr_phi_el(k)
     b_d = arr_b_d(k)
     i_eta = arr_i_eta(k)
     ip = arr_ip(k)
     ! PRINT *, 'i_eta: ',i_eta
     ! PRINT *, 'ip:    ',ip

!!$     ! not necessary
!!$     ! handle the shift in phi
!!$     delta_phi = 0.0_dp
!!$     IF (k .LT. arr_phi_count) THEN
!!$        IF (b_d .NE. 0.0_dp) THEN
!!$           delta_phi = ABS(delta_b / b_d)
!!$        ELSE
!!$           delta_phi = 1.0_dp
!!$        END IF
!!$        delta_phi = MIN( (phi_el-phi_sl)/10.0_dp , delta_phi)
!!$        delta_phi = - SIGN(1.0_dp,b_d)*delta_phi
!!$        phi_el = phi_el + delta_phi
!!$        arr_phi_el(k) = phi_el
!!$        arr_phi_sl(k+1) = phi_el
!!$     END IF

     !PRINT *, ip,phi_sl,phi_el,phi_el-phi_sl,delta_phi

     ! make the phi between phi_sl and phi_el
     ! last one is only recorded for the last sub-intervall
     IF (phi_place_mode .EQ. 1) THEN
        m = 1 ! or only one point 
     ELSE
        m = CEILING((phi_el-phi_sl)/(hphi))
        m = MAX(phi_split_min,m - MOD(m+1,2))
     END IF
     dpl = (phi_el-phi_sl) / DBLE(m+1) ! local delta

     IF (k .EQ. arr_phi_count) m = m+1 ! for the last one

     DO n = 0, m
        philast = phi_sl + n*dpl
        CALL set_new(phi_p,philast)
     END DO

     ! store the phi's belonging to eta
     IF (k .LT. arr_phi_count) THEN
        phi_count = phi_count + m + 1
        IF (ip .EQ. 1) THEN ! towards minimum
           phi_eta_ind(i_eta,1) = phi_count
        ELSE
           phi_eta_ind(i_eta,2) = phi_count
        END IF
     END IF
  END DO

  ! extract information to phiarr which is used in the caller
  CALL extract_array(phi_p,phiarr)

  ! finish with phi_eta_ind
  ub = UBOUND(phiarr,1)
  DO n = 0,u_eta
     IF (phi_eta_ind(n,2) .LT. 0) phi_eta_ind(n,2) = ub
  END DO

!!$  IF (fieldpropagator%tag .EQ. 82) THEN
!!$     OPEN(123,file='phiarr.dat')
!!$     DO i=LBOUND(phiarr,1),UBOUND(phiarr,1)
!!$        WRITE(123,*) i,phiarr(i)
!!$     ENDDO
!!$     CLOSE(123)
!!$
!!$     OPEN(123,file='phiind.dat')
!!$     DO i=0,u_eta
!!$        IF (phi_eta_ind(i,1) .GT. 0 .AND. phi_eta_ind(i,1) .LT. u_eta) THEN
!!$           WRITE(123,*) eta_m1(i),phiarr(phi_eta_ind(i,1))
!!$        END IF
!!$        IF (phi_eta_ind(i,2) .NE. 0.AND. phi_eta_ind(i,2) .LT. u_eta ) THEN
!!$           WRITE(123,*) eta_m1(i),phiarr(phi_eta_ind(i,2))
!!$        END IF
!!$     ENDDO
!!$     CLOSE(123)
!!$  END IF

  
  ! get rid of pointers
  CALL delete_all(phi_p)
  
  ! deallocation of local arrays
  IF (ALLOCATED(phi))  DEALLOCATE(phi)
  IF (ALLOCATED(barr)) DEALLOCATE(barr)
  IF (ALLOCATED(arr_phi_sl))  DEALLOCATE(arr_phi_sl)
  IF (ALLOCATED(arr_phi_el))  DEALLOCATE(arr_phi_el)
  IF (ALLOCATED(arr_b_d))  DEALLOCATE(arr_b_d)
  IF (ALLOCATED(arr_i_eta))  DEALLOCATE(arr_i_eta)
  IF (ALLOCATED(arr_ip))  DEALLOCATE(arr_ip)
  
  RETURN
END SUBROUTINE phi_placer

SUBROUTINE phi_divider(u_eta,phi_eta_ind)                  !<-in Winny
  ! helper routine who provides new phi-values according to the 
  ! array phi_divide
  USE device_mod
  USE magnetics_mod, ONLY : extract_array,set_new,delete_all,dnumber_struct
  USE flint_mod, ONLY : phiarr,phi_divide
  use nrtype, only : dp
  
  IMPLICIT NONE

  INTEGER,       INTENT(in) :: u_eta
  INTEGER,       INTENT(inout) :: phi_eta_ind(0:u_eta,2)

  INTEGER :: ub
  INTEGER :: iphi,idiv,cphi,cdiv,ieta
  INTEGER, ALLOCATABLE :: phi_eta_ind_loc(:,:)
  REAL(kind=dp) :: hphi
  REAL(kind=dp), ALLOCATABLE :: phi(:)
  ! pointer used for "unknown length"
  TYPE(dnumber_struct), POINTER :: phi_p => NULL()

  ALLOCATE(phi_eta_ind_loc(0:u_eta,2))
  phi_eta_ind_loc = phi_eta_ind

  ub   = UBOUND(fieldpropagator%coords%x2,1)
  ALLOCATE(phi(0:ub))
  phi    = fieldpropagator%coords%x2
  IF (ALLOCATED(phiarr)) DEALLOCATE(phiarr)

  ! Set the first value
  CALL set_new(phi_p,phi(0))
  cphi = 0 ! counter
  philoop: DO iphi = 1,ub
     cdiv = phi_divide(iphi)
     hphi = phi(iphi) - phi(iphi-1)
     divloop: DO idiv = 1,cdiv
        cphi = cphi + 1
        CALL set_new(phi_p, phi(iphi-1)+DBLE(idiv)*hphi/DBLE(cdiv))
     END DO divloop
     etaloop: DO ieta = 0,u_eta ! make replacements in phi_eta_ind
        IF (phi_eta_ind_loc(ieta,1) .EQ. iphi) THEN
           phi_eta_ind(ieta,1) = cphi
        END IF
        IF (phi_eta_ind_loc(ieta,2) .EQ. iphi) THEN
           phi_eta_ind(ieta,2) = cphi
        END IF
     END DO etaloop
  END DO philoop
  ! extract information to phiarr which is used in the caller
  CALL extract_array(phi_p,phiarr)

!!$  OPEN(unit=9999,file='divide.dat')
!!$  DO iphi = 1,ub
!!$     WRITE (9999,*) phi_divide(iphi)
!!$  END DO
!!$  CLOSE(unit=9999)
!!$  !stop

  ! get rid of pointers
  CALL delete_all(phi_p)
  
  ! deallocation of local arrays
  IF (ALLOCATED(phi))  DEALLOCATE(phi)
  IF (ALLOCATED(phi_eta_ind_loc)) DEALLOCATE(phi_eta_ind_loc)

END SUBROUTINE phi_divider

SUBROUTINE lagrange_coefs5(u,up,cu)
  use nrtype, only : dp

  IMPLICIT NONE

  INTEGER,          PARAMETER  :: mp=6
  real(kind=dp)                :: u
  real(kind=dp), DIMENSION(mp) :: up,cu

  cu(1) = (u - up(2))/(up(1) - up(2))        &
       * (u - up(3))/(up(1) - up(3))         &
       * (u - up(4))/(up(1) - up(4))         &
       * (u - up(5))/(up(1) - up(5))         &
       * (u - up(6))/(up(1) - up(6))
  cu(2) = (u - up(1))/(up(2) - up(1))        &
       * (u - up(3))/(up(2) - up(3))         &
       * (u - up(4))/(up(2) - up(4))         &
       * (u - up(5))/(up(2) - up(5))         &
       * (u - up(6))/(up(2) - up(6))
  cu(3) = (u - up(1))/(up(3) - up(1))        &
       * (u - up(2))/(up(3) - up(2))         &
       * (u - up(4))/(up(3) - up(4))         &
       * (u - up(5))/(up(3) - up(5))         &
       * (u - up(6))/(up(3) - up(6))
  cu(4) = (u - up(1))/(up(4) - up(1))        &
       * (u - up(2))/(up(4) - up(2))         &
       * (u - up(3))/(up(4) - up(3))         &
       * (u - up(5))/(up(4) - up(5))         &
       * (u - up(6))/(up(4) - up(6))
  cu(5) = (u - up(1))/(up(5) - up(1))        &
       * (u - up(2))/(up(5) - up(2))         &
       * (u - up(3))/(up(5) - up(3))         &
       * (u - up(4))/(up(5) - up(4))         &
       * (u - up(6))/(up(5) - up(6))
  cu(6) = (u - up(1))/(up(6) - up(1))        &
       * (u - up(2))/(up(6) - up(2))         &
       * (u - up(3))/(up(6) - up(3))         &
       * (u - up(4))/(up(6) - up(4))         &
       * (u - up(5))/(up(6) - up(5))
  !
  RETURN
END SUBROUTINE lagrange_coefs5


SUBROUTINE sort_theta()
  
  USE device_mod
  USE magnetics_mod
  use nrtype, only : dp

  IMPLICIT NONE

  TYPE(fieldperiod_struct), POINTER :: p_t,p_min,p_max

  REAL(kind=dp) :: t,t_min,t_max,t_min_last
  INTEGER :: found
  t_min = 1.0d100
  t_max = -1.9d100

  fieldperiod => fieldline%ch_fir
  ! find t_min and t_max
  findfirst: DO
     t = fieldperiod%theta_b
     IF (t .LT. t_min) THEN
        t_min = t
        p_min => fieldperiod
     END IF
     IF (t .GT. t_max) THEN
        t_max = t
        p_max => fieldperiod
     END IF
     NULLIFY(fieldperiod%prev_theta)
     NULLIFY(fieldperiod%next_theta)
     IF (.NOT. ASSOCIATED(fieldperiod%next)) EXIT
     fieldperiod => fieldperiod%next
  END DO findfirst

  ! 
  fieldperiod => fieldline%ch_fir
  p_t => p_min

  found = 1
  all: DO 
     found = 0
     fieldperiod => fieldline%ch_fir
     t_min_last = t_min
     t_min = t_max
     find: DO
        IF (.NOT. ASSOCIATED(fieldperiod%prev_theta)) THEN
           t = fieldperiod%theta_b
           IF (t .GT. t_min_last .AND. t .LT. t_min .AND. t .LT. t_max) THEN
              found = 1
              t_min = t
              p_min => fieldperiod
           END IF
        END IF
        IF (.NOT. ASSOCIATED(fieldperiod%next)) EXIT find
        fieldperiod => fieldperiod%next
     END DO find
     IF (found .EQ. 0) EXIT     
     p_t%next_theta   => p_min
     p_min%prev_theta => p_t
     p_t => p_min
  END DO all
  p_t%next_theta   => p_max
  p_max%prev_theta => p_t
  

  NULLIFY(p_t)
  NULLIFY(p_min)
  NULLIFY(p_max)
  RETURN
END SUBROUTINE sort_theta

SUBROUTINE write_volume_data(n_r,n_z,n_phi,fname)
  ! this is for christian
  ! I do not know whether I use the correct syntax for writing
  ! file contains 
  !   some constants
  !   then vectors r,z,phi
  !   then array bmod(i_z,i_r) repeated n_phi-times

  USE device_mod
  USE binarysplit_mod, ONLY : linspace
  !! Modifications by Andreas F. Martitsch (09.03.2014)
  ! Collection of subroutines (mag.f90) converted to a module.
  ! This allows one to make use of generic interfaces 
  ! for module procedures (e.g., for a different number of
  ! arguments).
  ! Note: This requires changes in "flint_prepare" and 
  ! "write_volume_data" (both flint.f90), and "rhs_kin"
  ! and "magdata_for_particles". 
  USE mag_sub, ONLY: mag
  !! End Modifications by Andreas F. Martitsch (09.03.2014)
  use math_constants, only : pi
  use nrtype, only : dp
  
  IMPLICIT NONE
  
  ! parameter list
  INTEGER, INTENT(in) :: n_r,n_z,n_phi
  CHARACTER(len=*), INTENT(in) :: fname

  ! local
  CHARACTER(len=100) :: name
  CHARACTER(len=30)  :: c_r,c_z,c_phi,f,f_f,f_r,f_z,f_phi
  LOGICAL :: opened
  INTEGER :: uw
  INTEGER :: nfp
  INTEGER :: i_r,i_z,i_phi
  REAL(kind=dp) :: r_max,r_min,z_max,z_min,phi_min,phi_max
  REAL(kind=dp) :: delta_r = 5.0_dp
  REAL(kind=dp) :: delta_z = 5.0_dp
  
  REAL(kind=dp), ALLOCATABLE :: r(:),z(:),phi(:)
  REAL(kind=dp), ALLOCATABLE :: b(:,:)

  ! to call mag
  REAL(kind=dp)                 :: bmod,sqrtg
  REAL(kind=dp), DIMENSION(3)   :: x,bder,hcovar,hctrvr
  REAL(kind=dp), DIMENSION(3,3) :: hcoder,hctder


  ! get data
  name = device%name
  nfp = device%nfp
  r_max = surface%r_max
  r_min = surface%r_min
  z_max = surface%z_max
  z_min = surface%z_min
  phi_min = 0.0_dp
  phi_max = 2*pi/DBLE(nfp)

  PRINT *, '-----------------------------------'
  PRINT *, 'This is program write_volume_data'
  PRINT *, 'name         ',TRIM(name)
  PRINT *, 'nfp          ',nfp
  PRINT *, 'r_min, r_max ',r_min,r_max
  PRINT *, 'z_min, z_max ',z_min,z_max
  PRINT *, 'filename     ',fname
  PRINT *, '-----------------------------------'

  ! formats
  WRITE(c_r,*) n_r
  WRITE(c_z,*) n_z
  WRITE(c_phi,*) n_phi

  WRITE(f_f,*) 'e15.5E3'
  WRITE(f,*) '(1x,',TRIM(ADJUSTL(f_f)),')'
  WRITE(f_r,*) '(1x,',TRIM(ADJUSTL(c_r)),TRIM(ADJUSTL(f_f)),')'
  WRITE(f_z,*) '(1x,',TRIM(ADJUSTL(c_z)),TRIM(ADJUSTL(f_f)),')'
  WRITE(f_phi,*) '(1x,',TRIM(ADJUSTL(c_phi)),TRIM(ADJUSTL(f_f)),')'


  ! create r, z, phi vectors
  CALL linspace(r_min-delta_r,r_max+delta_r,n_r,r)
  CALL linspace(z_min-delta_z,z_max+delta_z,n_z,z)
  CALL linspace(phi_min,phi_max,n_phi,phi)
  ! create b-vector
  IF (ALLOCATED(b)) DEALLOCATE(b)
  ALLOCATE(b(0:n_z-1,0:n_r-1))

  ! find free unit
  uw = 100
  DO
     INQUIRE(unit=uw,opened=opened)
     IF(.NOT. opened) EXIT
     uw = uw + 100
  END DO
  OPEN(unit=uw,file=fname)

  ! write
  WRITE (uw,*) '#v version = 1.0'
  WRITE (uw,*) '#e name = ',TRIM(name)
  WRITE (uw,*) '#s nfp = ',nfp
  WRITE (uw,*) '#s n_r = ',n_r
  WRITE (uw,*) '#s n_z = ',n_z
  WRITE (uw,*) '#s n_phi = ',n_phi
  WRITE (uw,*) '#s r_min = ',r_min-delta_r
  WRITE (uw,*) '#s r_max = ',r_max+delta_r
  WRITE (uw,*) '#s z_min = ',z_min-delta_z
  WRITE (uw,*) '#s z_max = ',z_max+delta_z
  WRITE (uw,*) '#s phi_min = ',phi_min
  WRITE (uw,*) '#s phi_max = ',phi_max

  WRITE (uw,*) '#g r'
  WRITE (uw,f_r) r
  WRITE (uw,*) '#g z'
  WRITE (uw,f_z) z
  WRITE (uw,*) '#g phi'
  WRITE (uw,f_phi) phi
  
  WRITE (uw,*) '#a bfield'
  WRITE (uw,*) '#c (n_z, n_r, n_phi)'
  WRITE (uw,*) '#k modb'
  phi_loop: DO i_phi = 0, n_phi - 1
     r_loop : DO i_r = 0, n_r - 1
        z_loop : DO i_z = 0, n_z -1
           x(1)=r(i_r)  
           x(2)=phi(i_phi)
           x(3)=z(i_z)  
           CALL mag(x,bmod,sqrtg,bder,hcovar,hctrvr,hcoder,hctder)
           b(i_z,i_r) = bmod
        END DO z_loop
     END DO r_loop
     WRITE (uw,f_r) TRANSPOSE(b)
  END DO phi_loop
  CLOSE(unit=uw)

  DEALLOCATE(r)
  DEALLOCATE(z)
  DEALLOCATE(phi)
  DEALLOCATE(b)

  RETURN
END SUBROUTINE write_volume_data

SUBROUTINE write_surface_data(fname)
  ! this is for christian
  ! the file produced here should be readable for you
  USE device_mod
  use nrtype, only : dp
  IMPLICIT NONE

  ! parameter list
  CHARACTER(len=*), INTENT(in) :: fname
  CHARACTER(len=30)  :: c_the,c_phi,f,f_f,f_the,f_phi

  ! local
  LOGICAL :: opened
  CHARACTER(len=100) :: name
  INTEGER :: uw
  INTEGER :: n_the,n_phi,i_the
  INTEGER :: ic,last_tag,lb,ub
  REAL(kind=dp) :: bmod0
  REAL(kind=dp), ALLOCATABLE, DIMENSION(:)   :: the_arr
  REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:) :: r_arr,phi_arr,z_arr,bmod_arr

  name = device%name
  bmod0 = surface%bmod0
  ! sort periods from theta_min to theta_max
  CALL sort_theta

  ! find n_the
  ! goes to the first theta
  fieldperiod => fieldline%ch_fir
  DO 
     IF (.NOT. ASSOCIATED(fieldperiod%prev_theta)) EXIT
     fieldperiod => fieldperiod%prev_theta
  END DO
  ! and continues then to the last
  n_the = 0
  DO 
     n_the = n_the + 1
     IF (.NOT. ASSOCIATED(fieldperiod%next_theta)) EXIT
     fieldperiod => fieldperiod%next_theta
  END DO

  ! check for n_phi (is a little bit complicated)
  n_phi = 0
  ic = 0
  fieldperiod => fieldline%ch_fir
  fieldpropagator => fieldperiod%ch_fir
  last_tag = fieldperiod%ch_las%tag
  DO
     ic = ic + 1
     n_phi = n_phi + SIZE(fieldpropagator%coords%x1,1)
     IF (.NOT. ASSOCIATED(fieldpropagator%next)) EXIT
     IF (fieldpropagator%tag .EQ. last_tag) EXIT
     fieldpropagator => fieldpropagator%next
  END DO
  n_phi = n_phi - ic + 1
  
  ! allocation
  ALLOCATE(the_arr(0:n_the-1))
  ALLOCATE(r_arr(0:n_phi-1,0:n_the-1))
  ALLOCATE(phi_arr(0:n_phi-1,0:n_the-1))
  ALLOCATE(z_arr(0:n_phi-1,0:n_the-1))
  ALLOCATE(bmod_arr(0:n_phi-1,0:n_the-1))

  ! fill arrays
  ! goes to the first theta
  fieldperiod => fieldline%ch_fir
  DO 
     IF (.NOT. ASSOCIATED(fieldperiod%prev_theta)) EXIT
     fieldperiod => fieldperiod%prev_theta
  END DO
  ! and continues then to the last
  i_the = -1
  periods :DO 
     i_the = i_the + 1
     fieldpropagator => fieldperiod%ch_fir
     last_tag = fieldperiod%ch_las%tag
     lb = 0
     the_arr(i_the) = fieldperiod%theta_b
     props : DO
        ub = UBOUND(fieldpropagator%coords%x1,1) + lb
        r_arr(lb:ub,i_the)    = fieldpropagator%coords%x1
        phi_arr(lb:ub,i_the)  = fieldpropagator%coords%x2
        z_arr(lb:ub,i_the)    = fieldpropagator%coords%x3
        bmod_arr(lb:ub,i_the) = fieldpropagator%mdata%bhat * bmod0
        IF (.NOT. ASSOCIATED(fieldpropagator%next)) EXIT props
        IF (fieldpropagator%tag .EQ. last_tag) EXIT props
        fieldpropagator => fieldpropagator%next
        lb = ub
     END DO props
     IF (.NOT. ASSOCIATED(fieldperiod%next_theta)) EXIT
     fieldperiod => fieldperiod%next_theta
  END DO periods

  ! now do the writing
  ! formats
  WRITE(c_phi,*) n_phi
  WRITE(c_the,*) n_the
  WRITE(f_f,*) 'e17.7E3'
  WRITE(f,*) '(1x,',TRIM(ADJUSTL(f_f)),')'
  WRITE(f_phi,*) '(1x,',TRIM(ADJUSTL(c_phi)),TRIM(ADJUSTL(f_f)),')'
  WRITE(f_the,*) '(1x,',TRIM(ADJUSTL(c_the)),TRIM(ADJUSTL(f_f)),')'


  ! find free unit
  uw = 100
  DO
     INQUIRE(unit=uw,opened=opened)
     IF(.NOT. opened) EXIT
     uw = uw + 100
  END DO
  !
  OPEN(unit=uw,file=fname)
  WRITE (uw,*) '#v version = 1.0'
  WRITE (uw,*) '#e name = ',TRIM(name)
  WRITE (uw,*) '#s n_phi = ',n_phi
  WRITE (uw,*) '#s n_the = ',n_the
  WRITE (uw,*) '#s r0 = ',device%r0

  WRITE (uw,*) '#g theta'
  WRITE (uw,f_the) the_arr
  
  WRITE (uw,*) '#a geom'
  WRITE (uw,*) '#c (n_the, n_phi)'
  WRITE (uw,*) '#c r [cm], phi [-], z [cm]'
  WRITE (uw,*) '#k r'
  WRITE (uw,f_phi) r_arr
  WRITE (uw,*) '#k phi'
  WRITE (uw,f_phi) phi_arr
  WRITE (uw,*) '#k r'
  WRITE (uw,f_phi) z_arr

  WRITE (uw,*) '#a bfield'
  WRITE (uw,*) '#c (n_the, n_phi)'
  WRITE (uw,*) '#k modb'
  WRITE (uw,f_phi) bmod_arr
  CLOSE(unit=uw)
  
  OPEN(unit=uw,file='r_arr.dat')
  WRITE (uw,f_phi) r_arr
  CLOSE(unit=uw)
  OPEN(unit=uw,file='phi_arr.dat')
  WRITE (uw,f_phi) phi_arr
  CLOSE(unit=uw)
  OPEN(unit=uw,file='z_arr.dat')
  WRITE (uw,f_phi) z_arr
  CLOSE(unit=uw)
  OPEN(unit=uw,file='the_arr.dat')
  WRITE (uw,f_the) the_arr
  CLOSE(unit=uw)


  ! deallocation
  DEALLOCATE(the_arr)
  DEALLOCATE(r_arr)
  DEALLOCATE(phi_arr)
  DEALLOCATE(z_arr)
  DEALLOCATE(bmod_arr)


  RETURN
END SUBROUTINE write_surface_data
