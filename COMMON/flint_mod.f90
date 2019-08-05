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
  LOGICAL :: bsfunc_lambda_loc_res
  DOUBLE PRECISION :: eta_savemem_dist1, eta_savemem_dist2, eta_savemem_sigma_mult

contains

  !> this routine is called before flint and does the setup for magnetics
  !> it creates the
  !>  device (stevvo, ....)
  !>  surface
  !>  fieldline and its children, grandchildren, grandgrandchildren
  !>    fieldperiod
  !>    fieldpropagator
  !>    fieldripple
  SUBROUTINE flint_prepare(phimi,rbeg,zbeg,nstep,nperiod,bin_split_mode,eta_s_lim)

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

    IMPLICIT NONE
    INTEGER, PARAMETER :: dp = KIND(1.0d0)
    ! parameter
    REAL(kind=dp), PARAMETER :: pi=3.14159265358979d0
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

END MODULE flint_mod
