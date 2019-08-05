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

  SUBROUTINE flint_prepare_2(bin_split_mode,eta_s_lim)
    USE collisionality_mod, ONLY : collpar_min
    USE  mag_interface_mod, ONLY : ripple_eta_magnetics

    IMPLICIT NONE
    INTEGER, PARAMETER :: dp = KIND(1.0d0)
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

    IMPLICIT NONE
    INTEGER, PARAMETER :: dp = KIND(1.0d0)
    REAL(kind=dp), PARAMETER :: pi=3.14159265358979d0

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
    IMPLICIT NONE
    INTEGER, PARAMETER :: dp = KIND(1.0d0)

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

  SUBROUTINE sort_theta()

    USE device_mod
    USE magnetics_mod

    IMPLICIT NONE

    TYPE(fieldperiod_struct), POINTER :: p_t,p_min,p_max
    INTEGER, PARAMETER :: dp = KIND(1.0d0)
    REAL(kind=dp) :: t,t_min,t_max,t_min_last
    INTEGER :: found

    t_min = 1.0d100
    t_max = -1.9d100

    p_min => fieldline%ch_fir ! initialization, should not be necessary?
    p_max => fieldline%ch_fir ! initialization, should not be necessary?

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

END MODULE flint_mod
