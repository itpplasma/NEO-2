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

    call set_collpar()

  END SUBROUTINE flint_prepare

  subroutine set_collpar()
    use mpiprovider_module, only : mpro

    use collisionality_mod, only : collpar, conl_over_mfp, collpar_min, &
      & collpar_max, lsw_multispecies, num_spec, conl_over_mfp_spec, collpar_spec
    use device_mod, only : device

    implicit none

    integer, parameter :: dp = kind(1.0d0)
    real(kind=dp), parameter :: pi=3.14159265358979d0

    ! species index
    integer :: ispec

    if (.NOT. lsw_multispecies) then
       ! negative input for conl_over_mfp should provide collpar directly
       if (conl_over_mfp .gt. 0.0d0) then
          collpar=4.d0/(2.d0*pi*device%r0)*conl_over_mfp
       else
          collpar=-conl_over_mfp
       end if
       collpar_max = collpar
       collpar_min = collpar
    else
       if (allocated(collpar_spec)) deallocate(collpar_spec)
       allocate(collpar_spec(0:num_spec-1))
       ! negative input for conl_over_mfp_spec should provide collpar_spec directly
       if (all(conl_over_mfp_spec .gt. 0.0d0)) then
          collpar_spec=4.d0/(2.d0*pi*device%r0)*conl_over_mfp_spec
       else if (all(conl_over_mfp_spec .le. 0.0d0)) then
          collpar_spec=-conl_over_mfp_spec
       else
          print *,"set_collpar: Values of collisionality &
               &parameter conl_over_spec not consistent!"
          print *,"set_collpar: All entries of conl_over_spec &
               &must be >0 or <=0"
          stop
       end if
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
    end if
  end subroutine set_collpar

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

  SUBROUTINE phi_placer(phi_place_mode,phi_split_min,u_eta,eta_m1, &
       phi_eta_ind,hphi_mult,phi_placer_status)
    ! helping routine which places phi-values according to eta_values
    !  eta_m1 stands for 1/eta
    USE device_mod
    USE magnetics_mod, ONLY : extract_array,set_new,delete_all,dnumber_struct
    USE plagrange_mod

    IMPLICIT NONE
    INTEGER, PARAMETER :: dp = KIND(1.0d0)
    REAL(kind=dp), PARAMETER :: twopi = 6.28318530717959_dp

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

    IMPLICIT NONE
    INTEGER, PARAMETER :: dp = KIND(1.0d0)

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
  !
    IMPLICIT NONE
    !
    INTEGER,          PARAMETER     :: mp=6
    DOUBLE PRECISION                :: u
    DOUBLE PRECISION, DIMENSION(mp) :: up,cu
    !
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

END MODULE flint_mod
