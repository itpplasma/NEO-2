MODULE magnetics_mod

  USE binarysplit_mod
  !************************************
  ! HDF5
  !************************************
  USE hdf5_tools
  USE hdf5_tools_f2003
  
  IMPLICIT NONE
  
  ! double precision
  INTEGER, PARAMETER, PRIVATE :: dp = KIND(1.0d0)
  
  ! talk or be silent
  LOGICAL, PUBLIC :: mag_talk =  .TRUE. 
  LOGICAL, PUBLIC :: mag_infotalk =  .TRUE. 
  LOGICAL, PUBLIC :: mag_write_hdf5 =  .FALSE. 

  ! handling of binarysplit
  LOGICAL, PUBLIC :: mag_split_ripple = .TRUE.

  ! internal counters
  INTEGER, PRIVATE :: device_tag = 0
  INTEGER, PRIVATE :: surface_tag = 0
  INTEGER, PRIVATE :: fieldline_tag = 0
  INTEGER, PRIVATE :: fieldperiod_tag = 0
  INTEGER, PRIVATE :: fieldpropagator_tag = 0
  INTEGER, PRIVATE :: fieldripple_tag = 0

  ! internal names
  CHARACTER(len=100), PRIVATE :: fieldperiod_name = 'period.dat'
  CHARACTER(len=100), PRIVATE :: fieldpropagator_name = 'propagator.dat'
  CHARACTER(len=100), PRIVATE :: fieldripple_name = 'ripple.dat'

  
  CHARACTER(len=100), PRIVATE :: h5_device_name          = 'device'
  CHARACTER(len=100), PRIVATE :: h5_surface_name         = 'surface'
  CHARACTER(len=100), PRIVATE :: h5_fieldline_name       = 'fieldline'
  CHARACTER(len=100), PRIVATE :: h5_fieldperiod_name     = 'fieldperiod'
  CHARACTER(len=100), PRIVATE :: h5_fieldpropagator_name = 'fieldpropagator'
  CHARACTER(len=100), PRIVATE :: h5_fieldripple_name     = 'fieldripple'
  CHARACTER(len=100), PUBLIC  :: h5_magnetics_file_name  = 'magnetics.h5'

  ! internal constants
  REAL(kind=dp), PARAMETER, PRIVATE :: pi=3.14159265358979_dp

  ! ---------------------------------------------------------------------------
  ! device_struct
  !> \brief Defines a magnetic fusion device (maybe even a certain configuration).
  !>
  !> Contains pointers to the first, last and actual surface.
  !> Has also a name.
  TYPE, PUBLIC :: device_struct
     TYPE(surface_struct),              POINTER     :: ch_act     => NULL()
     TYPE(surface_struct),              POINTER     :: ch_fir     => NULL()
     TYPE(surface_struct),              POINTER     :: ch_las     => NULL()
     INTEGER                                        :: tag = 0
     ! ------------------------------------------------------------------------
     CHARACTER(len=100)                             :: name = 'undefined'
     REAL(kind=dp)                                  :: r0   = 0.0_dp
     REAL(kind=dp)                                  :: z0   = 0.0_dp
     INTEGER                                        :: nfp  = 0
  END TYPE device_struct
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  ! surface_struct
  !> \brief Type to describe a magnetic surface.
  !>
  !> Contains a pointer to the device it belongs to.
  !> Contains pointers to the previous and next (in radial direction?)
  !> surfaces.
  !> Contains pointers to the actual, first and last fieldline.
  !> Contains some information about the fieldline.
  TYPE, PUBLIC :: surface_struct     
     TYPE(device_struct),               POINTER     :: parent     => NULL()
     TYPE(surface_struct),              POINTER     :: prev       => NULL()
     TYPE(surface_struct),              POINTER     :: next       => NULL()
     TYPE(fieldline_struct),            POINTER     :: ch_act     => NULL()
     TYPE(fieldline_struct),            POINTER     :: ch_fir     => NULL()
     TYPE(fieldline_struct),            POINTER     :: ch_las     => NULL()
     INTEGER                                        :: tag = 0
     ! ------------------------------------------------------------------------
     REAL(kind=dp)                                  :: bmod0     = 0.0_dp
     REAL(kind=dp)                                  :: aiota     = 0.0_dp
     REAL(kind=dp)                                  :: r_min     = 0.0_dp
     REAL(kind=dp)                                  :: r_max     = 0.0_dp
     REAL(kind=dp)                                  :: z_min     = 0.0_dp
     REAL(kind=dp)                                  :: z_max     = 0.0_dp
     REAL(kind=dp)                                  :: b_abs_min = 0.0_dp
     REAL(kind=dp)                                  :: b_abs_max = 0.0_dp
     INTEGER                                        :: nperiod   = 0
     INTEGER                                        :: nstep     = 0
     INTEGER                                        :: ndim      = 0
  END TYPE surface_struct
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  ! fieldline_struct
  TYPE, PUBLIC :: fieldline_struct
     TYPE(surface_struct),              POINTER     :: parent     => NULL()
     TYPE(fieldline_struct),            POINTER     :: prev       => NULL()
     TYPE(fieldline_struct),            POINTER     :: next       => NULL()
     TYPE(fieldperiod_struct),          POINTER     :: ch_act     => NULL()
     TYPE(fieldperiod_struct),          POINTER     :: ch_fir     => NULL()
     TYPE(fieldperiod_struct),          POINTER     :: ch_las     => NULL()
     INTEGER                                        :: tag = 0
     ! ------------------------------------------------------------------------
     INTEGER                                        :: abs_max_ptag = 0
     INTEGER                                        :: abs_min_ptag = 0
     REAL(kind=dp)                                  :: b_abs_min = 0.0_dp
     REAL(kind=dp)                                  :: b_abs_max = 0.0_dp
     REAL(kind=dp)                  :: xstart(3) = (/0.0_dp,0.0_dp,0.0_dp/)
  ! ---------------------------------------------------------------------------
  END TYPE fieldline_struct

  ! ---------------------------------------------------------------------------
  ! fieldperiod_struct
  TYPE, PUBLIC :: fieldperiod_struct
     TYPE(fieldline_struct),            POINTER     :: parent     => NULL()
     TYPE(fieldperiod_struct),          POINTER     :: prev       => NULL()
     TYPE(fieldperiod_struct),          POINTER     :: next       => NULL()
     TYPE(fieldperiod_struct),          POINTER     :: prev_theta => NULL()
     TYPE(fieldperiod_struct),          POINTER     :: next_theta => NULL()
     TYPE(fieldpropagator_struct),      POINTER     :: ch_act     => NULL()
     TYPE(fieldpropagator_struct),      POINTER     :: ch_fir     => NULL()
     TYPE(fieldpropagator_struct),      POINTER     :: ch_las     => NULL()
     TYPE(fieldpropagator_struct),      POINTER     :: ch_ext     => NULL()
     TYPE(inumber_struct),              POINTER     :: tag_child  => NULL()
     INTEGER                                        :: tag = 0
     INTEGER                                        :: extra = 0
     ! ------------------------------------------------------------------------
     REAL(kind=dp)                                  :: phi_l   = 0.0_dp
     REAL(kind=dp)                                  :: phi_r   = 0.0_dp
     REAL(kind=dp)                                  :: theta_b = 0.0_dp
     ! ------------------------------------------------------------------------
     REAL(kind=dp),        ALLOCATABLE              :: phi_ext(:)
     REAL(kind=dp),        ALLOCATABLE              :: bhat_ext(:)
     REAL(kind=dp),        ALLOCATABLE              :: dbp_ext(:)
     REAL(kind=dp),        ALLOCATABLE              :: d2bp_ext(:)
     INTEGER,              ALLOCATABLE              :: minmax(:)
     REAL(kind=dp),        ALLOCATABLE              :: width_left(:)
     REAL(kind=dp),        ALLOCATABLE              :: width_right(:)
     ! ------------------------------------------------------------------------
     TYPE(coordinates_struct),          POINTER     :: coords     => NULL()
     TYPE(magneticdata_struct),         POINTER     :: mdata      => NULL()
  ! ---------------------------------------------------------------------------
  END TYPE fieldperiod_struct

  ! ---------------------------------------------------------------------------
  ! fieldpropagator_struct
  TYPE, PUBLIC :: fieldpropagator_struct
     TYPE(fieldperiod_struct),          POINTER     :: parent     => NULL()
     TYPE(fieldpropagator_struct),      POINTER     :: prev       => NULL()
     TYPE(fieldpropagator_struct),      POINTER     :: next       => NULL()
     TYPE(fieldripple_struct),          POINTER     :: ch_act     => NULL()
     INTEGER                                        :: ch_tag
     INTEGER                                        :: tag = 0
     ! ------------------------------------------------------------------------
     REAL(kind=dp)                                  :: phi_l   = 0.0_dp
     REAL(kind=dp)                                  :: phi_r   = 0.0_dp
     REAL(kind=dp)                                  :: b_l     = 0.0_dp
     REAL(kind=dp)                                  :: b_r     = 0.0_dp
     REAL(kind=dp)                                  :: phi_min = 0.0_dp
     REAL(kind=dp)                                  :: b_min   = 0.0_dp
     INTEGER                                        :: i_min   = 0
     INTEGER                                        :: has_min = 0
!!$     REAL(kind=dp)                                  :: width   = 0.0_dp
!!$     REAL(kind=dp)                                  :: dist_l  = 0.0_dp
!!$     REAL(kind=dp)                                  :: dist_r  = 0.0_dp
     INTEGER,                           ALLOCATABLE :: phi_eta_ind(:,:)
     TYPE(coordinates_struct),          POINTER     :: coords     => NULL()
     TYPE(magneticdata_struct),         POINTER     :: mdata      => NULL()
     ! ------------------------------------------------------------------------
!!$     REAL(kind=dp),        ALLOCATABLE              :: eta(:)
!!$     INTEGER                                        :: bin_split_mode = 0
!!$     TYPE(dnumber_struct), POINTER                  :: eta_x0     => NULL()
!!$     TYPE(dnumber_struct), POINTER                  :: eta_s      => NULL()
!!$     TYPE(binarysplit)                              :: eta_bs
  ! ---------------------------------------------------------------------------
  END TYPE fieldpropagator_struct

  ! ---------------------------------------------------------------------------
  ! fieldripple_struct
  TYPE, PUBLIC :: fieldripple_struct
     TYPE(fieldpropagator_struct),      POINTER     :: parent     => NULL()
     TYPE(fieldpropagator_struct),      POINTER     :: pa_fir     => NULL()
     TYPE(fieldpropagator_struct),      POINTER     :: pa_las     => NULL()
     TYPE(fieldripple_struct),          POINTER     :: prev       => NULL()
     TYPE(fieldripple_struct),          POINTER     :: next       => NULL()
     INTEGER                                        :: tag = 0
     ! ------------------------------------------------------------------------
     REAL(kind=dp)                                  :: b_max_l    = 0.0_dp
     REAL(kind=dp)                                  :: b_max_r    = 0.0_dp
     REAL(kind=dp)                                  :: b_min      = 0.0_dp
     REAL(kind=dp)                                  :: d2bp_max_l = 0.0_dp
     REAL(kind=dp)                                  :: d2bp_max_r = 0.0_dp
     REAL(kind=dp)                                  :: d2bp_min   = 0.0_dp
     REAL(kind=dp)                                  :: width      = 0.0_dp
     REAL(kind=dp)                                  :: width_l    = 0.0_dp ! left of min
     REAL(kind=dp)                                  :: width_r    = 0.0_dp
     ! ------------------------------------------------------------------------
     REAL(kind=dp),        ALLOCATABLE              :: phi_inflection(:)
     REAL(kind=dp),        ALLOCATABLE              :: b_inflection(:)
     REAL(kind=dp),        ALLOCATABLE              :: dbdp_inflection(:)
     ! ------------------------------------------------------------------------
     REAL(kind=dp),        ALLOCATABLE              :: eta(:)
     REAL(kind=dp),        ALLOCATABLE              :: eta_loc(:)
     LOGICAL                                        :: shielding_ll
     LOGICAL                                        :: shielding_lr
     LOGICAL                                        :: shielding_rr
     LOGICAL                                        :: shielding_rl
     INTEGER                                        :: bin_split_mode = 0
     TYPE(dnumber_struct), POINTER                  :: eta_x0     => NULL()
     TYPE(dnumber_struct), POINTER                  :: eta_s      => NULL()
     TYPE(dnumber_struct), POINTER                  :: eta_cl     => NULL()
     TYPE(dnumber_struct), POINTER                  :: eta_shield => NULL()
     TYPE(dnumber_struct), POINTER                  :: eta_type   => NULL()
     TYPE(binarysplit)                              :: eta_bs
     TYPE(binarysplit)                              :: eta_bs_loc
     ! ------------------------------------------------------------------------
     REAL(kind=dp)                                  :: eta_boundary_left = 0.0_dp
     REAL(kind=dp)                                  :: eta_boundary_modification_left = 0.0_dp
     INTEGER                                        :: eta_boundary_index_left = 0
     REAL(kind=dp)                                  :: eta_boundary_right = 0.0_dp
     REAL(kind=dp)                                  :: eta_boundary_modification_right = 0.0_dp
     INTEGER                                        :: eta_boundary_index_right = 0
  ! ---------------------------------------------------------------------------
  END TYPE fieldripple_struct

  ! ---------------------------------------------------------------------------
  ! inumber_struct
  !> \brief Linked list element for storing integer values.
  ! ---------------------------------------------------------------------------
  TYPE, PUBLIC :: inumber_struct
     TYPE(inumber_struct), POINTER                  :: prev       => NULL()
     TYPE(inumber_struct), POINTER                  :: next       => NULL()
     ! ------------------------------------------------------------------------
     INTEGER                                        :: i = 0
  END TYPE inumber_struct
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  ! dnumber_struct
  !> \brief Linked list element for storing real values.
  ! ---------------------------------------------------------------------------
  TYPE, PUBLIC :: dnumber_struct
     TYPE(dnumber_struct), POINTER                  :: prev       => NULL()
     TYPE(dnumber_struct), POINTER                  :: next       => NULL()
     ! ------------------------------------------------------------------------
     REAL(kind=dp)                                  :: d = 0.0_dp
  END TYPE dnumber_struct
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  ! coordinates_struct
  !> \brief 3D coordinates of multiple points in one structure.
  !>
  !> Stores the three dimensional coordinates of points in three arrays,
  !> i.e. one for each dimension. The arrays are allocatable.
  TYPE, PUBLIC :: coordinates_struct
     REAL(kind=dp), ALLOCATABLE                     :: x1(:)
     REAL(kind=dp), ALLOCATABLE                     :: x2(:)
     REAL(kind=dp), ALLOCATABLE                     :: x3(:)
  END TYPE coordinates_struct
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  ! magneticdata_struct
  ! ---------------------------------------------------------------------------
  TYPE, PUBLIC :: magneticdata_struct
     REAL(kind=dp), ALLOCATABLE                     :: bhat(:)
     REAL(kind=dp), ALLOCATABLE                     :: geodcu(:)
     REAL(kind=dp), ALLOCATABLE                     :: h_phi(:)
     REAL(kind=dp), ALLOCATABLE                     :: dlogbdphi(:)
     !! Modifications by Andreas F. Martitsch (13.03.2014)
     ! Optional output (necessary for modeling the magnetic rotation)
     REAL(kind=dp), ALLOCATABLE                     :: dbcovar_s_hat_dphi(:)
     REAL(kind=dp), ALLOCATABLE                     :: bcovar_s_hat(:)
     REAL(kind=dp), ALLOCATABLE                     :: dlogbds(:)
     !! End Modifications by Andreas F. Martitsch (13.03.2014)
     REAL(kind=dp), ALLOCATABLE                     :: ybeg(:)
     REAL(kind=dp), ALLOCATABLE                     :: yend(:)
  END TYPE magneticdata_struct
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  ! interfaces for routines
  ! ---------------------------------------------------------------------------
  
  ! ---------------------------------------------------------------------------
  !! Modifications by Andreas F. Martitsch (11.06.2014)
  ! Optional output (necessary for modeling the magnetic rotation)
  ! allocation
  PUBLIC  set_magnetics_data
  PRIVATE  set_mag_data, set_mag_data_prop, set_mag_data_prop_de, &
       set_mag_data_prop2
  INTERFACE set_magnetics_data
     MODULE PROCEDURE set_mag_data, set_mag_data_prop, set_mag_data_prop_de, &
          set_mag_data_prop2
  END INTERFACE
  !! End Modifications by Andreas F. Martitsch (11.06.2014)
  ! ---------------------------------------------------------------------------
  
  ! ---------------------------------------------------------------------------
  ! constructor
  PUBLIC construct_magnetics
  PRIVATE                             &
       construct_mag_device,          &
       construct_mag_surface,         &
       construct_mag_fieldline,       &
       construct_mag_fieldperiod,     &
       construct_mag_fieldpropagator, &
       construct_mag_fieldripple
  INTERFACE construct_magnetics
     MODULE PROCEDURE                       &
          construct_mag_device,             &
          construct_mag_surface,            &
          construct_mag_fieldline,          &
          construct_mag_fieldperiod,        &
          construct_mag_fieldpropagator,    &
          construct_mag_fieldripple
  END INTERFACE
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  ! destructor
  PUBLIC destruct_magnetics
  PRIVATE                            &
       destruct_mag_device,          &
       destruct_mag_surface,         &
       destruct_mag_fieldline,       &
       destruct_mag_fieldperiod,     &
       destruct_mag_fieldpropagator, &
       destruct_mag_fieldripple
  INTERFACE destruct_magnetics
     MODULE PROCEDURE                       &
          destruct_mag_device,              &
          destruct_mag_surface,             &
          destruct_mag_fieldline,           &
          destruct_mag_fieldperiod,         &
          destruct_mag_fieldpropagator,     &
          destruct_mag_fieldripple
  END INTERFACE
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  ! info
  PUBLIC h5_magnetics
  PRIVATE                             &
       h5_mag_general_d1a,            &
       h5_mag_device,                 &
       h5_mag_surface,                &
       h5_mag_fieldline,              &
       h5_mag_fieldperiod,            &
       h5_mag_fieldpropagator,        &
       h5_mag_fieldripple
  INTERFACE h5_magnetics
     MODULE PROCEDURE                      &
          h5_mag_general_d1a,              &
          h5_mag_device,                   &
          h5_mag_surface,                  &
          h5_mag_fieldline,                &
          h5_mag_fieldperiod,              &
          h5_mag_fieldpropagator,          &
          h5_mag_fieldripple
  END INTERFACE 
  ! ---------------------------------------------------------------------------



  ! ---------------------------------------------------------------------------
  ! info
  PUBLIC info_magnetics
  PRIVATE                               &
       info_mag_device,                 &
       info_mag_surface,                &
       info_mag_fieldline,              &
       info_mag_fieldpropagator,        &
       info_mag_fieldperiod,            &
       info_mag_fieldripple
  INTERFACE info_magnetics
     MODULE PROCEDURE                        &
          info_mag_device,                   &
          info_mag_surface,                  &
          info_mag_fieldline,                &
          info_mag_fieldpropagator,          &
          info_mag_fieldperiod,              &
          info_mag_fieldripple
  END INTERFACE 
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  ! plot
  PUBLIC plot_magnetics
  PRIVATE                               &
       plot_mag_fieldpropagator,        &
       plot_mag_fieldpropagator_tag,    &
       plot_mag_fieldperiod,            &
       plot_mag_fieldripple
  INTERFACE plot_magnetics
     MODULE PROCEDURE                        &
          plot_mag_fieldpropagator,          &
          plot_mag_fieldpropagator_tag,      &
          plot_mag_fieldperiod,              &
          plot_mag_fieldripple
  END INTERFACE 
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  ! things related to simple pointers (double)
  PUBLIC create_new
  PRIVATE create_new_d
  INTERFACE create_new
     MODULE PROCEDURE create_new_d
  END INTERFACE

  PUBLIC create_before
  PRIVATE create_before_d
  INTERFACE create_before
     MODULE PROCEDURE create_before_d
  END INTERFACE

  PUBLIC set_new
  PRIVATE set_new_d
  INTERFACE set_new
     MODULE PROCEDURE set_new_d
  END INTERFACE

  PUBLIC set_before
  PRIVATE set_before_d
  INTERFACE set_before
     MODULE PROCEDURE set_before_d
  END INTERFACE

  PUBLIC delete_one
  PRIVATE delete_one_d
  INTERFACE delete_one
     MODULE PROCEDURE delete_one_d
  END INTERFACE

  PUBLIC delete_all
  PRIVATE delete_all_d
  INTERFACE delete_all
     MODULE PROCEDURE delete_all_d
  END INTERFACE

  PUBLIC goto_first
  PRIVATE goto_first_d
  INTERFACE goto_first
     MODULE PROCEDURE goto_first_d
  END INTERFACE

  PUBLIC goto_last
  PRIVATE goto_last_d
  INTERFACE goto_last
     MODULE PROCEDURE goto_last_d
  END INTERFACE

  PUBLIC extract_array
  PRIVATE extract_array_d1
  INTERFACE extract_array
     MODULE PROCEDURE extract_array_d1
  END INTERFACE

  PUBLIC sort
  PRIVATE sort_d1,sort_d2
  INTERFACE sort
     MODULE PROCEDURE sort_d1,sort_d2
  END INTERFACE

  PUBLIC disp
  PRIVATE disp_d1,disp_d2
  INTERFACE disp
     MODULE PROCEDURE disp_d1,disp_d2
  END INTERFACE
  ! ---------------------------------------------------------------------------


CONTAINS
  
  ! ---------------------------------------------------------------------------
  ! setting of data
  ! 

  ! ---------------------------------------------------------------------------
  ! deallocation, allocation, data
  SUBROUTINE set_mag_data(store,value)
    REAL(kind=dp), ALLOCATABLE, INTENT(inout) :: store(:)
    REAL(kind=dp), OPTIONAL,    INTENT(in)    :: value(0:)
    
    IF (ALLOCATED(store)) DEALLOCATE(store)
    IF (PRESENT(value)) THEN
       ALLOCATE(store(0:UBOUND(value,1)))
       store = value
    END IF
  END SUBROUTINE set_mag_data
  ! ---------------------------------------------------------------------------
  SUBROUTINE set_mag_data_prop(fieldpropagator, &
       x1,x2,x3,                                &
       bhat,geodcu,h_phi,dlogbdphi,             &
       ybeg,yend)
    TYPE(fieldpropagator_struct), POINTER :: fieldpropagator
    REAL(kind=dp), INTENT(in) :: x1(0:)
    REAL(kind=dp), INTENT(in) :: x2(0:)
    REAL(kind=dp), INTENT(in) :: x3(0:)
    REAL(kind=dp), INTENT(in) :: bhat(0:)
    REAL(kind=dp), INTENT(in) :: geodcu(0:)
    REAL(kind=dp), INTENT(in) :: h_phi(0:)
    REAL(kind=dp), INTENT(in) :: dlogbdphi(0:)
    REAL(kind=dp), INTENT(in) :: ybeg(0:)
    REAL(kind=dp), INTENT(in) :: yend(0:)
    CALL set_magnetics_data(fieldpropagator%coords%x1,x1)
    CALL set_magnetics_data(fieldpropagator%coords%x2,x2)
    CALL set_magnetics_data(fieldpropagator%coords%x3,x3)
    CALL set_magnetics_data(fieldpropagator%mdata%bhat,bhat)
    CALL set_magnetics_data(fieldpropagator%mdata%geodcu,geodcu)
    CALL set_magnetics_data(fieldpropagator%mdata%h_phi,h_phi)
    CALL set_magnetics_data(fieldpropagator%mdata%dlogbdphi,dlogbdphi)
    CALL set_magnetics_data(fieldpropagator%mdata%ybeg,ybeg)
    CALL set_magnetics_data(fieldpropagator%mdata%yend,yend)
  END SUBROUTINE set_mag_data_prop
  ! ---------------------------------------------------------------------------
  !! Modifications by Andreas F. Martitsch (13.03.2014)
  ! Optional output (necessary for modeling the magnetic rotation)
  SUBROUTINE set_mag_data_prop2(fieldpropagator, &
       x1,x2,x3,                                 &
       bhat,geodcu,h_phi,dlogbdphi,              &
       ybeg,yend,                                &
       dbcovar_s_hat_dphi,bcovar_s_hat,dlogbds   &
       )
    TYPE(fieldpropagator_struct), POINTER :: fieldpropagator
    REAL(kind=dp), INTENT(in) :: x1(0:)
    REAL(kind=dp), INTENT(in) :: x2(0:)
    REAL(kind=dp), INTENT(in) :: x3(0:)
    REAL(kind=dp), INTENT(in) :: bhat(0:)
    REAL(kind=dp), INTENT(in) :: geodcu(0:)
    REAL(kind=dp), INTENT(in) :: h_phi(0:)
    REAL(kind=dp), INTENT(in) :: dlogbdphi(0:)
    REAL(kind=dp), INTENT(in) :: ybeg(0:)
    REAL(kind=dp), INTENT(in) :: yend(0:)
    REAL(kind=dp), INTENT(in) :: dbcovar_s_hat_dphi(0:)
    REAL(kind=dp), INTENT(in) :: bcovar_s_hat(0:)
    REAL(kind=dp), INTENT(in) :: dlogbds(0:)
    CALL set_magnetics_data(fieldpropagator,x1,x2,x3,bhat,&
         geodcu,h_phi,dlogbdphi,ybeg,yend)
    !
    ! These are the additional entries:
    CALL set_magnetics_data(fieldpropagator%mdata%dbcovar_s_hat_dphi,&
         dbcovar_s_hat_dphi)
    CALL set_magnetics_data(fieldpropagator%mdata%bcovar_s_hat,bcovar_s_hat)
    CALL set_magnetics_data(fieldpropagator%mdata%dlogbds,dlogbds)
    !
  END SUBROUTINE set_mag_data_prop2
  !! End Modifications by Andreas F. Martitsch (13.03.2014)
  ! ---------------------------------------------------------------------------
  SUBROUTINE set_mag_data_prop_de(fieldpropagator,opt_in)
    TYPE(fieldpropagator_struct), POINTER :: fieldpropagator
    CHARACTER(len=3), INTENT(in), OPTIONAL :: opt_in

    CHARACTER(len=3) :: opt
    
    opt = 'all'
    IF (PRESENT(opt_in)) opt = opt_in

    IF (opt .EQ. 'all') THEN
       CALL set_magnetics_data(fieldpropagator%coords%x1)
       CALL set_magnetics_data(fieldpropagator%coords%x2)
       CALL set_magnetics_data(fieldpropagator%coords%x3)
       CALL set_magnetics_data(fieldpropagator%mdata%bhat)
       CALL set_magnetics_data(fieldpropagator%mdata%geodcu)
       CALL set_magnetics_data(fieldpropagator%mdata%h_phi)
       CALL set_magnetics_data(fieldpropagator%mdata%dlogbdphi)
       !! Modifications by Andreas F. Martitsch (11.06.2014)
       ! Optional output (necessary for modeling the magnetic rotation) ->
       ! Deallocate the additional entries (done, only if arrays allocated)
       CALL set_magnetics_data(fieldpropagator%mdata%dbcovar_s_hat_dphi)
       CALL set_magnetics_data(fieldpropagator%mdata%bcovar_s_hat)
       CALL set_magnetics_data(fieldpropagator%mdata%dlogbds)
       !! End Modifications by Andreas F. Martitsch (11.06.2014)
       CALL set_magnetics_data(fieldpropagator%mdata%ybeg)
       CALL set_magnetics_data(fieldpropagator%mdata%yend)
    ELSE IF (opt .EQ. 'sav') THEN
       ! x2, bhat, ybeg, yend stay
       CALL set_magnetics_data(fieldpropagator%coords%x1)
       CALL set_magnetics_data(fieldpropagator%coords%x3)
       CALL set_magnetics_data(fieldpropagator%mdata%geodcu)
       CALL set_magnetics_data(fieldpropagator%mdata%h_phi)
       CALL set_magnetics_data(fieldpropagator%mdata%dlogbdphi)
       !! Modifications by Andreas F. Martitsch (11.06.2014)
       ! Optional output (necessary for modeling the magnetic rotation) ->
       ! Deallocate the additional entries (done, only if arrays allocated)
       CALL set_magnetics_data(fieldpropagator%mdata%dbcovar_s_hat_dphi)
       CALL set_magnetics_data(fieldpropagator%mdata%bcovar_s_hat)
       CALL set_magnetics_data(fieldpropagator%mdata%dlogbds)       
       !! End Modifications by Andreas F. Martitsch (11.06.2014)
    ELSE
       PRINT *, 'ERROR: deallocate option in set_mag_data not known: ',opt
    END IF
  END SUBROUTINE set_mag_data_prop_de
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  !> \brief Constructor for device_struct.
  ! 
  SUBROUTINE construct_mag_device(device)
    TYPE(device_struct),     POINTER          :: device
    ! memory for device_struct
    ALLOCATE(device)
    ! tag
    device_tag = device_tag + 1
    device%tag = device_tag

    IF (mag_talk) PRINT *, 'magnetics: device created: ', device%tag
    RETURN
  END SUBROUTINE construct_mag_device
  ! ---------------------------------------------------------------------------
  !> \brief  Destructor for device_struct.
  SUBROUTINE destruct_mag_device(device)
    TYPE(device_struct),  POINTER :: device
    TYPE(surface_struct), POINTER :: surface
    INTEGER                       :: my_tag = 0
    my_tag = device%tag
    DO WHILE (ASSOCIATED(device%ch_act))
       surface => device%ch_act
       CALL destruct_magnetics(surface)
    END DO
    ! nullify children
    NULLIFY(device%ch_act)
    NULLIFY(device%ch_fir)
    NULLIFY(device%ch_las)
    ! additional deallocation
    ! nothing to be done for device
    ! final deallocation   
    DEALLOCATE(device)
    NULLIFY(device)
    IF (mag_talk) PRINT *, 'magnetics: device deleted: ',my_tag
  END SUBROUTINE destruct_mag_device
  ! end constructor and destructor for device_struct
  ! ---------------------------------------------------------------------------


  ! ---------------------------------------------------------------------------
  !> \brief constructor for surface_struct
  !
  SUBROUTINE construct_mag_surface(device,surface)
    TYPE(device_struct),     POINTER                 :: device
    TYPE(surface_struct),    POINTER                 :: surface
    INTEGER                                          :: my_tag = 0
    
    ! memory for surf
    ALLOCATE(surface)
    ! tag
    surface_tag = surface_tag + 1
    surface%tag = surface_tag
    my_tag      = surface_tag
    ! connect
    IF (ASSOCIATED(device%ch_act)) THEN
       device%ch_act%next => surface
       surface%prev       => device%ch_act
    ELSE
       device%ch_fir      => surface
    END IF
    device%ch_act  => surface
    device%ch_las  => surface
    surface%parent => device
    !
    IF (mag_talk) PRINT *, 'magnetics: surface added: ',my_tag, &
         ' parent: ',surface%parent%tag
    RETURN
  END SUBROUTINE construct_mag_surface
  ! ---------------------------------------------------------------------------
  !> \brief Destructor for surface_struct
  SUBROUTINE destruct_mag_surface(surface)
    TYPE(surface_struct),   POINTER :: surface
    TYPE(fieldline_struct), POINTER :: fieldline
    INTEGER                         :: my_tag = 0
    my_tag = surface%tag
    DO WHILE (ASSOCIATED(surface%ch_act))
       fieldline => surface%ch_act
       CALL destruct_magnetics(fieldline)
    END DO
    NULLIFY(surface%ch_act)
    NULLIFY(surface%ch_fir)
    NULLIFY(surface%ch_las)

    IF ( ASSOCIATED(surface%prev) .AND. &
         ASSOCIATED(surface%next) ) THEN
       surface%parent%ch_act => surface%next
       surface%prev%next     => surface%next
       surface%next%prev     => surface%prev
    ELSEIF ( ASSOCIATED(surface%prev) .AND. &
         .NOT. ASSOCIATED(surface%next) ) THEN
       surface%parent%ch_act => surface%prev
       surface%parent%ch_las => surface%prev
       NULLIFY(surface%prev%next)
    ELSEIF ( .NOT. ASSOCIATED(surface%prev) .AND. &
         ASSOCIATED(surface%next) ) THEN
       surface%parent%ch_act => surface%next
       surface%parent%ch_fir => surface%next
       NULLIFY(surface%next%prev)
    ELSEIF ( .NOT. ASSOCIATED(surface%prev) .AND. &
         .NOT. ASSOCIATED(surface%next) ) THEN
       NULLIFY(surface%parent%ch_act)
       NULLIFY(surface%parent%ch_fir)
       NULLIFY(surface%parent%ch_las)
    END IF
    ! additional deallocation
    ! nothing to be done for surface
    ! final deallocation
    NULLIFY(surface%prev)
    NULLIFY(surface%next)
    NULLIFY(surface%parent)
    DEALLOCATE(surface)
    NULLIFY(surface)
    ! message
    IF (mag_talk) PRINT *, 'magnetics: surface deleted: ',my_tag
  END SUBROUTINE destruct_mag_surface
  ! end constructor and destructor for surface_struct
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  !> \brief Constructor for fieldline_struct.
  ! 
  SUBROUTINE construct_mag_fieldline(surface,fieldline)
    TYPE(surface_struct),     POINTER                 :: surface
    TYPE(fieldline_struct),   POINTER                 :: fieldline
    INTEGER                                           :: my_tag = 0
    
    ! memory for fieldline
    ALLOCATE(fieldline)
    ! tag
    fieldline_tag = fieldline_tag + 1
    fieldline%tag = fieldline_tag
    my_tag        = fieldline_tag
    ! connect
    IF (ASSOCIATED(surface%ch_act)) THEN
       surface%ch_act%next => fieldline
       fieldline%prev      => surface%ch_act
    ELSE
       surface%ch_fir      => fieldline
    END IF
    surface%ch_act   => fieldline
    surface%ch_las   => fieldline
    fieldline%parent => surface
    ! additional quantities
    !
    ! end additional quantities    
    IF (mag_talk) PRINT *, 'magnetics: fieldline added: ',my_tag, &
         ' parent: ',fieldline%parent%tag
    RETURN
  END SUBROUTINE construct_mag_fieldline
  !> \brief Destructor for fieldline_struct.
  SUBROUTINE destruct_mag_fieldline(fieldline)
    TYPE(fieldline_struct),   POINTER :: fieldline
    TYPE(fieldperiod_struct), POINTER :: fieldperiod
    INTEGER                           :: my_tag = 0
    my_tag = fieldline%tag
    DO WHILE (ASSOCIATED(fieldline%ch_act))
       fieldperiod => fieldline%ch_act
       CALL destruct_magnetics(fieldperiod)
    END DO
    NULLIFY(fieldline%ch_act)
    NULLIFY(fieldline%ch_fir)
    NULLIFY(fieldline%ch_las)
    IF ( ASSOCIATED(fieldline%prev) .AND. &
         ASSOCIATED(fieldline%next) ) THEN
       fieldline%parent%ch_act => fieldline%next
       fieldline%prev%next     => fieldline%next
       fieldline%next%prev     => fieldline%prev
    ELSEIF ( ASSOCIATED(fieldline%prev) .AND. &
         .NOT. ASSOCIATED(fieldline%next) ) THEN
       fieldline%parent%ch_act => fieldline%prev
       fieldline%parent%ch_las => fieldline%prev
       NULLIFY(fieldline%prev%next)
    ELSEIF ( .NOT. ASSOCIATED(fieldline%prev) .AND. &
         ASSOCIATED(fieldline%next) ) THEN
       fieldline%parent%ch_act => fieldline%next
       fieldline%parent%ch_fir => fieldline%next
       NULLIFY(fieldline%next%prev)
    ELSEIF ( .NOT. ASSOCIATED(fieldline%prev) .AND. &
         .NOT. ASSOCIATED(fieldline%next) ) THEN
       NULLIFY(fieldline%parent%ch_act)
       NULLIFY(fieldline%parent%ch_fir)
       NULLIFY(fieldline%parent%ch_las)
    END IF
    ! additional deallocation
    !
    ! final deallocation
    NULLIFY(fieldline%prev)
    NULLIFY(fieldline%next)
    NULLIFY(fieldline%parent)
    DEALLOCATE(fieldline)
    NULLIFY(fieldline)
    ! message
    IF (mag_talk) PRINT *, 'magnetics: fieldline deleted: ',my_tag
  END SUBROUTINE destruct_mag_fieldline
  ! end constructor and destructor for fieldline_struct
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  !> \brief Constructor for fieldperiod_struct.
  ! 
  SUBROUTINE construct_mag_fieldperiod(fieldline,fieldperiod,opt_direction)
    TYPE(fieldline_struct),     POINTER               :: fieldline
    TYPE(fieldperiod_struct),   POINTER               :: fieldperiod
    INTEGER                                           :: my_tag = 0
    INTEGER,                    OPTIONAL              :: opt_direction
    INTEGER                                           :: direction
    ! optional variables
    IF ( PRESENT(opt_direction) ) THEN
       direction = opt_direction
    ELSE
       direction = 1
    END IF
    ! memory for fieldperiod
    ALLOCATE(fieldperiod)
    ! tag
    fieldperiod_tag = fieldperiod_tag + 1
    fieldperiod%tag = fieldperiod_tag
    my_tag          = fieldperiod_tag
    ! connect
    IF (direction .EQ. 1) THEN  ! normal behaviour forward
       IF (ASSOCIATED(fieldline%ch_act)) THEN
          fieldline%ch_act%next => fieldperiod
          fieldperiod%prev      => fieldline%ch_act
       ELSE
          fieldline%ch_fir      => fieldperiod
       END IF
       fieldline%ch_act   => fieldperiod
       fieldline%ch_las   => fieldperiod
       fieldperiod%parent => fieldline
    ELSEIF (direction .EQ. -1) THEN  ! backward
       IF (ASSOCIATED(fieldline%ch_act)) THEN
          fieldline%ch_act%prev => fieldperiod
          fieldperiod%next      => fieldline%ch_act
          fieldline%ch_fir      => fieldperiod
          fieldline%ch_act      => fieldperiod
          fieldperiod%parent    => fieldline
       ELSE
          PRINT *, 'Error: when going backward in construction of fieldperiods'
          PRINT *, '       the first fieldperiod must exist'
          STOP
       END IF
    END IF
    ! additional quantities
    !
    ! end additional quantities
    IF (mag_talk) PRINT *, 'magnetics: fieldperiod added, tag: ',my_tag, &
         ' parent: ',fieldperiod%parent%tag
    RETURN
  END SUBROUTINE construct_mag_fieldperiod
  !> \brief Destructor for fieldperiod_struct.
  SUBROUTINE destruct_mag_fieldperiod(fieldperiod)
    TYPE(fieldperiod_struct),     POINTER :: fieldperiod
    TYPE(fieldpropagator_struct), POINTER :: fieldpropagator
    TYPE(inumber_struct),         POINTER :: tag_child
    INTEGER                               :: my_tag = 0
    my_tag = fieldperiod%tag
    ! delete extra propagator (if it exists)
    IF (ASSOCIATED(fieldperiod%ch_ext)) THEN
       fieldpropagator => fieldperiod%ch_ext
       CALL destruct_magnetics(fieldpropagator,'ext')
       NULLIFY(fieldperiod%ch_ext)
    END IF
    ! delete childs
    IF (ASSOCIATED(fieldperiod%ch_act)) THEN
       ! find ch_las propagator which belongs to this period
       DO WHILE (ASSOCIATED(fieldperiod%tag_child%next))
          fieldperiod%tag_child => fieldperiod%tag_child%next
       END DO ! fieldperiod%tag_child points now to the ch_las one
       ! now go back
       DO
          ! search for the right tag
          DO WHILE (fieldperiod%ch_act%tag .LT. fieldperiod%tag_child%i)
             fieldperiod%ch_act => fieldperiod%ch_act%next
          END DO
          fieldpropagator => fieldperiod%ch_act
          CALL destruct_magnetics(fieldpropagator)   
          IF (.NOT. ASSOCIATED(fieldperiod%tag_child%prev)) EXIT
          fieldperiod%tag_child => fieldperiod%tag_child%prev
       END DO
       ! now clean tag_child PROBLEM
       DO WHILE (ASSOCIATED(fieldperiod%tag_child%next))
          fieldperiod%tag_child => fieldperiod%tag_child%next
       END DO ! fieldperiod%tag_child points now to the last one
       DO ! until all tag_child entries are removed 
          tag_child => fieldperiod%tag_child
          IF (ASSOCIATED(fieldperiod%tag_child%prev)) THEN
             NULLIFY(fieldperiod%tag_child%prev%next)
             fieldperiod%tag_child => fieldperiod%tag_child%prev
             DEALLOCATE(tag_child)
             NULLIFY(tag_child)
          ELSE
             DEALLOCATE(tag_child)
             NULLIFY(tag_child)
             EXIT
          END IF
       END DO
    END IF
    ! nullify children
    NULLIFY(fieldperiod%ch_act)
    NULLIFY(fieldperiod%ch_fir)
    NULLIFY(fieldperiod%ch_las)
    ! delete connections
    IF ( ASSOCIATED(fieldperiod%prev) .AND. &
         ASSOCIATED(fieldperiod%next) ) THEN
       fieldperiod%parent%ch_act => fieldperiod%next
       fieldperiod%prev%next => fieldperiod%next
       fieldperiod%next%prev => fieldperiod%prev
    ELSEIF ( ASSOCIATED(fieldperiod%prev) .AND. &
         .NOT. ASSOCIATED(fieldperiod%next) ) THEN
       fieldperiod%parent%ch_act => fieldperiod%prev
       fieldperiod%parent%ch_las => fieldperiod%prev
       NULLIFY(fieldperiod%prev%next)
    ELSEIF ( .NOT. ASSOCIATED(fieldperiod%prev) .AND. &
         ASSOCIATED(fieldperiod%next) ) THEN
       fieldperiod%parent%ch_act => fieldperiod%next
       fieldperiod%parent%ch_fir => fieldperiod%next
       NULLIFY(fieldperiod%next%prev)
    ELSEIF ( .NOT. ASSOCIATED(fieldperiod%prev) .AND. &
         .NOT. ASSOCIATED(fieldperiod%next) ) THEN
       NULLIFY(fieldperiod%parent%ch_act)
       NULLIFY(fieldperiod%parent%ch_fir)
       NULLIFY(fieldperiod%parent%ch_las)
    END IF
    ! additional deallocation
    IF (ALLOCATED(fieldperiod%phi_ext))  DEALLOCATE(fieldperiod%phi_ext)
    IF (ALLOCATED(fieldperiod%bhat_ext)) DEALLOCATE(fieldperiod%bhat_ext)
    IF (ALLOCATED(fieldperiod%dbp_ext))  DEALLOCATE(fieldperiod%dbp_ext)
    IF (ALLOCATED(fieldperiod%d2bp_ext)) DEALLOCATE(fieldperiod%d2bp_ext)
    IF (ALLOCATED(fieldperiod%minmax))   DEALLOCATE(fieldperiod%minmax)
    IF (ALLOCATED(fieldperiod%width_left)) DEALLOCATE(fieldperiod%width_left)
    IF (ALLOCATED(fieldperiod%width_right)) DEALLOCATE(fieldperiod%width_right)

    !CALL set_magnetics_data(fieldpropagator)
    IF (ALLOCATED(fieldperiod%coords%x1)) DEALLOCATE(fieldperiod%coords%x1)
    IF (ALLOCATED(fieldperiod%coords%x2)) DEALLOCATE(fieldperiod%coords%x2)
    IF (ALLOCATED(fieldperiod%coords%x3)) DEALLOCATE(fieldperiod%coords%x3)
    IF (ALLOCATED(fieldperiod%mdata%bhat)) DEALLOCATE(fieldperiod%mdata%bhat)
    IF (ALLOCATED(fieldperiod%mdata%geodcu)) DEALLOCATE(fieldperiod%mdata%geodcu)
    IF (ALLOCATED(fieldperiod%mdata%h_phi)) DEALLOCATE(fieldperiod%mdata%h_phi)
    IF (ALLOCATED(fieldperiod%mdata%dlogbdphi)) DEALLOCATE(fieldperiod%mdata%dlogbdphi)
    !! Modifications by Andreas F. Martitsch (11.06.2014)
    ! Optional output (necessary for modeling the magnetic rotation) ->
    ! Deallocate the additional entries (done, only if arrays allocated)
    IF (ALLOCATED(fieldperiod%mdata%dbcovar_s_hat_dphi)) &
         DEALLOCATE(fieldperiod%mdata%dbcovar_s_hat_dphi)
    IF (ALLOCATED(fieldperiod%mdata%bcovar_s_hat)) &
         DEALLOCATE(fieldperiod%mdata%bcovar_s_hat)
    IF (ALLOCATED(fieldperiod%mdata%dlogbds)) &
         DEALLOCATE(fieldperiod%mdata%dlogbds)
    !! End Modifications by Andreas F. Martitsch (11.06.2014)
    IF (ALLOCATED(fieldperiod%mdata%ybeg)) DEALLOCATE(fieldperiod%mdata%ybeg)
    IF (ALLOCATED(fieldperiod%mdata%yend)) DEALLOCATE(fieldperiod%mdata%yend)
    DEALLOCATE(fieldperiod%coords)
    DEALLOCATE(fieldperiod%mdata)
    NULLIFY(fieldperiod%coords)
    NULLIFY(fieldperiod%mdata)

    !
    ! final deallocation
    NULLIFY(fieldperiod%parent)
    NULLIFY(fieldperiod%prev)
    NULLIFY(fieldperiod%next)
    NULLIFY(fieldperiod%prev_theta)
    NULLIFY(fieldperiod%next_theta)
    NULLIFY(fieldperiod%tag_child)
    DEALLOCATE(fieldperiod)
    NULLIFY(fieldperiod)
    ! message
    IF (mag_talk) PRINT *, 'magnetics: fieldperiod deleted: ',my_tag
  END SUBROUTINE destruct_mag_fieldperiod
  ! end constructor and destructor for fieldperiod_struct
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  !> \brief Constructor for fieldpropagator_struct.
  ! 
  SUBROUTINE construct_mag_fieldpropagator(fieldperiod,fieldpropagator,action)
    TYPE(fieldperiod_struct),     POINTER               :: fieldperiod
    TYPE(fieldpropagator_struct), POINTER               :: fieldpropagator
    CHARACTER(len=3),             OPTIONAL              :: action
    TYPE(inumber_struct),         POINTER               :: tag_child
    INTEGER                                             :: my_tag = 0
    CHARACTER(len=3)                                    :: act

    IF (PRESENT(action)) THEN 
       act = action
    ELSE
       act = 'non'
    END IF
    ! memory for fieldpropagator
    ALLOCATE(fieldpropagator)
    ! tag
    fieldpropagator_tag = fieldpropagator_tag + 1
    fieldpropagator%tag = fieldpropagator_tag
    my_tag              = fieldpropagator_tag

    IF (act .EQ. 'non') THEN
       ALLOCATE(tag_child) ! new memory for a number
       tag_child%i = my_tag
       ! print *, 'setup propagator'
       ! connect
       IF (ASSOCIATED(fieldperiod%ch_act)) THEN
          ! print *, 'fieldperiod%ch_fir already defined'
          fieldperiod%ch_act%next    => fieldpropagator
          fieldpropagator%prev       => fieldperiod%ch_act
          fieldperiod%tag_child%next => tag_child
          tag_child%prev             => fieldperiod%tag_child
          fieldperiod%tag_child      => tag_child
       ELSE
          ! print *, 'fieldperiod%ch_fir not defined'
          IF (ASSOCIATED(fieldperiod%prev)) THEN 
             ! connect to propagator in previous period
             fieldperiod%prev%ch_las%next => fieldpropagator 
             fieldpropagator%prev         => fieldperiod%prev%ch_las
          END IF
          fieldperiod%ch_fir    => fieldpropagator
          fieldperiod%tag_child => tag_child
       END IF
       fieldperiod%ch_act     => fieldpropagator
       fieldperiod%ch_las     => fieldpropagator
       fieldpropagator%parent => fieldperiod
       NULLIFY(tag_child) ! destroy pointer
    ELSE IF (act .EQ. 'ext') THEN
       ! extra propagator at beginning of end
       ! not linked to other propagators through prev and next
       fieldpropagator%parent => fieldperiod
       fieldperiod%ch_ext     => fieldpropagator       
    END IF
    ! additional quantities
    !
    ! end additional quantities    
    IF (mag_talk) PRINT *, 'magnetics: fieldpropagator added: ',my_tag, &
         ' parent: ',fieldpropagator%parent%tag
    RETURN
  END SUBROUTINE construct_mag_fieldpropagator
  !> \brief Destructor for fieldpropagator_struct.
  SUBROUTINE destruct_mag_fieldpropagator(fieldpropagator,action)
    TYPE(fieldpropagator_struct),   POINTER  :: fieldpropagator
    TYPE(fieldripple_struct),       POINTER  :: fieldripple
    CHARACTER(len=3),               OPTIONAL :: action
    INTEGER                                  :: my_tag = 0
    INTEGER                                  :: tag_prev,tag_next
    CHARACTER(len=3)                         :: act

    IF (PRESENT(action)) THEN
       act = action
    ELSE
       act = 'non'
    END IF

    ! destroy ripples first
    my_tag = fieldpropagator%tag
    IF (ASSOCIATED(fieldpropagator%ch_act)) THEN
       fieldripple => fieldpropagator%ch_act
       IF (ASSOCIATED(fieldpropagator%prev)) THEN
          tag_prev = fieldpropagator%prev%ch_act%tag
       ELSE IF (ASSOCIATED(fieldpropagator%parent%ch_ext)) THEN
          IF (ASSOCIATED(fieldpropagator%parent%ch_ext%ch_act)) THEN
             tag_prev = fieldpropagator%parent%ch_ext%ch_act%tag
          ELSE
             tag_prev = fieldripple%tag + 1
          END IF
       ELSE
          tag_prev = fieldripple%tag + 1
       END IF

       IF (ASSOCIATED(fieldpropagator%next)) THEN
          tag_next = fieldpropagator%next%ch_act%tag
       ELSE IF (ASSOCIATED(fieldpropagator%parent%ch_ext)) THEN
          IF (ASSOCIATED(fieldpropagator%parent%ch_ext%ch_act)) THEN
             tag_next = fieldpropagator%parent%ch_ext%ch_act%tag
          ELSE
             tag_next = fieldripple%tag + 1
          END IF
       ELSE
          tag_next = fieldripple%tag + 1
       END IF
       
       IF (fieldripple%tag .NE. tag_prev .AND.  &
            fieldripple%tag .NE. tag_next) THEN
          CALL destruct_magnetics(fieldripple)
       END IF
    END IF

    IF ( ASSOCIATED(fieldpropagator%prev) .AND. &
         ASSOCIATED(fieldpropagator%next) ) THEN
       fieldpropagator%parent%ch_act => fieldpropagator%next
       fieldpropagator%prev%next     => fieldpropagator%next
       fieldpropagator%next%prev     => fieldpropagator%prev
    ELSEIF ( ASSOCIATED(fieldpropagator%prev) .AND. & 
         .NOT. ASSOCIATED(fieldpropagator%next) ) THEN
       fieldpropagator%parent%ch_act => fieldpropagator%prev
       fieldpropagator%parent%ch_las => fieldpropagator%prev
       NULLIFY(fieldpropagator%prev%next)
    ELSEIF ( .NOT. ASSOCIATED(fieldpropagator%prev) .AND. &
         ASSOCIATED(fieldpropagator%next) ) THEN
       fieldpropagator%parent%ch_act => fieldpropagator%next
       fieldpropagator%parent%ch_fir => fieldpropagator%next
       NULLIFY(fieldpropagator%next%prev)
    ELSEIF ( .NOT. ASSOCIATED(fieldpropagator%prev) .AND. &
         .NOT. ASSOCIATED(fieldpropagator%next) ) THEN
       IF (act .EQ. 'non') THEN
          NULLIFY(fieldpropagator%parent%ch_act)
          NULLIFY(fieldpropagator%parent%ch_fir)
          NULLIFY(fieldpropagator%parent%ch_las)
       END IF
    END IF
    NULLIFY(fieldpropagator%ch_act)
    ! additional deallocation
    CALL set_magnetics_data(fieldpropagator)
    DEALLOCATE(fieldpropagator%coords)
    DEALLOCATE(fieldpropagator%mdata)
    NULLIFY(fieldpropagator%coords)
    NULLIFY(fieldpropagator%mdata)
    IF (ALLOCATED(fieldpropagator%phi_eta_ind)) &
         DEALLOCATE(fieldpropagator%phi_eta_ind)
!!$    ! additional deallocation (eta)
!!$    IF (ALLOCATED(fieldpropagator%eta)) DEALLOCATE(fieldpropagator%eta)
!!$    IF (ASSOCIATED(fieldpropagator%eta_x0)) THEN
!!$       CALL delete_all(fieldpropagator%eta_x0)
!!$       NULLIFY(fieldpropagator%eta_x0)
!!$    END IF
!!$    IF (ASSOCIATED(fieldpropagator%eta_s)) THEN
!!$       CALL delete_all(fieldpropagator%eta_s)
!!$       NULLIFY(fieldpropagator%eta_s)
!!$    END IF
!!$    ! binarysplit
!!$    CALL deconstruct_binarysplit(fieldpropagator%eta_bs)
    ! final deallocation
    NULLIFY(fieldpropagator%prev)
    NULLIFY(fieldpropagator%next)
    NULLIFY(fieldpropagator%parent)
    DEALLOCATE(fieldpropagator)
    NULLIFY(fieldpropagator)
    ! message
    IF (mag_talk) PRINT *, 'magnetics: fieldpropagator deleted: ',my_tag
  END SUBROUTINE destruct_mag_fieldpropagator
  ! end constructor and destructor for fieldpropagator_struct
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  !> \brief Constructor for fieldripple_struct.
  ! 
  SUBROUTINE construct_mag_fieldripple(fieldpropagator,fieldripple,action)
    TYPE(fieldpropagator_struct), POINTER               :: fieldpropagator
    TYPE(fieldripple_struct),     POINTER               :: fieldripple
    CHARACTER(len=3),             OPTIONAL              :: action
    INTEGER                                             :: my_tag = 0
    CHARACTER(len=3)                                    :: act

    IF (PRESENT(action)) THEN 
       act = action
    ELSE
       act = 'non'
    END IF

    IF (act .EQ. 'non') THEN
       ! memory for fieldripple
       ALLOCATE(fieldripple)
       ! tag
       fieldripple_tag = fieldripple_tag + 1
       fieldripple%tag = fieldripple_tag
       my_tag          = fieldripple_tag
       ! connect
       IF (ASSOCIATED(fieldpropagator%prev)) THEN
          IF (ASSOCIATED(fieldpropagator%prev%ch_act)) THEN 
             ! a prev ripple exists
             fieldpropagator%prev%ch_act%next => fieldripple
             fieldripple%prev                 => fieldpropagator%prev%ch_act
          END IF
       END IF
       fieldpropagator%ch_act       => fieldripple
       fieldripple%parent           => fieldpropagator
       fieldripple%parent%ch_tag    =  my_tag
    ELSE IF (act .EQ. 'add') THEN
       IF (ASSOCIATED(fieldpropagator%prev)) THEN 
          fieldpropagator%ch_act    => fieldpropagator%prev%ch_act
          my_tag                    =  fieldpropagator%prev%ch_act%tag
          fieldripple               => fieldpropagator%prev%ch_act
       ELSE IF (ASSOCIATED(fieldpropagator%parent%ch_ext%ch_act)) THEN
          fieldpropagator%ch_act    => fieldpropagator%parent%ch_ext%ch_act
          my_tag                    =  fieldpropagator%parent%ch_ext%ch_act%tag
          fieldripple               => fieldpropagator%parent%ch_ext%ch_act
       ELSE ! at the end
          fieldpropagator%ch_act    => fieldpropagator%parent%ch_las%ch_act
          my_tag                    =  fieldpropagator%parent%ch_las%ch_act%tag
          fieldripple               => fieldpropagator%parent%ch_las%ch_act
       END IF
       fieldripple%parent           => fieldpropagator
       fieldripple%parent%ch_tag    =  my_tag                 
    END IF
    ! additional quantities
    !
    ! end additional quantities
    IF (mag_talk) PRINT *, 'magnetics: fieldripple added: ',my_tag, &
         ' parent: ',fieldripple%parent%tag
    RETURN
  END SUBROUTINE construct_mag_fieldripple
  !> \brief Destructor for fieldripple_struct.
  SUBROUTINE destruct_mag_fieldripple(fieldripple)
    TYPE(fieldripple_struct), POINTER :: fieldripple
    INTEGER                           :: my_tag = 0
    my_tag = fieldripple%tag

    IF ( ASSOCIATED(fieldripple%prev) .AND. &
         ASSOCIATED(fieldripple%next) ) THEN
       fieldripple%parent%ch_act => fieldripple%next
       fieldripple%prev%next     => fieldripple%next
       fieldripple%next%prev     => fieldripple%prev
    ELSEIF ( ASSOCIATED(fieldripple%prev) .AND. &
         .NOT. ASSOCIATED(fieldripple%next) ) THEN
       fieldripple%parent%ch_act => fieldripple%prev
       NULLIFY(fieldripple%prev%next)
    ELSEIF ( .NOT. ASSOCIATED(fieldripple%prev) .AND. &
         ASSOCIATED(fieldripple%next) ) THEN
       fieldripple%parent%ch_act => fieldripple%next
       NULLIFY(fieldripple%next%prev)
    ELSEIF ( .NOT. ASSOCIATED(fieldripple%prev) .AND. &
         .NOT. ASSOCIATED(fieldripple%next) ) THEN
       NULLIFY(fieldripple%parent%ch_act)
    END IF
    ! additional deallocation (eta)
    IF (ALLOCATED(fieldripple%eta))     DEALLOCATE(fieldripple%eta)
    IF (ALLOCATED(fieldripple%eta_loc)) DEALLOCATE(fieldripple%eta_loc)
    IF (ASSOCIATED(fieldripple%eta_x0)) THEN
       CALL delete_all(fieldripple%eta_x0)
       NULLIFY(fieldripple%eta_x0)
    END IF
    IF (ASSOCIATED(fieldripple%eta_s)) THEN
       CALL delete_all(fieldripple%eta_s)
       NULLIFY(fieldripple%eta_s)
    END IF
    IF (ASSOCIATED(fieldripple%eta_cl)) THEN
       CALL delete_all(fieldripple%eta_cl)
       NULLIFY(fieldripple%eta_cl)
    END IF
    IF (ASSOCIATED(fieldripple%eta_shield)) THEN
       CALL delete_all(fieldripple%eta_shield)
       NULLIFY(fieldripple%eta_shield)
    END IF
    IF (ASSOCIATED(fieldripple%eta_type)) THEN
       CALL delete_all(fieldripple%eta_type)
       NULLIFY(fieldripple%eta_type)
    END IF
    ! additional deallocation inflection
    IF (ALLOCATED(fieldripple%phi_inflection))  DEALLOCATE(fieldripple%phi_inflection)
    IF (ALLOCATED(fieldripple%b_inflection))    DEALLOCATE(fieldripple%b_inflection)
    IF (ALLOCATED(fieldripple%dbdp_inflection)) DEALLOCATE(fieldripple%dbdp_inflection)
    ! binarysplit
    CALL deconstruct_binarysplit(fieldripple%eta_bs)
    CALL deconstruct_binarysplit(fieldripple%eta_bs_loc)
    ! final deallocation
    NULLIFY(fieldripple%parent)
    NULLIFY(fieldripple%prev)
    NULLIFY(fieldripple%next)    
    DEALLOCATE(fieldripple)
    NULLIFY(fieldripple)
    ! message
    IF (mag_talk) PRINT *, 'magnetics: fieldripple deleted: ',my_tag
  END SUBROUTINE destruct_mag_fieldripple
  ! end constructor and destructor for fieldripple_struct
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  !> \brief Screen output of basic information for device_struct.
  SUBROUTINE info_mag_device(device)
    TYPE(device_struct), POINTER :: device
    IF (mag_infotalk .AND. ASSOCIATED(device)) THEN
       PRINT *, '-----------------------------------------------------------'
       PRINT *, ' device with tag ',device%tag
       PRINT *, '  name           : ', device%name
       PRINT *, '  r0             : ', device%r0  
       PRINT *, '  z0             : ', device%z0    
       PRINT *, '  nfp            : ', device%nfp    
       IF (ASSOCIATED(device%ch_act)) THEN
          PRINT *, '  childs (f,a,l) : ', &
               device%ch_fir%tag,device%ch_act%tag,device%ch_las%tag
       ELSE
          PRINT *, '  childs (f,a,l) : ', 'none'
       END IF
       PRINT *, '-----------------------------------------------------------'       
    END IF
    IF (mag_talk .AND. .NOT. ASSOCIATED(device)) &
       PRINT *, 'magnetics info: device not associated'
    RETURN
  END SUBROUTINE info_mag_device
  ! ---------------------------------------------------------------------------
  !> \brief Screen output of basic information for surface_struct.
  SUBROUTINE info_mag_surface(surface)
    TYPE(surface_struct), POINTER :: surface
    IF (mag_infotalk .AND. ASSOCIATED(surface)) THEN
       PRINT *, '-----------------------------------------------------------'
       PRINT *, ' surface with tag ',surface%tag,' and parent ',surface%parent%tag
       PRINT *, '  bmod0          : ', surface%bmod0
       PRINT *, '  aiota          : ', surface%aiota
       PRINT *, '  nperiod        : ', surface%nperiod  
       PRINT *, '  nstep          : ', surface%nstep    
       PRINT *, '  ndim           : ', surface%ndim    
       PRINT *, '  b_abs_min      : ', surface%b_abs_min
       PRINT *, '  b_abs_max      : ', surface%b_abs_max
       IF (ASSOCIATED(surface%ch_act)) THEN
          PRINT *, '  childs (f,a,l) : ', &
               surface%ch_fir%tag,surface%ch_act%tag,surface%ch_las%tag
       ELSE
          PRINT *, '  childs (f,a,l) : ', 'none'
       END IF
       PRINT *, '  sister (f,a,l) : ', &
            surface%parent%ch_fir%tag, &
            surface%parent%ch_act%tag, &
            surface%parent%ch_las%tag
       PRINT *, '-----------------------------------------------------------'       
    END IF
    IF (mag_talk .AND. .NOT. ASSOCIATED(surface)) &
       PRINT *, 'magnetics info: surface not associated'
    RETURN
  END SUBROUTINE info_mag_surface
  ! ---------------------------------------------------------------------------
  !> \brief Screen output of basic information for fieldline_struct.
  SUBROUTINE info_mag_fieldline(fieldline)
    TYPE(fieldline_struct), POINTER :: fieldline
    IF (mag_infotalk .AND. ASSOCIATED(fieldline)) THEN
       PRINT *, '-----------------------------------------------------------'
       PRINT *, ' fieldline with tag ',fieldline%tag,' and parent ',fieldline%parent%tag
       PRINT *, '  xstart(1)      : ', fieldline%xstart(1)
       PRINT *, '  xstart(2)      : ', fieldline%xstart(2)
       PRINT *, '  xstart(3)      : ', fieldline%xstart(3)
       PRINT *, '  b_abs_min      : ', fieldline%b_abs_min
       PRINT *, '  b_abs_max      : ', fieldline%b_abs_max
       IF (ASSOCIATED(fieldline%ch_act)) THEN
          PRINT *, '  childs (f,a,l) : ', &
               fieldline%ch_fir%tag,fieldline%ch_act%tag,fieldline%ch_las%tag
       ELSE
          PRINT *, '  childs (f,a,l) : ', 'none'
       END IF
       PRINT *, '  sister (f,a,l) : ', &
            fieldline%parent%ch_fir%tag, &
            fieldline%parent%ch_act%tag, &
            fieldline%parent%ch_las%tag
       PRINT *, '-----------------------------------------------------------'       
    END IF
    IF (mag_talk .AND. .NOT. ASSOCIATED(fieldline)) &
       PRINT *, 'magnetics info: fieldline not associated'
    RETURN
  END SUBROUTINE info_mag_fieldline
  ! ---------------------------------------------------------------------------
  !> \brief Screen output of basic information for fieldperiod_struct.
  SUBROUTINE info_mag_fieldperiod(fieldperiod)
    TYPE(fieldperiod_struct), POINTER :: fieldperiod
    IF (mag_infotalk .AND. ASSOCIATED(fieldperiod)) THEN
       PRINT *, '-----------------------------------------------------------'
       PRINT *, ' fieldperiod with tag ',fieldperiod%tag, &
            ' and parent ',fieldperiod%parent%tag
       PRINT *, '  phi_l          : ', fieldperiod%phi_l
       PRINT *, '  phi_r          : ', fieldperiod%phi_r
       PRINT *, '  theta_b        : ', fieldperiod%theta_b
       IF (ASSOCIATED(fieldperiod%ch_act)) THEN
          PRINT *, '  childs (f,a,l) : ', &
               fieldperiod%ch_fir%tag,fieldperiod%ch_act%tag,fieldperiod%ch_las%tag
       ELSE
          PRINT *, '  childs (f,a,l) : ', 'none'
       END IF
       IF (ASSOCIATED(fieldperiod%ch_ext)) THEN
          PRINT *, '  childs (e)     : ', fieldperiod%ch_ext%tag
       END IF
       PRINT *, '  sister (f,a,l) : ', &
            fieldperiod%parent%ch_fir%tag, &
            fieldperiod%parent%ch_act%tag, &
            fieldperiod%parent%ch_las%tag
       ! extra information
       IF (ALLOCATED(fieldperiod%phi_ext)) THEN
          PRINT *, '  phi_ext        : ',fieldperiod%phi_ext
       ELSE
          PRINT *, '  not allocated phi_ext'
       END IF
       IF (ALLOCATED(fieldperiod%bhat_ext)) THEN
          PRINT *, '  bhat_ext       : ',fieldperiod%bhat_ext
       ELSE
          PRINT *, '  not allocated bhat_ext'
       END IF
       IF (ALLOCATED(fieldperiod%dbp_ext)) THEN
          PRINT *, '  dbp_ext        : ',fieldperiod%dbp_ext
       ELSE
          PRINT *, '  not allocated dbp_ext'
       END IF
       IF (ALLOCATED(fieldperiod%d2bp_ext)) THEN
          PRINT *, '  d2bp_ext       : ',fieldperiod%d2bp_ext
       ELSE
          PRINT *, '  not allocated d2bp_ext'
       END IF
       IF (ALLOCATED(fieldperiod%minmax)) THEN
          PRINT *, '  minmax         : ',fieldperiod%minmax
       ELSE
          PRINT *, '  not allocated minmax'
       END IF
       PRINT *, '-----------------------------------------------------------'       
    END IF
    IF (mag_talk .AND. .NOT. ASSOCIATED(fieldperiod)) &
       PRINT *, 'magnetics info: fieldperiod not associated'
    RETURN
  END SUBROUTINE info_mag_fieldperiod
  ! ---------------------------------------------------------------------------
  !> \brief Screen output of basic information for fieldpropagator_struct.
  SUBROUTINE info_mag_fieldpropagator(fieldpropagator)
    TYPE(fieldpropagator_struct), POINTER :: fieldpropagator
    INTEGER :: ic
    ! problem one could add sisters belonging to the same ripple
    IF (mag_infotalk .AND. ASSOCIATED(fieldpropagator)) THEN
       PRINT *, '-----------------------------------------------------------'
       PRINT *, ' fieldpropagator with tag ',fieldpropagator%tag, &
            ' and parent ',fieldpropagator%parent%tag
       PRINT *, '  phi_l          : ', fieldpropagator%phi_l
       PRINT *, '  phi_r          : ', fieldpropagator%phi_r
       PRINT *, '  phi_min        : ', fieldpropagator%phi_min
       PRINT *, '  b_l            : ', fieldpropagator%b_l
       PRINT *, '  b_r            : ', fieldpropagator%b_r
       PRINT *, '  b_min          : ', fieldpropagator%b_min
       PRINT *, '  i_min          : ', fieldpropagator%i_min
       PRINT *, '  has_min        : ', fieldpropagator%has_min
       IF (ASSOCIATED(fieldpropagator%ch_act)) THEN
          PRINT *, '  child          : ', fieldpropagator%ch_act%tag
       ELSE
          PRINT *, '  child          : ', 'none'
       END IF
       PRINT *, '  child_tag      : ', fieldpropagator%ch_tag
       
       PRINT *, '  sister (f,a,l) : ', &
            fieldpropagator%parent%ch_fir%tag, &
            fieldpropagator%parent%ch_act%tag, &
            fieldpropagator%parent%ch_las%tag
       PRINT *, '-----------------------------------------------------------'    
!!$       IF (.NOT. mag_split_ripple) THEN
!!$          PRINT *, '  bin_split_mode : ', fieldpropagator%bin_split_mode
!!$          IF (ALLOCATED(fieldpropagator%eta)) THEN
!!$             PRINT *, '  eta allocated  : ', & 
!!$                  LBOUND(fieldpropagator%eta),UBOUND(fieldpropagator%eta)
!!$          ELSE
!!$             PRINT *, '  eta allocated  : ', 'none'
!!$          END IF
!!$          IF (ASSOCIATED(fieldpropagator%eta_x0)) THEN
!!$             CALL goto_first(fieldpropagator%eta_x0)
!!$             CALL goto_last(fieldpropagator%eta_x0,ic)       
!!$             PRINT *, '  eta_x0 assoc.  : ', ic
!!$          ELSE
!!$             PRINT *, '  eta_x0 assoc.  : ', 'none'
!!$          END IF
!!$       END IF
       PRINT *, '-----------------------------------------------------------'    
       IF (ALLOCATED(fieldpropagator%phi_eta_ind)) THEN
          PRINT *, '  phi_eta_ind    : ', &
               LBOUND(fieldpropagator%phi_eta_ind,1), &
               UBOUND(fieldpropagator%phi_eta_ind,1)
       ELSE
          PRINT *, '  phi_eta_ind    : ', 'not allocated'
       END IF
       IF (ASSOCIATED(fieldpropagator%coords)) THEN
          PRINT *, '  coords         : ', 'associated'
          IF (ALLOCATED(fieldpropagator%coords%x1)) THEN
             PRINT *, '   allocated     : ', &
                  LBOUND(fieldpropagator%coords%x1,1), &
                  UBOUND(fieldpropagator%coords%x1,1)
          ELSE
             PRINT *, '   allocated     : ', 'none'
          END IF
       ELSE
          PRINT *, '  coords         : ', 'none'
       END IF
       IF (ASSOCIATED(fieldpropagator%mdata)) THEN
          PRINT *, '  mdata          : ', 'associated'
          IF (ALLOCATED(fieldpropagator%mdata%bhat)) THEN
             PRINT *, '   allocated     : ', &
                  LBOUND(fieldpropagator%mdata%bhat,1), &
                  UBOUND(fieldpropagator%mdata%bhat,1)
          ELSE
             PRINT *, '   allocated     : ', 'none'
          END IF
       ELSE
          PRINT *, '  coords         : ', 'none'
       END IF
    END IF
    IF (mag_talk .AND. .NOT. ASSOCIATED(fieldpropagator)) &
       PRINT *, 'magnetics info: fieldpropagator not associated'
    RETURN
  END SUBROUTINE info_mag_fieldpropagator
  ! ---------------------------------------------------------------------------
  !> \brief Screen output of basic information for fieldripple_struct.
  SUBROUTINE info_mag_fieldripple(fieldripple)
    TYPE(fieldripple_struct), POINTER :: fieldripple
    INTEGER :: ic
    ! problem one could add here number of props belonging to this ripple
    IF (mag_infotalk .AND. ASSOCIATED(fieldripple)) THEN
       PRINT *, '-----------------------------------------------------------'
       PRINT *, ' fieldripple with tag ',fieldripple%tag, &
            ' and parent ', fieldripple%parent%tag
       IF ( ASSOCIATED(fieldripple%pa_fir) ) THEN
          PRINT *, '  phi_l          : ', fieldripple%pa_fir%phi_l, &
               '(Tag: ',fieldripple%pa_fir%tag,')'
       END IF
       IF ( ASSOCIATED(fieldripple%pa_las) ) THEN
          PRINT *, '  phi_r          : ', fieldripple%pa_las%phi_r, &
               '(Tag: ',fieldripple%pa_las%tag,')'
       END IF
       PRINT *, '  b_max_l        : ', fieldripple%b_max_l
       PRINT *, '  b_max_r        : ', fieldripple%b_max_r
       PRINT *, '  b_min          : ', fieldripple%b_min
       PRINT *, '  d2bp_max_l     : ', fieldripple%d2bp_max_l
       PRINT *, '  d2bp_max_r     : ', fieldripple%d2bp_max_r
       PRINT *, '  d2bp_min       : ', fieldripple%d2bp_min  
       PRINT *, '  width          : ', fieldripple%width
       PRINT *, '  width_l        : ', fieldripple%width_l
       PRINT *, '  width_r        : ', fieldripple%width_r
       IF (mag_split_ripple) THEN
          PRINT *, '  bin_split_mode : ', fieldripple%bin_split_mode
          IF (ALLOCATED(fieldripple%eta)) THEN
             PRINT *, '  eta allocated  : ', & 
                  LBOUND(fieldripple%eta),UBOUND(fieldripple%eta)
          ELSE
             PRINT *, '  eta allocated  : ', 'none'
          END IF
          IF (ASSOCIATED(fieldripple%eta_x0)) THEN
             CALL goto_first(fieldripple%eta_x0)
             CALL goto_last(fieldripple%eta_x0,ic)       
             PRINT *, '  eta_x0 assoc.  : ', ic
          ELSE
             PRINT *, '  eta_x0 assoc.  : ', 'none'
          END IF
       END IF
       IF (ALLOCATED(fieldripple%phi_inflection)) THEN
          PRINT *, '  phi_inflection  : ',fieldripple%phi_inflection
       ELSE
          PRINT *, '  not allocated phi_inflection'
       END IF
       IF (ALLOCATED(fieldripple%b_inflection)) THEN
          PRINT *, '  b_inflection    : ',fieldripple%b_inflection
       ELSE
          PRINT *, '  not allocated b_inflection'
       END IF
       IF (ALLOCATED(fieldripple%dbdp_inflection)) THEN
          PRINT *, '  dbdp_inflection : ',fieldripple%dbdp_inflection
       ELSE
          PRINT *, '  not allocated dbdp_inflection'
       END IF
       PRINT *, '-----------------------------------------------------------'       
    END IF
    IF (mag_talk .AND. .NOT. ASSOCIATED(fieldripple)) &
       PRINT *, 'magnetics info: fieldripple not associated'
    RETURN
  END SUBROUTINE info_mag_fieldripple
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  SUBROUTINE h5_mag_general_d1a(name,var,groupname_1_opt,groupname_2_opt)
    character(len=*) :: name
    real(kind=dp), dimension(:), allocatable :: var
    !class(*) :: var
    character(len=*), optional :: groupname_1_opt,groupname_2_opt
    character(len=100) :: groupname_1,groupname_2
    integer(HID_T) :: h5_file_id,h5_group_id,h5_group_1_id,h5_group_2_id

    IF ( mag_write_hdf5 ) THEN

       if ( present(groupname_1_opt) ) then
          groupname_1 = groupname_1_opt
       else
          groupname_1 = 'general'
       end if

       call h5_open_rw(h5_magnetics_file_name, h5_file_id)

       ! open group name for device
       if ( h5_exists(h5_file_id, groupname_1) ) then
          CALL h5_open_group(h5_file_id, groupname_1, h5_group_1_id)
       else
          CALL h5_define_group(h5_file_id, groupname_1, h5_group_1_id)
       end if
       h5_group_id = h5_group_1_id

       if ( present(groupname_2_opt) ) then
          groupname_2 = groupname_2_opt
          if ( h5_exists(h5_group_1_id, groupname_2) ) then
             CALL h5_open_group(h5_group_1_id, groupname_2, h5_group_2_id)
          else
             CALL h5_define_group(h5_group_1_id, groupname_2, h5_group_2_id)
          end if
          h5_group_id = h5_group_2_id
       end if

       CALL h5_add(h5_group_id, name,  var)

       if ( present(groupname_2_opt) ) then
          call h5_close_group(h5_group_2_id)
       end if

       call h5_close_group(h5_group_1_id)

       call h5_close(h5_file_id)
    end IF
  end SUBROUTINE h5_mag_general_d1a
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  SUBROUTINE h5_mag_device(device,h5_file_id_opt,one_opt)
    TYPE(device_struct), POINTER :: device
    integer(HID_T), INTENT(inout), optional :: h5_file_id_opt
    logical, intent(in), optional :: one_opt

    integer(HID_T) :: h5_file_id
    integer(HID_T) :: h5_category_id
    TYPE(device_struct), POINTER :: h5_device

    integer(HID_T) :: h5_group_id
    CHARACTER(len=100) :: group_name
    logical :: one
    logical :: close_file
    INTEGER :: ios

    IF (mag_write_hdf5 .AND. ASSOCIATED(device)) THEN

       if ( present(h5_file_id_opt) ) then
          h5_file_id =h5_file_id_opt
       else
          h5_file_id = 0
       end if
       if ( present(one_opt) ) then
          one = one_opt
       else
          one = .false.
       end if

       h5_device => device

       ! deletes file at first call, then opens it with rw access
       if (.not. h5_isvalid(h5_file_id)) then
          ! inquire (file=h5_magnetics_file_name, exist=f_exists)
          ! how does delete of file now work in modern fortran
          !OPEN(unit=1234, iostat=ios, file=h5_magnetics_file_name, status='old')
          !IF (ios .EQ. 0) CLOSE(unit=1234, status='delete')
          close_file = .true.
          call h5_open_rw(h5_magnetics_file_name, h5_file_id)
       else
          close_file = .false.
       end if

       ! open group name for device
       if ( h5_exists(h5_file_id, h5_device_name) ) then
          CALL h5_open_group(h5_file_id, h5_device_name, h5_category_id)
       else
          CALL h5_define_group(h5_file_id, h5_device_name, h5_category_id)
       end if
       ! group name is constructed from tag
       WRITE(group_name,'(I6)') h5_device%tag
       WRITE(group_name,'(100A)') TRIM(ADJUSTL(group_name))
       if ( h5_exists(h5_category_id, group_name) ) then
          CALL h5_open_group(h5_category_id, group_name, h5_group_id)
       else
          CALL h5_define_group(h5_category_id, group_name, h5_group_id)
       end if

       CALL h5_add(h5_group_id, 'tag',  h5_device%tag)
       CALL h5_add(h5_group_id, 'name', h5_device%name)
       CALL h5_add(h5_group_id, 'r0',   h5_device%r0)
       CALL h5_add(h5_group_id, 'z0',   h5_device%z0)
       CALL h5_add(h5_group_id, 'nfp',  h5_device%nfp)
       call h5_close_group(h5_group_id)

       ! close group for devices
       call h5_close_group(h5_category_id)

       if ( .not.one ) then
          ! surface
          call h5_magnetics(h5_device%ch_act,h5_file_id)
          ! fieldline
          call h5_magnetics(h5_device%ch_act%ch_act,h5_file_id)
          ! fieldperiod
          call h5_magnetics(h5_device%ch_act%ch_act%ch_act,h5_file_id)
          ! fieldpropagator
          call h5_magnetics(h5_device%ch_act%ch_act%ch_fir%ch_fir,h5_file_id)
          ! fieldripple
          call h5_magnetics(h5_device%ch_act%ch_act%ch_fir%ch_fir%ch_act,h5_file_id)
       end if

       ! close hdf-fie
       if (close_file) call h5_close(h5_file_id)

    end IF
  end SUBROUTINE h5_mag_device
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  SUBROUTINE h5_mag_surface(surface,h5_file_id_opt,one_opt)
    TYPE(surface_struct), POINTER :: surface
    integer(HID_T), INTENT(inout), optional :: h5_file_id_opt
    logical, intent(in), optional :: one_opt

    integer(HID_T) :: h5_file_id
    integer(HID_T) :: h5_category_id
    TYPE(surface_struct), POINTER :: h5_surface

    integer(HID_T) :: h5_group_id
    CHARACTER(len=100) :: group_name
    logical :: one
    logical :: close_file
    INTEGER :: ios

    IF (mag_write_hdf5 .AND. ASSOCIATED(surface)) THEN

       if ( present(h5_file_id_opt) ) then
          h5_file_id =h5_file_id_opt
       else
          h5_file_id = 0
       end if
       if ( present(one_opt) ) then
          one = one_opt
       else
          one = .false.
       end if

       if (one) then
          h5_surface => surface
       else
          h5_surface => surface%parent%ch_fir
       end if

       ! deletes file at first call, then opens it with rw access
       if (.not. h5_isvalid(h5_file_id)) then
          ! inquire (file=h5_magnetics_file_name, exist=f_exists)
          ! how does delete of file now work in modern fortran
          !OPEN(unit=1234, iostat=ios, file=h5_magnetics_file_name, status='old')
          !IF (ios .EQ. 0) CLOSE(unit=1234, status='delete')
          close_file = .true.
          call h5_open_rw(h5_magnetics_file_name, h5_file_id)
       else
          close_file = .false.
       end if

       ! open group name for surface
       if ( h5_exists(h5_file_id, h5_surface_name) ) then
          CALL h5_open_group(h5_file_id, h5_surface_name, h5_category_id)
       else
          CALL h5_define_group(h5_file_id, h5_surface_name, h5_category_id)
       end if

       allsurfaces: do
          ! group name is constructed from tag
          WRITE(group_name,'(I6)') h5_surface%tag
          WRITE(group_name,'(100A)') TRIM(ADJUSTL(group_name))
          if ( h5_exists(h5_category_id, group_name) ) then
             CALL h5_open_group(h5_category_id, group_name, h5_group_id)
          else
             CALL h5_define_group(h5_category_id, group_name, h5_group_id)
          end if

          ! ------------------------------------------------------------------------
          CALL h5_add(h5_group_id, 'tag',    h5_surface%tag)
          CALL h5_add(h5_group_id, 'parent', h5_surface%parent%tag)
          CALL h5_add(h5_group_id, 'ch_act', h5_surface%ch_act%tag)
          CALL h5_add(h5_group_id, 'ch_fir', h5_surface%ch_fir%tag)
          CALL h5_add(h5_group_id, 'ch_las', h5_surface%ch_las%tag)
          ! ------------------------------------------------------------------------
          CALL h5_add(h5_group_id, 'bmod0',     h5_surface%bmod0, comment='reference magnetic field in Tesla', unit='T')
          CALL h5_add(h5_group_id, 'aiota',     h5_surface%aiota)
          CALL h5_add(h5_group_id, 'r_min',     h5_surface%r_min)
          CALL h5_add(h5_group_id, 'r_max',     h5_surface%r_max)
          CALL h5_add(h5_group_id, 'z_min',     h5_surface%z_min)
          CALL h5_add(h5_group_id, 'z_max',     h5_surface%z_max)
          CALL h5_add(h5_group_id, 'b_abs_min', h5_surface%b_abs_min)
          CALL h5_add(h5_group_id, 'b_abs_max', h5_surface%b_abs_max)
          CALL h5_add(h5_group_id, 'nperiod',   h5_surface%nperiod)
          CALL h5_add(h5_group_id, 'nstep',     h5_surface%nstep)
          CALL h5_add(h5_group_id, 'ndim',      h5_surface%ndim)
          ! ------------------------------------------------------------------------
          ! close group
          CALL h5_close_group(h5_group_id)

          ! got to next ripple or exit
          if ( .not.one .and. associated(h5_surface%next) ) then
             h5_surface => h5_surface%next
          else
             exit allsurfaces
          end if
       end do allsurfaces


       ! close group for surface
       call h5_close_group(h5_category_id)

       ! close hdf-fie
       if (close_file) call h5_close(h5_file_id)

    end IF
  end SUBROUTINE h5_mag_surface
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  SUBROUTINE h5_mag_fieldline(fieldline,h5_file_id_opt,one_opt)
    TYPE(fieldline_struct), POINTER :: fieldline
    integer(HID_T), INTENT(inout), optional :: h5_file_id_opt
    logical, intent(in), optional :: one_opt

    integer(HID_T) :: h5_file_id
    integer(HID_T) :: h5_category_id
    TYPE(fieldline_struct), POINTER :: h5_fieldline

    integer(HID_T) :: h5_group_id
    CHARACTER(len=100) :: group_name
    logical :: one
    logical :: close_file
    INTEGER :: ios

    IF (mag_write_hdf5 .AND. ASSOCIATED(fieldline)) THEN

       if ( present(h5_file_id_opt) ) then
          h5_file_id =h5_file_id_opt
       else
          h5_file_id = 0
       end if
       if ( present(one_opt) ) then
          one = one_opt
       else
          one = .false.
       end if

       if (one) then
          h5_fieldline => fieldline
       else
          h5_fieldline => fieldline%parent%ch_fir
       end if

       ! deletes file at first call, then opens it with rw access
       if (.not. h5_isvalid(h5_file_id)) then
          ! inquire (file=h5_magnetics_file_name, exist=f_exists)
          ! how does delete of file now work in modern fortran
          !OPEN(unit=1234, iostat=ios, file=h5_magnetics_file_name, status='old')
          !IF (ios .EQ. 0) CLOSE(unit=1234, status='delete')
          close_file = .true.
          call h5_open_rw(h5_magnetics_file_name, h5_file_id)
       else
          close_file = .false.
       end if

       ! open group name for fieldline
       if ( h5_exists(h5_file_id, h5_fieldline_name) ) then
          CALL h5_open_group(h5_file_id, h5_fieldline_name, h5_category_id)
       else
          CALL h5_define_group(h5_file_id, h5_fieldline_name, h5_category_id)
       end if

       allfieldlines: do
          ! group name is constructed from tag
          WRITE(group_name,'(I6)') h5_fieldline%tag
          WRITE(group_name,'(100A)') TRIM(ADJUSTL(group_name))
          if ( h5_exists(h5_category_id, group_name) ) then
             CALL h5_open_group(h5_category_id, group_name, h5_group_id)
          else
             CALL h5_define_group(h5_category_id, group_name, h5_group_id)
          end if

          ! ------------------------------------------------------------------------
          CALL h5_add(h5_group_id, 'tag',    h5_fieldline%tag)
          CALL h5_add(h5_group_id, 'parent', h5_fieldline%parent%tag)
          CALL h5_add(h5_group_id, 'ch_act', h5_fieldline%ch_act%tag)
          CALL h5_add(h5_group_id, 'ch_fir', h5_fieldline%ch_fir%tag)
          CALL h5_add(h5_group_id, 'ch_las', h5_fieldline%ch_las%tag)
          ! ------------------------------------------------------------------------
          CALL h5_add(h5_group_id, 'abs_max_ptag', h5_fieldline%abs_max_ptag)
          CALL h5_add(h5_group_id, 'abs_min_ptag', h5_fieldline%abs_min_ptag)
          CALL h5_add(h5_group_id, 'b_abs_min',    h5_fieldline%b_abs_min)
          CALL h5_add(h5_group_id, 'b_abs_max',    h5_fieldline%b_abs_max)
          CALL h5_add(h5_group_id, 'xstart',       h5_fieldline%xstart, &
               lbound(h5_fieldline%xstart),ubound(h5_fieldline%xstart))
         ! ------------------------------------------------------------------------
          ! close group
          CALL h5_close_group(h5_group_id)

          ! got to next fieldline or exit
          if ( .not.one .and. associated(h5_fieldline%next) ) then
             h5_fieldline => h5_fieldline%next
          else
             exit allfieldlines
          end if
       end do allfieldlines


       ! close group for fieldline
       call h5_close_group(h5_category_id)

       ! close hdf-fie
       if (close_file) call h5_close(h5_file_id)

    end IF
  end SUBROUTINE h5_mag_fieldline
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  SUBROUTINE h5_mag_fieldperiod(fieldperiod,h5_file_id_opt,one_opt)
    TYPE(fieldperiod_struct), POINTER :: fieldperiod
    integer(HID_T), INTENT(inout), optional :: h5_file_id_opt
    logical, intent(in), optional :: one_opt

    integer(HID_T) :: h5_file_id
    integer(HID_T) :: h5_category_id
    TYPE(fieldperiod_struct), POINTER :: h5_fieldperiod

    integer(HID_T) :: h5_group_id
    CHARACTER(len=100) :: group_name
    logical :: one
    logical :: close_file
    INTEGER :: ios

    IF (mag_write_hdf5 .AND. ASSOCIATED(fieldperiod)) THEN

       if ( present(h5_file_id_opt) ) then
          h5_file_id =h5_file_id_opt
       else
          h5_file_id = 0
       end if
       if ( present(one_opt) ) then
          one = one_opt
       else
          one = .false.
       end if

       if (one) then
          h5_fieldperiod => fieldperiod
       else
          h5_fieldperiod => fieldperiod%parent%ch_fir
       end if

       ! deletes file at first call, then opens it with rw access
       if (.not. h5_isvalid(h5_file_id)) then
          ! inquire (file=h5_magnetics_file_name, exist=f_exists)
          ! how does delete of file now work in modern fortran
          !OPEN(unit=1234, iostat=ios, file=h5_magnetics_file_name, status='old')
          !IF (ios .EQ. 0) CLOSE(unit=1234, status='delete')
          close_file = .true.
          call h5_open_rw(h5_magnetics_file_name, h5_file_id)
       else
          close_file = .false.
       end if

       ! open group name for fieldperiod
       if ( h5_exists(h5_file_id, h5_fieldperiod_name) ) then
          CALL h5_open_group(h5_file_id, h5_fieldperiod_name, h5_category_id)
       else
          CALL h5_define_group(h5_file_id, h5_fieldperiod_name, h5_category_id)
       end if

       allfieldperiods: do
          ! group name is constructed from tag
          WRITE(group_name,'(I6)') h5_fieldperiod%tag
          WRITE(group_name,'(100A)') TRIM(ADJUSTL(group_name))
          if ( h5_exists(h5_category_id, group_name) ) then
             CALL h5_open_group(h5_category_id, group_name, h5_group_id)
          else
             CALL h5_define_group(h5_category_id, group_name, h5_group_id)
          end if
          
          ! ------------------------------------------------------------------------
          CALL h5_add(h5_group_id, 'tag',        h5_fieldperiod%tag)
          CALL h5_add(h5_group_id, 'parent',     h5_fieldperiod%parent%tag)

          if ( associated(h5_fieldperiod%prev_theta) ) then
             CALL h5_add(h5_group_id, 'prev_theta', h5_fieldperiod%prev_theta%tag)
          else
             CALL h5_add(h5_group_id, 'prev_theta', 0)
          end if
          if ( associated(h5_fieldperiod%next_theta) ) then
             CALL h5_add(h5_group_id, 'next_theta', h5_fieldperiod%next_theta%tag)
          else
             CALL h5_add(h5_group_id, 'next_theta', 0)
          end if
          CALL h5_add(h5_group_id, 'ch_act',     h5_fieldperiod%ch_act%tag)
          CALL h5_add(h5_group_id, 'ch_fir',     h5_fieldperiod%ch_fir%tag)
          CALL h5_add(h5_group_id, 'ch_las',     h5_fieldperiod%ch_las%tag)
          if ( associated(h5_fieldperiod%ch_ext) ) then
             CALL h5_add(h5_group_id, 'ch_ext', h5_fieldperiod%ch_ext%tag)
          else
             CALL h5_add(h5_group_id, 'ch_ext', 0)
          end if
          CALL h5_add(h5_group_id, 'extra',      h5_fieldperiod%extra)
          ! inumber_struct :: tag_child
          ! ------------------------------------------------------------------------
          CALL h5_add(h5_group_id, 'phi_l',   h5_fieldperiod%phi_l)
          CALL h5_add(h5_group_id, 'phi_r',   h5_fieldperiod%phi_r)
          CALL h5_add(h5_group_id, 'theta_b', h5_fieldperiod%theta_b)
         ! ------------------------------------------------------------------------
          CALL h5_add(h5_group_id, 'phi_ext',     h5_fieldperiod%phi_ext)
          CALL h5_add(h5_group_id, 'bhat_ext',    h5_fieldperiod%bhat_ext)
          CALL h5_add(h5_group_id, 'dbp_ext',     h5_fieldperiod%dbp_ext)
          CALL h5_add(h5_group_id, 'd2bp_ext',    h5_fieldperiod%d2bp_ext)
          CALL h5_add(h5_group_id, 'minmax',      h5_fieldperiod%minmax)
          CALL h5_add(h5_group_id, 'width_left',  h5_fieldperiod%width_left)
          CALL h5_add(h5_group_id, 'width_right', h5_fieldperiod%width_right)
          ! ------------------------------------------------------------------------
          CALL h5_add(h5_group_id, 'x1',        h5_fieldperiod%coords%x1)
          CALL h5_add(h5_group_id, 'x2',        h5_fieldperiod%coords%x2)
          CALL h5_add(h5_group_id, 'x3',        h5_fieldperiod%coords%x3)
          CALL h5_add(h5_group_id, 'bhat',      h5_fieldperiod%mdata%bhat)
          CALL h5_add(h5_group_id, 'geodcu',    h5_fieldperiod%mdata%geodcu)
          CALL h5_add(h5_group_id, 'h_phi',     h5_fieldperiod%mdata%h_phi)
          CALL h5_add(h5_group_id, 'dlogbdphi', h5_fieldperiod%mdata%dlogbdphi)
          CALL h5_add(h5_group_id, 'ybeg',      h5_fieldperiod%mdata%ybeg)
          CALL h5_add(h5_group_id, 'yend',      h5_fieldperiod%mdata%yend)
          ! ------------------------------------------------------------------------
          ! close group
          CALL h5_close_group(h5_group_id)
          ! got to next fieldperiod or exit
          if ( .not.one .and. associated(h5_fieldperiod%next) ) then
             if ( h5_fieldperiod%next%tag .eq. fieldperiod%parent%ch_fir%tag) exit allfieldperiods
             h5_fieldperiod => h5_fieldperiod%next
          else
             exit allfieldperiods
          end if
       end do allfieldperiods

       ! close group for fieldperiod
       call h5_close_group(h5_category_id)

       ! close hdf-fie
       if (close_file) call h5_close(h5_file_id)

    end IF
  end SUBROUTINE h5_mag_fieldperiod
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  SUBROUTINE h5_mag_fieldpropagator(fieldpropagator,h5_file_id_opt,one_opt)
    TYPE(fieldpropagator_struct), POINTER :: fieldpropagator
    integer(HID_T), INTENT(inout), optional :: h5_file_id_opt
    logical, intent(in), optional :: one_opt

    integer(HID_T) :: h5_file_id
    integer(HID_T) :: h5_category_id
    TYPE(fieldpropagator_struct), POINTER :: h5_fieldpropagator

    integer(HID_T) :: h5_group_id
    CHARACTER(len=100) :: group_name
    logical :: one
    logical :: close_file
    INTEGER :: ios

    IF (mag_write_hdf5 .AND. ASSOCIATED(fieldpropagator)) THEN

       if ( present(h5_file_id_opt) ) then
          h5_file_id =h5_file_id_opt
       else
          h5_file_id = 0
       end if
       if ( present(one_opt) ) then
          one = one_opt
       else
          one = .false.
       end if

       if (one) then
          h5_fieldpropagator => fieldpropagator
       else
          h5_fieldpropagator => fieldpropagator%parent%parent%ch_fir%ch_fir
       end if

       ! deletes file at first call, then opens it with rw access
       if (.not. h5_isvalid(h5_file_id)) then
          ! inquire (file=h5_magnetics_file_name, exist=f_exists)
          ! how does delete of file now work in modern fortran
          !OPEN(unit=1234, iostat=ios, file=h5_magnetics_file_name, status='old')
          !IF (ios .EQ. 0) CLOSE(unit=1234, status='delete')
          close_file = .true.
          call h5_open_rw(h5_magnetics_file_name, h5_file_id)
       else
          close_file = .false.
       end if

       ! open group name for fieldpropagator
       if ( h5_exists(h5_file_id, h5_fieldpropagator_name) ) then
          CALL h5_open_group(h5_file_id, h5_fieldpropagator_name, h5_category_id)
       else
          CALL h5_define_group(h5_file_id, h5_fieldpropagator_name, h5_category_id)
       end if

       allfieldpropagators: do
          ! group name is constructed from tag
          WRITE(group_name,'(I6)') h5_fieldpropagator%tag
          WRITE(group_name,'(100A)') TRIM(ADJUSTL(group_name))
          if ( h5_exists(h5_category_id, group_name) ) then
             CALL h5_open_group(h5_category_id, group_name, h5_group_id)
          else
             CALL h5_define_group(h5_category_id, group_name, h5_group_id)
          end if
          
          ! ------------------------------------------------------------------------
          CALL h5_add(h5_group_id, 'tag',        h5_fieldpropagator%tag)
          CALL h5_add(h5_group_id, 'parent',     h5_fieldpropagator%parent%tag)
          CALL h5_add(h5_group_id, 'ch_tag',     h5_fieldpropagator%ch_tag)
          ! ------------------------------------------------------------------------
          CALL h5_add(h5_group_id, 'phi_l',   h5_fieldpropagator%phi_l)
          CALL h5_add(h5_group_id, 'phi_r',   h5_fieldpropagator%phi_r)
          CALL h5_add(h5_group_id, 'b_l',     h5_fieldpropagator%b_l)
          CALL h5_add(h5_group_id, 'b_r',     h5_fieldpropagator%b_r)
          CALL h5_add(h5_group_id, 'phi_min', h5_fieldpropagator%phi_min)
          CALL h5_add(h5_group_id, 'b_min',   h5_fieldpropagator%b_min)
          CALL h5_add(h5_group_id, 'i_min',   h5_fieldpropagator%i_min)
          CALL h5_add(h5_group_id, 'has_min', h5_fieldpropagator%has_min)
          ! ------------------------------------------------------------------------
          CALL h5_add(h5_group_id, 'phi_eta_ind', h5_fieldpropagator%phi_eta_ind)
          ! ------------------------------------------------------------------------
          CALL h5_add(h5_group_id, 'x1',        h5_fieldpropagator%coords%x1)
          CALL h5_add(h5_group_id, 'x2',        h5_fieldpropagator%coords%x2)
          CALL h5_add(h5_group_id, 'x3',        h5_fieldpropagator%coords%x3)
          CALL h5_add(h5_group_id, 'bhat',      h5_fieldpropagator%mdata%bhat)
          CALL h5_add(h5_group_id, 'geodcu',    h5_fieldpropagator%mdata%geodcu)
          CALL h5_add(h5_group_id, 'h_phi',     h5_fieldpropagator%mdata%h_phi)
          CALL h5_add(h5_group_id, 'dlogbdphi', h5_fieldpropagator%mdata%dlogbdphi)
          CALL h5_add(h5_group_id, 'ybeg',      h5_fieldpropagator%mdata%ybeg)
          CALL h5_add(h5_group_id, 'yend',      h5_fieldpropagator%mdata%yend)
          ! ------------------------------------------------------------------------
          ! close group
          CALL h5_close_group(h5_group_id)
          ! got to next fieldperiod or exit
          if ( .not.one .and. associated(h5_fieldpropagator%next) ) then
             if ( h5_fieldpropagator%next%tag .eq. fieldpropagator%parent%parent%ch_fir%ch_fir%tag) exit allfieldpropagators
             h5_fieldpropagator => h5_fieldpropagator%next
          else
             exit allfieldpropagators
          end if
       end do allfieldpropagators

       ! close group for fieldpropagator
       call h5_close_group(h5_category_id)

       ! close hdf-fie
       if (close_file) call h5_close(h5_file_id)

    end IF
  end SUBROUTINE h5_mag_fieldpropagator
  ! ---------------------------------------------------------------------------


  ! ---------------------------------------------------------------------------
  SUBROUTINE h5_mag_fieldripple(fieldripple,h5_file_id_opt,one_opt)
    TYPE(fieldripple_struct), POINTER :: fieldripple
    integer(HID_T), INTENT(inout), optional :: h5_file_id_opt
    logical, intent(in), optional :: one_opt

    integer(HID_T) :: h5_file_id
    integer(HID_T) :: h5_category_id
    TYPE(fieldripple_struct), POINTER :: h5_fieldripple
    TYPE(fieldpropagator_struct), POINTER :: fieldpropagator

    integer(HID_T) :: h5_group_id
    CHARACTER(len=100) :: group_name
    logical :: one
    logical :: close_file
    INTEGER :: ios

    REAL(kind=dp), allocatable, dimension(:) :: dummy_arr

    IF (mag_write_hdf5 .AND. ASSOCIATED(fieldripple)) THEN

       if ( present(h5_file_id_opt) ) then
          h5_file_id =h5_file_id_opt
       else
          h5_file_id = 0
       end if

       if ( present(one_opt) ) then
          one = one_opt
       else
          one = .false.
       end if

       if (one) then
          h5_fieldripple => fieldripple
       else
          h5_fieldripple => fieldripple%parent%parent%parent%ch_fir%ch_fir%ch_act
       end if

       ! deletes file at first call, then opens it with rw access
       if (.not. h5_isvalid(h5_file_id)) then
          ! inquire (file=h5_magnetics_file_name, exist=f_exists)
          ! how does delete of file now work in modern fortran
          !OPEN(unit=1234, iostat=ios, file=h5_magnetics_file_name, status='old')
          !IF (ios .EQ. 0) CLOSE(unit=1234, status='delete')
          call h5_open_rw(h5_magnetics_file_name, h5_file_id)
          close_file = .true.
       else
          close_file = .false.
       end if
       

       ! open group name for all fieldripples
       if ( h5_exists(h5_file_id, h5_fieldripple_name) ) then
          CALL h5_open_group(h5_file_id, h5_fieldripple_name, h5_category_id)
       else
          CALL h5_define_group(h5_file_id, h5_fieldripple_name, h5_category_id)
       end if


       allripples: do
          ! group name is constructed from tag
          WRITE(group_name,'(I6)') h5_fieldripple%tag
          WRITE(group_name,'(100A)') TRIM(ADJUSTL(group_name))
          if ( h5_exists(h5_category_id, group_name) ) then
             CALL h5_open_group(h5_category_id, group_name, h5_group_id)
          else
             CALL h5_define_group(h5_category_id, group_name, h5_group_id)
          end if

          ! add to group
          CALL h5_add(h5_group_id, 'parent', h5_fieldripple%parent%tag)
          CALL h5_add(h5_group_id, 'pa_fir', h5_fieldripple%pa_fir%tag)
          CALL h5_add(h5_group_id, 'pa_las', h5_fieldripple%pa_las%tag)
          CALL h5_add(h5_group_id, 'tag',    h5_fieldripple%tag)
          ! ------------------------------------------------------------------------
          CALL h5_add(h5_group_id, 'b_max_l',    h5_fieldripple%b_max_l)
          CALL h5_add(h5_group_id, 'b_max_r',    h5_fieldripple%b_max_r)
          CALL h5_add(h5_group_id, 'b_min',      h5_fieldripple%b_min)
          CALL h5_add(h5_group_id, 'd2bp_max_l', h5_fieldripple%d2bp_max_l)
          CALL h5_add(h5_group_id, 'd2bp_max_r', h5_fieldripple%d2bp_max_r)
          CALL h5_add(h5_group_id, 'd2bp_min',   h5_fieldripple%d2bp_min)
          CALL h5_add(h5_group_id, 'width',      h5_fieldripple%width)
          CALL h5_add(h5_group_id, 'width_l',    h5_fieldripple%width_l)
          CALL h5_add(h5_group_id, 'width_r',    h5_fieldripple%width_r)
          ! ------------------------------------------------------------------------
          CALL h5_add(h5_group_id, 'phi_inflection',  h5_fieldripple%phi_inflection)
          CALL h5_add(h5_group_id, 'b_inflection',    h5_fieldripple%b_inflection)
          CALL h5_add(h5_group_id, 'dbdp_inflection', h5_fieldripple%dbdp_inflection)         
          ! ------------------------------------------------------------------------
          CALL h5_add(h5_group_id, 'eta',            h5_fieldripple%eta)
          CALL h5_add(h5_group_id, 'eta_loc',        h5_fieldripple%eta_loc)
          CALL h5_add(h5_group_id, 'shielding_ll',   h5_fieldripple%shielding_ll)
          CALL h5_add(h5_group_id, 'shielding_lr',   h5_fieldripple%shielding_lr)
          CALL h5_add(h5_group_id, 'shielding_rr',   h5_fieldripple%shielding_rr)
          CALL h5_add(h5_group_id, 'shielding_rl',   h5_fieldripple%shielding_rl)
          CALL h5_add(h5_group_id, 'bin_split_mode', h5_fieldripple%bin_split_mode)
          ! ------------------------------------------------------------------------
          CALL extract_array( h5_fieldripple%eta_x0, dummy_arr, 1 )
          CALL h5_add(h5_group_id, 'eta_x0',     dummy_arr)
          CALL extract_array( h5_fieldripple%eta_s, dummy_arr, 1 )
          CALL h5_add(h5_group_id, 'eta_s',      dummy_arr)
          CALL extract_array( h5_fieldripple%eta_cl, dummy_arr, 1 )
          CALL h5_add(h5_group_id, 'eta_cl',     dummy_arr)
          CALL extract_array( h5_fieldripple%eta_shield, dummy_arr, 1 )
          CALL h5_add(h5_group_id, 'eta_shield', dummy_arr)
          CALL extract_array( h5_fieldripple%eta_type, dummy_arr, 1 )
          CALL h5_add(h5_group_id, 'eta_type',   dummy_arr)
          if ( allocated(dummy_arr) ) deallocate( dummy_arr )
          ! ------------------------------------------------------------------------
          CALL h5_add(h5_group_id, 'eta_boundary_left',               h5_fieldripple%eta_boundary_left)
          CALL h5_add(h5_group_id, 'eta_boundary_modification_left',  h5_fieldripple%eta_boundary_modification_left)
          CALL h5_add(h5_group_id, 'eta_boundary_index_left',         h5_fieldripple%eta_boundary_index_left)
          CALL h5_add(h5_group_id, 'eta_boundary_right',              h5_fieldripple%eta_boundary_right)
          CALL h5_add(h5_group_id, 'eta_boundary_modification_right', h5_fieldripple%eta_boundary_modification_right)
          CALL h5_add(h5_group_id, 'eta_boundary_index_right',        h5_fieldripple%eta_boundary_index_right)
          ! ------------------------------------------------------------------------
          ! additional stuff
          CALL h5_add(h5_group_id, 'phi_l', h5_fieldripple%pa_fir%phi_l)
          CALL h5_add(h5_group_id, 'phi_r', h5_fieldripple%pa_las%phi_r)
          fieldpropagator => h5_fieldripple%pa_fir
          b_min_search: do
             if ( fieldpropagator%has_min .ne. 0 ) then
                CALL h5_add(h5_group_id, 'phi_min', fieldpropagator%phi_min)
                exit b_min_search
             end if
             if ( associated(fieldpropagator%next) ) then
                fieldpropagator => fieldpropagator%next
             else
                exit b_min_search
             end if
          end do b_min_search
          ! ------------------------------------------------------------------------
          ! close group
          CALL h5_close_group(h5_group_id)

          ! got to next ripple or exit
          if ( .not.one .and. associated(h5_fieldripple%next) ) then
             h5_fieldripple => h5_fieldripple%next
          else
             exit allripples
          end if
       end do allripples

       ! close group for fieldripples
       call h5_close_group(h5_category_id)
       ! close hdf-file
       if (close_file) call h5_close(h5_file_id)

    end IF

!!$    ! problem one could add here number of props belonging to this ripple
!!$    IF (mag_infotalk .AND. ASSOCIATED(fieldripple)) THEN
!!$       PRINT *, ' fieldripple with tag ',fieldripple%tag, &
!!$            ' and parent ', fieldripple%parent%tag
!!$       IF ( ASSOCIATED(fieldripple%pa_fir) ) THEN
!!$          PRINT *, '  phi_l          : ', fieldripple%pa_fir%phi_l, &
!!$               '(Tag: ',fieldripple%pa_fir%tag,')'
!!$       END IF
!!$       IF ( ASSOCIATED(fieldripple%pa_las) ) THEN
!!$          PRINT *, '  phi_r          : ', fieldripple%pa_las%phi_r, &
!!$               '(Tag: ',fieldripple%pa_las%tag,')'
!!$       END IF
!!$       PRINT *, '  d2bp_max_l     : ', fieldripple%d2bp_max_l
!!$       PRINT *, '  d2bp_max_r     : ', fieldripple%d2bp_max_r
!!$       PRINT *, '  d2bp_min       : ', fieldripple%d2bp_min  
!!$       PRINT *, '  width          : ', fieldripple%width
!!$       PRINT *, '  width_l        : ', fieldripple%width_l
!!$       PRINT *, '  width_r        : ', fieldripple%width_r
!!$       IF (mag_split_ripple) THEN
!!$          PRINT *, '  bin_split_mode : ', fieldripple%bin_split_mode
!!$          IF (ALLOCATED(fieldripple%eta)) THEN
!!$             PRINT *, '  eta allocated  : ', & 
!!$                  LBOUND(fieldripple%eta),UBOUND(fieldripple%eta)
!!$          ELSE
!!$             PRINT *, '  eta allocated  : ', 'none'
!!$          END IF
!!$          IF (ASSOCIATED(fieldripple%eta_x0)) THEN
!!$             CALL goto_first(fieldripple%eta_x0)
!!$             CALL goto_last(fieldripple%eta_x0,ic)       
!!$             PRINT *, '  eta_x0 assoc.  : ', ic
!!$          ELSE
!!$             PRINT *, '  eta_x0 assoc.  : ', 'none'
!!$          END IF
!!$       END IF
!!$       IF (ALLOCATED(fieldripple%phi_inflection)) THEN
!!$          PRINT *, '  phi_inflection  : ',fieldripple%phi_inflection
!!$       ELSE
!!$          PRINT *, '  not allocated phi_inflection'
!!$       END IF
!!$       IF (ALLOCATED(fieldripple%b_inflection)) THEN
!!$          PRINT *, '  b_inflection    : ',fieldripple%b_inflection
!!$       ELSE
!!$          PRINT *, '  not allocated b_inflection'
!!$       END IF
!!$       IF (ALLOCATED(fieldripple%dbdp_inflection)) THEN
!!$          PRINT *, '  dbdp_inflection : ',fieldripple%dbdp_inflection
!!$       ELSE
!!$          PRINT *, '  not allocated dbdp_inflection'
!!$       END IF
!!$       PRINT *, '-----------------------------------------------------------'       
!!$    END IF
!!$    IF (mag_talk .AND. .NOT. ASSOCIATED(fieldripple)) &
!!$         PRINT *, 'magnetics info: fieldripple not associated'

    RETURN
  END SUBROUTINE h5_mag_fieldripple
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  ! plotting
  SUBROUTINE plot_mag_fieldperiod(fieldperiod,tags_in,tage_in,name_in)
    TYPE(fieldperiod_struct), POINTER :: fieldperiod
    INTEGER, OPTIONAL, INTENT(in) :: tags_in,tage_in
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: name_in

    TYPE(fieldpropagator_struct), POINTER :: fieldpropagator
    INTEGER :: tags,tage,unit,ic
    LOGICAL :: opened
    CHARACTER(len=100) :: name

    IF (.NOT. ASSOCIATED(fieldperiod)) THEN
       IF (mag_talk) PRINT *, 'magnetics plot: fieldperiod not associated '
       RETURN
    END IF
    IF (.NOT. PRESENT(tags_in)) tags = fieldperiod%tag
    IF (      PRESENT(tags_in)) tags = tags_in
    IF (.NOT. PRESENT(tage_in)) tage = fieldperiod%tag
    IF (      PRESENT(tage_in)) tage = tage_in
    IF (.NOT. PRESENT(name_in)) name = fieldperiod_name
    IF (      PRESENT(name_in)) name = name_in

    ! find free unit
    unit = 100
    DO
       INQUIRE(unit=unit,opened=opened)
       IF(.NOT. opened) EXIT
       unit = unit + 1
    END DO

    ! walk through fieldperiods
    fieldperiod => fieldperiod%parent%ch_fir
    ic = 1 ! mark the first propagator which is plotted
    DO 
       ! relevant periods
       IF (fieldperiod%tag .GE. tags .AND. fieldperiod%tag .LE. tage) THEN
          IF (ASSOCIATED(fieldperiod%ch_fir)) THEN
             ! walk through propagators of period
             fieldpropagator => fieldperiod%ch_fir ! first
             DO
                IF (ic .EQ. 1) THEN
                   CALL plot_magnetics(fieldpropagator,unit,name)
                   ic = ic + 1
                ELSE
                   CALL plot_magnetics(fieldpropagator,unit)
                END IF
                ! go to the next propagator or exit
                IF (ASSOCIATED(fieldpropagator%next)) THEN
                   IF (fieldpropagator%next%tag .LE. fieldperiod%ch_las%tag) THEN
                      fieldpropagator => fieldpropagator%next
                   ELSE
                      EXIT
                   END IF
                ELSE
                   EXIT
                END IF
             END DO
          END IF
       ELSE IF(fieldperiod%tag .GT. tage) THEN
          ! periods are processed
          EXIT
       END IF
       ! next period or exit
       IF (ASSOCIATED(fieldperiod%next)) THEN
          fieldperiod => fieldperiod%next
       ELSE
          IF (mag_talk) PRINT *, &
               'magnetics plot: fieldperiod tag not found: ',tage
          EXIT
       END IF
    END DO
    ! final closing
    INQUIRE(unit=unit,opened=opened)
    IF (opened) CLOSE(unit)
    !
    RETURN
  END SUBROUTINE plot_mag_fieldperiod
  ! ---------------------------------------------------------------------------
  SUBROUTINE plot_mag_fieldpropagator_tag(fieldpropagator,tags,tage,name_in)
    TYPE(fieldpropagator_struct), POINTER    :: fieldpropagator
    INTEGER,                      INTENT(in) :: tags,tage
    CHARACTER(len=*), OPTIONAL,   INTENT(in) :: name_in

    INTEGER :: i,unit
    LOGICAL :: opened
    CHARACTER(len=100) :: name
    
    IF (.NOT. PRESENT(name_in)) name = fieldpropagator_name
    IF (      PRESENT(name_in)) name = name_in

    ! find free unit
    unit = 100
    DO
       INQUIRE(unit=unit,opened=opened)
       IF(.NOT. opened) EXIT
       unit = unit + 1
    END DO

    ! problem: extra propagatoprs make problems here

    ! find propagator with starting tag
    fieldpropagator => fieldpropagator%parent%parent%ch_fir%ch_fir
    DO 
       IF (fieldpropagator%tag .EQ. tags) EXIT
       IF (.NOT. ASSOCIATED(fieldpropagator%next)) THEN
          IF (mag_talk) PRINT *, 'magnetics plot: tag not found: ',tags
          RETURN
       END IF
       fieldpropagator => fieldpropagator%next
    END DO

    ! loop for all fieldpropagators
    i = 0
    DO 
       i = i + 1
       IF (i .EQ. 1) THEN
          CALL plot_magnetics(fieldpropagator,unit,name)
       ELSE
          CALL plot_magnetics(fieldpropagator,unit)
       END IF
       IF (.NOT. ASSOCIATED(fieldpropagator%next)) THEN
          IF (mag_talk) PRINT *, 'magnetics plot: tag not found: ',tage
          EXIT
       END IF
       fieldpropagator => fieldpropagator%next
       IF (fieldpropagator%tag .GT. tage) EXIT
    END DO
    ! final closing
    INQUIRE(unit=unit,opened=opened)
    IF (opened) CLOSE(unit)
    RETURN
  END SUBROUTINE plot_mag_fieldpropagator_tag
  ! ---------------------------------------------------------------------------
  SUBROUTINE plot_mag_fieldpropagator(fieldpropagator,unit,name_in)
    TYPE(fieldpropagator_struct), POINTER    :: fieldpropagator
    INTEGER, INTENT(inout) :: unit
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: name_in
    !
    CHARACTER(len=100) :: name
    LOGICAL :: opened
    INTEGER :: i,is,ie

    IF (.NOT. PRESENT(name_in)) name = fieldpropagator_name
    IF (      PRESENT(name_in)) name = name_in

    INQUIRE(unit=unit,opened=opened)
    IF (ASSOCIATED(fieldpropagator)) THEN
       IF (ASSOCIATED(fieldpropagator%coords) .AND. &
            ASSOCIATED(fieldpropagator%mdata)) THEN
          IF (.NOT. opened) THEN
             OPEN(unit=unit,file=TRIM(name))
             is = 0
          ELSE
             is = 1
          END IF
          ie = UBOUND(fieldpropagator%coords%x1,1)
          DO i = is,ie
             !! Modifications by Andreas F. Martitsch (11.06.2014)
             ! Optional output (necessary for modeling the magnetic rotation)
             !
             IF ( ALLOCATED(fieldpropagator%mdata%dbcovar_s_hat_dphi) .AND. &
                  ALLOCATED(fieldpropagator%mdata%bcovar_s_hat)       .AND. &
                  ALLOCATED(fieldpropagator%mdata%dlogbds) ) THEN
                !WRITE(unit,'(1000(1x,e15.8))')               &
                WRITE(unit,*)                                 &
                     fieldpropagator%coords%x1(i),            &
                     fieldpropagator%coords%x2(i),            &
                     fieldpropagator%coords%x3(i),            &
                     fieldpropagator%mdata%bhat(i),           &
                     fieldpropagator%mdata%geodcu(i),         &
                     fieldpropagator%mdata%h_phi(i),          &
                     fieldpropagator%mdata%dlogbdphi(i),      &
                     fieldpropagator%mdata%dbcovar_s_hat_dphi,&
                     fieldpropagator%mdata%bcovar_s_hat,      &
                     fieldpropagator%mdata%dlogbds
             ELSE ! This is the old version
                !WRITE(unit,'(1000(1x,e15.8))')               &
                WRITE(unit,*)                                 &
                     fieldpropagator%coords%x1(i),            &
                     fieldpropagator%coords%x2(i),            &
                     fieldpropagator%coords%x3(i),            &
                     fieldpropagator%mdata%bhat(i),           &
                     fieldpropagator%mdata%geodcu(i),         &
                     fieldpropagator%mdata%h_phi(i),          &
                     fieldpropagator%mdata%dlogbdphi(i)
             END IF
             !
             !! End Modifications by Andreas F. Martitsch (11.06.2014)
          END DO
       ELSE
          IF (mag_talk) PRINT *, &
               'magnetics plot: fieldpropagator data not associated: ', &
               fieldpropagator%tag
       END IF
    ELSE
       IF (mag_talk) PRINT *, &
            'magnetics plot: fieldpropagator not associated '
       
    END IF
  END SUBROUTINE plot_mag_fieldpropagator
  ! ---------------------------------------------------------------------------
  SUBROUTINE plot_mag_fieldripple(fieldripple,tags_in,tage_in,name_in)
    TYPE(fieldripple_struct), POINTER :: fieldripple
    INTEGER, OPTIONAL, INTENT(in) :: tags_in,tage_in
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: name_in

    TYPE(fieldperiod_struct), POINTER     :: fieldperiod
    TYPE(fieldpropagator_struct), POINTER :: fieldpropagator
    INTEGER :: tags,tage,unit,ic
    INTEGER :: ripple_tag,ext_first
    LOGICAL :: opened
    CHARACTER(len=100) :: name

    IF (.NOT. ASSOCIATED(fieldripple)) THEN
       IF (mag_talk) PRINT *, 'magnetics plot: fieldripple not associated '
       RETURN
    END IF
    IF (.NOT. PRESENT(tags_in)) tags = fieldripple%tag
    IF (      PRESENT(tags_in)) tags = tags_in
    IF (.NOT. PRESENT(tage_in)) tage = fieldripple%tag
    IF (      PRESENT(tage_in)) tage = tage_in
    IF (.NOT. PRESENT(name_in)) name = fieldripple_name
    IF (      PRESENT(name_in)) name = name_in

    ! find free unit
    unit = 100
    DO
       INQUIRE(unit=unit,opened=opened)
       IF(.NOT. opened) EXIT
       unit = unit + 1
    END DO

    ! walk through fieldripples
    fieldripple => fieldripple%parent%parent%parent%ch_fir%ch_fir%ch_act ! first
    ic = 1 ! mark the first propagator which is plotted
    ripple: DO 
       ! relevant ripple
       ripple_tag = fieldripple%tag
       IF (ripple_tag .GE. tags .AND. ripple_tag .LE. tage) THEN
          ! walk through propagators of period
          fieldpropagator => fieldripple%parent ! first
          fieldperiod     => fieldpropagator%parent
          ! extra propagator
          ext_first = 0
          IF (ASSOCIATED(fieldperiod%ch_ext)) THEN
             IF (fieldperiod%ch_ext%ch_tag .EQ. ripple_tag) THEN
                ext_first = 1
             END IF
          END IF
          ! find first propagator belonging to ripple
          fieldpropagator => fieldperiod%parent%ch_fir%ch_fir
          firstprop: DO 
             IF (fieldpropagator%ch_tag .EQ. ripple_tag) EXIT
             IF (ASSOCIATED(fieldpropagator%next)) THEN
                fieldpropagator => fieldpropagator%next
             ELSE
                IF (mag_talk) PRINT *, 'magnetics plot: propagator not found'
                EXIT ripple
             END IF
          END DO firstprop
          
          ! walk through props belonging to ripple
          prop: DO
             ! handle the first extra prop if necessary
             IF (ext_first .EQ. 1) THEN
                IF (fieldperiod%ch_ext%tag .LT. fieldpropagator%tag) THEN
                   CALL plot_magnetics(fieldperiod%ch_ext,unit,name)
                   ic = ic + 1
                END IF
             END IF
             IF (ic .EQ. 1) THEN
                CALL plot_magnetics(fieldpropagator,unit,name)
                ic = ic + 1
             ELSE
                CALL plot_magnetics(fieldpropagator,unit)
             END IF
             ! handle the last extra prop if necessary
             IF (ext_first .EQ. 1) THEN
                IF (fieldperiod%ch_ext%tag .GT. fieldpropagator%tag) THEN
                   IF (ic .EQ. 1) THEN
                      CALL plot_magnetics(fieldperiod%ch_ext,unit,name)
                      ic = ic + 1
                   ELSE
                      CALL plot_magnetics(fieldperiod%ch_ext,unit)
                   END IF
                END IF
             END IF
             ! go to the next propagator or exit
             IF (ASSOCIATED(fieldpropagator%next)) THEN
                IF (fieldpropagator%next%ch_tag .EQ. ripple_tag) THEN
                   fieldpropagator => fieldpropagator%next
                ELSE
                   EXIT prop
                END IF
             ELSE
                EXIT prop
             END IF
          END DO prop
       ELSE IF(fieldripple%tag .GT. tage) THEN
          ! ripples are processed
          EXIT
       END IF
       ! next ripple or exit
       IF (ASSOCIATED(fieldripple%next)) THEN
          fieldripple => fieldripple%next
       ELSE
          IF (mag_talk) PRINT *, &
               'magnetics plot: fieldripple tag not found: ',tage
          EXIT
       END IF
    END DO ripple
    ! final closing
    INQUIRE(unit=unit,opened=opened)
    IF (opened) CLOSE(unit)
    !
    RETURN
  END SUBROUTINE plot_mag_fieldripple
  ! end plot

  ! ---------------------------------------------------------------------------
  ! simple pointers
  SUBROUTINE disp_d1(dnumber)
    TYPE(dnumber_struct), POINTER :: dnumber
    IF (ASSOCIATED(dnumber)) THEN
       CALL goto_first(dnumber)
       DO
          PRINT *, dnumber%d
          IF (.NOT. ASSOCIATED(dnumber%next)) EXIT
          dnumber => dnumber%next
       END DO
       CALL goto_first(dnumber)    
    END IF
  END SUBROUTINE disp_d1
  ! ---------------------------------------------------------------------------
  SUBROUTINE disp_d2(dnumber,dnumber2)
    TYPE(dnumber_struct), POINTER :: dnumber,dnumber2
    IF (ASSOCIATED(dnumber) .AND. ASSOCIATED(dnumber2)) THEN
       CALL goto_first(dnumber)
       CALL goto_first(dnumber2)
       DO
          PRINT *, dnumber%d,dnumber2%d
          IF (.NOT. ASSOCIATED(dnumber%next) .OR. &
               .NOT. ASSOCIATED(dnumber2%next) ) EXIT
          dnumber  => dnumber%next
          dnumber2 => dnumber2%next
       END DO
       CALL goto_first(dnumber)
       CALL goto_first(dnumber2)
    END IF
  END SUBROUTINE disp_d2
  ! ---------------------------------------------------------------------------
  SUBROUTINE create_new_d(dnumber)
    TYPE(dnumber_struct), POINTER :: dnumber
    TYPE(dnumber_struct), POINTER :: nnumber => NULL()
    ALLOCATE(nnumber)

    nnumber%prev => NULL()
    nnumber%next => NULL()

    !PRINT * ,'ASSOCIATED(dnumber) ',ASSOCIATED(dnumber)
    !PRINT * ,'ASSOCIATED(nnumber) ',ASSOCIATED(nnumber)

    IF (ASSOCIATED(dnumber)) THEN
       
       !PRINT * ,'ASSOCIATED(dnumber%next) ',ASSOCIATED(dnumber%next)
       IF (ASSOCIATED(dnumber%next)) THEN
          !PRINT *, '1'
          !PRINT * ,'ASSOCIATED(nnumber) ',ASSOCIATED(nnumber)
          !PRINT * ,'ASSOCIATED(dnumber%next%prev) ',ASSOCIATED(dnumber%next%prev)
          !dnumber%next%prev => nnumber
          !dnumber => nnumber
          !PRINT *, '2'
          nnumber%prev      => dnumber
          !PRINT *, '3'
          nnumber%next      => dnumber%next
          !PRINT *, '4'
          dnumber%next%prev => nnumber
          !PRINT *, '5'
          dnumber%next => nnumber
       ELSE
          !PRINT *, '6'
          nnumber%prev      => dnumber
          !PRINT *, '6a'
          nnumber%next      => NULL()
          !PRINT *, '7'
          dnumber%next => nnumber
          !PRINT *, '8'
       END IF
       !PRINT *, '9'
       dnumber => dnumber%next
       !PRINT *, '10'
    ELSE
       !PRINT *, '11'
       dnumber => nnumber
       !PRINT *, '12'
    END IF
    !PRINT *, '13'
    NULLIFY(nnumber)
    !PRINT *, '14'
  END SUBROUTINE create_new_d
  ! ---------------------------------------------------------------------------
  SUBROUTINE create_before_d(dnumber)
    TYPE(dnumber_struct), POINTER :: dnumber
    TYPE(dnumber_struct), POINTER :: nnumber => NULL()
    ALLOCATE(nnumber)
    
    nnumber%prev => NULL()
    nnumber%next => NULL()
    
    IF (ASSOCIATED(dnumber)) THEN
       IF (ASSOCIATED(dnumber%prev)) THEN
          nnumber%prev      => dnumber%prev
          nnumber%next      => dnumber
          
          dnumber%prev%next => nnumber
          dnumber%prev      => nnumber
       ELSE
          nnumber%next      => dnumber
          dnumber%prev      => nnumber
       END IF
    END IF
    dnumber => nnumber
    NULLIFY(nnumber)
  END SUBROUTINE create_before_d
  ! ---------------------------------------------------------------------------
  SUBROUTINE set_new_d(dnumber,d)
    TYPE(dnumber_struct), POINTER :: dnumber
    REAL(kind=dp), INTENT(in) :: d
    CALL create_new(dnumber)
    !PRINT *, 'SET 1'
    !PRINT *, 'ass ',associated(dnumber)
    dnumber%d = d
    !PRINT *, 'SET 2'
  END SUBROUTINE set_new_d
  ! ---------------------------------------------------------------------------
  SUBROUTINE set_before_d(dnumber,d)
    TYPE(dnumber_struct), POINTER :: dnumber
    REAL(kind=dp), INTENT(in) :: d
    CALL create_before(dnumber)
    dnumber%d = d
  END SUBROUTINE set_before_d
  ! ---------------------------------------------------------------------------
  SUBROUTINE delete_one_d(dnumber)
    TYPE(dnumber_struct), POINTER :: dnumber
    TYPE(dnumber_struct), POINTER :: nnumber => NULL()
    IF (ASSOCIATED(dnumber)) THEN
       nnumber => dnumber
       IF ( ASSOCIATED(dnumber%prev) .AND. &
            ASSOCIATED(dnumber%next) ) THEN
          dnumber%prev%next     => dnumber%next
          dnumber%next%prev     => dnumber%prev
          dnumber               => dnumber%next
       ELSEIF ( ASSOCIATED(dnumber%prev) .AND. &
            .NOT. ASSOCIATED(dnumber%next) ) THEN
          NULLIFY(dnumber%prev%next)
          dnumber => dnumber%prev
       ELSEIF ( .NOT. ASSOCIATED(dnumber%prev) .AND. &
            ASSOCIATED(dnumber%next) ) THEN
          NULLIFY(dnumber%next%prev)
          dnumber => dnumber%next
       ELSEIF ( .NOT. ASSOCIATED(dnumber%prev) .AND. &
            .NOT. ASSOCIATED(dnumber%next) ) THEN
          NULLIFY(dnumber)
       END IF
       DEALLOCATE(nnumber)
       NULLIFY(nnumber)
    END IF
  END SUBROUTINE delete_one_d
  ! ---------------------------------------------------------------------------
  SUBROUTINE delete_all_d(dnumber)
    TYPE(dnumber_struct), POINTER :: dnumber
    !PRINT *, 'delete_all'
    IF (ASSOCIATED(dnumber)) THEN
       !PRINT *, 'A'
       CALL goto_last(dnumber)
       !PRINT *, 'B'
       DO
          !PRINT *, 'C'
          CALL delete_one(dnumber)
          !PRINT *, 'D'
          IF (.NOT. ASSOCIATED(dnumber)) EXIT
          IF (.NOT. ASSOCIATED(dnumber%prev)) EXIT
          !PRINT *, 'E'
          dnumber => dnumber%prev
          !PRINT *, 'F'
       END DO
    END IF
    !PRINT *, 'G'
    IF (ASSOCIATED(dnumber)) THEN
       !PRINT *, 'H'
       CALL delete_one(dnumber)
       !PRINT *, 'I'
    END IF
    !PRINT *, 'J'
  END SUBROUTINE delete_all_d
  ! ---------------------------------------------------------------------------
  SUBROUTINE goto_first_d(dnumber)
    TYPE(dnumber_struct), POINTER :: dnumber
    IF (ASSOCIATED(dnumber)) THEN
       DO
          IF (.NOT. ASSOCIATED(dnumber%prev)) EXIT
          dnumber => dnumber%prev
       END DO
    END IF
  END SUBROUTINE goto_first_d
  ! ---------------------------------------------------------------------------
  SUBROUTINE goto_last_d(dnumber,ic)
    TYPE(dnumber_struct), POINTER :: dnumber
    INTEGER, OPTIONAL :: ic
    INTEGER           :: i
    i = 0
    IF (ASSOCIATED(dnumber)) THEN
       i = 1
       DO
          IF (.NOT. ASSOCIATED(dnumber%next)) EXIT
          dnumber => dnumber%next
          i = i + 1 
       END DO
    END IF
    IF (PRESENT(ic)) ic = i
  END SUBROUTINE goto_last_d
  ! ---------------------------------------------------------------------------
  SUBROUTINE extract_array_d1(dnumber,arr,lb_in)
    TYPE(dnumber_struct), POINTER :: dnumber
    REAL(kind=dp), ALLOCATABLE, INTENT(inout) :: arr(:)
    INTEGER, OPTIONAL, INTENT(in) :: lb_in
    INTEGER :: ic,i,lb
    lb = 0
    IF (PRESENT(lb_in)) lb=lb_in

    IF (ALLOCATED(arr)) DEALLOCATE(arr)
    IF (ASSOCIATED(dnumber)) THEN
       CALL goto_first(dnumber)
       CALL goto_last(dnumber,ic)
       ALLOCATE(arr(lb:lb+ic-1))
       CALL goto_first(dnumber)
       i = lb
       DO 
          arr(i) = dnumber%d
          IF (.NOT. ASSOCIATED(dnumber%next)) EXIT
          dnumber => dnumber%next
          i = i + 1
       END DO
    END IF
  END SUBROUTINE extract_array_d1
  ! ---------------------------------------------------------------------------
  SUBROUTINE sort_d1(dnumber)
    TYPE(dnumber_struct), POINTER :: dnumber
    TYPE(dnumber_struct), POINTER :: w => NULL()
    TYPE(dnumber_struct), POINTER :: m => NULL()

    INTEGER       :: i
    REAL(kind=dp) :: min_val

    IF (ASSOCIATED(dnumber)) THEN
       CALL goto_first(dnumber)
       all: DO
          i = -1
          min_val = 1.0d199
          w => dnumber
          search: DO 
             IF (w%d .LE. min_val) THEN
                min_val = w%d
                m => w
                i = i + 1
             END IF
             IF (.NOT. ASSOCIATED(w%next)) EXIT search
             w => w%next
          END DO search
          ! put in new place
          IF (i .GT. 0) THEN
             CALL delete_one(m)
             CALL set_before(dnumber,min_val)
          END IF
          IF (.NOT. ASSOCIATED(dnumber%next)) EXIT all
          dnumber => dnumber%next
       END DO all
       NULLIFY(w)
       NULLIFY(m)
    END IF
  END SUBROUTINE sort_d1
  ! ---------------------------------------------------------------------------
  SUBROUTINE sort_d2(dnumber,dnumber2)
    TYPE(dnumber_struct), POINTER :: dnumber,dnumber2
    TYPE(dnumber_struct), POINTER :: w => NULL()
    TYPE(dnumber_struct), POINTER :: m => NULL()
    TYPE(dnumber_struct), POINTER :: w2 => NULL()
    TYPE(dnumber_struct), POINTER :: m2 => NULL()

    INTEGER       :: i
    REAL(kind=dp) :: min_val,min_val2

    min_val = 1.234e5
    min_val2 = 1.234e5

    IF (ASSOCIATED(dnumber) .AND. ASSOCIATED(dnumber2)) THEN
       CALL goto_first(dnumber )
       CALL goto_first(dnumber2)
       all: DO
          i = -1
          min_val = 1.0d199
          w  => dnumber
          w2 => dnumber2
          search: DO 
             IF (w%d .LE. min_val) THEN
                min_val  = w%d
                min_val2 = w2%d
                m  => w
                m2 => w2
                i = i + 1
             END IF
             IF (.NOT. ASSOCIATED(w%next)) EXIT search
             w  => w%next
             w2 => w2%next
          END DO search
          ! put in new place
          IF (i .GT. 0) THEN
             CALL delete_one(m )
             CALL delete_one(m2)
             CALL set_before(dnumber, min_val )
             CALL set_before(dnumber2,min_val2)
          END IF
          IF (.NOT. ASSOCIATED(dnumber%next)) EXIT all
          dnumber  => dnumber%next
          dnumber2 => dnumber2%next
       END DO all
       NULLIFY(w)
       NULLIFY(w2)
       NULLIFY(m)
       NULLIFY(m2)
    END IF
  END SUBROUTINE sort_d2
  ! end simple pointers
  ! ---------------------------------------------------------------------------


END MODULE magnetics_mod
