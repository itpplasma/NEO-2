!! Modification by Andreas F. Martitsch (14.07.2015)
MODULE ntv_mod
  !
  ! Module containing additional data from/for the NTV version
  ! of ripple_solver
  !
  ! Input (switches, mach number, normalized magnetic drift frequency,
  ! boozer_s, collisionality,...) provided by neo2.in
  !
  ! Used within programs:
  !   neo2 (main)
  !   ripple_solver
  !   diag_propagator_res
  !   neo_magfie_perturbation
  !
  ! module containing numerical constants
  USE neo_precision
  !
  IMPLICIT NONE
  !
  ! INPUT
  ! switch: turn on(=1)/off(=0) ntv mode (not used at the moment)
  INTEGER, PUBLIC :: isw_ntv_mode
  ! switch: 0=compute qflux only for the symmetric case; 1=do all computations
  INTEGER, PUBLIC :: isw_qflux_NA
  ! switch for rippler_solver versions
  ! (1=preconditioned; 2=Arnoldi Order 1; 3=Arnoldi Order 2)
  INTEGER, PUBLIC :: isw_ripple_solver
  ! name of perturbation file
  CHARACTER(len=100), PUBLIC :: in_file_pert
  ! toroidal mach number over R_major (Mt/R), Larmor radius associated with
  ! $B_{00}^{Booz}$ (rho_L_loc) times B
  REAL(kind=dp), PUBLIC :: MtOvR, B_rho_L_loc
  !
  ! OUTPUT
  ! value of the average ripple of the perturbation field
  REAL(kind=dp), PUBLIC :: eps_M_2_val
  ! value of the flux surface average of $g_{\varphi\varphi}$
  ! for symmetry flux coordinates
  REAL(kind=dp), PUBLIC :: av_gphph_val
  ! value of the flux surface average of $\frac{1}{B}$
  REAL(kind=dp), PUBLIC :: av_inv_bhat_val
  !
  ! LOCAL DEFINITIONS
  ! storage array for qflux_symm (If isw_qflux_symm=0, this quantity stores the
  ! qflux-matrix for the symmetric field. If isw_qflux_symm=1, ripple_solver
  ! returns to the calling routine after the computation of the qflux-matrix
  ! for the symmetric field and array for qflux_symm is not allocated!)
  REAL(kind=dp), DIMENSION(:,:),  ALLOCATABLE, PUBLIC :: qflux_symm
  ! starting point of field line for cylindrical coordinates
  ! (used for normalizations)
  REAL(kind=dp), PUBLIC :: xstart_cyl(3)
  ! toroidal mode number of the perturbation field
  INTEGER, PUBLIC :: m_phi
  !
  PUBLIC write_ntv_output
  PRIVATE write_ntv_output_a
  INTERFACE write_ntv_output
     MODULE PROCEDURE write_ntv_output_a
  END INTERFACE write_ntv_output
  !
CONTAINS
  !
  SUBROUTINE write_ntv_output_a(isw_qflux_NA_in,qflux_NA_in,ind_map_in,&
       beta_out_in,y_in,aiota_loc_in,rt0_in,avnabpsi_in)
    !
    USE mag_interface_mod, ONLY : mag_coordinates, &
         boozer_s, boozer_theta_beg, boozer_phi_beg
    USE partpa_mod,  ONLY : bmod0
    USE mag_sub, ONLY: mag
    USE neo_magfie_mod, ONLY: boozer_curr_pol_hat, boozer_psi_pr_hat
    USE collisionality_mod, ONLY : collpar
    !
    ! input:
    INTEGER, INTENT(in) :: isw_qflux_NA_in
    REAL(kind=dp), DIMENSION(:,:), INTENT(in) :: qflux_NA_in
    INTEGER, DIMENSION(3), INTENT(in) :: ind_map_in
    REAL(kind=dp), DIMENSION(3), INTENT(in) :: beta_out_in
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: y_in
    REAL(kind=dp), INTENT(in) :: aiota_loc_in, rt0_in, avnabpsi_in
    ! local definitions:
    ! D11 and D12 for the non-axisymmetric problem
    ! normalized with the plateau coefficient
    REAL(kind=dp) :: D11_NA_Dpl, D12_NA_Dpl
    ! D31 and D32 for the axisymmetric problem
    ! normalized with analytical value of D31
    REAL(kind=dp) :: D31_AX_D31ref, D32_AX_D31ref, k_cof
    ! D13 for the non-axisymmetric problem
    ! normalized with analytical value of D31
    REAL(kind=dp) :: D13_NA_D31ref
    ! indices, file id
    LOGICAL :: opened
    INTEGER :: i_p, j_p, uw
    ! conversion factors for normalization
    REAL(kind=dp) :: fac1, fac2, fac3
    ! normalized co-variant phi-component of B,  $\sqrt{g}B^\vartheta$
    REAL(kind=dp) :: bcovar_phi_hat, sqrtg_bctrvr_tht
    ! Only used to call mag for normalizations
    ! related to D31 and D32
    REAL(kind=dp)                 :: bmod_tmp,sqrtg_tmp
    REAL(kind=dp), DIMENSION(3)   :: x_tmp,bder_tmp,hcovar_tmp,hctrvr_tmp
    REAL(kind=dp), DIMENSION(3,3) :: hcoder_tmp,hctder_tmp
    ! Physical output (Mach number, collisionality, $\sqrt{g}B^\varphi$)
    REAL(kind=dp) :: Mt_val, nu_star, sqrtg_bctrvr_phi
    ! Physical output ($B_\varphi$,$B_\vartheta$,\langle{B^2}\rangle)
    REAL(kind=dp) :: bcovar_phi, bcovar_tht, avbhat2, avb2
    !
    ! computation of the normalization for D31 and D32 (-> D31_ref)
    IF (mag_coordinates .EQ. 0) THEN
       ! cylindrical coordinates
       x_tmp = xstart_cyl
       CALL mag(x_tmp,bmod_tmp,sqrtg_tmp,bder_tmp,hcovar_tmp,&
            hctrvr_tmp,hcoder_tmp,hctder_tmp)
       ! normalized co-variant phi-component of B
       ! (Note! co-variant phi-component of B is the
       ! same for cylindrical coordinates and symmetry flux
       ! coodrinates --> no conversion needed)
       bcovar_phi_hat = hcovar_tmp(2)*(bmod_tmp/bmod0)
       ! restore value of $\sqrt{g}B^\vartheta$ for
       ! symmetry flux coordinates from quantities
       ! given in cylindircal coordinates
       sqrtg_bctrvr_tht = avnabpsi_in*sqrtg_tmp*(hctrvr_tmp(3)*bmod_tmp*1.0e4_dp)
    ELSE
       ! boozer coordinates
       x_tmp = (/boozer_s,boozer_phi_beg,boozer_theta_beg/)
       CALL mag(x_tmp,bmod_tmp,sqrtg_tmp,bder_tmp,hcovar_tmp,&
            hctrvr_tmp,hcoder_tmp,hctder_tmp)
       !
       ! normalized co-variant phi-component of B
       ! (Note! There is no difference between the co-variant
       ! phi-component of B for Boozer coordinates and those for
       ! symmetry flux coordinates, which were used for the
       ! computation of the normalization of D31.)
       bcovar_phi_hat = hcovar_tmp(2)*(bmod_tmp/bmod0)
       ! actually this is the same as:
       ! bcovar_phi_hat = boozer_curr_pol_hat
       !
       ! restore value of $\sqrt{g}B^\vartheta$ for
       ! symmetry flux coordinates from the quantities
       ! given in Boozer coordinates
       sqrtg_bctrvr_tht = avnabpsi_in*aiota_loc_in*boozer_psi_pr_hat*(bmod0*1.0e4_dp)
       ! this is the same up to a minus sign resulting from the
       ! definition of sqrtg_tmp (right-handed system)
       !sqrtg_bctrvr_tht = avnabpsi_in*sqrtg_tmp*(hctrvr_tmp(3)*bmod_tmp*1.0d4)                
    END IF
    fac3 = - (2.0_dp/(bmod0*1.0e4_dp)) * beta_out_in(3) * beta_out_in(1) / y_in(6)
    !
    ! extra output for NTV computations
    IF(isw_qflux_NA_in .EQ. 1) THEN
       ! 1) axisymmetric and non-axisymmetric solution
       ! have been computed
       ! a) normalized diffusion coefficients for the
       ! non-axisymmetric case
       !
       ! normalization factor from plateau coefficient
       fac1=16.0_dp*rt0_in*aiota_loc_in/PI
       ! normalization factor from gamma matrices
       fac2= - beta_out_in(1) * beta_out_in(1) / y_in(6)
       ! convert indices for the gamma matrices according 
       ! to the paper Kernbichler(2008)
       i_p = ind_map_in(1)
       j_p = ind_map_in(1)
       D11_NA_Dpl=fac1*fac2*qflux_NA_in(i_p,j_p)
       i_p = ind_map_in(1)
       j_p = ind_map_in(2)
       D12_NA_Dpl=fac1*fac2*qflux_NA_in(i_p,j_p)
       !
       ! $D_{13}^{\rm NA}$ normalized with D31_ref
       i_p = ind_map_in(1)
       j_p = ind_map_in(3)
       D13_NA_D31ref = (fac3*sqrtg_bctrvr_tht/bcovar_phi_hat)*qflux_NA_in(i_p,j_p)
       !
       ! normalized diffusion coefficients for the
       ! axisymmetric case
       IF (ALLOCATED(qflux_symm)) THEN
          i_p = ind_map_in(3)
          j_p = ind_map_in(1)
          D31_AX_D31ref = (fac3*sqrtg_bctrvr_tht/bcovar_phi_hat)*qflux_symm(i_p,j_p)
          i_p = ind_map_in(3)
          j_p = ind_map_in(2)
          D32_AX_D31ref = (fac3*sqrtg_bctrvr_tht/bcovar_phi_hat)*qflux_symm(i_p,j_p)
          k_cof = (2.5_dp-D32_AX_D31ref/D31_AX_D31ref)
       ELSE
          ! in case of bounce-averaged model the solution
          ! for the axisymmetric problem is not computed 
          D31_AX_D31ref = 0.0_dp
          D32_AX_D31ref = 0.0_dp
          k_cof = 0.0_dp
       END IF
       !
    ELSE
       ! only axisymmteric solution has been computed
       ! (non-axisymmteric coefficients set to zero)
       D11_NA_Dpl=0.0d0
       D12_NA_Dpl=0.0d0
       D13_NA_D31ref=0.0d0
       ! in this case "qflux" is stored within the actual propagator "prop_a"
       i_p = ind_map_in(3)
       j_p = ind_map_in(1)
       D31_AX_D31ref = (fac3*sqrtg_bctrvr_tht/bcovar_phi_hat)*qflux_NA_in(i_p,j_p)
       i_p = ind_map_in(3)
       j_p = ind_map_in(2)
       D32_AX_D31ref = (fac3*sqrtg_bctrvr_tht/bcovar_phi_hat)*qflux_NA_in(i_p,j_p)
       k_cof = (2.5_dp-D32_AX_D31ref/D31_AX_D31ref)
    END IF
    !
    !PRINT *,D11_NA_Dpl, D12_NA_Dpl, D13_NA_D31ref
    !PRINT *,D31_AX_D31ref, D32_AX_D31ref
    !
    ! write output:
    !
    ! find free unit
    uw = 100
    DO
       INQUIRE(unit=uw,opened=opened)
       IF(.NOT. opened) EXIT
       uw = uw + 100
    END DO
    !
    !PRINT *,'1'
    ! physical output
    Mt_val=MtOvR*rt0_in
    nu_star=collpar*rt0_in/aiota_loc_in
    sqrtg_bctrvr_phi=sqrtg_bctrvr_tht/aiota_loc_in
    bcovar_phi=hcovar_tmp(2)*(bmod_tmp*1.0e4_dp)
    bcovar_tht=hcovar_tmp(3)*(bmod_tmp*1.0e4_dp)
    avbhat2=y_in(9)/y_in(6)
    avb2=avbhat2*((bmod_tmp*1.0e4_dp)**2)
    !
    !PRINT *,'2'
    OPEN(uw,file='ntv_out.dat',status='replace')
    WRITE (uw,'(1000(1x,e18.5))') &
         boozer_s, Mt_val, nu_star, B_rho_L_loc, &
         D31_AX_D31ref, D32_AX_D31ref, k_cof, &
         D11_NA_Dpl, D12_NA_Dpl, D13_NA_D31ref, &
         aiota_loc_in, rt0_in, (bmod_tmp*1.0e4_dp), &
         boozer_psi_pr_hat, avnabpsi_in, &
         sqrtg_bctrvr_tht, sqrtg_bctrvr_phi, bcovar_tht, bcovar_phi, &
         DBLE(m_phi), avbhat2, av_inv_bhat_val, eps_M_2_val, av_gphph_val
    CLOSE(uw)
    !PRINT *,'3'
    !
  END SUBROUTINE write_ntv_output_a
  !
END MODULE ntv_mod
!! End Modification by Andreas F. Martitsch (14.07.2015)

MODULE propagator_mod
  ! Module to handle Propagators for Neo2
  !
  ! External programs:
  !   subroutine ripple_solver
  !   subroutine join_ripples
  !
  ! Called through programs:
  !   subroutine propagator_solver
  !
  ! ATTENTION:
  !  
  !  If physical content is changed then pertinent changes have to
  !  be made in several places!
  !
  ! External quantities: <default>
  !  prop_timing      0: no timing; <1: timing of module>
  !  prop_diagnostic  0: no diagnostic; <1: normal>; 2: extented 
  !  prop_binary      0: normal joining; <1: binary joining> 
  !  
  ! TODO: Make it save for subsequent calls to flint
  !
  ! Winfried Kernbichler 19.08.2004
  !
  ! ---------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------
  ! mnemonics
  !  
  ! "_p" - means "plus"  - particles with positive velocity (left  to right)
  ! "_m" - means "minus" - particles with negative velocity (right to left )
  !
  ! ---------------------------------------------------------------------------
  ! dimensions
  !
  ! npart          - number of particles
  ! npass_l        - number of passing particles on the left  side
  ! npass_r        - number of passing particles on the right side 
  ! npart_halfband
  !
  ! ---------------------------------------------------------------------------
  ! flux and current (from the outside)
  !
  ! flux_p(npass_l) - particles enter with positive velocity from the left side 
  ! flux_m(npass_r) - particles enter with negative velocity from the left side
  ! curr_p(npass_l) - particles enter with positive velocity from the left side 
  ! curr_m(npass_r) - particles enter with negative velocity from the left side
  !
  ! ---------------------------------------------------------------------------
  ! sources (internal)
  !  this is a distribution of particles generated by
  !  internal sources in the ripple which are leaving the ripple
  !
  ! source_p_*(npass_r) - positive velocity leaving on right side
  ! source_m_*(npass_l) - negative velocity leaving on left  side
  !
  !  * can be g (gradient driven) or _e (parallel electric field driven)
  !
  ! ---------------------------------------------------------------------------
  ! matrices
  !
  ! amat_p_m(exit,entry) - enter with positive velocity (entry left) and
  !                        leave with negative velocity (exit  left)
  !
  ! amat_p_m(npass_l,npass_l)
  ! amat_m_p(npass_r,npass_r)
  ! amat_p_p(npass_r,npass_l)
  ! amat_m_m(npass_l,npass_r)
  !
  ! ---------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------


  USE binarysplit_mod

  IMPLICIT NONE
  
  ! ---------------------------------------------------------------------------
  ! private variables
  INTEGER, PARAMETER, PRIVATE :: dp = KIND(1.0D0)
  REAL(kind=dp),      PRIVATE :: time_o
  REAL(kind=dp),      PRIVATE :: time_co,time_jp,time_ja,time_jf,time_so
  REAL(kind=dp),      PRIVATE :: stime_co,stime_jp,stime_ja,stime_jf,stime_so
  REAL(kind=dp),      PRIVATE :: stime_tot = 0.0_dp, time_tot,time_tot_o
  REAL(kind=dp), PARAMETER, PRIVATE :: pi=3.14159265358979d0

  CHARACTER(len=9),   PARAMETER, PRIVATE         :: prop_format = 'formatted'
  CHARACTER(len=4),   PARAMETER, PRIVATE         :: prop_cext = 'prop'
  CHARACTER(len=6),   PARAMETER, PRIVATE         :: prop_cperiod = 'period'
  CHARACTER(len=10),  PARAMETER, PRIVATE         :: prop_cpropagator = 'propagator'
  CHARACTER(len=8),   PARAMETER, PRIVATE         :: prop_cboundary = 'boundary'
  CHARACTER(len=12),  PARAMETER, PRIVATE         :: prop_ctaginfo = 'taginfo.prop'
  CHARACTER(len=16),  PARAMETER, PRIVATE         :: prop_cresult = 'reconstruct'
  CHARACTER(len=100),            PRIVATE         :: prop_cfilename
  INTEGER,                       PRIVATE         :: prop_unit = 150
  INTEGER,                       PRIVATE         :: prop_first_tag = 0
  INTEGER,                       PRIVATE         :: prop_last_tag = 0



  ! ---------------------------------------------------------------------------


  ! public variables (input file)
  INTEGER,            PUBLIC  :: prop_diagphys 
  INTEGER,            PUBLIC  :: prop_overwrite 
  INTEGER,            PUBLIC  :: prop_diagnostic = 1
  INTEGER,            PUBLIC  :: prop_binary = 0
  INTEGER,            PUBLIC  :: prop_timing = 1
  INTEGER,            PUBLIC  :: prop_join_ends = 0
  INTEGER,            PUBLIC  :: prop_fluxsplitmode = 1
  INTEGER,            PUBLIC  :: prop_write = 0
  INTEGER,            PUBLIC  :: prop_reconstruct = 0
  INTEGER,            PUBLIC  :: prop_ripple_plot = 0
  ! usage for communication purposes
  INTEGER,            PUBLIC  :: prop_count_call = 0
  INTEGER,            PUBLIC  :: prop_ibegperiod = 1
  INTEGER,            PUBLIC  :: prop_modifyold = 1
  !
  ! fluxes after reconstruction
  REAL(kind=dp), ALLOCATABLE, PUBLIC  :: flux_mr(:,:),flux_pl(:,:)
  REAL(kind=dp),              PUBLIC  :: eta_modboundary_l,eta_modboundary_r
  INTEGER,                    PUBLIC  :: sw_first_prop,sw_last_prop
  !
  ! new switch for reconstruction of levels
  INTEGER,            PUBLIC  :: prop_reconstruct_levels = 0
  


  ! ---------------------------------------------------------------------------
  ! type declaration for propagator physics
  !  with gradients _g and parallel electric field _e
  PUBLIC prop_qe
  TYPE prop_qe
     INTEGER                                    :: npart
     INTEGER                                    :: npass_l
     INTEGER                                    :: npass_r
!->out     INTEGER                                    :: npart_halfband
     INTEGER                                    :: nvelocity                !<-in
     REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: amat_p_p
     REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: amat_m_m
     REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: amat_p_m
     REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: amat_m_p
!->out     REAL(kind=dp), DIMENSION(:),   ALLOCATABLE :: source_p_g
!->out     REAL(kind=dp), DIMENSION(:),   ALLOCATABLE :: source_m_g
!->out     REAL(kind=dp), DIMENSION(:),   ALLOCATABLE :: source_p_e
!->out     REAL(kind=dp), DIMENSION(:),   ALLOCATABLE :: source_m_e
     REAL(kind=dp), DIMENSION(:,:),   ALLOCATABLE :: source_p              !<-in
     REAL(kind=dp), DIMENSION(:,:),   ALLOCATABLE :: source_m              !<-in
!->out     REAL(kind=dp), DIMENSION(:),   ALLOCATABLE :: flux_p
!->out     REAL(kind=dp), DIMENSION(:),   ALLOCATABLE :: flux_m
!->out     REAL(kind=dp), DIMENSION(:),   ALLOCATABLE :: curr_p
!->out     REAL(kind=dp), DIMENSION(:),   ALLOCATABLE :: curr_m
     REAL(kind=dp), DIMENSION(:,:),   ALLOCATABLE :: flux_p                !<-in
     REAL(kind=dp), DIMENSION(:,:),   ALLOCATABLE :: flux_m                !<-in
!->out     REAL(kind=dp)                              :: qflux_g
!->out     REAL(kind=dp)                              :: qflux_e
!->out     REAL(kind=dp)                              :: qcurr_g
!->out     REAL(kind=dp)                              :: qcurr_e
     REAL(kind=dp), DIMENSION(:,:),       ALLOCATABLE :: qflux             !<-in
     !
     REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: cmat
     !
     ! eta at left and right boundary
     REAL(kind=dp), DIMENSION(:),   ALLOCATABLE :: eta_l
     REAL(kind=dp), DIMENSION(:),   ALLOCATABLE :: eta_r
     REAL(kind=dp)                              :: eta_boundary_l
     REAL(kind=dp)                              :: eta_boundary_r
     !
     ! working field for intermediate storage during reallocation
     REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: w
  END TYPE prop_qe
  
  PUBLIC prop_boundary
  TYPE prop_boundary
          REAL(kind=dp)  :: fieldpropagator_tag_left
          REAL(kind=dp)  :: fieldpropagator_tag_right
          REAL(kind=dp)  :: fieldperiod_tag_left
          REAL(kind=dp)  :: fieldperiod_tag_right
          REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: c_forward
          REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: c_backward
  END TYPE prop_boundary

  ! ---------------------------------------------------------------------------
  ! type declaration for propagator
  PUBLIC propagator
  TYPE propagator
     ! internal service pointers: Links to previous and next propagator
     TYPE(propagator),              POINTER     :: prev => NULL()
     TYPE(propagator),              POINTER     :: next => NULL()
     ! internal service quantity: Current usage of propagator
     !  not allocated   -3
     !  allocated       -2
     !  used            -1
     !  containing sums >= 0
     !   0: 2^0 - 1 period
     !   1: 2^1 - 2 periods
     !   2: 2^2 - 4 periods ....
     INTEGER                                    :: nr_joined = -3
     ! physical quantities
     ! have to be changed if pyhsics changes!
     !  see also: deallocate_propagator_cont    (deallocation)
     !            propagator_solver_int         (propagator solver)
     !            ripple_solver_int             (ripple solver)
     !            join_ripples_int              (joining ripples)
     !          
     INTEGER                                    :: bin_split_mode
     TYPE(binarysplit)                          :: eta_bs_l
     TYPE(binarysplit)                          :: eta_bs_r

     ! tag of fieldpropagator (start, end)
     INTEGER                                    :: fieldpropagator_tag_s
     INTEGER                                    :: fieldpropagator_tag_e
     INTEGER                                    :: fieldperiod_tag_s
     INTEGER                                    :: fieldperiod_tag_e
     ! y values from RK equation solver at end of propagator
     REAL(kind=dp), ALLOCATABLE                 :: y(:)
     ! phi value at left and right side
     REAL(kind=dp)                              :: phi_l
     REAL(kind=dp)                              :: phi_r
     ! real propagator content
     TYPE(prop_qe)                              :: p
  END TYPE propagator

  ! ---------------------------------------------------------------------------
  ! private propagators
  TYPE(propagator), POINTER, PRIVATE            :: prop_a ! actual for diag
  TYPE(propagator), POINTER, PRIVATE            :: prop_r ! root
  TYPE(propagator), POINTER, PUBLIC             :: prop_c ! current
  TYPE(propagator), POINTER, PRIVATE            :: prop_n ! new
  TYPE(propagator), POINTER, PRIVATE            :: prop_l ! last
  
  ! ---------------------------------------------------------------------------
  ! public routines for propagator handling
  PUBLIC  propagator_solver
  PRIVATE propagator_solver_loc
  INTERFACE propagator_solver
     MODULE PROCEDURE propagator_solver_loc
  END INTERFACE

  PUBLIC write_propagator_content
  PRIVATE write_propagator_cont
  INTERFACE write_propagator_content
     MODULE PROCEDURE write_propagator_cont
  END INTERFACE

  PUBLIC reconstruct_prop_dist
  PRIVATE reconstruct_propagator_dist,reconstruct_propagator_dist_1
  INTERFACE reconstruct_prop_dist
     MODULE PROCEDURE reconstruct_propagator_dist,reconstruct_propagator_dist_1
  END INTERFACE

  PUBLIC write_prop_bound_content
  PRIVATE write_prop_bound_cont
  INTERFACE write_prop_bound_content
     MODULE PROCEDURE write_prop_bound_cont
  END INTERFACE

  PUBLIC read_propagator_content
  PRIVATE read_propagator_cont
  INTERFACE read_propagator_content
     MODULE PROCEDURE read_propagator_cont
  END INTERFACE

  PUBLIC read_prop_bound_content
  PRIVATE read_prop_bound_cont
  INTERFACE read_prop_bound_content
     MODULE PROCEDURE read_prop_bound_cont
  END INTERFACE

  PUBLIC read_prop_recon_content
  PRIVATE read_prop_recon_cont
  INTERFACE read_prop_recon_content
     MODULE PROCEDURE read_prop_recon_cont
  END INTERFACE
  ! ---------------------------------------------------------------------------
  ! Private helpers
  PRIVATE unit_propagator
  PRIVATE unit_prop
  INTERFACE unit_propagator
     MODULE PROCEDURE unit_prop
  END INTERFACE

  PRIVATE filename_propagator
  PRIVATE filename_prop
  INTERFACE filename_propagator
     MODULE PROCEDURE filename_prop
  END INTERFACE

  ! ---------------------------------------------------------------------------
  ! Private routines for propagator handling
  PRIVATE construct_propagator
  PRIVATE construct_prop
  INTERFACE construct_propagator
     MODULE PROCEDURE construct_prop
  END INTERFACE
  ! ---------------------------------------------------------------------------
  PRIVATE destruct_propagator 
  PRIVATE destruct_prop
  INTERFACE destruct_propagator
     MODULE PROCEDURE destruct_prop
  END INTERFACE
  ! ---------------------------------------------------------------------------
  PRIVATE destruct_all_propagators 
  PRIVATE destruct_all_prop
  INTERFACE destruct_all_propagators
     MODULE PROCEDURE destruct_all_prop
  END INTERFACE
  ! ---------------------------------------------------------------------------
  PRIVATE deallocate_propagator_content
  PRIVATE deallocate_propagator_cont
  INTERFACE deallocate_propagator_content
     MODULE PROCEDURE deallocate_propagator_cont
  END INTERFACE
  ! ---------------------------------------------------------------------------
  !PUBLIC ASSIGNMENT(=)
  PUBLIC assign_propagator_content
  PRIVATE assign_propagator_cont,assign_propagator_cont_qe
  !INTERFACE ASSIGNMENT(=)
  INTERFACE assign_propagator_content
     MODULE PROCEDURE assign_propagator_cont,assign_propagator_cont_qe
  END INTERFACE
  ! ---------------------------------------------------------------------------
  PRIVATE diag_propagator_content
  PRIVATE diag_propagator_cont
  INTERFACE diag_propagator_content
     MODULE PROCEDURE diag_propagator_cont
  END INTERFACE
  ! ---------------------------------------------------------------------------
  PRIVATE diag_propagator_result
  PRIVATE diag_propagator_res
  INTERFACE diag_propagator_result
     MODULE PROCEDURE diag_propagator_res
  END INTERFACE
  ! ---------------------------------------------------------------------------
  PRIVATE diag_propagator_distrf
  PRIVATE diag_propagator_dis
  INTERFACE diag_propagator_distrf
     MODULE PROCEDURE diag_propagator_dis
  END INTERFACE
 
  ! ---------------------------------------------------------------------------
  ! Private interface to the ripple_solver and join_ripples subroutine
  PRIVATE ripple_solver_interface
  PRIVATE ripple_solver_int
  INTERFACE ripple_solver_interface
     MODULE PROCEDURE ripple_solver_int
  END INTERFACE
  ! ---------------------------------------------------------------------------
  PRIVATE plot_distrf_interface
  PRIVATE plot_distrf_int
  INTERFACE plot_distrf_interface
     MODULE PROCEDURE plot_distrf_int
  END INTERFACE
  ! ---------------------------------------------------------------------------
  PRIVATE join_ripples_interface
  PRIVATE join_ripples_int
  INTERFACE join_ripples_interface
     MODULE PROCEDURE join_ripples_int
  END INTERFACE
  ! ---------------------------------------------------------------------------

CONTAINS

  ! ---------------------------------------------------------------------------
  SUBROUTINE construct_prop(before_in)

    INTEGER, OPTIONAL :: before_in
    INTEGER           :: before

    IF (PRESENT(before_in)) THEN
       before = before_in
    ELSE
       before = 0
    END IF


    !
    ! Subroutine to setup a new propagator
    !
    ! At the moment prop_c has to be placed at the end and the
    ! new prop_c is associated with the new end of the linked list!
    ! 
    ! Allocate memory for a new node.
    ALLOCATE( prop_n )
    ! At the beginning of the loop root is not associated to some
    ! variable, so it has to be associated with the newnode.
    IF ( .NOT. ASSOCIATED( prop_r ) ) prop_r => prop_n
    ! If current is already associated its pointer next has to point
    ! towards the newly created node.
    ! The pointer prev of newnode has to point towards current.

    IF (before .EQ. 0) THEN
       IF ( ASSOCIATED( prop_c) ) THEN   
          prop_c%next => prop_n
          prop_n%prev => prop_c
       END IF
    ELSE
       IF ( ASSOCIATED( prop_c) ) THEN   
          prop_c%prev => prop_n
          prop_n%next => prop_c
       END IF
    END IF
    ! Now the current node is associated with the newly created node.
    prop_c => prop_n
    ! The last node is associated with the current node.
    prop_l => prop_c
    ! When newnode is nullified, the association between newnode and the
    ! associated memory is deleted but the memory is not deallocated.
    ! At this moment the pointer current points towards this memory.
    NULLIFY( prop_n )
    ! indicated that it is not being used at the moment
    prop_c%nr_joined = -2 
    !
  END SUBROUTINE construct_prop

  ! ---------------------------------------------------------------------------
  SUBROUTINE destruct_prop
    !
    ! Removes the current propagator prop_c 
    !
    ! Remove content
    CALL deallocate_propagator_content
    IF (.NOT. ASSOCIATED(prop_c%prev) .AND. ASSOCIATED(prop_c%next)) THEN
       ! Root is removed
       prop_r => prop_c%next
       prop_c%next%prev => NULL()
       DEALLOCATE( prop_c )
       prop_c => prop_r
    ELSEIF (ASSOCIATED(prop_c%prev) .AND. ASSOCIATED(prop_c%next)) THEN
       ! Middle is removed
       prop_n => prop_c
       prop_c%prev%next => prop_c%next
       prop_c%next%prev => prop_c%prev
       prop_c => prop_c%prev
       DEALLOCATE( prop_n )
       NULLIFY( prop_n )
    ELSEIF (ASSOCIATED(prop_c%prev) .AND. .NOT. ASSOCIATED(prop_c%next)) THEN
       ! Last is removed
       prop_l => prop_c%prev
       prop_c%prev%next => NULL()
       DEALLOCATE( prop_c )
       prop_c => prop_l
    ELSEIF (.NOT. ASSOCIATED(prop_c%prev) .AND. .NOT. ASSOCIATED(prop_c%next)) THEN
       ! Root and last are removed - nothing remains
       NULLIFY( prop_r )
       NULLIFY( prop_l )
       DEALLOCATE( prop_c )
       NULLIFY( prop_c )
    END IF
    !
  END SUBROUTINE destruct_prop

  ! ---------------------------------------------------------------------------
  SUBROUTINE destruct_all_prop
    !
    ! Subroutine to remove all propagators starting from the 
    ! last one (prop_l)
    prop_c => prop_l
    DO
       CALL destruct_propagator
       IF (.NOT. ASSOCIATED(prop_c)) EXIT
    END DO
    !
  END SUBROUTINE destruct_all_prop

  ! ---------------------------------------------------------------------------
  SUBROUTINE deallocate_propagator_cont()
    !
    ! Deallocates all parts of current propagator prop_c
    !
    ! Has to be changed if physical content of propagator changes!
    ! (see type declaration of propagator)
    !
    ! 2-D quantities
    IF (ALLOCATED(prop_c%p%amat_p_p)) DEALLOCATE(prop_c%p%amat_p_p)
    IF (ALLOCATED(prop_c%p%amat_m_m)) DEALLOCATE(prop_c%p%amat_m_m)
    IF (ALLOCATED(prop_c%p%amat_p_m)) DEALLOCATE(prop_c%p%amat_p_m)
    IF (ALLOCATED(prop_c%p%amat_m_p)) DEALLOCATE(prop_c%p%amat_m_p)
    !
    IF (ALLOCATED(prop_c%p%cmat)) DEALLOCATE(prop_c%p%cmat)
!->out    ! 1-D quantities
!->out    IF (ALLOCATED(prop_c%p%source_p_g)) DEALLOCATE(prop_c%p%source_p_g)
!->out    IF (ALLOCATED(prop_c%p%source_m_g)) DEALLOCATE(prop_c%p%source_m_g)
!->out    IF (ALLOCATED(prop_c%p%source_p_e)) DEALLOCATE(prop_c%p%source_p_e)
!->out    IF (ALLOCATED(prop_c%p%source_m_e)) DEALLOCATE(prop_c%p%source_m_e)
    IF (ALLOCATED(prop_c%p%source_p)) DEALLOCATE(prop_c%p%source_p)        !<-in
    IF (ALLOCATED(prop_c%p%source_m)) DEALLOCATE(prop_c%p%source_m)        !<-in

    IF (ALLOCATED(prop_c%p%flux_p))   DEALLOCATE(prop_c%p%flux_p)
    IF (ALLOCATED(prop_c%p%flux_m))   DEALLOCATE(prop_c%p%flux_m)

!->out    IF (ALLOCATED(prop_c%p%curr_p))   DEALLOCATE(prop_c%p%curr_p)
!->out    IF (ALLOCATED(prop_c%p%curr_m))   DEALLOCATE(prop_c%p%curr_m)

    IF (ALLOCATED(prop_c%p%qflux))   DEALLOCATE(prop_c%p%qflux)            !<-in

    ! 1-D quantities                                                       !<-in
    IF (ALLOCATED(prop_c%p%eta_l))   DEALLOCATE(prop_c%p%eta_l)
    IF (ALLOCATED(prop_c%p%eta_r))   DEALLOCATE(prop_c%p%eta_r)

    ! working array
    IF (ALLOCATED(prop_c%p%w)) DEALLOCATE(prop_c%p%w)
    ! y-vector (magnetics)
    IF (ALLOCATED(prop_c%y)) DEALLOCATE(prop_c%y)
    ! binarysplit
    CALL deconstruct_binarysplit(prop_c%eta_bs_l)
    CALL deconstruct_binarysplit(prop_c%eta_bs_r)
    !
  END SUBROUTINE deallocate_propagator_cont

  ! ---------------------------------------------------------------------------
  SUBROUTINE assign_propagator_cont(n,o)
    TYPE(propagator), POINTER  :: n
    TYPE(propagator), POINTER  :: o
    
    n%nr_joined             = o%nr_joined
    n%bin_split_mode        = o%bin_split_mode
    IF (o%bin_split_mode .EQ. 1) THEN
       n%eta_bs_l              = o%eta_bs_l
       n%eta_bs_r              = o%eta_bs_r
    END IF
    n%fieldpropagator_tag_s = o%fieldpropagator_tag_s
    n%fieldpropagator_tag_e = o%fieldpropagator_tag_e
    n%fieldperiod_tag_s     = o%fieldperiod_tag_s
    n%fieldperiod_tag_e     = o%fieldperiod_tag_e
    n%phi_l                 = o%phi_l
    n%phi_r                 = o%phi_r

    IF(ALLOCATED(o%y)) THEN
       ALLOCATE(n%y(LBOUND(o%y,1):UBOUND(o%y,1)))
       n%y                  = o%y
    END IF

    !n%p                     = o%p
    CALL assign_propagator_content(n%p,o%p)
  END SUBROUTINE assign_propagator_cont

  SUBROUTINE assign_propagator_cont_qe(pn,po)
    TYPE(prop_qe), INTENT(inout) :: pn
    TYPE(prop_qe), INTENT(in)    :: po

    pn%npart          = po%npart
    pn%npass_l        = po%npass_l
    pn%npass_r        = po%npass_r
!->out    pn%npart_halfband = po%npart_halfband
    pn%nvelocity       = po%nvelocity                                        !<-in

!->out    pn%qflux_g        = po%qflux_g
!->out    pn%qflux_e        = po%qflux_e
!->out    pn%qcurr_g        = po%qcurr_g
!->out    pn%qcurr_e        = po%qcurr_e
    IF (ALLOCATED(po%qflux)) THEN                                          !<-in
       ALLOCATE(pn%qflux(SIZE(po%qflux,1),SIZE(po%qflux,2)))               !<-in
       pn%qflux    = po%qflux                                              !<-in
    END IF                                                                 !<-in


    pn%eta_boundary_l = po%eta_boundary_l
    pn%eta_boundary_r = po%eta_boundary_r

    IF (ALLOCATED(po%amat_p_p)) THEN
       ALLOCATE(pn%amat_p_p(SIZE(po%amat_p_p,1),SIZE(po%amat_p_p,2)))
       pn%amat_p_p    = po%amat_p_p
    END IF
    IF (ALLOCATED(po%amat_m_m)) THEN
       ALLOCATE(pn%amat_m_m(SIZE(po%amat_m_m,1),SIZE(po%amat_m_m,2)))
       pn%amat_m_m    = po%amat_m_m
    END IF
    IF (ALLOCATED(po%amat_p_m)) THEN
       ALLOCATE(pn%amat_p_m(SIZE(po%amat_p_m,1),SIZE(po%amat_p_m,2)))
       pn%amat_p_m    = po%amat_p_m
    END IF
    IF (ALLOCATED(po%amat_m_p)) THEN
       ALLOCATE(pn%amat_m_p(SIZE(po%amat_m_p,1),SIZE(po%amat_m_p,2)))
       pn%amat_m_p    = po%amat_m_p
    END IF

    IF (ALLOCATED(po%cmat)) THEN
       ALLOCATE(pn%cmat(SIZE(po%cmat,1),SIZE(po%cmat,2)))
       pn%cmat        = po%cmat
    END IF
    
    IF (ALLOCATED(po%w)) THEN
       ALLOCATE(pn%w(SIZE(po%w,1),SIZE(po%w,2)))
       pn%w           = po%w
    END IF
    
!->out    IF (ALLOCATED(po%source_p_g)) THEN
!->out       ALLOCATE(pn%source_p_g(SIZE(po%source_p_g,1)))
!->out       pn%source_p_g  = po%source_p_g
!->out    END IF
!->out    IF (ALLOCATED(po%source_m_g)) THEN
!->out       ALLOCATE(pn%source_m_g(SIZE(po%source_m_g,1)))
!->out       pn%source_m_g  = po%source_m_g
!->out    END IF
!->out    IF (ALLOCATED(po%source_p_e)) THEN
!->out       ALLOCATE(pn%source_p_e(SIZE(po%source_p_e,1)))
!->out       pn%source_p_e  = po%source_p_e
!->out    END IF
!->out    IF (ALLOCATED(po%source_m_e)) THEN
!->out       ALLOCATE(pn%source_m_e(SIZE(po%source_m_e,1)))
!->out       pn%source_m_e  = po%source_m_e
!->out    END IF
    IF (ALLOCATED(po%source_p)) THEN                                       !<-in
       ALLOCATE(pn%source_p(SIZE(po%source_p,1),SIZE(po%source_p,2)))      !<-in
       pn%source_p  = po%source_p                                          !<-in
    END IF                                                                 !<-in
    IF (ALLOCATED(po%source_m)) THEN                                       !<-in
       ALLOCATE(pn%source_m(SIZE(po%source_m,1),SIZE(po%source_m,2)))      !<-in
       pn%source_m  = po%source_m                                          !<-in
    END IF                                                                 !<-in
    
    IF (ALLOCATED(po%flux_p)) THEN
!->out       ALLOCATE(pn%flux_p(SIZE(po%flux_p,1)))
       ALLOCATE(pn%flux_p(SIZE(po%flux_p,1),SIZE(po%flux_p,2)))            !<-in
       pn%flux_p     = po%flux_p
    END IF
    IF (ALLOCATED(po%flux_m)) THEN
!->out       ALLOCATE(pn%flux_m(SIZE(po%flux_m,1)))
       ALLOCATE(pn%flux_m(SIZE(po%flux_m,1),SIZE(po%flux_m,2)))            !<-in
       pn%flux_m     = po%flux_m
    END IF
!->out    IF (ALLOCATED(po%curr_p)) THEN
!->out       ALLOCATE(pn%curr_p(SIZE(po%curr_p,1)))
!->out       pn%curr_p     = po%curr_p
!->out    END IF
!->out    IF (ALLOCATED(po%curr_m)) THEN
!->out       ALLOCATE(pn%curr_m(SIZE(po%curr_m,1)))
!->out       pn%curr_m     = po%curr_m
!->out    END IF

    IF (ALLOCATED(po%eta_l)) THEN
       ALLOCATE(pn%eta_l(LBOUND(po%eta_l,1):UBOUND(po%eta_l,1)))
       pn%eta_l      = po%eta_l
    END IF
    IF (ALLOCATED(po%eta_r)) THEN
       ALLOCATE(pn%eta_r(LBOUND(po%eta_r,1):UBOUND(po%eta_r,1)))
       pn%eta_r      = po%eta_r
    END IF


  END SUBROUTINE assign_propagator_cont_qe

  ! ---------------------------------------------------------------------------
  SUBROUTINE diag_propagator_cont()
    !
    USE magnetics_mod
    USE device_mod
    ! physical diagnostic
    ! all quantities have to be taken from prop_a
    INTEGER :: uw
    INTEGER :: i,proptag
    LOGICAL :: opened
    CHARACTER(len=10) :: cname,ctag_s,ctag_e
    CHARACTER(len=22) :: cadd,tags

    TYPE(fieldpropagator_struct), POINTER :: plotpropagator
    
    ! find free unit
    uw = 100
    DO
       INQUIRE(unit=uw,opened=opened)
       IF(.NOT. opened) EXIT
       uw = uw + 100
    END DO

    WRITE(ctag_s,*) prop_a%fieldpropagator_tag_s
    WRITE(ctag_e,*) prop_a%fieldpropagator_tag_e
    IF (prop_a%fieldpropagator_tag_s .EQ. prop_a%fieldpropagator_tag_e) THEN
       cadd = '_'//TRIM(ADJUSTL(ctag_s))//'.dat'
    ELSE
       cadd = '_'//TRIM(ADJUSTL(ctag_s))//'_'//TRIM(ADJUSTL(ctag_e))//'.dat'
    END IF
    tags = cadd

    IF (prop_overwrite .EQ. 1) THEN
       cadd = '.dat'
    END IF
    
    cname = 'amat_p_p'
    OPEN(unit=uw,file=TRIM(ADJUSTL(cname))//TRIM(ADJUSTL(cadd)))
    DO i = 1,SIZE(prop_a%p%amat_p_p,1)
       WRITE(uw,*) prop_a%p%amat_p_p(i,:)
    END DO
    CLOSE(uw)

    cname = 'amat_m_m'
    OPEN(unit=uw,file=TRIM(ADJUSTL(cname))//TRIM(ADJUSTL(cadd)))
    DO i = 1,SIZE(prop_a%p%amat_m_m,1)
       WRITE(uw,*) prop_a%p%amat_m_m(i,:)
    END DO
    CLOSE(uw)

    cname = 'amat_p_m'
    OPEN(unit=uw,file=TRIM(ADJUSTL(cname))//TRIM(ADJUSTL(cadd)))
    DO i = 1,SIZE(prop_a%p%amat_p_m,1)
       WRITE(uw,*) prop_a%p%amat_p_m(i,:)
    END DO
    CLOSE(uw)

    cname = 'amat_m_p'
    OPEN(unit=uw,file=TRIM(ADJUSTL(cname))//TRIM(ADJUSTL(cadd)))
    DO i = 1,SIZE(prop_a%p%amat_m_p,1)
       WRITE(uw,*) prop_a%p%amat_m_p(i,:)
    END DO
    CLOSE(uw)

    cname = 'source_p_g'
    OPEN(unit=uw,file=TRIM(ADJUSTL(cname))//TRIM(ADJUSTL(cadd)))
!->out    WRITE(uw,*) prop_a%p%source_p_g
    WRITE(uw,*) prop_a%p%source_p(:,1)                                     !<-in
    CLOSE(uw)
    cname = 'source_p_e'
    OPEN(unit=uw,file=TRIM(ADJUSTL(cname))//TRIM(ADJUSTL(cadd)))
!->out    WRITE(uw,*) prop_a%p%source_p_e
    WRITE(uw,*) prop_a%p%source_p(:,2)                                     !<-in
    CLOSE(uw)
    cname = 'source_m_g'
    OPEN(unit=uw,file=TRIM(ADJUSTL(cname))//TRIM(ADJUSTL(cadd)))
!->out    WRITE(uw,*) prop_a%p%source_m_g
    WRITE(uw,*) prop_a%p%source_m(:,1)                                     !<-in
    CLOSE(uw)
    cname = 'source_m_e'
    OPEN(unit=uw,file=TRIM(ADJUSTL(cname))//TRIM(ADJUSTL(cadd)))
!->out    WRITE(uw,*) prop_a%p%source_m_e
    WRITE(uw,*) prop_a%p%source_m(:,2)                                     !<-in
    CLOSE(uw)

    cname = 'flux_p'
    OPEN(unit=uw,file=TRIM(ADJUSTL(cname))//TRIM(ADJUSTL(cadd)))
!->out    WRITE(uw,*) prop_a%p%flux_p
    WRITE(uw,*) prop_a%p%flux_p(1,:)                                       !<-in
    CLOSE(uw)    
    cname = 'flux_m'
    OPEN(unit=uw,file=TRIM(ADJUSTL(cname))//TRIM(ADJUSTL(cadd)))
!->out    WRITE(uw,*) prop_a%p%flux_m
    WRITE(uw,*) prop_a%p%flux_m(1,:)                                       !<-in
    CLOSE(uw)
    cname = 'curr_p'
    OPEN(unit=uw,file=TRIM(ADJUSTL(cname))//TRIM(ADJUSTL(cadd)))
!->out    WRITE(uw,*) prop_a%p%curr_p
    WRITE(uw,*) prop_a%p%flux_p(2,:)                                       !<-in
    CLOSE(uw)    
    cname = 'curr_m'
    OPEN(unit=uw,file=TRIM(ADJUSTL(cname))//TRIM(ADJUSTL(cadd)))
!->out    WRITE(uw,*) prop_a%p%curr_m
    WRITE(uw,*) prop_a%p%flux_m(2,:)                                       !<-in
    CLOSE(uw)
    
    cname = 'eta'
    OPEN(unit=uw,file=TRIM(ADJUSTL(cname))//TRIM(ADJUSTL(cadd)))
    IF (prop_a%bin_split_mode .EQ. 1) THEN
       WRITE(uw,*) prop_a%p%eta_boundary_l,    &
            prop_a%p%eta_boundary_r,           &
            prop_a%eta_bs_l%x(1:)
    ELSE
       WRITE(uw,*) prop_a%p%eta_boundary_l,    &
            prop_a%p%eta_boundary_r,           &
            fieldpropagator%ch_act%eta(1:)
    END IF
    CLOSE(uw)
    
    ! summary
    cname = 'summary'
    OPEN(unit=uw,file=TRIM(ADJUSTL(cname))//TRIM(ADJUSTL(cadd)))
    WRITE(uw,*) 'bin_split_mode  ',prop_a%bin_split_mode
    WRITE(uw,*) 'tag_start       ',prop_a%fieldpropagator_tag_s
    WRITE(uw,*) 'tag_end         ',prop_a%fieldpropagator_tag_e
    WRITE(uw,*) 'tag_period_start',prop_a%fieldperiod_tag_s
    WRITE(uw,*) 'tag_period_end  ',prop_a%fieldperiod_tag_e
    WRITE(uw,*) ' '
    WRITE(uw,*) 'npart           ',prop_a%p%npart
    WRITE(uw,*) 'npass_l         ',prop_a%p%npass_l
    WRITE(uw,*) 'npass_r         ',prop_a%p%npass_r
!->out    WRITE(uw,*) 'npart_halfband  ',prop_a%p%npart_halfband
    WRITE(uw,*) 'nvelocity        ',prop_a%p%nvelocity                       !<-in
!->out    WRITE(uw,*) 'qflux_g         ',prop_a%p%qflux_g
!->out    WRITE(uw,*) 'qflux_e         ',prop_a%p%qflux_e
!->out    WRITE(uw,*) 'qcurr_g         ',prop_a%p%qcurr_g
!->out    WRITE(uw,*) 'qcurr_e         ',prop_a%p%qcurr_e
    WRITE(uw,*) 'qflux_g         ',prop_a%p%qflux(1,1)                     !<-in
    WRITE(uw,*) 'qflux_e         ',prop_a%p%qflux(1,2)                     !<-in
    WRITE(uw,*) 'qcurr_g         ',prop_a%p%qflux(2,1)                     !<-in
    WRITE(uw,*) 'qcurr_e         ',prop_a%p%qflux(2,2)                     !<-in
    WRITE(uw,*) 'eta_boundary_l  ',prop_a%p%eta_boundary_l
    WRITE(uw,*) 'eta_boundary_r  ',prop_a%p%eta_boundary_r
    CLOSE(uw)
    
    PRINT *, 'Physical Output written on files *'//TRIM(ADJUSTL(cadd))
    PRINT *, 'Involved propagator tags:         '  &
         //TRIM(ADJUSTL(ctag_s))//' '//TRIM(ADJUSTL(ctag_e))

    ! output only for original propagators (not joined)
    IF (prop_a%fieldpropagator_tag_s .EQ. prop_a%fieldpropagator_tag_e) THEN
       proptag = prop_a%fieldpropagator_tag_s
       CALL info_magnetics(fieldpropagator)
       cname = 'propm'
       OPEN(unit=uw,file=TRIM(ADJUSTL(cname))//TRIM(ADJUSTL(cadd)))

       plotpropagator => fieldpropagator
       CALL plot_magnetics(plotpropagator,proptag,proptag,  &
            TRIM(ADJUSTL(cname))//TRIM(ADJUSTL(cadd)) &
            )
        
       cname = 'etam'
       OPEN(unit=uw,status='replace',                       &
            file=TRIM(ADJUSTL(cname))//TRIM(ADJUSTL(cadd))  &
            )
       DO i = 1,UBOUND(fieldpropagator%ch_act%eta,1)
          WRITE(uw,'(3i6,1000(1x,e15.8))')              &
               i,fieldpropagator%phi_eta_ind(i,1),      &
               fieldpropagator%phi_eta_ind(i,2),        &
               fieldpropagator%ch_act%eta(i),           &
               1.0_dp - fieldpropagator%ch_act%eta(i) * & 
               fieldpropagator%mdata%bhat(              &
               fieldpropagator%phi_eta_ind(i,1)         &
               ),                                       &
               1.0_dp - fieldpropagator%ch_act%eta(i) * &
               fieldpropagator%mdata%bhat(              &
               fieldpropagator%phi_eta_ind(i,2)         &
               )
       END DO
       CLOSE(unit=uw)
       
    END IF
    
    !PAUSE

  END SUBROUTINE diag_propagator_cont

  SUBROUTINE diag_propagator_res(iend)
    ! writes out eps_eff and related stuff
    ! ATTENTION hxeta does not exist any more

    USE device_mod
    USE collisionality_mod, ONLY : collpar, conl_over_mfp, &
         isw_lorentz, isw_integral, isw_energy, isw_axisymm, y_axi_averages
    USE rkstep_mod, ONLY : asource,anumm,ailmm,lag,leg
    USE collop, ONLY : z_eff
    USE mag_interface_mod, ONLY : magnetic_device,mag_magfield,&
         mag_coordinates,boozer_s,boozer_theta_beg,boozer_phi_beg
    !! Modification by Andreas F. Martitsch (14.07.2015)
    ! Extra output/input for NTV computations
    USE ntv_mod, ONLY : isw_qflux_NA, write_ntv_output
    !! End Modification by Andreas F. Martitsch (14.07.2015)

    INTEGER, INTENT(in) :: iend

    LOGICAL :: opened
    INTEGER :: uw

    REAL(kind=dp) :: g_bs                                                 !<-GBS

    REAL(kind=dp) :: transport_factor
    REAL(kind=dp) :: qflux_g,qcurr_g
    REAL(kind=dp) :: qflux_e,qcurr_e
    REAL(kind=dp) :: dmono_over_dplateau,epseff3_2
    REAL(kind=dp) :: alambda_b,alambda_bb
    REAL(kind=dp) :: gamma_E
    REAL(kind=dp) :: aiota_loc,rt0
    REAL(kind=dp) :: phi

    REAL(kind=dp), ALLOCATABLE :: y(:)

    ! Declarations for final output
    REAL(kind=dp) :: gamma_fco(3,3)
    REAL(kind=dp) :: gamma_out(3,3)
    REAL(kind=dp) :: beta_out(3)
    REAL(kind=dp) :: avnabpsi,avbhat2,dl1obhat

    INTEGER, PARAMETER, DIMENSION(3) :: ind_map = (/1,3,2/)
    INTEGER :: i, i_p, j, j_p
    INTEGER :: full_version

    ! taken from Sergei
!->out    qflux_g = prop_a%p%qflux_g
!->out    qcurr_g = prop_a%p%qcurr_g
!->out    qflux_e = prop_a%p%qflux_e
!->out    qcurr_e = prop_a%p%qcurr_e
    qflux_g = prop_a%p%qflux(1,1)                                          !<-in
    qcurr_g = prop_a%p%qflux(2,1)                                          !<-in
    qflux_e = prop_a%p%qflux(1,2)                                          !<-in
    qcurr_e = prop_a%p%qflux(2,2)                                          !<-in
    IF ( (magnetic_device .EQ. 0 .AND. isw_axisymm .EQ. 1) .OR. mag_magfield .EQ. 0 ) THEN
       ALLOCATE(y(SIZE(y_axi_averages,1)))
       y = y_axi_averages
    ELSE
       ALLOCATE(y(SIZE(prop_a%y,1)))
       y = prop_a%y
    END IF
    aiota_loc = surface%aiota
    rt0 = device%r0
    phi = prop_a%phi_r

!!$    PRINT *, 'y ',y
!!$    PRINT *, 'aiota_loc ',aiota_loc 
!!$    PRINT *, 'rt0 ',rt0
!!$    PRINT *, 'phi ',phi
!!$    PRINT *, 'qflux_g ',qflux_g
    transport_factor = qflux_g*y(6)*(y(14)/(y(7)*y(13)))**2
!!$    PRINT *, 'transport_factor ',transport_factor
!!$    PAUSE
    dmono_over_dplateau=-2.d0*SQRT(2.d0)/pi*rt0*aiota_loc*transport_factor

    epseff3_2=-(9.d0*pi/(16.d0*SQRT(2.d0)))*collpar*rt0**2     &
         *transport_factor 

    alambda_b=-0.75d0 * qcurr_g *y(6)/(y(7)*y(9))
    alambda_bb=alambda_b * y(9) / y(6)
    g_bs=2.d0*qcurr_g/(y(7)*(collpar*qcurr_e/y(9)-8.d0/3.d0))             !<-GBS
    ! gamma_E from the Spitzer-Haerm paper 
    ! to get the values there, the result has to be multiplied by Zeff
    ! sigma / sigma_lorentz(zeff=1)
!    gamma_E = (3.d0*SQRT(pi)*qcurr_e*collpar/(32.d0*y(9)))  !***change19.09.07
    gamma_E = (3.d0*pi*qcurr_e*collpar/(32.d0*y(9)))
    !PRINT *, 'I am in diag_propagator_result'
    !PRINT *, 'iota = ',aiota_loc, 'qflux_g = ',qflux_g,' qcurr_g = ',qcurr_g
    
    ! find free unit
    uw = 100
    DO
       INQUIRE(unit=uw,opened=opened)
       IF(.NOT. opened) EXIT
       uw = uw + 100
    END DO
    OPEN(uw,file='evolve.dat',position='append')
!!$    WRITE (uw,'(1000(1x,e12.5))')                                   &
!!$         REAL(phi),REAL(y(1:2)),REAL(aiota_loc),                    &
!!$         REAL(dmono_over_dplateau),REAL(epseff3_2),REAL(alambda_b), &
!!$         REAL(qflux_g),REAL(qflux_e),REAL(qcurr_g),REAL(qcurr_e)    &
!!$         ,REAL(alambda_bb),REAL(3.d0*SQRT(pi)*qcurr_e*collpar/(32.d0*y(9)))*1.d0 &
!!$         ,REAL(g_bs)    &                                              !<-GBS
!!$         ,REAL(device%r0),REAL(surface%bmod0)
    WRITE (uw,'(1000(1x,e18.5))')                                   &
         (phi),(y(1:2)),(aiota_loc),                    &
         (dmono_over_dplateau),(epseff3_2),(alambda_b), &
         (qflux_g),(qflux_e),(qcurr_g),(qcurr_e),    &
         (alambda_bb),(gamma_E), &
         (g_bs),    &                                              !<-GBS
         (device%r0),(surface%bmod0)  !, &
         !y(6),y(7),y(9),y(13),y(14)
    ! WRITE (uw,*)                                   &
    !     y
    CLOSE(uw)

    ! Final output
    IF ( iend .EQ. 1) THEN
       full_version = 2
       avnabpsi = y(7) / y(6)
       avbhat2 = y(9) / y(6)
       dl1obhat = y(6)
!       beta_out = (/ 1.0_dp/avnabpsi, 1.0_dp/avnabpsi, 1.0_dp /) ***change: 19.09.07
       beta_out = (/ y(14)/y(13)/avnabpsi, y(14)/y(13)/avnabpsi, y(13)/y(14) /)
    
       DO i = 1,3
          i_p = ind_map(i)
          DO j = 1,3
             j_p = ind_map(j)
!             gamma_fco(i,j) = prop_a%p%qflux(i_p,j_p) / y(9)!***change: 19.09.07
             gamma_fco(i,j) = - prop_a%p%qflux(i_p,j_p) / y(6)
             gamma_out(i,j) = gamma_fco(i,j) * beta_out(i) * beta_out(j)
          END DO
       END DO
       !
       OPEN(uw,file='fulltransp.dat',status='replace')
       WRITE (uw,'(6(1x,i4),1000(1x,e18.5))')  &
            full_version, &
            isw_lorentz, isw_integral, isw_energy, lag, leg, &
            conl_over_mfp, collpar, z_eff, &
            avnabpsi, avbhat2, dl1obhat, &
            gamma_out
       CLOSE(uw)
       !
       OPEN(uw,file='efinal.dat',status='replace')
       WRITE (uw,'(1000(1x,e18.5))')                                   &
            (phi),(y(1:2)),(aiota_loc),                    &
            (dmono_over_dplateau),(epseff3_2),(alambda_b), &
            (qflux_g),(qflux_e),(qcurr_g),(qcurr_e),    &
            (alambda_bb),(gamma_E), &
            (g_bs),    &
            (device%r0),(surface%bmod0), &
            y(6),y(7),y(9),y(13),y(14)
       CLOSE(uw)
       !
       OPEN(uw,file='sigma_alex.dat')
!       write(uw,*) 1.5d0*sqrt(pi)**3*collpar*(device%r0), &
       WRITE(uw,*) 0.75d0*SQRT(pi)**3*collpar*(device%r0)/aiota_loc, &
                   -3.d0*pi/32.d0*collpar*gamma_out(3,3)
       CLOSE(uw)
       !
       !! Modification by Andreas F. Martitsch (14.07.2015)
       !! Extra output for NTV computations
       CALL write_ntv_output(isw_qflux_NA,prop_a%p%qflux,ind_map,&
            beta_out,y,aiota_loc,rt0,avnabpsi)
       !! End Modification by Andreas F. Martitsch (14.07.2015)
!
    END IF
!    
    DEALLOCATE(y)
  END SUBROUTINE diag_propagator_res
  ! ---------------------------------------------------------------------------
  SUBROUTINE diag_propagator_dis

    CALL plot_distrf_interface

  END SUBROUTINE diag_propagator_dis
  ! ---------------------------------------------------------------------------
  SUBROUTINE propagator_solver_loc(iend,iendperiod,bin_split_mode,            &
       eta_ori,ierr_solv,ierr_join                                            &
       )
!  SUBROUTINE propagator_solver_loc(iend,iendperiod,bin_split_mode,            &
!       eta_ori,eta_bs,ierr_solv,ierr_join                                     &
!       )
    !
    USE device_mod
    USE flint_mod , ONLY : phi_split_mode,phi_place_mode,  &
         phi_split_min,hphi_mult,max_solver_try
    USE collisionality_mod, ONLY : isw_axisymm
    USE mag_interface_mod, ONLY : magnetic_device,mag_magfield
             

    ! parameter list
    INTEGER,                      INTENT(in)  :: iend
    INTEGER,                      INTENT(in)  :: iendperiod
    INTEGER,                      INTENT(in)  :: bin_split_mode
    REAL(kind=dp), DIMENSION(0:), INTENT(in)  :: eta_ori
    ! TYPE(binarysplit),            INTENT(in)  :: eta_bs
    INTEGER,                      INTENT(out) :: ierr_solv
    INTEGER,                      INTENT(out) :: ierr_join
    ! 
    ! local quantity

    TYPE(binarysplit)                           :: eta_bs

    INTEGER                                     :: k,sy
    INTEGER                                     :: prop_npart
    INTEGER                                     :: i_joined
    INTEGER                                     :: count_solv
    INTEGER                                     :: phi_split_mode_ori   !<-in Winny
    REAL(kind=dp)                               :: mult_solv

    !
    ! initialize
    ierr_solv = 0
    ierr_join = 0
    sw_last_prop = iend

    fieldripple => fieldpropagator%ch_act
    eta_bs = fieldripple%eta_bs

    IF (bin_split_mode .EQ. 0) THEN
       prop_npart = UBOUND(eta_ori,1)
    ELSE
       prop_npart = UBOUND(eta_bs%x_split,1)
    END IF

    IF (prop_write .NE. 0) THEN
       prop_binary = 0
    END IF
    !
    IF (prop_timing .EQ. 1) CALL CPU_TIME(time_tot_o)
    ! Counter
    prop_count_call = prop_count_call + 1

       ! First call
       sw_first_prop = 0
       IF (prop_count_call .EQ. 1) THEN
          sw_first_prop = 1
          IF (prop_timing .EQ. 1) THEN
             stime_co  = 0.0_dp
             stime_jp  = 0.0_dp
             stime_ja  = 0.0_dp
             stime_jf  = 0.0_dp
             stime_so  = 0.0_dp
          END IF
          IF (prop_diagnostic .GE. 1) THEN
             PRINT *, '               construct propagator'
          END IF
          ! now construct the new propagator
          IF (prop_timing .EQ. 1) CALL CPU_TIME(time_o)
          CALL construct_propagator
          prop_c%nr_joined = -1
          IF (prop_timing .EQ. 1) THEN
             CALL CPU_TIME(time_co)
             stime_co = stime_co + time_co - time_o
          END IF
          prop_c%nr_joined = -1
       END IF

    IF ( (magnetic_device .EQ. 0 .AND. isw_axisymm .EQ. 1) .OR. mag_magfield .EQ. 0 ) THEN ! tokamak case

       ! ---------------------------------------------------------------------------
       ! call to ripple_solver for a propagator
       prop_c%bin_split_mode = bin_split_mode
       prop_c%p%npart = prop_npart
       prop_c%fieldpropagator_tag_s = fieldpropagator%tag
       prop_c%fieldpropagator_tag_e = fieldpropagator%tag
       prop_c%fieldperiod_tag_s     = fieldpropagator%parent%tag
       prop_c%fieldperiod_tag_e     = fieldpropagator%parent%tag
       prop_c%phi_l = fieldpropagator%phi_l
       prop_c%phi_r = fieldpropagator%phi_r

       ! put in eta_information
       IF (bin_split_mode .EQ. 1) THEN
          prop_c%eta_bs_l = eta_bs
          prop_c%eta_bs_r = eta_bs
       ENDIF
       IF (ALLOCATED(prop_c%p%eta_l)) DEALLOCATE(prop_c%p%eta_l)
       ALLOCATE(prop_c%p%eta_l( &
            LBOUND(fieldpropagator%ch_act%eta,1):UBOUND(fieldpropagator%ch_act%eta,1) &
            ))
       prop_c%p%eta_l  = fieldpropagator%ch_act%eta
       IF (ALLOCATED(prop_c%p%eta_r)) DEALLOCATE(prop_c%p%eta_r)
       ALLOCATE(prop_c%p%eta_r( &
            LBOUND(fieldpropagator%ch_act%eta,1):UBOUND(fieldpropagator%ch_act%eta,1) &
            ))
       prop_c%p%eta_r  = fieldpropagator%ch_act%eta

       IF (prop_diagnostic .GE. 2) THEN
          PRINT *, '               in ripple_solver ',prop_count_call
          PRINT *, '         ibegperiod, iendperiod ',prop_ibegperiod,iendperiod
       END IF
       IF (prop_timing .EQ. 1) CALL CPU_TIME(time_o)

       ! try the ripple_solver
       count_solv = 0
       mult_solv  = hphi_mult
       reduce_hphi_axisym: DO
          count_solv = count_solv + 1
          ierr_solv = 0
          CALL ripple_solver_interface(ierr_solv)
          IF (ierr_solv .EQ. 0) EXIT reduce_hphi_axisym
          IF (ierr_solv .NE. 0 .AND. ierr_solv .NE. 3) THEN
             PRINT *, 'Error in ripple_solver: ',ierr_solv
             PRINT *, ' I give up'
             STOP
          END IF

          IF (ierr_solv .EQ. 3 .AND. count_solv .LT. max_solver_try) THEN
             mult_solv = mult_solv * 0.5_dp
             phi_split_mode_ori = phi_split_mode                  !<-in Winny
             phi_split_mode = 3                                   !<-in Winny
             CALL modify_propagator(phi_split_mode,phi_place_mode,phi_split_min, &
                  UBOUND(prop_c%p%eta_l,1),prop_c%p%eta_l,mult_solv,count_solv)
             phi_split_mode = phi_split_mode_ori                  !<-in Winny
             PRINT *, 'Error in ripple_solver: ',ierr_solv
             PRINT *, ' I try it again ',count_solv+1
!!$          IF (count_solv .GT. 1) THEN
!!$             PRINT *, 'PAUSE - MODE'
!!$             PAUSE
!!$          END IF
          ELSE IF (ierr_solv .EQ. 3 .AND. count_solv .GE. max_solver_try) THEN
             PRINT *, 'Error in ripple_solver: ',ierr_solv
             PRINT *, ' I give up'
             STOP
          END IF
       END DO reduce_hphi_axisym
       ! write the end value for y of the field propagator into prop_c
       ! this is for computing physical output
       !sy = SIZE(fieldpropagator%mdata%yend,1)
       sy = SIZE(fieldpropagator%parent%mdata%yend,1) ! WINNY YEND
       IF (ALLOCATED(prop_c%y)) DEALLOCATE(prop_c%y)
       ALLOCATE(prop_c%y(sy))
       !prop_c%y = fieldpropagator%mdata%yend
       prop_c%y = fieldpropagator%parent%mdata%yend ! WINNY YEND
       ! write eta at boundaries - replaced by modified stuff
       !prop_c%p%eta_boundary_l = 1.0_dp / fieldpropagator%b_l
       !prop_c%p%eta_boundary_r = 1.0_dp / fieldpropagator%b_r
       prop_c%p%eta_boundary_l = eta_modboundary_l
       prop_c%p%eta_boundary_r = eta_modboundary_r
       ! link actual to current (mainly for results and diagnostic)
       prop_a => prop_c
       ! writing of propagators
       !IF (prop_write .EQ. 2) THEN
       !   CALL write_propagator_content(prop_a,3)
       !   IF (prop_first_tag .EQ. 0) prop_first_tag = fieldpropagator%tag
       !   prop_last_tag = fieldpropagator%tag
       !END IF
       ! ---------------------------------------------------------------------------
       CALL diag_propagator_result(iend)
       !! Modification by Andreas F. Martitsch (16.07.2015)
       ! not needed for NTV computations
       ! (furthermore, several uninitialised value(s) are used - valgrind)
       !CALL diag_propagator_distrf
       !! End Modification by Andreas F. Martitsch (16.07.2015)
    ELSE ! this is the general case
   
       IF (prop_reconstruct .EQ. 0) THEN
          ! Begin of period (no joining of all propagators
          IF (prop_ibegperiod .EQ. 1 .AND. prop_write .NE. 2) THEN
             IF (prop_c%nr_joined .EQ. -1 .AND. .NOT. ASSOCIATED(prop_c%next)) THEN
                IF (prop_diagnostic .GE. 1) THEN
                   PRINT *, '               construct propagator'
                END IF
                IF (prop_timing .EQ. 1) CALL CPU_TIME(time_o)
                CALL construct_propagator
                IF (prop_timing .EQ. 1) THEN
                   CALL CPU_TIME(time_co)
                   stime_co = stime_co + time_co - time_o
                END IF
                prop_c%nr_joined = -1
                prop_c => prop_c%prev
             END IF
             prop_a => prop_c
          END IF
       END IF

       IF (prop_reconstruct .EQ. 2) THEN
          CALL read_prop_recon_content(fieldpropagator%tag)
       END IF
       ! ---------------------------------------------------------------------------
       ! call to ripple_solver for a propagator
       prop_c%bin_split_mode = bin_split_mode
       prop_c%p%npart = prop_npart
       prop_c%fieldpropagator_tag_s = fieldpropagator%tag
       prop_c%fieldpropagator_tag_e = fieldpropagator%tag
       prop_c%fieldperiod_tag_s     = fieldpropagator%parent%tag
       prop_c%fieldperiod_tag_e     = fieldpropagator%parent%tag
       prop_c%phi_l = fieldpropagator%phi_l
       prop_c%phi_r = fieldpropagator%phi_r

       ! put in eta_information
       IF (bin_split_mode .EQ. 1) THEN
          prop_c%eta_bs_l = eta_bs
          prop_c%eta_bs_r = eta_bs
       ENDIF
       IF (ALLOCATED(prop_c%p%eta_l)) DEALLOCATE(prop_c%p%eta_l)
       ALLOCATE(prop_c%p%eta_l( &
            LBOUND(fieldpropagator%ch_act%eta,1):UBOUND(fieldpropagator%ch_act%eta,1) &
            ))
       prop_c%p%eta_l  = fieldpropagator%ch_act%eta
       IF (ALLOCATED(prop_c%p%eta_r)) DEALLOCATE(prop_c%p%eta_r)
       ALLOCATE(prop_c%p%eta_r( &
            LBOUND(fieldpropagator%ch_act%eta,1):UBOUND(fieldpropagator%ch_act%eta,1) &
            ))
       prop_c%p%eta_r  = fieldpropagator%ch_act%eta

       IF (prop_diagnostic .GE. 2) THEN
          PRINT *, '               in ripple_solver ',prop_count_call
          PRINT *, '         ibegperiod, iendperiod ',prop_ibegperiod,iendperiod
       END IF
       IF (prop_timing .EQ. 1) CALL CPU_TIME(time_o)

       ! try the ripple_solver
       count_solv = 0
       mult_solv  = hphi_mult
       reduce_hphi: DO
          count_solv = count_solv + 1
          ierr_solv = 0
          CALL ripple_solver_interface(ierr_solv)
          IF (ierr_solv .EQ. 0) EXIT reduce_hphi
          IF (ierr_solv .NE. 0 .AND. ierr_solv .NE. 3) THEN
             PRINT *, 'Error in ripple_solver: ',ierr_solv
             PRINT *, ' I give up'
             STOP
          END IF

          IF (ierr_solv .EQ. 3 .AND. count_solv .LT. max_solver_try) THEN
             mult_solv = mult_solv * 0.5_dp
             phi_split_mode_ori = phi_split_mode                  !<-in Winny
             phi_split_mode = 3                                   !<-in Winny
             CALL modify_propagator(phi_split_mode,phi_place_mode,phi_split_min, &
                  UBOUND(prop_c%p%eta_l,1),prop_c%p%eta_l,mult_solv,count_solv)
             phi_split_mode = phi_split_mode_ori                  !<-in Winny
             PRINT *, 'Error in ripple_solver: ',ierr_solv
             PRINT *, ' I try it again ',count_solv+1
!!$          IF (count_solv .GT. 1) THEN
!!$             PRINT *, 'PAUSE - MODE'
!!$             PAUSE
!!$          END IF
          ELSE IF (ierr_solv .EQ. 3 .AND. count_solv .GE. max_solver_try) THEN
             PRINT *, 'Error in ripple_solver: ',ierr_solv
             PRINT *, ' I give up'
             STOP
          END IF
       END DO reduce_hphi

       IF (prop_reconstruct .EQ. 2) RETURN

       ! write the end value for y of the field propagator into prop_c
       ! this is for computing physical output
       !sy = SIZE(fieldpropagator%mdata%yend,1)
       sy = SIZE(fieldpropagator%parent%mdata%yend,1) ! WINNY YEND
       IF (ALLOCATED(prop_c%y)) DEALLOCATE(prop_c%y)
       ALLOCATE(prop_c%y(sy))
       !prop_c%y = fieldpropagator%mdata%yend
       prop_c%y = fieldpropagator%parent%mdata%yend ! WINNY YEND
       ! write eta at boundaries - replaced by modified stuff
       !prop_c%p%eta_boundary_l = 1.0_dp / fieldpropagator%b_l
       !prop_c%p%eta_boundary_r = 1.0_dp / fieldpropagator%b_r
       prop_c%p%eta_boundary_l = eta_modboundary_l
       prop_c%p%eta_boundary_r = eta_modboundary_r
       ! link actual to current (mainly for results and diagnostic)
       prop_a => prop_c
       ! writing of propagators
       IF (prop_write .EQ. 2) THEN
          CALL write_propagator_content(prop_a,3)
          IF (prop_first_tag .EQ. 0) prop_first_tag = fieldpropagator%tag
          prop_last_tag = fieldpropagator%tag
       END IF
       ! ---------------------------------------------------------------------------
       IF (prop_diagphys .EQ. 1) THEN
          CALL diag_propagator_content
       END IF
       ! ---------------------------------------------------------------------------
       IF (prop_timing .EQ. 1) THEN
          CALL CPU_TIME(time_so)
          stime_so = stime_so + time_so - time_o
       END IF

       ! make a new propagator or join
       i_joined = 0
       ! ---------------------------------------------------------------------------
       IF (prop_ibegperiod .EQ. 1 .AND. iendperiod .EQ. 1 .AND. prop_write .NE. 2) THEN 
          i_joined = 1
          prop_a => prop_c
          prop_ibegperiod = 1
          prop_c%nr_joined = 0
          prop_c => prop_c%next
          IF (.NOT. ASSOCIATED(prop_c%next)) THEN
             CALL construct_propagator
             prop_c%nr_joined = -3
             prop_c => prop_c%prev
          END IF
       ELSEIF (prop_count_call .EQ. 1 .AND. prop_write .EQ. 2) THEN
          CALL construct_propagator
          prop_c%nr_joined = -1
          ! prop_c => prop_c%prev
          prop_a => prop_c       
       ELSE
          IF (prop_ibegperiod .EQ. 0 .OR. iendperiod .EQ. 1 .OR. prop_write .EQ. 2) THEN 
             IF (prop_diagnostic .GE. 2) THEN
                PRINT *, '               join_ripples within period'
             END IF
             IF (prop_timing .EQ. 1) CALL CPU_TIME(time_o)
             CALL join_ripples_interface(ierr_join,'inter',0)
             i_joined = 1
             IF (prop_timing .EQ. 1) THEN
                CALL CPU_TIME(time_jp)
                stime_jp = stime_jp + time_jp - time_o
             END IF
             prop_a => prop_c%prev
             ! writing of propagators (joined within period)
             IF (prop_write .EQ. 2) THEN 
                CALL write_prop_bound_content(prop_c%prev,prop_c,3)
                CALL write_propagator_content(prop_c%prev,3,2) ! reduced for joined
             END IF
             IF (ALLOCATED(prop_c%prev%p%cmat)) DEALLOCATE(prop_c%prev%p%cmat)
             IF (ALLOCATED(prop_c%p%cmat)) DEALLOCATE(prop_c%p%cmat)

          END IF
          IF (prop_ibegperiod .EQ. 1 .AND. prop_write .NE. 2) THEN
             prop_c => prop_c%next
          END IF
          IF (prop_ibegperiod .EQ. 1) THEN
             prop_ibegperiod = 0
          END IF
          IF (iendperiod .EQ. 1) THEN
             ! writing of periods
             ! Winny: This is the place to write things to the file
             ! Things at this point are in prop_a => prop_c%prev
             ! cmat is computed when you join periods, so it can not be written here
             IF (prop_write .EQ. 1) THEN 
                CALL write_propagator_content(prop_a,1)
                IF (prop_first_tag .EQ. 0) prop_first_tag = fieldpropagator%parent%tag
                prop_last_tag = fieldpropagator%parent%tag
             END IF
             ! Winny
             prop_ibegperiod = 1
             prop_c%nr_joined = -1
             IF (prop_write .NE. 2) THEN
                prop_c%prev%nr_joined = 0
             END IF
          END IF
       END IF
       ! ---------------------------------------------------------------------------


       ! ---------------------------------------------------------------------------
       ! Check for advanced joining
       IF (iendperiod .EQ. 1 .AND. prop_write .NE. 2) THEN
          prop_c => prop_c%prev
          ! Binary joining
          IF (prop_binary .EQ. 1) THEN
             k = 1
             DO
                IF ( .NOT. ASSOCIATED(prop_c%prev) ) EXIT
                IF ( prop_c%nr_joined .EQ. prop_c%prev%nr_joined ) THEN
                   k = k + 1
                   IF (prop_diagnostic .GE. 1) THEN
                      PRINT *, '               join_ripples advanced binary',k
                   END IF
                   IF (prop_timing .EQ. 1) CALL CPU_TIME(time_o)
                   CALL join_ripples_interface(ierr_join)
                   IF (ALLOCATED(prop_c%prev%p%cmat)) DEALLOCATE(prop_c%prev%p%cmat)
                   IF (ALLOCATED(prop_c%p%cmat)) DEALLOCATE(prop_c%p%cmat)
                   i_joined = 1
                   IF (prop_timing .EQ. 1) THEN
                      CALL CPU_TIME(time_ja)
                      stime_ja = stime_ja + time_ja - time_o
                   END IF
                   prop_c%nr_joined = -1
                   prop_c => prop_c%prev
                   prop_c%nr_joined = k
                   prop_a => prop_c
                ELSE
                   EXIT
                END IF
             END DO
          ELSE
             IF ( ASSOCIATED(prop_c%prev) ) THEN
                k = 2
                IF (prop_diagnostic .GE. 1) THEN
                   PRINT *, '               join_ripples advanced period',k
                END IF
                IF (prop_timing .EQ. 1) CALL CPU_TIME(time_o)
                CALL join_ripples_interface(ierr_join,'inter',0)
                ! Winny 
                ! here the periods are joined, so cmat should be available
                ! but cmat is not stored, so this has to be handled
                ! PRINT *, 'Joining of periods'
                ! PRINT *,fieldpropagator%parent%prev%tag
                ! PRINT *,fieldpropagator%parent%tag
                IF (prop_write .EQ. 1) CALL write_prop_bound_content(prop_c%prev,prop_c,1)
                IF (ALLOCATED(prop_c%prev%p%cmat)) DEALLOCATE(prop_c%prev%p%cmat)
                IF (ALLOCATED(prop_c%p%cmat)) DEALLOCATE(prop_c%p%cmat)
                IF (prop_write .EQ. 1) CALL write_propagator_content(prop_c%prev,1,2)
                ! Winny
                i_joined = 1
                IF (prop_timing .EQ. 1) THEN
                   CALL CPU_TIME(time_ja)
                   stime_ja = stime_ja + time_ja - time_o
                END IF
                prop_c%nr_joined = -1
                prop_c => prop_c%prev
                prop_c%nr_joined = 1
                prop_a => prop_c
             END IF
          END IF
          prop_c => prop_c%next
       END IF

       ! Final joining
       IF (iend .EQ. 1) THEN
          k = 1
          DO 
             prop_c => prop_c%prev
             IF ( ASSOCIATED(prop_c%prev) ) THEN
                k = k + 1
                IF (prop_diagnostic .GE. 1) THEN
                   PRINT *, '               final joining ',k
                END IF
                IF (prop_timing .EQ. 1) CALL CPU_TIME(time_o)
                CALL join_ripples_interface(ierr_join)
                i_joined = 1 
                IF (prop_timing .EQ. 1) THEN
                   CALL CPU_TIME(time_jf)
                   stime_jf = stime_jf + time_jf - time_o
                END IF
             ELSE
                ! join ends
                IF (prop_join_ends .EQ. 1) THEN
                   i_joined = 1
                   CALL construct_propagator()
                   !prop_c = prop_c%prev
                   CALL assign_propagator_content(prop_c,prop_c%prev)
                   CALL join_ripples_interface(ierr_join,'final')
                   prop_c => prop_c%prev
                END IF
                !
                prop_a => prop_c
                IF (prop_write .EQ. 1) THEN
                   ! final joining
                   CALL write_propagator_content(prop_a,2)
                ELSEIF (prop_write .EQ. 2) THEN
                   ! final joining
                   CALL write_propagator_content(prop_a,4)
                END IF
                IF (prop_write .EQ. 1 .OR. prop_write .EQ. 2) THEN
                   ! taginfo
                   CALL unit_propagator
                   OPEN(unit=prop_unit,file=prop_ctaginfo,status='replace', &
                        form=prop_format,action='write')
                   WRITE(prop_unit,*) prop_write
                   WRITE(prop_unit,*) prop_first_tag
                   WRITE(prop_unit,*) prop_last_tag
                   CLOSE(unit=prop_unit)
                END IF
                EXIT
             END IF
          END DO
       END IF
       !
       IF (prop_timing .EQ. 1) THEN
          CALL CPU_TIME(time_tot)
          stime_tot = stime_tot + time_tot - time_tot_o
       END IF
       ! diagnostic
       IF (i_joined .EQ. 1 .AND. prop_diagphys .EQ. 1) THEN
          CALL diag_propagator_content 
       END IF
       ! physical output
       IF (i_joined .EQ. 1 .AND. iendperiod .EQ. 1) THEN
          CALL diag_propagator_result(iend)
       END IF

       IF (i_joined .EQ. 1 .AND. iend .EQ. 1) THEN
          CALL diag_propagator_distrf 
       END IF
    END IF

    IF (iend .EQ. 1 .AND. prop_timing .EQ. 1) THEN
       PRINT *, ' '
       PRINT *, 'Time for construction      :', stime_co 
       PRINT *, 'Time for solver            :', stime_so
       PRINT *, 'Time for joining in period :', stime_jp 
       PRINT *, 'Time for advanced joining  :', stime_ja 
       PRINT *, 'Time for final joining     :', stime_jf
       PRINT *, ' '
       PRINT *, 'Total Time for joining     :', stime_jp + stime_ja + stime_jf
       PRINT *, 'Total Time in module       :', stime_tot 
       PRINT *, ' '
    END IF

    CALL deconstruct_binarysplit(eta_bs)

    !
  END SUBROUTINE propagator_solver_loc
  
  ! ---------------------------------------------------------------------------
  SUBROUTINE ripple_solver_int(ierr)
    !
    ! Interface routine for subroutine ripple_solver
    ! Connects to the outside
    !
    ! Has to be changed if physical content of propagator changes!
    ! (see type declaration of propagator)
    !
    ! Call to external program ripple_solver
    !
    !! Modification by Andreas F. Martitsch (27.07.2015)
    ! Switch for ripple_solver version
    USE ntv_mod, ONLY : isw_ripple_solver
    !! End Modification by Andreas F. Martitsch (27.07.2015)
    !
    INTEGER :: ierr
    
    INTERFACE ripple_solver
       SUBROUTINE ripple_solver(                                 &
!->out            npass_l,npass_r,npart_halfband,                      &
            npass_l,npass_r,nvelocity,                            &         !<-in
            amat_plus_plus,amat_minus_minus,                     &
            amat_plus_minus,amat_minus_plus,                     &
!->out            source_p_g,source_m_g,source_p_e,source_m_e,         &
            source_p,source_m,                                   &         !<-in
            flux_p,flux_m,                                       &
!->out            curr_p,curr_m,                                       &
!->out            qflux_g,qflux_e,                                     &
!->out            qcurr_g,qcurr_e,                                     &
            qflux,                                               &         !<-in
            ierr                                                 &
            )
         INTEGER, PARAMETER                                        :: dp = KIND(1.0d0)

         INTEGER,                                    INTENT(out)   :: npass_l
         INTEGER,                                    INTENT(out)   :: npass_r
!->out         INTEGER,                                    INTENT(out)   :: npart_halfband
         INTEGER,   INTENT(out)   :: nvelocity                              !<-in
         REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: amat_plus_plus
         REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: amat_minus_minus
         REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: amat_plus_minus
         REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: amat_minus_plus
!->out         REAL(kind=dp), DIMENSION(:),   ALLOCATABLE, INTENT(inout) :: source_p_g
!->out         REAL(kind=dp), DIMENSION(:),   ALLOCATABLE, INTENT(inout) :: source_m_g
!->out         REAL(kind=dp), DIMENSION(:),   ALLOCATABLE, INTENT(inout) :: source_p_e
!->out         REAL(kind=dp), DIMENSION(:),   ALLOCATABLE, INTENT(inout) :: source_m_e
!->out         REAL(kind=dp), DIMENSION(:),   ALLOCATABLE, INTENT(inout) :: flux_p
!->out         REAL(kind=dp), DIMENSION(:),   ALLOCATABLE, INTENT(inout) :: flux_m
!->out         REAL(kind=dp), DIMENSION(:),   ALLOCATABLE, INTENT(inout) :: curr_p
!->out         REAL(kind=dp), DIMENSION(:),   ALLOCATABLE, INTENT(inout) :: curr_m
         REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(inout) ::         &
                                       source_p,source_m,flux_p,flux_m     !<-in
!->out         REAL(kind=dp),                              INTENT(out)   :: qflux_g
!->out         REAL(kind=dp),                              INTENT(out)   :: qflux_e
!->out         REAL(kind=dp),                              INTENT(out)   :: qcurr_g
!->out         REAL(kind=dp),                              INTENT(out)   :: qcurr_e
         REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(out) :: qflux  !<-in
         INTEGER,                                    INTENT(out)   :: ierr
       END SUBROUTINE ripple_solver
    END INTERFACE ripple_solver
    !! Modification by Andreas F. Martitsch (27.07.2015)
    ! Interface Arnoldi Solver Order 1
    INTERFACE ripple_solver_ArnoldiO1
       SUBROUTINE ripple_solver_ArnoldiO1(                       &
!->out            npass_l,npass_r,npart_halfband,                      &
            npass_l,npass_r,nvelocity,                            &         !<-in
            amat_plus_plus,amat_minus_minus,                     &
            amat_plus_minus,amat_minus_plus,                     &
!->out            source_p_g,source_m_g,source_p_e,source_m_e,         &
            source_p,source_m,                                   &         !<-in
            flux_p,flux_m,                                       &
!->out            curr_p,curr_m,                                       &
!->out            qflux_g,qflux_e,                                     &
!->out            qcurr_g,qcurr_e,                                     &
            qflux,                                               &         !<-in
            ierr                                                 &
            )
         INTEGER, PARAMETER                                        :: dp = KIND(1.0d0)

         INTEGER,                                    INTENT(out)   :: npass_l
         INTEGER,                                    INTENT(out)   :: npass_r
!->out         INTEGER,                                    INTENT(out)   :: npart_halfband
         INTEGER,   INTENT(out)   :: nvelocity                              !<-in
         REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: amat_plus_plus
         REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: amat_minus_minus
         REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: amat_plus_minus
         REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: amat_minus_plus
!->out         REAL(kind=dp), DIMENSION(:),   ALLOCATABLE, INTENT(inout) :: source_p_g
!->out         REAL(kind=dp), DIMENSION(:),   ALLOCATABLE, INTENT(inout) :: source_m_g
!->out         REAL(kind=dp), DIMENSION(:),   ALLOCATABLE, INTENT(inout) :: source_p_e
!->out         REAL(kind=dp), DIMENSION(:),   ALLOCATABLE, INTENT(inout) :: source_m_e
!->out         REAL(kind=dp), DIMENSION(:),   ALLOCATABLE, INTENT(inout) :: flux_p
!->out         REAL(kind=dp), DIMENSION(:),   ALLOCATABLE, INTENT(inout) :: flux_m
!->out         REAL(kind=dp), DIMENSION(:),   ALLOCATABLE, INTENT(inout) :: curr_p
!->out         REAL(kind=dp), DIMENSION(:),   ALLOCATABLE, INTENT(inout) :: curr_m
         REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(inout) ::         &
                                       source_p,source_m,flux_p,flux_m     !<-in
!->out         REAL(kind=dp),                              INTENT(out)   :: qflux_g
!->out         REAL(kind=dp),                              INTENT(out)   :: qflux_e
!->out         REAL(kind=dp),                              INTENT(out)   :: qcurr_g
!->out         REAL(kind=dp),                              INTENT(out)   :: qcurr_e
         REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(out) :: qflux  !<-in
         INTEGER,                                    INTENT(out)   :: ierr
       END SUBROUTINE ripple_solver_ArnoldiO1
    END INTERFACE ripple_solver_ArnoldiO1
    ! Interface Arnoldi Solver Order 2
    INTERFACE ripple_solver_ArnoldiO2
       SUBROUTINE ripple_solver_ArnoldiO2(                       &
!->out            npass_l,npass_r,npart_halfband,                      &
            npass_l,npass_r,nvelocity,                            &         !<-in
            amat_plus_plus,amat_minus_minus,                     &
            amat_plus_minus,amat_minus_plus,                     &
!->out            source_p_g,source_m_g,source_p_e,source_m_e,         &
            source_p,source_m,                                   &         !<-in
            flux_p,flux_m,                                       &
!->out            curr_p,curr_m,                                       &
!->out            qflux_g,qflux_e,                                     &
!->out            qcurr_g,qcurr_e,                                     &
            qflux,                                               &         !<-in
            ierr                                                 &
            )
         INTEGER, PARAMETER                                        :: dp = KIND(1.0d0)

         INTEGER,                                    INTENT(out)   :: npass_l
         INTEGER,                                    INTENT(out)   :: npass_r
!->out         INTEGER,                                    INTENT(out)   :: npart_halfband
         INTEGER,   INTENT(out)   :: nvelocity                              !<-in
         REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: amat_plus_plus
         REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: amat_minus_minus
         REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: amat_plus_minus
         REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: amat_minus_plus
!->out         REAL(kind=dp), DIMENSION(:),   ALLOCATABLE, INTENT(inout) :: source_p_g
!->out         REAL(kind=dp), DIMENSION(:),   ALLOCATABLE, INTENT(inout) :: source_m_g
!->out         REAL(kind=dp), DIMENSION(:),   ALLOCATABLE, INTENT(inout) :: source_p_e
!->out         REAL(kind=dp), DIMENSION(:),   ALLOCATABLE, INTENT(inout) :: source_m_e
!->out         REAL(kind=dp), DIMENSION(:),   ALLOCATABLE, INTENT(inout) :: flux_p
!->out         REAL(kind=dp), DIMENSION(:),   ALLOCATABLE, INTENT(inout) :: flux_m
!->out         REAL(kind=dp), DIMENSION(:),   ALLOCATABLE, INTENT(inout) :: curr_p
!->out         REAL(kind=dp), DIMENSION(:),   ALLOCATABLE, INTENT(inout) :: curr_m
         REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(inout) ::         &
                                       source_p,source_m,flux_p,flux_m     !<-in
!->out         REAL(kind=dp),                              INTENT(out)   :: qflux_g
!->out         REAL(kind=dp),                              INTENT(out)   :: qflux_e
!->out         REAL(kind=dp),                              INTENT(out)   :: qcurr_g
!->out         REAL(kind=dp),                              INTENT(out)   :: qcurr_e
         REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(out) :: qflux  !<-in
         INTEGER,                                    INTENT(out)   :: ierr
       END SUBROUTINE ripple_solver_ArnoldiO2
    END INTERFACE ripple_solver_ArnoldiO2
    !! End Modification by Andreas F. Martitsch (27.07.2015)

    !! Modification by Andreas F. Martitsch (27.07.2015)
    ! Select ripple_solver version
    IF (isw_ripple_solver .EQ. 1) THEN
       CALL ripple_solver(                                                  &
!->out         prop_c%p%npass_l,prop_c%p%npass_r,prop_c%p%npart_halfband,      &
            prop_c%p%npass_l,prop_c%p%npass_r,prop_c%p%nvelocity,            & !<-in
            prop_c%p%amat_p_p, prop_c%p%amat_m_m,                           &
            prop_c%p%amat_p_m, prop_c%p%amat_m_p,                           &
!->out         prop_c%p%source_p_g, prop_c%p%source_m_g,                       &
!->out         prop_c%p%source_p_e, prop_c%p%source_m_e,                       &
            prop_c%p%source_p, prop_c%p%source_m,                           & !<-in
            prop_c%p%flux_p, prop_c%p%flux_m,                               & 
!->out         prop_c%p%curr_p, prop_c%p%curr_m,                               &
!->out         prop_c%p%qflux_g, prop_c%p%qflux_e,                             &
!->out         prop_c%p%qcurr_g, prop_c%p%qcurr_e,                             &
            prop_c%p%qflux,                                                 & !<-in
            ierr                                                            &
            )
    ELSEIF (isw_ripple_solver .EQ. 2) THEN
       CALL ripple_solver_ArnoldiO1(                                        &
!->out         prop_c%p%npass_l,prop_c%p%npass_r,prop_c%p%npart_halfband,      &
            prop_c%p%npass_l,prop_c%p%npass_r,prop_c%p%nvelocity,            & !<-in
            prop_c%p%amat_p_p, prop_c%p%amat_m_m,                           &
            prop_c%p%amat_p_m, prop_c%p%amat_m_p,                           &
!->out         prop_c%p%source_p_g, prop_c%p%source_m_g,                       &
!->out         prop_c%p%source_p_e, prop_c%p%source_m_e,                       &
            prop_c%p%source_p, prop_c%p%source_m,                           & !<-in
            prop_c%p%flux_p, prop_c%p%flux_m,                               & 
!->out         prop_c%p%curr_p, prop_c%p%curr_m,                               &
!->out         prop_c%p%qflux_g, prop_c%p%qflux_e,                             &
!->out         prop_c%p%qcurr_g, prop_c%p%qcurr_e,                             &
            prop_c%p%qflux,                                                 & !<-in
            ierr                                                            &
            )
    ELSEIF (isw_ripple_solver .EQ. 3) THEN
       CALL ripple_solver_ArnoldiO2(                                        &
!->out         prop_c%p%npass_l,prop_c%p%npass_r,prop_c%p%npart_halfband,      &
            prop_c%p%npass_l,prop_c%p%npass_r,prop_c%p%nvelocity,            & !<-in
            prop_c%p%amat_p_p, prop_c%p%amat_m_m,                           &
            prop_c%p%amat_p_m, prop_c%p%amat_m_p,                           &
!->out         prop_c%p%source_p_g, prop_c%p%source_m_g,                       &
!->out         prop_c%p%source_p_e, prop_c%p%source_m_e,                       &
            prop_c%p%source_p, prop_c%p%source_m,                           & !<-in
            prop_c%p%flux_p, prop_c%p%flux_m,                               & 
!->out         prop_c%p%curr_p, prop_c%p%curr_m,                               &
!->out         prop_c%p%qflux_g, prop_c%p%qflux_e,                             &
!->out         prop_c%p%qcurr_g, prop_c%p%qcurr_e,                             &
            prop_c%p%qflux,                                                 & !<-in
            ierr                                                            &
            )
    ELSE
       STOP "Undefined version of ripple_solver selected (isw_ripple_solver)!"
    ENDIF
    !! End Modification by Andreas F. Martitsch (27.07.2015)

!!$    PRINT *, '----------------------------------------------------------'
!!$    PRINT *, 'SOLVER INTER : npass_l,npass_r:   ',  prop_c%p%npass_l,prop_c%p%npass_r
!->out    prop_c%p%npass_l = SIZE(prop_c%p%amat_p_m,1)
!->out    prop_c%p%npass_r = SIZE(prop_c%p%amat_m_p,1)
    prop_c%p%npass_l = SIZE(prop_c%p%amat_p_m,1)/(prop_c%p%nvelocity+1)     !<-in
    prop_c%p%npass_r = SIZE(prop_c%p%amat_m_p,1)/(prop_c%p%nvelocity+1)     !<-in
!!$    PRINT *, 'SOLVER INTER : npass_l,npass_r:   ',  prop_c%p%npass_l,prop_c%p%npass_r
!!$    PRINT *, 'SOLVER INTER : prop_c%p%amat_p_m: ',  &
!!$         SIZE(prop_c%p%amat_p_m,1),SIZE(prop_c%p%amat_p_m,2)
!!$    PRINT *, 'SOLVER INTER : prop_c%p%amat_m_p: ',  &
!!$         SIZE(prop_c%p%amat_m_p,1),SIZE(prop_c%p%amat_m_p,2)
!!$    PRINT *, '----------------------------------------------------------'
!!$    PAUSE
    


    !
  END SUBROUTINE ripple_solver_int

  SUBROUTINE plot_distrf_int

    !! Modification by Andreas F. Martitsch (17.07.2015)
    ! interface name ambiguous
    !INTERFACE ripple_solver
    INTERFACE plot_distrf
    !! End Modification by Andreas F. Martitsch (17.07.2015)
       SUBROUTINE plot_distrf(source_p,source_m,eta_l,eta_r,eta_boundary_l,eta_boundary_r)
         INTEGER, PARAMETER :: dp = KIND(1.0d0)
         
         REAL(kind=dp), DIMENSION(:,:),   ALLOCATABLE, INTENT(inout) :: source_p
         REAL(kind=dp), DIMENSION(:,:),   ALLOCATABLE, INTENT(inout) :: source_m
         REAL(kind=dp), DIMENSION(:),     ALLOCATABLE, INTENT(inout) :: eta_l
         REAL(kind=dp), DIMENSION(:),     ALLOCATABLE, INTENT(inout) :: eta_r
         REAL(kind=dp)                               , INTENT(inout) :: eta_boundary_l
         REAL(kind=dp)                               , INTENT(inout) :: eta_boundary_r
       END SUBROUTINE plot_distrf
    END INTERFACE

    CALL plot_distrf(prop_a%p%source_p,prop_a%p%source_m,prop_a%p%eta_l,prop_a%p%eta_r, &
         prop_a%p%eta_boundary_l,prop_a%p%eta_boundary_r)

  END SUBROUTINE plot_distrf_int

  SUBROUTINE join_ripples_int(ierr,cstat_in,deall_in)
    !  
    ! Interface for subroutine join_ripples
    ! Joines current propagator prop_c with previous propagator prop_c%prev
    ! and puts result into the previous one.
    !
    ! Has to be changed if physical content of propagator changes!
    ! (see type declaration of propagator)
    !
    ! Call to external program join_ripples
    USE binarysplit_mod
    
    INTEGER, INTENT(out) :: ierr
    CHARACTER(len=5), OPTIONAL, INTENT(in) :: cstat_in
    INTEGER, OPTIONAL, INTENT(in) :: deall_in

    CHARACTER(len=5)                     :: cstat
    
    TYPE(propagator), POINTER            :: o
    TYPE(propagator), POINTER            :: n
    
    TYPE(binarysplit)                    :: loc_bs_1a,loc_bs_2a
    TYPE(binarysplit)                    :: loc_bs_1b,loc_bs_2b
    
    INTEGER :: i
    INCLUDE 'longint.f90'
    INTEGER(kind=longint), DIMENSION(:,:), ALLOCATABLE :: bin1,bin2
    INTEGER :: deall

    REAL(kind=dp),   DIMENSION(:,:), ALLOCATABLE :: cmat_help
    ierr = 0
    
    o => prop_c%prev
    n => prop_c

    IF (PRESENT(cstat_in)) THEN
       cstat = cstat_in
    ELSE
       cstat = 'inter'
    END IF
    IF (PRESENT(deall_in)) THEN
       deall = deall_in
    ELSE
       deall = 1
    END IF


    ! allocate cmat (c_forward and c_backward)
    ! c_forward
    IF (ALLOCATED(o%p%cmat)) DEALLOCATE(o%p%cmat)
    ALLOCATE( o%p%cmat(o%p%npass_r,o%p%npass_r) )
    o%p%cmat = 0.0_dp
    DO i = 1,o%p%npass_r
       o%p%cmat(i,i) = 1.0_dp
    END DO
    ! c_backward
    IF (ALLOCATED(n%p%cmat)) DEALLOCATE(n%p%cmat)
    ALLOCATE( n%p%cmat(n%p%npass_l,n%p%npass_l) )
    n%p%cmat = 0.0_dp
    DO i = 1,n%p%npass_l
       n%p%cmat(i,i) = 1.0_dp
    END DO
    ! 

    IF (prop_diagnostic .GE. 1) THEN
       PRINT *, ' '
       PRINT *, 'I am in join_ripples'
       PRINT *, ' '
       PRINT *, ' bin_split_mode old new   ', o%bin_split_mode,n%bin_split_mode
       PRINT *, ' ubound(x_split)old new l ', UBOUND(o%eta_bs_l%x_split,1), &
            UBOUND(n%eta_bs_l%x_split,1)
       PRINT *, ' ubound(x_split)old new r ', UBOUND(o%eta_bs_r%x_split,1), &
            UBOUND(n%eta_bs_r%x_split,1)
       PRINT *, ' npart old new            ', o%p%npart,n%p%npart
       PRINT *, ' npass_l old new          ', o%p%npass_l,n%p%npass_l
       PRINT *, ' npass_r old new          ', o%p%npass_r,n%p%npass_r
       PRINT *, ' '
    END IF
        
    IF (o%bin_split_mode .NE. n%bin_split_mode) THEN
       PRINT *, 'ERROR: A difference of bin_split_mode of 1 and 2'
       PRINT *, '       cannot be handled!'
       STOP
    END IF

    ! now we have to modifiy c_forward and c_backward
    IF (o%bin_split_mode .EQ. 1) THEN

       ! we first fix it forward
       prop_modifyold = 1       
       ! first: remove those splits in old which are not in new
       IF (prop_diagnostic .GE. 3) THEN
          PRINT *, 'JOIN binarysplit action - forward'
       END IF
       ! look for splits which are in o%eta_bs_r and not in n%eta_bs_l
       CALL compare_binarysplit(o%eta_bs_r,n%eta_bs_l,bin1,'diff')
       ! do the joining of levels
       ! use bin1 to remove them from o%eta_bs_r
       CALL join_binarysplit(loc_bs_1a,o%eta_bs_r,bin1)
       ! now loc_bs_1a is the modified o%eta_bs_r

       ! second: create those splits in old which are not in new
       IF (prop_diagnostic .GE. 3) THEN
          PRINT *, 'SPLIT binarysplit action - forward'
       END IF
       ! 
       CALL compare_binarysplit(n%eta_bs_l,loc_bs_1a,bin2,'diff')
       ! do the splitting of levels
       CALL dosplit_binarysplit(loc_bs_2a,loc_bs_1a,bin2)
       ! now loc_bs_2a is the final modified o%eta_bs_r

       ! remove unnecessary things
       IF (ALLOCATED(bin1)) DEALLOCATE(bin1)
       IF (ALLOCATED(bin2)) DEALLOCATE(bin2)

       ! now we fix it backward
       prop_modifyold = 0       
       ! first: remove those splits in new which are not in old
       IF (prop_diagnostic .GE. 3) THEN
          PRINT *, 'JOIN binarysplit action - backward'
       END IF
       CALL compare_binarysplit(n%eta_bs_l,o%eta_bs_r,bin1,'diff')
       ! do the joining of levels
       CALL join_binarysplit(loc_bs_1b,n%eta_bs_l,bin1)
       ! now loc_bs_1b is the modified n%eta_bs_l

       ! second: create those splits in new which are not in old
       IF (prop_diagnostic .GE. 3) THEN
          PRINT *, 'SPLIT binarysplit action - backward'
       END IF
       CALL compare_binarysplit(o%eta_bs_r,loc_bs_1b,bin2,'diff')
       ! do the splitting of levels
       CALL dosplit_binarysplit(loc_bs_2b,loc_bs_1b,bin2)
       ! now loc_bs_2b is the final modified n%eta_bs_l - not needed

       ! put the eta_information on right side of propagator
       ! o%eta_bs_r = loc_bs_2a 
       ! CALL get_binarysplit(loc_bs_2a,o%p%eta_r,'x')
       ! put the eta_information on left side of propagator
       ! n%eta_bs_l = loc_bs_2b 
       ! CALL get_binarysplit(loc_bs_2b,n%p%eta_l,'x')

       ! remove unnecessary things
       IF (ALLOCATED(bin1)) DEALLOCATE(bin1)
       IF (ALLOCATED(bin2)) DEALLOCATE(bin2)
       
       IF (cstat .EQ. 'final') THEN
       IF (prop_diagnostic .GE. 2) THEN     
          OPEN(unit=1000,file='c_forward.dat')
          DO i = 1,SIZE(o%p%cmat,1)
             WRITE(1000,'(1000e14.5)') o%p%cmat(i,:)
          END DO
          CLOSE(unit=1000)
          OPEN(unit=1000,file='c_backward.dat')
          DO i = 1,SIZE(n%p%cmat,1)
             WRITE(1000,'(1000e14.5)') n%p%cmat(i,:)
          END DO
          CLOSE(unit=1000)

          ALLOCATE(cmat_help( n%p%npass_l,n%p%npass_l ))
          OPEN(unit=1000,file='forward_backward.dat')
          cmat_help = MATMUL(o%p%cmat,n%p%cmat)
          DO i = 1,n%p%npass_l
             WRITE(1000,'(1000e14.5)') cmat_help(i,:)
          END DO
          CLOSE(unit=1000)
          DEALLOCATE(cmat_help)

          ALLOCATE(cmat_help( o%p%npass_r,o%p%npass_r ))
          OPEN(unit=1000,file='backward_forward.dat')
          cmat_help = MATMUL(n%p%cmat,o%p%cmat)
          DO i = 1,o%p%npass_r
             WRITE(1000,'(1000e14.5)') cmat_help(i,:)
          END DO
          CLOSE(unit=1000)
          DEALLOCATE(cmat_help)

          PRINT *,  'size of c_forward:    ', SIZE(o%p%cmat,1),SIZE(o%p%cmat,2)
          PRINT *,  ' should be            ', n%p%npass_l, o%p%npass_r
          PRINT *,  'size of c_backward:   ', SIZE(n%p%cmat,1),SIZE(n%p%cmat,2)
          PRINT *,  ' should be            ', o%p%npass_r, n%p%npass_l

          PRINT *,  'old propagator tag:   ', o%fieldpropagator_tag_s,o%fieldpropagator_tag_e
          PRINT *,  'size of o%p%amat_m_p: ', SIZE(o%p%amat_m_p,1),SIZE(o%p%amat_m_p,2)
!->out          PRINT *,  ' should be            ', o%p%npass_r, o%p%npass_r
          PRINT *,  ' should be            ', o%p%npass_r*(o%p%nvelocity+1), o%p%npass_r*(o%p%nvelocity+1) !<-in

          PRINT *,  'new propagator tag:   ', n%fieldpropagator_tag_s,n%fieldpropagator_tag_e
          PRINT *,  'size of n%p%amat_p_m: ', SIZE(n%p%amat_p_m,1),SIZE(n%p%amat_p_m,2)
!->out          PRINT *,  ' should be            ', n%p%npass_l, n%p%npass_l
          PRINT *,  ' should be            ', n%p%npass_l*(n%p%nvelocity+1), n%p%npass_l*(n%p%nvelocity+1)  !<-in
          

          PRINT *, 'Difference in forward'
          DO i = 0,UBOUND(loc_bs_2a%x_split,1)
             IF (loc_bs_2a%x_split(i) .NE. n%eta_bs_l%x_split(i)) THEN
                PRINT *, 'loc_bs not equal: ',i,loc_bs_2a%x_split(i),n%eta_bs_l%x_split(i)
             END IF
          END DO

          PRINT *, 'Difference in backward'
          DO i = 0,UBOUND(loc_bs_2b%x_split,1)
             IF (loc_bs_2b%x_split(i) .NE. o%eta_bs_r%x_split(i)) THEN
                PRINT *, 'loc_bs not equal: ',i,loc_bs_2b%x_split(i),o%eta_bs_r%x_split(i)
             END IF
          END DO

          CALL compare_binarysplit(loc_bs_2a,n%eta_bs_l,bin1,'diff')
          !CALL printbin_binarysplit(bin1)
          PRINT *, 'Count difference forward:  ', COUNT(bin1 .NE. 0)
          IF (ALLOCATED(bin1)) DEALLOCATE(bin1)
          CALL compare_binarysplit(loc_bs_2b,o%eta_bs_r,bin1,'diff')
          !CALL printbin_binarysplit(bin1)
          PRINT *, 'Count difference backward:   ', COUNT(bin1 .NE. 0)
          IF (ALLOCATED(bin1)) DEALLOCATE(bin1)

          PRINT *, 'c_forward.dat and c_backward.dat written'
          !PAUSE
       END IF       
       END IF
       CALL deconstruct_binarysplit(loc_bs_1a)
       CALL deconstruct_binarysplit(loc_bs_2a)
       CALL deconstruct_binarysplit(loc_bs_1b)
       CALL deconstruct_binarysplit(loc_bs_2b)
    END IF



    
    IF (prop_diagnostic .GE. 2) THEN
       PRINT *, ' '
       PRINT *, ' I am now after binarysplit joining amd splitting'
       PRINT *, ' npart old new          ', o%p%npart,n%p%npart
       PRINT *, ' npass_l old new        ', o%p%npass_l,n%p%npass_l
       PRINT *, ' npass_r old new        ', o%p%npass_r,n%p%npass_r
       PRINT *, ' '    
    END IF

    ! here the two array c_forward and c_backward should be finished
    ! and passed on to the join_ripples program
    CALL join_ripples_nn(ierr,cstat)
    IF (ierr .NE. 0) THEN
       PRINT *, 'Error from Joining, ierr=',ierr
       STOP
    END IF

    IF (cstat .EQ. 'inter') THEN
       o%fieldperiod_tag_e = n%fieldperiod_tag_e
       o%fieldpropagator_tag_e = n%fieldpropagator_tag_e
    END IF


    ! final cleaning
    IF (deall .EQ. 1) THEN
       IF (ALLOCATED(o%p%cmat)) DEALLOCATE(o%p%cmat)
       IF (ALLOCATED(n%p%cmat)) DEALLOCATE(n%p%cmat)
    END IF

    NULLIFY(o)
    NULLIFY(n)
    
    !
    RETURN
  END SUBROUTINE join_ripples_int
  ! ---------------------------------------------------------------------------

  SUBROUTINE write_propagator_cont(o,prop_type,prop_showall_in)
    TYPE(propagator), POINTER  :: o

    INTEGER, INTENT(in) :: prop_type
    INTEGER, INTENT(in),OPTIONAL  :: prop_showall_in

    INTEGER :: prop_bound
    INTEGER :: prop_start
    INTEGER :: prop_end
    INTEGER :: prop_showall

    IF (PRESENT(prop_showall_in)) THEN
       prop_showall = prop_showall_in
    ELSE
       prop_showall = 1
    END IF

    prop_bound = 0
    IF (prop_type .EQ. 1) THEN
       prop_start = o%fieldperiod_tag_s
       prop_end   = o%fieldperiod_tag_e
    ELSEIF (prop_type .EQ. 2 .OR. prop_type .EQ. 4) THEN
       prop_start = 0
       prop_end   = 0
    ELSEIF (prop_type .EQ. 3) THEN ! final
       prop_start = o%fieldpropagator_tag_s
       prop_end   = o%fieldpropagator_tag_e
    ELSE
       PRINT *, 'Propagator Writing: prop_type not implemented: ',prop_type
       RETURN
    END IF

    CALL filename_propagator(prop_type,prop_bound,prop_start,prop_end) 
    CALL unit_propagator
    OPEN(unit=prop_unit,file=prop_cfilename,status='replace', &
         form=prop_format,action='write')
    ! tags
    WRITE(prop_unit,*) prop_start
    WRITE(prop_unit,*) prop_end
    
    ! info
    IF (prop_showall .EQ. 1) THEN
       WRITE(prop_unit,*) o%nr_joined
       WRITE(prop_unit,*) o%fieldpropagator_tag_s
       WRITE(prop_unit,*) o%fieldpropagator_tag_e
       WRITE(prop_unit,*) o%fieldperiod_tag_s
       WRITE(prop_unit,*) o%fieldperiod_tag_e
       IF (ALLOCATED(o%y)) THEN
          WRITE(prop_unit,*) LBOUND(o%y,1),UBOUND(o%y,1)
          WRITE(prop_unit,*) o%y
       ELSE
          WRITE(prop_unit,*) 0,0
       END IF
       WRITE(prop_unit,*) o%phi_l
       WRITE(prop_unit,*) o%phi_r
    END IF

    ! Binarysplit stuff is not dumped
    IF (prop_showall .EQ. 0) THEN
       WRITE(prop_unit,*) o%bin_split_mode
    END IF

    ! sizes
    IF (prop_showall .GE. 1) THEN
       WRITE(prop_unit,*) o%p%npart
       WRITE(prop_unit,*) o%p%npass_l
       WRITE(prop_unit,*) o%p%npass_r
       WRITE(prop_unit,*) o%p%nvelocity
    END IF
    ! amat_p_p
    IF (prop_showall .GE. 1) THEN
       IF (ALLOCATED(o%p%amat_p_p)) THEN
          WRITE(prop_unit,*) LBOUND(o%p%amat_p_p,1),UBOUND(o%p%amat_p_p,1)
          WRITE(prop_unit,*) LBOUND(o%p%amat_p_p,2),UBOUND(o%p%amat_p_p,2)
          WRITE(prop_unit,*) o%p%amat_p_p
       ELSE
          WRITE(prop_unit,*) 0,0
          WRITE(prop_unit,*) 0,0
       END IF
    END IF
    ! amat_m_m
    IF (prop_showall .EQ. 1) THEN
       IF (ALLOCATED(o%p%amat_m_m)) THEN
          WRITE(prop_unit,*) LBOUND(o%p%amat_m_m,1),UBOUND(o%p%amat_m_m,1)
          WRITE(prop_unit,*) LBOUND(o%p%amat_m_m,2),UBOUND(o%p%amat_m_m,2)
          WRITE(prop_unit,*) o%p%amat_m_m
       ELSE
          WRITE(prop_unit,*) 0,0
          WRITE(prop_unit,*) 0,0
       END IF
    END IF
    ! amat_p_m
    IF (prop_showall .EQ. 1) THEN
       IF (ALLOCATED(o%p%amat_p_m)) THEN
          WRITE(prop_unit,*) LBOUND(o%p%amat_p_m,1),UBOUND(o%p%amat_p_m,1)
          WRITE(prop_unit,*) LBOUND(o%p%amat_p_m,2),UBOUND(o%p%amat_p_m,2)
          WRITE(prop_unit,*) o%p%amat_p_m
       ELSE
          WRITE(prop_unit,*) 0,0
          WRITE(prop_unit,*) 0,0
       END IF
    END IF
    ! amat_m_p
    IF (prop_showall .GE. 1) THEN
       IF (ALLOCATED(o%p%amat_m_p)) THEN
          WRITE(prop_unit,*) LBOUND(o%p%amat_m_p,1),UBOUND(o%p%amat_m_p,1)
          WRITE(prop_unit,*) LBOUND(o%p%amat_m_p,2),UBOUND(o%p%amat_m_p,2)
          WRITE(prop_unit,*) o%p%amat_m_p
       ELSE
          WRITE(prop_unit,*) 0,0
          WRITE(prop_unit,*) 0,0
       END IF
    END IF
    
    ! source_p
    IF (prop_showall .GE. 1) THEN
       IF (ALLOCATED(o%p%source_p)) THEN
          WRITE(prop_unit,*) LBOUND(o%p%source_p,1),UBOUND(o%p%source_p,1)
          WRITE(prop_unit,*) LBOUND(o%p%source_p,2),UBOUND(o%p%source_p,2)
          WRITE(prop_unit,*) o%p%source_p
       ELSE
          WRITE(prop_unit,*) 0,0
          WRITE(prop_unit,*) 0,0
       END IF
    END IF
    ! source_m
    IF (prop_showall .EQ. 1) THEN
       IF (ALLOCATED(o%p%source_m)) THEN
          WRITE(prop_unit,*) LBOUND(o%p%source_m,1),UBOUND(o%p%source_m,1)
          WRITE(prop_unit,*) LBOUND(o%p%source_m,2),UBOUND(o%p%source_m,2)
          WRITE(prop_unit,*) o%p%source_m
       ELSE
          WRITE(prop_unit,*) 0,0
          WRITE(prop_unit,*) 0,0
       END IF
    END IF

    
    ! flux_p
    IF (prop_showall .EQ. 1) THEN
       IF (ALLOCATED(o%p%flux_p)) THEN
          WRITE(prop_unit,*) LBOUND(o%p%flux_p,1),UBOUND(o%p%flux_p,1)
          WRITE(prop_unit,*) LBOUND(o%p%flux_p,2),UBOUND(o%p%flux_p,2)
          WRITE(prop_unit,*) o%p%flux_p
       ELSE
          WRITE(prop_unit,*) 0,0
          WRITE(prop_unit,*) 0,0
       END IF
    END IF
    ! flux_m
    IF (prop_showall .EQ. 1) THEN
       IF (ALLOCATED(o%p%flux_m)) THEN
          WRITE(prop_unit,*) LBOUND(o%p%flux_m,1),UBOUND(o%p%flux_m,1)
          WRITE(prop_unit,*) LBOUND(o%p%flux_m,2),UBOUND(o%p%flux_m,2)
          WRITE(prop_unit,*) o%p%flux_m
       ELSE
          WRITE(prop_unit,*) 0,0
          WRITE(prop_unit,*) 0,0
       END IF
    END IF

    ! qflux
    IF (prop_showall .EQ. 1) THEN
       WRITE(prop_unit,*) o%p%qflux
    END IF

    ! eta
    IF (prop_showall .EQ. 1) THEN
       IF (ALLOCATED(o%p%eta_l)) THEN
          WRITE(prop_unit,*) LBOUND(o%p%eta_l,1),UBOUND(o%p%eta_l,1)
          WRITE(prop_unit,*) o%p%eta_l
       ELSE
          WRITE(prop_unit,*) 0,0
       END IF
       IF (ALLOCATED(o%p%eta_r)) THEN
          WRITE(prop_unit,*) LBOUND(o%p%eta_r,1),UBOUND(o%p%eta_r,1)
          WRITE(prop_unit,*) o%p%eta_r
       ELSE
          WRITE(prop_unit,*) 0,0
       END IF
       WRITE(prop_unit,*) o%p%eta_boundary_l
       WRITE(prop_unit,*) o%p%eta_boundary_r
    END IF
    
    CLOSE(unit=prop_unit)


  END SUBROUTINE write_propagator_cont
  ! ---------------------------------------------------------------------------

  SUBROUTINE write_prop_bound_cont(o,n,prop_type)
    TYPE(propagator), POINTER  :: o,n

    INTEGER, INTENT(in) :: prop_type
    INTEGER :: prop_bound
    INTEGER :: prop_left
    INTEGER :: prop_right

    prop_bound = 1
    IF (prop_type .EQ. 1) THEN
       prop_right = n%fieldperiod_tag_s
       !prop_left  = o%fieldperiod_tag_e
       prop_left  = prop_right - 1
    ELSEIF (prop_type .EQ. 3) THEN
       prop_right = n%fieldpropagator_tag_s
       !prop_left  = o%fieldperiod_tag_e
       prop_left  = prop_right - 1
    ELSE
       PRINT *, 'Propagator Writing: prop_type not implemented: ',prop_type
       RETURN
    END IF

    CALL filename_propagator(prop_type,prop_bound,prop_left,prop_right) 
    CALL unit_propagator
    OPEN(unit=prop_unit,file=prop_cfilename,status='replace', &
         form=prop_format,action='write')
    ! tags
    WRITE(prop_unit,*) n%fieldpropagator_tag_s - 1
    WRITE(prop_unit,*) n%fieldpropagator_tag_s
    WRITE(prop_unit,*) n%fieldperiod_tag_s - 1
    WRITE(prop_unit,*) n%fieldperiod_tag_s
    ! forward
    IF (ALLOCATED(o%p%cmat)) THEN
       WRITE(prop_unit,*) LBOUND(o%p%cmat,1),UBOUND(o%p%cmat,1)
       WRITE(prop_unit,*) LBOUND(o%p%cmat,2),UBOUND(o%p%cmat,2)
       WRITE(prop_unit,*) o%p%cmat
    ELSE
       WRITE(prop_unit,*) 0,0
       WRITE(prop_unit,*) 0,0
    END IF
    ! backward
    IF (ALLOCATED(n%p%cmat)) THEN
       WRITE(prop_unit,*) LBOUND(n%p%cmat,1),UBOUND(n%p%cmat,1)
       WRITE(prop_unit,*) LBOUND(n%p%cmat,2),UBOUND(n%p%cmat,2)
       WRITE(prop_unit,*) n%p%cmat    
    ELSE
       WRITE(prop_unit,*) 0,0
       WRITE(prop_unit,*) 0,0
    END IF

    CLOSE(unit=prop_unit)
    
  END SUBROUTINE write_prop_bound_cont
  ! ---------------------------------------------------------------------------

  SUBROUTINE read_propagator_cont(o,prop_type,prop_start,prop_end,prop_showall_in)
    TYPE(propagator), POINTER  :: o
    INTEGER, INTENT(in), OPTIONAL :: prop_showall_in

    INTEGER, INTENT(in) :: prop_type
    INTEGER, INTENT(in) :: prop_start
    INTEGER, INTENT(in) :: prop_end

    INTEGER :: prop_bound
    INTEGER :: lb1,ub1
    INTEGER :: lb2,ub2
    INTEGER :: dummy

    INTEGER :: prop_showall

    IF (PRESENT(prop_showall_in)) THEN
       prop_showall = prop_showall_in
    ELSE
       prop_showall = 1
    END IF

    prop_bound = 0

    CALL filename_propagator(prop_type,prop_bound,prop_start,prop_end)
    CALL unit_propagator

    !PRINT *, prop_cfilename

    OPEN(unit=prop_unit,file=prop_cfilename,status='old', &
         form=prop_format,action='read')

    ! tags
    READ(prop_unit,*) dummy
    READ(prop_unit,*) dummy
    
    ! info
    IF (prop_showall .EQ. 1) THEN
       READ(prop_unit,*) o%nr_joined
       READ(prop_unit,*) o%fieldpropagator_tag_s
       READ(prop_unit,*) o%fieldpropagator_tag_e
       READ(prop_unit,*) o%fieldperiod_tag_s
       READ(prop_unit,*) o%fieldperiod_tag_e

       READ(prop_unit,*) lb1,ub1
       IF (ub1 .GT. 0) THEN
          IF (ALLOCATED(o%y)) DEALLOCATE(o%y)
          ALLOCATE(o%y(lb1:ub1))
          READ(prop_unit,*) o%y
       END IF
       READ(prop_unit,*) o%phi_l
       READ(prop_unit,*) o%phi_r
    END IF

    ! Binarysplit stuff is not dumped
    IF (prop_showall .EQ. 0) THEN
       READ(prop_unit,*) o%bin_split_mode
    END IF

    ! sizes
    IF (prop_showall .GE. 1) THEN
       READ(prop_unit,*) o%p%npart
       READ(prop_unit,*) o%p%npass_l
       READ(prop_unit,*) o%p%npass_r
       READ(prop_unit,*) o%p%nvelocity
    END IF

    ! amat_p_p
    IF (prop_showall .GE. 1) THEN
       READ(prop_unit,*) lb1,ub1
       READ(prop_unit,*) lb2,ub2
       IF (ub1 .GT. 0 .AND. ub2 .GT. 0) THEN
          IF (ALLOCATED(o%p%amat_p_p)) DEALLOCATE(o%p%amat_p_p)
          ALLOCATE(o%p%amat_p_p(lb1:ub1,lb2:ub2))
          READ(prop_unit,*) o%p%amat_p_p
       END IF
    END IF
    ! amat_m_m
    IF (prop_showall .EQ. 1) THEN
       READ(prop_unit,*) lb1,ub1
       READ(prop_unit,*) lb2,ub2
       IF (ub1 .GT. 0 .AND. ub2 .GT. 0) THEN
          IF (ALLOCATED(o%p%amat_m_m)) DEALLOCATE(o%p%amat_m_m)
          ALLOCATE(o%p%amat_m_m(lb1:ub1,lb2:ub2))
          READ(prop_unit,*) o%p%amat_m_m
       END IF
    END IF
    ! amat_p_m
    IF (prop_showall .EQ. 1) THEN
       READ(prop_unit,*) lb1,ub1
       READ(prop_unit,*) lb2,ub2
       IF (ub1 .GT. 0 .AND. ub2 .GT. 0) THEN
          IF (ALLOCATED(o%p%amat_p_m)) DEALLOCATE(o%p%amat_p_m)
          ALLOCATE(o%p%amat_p_m(lb1:ub1,lb2:ub2))
          READ(prop_unit,*) o%p%amat_p_m
       END IF
    END IF
    ! amat_m_p
    IF (prop_showall .GE. 1) THEN
       READ(prop_unit,*) lb1,ub1
       READ(prop_unit,*) lb2,ub2
       IF (ub1 .GT. 0 .AND. ub2 .GT. 0) THEN
          IF (ALLOCATED(o%p%amat_m_p)) DEALLOCATE(o%p%amat_m_p)
          ALLOCATE(o%p%amat_m_p(lb1:ub1,lb2:ub2))
          READ(prop_unit,*) o%p%amat_m_p
       END IF
    END IF
    
    ! source_p
    IF (prop_showall .GE. 1) THEN
       READ(prop_unit,*) lb1,ub1
       READ(prop_unit,*) lb2,ub2
       IF (ub1 .GT. 0 .AND. ub2 .GT. 0) THEN
          IF (ALLOCATED(o%p%source_p)) DEALLOCATE(o%p%source_p)
          ALLOCATE(o%p%source_p(lb1:ub1,lb2:ub2))
          READ(prop_unit,*) o%p%source_p
       END IF
    END IF
    ! source_m
    IF (prop_showall .EQ. 1) THEN
       READ(prop_unit,*) lb1,ub1
       READ(prop_unit,*) lb2,ub2
       IF (ub1 .GT. 0 .AND. ub2 .GT. 0) THEN
          IF (ALLOCATED(o%p%source_m)) DEALLOCATE(o%p%source_m)
          ALLOCATE(o%p%source_m(lb1:ub1,lb2:ub2))
          READ(prop_unit,*) o%p%source_m
       END IF
    END IF

    ! flux_p
    IF (prop_showall .EQ. 1) THEN
       READ(prop_unit,*) lb1,ub1
       READ(prop_unit,*) lb2,ub2
       IF (ub1 .GT. 0 .AND. ub2 .GT. 0) THEN
          IF (ALLOCATED(o%p%flux_p)) DEALLOCATE(o%p%flux_p)
          ALLOCATE(o%p%flux_p(lb1:ub1,lb2:ub2))
          READ(prop_unit,*) o%p%flux_p
       END IF
    END IF
    ! flux_m
    IF (prop_showall .EQ. 1) THEN
       READ(prop_unit,*) lb1,ub1
       READ(prop_unit,*) lb2,ub2
       IF (ub1 .GT. 0 .AND. ub2 .GT. 0) THEN
          IF (ALLOCATED(o%p%flux_m)) DEALLOCATE(o%p%flux_m)
          ALLOCATE(o%p%flux_m(lb1:ub1,lb2:ub2))
          READ(prop_unit,*) o%p%flux_m
       END IF
    END IF

    ! qflux
    IF (prop_showall .EQ. 1) THEN
       READ(prop_unit,*) o%p%qflux
    END IF

    ! eta
    IF (prop_showall .EQ. 1) THEN
       READ(prop_unit,*) lb1,ub1
       IF (ub1 .GT. 0) THEN
          IF (ALLOCATED(o%p%eta_l)) DEALLOCATE(o%p%eta_l)
          ALLOCATE(o%p%eta_l(lb1:ub1))
          READ(prop_unit,*) o%p%eta_l
       END IF
       READ(prop_unit,*) lb1,ub1
       IF (ub1 .GT. 0) THEN
          IF (ALLOCATED(o%p%eta_r)) DEALLOCATE(o%p%eta_r)
          ALLOCATE(o%p%eta_r(lb1:ub1))
          READ(prop_unit,*) o%p%eta_r
       END IF
       READ(prop_unit,*) o%p%eta_boundary_l
       READ(prop_unit,*) o%p%eta_boundary_r
    END IF
    
    CLOSE(unit=prop_unit)

  END SUBROUTINE read_propagator_cont
  ! ---------------------------------------------------------------------------

  SUBROUTINE read_prop_bound_cont(b,prop_type,prop_left,prop_right)
    TYPE(prop_boundary)  :: b
    
    INTEGER, INTENT(in) :: prop_type
    INTEGER, INTENT(in) :: prop_left
    INTEGER, INTENT(in) :: prop_right

    INTEGER :: prop_bound
    INTEGER :: lb1,ub1
    INTEGER :: lb2,ub2

    prop_bound = 1

    CALL filename_propagator(prop_type,prop_bound,prop_left,prop_right)
    CALL unit_propagator
    OPEN(unit=prop_unit,file=prop_cfilename,status='old', &
         form=prop_format,action='read')
    ! tags
    READ(prop_unit,*) b%fieldpropagator_tag_left
    READ(prop_unit,*) b%fieldpropagator_tag_right
    READ(prop_unit,*) b%fieldperiod_tag_left
    READ(prop_unit,*) b%fieldperiod_tag_right
    ! forward
    READ(prop_unit,*) lb1,ub1
    READ(prop_unit,*) lb2,ub2
    IF (ub1 .GT. 0 .AND. ub2 .GT. 0) THEN
       IF (ALLOCATED(b%c_forward)) DEALLOCATE(b%c_forward)
       ALLOCATE(b%c_forward(lb1:ub1,lb2:ub2))
       READ(prop_unit,*) b%c_forward
    END IF
    ! backward
    READ(prop_unit,*) lb1,ub1
    READ(prop_unit,*) lb2,ub2
    IF (ub1 .GT. 0 .AND. ub2 .GT. 0) THEN
       IF (ALLOCATED(b%c_backward)) DEALLOCATE(b%c_backward)
       ALLOCATE(b%c_backward(lb1:ub1,lb2:ub2))
       READ(prop_unit,*) b%c_backward
    END IF
   
    CLOSE(unit=prop_unit)
  
  END SUBROUTINE read_prop_bound_cont
  ! ---------------------------------------------------------------------------

  SUBROUTINE read_prop_recon_cont(tag)
    INTEGER, INTENT(in) :: tag
    INTEGER :: dummy,lb1,ub1,lb2,ub2

    CALL filename_propagator(5,0,tag,tag)
    OPEN(unit=prop_unit,file=prop_cfilename,status='old', &
         form=prop_format,action='read')
    READ(prop_unit,*) dummy
    ! flux_mr for tag
    READ(prop_unit,*) lb1,ub1
    READ(prop_unit,*) lb2,ub2
    IF (ALLOCATED(flux_mr)) DEALLOCATE(flux_mr)
    ALLOCATE(flux_mr(lb1:ub1,lb2:ub2))
    READ(prop_unit,*) flux_mr
    ! flux_pl for tag
    READ(prop_unit,*) lb1,ub1
    READ(prop_unit,*) lb2,ub2
    IF (ALLOCATED(flux_pl)) DEALLOCATE(flux_pl)
    ALLOCATE(flux_pl(lb1:ub1,lb2:ub2))
    READ(prop_unit,*) flux_pl
    CLOSE(unit=prop_unit)

  END SUBROUTINE read_prop_recon_cont

  ! ---------------------------------------------------------------------------
  SUBROUTINE reconstruct_propagator_dist

    TYPE(propagator), POINTER            :: l  ! left
    TYPE(propagator), POINTER            :: r  ! right
    TYPE(prop_boundary)  :: b

    INTEGER :: prop_type
    INTEGER :: prop_start
    INTEGER :: prop_end
    INTEGER :: prop_left
    INTEGER :: prop_right
    INTEGER :: prop_showall

    INTEGER :: N,i
    INTEGER :: lb1,ub1,lb2,ub2

    INTEGER :: ierr_join

    ! from the finally joined propagator
    REAL(kind=dp), DIMENSION(:,:),   ALLOCATABLE :: source_p_0               
    REAL(kind=dp), DIMENSION(:,:),   ALLOCATABLE :: source_m_N  
    ! 
    REAL(kind=dp), DIMENSION(:,:),   ALLOCATABLE :: source_p_N               
    REAL(kind=dp), DIMENSION(:,:),   ALLOCATABLE :: source_m_N1  

    ! read the information about tags
    CALL unit_propagator
    OPEN(unit=prop_unit,file=prop_ctaginfo,status='old', &
         form=prop_format,action='read')
    READ(prop_unit,*) prop_write
    READ(prop_unit,*) prop_first_tag
    READ(prop_unit,*) prop_last_tag
    CLOSE(unit=prop_unit)
    
    IF (prop_write .EQ. 1) THEN
       prop_type = 1
    ELSEIF (prop_write .EQ. 2) THEN
       prop_type = 3
    ELSE
       PRINT *, 'No reconstruction possible'
       PRINT *, 'Set prop_write to 1 (period) or 2 (propagator)'
       STOP
    END IF

    ! read the final joined propagator (from join_ends)
    ! and keep only the fluxes (which are in source)
    ! (only one propagator exists at this moment)
    CALL construct_propagator
    l => prop_c
    prop_start = 0
    prop_end = 0
    prop_showall = 1
    CALL read_propagator_content(l,prop_type,prop_start,prop_end,prop_showall)
    lb1 = LBOUND(l%p%source_p,1)
    ub1 = UBOUND(l%p%source_p,1)
    lb2 = LBOUND(l%p%source_p,2)
    ub2 = UBOUND(l%p%source_p,2)
    ALLOCATE(source_p_0(lb1:ub1,lb2:ub2))
    source_p_0 = l%p%source_p
    lb1 = LBOUND(l%p%source_m,1)
    ub1 = UBOUND(l%p%source_m,1)
    lb2 = LBOUND(l%p%source_m,2)
    ub2 = UBOUND(l%p%source_m,2)
    ALLOCATE(source_m_N(lb1:ub1,lb2:ub2))
    source_m_N = l%p%source_m
    NULLIFY(l)

    !PRINT *, 'source_p_0', SIZE(source_p_0,1),SIZE(source_p_0,2)
    !PRINT *, 'source_m_N', SIZE(source_m_N,1),SIZE(source_m_N,2)

    ! write the results - starting point
    CALL filename_propagator(5,0,prop_last_tag,prop_last_tag)
    CALL unit_propagator
    OPEN(unit=prop_unit,file=prop_cfilename,status='replace', &
         form=prop_format,action='write')
    WRITE(prop_unit,*) prop_last_tag
    WRITE(prop_unit,*) LBOUND(source_m_N,1),UBOUND(source_m_N,1)
    WRITE(prop_unit,*) LBOUND(source_m_N,2),UBOUND(source_m_N,2)
    WRITE(prop_unit,*) source_m_N
    CLOSE(unit=prop_unit)
    
    ! propagators
    r => prop_c ! the first one is the right one
    CALL construct_propagator(1)
    l => prop_c ! the second one is the left one

    DO N = prop_last_tag, prop_first_tag + 1, -1

       PRINT *, N,prop_last_tag, prop_first_tag
       
       ! now read the propagator N into the right (r)
       prop_start = N
       prop_end   = N
       prop_showall = 1
       CALL read_propagator_content(r,prop_type,prop_start,prop_end,prop_showall)
    
       ! read the one from the first to (N-1) into the left (l) propagator
       prop_start = prop_first_tag
       prop_end     = N - 1
       prop_showall = 2
       IF (prop_start .EQ. prop_end) prop_showall = 1
       CALL read_propagator_content(l,prop_type,prop_start,prop_end,prop_showall)

       ! read the boundary between those two
       prop_left  = N - 1
       prop_right = N
       CALL read_prop_bound_cont(b,prop_type,prop_left,prop_right)

       ! do the computation with left (l), right (r), boundary (b)
       CALL reconstruct_prop_dist(l,r,b, &
            source_p_0,source_m_N,source_m_N1,source_p_N)
       
       ! output of new results
       CALL filename_propagator(5,0,N,N)
       CALL unit_propagator
       OPEN(unit=prop_unit,file=prop_cfilename,status='old', &
            form=prop_format,action='write',position='append')
       ! source_p for N
       WRITE(prop_unit,*) LBOUND(source_p_N,1),UBOUND(source_p_N,1)
       WRITE(prop_unit,*) LBOUND(source_p_N,2),UBOUND(source_p_N,2)
       WRITE(prop_unit,*) source_p_N
       CLOSE(prop_unit)

       CALL filename_propagator(5,0,N-1,N-1)
       CALL unit_propagator
       OPEN(unit=prop_unit,file=prop_cfilename,status='replace', &
            form=prop_format,action='write')
       WRITE(prop_unit,*) N - 1
       ! source_m for N-1       
       WRITE(prop_unit,*) LBOUND(source_m_N1,1),UBOUND(source_m_N1,1)
       WRITE(prop_unit,*) LBOUND(source_m_N1,2),UBOUND(source_m_N1,2)
       WRITE(prop_unit,*) source_m_N1
       CLOSE(unit=prop_unit)
      
       ! now make source_m or N-1 the new starting value
       DEALLOCATE(source_m_N)
!       ALLOCATE(source_m_N(lb1:ub1,lb2:ub2))
       ALLOCATE(source_m_N(SIZE(source_m_N1,1),SIZE(source_m_N1,2))) !<-SERGEI
       source_m_N = source_m_N1
       ! and continue with the backward recurance
    END DO

    ! write the results - final point (from join_ends)
    CALL filename_propagator(5,0,prop_first_tag,prop_first_tag)
    CALL unit_propagator
    OPEN(unit=prop_unit,file=prop_cfilename,status='old', &
         form=prop_format,action='write',position='append')
    WRITE(prop_unit,*) LBOUND(source_p_0,1),UBOUND(source_p_0,1)
    WRITE(prop_unit,*) LBOUND(source_p_0,2),UBOUND(source_p_0,2)
    WRITE(prop_unit,*) source_p_0
    CLOSE(unit=prop_unit)



!!$    ! for joining the cmat must go into the propagator
!!$    ! forward goes to the left - l(eft)
!!$    lb1 = LBOUND(b%c_forward,1)
!!$    ub1 = UBOUND(b%c_forward,1)
!!$    lb2 = LBOUND(b%c_forward,2)
!!$    ub2 = UBOUND(b%c_forward,2)
!!$    IF (ALLOCATED(l%p%cmat)) DEALLOCATE(l%p%cmat)
!!$    ALLOCATE(l%p%cmat(lb1:ub1,lb2:ub2))
!!$    l%p%cmat = b%c_forward
!!$    ! backward goes to the right - r(ight)
!!$    lb1 = LBOUND(b%c_backward,1)
!!$    ub1 = UBOUND(b%c_backward,1)
!!$    lb2 = LBOUND(b%c_backward,2)
!!$    ub2 = UBOUND(b%c_backward,2)
!!$    IF (ALLOCATED(r%p%cmat)) DEALLOCATE(r%p%cmat)
!!$    ALLOCATE(r%p%cmat(lb1:ub1,lb2:ub2))
!!$    r%p%cmat = b%c_backward
!!$
!!$    ! when you do joining the result is in prop_c
!!$    prop_c => prop_c%next
!!$    CALL join_ripples(ierr_join)




  END SUBROUTINE reconstruct_propagator_dist

  ! ---------------------------------------------------------------------------
  SUBROUTINE reconstruct_propagator_dist_1(l,r,b, &
       source_p_0,source_m_N,source_m_N1,source_p_N)
!
    USE lapack_band
!
    TYPE(propagator), POINTER            :: l  ! left
    TYPE(propagator), POINTER            :: r  ! right
    TYPE(prop_boundary)  :: b

    INTEGER :: lb1,ub1,lb2,ub2
    INTEGER :: nvel,ndim,ndim1,m,k,k1,info,i


    ! from the finally joined propagator
    REAL(kind=dp), DIMENSION(:,:),   ALLOCATABLE :: source_p_0               
    ! from the last iteration
    REAL(kind=dp), DIMENSION(:,:),   ALLOCATABLE :: source_m_N  
    ! new results
    REAL(kind=dp), DIMENSION(:,:),   ALLOCATABLE :: source_m_N1  
    REAL(kind=dp), DIMENSION(:,:),   ALLOCATABLE :: source_p_N               

    REAL(kind=dp), DIMENSION(:,:),   ALLOCATABLE :: amat,bvec_lapack
    REAL(kind=dp), DIMENSION(:,:),   ALLOCATABLE :: a_mp,a_pm,q_p,q_m
    INTEGER,       DIMENSION(:),     ALLOCATABLE :: ipivot

    IF (ALLOCATED(source_m_N1)) DEALLOCATE(source_m_N1)
    IF (ALLOCATED(source_p_N))  DEALLOCATE(source_p_N)

    ! one has
    !  source_m_N
    !  source_p_0

    ! so what one needs is
    ! source_m_N1     - N1 means N-1
    ! source_p_N

    ! there are two pointers
    ! l : jpoined from the beginning to N-1 
    !     (with limited information) A_p_p,A_m_p,q_p and the sizes)
    ! r : the actual propagator at N 
    !     (with all information) 

    ! there is also the boundary
    ! b : with
    ! b%c_forward
    ! b%c_backward

    ! SERGEI - here comes all the stuff
  nvel = r%p%nvelocity
!
  ndim=r%p%npass_l*(nvel+1)
  ndim1=l%p%npass_r*(nvel+1)
!
  ALLOCATE(amat(ndim,ndim),bvec_lapack(ndim,3),ipivot(ndim))
  ALLOCATE(a_mp(ndim,ndim1),a_pm(ndim1,ndim),q_p(ndim,3),q_m(ndim1,3))
!
  DO m=0,nvel
    k=m*r%p%npass_l
    k1=m*l%p%npass_r
    a_mp(k+1:k+r%p%npass_l,:)                                       &
           =MATMUL(b%c_forward,l%p%amat_m_p(k1+1:k1+l%p%npass_r,:))
    a_pm(k1+1:k1+l%p%npass_r,:)                                     &
           =MATMUL(b%c_backward,r%p%amat_p_m(k+1:k+r%p%npass_l,:))
    q_p(k+1:k+r%p%npass_l,:)                                        &
           =MATMUL(b%c_forward,l%p%source_p(k1+1:k1+l%p%npass_r,:)  &
           +MATMUL(l%p%amat_p_p(k1+1:k1+l%p%npass_r,:),source_p_0))
    q_m(k1+1:k1+l%p%npass_r,:)                                      &
           =MATMUL(b%c_backward,r%p%source_m(k+1:k+r%p%npass_l,:)   &
           +MATMUL(r%p%amat_m_m(k+1:k+r%p%npass_l,:),source_m_N))
  ENDDO
!
! Now the set of equations has the aligned dimensions and is of the form:
! $f^+ = A^{-+} f^- + q^+$
! $f^- = A^{+-} f^+ + q^-$
!
  amat=0.d0
  DO i=1,ndim
    amat(i,i)=1.d0
  ENDDO
  amat=amat-MATMUL(a_mp,a_pm)
  bvec_lapack=q_p+MATMUL(a_mp,q_m)
!
  CALL gbsv(ndim,ndim,amat,ipivot,bvec_lapack,info)
!
  IF(info.NE.0) THEN
    PRINT *,'gbsv error ',info,' in reconstruct_propagator_dist'
    RETURN
  ENDIF
!
  ALLOCATE(source_p_N(ndim,3),source_m_N1(ndim1,3))
!
  source_p_N=bvec_lapack
  source_m_N1=q_m+MATMUL(a_pm,bvec_lapack)
!
  DEALLOCATE(amat,bvec_lapack,ipivot,a_mp,a_pm,q_p,q_m)
!

  END SUBROUTINE reconstruct_propagator_dist_1
  ! ---------------------------------------------------------------------------

  SUBROUTINE unit_prop
    LOGICAL :: opened
    DO
       INQUIRE(unit=prop_unit,opened=opened)
       IF (.NOT. opened) EXIT
       prop_unit = prop_unit + 1
    END DO
  END SUBROUTINE unit_prop
  ! ---------------------------------------------------------------------------

  SUBROUTINE filename_prop(prop_type,prop_bound,prop_start,prop_end)

    INTEGER, INTENT(in) :: prop_type
    INTEGER, INTENT(in) :: prop_bound
    INTEGER, INTENT(in) :: prop_start
    INTEGER, INTENT(in) :: prop_end

    CHARACTER(len=15) :: ctag1,ctag2

    WRITE(ctag1,*) prop_start
    WRITE(ctag2,*) prop_end

    ! choose basename
    IF (prop_type .EQ. 1 .OR. prop_type .EQ. 2) THEN ! period
       prop_cfilename = prop_cperiod
    ELSEIF (prop_type .EQ. 3 .OR. prop_type .EQ. 4) THEN ! propagator
       prop_cfilename = prop_cpropagator
    ELSEIF (prop_type .EQ. 5) THEN ! propagator
       prop_cfilename = prop_cresult
    END IF

    ! add boundary to name
    IF (prop_bound .EQ. 1) THEN
       WRITE(prop_cfilename,'(100A)') &
         TRIM(ADJUSTL(prop_cfilename)),'_', &
         TRIM(ADJUSTL(prop_cboundary))
    END IF

    ! add numbers and extension
    WRITE(prop_cfilename,'(100A)') &
         TRIM(ADJUSTL(prop_cfilename)),'_', &
         TRIM(ADJUSTL(ctag1)),'_', &
         TRIM(ADJUSTL(ctag2)),'.', &
         TRIM(ADJUSTL(prop_cext))


  END SUBROUTINE filename_prop
  ! ---------------------------------------------------------------------------



END MODULE propagator_mod
