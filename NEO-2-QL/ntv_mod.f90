!> Module containing additional data from/for the NTV version
!> of ripple_solver
!>
!> Input (switches, mach number, normalized magnetic drift frequency,
!> boozer_s, collisionality,...) provided by neo2.in
!>
!> Used within programs:
!>   neo2 (main)
!>   ripple_solver
!>   diag_propagator_res
!>   neo_magfie_perturbation
MODULE ntv_mod
  ! module containing numerical constants
  USE neo_precision
  !
  IMPLICIT NONE
  !
  ! Define physical constants (cgs-units)
  REAL(kind=dp), PARAMETER, PUBLIC :: c=2.9979e10_dp       ! speed of light
  REAL(kind=dp), PARAMETER, PUBLIC :: e=4.8032e-10_dp      ! elementary charge
  REAL(kind=dp), PARAMETER, PUBLIC :: u=1.660539040e-24_dp ! atomic mass unit

  real(kind=dp), parameter, public :: epsilon_transport_coefficients = 1.0e-3
  real(kind=dp), parameter, public :: epsilon_particle_flux = 1.0e-6
  !
  ! INPUT
  !> switch: turn on(=1)/off(=0) ntv mode (not used at the moment)
  INTEGER, PUBLIC :: isw_ntv_mode
  !> switch: 0=compute qflux only for the symmetric case; 1=do all computations
  INTEGER, PUBLIC :: isw_qflux_NA
  !> switch for rippler_solver versions
  !> (1=preconditioned; 2=Arnoldi Order 1; 3=Arnoldi Order 2)
  INTEGER, PUBLIC :: isw_ripple_solver
  !> name of perturbation file
  CHARACTER(len=100), PUBLIC :: in_file_pert
  !> toroidal mach number over R_major (Mt/R)
  REAL(kind=dp), PUBLIC :: MtOvR
  !> Larmor radius associated with $B_{00}^{Booz}$ (rho_L_loc) times B
  REAL(kind=dp), PUBLIC :: B_rho_L_loc

  ! ADDITIONAL INPUT FOR MULTI-SPECIES COMPUTATIONS (neo2.in)
  !
  !> switch: turn on(=1)/off(=0) computation of E_r
  INTEGER, PUBLIC :: isw_calc_Er
  !> switch: turn on(=1)/off(=0) computation of magnetic drift
  INTEGER, PUBLIC :: isw_calc_MagDrift
  !> species tag of measured toroidal rotation frequency
  !> (e.g., toroidal rotation frequency measured for main ion species)
  !> -> species velocity, density and temperature must be known at each
  !>   flux surface
  INTEGER, PUBLIC :: species_tag_Vphi
  !> isw_Vphi_loc=0: value of "Vphi" corresponds to flux surface average (<V_\varphi>)
  !> isw_Vphi_loc=1: value of "Vphi" is specified locally for given (R,Z)-position
  !> isw_Vphi_loc=2: value of "Vphi" is specified for given \vartheta_B postion?
  INTEGER, PUBLIC :: isw_Vphi_loc
  !> toroidal (= geometric angle) rotation frequency of species i (e.g., main ion species)
  REAL(kind=dp), PUBLIC :: Vphi
  !> only used for isw_Vphi_loc=1: (R,Z)-position of V_\varphi
  REAL(kind=dp), PUBLIC :: R_Vphi, Z_Vphi
  !> only used for isw_Vphi_loc=2: \vartheta_B postion of V_\varphi
  REAL(kind=dp), PUBLIC :: boozer_theta_Vphi
  ! Radial derivatives (w.r.t. boozer_s) of plasma parameter profiles (for all species)
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: dn_spec_ov_ds ! species density [cm^-3]
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: dT_spec_ov_ds ! species temperature [erg]
  !
  ! ADDITIONAL OUTPUT FOR MULTI-SPECIES COMPUTATIONS
  !
  !> species toroidal Mach numbers over R_major (= Mt/R)
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: MtOvR_spec
  !> species hatOmegaB_ref (= rho_L_loc*B)
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: B_rho_L_loc_spec 
  !> species poloidal variation of Vtor: poloidal resolution
  INTEGER, PARAMETER :: num_thtB_VphiProf=101
  !
  ! OUTPUT
  !> value of the average ripple of the perturbation field
  REAL(kind=dp), PUBLIC :: eps_M_2_val
  !> value of the flux surface average of $g_{\varphi\varphi}$
  !> for symmetry flux coordinates
  REAL(kind=dp), PUBLIC :: av_gphph_val
  !> value of the flux surface average of $\frac{1}{B}$
  REAL(kind=dp), PUBLIC :: av_inv_bhat_val
  !
  ! LOCAL DEFINITIONS
  !! Modification by Andreas F. Martitsch (23.08.2015)
  ! NEO-2 can treat now multiple species -> qflux is now a 4D array
  ! (at the moment these arrays cannot be handled correctly using the
  ! propagator structure -> global variables used):
  ! storage array for qflux_symm
  REAL(kind=dp), DIMENSION(:,:,:,:),  ALLOCATABLE, PUBLIC :: qflux_symm_allspec
  ! storage array for qflux_ntv
  REAL(kind=dp), DIMENSION(:,:,:,:),  ALLOCATABLE, PUBLIC :: qflux_ntv_allspec
  !! End Modification by Andreas F. Martitsch (23.08.2015)
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
  ! Local copy of y-vector (see definition in rhs_kin.f90)
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE, PUBLIC :: y_ntv_mod
  !
  !
  PUBLIC write_ntv_output
  PRIVATE write_ntv_output_a
  INTERFACE write_ntv_output
     MODULE PROCEDURE write_ntv_output_a
  END INTERFACE write_ntv_output
  !
  PUBLIC write_multispec_output
  PRIVATE write_multispec_output_a
  INTERFACE write_multispec_output
     MODULE PROCEDURE write_multispec_output_a
  END INTERFACE write_multispec_output
  !
  !
  PUBLIC compute_Dij_norm
  PRIVATE compute_Dij_norm_a
  INTERFACE compute_Dij_norm
     MODULE PROCEDURE compute_Dij_norm_a
  END INTERFACE compute_Dij_norm
  !
  PUBLIC compute_Dijab
  PRIVATE compute_Dijab_a
  INTERFACE compute_Dijab
     MODULE PROCEDURE compute_Dijab_a
  END INTERFACE compute_Dijab
  !
  PUBLIC compute_A3norm
  PRIVATE compute_A3norm_a
  INTERFACE compute_A3norm
     MODULE PROCEDURE compute_A3norm_a
  END INTERFACE compute_A3norm
  !
  PUBLIC compute_Er_and_A3norm
  PRIVATE compute_Er_and_A3norm_a
  INTERFACE compute_Er_and_A3norm
     MODULE PROCEDURE compute_Er_and_A3norm_a
  END INTERFACE compute_Er_and_A3norm
  !
  PUBLIC get_Er
  PRIVATE get_Er_a, get_Er_b
  INTERFACE get_Er
     MODULE PROCEDURE get_Er_a, get_Er_b
  END INTERFACE get_Er
  !
  PUBLIC get_B_rho_L_loc
  PRIVATE get_B_rho_L_loc_a
  INTERFACE get_B_rho_L_loc
     MODULE PROCEDURE get_B_rho_L_loc_a
  END INTERFACE get_B_rho_L_loc
  !
  PUBLIC compute_VthtB_and_VphiB
  PRIVATE compute_VthtB_and_VphiB_a, compute_VthtB_and_VphiB_b
  INTERFACE compute_VthtB_and_VphiB
     MODULE PROCEDURE compute_VthtB_and_VphiB_a, compute_VthtB_and_VphiB_b
  END INTERFACE compute_VthtB_and_VphiB
  !
  PUBLIC compute_Gamma
  PRIVATE compute_Gamma_a, compute_Gamma_b
  INTERFACE compute_Gamma
     MODULE PROCEDURE compute_Gamma_a, compute_Gamma_b
  END INTERFACE compute_Gamma
  !
  PUBLIC compute_Qflux
  PRIVATE compute_Qflux_a, compute_Qflux_b
  INTERFACE compute_Qflux
     MODULE PROCEDURE compute_Qflux_a, compute_Qflux_b
  END INTERFACE compute_Qflux
  !
  PUBLIC compute_ParFlow
  PRIVATE compute_ParFlow_a, compute_ParFlow_b
  INTERFACE compute_ParFlow
     MODULE PROCEDURE compute_ParFlow_a, compute_ParFlow_b
  END INTERFACE compute_ParFlow
  !
  PUBLIC compute_TphiNA
  PRIVATE compute_TphiNA_a
  INTERFACE compute_TphiNA
     MODULE PROCEDURE compute_TphiNA_a
  END INTERFACE compute_TphiNA
  !
CONTAINS
  !
  SUBROUTINE compute_Dij_norm_a(qflux_NA_in,qflux_AX_in,Dij_NA,Dij_AX)
    !
    USE neo_precision
    USE neo_control, ONLY: lab_swi
    USE device_mod, ONLY : device, surface
    USE mag_interface_mod, ONLY : mag_coordinates, &
         boozer_s, boozer_theta_beg, boozer_phi_beg
    USE partpa_mod,  ONLY : bmod0
    USE mag_sub, ONLY: mag
    use neo_magfie, only : boozer_curr_pol_hat, boozer_psi_pr_hat
    !
    ! ---------------------------------------------------------------!
    ! input:
    ! ---------------------------------------------------------------!
    REAL(kind=dp), DIMENSION(3,3), INTENT(in) :: qflux_NA_in
    REAL(kind=dp), DIMENSION(3,3), INTENT(in) :: qflux_AX_in
    ! ---------------------------------------------------------------!
    ! output:
    ! ---------------------------------------------------------------!
    REAL(kind=dp), DIMENSION(3,3), INTENT(out) :: Dij_NA
    REAL(kind=dp), DIMENSION(3,3), INTENT(out) :: Dij_AX
    ! ---------------------------------------------------------------!
    ! local definitions:
    ! ---------------------------------------------------------------!
    INTEGER, DIMENSION(3), PARAMETER :: ind_map = (/1,3,2/)
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: y
    REAL(kind=dp), DIMENSION(3) :: beta_out
    REAL(kind=dp) :: aiota_loc, rt0, avnabpsi
    ! ---------------------------------------------------------------!
    ! D11 and D12 for the non-axisymmetric problem
    ! normalized with the plateau coefficient
    REAL(kind=dp) :: D11_NA_Dpl, D12_NA_Dpl
    ! D13 for the non-axisymmetric problem
    ! normalized with analytical value of D31
    REAL(kind=dp) :: D13_NA_D31ref
    ! D21 and D22 for the non-axisymmetric problem
    ! normalized with the plateau coefficient
    REAL(kind=dp) :: D21_NA_Dpl, D22_NA_Dpl
    ! D23 for the non-axisymmetric problem
    ! normalized with analytical value of D31
    REAL(kind=dp) :: D23_NA_D31ref
    ! D31 and D32 for the non-axisymmetric problem
    ! normalized with analytical value of D31
    REAL(kind=dp) :: D31_NA_D31ref, D32_NA_D31ref
    ! D33 for the non-axisymmetric problem
    ! normalized with coefficient from Lorentz model
    REAL(kind=dp) :: D33_NA_norm
    ! ---------------------------------------------------------------!
    ! D11 and D12 for the axisymmetric problem
    ! normalized with the plateau coefficient
    REAL(kind=dp) :: D11_AX_Dpl, D12_AX_Dpl
    ! D13 for the axisymmetric problem
    ! normalized with analytical value of D31
    REAL(kind=dp) :: D13_AX_D31ref
    ! D21 and D22 for the axisymmetric problem
    ! normalized with the plateau coefficient
    REAL(kind=dp) :: D21_AX_Dpl, D22_AX_Dpl
    ! D23 for the axisymmetric problem
    ! normalized with analytical value of D31
    REAL(kind=dp) :: D23_AX_D31ref
    ! D31 and D32 for the axisymmetric problem
    ! normalized with analytical value of D31
    REAL(kind=dp) :: D31_AX_D31ref, D32_AX_D31ref
    ! D33 for the axisymmetric problem
    ! normalized with coefficient from Lorentz model
    REAL(kind=dp) :: D33_AX_norm
    ! ---------------------------------------------------------------!
    ! indices
    INTEGER :: i_p, j_p
    ! conversion factors for normalization
    REAL(kind=dp) :: fac1, fac2, fac3, fac4
    ! normalized co-variant phi-component of B,  $\sqrt{g}B^\vartheta$
    REAL(kind=dp) :: bcovar_phi_hat, sqrtg_bctrvr_tht
    ! Only used to call mag for normalizations
    ! related to D31 and D32
    REAL(kind=dp)                 :: bmod_tmp,sqrtg_tmp
    REAL(kind=dp), DIMENSION(3)   :: x_tmp,bder_tmp,hcovar_tmp,hctrvr_tmp
    REAL(kind=dp), DIMENSION(3,3) :: hcoder_tmp,hctder_tmp
    !
    ! copy y-vector (see definition in rhs_kin.f90)
    ALLOCATE(y(SIZE(y_ntv_mod,1)))
    y = y_ntv_mod
    !
    ! R0, aiota, avnabpsi, beta_out:
    rt0 = device%r0
    aiota_loc = surface%aiota
    avnabpsi = y(7) / y(6)
    beta_out = (/ y(14)/y(13)/avnabpsi, y(14)/y(13)/avnabpsi, y(13)/y(14) /)
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
       sqrtg_bctrvr_tht = avnabpsi*sqrtg_tmp*(hctrvr_tmp(3)*bmod_tmp*1.0e4_dp)
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
       IF (lab_swi .EQ. 10) THEN ! ASDEX-U (E. Strumberger)
          ! this is the same up to a minus sign resulting from the
          ! definition of sqrtg_tmp (left-handed system)
          sqrtg_bctrvr_tht = avnabpsi*sqrtg_tmp*(hctrvr_tmp(3)*bmod_tmp*1.0d4) 
       ELSE
          ! restore value of $\sqrt{g}B^\vartheta$ for
          ! symmetry flux coordinates from the quantities
          ! given in Boozer coordinates (only valid for right-handed system)
          sqrtg_bctrvr_tht = avnabpsi*aiota_loc*boozer_psi_pr_hat*(bmod0*1.0e4_dp)
       END IF               
    END IF
    !
    ! normalization factor from plateau coefficient
    fac1=16.0_dp*rt0*aiota_loc/PI
    ! normalization factor from gamma matrices
    fac2= - beta_out(1) * beta_out(1) / y(6)
    !PRINT *,fac1,fac2
    fac3 = - (2.0_dp/(bmod0*1.0e4_dp)) * beta_out(3) * beta_out(1) / y(6)
    fac4 = (3.0_dp*PI/16.0_dp) * beta_out(3) * beta_out(3) / y(6)
    !PRINT *,fac3, fac4
    !
    ! conversion according to definitions in Kasilov et al (2014) 
    IF(isw_qflux_NA .EQ. 1) THEN
       ! 1) axisymmetric and non-axisymmetric solution
       ! have been computed
       ! a) normalized diffusion coefficients for the
       ! non-axisymmetric case
       !
       ! $D_{11}^{\rm NA}$ and $D_{12}^{\rm NA}$ normalized with D_p
       i_p = ind_map(1)
       j_p = ind_map(1)
       D11_NA_Dpl=fac1*fac2*qflux_NA_in(i_p,j_p)
       i_p = ind_map(1)
       j_p = ind_map(2)
       D12_NA_Dpl=fac1*fac2*qflux_NA_in(i_p,j_p)
       ! $D_{13}^{\rm NA}$ normalized with D31_ref
       i_p = ind_map(1)
       j_p = ind_map(3)
       D13_NA_D31ref = (fac3*sqrtg_bctrvr_tht/bcovar_phi_hat)*qflux_NA_in(i_p,j_p)
       ! $D_{21}^{\rm NA}$ and $D_{22}^{\rm NA}$ normalized with D_p
       i_p = ind_map(2)
       j_p = ind_map(1)
       D21_NA_Dpl=fac1*fac2*qflux_NA_in(i_p,j_p)
       i_p = ind_map(2)
       j_p = ind_map(2)
       D22_NA_Dpl=fac1*fac2*qflux_NA_in(i_p,j_p)
       ! $D_{23}^{\rm NA}$ normalized with D31_ref
       i_p = ind_map(2)
       j_p = ind_map(3)
       D23_NA_D31ref = (fac3*sqrtg_bctrvr_tht/bcovar_phi_hat)*qflux_NA_in(i_p,j_p)
       ! $D_{31}^{\rm NA}$ and $D_{32}^{\rm NA}$ normalized with D31_ref
       i_p = ind_map(3)
       j_p = ind_map(1)
       D31_NA_D31ref = (fac3*sqrtg_bctrvr_tht/bcovar_phi_hat)*qflux_NA_in(i_p,j_p)
       i_p = ind_map(3)
       j_p = ind_map(2)
       D32_NA_D31ref = (fac3*sqrtg_bctrvr_tht/bcovar_phi_hat)*qflux_NA_in(i_p,j_p)
       ! $D_{33}^{\rm NA}$ normalized (experimental!!!)
       i_p = ind_map(3)
       j_p = ind_map(3)
       D33_NA_norm = fac4*qflux_NA_in(i_p,j_p)
       !
    ELSE
       ! 2) only axisymmteric solution has been computed
       ! a) non-axisymmteric coefficients set to zero
       !
       D11_NA_Dpl=0.0d0
       D12_NA_Dpl=0.0d0
       D13_NA_D31ref=0.0d0
       !
       D21_NA_Dpl=0.0d0
       D22_NA_Dpl=0.0d0
       D23_NA_D31ref=0.0d0
       !
       D31_NA_D31ref=0.0d0
       D32_NA_D31ref=0.0d0
       D33_NA_norm=0.0d0
       !
    END IF
    !
    ! b) normalized diffusion coefficients for the
    ! axisymmetric case
    !
    ! $D_{11}^{\rm NA}$ and $D_{12}^{\rm NA}$ normalized with D_p (experimental!)
    i_p = ind_map(1)
    j_p = ind_map(1)
    D11_AX_Dpl=fac1*fac2*qflux_AX_in(i_p,j_p)
    i_p = ind_map(1)
    j_p = ind_map(2)
    D12_AX_Dpl=fac1*fac2*qflux_AX_in(i_p,j_p)
    ! $D_{13}^{\rm NA}$ normalized with D31_ref
    i_p = ind_map(1)
    j_p = ind_map(3)
    D13_AX_D31ref = (fac3*sqrtg_bctrvr_tht/bcovar_phi_hat)*qflux_AX_in(i_p,j_p)
    ! $D_{21}^{\rm NA}$ and $D_{22}^{\rm NA}$ normalized with D_p
    i_p = ind_map(2)
    j_p = ind_map(1)
    D21_AX_Dpl=fac1*fac2*qflux_AX_in(i_p,j_p)
    i_p = ind_map(2)
    j_p = ind_map(2)
    D22_AX_Dpl=fac1*fac2*qflux_AX_in(i_p,j_p)
    ! $D_{23}^{\rm NA}$ normalized with D31_ref
    i_p = ind_map(2)
    j_p = ind_map(3)
    D23_AX_D31ref = (fac3*sqrtg_bctrvr_tht/bcovar_phi_hat)*qflux_AX_in(i_p,j_p)
    ! $D_{31}^{\rm AX}$ and $D_{32}^{\rm AX}$ normalized with D31_ref
    i_p = ind_map(3)
    j_p = ind_map(1)
    D31_AX_D31ref = (fac3*sqrtg_bctrvr_tht/bcovar_phi_hat)*qflux_AX_in(i_p,j_p)
    i_p = ind_map(3)
    j_p = ind_map(2)
    D32_AX_D31ref = (fac3*sqrtg_bctrvr_tht/bcovar_phi_hat)*qflux_AX_in(i_p,j_p)
    ! $D_{33}^{\rm AX}$ normalized (experimental!!!)
    i_p = ind_map(3)
    j_p = ind_map(3)
    D33_AX_norm = fac4*qflux_AX_in(i_p,j_p)
    !
    ! return normalized diffusion tensors
    ! a) non-axisymmteric coefficients
    Dij_NA(1,:) = (/ D11_NA_Dpl, D12_NA_Dpl, D13_NA_D31ref /)
    Dij_NA(2,:) = (/ D21_NA_Dpl, D22_NA_Dpl, D23_NA_D31ref /)
    Dij_NA(3,:) = (/ D31_NA_D31ref, D32_NA_D31ref, D33_NA_norm /)
    ! b) axisymmteric coefficients
    Dij_AX(1,:) = (/ D11_AX_Dpl, D12_AX_Dpl, D13_AX_D31ref /)
    Dij_AX(2,:) = (/ D21_AX_Dpl, D22_AX_Dpl, D23_AX_D31ref /)
    Dij_AX(3,:) = (/ D31_AX_D31ref, D32_AX_D31ref, D33_AX_norm /)
    !
  END SUBROUTINE compute_Dij_norm_a
  !
  SUBROUTINE compute_Dijab_a(qflux_ab_NA_in, qflux_ab_AX_in, &
       row_ind_spec, col_ind_spec, row_ind_ptr, col_ind_ptr, &
       Dijab_NA, Dijab_AX, Dijab_norm_NA, Dijab_norm_AX)
    !
    USE neo_precision
    USE neo_control, ONLY: lab_swi
    USE device_mod, ONLY : device, surface
    USE mag_interface_mod, ONLY : mag_coordinates, &
         boozer_s, boozer_theta_beg, boozer_phi_beg
    USE partpa_mod,  ONLY : bmod0
    USE mag_sub, ONLY: mag
    use neo_magfie, only : boozer_curr_pol_hat, boozer_psi_pr_hat
    USE collisionality_mod, ONLY : num_spec, species_tag, &
         z_spec, m_spec, n_spec, T_spec, collpar_spec
    !
    ! ---------------------------------------------------------------!
    ! input:
    ! ---------------------------------------------------------------!
    REAL(kind=dp), DIMENSION(:,:,:,:), INTENT(in) :: qflux_ab_NA_in
    REAL(kind=dp), DIMENSION(:,:,:,:), INTENT(in) :: qflux_ab_AX_in
    ! ---------------------------------------------------------------!
    ! output:
    ! ---------------------------------------------------------------!
    ! row- and column-indices (=species_tag) of
    ! diffusion tensor elements (e.g, D11, D12, ...)
    ! -> map to all species of the given profile (global)
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(inout) :: row_ind_spec, col_ind_spec
    ! local row- and column-indices of diffusion tensor elements
    ! (e.g, D11, D12, ...)
    ! -> map to all species of the selected radial point (local: 0:num_spec-1)
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(inout) :: row_ind_ptr, col_ind_ptr
    ! species diffusion coefficients:
    REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: Dijab_NA
    REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: Dijab_AX
    REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: Dijab_norm_NA
    REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: Dijab_norm_AX
    ! ---------------------------------------------------------------!
    ! local definitions:
    ! ---------------------------------------------------------------!
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: y
    REAL(kind=dp) :: aiota_loc, rt0, avnabpsi
    ! ---------------------------------------------------------------!
    INTEGER :: ispec_row, ispec_col, ispec_ctr
    REAL(kind=dp) :: ma, mb, m0
    REAL(kind=dp) :: vta, vtb, vt0
    REAL(kind=dp) :: za, zb, z0
    REAL(kind=dp) :: Dp00, D31ref00, D33L00_Zeff, rho0
    REAL(kind=dp) :: Dpab_ov_Dp00
    REAL(kind=dp) :: D31refab_ov_D31ref00
    REAL(kind=dp) :: D13refab_ov_D31ref00
    REAL(kind=dp) :: D33norm_fac_a0
    REAL(kind=dp), DIMENSION(3,3) :: Dij_AX, Dij_NA
    REAL(kind=dp), DIMENSION(:,:,:,:), ALLOCATABLE :: qflux_ab_NA_tmp
    REAL(kind=dp), DIMENSION(:,:,:,:), ALLOCATABLE :: qflux_ab_AX_tmp
    ! ---------------------------------------------------------------!
    ! normalized co-variant phi-component of B,  $\sqrt{g}B^\vartheta$
    REAL(kind=dp) :: bcovar_phi_hat, sqrtg_bctrvr_tht
    ! Only used to call mag for normalizations
    ! related to D31 and D32
    REAL(kind=dp)                 :: bmod_tmp,sqrtg_tmp
    REAL(kind=dp), DIMENSION(3)   :: x_tmp,bder_tmp,hcovar_tmp,hctrvr_tmp
    REAL(kind=dp), DIMENSION(3,3) :: hcoder_tmp,hctder_tmp
    ! ---------------------------------------------------------------!
    !
    ! copy y-vector (see definition in rhs_kin.f90)
    ALLOCATE(y(SIZE(y_ntv_mod,1)))
    y = y_ntv_mod
    !
    ! R0, aiota, avnabpsi :
    rt0 = device%r0
    aiota_loc = surface%aiota
    avnabpsi = y(7) / y(6)
    !
    ! computation of the normalization for D31 and D32 (-> D31_ref)
    IF (mag_coordinates .EQ. 0) THEN
       ! cylindrical coordinates
       x_tmp = xstart_cyl
       CALL mag(x_tmp,bmod_tmp,sqrtg_tmp,bder_tmp,hcovar_tmp,&
            hctrvr_tmp,hcoder_tmp,hctder_tmp)
       !
       ! normalized co-variant phi-component of B
       ! (Note! co-variant phi-component of B is the
       ! same for cylindrical coordinates and symmetry flux
       ! coodrinates --> no conversion needed)
       bcovar_phi_hat = hcovar_tmp(2)*(bmod_tmp/bmod0)
       !
       ! restore value of $\sqrt{g}B^\vartheta$ for
       ! symmetry flux coordinates from quantities
       ! given in cylindircal coordinates
       sqrtg_bctrvr_tht = avnabpsi*sqrtg_tmp*(hctrvr_tmp(3)*bmod_tmp*1.0e4_dp)
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
       IF (lab_swi .EQ. 10) THEN ! ASDEX-U (E. Strumberger)
          ! this is the same up to a minus sign resulting from the
          ! definition of sqrtg_tmp (left-handed system)
          sqrtg_bctrvr_tht = avnabpsi*sqrtg_tmp*(hctrvr_tmp(3)*bmod_tmp*1.0d4) 
       ELSE
          ! restore value of $\sqrt{g}B^\vartheta$ for
          ! symmetry flux coordinates from the quantities
          ! given in Boozer coordinates (only valid for right-handed system)
          sqrtg_bctrvr_tht = avnabpsi*aiota_loc*boozer_psi_pr_hat*(bmod0*1.0e4_dp)
       END IF
    END IF
    !
    ! Allocate temporary storage arrays for qflux_ab
    IF (ALLOCATED(qflux_ab_NA_tmp)) DEALLOCATE(qflux_ab_NA_tmp)
    ALLOCATE(qflux_ab_NA_tmp(1:3,1:3,0:num_spec-1,0:num_spec-1))
    qflux_ab_NA_tmp = qflux_ab_NA_in
    !
    IF (ALLOCATED(qflux_ab_AX_tmp)) DEALLOCATE(qflux_ab_AX_tmp)
    ALLOCATE(qflux_ab_AX_tmp(1:3,1:3,0:num_spec-1,0:num_spec-1))
    qflux_ab_AX_tmp = qflux_ab_AX_in
    !PRINT *,qflux_ab_AX_tmp(1,1,0,0),qflux_ab_AX_tmp(1,1,0,1)
    !PRINT *,qflux_ab_AX_tmp(1,1,1,0),qflux_ab_AX_tmp(1,1,1,1)
    !
    ! Allocate row- and column-indices (=species_tag)
    ! -> map to all species of the given profile (global)
    IF (ALLOCATED(row_ind_spec)) DEALLOCATE(row_ind_spec)
    ALLOCATE(row_ind_spec(0:num_spec**2-1))
    !
    IF (ALLOCATED(col_ind_spec)) DEALLOCATE(col_ind_spec)
    ALLOCATE(col_ind_spec(0:num_spec**2-1))
    !
    ! Allocate local row- and column-indices
    ! -> map to all species of the selected radial point (local: 0:num_spec-1)
    IF (ALLOCATED(row_ind_ptr)) DEALLOCATE(row_ind_ptr)
    ALLOCATE(row_ind_ptr(0:num_spec**2-1))
    !
    IF (ALLOCATED(col_ind_ptr)) DEALLOCATE(col_ind_ptr)
    ALLOCATE(col_ind_ptr(0:num_spec**2-1))
    !
    ! Allocate species diffusion coefficients
    IF (ALLOCATED(Dijab_NA)) DEALLOCATE(Dijab_NA)
    ALLOCATE(Dijab_NA(1:3,1:3,0:num_spec**2-1))
    !
    IF (ALLOCATED(Dijab_AX)) DEALLOCATE(Dijab_AX)
    ALLOCATE(Dijab_AX(1:3,1:3,0:num_spec**2-1))
    !
    IF (ALLOCATED(Dijab_norm_NA)) DEALLOCATE(Dijab_norm_NA)
    ALLOCATE(Dijab_norm_NA(1:3,1:3,0:num_spec**2-1))
    !
    IF (ALLOCATED(Dijab_norm_AX)) DEALLOCATE(Dijab_norm_AX)
    ALLOCATE(Dijab_norm_AX(1:3,1:3,0:num_spec**2-1))
    !
    ! reference mass, thermal velocity and charge number (=electrons)
    m0 = m_spec(0)
    vt0 = SQRT(2*T_spec(0)/m0)
    z0 = z_spec(0)
    rho0 = vt0*m0*c/(z0*e*bmod0*1.0e4_dp) ! [cm] (cgs-units)
    !
    ! Compute D_p, D_{31,ref} and D_{33]^L for 0-species (=electrons)
    Dp00 = PI*vt0*(rho0**2)/(16.0_dp*aiota_loc*rt0)
    ! Modification by Andreas F. Martitsch (27.06.2017)
    !-> old: Expression given below is only valid for
    !-> right-handed coordinate system!!!
    !-> However, normalizations used in 'compute_Dij_norm_a' are
    !-> generally valid and, therfore, this is inconsistent. 
    !D31ref00 = vt0 * rho0 * bcovar_phi_hat * (bmod0*1.0e4_dp) / &
    !     (2.0_dp * aiota_loc * boozer_psi_pr_hat * avnabpsi)
    !-> new:
    D31ref00 = vt0 * rho0 * bcovar_phi_hat * ((bmod0*1.0e4_dp)**2) / &
         (2.0_dp * sqrtg_bctrvr_tht)
    ! End Modification by Andreas F. Martitsch (27.06.2017)
    D33L00_Zeff = (-16.0_dp/(3.0_dp*PI)) * ((bmod0*1.0e4_dp)**2) * &
         vt0 * (2.0_dp/collpar_spec(0))
    !
    ! Computation of normalized diffusion coefficients
    ! (w.r.t. 0-species (=electrons))
    !
    DO ispec_row = 0,num_spec-1
       !
       ma = m_spec(ispec_row)
       vta = SQRT(2*T_spec(ispec_row)/ma)
       za = z_spec(ispec_row)
       !
       DO ispec_col = 0,num_spec-1
          !
          ispec_ctr = ispec_row * num_spec + ispec_col
          !
          CALL compute_Dij_norm_a(qflux_ab_NA_tmp(:,:,ispec_row,ispec_col),&
               qflux_ab_AX_tmp(:,:,ispec_row,ispec_col), Dij_NA, Dij_AX)
          !
          ! re-normalize w.r.t. 0-species (=electrons)
          mb = m_spec(ispec_col)
          vtb = SQRT(2*T_spec(ispec_col)/mb)
          zb = z_spec(ispec_col)
          !
          Dpab_ov_Dp00 = (ma*mb/(m0**2)) * ((vta**2)*vtb/(vt0**3)) * &
               ((z0**2)/(za*zb))
          D31refab_ov_D31ref00 = (mb/m0) * (vta*vtb/(vt0**2)) * (z0/zb)
          D13refab_ov_D31ref00 = (ma/m0) * ((vta**2)/(vt0**2)) * (z0/za)
          D33norm_fac_a0 = (collpar_spec(0)/2.0_dp)*(vta/vt0)
          !
          ! row- and column-indices (=species) of diffusion coefficients
          ! -> map to all species of the given profile (global)
          row_ind_spec(ispec_ctr) = species_tag(ispec_row)
          col_ind_spec(ispec_ctr) = species_tag(ispec_col)
          !
          ! local row- and column-indices of diffusion coefficients
          ! -> map to all species of the selected radial point (0:num_spec-1)
          row_ind_ptr(ispec_ctr) = ispec_row
          col_ind_ptr(ispec_ctr) = ispec_col
          !
          ! non-axisymmetric coefficients (normalized w.r.t. 0-species)
          Dijab_norm_NA(1,1,ispec_ctr) = Dij_NA(1,1)*Dpab_ov_Dp00
          Dijab_norm_NA(1,2,ispec_ctr) = Dij_NA(1,2)*Dpab_ov_Dp00
          Dijab_norm_NA(1,3,ispec_ctr) = Dij_NA(1,3)*D13refab_ov_D31ref00
          !
          Dijab_norm_NA(2,1,ispec_ctr) = Dij_NA(2,1)*Dpab_ov_Dp00
          Dijab_norm_NA(2,2,ispec_ctr) = Dij_NA(2,2)*Dpab_ov_Dp00
          Dijab_norm_NA(2,3,ispec_ctr) = Dij_NA(2,3)*D13refab_ov_D31ref00
          !
          Dijab_norm_NA(3,1,ispec_ctr) = Dij_NA(3,1)*D31refab_ov_D31ref00
          Dijab_norm_NA(3,2,ispec_ctr) = Dij_NA(3,2)*D31refab_ov_D31ref00
          Dijab_norm_NA(3,3,ispec_ctr) = Dij_NA(3,3)*D33norm_fac_a0
          !
          ! axisymmetric coefficients (normalized w.r.t. 0-species)
          Dijab_norm_AX(1,1,ispec_ctr) = Dij_AX(1,1)*Dpab_ov_Dp00
          Dijab_norm_AX(1,2,ispec_ctr) = Dij_AX(1,2)*Dpab_ov_Dp00
          Dijab_norm_AX(1,3,ispec_ctr) = Dij_AX(1,3)*D13refab_ov_D31ref00
          !
          Dijab_norm_AX(2,1,ispec_ctr) = Dij_AX(2,1)*Dpab_ov_Dp00
          Dijab_norm_AX(2,2,ispec_ctr) = Dij_AX(2,2)*Dpab_ov_Dp00
          Dijab_norm_AX(2,3,ispec_ctr) = Dij_AX(2,3)*D13refab_ov_D31ref00
          !
          Dijab_norm_AX(3,1,ispec_ctr) = Dij_AX(3,1)*D31refab_ov_D31ref00
          Dijab_norm_AX(3,2,ispec_ctr) = Dij_AX(3,2)*D31refab_ov_D31ref00
          Dijab_norm_AX(3,3,ispec_ctr) = Dij_AX(3,3)*D33norm_fac_a0
          !
          ! non-axisymmetric coefficients (dimensional, cgs-units)
          Dijab_NA(1,1,ispec_ctr) = Dijab_norm_NA(1,1,ispec_ctr) * Dp00
          Dijab_NA(1,2,ispec_ctr) = Dijab_norm_NA(1,2,ispec_ctr) * Dp00
          Dijab_NA(1,3,ispec_ctr) = Dijab_norm_NA(1,3,ispec_ctr) * D31ref00
          !
          Dijab_NA(2,1,ispec_ctr) = Dijab_norm_NA(2,1,ispec_ctr) * Dp00
          Dijab_NA(2,2,ispec_ctr) = Dijab_norm_NA(2,2,ispec_ctr) * Dp00
          Dijab_NA(2,3,ispec_ctr) = Dijab_norm_NA(2,3,ispec_ctr) * D31ref00
          !
          Dijab_NA(3,1,ispec_ctr) = Dijab_norm_NA(3,1,ispec_ctr) * D31ref00
          Dijab_NA(3,2,ispec_ctr) = Dijab_norm_NA(3,2,ispec_ctr) * D31ref00
          Dijab_NA(3,3,ispec_ctr) = Dijab_norm_NA(3,3,ispec_ctr) * D33L00_Zeff
          !
          ! axisymmetric coefficients (dimensional, cgs-units)
          Dijab_AX(1,1,ispec_ctr) = Dijab_norm_AX(1,1,ispec_ctr) * Dp00
          Dijab_AX(1,2,ispec_ctr) = Dijab_norm_AX(1,2,ispec_ctr) * Dp00
          Dijab_AX(1,3,ispec_ctr) = Dijab_norm_AX(1,3,ispec_ctr) * D31ref00
          !
          Dijab_AX(2,1,ispec_ctr) = Dijab_norm_AX(2,1,ispec_ctr) * Dp00
          Dijab_AX(2,2,ispec_ctr) = Dijab_norm_AX(2,2,ispec_ctr) * Dp00
          Dijab_AX(2,3,ispec_ctr) = Dijab_norm_AX(2,3,ispec_ctr) * D31ref00
          !
          Dijab_AX(3,1,ispec_ctr) = Dijab_norm_AX(3,1,ispec_ctr) * D31ref00
          Dijab_AX(3,2,ispec_ctr) = Dijab_norm_AX(3,2,ispec_ctr) * D31ref00
          Dijab_AX(3,3,ispec_ctr) = Dijab_norm_AX(3,3,ispec_ctr) * D33L00_Zeff
          !
       END DO
       !
    END DO
    !
  END SUBROUTINE compute_Dijab_a
  !
  SUBROUTINE write_ntv_output_a(qflux_in)
    !
    USE neo_precision
    USE neo_control, ONLY: lab_swi
    USE device_mod, ONLY : device, surface
    USE mag_interface_mod, ONLY : mag_coordinates, &
         boozer_s, boozer_theta_beg, boozer_phi_beg
    USE partpa_mod,  ONLY : bmod0
    USE mag_sub, ONLY: mag
    use neo_magfie, only : boozer_curr_pol_hat, boozer_psi_pr_hat
    USE collisionality_mod, ONLY : collpar
    !
    ! ---------------------------------------------------------------!
    ! input:
    ! ---------------------------------------------------------------!
    REAL(kind=dp), DIMENSION(:,:), INTENT(in) :: qflux_in
    ! ---------------------------------------------------------------!
    ! local definitions:
    ! ---------------------------------------------------------------!
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: y
    REAL(kind=dp) :: aiota_loc, rt0, avnabpsi
    ! ---------------------------------------------------------------!
    REAL(kind=dp), DIMENSION(3,3) :: Dij_AX, Dij_NA, dummy_33
    ! D11 and D12 for the non-axisymmetric problem
    ! normalized with the plateau coefficient
    REAL(kind=dp) :: D11_NA_Dpl, D12_NA_Dpl    
    ! D13 for the non-axisymmetric problem
    ! normalized with analytical value of D31
    REAL(kind=dp) :: D13_NA_D31ref
    ! ---------------------------------------------------------------!
    ! D11 and D12 for the axisymmetric problem
    ! normalized with the plateau coefficient
    REAL(kind=dp) :: D11_AX_Dpl, D12_AX_Dpl
    ! D31 and D32 for the axisymmetric problem
    ! normalized with analytical value of D31
    REAL(kind=dp) :: D31_AX_D31ref, D32_AX_D31ref, k_cof
    ! ---------------------------------------------------------------!
    ! indices, file id
    LOGICAL :: opened
    INTEGER :: i_p, j_p, uw    
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
    REAL(kind=dp) :: bcovar_phi, bcovar_tht, avbhat2, avb2, avbhat
    CHARACTER(len=30) :: file_name
    !
    ! copy y-vector (see definition in rhs_kin.f90)
    ALLOCATE(y(SIZE(y_ntv_mod,1)))
    y = y_ntv_mod
    !
    ! R0, aiota, avnabpsi :
    rt0 = device%r0
    aiota_loc = surface%aiota
    avnabpsi = y(7) / y(6)
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
       sqrtg_bctrvr_tht = avnabpsi*sqrtg_tmp*(hctrvr_tmp(3)*bmod_tmp*1.0e4_dp)
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
       IF (lab_swi .EQ. 10) THEN ! ASDEX-U (E. Strumberger)
          ! this is the same up to a minus sign resulting from the
          ! definition of sqrtg_tmp (left-handed system)
          sqrtg_bctrvr_tht = avnabpsi*sqrtg_tmp*(hctrvr_tmp(3)*bmod_tmp*1.0d4) 
       ELSE
          ! restore value of $\sqrt{g}B^\vartheta$ for
          ! symmetry flux coordinates from the quantities
          ! given in Boozer coordinates (only valid for right-handed system)
          sqrtg_bctrvr_tht = avnabpsi*aiota_loc*boozer_psi_pr_hat*(bmod0*1.0e4_dp)
       END IF
    END IF
    !
    ! extra output for NTV computations
    IF(isw_qflux_NA .EQ. 1) THEN
       IF (ALLOCATED(qflux_symm)) THEN
          !
          CALL compute_Dij_norm_a(qflux_in, qflux_symm, Dij_NA, Dij_AX)
          !
          ! non-axisymmetric coefficients
          D11_NA_Dpl = Dij_NA(1,1)
          D12_NA_Dpl = Dij_NA(1,2)
          D13_NA_D31ref = Dij_NA(1,3)
          !
          ! axisymmetric solution
          D11_AX_Dpl = Dij_AX(1,1)
          D12_AX_Dpl = Dij_AX(1,2)
          D31_AX_D31ref = Dij_AX(3,1)
          D32_AX_D31ref = Dij_AX(3,2)
          k_cof = (2.5_dp-D32_AX_D31ref/D31_AX_D31ref)
          !
       ELSE
          !
          dummy_33 = 0.0_dp
          CALL compute_Dij_norm_a(qflux_in, dummy_33, Dij_NA, Dij_AX)
          !
          ! non-axisymmetric coefficients
          D11_NA_Dpl = Dij_NA(1,1)
          D12_NA_Dpl = Dij_NA(1,2)
          D13_NA_D31ref = Dij_NA(1,3)
          !
          ! in case of bounce-averaged model the solution
          ! for the axisymmetric problem is not computed
          D11_AX_Dpl = 0.0_dp
          D12_AX_Dpl = 0.0_dp
          D31_AX_D31ref = 0.0_dp
          D32_AX_D31ref = 0.0_dp
          k_cof = 0.0_dp
          !
       END IF
       !
    ELSE
       !
       ! in this case "qflux_symm" is stored within the actual propagator "prop_a"
       dummy_33 = 0.0_dp
       CALL compute_Dij_norm_a(dummy_33, qflux_in, Dij_NA, Dij_AX)
       !
       ! only axisymmteric solution has been computed
       ! (non-axisymmteric coefficients set to zero)
       D11_NA_Dpl=0.0_dp
       D12_NA_Dpl=0.0_dp
       D13_NA_D31ref=0.0_dp
       !
       ! axisymmetric solution
       D11_AX_Dpl = Dij_AX(1,1)
       D12_AX_Dpl = Dij_AX(1,2)
       D31_AX_D31ref = Dij_AX(3,1)
       D32_AX_D31ref = Dij_AX(3,2)
       k_cof = (2.5_dp-D32_AX_D31ref/D31_AX_D31ref)
       !
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
    ! physical output
    Mt_val=MtOvR*rt0
    nu_star=collpar*rt0/aiota_loc
    sqrtg_bctrvr_phi=sqrtg_bctrvr_tht/aiota_loc
    bcovar_phi=hcovar_tmp(2)*(bmod_tmp*1.0e4_dp)
    bcovar_tht=hcovar_tmp(3)*(bmod_tmp*1.0e4_dp)
    avbhat2=y(9)/y(6)
    avb2=avbhat2*((bmod0*1.0e4_dp)**2)
    avbhat=y(14)/y(13)
    !
    !PRINT *,'avbhat: ',avbhat
    !PRINT *,'avnabpsi: ',avnabpsi
    !PRINT *,'$\int \rd s / B$: ',y(6)
    !PRINT *,'R: ',rt0
    !PRINT *,'q: ',1/aiota_loc
    !
    file_name='ntv_out.dat'
    OPEN(uw,file=TRIM(ADJUSTL(file_name)),status='replace')
    WRITE (uw,'(1000(1x,e18.5))') &
         boozer_s, Mt_val, nu_star, B_rho_L_loc, &
         D31_AX_D31ref, D32_AX_D31ref, k_cof, &
         D11_NA_Dpl, D12_NA_Dpl, D13_NA_D31ref, &
         aiota_loc, rt0, (bmod0*1.0e4_dp), &
         boozer_psi_pr_hat, avnabpsi, &
         sqrtg_bctrvr_tht, sqrtg_bctrvr_phi, bcovar_tht, bcovar_phi, &
         DBLE(m_phi), avbhat2, av_inv_bhat_val, eps_M_2_val, &
         av_gphph_val, avbhat, D11_AX_Dpl, D12_AX_Dpl
    CLOSE(uw)
    !
  END SUBROUTINE write_ntv_output_a

  !> \brief Write multispecies output to hdf5 file.
  !>
  !> This subroutine makes sure multispecies output quantities are
  !> computed, and writes them to a hdf5 file.
  !>
  !> input:
  !> ------
  !> none
  !>
  !> output:
  !> -------
  !> none
  !>
  !> side effects:
  !> -------------
  !> - writes file to disk in current folder
  !> - (?)
  SUBROUTINE write_multispec_output_a()
    !
    USE neo_precision
    USE neo_control, ONLY: lab_swi
    USE device_mod, ONLY : device, surface
    USE mag_interface_mod, ONLY : mag_coordinates, &
         boozer_s, boozer_theta_beg, boozer_phi_beg
    USE partpa_mod,  ONLY : bmod0
    USE mag_sub, ONLY: mag
    use neo_magfie, only : boozer_curr_pol_hat, boozer_psi_pr_hat, &
         boozer_curr_tor_hat_s, boozer_curr_pol_hat_s, &
      & boozer_iota_s
    USE collisionality_mod, ONLY : collpar, num_spec, species_tag, &
         z_spec, m_spec, n_spec, T_spec, collpar_spec, isw_coul_log
    ! Output stored as HDF5-file
    USE hdf5_tools
    !
    ! ---------------------------------------------------------------!
    ! local definitions:
    ! ---------------------------------------------------------------!
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: y
    REAL(kind=dp) :: aiota_loc, rt0, avnabpsi
    ! ---------------------------------------------------------------!
    INTEGER :: ispec_row, ispec_col, ispec_ctr
    REAL(kind=dp) :: ma, mb, m0
    REAL(kind=dp) :: vta, vtb, vt0
    REAL(kind=dp) :: za, zb, z0
    REAL(kind=dp) :: Dp00, D31ref00, D33L00_Zeff, rho0
    REAL(kind=dp) :: Dpab_ov_Dp00
    REAL(kind=dp) :: D31refab_ov_D31ref00
    REAL(kind=dp) :: D13refab_ov_D31ref00
    REAL(kind=dp) :: D33norm_fac_a0
    REAL(kind=dp), DIMENSION(3,3) :: Dij_AX, Dij_NA
    ! ---------------------------------------------------------------!
    ! row- and column-indices (=species_tag) of
    ! diffusion tensor elements (e.g, D11, D12, ...)
    ! -> map to all species of the given profile (global)
    INTEGER, DIMENSION(:), ALLOCATABLE :: row_ind_spec, col_ind_spec
    ! local row- and column-indices of diffusion tensor elements
    ! (e.g, D11, D12, ...)
    ! -> map to all species of the selected radial point (local: 0:num_spec-1)
    INTEGER, DIMENSION(:), ALLOCATABLE :: row_ind_ptr, col_ind_ptr
    ! ---------------------------------------------------------------!
    ! D11 and D12 for the non-axisymmetric problem
    ! normalized with the plateau coefficient
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: D11_NA_Dpl, D12_NA_Dpl    
    ! D13 for the non-axisymmetric problem
    ! normalized with analytical value of D31
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: D13_NA_D31ref
    ! D21 and D22 for the non-axisymmetric problem
    ! normalized with the plateau coefficient
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: D21_NA_Dpl, D22_NA_Dpl
    ! D23 for the non-axisymmetric problem
    ! normalized with analytical value of D31
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: D23_NA_D31ref
    ! D31 and D32 for the non-axisymmetric problem
    ! normalized with analytical value of D31
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: D31_NA_D31ref, D32_NA_D31ref
    ! D33 for the non-axisymmetric problem
    ! normalized with coefficient from Lorentz model
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: D33_NA_norm
    ! ---------------------------------------------------------------!
    ! D11 and D12 for the non-axisymmetric problem (cgs-units)
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: D11_NA, D12_NA   
    ! D13 for the non-axisymmetric problem (cgs-units)
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: D13_NA
    ! D21 and D22 for the non-axisymmetric problem (cgs-units)
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: D21_NA, D22_NA
    ! D23 for the non-axisymmetric problem (cgs-units)
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: D23_NA
    ! D31 and D32 for the non-axisymmetric problem (cgs-units)
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: D31_NA, D32_NA
    ! D33 for the non-axisymmetric problem (cgs-units)
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: D33_NA
    ! ---------------------------------------------------------------!
    ! D11 and D12 for the axisymmetric problem
    ! normalized with the plateau coefficient
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: D11_AX_Dpl, D12_AX_Dpl
    ! D13 for the axisymmetric problem
    ! normalized with analytical value of D31
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: D13_AX_D31ref
    ! D21 and D22 for the axisymmetric problem
    ! normalized with the plateau coefficient
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: D21_AX_Dpl, D22_AX_Dpl
    ! D23 for the axisymmetric problem
    ! normalized with analytical value of D31
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: D23_AX_D31ref
    ! D31 and D32 for the axisymmetric problem
    ! normalized with analytical value of D31
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: D31_AX_D31ref, D32_AX_D31ref
    ! D33 for the axisymmetric problem
    ! normalized with coefficient from Lorentz model
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: D33_AX_norm
    ! ---------------------------------------------------------------!
    ! D11 and D12 for the axisymmetric problem (cgs-units)
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: D11_AX, D12_AX   
    ! D13 for the axisymmetric problem (cgs-units)
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: D13_AX
    ! D21 and D22 for the axisymmetric problem (cgs-units)
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: D21_AX, D22_AX
    ! D23 for the axisymmetric problem (cgs-units)
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: D23_AX
    ! D31 and D32 for the axisymmetric problem (cgs-units)
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: D31_AX, D32_AX
    ! D33 for the axisymmetric problem (cgs-units)
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: D33_AX
    ! ---------------------------------------------------------------!
    REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: Dijab_NA
    REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: Dijab_AX
    REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: Dijab_norm_NA
    REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: Dijab_norm_AX
    REAL(kind=dp), DIMENSION(:,:,:,:), ALLOCATABLE :: dummy_qflux
    ! ---------------------------------------------------------------!
    ! species poloidal and toroidal rotation velocities
    ! (defined w.r.t Boozer angles + flux surface averaged):
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: VthtB_spec
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: VphiB_spec
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: VthtB_Ware_spec
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: VphiB_Ware_spec
    ! species poloidal variation of Vtor and Vpol
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: R_Vphi_prof, Z_Vphi_prof
    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: Vphi_prof_spec
    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: Vtht_prof_spec
    ! species poloidal variation of Vtor and Vpol without E_{\parallel}-contribution
    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: Vphi_prof_woWare_spec
    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: Vtht_prof_woWare_spec
    ! ---------------------------------------------------------------!
    ! species particle flux densities (AX + NA) and NTV torque densities
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: Gamma_AX_spec, Gamma_AX_Ware_spec
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: Gamma_NA_spec, Gamma_NA_Ware_spec
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: TphiNA_spec, TphiNA_Ware_spec
    REAL(kind=dp) :: TphiNA_tot, TphiNA_Ware_tot
    ! species heat flux densities (AX + NA)
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: Qflux_AX_spec, Qflux_AX_Ware_spec
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: Qflux_NA_spec, Qflux_NA_Ware_spec
    ! species parallel flow (AX + NA)
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: ParFlow_AX_spec, ParFlow_AX_Ware_spec
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: ParFlow_NA_spec, ParFlow_NA_Ware_spec
    ! ---------------------------------------------------------------!
    ! radial electric field Er, <E_par*B>/<B^2>
    REAL(kind=dp) :: Er, avEparB_ov_avb2
    ! normalized co-variant phi-component of B,  $\sqrt{g}B^\vartheta$
    REAL(kind=dp) :: bcovar_phi_hat, sqrtg_bctrvr_tht
    ! Only used to call mag for normalizations
    ! related to D31 and D32
    REAL(kind=dp)                 :: bmod_tmp,sqrtg_tmp
    REAL(kind=dp), DIMENSION(3)   :: x_tmp,bder_tmp,hcovar_tmp,hctrvr_tmp
    REAL(kind=dp), DIMENSION(3,3) :: hcoder_tmp,hctder_tmp
    ! Physical output (collisionality, $\sqrt{g}B^\varphi$)
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: nu_star_spec
    REAL(kind=dp) :: sqrtg_bctrvr_phi
    ! Physical output ($B_\varphi$,$B_\vartheta$,\langle{B^2}\rangle)
    REAL(kind=dp) :: bcovar_phi, bcovar_tht, avbhat2, avb2, avbhat
    REAL(kind=dp) :: dbcovar_theta_hat_ds, dbcovar_phi_hat_ds
    ! HDF5 file id
    INTEGER(HID_T)   :: h5id_multispec
    ! HDF5 file name
    CHARACTER(len=30) :: file_name
    !
    ! copy y-vector (see definition in rhs_kin.f90)
    ALLOCATE(y(SIZE(y_ntv_mod,1)))
    y = y_ntv_mod
    !
    ! R0, aiota, avnabpsi :
    rt0 = device%r0
    aiota_loc = surface%aiota
    avnabpsi = y(7) / y(6)
    !
    ! reference mass, thermal velocity and charge number (=electrons)
    m0 = m_spec(0)
    vt0 = SQRT(2*T_spec(0)/m0)
    z0 = z_spec(0)
    rho0 = vt0*m0*c/(z0*e*bmod0*1.0e4_dp) ! [cm] (cgs-units)
    !
    ! computation of the normalization for D31 and D32 (-> D31_ref)
    IF (mag_coordinates .EQ. 0) THEN
       ! cylindrical coordinates
       x_tmp = xstart_cyl
       CALL mag(x_tmp,bmod_tmp,sqrtg_tmp,bder_tmp,hcovar_tmp,&
            hctrvr_tmp,hcoder_tmp,hctder_tmp)
       !
       ! normalized co-variant phi-component of B
       ! (Note! co-variant phi-component of B is the
       ! same for cylindrical coordinates and symmetry flux
       ! coodrinates --> no conversion needed)
       bcovar_phi_hat = hcovar_tmp(2)*(bmod_tmp/bmod0)
       !
       dbcovar_theta_hat_ds = 0.0_dp ! not available at the moment !!!
       dbcovar_phi_hat_ds = 0.0_dp ! not available at the moment !!!
       !
       ! restore value of $\sqrt{g}B^\vartheta$ for
       ! symmetry flux coordinates from quantities
       ! given in cylindircal coordinates
       sqrtg_bctrvr_tht = avnabpsi*sqrtg_tmp*(hctrvr_tmp(3)*bmod_tmp*1.0e4_dp)
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
       dbcovar_theta_hat_ds = boozer_curr_tor_hat_s
       dbcovar_phi_hat_ds = boozer_curr_pol_hat_s
       !
       IF (lab_swi .EQ. 10) THEN ! ASDEX-U (E. Strumberger)
          ! this is the same up to a minus sign resulting from the
          ! definition of sqrtg_tmp (left-handed system)
          sqrtg_bctrvr_tht = avnabpsi*sqrtg_tmp*(hctrvr_tmp(3)*bmod_tmp*1.0d4) 
       ELSE
          ! restore value of $\sqrt{g}B^\vartheta$ for
          ! symmetry flux coordinates from the quantities
          ! given in Boozer coordinates (only valid for right-handed system)
          sqrtg_bctrvr_tht = avnabpsi*aiota_loc*boozer_psi_pr_hat*(bmod0*1.0e4_dp)
       END IF
    END IF
    !
    ! Allocate row- and column-indices (=species_tag)
    ! -> map to all species of the given profile (global)
    IF (ALLOCATED(row_ind_spec)) DEALLOCATE(row_ind_spec)
    ALLOCATE(row_ind_spec(0:num_spec**2-1))
    !
    IF (ALLOCATED(col_ind_spec)) DEALLOCATE(col_ind_spec)
    ALLOCATE(col_ind_spec(0:num_spec**2-1))
    !
    ! Allocate local row- and column-indices
    ! -> map to all species of the selected radial point (local: 0:num_spec-1)
    IF (ALLOCATED(row_ind_ptr)) DEALLOCATE(row_ind_ptr)
    ALLOCATE(row_ind_ptr(0:num_spec**2-1))
    !
    IF (ALLOCATED(col_ind_ptr)) DEALLOCATE(col_ind_ptr)
    ALLOCATE(col_ind_ptr(0:num_spec**2-1))
    !
    ! Allocate storage arrays (axisymmetric solution)
    IF (ALLOCATED(D11_AX_Dpl)) DEALLOCATE(D11_AX_Dpl)
    ALLOCATE(D11_AX_Dpl(0:num_spec**2-1))
    IF (ALLOCATED(D12_AX_Dpl)) DEALLOCATE(D12_AX_Dpl)
    ALLOCATE(D12_AX_Dpl(0:num_spec**2-1))
    IF (ALLOCATED(D13_AX_D31ref)) DEALLOCATE(D13_AX_D31ref)
    ALLOCATE(D13_AX_D31ref(0:num_spec**2-1))
    !
    IF (ALLOCATED(D21_AX_Dpl)) DEALLOCATE(D21_AX_Dpl)
    ALLOCATE(D21_AX_Dpl(0:num_spec**2-1))
    IF (ALLOCATED(D22_AX_Dpl)) DEALLOCATE(D22_AX_Dpl)
    ALLOCATE(D22_AX_Dpl(0:num_spec**2-1))
    IF (ALLOCATED(D23_AX_D31ref)) DEALLOCATE(D23_AX_D31ref)
    ALLOCATE(D23_AX_D31ref(0:num_spec**2-1))
    !
    IF (ALLOCATED(D31_AX_D31ref)) DEALLOCATE(D31_AX_D31ref)
    ALLOCATE(D31_AX_D31ref(0:num_spec**2-1))
    IF (ALLOCATED(D32_AX_D31ref)) DEALLOCATE(D32_AX_D31ref)
    ALLOCATE(D32_AX_D31ref(0:num_spec**2-1))
    IF (ALLOCATED(D33_AX_norm)) DEALLOCATE(D33_AX_norm)
    ALLOCATE(D33_AX_norm(0:num_spec**2-1))
    !
    ! Allocate storage arrays (non-axisymmetric solution)
    IF (ALLOCATED(D11_NA_Dpl)) DEALLOCATE(D11_NA_Dpl)
    ALLOCATE(D11_NA_Dpl(0:num_spec**2-1))
    IF (ALLOCATED(D12_NA_Dpl)) DEALLOCATE(D12_NA_Dpl)
    ALLOCATE(D12_NA_Dpl(0:num_spec**2-1))
    IF (ALLOCATED(D13_NA_D31ref)) DEALLOCATE(D13_NA_D31ref)
    ALLOCATE(D13_NA_D31ref(0:num_spec**2-1))
    !
    IF (ALLOCATED(D21_NA_Dpl)) DEALLOCATE(D21_NA_Dpl)
    ALLOCATE(D21_NA_Dpl(0:num_spec**2-1))
    IF (ALLOCATED(D22_NA_Dpl)) DEALLOCATE(D22_NA_Dpl)
    ALLOCATE(D22_NA_Dpl(0:num_spec**2-1))
    IF (ALLOCATED(D23_NA_D31ref)) DEALLOCATE(D23_NA_D31ref)
    ALLOCATE(D23_NA_D31ref(0:num_spec**2-1))
    !
    IF (ALLOCATED(D31_NA_D31ref)) DEALLOCATE(D31_NA_D31ref)
    ALLOCATE(D31_NA_D31ref(0:num_spec**2-1))
    IF (ALLOCATED(D32_NA_D31ref)) DEALLOCATE(D32_NA_D31ref)
    ALLOCATE(D32_NA_D31ref(0:num_spec**2-1))
    IF (ALLOCATED(D33_NA_norm)) DEALLOCATE(D33_NA_norm)
    ALLOCATE(D33_NA_norm(0:num_spec**2-1))
    !
    ! Allocate storage arrays (axisymmetric solution, cgs-units)
    IF (ALLOCATED(D11_AX)) DEALLOCATE(D11_AX)
    ALLOCATE(D11_AX(0:num_spec**2-1))
    IF (ALLOCATED(D12_AX)) DEALLOCATE(D12_AX)
    ALLOCATE(D12_AX(0:num_spec**2-1))
    IF (ALLOCATED(D13_AX)) DEALLOCATE(D13_AX)
    ALLOCATE(D13_AX(0:num_spec**2-1))
    !
    IF (ALLOCATED(D21_AX)) DEALLOCATE(D21_AX)
    ALLOCATE(D21_AX(0:num_spec**2-1))
    IF (ALLOCATED(D22_AX)) DEALLOCATE(D22_AX)
    ALLOCATE(D22_AX(0:num_spec**2-1))
    IF (ALLOCATED(D23_AX)) DEALLOCATE(D23_AX)
    ALLOCATE(D23_AX(0:num_spec**2-1))
    !
    IF (ALLOCATED(D31_AX)) DEALLOCATE(D31_AX)
    ALLOCATE(D31_AX(0:num_spec**2-1))
    IF (ALLOCATED(D32_AX)) DEALLOCATE(D32_AX)
    ALLOCATE(D32_AX(0:num_spec**2-1))
    IF (ALLOCATED(D33_AX)) DEALLOCATE(D33_AX)
    ALLOCATE(D33_AX(0:num_spec**2-1))
    !
    ! Allocate storage arrays (non-axisymmetric solution, cgs-units)
    IF (ALLOCATED(D11_NA)) DEALLOCATE(D11_NA)
    ALLOCATE(D11_NA(0:num_spec**2-1))
    IF (ALLOCATED(D12_NA)) DEALLOCATE(D12_NA)
    ALLOCATE(D12_NA(0:num_spec**2-1))
    IF (ALLOCATED(D13_NA)) DEALLOCATE(D13_NA)
    ALLOCATE(D13_NA(0:num_spec**2-1))
    !
    IF (ALLOCATED(D21_NA)) DEALLOCATE(D21_NA)
    ALLOCATE(D21_NA(0:num_spec**2-1))
    IF (ALLOCATED(D22_NA)) DEALLOCATE(D22_NA)
    ALLOCATE(D22_NA(0:num_spec**2-1))
    IF (ALLOCATED(D23_NA)) DEALLOCATE(D23_NA)
    ALLOCATE(D23_NA(0:num_spec**2-1))
    !
    IF (ALLOCATED(D31_NA)) DEALLOCATE(D31_NA)
    ALLOCATE(D31_NA(0:num_spec**2-1))
    IF (ALLOCATED(D32_NA)) DEALLOCATE(D32_NA)
    ALLOCATE(D32_NA(0:num_spec**2-1))
    IF (ALLOCATED(D33_NA)) DEALLOCATE(D33_NA)
    ALLOCATE(D33_NA(0:num_spec**2-1))
    !
    ! Computation of normalized diffusion coefficients
    ! (w.r.t. 0-species (=electrons))
    IF(isw_qflux_NA .EQ. 1) THEN
       !
       IF (.NOT. ALLOCATED(qflux_ntv_allspec)) &
            STOP "Non-axisymmetric solution does not exist!"
       !
       IF (ALLOCATED(qflux_symm_allspec)) THEN
          !
          CALL compute_Dijab(qflux_ntv_allspec, qflux_symm_allspec, &
               row_ind_spec, col_ind_spec, row_ind_ptr, col_ind_ptr, &
               Dijab_NA, Dijab_AX, Dijab_norm_NA, Dijab_norm_AX)
          !
          ! non-axisymmetric coefficients (normalized w.r.t. 0-species)
          D11_NA_Dpl = Dijab_norm_NA(1,1,:)
          D12_NA_Dpl = Dijab_norm_NA(1,2,:)
          D13_NA_D31ref = Dijab_norm_NA(1,3,:)
          !
          D21_NA_Dpl = Dijab_norm_NA(2,1,:)
          D22_NA_Dpl = Dijab_norm_NA(2,2,:)
          D23_NA_D31ref = Dijab_norm_NA(2,3,:)
          !
          D31_NA_D31ref = Dijab_norm_NA(3,1,:)
          D32_NA_D31ref = Dijab_norm_NA(3,2,:)
          D33_NA_norm = Dijab_norm_NA(3,3,:)
          !
          ! axisymmetric coefficients (normalized w.r.t. 0-species)
          D11_AX_Dpl = Dijab_norm_AX(1,1,:)
          D12_AX_Dpl = Dijab_norm_AX(1,2,:)
          D13_AX_D31ref = Dijab_norm_AX(1,3,:)
          !
          D21_AX_Dpl = Dijab_norm_AX(2,1,:)
          D22_AX_Dpl = Dijab_norm_AX(2,2,:)
          D23_AX_D31ref = Dijab_norm_AX(2,3,:)
          !
          D31_AX_D31ref = Dijab_norm_AX(3,1,:)
          D32_AX_D31ref = Dijab_norm_AX(3,2,:)
          D33_AX_norm = Dijab_norm_AX(3,3,:)
          !
          ! non-axisymmetric coefficients (cgs-units)
          D11_NA = Dijab_NA(1,1,:)
          D12_NA = Dijab_NA(1,2,:)
          D13_NA = Dijab_NA(1,3,:)
          !
          D21_NA = Dijab_NA(2,1,:)
          D22_NA = Dijab_NA(2,2,:)
          D23_NA = Dijab_NA(2,3,:)
          !
          D31_NA = Dijab_NA(3,1,:)
          D32_NA = Dijab_NA(3,2,:)
          D33_NA = Dijab_NA(3,3,:)
          !
          ! axisymmetric coefficients (cgs-units)
          D11_AX = Dijab_AX(1,1,:)
          D12_AX = Dijab_AX(1,2,:)
          D13_AX = Dijab_AX(1,3,:)
          !
          D21_AX = Dijab_AX(2,1,:)
          D22_AX = Dijab_AX(2,2,:)
          D23_AX = Dijab_AX(2,3,:)
          !
          D31_AX = Dijab_AX(3,1,:)
          D32_AX = Dijab_AX(3,2,:)
          D33_AX = Dijab_AX(3,3,:)
          !
       ELSE
          !
          IF (ALLOCATED(dummy_qflux)) DEALLOCATE(dummy_qflux)
          ALLOCATE(dummy_qflux(1:3,1:3,0:num_spec,0:num_spec))
          dummy_qflux = 0.0_dp
          !
          CALL compute_Dijab(qflux_ntv_allspec, dummy_qflux, &
               row_ind_spec, col_ind_spec, row_ind_ptr, col_ind_ptr, &
               Dijab_NA, Dijab_AX, Dijab_norm_NA, Dijab_norm_AX)
          !
          DEALLOCATE(dummy_qflux)
          !
          ! non-axisymmetric coefficients (normalized w.r.t. 0-species)
          D11_NA_Dpl = Dijab_norm_NA(1,1,:)
          D12_NA_Dpl = Dijab_norm_NA(1,2,:)
          D13_NA_D31ref = Dijab_norm_NA(1,3,:)
          !
          D21_NA_Dpl = Dijab_norm_NA(2,1,:)
          D22_NA_Dpl = Dijab_norm_NA(2,2,:)
          D23_NA_D31ref = Dijab_norm_NA(2,3,:)
          !
          D31_NA_D31ref = Dijab_norm_NA(3,1,:)
          D32_NA_D31ref = Dijab_norm_NA(3,2,:)
          D33_NA_norm = Dijab_norm_NA(3,3,:)
          !
          ! in case of bounce-averaged model the solution
          ! for the axisymmetric problem is not computed
          D11_AX_Dpl = 0.0_dp
          D12_AX_Dpl = 0.0_dp
          D13_AX_D31ref = 0.0_dp
          !
          D21_AX_Dpl = 0.0_dp
          D22_AX_Dpl = 0.0_dp
          D23_AX_D31ref = 0.0_dp
          !
          D31_AX_D31ref = 0.0_dp
          D32_AX_D31ref = 0.0_dp
          D33_AX_norm = 0.0_dp
          !
          ! non-axisymmetric coefficients (cgs-units)
          D11_NA = Dijab_NA(1,1,:)
          D12_NA = Dijab_NA(1,2,:)
          D13_NA = Dijab_NA(1,3,:)
          !
          D21_NA = Dijab_NA(2,1,:)
          D22_NA = Dijab_NA(2,2,:)
          D23_NA = Dijab_NA(2,3,:)
          !
          D31_NA = Dijab_NA(3,1,:)
          D32_NA = Dijab_NA(3,2,:)
          D33_NA = Dijab_NA(3,3,:)
          !
          ! axisymmetric coefficients (cgs-units)
          D11_AX = 0.0_dp
          D12_AX = 0.0_dp
          D13_AX = 0.0_dp
          !
          D21_AX = 0.0_dp
          D22_AX = 0.0_dp
          D23_AX = 0.0_dp
          !
          D31_AX = 0.0_dp
          D32_AX = 0.0_dp
          D33_AX = 0.0_dp
          !
       END IF
       !
    ELSE
       !
       IF (.NOT. ALLOCATED(qflux_symm_allspec)) &
            STOP "Axisymmetric solution does not exist!"
       !
       IF (ALLOCATED(dummy_qflux)) DEALLOCATE(dummy_qflux)
       ALLOCATE(dummy_qflux(1:3,1:3,0:num_spec,0:num_spec))
       dummy_qflux = 0.0_dp
       !
       CALL compute_Dijab(dummy_qflux, qflux_symm_allspec, &
            row_ind_spec, col_ind_spec, row_ind_ptr, col_ind_ptr, &
            Dijab_NA, Dijab_AX, Dijab_norm_NA, Dijab_norm_AX)
       !PRINT *,qflux_symm_allspec(1,1,0,0),qflux_symm_allspec(1,1,0,1)
       !PRINT *,qflux_symm_allspec(1,1,1,0),qflux_symm_allspec(1,1,1,1)
       !
       DEALLOCATE(dummy_qflux)
       !
       ! only axisymmteric solution has been computed
       ! (non-axisymmteric coefficients set to zero)
       D11_NA_Dpl = 0.0_dp
       D12_NA_Dpl = 0.0_dp
       D13_NA_D31ref = 0.0_dp
       !
       D21_NA_Dpl = 0.0_dp
       D22_NA_Dpl = 0.0_dp
       D23_NA_D31ref = 0.0_dp
       !
       D31_NA_D31ref = 0.0_dp
       D32_NA_D31ref = 0.0_dp
       D33_NA_norm = 0.0_dp
       !
       ! axisymmetric coefficients (normalized w.r.t. 0-species)
       D11_AX_Dpl = Dijab_norm_AX(1,1,:)
       D12_AX_Dpl = Dijab_norm_AX(1,2,:)
       D13_AX_D31ref = Dijab_norm_AX(1,3,:)
       !
       D21_AX_Dpl = Dijab_norm_AX(2,1,:)
       D22_AX_Dpl = Dijab_norm_AX(2,2,:)
       D23_AX_D31ref = Dijab_norm_AX(2,3,:)
       !
       D31_AX_D31ref = Dijab_norm_AX(3,1,:)
       D32_AX_D31ref = Dijab_norm_AX(3,2,:)
       D33_AX_norm = Dijab_norm_AX(3,3,:)
       !
       ! non-axisymmetric coefficients (cgs-units)
       D11_NA = 0.0_dp
       D12_NA = 0.0_dp
       D13_NA = 0.0_dp
       !
       D21_NA = 0.0_dp
       D22_NA = 0.0_dp
       D23_NA = 0.0_dp
       !
       D31_NA = 0.0_dp
       D32_NA = 0.0_dp
       D33_NA = 0.0_dp
       !
       ! axisymmetric coefficients (cgs-units)
       D11_AX = Dijab_AX(1,1,:)
       D12_AX = Dijab_AX(1,2,:)
       D13_AX = Dijab_AX(1,3,:)
       !
       D21_AX = Dijab_AX(2,1,:)
       D22_AX = Dijab_AX(2,2,:)
       D23_AX = Dijab_AX(2,3,:)
       !
       D31_AX = Dijab_AX(3,1,:)
       D32_AX = Dijab_AX(3,2,:)
       D33_AX = Dijab_AX(3,3,:)
       !
    END IF
    !
    ! Computation of dimensional diffusion coefficients (cgs-units)
    !
    ! Compute D_p, D_{31,ref} and D_{33]^L for 0-species (=electrons)
    Dp00 = PI*vt0*(rho0**2)/(16.0_dp*aiota_loc*rt0)
    D31ref00 = vt0 * rho0 * bcovar_phi_hat * (bmod0*1.0e4_dp) / &
         (2.0_dp * aiota_loc * boozer_psi_pr_hat * avnabpsi)
    D33L00_Zeff = (-16.0_dp/(3.0_dp*PI)) * ((bmod0*1.0e4_dp)**2) * &
         vt0 * (2.0_dp/collpar_spec(0))
    !
    ! Write output to HDF5-file
    !
    ! compute nu^\star
    IF (ALLOCATED(nu_star_spec)) DEALLOCATE(nu_star_spec)
    ALLOCATE(nu_star_spec(0:num_spec-1))
    nu_star_spec = collpar_spec*rt0/aiota_loc
    !
    ! compute additional B-field quantities
    sqrtg_bctrvr_phi=sqrtg_bctrvr_tht/aiota_loc
    bcovar_phi=hcovar_tmp(2)*(bmod_tmp*1.0e4_dp)
    bcovar_tht=hcovar_tmp(3)*(bmod_tmp*1.0e4_dp)
    avbhat=y(14)/y(13)
    avbhat2=y(9)/y(6)
    avb2=avbhat2*((bmod0*1.0e4_dp)**2)
    !
    ! compute radial electric field and species Mach numbers
    IF (isw_calc_Er .EQ. 1) THEN
       IF (num_spec .EQ. 1) THEN
          !
          PRINT *,'------------------------------'
          CALL get_Er(qflux_symm_allspec, Er)
          PRINT *,'Er: ', Er
          PRINT *,'------------------------------'
          !
          CALL compute_VthtB_and_VphiB(row_ind_ptr, col_ind_ptr, &
               D31_AX, D32_AX, Er, VthtB_spec, VphiB_spec)
          !
          CALL compute_Vphi_profile(row_ind_ptr, col_ind_ptr, &
               D31_AX, D32_AX, Er, R_Vphi_prof, Z_Vphi_prof, &
               Vphi_prof_spec, Vtht_prof_spec)
          !
          CALL compute_Gamma(row_ind_ptr, col_ind_ptr, &
               D11_AX, D12_AX, Er, Gamma_AX_spec)
          CALL compute_Qflux(row_ind_ptr, col_ind_ptr, &
               D21_AX, D22_AX, Er, Qflux_AX_spec)
          CALL compute_ParFlow(row_ind_ptr, col_ind_ptr, &
               D31_AX, D32_AX, Er, ParFlow_AX_spec)
          CALL compute_Gamma(row_ind_ptr, col_ind_ptr, &
               D11_NA, D12_NA, Er, Gamma_NA_spec)
          CALL compute_TphiNA(Gamma_NA_spec, TphiNA_spec, TphiNA_tot)
          CALL compute_Qflux(row_ind_ptr, col_ind_ptr, &
               D21_NA, D22_NA, Er, Qflux_NA_spec)
          CALL compute_ParFlow(row_ind_ptr, col_ind_ptr, &
               D31_NA, D32_NA, Er, ParFlow_NA_spec)
          !
       ELSE
          !
          PRINT *,'------------------------------'
          CALL get_Er(qflux_symm_allspec, Er)
          CALL compute_Vphi_profile(row_ind_ptr, col_ind_ptr, &
               D31_AX, D32_AX, Er, R_Vphi_prof, Z_Vphi_prof, &
               Vphi_prof_woWare_spec, Vtht_prof_woWare_spec)
          PRINT *,'Er: ', Er
          PRINT *,'------------------------------'
          CALL compute_Er_and_A3norm(row_ind_ptr, col_ind_ptr, D31_AX, D32_AX, &
               D33_AX, Er, avEparB_ov_avb2)
          PRINT *,'Er, avEparB_ov_avb2: ', Er, avEparB_ov_avb2
          PRINT *,'------------------------------'
          CALL get_Er(qflux_symm_allspec, Er, avEparB_ov_avb2)
          PRINT *,'Er, avEparB_ov_avb2: ', Er, avEparB_ov_avb2
          PRINT *,'------------------------------'
          !
          CALL compute_VthtB_and_VphiB_b(row_ind_ptr, col_ind_ptr, &
               D31_AX, D32_AX, D33_AX, Er, avEparB_ov_avb2, &
               VthtB_spec, VphiB_spec, VthtB_Ware_spec, VphiB_Ware_spec)
          !
          CALL compute_Vphi_profile(row_ind_ptr, col_ind_ptr, &
               & D31_AX, D32_AX, Er, R_Vphi_prof, Z_Vphi_prof, &
               & Vphi_prof_spec, Vtht_prof_spec, D33_AX, avEparB_ov_avb2)
          !
          CALL compute_Gamma(row_ind_ptr, col_ind_ptr, &
               D11_AX, D12_AX, D13_AX, Er, avEparB_ov_avb2, &
               Gamma_AX_spec, Gamma_AX_Ware_spec)
          CALL compute_Qflux(row_ind_ptr, col_ind_ptr, &
               D21_AX, D22_AX, D23_AX, Er, avEparB_ov_avb2, &
               Qflux_AX_spec, Qflux_AX_Ware_spec)
          CALL compute_ParFlow(row_ind_ptr, col_ind_ptr, &
               D31_AX, D32_AX, D33_AX, Er, avEparB_ov_avb2, &
               ParFlow_AX_spec, ParFlow_AX_Ware_spec)
          CALL compute_Gamma(row_ind_ptr, col_ind_ptr, &
               D11_NA, D12_NA, D13_NA, Er, avEparB_ov_avb2, &
               Gamma_NA_spec, Gamma_NA_Ware_spec)
          CALL compute_TphiNA(Gamma_NA_spec, TphiNA_spec, TphiNA_tot)
          CALL compute_TphiNA(Gamma_NA_Ware_spec, TphiNA_Ware_spec, TphiNA_Ware_tot)
          CALL compute_Qflux(row_ind_ptr, col_ind_ptr, &
               D21_NA, D22_NA, D23_NA, Er, avEparB_ov_avb2, &
               Qflux_NA_spec, Qflux_NA_Ware_spec)
          CALL compute_ParFlow(row_ind_ptr, col_ind_ptr, &
               D31_NA, D32_NA, D33_NA, Er, avEparB_ov_avb2, &
               ParFlow_NA_spec, ParFlow_NA_Ware_spec)
          !
       END IF
    END IF
    !
    ! initialize HDF5 file
    file_name='neo2_multispecies_out.h5'
    CALL h5_create(TRIM(ADJUSTL(file_name)), h5id_multispec)
    !
    ! add B-field used for normalizations and further computations
    CALL h5_add(h5id_multispec, 'boozer_s', boozer_s, comment='', unit='1')
    CALL h5_add(h5id_multispec, 'aiota', aiota_loc, comment='toroidal transform, 1/q', unit='1')
    call h5_add(h5id_multispec, 'diota_ds', boozer_iota_s, &
      & comment='as calculated from spline interpolation', unit='1')
    CALL h5_add(h5id_multispec, 'R0', rt0, comment='major radius', unit='cm')
    CALL h5_add(h5id_multispec, 'Bref', (bmod0*1.0e4_dp), comment='reference magnetic field in gauss', unit='G')
    CALL h5_add(h5id_multispec, 'psi_pr_hat', boozer_psi_pr_hat)
    CALL h5_add(h5id_multispec, 'avnabpsi', avnabpsi)
    CALL h5_add(h5id_multispec, 'sqrtg_bctrvr_tht', sqrtg_bctrvr_tht)
    CALL h5_add(h5id_multispec, 'sqrtg_bctrvr_phi', sqrtg_bctrvr_phi)
    CALL h5_add(h5id_multispec, 'bcovar_tht', bcovar_tht)
    CALL h5_add(h5id_multispec, 'bcovar_phi', bcovar_phi)
    CALL h5_add(h5id_multispec, 'dbcovar_theta_ds', dbcovar_theta_hat_ds*bmod0*1.0e4_dp)
    CALL h5_add(h5id_multispec, 'dbcovar_phi_ds', dbcovar_phi_hat_ds*bmod0*1.0e4_dp)
    CALL h5_add(h5id_multispec, 'avbhat', avbhat)
    CALL h5_add(h5id_multispec, 'avbhat2', avbhat2)
    CALL h5_add(h5id_multispec, 'av_inv_bhat', av_inv_bhat_val)
    CALL h5_add(h5id_multispec, 'av_gphph', av_gphph_val)
    CALL h5_add(h5id_multispec, 'm_phi', DBLE(m_phi))
    CALL h5_add(h5id_multispec, 'eps_M_2', eps_M_2_val)
    !
    ! add species tags, charge numbers, mass, temperatures,
    ! collpar, ...
    CALL h5_add(h5id_multispec, 'num_spec', num_spec, comment='number of species', unit='1')
    CALL h5_add(h5id_multispec, 'isw_coul_log', isw_coul_log)
    CALL h5_add(h5id_multispec, 'species_tag', species_tag, &
         LBOUND(species_tag), UBOUND(species_tag))
    CALL h5_add(h5id_multispec, 'z_spec', z_spec, LBOUND(z_spec), UBOUND(z_spec), comment='charges of the species', unit='1')
    CALL h5_add(h5id_multispec, 'm_spec', m_spec, LBOUND(m_spec), UBOUND(m_spec), comment='mass of the species', unit='g')
    CALL h5_add(h5id_multispec, 'n_spec', n_spec, LBOUND(n_spec), UBOUND(n_spec), comment='density of the species', unit='1/cm^3')
    CALL h5_add(h5id_multispec, 'T_spec', T_spec, LBOUND(T_spec), UBOUND(T_spec), comment='temperature of the species', unit='erg')
    CALL h5_add(h5id_multispec, 'collpar_spec', collpar_spec, &
         LBOUND(collpar_spec), UBOUND(collpar_spec))
    CALL h5_add(h5id_multispec, 'nu_star_spec', nu_star_spec, &
         LBOUND(nu_star_spec), UBOUND(nu_star_spec))
    !
    ! add D_p, D_{31,ref} and D_{33]^L for 0-species (=electrons)
    CALL h5_add(h5id_multispec, 'Dp0', Dp00)
    CALL h5_add(h5id_multispec, 'D31ref0', D31ref00)
    CALL h5_add(h5id_multispec, 'D33L0_Zeff', D33L00_Zeff)
    !
    ! add row- and column-indices of diffusion coefficients
    CALL h5_add(h5id_multispec, 'row_ind_spec', row_ind_spec, LBOUND(row_ind_spec), UBOUND(row_ind_spec))
    CALL h5_add(h5id_multispec, 'col_ind_spec', col_ind_spec, LBOUND(col_ind_spec), UBOUND(col_ind_spec))
    !
    ! add normalized diffusion coefficients
    ! (w.r.t. 0-species (=electrons))
    !-> axisymmetric solution
    CALL h5_add(h5id_multispec, 'D11_AX_Dpl', D11_AX_Dpl, &
         LBOUND(D11_AX_Dpl), UBOUND(D11_AX_Dpl))
    CALL h5_add(h5id_multispec, 'D12_AX_Dpl', D12_AX_Dpl, &
         LBOUND(D12_AX_Dpl), UBOUND(D12_AX_Dpl))
    CALL h5_add(h5id_multispec, 'D13_AX_D31ref', D13_AX_D31ref, &
         LBOUND(D13_AX_D31ref), UBOUND(D13_AX_D31ref))
    !
    CALL h5_add(h5id_multispec, 'D21_AX_Dpl', D21_AX_Dpl, &
         LBOUND(D21_AX_Dpl), UBOUND(D21_AX_Dpl))
    CALL h5_add(h5id_multispec, 'D22_AX_Dpl', D22_AX_Dpl, &
         LBOUND(D22_AX_Dpl), UBOUND(D22_AX_Dpl))
    CALL h5_add(h5id_multispec, 'D23_AX_D31ref', D23_AX_D31ref, &
         LBOUND(D23_AX_D31ref), UBOUND(D23_AX_D31ref))
    !
    CALL h5_add(h5id_multispec, 'D31_AX_D31ref', D31_AX_D31ref, &
         LBOUND(D31_AX_D31ref), UBOUND(D31_AX_D31ref))
    CALL h5_add(h5id_multispec, 'D32_AX_D31ref', D32_AX_D31ref, &
         LBOUND(D32_AX_D31ref), UBOUND(D32_AX_D31ref))
    CALL h5_add(h5id_multispec, 'D33_AX_norm', D33_AX_norm, &
         LBOUND(D33_AX_norm), UBOUND(D33_AX_norm))
    !-> non-axisymmetric solution
    CALL h5_add(h5id_multispec, 'D11_NA_Dpl', D11_NA_Dpl, &
         LBOUND(D11_NA_Dpl), UBOUND(D11_NA_Dpl))
    CALL h5_add(h5id_multispec, 'D12_NA_Dpl', D12_NA_Dpl, &
         LBOUND(D12_NA_Dpl), UBOUND(D12_NA_Dpl))
    CALL h5_add(h5id_multispec, 'D13_NA_D31ref', D13_NA_D31ref, &
         LBOUND(D13_NA_D31ref), UBOUND(D13_NA_D31ref))
    !
    CALL h5_add(h5id_multispec, 'D21_NA_Dpl', D21_NA_Dpl, &
         LBOUND(D21_NA_Dpl), UBOUND(D21_NA_Dpl))
    CALL h5_add(h5id_multispec, 'D22_NA_Dpl', D22_NA_Dpl, &
         LBOUND(D22_NA_Dpl), UBOUND(D22_NA_Dpl))
    CALL h5_add(h5id_multispec, 'D23_NA_D31ref', D23_NA_D31ref, &
         LBOUND(D23_NA_D31ref), UBOUND(D23_NA_D31ref))
    !
    CALL h5_add(h5id_multispec, 'D31_NA_D31ref', D31_NA_D31ref, &
         LBOUND(D31_NA_D31ref), UBOUND(D31_NA_D31ref))
    CALL h5_add(h5id_multispec, 'D32_NA_D31ref', D32_NA_D31ref, &
         LBOUND(D32_NA_D31ref), UBOUND(D32_NA_D31ref))
    CALL h5_add(h5id_multispec, 'D33_NA_norm', D33_NA_norm, &
         LBOUND(D33_NA_norm), UBOUND(D33_NA_norm))

    ! add dimensional diffusion coefficients (cgs-units)
    !-> axisymmetric solution
    CALL h5_add(h5id_multispec, 'D11_AX', D11_AX, LBOUND(D11_AX), UBOUND(D11_AX), &
      & comment='dimensional diffusion coefficient for axisymmetric solution', unit='cm^2/s')
    CALL h5_add(h5id_multispec, 'D12_AX', D12_AX, LBOUND(D12_AX), UBOUND(D12_AX), &
      & comment='dimensional diffusion coefficient for axisymmetric solution', unit='cm^2/s')
    CALL h5_add(h5id_multispec, 'D13_AX', D13_AX, LBOUND(D13_AX), UBOUND(D13_AX), &
      & comment='dimensional diffusion coefficient for axisymmetric solution', unit='cm^2/s')

    CALL h5_add(h5id_multispec, 'D21_AX', D21_AX, LBOUND(D21_AX), UBOUND(D21_AX), &
      & comment='dimensional diffusion coefficient for axisymmetric solution', unit='cm^2/s')
    CALL h5_add(h5id_multispec, 'D22_AX', D22_AX, LBOUND(D22_AX), UBOUND(D22_AX), &
      & comment='dimensional diffusion coefficient for axisymmetric solution', unit='cm^2/s')
    CALL h5_add(h5id_multispec, 'D23_AX', D23_AX, LBOUND(D23_AX), UBOUND(D23_AX), &
      & comment='dimensional diffusion coefficient for axisymmetric solution', unit='cm^2/s')

    CALL h5_add(h5id_multispec, 'D31_AX', D31_AX, LBOUND(D31_AX), UBOUND(D31_AX), &
      & comment='dimensional diffusion coefficient for axisymmetric solution', unit='cm^2/s')
    CALL h5_add(h5id_multispec, 'D32_AX', D32_AX, LBOUND(D32_AX), UBOUND(D32_AX), &
      & comment='dimensional diffusion coefficient for axisymmetric solution', unit='cm^2/s')
    CALL h5_add(h5id_multispec, 'D33_AX', D33_AX, LBOUND(D33_AX), UBOUND(D33_AX), &
      & comment='dimensional diffusion coefficient for axisymmetric solution', unit='cm^2/s')
    !-> non-axisymmetric solution
    CALL h5_add(h5id_multispec, 'D11_NA', D11_NA, LBOUND(D11_NA), UBOUND(D11_NA), &
      & comment='dimensional diffusion coefficient for non-axisymmetric solution', unit='cm^2/s')
    CALL h5_add(h5id_multispec, 'D12_NA', D12_NA, LBOUND(D12_NA), UBOUND(D12_NA), &
      & comment='dimensional diffusion coefficient for non-axisymmetric solution', unit='cm^2/s')
    CALL h5_add(h5id_multispec, 'D13_NA', D13_NA, LBOUND(D13_NA), UBOUND(D13_NA), &
      & comment='dimensional diffusion coefficient for non-axisymmetric solution', unit='cm^2/s')

    CALL h5_add(h5id_multispec, 'D21_NA', D21_NA, LBOUND(D21_NA), UBOUND(D21_NA), &
      & comment='dimensional diffusion coefficient for non-axisymmetric solution', unit='cm^2/s')
    CALL h5_add(h5id_multispec, 'D22_NA', D22_NA, LBOUND(D22_NA), UBOUND(D22_NA), &
      & comment='dimensional diffusion coefficient for non-axisymmetric solution', unit='cm^2/s')
    CALL h5_add(h5id_multispec, 'D23_NA', D23_NA, LBOUND(D23_NA), UBOUND(D23_NA), &
      & comment='dimensional diffusion coefficient for non-axisymmetric solution', unit='cm^2/s')

    CALL h5_add(h5id_multispec, 'D31_NA', D31_NA, LBOUND(D31_NA), UBOUND(D31_NA), &
      & comment='dimensional diffusion coefficient for non-axisymmetric solution', unit='cm^2/s')
    CALL h5_add(h5id_multispec, 'D32_NA', D32_NA, LBOUND(D32_NA), UBOUND(D32_NA), &
      & comment='dimensional diffusion coefficient for non-axisymmetric solution', unit='cm^2/s')
    CALL h5_add(h5id_multispec, 'D33_NA', D33_NA, LBOUND(D33_NA), UBOUND(D33_NA), &
      & comment='dimensional diffusion coefficient for non-axisymmetric solution', unit='cm^2/s')

    ! add radial electric field and species Mach numbers
    IF (isw_calc_Er .EQ. 1) THEN
       !
       CALL h5_add(h5id_multispec, 'Er', Er)
       CALL h5_add(h5id_multispec, 'MtOvR', MtOvR_spec, &
            LBOUND(MtOvR_spec), UBOUND(MtOvR_spec))
       !
       CALL h5_add(h5id_multispec, 'VthtB_spec', VthtB_spec, &
            LBOUND(VthtB_spec), UBOUND(VthtB_spec))
       CALL h5_add(h5id_multispec, 'VphiB_spec', VphiB_spec, &
            LBOUND(VphiB_spec), UBOUND(VphiB_spec))
       !
       CALL h5_add(h5id_multispec, 'R_Vphi_prof', R_Vphi_prof, &
            LBOUND(R_Vphi_prof), UBOUND(R_Vphi_prof))
       CALL h5_add(h5id_multispec, 'Z_Vphi_prof', Z_Vphi_prof, &
            LBOUND(Z_Vphi_prof), UBOUND(Z_Vphi_prof))
       CALL h5_add(h5id_multispec, 'Vphi_prof_spec', Vphi_prof_spec, &
            LBOUND(Vphi_prof_spec), UBOUND(Vphi_prof_spec))
       CALL h5_add(h5id_multispec, 'Vtht_prof_spec', Vtht_prof_spec, &
            LBOUND(Vtht_prof_spec), UBOUND(Vtht_prof_spec))
       !
       CALL h5_add(h5id_multispec, 'Gamma_AX_spec', Gamma_AX_spec, &
            LBOUND(Gamma_AX_spec), UBOUND(Gamma_AX_spec))
       CALL h5_add(h5id_multispec, 'Qflux_AX_spec', Qflux_AX_spec, &
            LBOUND(Qflux_AX_spec), UBOUND(Qflux_AX_spec))
       CALL h5_add(h5id_multispec, 'ParFlow_AX_spec', ParFlow_AX_spec, &
            LBOUND(ParFlow_AX_spec), UBOUND(ParFlow_AX_spec))
       !
       CALL h5_add(h5id_multispec, 'Gamma_NA_spec', Gamma_NA_spec, &
            LBOUND(Gamma_NA_spec), UBOUND(Gamma_NA_spec))
       CALL h5_add(h5id_multispec, 'Qflux_NA_spec', Qflux_NA_spec, &
            LBOUND(Qflux_NA_spec), UBOUND(Qflux_NA_spec))
       CALL h5_add(h5id_multispec, 'ParFlow_NA_spec', ParFlow_NA_spec, &
            LBOUND(ParFlow_NA_spec), UBOUND(ParFlow_NA_spec))
       !
       CALL h5_add(h5id_multispec, 'TphiNA_spec', TphiNA_spec, &
            LBOUND(TphiNA_spec), UBOUND(TphiNA_spec))
       CALL h5_add(h5id_multispec, 'TphiNA_tot', TphiNA_tot)
       !
       IF (num_spec .GT. 1) THEN
          CALL h5_add(h5id_multispec, 'avEparB_ov_avb2', avEparB_ov_avb2)
          !
          CALL h5_add(h5id_multispec, 'VthtB_Ware_spec', VthtB_Ware_spec, &
               LBOUND(VthtB_Ware_spec), UBOUND(VthtB_Ware_spec))
          CALL h5_add(h5id_multispec, 'VphiB_Ware_spec', VphiB_Ware_spec, &
               LBOUND(VphiB_Ware_spec), UBOUND(VphiB_Ware_spec))
          !
          CALL h5_add(h5id_multispec, 'Vphi_prof_woWare_spec', Vphi_prof_woWare_spec, &
               LBOUND(Vphi_prof_woWare_spec), UBOUND(Vphi_prof_woWare_spec))
          CALL h5_add(h5id_multispec, 'Vtht_prof_woWare_spec', Vtht_prof_woWare_spec, &
               LBOUND(Vtht_prof_woWare_spec), UBOUND(Vtht_prof_woWare_spec))
          !
          CALL h5_add(h5id_multispec, 'Gamma_AX_Ware_spec', Gamma_AX_Ware_spec, &
               LBOUND(Gamma_AX_Ware_spec), UBOUND(Gamma_AX_Ware_spec))
          CALL h5_add(h5id_multispec, 'Qflux_AX_Ware_spec', Qflux_AX_Ware_spec, &
               LBOUND(Qflux_AX_Ware_spec), UBOUND(Qflux_AX_Ware_spec))
          CALL h5_add(h5id_multispec, 'ParFlow_AX_Ware_spec', ParFlow_AX_Ware_spec, &
               LBOUND(ParFlow_AX_Ware_spec), UBOUND(ParFlow_AX_Ware_spec))
          !
          CALL h5_add(h5id_multispec, 'Gamma_NA_Ware_spec', Gamma_NA_Ware_spec, &
               LBOUND(Gamma_NA_Ware_spec), UBOUND(Gamma_NA_Ware_spec))
          CALL h5_add(h5id_multispec, 'Qflux_NA_Ware_spec', Qflux_NA_Ware_spec, &
               LBOUND(Qflux_NA_Ware_spec), UBOUND(Qflux_NA_Ware_spec))
          CALL h5_add(h5id_multispec, 'ParFlow_NA_Ware_spec', ParFlow_NA_Ware_spec, &
               LBOUND(ParFlow_NA_Ware_spec), UBOUND(ParFlow_NA_Ware_spec))
          !
          CALL h5_add(h5id_multispec, 'TphiNA_Ware_spec', TphiNA_Ware_spec, &
               LBOUND(TphiNA_Ware_spec), UBOUND(TphiNA_Ware_spec))
          CALL h5_add(h5id_multispec, 'TphiNA_Ware_tot', TphiNA_Ware_tot)
       END IF
    END IF
    !
    CALL h5_close(h5id_multispec)

    if (.not. check_coefficients(.true.)) then
      write(*,*) 'WARNING: sanity checks of the D1-_AX coefficients failed.'
    end if
    if (.not. check_ambipolarity_particle_flux(.true.)) then
      write(*,*) 'WARNING: sanity check of ambipolarity of particle flux failed.'
    end if

  contains

    !> \brief Perform some sanity checks on the coefficients.
    function check_coefficients(verbose) result(passed)
      logical :: passed
      logical, intent(in) :: verbose

      passed = check_ambipolarity_conditions(verbose) &
        & .and. check_independence_radial_electric_field_condition(verbose)
    end function check_coefficients

    !> \brief Checks that come from the ambipolarity condition.
    function check_ambipolarity_conditions(verbose) result(passed)
      logical :: passed
      logical, intent(in) :: verbose

      passed = check_ambipolarity_condition_density(verbose) &
        & .and. check_ambipolarity_condition_temperature(verbose) &
        & .and. check_ambipolarity_condition_from_parallel_field(verbose) &
        & .and. check_ambipolarity_condition_from_radial_field(verbose)
    end function check_ambipolarity_conditions

    !> \brief Check from the ambpolarity condition involving the density.
    function check_ambipolarity_condition_density(verbose) result(passed)
      use collisionality_mod, only : num_spec, z_spec, n_spec

      implicit none

      logical :: passed
      logical, intent(in) :: verbose

      real(kind=dp) :: sum_d11, sum_abs_d11, summand
      integer :: k,l

      sum_d11 = 0.0
      sum_abs_d11 = 0.0

      do k = 1,num_spec
        do l = 1,num_spec
          ! \note dn_spec_ov_ds uses zero based index.
          summand = z_spec(k-1) * n_spec(k-1) &
            & * D11_AX(return_linear_species_index(k, l)) * dn_spec_ov_ds(l-1) / n_spec(l-1)
          sum_d11 = sum_d11 + summand
          sum_abs_d11 = sum_abs_d11 + abs(summand)
        end do
      end do

      if (abs(sum_d11 / sum_abs_d11) < epsilon_transport_coefficients) then
        passed = .true.
      else if(verbose) then
        write(*,*) 'WARNING: sanity check check_ambipolarity_condition_density failed.'
        write(*,*) '  Relative error: ', abs(sum_d11 / sum_abs_d11)
      end if
    end function check_ambipolarity_condition_density

    !> \brief Check from the ambpolarity condition involving the temperature.
    function check_ambipolarity_condition_temperature(verbose) result(passed)
      use collisionality_mod, only : num_spec, z_spec, T_spec, n_spec

      implicit none

      logical :: passed
      logical, intent(in) :: verbose

      real(kind=dp) :: sum_d112, sum_abs_d112
      integer :: k,l

      sum_d112 = 0.0
      sum_abs_d112 = 0.0

      do k = 1,num_spec
        do l = 1,num_spec
          ! \note dT_spec_ov_ds uses zero based index.
          sum_d112 = sum_d112 &
            & + z_spec(k-1) * n_spec(k-1) * dT_spec_ov_ds(l-1) / T_spec(l-1) &
            & * (3*D11_AX(return_linear_species_index(k, l))/2 - D12_AX(return_linear_species_index(k, l)) )
          sum_abs_d112 = sum_abs_d112 &
            & + abs(z_spec(k-1) * n_spec(k-1) * dT_spec_ov_ds(l-1) / T_spec(l-1) &
            & * (3*D11_AX(return_linear_species_index(k, l))/2 - D12_AX(return_linear_species_index(k, l)) ))
        end do
      end do

      if (abs(sum_d112 / sum_abs_d112) < epsilon_transport_coefficients) then
        passed = .true.
      else if(verbose) then
        write(*,*) 'WARNING: sanity check check_ambipolarity_condition_temperature failed.'
        write(*,*) '  Relative error: ', abs(sum_d112 / sum_abs_d112)
      end if
    end function check_ambipolarity_condition_temperature

    !> \brief Check from the ambpolarity condition from the term involving the parallel electric field.
    function check_ambipolarity_condition_from_parallel_field(verbose) result(passed)
      use collisionality_mod, only : num_spec, z_spec, T_spec, n_spec

      implicit none

      logical :: passed
      logical, intent(in) :: verbose

      real(kind=dp) :: sum_d13, sum_abs_d13
      integer :: k,l

      sum_d13 = 0.0
      sum_abs_d13 = 0.0

      do k = 1,num_spec
        do l = 1,num_spec
          sum_d13 = sum_d13 + z_spec(k-1) * n_spec(k-1) &
            & * D13_AX(return_linear_species_index(k, l)) * z_spec(l-1) / T_spec(l-1)
          sum_abs_d13 = sum_abs_d13 &
            & + abs(z_spec(k-1) * n_spec(k-1) &
            &  * D13_AX(return_linear_species_index(k, l)) * z_spec(l-1) / T_spec(l-1))
        end do
      end do

      if (abs(sum_d13 / sum_abs_d13) < epsilon_transport_coefficients) then
        passed = .true.
      else if(verbose) then
        write(*,*) 'WARNING: sanity check check_ambipolarity_condition_from_parallel_field failed.'
        write(*,*) '  Relative error: ', abs(sum_d13 / sum_abs_d13)
      end if
    end function check_ambipolarity_condition_from_parallel_field

    !> \brief Check from the ambpolarity condition from the term involving the radial electric field.
    function check_ambipolarity_condition_from_radial_field(verbose) result(passed)
      use collisionality_mod, only : num_spec, z_spec, T_spec, n_spec

      implicit none

      logical :: passed
      logical, intent(in) :: verbose

      real(kind=dp) :: sum_d11, sum_abs_d11
      integer :: k,l

      sum_d11 = 0.0
      sum_abs_d11 = 0.0

      do k = 1,num_spec
        do l = 1,num_spec
          sum_d11 = sum_d11 + z_spec(k-1) * n_spec(k-1) &
            & * D11_AX(return_linear_species_index(k, l)) * z_spec(l-1) / T_spec(l-1)
          sum_abs_d11 = sum_abs_d11 &
            & + abs(z_spec(k-1) * n_spec(k-1) &
            &  * D11_AX(return_linear_species_index(k, l)) * z_spec(l-1) / T_spec(l-1))
        end do
      end do

      if (abs(sum_d11 / sum_abs_d11) < epsilon_transport_coefficients) then
        passed = .true.
      else if(verbose) then
        write(*,*) 'WARNING: sanity check check_ambipolarity_condition_from_radial_field failed.'
        write(*,*) '  Relative error: ', abs(sum_d11 / sum_abs_d11)
      end if
    end function check_ambipolarity_condition_from_radial_field

    !> \brief Checks that come from the condition that particle flux is independent off the radial electric field.
    !>
    !> This function performs checks that come from the assumption that
    !> the particle flux does not depend on the radial electric field,
    !> which is valid when the centrifugal forces can be neglected.
    function check_independence_radial_electric_field_condition(verbose) result(passed)
      use collisionality_mod, only : num_spec, z_spec, T_spec

      implicit none

      logical :: passed
      logical, intent(in) :: verbose

      real(kind=dp), dimension(:), allocatable :: d11_alpha, d11_abs_alpha
      integer :: k,l

      passed = .false.

      allocate(d11_alpha(1:num_spec))
      d11_alpha = 0.0
      allocate(d11_abs_alpha(1:num_spec))
      d11_abs_alpha = 0.0

      do k = 1,num_spec
        do l = 1,num_spec
          d11_alpha(k) = d11_alpha(k) &
            & + D11_AX(return_linear_species_index(k, l)) * z_spec(l-1) / T_spec(l-1)
          d11_abs_alpha(k) = d11_abs_alpha(k) &
            & + abs(D11_AX(return_linear_species_index(k, l)) * z_spec(l-1) / T_spec(l-1))
        end do
      end do

      if (all(abs(d11_alpha / d11_abs_alpha) < epsilon_transport_coefficients)) then
        passed = .true.
      else if(verbose) then
        write(*,*) 'WARNING: sanity check check_independence_radial_electric_field_condition failed.'
        write(*,*) '  Relative error: ', abs(d11_alpha / d11_abs_alpha)
      end if

      if (allocated(d11_alpha)) deallocate(d11_alpha)

    end function check_independence_radial_electric_field_condition

    function return_linear_species_index(k, l) result(ind)
      use collisionality_mod, only : num_spec

      implicit none

      integer, intent(in) :: k, l

      integer :: ind

      ind = (k-1)*num_spec + (l-1)
    end function return_linear_species_index

    function check_ambipolarity_particle_flux(verbose) result(passed)
      use collisionality_mod, only : num_spec, z_spec

      implicit none

      logical :: passed
      logical, intent(in) :: verbose

      real(kind=dp) :: sum_fluxes, sum_abs_fluxes

      integer :: k

      sum_fluxes = 0.0
      sum_abs_fluxes = 0.0
      passed = .false.

      if (isw_calc_Er == 1) then
        do k = 1, num_spec
          sum_fluxes = sum_fluxes + Gamma_AX_spec(k-1)*z_spec(k-1)
          sum_abs_fluxes = sum_abs_fluxes + abs(Gamma_AX_spec(k-1)*z_spec(k-1))
        end do

        if (abs(sum_fluxes/sum_abs_fluxes) < epsilon_particle_flux) then
          passed = .true.
        else if (verbose) then
          write(*,*) 'WARNING: particle flux not ambipolar to relative accuracy ', epsilon_particle_flux
          write(*,*) '  sum is: ', sum_fluxes, ' relative sum is: ', sum_fluxes/sum_abs_fluxes
        end if
      else
        write(*,*) 'WARNING: particle flux ambipolarity could not be checked,'
        write(*,*) '  as isw_calc_er is 0'
      end if
    end function check_ambipolarity_particle_flux
  END SUBROUTINE write_multispec_output_a
  !
  SUBROUTINE compute_Er(row_ind_ptr, col_ind_ptr, D31AX_spec, D32AX_spec, &
       & Er, D33AX_spec_in, avEparB_ov_avb2_in)
    !
    USE neo_precision
    USE neo_control, ONLY: lab_swi
    USE device_mod, ONLY : surface
    USE mag_interface_mod, ONLY : mag_coordinates, &
         boozer_s, boozer_theta_beg, boozer_phi_beg
    USE partpa_mod,  ONLY : bmod0
    USE mag_sub, ONLY: mag
    use neo_magfie, only : boozer_curr_pol_hat, boozer_psi_pr_hat, &
         compute_Gsymm, calc_thetaB_RZloc
    USE collisionality_mod, ONLY : num_spec, species_tag, &
         z_spec, m_spec, n_spec, T_spec
    !
    ! ---------------------------------------------------------------!
    ! input:
    ! ---------------------------------------------------------------!
    ! row- and column indices (=species) of diffusion coefficients
    INTEGER, DIMENSION(:), INTENT(in) :: row_ind_ptr, col_ind_ptr
    ! species parallel flow
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: D31AX_spec, D32AX_spec
    REAL(kind=dp), DIMENSION(:), INTENT(in), optional :: D33AX_spec_in
    ! drive A_3 normalized (=inductive electric field)
    REAL(kind=dp), INTENT(in), optional :: avEparB_ov_avb2_in
    ! ---------------------------------------------------------------!
    ! output:
    ! ---------------------------------------------------------------!
    ! radial electric field (w.r.t. effective radius)
    REAL(kind=dp), INTENT(out) :: Er
    ! ---------------------------------------------------------------!
    ! local:
    ! ---------------------------------------------------------------!
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: y
    REAL(kind=dp) :: aiota_loc, avnabpsi
    ! ---------------------------------------------------------------!
    ! co- and contra-variant B-field components using
    ! r_eff as a flux-surface label
    REAL(kind=dp) :: bcovar_tht, bcovar_phi
    REAL(kind=dp) :: sqrtg_bctrvr_phi, boozer_psi_pr
    REAL(kind=dp) :: avbhat2, avb2
    ! normalized co-variant phi-component of B,  $\sqrt{g}B^\vartheta$
    REAL(kind=dp) :: bcovar_phi_hat, sqrtg_bctrvr_tht
    ! Only used to call mag for normalizations
    ! related to D31 and D32
    REAL(kind=dp)                 :: bmod_tmp,sqrtg_tmp
    REAL(kind=dp), DIMENSION(3)   :: x_tmp,bder_tmp,hcovar_tmp,hctrvr_tmp
    REAL(kind=dp), DIMENSION(3,3) :: hcoder_tmp,hctder_tmp
    ! Boozer coordinates of local Vphi value (only used for isw_Vphi_loc=1)
    REAL(kind=dp)                 :: thetaB
    REAL(kind=dp), DIMENSION(3)   :: x_start
    ! transformation function Boozer coord. -> Symm. flux coord.
    REAL(kind=dp)                 :: G_symm, G_symm_tb, G_symm_pb
    ! ---------------------------------------------------------------!
    ! temperature, pressure and density of species i 
    ! (i = species of measured toroidal rotation frequency)
    REAL(kind=dp) :: T_ions, n_ions, p_ions, z_ions
    REAL(kind=dp) :: dT_ions_ov_dr, dn_ions_ov_dr, dp_ions_ov_dr
    ! loop indices, temporary variables
    INTEGER :: ispec_ctr, spec_i, irow_spec, icol_spec, num_ctr
    REAL(kind=dp) :: denom_Er, nom_Er, fac1
    REAL(kind=dp) :: denom_Er_1, nom_Er_1, nom_Er_2, nom_Er_3, nom_Er_4, nom_Er_5

    ! Variables for the optional parameters.
    real(kind=dp), dimension(:) , allocatable:: D33AX_spec
    real(kind=dp) :: avEparB_ov_avb2

    allocate(D33AX_spec(lbound(row_ind_ptr,1):ubound(row_ind_ptr,1)))
    if (present(D33AX_spec_in)) then
      D33AX_spec = D33AX_spec_in
    else
      D33AX_spec = 0.0
    end if
    if (present(avEparB_ov_avb2_in)) then
      avEparB_ov_avb2 = avEparB_ov_avb2_in
    else
      avEparB_ov_avb2 = 0.0
    end if
    ! ---------------------------------------------------------------!
    !
    ! copy y-vector (see definition in rhs_kin.f90)
    ALLOCATE(y(SIZE(y_ntv_mod,1)))
    y = y_ntv_mod
    !
    ! aiota, avnabpsi :
    aiota_loc = surface%aiota
    avnabpsi = y(7) / y(6)
    avbhat2 = y(9) / y(6)

    spec_i = -1

    denom_Er = 0.0
    nom_Er = 0.0

    ! computation of the normalization for D31 and D32 (-> D31_ref)
    IF (mag_coordinates .EQ. 0) THEN
       ! cylindrical coordinates
       PRINT *,"ntv_mod.f90: Computation of radial electric field &
            &only implemented for Boozer coordinates at the moment!"
       STOP
       !
       x_tmp = xstart_cyl
       CALL mag(x_tmp,bmod_tmp,sqrtg_tmp,bder_tmp,hcovar_tmp,&
            hctrvr_tmp,hcoder_tmp,hctder_tmp)
       !
       ! normalized co-variant phi-component of B
       ! (Note! co-variant phi-component of B is the
       ! same for cylindrical coordinates and symmetry flux
       ! coodrinates --> no conversion needed)
       bcovar_phi_hat = hcovar_tmp(2)*(bmod_tmp/bmod0)
       !
       ! restore value of $\sqrt{g}B^\vartheta$ for
       ! symmetry flux coordinates from quantities
       ! given in cylindircal coordinates
       sqrtg_bctrvr_tht = avnabpsi*sqrtg_tmp*(hctrvr_tmp(3)*bmod_tmp*1.0e4_dp)
    ELSE
       ! boozer coordinates
       !
       IF (isw_Vphi_loc .EQ. 0) THEN
          ! isw_Vphi_loc=0: value of "Vphi" corresponds to
          ! flux surface average (<V_\varphi>)
          x_tmp = (/boozer_s,boozer_phi_beg,boozer_theta_beg/)
          CALL mag(x_tmp,bmod_tmp,sqrtg_tmp,bder_tmp,hcovar_tmp,&
               hctrvr_tmp,hcoder_tmp,hctder_tmp)
       ELSE IF (isw_Vphi_loc .EQ. 1) THEN
          ! Caution: This branch is not tested!
          ! isw_Vphi_loc=1: value of "Vphi" is specified
          ! locally for given (R,Z)-position
          x_start = (/boozer_s,boozer_phi_beg,0.0_dp/)
          CALL calc_thetaB_RZloc(R_Vphi, Z_Vphi, x_start, thetaB)
          x_tmp = (/boozer_s,boozer_phi_beg,thetaB/)
          CALL mag(x_tmp,bmod_tmp,sqrtg_tmp,bder_tmp,hcovar_tmp,&
               hctrvr_tmp,hcoder_tmp,hctder_tmp)
          CALL compute_Gsymm( x_tmp, G_symm, G_symm_tb, G_symm_pb )
       ELSE IF (isw_Vphi_loc .EQ. 2) THEN
          ! isw_Vphi_loc=2: value of "Vphi" is specified
          ! locally for given \vartheta_B position
          x_tmp = (/boozer_s,boozer_phi_beg,boozer_theta_Vphi/)
          CALL mag(x_tmp,bmod_tmp,sqrtg_tmp,bder_tmp,hcovar_tmp,&
               hctrvr_tmp,hcoder_tmp,hctder_tmp)
          CALL compute_Gsymm( x_tmp, G_symm, G_symm_tb, G_symm_pb )
       ELSE
          PRINT *,"ntv_mod.f90: Undefined state of switch isw_Vphi_loc (= 0 / 1 / 2)!"
          STOP
       END IF
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
       IF (lab_swi .EQ. 10) THEN ! ASDEX-U (E. Strumberger)
          ! this is the same up to a minus sign resulting from the
          ! definition of sqrtg_tmp (left-handed system)
          sqrtg_bctrvr_tht = avnabpsi*sqrtg_tmp*(hctrvr_tmp(3)*bmod_tmp*1.0d4) 
       ELSE
          ! restore value of $\sqrt{g}B^\vartheta$ for
          ! symmetry flux coordinates from the quantities
          ! given in Boozer coordinates (only valid for right-handed system)
          sqrtg_bctrvr_tht = avnabpsi*aiota_loc*boozer_psi_pr_hat*(bmod0*1.0e4_dp)
       END IF
    END IF
    !
    ! compute additional B-field quantities
    sqrtg_bctrvr_phi=sqrtg_bctrvr_tht/aiota_loc
    bcovar_phi=hcovar_tmp(2)*(bmod_tmp*1.0e4_dp)
    bcovar_tht=hcovar_tmp(3)*(bmod_tmp*1.0e4_dp)
    boozer_psi_pr = boozer_psi_pr_hat*(bmod0*1.0e4_dp)
    avb2 = avbhat2*((bmod0*1.0e4_dp)**2)
    !
    ! detect species of measured V_phi
    num_ctr = 0
    DO ispec_ctr = 0,num_spec-1
       IF (species_tag(ispec_ctr) .EQ. species_tag_Vphi) THEN
          spec_i = ispec_ctr
          num_ctr = num_ctr + 1
       END IF
    END DO
    IF (num_ctr .EQ. 0) THEN
       PRINT *,"ntv_mod.f90: Subroutine compute_Er - &
            &Species of measured V_phi not found!"
       STOP
    ELSEIF (num_ctr .GT. 1) THEN
       PRINT *,"ntv_mod.f90: Subroutine compute_Er - &
            &Multiple definition of species V_phi!"
       STOP
    END IF
    !
    z_ions = z_spec(spec_i)
    T_ions = T_spec(spec_i)
    dT_ions_ov_dr = dT_spec_ov_ds(spec_i) * avnabpsi
    n_ions = n_spec(spec_i)
    dn_ions_ov_dr = dn_spec_ov_ds(spec_i) * avnabpsi
    p_ions = n_ions*T_ions
    dp_ions_ov_dr = T_ions * dn_ions_ov_dr + n_ions * dT_ions_ov_dr
    !
    IF (isw_Vphi_loc .EQ. 0) THEN
       denom_Er = c *  bcovar_tht / sqrtg_bctrvr_phi
       !PRINT *,'denom_Er - 1st term: ',denom_Er
       DO ispec_ctr = LBOUND(row_ind_ptr,1),UBOUND(row_ind_ptr,1)
          irow_spec = row_ind_ptr(ispec_ctr)
          IF (irow_spec .NE. spec_i) CYCLE
          !
          icol_spec = col_ind_ptr(ispec_ctr)
          denom_Er_1 = &
               D31AX_spec(ispec_ctr) * (z_spec(icol_spec)*e) / T_spec(icol_spec)
          denom_Er = denom_Er + denom_Er_1
          !PRINT *, 'denom_Er - 2nd term: ', denom_Er_1
          !
       END DO
    ELSE IF (isw_Vphi_loc .GE. 1) THEN
       fac1 = (hctrvr_tmp(2)*bmod_tmp*1.0e4_dp) * &
            (1.0_dp+TWOPI*aiota_loc*boozer_psi_pr*G_symm_tb) / avb2
       denom_Er = (c / (avnabpsi*aiota_loc*boozer_psi_pr)) + &
            fac1 * (-c*bcovar_phi/sqrtg_bctrvr_tht)
       !PRINT *,'denom_Er - 1st term: ',denom_Er
       DO ispec_ctr = LBOUND(row_ind_ptr,1),UBOUND(row_ind_ptr,1)
          irow_spec = row_ind_ptr(ispec_ctr)
          IF (irow_spec .NE. spec_i) CYCLE
          !
          icol_spec = col_ind_ptr(ispec_ctr)
          denom_Er_1 = &
               D31AX_spec(ispec_ctr) * (z_spec(icol_spec)*e) / T_spec(icol_spec)
          denom_Er = denom_Er + fac1 * denom_Er_1
          !PRINT *, 'denom_Er - 2nd term: ', denom_Er_1
          !
       END DO
    END IF
    !
    IF (isw_Vphi_loc .EQ. 0) THEN
       nom_Er_1 = Vphi * (aiota_loc * bcovar_tht + bcovar_phi)
       nom_Er_2 = (c * T_ions * bcovar_tht / (z_ions * e * sqrtg_bctrvr_phi)) * &
            (dp_ions_ov_dr / p_ions)
       nom_Er = nom_Er_1 + nom_Er_2
       !PRINT *,'Er - 1st term: ', nom_Er_1 / denom_Er
       !PRINT *,Vphi, aiota_loc, bcovar_tht, bcovar_phi, denom_Er
       DO ispec_ctr = LBOUND(row_ind_ptr,1),UBOUND(row_ind_ptr,1)
          irow_spec = row_ind_ptr(ispec_ctr)
          IF (irow_spec .NE. spec_i) CYCLE
          !
          icol_spec = col_ind_ptr(ispec_ctr)
          nom_Er_3 = avnabpsi * D31AX_spec(ispec_ctr) * &
               (dn_spec_ov_ds(icol_spec) / n_spec(icol_spec) + &
               dT_spec_ov_ds(icol_spec) / T_spec(icol_spec))
          nom_Er_4 = avnabpsi * (dT_spec_ov_ds(icol_spec) / T_spec(icol_spec)) * &
               (D32AX_spec(ispec_ctr) - 2.5_dp * D31AX_spec(ispec_ctr))
          nom_Er_5 = D33AX_spec(ispec_ctr) * avEparB_ov_avb2 * &
               (z_spec(icol_spec)*e) / T_spec(icol_spec)
          nom_Er = nom_Er + nom_Er_3 + nom_Er_4 + nom_Er_5
          !PRINT *,'Er - 2nd term: ', nom_Er_3 / denom_Er
          !PRINT *,'Er - 3rd term: ', nom_Er_4 / denom_Er
          !PRINT *,'Er - 4th term: ', nom_Er_5 / denom_Er
          !
       END DO
    ELSE IF (isw_Vphi_loc .GE. 1) THEN
       nom_Er_1 = ((c*T_ions/(z_ions*e)) / (avnabpsi*aiota_loc*boozer_psi_pr)) * &
            (dp_ions_ov_dr / p_ions)
       nom_Er_2 = (c*T_ions/(z_ions*e)) * (-fac1*bcovar_phi/sqrtg_bctrvr_tht) * &
            (dp_ions_ov_dr / p_ions)
       nom_Er = Vphi + nom_Er_1 + nom_Er_2
       !PRINT *,'Er - 1st term: ', nom_Er_1 / denom_Er
       !PRINT *,Vphi, aiota_loc, bcovar_tht, bcovar_phi, denom_Er
       DO ispec_ctr = LBOUND(row_ind_ptr,1),UBOUND(row_ind_ptr,1)
          irow_spec = row_ind_ptr(ispec_ctr)
          IF (irow_spec .NE. spec_i) CYCLE
          !
          icol_spec = col_ind_ptr(ispec_ctr)
          nom_Er_3 = avnabpsi * D31AX_spec(ispec_ctr) * &
               (dn_spec_ov_ds(icol_spec) / n_spec(icol_spec) + &
               dT_spec_ov_ds(icol_spec) / T_spec(icol_spec))
          nom_Er_4 = avnabpsi * (dT_spec_ov_ds(icol_spec) / T_spec(icol_spec)) * &
               (D32AX_spec(ispec_ctr) - 2.5_dp * D31AX_spec(ispec_ctr))
          nom_Er_5 = D33AX_spec(ispec_ctr) * avEparB_ov_avb2 * &
               (z_spec(icol_spec)*e) / T_spec(icol_spec)
          nom_Er = nom_Er + fac1 * (nom_Er_3 + nom_Er_4 + nom_Er_5)
          !PRINT *,'Er - 2nd term: ', nom_Er_3 / denom_Er
          !PRINT *,'Er - 3rd term: ', nom_Er_4 / denom_Er
          !PRINT *,'Er - 4th term: ', nom_Er_5 / denom_Er
          !
       END DO
    END IF
    !
    ! compute radial electric field
    Er = nom_Er / denom_Er
    !
    ! compute species Mach numbers
    IF (ALLOCATED(MtOvR_spec)) DEALLOCATE(MtOvR_spec)
    ALLOCATE(MtOvR_spec(0:num_spec-1))
    MtOvR_spec = (c * Er / (aiota_loc * sqrtg_bctrvr_phi)) / &
         SQRT(2.0_dp * T_spec / m_spec)

    if (allocated(D33AX_spec)) deallocate(D33AX_spec)
  END SUBROUTINE compute_Er
  !
  SUBROUTINE compute_A3norm_a(row_ind_ptr, col_ind_ptr, D31AX_spec, D32AX_spec, &
       D33AX_spec, Er, avEparB_ov_avb2)
    !
    USE neo_precision
    USE neo_control, ONLY: lab_swi
    USE device_mod, ONLY : surface
    USE mag_interface_mod, ONLY : mag_coordinates, &
         boozer_s, boozer_theta_beg, boozer_phi_beg
    USE partpa_mod,  ONLY : bmod0
    USE mag_sub, ONLY: mag
    use neo_magfie, only : boozer_curr_pol_hat, boozer_psi_pr_hat, &
         boozer_curr_pol_hat_s, boozer_curr_tor_hat_s, boozer_isqrg
    USE collisionality_mod, ONLY : num_spec, species_tag, &
         z_spec, m_spec, n_spec, T_spec
    !
    ! ---------------------------------------------------------------!
    ! input:
    ! ---------------------------------------------------------------!
    ! row- and column indices (=species) of diffusion coefficients
    INTEGER, DIMENSION(:), INTENT(in) :: row_ind_ptr, col_ind_ptr
    ! species parallel flow
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: D31AX_spec, D32AX_spec
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: D33AX_spec
    ! radial electric field (w.r.t. effective radius)
    REAL(kind=dp), INTENT(in) :: Er
    ! ---------------------------------------------------------------!
    ! output:
    ! ---------------------------------------------------------------!
    ! drive A_3 normalized (=inductive electric field)
    REAL(kind=dp), INTENT(out) :: avEparB_ov_avb2
    ! ---------------------------------------------------------------!
    ! local:
    ! ---------------------------------------------------------------!
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: y
    REAL(kind=dp) :: aiota_loc, avnabpsi
    ! ---------------------------------------------------------------!
    ! co- and contra-variant B-field components using
    ! r_eff as a flux-surface label
    REAL(kind=dp) :: bcovar_tht, bcovar_phi
    REAL(kind=dp) :: dbcovar_theta_ds, dbcovar_phi_ds
    REAL(kind=dp) :: sqrtg_bctrvr_phi, boozer_psi_pr
    REAL(kind=dp) :: avbhat2, avb2
    ! normalized co-variant phi-component of B,  $\sqrt{g}B^\vartheta$
    REAL(kind=dp) :: bcovar_phi_hat, sqrtg_bctrvr_tht
    REAL(kind=dp) :: dbcovar_theta_hat_ds, dbcovar_phi_hat_ds
    ! Only used to call mag for normalizations
    ! related to D31 and D32
    REAL(kind=dp)                 :: bmod_tmp,sqrtg_tmp
    REAL(kind=dp), DIMENSION(3)   :: x_tmp,bder_tmp,hcovar_tmp,hctrvr_tmp
    REAL(kind=dp), DIMENSION(3,3) :: hcoder_tmp,hctder_tmp
    ! Boozer coordinates of local Vphi value (only used for isw_Vphi_loc=1)
    REAL(kind=dp)                 :: thetaB
    ! transformation function Boozer coord. -> Symm. flux coord.
    REAL(kind=dp)                 :: G_symm, G_symm_tb, G_symm_pb
    ! ---------------------------------------------------------------!
    ! loop indices, temporary variables
    INTEGER :: ispec_ctr, spec_i, irow_spec, icol_spec, num_ctr
    REAL(kind=dp) :: denom_Epar, nom_Epar
    REAL(kind=dp) :: denom_Epar_1, nom_Epar_1
    REAL(kind=dp) :: av_jpar_tot_B
    REAL(kind=dp) :: A1_b, A2_b
    ! ---------------------------------------------------------------!
    !
    ! copy y-vector (see definition in rhs_kin.f90)
    ALLOCATE(y(SIZE(y_ntv_mod,1)))
    y = y_ntv_mod
    !
    ! aiota, avnabpsi, avbhat2:
    aiota_loc = surface%aiota
    avnabpsi = y(7) / y(6)
    avbhat2 = y(9) / y(6)
    !
    ! computation of the normalization for D31 and D32 (-> D31_ref)
    ! + total parallel current (at the moment only available for Boozer coordinates)
    IF (mag_coordinates .EQ. 0) THEN
       ! cylindrical coordinates
       !
       PRINT *,"ntv_mod.f90: Computation of inductive electric field &
            &only implemented for Boozer coordinates at the moment!"
       STOP
       !
       x_tmp = xstart_cyl
       CALL mag(x_tmp,bmod_tmp,sqrtg_tmp,bder_tmp,hcovar_tmp,&
            hctrvr_tmp,hcoder_tmp,hctder_tmp)
       !
       ! normalized co-variant phi-component of B
       ! (Note! co-variant phi-component of B is the
       ! same for cylindrical coordinates and symmetry flux
       ! coodrinates --> no conversion needed)
       bcovar_phi_hat = hcovar_tmp(2)*(bmod_tmp/bmod0)
       !
       ! restore value of $\sqrt{g}B^\vartheta$ for
       ! symmetry flux coordinates from quantities
       ! given in cylindircal coordinates
       sqrtg_bctrvr_tht = avnabpsi*sqrtg_tmp*(hctrvr_tmp(3)*bmod_tmp*1.0e4_dp)
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
       dbcovar_theta_hat_ds = boozer_curr_tor_hat_s
       dbcovar_phi_hat_ds = boozer_curr_pol_hat_s
       !
       IF (lab_swi .EQ. 10) THEN ! ASDEX-U (E. Strumberger)
          ! this is the same up to a minus sign resulting from the
          ! definition of sqrtg_tmp (left-handed system)
          sqrtg_bctrvr_tht = avnabpsi*sqrtg_tmp*(hctrvr_tmp(3)*bmod_tmp*1.0d4) 
       ELSE
          ! restore value of $\sqrt{g}B^\vartheta$ for
          ! symmetry flux coordinates from the quantities
          ! given in Boozer coordinates (only valid for right-handed system)
          sqrtg_bctrvr_tht = avnabpsi*aiota_loc*boozer_psi_pr_hat*(bmod0*1.0e4_dp)
       END IF
    END IF
    !
    ! compute additional B-field quantities
    sqrtg_bctrvr_phi = sqrtg_bctrvr_tht/aiota_loc
    bcovar_phi = hcovar_tmp(2)*(bmod_tmp*1.0e4_dp)
    bcovar_tht = hcovar_tmp(3)*(bmod_tmp*1.0e4_dp)
    dbcovar_phi_ds = dbcovar_phi_hat_ds*(bmod0*1.0e4_dp)
    dbcovar_theta_ds = dbcovar_theta_hat_ds*(bmod0*1.0e4_dp)
    boozer_psi_pr = boozer_psi_pr_hat*(bmod0*1.0e4_dp)
    avb2 = avbhat2*((bmod0*1.0e4_dp)**2)
    !
    ! compute <jE_par*B>/(<E_par*B>/<B^2>)
    denom_Epar = 0.0_dp
    DO ispec_ctr = LBOUND(row_ind_ptr,1),UBOUND(row_ind_ptr,1) ! sum over all species
       irow_spec = row_ind_ptr(ispec_ctr)
       icol_spec = col_ind_ptr(ispec_ctr)
       denom_Epar_1 = - D33AX_spec(ispec_ctr) * &
            (z_spec(irow_spec)*e) * n_spec(irow_spec) * &
            (z_spec(icol_spec)*e) / T_spec(icol_spec)
       denom_Epar = denom_Epar + denom_Epar_1
       !
    END DO
    PRINT *,'<jE_par*B>/(<E_par*B>/<B^2>): ',denom_Epar
    !
    ! compute <jB_par*B>
    nom_Epar = 0.0_dp
    DO ispec_ctr = LBOUND(row_ind_ptr,1),UBOUND(row_ind_ptr,1) ! sum over all species
       !
       irow_spec = row_ind_ptr(ispec_ctr)
       icol_spec = col_ind_ptr(ispec_ctr)
       !
       A2_b = avnabpsi * (dT_spec_ov_ds(icol_spec) / T_spec(icol_spec))
       A1_b = avnabpsi * (dn_spec_ov_ds(icol_spec) / n_spec(icol_spec)) - &
            1.5_dp * A2_b - Er * (z_spec(icol_spec)*e) / T_spec(icol_spec)
       !
       nom_Epar_1 = - (z_spec(irow_spec)*e) * n_spec(irow_spec) * &
            (D31AX_spec(ispec_ctr)*A1_b + D32AX_spec(ispec_ctr)*A2_b)
       nom_Epar = nom_Epar + nom_Epar_1
       !
    END DO
    PRINT *,'<jB_par*B>: ',nom_Epar
    !
    ! compute total parallel current
    av_jpar_tot_B = ( (avb2*c) / (4*PI) ) * SIGN(1.0_dp,boozer_isqrg) * &
         (bcovar_phi*dbcovar_theta_ds - bcovar_tht*dbcovar_phi_ds) / &
         ABS(boozer_psi_pr*(bcovar_phi + aiota_loc * bcovar_tht))
    PRINT *,'av_jpar_tot_B: ',av_jpar_tot_B
    !
    !PRINT *,(avb2 / (4*PI*boozer_psi_pr)) * c
    !PRINT *,avb2, 4*PI, boozer_psi_pr, c
    !PRINT *,(bcovar_phi*dbcovar_theta_ds - bcovar_tht*dbcovar_phi_ds)
    !PRINT *,(bcovar_phi + aiota_loc * bcovar_tht)
    !
    ! compute inductive electric field
    nom_Epar = av_jpar_tot_B - nom_Epar
    avEparB_ov_avb2 = nom_Epar / denom_Epar
    !
  END SUBROUTINE compute_A3norm_a
  !
  SUBROUTINE compute_Er_and_A3norm_a(row_ind_ptr, col_ind_ptr, D31AX_spec, D32AX_spec, &
       D33AX_spec, Er, avEparB_ov_avb2)
    !
    USE neo_precision
    USE neo_control, ONLY: lab_swi
    USE device_mod, ONLY : surface
    USE mag_interface_mod, ONLY : mag_coordinates, &
         boozer_s, boozer_theta_beg, boozer_phi_beg
    USE partpa_mod,  ONLY : bmod0
    USE mag_sub, ONLY: mag
    use neo_magfie, only : boozer_curr_pol_hat, boozer_psi_pr_hat, &
         boozer_curr_pol_hat_s, boozer_curr_tor_hat_s, &
         compute_Gsymm, calc_thetaB_RZloc, boozer_isqrg
    USE collisionality_mod, ONLY : num_spec, species_tag, &
         z_spec, m_spec, n_spec, T_spec
    !
    ! ---------------------------------------------------------------!
    ! input:
    ! ---------------------------------------------------------------!
    ! row- and column indices (=species) of diffusion coefficients
    INTEGER, DIMENSION(:), INTENT(in) :: row_ind_ptr, col_ind_ptr
    ! species parallel flow
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: D31AX_spec, D32AX_spec
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: D33AX_spec
    ! ---------------------------------------------------------------!
    ! output:
    ! ---------------------------------------------------------------!
    ! radial electric field (w.r.t. effective radius) +
    ! drive A_3 normalized (=inductive electric field)
    REAL(kind=dp), INTENT(out) :: Er, avEparB_ov_avb2
    ! ---------------------------------------------------------------!
    ! local:
    ! ---------------------------------------------------------------!
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: y
    REAL(kind=dp) :: aiota_loc, avnabpsi
    ! ---------------------------------------------------------------!
    ! co- and contra-variant B-field components using
    ! r_eff as a flux-surface label
    REAL(kind=dp) :: bcovar_tht, bcovar_phi
    REAL(kind=dp) :: dbcovar_theta_ds, dbcovar_phi_ds
    REAL(kind=dp) :: sqrtg_bctrvr_phi, boozer_psi_pr
    REAL(kind=dp) :: avbhat2, avb2
    ! normalized co-variant phi-component of B,  $\sqrt{g}B^\vartheta$
    REAL(kind=dp) :: bcovar_phi_hat, sqrtg_bctrvr_tht
    REAL(kind=dp) :: dbcovar_theta_hat_ds, dbcovar_phi_hat_ds
    ! Only used to call mag for normalizations
    ! related to D31 and D32
    REAL(kind=dp)                 :: bmod_tmp,sqrtg_tmp
    REAL(kind=dp), DIMENSION(3)   :: x_tmp,bder_tmp,hcovar_tmp,hctrvr_tmp
    REAL(kind=dp), DIMENSION(3,3) :: hcoder_tmp,hctder_tmp
    ! Boozer coordinates of local Vphi value (only used for isw_Vphi_loc=1)
    REAL(kind=dp)                 :: thetaB
    REAL(kind=dp), DIMENSION(3)   :: x_start
    ! transformation function Boozer coord. -> Symm. flux coord.
    REAL(kind=dp)                 :: G_symm, G_symm_tb, G_symm_pb
    ! ---------------------------------------------------------------!
    ! temperature, pressure and density of species i 
    ! (i = species of measured toroidal rotation frequency)
    REAL(kind=dp) :: T_ions, n_ions, p_ions, z_ions
    REAL(kind=dp) :: dT_ions_ov_dr, dn_ions_ov_dr, dp_ions_ov_dr
    ! loop indices, temporary variables
    INTEGER :: ispec_ctr, spec_i, irow_spec, icol_spec, num_ctr
    REAL(kind=dp) :: denom_Er_a, denom_Er_c, fac1
    REAL(kind=dp) :: denom_Er_a_1, denom_Er_c_1
    REAL(kind=dp) :: nom_Er
    REAL(kind=dp) :: nom_Er_1, nom_Er_2, nom_Er_3, nom_Er_4, nom_Er_5
    REAL(kind=dp) :: denom_Epar_b, denom_Epar_d, nom_Epar
    REAL(kind=dp) :: denom_Epar_b_1, denom_Epar_d_1, nom_Epar_1
    REAL(kind=dp) :: av_jpar_tot_B
    REAL(kind=dp) :: dlogT_ov_dr, dlogn_ov_dr
    ! ---------------------------------------------------------------!
    !
    ! copy y-vector (see definition in rhs_kin.f90)
    ALLOCATE(y(SIZE(y_ntv_mod,1)))
    y = y_ntv_mod
    !
    ! aiota, avnabpsi, avbhat2:
    aiota_loc = surface%aiota
    avnabpsi = y(7) / y(6)
    avbhat2 = y(9) / y(6)

    spec_i = -1

    denom_Er_a = 0.0
    denom_Er_c = 0.0
    nom_Er = 0.0
    denom_Epar_b = 0.0
    denom_Epar_d = 0.0
    nom_Epar = 0.0

    ! computation of the normalization for D31 and D32 (-> D31_ref)
    ! + total parallel current (at the moment only available for Boozer coordinates)
    IF (mag_coordinates .EQ. 0) THEN
       ! cylindrical coordinates
       !
       PRINT *,"ntv_mod.f90: Computation of Er and inductive electric field &
            &only implemented for Boozer coordinates at the moment!"
       STOP
       !
       x_tmp = xstart_cyl
       CALL mag(x_tmp,bmod_tmp,sqrtg_tmp,bder_tmp,hcovar_tmp,&
            hctrvr_tmp,hcoder_tmp,hctder_tmp)
       !
       ! normalized co-variant phi-component of B
       ! (Note! co-variant phi-component of B is the
       ! same for cylindrical coordinates and symmetry flux
       ! coodrinates --> no conversion needed)
       bcovar_phi_hat = hcovar_tmp(2)*(bmod_tmp/bmod0)
       !
       ! restore value of $\sqrt{g}B^\vartheta$ for
       ! symmetry flux coordinates from quantities
       ! given in cylindircal coordinates
       sqrtg_bctrvr_tht = avnabpsi*sqrtg_tmp*(hctrvr_tmp(3)*bmod_tmp*1.0e4_dp)
    ELSE
       ! boozer coordinates
       !
       IF (isw_Vphi_loc .EQ. 0) THEN
          ! isw_Vphi_loc=0: value of "Vphi" corresponds to
          ! flux surface average (<V_\varphi>)
          x_tmp = (/boozer_s,boozer_phi_beg,boozer_theta_beg/)
          CALL mag(x_tmp,bmod_tmp,sqrtg_tmp,bder_tmp,hcovar_tmp,&
               hctrvr_tmp,hcoder_tmp,hctder_tmp)
       ELSE IF (isw_Vphi_loc .EQ. 1) THEN
          ! Caution: This branch is not tested!
          ! isw_Vphi_loc=1: value of "Vphi" is specified
          ! locally for given (R,Z)-position
          x_start = (/boozer_s,boozer_phi_beg,0.0_dp/)
          CALL calc_thetaB_RZloc(R_Vphi, Z_Vphi, x_start, thetaB)
          x_tmp = (/boozer_s,boozer_phi_beg,thetaB/)
          CALL mag(x_tmp,bmod_tmp,sqrtg_tmp,bder_tmp,hcovar_tmp,&
               hctrvr_tmp,hcoder_tmp,hctder_tmp)
          CALL compute_Gsymm( x_tmp, G_symm, G_symm_tb, G_symm_pb )
       ELSE IF (isw_Vphi_loc .EQ. 2) THEN
          ! isw_Vphi_loc=2: value of "Vphi" is specified
          ! locally for given \vartheta_B position
          x_tmp = (/boozer_s,boozer_phi_beg,boozer_theta_Vphi/)
          CALL mag(x_tmp,bmod_tmp,sqrtg_tmp,bder_tmp,hcovar_tmp,&
               hctrvr_tmp,hcoder_tmp,hctder_tmp)
          CALL compute_Gsymm( x_tmp, G_symm, G_symm_tb, G_symm_pb )
       ELSE
          PRINT *,"ntv_mod.f90: Undefined state of switch isw_Vphi_loc (= 0 / 1 / 2)!"
          STOP
       END IF
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
       dbcovar_theta_hat_ds = boozer_curr_tor_hat_s
       dbcovar_phi_hat_ds = boozer_curr_pol_hat_s
       !
       IF (lab_swi .EQ. 10) THEN ! ASDEX-U (E. Strumberger)
          ! this is the same up to a minus sign resulting from the
          ! definition of sqrtg_tmp (left-handed system)
          sqrtg_bctrvr_tht = avnabpsi*sqrtg_tmp*(hctrvr_tmp(3)*bmod_tmp*1.0d4) 
       ELSE
          ! restore value of $\sqrt{g}B^\vartheta$ for
          ! symmetry flux coordinates from the quantities
          ! given in Boozer coordinates (only valid for right-handed system)
          sqrtg_bctrvr_tht = avnabpsi*aiota_loc*boozer_psi_pr_hat*(bmod0*1.0e4_dp)
       END IF
    END IF
    !
    ! compute additional B-field quantities
    sqrtg_bctrvr_phi = sqrtg_bctrvr_tht/aiota_loc
    bcovar_phi = hcovar_tmp(2)*(bmod_tmp*1.0e4_dp)
    bcovar_tht = hcovar_tmp(3)*(bmod_tmp*1.0e4_dp)
    dbcovar_phi_ds = dbcovar_phi_hat_ds*(bmod0*1.0e4_dp)
    dbcovar_theta_ds = dbcovar_theta_hat_ds*(bmod0*1.0e4_dp)
    boozer_psi_pr = boozer_psi_pr_hat*(bmod0*1.0e4_dp)
    avb2 = avbhat2*((bmod0*1.0e4_dp)**2)
    !
    ! detect species of measured V_phi
    num_ctr = 0
    DO ispec_ctr = 0,num_spec-1
       IF (species_tag(ispec_ctr) .EQ. species_tag_Vphi) THEN
          spec_i = ispec_ctr
          num_ctr = num_ctr + 1
       END IF
    END DO
    IF (num_ctr .EQ. 0) THEN
       PRINT *,"ntv_mod.f90: Subroutine compute_Er - &
            &Species of measured V_phi not found!"
       STOP
    ELSEIF (num_ctr .GT. 1) THEN
       PRINT *,"ntv_mod.f90: Subroutine compute_Er - &
            &Multiple definition of species V_phi!"
       STOP
    END IF
    !
    z_ions = z_spec(spec_i)
    T_ions = T_spec(spec_i)
    dT_ions_ov_dr = dT_spec_ov_ds(spec_i) * avnabpsi
    n_ions = n_spec(spec_i)
    dn_ions_ov_dr = dn_spec_ov_ds(spec_i) * avnabpsi
    p_ions = n_ions*T_ions
    dp_ions_ov_dr = T_ions * dn_ions_ov_dr + n_ions * dT_ions_ov_dr
    !
    ! compute <jE_par*B>/(<E_par*B>/<B^2>) and denom_Er
    IF (isw_Vphi_loc .EQ. 0) THEN
       denom_Er_a = c *  bcovar_tht / sqrtg_bctrvr_phi
       denom_Er_c = 0.0_dp
       denom_Epar_b = 0.0_dp
       denom_Epar_d = 0.0_dp
       DO ispec_ctr = LBOUND(row_ind_ptr,1),UBOUND(row_ind_ptr,1) ! sum over all species
          irow_spec = row_ind_ptr(ispec_ctr)
          icol_spec = col_ind_ptr(ispec_ctr)
          !
          IF (irow_spec .EQ. spec_i) THEN
             denom_Er_a_1 = &
                  D31AX_spec(ispec_ctr) * (z_spec(icol_spec)*e) / T_spec(icol_spec)
             denom_Er_a = denom_Er_a + denom_Er_a_1
             !
             denom_Epar_b_1 = &
                  - D33AX_spec(ispec_ctr) * (z_spec(icol_spec)*e) / T_spec(icol_spec)
             denom_Epar_b = denom_Epar_b + denom_Epar_b_1
          END IF
          !
          denom_Er_c_1 = D31AX_spec(ispec_ctr) * &
               (z_spec(irow_spec)*e) * n_spec(irow_spec) * &
               (z_spec(icol_spec)*e) / T_spec(icol_spec)
          denom_Er_c = denom_Er_c + denom_Er_c_1
          !
          denom_Epar_d_1 = - D33AX_spec(ispec_ctr) * &
               (z_spec(irow_spec)*e) * n_spec(irow_spec) * &
               (z_spec(icol_spec)*e) / T_spec(icol_spec)
          denom_Epar_d = denom_Epar_d + denom_Epar_d_1
          !
       END DO
    ELSE IF (isw_Vphi_loc .GE. 1) THEN
       fac1 = (hctrvr_tmp(2)*bmod_tmp*1.0e4_dp) * &
            (1.0_dp+TWOPI*aiota_loc*boozer_psi_pr*G_symm_tb) / avb2
       denom_Er_a = (c / (avnabpsi*aiota_loc*boozer_psi_pr)) + &
            fac1 * (-c*bcovar_phi/sqrtg_bctrvr_tht)
       denom_Er_c = 0.0_dp
       denom_Epar_b = 0.0_dp
       denom_Epar_d = 0.0_dp
       DO ispec_ctr = LBOUND(row_ind_ptr,1),UBOUND(row_ind_ptr,1)
          irow_spec = row_ind_ptr(ispec_ctr)
          icol_spec = col_ind_ptr(ispec_ctr)
          !
          IF (irow_spec .EQ. spec_i) THEN
             denom_Er_a_1 = &
                  D31AX_spec(ispec_ctr) * (z_spec(icol_spec)*e) / T_spec(icol_spec)
             denom_Er_a = denom_Er_a + fac1 * denom_Er_a_1
             !
             denom_Epar_b_1 = &
                  - D33AX_spec(ispec_ctr) * (z_spec(icol_spec)*e) / T_spec(icol_spec)
             denom_Epar_b = denom_Epar_b + fac1 * denom_Epar_b_1
          END IF
          !
          denom_Er_c_1 = D31AX_spec(ispec_ctr) * &
               (z_spec(irow_spec)*e) * n_spec(irow_spec) * &
               (z_spec(icol_spec)*e) / T_spec(icol_spec)
          denom_Er_c = denom_Er_c + denom_Er_c_1
          !
          denom_Epar_d_1 = - D33AX_spec(ispec_ctr) * &
               (z_spec(irow_spec)*e) * n_spec(irow_spec) * &
               (z_spec(icol_spec)*e) / T_spec(icol_spec)
          denom_Epar_d = denom_Epar_d + denom_Epar_d_1
          !
       END DO
    END IF
    !
    ! compute <jB_par*B> and nom_Er
    IF (isw_Vphi_loc .EQ. 0) THEN
       nom_Er_1 = Vphi * (aiota_loc * bcovar_tht + bcovar_phi)
       nom_Er_2 = (c * T_ions * bcovar_tht / (z_ions * e * sqrtg_bctrvr_phi)) * &
            (dp_ions_ov_dr / p_ions)
       nom_Er = nom_Er_1 + nom_Er_2
       nom_Epar = 0.0_dp
       DO ispec_ctr = LBOUND(row_ind_ptr,1),UBOUND(row_ind_ptr,1) ! sum over all species
          !
          irow_spec = row_ind_ptr(ispec_ctr)
          icol_spec = col_ind_ptr(ispec_ctr)
          !
          IF (irow_spec .EQ. spec_i) THEN
             nom_Er_3 = avnabpsi * D31AX_spec(ispec_ctr) * &
                  (dn_spec_ov_ds(icol_spec) / n_spec(icol_spec) + &
                  dT_spec_ov_ds(icol_spec) / T_spec(icol_spec))
             nom_Er_4 = avnabpsi * (dT_spec_ov_ds(icol_spec) / T_spec(icol_spec)) * &
                  (D32AX_spec(ispec_ctr) - 2.5_dp * D31AX_spec(ispec_ctr))
             nom_Er = nom_Er + nom_Er_3 + nom_Er_4
          END IF
          !
          dlogT_ov_dr = avnabpsi * (dT_spec_ov_ds(icol_spec) / T_spec(icol_spec))
          dlogn_ov_dr = avnabpsi * (dn_spec_ov_ds(icol_spec) / n_spec(icol_spec))
          !
          nom_Epar_1 = (z_spec(irow_spec)*e) * n_spec(irow_spec) * &
               (D31AX_spec(ispec_ctr)*(dlogn_ov_dr+dlogT_ov_dr) + &
               (D32AX_spec(ispec_ctr)-2.5_dp*D31AX_spec(ispec_ctr)) * &
               dlogT_ov_dr)
          nom_Epar = nom_Epar + nom_Epar_1
          !
       END DO
    ELSE IF (isw_Vphi_loc .GE. 1) THEN
       nom_Er_1 = ((c*T_ions/(z_ions*e)) / (avnabpsi*aiota_loc*boozer_psi_pr)) * &
            (dp_ions_ov_dr / p_ions)
       nom_Er_2 = (c*T_ions/(z_ions*e)) * (-fac1*bcovar_phi/sqrtg_bctrvr_tht) * &
            (dp_ions_ov_dr / p_ions)
       nom_Er = Vphi + nom_Er_1 + nom_Er_2
       nom_Epar = 0.0_dp
       DO ispec_ctr = LBOUND(row_ind_ptr,1),UBOUND(row_ind_ptr,1)
          !
          irow_spec = row_ind_ptr(ispec_ctr)
          icol_spec = col_ind_ptr(ispec_ctr)
          !
          IF (irow_spec .EQ. spec_i) THEN
             nom_Er_3 = avnabpsi * D31AX_spec(ispec_ctr) * &
                  (dn_spec_ov_ds(icol_spec) / n_spec(icol_spec) + &
                  dT_spec_ov_ds(icol_spec) / T_spec(icol_spec))
             nom_Er_4 = avnabpsi * (dT_spec_ov_ds(icol_spec) / T_spec(icol_spec)) * &
                  (D32AX_spec(ispec_ctr) - 2.5_dp * D31AX_spec(ispec_ctr))
             nom_Er = nom_Er + fac1 * (nom_Er_3 + nom_Er_4)
          END IF
          !
          dlogT_ov_dr = avnabpsi * (dT_spec_ov_ds(icol_spec) / T_spec(icol_spec))
          dlogn_ov_dr = avnabpsi * (dn_spec_ov_ds(icol_spec) / n_spec(icol_spec))
          !
          nom_Epar_1 = (z_spec(irow_spec)*e) * n_spec(irow_spec) * &
               (D31AX_spec(ispec_ctr)*(dlogn_ov_dr+dlogT_ov_dr) + &
               (D32AX_spec(ispec_ctr)-2.5_dp*D31AX_spec(ispec_ctr)) * &
               dlogT_ov_dr)
          nom_Epar = nom_Epar + nom_Epar_1
          !
       END DO
    END IF
    !PRINT *,'Er (simple): ',nom_Er/denom_Er_a
    !
    ! compute total parallel current
    av_jpar_tot_B = ( (avb2*c) / (4*PI) ) * SIGN(1.0_dp,boozer_isqrg) * &
         (bcovar_phi*dbcovar_theta_ds - bcovar_tht*dbcovar_phi_ds) / &
         ABS(boozer_psi_pr*(bcovar_phi + aiota_loc * bcovar_tht))
    !PRINT *,'av_jpar_tot_B: ',av_jpar_tot_B
    !
    !PRINT *,(avb2 / (4*PI*boozer_psi_pr)) * c
    !PRINT *,avb2, 4*PI, boozer_psi_pr, c
    !PRINT *,(bcovar_phi*dbcovar_theta_ds - bcovar_tht*dbcovar_phi_ds)
    !PRINT *,(bcovar_phi + aiota_loc * bcovar_tht)
    !    
    nom_Epar = av_jpar_tot_B + nom_Epar
    !
    ! compute inductive electric field and radial electric field
    avEparB_ov_avb2 = (nom_Er*denom_Er_c - nom_Epar*denom_Er_a) / &
         (denom_Epar_b*denom_Er_c - denom_Epar_d*denom_Er_a)
    Er = (nom_Er-denom_Epar_b*avEparB_ov_avb2)/denom_Er_a
    !
    ! compute species Mach numbers
    IF (ALLOCATED(MtOvR_spec)) DEALLOCATE(MtOvR_spec)
    ALLOCATE(MtOvR_spec(0:num_spec-1))
    MtOvR_spec = (c * Er / (aiota_loc * sqrtg_bctrvr_phi)) / &
         SQRT(2.0_dp * T_spec / m_spec)
    !
  END SUBROUTINE compute_Er_and_A3norm_a
  !
  SUBROUTINE get_Er_a(qflux_ab_AX_in, Er)
    !
    USE neo_precision
    USE collisionality_mod, ONLY : num_spec
    !
    ! ---------------------------------------------------------------!
    ! input:
    ! ---------------------------------------------------------------!
    REAL(kind=dp), DIMENSION(:,:,:,:), INTENT(in) :: qflux_ab_AX_in
    ! ---------------------------------------------------------------!
    ! output:
    ! ---------------------------------------------------------------!
    ! radial electric field (w.r.t. effective radius)
    REAL(kind=dp), INTENT(out) :: Er
    ! ---------------------------------------------------------------!
    ! local definitions:
    ! ---------------------------------------------------------------!
    ! row- and column-indices (=species_tag) of
    ! diffusion tensor elements (e.g, D11, D12, ...)
    ! -> map to all species of the given profile (global)
    INTEGER, DIMENSION(:), ALLOCATABLE :: row_ind_spec, col_ind_spec
    ! local row- and column-indices of diffusion tensor elements
    ! (e.g, D11, D12, ...)
    ! -> map to all species of the selected radial point (local: 0:num_spec-1)
    INTEGER, DIMENSION(:), ALLOCATABLE :: row_ind_ptr, col_ind_ptr
    ! species diffusion coefficients:
    REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: Dijab_NA
    REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: Dijab_AX
    REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: Dijab_norm_NA
    REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: Dijab_norm_AX
    REAL(kind=dp), DIMENSION(:,:,:,:), ALLOCATABLE :: qflux_ab_NA_tmp
    ! ---------------------------------------------------------------!
    !
    ! Allocate dummy storage arrays for qflux_ab_NA
    IF (ALLOCATED(qflux_ab_NA_tmp)) DEALLOCATE(qflux_ab_NA_tmp)
    ALLOCATE(qflux_ab_NA_tmp(1:3,1:3,0:num_spec-1,0:num_spec-1))
    qflux_ab_NA_tmp = 0.0_dp
    !
    ! Compute diffusion coefficients
    CALL compute_Dijab(qflux_ab_NA_tmp, qflux_ab_AX_in, &
         row_ind_spec, col_ind_spec, row_ind_ptr, col_ind_ptr, &
         Dijab_NA, Dijab_AX, Dijab_norm_NA, Dijab_norm_AX)
    !
    ! Compute Er
    CALL compute_Er(row_ind_ptr, col_ind_ptr, Dijab_AX(3,1,:), Dijab_AX(3,2,:), Er)
    !
  END SUBROUTINE get_Er_a
  !
  SUBROUTINE get_Er_b(qflux_ab_AX_in, Er, avEparB_ov_avb2)
    !
    USE neo_precision
    USE collisionality_mod, ONLY : num_spec
    !
    ! ---------------------------------------------------------------!
    ! input:
    ! ---------------------------------------------------------------!
    REAL(kind=dp), DIMENSION(:,:,:,:), INTENT(in) :: qflux_ab_AX_in
    ! ---------------------------------------------------------------!
    ! output:
    ! ---------------------------------------------------------------!
    ! radial electric field (w.r.t. effective radius)
    REAL(kind=dp), INTENT(out) :: Er, avEparB_ov_avb2
    ! ---------------------------------------------------------------!
    ! local definitions:
    ! ---------------------------------------------------------------!
    ! relative accuracy, loop indices, counters, ...
    REAL(kind=dp), PARAMETER :: epserr_iter = 1.0e-5_dp
    INTEGER, PARAMETER :: k_max = 100
    LOGICAL :: break_cond
    INTEGER :: k
    REAL(kind=dp) :: Er_prev, avEparB_ov_avb2_prev
    ! row- and column-indices (=species_tag) of
    ! diffusion tensor elements (e.g, D11, D12, ...)
    ! -> map to all species of the given profile (global)
    INTEGER, DIMENSION(:), ALLOCATABLE :: row_ind_spec, col_ind_spec
    ! local row- and column-indices of diffusion tensor elements
    ! (e.g, D11, D12, ...)
    ! -> map to all species of the selected radial point (local: 0:num_spec-1)
    INTEGER, DIMENSION(:), ALLOCATABLE :: row_ind_ptr, col_ind_ptr
    ! species diffusion coefficients:
    REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: Dijab_NA
    REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: Dijab_AX
    REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: Dijab_norm_NA
    REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: Dijab_norm_AX
    REAL(kind=dp), DIMENSION(:,:,:,:), ALLOCATABLE :: qflux_ab_NA_tmp
    ! ---------------------------------------------------------------!
    !
    ! Allocate dummy storage arrays for qflux_ab_NA
    IF (ALLOCATED(qflux_ab_NA_tmp)) DEALLOCATE(qflux_ab_NA_tmp)
    ALLOCATE(qflux_ab_NA_tmp(1:3,1:3,0:num_spec-1,0:num_spec-1))
    qflux_ab_NA_tmp = 0.0_dp
    !
    ! Compute diffusion coefficients
    CALL compute_Dijab(qflux_ab_NA_tmp, qflux_ab_AX_in, &
         row_ind_spec, col_ind_spec, row_ind_ptr, col_ind_ptr, &
         Dijab_NA, Dijab_AX, Dijab_norm_NA, Dijab_norm_AX)
    !
    ! Compute Er (without account of A_3)
    CALL compute_Er(row_ind_ptr, col_ind_ptr, Dijab_AX(3,1,:), Dijab_AX(3,2,:), Er)
    !
    ! Compute <E_par*B>/<B^2>
    CALL compute_A3norm(row_ind_ptr, col_ind_ptr, Dijab_AX(3,1,:), &
         Dijab_AX(3,2,:), Dijab_AX(3,3,:), Er, avEparB_ov_avb2)
    !PAUSE
    !
    ! Determine Er and <E_par*B>/<B^2> iteratively
    break_cond = .FALSE.
    k = 0
    DO WHILE (.NOT. break_cond)
       !
       k = k + 1
       IF (k .GT. k_max) THEN
          PRINT *,"propagator.f90: Subroutine get_Er - &
               &maximum number of iterations reached!"
          STOP
       END IF
       !
       Er_prev = Er
       avEparB_ov_avb2_prev = avEparB_ov_avb2
       !
       ! Compute Er (with account of A_3)
       CALL compute_Er(row_ind_ptr, col_ind_ptr, Dijab_AX(3,1,:), Dijab_AX(3,2,:), &
            Er, Dijab_AX(3,3,:), avEparB_ov_avb2_prev)
       !
       ! Compute <E_par*B>/<B^2>
       CALL compute_A3norm(row_ind_ptr, col_ind_ptr, Dijab_AX(3,1,:), &
            Dijab_AX(3,2,:), Dijab_AX(3,3,:), Er, avEparB_ov_avb2)
       !PAUSE
       !
       break_cond = (ABS(Er-Er_prev) .LT. ABS(Er*epserr_iter)) .AND. &
            (ABS(avEparB_ov_avb2 - avEparB_ov_avb2_prev) .LT. &
            ABS(avEparB_ov_avb2*epserr_iter))
       !
    END DO
    !
  END SUBROUTINE get_Er_b
  !
  SUBROUTINE get_B_rho_L_loc_a()
    !
    USE neo_precision, ONLY : dp
    USE collisionality_mod, ONLY : num_spec, z_spec, m_spec, T_spec
    !
    ! local indices
    INTEGER :: ind_spec
    !
    ! compute species hatOmegaB_ref
    IF (ALLOCATED(B_rho_L_loc_spec)) DEALLOCATE(B_rho_L_loc_spec)
    ALLOCATE(B_rho_L_loc_spec(0:num_spec-1))
    DO ind_spec = 0,num_spec-1
       B_rho_L_loc_spec(ind_spec) = &
            c * SQRT(2.0_dp * m_spec(ind_spec) * T_spec(ind_spec)) / &
            (z_spec(ind_spec) * e)
    END DO
    !
  END SUBROUTINE get_B_rho_L_loc_a
  !
  SUBROUTINE compute_Vphi_profile(row_ind_ptr, col_ind_ptr, &
       & D31AX_spec, D32AX_spec, Er, R_Vphi_prof, Z_Vphi_prof, &
       & Vphi_prof_spec, Vtht_prof_spec, D33AX_spec_in, avEparB_ov_avb2_in)
    !
    USE neo_precision
    USE neo_control, ONLY: lab_swi
    USE device_mod, ONLY : surface
    USE mag_interface_mod, ONLY : mag_coordinates, &
         boozer_s, boozer_theta_beg, boozer_phi_beg
    USE partpa_mod,  ONLY : bmod0
    USE mag_sub, ONLY: mag
    use neo_magfie, only : boozer_curr_pol_hat, boozer_psi_pr_hat, &
         compute_Gsymm, compute_RZ
    USE collisionality_mod, ONLY : num_spec, species_tag, &
         z_spec, m_spec, n_spec, T_spec
    !
    ! ---------------------------------------------------------------!
    ! input:
    ! ---------------------------------------------------------------!
    ! row- and column indices (=species) of diffusion coefficients
    INTEGER, DIMENSION(:), INTENT(in) :: row_ind_ptr, col_ind_ptr
    ! species parallel flow
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: D31AX_spec, D32AX_spec
    real(kind=dp), dimension(:), intent(in), optional :: D33AX_spec_in
    ! radial electric field (w.r.t. effective radius) +
    ! drive A_3 normalized (=inductive electric field)
    REAL(kind=dp), INTENT(in) :: Er
    real(kind=dp), intent(in), optional :: avEparB_ov_avb2_in
    ! ---------------------------------------------------------------!
    ! output:
    ! ---------------------------------------------------------------!
    ! species poloidal variation of toroidal rotation frequency
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(inout) :: R_Vphi_prof
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(inout) :: Z_Vphi_prof
    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: Vphi_prof_spec
    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: Vtht_prof_spec
    ! ---------------------------------------------------------------!
    ! local:
    ! ---------------------------------------------------------------!
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: y
    REAL(kind=dp) :: aiota_loc, avnabpsi
    ! ---------------------------------------------------------------!
    ! co- and contra-variant B-field components using
    ! r_eff as a flux-surface label
    REAL(kind=dp) :: bcovar_tht, bcovar_phi
    REAL(kind=dp) :: sqrtg_bctrvr_phi, boozer_psi_pr
    REAL(kind=dp) :: avbhat2, avb2
    ! normalized co-variant phi-component of B,  $\sqrt{g}B^\vartheta$
    REAL(kind=dp) :: bcovar_phi_hat, sqrtg_bctrvr_tht
    ! Only used to call mag for normalizations
    ! related to D31 and D32
    REAL(kind=dp)                 :: bmod_tmp,sqrtg_tmp
    REAL(kind=dp), DIMENSION(3)   :: x_tmp,bder_tmp,hcovar_tmp,hctrvr_tmp
    REAL(kind=dp), DIMENSION(3,3) :: hcoder_tmp,hctder_tmp
    ! transformation function Boozer coord. -> Symm. flux coord.
    REAL(kind=dp)                 :: G_symm, G_symm_tb, G_symm_pb
    ! cylindrical coordinates
    REAL(kind=dp)                 :: R, R_tb, Z, Z_tb
    ! ---------------------------------------------------------------!
    ! loop indices, temporary variables
    INTEGER :: ispec_ctr, spec_i, irow_spec, icol_spec, irow_ctr
    REAL(kind=dp) :: fac1
    REAL(kind=dp) :: denom_Er_1, nom_Er_3, nom_Er_4, nom_Er_5
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: denom_Er, nom_Er
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: nom_Er_1, nom_Er_2
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: coef_K, coef_K_Er, coef_K_grad
    ! ---------------------------------------------------------------!
    ! poloidal dependence of Vtor: loop indices
    INTEGER:: ind_thtB
    REAL(kind=dp) :: thetaB

    ! Variables for the optional parameters.
    real(kind=dp), dimension(:) , allocatable:: D33AX_spec
    real(kind=dp) :: avEparB_ov_avb2
    ! ---------------------------------------------------------------!

    allocate(D33AX_spec(lbound(row_ind_ptr,1):ubound(row_ind_ptr,1)))
    if (present(D33AX_spec_in)) then
      D33AX_spec = D33AX_spec_in
    else
      D33AX_spec = 0.0
    end if
    if (present(avEparB_ov_avb2_in)) then
      avEparB_ov_avb2 = avEparB_ov_avb2_in
    else
      avEparB_ov_avb2 = 0.0
    end if

    ! copy y-vector (see definition in rhs_kin.f90)
    ALLOCATE(y(SIZE(y_ntv_mod,1)))
    y = y_ntv_mod
    !
    ! aiota, avnabpsi :
    aiota_loc = surface%aiota
    avnabpsi = y(7) / y(6)
    avbhat2 = y(9) / y(6)
    !
    ! computation of the normalization for D31 and D32 (-> D31_ref)
    IF (mag_coordinates .EQ. 0) THEN
       ! cylindrical coordinates
       PRINT *,"ntv_mod.f90: Evaluation of poloidal variation of Vtor&
            &only implemented for Boozer coordinates at the moment!"
       STOP
       !
       x_tmp = xstart_cyl
       CALL mag(x_tmp,bmod_tmp,sqrtg_tmp,bder_tmp,hcovar_tmp,&
            hctrvr_tmp,hcoder_tmp,hctder_tmp)
       !
       ! normalized co-variant phi-component of B
       ! (Note! co-variant phi-component of B is the
       ! same for cylindrical coordinates and symmetry flux
       ! coodrinates --> no conversion needed)
       bcovar_phi_hat = hcovar_tmp(2)*(bmod_tmp/bmod0)
       !
       ! restore value of $\sqrt{g}B^\vartheta$ for
       ! symmetry flux coordinates from quantities
       ! given in cylindircal coordinates
       sqrtg_bctrvr_tht = avnabpsi*sqrtg_tmp*(hctrvr_tmp(3)*bmod_tmp*1.0e4_dp)
    ELSE
       ! boozer coordinates
       !
       IF (isw_Vphi_loc.GE.0 .AND. isw_Vphi_loc.LE.2) THEN
          x_tmp = (/boozer_s,boozer_phi_beg,boozer_theta_beg/)
          CALL mag(x_tmp,bmod_tmp,sqrtg_tmp,bder_tmp,hcovar_tmp,&
               hctrvr_tmp,hcoder_tmp,hctder_tmp)
       ELSE
          PRINT *,"ntv_mod.f90: Undefined state of switch isw_Vphi_loc (= 0 / 1 / 2)!"
          STOP
       END IF
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
       IF (lab_swi .EQ. 10) THEN ! ASDEX-U (E. Strumberger)
          ! this is the same up to a minus sign resulting from the
          ! definition of sqrtg_tmp (left-handed system)
          sqrtg_bctrvr_tht = avnabpsi*sqrtg_tmp*(hctrvr_tmp(3)*bmod_tmp*1.0d4) 
       ELSE
          ! restore value of $\sqrt{g}B^\vartheta$ for
          ! symmetry flux coordinates from the quantities
          ! given in Boozer coordinates (only valid for right-handed system)
          sqrtg_bctrvr_tht = avnabpsi*aiota_loc*boozer_psi_pr_hat*(bmod0*1.0e4_dp)
       END IF
    END IF
    !
    IF (ALLOCATED(R_Vphi_prof)) DEALLOCATE(R_Vphi_prof)
    ALLOCATE(R_Vphi_prof(1:num_thtB_VphiProf))
    R_Vphi_prof = 0.0_dp
    IF (ALLOCATED(Z_Vphi_prof)) DEALLOCATE(Z_Vphi_prof)
    ALLOCATE(Z_Vphi_prof(1:num_thtB_VphiProf))
    Z_Vphi_prof = 0.0_dp
    IF (ALLOCATED(Vphi_prof_spec)) DEALLOCATE(Vphi_prof_spec)
    ALLOCATE(Vphi_prof_spec(0:num_spec-1,1:num_thtB_VphiProf))
    Vphi_prof_spec = 0.0_dp
    !
    IF (ALLOCATED(denom_Er)) DEALLOCATE(denom_Er)
    ALLOCATE(denom_Er(0:num_spec-1))
    denom_Er = 0.0_dp
    IF (ALLOCATED(nom_Er)) DEALLOCATE(nom_Er)
    ALLOCATE(nom_Er(0:num_spec-1))
    nom_Er = 0.0_dp
    IF (ALLOCATED(nom_Er_1)) DEALLOCATE(nom_Er_1)
    ALLOCATE(nom_Er_1(0:num_spec-1))
    nom_Er_1 = 0.0_dp
    IF (ALLOCATED(nom_Er_2)) DEALLOCATE(nom_Er_2)
    ALLOCATE(nom_Er_2(0:num_spec-1))
    nom_Er_2 = 0.0_dp
    !
    IF (ALLOCATED(Vtht_prof_spec)) DEALLOCATE(Vtht_prof_spec)
    ALLOCATE(Vtht_prof_spec(0:num_spec-1,1:num_thtB_VphiProf))
    Vtht_prof_spec = 0.0_dp
    !
    IF (ALLOCATED(coef_K)) DEALLOCATE(coef_K)
    ALLOCATE(coef_K(0:num_spec-1))
    coef_K = 0.0_dp
    IF (ALLOCATED(coef_K_Er)) DEALLOCATE(coef_K_Er)
    ALLOCATE(coef_K_Er(0:num_spec-1))
    coef_K_Er = 0.0_dp
    IF (ALLOCATED(coef_K_grad)) DEALLOCATE(coef_K_grad)
    ALLOCATE(coef_K_grad(0:num_spec-1))
    coef_K_grad = 0.0_dp
    !
    DO ind_thtB = 1,num_thtB_VphiProf
       !
       thetaB = DBLE(ind_thtB-1)*TWOPI/DBLE(num_thtB_VphiProf-1)
       x_tmp = (/boozer_s,boozer_phi_beg,thetaB/)
       !PRINT *,ind_thtB, x_tmp
       !
       ! TESTING
       !x_tmp = (/boozer_s,boozer_phi_beg,boozer_theta_Vphi/)
       !
       CALL mag(x_tmp,bmod_tmp,sqrtg_tmp,bder_tmp,hcovar_tmp,&
            hctrvr_tmp,hcoder_tmp,hctder_tmp)
       CALL compute_Gsymm( x_tmp, G_symm, G_symm_tb, G_symm_pb )
       CALL compute_RZ( x_tmp, R, R_tb, Z, Z_tb )
       R_Vphi_prof(ind_thtB) = R
       Z_Vphi_prof(ind_thtB) = Z
       !
       ! compute additional B-field quantities
       sqrtg_bctrvr_phi=sqrtg_bctrvr_tht/aiota_loc
       bcovar_phi=hcovar_tmp(2)*(bmod_tmp*1.0e4_dp)
       bcovar_tht=hcovar_tmp(3)*(bmod_tmp*1.0e4_dp)
       boozer_psi_pr = boozer_psi_pr_hat*(bmod0*1.0e4_dp)
       avb2 = avbhat2*((bmod0*1.0e4_dp)**2)
       !
       coef_K_Er = -c*bcovar_phi/sqrtg_bctrvr_tht
       fac1 = (hctrvr_tmp(2)*bmod_tmp*1.0e4_dp) * &
            (1.0_dp+TWOPI*aiota_loc*boozer_psi_pr*G_symm_tb) / avb2
       denom_Er = (c / (avnabpsi*aiota_loc*boozer_psi_pr)) + &
            fac1 * (-c*bcovar_phi/sqrtg_bctrvr_tht)
       !PRINT *,'denom_Er - 1st term: ',denom_Er
       DO ispec_ctr = LBOUND(row_ind_ptr,1),UBOUND(row_ind_ptr,1)
          irow_spec = row_ind_ptr(ispec_ctr)
          icol_spec = col_ind_ptr(ispec_ctr)
          denom_Er_1 = &
               D31AX_spec(ispec_ctr) * (z_spec(icol_spec)*e) / T_spec(icol_spec)
          denom_Er(irow_spec) = denom_Er(irow_spec) + fac1 * denom_Er_1
          !PRINT *, 'denom_Er - 2nd term: ', denom_Er_1
          coef_K_Er(irow_spec) = coef_K_Er(irow_spec) + denom_Er_1
          !
       END DO
       !
       coef_K_grad = 0.0_dp
       nom_Er_1 = 0.0_dp
       nom_Er_2 = 0.0_dp
       DO irow_ctr = 0,num_spec-1
          nom_Er_1(irow_ctr) = ((c*T_spec(irow_ctr)/(z_spec(irow_ctr)*e)) / &
               (avnabpsi*aiota_loc*boozer_psi_pr)) * &
               (dn_spec_ov_ds(irow_ctr) / n_spec(irow_ctr) + &
               dT_spec_ov_ds(irow_ctr) / T_spec(irow_ctr)) * avnabpsi
          nom_Er_2(irow_ctr) = (c*T_spec(irow_ctr)/(z_spec(irow_ctr)*e)) * &
               (-fac1*bcovar_phi/sqrtg_bctrvr_tht) * &
               (dn_spec_ov_ds(irow_ctr) / n_spec(irow_ctr) + &
               dT_spec_ov_ds(irow_ctr) / T_spec(irow_ctr)) * avnabpsi
          coef_K_grad(irow_ctr) = nom_Er_2(irow_ctr) / (-fac1)
       END DO
       nom_Er = nom_Er_1 + nom_Er_2
       DO ispec_ctr = LBOUND(row_ind_ptr,1),UBOUND(row_ind_ptr,1)
          irow_spec = row_ind_ptr(ispec_ctr)
          icol_spec = col_ind_ptr(ispec_ctr)
          !
          nom_Er_3 = avnabpsi * D31AX_spec(ispec_ctr) * &
               (dn_spec_ov_ds(icol_spec) / n_spec(icol_spec) + &
               dT_spec_ov_ds(icol_spec) / T_spec(icol_spec))
          nom_Er_4 = avnabpsi * (dT_spec_ov_ds(icol_spec) / T_spec(icol_spec)) * &
               (D32AX_spec(ispec_ctr) - 2.5_dp * D31AX_spec(ispec_ctr))
          nom_Er_5 = D33AX_spec(ispec_ctr) * avEparB_ov_avb2 * &
               (z_spec(icol_spec)*e) / T_spec(icol_spec)
          nom_Er(irow_spec) = nom_Er(irow_spec) + fac1 * (nom_Er_3 + nom_Er_4 + nom_Er_5)
          !PRINT *,'Er - 2nd term: ', nom_Er_3 / denom_Er
          !PRINT *,'Er - 3rd term: ', nom_Er_4 / denom_Er
          coef_K_grad(irow_spec) = coef_K_grad(irow_spec) - (nom_Er_3 + nom_Er_4 + nom_Er_5)
          !
       END DO
       !
       ! compute species poloidal dependence of Vtor
       Vphi_prof_spec(:,ind_thtB) = denom_Er * Er - nom_Er
       !
       ! compute species poloidal (Boozer) rotation frequency
       Vtht_prof_spec(:,ind_thtB) = (hctrvr_tmp(3) * bmod_tmp*1.0e4_dp / avb2) * &
            (coef_K_Er * Er + coef_K_grad)
       !
    END DO ! end loop over theta-grid
    !
    IF (ALLOCATED(denom_Er)) DEALLOCATE(denom_Er)
    IF (ALLOCATED(nom_Er)) DEALLOCATE(nom_Er)
    IF (ALLOCATED(nom_Er_1)) DEALLOCATE(nom_Er_1)
    IF (ALLOCATED(nom_Er_2)) DEALLOCATE(nom_Er_2)
    IF (ALLOCATED(coef_K)) DEALLOCATE(coef_K)
    IF (ALLOCATED(coef_K_Er)) DEALLOCATE(coef_K_Er)
    IF (ALLOCATED(coef_K_grad)) DEALLOCATE(coef_K_grad)
    if (allocated(D33AX_spec)) deallocate(D33AX_spec)
  END SUBROUTINE compute_Vphi_profile
  !
  SUBROUTINE compute_VthtB_and_VphiB_a(row_ind_ptr, col_ind_ptr, &
       D31AX_spec, D32AX_spec, Er, VthtB_spec, VphiB_spec)
    !
    USE neo_precision
    USE neo_control, ONLY: lab_swi
    USE device_mod, ONLY : surface
    USE mag_interface_mod, ONLY : mag_coordinates, &
         boozer_s, boozer_theta_beg, boozer_phi_beg
    USE partpa_mod,  ONLY : bmod0
    USE mag_sub, ONLY: mag
    use neo_magfie, only : boozer_curr_pol_hat, boozer_psi_pr_hat
    USE collisionality_mod, ONLY : num_spec, species_tag, &
         z_spec, m_spec, n_spec, T_spec
    !
    ! ---------------------------------------------------------------!
    ! input:
    ! ---------------------------------------------------------------!
    ! row- and column indices (=species) of diffusion coefficients
    INTEGER, DIMENSION(:), INTENT(in) :: row_ind_ptr, col_ind_ptr
    ! species parallel flow
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: D31AX_spec, D32AX_spec
    ! radial electric field (w.r.t. effective radius)
    REAL(kind=dp), INTENT(in) :: Er
    ! ---------------------------------------------------------------!
    ! output:
    ! ---------------------------------------------------------------!
    ! species poloidal and toroidal rotation velocities:
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(inout) :: VthtB_spec
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(inout) :: VphiB_spec
    ! ---------------------------------------------------------------!
    ! local:
    ! ---------------------------------------------------------------!
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: y
    REAL(kind=dp) :: aiota_loc, avnabpsi
    ! ---------------------------------------------------------------!
    ! co- and contra-variant B-field components using
    ! r_eff as a flux-surface label
    REAL(kind=dp) :: bcovar_tht, bcovar_phi
    REAL(kind=dp) :: sqrtg_bctrvr_phi, boozer_psi_pr
    REAL(kind=dp) :: sqrtg_b2
    ! normalized co-variant phi-component of B,  $\sqrt{g}B^\vartheta$
    REAL(kind=dp) :: bcovar_phi_hat, sqrtg_bctrvr_tht
    ! Only used to call mag for normalizations
    ! related to D31 and D32
    REAL(kind=dp)                 :: bmod_tmp,sqrtg_tmp
    REAL(kind=dp), DIMENSION(3)   :: x_tmp,bder_tmp,hcovar_tmp,hctrvr_tmp
    REAL(kind=dp), DIMENSION(3,3) :: hcoder_tmp,hctder_tmp
    ! ---------------------------------------------------------------!
    ! loop indices, temporary variables
    INTEGER :: ispec_ctr, irow_spec, icol_spec
    REAL(kind=dp) :: pressure_a, dpressure_ov_dr_a
    REAL(kind=dp) :: Vtht_a_perp, Vphi_a_perp
    REAL(kind=dp) :: fac_Vtht, fac_Vphi
    REAL(kind=dp) :: Ta_ov_ea
    REAL(kind=dp) :: A1_b, A2_b, VparB_a
    ! ---------------------------------------------------------------!
    !
    ! copy y-vector (see definition in rhs_kin.f90)
    ALLOCATE(y(SIZE(y_ntv_mod,1)))
    y = y_ntv_mod
    !
    ! aiota, avnabpsi :
    aiota_loc = surface%aiota
    avnabpsi = y(7) / y(6)
    !
    ! computation of the normalization for D31 and D32 (-> D31_ref)
    IF (mag_coordinates .EQ. 0) THEN
       ! cylindrical coordinates
       x_tmp = xstart_cyl
       CALL mag(x_tmp,bmod_tmp,sqrtg_tmp,bder_tmp,hcovar_tmp,&
            hctrvr_tmp,hcoder_tmp,hctder_tmp)
       !
       ! normalized co-variant phi-component of B
       ! (Note! co-variant phi-component of B is the
       ! same for cylindrical coordinates and symmetry flux
       ! coodrinates --> no conversion needed)
       bcovar_phi_hat = hcovar_tmp(2)*(bmod_tmp/bmod0)
       !
       ! restore value of $\sqrt{g}B^\vartheta$ for
       ! symmetry flux coordinates from quantities
       ! given in cylindircal coordinates
       sqrtg_bctrvr_tht = avnabpsi*sqrtg_tmp*(hctrvr_tmp(3)*bmod_tmp*1.0e4_dp)
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
       IF (lab_swi .EQ. 10) THEN ! ASDEX-U (E. Strumberger)
          ! this is the same up to a minus sign resulting from the
          ! definition of sqrtg_tmp (left-handed system)
          sqrtg_bctrvr_tht = avnabpsi*sqrtg_tmp*(hctrvr_tmp(3)*bmod_tmp*1.0d4) 
       ELSE
          ! restore value of $\sqrt{g}B^\vartheta$ for
          ! symmetry flux coordinates from the quantities
          ! given in Boozer coordinates (only valid for right-handed system)
          sqrtg_bctrvr_tht = avnabpsi*aiota_loc*boozer_psi_pr_hat*(bmod0*1.0e4_dp)
       END IF
    END IF
    !
    ! compute additional B-field quantities
    sqrtg_bctrvr_phi = sqrtg_bctrvr_tht/aiota_loc
    bcovar_phi = hcovar_tmp(2)*(bmod_tmp*1.0e4_dp)
    bcovar_tht = hcovar_tmp(3)*(bmod_tmp*1.0e4_dp)
    boozer_psi_pr = boozer_psi_pr_hat*(bmod0*1.0e4_dp)
    sqrtg_b2 = avnabpsi * boozer_psi_pr * &
         (aiota_loc * bcovar_tht + bcovar_phi)
    !
    ! allocate species poloidal and toroidal rotation velocities
    IF(ALLOCATED(VthtB_spec)) DEALLOCATE(VthtB_spec)
    ALLOCATE(VthtB_spec(0:num_spec-1))
    VthtB_spec = 0.0_dp
    IF(ALLOCATED(VphiB_spec)) DEALLOCATE(VphiB_spec)
    ALLOCATE(VphiB_spec(0:num_spec-1))
    VphiB_spec = 0.0_dp
    !
    ! compute perpendicular components of VthtB_spec and VphiB_spec
    fac_Vtht = c * bcovar_phi / sqrtg_b2
    fac_Vphi = -c * bcovar_tht / sqrtg_b2
    DO irow_spec = 0,num_spec-1
       !
       pressure_a = n_spec(irow_spec) * T_spec(irow_spec)
       dpressure_ov_dr_a = avnabpsi * &
            (T_spec(irow_spec) * dn_spec_ov_ds(irow_spec) + &
            n_spec(irow_spec) * dT_spec_ov_ds(irow_spec))
       Ta_ov_ea = T_spec(irow_spec) / (z_spec(irow_spec) * e)
       !
       Vtht_a_perp = fac_Vtht * Ta_ov_ea * &
            (dpressure_ov_dr_a / pressure_a - Er / Ta_ov_ea)
       VthtB_spec(irow_spec) = Vtht_a_perp
       !
       Vphi_a_perp = fac_Vphi * Ta_ov_ea * &
            (dpressure_ov_dr_a / pressure_a - Er / Ta_ov_ea)
       VphiB_spec(irow_spec) = Vphi_a_perp
       !
    END DO
    !
    ! compute parallel components of VthtB_spec and VphiB_spec
    DO ispec_ctr = LBOUND(row_ind_ptr,1),UBOUND(row_ind_ptr,1)
       !
       irow_spec = row_ind_ptr(ispec_ctr)
       icol_spec = col_ind_ptr(ispec_ctr)
       !
       A2_b = avnabpsi * (dT_spec_ov_ds(icol_spec) / T_spec(icol_spec))
       A1_b = avnabpsi * (dn_spec_ov_ds(icol_spec) / n_spec(icol_spec)) - &
            1.5_dp * A2_b - Er * (z_spec(icol_spec)*e) / T_spec(icol_spec)
       !
       VparB_a = -(D31AX_spec(ispec_ctr) * A1_b + D32AX_spec(ispec_ctr) * A2_b)
       !
       VthtB_spec(irow_spec) = VthtB_spec(irow_spec) + &
            VparB_a * aiota_loc / (aiota_loc * bcovar_tht + bcovar_phi)
       !
       VphiB_spec(irow_spec) = VphiB_spec(irow_spec) + &
            VparB_a / (aiota_loc * bcovar_tht + bcovar_phi)
       !
    END DO
    !
  END SUBROUTINE compute_VthtB_and_VphiB_a
  !> \note Makes use of compute_VthtB_and_VphiB_a, to compute the
  !>   corresponding quantities.
  SUBROUTINE compute_VthtB_and_VphiB_b(row_ind_ptr, col_ind_ptr, &
       D31AX_spec, D32AX_spec, D33AX_spec, Er, avEparB_ov_avb2, &
       VthtB_spec, VphiB_spec, VthtB_Ware_spec, VphiB_Ware_spec)
    !
    USE neo_precision
    USE neo_control, ONLY: lab_swi
    USE device_mod, ONLY : surface
    USE mag_interface_mod, ONLY : mag_coordinates, &
         boozer_s, boozer_theta_beg, boozer_phi_beg
    USE partpa_mod,  ONLY : bmod0
    USE mag_sub, ONLY: mag
    use neo_magfie, only : boozer_curr_pol_hat, boozer_psi_pr_hat
    USE collisionality_mod, ONLY : num_spec, species_tag, &
         z_spec, m_spec, n_spec, T_spec
    !
    ! ---------------------------------------------------------------!
    ! input:
    ! ---------------------------------------------------------------!
    ! row- and column indices (=species) of diffusion coefficients
    INTEGER, DIMENSION(:), INTENT(in) :: row_ind_ptr, col_ind_ptr
    ! species parallel flow
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: D31AX_spec, D32AX_spec
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: D33AX_spec
    ! radial electric field (w.r.t. effective radius)
    ! + drive A3 normalized (=inductive electric field)
    REAL(kind=dp), INTENT(in) :: Er, avEparB_ov_avb2
    ! ---------------------------------------------------------------!
    ! output:
    ! ---------------------------------------------------------------!
    ! species poloidal and toroidal rotation velocities:
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(inout) :: VthtB_spec
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(inout) :: VphiB_spec
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(inout) :: VthtB_Ware_spec
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(inout) :: VphiB_Ware_spec
    ! ---------------------------------------------------------------!
    ! local:
    ! ---------------------------------------------------------------!
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: y
    REAL(kind=dp) :: aiota_loc, avnabpsi
    ! ---------------------------------------------------------------!
    ! co- and contra-variant B-field components using
    ! r_eff as a flux-surface label
    REAL(kind=dp) :: bcovar_tht, bcovar_phi
    REAL(kind=dp) :: sqrtg_bctrvr_phi, boozer_psi_pr
    REAL(kind=dp) :: sqrtg_b2
    ! normalized co-variant phi-component of B,  $\sqrt{g}B^\vartheta$
    REAL(kind=dp) :: bcovar_phi_hat, sqrtg_bctrvr_tht
    ! Only used to call mag for normalizations
    ! related to D31 and D32
    REAL(kind=dp)                 :: bmod_tmp,sqrtg_tmp
    REAL(kind=dp), DIMENSION(3)   :: x_tmp,bder_tmp,hcovar_tmp,hctrvr_tmp
    REAL(kind=dp), DIMENSION(3,3) :: hcoder_tmp,hctder_tmp
    ! ---------------------------------------------------------------!
    ! loop indices, temporary variables
    INTEGER :: ispec_ctr, irow_spec, icol_spec
    REAL(kind=dp) :: A3_b, VparB_a
    ! ---------------------------------------------------------------!
    !
    ! copy y-vector (see definition in rhs_kin.f90)
    ALLOCATE(y(SIZE(y_ntv_mod,1)))
    y = y_ntv_mod
    !
    ! aiota, avnabpsi :
    aiota_loc = surface%aiota
    avnabpsi = y(7) / y(6)
    !
    ! computation of the normalization for D31 and D32 (-> D31_ref)
    IF (mag_coordinates .EQ. 0) THEN
       ! cylindrical coordinates
       x_tmp = xstart_cyl
       CALL mag(x_tmp,bmod_tmp,sqrtg_tmp,bder_tmp,hcovar_tmp,&
            hctrvr_tmp,hcoder_tmp,hctder_tmp)
       !
       ! normalized co-variant phi-component of B
       ! (Note! co-variant phi-component of B is the
       ! same for cylindrical coordinates and symmetry flux
       ! coodrinates --> no conversion needed)
       bcovar_phi_hat = hcovar_tmp(2)*(bmod_tmp/bmod0)
       !
       ! restore value of $\sqrt{g}B^\vartheta$ for
       ! symmetry flux coordinates from quantities
       ! given in cylindircal coordinates
       sqrtg_bctrvr_tht = avnabpsi*sqrtg_tmp*(hctrvr_tmp(3)*bmod_tmp*1.0e4_dp)
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
       IF (lab_swi .EQ. 10) THEN ! ASDEX-U (E. Strumberger)
          ! this is the same up to a minus sign resulting from the
          ! definition of sqrtg_tmp (left-handed system)
          sqrtg_bctrvr_tht = avnabpsi*sqrtg_tmp*(hctrvr_tmp(3)*bmod_tmp*1.0d4) 
       ELSE
          ! restore value of $\sqrt{g}B^\vartheta$ for
          ! symmetry flux coordinates from the quantities
          ! given in Boozer coordinates (only valid for right-handed system)
          sqrtg_bctrvr_tht = avnabpsi*aiota_loc*boozer_psi_pr_hat*(bmod0*1.0e4_dp)
       END IF
    END IF
    !
    ! compute additional B-field quantities
    sqrtg_bctrvr_phi = sqrtg_bctrvr_tht/aiota_loc
    bcovar_phi = hcovar_tmp(2)*(bmod_tmp*1.0e4_dp)
    bcovar_tht = hcovar_tmp(3)*(bmod_tmp*1.0e4_dp)
    boozer_psi_pr = boozer_psi_pr_hat*(bmod0*1.0e4_dp)
    sqrtg_b2 = avnabpsi * boozer_psi_pr * &
         (aiota_loc * bcovar_tht + bcovar_phi)
    !
    ! compute species poloidal and toroidal rotation velocities
    ! without account of inductive electric field
    CALL compute_VthtB_and_VphiB_a(row_ind_ptr, col_ind_ptr, &
         D31AX_spec, D32AX_spec, Er, VthtB_spec, VphiB_spec)
    !
    ! allocate species poloidal and toroidal rotation velocities (Ware pinch contribution)
    IF(ALLOCATED(VthtB_Ware_spec)) DEALLOCATE(VthtB_Ware_spec)
    ALLOCATE(VthtB_Ware_spec(0:num_spec-1))
    VthtB_Ware_spec = 0.0_dp
    IF(ALLOCATED(VphiB_Ware_spec)) DEALLOCATE(VphiB_Ware_spec)
    ALLOCATE(VphiB_Ware_spec(0:num_spec-1))
    VphiB_Ware_spec = 0.0_dp
    !
    ! add contribution from inductive electric field to
    ! VthtB_spec and VphiB_spec    
    DO ispec_ctr = LBOUND(row_ind_ptr,1),UBOUND(row_ind_ptr,1)
       !
       irow_spec = row_ind_ptr(ispec_ctr)
       icol_spec = col_ind_ptr(ispec_ctr)
       !
       A3_b = avEparB_ov_avb2 * (z_spec(icol_spec)*e) / T_spec(icol_spec)
       !
       VparB_a = -D33AX_spec(ispec_ctr) * A3_b
       !
       VthtB_Ware_spec(irow_spec) = VthtB_Ware_spec(irow_spec) + &
            VparB_a * aiota_loc / (aiota_loc * bcovar_tht + bcovar_phi)
       VthtB_spec(irow_spec) = VthtB_spec(irow_spec) + &
            VparB_a * aiota_loc / (aiota_loc * bcovar_tht + bcovar_phi)
       !
       VphiB_Ware_spec(irow_spec) = VphiB_Ware_spec(irow_spec) + &
            VparB_a / (aiota_loc * bcovar_tht + bcovar_phi)
       VphiB_spec(irow_spec) = VphiB_spec(irow_spec) + &
            VparB_a / (aiota_loc * bcovar_tht + bcovar_phi)
       !
    END DO
    !
  END SUBROUTINE compute_VthtB_and_VphiB_b
  !
  SUBROUTINE compute_Gamma_a(row_ind_ptr, col_ind_ptr, &
       D11_spec, D12_spec, Er, Gamma_spec)
    !
    USE neo_precision
    USE collisionality_mod, ONLY : num_spec, species_tag, &
         z_spec, m_spec, n_spec, T_spec
    !
    ! ---------------------------------------------------------------!
    ! input:
    ! ---------------------------------------------------------------!
    ! row- and column indices (=species) of diffusion coefficients
    INTEGER, DIMENSION(:), INTENT(in) :: row_ind_ptr, col_ind_ptr
    ! species diffusion coefficients
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: D11_spec, D12_spec
    ! radial electric field (w.r.t. effective radius)
    REAL(kind=dp), INTENT(in) :: Er
    ! ---------------------------------------------------------------!
    ! output:
    ! ---------------------------------------------------------------!
    ! species particle flux density:
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(inout) :: Gamma_spec
    ! ---------------------------------------------------------------!
    ! local:
    ! ---------------------------------------------------------------!
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: y
    REAL(kind=dp) :: avnabpsi
    ! ---------------------------------------------------------------!
    ! loop indices, temporary variables
    INTEGER :: ispec_ctr, irow_spec, icol_spec
    REAL(kind=dp) :: A1_b, A2_b, flux_a
    ! ---------------------------------------------------------------!
    !
    ! copy y-vector (see definition in rhs_kin.f90)
    ALLOCATE(y(SIZE(y_ntv_mod,1)))
    y = y_ntv_mod
    !
    ! avnabpsi :
    avnabpsi = y(7) / y(6)
    !
    ! allocate species particle flux density
    IF(ALLOCATED(Gamma_spec)) DEALLOCATE(Gamma_spec)
    ALLOCATE(Gamma_spec(0:num_spec-1))
    Gamma_spec = 0.0_dp
    !
    ! compute species particle flux density
    DO ispec_ctr = LBOUND(row_ind_ptr,1),UBOUND(row_ind_ptr,1)
       !
       irow_spec = row_ind_ptr(ispec_ctr)
       icol_spec = col_ind_ptr(ispec_ctr)
       !
       A2_b = avnabpsi * (dT_spec_ov_ds(icol_spec) / T_spec(icol_spec))
       A1_b = avnabpsi * (dn_spec_ov_ds(icol_spec) / n_spec(icol_spec)) - &
            1.5_dp * A2_b - Er * (z_spec(icol_spec)*e) / T_spec(icol_spec)
       !
       flux_a = -n_spec(irow_spec) * &
            (D11_spec(ispec_ctr) * A1_b + D12_spec(ispec_ctr) * A2_b)
       !
       Gamma_spec(irow_spec) = Gamma_spec(irow_spec) + flux_a
       !
    END DO
    !
  END SUBROUTINE compute_Gamma_a
  !> \note Makes use of compute_Gamma_a, to compute the corresponding
  !>   quantities.
  SUBROUTINE compute_Gamma_b(row_ind_ptr, col_ind_ptr, &
       D11_spec, D12_spec, D13_spec, Er, avEparB_ov_avb2, &
       Gamma_spec, Gamma_Ware_spec)
    !
    USE neo_precision
    USE collisionality_mod, ONLY : num_spec, species_tag, &
         z_spec, m_spec, n_spec, T_spec
    !
    ! ---------------------------------------------------------------!
    ! input:
    ! ---------------------------------------------------------------!
    ! row- and column indices (=species) of diffusion coefficients
    INTEGER, DIMENSION(:), INTENT(in) :: row_ind_ptr, col_ind_ptr
    ! species diffusion coefficients
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: D11_spec, D12_spec
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: D13_spec
    ! radial electric field (w.r.t. effective radius)
    ! + drive A3 normalized (=inductive electric field)
    REAL(kind=dp), INTENT(in) :: Er, avEparB_ov_avb2
    ! ---------------------------------------------------------------!
    ! output:
    ! ---------------------------------------------------------------!
    ! species particle flux density:
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(inout) :: Gamma_spec
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(inout) :: Gamma_Ware_spec
    ! ---------------------------------------------------------------!
    ! local:
    ! ---------------------------------------------------------------!
    ! loop indices, temporary variables
    INTEGER :: ispec_ctr, irow_spec, icol_spec
    REAL(kind=dp) :: A3_b, flux_a
    ! ---------------------------------------------------------------!
    !
    ! allocate species particle flux density (Ware pinch contribution)
    IF(ALLOCATED(Gamma_Ware_spec)) DEALLOCATE(Gamma_Ware_spec)
    ALLOCATE(Gamma_Ware_spec(0:num_spec-1))
    Gamma_Ware_spec = 0.0_dp
    !
    ! compute species particle flux density
    ! without account of inductive electric field
    CALL compute_Gamma_a(row_ind_ptr, col_ind_ptr, &
         D11_spec, D12_spec, Er, Gamma_spec)
    !
    ! add contribution from inductive electric field to Gamma_spec
    DO ispec_ctr = LBOUND(row_ind_ptr,1),UBOUND(row_ind_ptr,1)
       !
       irow_spec = row_ind_ptr(ispec_ctr)
       icol_spec = col_ind_ptr(ispec_ctr)
       !
       A3_b = avEparB_ov_avb2 * (z_spec(icol_spec)*e) / T_spec(icol_spec)
       !
       flux_a = -n_spec(irow_spec) * D13_spec(ispec_ctr) * A3_b
       !
       Gamma_Ware_spec(irow_spec) = Gamma_Ware_spec(irow_spec) + flux_a
       Gamma_spec(irow_spec) = Gamma_spec(irow_spec) + flux_a
       !
    END DO
    !
  END SUBROUTINE compute_Gamma_b
  !
  SUBROUTINE compute_Qflux_a(row_ind_ptr, col_ind_ptr, &
       D21_spec, D22_spec, Er, Qflux_spec)
    !
    USE neo_precision
    USE collisionality_mod, ONLY : num_spec, species_tag, &
         z_spec, m_spec, n_spec, T_spec
    !
    ! ---------------------------------------------------------------!
    ! input:
    ! ---------------------------------------------------------------!
    ! row- and column indices (=species) of diffusion coefficients
    INTEGER, DIMENSION(:), INTENT(in) :: row_ind_ptr, col_ind_ptr
    ! species diffusion coefficients
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: D21_spec, D22_spec
    ! radial electric field (w.r.t. effective radius)
    REAL(kind=dp), INTENT(in) :: Er
    ! ---------------------------------------------------------------!
    ! output:
    ! ---------------------------------------------------------------!
    ! species heat flux density:
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(inout) :: Qflux_spec
    ! ---------------------------------------------------------------!
    ! local:
    ! ---------------------------------------------------------------!
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: y
    REAL(kind=dp) :: avnabpsi
    ! ---------------------------------------------------------------!
    ! loop indices, temporary variables
    INTEGER :: ispec_ctr, irow_spec, icol_spec
    REAL(kind=dp) :: A1_b, A2_b, flux_a
    ! ---------------------------------------------------------------!
    !
    ! copy y-vector (see definition in rhs_kin.f90)
    ALLOCATE(y(SIZE(y_ntv_mod,1)))
    y = y_ntv_mod
    !
    ! avnabpsi :
    avnabpsi = y(7) / y(6)
    !
    ! allocate species heat flux density
    IF(ALLOCATED(Qflux_spec)) DEALLOCATE(Qflux_spec)
    ALLOCATE(Qflux_spec(0:num_spec-1))
    Qflux_spec = 0.0_dp
    !
    ! compute species heat flux density
    DO ispec_ctr = LBOUND(row_ind_ptr,1),UBOUND(row_ind_ptr,1)
       !
       irow_spec = row_ind_ptr(ispec_ctr)
       icol_spec = col_ind_ptr(ispec_ctr)
       !
       A2_b = avnabpsi * (dT_spec_ov_ds(icol_spec) / T_spec(icol_spec))
       A1_b = avnabpsi * (dn_spec_ov_ds(icol_spec) / n_spec(icol_spec)) - &
            1.5_dp * A2_b - Er * (z_spec(icol_spec)*e) / T_spec(icol_spec)
       !
       flux_a = -n_spec(irow_spec) * &
            (D21_spec(ispec_ctr) * A1_b + D22_spec(ispec_ctr) * A2_b)
       !
       Qflux_spec(irow_spec) = Qflux_spec(irow_spec) + flux_a
       !
    END DO
    !
    Qflux_spec = Qflux_spec * T_spec
    !
  END SUBROUTINE compute_Qflux_a
  !> \note Makes use of compute_Qflux_a, to compute the corresponding
  !>   quantities.
  SUBROUTINE compute_Qflux_b(row_ind_ptr, col_ind_ptr, &
       D21_spec, D22_spec, D23_spec, Er, avEparB_ov_avb2, &
       Qflux_spec, Qflux_Ware_spec)
    !
    USE neo_precision
    USE collisionality_mod, ONLY : num_spec, species_tag, &
         z_spec, m_spec, n_spec, T_spec
    !
    ! ---------------------------------------------------------------!
    ! input:
    ! ---------------------------------------------------------------!
    ! row- and column indices (=species) of diffusion coefficients
    INTEGER, DIMENSION(:), INTENT(in) :: row_ind_ptr, col_ind_ptr
    ! species diffusion coefficients
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: D21_spec, D22_spec
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: D23_spec
    ! radial electric field (w.r.t. effective radius)
    ! + drive A3 normalized (=inductive electric field)
    REAL(kind=dp), INTENT(in) :: Er, avEparB_ov_avb2
    ! ---------------------------------------------------------------!
    ! output:
    ! ---------------------------------------------------------------!
    ! species heat flux density:
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(inout) :: Qflux_spec
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(inout) :: Qflux_Ware_spec
    ! ---------------------------------------------------------------!
    ! local:
    ! ---------------------------------------------------------------!
    ! loop indices, temporary variables
    INTEGER :: ispec_ctr, irow_spec, icol_spec
    REAL(kind=dp) :: A3_b, flux_a
    ! ---------------------------------------------------------------!
    !
    ! allocate species heat flux density (Ware pinch contribution)
    IF(ALLOCATED(Qflux_Ware_spec)) DEALLOCATE(Qflux_Ware_spec)
    ALLOCATE(Qflux_Ware_spec(0:num_spec-1))
    Qflux_Ware_spec = 0.0_dp
    !
    ! compute species heat flux density
    ! without account of inductive electric field
    CALL compute_Qflux_a(row_ind_ptr, col_ind_ptr, &
         D21_spec, D22_spec, Er, Qflux_spec)
    !
    ! add contribution from inductive electric field to Qflux_spec
    DO ispec_ctr = LBOUND(row_ind_ptr,1),UBOUND(row_ind_ptr,1)
       !
       irow_spec = row_ind_ptr(ispec_ctr)
       icol_spec = col_ind_ptr(ispec_ctr)
       !
       A3_b = avEparB_ov_avb2 * (z_spec(icol_spec)*e) / T_spec(icol_spec)
       !
       flux_a = -n_spec(irow_spec) * D23_spec(ispec_ctr) * A3_b
       !
       Qflux_Ware_spec(irow_spec) = Qflux_Ware_spec(irow_spec) + flux_a
       !
    END DO
    !
    Qflux_Ware_spec = Qflux_Ware_spec * T_spec
    Qflux_spec = Qflux_spec + Qflux_Ware_spec
    !
  END SUBROUTINE compute_Qflux_b
  !
  SUBROUTINE compute_ParFlow_a(row_ind_ptr, col_ind_ptr, &
       D31_spec, D32_spec, Er, ParFlow_spec)
    !
    USE neo_precision
    USE collisionality_mod, ONLY : num_spec, species_tag, &
         z_spec, m_spec, n_spec, T_spec
    !
    ! ---------------------------------------------------------------!
    ! input:
    ! ---------------------------------------------------------------!
    ! row- and column indices (=species) of diffusion coefficients
    INTEGER, DIMENSION(:), INTENT(in) :: row_ind_ptr, col_ind_ptr
    ! species diffusion coefficients
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: D31_spec, D32_spec
    ! radial electric field (w.r.t. effective radius)
    REAL(kind=dp), INTENT(in) :: Er
    ! ---------------------------------------------------------------!
    ! output:
    ! ---------------------------------------------------------------!
    ! species parallel flow:
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(inout) :: ParFlow_spec
    ! ---------------------------------------------------------------!
    ! local:
    ! ---------------------------------------------------------------!
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: y
    REAL(kind=dp) :: avnabpsi
    ! ---------------------------------------------------------------!
    ! loop indices, temporary variables
    INTEGER :: ispec_ctr, irow_spec, icol_spec
    REAL(kind=dp) :: A1_b, A2_b, flux_a
    ! ---------------------------------------------------------------!
    !
    ! copy y-vector (see definition in rhs_kin.f90)
    ALLOCATE(y(SIZE(y_ntv_mod,1)))
    y = y_ntv_mod
    !
    ! avnabpsi :
    avnabpsi = y(7) / y(6)
    !
    ! allocate species parallel flow
    IF(ALLOCATED(ParFlow_spec)) DEALLOCATE(ParFlow_spec)
    ALLOCATE(ParFlow_spec(0:num_spec-1))
    ParFlow_spec = 0.0_dp
    !
    ! compute species parallel flow
    DO ispec_ctr = LBOUND(row_ind_ptr,1),UBOUND(row_ind_ptr,1)
       !
       irow_spec = row_ind_ptr(ispec_ctr)
       icol_spec = col_ind_ptr(ispec_ctr)
       !
       A2_b = avnabpsi * (dT_spec_ov_ds(icol_spec) / T_spec(icol_spec))
       A1_b = avnabpsi * (dn_spec_ov_ds(icol_spec) / n_spec(icol_spec)) - &
            1.5_dp * A2_b - Er * (z_spec(icol_spec)*e) / T_spec(icol_spec)
       !
       flux_a = -n_spec(irow_spec) * &
            (D31_spec(ispec_ctr) * A1_b + D32_spec(ispec_ctr) * A2_b)
       !
       ParFlow_spec(irow_spec) = ParFlow_spec(irow_spec) + flux_a
       !
    END DO
    !
    ParFlow_spec = ParFlow_spec / n_spec
    !
  END SUBROUTINE compute_ParFlow_a
  !> \note Makes use of compute_ParFlow_a, to compute the corresponding
  !>   quantities.
  SUBROUTINE compute_ParFlow_b(row_ind_ptr, col_ind_ptr, &
       D31_spec, D32_spec, D33_spec, Er, avEparB_ov_avb2, &
       ParFlow_spec, ParFlow_Ware_spec)
    !
    USE neo_precision
    USE collisionality_mod, ONLY : num_spec, species_tag, &
         z_spec, m_spec, n_spec, T_spec
    !
    ! ---------------------------------------------------------------!
    ! input:
    ! ---------------------------------------------------------------!
    ! row- and column indices (=species) of diffusion coefficients
    INTEGER, DIMENSION(:), INTENT(in) :: row_ind_ptr, col_ind_ptr
    ! species diffusion coefficients
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: D31_spec, D32_spec
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: D33_spec
    ! radial electric field (w.r.t. effective radius)
    ! + drive A3 normalized (=inductive electric field)
    REAL(kind=dp), INTENT(in) :: Er, avEparB_ov_avb2
    ! ---------------------------------------------------------------!
    ! output:
    ! ---------------------------------------------------------------!
    ! species parallel flow:
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(inout) :: ParFlow_spec
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(inout) :: ParFlow_Ware_spec
    ! ---------------------------------------------------------------!
    ! local:
    ! ---------------------------------------------------------------!
    ! loop indices, temporary variables
    INTEGER :: ispec_ctr, irow_spec, icol_spec
    REAL(kind=dp) :: A3_b, flux_a
    ! ---------------------------------------------------------------!
    !
    ! allocate species parallel flow (Ware pinch contribution)
    IF(ALLOCATED(ParFlow_Ware_spec)) DEALLOCATE(ParFlow_Ware_spec)
    ALLOCATE(ParFlow_Ware_spec(0:num_spec-1))
    ParFlow_Ware_spec = 0.0_dp
    !
    ! compute species parallel flow
    ! without account of inductive electric field
    CALL compute_ParFlow_a(row_ind_ptr, col_ind_ptr, &
         D31_spec, D32_spec, Er, ParFlow_spec)
    !
    ! add contribution from inductive electric field to ParFlow_spec
    DO ispec_ctr = LBOUND(row_ind_ptr,1),UBOUND(row_ind_ptr,1)
       !
       irow_spec = row_ind_ptr(ispec_ctr)
       icol_spec = col_ind_ptr(ispec_ctr)
       !
       A3_b = avEparB_ov_avb2 * (z_spec(icol_spec)*e) / T_spec(icol_spec)
       !
       flux_a = -n_spec(irow_spec) * D33_spec(ispec_ctr) * A3_b
       !
       ParFlow_Ware_spec(irow_spec) = ParFlow_Ware_spec(irow_spec) + flux_a
       !
    END DO
    !
    ParFlow_Ware_spec = ParFlow_Ware_spec / n_spec
    ParFlow_spec = ParFlow_spec + ParFlow_Ware_spec
    !
  END SUBROUTINE compute_ParFlow_b
  !
  SUBROUTINE compute_TphiNA_a(Gamma_NA_spec, TphiNA_spec, TphiNA_tot)
    !
    USE neo_precision
    USE neo_control, ONLY: lab_swi
    USE device_mod, ONLY : surface
    USE mag_interface_mod, ONLY : mag_coordinates, &
         boozer_s, boozer_theta_beg, boozer_phi_beg
    USE partpa_mod,  ONLY : bmod0
    USE mag_sub, ONLY: mag
    use neo_magfie, only : boozer_curr_pol_hat, boozer_psi_pr_hat
    USE collisionality_mod, ONLY : num_spec, z_spec
    !
    ! ---------------------------------------------------------------!
    ! input:
    ! ---------------------------------------------------------------!
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: Gamma_NA_spec
    ! ---------------------------------------------------------------!
    ! output:
    ! ---------------------------------------------------------------!
    ! species NTV torque density and total NTV torque density:
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(inout) :: TphiNA_spec
    REAL(kind=dp) :: TphiNA_tot
    ! ---------------------------------------------------------------!
    ! local:
    ! ---------------------------------------------------------------!
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: y
    REAL(kind=dp) :: aiota_loc, avnabpsi
    ! ---------------------------------------------------------------!
    ! normalized co-variant phi-component of B,  $\sqrt{g}B^\vartheta$
    REAL(kind=dp) :: bcovar_phi_hat, sqrtg_bctrvr_tht
    ! Only used to call mag for normalizations
    ! related to D31 and D32
    REAL(kind=dp)                 :: bmod_tmp,sqrtg_tmp
    REAL(kind=dp), DIMENSION(3)   :: x_tmp,bder_tmp,hcovar_tmp,hctrvr_tmp
    REAL(kind=dp), DIMENSION(3,3) :: hcoder_tmp,hctder_tmp
    ! ---------------------------------------------------------------!
    ! loop indices, temporary variables
    INTEGER :: irow_spec
    REAL(kind=dp) :: fac_flux_force
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: Gamma_NA_spec_tmp
    ! ---------------------------------------------------------------!
    !
    ! allocate temporary storage array
    IF(ALLOCATED(Gamma_NA_spec_tmp)) DEALLOCATE(Gamma_NA_spec_tmp)
    ALLOCATE(Gamma_NA_spec_tmp(0:num_spec-1))
    Gamma_NA_spec_tmp = Gamma_NA_spec
    !
    ! copy y-vector (see definition in rhs_kin.f90)
    ALLOCATE(y(SIZE(y_ntv_mod,1)))
    y = y_ntv_mod
    !
    ! aiota, avnabpsi :
    aiota_loc = surface%aiota
    avnabpsi = y(7) / y(6)
    !
    ! computation of the normalization for D31 and D32 (-> D31_ref)
    IF (mag_coordinates .EQ. 0) THEN
       ! cylindrical coordinates
       x_tmp = xstart_cyl
       CALL mag(x_tmp,bmod_tmp,sqrtg_tmp,bder_tmp,hcovar_tmp,&
            hctrvr_tmp,hcoder_tmp,hctder_tmp)
       !
       ! normalized co-variant phi-component of B
       ! (Note! co-variant phi-component of B is the
       ! same for cylindrical coordinates and symmetry flux
       ! coodrinates --> no conversion needed)
       bcovar_phi_hat = hcovar_tmp(2)*(bmod_tmp/bmod0)
       !
       ! restore value of $\sqrt{g}B^\vartheta$ for
       ! symmetry flux coordinates from quantities
       ! given in cylindircal coordinates
       sqrtg_bctrvr_tht = avnabpsi*sqrtg_tmp*(hctrvr_tmp(3)*bmod_tmp*1.0e4_dp)
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
       IF (lab_swi .EQ. 10) THEN ! ASDEX-U (E. Strumberger)
          ! this is the same up to a minus sign resulting from the
          ! definition of sqrtg_tmp (left-handed system)
          sqrtg_bctrvr_tht = avnabpsi*sqrtg_tmp*(hctrvr_tmp(3)*bmod_tmp*1.0d4) 
       ELSE
          ! restore value of $\sqrt{g}B^\vartheta$ for
          ! symmetry flux coordinates from the quantities
          ! given in Boozer coordinates (only valid for right-handed system)
          sqrtg_bctrvr_tht = avnabpsi*aiota_loc*boozer_psi_pr_hat*(bmod0*1.0e4_dp)
       END IF
    END IF
    !
    ! allocate species NTV torque density
    IF(ALLOCATED(TphiNA_spec)) DEALLOCATE(TphiNA_spec)
    ALLOCATE(TphiNA_spec(0:num_spec-1))
    TphiNA_spec = 0.0_dp
    !
    ! compute species NTV torque density and total NTV torque
    fac_flux_force = -sqrtg_bctrvr_tht / c
    TphiNA_tot = 0.0_dp
    DO irow_spec = 0,num_spec-1
       TphiNA_spec(irow_spec) = &
            fac_flux_force * (z_spec(irow_spec)*e) * Gamma_NA_spec_tmp(irow_spec)
       TphiNA_tot = TphiNA_tot + TphiNA_spec(irow_spec)
    END DO
    !
  END SUBROUTINE compute_TphiNA_a
  !
END MODULE ntv_mod
