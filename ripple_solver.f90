!
  MODULE ntv_eqmat_mod
    INTEGER                                       :: nz_symm,nz_asymm,nz_regper
    INTEGER                                       :: nz_per_pos,nz_per_neg
    INTEGER,          DIMENSION(:),   ALLOCATABLE :: irow_symm,icol_symm
    INTEGER,          DIMENSION(:),   ALLOCATABLE :: irow_regper,icol_regper
    INTEGER,          DIMENSION(:),   ALLOCATABLE :: irow_asymm,icol_asymm
    INTEGER,          DIMENSION(:),   ALLOCATABLE :: irow_per_pos,icol_per_pos
    INTEGER,          DIMENSION(:),   ALLOCATABLE :: irow_per_neg,icol_per_neg
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: amat_symm
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: amat_regper
    DOUBLE COMPLEX,   DIMENSION(:),   ALLOCATABLE :: amat_asymm
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: f0_coll,f0_ttmp
  END MODULE ntv_eqmat_mod
!
!Sergei 20.07.2006 : modification of boundary layer is done now locally,
!                    filtering of the magnetic field maxima has been removed,
!                    old unused lines which were commented have been removed

SUBROUTINE ripple_solver(                                 &
     npass_l,npass_r,nvelocity,                            &
     amat_plus_plus,amat_minus_minus,                     &
     amat_plus_minus,amat_minus_plus,                     &
     source_p,source_m,                                   &
     flux_p,flux_m,                                       &
     qflux,                                               &
     ierr                                                 &
     )

  USE device_mod
  USE flint_mod, ONLY : phi_divide                      !<-in Winny
  USE size_mod, ONLY : ndim0
  USE collisionality_mod, ONLY : collpar,conl_over_mfp,isw_lorentz, &
       isw_energy,isw_integral,isw_axisymm, & !<-in Winny
       isw_momentum,vel_distri_swi,vel_num,vel_max, &
       nvel,vel_array
! collpar - the same as $\kappa$ - inverse mean-free path times 4
  !! Modifications by Andreas F. Martitsch (01.04.2015)
  !
  ! DEFINITION OF INPUT QUANTITIES:
  !
  ! a) collpar (also called $\kappa$ in the write-up): collpar = $\frac{4}{l_c}$,
  ! where the mean free path $l_c = 2 v_{Ta} \tau_{aa}$. (Note: The definition of
  ! $l_c$ differs from the one used in ntv_booz (Version Nov 2013) by a factor two, 
  ! which entered the part related to NTV, e.g., $\omega_{mm'}$. This scaling is
  ! fixed now according to the benchmarks published in PoP 21, 092506 (2014).)
  ! The quantites $v_{Ta}$ and $\tau_{aa}$ are given in cgs-units:
  ! -> $v_{Ta} = \sqrt{2T/m}$
  ! -> $\tau_{aa} = \frac{3 m_a^2 v_{Ta}^3}{16 \sqrt{\pi} n_a e_a^4 \Lambda_a}$
  !
  ! b) Mach number $M_t$: The electric rotation frequency ($\Omega_{tE}$) is
  ! computed from the toroidal Mach number, $\Omega_{tE} = M_t v_{Ta} / R$.
  ! Here $R$ denotes the major radius. The normalized toroidal rotation frequency
  ! becomes then: hatOmegaE = $\bar{\Omega}_{tE}$ = $\frac{M_t}{\kappa R}$
  !
  !! End Modifications by Andreas F. Martitsch (01.04.2015)
  USE lapack_band
  USE rkstep_mod
  USE polleg_mod
  ! USE binarysplit_mod, ONLY : bsfunc_reconstruct_levels !<- LOCAL
  ! new switch in propagator_mod
  USE propagator_mod,ONLY : prop_ripple_plot,prop_reconstruct,flux_mr,flux_pl, &
                            eta_modboundary_l,eta_modboundary_r,               &
                            sw_first_prop,sw_last_prop,                        &
                            prop_reconstruct_levels
  USE sparse_mod, ONLY : sparse_talk,sparse_solve_method,sparse_solve, &
       column_full2pointer,remap_rc,sparse_solver_test
  USE mag_interface_mod, ONLY: average_bhat,average_one_over_bhat,             &
                               surface_boozer_B00,travis_convfac,              &
                               mag_magfield
  USE ntv_eqmat_mod, ONLY : nz_symm,nz_asymm,nz_per_pos,nz_per_neg,            &
                            irow_symm,icol_symm,amat_symm,                     &
                            irow_per_pos,icol_per_pos,                         &
                            irow_per_neg,icol_per_neg,                         &
                            irow_asymm,icol_asymm,amat_asymm,                  &
                            f0_coll,f0_ttmp,                                   &
                            nz_regper,irow_regper,icol_regper,amat_regper
  USE partpa_mod, ONLY : bmod0
  USE development
  !! Modifications by Andreas F. Martitsch (12.03.2014)
  ! boozer_curr_tor, boozer_curr_pol, boozer_psi_pr,
  ! boozer_sqrtg11 and boozer_isqrg are now converted
  ! to cgs-units within neo_magfie.
  ! This step requires changes within rhs_kin.f90 and
  ! ripple_solver.f90!
  USE neo_magfie_mod, ONLY: boozer_iota,boozer_curr_pol_hat,&
       boozer_curr_tor_hat,boozer_psi_pr_hat,boozer_curr_pol_hat_s,&
       boozer_curr_tor_hat_s
  !! End Modifications by Andreas F. Martitsch (12.03.2014)
  !! Modifications by Andreas F. Martitsch (12.06.2014)
  ! quantities of the perturbation field extracted
  ! from the Boozer file (e.g., toroidal mode number m_phi)
  USE neo_magfie_perturbation, ONLY: calc_bnoverb0_arr, calc_ntv_output
  !! End Modifications by Andreas F. Martitsch (12.06.2014)
  !! Modification by Andreas F. Martitsch (14.07.2015)
  ! Extra input for NTV computations
  USE ntv_mod, ONLY : isw_qflux_NA, MtOvR, B_rho_L_loc, &
       m_phi,  qflux_symm, eps_M_2_val, av_gphph_val, av_inv_bhat_val
  USE neo_precision, ONLY : PI
  !! End Modification by Andreas F. Martitsch (14.07.2015)
  
  IMPLICIT NONE
  INTEGER, PARAMETER :: dp = KIND(1.0d0)
  DOUBLE COMPLEX, PARAMETER :: imun=(0.d0,1.d0)

  ! parameter list
  INTEGER,                                    INTENT(out)   :: npass_l
  INTEGER,                                    INTENT(out)   :: npass_r
  INTEGER,                                    INTENT(out)   :: nvelocity
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: amat_plus_plus
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: amat_minus_minus
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: amat_plus_minus
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: amat_minus_plus
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: source_p
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: source_m
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: flux_p
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: flux_m
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(out)   :: qflux
  INTEGER,                                    INTENT(out)   :: ierr
  
  ! local stuff
  ! Winny's output
  ! INTEGER, PARAMETER :: solver_talk = 0 ! (0: silent, 1: output)

  INTEGER :: add_global_eta
  INTEGER :: ub_eta,ub_mag,ub_eta_loc
  INTEGER :: ub_eta_prev,ub_eta_next,npass_l_out,npass_r_out
  INTEGER :: npart
  INTEGER :: ibeg,iend
  INTEGER :: bin_split_mode
  INTEGER :: i, i_loc
  INTEGER, ALLOCATABLE :: eta_loc_ind(:), eta_glob_ind(:)
  REAL(kind=dp) :: phibeg,phiend
  REAL(kind=dp) :: xetami,xetama
  REAL(kind=dp) :: rt0
  REAL(kind=dp) :: b_prop_l,b_prop_r,b_prop_min
  REAL(kind=dp) :: b_max_l,b_max_r,b_min,width
  REAL(kind=dp) :: b_l,b_r,eta_l,eta_r
  REAL(kind=dp), ALLOCATABLE :: eta(:),eta_prev(:),eta_next(:)
  REAL(kind=dp), ALLOCATABLE :: eta_loc(:),eta_glob(:) !<- LOCAL

  ! Timing
  REAL(kind=dp) :: time_start,time_factorization,time_solver


  !------------------------------------------------------------------------
  ! SERGEI
  !------------------------------------------------------------------------
  ! additional definitions from Sergei
  ! you can also use "double precision" isntead of  "REAL(kind=dp)"
!
  LOGICAL :: precond
!
  INTEGER :: itotstep,ntotstep,npart_loc,info
  INTEGER :: ndim,nhalf,istep,nstep,nback,npoints,npassing,ioddeven
  INTEGER :: ipart,nhalfstep,idir,inhom,ipart1
  INTEGER :: istp1,istp2,ifullstep,iendperiod,i1,ndim2
  INTEGER :: npass_l_min,npass_r_min,ibeg_ren,iend_ren,iempty,ndummy
  INTEGER :: nstep_renorm,nrenorm,irenorm,k,i1min,i1max,imat_min,imat_max
  INTEGER :: npass_old,npass_new
  INTEGER :: kmax,kshift,npass_bound,ibottom
  INTEGER, DIMENSION(:), ALLOCATABLE :: ipivot,iminvec,imaxvec
  INTEGER, DIMENSION(:), ALLOCATABLE :: iend_renorm
!
  DOUBLE PRECISION                         :: hf,theta_beg,theta_end,aiota
  DOUBLE PRECISION                         :: phi,hneg,eta0
  DOUBLE PRECISION                         :: dlu,subsq,subsqmin
  DOUBLE PRECISION                         :: diflam,diflampow,coefdir
  DOUBLE PRECISION                         :: coefenu,coefenu_averb   !!!term[1]
  DOUBLE PRECISION :: alambd_save1,alambd_save2,alambd_save3
  DOUBLE PRECISION :: diflam_flux
  DOUBLE PRECISION :: hxeta,fun_bound_new,amin2ovb
  DOUBLE PRECISION :: a11,a12,a21,a22,determ,deleta_b,dum_phi,dum_a3
  DOUBLE PRECISION :: exl,exprenorm,accurfac,renormfac
  DOUBLE COMPLEX :: coef_cf
!
  DOUBLE PRECISION, DIMENSION(6)           :: alp,bet,gam,del
  DOUBLE PRECISION, DIMENSION(4)           :: old_fun
  DOUBLE PRECISION, DIMENSION(ndim0)       :: ybeg
  DOUBLE PRECISION, DIMENSION(4,4)         :: a_beg,b_beg,a_end,b_end
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: amat,bvec_lapack,deriv_coef
  DOUBLE COMPLEX,   DIMENSION(:,:), ALLOCATABLE :: amat_z,bvec_lapack_z
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: prod_p,prod_m,prod
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: fun_coef
  DOUBLE PRECISION, DIMENSION(:,:),   ALLOCATABLE :: g_plus_beg
  DOUBLE PRECISION, DIMENSION(:,:),   ALLOCATABLE :: f_plus_beg
  DOUBLE PRECISION, DIMENSION(:,:),   ALLOCATABLE :: g_minus_end
  DOUBLE PRECISION, DIMENSION(:,:),   ALLOCATABLE :: f_minus_end
  DOUBLE PRECISION, DIMENSION(:,:),   ALLOCATABLE :: f_plus,f_minus
  DOUBLE PRECISION, DIMENSION(:,:),   ALLOCATABLE :: beta
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: g_plus_end
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: f_plus_end
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: renorm_c
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: alpha
  DOUBLE PRECISION, DIMENSION(:,:),   ALLOCATABLE :: alambd
  DOUBLE COMPLEX,   DIMENSION(:,:),   ALLOCATABLE :: Vg_vp_over_B
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: redmat_forw_l,redmat_back_l
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: redmat_forw_r,redmat_back_r
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: delta_eta
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: enu_coef        !!!term[1]
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: enu_coef2       !!!NTV
  INTEGER :: km1,kp1,m,m1,kb,kb1,ke,ke1,mfactorial,nplp1,k_bound,k1
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE ::  source_pow !!!term[2]
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE ::  conv_pow !!!term[3]
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE ::  alampow !!!term[3]
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE ::  vrecurr !!!term[3]
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE ::  dellampow !!!term[3]
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE ::  dellampow2 !!!NTV
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE ::  convol_polpow !!!term[3]
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE ::  coefleg      !!!terms[2,3]
!
  INTEGER :: ntotsize,ntotsize_m,nts_r,nts_l,kk,kk1,isplit
!
  DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE :: derivs_plot,fun_write
  DOUBLE PRECISION, DIMENSION(:,:,:),   ALLOCATABLE :: alpha_plot,beta_plot
  DOUBLE PRECISION, DIMENSION(:,:,:),   ALLOCATABLE :: f_plus_plot,g_plus_plot
  DOUBLE PRECISION, DIMENSION(:,:,:),   ALLOCATABLE :: fun_lambda
  DOUBLE PRECISION, DIMENSION(:,:),     ALLOCATABLE :: fun_pl,fun_mr
  INTEGER :: iplot,nphiplot,iunit_phi,iunit_sizes,isig,nlam
  INTEGER :: iunit_dt_p,iunit_dt_m
  INTEGER :: iunit_sp_p,iunit_sp_m
  INTEGER :: iunit_et_p,iunit_et_m
  INTEGER, DIMENSION(0:3) :: iunitp,iunitm,iunit_lam
  DOUBLE PRECISION :: phiplot,delphiplot,deleta_b_min,facnorm_p,facnorm_m
  DOUBLE PRECISION :: boundlayer_ignore
  INTEGER :: ignore_lb,ignore_rb,ignore_lb_out,ignore_rb_out,modify_bl,modify_br
  DOUBLE PRECISION :: bhat_changed_l,bhat_changed_r
  DOUBLE PRECISION :: bhat_changed_l_out,bhat_changed_r_out
  DOUBLE PRECISION :: stiffpar,delphi_stfp,phi_stfp_prev,phinext
  DOUBLE PRECISION :: sign_of_bphi                                     !08.12.08
  INTEGER :: nsplit_stfp,npass_prev
  INTEGER, DIMENSION(:),   ALLOCATABLE :: irkstep_stfp
  INTEGER :: icounter
!
  CHARACTER(len=100) :: propname
  INTEGER :: n_2d_size,nrow,ncol,iopt,nz,nz_sq,nz_beg,npassing_prev,k_prev,mm
  INTEGER :: iter,niter,nphiequi,npassing_next
  INTEGER :: isw_regper,nz_coll,nz_ttmp,nz_coll_beg
  DOUBLE PRECISION :: delphim1,deloneovb,step_factor_p,step_factor_m
  DOUBLE PRECISION :: epserr_iter,deleta_factor
  DOUBLE COMPLEX   :: epserr_sink
  INTEGER,          DIMENSION(:),   ALLOCATABLE :: ind_start,irow,icol,ipcol
  INTEGER,          DIMENSION(:),   ALLOCATABLE :: irow_coll,icol_coll
  INTEGER,          DIMENSION(:),   ALLOCATABLE :: irow_ttmp,icol_ttmp
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: funsol_p,funsol_m
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: amat_coll,amat_ttmp
  DOUBLE COMPLEX,   DIMENSION(:),   ALLOCATABLE :: amat_sp,bvec_sp
  DOUBLE COMPLEX,   DIMENSION(:),   ALLOCATABLE :: bvec_iter,bvec_prev
  DOUBLE COMPLEX,   DIMENSION(:,:), ALLOCATABLE :: flux_vector,source_vector
  DOUBLE COMPLEX,   DIMENSION(:,:), ALLOCATABLE :: basevec_p
  INTEGER :: isw_lor,isw_ene,isw_intp
  INTEGER,          DIMENSION(:),       ALLOCATABLE :: npl
  DOUBLE PRECISION, DIMENSION(:,:,:),   ALLOCATABLE :: rhs_mat_fzero
  DOUBLE PRECISION, DIMENSION(:,:,:),   ALLOCATABLE :: rhs_mat_lorentz
  DOUBLE PRECISION, DIMENSION(:,:,:),   ALLOCATABLE :: ttmp_mat
  DOUBLE COMPLEX,   DIMENSION(:,:,:),   ALLOCATABLE :: q_rip
  DOUBLE COMPLEX,   DIMENSION(:,:),     ALLOCATABLE :: q_rip_1
  DOUBLE PRECISION, DIMENSION(:,:,:),   ALLOCATABLE :: rhs_mat_energ
  DOUBLE PRECISION, DIMENSION(:,:,:),   ALLOCATABLE :: rhs_mat_energ2     !NTV
  DOUBLE PRECISION, DIMENSION(:,:,:),   ALLOCATABLE :: pleg_bra,pleg_ket
  DOUBLE COMPLEX,   DIMENSION(:,:),     ALLOCATABLE :: convol_flux,convol_curr
  DOUBLE COMPLEX,   DIMENSION(:,:),     ALLOCATABLE :: convol_flux_0
  DOUBLE PRECISION, DIMENSION(:,:),     ALLOCATABLE :: scalprod_pleg
  DOUBLE PRECISION, DIMENSION(:,:),     ALLOCATABLE :: x1mm,x2mm
  DOUBLE COMPLEX,   DIMENSION(:), ALLOCATABLE :: scalprod
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: phi_mfl
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: bhat_mfl,h_phi_mfl
  DOUBLE COMPLEX,   DIMENSION(:), ALLOCATABLE :: geodcu_mfl
  DOUBLE COMPLEX,   DIMENSION(:), ALLOCATABLE :: geodcu_forw,geodcu_back
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: dlogbdphi_mfl
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: delt_pos,delt_neg
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: fact_pos_b,fact_neg_b
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: fact_pos_e,fact_neg_e
  INTEGER          :: nreal,ncomp
  DOUBLE COMPLEX   :: expforw,expbackw,perbou_pos,perbou_neg,rotfactor
  DOUBLE PRECISION :: a1b,a2b,hatOmegaE,hatOmegaB,denomjac
  !! Modifications by Andreas F. Martitsch (14.03.2014)
  ! Subsequent quantities are given now in cgs-units and they are 
  ! renormalized using bmod0 within neo_magfie: 
  !DOUBLE PRECISION :: bcovar_theta,bcovar_phi,dbcovar_theta_ds,dbcovar_phi_ds
  ! For this reason these variables are renamed:
  DOUBLE PRECISION :: bcovar_theta_hat,bcovar_phi_hat
  DOUBLE PRECISION :: dbcovar_theta_hat_ds,dbcovar_phi_hat_ds
  !! End Modifications by Andreas F. Martitsch (14.03.2014)
  DOUBLE PRECISION :: scalefac_kG
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: arr_real
  DOUBLE COMPLEX,   DIMENSION(:,:), ALLOCATABLE :: arr_comp
  !! Modifications by Andreas F. Martitsch (13.06.2014)
  ! Subsequent quantities (given now in cgs-units) are computed by 
  ! magdata_for_particles and stored within the fieldpropagator-structure.
  ! This step required changes within neo_magfie, magfie, mag,
  ! magdata_for_particles, mag_interface_mod, plagrange_mod,
  ! modify_propagator and magnetics_mod.
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: dlogbds_mfl
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: bcovar_s_hat_mfl
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: dbcovar_s_hat_dphi_mfl
  DOUBLE COMPLEX,   DIMENSION(:),   ALLOCATABLE :: bnoverb0,dbnoverb0_dphi_mfl
  ! For testing you can specify here an artificial perturbation field
  DOUBLE COMPLEX,   DIMENSION(:),   ALLOCATABLE :: bnoverb0_test,dbnoverb0_dphi_mfl_test
  ! amplitude
  DOUBLE COMPLEX :: bnoverb0_test_val=(1.0d-3,0.0d-0)
  ! poloidal mode number
  INTEGER :: m_theta = 0 
  !! End Modifications by Andreas F. Martitsch (13.06.2014)
  DOUBLE COMPLEX,   DIMENSION(:),   ALLOCATABLE :: ttmpfact
  LOGICAL :: colltest=.FALSE.
  LOGICAL :: ttmptest=.FALSE.
!  logical :: ttmptest=.true.
  LOGICAL :: nobounceaver=.TRUE.
!  logical :: nobounceaver=.false.
! 
  ! integer :: isw_axisymm=0 ! now in collisionality_mod
 
  niter=100000       !maximum number of integral part iterations
  epserr_iter=1.d-5 !5  !relative error of integral part iterations
  precond=.TRUE.     !preconditioner for integral part iterations
  isw_regper=1       !regulariization by periodic boundary condition
  epserr_sink=0.d0 !1.d-12 !sink for regularization, it is equal to
!                    $\nu_s/(\sqrt{2} v_T \kappa)$ where
!                    $\bu_s$ is sink rate, $v_T=\sqrt{T/m}$, and
!                    $\kappa$ is inverse m.f.p. times 4 ("collpar")
  sparse_solve_method = 3 !2 !2,3 - with and without iterative refinement, resp.
!
  !! Modifications by Andreas F. Martitsch (14.07.2015)
  ! normalized electric rotation frequency ($\hat{\Omega}_{tE}$)
  ! specified via toroidal Mach number over R_major (Mt/R)
  !conversion conl_over_mfp to collpar by a factor 2 wrong for
  !the full collision operator (therefore numbers for Mt in the paper were rescaled)
  !hatOmegaE=Mt*PI/(conl_over_mfp)
  !this shows the factor 2
  !(results are identical to the ones above; collpar=$\frac{2}{v_{Ta}\tau_{aa}}$)
  !hatOmegaE=2.0d0*Mt/(collpar*(device%r0))
  !this is the correct definition with collpar=4/$l_c$ ($l_c=2 v_{Ta} \tau_{aa}$)
  hatOmegaE=MtOvR/collpar
  ! normalized magnetic rotation frequency ($\hat{\Omega}_{tB}^{\rm ref}$)
  ! specified via Larmor radius associated with $B_{00}^{Booz}$ (rho_L_loc)
  !definition with conl_over_mfp suffers again from wrong conversion to
  !collpar
  !hatOmegaB=((device%r0*PI)/(2.0d0*boozer_psi_pr_hat))*&
  !     (B_rho_L_loc/(avbhat*conl_over_mfp))
  !correct definition with collpar=4/$l_c$ ($l_c=2 v_{Ta} \tau_{aa}$)
  hatOmegaB=B_rho_L_loc/(2.0d0*boozer_psi_pr_hat*(bmod0*1.0d4)*collpar)
  ! print normalized rotation frequencies
  PRINT *,'hatOmegaB,hatOmegaE: ',hatOmegaB,hatOmegaE
  !STOP
  !! End Modifications by Andreas F. Martitsch (14.07.2015)
!  
  !------------------------------------------------------------------------
  ! END SERGEI
  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
  ! Winny: sparse
  !------------------------------------------------------------------------
  ! if sparse_solve should talk
  !  sparse_talk = .TRUE. ! default .FALSE. - neo2.in - settings
  !
  ! sparse method - only 1 (SuperLU) implemented
  !  sparse_solve_method = 1 ! default 0 - neo2.in - settings
  ! if sparse_solve is called with (sparse_solve_method .eq. 0)
  !  program stops
  !  one should call a normal solver in this case
  !
  ! These are three possibilities to call sparse_solve
  !  CALL sparse_solve(nrow,ncol,nz,irow,icol,val,b)  ! full column index
  !  CALL sparse_solve(nrow,ncol,nz,irow,pcol,val,b)  ! column pointer
  !  CALL sparse_solve(A,b)                           ! full matrix
  ! results are returned in b
  !
  ! Input
  !  INTEGER :: nrow,ncol,nz,nrhs
  !  INTEGER, DIMENSION(:), ALLOCATABLE :: irow,icol,pcol
  !  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: val
  !  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: A
  !
  ! In/output
  !  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: b
  ! or
  !  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: b
  !
  ! Matrix - square
  !  ncol = nrow
  !  ALLOCATE( irow(nz) )     ! row index
  !  ALLOCATE( icol(nz) )     ! column index (column ordering)
  !  ALLOCATE( pcol(ncol+1) ) ! column pointer (other possibility)
  !  ALLOCATE( val(nz) )      ! values
  !  ALLOCATE( A(nrow,ncol) ) ! full matrix will be converted to sparse
  !
  ! rhs
  !  ALLOCATE( b(nrow) )
  ! or
  !  ALLOCATE( b(nrow,nrhs) )
  !
  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
  ! reconstruction work by Winny
  !------------------------------------------------------------------------
  ! values of reconstructed fluxes are put into fun_mr and fun_pl
  ! allocation in the 2nd-dimension is from 0 on
  ! (as requested by Sergei)
  IF (prop_reconstruct .EQ. 2) THEN
     prop_ripple_plot = 1
!     ALLOCATE(fun_mr(LBOUND(flux_mr,1):UBOUND(flux_mr,1),0:UBOUND(flux_mr,2)))
!     fun_mr = 0.0_dp
!     fun_mr(LBOUND(flux_mr,1):UBOUND(flux_mr,1),1:UBOUND(flux_mr,2)) = flux_mr
!     ALLOCATE(fun_pl(LBOUND(flux_pl,1):UBOUND(flux_pl,1),0:UBOUND(flux_pl,2)))
!     fun_pl = 0.0_dp
!     fun_pl(LBOUND(flux_pl,1):UBOUND(flux_pl,1),1:UBOUND(flux_pl,2)) = flux_pl
  END IF
  !------------------------------------------------------------------------
  ! end reconstruction work by Winny
  !------------------------------------------------------------------------

  IF (isw_momentum .EQ. 1) THEN ! Grid
     PRINT *, 'isw_momentum = ',isw_momentum,' not implemented in ripple solver!'
     PRINT *, 'I stop here'
     STOP
  END IF



  ierr=0
  !accurfac=3.d0 !10.d0
  accurfac=ripple_solver_accurfac ! in develpment
  renormfac=1.d5 !1.d3 !1d1
  deleta_b_min=1.d-12
  boundlayer_ignore=0.01d0
!
  iunit_phi=200
  iunit_sizes=201
  iunit_dt_p=300
  iunit_dt_m=301
  iunit_sp_p=302
  iunit_sp_m=303
  iunit_et_p=304
  iunit_et_m=305
!
  IF(isw_lorentz.EQ.1) THEN
    lag=0
    nvel=0
    legmax=2            !Caution: legmax cannot be less than 2 (for convol_flux)
    isw_lor=1
    isw_ene=0
    isw_intp=0
    anumm(0,0)=1.d0
    asource(0,1:3)=1.d0
    weightlag(1:3,0)=1.d0
!    nobounceaver=.false.
  ELSE
    legmax=MAX(leg,2)!Caution: legmax cannot be less than 2 (for convol_flux)
    isw_lor=1
    isw_ene=isw_energy
    isw_intp=isw_integral
  ENDIF
!
  nvelocity=nvel
!

iprintflag=1
  !------------------------------------------------------------------------
  ! here i get everything for you what you might need
  !

  ! eta related stuff
  ub_eta = UBOUND(fieldpropagator%ch_act%eta,1)
  npart  = ub_eta
  add_global_eta = 0
  ! IF (bsfunc_local_solver .EQ. 0 .or. bsfunc_local_solver .GT. 2) THEN  !<- SOLVER
  IF (prop_reconstruct_levels .EQ. 0) THEN
     ! This is the old stuff
     ! here eta = eta_glob = eta_loc
     ub_eta_loc = ub_eta   

     ! eta
     IF (ALLOCATED(eta)) DEALLOCATE(eta)
     ALLOCATE(eta(0:ub_eta))
     eta = fieldpropagator%ch_act%eta
     ! eta_glob = eta
     IF (ALLOCATED(eta_glob)) DEALLOCATE(eta_glob)
     ALLOCATE(eta_glob(0:ub_eta))
     eta_glob = fieldpropagator%ch_act%eta

     ! eta_loc = eta
     IF (ALLOCATED(eta_loc)) DEALLOCATE(eta_loc)
     ALLOCATE(eta_loc(0:ub_eta_loc))
     eta_loc = fieldpropagator%ch_act%eta

     ! index
     IF (ALLOCATED(eta_glob_ind)) DEALLOCATE(eta_glob_ind)
     ALLOCATE(eta_glob_ind(0:ub_eta))
     IF (ALLOCATED(eta_loc_ind)) DEALLOCATE(eta_loc_ind)
     ALLOCATE(eta_loc_ind(0:ub_eta_loc))
     DO i = 0,ub_eta
        eta_loc_ind(i) = i
        eta_glob_ind(i) = i
     END DO

  ! ELSEIF (bsfunc_local_solver .EQ. 1 .OR. bsfunc_local_solver .EQ. 2) THEN  !<- SOLVER
  ELSE
     ! This is the new stuff from Winny
     !
     ! local: important for the local propagator
     ! additional: additional stuff which should be reconstructed after solving the local stuff
     ! global = local+additional: all together, used also for allocation of output quantities
     ! 
     ! eta is the global eta
     ! In addition eta_loc (local) and eta_glob (global) exist
     !
     ! All procedures after the solver (e.g. joining do not know about this local 
     ! treatment of eta. So before one returns to the calling routine, all relevant
     ! external quantities have to be interpolated to match the global eta!
     !
     ! compared to eta_loc there are extra levels added in eta = eta_glob, so no levels
     ! are missing there
     !
     ! Now there exist the following new quantities
     !
     ! add_global_eta - integer: 
     !    1 - if additional levels exist; 0 - otherwise
     !
     ! eta_loc_ind(0:ub_eta_loc) - integer:
     !    index where the value of eta_loc can be found in eta_glob
     !
     ! eta_glob_ind(0:ub_eta) - integer:
     !    index where the value of eta_glob can be found in eta_loc
     !    or -1 if it does not exist in eta_loc
     !
     ! reconstruction has to be done, 
     !    if (add_global_eta .gt. 0) 
     ! 
     ! So the usage is now (input file):
     !
     !  prop_reconstruct_levels = 0 : no reconstruction
     !                            1 : reconstruction (if proper bsfunc_local_solver value)
     !
     !  bsfunc_local_solver = 0 : original stuff (no reconstruction possible)
     !                        1 : only local (no reconstruction possible)
     !                        2 : local + "absolute maximum" (reconstruction possible)
     !                        3 : local + "absolute maximum" + rest (reconstruction possible)
     !                        4 : not recommended (no reconstruction possible) 
     !     if reconstruction is not possible, prop_reconstruct_levels is automatically set to 0
     !
     !  bsfunc_sigma_mult .ne. 1.0d0 : then splitting is done with the original value of sigma and 
     !                                 then also with the multiplied value
     !                                 sigma * bsfunc_sigma_mult
     !
     !  isw_lorentz  = 1 : Lorentz operator only (default)
     !                 0 : Lorentz + either integral part or energy diffusion or both
     !  isw_integral = 0 : no integral part (default)
     !                 1 : integral part
     !  isw_energy   = 0 : no energy diffusion (default)
     !                 1 : energy diffusion

     ub_eta_loc = UBOUND(fieldpropagator%ch_act%eta_loc,1)
     ! eta
     IF (ALLOCATED(eta)) DEALLOCATE(eta)
     ALLOCATE(eta(0:ub_eta))
     eta = fieldpropagator%ch_act%eta

     ! eta_glob = eta
     IF (ALLOCATED(eta_glob)) DEALLOCATE(eta_glob)
     ALLOCATE(eta_glob(0:ub_eta))
     eta_glob = fieldpropagator%ch_act%eta

     ! eta_loc
     IF (ALLOCATED(eta_loc)) DEALLOCATE(eta_loc)
     ALLOCATE(eta_loc(0:ub_eta_loc))
     eta_loc = fieldpropagator%ch_act%eta_loc
     
     ! index
     IF (ALLOCATED(eta_glob_ind)) DEALLOCATE(eta_glob_ind)
     ALLOCATE(eta_glob_ind(0:ub_eta))
     IF (ALLOCATED(eta_loc_ind)) DEALLOCATE(eta_loc_ind)
     ALLOCATE(eta_loc_ind(0:ub_eta_loc))
     i_loc = 0
     DO i = 0,ub_eta
        IF ( eta_glob(i) .EQ. eta_loc(i_loc) ) THEN
           eta_loc_ind(i_loc) = i
           eta_glob_ind(i) = i_loc
           i_loc = i_loc + 1
        ELSE
           eta_glob_ind = -1
           add_global_eta = 1
        END IF
     END DO

  END IF  !<- SOLVER
  xetami = eta(0)
  xetama = eta(ub_eta)

  ! previous eta
  IF (ASSOCIATED(fieldpropagator%prev)) THEN
     fieldripple => fieldpropagator%prev%ch_act
  ELSE
     fieldripple => fieldpropagator%parent%ch_las%ch_act
  END IF
  ub_eta_prev = UBOUND(fieldripple%eta,1)
  IF (ALLOCATED(eta_prev)) DEALLOCATE(eta_prev)
  ALLOCATE(eta_prev(0:ub_eta_prev))
  eta_prev = fieldripple%eta
  ! and next eta
  IF (ASSOCIATED(fieldpropagator%next)) THEN
     fieldripple => fieldpropagator%next%ch_act
  ELSE
     fieldripple => fieldpropagator%parent%ch_fir%ch_act
  END IF
  ub_eta_next = UBOUND(fieldripple%eta,1)
  IF (ALLOCATED(eta_next)) DEALLOCATE(eta_next)
  ALLOCATE(eta_next(0:ub_eta_next))
  eta_next = fieldripple%eta
  ! fieldripple back to original
  fieldripple => fieldpropagator%ch_act

  ! bhat (eta) on the left and right side of propagator
  b_l = fieldpropagator%b_l
  b_r = fieldpropagator%b_r
  eta_l = 1.0_dp / b_l
  eta_r = 1.0_dp / b_r

  ! additional stuff
  ! device
  rt0 = fieldpropagator%parent%parent%parent%parent%r0
  ! fieldpropagator
  b_prop_l   = fieldpropagator%b_l
  b_prop_r   = fieldpropagator%b_r
  b_prop_min = fieldpropagator%b_min

  ! fieldripple
  bin_split_mode = fieldpropagator%ch_act%bin_split_mode
  b_max_l = fieldpropagator%ch_act%b_max_l
  b_max_r = fieldpropagator%ch_act%b_max_r
  b_min   = fieldpropagator%ch_act%b_min
  width   = fieldpropagator%ch_act%width


  DO i = 1,ub_eta
     IF(1.d0-b_l*eta(i)+10.d0*EPSILON(1.d0).GT.0.d0) npass_l = i
     IF(1.d0-b_r*eta(i)+10.d0*EPSILON(1.d0).GT.0.d0) npass_r = i
  END DO
  DO i=1,ub_eta_prev
    IF(1.d0-b_l*eta_prev(i)+10.d0*EPSILON(1.d0).GT.0.d0) npass_l_out = i
  ENDDO
  DO i=1,ub_eta_next
    IF(1.d0-b_r*eta_next(i)+10.d0*EPSILON(1.d0).GT.0.d0) npass_r_out = i
  ENDDO
!  
! Ignore the boundary layer if it is too narrow
ignore_lb=0
bhat_changed_l=0.d0
ignore_lb_out=0
bhat_changed_l_out=0.d0
modify_bl=0
ignore_rb=0
bhat_changed_r=0.d0
ignore_rb_out=0
bhat_changed_r_out=0.d0
modify_br=0
GOTO 10
!
! Left boundary:
!
! check own eta-levels
  IF(eta_l-eta(npass_l) .LT.                                         &
    (eta(npass_l)-eta(npass_l-1))*boundlayer_ignore) THEN
    ignore_lb=1
    bhat_changed_l=1.d0/eta(npass_l)+100.d0*EPSILON(1.d0)
  ELSE
    ignore_lb=0
    bhat_changed_l=0.d0
  ENDIF
!
! check outer eta-levels
  IF(eta_l-eta_prev(npass_l_out) .LT.                                &
    (eta_prev(npass_l_out)-eta_prev(npass_l_out-1))*boundlayer_ignore) THEN
    ignore_lb_out=1
    bhat_changed_l_out=1.d0/eta_prev(npass_l_out)+100.d0*EPSILON(1.d0)
  ELSE
    ignore_lb_out=0
    bhat_changed_l_out=0.d0
  ENDIF
!
! forbid bhat modification if a regular band is eliminated in addition to b.l.
  IF(1.d0-bhat_changed_l*eta_prev(npass_l_out-1)+10.d0*EPSILON(1.d0) &
     .LE.0.d0) THEN
    ignore_lb=0
    bhat_changed_l=0.d0
    PRINT *,'cannot ignore left boundary layer: jump over normal band'
  ENDIF
  IF(1.d0-bhat_changed_l_out*eta(npass_l-1)+10.d0*EPSILON(1.d0)      &
     .LE.0.d0) THEN
    ignore_lb_out=0
    bhat_changed_l_out=0.d0
    PRINT *,'cannot ignore right boundary layer: jump over normal band'
  ENDIF
!
! final value of modified bhat
  IF(ignore_lb.EQ.1 .OR. ignore_lb_out.EQ.1) THEN
    bhat_changed_l=MAX(bhat_changed_l,bhat_changed_l_out)
    modify_bl=1
PRINT *,'field at the left boundary modified'
  ELSE
    modify_bl=0
  ENDIF
!
! final decision on the boundary layer
  IF(modify_bl.EQ.1) THEN
    IF(1.d0-bhat_changed_l*eta(npass_l)+10.d0*EPSILON(1.d0).LE.0.d0) THEN
      ignore_lb=1
PRINT *,'left boundary layer ignored'
    ELSE
      ignore_lb=0
    ENDIF
  ENDIF
!
! Right boundary:
!
! check own eta-levels
  IF(eta_r-eta(npass_r) .LT.                                         &
    (eta(npass_r)-eta(npass_r-1))*boundlayer_ignore) THEN
    ignore_rb=1
    bhat_changed_r=1.d0/eta(npass_r)+100.d0*EPSILON(1.d0)
  ELSE
    ignore_rb=0
    bhat_changed_r=0.d0
  ENDIF
!
! check outer eta-levels
  IF(eta_r-eta_next(npass_r_out) .LT.                                &
    (eta_next(npass_r_out)-eta_next(npass_r_out-1))*boundlayer_ignore) THEN
    ignore_rb_out=1
    bhat_changed_r_out=1.d0/eta_next(npass_r_out)+100.d0*EPSILON(1.d0)
  ELSE
    ignore_rb_out=0
    bhat_changed_r_out=0.d0
  ENDIF
!
! forbid bhat modification if a regular band is eliminated in addition to b.l.
  IF(1.d0-bhat_changed_r*eta_next(npass_r_out-1)+10.d0*EPSILON(1.d0) &
     .LE.0.d0) THEN
    ignore_rb=0
    bhat_changed_r=0.d0
    PRINT *,'cannot ignore right boundary layer: jump over normal band'
  ENDIF
  IF(1.d0-bhat_changed_r_out*eta(npass_r-1)+10.d0*EPSILON(1.d0)      &
     .LE.0.d0) THEN
    ignore_rb_out=0
    bhat_changed_r_out=0.d0
    PRINT *,'cannot ignore left boundary layer: jump over normal band'
  ENDIF
!
! final value of modified bhat
  IF(ignore_rb.EQ.1 .OR. ignore_rb_out.EQ.1) THEN
    bhat_changed_r=MAX(bhat_changed_r,bhat_changed_r_out)
    modify_br=1
PRINT *,'field at the right boundary modified'
  ELSE
    modify_br=0
  ENDIF
!
! final decision on the boundary layer
  IF(modify_br.EQ.1) THEN
    IF(1.d0-bhat_changed_r*eta(npass_r)+10.d0*EPSILON(1.d0).LE.0.d0) THEN
      ignore_rb=1
PRINT *,'right boundary layer ignored'
    ELSE
      ignore_rb=0
    ENDIF
  ENDIF
10 CONTINUE
!
  ! place for boundary
  npass_l = npass_l + 1 - ignore_lb
  npass_r = npass_r + 1 - ignore_rb

  ! allocate and copy the magnetic stuff
  ub_mag = UBOUND(fieldpropagator%coords%x2,1)
  ibeg   = 0 - 2*modify_bl
  iend   = ub_mag + 2*modify_br
PRINT *,ub_mag,ibeg,iend

  IF (ALLOCATED(phi_divide)) DEALLOCATE(phi_divide)           !<-in Winny
  ALLOCATE(phi_divide(1:ub_mag))                              !<-in Winny
  phi_divide = 2                                              !<-in Winny

  IF (ALLOCATED(phi_mfl)) DEALLOCATE(phi_mfl)
  ALLOCATE(phi_mfl(ibeg:iend))
  phi_mfl(0:ub_mag) = fieldpropagator%coords%x2
  phibeg = phi_mfl(0)
  phiend = phi_mfl(ub_mag)

  IF (ALLOCATED(bhat_mfl)) DEALLOCATE(bhat_mfl)
  ALLOCATE(bhat_mfl(ibeg:iend))
  bhat_mfl(0:ub_mag) = fieldpropagator%mdata%bhat

  IF (ALLOCATED(geodcu_mfl)) DEALLOCATE(geodcu_mfl)
  ALLOCATE(geodcu_mfl(ibeg:iend))
  geodcu_mfl(0:ub_mag) = fieldpropagator%mdata%geodcu
!
  IF (ALLOCATED(geodcu_forw)) DEALLOCATE(geodcu_forw)
  ALLOCATE(geodcu_forw(ibeg:iend))
  IF (ALLOCATED(geodcu_back)) DEALLOCATE(geodcu_back)
  ALLOCATE(geodcu_back(ibeg:iend))

  IF (ALLOCATED(h_phi_mfl)) DEALLOCATE(h_phi_mfl)
  ALLOCATE(h_phi_mfl(ibeg:iend))
  h_phi_mfl(0:ub_mag) = fieldpropagator%mdata%h_phi
  sign_of_bphi= SIGN(1.d0,h_phi_mfl(0))                                !08.12.08
  h_phi_mfl(0:ub_mag)=h_phi_mfl(0:ub_mag)*sign_of_bphi                 !08.12.08
  geodcu_mfl(0:ub_mag)=geodcu_mfl(0:ub_mag)*sign_of_bphi               !08.12.08

  IF (ALLOCATED(dlogbdphi_mfl)) DEALLOCATE(dlogbdphi_mfl)
  ALLOCATE(dlogbdphi_mfl(ibeg:iend))
  dlogbdphi_mfl(0:ub_mag) = fieldpropagator%mdata%dlogbdphi

  !! Modifications by Andreas F. Martitsch (14.03.2014)
  ! Optional output (necessary for modeling the magnetic rotation)
  IF (ALLOCATED(dlogbds_mfl)) DEALLOCATE(dlogbds_mfl)
  ALLOCATE(dlogbds_mfl(ibeg:iend))
  dlogbds_mfl(0:ub_mag) = fieldpropagator%mdata%dlogbds
  !
  IF (ALLOCATED(bcovar_s_hat_mfl)) DEALLOCATE(bcovar_s_hat_mfl)
  ALLOCATE(bcovar_s_hat_mfl(ibeg:iend))
  bcovar_s_hat_mfl(0:ub_mag) = fieldpropagator%mdata%bcovar_s_hat
  !
  IF (ALLOCATED(dbcovar_s_hat_dphi_mfl)) DEALLOCATE(dbcovar_s_hat_dphi_mfl)
  ALLOCATE(dbcovar_s_hat_dphi_mfl(ibeg:iend))
  dbcovar_s_hat_dphi_mfl(0:ub_mag) = fieldpropagator%mdata%dbcovar_s_hat_dphi
  !! End Modifications by Andreas F. Martitsch (14.03.2014)
  
  !! Modifications by Andreas F. Martitsch (13.06.2014)
  ! Perturbation field extracted from Boozer file
  !IF (ALLOCATED(bnoverb0)) DEALLOCATE(bnoverb0)
  !ALLOCATE(bnoverb0(ibeg:iend))
  !bnoverb0(0:ub_mag) = fieldpropagator%mdata%bnoverb0
  ! Derivative pf the perturbation field extracted from Boozer file
  ! (direct computation without Lagrange interpolation preferred in
  ! order to avoid artefacts nearby the boundary)
  !IF (ALLOCATED(dbnoverb0_dphi_mfl)) DEALLOCATE(dbnoverb0_dphi_mfl)
  !ALLOCATE(dbnoverb0_dphi_mfl(ibeg:iend))
  !dbnoverb0_dphi_mfl(0:ub_mag) = fieldpropagator%mdata%dbnoverb0_dphi
  !! End Modifications by Andreas F. Martitsch (13.06.2014)

  IF(modify_bl.EQ.1) THEN
    phi_mfl(-2:-1)=phi_mfl(0)
    bhat_mfl(-2:-1)=bhat_changed_l
    geodcu_mfl(-2:-1)=geodcu_mfl(0)
    h_phi_mfl(-2:-1)=h_phi_mfl(0)
    dlogbdphi_mfl(-2:-1)=dlogbdphi_mfl(0)
    !! Modifications by Andreas F. Martitsch (14.03.2014)
    ! Optional output (necessary for modeling the magnetic rotation)
    dlogbds_mfl(-2:-1)=dlogbds_mfl(0)
    bcovar_s_hat_mfl(-2:-1)=bcovar_s_hat_mfl(0)
    dbcovar_s_hat_dphi_mfl(-2:-1)=dbcovar_s_hat_dphi_mfl(0)
    !! End Modifications by Andreas F. Martitsch (14.03.2014)
    !! Modifications by Andreas F. Martitsch (13.06.2014)
    ! Perturbation field extracted from Boozer file
    !bnoverb0(-2:-1)=bnoverb0(0)
    ! Derivative pf the perturbation field extracted from Boozer file
    ! (direct computation without Lagrange interpolation preferred in
    ! order to avoid artefacts nearby the boundary)
    !dbnoverb0_dphi_mfl(-2:-1)=dbnoverb0_dphi_mfl(0)
    !! End Modifications by Andreas F. Martitsch (13.06.2014)
  ENDIF
  IF(modify_br.EQ.1) THEN
    phi_mfl(iend-1:iend)=phi_mfl(ub_mag)
    bhat_mfl(iend-1:iend)=bhat_changed_r
    geodcu_mfl(iend-1:iend)=geodcu_mfl(ub_mag)
    h_phi_mfl(iend-1:iend)=h_phi_mfl(ub_mag)
    dlogbdphi_mfl(iend-1:iend)=dlogbdphi_mfl(ub_mag)
    !! Modifications by Andreas F. Martitsch (14.03.2014)
    ! Optional output (necessary for modeling the magnetic rotation)
    dlogbds_mfl(iend-1:iend)=dlogbds_mfl(ub_mag)
    bcovar_s_hat_mfl(iend-1:iend)=bcovar_s_hat_mfl(ub_mag)
    dbcovar_s_hat_dphi_mfl(iend-1:iend)=dbcovar_s_hat_dphi_mfl(ub_mag)
    !! End Modifications by Andreas F. Martitsch (14.03.2014)
    !! Modifications by Andreas F. Martitsch (13.06.2014)
    ! Perturbation field extracted from Boozer file
    !bnoverb0(iend-1:iend)=bnoverb0(ub_mag)
    ! Derivative pf the perturbation field extracted from Boozer file
    ! (direct computation without Lagrange interpolation preferred in
    ! order to avoid artefacts nearby the boundary)
    !dbnoverb0_dphi_mfl(iend-1:iend)=dbnoverb0_dphi_mfl(ub_mag)
    !! End Modifications by Andreas F. Martitsch (13.06.2014)
  ENDIF
  !! Modifications by Andreas F. Martitsch (25.08.2014)
  ! Computation of the perturbed quantities without 
  ! usage of interfaces (fieldpropagator-structure,...)
  IF (ALLOCATED(bnoverb0)) DEALLOCATE(bnoverb0)
  ALLOCATE(bnoverb0(ibeg:iend))
  !
  IF (ALLOCATED(dbnoverb0_dphi_mfl)) DEALLOCATE(dbnoverb0_dphi_mfl)
  ALLOCATE(dbnoverb0_dphi_mfl(ibeg:iend))
  !
  CALL calc_bnoverb0_arr(phi_mfl,ibeg,iend,bnoverb0,dbnoverb0_dphi_mfl)
  CALL calc_ntv_output(phi_mfl,bhat_mfl,bnoverb0,ibeg,iend,&
       eps_M_2_val,av_inv_bhat_val,av_gphph_val)
  !PRINT *,'eps_M_2: ',eps_M_2_val
  !STOP
  !! End Modifications by Andreas F. Martitsch (25.08.2014)
  eta_modboundary_l=1.d0/bhat_mfl(ibeg)
  eta_modboundary_r=1.d0/bhat_mfl(iend)
  ! allocation
  !  at the moment everything from 1:npart
  !  attention eta is from 0:npart-1
  !
  ! 2-D quantities
  nts_l=(lag+1)*npass_l
  nts_r=(lag+1)*npass_r
  IF (ALLOCATED(amat_plus_plus)) DEALLOCATE(amat_plus_plus)
  ALLOCATE(amat_plus_plus(nts_r,nts_l))
  IF (ALLOCATED(amat_minus_minus)) DEALLOCATE(amat_minus_minus)
  ALLOCATE(amat_minus_minus(nts_l,nts_r))
  IF (ALLOCATED(amat_plus_minus)) DEALLOCATE(amat_plus_minus)
  ALLOCATE(amat_plus_minus(nts_l,nts_l))
  IF (ALLOCATED(amat_minus_plus)) DEALLOCATE(amat_minus_plus)
  ALLOCATE(amat_minus_plus(nts_r,nts_r))
!
  IF (ALLOCATED(source_p)) DEALLOCATE(source_p)
  ALLOCATE(source_p(nts_r,3))
  IF (ALLOCATED(source_m)) DEALLOCATE(source_m)
  ALLOCATE(source_m(nts_l,3))
  !
  IF (ALLOCATED(flux_p)) DEALLOCATE(flux_p)
  ALLOCATE(flux_p(3,nts_l))
  IF (ALLOCATED(flux_m)) DEALLOCATE(flux_m)
  ALLOCATE(flux_m(3,nts_r))
  !
  IF (ALLOCATED(qflux)) DEALLOCATE(qflux)
  ALLOCATE(qflux(3,3))

     PRINT *, 'propagator tag         ', fieldpropagator%tag
!solver_talk=1
  IF (solver_talk .EQ. 1) THEN
     PRINT *, ' '
     PRINT *, 'I am in ripple_solver'
     PRINT *, ' '
     PRINT *, 'fieldpropagator tag    ', fieldpropagator%tag
     PRINT *, 'fieldripple tag        ', fieldpropagator%ch_act%tag
     PRINT *, ' fieldprop first last  ', fieldpropagator%ch_act%pa_fir%tag,fieldpropagator%ch_act%pa_las%tag
     PRINT *, ' b_prop left,right,min ', b_prop_l,b_prop_r,b_prop_min
     PRINT *, '                       ', fieldpropagator%mdata%bhat(0),fieldpropagator%mdata%bhat(ub_mag)
     PRINT *, ' b_max left,right,min  ', b_max_l,b_max_r,b_min
     PRINT *, ' width                 ', width
     PRINT *, 'bin_split_mode         ', bin_split_mode
     PRINT *, 'phibeg,phiend          ', phibeg,phiend
     PRINT *, '                       ', fieldpropagator%coords%x2(0),fieldpropagator%coords%x2(ub_mag)
     PRINT *, 'ub_mag                 ', ub_mag
     PRINT *, 'ignore_lb,ignore_rb    ', ignore_lb,ignore_rb
     PRINT *, 'ibeg,iend              ', ibeg,iend
     PRINT *, 'ub_eta                 ', ub_eta
     PRINT *, 'npart                  ', npart
     PRINT *, 'npass_l,npass_r        ', npass_l,npass_r
     PRINT *, 'eta_l,eta_r            ', eta_l,eta_r
     PRINT *, 'xetami,xetama          ', xetami,xetama
     PRINT *, 'rt0                    ', rt0
     PRINT *, 'collpar,conl_over_mfp  ', collpar,conl_over_mfp

     !! Modifications by Andreas F. Martitsch (14.03.2014)
     ! Optional output (necessary for modeling the magnetic rotation)
     OPEN(123,file='bhat_mfl.dat')
     DO i=ibeg,iend
        WRITE(123,*) phi_mfl(i),bhat_mfl(i),REAL(geodcu_mfl(i)),h_phi_mfl(i),dlogbdphi_mfl(i),&
             dlogbds_mfl(i),bcovar_s_hat_mfl(i),dbcovar_s_hat_dphi_mfl(i),&
             REAL(bnoverb0(i)),AIMAG(bnoverb0(i)),REAL(dbnoverb0_dphi_mfl(i)),AIMAG(dbnoverb0_dphi_mfl(i))
     ENDDO
     CLOSE(123)
     !! End Modifications by Andreas F. Martitsch (14.03.2014)

     OPEN(123,file='eta.dat')
     DO i=0,ub_eta
        WRITE(123,*) phi_mfl(ibeg),eta(i)
        WRITE(123,*) phi_mfl(iend),eta(i)
        WRITE(123,*) ' '
     ENDDO
     CLOSE(123)

!     PAUSE 'bmod and eta written'
  END IF
  !------------------------------------------------------------------------
  ! SERGEI
  !------------------------------------------------------------------------
  epserr_sink=(0.d0,0.d0)
!
! Check for axisymmetry:
!
  IF(isw_axisymm.EQ.1.AND.npass_l.NE.npass_r) THEN
    PRINT *,'ripple_solver: cannot run axisymmetric mode, sizes do not fit'
    ierr=1
    RETURN
  ENDIF
!
  iplot=prop_ripple_plot
!
! Preparation of coefficients for the kinetic equation solver
!
  ALLOCATE(deriv_coef(4,0:npart+1))
  ALLOCATE(fun_coef(4,0:npart+1))
  ALLOCATE(enu_coef(4,npart+1))                                    !!!term[1]
  ALLOCATE(enu_coef2(4,npart+1))                                   !!!NTV
  ALLOCATE(alambd(0:npart+3,ibeg:iend),Vg_vp_over_B(0:npart,ibeg:iend))
  ALLOCATE(scalprod_pleg(0:lag,0:legmax))                          !!!term[3]
  ALLOCATE(alampow(legmax+1,0:npart+1))                            !!!terms[2,3]
  ALLOCATE(vrecurr(0:legmax,0:3,1:npart+1))                        !!!terms[2,3]
  ALLOCATE(dellampow(4,1:npart+1))                                 !!!terms[1-3]
  ALLOCATE(dellampow2(4,1:npart+1))                                !!!NTV
  ALLOCATE(convol_polpow(0:legmax,1:npart+3))                      !!!terms[2,3]
  ALLOCATE(pleg_bra(0:legmax,1:npart+1,ibeg:iend))                 !!!terms[2,3]
  ALLOCATE(pleg_ket(0:legmax,1:npart+1,ibeg:iend))                 !!!terms[2,3]
  ALLOCATE(npl(ibeg:iend))
  ALLOCATE(rhs_mat_fzero(4,ibeg:iend,0:1))
  ALLOCATE(rhs_mat_lorentz(5,npart+1,ibeg:iend))
  ALLOCATE(ttmp_mat(5,npart+1,ibeg:iend))
  ALLOCATE(rhs_mat_energ(4,npart+1,ibeg:iend))
  ALLOCATE(rhs_mat_energ2(4,npart+1,ibeg:iend))              !NTV
  ALLOCATE(q_rip(npart+2,ibeg:iend,0:2))
  ALLOCATE(q_rip_1(npart+2,ibeg:iend))
  ALLOCATE(convol_flux(npart+1,ibeg:iend),convol_curr(npart+1,ibeg:iend))
  ALLOCATE(convol_flux_0(npart+1,ibeg:iend))
  ALLOCATE(ind_start(ibeg:iend))
  IF(iplot.EQ.1) ALLOCATE(derivs_plot(0:3,4,npart+1,ibeg:iend))
!
!tmp:
aiota=boozer_iota
!! Modifications by Andreas F. Martitsch (14.03.2014)
! boozer_curr_tor, boozer_curr_pol, boozer_psi_pr,
! boozer_sqrtg11 and boozer_isqrg are now converted
! to cgs-units within neo_magfie.
! This step requires changes within rhs_kin.f90 and
! ripple_solver.f90!
bcovar_theta_hat=boozer_curr_tor_hat
bcovar_phi_hat=boozer_curr_pol_hat
!print *,bcovar_theta_hat,bcovar_phi_hat
!pause
dbcovar_theta_hat_ds=boozer_curr_tor_hat_s
dbcovar_phi_hat_ds=boozer_curr_pol_hat_s
!print *,dbcovar_theta_hat_ds,dbcovar_phi_hat_ds
!pause
!
! These quantities are now available from the
! fieldpropagator-structure.
! Caution! They have also been renamed!
!ALLOCATE(dlogbds_mfl(ibeg:iend))
!ALLOCATE(bcovar_s(ibeg:iend),dbcovar_s_dtheta(ibeg:iend))
!ALLOCATE(bnoverb0(ibeg:iend),dbnoverb0_dtheta(ibeg:iend))
!dlogbds_mfl=0.d0
!bcovar_s=0.d0
!dbcovar_s_dtheta=0.d0
!! End Modifications by Andreas F. Martitsch (14.03.2014)
!
!! Modification by Andreas F. Martitsch (13.06.2014)
!! Block commented out (Old way to specify artificial
!! perturbation field + Derivative over the periodic
!! Boozer angle is now redefined to the derivative 
!! along the field line (phi_mfl))
!bnoverb0=(1.d0,0.d0)
! Read perturbation field amplitude from external
! input file (bnoverb0_Val) and include dependency on
! poloidal mode number (m_theta)
!bnoverb0=bnoverb0_Val*EXP(imun*m_theta*aiota*phi_mfl)
!dbnoverb0_dtheta=(0.d0,0.d0)
! Since perturbation field amplitude is now dependent on m_theta,
! also the derivative over theta gives a finite contribution
! (Note: theta denotes here the periodic Boozer angle and
! not the field-aligned one - Derivative along periodic Boozer angle)
!dbnoverb0_dtheta=imun*m_theta*bnoverb0
!! End Modification by Andreas F. Martitsch (13.06.2014)
!
!! Modification by Andreas F. Martitsch (13.06.2014)
!! Block commented out (For testing you can specify here
!! an artificial perturbation field)
ALLOCATE(bnoverb0_test(ibeg:iend))
ALLOCATE(dbnoverb0_dphi_mfl_test(ibeg:iend))
bnoverb0_test=bnoverb0_test_val*EXP(imun*(m_theta*aiota+m_phi)*phi_mfl)
dbnoverb0_dphi_mfl_test=imun*(m_theta*aiota+m_phi)*bnoverb0_test
!PRINT *,'m_phi: ',m_phi
!DO istep=ibeg,iend
!   PRINT *,'----'
!   PRINT *,bnoverb0(istep),bnoverb0_test(istep)
!   PRINT *,dbnoverb0_dphi_mfl(istep),dbnoverb0_dphi_mfl_test(istep)
!END DO
!bnoverb0=bnoverb0_test
!dbnoverb0_dphi_mfl=dbnoverb0_dphi_mfl_test
DEALLOCATE(bnoverb0_test)
DEALLOCATE(dbnoverb0_dphi_mfl_test)
!! End Modification by Andreas F. Martitsch (13.06.2014)
!
!! Modifications by Andreas F. Martitsch (12.03.2014)
! boozer_curr_tor, boozer_curr_pol, boozer_psi_pr,
! boozer_sqrtg11 and boozer_isqrg are now converted
! to cgs-units within neo_magfie.
! This step requires changes within rhs_kin.f90 and
! ripple_solver.f90!
!             scalefac_kG=$B_{ref}/(\iota \psi_{tor}^a)$
!This was the old version (conversion to cgs is done here):
!scalefac_kG=1d-4*bmod0/(aiota*boozer_psi_pr)
!Now the quantities are already converted within neo_magfie:
scalefac_kG=1.0d0/(aiota*boozer_psi_pr_hat)
!PRINT *,'scalefac_kG: ',scalefac_kG
!! Modifications by Andreas F. Martitsch (12.03.2014)
!
!! Modification by Andreas F. Martitsch (17.12.2013)
!! Block commented out since data is processed by another extra input file
!hatOmegaE=0.d0 !1.d0
!hatOmegaB=0.d0 !2.d2
!OPEN(4321,file='m_phi_Omegas.inp')
!READ(4321,*) m_phi
!READ(4321,*) hatOmegaE
!READ(4321,*) hatOmegaB
!READ(4321,*) m_theta
!CLOSE(4321)
!! End Modifications by Andreas F. Martitsch (17.12.2013)
rotfactor=imun*m_phi
!! Modification by Andreas F. Martitsch (17.12.2013)
!! Block commented out
! --> Electric field is already set to a constant
!hatOmegaE=hatOmegaE/conl_over_mfp
!hatOmegaB=hatOmegaB/conl_over_mfp
! --> Amplitude of the perturbation field is specified by
! --> the external input file
!bnoverb0=EXP(imun*m_theta*aiota*phi_mfl)
!dbnoverb0_dtheta=imun*m_theta*bnoverb0
!! End Modifications by Andreas F. Martitsch (17.12.2013)
!
!end tmp
!
! Compute coefficients of Legendre polynomials of the order 0,...,legmax:
  CALL polleg(legmax,coefleg)
! coefleg(l,k) - coefficient of $x^k$ of polynomial $P_l(x)$
!
  q_rip(:,:,0)=0.d0
  DO i=1,npart
    q_rip(i,:,2)=eta(i)-eta(i-1)
  ENDDO
!
  ndim=4
  ALLOCATE(amat(ndim,ndim),bvec_lapack(ndim,ndim),ipivot(ndim))
!
  npart_loc=0
  subsqmin=1.d5*EPSILON(1.d0)
!
  ALLOCATE(delt_pos(ibeg:iend),delt_neg(ibeg:iend))
  ALLOCATE(fact_pos_b(ibeg:iend),fact_neg_b(ibeg:iend))
  ALLOCATE(fact_pos_e(ibeg:iend),fact_neg_e(ibeg:iend))
!
  nreal=1
  ncomp=1
  ALLOCATE(arr_real(ibeg:iend,nreal),arr_comp(ibeg:iend,ncomp))
  arr_real(:,1)=h_phi_mfl
  arr_comp(:,1)=geodcu_mfl
!
  CALL rearrange_phideps(ibeg,iend,npart,ncomp,nreal,subsqmin,phi_divide, &
                         phi_mfl,bhat_mfl,arr_real,arr_comp,eta,          &
                         delt_pos,delt_neg,                               &
                         fact_pos_b,fact_neg_b,fact_pos_e,fact_neg_e)
!
  h_phi_mfl=arr_real(:,1)
  geodcu_mfl=arr_comp(:,1)
  DEALLOCATE(arr_real,arr_comp)
!
  IF(MAXVAL(phi_divide).GT.1) THEN
    ierr=3
    DEALLOCATE(deriv_coef,npl)
    DEALLOCATE(rhs_mat_lorentz,rhs_mat_energ)
    DEALLOCATE(fun_coef,ttmp_mat)
    DEALLOCATE(rhs_mat_energ2)        !NTV
    DEALLOCATE(q_rip,q_rip_1)
    DEALLOCATE(convol_flux,convol_curr,convol_flux_0)
    DEALLOCATE(pleg_bra,pleg_ket,scalprod_pleg)
    RETURN
  ENDIF
!
  DO istep=ibeg,iend
!
! semi-levels
!
    eta0=1.d0/bhat_mfl(istep)
!
    DO i=0,npart
      subsq=1.d0-bhat_mfl(istep)*eta(i)
      IF(subsq.GT.subsqmin) THEN
        npassing=i
        alambd(i,istep)=SQRT(subsq)
        Vg_vp_over_B(i,istep)=alambd(i,istep)*eta0/h_phi_mfl(istep)            &
!                             *(4.d0*eta0-eta(i))*geodcu_mfl(istep)/3.d0
                             *(4.d0*eta0-eta(i))/3.d0
      ELSE
        alambd(i,istep)=0.d0
        Vg_vp_over_B(i,istep)=0.d0
      ENDIF
    ENDDO
!
    alambd(npart+1:npart+3,istep)=0.d0
    alambd(npassing+1,istep)=0.d0
    alambd(npassing+2,istep)=-alambd(npassing,istep)
    alambd(npassing+3,istep)=-alambd(npassing-1,istep)
!
    npl(istep)=npassing
!
    npart_loc=MAX(npart_loc,npassing)
!
    IF(istep.EQ.ibeg) THEN 
      npass_l=npassing+1
    ELSEIF(istep.EQ.iend) THEN 
      npass_r=npassing+1
    ENDIF
!
  ENDDO
!
! compute starting index for 2D vectors
!
  ind_start(ibeg)=0
!
  DO istep=ibeg,iend-1
    ind_start(istep+1)=ind_start(istep)+2*(lag+1)*(npl(istep)+1)
  ENDDO
!
  n_2d_size=ind_start(iend)+2*(lag+1)*(npl(iend)+1)
!
!
  IF(ALLOCATED(alam_l)) DEALLOCATE(alam_l)
  IF(ALLOCATED(alam_r)) DEALLOCATE(alam_r)
  IF(ALLOCATED(delta_eta_l)) DEALLOCATE(delta_eta_l)
  IF(ALLOCATED(delta_eta_r)) DEALLOCATE(delta_eta_r)
  ALLOCATE(delta_eta(npart_loc),alam_l(npass_l),alam_r(npass_r))
  ALLOCATE(delta_eta_l(npass_l),delta_eta_r(npass_r))
  delta_eta=eta(1:npart_loc)-eta(0:npart_loc-1)
  DO i=1,npass_l-1
    alam_l(i)=SQRT(1.d0-0.5d0*(eta(i-1)+eta(i))*bhat_mfl(ibeg))
    delta_eta_l(i)=eta(i)-eta(i-1)
  ENDDO
  i=npass_l
  alam_l(i)=SQRT(1.d0-0.5d0*(eta(i-1)*bhat_mfl(ibeg)+1.d0))
  delta_eta_l(i)=1.d0/bhat_mfl(ibeg)-eta(i-1)
  DO i=1,npass_r-1
    alam_r(i)=SQRT(1.d0-0.5d0*(eta(i-1)+eta(i))*bhat_mfl(iend))
    delta_eta_r(i)=eta(i)-eta(i-1)
  ENDDO
  i=npass_r
  alam_r(i)=SQRT(1.d0-0.5d0*(eta(i-1)*bhat_mfl(iend)+1.d0))
  delta_eta_r(i)=1.d0/bhat_mfl(iend)-eta(i-1)
!
!
! Calculation of the ODE coefficients
!
  DO istep=ibeg,iend
!
    npassing=npl(istep)
    eta0=1.d0/bhat_mfl(istep)
    amin2ovb=-2.d0/bhat_mfl(istep)
    coefdir=0.5*collpar/h_phi_mfl(istep)
    coefenu_averb=0.5d0*collpar/h_phi_mfl(istep)           !!!term[1]
    coefenu=-coefenu_averb*2.d0/bhat_mfl(istep)             !!!term[1]
!
!
!-----------------------------
! begin terms[2,3]
!
! here dellampow(m,n)=$(\lambda_{n-1}-\lambda_{n})^m$
    dellampow(1,1:npassing+1)                                                &
         =alambd(0:npassing,istep)-alambd(1:npassing+1,istep)
    DO k=2,4
      km1=k-1
      dellampow(k,1:npassing+1)=dellampow(km1,1:npassing+1)                  &
                               *dellampow(1,1:npassing+1)
    ENDDO
!
    alampow(1,0:npassing+1)=alambd(0:npassing+1,istep)
    DO k=2,legmax+1
      km1=k-1
      alampow(k,0:npassing+1)=alambd(0:npassing+1,istep)                     &
                             *alampow(km1,0:npassing+1)
    ENDDO
!
    DO k=1,legmax+1
! Caution: 
! here power index is shifted (instead of term (k) -> term (k-1) is computed)
      km1=k-1
      vrecurr(km1,0,1:npassing+1)                                            &
         =(alampow(k,0:npassing)-alampow(k,1:npassing+1))/DBLE(k)
      DO m=1,3
        vrecurr(km1,m,1:npassing+1)                                          &
           =(alampow(k,0:npassing)*dellampow(m,1:npassing+1)                 &
           -DBLE(m)*alambd(1:npassing+1,istep)*vrecurr(km1,m-1,1:npassing+1))&
           /DBLE(k+m)
      ENDDO
! divide by factorial
      mfactorial=1
      DO m=1,3
        mfactorial=mfactorial*m
        vrecurr(km1,m,1:npassing+1)=vrecurr(km1,m,1:npassing+1)              &
                                   /DBLE(mfactorial)
      ENDDO
    ENDDO
!
! re-definition: here dellampow(m,n)=$(\lambda_{n-1}-\lambda_{n})^m/m!$
! (divided by factorial)
!
    mfactorial=1
    DO m=2,4
      mfactorial=mfactorial*m
      dellampow(m,1:npassing+1)=dellampow(m,1:npassing+1)/DBLE(mfactorial)
    ENDDO
!
! new stuff: NTV
!
    DO m=1,4
      dellampow2(m,1:npassing+1)=dellampow(m,1:npassing+1)              &
        *(alambd(1:npassing+1,istep)**2                                 &
        + 2.d0*alambd(1:npassing+1,istep)*dellampow(1,1:npassing+1)     &
              *dfloat(m)/dfloat(m+1)                                    &
        + dellampow(1,1:npassing+1)**2*dfloat(m)/dfloat(m+2))
    ENDDO
!
! end new stuff: NTV
!
! term[2] (Legendre polynomials) -  ket-vector
! Caution: ket-vector cooresponds do discretization of 
! P_l(lambda)/|lambda|, not P_l(lambda). Scalar prduct with bra-vector
! does not mean the product of Legendre polynomials!
! even powers:
    DO m=0,legmax,2
      pleg_ket(m,1:npassing+1,istep)=amin2ovb*coefleg(m,0)       &
           *(alampow(1,1:npassing+1)-alampow(1,0:npassing))
      DO k=2,m,2
        kp1=k+1
        pleg_ket(m,1:npassing+1,istep)                           &
           =pleg_ket(m,1:npassing+1,istep)+amin2ovb*coefleg(m,k) &
           *(alampow(kp1,1:npassing+1)-alampow(kp1,0:npassing))/DBLE(kp1)
      ENDDO
    ENDDO
! odd powers:
    DO m=1,legmax,2
      pleg_ket(m,1:npassing+1,istep)=amin2ovb*coefleg(m,1)       &
           *(alampow(2,1:npassing+1)-alampow(2,0:npassing))/2.d0
      DO k=3,m,2
        kp1=k+1
        pleg_ket(m,1:npassing+1,istep)                           &
           =pleg_ket(m,1:npassing+1,istep)+amin2ovb*coefleg(m,k) &
           *(alampow(kp1,1:npassing+1)-alampow(kp1,0:npassing))/DBLE(kp1)
      ENDDO
    ENDDO
!
    convol_polpow=0.d0
!
! end terms[2,3]
!---------------------------------------
!
    DO i=1,npassing+1
!
      i1min=MAX(0,i-2)
!
      kmax=5
!
      DO k=1,kmax
        i1=k-1+i1min
        diflam=alambd(i1,istep)-alambd(i,istep)
        diflampow=diflam
        alp(k)=(alambd(i,istep)+diflam/2.d0)*diflampow
        diflampow=diflam*diflampow
        bet(k)=(alambd(i,istep)/2.d0+diflam/3.d0)*diflampow
        diflampow=diflam*diflampow
        gam(k)=(alambd(i,istep)/6.d0+diflam/8.d0)*diflampow
        diflampow=diflam*diflampow
        del(k)=(alambd(i,istep)/24.d0+diflam/30.d0)*diflampow
      ENDDO
!
      DO k=1,4
        amat(k,1)=(alp(k+1)-alp(k))*amin2ovb
        amat(k,2)=(bet(k+1)-bet(k))*amin2ovb
        amat(k,3)=(gam(k+1)-gam(k))*amin2ovb
        amat(k,4)=(del(k+1)-del(k))*amin2ovb
      ENDDO
!
      IF(i.EQ.npassing) THEN
        amat(4,:)=-amat(4,:)
      ELSEIF(i.EQ.npassing+1) THEN
        amat(3,:)=-amat(3,:)
        amat(4,:)=-amat(4,:)
      ENDIF
!
      bvec_lapack=0.d0
      DO k=1,ndim
        bvec_lapack(k,k)=1.d0
      ENDDO
!
      CALL gbsv(ndim,ndim,amat,ipivot,bvec_lapack,info)
!
! bvec_lapack(j,k) - contribution to the derivative of the distribution 
! function $\hat f^\sigma$ of the order j-1=0,1,2,3 at the boundary 
! $\lambda=\lambda_i$ (at the level $\eta=\eta_i$) from the band i+k-2, 
! where k=1,2,3,4. If i=1 contributing bands are i+k-1=1,2,3,4 (shift up by 1).
! If i=npassing, sigma=-1 fluxes start contribute:
! contributions for k=1,2,3,4 come from fun(N-1),fun(N),fun(N+1),fun_b(N*1)
! If i=npassing+1
! contributions for k=1,2,3,4 come from fun(N),fun(N+1),fun_b(N+1),fun_b(N)
! Actual derivative can be obtained by summation of corresponding
! band-integrated fluxes, $f_{i+k-2}$, multiplied with these contributions
!
!
      IF(iplot.EQ.1) derivs_plot(0:3,1:4,i,istep)=bvec_lapack(1:4,1:4)
      deriv_coef(:,i)=bvec_lapack(2,:)*coefdir*MIN(eta(i),eta0)
      fun_coef(:,i)=bvec_lapack(1,:)*MIN(eta(i),eta0)
!
      enu_coef(:,i)=MATMUL(dellampow(:,i),bvec_lapack)*coefenu      !!!term[1]
      enu_coef2(:,i)=MATMUL(dellampow2(:,i),bvec_lapack)*coefenu    !!!NTV
!
      IF(i.EQ.1) THEN
        convol_polpow(0:legmax,i:i+3)=convol_polpow(0:legmax,i:i+3)          &
                 +MATMUL(vrecurr(0:legmax,0:3,i),bvec_lapack)       !!!term[3]
      ELSE
        convol_polpow(0:legmax,i-1:i+2)=convol_polpow(0:legmax,i-1:i+2)      &
                 +MATMUL(vrecurr(0:legmax,0:3,i),bvec_lapack)       !!!term[3]
      ENDIF
      IF(i.EQ.npassing+1) THEN
! distribution function at the trapped-passing boundary:
! f0=sum(fun(npassing:npassing+3)*rhs_mat_fzero(:,istep,0))
        rhs_mat_fzero(:,istep,0)=bvec_lapack(1,:)
      ENDIF
!
    ENDDO
!
! Eliminate stepping over the boundary:
!
    DO k=0,legmax,2
      convol_polpow(k,npassing)  =convol_polpow(k,npassing)                    &
                                 +convol_polpow(k,npassing+3)
      convol_polpow(k,npassing+1)=convol_polpow(k,npassing+1)                  &
                                 +convol_polpow(k,npassing+2)
    ENDDO
!
    DO k=1,legmax,2
      convol_polpow(k,npassing)  =convol_polpow(k,npassing)                    &
                                 -convol_polpow(k,npassing+3)
      convol_polpow(k,npassing+1)=convol_polpow(k,npassing+1)                  &
                                 -convol_polpow(k,npassing+2)
    ENDDO
!
! term[3] (Legendre polynomials) -  bra-vector
! numbering of levels (N=npassing): 1-N - $f_n$, N+1 - $f^b$, N+2 - $f^a$
! even powers:
    DO m=0,legmax,2
      pleg_bra(m,1:npassing+1,istep)=coefleg(m,0)                &
                                                *convol_polpow(0,1:npassing+1)
      DO k=2,m,2
        pleg_bra(m,1:npassing+1,istep)                           &
           =pleg_bra(m,1:npassing+1,istep)+coefleg(m,k)          &
                                                *convol_polpow(k,1:npassing+1)
      ENDDO
    ENDDO
! odd powers:
    DO m=1,legmax,2
      pleg_bra(m,1:npassing+1,istep)=coefleg(m,1)                &
                                                *convol_polpow(1,1:npassing+1)
      DO k=3,m,2
        pleg_bra(m,1:npassing+1,istep)                           &
           =pleg_bra(m,1:npassing+1,istep)+coefleg(m,k)          &
                                                *convol_polpow(k,1:npassing+1)
      ENDDO
    ENDDO
!
    pleg_bra(0:legmax,1:npassing+1,istep)                        &
           =pleg_bra(0:legmax,1:npassing+1,istep)*coefenu_averb
!
!    coef_cf=geodcu_mfl(istep)/bhat_mfl(istep)**2/h_phi_mfl(istep)
    coef_cf=(1.d0,0.d0)/bhat_mfl(istep)**2/h_phi_mfl(istep)
    convol_flux_0(1:npassing+1,istep)                                        &
           =(convol_polpow(0,1:npassing+1)+convol_polpow(2,1:npassing+1))    &
           *coef_cf
!
! levels
!
    rhs_mat_lorentz(5,1:npassing+1,istep)=0.d0
!
    rhs_mat_lorentz(1:4,1,istep)=deriv_coef(:,1)
    rhs_mat_lorentz(1:4,2,istep)=deriv_coef(:,2)-deriv_coef(:,1)
    rhs_mat_lorentz(1:4,3:npassing+1,istep)                      &
                  =-deriv_coef(:,2:npassing)
    rhs_mat_lorentz(2:5,3:npassing+1,istep)                      &
                 =rhs_mat_lorentz(2:5,3:npassing+1,istep)        &
                 +deriv_coef(:,3:npassing+1)
!
    ttmp_mat(5,1:npassing+1,istep)=0.d0
!
    ttmp_mat(1:4,1,istep)=fun_coef(:,1)
    ttmp_mat(1:4,2,istep)=fun_coef(:,2)-fun_coef(:,1)
    ttmp_mat(1:4,3:npassing+1,istep)=-fun_coef(:,2:npassing)
    ttmp_mat(2:5,3:npassing+1,istep)                             &
                 =ttmp_mat(2:5,3:npassing+1,istep)               &
                 +fun_coef(:,3:npassing+1)
!
    rhs_mat_fzero(:,istep,1)=deriv_coef(:,npassing+1)
!
! Change of the boundary layer width
!
!
!
! begin term[1]:
!
    rhs_mat_energ(:,1:npassing+1,istep)=enu_coef(:,1:npassing+1)
    rhs_mat_energ2(:,1:npassing+1,istep)=enu_coef2(:,1:npassing+1) !NTV
!
! end term[1]
!
    q_rip_1(1:npassing,istep)                                    &
          =Vg_vp_over_B(1:npassing,istep)-Vg_vp_over_B(0:npassing-1,istep)
    q_rip_1(npassing+1,istep)=-Vg_vp_over_B(npassing,istep)
    q_rip(npassing+1,istep,2)=eta0-eta(npassing)
    q_rip(1:npassing+1,istep,2)                                  &
          =q_rip(1:npassing+1,istep,2)                           &
          *bhat_mfl(istep)/h_phi_mfl(istep)
!
    convol_curr(1:npassing+1,istep)=bhat_mfl(istep)/h_phi_mfl(istep)
!
  ENDDO
!
  DEALLOCATE(amat,bvec_lapack,ipivot)
!
!
! Preparation of data for sparce solver
!
! The solution vector has the following structure:
! 1) vector is split in the blocks corresponding to spatial position -
!    "spatial blocks", sequence of blocks is ibeg:iend
! 2) each spatial block is split into "Laguerre harmonic blocks",
!    sequence of these blocks is 0:lag
! 3) Each Laguerre harmonic block contains 2*npassing+2 elements, where
!    npassing is the number of complete bands (i.e. number of levels) at the
!    given spatial point; this block is split in two parts:
! a) co-passing particles:
!    These occupy first npassing+1 elements of the Laguerre block, sequence
!    of these elements is direct, 1:npassing+1 - 1st element corresponds 
!    to f_1, i.e. the flux through the 1st band and last element - to
!    the flux through incoplete band - boundary layer
! b) counter-passing particles:
!    These occupy next npassing+1 elements of the Laguerre block, sequence
!    of these elements is reverse, npassing+1:1:-1 - element npassing+2
!    contains the counter-passing flux through the boundary layer, and
!    element 2*npassing+2 contains the counter-passing flux through the
!    first band.
!
!
! Matrix size:
!
  nrow=n_2d_size
  ncol=n_2d_size
!
! Compute vectors for convolution of fluxes and source vectors:
!
  ALLOCATE(flux_vector(3,n_2d_size),source_vector(n_2d_size,3))
  ALLOCATE(x1mm(0:lag,0:lag))
  ALLOCATE(x2mm(0:lag,0:lag))
!
  CALL lagxmm(lag,x1mm,x2mm)
!
  IF(isw_lorentz.EQ.1) THEN
    x1mm(0,0)=1.d0
    x2mm(0,0)=1.d0
  ENDIF
!
! Determine the size of arrays (number of non-zero elements):
!
  nz=0
  nz_coll=0
  nz_ttmp=0
  nz_regper=0
!
! Co-passing: sigma=1
!
  istep=ibeg
  npassing=npl(istep)
!
! entry:
!
  DO m=0,lag
    k=ind_start(istep)+2*(npassing+1)*m
!
    DO ipart=1,npassing+1
      nz=nz+1
!      irow(nz)=k+ipart
!      icol(nz)=k+ipart
!      amat_sp(nz)=1.d0
    ENDDO
!
    IF(isw_axisymm.EQ.1) THEN
! periodicity:
      k_prev=ind_start(iend)+2*(npassing+1)*m
!
      DO ipart=1,npassing+1
        nz=nz+1
!        irow(nz)=k+ipart
!        icol(nz)=k_prev+ipart
!        amat_sp(nz)=(-1.d0,0.d0)
      ENDDO
!
      IF(isw_regper.EQ.1.AND.m.LE.1) THEN
        DO ipart=1,npassing+1
!          if(ipart.le.npassing) then
!            deleta_factor=(eta(ipart)-eta(ipart-1))*bhat_mfl(istep)
!          else
!            deleta_factor=1.d0-eta(ipart-1)*bhat_mfl(istep)
!          endif
!
          DO ipart1=1,npassing+1
            nz_regper=nz_regper+1
!            irow_regper(nz_regper)=k+ipart
!            icol_regper(nz_regper)=k_prev+ipart1
!            amat_regper(nz_regper)=deleta_factor
            nz_regper=nz_regper+1
!            irow_regper(nz_regper)=k+ipart
!            icol_regper(nz_regper)=k+npassing+1+ipart1
!            amat_regper(nz_regper)=deleta_factor
          ENDDO
!
        ENDDO
      ENDIF
!
    ENDIF
!
  ENDDO
!
  DO istep=ibeg+1,iend
    npassing_prev=npl(istep-1)
    npassing=npl(istep)
!    delphim1=1.d0/delt_pos(istep)
!    deloneovb=0.5d0*delphim1*(1.d0/bhat_mfl(istep-1)-1.d0/bhat_mfl(istep))
!
    DO m=0,lag
      k_prev=ind_start(istep-1)+2*(npassing_prev+1)*m
      k=ind_start(istep)+2*(npassing+1)*m
!
! free flight:
!
      DO ipart=1,npassing+1
        nz=nz+1
!        irow(nz)=k+ipart
!        icol(nz)=k+ipart
!        amat_sp(nz)=delphim1
        IF(colltest) nz_ttmp=nz_ttmp+1
      ENDDO
!
      DO ipart=1,npassing
        nz=nz+1
!        irow(nz)=k+ipart
!        icol(nz)=k_prev+ipart
!        amat_sp(nz)=-delphim1
        IF(colltest) nz_ttmp=nz_ttmp+1
      ENDDO
!
      IF(npassing_prev.GE.npassing) THEN
        nz=nz+1
!        irow(nz)=k+npassing+1
!        icol(nz)=k_prev+npassing+1
!        amat_sp(nz)=-delphim1
        IF(colltest) nz_ttmp=nz_ttmp+1
      ENDIF
!
! mirroring:
!
      IF(npassing_prev.EQ.npassing) THEN
!
        DO kk=1,4
          nz=nz+1
!          irow(nz)=k+npassing+1
!          icol(nz)=k+npassing+kk-1
!          amat_sp(nz)=deloneovb*rhs_mat_fzero(kk,istep,0)
          IF(colltest) nz_ttmp=nz_ttmp+1
        ENDDO
!
        DO kk=1,4
          nz=nz+1
!          irow(nz)=k+npassing+1
!          icol(nz)=k_prev+npassing_prev+kk-1
!          amat_sp(nz)=deloneovb*rhs_mat_fzero(kk,istep-1,0)
          IF(colltest) nz_ttmp=nz_ttmp+1
        ENDDO
!
      ELSEIF(npassing_prev.GT.npassing) THEN
        nz=nz+1
!        irow(nz)=k_prev+npassing_prev+2
!        icol(nz)=k_prev+npassing_prev+1
!        amat_sp(nz)=-delphim1
         IF(colltest) nz_ttmp=nz_ttmp+1
      ENDIF
!
! collisions:
!
      IF(fact_pos_e(istep).NE.0.d0) THEN
!
! Lorentz operator:
!
        IF(isw_lor.EQ.1) THEN
!
          DO ipart=1,npassing+1
            DO kk=1,5
              DO mm=0,lag
                nz=nz+1
                nz_coll=nz_coll+1
!                irow(nz)=k+ipart
!                icol(nz)=k+max(0,ipart-3)+kk+2*(npassing+1)*(mm-m)
!                amat_sp(nz)=anumm(m,mm)*rhs_mat_lorentz(kk,ipart,istep) &
!                           *fact_pos_e(istep)
                IF(.NOT.colltest.AND.mm.EQ.m) THEN
                  nz_ttmp=nz_ttmp+1
                ENDIF
              ENDDO
            ENDDO
          ENDDO
!
! matching collisional flux through the boundary in backward Euler scheme:
!
          DO kk=1,4
            DO mm=0,lag
              nz=nz+1
              nz_coll=nz_coll+1
!              irow(nz)=k+npassing+1
!              icol(nz)=k+npassing-1+kk+2*(npassing+1)*(mm-m)
!              amat_sp(nz)=-0.5d0*anumm(m,mm)*rhs_mat_fzero(kk,istep,1)  &
!                         *fact_pos_e(istep)
            ENDDO
          ENDDO
!
          DO kk=1,4
            DO mm=0,lag
              nz=nz+1
              nz_coll=nz_coll+1
!              irow(nz)=k+npassing+1
!              icol(nz)=k_prev+npassing_prev-1+kk+2*(npassing_prev+1)*(mm-m)
!              amat_sp(nz)=0.5d0*anumm(m,mm)*rhs_mat_fzero(kk,istep-1,1)
!                         *fact_pos_b(istep-1)
            ENDDO
          ENDDO
!
        ENDIF
!
!        nz_beg=nz+1
!
! energy diffusion operator:
!
        IF(isw_ene.EQ.1) THEN
!
          DO ipart=1,npassing
            DO kk=1,4
              DO mm=0,lag
                nz=nz+1
                nz_coll=nz_coll+1
!                irow(nz)=k+ipart
!                icol(nz)=k+max(0,ipart-2)+kk+2*(npassing+1)*(mm-m)
!                amat_sp(nz)=denmm(m,mm)*6hs_mat_energ(kk,ipart,istep)
              ENDDO
            ENDDO
          ENDDO
!
          DO kk=1,4
            DO mm=0,lag
              nz=nz+1
              nz_coll=nz_coll+1
!              irow(nz)=k+npassing+1
!              icol(nz)=k+npassing-1+kk+2*(npassing+1)*(mm-m)
!              amat_sp(nz)=0.5d0*denmm(m,mm)                                  &
!                         *rhs_mat_energ(kk,npassing+1,istep)
            ENDDO
          ENDDO
!
          DO kk=1,4
            DO mm=0,lag
              nz=nz+1
              nz_coll=nz_coll+1
!              irow(nz)=k+npassing+1
!              icol(nz)=k_prev+npassing_prev-1+kk+2*(npassing_prev+1)*(mm-m)
!              amat_sp(nz)=0.5d0*denmm(m,mm)                                  &
!                         *rhs_mat_energ(kk,npassing_prev+1,istep-1)
            ENDDO
          ENDDO
!
        ENDIF
!
!        amat_sp(nz_beg:nz)=fact_pos_e(istep)*amat_sp(nz_beg:nz)
!
!        nz_beg=nz+1
!
! regularization sink:
!
        DO ipart=1,npassing
          DO kk=1,4
            DO mm=0,lag
              nz=nz+1
!              irow(nz)=k+ipart
!              icol(nz)=k+max(0,ipart-2)+kk+2*(npassing+1)*(mm-m)
!              amat_sp(nz)=x1mm(m,mm)*rhs_mat_energ(kk,ipart,istep)
            ENDDO
          ENDDO
        ENDDO
!
        DO kk=1,4
          DO mm=0,lag
            nz=nz+1
!            irow(nz)=k+npassing+1
!            icol(nz)=k+npassing-1+kk+2*(npassing+1)*(mm-m)
!            amat_sp(nz)=0.5d0*x1mm(m,mm)                                  &
!                       *rhs_mat_energ(kk,npassing+1,istep)
          ENDDO
        ENDDO
!
        DO kk=1,4
          DO mm=0,lag
            nz=nz+1
!            irow(nz)=k+npassing+1
!            icol(nz)=k_prev+npassing_prev-1+kk+2*(npassing_prev+1)*(mm-m)
!            amat_sp(nz)=0.5d0*x1mm(m,mm)                                  &
!                       *rhs_mat_energ(kk,npassing_prev+1,istep-1)
          ENDDO
        ENDDO
!
!        amat_sp(nz_beg:nz)=-epserr_sink*fact_pos_e(istep)*amat_sp(nz_beg:nz)
!
      ENDIF
!
    ENDDO
!
  ENDDO
!
! Counter-passing: sigma=-1
!
  istep=iend
  npassing=npl(istep)
!
! entry:
!
  DO m=0,lag
    k=ind_start(istep)+2*(npassing+1)*m+2*npassing+3
!
    DO ipart=1,npassing+1
      nz=nz+1
!      irow(nz)=k-ipart
!      icol(nz)=k-ipart
!      amat_sp(nz)=1.d0
    ENDDO
!
    IF(isw_axisymm.EQ.1) THEN
! periodicity:
      k_prev=ind_start(ibeg)+2*(npassing+1)*m+2*npassing+3
!
      DO ipart=1,npassing+1
        nz=nz+1
!        irow(nz)=k-ipart
!        icol(nz)=k_prev-ipart
!        amat_sp(nz)=(-1.d0,0.d0)
      ENDDO
!
      IF(isw_regper.EQ.1.AND.m.LE.1) THEN
        DO ipart=1,npassing+1
!          if(ipart.le.npassing) then
!            deleta_factor=(eta(ipart)-eta(ipart-1))*bhat_mfl(istep)
!          else
!            deleta_factor=1.d0-eta(ipart-1)*bhat_mfl(istep)
!          endif
!
          DO ipart1=1,npassing+1
            nz_regper=nz_regper+1
!            irow_regper(nz_regper)=k-ipart
!            icol_regper(nz_regper)=k_prev-ipart1
!            amat_regper(nz_regper)=deleta_factor
            nz_regper=nz_regper+1
!            irow_regper(nz_regper)=k-ipart
!            icol_regper(nz_regper)=k-npassing-1-ipart1
!            amat_regper(nz_regper)=deleta_factor
          ENDDO
!
        ENDDO
      ENDIF
!
    ENDIF
!
  ENDDO
!
  DO istep=ibeg,iend-1
    npassing_prev=npl(istep+1)
    npassing=npl(istep)
!    delphim1=1.d0/delt_neg(istep)
!    deloneovb=0.5d0*delphim1*(1.d0/bhat_mfl(istep+1)-1.d0/bhat_mfl(istep))
!
    DO m=0,lag
      k_prev=ind_start(istep+1)+2*(npassing_prev+1)*m+2*npassing_prev+3
      k=ind_start(istep)+2*(npassing+1)*m+2*npassing+3
!
! free flight:
!
      DO ipart=1,npassing+1
        nz=nz+1
!        irow(nz)=k-ipart
!        icol(nz)=k-ipart
!        amat_sp(nz)=delphim1
        IF(colltest) nz_ttmp=nz_ttmp+1
      ENDDO
!
      DO ipart=1,npassing
        nz=nz+1
!        irow(nz)=k-ipart
!        icol(nz)=k_prev-ipart
!        amat_sp(nz)=-delphim1
        IF(colltest) nz_ttmp=nz_ttmp+1
      ENDDO
!
      IF(npassing_prev.GE.npassing) THEN
        nz=nz+1
!        irow(nz)=k-npassing-1
!        icol(nz)=k_prev-npassing-1
!        amat_sp(nz)=-delphim1
        IF(colltest) nz_ttmp=nz_ttmp+1
      ENDIF
!
! mirroring:
!
      IF(npassing_prev.EQ.npassing) THEN
!
        DO kk=1,4
          nz=nz+1
!          irow(nz)=k-npassing-1
!          icol(nz)=k-npassing-kk+1
!          amat_sp(nz)=deloneovb*rhs_mat_fzero(kk,istep,0)
          IF(colltest) nz_ttmp=nz_ttmp+1
        ENDDO
!
        DO kk=1,4
          nz=nz+1
!          irow(nz)=k-npassing-1
!          icol(nz)=k_prev-npassing_prev-kk+1
!          amat_sp(nz)=deloneovb*rhs_mat_fzero(kk,istep+1,0)
          IF(colltest) nz_ttmp=nz_ttmp+1
        ENDDO
!
      ELSEIF(npassing_prev.GT.npassing) THEN
        nz=nz+1
!        irow(nz)=k_prev-npassing_prev-2
!        icol(nz)=k_prev-npassing_prev-1
!        amat_sp(nz)=-delphim1
        IF(colltest) nz_ttmp=nz_ttmp+1
      ENDIF
!
! collisions:
!
      IF(fact_neg_e(istep).NE.0.d0) THEN
!
! Lorentz operator:
!
        IF(isw_lor.EQ.1) THEN
!
          DO ipart=1,npassing+1
            DO kk=1,5
              DO mm=0,lag
                nz=nz+1
                nz_coll=nz_coll+1
!                irow(nz)=k-ipart
!                icol(nz)=k-max(0,ipart-3)-kk+2*(npassing+1)*(mm-m)
!                amat_sp(nz)=anumm(m,mm)*rhs_mat_lorentz(kk,ipart,istep) &
!                           *fact_neg_e(istep)
                IF(.NOT.colltest.AND.mm.EQ.m) THEN
                  nz_ttmp=nz_ttmp+1
                ENDIF
              ENDDO
            ENDDO
          ENDDO
!
! matching collisional flux through the boundary in backward Euler scheme:
!
          DO kk=1,4
            DO mm=0,lag
              nz=nz+1
              nz_coll=nz_coll+1
!              irow(nz)=k-npassing-1
!              icol(nz)=k-npassing+1-kk+2*(npassing+1)*(mm-m)
!              amat_sp(nz)=-0.5d0*anumm(m,mm)*rhs_mat_fzero(kk,istep,1) &
!                         *fact_neg_e(istep)
            ENDDO
          ENDDO
!
          DO kk=1,4
            DO mm=0,lag
              nz=nz+1
              nz_coll=nz_coll+1
!              irow(nz)=k-npassing-1
!              icol(nz)=k_prev-npassing_prev+1-kk+2*(npassing_prev+1)*(mm-m)
!              amat_sp(nz)=0.5d0*anumm(m,mm)*rhs_mat_fzero(kk,istep+1,1) &
!                         *fact_neg_b(istep+1)
            ENDDO
          ENDDO
!
        ENDIF
!
!        nz_beg=nz+1
!
! energy diffusion operator:
!
        IF(isw_ene.EQ.1) THEN
!
          DO ipart=1,npassing
            DO kk=1,4
              DO mm=0,lag
                nz=nz+1
                nz_coll=nz_coll+1
!                irow(nz)=k-ipart
!                icol(nz)=k-max(0,ipart-2)-kk+2*(npassing+1)*(mm-m)
!                amat_sp(nz)=denmm(m,mm)*rhs_mat_energ(kk,ipart,istep)
              ENDDO
            ENDDO
          ENDDO
!
          DO kk=1,4
            DO mm=0,lag
              nz=nz+1
              nz_coll=nz_coll+1
!              irow(nz)=k-npassing-1
!              icol(nz)=k-npassing+1-kk+2*(npassing+1)*(mm-m)
!              amat_sp(nz)=0.5d0*denmm(m,mm)                                  &
!                         *rhs_mat_energ(kk,npassing+1,istep)
            ENDDO
          ENDDO
!
          DO kk=1,4
            DO mm=0,lag
              nz=nz+1
              nz_coll=nz_coll+1
!              irow(nz)=k-npassing-1
!              icol(nz)=k_prev-npassing_prev+1-kk+2*(npassing_prev+1)*(mm-m)
!              amat_sp(nz)=0.5d0*denmm(m,mm)                                  &
!                         *rhs_mat_energ(kk,npassing_prev+1,istep+1)
            ENDDO
          ENDDO
!
        ENDIF
!
!        amat_sp(nz_beg:nz)=fact_neg_e(istep)*amat_sp(nz_beg:nz)
!
!        nz_beg=nz+1
!
! regularization sink:
!
        DO ipart=1,npassing
          DO kk=1,4
            DO mm=0,lag
              nz=nz+1
!              irow(nz)=k-ipart
!              icol(nz)=k-max(0,ipart-2)-kk+2*(npassing+1)*(mm-m)
!              amat_sp(nz)=x1mm(m,mm)*rhs_mat_energ(kk,ipart,istep)
            ENDDO
          ENDDO
        ENDDO
!
        DO kk=1,4
          DO mm=0,lag
            nz=nz+1
!            irow(nz)=k-npassing-1
!            icol(nz)=k-npassing+1-kk+2*(npassing+1)*(mm-m)
!            amat_sp(nz)=0.5d0*x1mm(m,mm)                                  &
!                       *rhs_mat_energ(kk,npassing+1,istep)
          ENDDO
        ENDDO
!
        DO kk=1,4
          DO mm=0,lag
            nz=nz+1
!            irow(nz)=k-npassing-1
!            icol(nz)=k_prev-npassing_prev+1-kk+2*(npassing_prev+1)*(mm-m)
!            amat_sp(nz)=0.5d0*x1mm(m,mm)                                  &
!                       *rhs_mat_energ(kk,npassing_prev+1,istep+1)
          ENDDO
        ENDDO
!
!        amat_sp(nz_beg:nz)=-epserr_sink*fact_neg_e(istep)*amat_sp(nz_beg:nz)
!
      ENDIF
!
    ENDDO
!
  ENDDO
!
!
  ALLOCATE(irow(nz),icol(nz),amat_sp(nz))
  ALLOCATE(irow_coll(nz_coll),icol_coll(nz_coll),amat_coll(nz_coll))
  ALLOCATE(irow_ttmp(nz_ttmp),icol_ttmp(nz_ttmp),amat_ttmp(nz_ttmp))
  ALLOCATE(irow_regper(nz_regper),icol_regper(nz_regper),amat_regper(nz_regper))
!
! Fill the arrays:
!
  nz=0
  nz_coll=0
  nz_ttmp=0
  nz_regper=0
!
! Co-passing: sigma=1
!
  istep=ibeg
  npassing=npl(istep)
!
! entry:
!
  DO m=0,lag
    k=ind_start(istep)+2*(npassing+1)*m
!
    DO ipart=1,npassing+1
      nz=nz+1
      irow(nz)=k+ipart
      icol(nz)=k+ipart
      amat_sp(nz)=1.d0
    ENDDO
!
    IF(isw_axisymm.EQ.1) THEN
! periodicity:
      k_prev=ind_start(iend)+2*(npassing+1)*m
!
      DO ipart=1,npassing+1
        nz=nz+1
        irow(nz)=k+ipart
        icol(nz)=k_prev+ipart
        amat_sp(nz)=(-1.d0,0.d0)
      ENDDO
!
      IF(isw_regper.EQ.1.AND.m.LE.1) THEN
        DO ipart=1,npassing+1
          IF(ipart.LE.npassing) THEN
            deleta_factor=(eta(ipart)-eta(ipart-1))*bhat_mfl(istep)
          ELSE
            deleta_factor=1.d0-eta(ipart-1)*bhat_mfl(istep)
          ENDIF
!
          DO ipart1=1,npassing+1
            nz_regper=nz_regper+1
            irow_regper(nz_regper)=k+ipart
            icol_regper(nz_regper)=k_prev+ipart1
            amat_regper(nz_regper)=deleta_factor
            nz_regper=nz_regper+1
            irow_regper(nz_regper)=k+ipart
            icol_regper(nz_regper)=k+npassing+1+ipart1
            amat_regper(nz_regper)=deleta_factor
          ENDDO
!
        ENDDO
      ENDIF
!
    ENDIF
!
  ENDDO
!
  DO istep=ibeg+1,iend
    npassing_prev=npl(istep-1)
    npassing=npl(istep)
    delphim1=1.d0/delt_pos(istep)
    deloneovb=0.5d0*delphim1*(1.d0/bhat_mfl(istep-1)-1.d0/bhat_mfl(istep))
!
    DO m=0,lag
      k_prev=ind_start(istep-1)+2*(npassing_prev+1)*m
      k=ind_start(istep)+2*(npassing+1)*m
!
! free flight:
!
      DO ipart=1,npassing+1
        nz=nz+1
        irow(nz)=k+ipart
        icol(nz)=k+ipart
        amat_sp(nz)=delphim1
        IF(colltest) THEN
          nz_ttmp=nz_ttmp+1
          irow_ttmp(nz_ttmp)=irow(nz)
          icol_ttmp(nz_ttmp)=icol(nz)
          !! Modification by Andreas F. Martitsch (17.07.2015)
          ! fixed warning: Possible change of value in conversion
          ! from COMPLEX(8) to REAL(8)
          amat_ttmp(nz_ttmp)=REAL(amat_sp(nz))
          !! End Modification by Andreas F. Martitsch (17.07.2015)                
        ENDIF
      ENDDO
!
      DO ipart=1,npassing
        nz=nz+1
        irow(nz)=k+ipart
        icol(nz)=k_prev+ipart
        amat_sp(nz)=-delphim1
        IF(colltest) THEN
          nz_ttmp=nz_ttmp+1
          irow_ttmp(nz_ttmp)=irow(nz)
          icol_ttmp(nz_ttmp)=icol(nz)
          !! Modification by Andreas F. Martitsch (17.07.2015)
          ! fixed warning: Possible change of value in conversion
          ! from COMPLEX(8) to REAL(8)
          amat_ttmp(nz_ttmp)=REAL(amat_sp(nz))
          !! End Modification by Andreas F. Martitsch (17.07.2015)            
        ENDIF
      ENDDO
!
      IF(npassing_prev.GE.npassing) THEN
        nz=nz+1
        irow(nz)=k+npassing+1
        icol(nz)=k_prev+npassing+1
        amat_sp(nz)=-delphim1
        IF(colltest) THEN
          nz_ttmp=nz_ttmp+1
          irow_ttmp(nz_ttmp)=irow(nz)
          icol_ttmp(nz_ttmp)=icol(nz)
          !! Modification by Andreas F. Martitsch (17.07.2015)
          ! fixed warning: Possible change of value in conversion
          ! from COMPLEX(8) to REAL(8)
          amat_ttmp(nz_ttmp)=REAL(amat_sp(nz))
          !! End Modification by Andreas F. Martitsch (17.07.2015)            
        ENDIF
      ENDIF
!
! mirroring:
!
      IF(npassing_prev.EQ.npassing) THEN
!
        DO kk=1,4
          nz=nz+1
          irow(nz)=k+npassing+1
          icol(nz)=k+npassing+kk-1
          amat_sp(nz)=deloneovb*rhs_mat_fzero(kk,istep,0)
          IF(colltest) THEN
            nz_ttmp=nz_ttmp+1
            irow_ttmp(nz_ttmp)=irow(nz)
            icol_ttmp(nz_ttmp)=icol(nz)
            !! Modification by Andreas F. Martitsch (17.07.2015)
            ! fixed warning: Possible change of value in conversion
            ! from COMPLEX(8) to REAL(8)
            amat_ttmp(nz_ttmp)=REAL(amat_sp(nz))
            !! End Modification by Andreas F. Martitsch (17.07.2015)            
          ENDIF
        ENDDO
!
        DO kk=1,4
          nz=nz+1
          irow(nz)=k+npassing+1
          icol(nz)=k_prev+npassing_prev+kk-1
          amat_sp(nz)=deloneovb*rhs_mat_fzero(kk,istep-1,0)
          IF(colltest) THEN
            nz_ttmp=nz_ttmp+1
            irow_ttmp(nz_ttmp)=irow(nz)
            icol_ttmp(nz_ttmp)=icol(nz)
            !! Modification by Andreas F. Martitsch (17.07.2015)
            ! fixed warning: Possible change of value in conversion
            ! from COMPLEX(8) to REAL(8)
            amat_ttmp(nz_ttmp)=REAL(amat_sp(nz))
            !! End Modification by Andreas F. Martitsch (17.07.2015)
          ENDIF
        ENDDO
!
      ELSEIF(npassing_prev.GT.npassing) THEN
        nz=nz+1
        irow(nz)=k_prev+npassing_prev+2
        icol(nz)=k_prev+npassing_prev+1
        amat_sp(nz)=-delphim1
        IF(colltest) THEN
          nz_ttmp=nz_ttmp+1
          irow_ttmp(nz_ttmp)=irow(nz)
          icol_ttmp(nz_ttmp)=icol(nz)
          !! Modification by Andreas F. Martitsch (17.07.2015)
          ! fixed warning: Possible change of value in conversion
          ! from COMPLEX(8) to REAL(8)
          amat_ttmp(nz_ttmp)=REAL(amat_sp(nz))
          !! End Modification by Andreas F. Martitsch (17.07.2015)
        ENDIF
      ENDIF
!
! collisions:
!
      IF(fact_pos_e(istep).NE.0.d0) THEN
!
! Lorentz operator:
!
        IF(isw_lor.EQ.1) THEN
!
          DO ipart=1,npassing+1
            DO kk=1,5
              DO mm=0,lag
                nz=nz+1
                irow(nz)=k+ipart
                icol(nz)=k+MAX(0,ipart-3)+kk+2*(npassing+1)*(mm-m)
                amat_sp(nz)=anumm(m,mm)*rhs_mat_lorentz(kk,ipart,istep) &
                           *fact_pos_e(istep)
                nz_coll=nz_coll+1
                irow_coll(nz_coll)=irow(nz)
                icol_coll(nz_coll)=icol(nz)
                amat_coll(nz_coll)=anumm(m,mm)*rhs_mat_lorentz(kk,ipart,istep)
                IF(.NOT.colltest.AND.mm.EQ.m) THEN
                  nz_ttmp=nz_ttmp+1
                  irow_ttmp(nz_ttmp)=irow(nz)
                  icol_ttmp(nz_ttmp)=icol(nz)
                  amat_ttmp(nz_ttmp)=-ttmp_mat(kk,ipart,istep)
                  IF(irow(nz).EQ.icol(nz)) THEN
                    amat_ttmp(nz_ttmp)=amat_ttmp(nz_ttmp)+1.d0
                  ENDIF
                ENDIF
              ENDDO
            ENDDO
          ENDDO
!
! matching collisional flux through the boundary in backward Euler scheme:
!
          DO kk=1,4
            DO mm=0,lag
              nz=nz+1
              irow(nz)=k+npassing+1
              icol(nz)=k+npassing-1+kk+2*(npassing+1)*(mm-m)
              amat_sp(nz)=-0.5d0*anumm(m,mm)*rhs_mat_fzero(kk,istep,1)  &
                         *fact_pos_e(istep)
              nz_coll=nz_coll+1
              irow_coll(nz_coll)=irow(nz)
              icol_coll(nz_coll)=icol(nz)
              amat_coll(nz_coll)=-0.5d0*anumm(m,mm)*rhs_mat_fzero(kk,istep,1)
            ENDDO
          ENDDO
!
          DO kk=1,4
            DO mm=0,lag
              nz=nz+1
              irow(nz)=k+npassing+1
              icol(nz)=k_prev+npassing_prev-1+kk+2*(npassing_prev+1)*(mm-m)
              amat_sp(nz)=0.5d0*anumm(m,mm)*rhs_mat_fzero(kk,istep-1,1) &
                         *fact_pos_b(istep-1)
              nz_coll=nz_coll+1
              irow_coll(nz_coll)=irow(nz)
              icol_coll(nz_coll)=icol(nz)
              amat_coll(nz_coll)=0.5d0*anumm(m,mm)*rhs_mat_fzero(kk,istep-1,1)
            ENDDO
          ENDDO
!
        ENDIF
!
        nz_beg=nz+1
        nz_coll_beg=nz_coll+1
!
! energy diffusion operator:
!
        IF(isw_ene.EQ.1) THEN
!
          DO ipart=1,npassing
            DO kk=1,4
              DO mm=0,lag
                nz=nz+1
                irow(nz)=k+ipart
                icol(nz)=k+MAX(0,ipart-2)+kk+2*(npassing+1)*(mm-m)
                amat_sp(nz)=denmm(m,mm)*rhs_mat_energ(kk,ipart,istep)
                nz_coll=nz_coll+1
              ENDDO
            ENDDO
          ENDDO
!
          DO kk=1,4
            DO mm=0,lag
              nz=nz+1
              irow(nz)=k+npassing+1
              icol(nz)=k+npassing-1+kk+2*(npassing+1)*(mm-m)
              amat_sp(nz)=0.5d0*denmm(m,mm)                                  &
                         *rhs_mat_energ(kk,npassing+1,istep)
              nz_coll=nz_coll+1
            ENDDO
          ENDDO
!
          DO kk=1,4
            DO mm=0,lag
              nz=nz+1
              irow(nz)=k+npassing+1
              icol(nz)=k_prev+npassing_prev-1+kk+2*(npassing_prev+1)*(mm-m)
              amat_sp(nz)=0.5d0*denmm(m,mm)                                  &
                         *rhs_mat_energ(kk,npassing_prev+1,istep-1)
              nz_coll=nz_coll+1
            ENDDO
          ENDDO
!
        ENDIF
!
        irow_coll(nz_coll_beg:nz_coll)=irow(nz_beg:nz)
        icol_coll(nz_coll_beg:nz_coll)=icol(nz_beg:nz)
        !! Modification by Andreas F. Martitsch (17.07.2015)
        ! fixed warning: Possible change of value in conversion
        ! from COMPLEX(8) to REAL(8)
        amat_coll(nz_coll_beg:nz_coll)=REAL(amat_sp(nz_beg:nz))
        !! End Modification by Andreas F. Martitsch (17.07.2015)
!
        amat_sp(nz_beg:nz)=fact_pos_e(istep)*amat_sp(nz_beg:nz)
!
        nz_beg=nz+1
!
! regularization sink:
!
        DO ipart=1,npassing
          DO kk=1,4
            DO mm=0,lag
              nz=nz+1
              irow(nz)=k+ipart
              icol(nz)=k+MAX(0,ipart-2)+kk+2*(npassing+1)*(mm-m)
              amat_sp(nz)=x1mm(m,mm)*rhs_mat_energ(kk,ipart,istep)
            ENDDO
          ENDDO
        ENDDO
!
        DO kk=1,4
          DO mm=0,lag
            nz=nz+1
            irow(nz)=k+npassing+1
            icol(nz)=k+npassing-1+kk+2*(npassing+1)*(mm-m)
            amat_sp(nz)=0.5d0*x1mm(m,mm)                                  &
                       *rhs_mat_energ(kk,npassing+1,istep)
          ENDDO
        ENDDO
!
        DO kk=1,4
          DO mm=0,lag
            nz=nz+1
            irow(nz)=k+npassing+1
            icol(nz)=k_prev+npassing_prev-1+kk+2*(npassing_prev+1)*(mm-m)
            amat_sp(nz)=0.5d0*x1mm(m,mm)                                  &
                       *rhs_mat_energ(kk,npassing_prev+1,istep-1)
          ENDDO
        ENDDO
!
        amat_sp(nz_beg:nz)=-epserr_sink*fact_pos_e(istep)*amat_sp(nz_beg:nz)
!
      ENDIF
!
    ENDDO
!
  ENDDO
!
! Counter-passing: sigma=-1
!
  istep=iend
  npassing=npl(istep)
!
! entry:
!
  DO m=0,lag
    k=ind_start(istep)+2*(npassing+1)*m+2*npassing+3
!
    DO ipart=1,npassing+1
      nz=nz+1
      irow(nz)=k-ipart
      icol(nz)=k-ipart
      amat_sp(nz)=1.d0
    ENDDO
!
    IF(isw_axisymm.EQ.1) THEN
! periodicity:
      k_prev=ind_start(ibeg)+2*(npassing+1)*m+2*npassing+3
!
      DO ipart=1,npassing+1
        nz=nz+1
        irow(nz)=k-ipart
        icol(nz)=k_prev-ipart
        amat_sp(nz)=(-1.d0,0.d0)
      ENDDO
!
      IF(isw_regper.EQ.1.AND.m.LE.1) THEN
        DO ipart=1,npassing+1
          IF(ipart.LE.npassing) THEN
            deleta_factor=(eta(ipart)-eta(ipart-1))*bhat_mfl(istep)
          ELSE
            deleta_factor=1.d0-eta(ipart-1)*bhat_mfl(istep)
          ENDIF
!
          DO ipart1=1,npassing+1
            nz_regper=nz_regper+1
            irow_regper(nz_regper)=k-ipart
            icol_regper(nz_regper)=k_prev-ipart1
            amat_regper(nz_regper)=deleta_factor
            nz_regper=nz_regper+1
            irow_regper(nz_regper)=k-ipart
            icol_regper(nz_regper)=k-npassing-1-ipart1
            amat_regper(nz_regper)=deleta_factor
          ENDDO
!
        ENDDO
      ENDIF
!
    ENDIF
!
  ENDDO
!
  DO istep=ibeg,iend-1
    npassing_prev=npl(istep+1)
    npassing=npl(istep)
    delphim1=1.d0/delt_neg(istep)
    deloneovb=0.5d0*delphim1*(1.d0/bhat_mfl(istep+1)-1.d0/bhat_mfl(istep))
!
    DO m=0,lag
      k_prev=ind_start(istep+1)+2*(npassing_prev+1)*m+2*npassing_prev+3
      k=ind_start(istep)+2*(npassing+1)*m+2*npassing+3
!
! free flight:
!
      DO ipart=1,npassing+1
        nz=nz+1
        irow(nz)=k-ipart
        icol(nz)=k-ipart
        amat_sp(nz)=delphim1
        IF(colltest) THEN
          nz_ttmp=nz_ttmp+1
          irow_ttmp(nz_ttmp)=irow(nz)
          icol_ttmp(nz_ttmp)=icol(nz)
          !! Modification by Andreas F. Martitsch (17.07.2015)
          ! fixed warning: Possible change of value in conversion
          ! from COMPLEX(8) to REAL(8)
          amat_ttmp(nz_ttmp)=REAL(amat_sp(nz))
          !! End Modification by Andreas F. Martitsch (17.07.2015)
        ENDIF
      ENDDO
!
      DO ipart=1,npassing
        nz=nz+1
        irow(nz)=k-ipart
        icol(nz)=k_prev-ipart
        amat_sp(nz)=-delphim1
        IF(colltest) THEN
          nz_ttmp=nz_ttmp+1
          irow_ttmp(nz_ttmp)=irow(nz)
          icol_ttmp(nz_ttmp)=icol(nz)
          !! Modification by Andreas F. Martitsch (17.07.2015)
          ! fixed warning: Possible change of value in conversion
          ! from COMPLEX(8) to REAL(8)
          amat_ttmp(nz_ttmp)=REAL(amat_sp(nz))
          !! End Modification by Andreas F. Martitsch (17.07.2015)
        ENDIF
      ENDDO
!
      IF(npassing_prev.GE.npassing) THEN
        nz=nz+1
        irow(nz)=k-npassing-1
        icol(nz)=k_prev-npassing-1
        amat_sp(nz)=-delphim1
        IF(colltest) THEN
          nz_ttmp=nz_ttmp+1
          irow_ttmp(nz_ttmp)=irow(nz)
          icol_ttmp(nz_ttmp)=icol(nz)
          !! Modification by Andreas F. Martitsch (17.07.2015)
          ! fixed warning: Possible change of value in conversion
          ! from COMPLEX(8) to REAL(8)
          amat_ttmp(nz_ttmp)=REAL(amat_sp(nz))
          !! End Modification by Andreas F. Martitsch (17.07.2015)
        ENDIF
      ENDIF
!
! mirroring:
!
      IF(npassing_prev.EQ.npassing) THEN
!
        DO kk=1,4
          nz=nz+1
          irow(nz)=k-npassing-1
          icol(nz)=k-npassing-kk+1
          amat_sp(nz)=deloneovb*rhs_mat_fzero(kk,istep,0)
          IF(colltest) THEN
            nz_ttmp=nz_ttmp+1
            irow_ttmp(nz_ttmp)=irow(nz)
            icol_ttmp(nz_ttmp)=icol(nz)
            !! Modification by Andreas F. Martitsch (17.07.2015)
            ! fixed warning: Possible change of value in conversion
            ! from COMPLEX(8) to REAL(8)
            amat_ttmp(nz_ttmp)=REAL(amat_sp(nz))
            !! End Modification by Andreas F. Martitsch (17.07.2015)
          ENDIF
        ENDDO
!
        DO kk=1,4
          nz=nz+1
          irow(nz)=k-npassing-1
          icol(nz)=k_prev-npassing_prev-kk+1
          amat_sp(nz)=deloneovb*rhs_mat_fzero(kk,istep+1,0)
          IF(colltest) THEN
            nz_ttmp=nz_ttmp+1
            irow_ttmp(nz_ttmp)=irow(nz)
            icol_ttmp(nz_ttmp)=icol(nz)
            !! Modification by Andreas F. Martitsch (17.07.2015)
            ! fixed warning: Possible change of value in conversion
            ! from COMPLEX(8) to REAL(8)
            amat_ttmp(nz_ttmp)=REAL(amat_sp(nz))
            !! End Modification by Andreas F. Martitsch (17.07.2015)
          ENDIF
        ENDDO
!
      ELSEIF(npassing_prev.GT.npassing) THEN
        nz=nz+1
        irow(nz)=k_prev-npassing_prev-2
        icol(nz)=k_prev-npassing_prev-1
        amat_sp(nz)=-delphim1
        IF(colltest) THEN
          nz_ttmp=nz_ttmp+1
          irow_ttmp(nz_ttmp)=irow(nz)
          icol_ttmp(nz_ttmp)=icol(nz)
          !! Modification by Andreas F. Martitsch (17.07.2015)
          ! fixed warning: Possible change of value in conversion
          ! from COMPLEX(8) to REAL(8)
          amat_ttmp(nz_ttmp)=REAL(amat_sp(nz))
          !! End Modification by Andreas F. Martitsch (17.07.2015)
        ENDIF
      ENDIF
!
! collisions:
!
      IF(fact_neg_e(istep).NE.0.d0) THEN
!
! Lorentz operator:
!
        IF(isw_lor.EQ.1) THEN
!
          DO ipart=1,npassing+1
            DO kk=1,5
              DO mm=0,lag
                nz=nz+1
                irow(nz)=k-ipart
                icol(nz)=k-MAX(0,ipart-3)-kk+2*(npassing+1)*(mm-m)
                amat_sp(nz)=anumm(m,mm)*rhs_mat_lorentz(kk,ipart,istep) &
                           *fact_neg_e(istep)
                nz_coll=nz_coll+1
                irow_coll(nz_coll)=irow(nz)
                icol_coll(nz_coll)=icol(nz)
                amat_coll(nz_coll)=anumm(m,mm)*rhs_mat_lorentz(kk,ipart,istep)
                IF(.NOT.colltest.AND.mm.EQ.m) THEN
                  nz_ttmp=nz_ttmp+1
                  irow_ttmp(nz_ttmp)=irow(nz)
                  icol_ttmp(nz_ttmp)=icol(nz)
                  amat_ttmp(nz_ttmp)=ttmp_mat(kk,ipart,istep)
                  IF(irow(nz).EQ.icol(nz)) THEN
                    amat_ttmp(nz_ttmp)=amat_ttmp(nz_ttmp)-1.d0
                  ENDIF
                ENDIF
              ENDDO
            ENDDO
          ENDDO
!
! matching collisional flux through the boundary in backward Euler scheme:
!
          DO kk=1,4
            DO mm=0,lag
              nz=nz+1
              irow(nz)=k-npassing-1
              icol(nz)=k-npassing+1-kk+2*(npassing+1)*(mm-m)
              amat_sp(nz)=-0.5d0*anumm(m,mm)*rhs_mat_fzero(kk,istep,1)  &
                         *fact_neg_e(istep)
              nz_coll=nz_coll+1
              irow_coll(nz_coll)=irow(nz)
              icol_coll(nz_coll)=icol(nz)
              amat_coll(nz_coll)=-0.5d0*anumm(m,mm)*rhs_mat_fzero(kk,istep,1)
            ENDDO
          ENDDO
!
          DO kk=1,4
            DO mm=0,lag
              nz=nz+1
              irow(nz)=k-npassing-1
              icol(nz)=k_prev-npassing_prev+1-kk+2*(npassing_prev+1)*(mm-m)
              amat_sp(nz)=0.5d0*anumm(m,mm)*rhs_mat_fzero(kk,istep+1,1) &
                         *fact_neg_b(istep+1)
              nz_coll=nz_coll+1
              irow_coll(nz_coll)=irow(nz)
              icol_coll(nz_coll)=icol(nz)
              amat_coll(nz_coll)=0.5d0*anumm(m,mm)*rhs_mat_fzero(kk,istep+1,1)
            ENDDO
          ENDDO
!
        ENDIF
!
        nz_beg=nz+1
        nz_coll_beg=nz_coll+1
!
! energy diffusion operator:
!
        IF(isw_ene.EQ.1) THEN
!
          DO ipart=1,npassing
            DO kk=1,4
              DO mm=0,lag
                nz=nz+1
                irow(nz)=k-ipart
                icol(nz)=k-MAX(0,ipart-2)-kk+2*(npassing+1)*(mm-m)
                amat_sp(nz)=denmm(m,mm)*rhs_mat_energ(kk,ipart,istep)
                nz_coll=nz_coll+1
              ENDDO
            ENDDO
          ENDDO
!
          DO kk=1,4
            DO mm=0,lag
              nz=nz+1
              irow(nz)=k-npassing-1
              icol(nz)=k-npassing+1-kk+2*(npassing+1)*(mm-m)
              amat_sp(nz)=0.5d0*denmm(m,mm)                                  &
                         *rhs_mat_energ(kk,npassing+1,istep)
              nz_coll=nz_coll+1
            ENDDO
          ENDDO
!
          DO kk=1,4
            DO mm=0,lag
              nz=nz+1
              irow(nz)=k-npassing-1
              icol(nz)=k_prev-npassing_prev+1-kk+2*(npassing_prev+1)*(mm-m)
              amat_sp(nz)=0.5d0*denmm(m,mm)                                  &
                         *rhs_mat_energ(kk,npassing_prev+1,istep+1)
              nz_coll=nz_coll+1
            ENDDO
          ENDDO
!
        ENDIF
!
        irow_coll(nz_coll_beg:nz_coll)=irow(nz_beg:nz)
        icol_coll(nz_coll_beg:nz_coll)=icol(nz_beg:nz)
        !! Modification by Andreas F. Martitsch (17.07.2015)
        ! fixed warning: Possible change of value in conversion
        ! from COMPLEX(8) to REAL(8)
        amat_coll(nz_coll_beg:nz_coll)=REAL(amat_sp(nz_beg:nz))
        !! End Modification by Andreas F. Martitsch (17.07.2015)
!
        amat_sp(nz_beg:nz)=fact_neg_e(istep)*amat_sp(nz_beg:nz)
!
        nz_beg=nz+1
!
! regularization sink:
!
        DO ipart=1,npassing
          DO kk=1,4
            DO mm=0,lag
              nz=nz+1
              irow(nz)=k-ipart
              icol(nz)=k-MAX(0,ipart-2)-kk+2*(npassing+1)*(mm-m)
              amat_sp(nz)=x1mm(m,mm)*rhs_mat_energ(kk,ipart,istep)
            ENDDO
          ENDDO
        ENDDO
!
        DO kk=1,4
          DO mm=0,lag
            nz=nz+1
            irow(nz)=k-npassing-1
            icol(nz)=k-npassing+1-kk+2*(npassing+1)*(mm-m)
            amat_sp(nz)=0.5d0*x1mm(m,mm)                                  &
                       *rhs_mat_energ(kk,npassing+1,istep)
          ENDDO
        ENDDO
!
        DO kk=1,4
          DO mm=0,lag
            nz=nz+1
            irow(nz)=k-npassing-1
            icol(nz)=k_prev-npassing_prev+1-kk+2*(npassing_prev+1)*(mm-m)
            amat_sp(nz)=0.5d0*x1mm(m,mm)                                  &
                       *rhs_mat_energ(kk,npassing_prev+1,istep+1)
          ENDDO
        ENDDO
!
        amat_sp(nz_beg:nz)=-epserr_sink*fact_neg_e(istep)*amat_sp(nz_beg:nz)
!
      ENDIF
!
    ENDDO
!
  ENDDO
!
! Save the symmetric matrix:
!
  nz_symm=nz
  ALLOCATE(irow_symm(nz_symm),icol_symm(nz_symm),amat_symm(nz_symm))
  irow_symm=irow
  icol_symm=icol
  !! Modification by Andreas F. Martitsch (17.07.2015)
  ! fixed warning: Possible change of value in conversion
  ! from COMPLEX(8) to REAL(8)
  amat_symm=REAL(amat_sp)
  !! End Modification by Andreas F. Martitsch (17.07.2015)
!
! End save symmetric matrix
!
  DEALLOCATE(irow,icol,amat_sp)
!
! Nox-axisymmetric matrices:
!
! Periodicity for Co-passing, sigma=1
!
! Determine the size of arrays:
!
  nz=nz_symm
!
  istep=ibeg
  npassing=npl(istep)
!
  DO m=0,lag
    k=ind_start(istep)+2*(npassing+1)*m
    k_prev=ind_start(iend)+2*(npassing+1)*m
!
    DO ipart=1,npassing+1
      nz=nz+1
!      irow_per_pos(nz)=k+ipart
!      icol_per_pos(nz)=k_prev+ipart
    ENDDO
!
  ENDDO
!
  nz_per_pos=nz
!
  ALLOCATE(irow_per_pos(nz_symm+1:nz_per_pos))
  ALLOCATE(icol_per_pos(nz_symm+1:nz_per_pos))
!
! Fill the arrays:
!
  nz=nz_symm
!
  istep=ibeg
  npassing=npl(istep)
!
  DO m=0,lag
    k=ind_start(istep)+2*(npassing+1)*m
    k_prev=ind_start(iend)+2*(npassing+1)*m
!
    DO ipart=1,npassing+1
      nz=nz+1
      irow_per_pos(nz)=k+ipart
      icol_per_pos(nz)=k_prev+ipart
    ENDDO
!
  ENDDO
!
! End periodicity for Co-passing, sigma=1
!
! Periodicity for Counter-passing, sigma=-1
!
! Determine the size of arrays:
!
  nz=nz_per_pos
!
  istep=iend
  npassing=npl(istep)
!
  DO m=0,lag
    k=ind_start(istep)+2*(npassing+1)*m+2*npassing+3
    k_prev=ind_start(ibeg)+2*(npassing+1)*m+2*npassing+3
!
    DO ipart=1,npassing+1
      nz=nz+1
!      irow_per_neg(nz)=k-ipart
!      icol_per_neg(nz)=k_prev-ipart
    ENDDO
!
  ENDDO
!
  nz_per_neg=nz
!
  ALLOCATE(irow_per_neg(nz_per_pos+1:nz_per_neg))
  ALLOCATE(icol_per_neg(nz_per_pos+1:nz_per_neg))
!
! Fill the arrays:
!
  nz=nz_per_pos
!
  istep=iend
  npassing=npl(istep)
!
  DO m=0,lag
    k=ind_start(istep)+2*(npassing+1)*m+2*npassing+3
    k_prev=ind_start(ibeg)+2*(npassing+1)*m+2*npassing+3
!
    DO ipart=1,npassing+1
      nz=nz+1
      irow_per_neg(nz)=k-ipart
      icol_per_neg(nz)=k_prev-ipart
    ENDDO
!
  ENDDO
!
! End periodicity for Counter-passing, sigma=-1
!
! Rotation matrix:
!
! Determine the size of arrays:
!
  nz=nz_per_neg
!
! Co-passing: sigma=1
!
  DO istep=ibeg+1,iend
    npassing_prev=npl(istep-1)
    npassing=npl(istep)
!
    DO m=0,lag
      k_prev=ind_start(istep-1)+2*(npassing_prev+1)*m
      k=ind_start(istep)+2*(npassing+1)*m
!
      IF(fact_pos_e(istep).NE.0.d0) THEN
!        nz_beg=nz+1
!
! Toroidal rotation:
!
        DO ipart=1,npassing
          DO kk=1,4
            DO mm=0,lag
              nz=nz+1
!              irow_asymm(nz)=k+ipart
!              icol_asymm(nz)=k+max(0,ipart-2)+kk+2*(npassing+1)*(mm-m)
!              amat_asymm(nz)=x1mm(m,mm)*rhs_mat_energ(kk,ipart,istep)
            ENDDO
          ENDDO
        ENDDO
!
        DO kk=1,4
          DO mm=0,lag
            nz=nz+1
!            irow_asymm(nz)=k+npassing+1
!            icol_asymm(nz)=k+npassing-1+kk+2*(npassing+1)*(mm-m)
!            amat_asymm(nz)=0.5d0*x1mm(m,mm)                                  &
!                       *rhs_mat_energ(kk,npassing+1,istep)
          ENDDO
        ENDDO
!
        DO kk=1,4
          DO mm=0,lag
            nz=nz+1
!            irow_asymm(nz)=k+npassing+1
!            icol_asymm(nz)=k_prev+npassing_prev-1+kk+2*(npassing_prev+1)*(mm-m)
!            amat_asymm(nz)=0.5d0*x1mm(m,mm)                                  &
!                       *rhs_mat_energ(kk,npassing_prev+1,istep-1)
          ENDDO
        ENDDO
!
!        amat_asymm(nz_beg:nz)=-fact_pos_e(istep)*amat_asymm(nz_beg:nz)
!
      ENDIF
!
    ENDDO
!
  ENDDO
!
! Counter-passing: sigma=-1
!
  DO istep=ibeg,iend-1
    npassing_prev=npl(istep+1)
    npassing=npl(istep)
!
    DO m=0,lag
      k_prev=ind_start(istep+1)+2*(npassing_prev+1)*m+2*npassing_prev+3
      k=ind_start(istep)+2*(npassing+1)*m+2*npassing+3
!
      IF(fact_neg_e(istep).NE.0.d0) THEN
!        nz_beg=nz+1
!
! Toroidal rotation:
!
        DO ipart=1,npassing
          DO kk=1,4
            DO mm=0,lag
              nz=nz+1
!              irow_asymm(nz)=k-ipart
!              icol_asymm(nz)=k-max(0,ipart-2)-kk+2*(npassing+1)*(mm-m)
!              amat_asymm(nz)=x1mm(m,mm)*rhs_mat_energ(kk,ipart,istep)
            ENDDO
          ENDDO
        ENDDO
!
        DO kk=1,4
          DO mm=0,lag
            nz=nz+1
!            irow_asymm(nz)=k-npassing-1
!            icol_asymm(nz)=k-npassing+1-kk+2*(npassing+1)*(mm-m)
!            amat_asymm(nz)=0.5d0*x1mm(m,mm)                                  &
!                       *rhs_mat_energ(kk,npassing+1,istep)
          ENDDO
        ENDDO
!
        DO kk=1,4
          DO mm=0,lag
            nz=nz+1
!            irow_asymm(nz)=k-npassing-1
!            icol_asymm(nz)=k_prev-npassing_prev+1-kk+2*(npassing_prev+1)*(mm-m)
!            amat_asymm(nz)=0.5d0*x1mm(m,mm)                                  &
!                       *rhs_mat_energ(kk,npassing_prev+1,istep+1)
          ENDDO
        ENDDO
!
!        amat_asymm(nz_beg:nz)=-fact_neg_e(istep)*amat_asymm(nz_beg:nz)
!
      ENDIF
!
    ENDDO
!
  ENDDO
!
  nz_asymm=nz
!
  ALLOCATE(irow_asymm(nz_per_neg+1:nz_asymm))
  ALLOCATE(icol_asymm(nz_per_neg+1:nz_asymm))
  ALLOCATE(amat_asymm(nz_per_neg+1:nz_asymm))
!
! Fill the arrays:
!
! Notes on normalization: 
! matrix rhs_mat_energ corresponds to discretization over $\eta$
! of the following function, $-\kappa f/(2 h^\varphi |\lambda|)$
! matrix rhs_mat_energ2 corresponds to discretization over $\eta$
! of the following function, $-\kappa f |\lambda|/(2 h^\varphi)$
! where $\kappa$=collpar
!
!
  !! Modifications by Andreas F. Martitsch (01.04.2015)
  ! Definition of $\kappa$ specified in the comment above originates
  ! from the ntv_booz (Version Nov 2013), which differs by a factor 2
  ! from the quantity collpar that is used here. Therefore, rhs_mat_energ
  ! and rhs_mat_energ2 are corrected by a factor two (no 2 in the denominator)
  !! End Modifications by Andreas F. Martitsch (01.04.2015)
!
  denomjac=aiota*bcovar_theta_hat+bcovar_phi_hat
!
  nz=nz_per_neg
!
! Co-passing: sigma=1
!
  DO istep=ibeg+1,iend
!
!    a1b=bcovar_s(istep)*dbhat_dtheta(istep)/(bhat_mfl(istep)*denomjac)   &
!       -dbhat_ds(istep)/(bhat_mfl(istep)*aiota)
    !! Modifications by Andreas F. Martitsch (14.03.2014)
    ! Optional output (necessary for modeling the magnetic rotation)
    a1b=(bcovar_s_hat_mfl(istep)*dlogbdphi_mfl(istep)/denomjac                   &
       - dlogbds_mfl(istep))/aiota 
    a2b=a1b+2.d0*(dbcovar_theta_hat_ds+dbcovar_phi_hat_ds/aiota                  &
       -          dbcovar_s_hat_dphi_mfl(istep)/aiota)/denomjac
    !! End Modifications by Andreas F. Martitsch (14.03.2014)
!
    npassing_prev=npl(istep-1)
    npassing=npl(istep)
!
    DO m=0,lag
      k_prev=ind_start(istep-1)+2*(npassing_prev+1)*m
      k=ind_start(istep)+2*(npassing+1)*m
!
      IF(fact_pos_e(istep).NE.0.d0) THEN
        nz_beg=nz+1
!
! Toroidal rotation:
!
        DO ipart=1,npassing
          DO kk=1,4
            DO mm=0,lag
              nz=nz+1
              irow_asymm(nz)=k+ipart
              icol_asymm(nz)=k+MAX(0,ipart-2)+kk+2*(npassing+1)*(mm-m)
              !! Modifications by Andreas F. Martitsch (01.04.2015)
              ! rhs_mat_energ and rhs_mat_energ2 are corrected by a factor two
              ! (see comments related to "DEFINITION OF INPUT QUANTITIES" and
              ! "Fill the arrays")
              amat_asymm(nz)=(hatOmegaE*x1mm(m,mm)+hatOmegaB*a1b*x2mm(m,mm)) &
                            * (rhs_mat_energ(kk,ipart,istep)*2.0d0)          &
                            + hatOmegaB*a2b*x2mm(m,mm)                       &
                            * (rhs_mat_energ2(kk,ipart,istep)*2.0d0)
              !! End Modifications by Andreas F. Martitsch (01.04.2015)
            ENDDO
          ENDDO
        ENDDO
!
        DO kk=1,4
          DO mm=0,lag
            nz=nz+1
            irow_asymm(nz)=k+npassing+1
            icol_asymm(nz)=k+npassing-1+kk+2*(npassing+1)*(mm-m)
            !! Modifications by Andreas F. Martitsch (01.04.2015)
            ! rhs_mat_energ and rhs_mat_energ2 are corrected by a factor two
            ! (see comments related to "DEFINITION OF INPUT QUANTITIES" and
            ! "Fill the arrays")
            amat_asymm(nz)=0.5d0*(                                           &
                           (hatOmegaE*x1mm(m,mm)+hatOmegaB*a1b*x2mm(m,mm))   &
                          * (rhs_mat_energ(kk,npassing+1,istep)*2.0d0)       &
                          + hatOmegaB*a2b*x2mm(m,mm)                         &
                          * (rhs_mat_energ2(kk,npassing+1,istep)*2.0d0)      &
                                 )
            !! End Modifications by Andreas F. Martitsch (01.04.2015)
          ENDDO
        ENDDO
!
        DO kk=1,4
          DO mm=0,lag
            nz=nz+1
            irow_asymm(nz)=k+npassing+1
            icol_asymm(nz)=k_prev+npassing_prev-1+kk+2*(npassing_prev+1)*(mm-m)
            !! Modifications by Andreas F. Martitsch (01.04.2015)
            ! rhs_mat_energ and rhs_mat_energ2 are corrected by a factor two
            ! (see comments related to "DEFINITION OF INPUT QUANTITIES" and
            ! "Fill the arrays")
            amat_asymm(nz)=0.5d0*(                                             &
                           (hatOmegaE*x1mm(m,mm)+hatOmegaB*a1b*x2mm(m,mm))     &
                          * (rhs_mat_energ(kk,npassing_prev+1,istep-1)*2.0d0)  &
                          + hatOmegaB*a2b*x2mm(m,mm)                           &
                          * (rhs_mat_energ2(kk,npassing_prev+1,istep-1)*2.0d0) &
                                 )
            !! End Modifications by Andreas F. Martitsch (01.04.2015)
          ENDDO
        ENDDO
!
        amat_asymm(nz_beg:nz)=-fact_pos_e(istep)*amat_asymm(nz_beg:nz)
!
      ENDIF
!
    ENDDO
!
  ENDDO
!
! Counter-passing: sigma=-1
!
  DO istep=ibeg,iend-1
!
!    a1b=bcovar_s(istep)*dbhat_dtheta(istep)/(bhat_mfl(istep)*denomjac)   &
!       -dbhat_ds(istep)/(bhat_mfl(istep)*aiota)
     !! Modifications by Andreas F. Martitsch (14.03.2014)
     ! Optional output (necessary for modeling the magnetic rotation)
     a1b=(bcovar_s_hat_mfl(istep)*dlogbdphi_mfl(istep)/denomjac                   &
          - dlogbds_mfl(istep))/aiota 
     a2b=a1b+2.d0*(dbcovar_theta_hat_ds+dbcovar_phi_hat_ds/aiota                  &
          -          dbcovar_s_hat_dphi_mfl(istep)/aiota)/denomjac
     !! End Modifications by Andreas F. Martitsch (14.03.2014)
!
    npassing_prev=npl(istep+1)
    npassing=npl(istep)
!
    DO m=0,lag
      k_prev=ind_start(istep+1)+2*(npassing_prev+1)*m+2*npassing_prev+3
      k=ind_start(istep)+2*(npassing+1)*m+2*npassing+3
!
      IF(fact_neg_e(istep).NE.0.d0) THEN
        nz_beg=nz+1
!
! Toroidal rotation:
!
        DO ipart=1,npassing
          DO kk=1,4
            DO mm=0,lag
              nz=nz+1
              irow_asymm(nz)=k-ipart
              icol_asymm(nz)=k-MAX(0,ipart-2)-kk+2*(npassing+1)*(mm-m)
              !! Modifications by Andreas F. Martitsch (01.04.2015)
              ! rhs_mat_energ and rhs_mat_energ2 are corrected by a factor two
              ! (see comments related to "DEFINITION OF INPUT QUANTITIES" and
              ! "Fill the arrays")
              amat_asymm(nz)=(hatOmegaE*x1mm(m,mm)+hatOmegaB*a1b*x2mm(m,mm)) &
                            * (rhs_mat_energ(kk,ipart,istep)*2.0d0)          &
                            + hatOmegaB*a2b*x2mm(m,mm)                       &
                            * (rhs_mat_energ2(kk,ipart,istep)*2.0d0)
              !! End Modifications by Andreas F. Martitsch (01.04.2015)
            ENDDO
          ENDDO
        ENDDO
!
        DO kk=1,4
          DO mm=0,lag
            nz=nz+1
            irow_asymm(nz)=k-npassing-1
            icol_asymm(nz)=k-npassing+1-kk+2*(npassing+1)*(mm-m)
            !! Modifications by Andreas F. Martitsch (01.04.2015)
            ! rhs_mat_energ and rhs_mat_energ2 are corrected by a factor two
            ! (see comments related to "DEFINITION OF INPUT QUANTITIES" and
            ! "Fill the arrays")
            amat_asymm(nz)=0.5d0*(                                           &
                           (hatOmegaE*x1mm(m,mm)+hatOmegaB*a1b*x2mm(m,mm))   &
                          * (rhs_mat_energ(kk,npassing+1,istep)*2.0d0)       &
                          + hatOmegaB*a2b*x2mm(m,mm)                         &
                          * (rhs_mat_energ2(kk,npassing+1,istep)*2.0d0)      &
                                 )
            !! End Modifications by Andreas F. Martitsch (01.04.2015)
          ENDDO
        ENDDO
!
        DO kk=1,4
          DO mm=0,lag
            nz=nz+1
            irow_asymm(nz)=k-npassing-1
            icol_asymm(nz)=k_prev-npassing_prev+1-kk+2*(npassing_prev+1)*(mm-m)
            !! Modifications by Andreas F. Martitsch (01.04.2015)
            ! rhs_mat_energ and rhs_mat_energ2 are corrected by a factor two
            ! (see comments related to "DEFINITION OF INPUT QUANTITIES" and
            ! "Fill the arrays")
            amat_asymm(nz)=0.5d0*(                                             &
                           (hatOmegaE*x1mm(m,mm)+hatOmegaB*a1b*x2mm(m,mm))     &
                          * (rhs_mat_energ(kk,npassing_prev+1,istep+1)*2.0d0)  &
                          + hatOmegaB*a2b*x2mm(m,mm)                           &
                          * (rhs_mat_energ2(kk,npassing_prev+1,istep+1)*2.0d0) &
                                 )
            !! End Modifications by Andreas F. Martitsch (01.04.2015)
          ENDDO
        ENDDO
!
        amat_asymm(nz_beg:nz)=-fact_neg_e(istep)*amat_asymm(nz_beg:nz)
!
      ENDIF
!
    ENDDO
!
  ENDDO
!
  DEALLOCATE(x1mm,x2mm)
!
!! Modifications by Andreas F. Martitsch (02.09.2014)
! For the computation of hatOmegaE from the profile
! one must know the solution of the axisymmetric
! equation set.
  ALLOCATE(ipcol(ncol),bvec_sp(ncol))
!
! Solve the axisymmetric equation set:
!
  IF(ALLOCATED(qflux_symm)) DEALLOCATE(qflux_symm)
  IF(nobounceaver) THEN
!
    nz=nz_symm+nz_regper
!
    ALLOCATE(irow(nz),icol(nz),amat_sp(nz))
!
    irow(1:nz_symm)=irow_symm
    icol(1:nz_symm)=icol_symm
    amat_sp(1:nz_symm)=amat_symm
    IF(nz_regper.GT.0) THEN
      irow(nz_symm+1:nz)=irow_regper
      icol(nz_symm+1:nz)=icol_regper
      amat_sp(nz_symm+1:nz)=amat_regper
    ENDIF
!
    !! Modifications by Andreas F. Martitsch (28.08.2014)
    ! geodesic curvature for the axisymmetric field computed by
    ! external routines
    !geodcu_forw=geodcu_mfl
    !geodcu_back=geodcu_mfl
    ! computation of geodesic curvature according to
    ! $\|{\nabla}s\| k_{G0} = - \frac{B_\ph}{\iota B_{\tht}+B_{\ph}}
    ! \frac{B_{\rm ref}}{\psi_{\rm tor}^{a}} \difp{B_0}{\tht}$
    denomjac=-scalefac_kG*bcovar_phi_hat/(aiota*bcovar_theta_hat+bcovar_phi_hat)
    ! geodcu_forw used for computation of q_rip(1:npassing+1,istep,1),
    ! which in turn enters the source_vector
    geodcu_forw=denomjac*dlogbdphi_mfl*bhat_mfl
    ! geodcu_back used for the computation of convol_flux, which enters q_flux
    ! via flux_vector(1,:) and flux_vector(3,:)
    !--> Computation of D31/D32 not affected ( flux_vector(2,:) determined by convol_curr )
    geodcu_back=geodcu_forw
    !! End Modifications by Andreas F. Martitsch (28.08.2014)
!
    CALL source_flux
!
    CALL solve_eqs(.TRUE.)
!
    DEALLOCATE(irow,icol,amat_sp)
!
! Caution!!! factor 2 is not needed!!!
    qflux=2.d0*qflux
    OPEN(1234,file='qflux_symm.dat')
    WRITE(1234,*) conl_over_mfp/boozer_iota,-qflux(2,2:3),-qflux(3,2:3)
    CLOSE(1234)

    !! Modifications by Andreas F. Martitsch (28.07.2014)
    ! Save here the solution of the axisymmetric case (if non-axisymmetric
    ! solution is also computed). Otherwise exit here ripple_solver and
    ! return to calling routine propagator_solver.
    IF(isw_qflux_NA .EQ. 0) THEN
       ! compute qflux only for the axisymmetric equation set
       ! --> qflux already available from prop_a%p%qflux
       RETURN
    ELSE IF(isw_qflux_NA .EQ. 1) THEN
       ! save qflux for the axisymmetric equation set
       ! and proceed with the solution of the non-axisymmetric
       ! equation set (stored within prop_a%p%qflux)
       ALLOCATE(qflux_symm(3,3))
       qflux_symm=qflux
    ELSE
       STOP "Invalid input for isw_qflux_symm (0/1)!"
    END IF
    !! End Modifications by Andreas F. Martitsch (28.07.2014)
!
    ALLOCATE(f0_coll(n_2d_size,3),f0_ttmp(n_2d_size,3))
    f0_coll=0.d0
    f0_ttmp=0.d0
!
! Here f0_coll=$-\frac{1}{h^\varphi}\sum_{m^\prime}\hat L_{mm^\prime}^c
!                bar f_{m^\prime}^{\sigma (k)}$ :
!
    DO nz=1,nz_coll
      !! Modification by Andreas F. Martitsch (17.07.2015)
      ! fixed warning: Possible change of value in conversion
      ! from COMPLEX(8) to REAL(8)
      f0_coll(irow_coll(nz),:)=f0_coll(irow_coll(nz),:)+amat_coll(nz)  &
                              *REAL(source_vector(icol_coll(nz),:))
      !! End Modification by Andreas F. Martitsch (17.07.2015)
    ENDDO
!
    IF(isw_intp.EQ.1) THEN
      ALLOCATE(bvec_iter(ncol),bvec_prev(ncol))
!
      DO i=1,3
        bvec_prev=source_vector(:,i)
!
        CALL integral_part(bvec_prev,bvec_iter)
!
        !! Modification by Andreas F. Martitsch (17.07.2015)
        ! fixed warning: Possible change of value in conversion
        ! from COMPLEX(8) to REAL(8)
        f0_coll(:,i)=f0_coll(:,i)-REAL(bvec_iter)
        !! End Modification by Andreas F. Martitsch (17.07.2015)
      ENDDO
!
      DEALLOCATE(bvec_iter,bvec_prev)
    ENDIF
!
! Here f0_ttmp=$-\sigma \eta \difp{}{\eta} f_{m^\prime}^{\sigma (k)}$ :
!
    DO nz=1,nz_ttmp
      !! Modification by Andreas F. Martitsch (17.07.2015)
      ! fixed warning: Possible change of value in conversion
      ! from COMPLEX(8) to REAL(8)
      f0_ttmp(irow_ttmp(nz),:)=f0_ttmp(irow_ttmp(nz),:)+amat_ttmp(nz)  &
                              *REAL(source_vector(icol_ttmp(nz),:))
      !! End Modification by Andreas F. Martitsch (17.07.2015)
    ENDDO
!
    IF(colltest) THEN
!
      CALL source_flux
!
      istep=(ibeg+iend)/3
      npassing=npl(istep)
!
      DO m=0,lag
        k=ind_start(istep)+2*(npassing+1)*m
        DO i=k+1,k+2*(npassing+1)
          WRITE(7000+m,*) sngl(f0_coll(i,:)+f0_ttmp(i,:))
          WRITE(8000+m,*) sngl(REAL(source_vector(i,:)))
        ENDDO
      ENDDO
    ENDIF
!
    IF(ttmptest) THEN
!
! Plot the mirroring force $-\lambda \eta \difp{}{\eta} f_{m^\prime}^{\sigma (k)}$ 
! as function of $\lambda$ :
!
      istep=(ibeg+iend)/3
!
      CALL plotsource(9000,f0_ttmp)
!
    ENDIF
!
  ENDIF
!
! End solve the axisymmetric equation set
!
!! End Modifications by Andreas F. Martitsch (02.09.2014)
!
! Solve the non-axisymmetric equation set:
!
  !! Modifications by Andreas F. Martitsch (13.06.2014)
  ! quantities of the perturbation field are extracted
  ! from the Boozer file
  ! (-> no further need to specify them manually +
  ! tests using an artificial perturbation field can
  ! be done inside tmp section (see above))  
  !bnoverb0=bnoverb0*EXP(imun*m_phi*phi_mfl)
  !dbnoverb0_dtheta=dbnoverb0_dtheta*EXP(imun*m_phi*phi_mfl)
  !! End Modifications by Andreas F. Martitsch (13.06.2014)
  denomjac=aiota*bcovar_theta_hat+bcovar_phi_hat
!
  geodcu_back=scalefac_kG*imun*m_phi*bnoverb0*bhat_mfl
  geodcu_forw=geodcu_back
!
  !! Modifications by Andreas F. Martitsch (10.12.2014)
  ! Ware pinch effect
  DO istep=ibeg,iend
     q_rip(:,istep,2)=-2.d0*q_rip(:,istep,2)*bnoverb0(istep)
  END DO
  !! End Modifications by Andreas F. Martitsch (10.12.2014)
!  
  IF(nobounceaver) THEN
    !! Modifications by Andreas F. Martitsch (18.08.2014)
    ! derivative along the periodic Boozer angle theta has
    ! been redefined to a derivative along the field line (phi_mfl)
    ! (changes concerning dbnoverb0_dtheta are required!)
    !geodcu_forw=geodcu_forw-scalefac_kG*(bcovar_phi_hat/denomjac)    &
    !           *(aiota*dbnoverb0_dtheta+bnoverb0*(imun*m_phi     &
    !           - dlogbdphi_mfl*aiota))*bhat_mfl
    !geodcu_forw=geodcu_forw-scalefac_kG*(bcovar_phi_hat/denomjac)&
    !           *(dbnoverb0_dphi_mfl-bnoverb0*dlogbdphi_mfl*aiota)*bhat_mfl
    ! incorrect factor aiota removed from last term 
    geodcu_forw=geodcu_forw-scalefac_kG*(bcovar_phi_hat/denomjac)&
         *(dbnoverb0_dphi_mfl-bnoverb0*dlogbdphi_mfl)*bhat_mfl
    !! End Modifications by Andreas F. Martitsch (28.08.2014)    
    !PRINT *,'geodcu_forw: ',geodcu_forw
  ENDIF
!
  CALL source_flux
!
  IF(nobounceaver) THEN
    ALLOCATE(ttmpfact(ibeg:iend))
    !! Modifications by Andreas F. Martitsch (13.06.2014)
    ! derivative along the periodic Boozer angle theta has
    ! been redefined to a derivative along the field line (phi_mfl)
    ! (changes concerning dbnoverb0_dtheta are required!)
    !ttmpfact=aiota*dbnoverb0_dtheta+imun*m_phi*bnoverb0
    ttmpfact=dbnoverb0_dphi_mfl
    !! End Modifications by Andreas F. Martitsch (13.06.2014)
!
    CALL add_f01_source
!
    DEALLOCATE(ttmpfact)
  ENDIF
!
  expforw=EXP(imun*m_phi*(phi_mfl(iend)-phi_mfl(ibeg)))
  expbackw=(1.d0,0.d0)/expforw
  perbou_pos=(1.d0,0.d0)-expbackw
  perbou_neg=(1.d0,0.d0)-expforw
!
  nz=nz_asymm
!
  ALLOCATE(irow(nz),icol(nz),amat_sp(nz))
!
  irow(1:nz_symm)=irow_symm
  icol(1:nz_symm)=icol_symm
  amat_sp(1:nz_symm)=amat_symm
!
  irow(nz_symm+1:nz_per_pos)=irow_per_pos(nz_symm+1:nz_per_pos)
  icol(nz_symm+1:nz_per_pos)=icol_per_pos(nz_symm+1:nz_per_pos)
  amat_sp(nz_symm+1:nz_per_pos)=perbou_pos
!
  irow(nz_per_pos+1:nz_per_neg)=irow_per_neg(nz_per_pos+1:nz_per_neg)
  icol(nz_per_pos+1:nz_per_neg)=icol_per_neg(nz_per_pos+1:nz_per_neg)
  amat_sp(nz_per_pos+1:nz_per_neg)=perbou_neg
!
  irow(nz_per_neg+1:nz_asymm)=irow_asymm(nz_per_neg+1:nz_asymm)
  icol(nz_per_neg+1:nz_asymm)=icol_asymm(nz_per_neg+1:nz_asymm)
  amat_sp(nz_per_neg+1:nz_asymm)=amat_asymm(nz_per_neg+1:nz_asymm)*rotfactor
!
  CALL solve_eqs(.TRUE.)
!istep=(ibeg+iend)/2
!call plotsource(10000,real(source_vector))
!call plotsource(11000,dimag(source_vector))
!
  OPEN(1234,file='qflux_ntv.dat')
  WRITE(1234,*) conl_over_mfp/boozer_iota,-qflux(1,1:3:2),-qflux(3,1:3:2)
  CLOSE(1234)
!! Modifications by Andreas F. Martitsch (17.12.2013)
!! Create standard output (efinal.dat,fulltransp.dat)
RETURN
!! or stop here
!! STOP
!! End Modifications by Andreas F. Martitsch (17.12.2013)
!
  DO m=0,lag
    DO kk=1,3
      k=ind_start(ibeg)+2*npass_l*m
      !! Modification by Andreas F. Martitsch (17.07.2015)
      ! fixed warning: Possible change of value in conversion
      ! from COMPLEX(8) to REAL(8)
      source_m(npass_l*m+1:npass_l*m+npass_l,kk) &
          =REAL(source_vector(k+2*npass_l:k+npass_l+1:-1,kk))
      !! End Modification by Andreas F. Martitsch (17.07.2015)
      k=ind_start(iend)+2*npass_r*m
      !! Modification by Andreas F. Martitsch (17.07.2015)
      ! fixed warning: Possible change of value in conversion
      ! from COMPLEX(8) to REAL(8)
      source_p(npass_r*m+1:npass_r*m+npass_r,kk) &
          =REAL(source_vector(k+1:k+npass_r,kk))
      !! End Modification by Andreas F. Martitsch (17.07.2015)
    ENDDO
  ENDDO
!
  IF(iplot.EQ.1.OR.isw_axisymm.EQ.1) THEN
!
    flux_p=0.d0
    flux_m=0.d0
    amat_plus_plus=0.d0
    amat_plus_minus=0.d0
    amat_minus_minus=0.d0
    amat_minus_plus=0.d0
!
    DEALLOCATE(flux_vector,source_vector,irow,icol,amat_sp,ipcol,bvec_sp)
    IF(isw_intp.EQ.1) DEALLOCATE(bvec_iter,bvec_prev)
    DEALLOCATE(deriv_coef,enu_coef,alambd,Vg_vp_over_B,scalprod_pleg)
    DEALLOCATE(alampow,vrecurr,dellampow,convol_polpow,pleg_bra,pleg_ket)
    DEALLOCATE(enu_coef2,dellampow2,rhs_mat_energ2)
    DEALLOCATE(npl,rhs_mat_fzero,rhs_mat_lorentz,rhs_mat_energ,q_rip,q_rip_1)
    DEALLOCATE(fun_coef,ttmp_mat)
    DEALLOCATE(convol_flux,convol_curr,ind_start,convol_flux_0)
    DEALLOCATE(phi_mfl,bhat_mfl,geodcu_mfl,h_phi_mfl,dlogbdphi_mfl,eta)
    !! Modifications by Andreas F. Martitsch (14.03.2014)
    ! Optional output (necessary for modeling the magnetic rotation)
    ! --> Deallocate the arrays
    DEALLOCATE(dlogbds_mfl,bcovar_s_hat_mfl,dbcovar_s_hat_dphi_mfl)
    !! End Modifications by Andreas F. Martitsch (14.03.2014)
    DEALLOCATE(delt_pos,delt_neg,fact_pos_b,fact_neg_b,fact_pos_e,fact_neg_e)
    IF (ALLOCATED(eta_prev)) DEALLOCATE(eta_prev)
    IF (ALLOCATED(eta_next)) DEALLOCATE(eta_next)
    IF (ALLOCATED(amat_z)) DEALLOCATE(amat_z,bvec_lapack_z,ipivot)
    IF (ALLOCATED(scalprod)) DEALLOCATE(scalprod,basevec_p)
!
    CALL CPU_TIME(time_solver)
    PRINT *,'solving completed       ',time_solver - time_factorization,' sec'
!
    RETURN
!
  ENDIF
!
! Calculation of propagators:
!
  DO m=0,lag
    k=ind_start(ibeg)+2*npass_l*m
    DO i=1,npass_l
      bvec_sp=0.d0
      bvec_sp(k+i)=1.d0
!
      CALL sparse_solve(nrow,ncol,nz,irow(1:nz),ipcol,amat_sp(1:nz),   &
                        bvec_sp,iopt)
!
! integral part:
!
      IF(isw_intp.EQ.1) THEN
!
        CALL iterate(bvec_sp)
!
        IF(precond) THEN
          bvec_lapack_z(:,1)=scalprod
!
          CALL zgesv(ndim,1,amat_z,ndim,ipivot,bvec_lapack_z,ndim,info)
!
          DO mm=0,lag
            bvec_sp=bvec_sp+bvec_lapack_z(mm+1,1)*basevec_p(:,mm)
          ENDDO
!
        ENDIF
!
      ENDIF
!
      DO mm=0,lag
        kk=ind_start(iend)+2*npass_r*mm
        !! Modification by Andreas F. Martitsch (17.07.2015)
        ! fixed warning: Possible change of value in conversion
        ! from COMPLEX(8) to REAL(8)
        amat_plus_plus(npass_r*mm+1:npass_r*mm+npass_r,npass_l*m+i)    &
                     =REAL(bvec_sp(kk+1:kk+npass_r))
        !! End Modification by Andreas F. Martitsch (17.07.2015)
      ENDDO
      DO mm=0,lag
        kk=ind_start(ibeg)+2*npass_l*mm
        !! Modification by Andreas F. Martitsch (17.07.2015)
        ! fixed warning: Possible change of value in conversion
        ! from COMPLEX(8) to REAL(8)
        amat_plus_minus(npass_l*mm+1:npass_l*mm+npass_l,npass_l*m+i)   &
                     =REAL(bvec_sp(kk+2*npass_l:kk+npass_l+1:-1))
        !! End Modification by Andreas F. Martitsch (17.07.2015)
      ENDDO
      !! Modification by Andreas F. Martitsch (17.07.2015)
      ! fixed warning: Possible change of value in conversion
      ! from COMPLEX(8) to REAL(8)
      flux_p(:,npass_l*m+i)=REAL(MATMUL(flux_vector,bvec_sp(:)))
      !! End Modification by Andreas F. Martitsch (17.07.2015)
    ENDDO
  ENDDO
!
  DO m=0,lag
    k=ind_start(iend)+2*npass_r*m
    DO i=1,npass_r
      bvec_sp=0.d0
      bvec_sp(k+2*npass_r+1-i)=1.d0
!
      CALL sparse_solve(nrow,ncol,nz,irow(1:nz),ipcol,amat_sp(1:nz),   &
                        bvec_sp,iopt)
!
! integral part:
!
      IF(isw_intp.EQ.1) THEN
!
        CALL iterate(bvec_sp)
!
        IF(precond) THEN
          bvec_lapack_z(:,1)=scalprod
!
          CALL zgesv(ndim,1,amat_z,ndim,ipivot,bvec_lapack_z,ndim,info)
!
          DO mm=0,lag
            bvec_sp=bvec_sp+bvec_lapack_z(mm+1,1)*basevec_p(:,mm)
          ENDDO
!
        ENDIF
!
      ENDIF
!
      DO mm=0,lag
        kk=ind_start(iend)+2*npass_r*mm
        !! Modification by Andreas F. Martitsch (17.07.2015)
        ! fixed warning: Possible change of value in conversion
        ! from COMPLEX(8) to REAL(8)
        amat_minus_plus(npass_r*mm+1:npass_r*mm+npass_r,npass_r*m+i)   &
                     =REAL(bvec_sp(kk+1:kk+npass_r))
        !! End Modification by Andreas F. Martitsch (17.07.2015)) 
      ENDDO
      DO mm=0,lag
        kk=ind_start(ibeg)+2*npass_l*mm
        !! Modification by Andreas F. Martitsch (17.07.2015)
        ! fixed warning: Possible change of value in conversion
        ! from COMPLEX(8) to REAL(8) 
        amat_minus_minus(npass_l*mm+1:npass_l*mm+npass_l,npass_r*m+i)  &
                     =REAL(bvec_sp(kk+2*npass_l:kk+npass_l+1:-1))
        !! End Modification by Andreas F. Martitsch (17.07.2015)) 
      ENDDO
      !! Modification by Andreas F. Martitsch (17.07.2015)
      ! fixed warning: Possible change of value in conversion
      ! from COMPLEX(8) to REAL(8) 
      flux_m(:,npass_r*m+i)=REAL(MATMUL(flux_vector,bvec_sp(:)))
      !! End Modification by Andreas F. Martitsch (17.07.2015))
    ENDDO
  ENDDO
!
  iopt=3
!
  CALL sparse_solve(nrow,ncol,nz,irow(1:nz),ipcol,amat_sp(1:nz),bvec_sp,iopt)
!
  DEALLOCATE(flux_vector,source_vector,irow,icol,amat_sp,ipcol,bvec_sp)
  IF(isw_intp.EQ.1) THEN 
    DEALLOCATE(bvec_iter,bvec_prev)
    IF(precond) DEALLOCATE(amat_z,bvec_lapack_z,ipivot,scalprod,basevec_p)
  ENDIF
!
!
!
  CALL CPU_TIME(time_solver)
  PRINT *,'solving completed       ',time_solver - time_factorization,' sec'
!
!
!
!
!
  nplp1=npart_loc+1
  ntotsize=(lag+1)*nplp1
!
PRINT *,'npart_loc = ',npart_loc,' npass_l = ',npass_l,' npass_r = ',npass_r
!
  ndim=ntotsize
!
IF(fieldpropagator%tag.EQ.2) THEN
OPEN(112,file='onsager.dat')
OPEN(113,file='levhist.dat')
OPEN(111,file='qfluxhist.dat')
ELSE
OPEN(112,file='onsager.dat',position='append')
OPEN(113,file='levhist.dat',position='append')
OPEN(111,file='qfluxhist.dat',position='append')
ENDIF
facnorm_m=1.d0
facnorm_p=1.d0
!WRITE(112,*) fieldpropagator%tag,qflux(1:2,1),qflux(1:2,2),npart_loc &
!             ,b_max_l,b_max_r,facnorm_p,facnorm_m,ignore_boundary_layer_new &
!             ,ignore_boundary_layer,ifilter_l,ifilter_r
WRITE(112,*) fieldpropagator%tag,qflux(1:2,1),qflux(1:2,2),npart_loc &
             ,b_max_l,b_max_r,facnorm_p,facnorm_m,ignore_lb &
             ,ignore_rb,modify_bl,modify_br
!             ,ignore_rb,ifilter_l,ifilter_r
WRITE(111,*) fieldpropagator%tag,qflux
DO i=0,npart_loc
WRITE(113,*) fieldpropagator%tag,eta(i)
ENDDO
CLOSE(112)
CLOSE(113)
CLOSE(111)
PRINT *,qflux(1:2,1),qflux(1:2,2)
!!PAUSE
!
GOTO 1
!
OPEN(111,file='flux_p.dat')
WRITE(111,'(4(1x,e12.5))') (alam_l(i),flux_p(:,i) ,i=1,npass_l)
CLOSE(111)
OPEN(111,file='flux_m.dat')
WRITE(111,'(4(1x,e12.5))') (alam_r(i),flux_m(:,i) ,i=1,npass_r)
CLOSE(111)
OPEN(111,file='source_p.dat')
WRITE(111,'(4(1x,e12.5))') (alam_r(i),source_p(i,:)/delta_eta_r(i) ,i=1,npass_r)
CLOSE(111)
OPEN(111,file='source_m.dat')
WRITE(111,'(4(1x,e12.5))') (alam_l(i),source_m(i,:)/delta_eta_l(i) ,i=1,npass_l)
CLOSE(111)
!amat_plus_plus=0.d0
!amat_plus_minus=0.d0
!amat_minus_plus=0.d0
!amat_minus_minus=0.d0
!do i=1,min(npass_l,npass_r)
!amat_plus_plus(i,i)=1.d0
!amat_minus_minus(i,i)=1.d0
!enddo
!do i=npass_r+1,npass_l
!amat_plus_minus(i,i)=0.d0 !1.d0
!enddo
!do i=npass_l+1,npass_r
!amat_minus_plus(i,i)=0.d0 !1.d0
!enddo
!
OPEN(111,file='amat_p_p.dat')
DO i=1,npass_r
WRITE(111,*) amat_plus_plus(i,1:npass_l)
!WRITE(111,*) 0.5d0*(eta(i)+eta(i-1)),amat_plus_plus(i,1:npass_l)
ENDDO
CLOSE(111)
OPEN(111,file='amat_p_m.dat')
DO i=1,npass_l
WRITE(111,*) amat_plus_minus(i,1:npass_l)
!WRITE(111,*) 0.5d0*(eta(i)+eta(i-1)),amat_plus_minus(i,1:npass_l)
ENDDO
CLOSE(111)
OPEN(111,file='amat_m_p.dat')
DO i=1,npass_r
WRITE(111,*) amat_minus_plus(i,1:npass_r)
!WRITE(111,*) 0.5d0*(eta(i)+eta(i-1)),amat_minus_plus(i,1:npass_r)
ENDDO
CLOSE(111)
OPEN(111,file='amat_m_m.dat')
DO i=1,npass_l
WRITE(111,*) amat_minus_minus(i,1:npass_r)
!WRITE(111,*) 0.5d0*(eta(i)+eta(i-1)),amat_minus_minus(i,1:npass_r)
ENDDO
CLOSE(111)
!open(111,file='lambda_l.dat')
!do i=1,npass_l
!write (111,*) -alam_l(i),delta_eta_l(i)
!enddo
!do i=npass_l,1,-1
!write (111,*) alam_l(i),delta_eta_l(i)
!enddo
!close(111)
!PAUSE 'written'
1 CONTINUE
!
PRINT *,' '
  DEALLOCATE(deriv_coef,enu_coef,alambd,Vg_vp_over_B,scalprod_pleg)
  DEALLOCATE(alampow,vrecurr,dellampow,convol_polpow,pleg_bra,pleg_ket)
  DEALLOCATE(enu_coef2,dellampow2,rhs_mat_energ2)
  DEALLOCATE(npl,rhs_mat_fzero,rhs_mat_lorentz,rhs_mat_energ,q_rip,q_rip_1)
  DEALLOCATE(fun_coef,ttmp_mat)
  DEALLOCATE(convol_flux,convol_curr,ind_start,convol_flux_0)
  DEALLOCATE(delt_pos,delt_neg,fact_pos_b,fact_neg_b,fact_pos_e,fact_neg_e)
!
  !------------------------------------------------------------------------
  ! END SERGEI
  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
  ! final deallocation of initial things provided by winny
  IF (ALLOCATED(phi_mfl)) DEALLOCATE(phi_mfl)
  IF (ALLOCATED(bhat_mfl)) DEALLOCATE(bhat_mfl)
  IF (ALLOCATED(geodcu_mfl)) DEALLOCATE(geodcu_mfl)
  IF (ALLOCATED(h_phi_mfl)) DEALLOCATE(h_phi_mfl)
  IF (ALLOCATED(dlogbdphi_mfl)) DEALLOCATE(dlogbdphi_mfl)
  !! Modifications by Andreas F. Martitsch (14.03.2014)
  ! Optional output (necessary for modeling the magnetic rotation)
  ! --> Deallocate the arrays
  IF (ALLOCATED(dlogbds_mfl)) DEALLOCATE(dlogbds_mfl)
  IF (ALLOCATED(bcovar_s_hat_mfl)) DEALLOCATE(bcovar_s_hat_mfl)
  IF (ALLOCATED(dbcovar_s_hat_dphi_mfl)) DEALLOCATE(dbcovar_s_hat_dphi_mfl)
  !! End Modifications by Andreas F. Martitsch (14.03.2014)
  IF (ALLOCATED(eta)) DEALLOCATE(eta)
  IF (ALLOCATED(eta_prev)) DEALLOCATE(eta_prev)
  IF (ALLOCATED(eta_next)) DEALLOCATE(eta_next)

  !------------------------------------------------------------------------
  RETURN
!
  CONTAINS
!
!------------------------------------------------------------------------
!
    SUBROUTINE plotsource(iunit_base,sourcevec_tmp)
!
    INTEGER :: iunit_base
    DOUBLE PRECISION, DIMENSION(n_2d_size,3) :: sourcevec_tmp
!
    npassing=npl(istep)
    delta_eta=eta(1:npassing)-eta(0:npassing-1)
    eta0=1.d0/bhat_mfl(istep)
    DO m=0,lag
      k=ind_start(istep)+2*(npassing+1)*m
      DO i=1,npassing+1
        IF(i.LE.npassing) THEN
          alambd_save1=0.5d0*(alambd(i,istep)+alambd(i-1,istep))
          WRITE(iunit_base+m,*) -alambd_save1,alambd_save1 &
                *sngl(sourcevec_tmp(k+2*(npassing+1)-i+1,:))/delta_eta(i)
        ELSE
          alambd_save1=0.5d0*alambd(i-1,istep)
          WRITE(iunit_base+m,*) -alambd_save1,alambd_save1 &
                *sngl(sourcevec_tmp(k+2*(npassing+1)-i+1,:))/(eta0-eta(i-1))
        ENDIF
      ENDDO
      DO i=npassing+1,1,-1
        IF(i.LE.npassing) THEN
          alambd_save1=0.5d0*(alambd(i,istep)+alambd(i-1,istep))
          WRITE(iunit_base+m,*) alambd_save1,alambd_save1  &
                *sngl(sourcevec_tmp(k+i,:))/delta_eta(i)
        ELSE
          alambd_save1=0.5d0*alambd(i-1,istep)
          WRITE(iunit_base+m,*) alambd_save1,alambd_save1  &
               *sngl(sourcevec_tmp(k+i,:))/(eta0-eta(i-1))
        ENDIF
      ENDDO
    ENDDO
!
    END SUBROUTINE plotsource
!
!------------------------------------------------------------------------
!
    SUBROUTINE solve_eqs(clean)
!
! Solve the linear equation set:
!
    IMPLICIT NONE    
    LOGICAL :: clean
!
    IF(isw_intp.EQ.1) ALLOCATE(bvec_iter(ncol),bvec_prev(ncol))
!
    CALL  remap_rc(nz,nz_sq,irow,icol,amat_sp)
    !! Modifications by Andreas F. Martitsch (21.05.2014)
    ! Write sparse matrix into a file (-> used for diagnosis - test iterative solver from matlab)
    !OPEN(unit=120514,file='kinetic_eq_linearized.dat')
    !DO m=1,nz_sq
    !   WRITE(unit=120514,fmt=*) irow(m),icol(m),dimag(amat_sp(m)),dreal(amat_sp(m))
    !END DO
    !CLOSE(unit=120514)
    !OPEN(unit=120514,file='kinetic_eq_sources.dat')
    !DO m=1,nrow
    !   WRITE(unit=120514,fmt=*) dimag(source_vector(m,1:3)), dreal(source_vector(m,1:3))
    !END DO
    !CLOSE(unit=120514)
    !STOP 'Exit solve_eqs ...'
    !! End Modifications by Andreas F. Martitsch (21.05.2014)
!
    PRINT *,'system size = ',n_2d_size
    PRINT *,'non-zeros before and after truncation = ',nz,nz_sq
    nz=nz_sq
    !PRINT *,icol(1:nz)
    !PRINT *,ipcol
!
    CALL column_full2pointer(icol(1:nz),ipcol)
!
! There are now three different calls sparse_solve
!   iopt = 1 ; factorization
!   iopt = 2 ; solve
!   iopt = 3 ; free memory
! without iopt or with iopt = 0: old behaviour (factorize,solve,free)
!
! factorization:
!
    bvec_sp=0.d0
    !PRINT *,bvec_sp
    iopt=1
!
    CALL CPU_TIME(time_start)
    CALL sparse_solve(nrow,ncol,nz,irow(1:nz),ipcol,amat_sp(1:nz),       &
                      bvec_sp,iopt)
!
    CALL CPU_TIME(time_factorization)
    PRINT *,'factorization completed ',time_factorization - time_start,' sec'
!
    iopt=2
!
! Solution of inhomogeneus equation (account of sources):
!
    CALL sparse_solve(nrow,ncol,nz,irow(1:nz),ipcol,amat_sp(1:nz),       &
                      source_vector(:,1:3),iopt)
!
! integral part:
!
    IF(isw_intp.EQ.1) THEN
!
      IF(precond) THEN
        ndim=lag+1
        ALLOCATE(amat_z(ndim,ndim),bvec_lapack_z(ndim,3),ipivot(ndim))
        ALLOCATE(scalprod(0:lag),basevec_p(n_2d_size,0:lag))
!
        amat_z=0.d0
!
        DO m=0,lag
          amat_z(m+1,m+1)=amat_z(m+1,m+1)+1.d0
!
          CALL basevec_precond(m,basevec_p(:,m))
!
          PRINT *,'basevec',m,':'
!
          CALL iterate(basevec_p(:,m))
!
          amat_z(:,m+1)=amat_z(:,m+1)-scalprod
        ENDDO
!
      ENDIF
!
      DO k=1,3
!
        PRINT *,'source',k,':'
!
        CALL iterate(source_vector(:,k))
!
        IF(precond) THEN
          bvec_lapack_z(:,k)=scalprod
        ENDIF
!
      ENDDO
!
      IF(precond) THEN
!
        CALL zgesv(ndim,3,amat_z,ndim,ipivot,bvec_lapack_z,ndim,info)
!
        DO m=0,lag
!
          DO k=1,3
            source_vector(:,k)=source_vector(:,k)                   &
                              +bvec_lapack_z(m+1,k)*basevec_p(:,m)
          ENDDO
!
        ENDDO
!
        DEALLOCATE(bvec_lapack_z)
        ALLOCATE(bvec_lapack_z(ndim,1))
      ENDIF
!
    ENDIF
!
! Plotting:
!
    IF(iplot.EQ.1) THEN
      nphiplot=200    !serves as an upper limit for data output
      delphiplot=MAXVAL(phi_mfl(ibeg+1:iend)-phi_mfl(ibeg:iend-1))
      IF(delphiplot.GT.EPSILON(1.d0)) THEN
        !! Modification by Andreas F. Martitsch (17.07.2015)
        ! fixed warning: Possible change of value in conversion
        ! from REAL(8) to INTEGER(4)
        nphiequi=INT((phi_mfl(iend)-phi_mfl(ibeg))/delphiplot)         
        !! End Modification by Andreas F. Martitsch (17.07.2015)
        nphiequi=MAX(1,nphiequi)
      ELSE
        nphiequi=1
      ENDIF
      delphiplot=(phi_mfl(iend)-phi_mfl(ibeg))/MIN(nphiplot,nphiequi)
      nplp1=npart_loc+1
      ALLOCATE(fun_write(0:lag,0:3,0:nplp1,3))
      icounter=0
      phiplot=phi_mfl(ibeg)-1.d0
!
      WRITE(propname,*) fieldpropagator%tag
      OPEN(iunit_phi,file='phi_mesh.'                        &
           //TRIM(ADJUSTL(propname))//'.dat')
      OPEN(iunit_dt_p,form='unformatted',file='dentf_p.'     &
           //TRIM(ADJUSTL(propname))//'.dat')
      OPEN(iunit_dt_m,form='unformatted',file='dentf_m.'     &
           //TRIM(ADJUSTL(propname))//'.dat')
      OPEN(iunit_sp_p,form='unformatted',file='spitf_p.'     &
           //TRIM(ADJUSTL(propname))//'.dat')
      OPEN(iunit_sp_m,form='unformatted',file='spitf_m.'     &
           //TRIM(ADJUSTL(propname))//'.dat')
      OPEN(iunit_et_p,form='unformatted',file='enetf_p.'     &
           //TRIM(ADJUSTL(propname))//'.dat')
      OPEN(iunit_et_m,form='unformatted',file='enetf_m.'     &
           //TRIM(ADJUSTL(propname))//'.dat')
!
      DO istep=ibeg,iend
        IF(phi_mfl(istep).LT.phiplot.AND.istep.NE.iend) CYCLE
        icounter=icounter+1
        phiplot=phi_mfl(istep)+delphiplot
        npassing=npl(istep)
        eta0=1.d0/bhat_mfl(istep)
        WRITE (iunit_phi,*) phi_mfl(istep),npassing,bhat_mfl(istep)
!
        fun_write=0.d0
        DO m=0,lag
          k=ind_start(istep)+2*(npassing+1)*m
          !! Modification by Andreas F. Martitsch (17.07.2015)
          ! fixed warning: Possible change of value in conversion
          ! from COMPLEX(8) to REAL(8)
          fun_write(m,:,1,:)=REAL(MATMUL(derivs_plot(:,:,1,istep),        &
                                  source_vector(k+1:k+4,:)))
          !! End Modification by Andreas F. Martitsch (17.07.2015) 
          DO i=2,npassing+1
            !! Modification by Andreas F. Martitsch (17.07.2015)
            ! fixed warning: Possible change of value in conversion
            ! from COMPLEX(8) to REAL(8)
            fun_write(m,:,i,:)=REAL(MATMUL(derivs_plot(:,:,i,istep),      &
                                    source_vector(k+i-1:k+i+2,:)))
            !! End Modification by Andreas F. Martitsch (17.07.2015) 
          ENDDO
        ENDDO
        WRITE(iunit_dt_p) fun_write(:,:,:,1)
        WRITE(iunit_sp_p) fun_write(:,:,:,2)/surface_boozer_B00
        WRITE(iunit_et_p) fun_write(:,:,:,3)
!
        fun_write=0.d0
        DO m=0,lag
          k=ind_start(istep)+2*(npassing+1)*(m+1)
          !! Modification by Andreas F. Martitsch (17.07.2015)
          ! fixed warning: Possible change of value in conversion
          ! from COMPLEX(8) to REAL(8)
          fun_write(m,:,1,:)=REAL(MATMUL(derivs_plot(:,:,1,istep),        &
                                  source_vector(k:k-3:-1,:)))
          !! End Modification by Andreas F. Martitsch (17.07.2015) 
          DO i=2,npassing+1
            !! Modification by Andreas F. Martitsch (17.07.2015)
            ! fixed warning: Possible change of value in conversion
            ! from COMPLEX(8) to REAL(8)
            fun_write(m,:,i,:)=REAL(MATMUL(derivs_plot(:,:,i,istep),      &
                                    source_vector(k-i+2:k-i-1:-1,:)))            
            !! End Modification by Andreas F. Martitsch (17.07.2015)
          ENDDO
        ENDDO
        WRITE(iunit_dt_m) fun_write(:,:,:,1)
        WRITE(iunit_sp_m) fun_write(:,:,:,2)/surface_boozer_B00
        WRITE(iunit_et_m) fun_write(:,:,:,3)
!
      ENDDO
!
      CLOSE(iunit_phi)
      CLOSE(iunit_dt_p)
      CLOSE(iunit_dt_m)
      CLOSE(iunit_sp_p)
      CLOSE(iunit_sp_m)
      CLOSE(iunit_et_p)
      CLOSE(iunit_et_m)
      OPEN(iunit_sizes,file='sizeplot_etalev.'               &
           //TRIM(ADJUSTL(propname))//'.dat')
      WRITE(iunit_sizes,*) lag,nplp1,icounter,collpar,travis_convfac
      WRITE(iunit_sizes,*) eta(0:nplp1)
      CLOSE(iunit_sizes)
!
      DEALLOCATE(fun_write)
!
    ENDIF
!
    qflux=0.5d0*REAL(MATMUL(CONJG(flux_vector),source_vector))
!
    IF(clean) THEN
!
      iopt=3
!
      CALL sparse_solve(nrow,ncol,nz,irow(1:nz),ipcol,amat_sp(1:nz),bvec_sp,iopt)
!
      IF(isw_intp.EQ.1) THEN 
        DEALLOCATE(bvec_iter,bvec_prev)
        IF(precond) DEALLOCATE(amat_z,bvec_lapack_z,ipivot,scalprod,basevec_p)
      ENDIF
!
    ENDIF
!
    END SUBROUTINE solve_eqs
!
!------------------------------------------------------------------------
!
    SUBROUTINE source_flux
!
    DO istep=ibeg,iend
      npassing=npl(istep)
      q_rip(1:npassing+1,istep,1)=q_rip_1(1:npassing+1,istep)           &
                                 *geodcu_forw(istep)
      convol_flux(1:npassing+1,istep)=convol_flux_0(1:npassing+1,istep) &
                                     *geodcu_back(istep)
    ENDDO
!
    flux_vector=0.d0
    source_vector=0.d0
!
    DO istep=ibeg,iend
!
      ioddeven=MOD(istep-ibeg,2) !0 for full RK step, 1 for half RK step
!
      IF(istep.EQ.ibeg) THEN
        step_factor_p=delt_pos(ibeg+1)*fact_pos_b(ibeg)/3.d0
        step_factor_m=delt_neg(ibeg)*fact_neg_e(ibeg)/3.d0
      ELSEIF(istep.EQ.iend) THEN
        step_factor_p=delt_pos(iend)*fact_pos_e(iend)/3.d0
        step_factor_m=delt_neg(iend-1)*fact_neg_b(iend)/3.d0
      ELSEIF(ioddeven.EQ.1) THEN
        step_factor_p=(delt_pos(istep+1)*fact_pos_b(istep)        &
                     + delt_pos(istep)*fact_pos_e(istep))/1.5d0
        step_factor_m=(delt_neg(istep-1)*fact_neg_b(istep)       &
                     + delt_neg(istep)*fact_neg_e(istep))/1.5d0
      ELSE
        step_factor_p=(delt_pos(istep+1)*fact_pos_b(istep)        &
                     + delt_pos(istep)*fact_pos_e(istep))/3.d0
        step_factor_m=(delt_neg(istep-1)*fact_neg_b(istep)       &
                     + delt_neg(istep)*fact_neg_e(istep))/3.d0
      ENDIF
!
      npassing=npl(istep)
!
      DO m=0,lag
        k=ind_start(istep)+2*(npassing+1)*m
!
        flux_vector(1,k+1:k+npassing+1) =                                     &
              step_factor_p*weightlag(1,m)*convol_flux(1:npassing+1,istep)
        flux_vector(1,k+npassing+2:k+2*npassing+2)=                           &
              step_factor_m*weightlag(1,m)*convol_flux(npassing+1:1:-1,istep)
!
        flux_vector(2,k+1:k+npassing+1) =                                     &
              step_factor_p*weightlag(2,m)*convol_curr(1:npassing+1,istep)
        flux_vector(2,k+npassing+2:k+2*npassing+2)=                           &
             -step_factor_m*weightlag(2,m)*convol_curr(npassing+1:1:-1,istep)
!
        flux_vector(3,k+1:k+npassing+1) =                                     &
              step_factor_p*weightlag(3,m)*convol_flux(1:npassing+1,istep)
        flux_vector(3,k+npassing+2:k+2*npassing+2) =                          &
              step_factor_m*weightlag(3,m)*convol_flux(npassing+1:1:-1,istep)
!
        IF(istep.GT.ibeg) THEN
          npassing_prev=npl(istep-1)
          k_prev=ind_start(istep-1)+2*(npassing_prev+1)*m
          IF(ioddeven.EQ.1) THEN
            npassing_next=npl(istep+1)
            source_vector(k+1:k+npassing+1,1)                                &
                 =source_vector(k+1:k+npassing+1,1)                          &
                 +asource(m,1)/1.5d0*q_rip(1:npassing+1,istep,1)             &
                 *fact_pos_e(istep)
            source_vector(k+1:k+npassing+1,2)                                &
                 =source_vector(k+1:k+npassing+1,2)                          &
                 +asource(m,2)/1.5d0*q_rip(1:npassing+1,istep,2)             &
                 *fact_pos_e(istep)
            source_vector(k+1:k+npassing+1,3)                                &
                 =source_vector(k+1:k+npassing+1,3)                          &
                 +asource(m,3)/1.5d0*q_rip(1:npassing+1,istep,1)             &
                 *fact_pos_e(istep)
            source_vector(k+1:k+npassing_prev+1,1)                           &
                 =source_vector(k+1:k+npassing_prev+1,1)                     &
                 +asource(m,1)/2.4d0*q_rip(1:npassing_prev+1,istep-1,1)      &
                 *fact_pos_b(istep-1)
            source_vector(k+1:k+npassing_prev+1,2)                           &
                 =source_vector(k+1:k+npassing_prev+1,2)                     &
                 +asource(m,2)/2.4d0*q_rip(1:npassing_prev+1,istep-1,2)      &
                 *fact_pos_b(istep-1)
            source_vector(k+1:k+npassing_prev+1,3)                           &
                 =source_vector(k+1:k+npassing_prev+1,3)                     &
                 +asource(m,3)/2.4d0*q_rip(1:npassing_prev+1,istep-1,1)      &
                 *fact_pos_b(istep-1)
            source_vector(k+1:k+npassing_next+1,1)                           &
                 =source_vector(k+1:k+npassing_next+1,1)                     &
                 -asource(m,1)/12d0*q_rip(1:npassing_next+1,istep+1,1)       &
                 *fact_pos_e(istep+1)
            source_vector(k+1:k+npassing_next+1,2)                           &
                 =source_vector(k+1:k+npassing_next+1,2)                     &
                 -asource(m,2)/12d0*q_rip(1:npassing_next+1,istep+1,2)       &
                 *fact_pos_e(istep+1)
            source_vector(k+1:k+npassing_next+1,3)                           &
                 =source_vector(k+1:k+npassing_next+1,3)                     &
                 -asource(m,3)/12d0*q_rip(1:npassing_next+1,istep+1,1)       &
                 *fact_pos_e(istep+1)
          ELSE
            npassing_next=npl(istep-2)
            source_vector(k+1:k+npassing+1,1)                                &
                 =source_vector(k+1:k+npassing+1,1)                          &
                 +asource(m,1)/2.4d0*q_rip(1:npassing+1,istep,1)             &
                 *fact_pos_e(istep)
            source_vector(k+1:k+npassing+1,2)                                &
                 =source_vector(k+1:k+npassing+1,2)                          &
                 +asource(m,2)/2.4d0*q_rip(1:npassing+1,istep,2)             &
                 *fact_pos_e(istep)
            source_vector(k+1:k+npassing+1,3)                                &
                 =source_vector(k+1:k+npassing+1,3)                          &
                 +asource(m,3)/2.4d0*q_rip(1:npassing+1,istep,1)             &
                 *fact_pos_e(istep)
            IF(npassing_prev.LE.npassing) THEN
              source_vector(k+1:k+npassing_prev+1,1)                         &
                   =source_vector(k+1:k+npassing_prev+1,1)                   &
                   +asource(m,1)/1.5d0*q_rip(1:npassing_prev+1,istep-1,1)    &
                   *fact_pos_b(istep-1)
              source_vector(k+1:k+npassing_prev+1,2)                         &
                   =source_vector(k+1:k+npassing_prev+1,2)                   &
                   +asource(m,2)/1.5d0*q_rip(1:npassing_prev+1,istep-1,2)    &
                   *fact_pos_b(istep-1)
              source_vector(k+1:k+npassing_prev+1,3)                         &
                   =source_vector(k+1:k+npassing_prev+1,3)                   &
                   +asource(m,3)/1.5d0*q_rip(1:npassing_prev+1,istep-1,1)    &
                   *fact_pos_b(istep-1)
            ELSE
              source_vector(k+1:k+npassing+1,1)                              &
                   =source_vector(k+1:k+npassing+1,1)                        &
                   +asource(m,1)/1.5d0*q_rip(1:npassing+1,istep-1,1)         &
                   *fact_pos_b(istep-1)
              source_vector(k_prev+npassing_prev+2,1)                        &
                     =source_vector(k_prev+npassing_prev+2,1)                &
                     +asource(m,1)/1.5d0*q_rip(npassing_prev+1,istep-1,1)    &
                     *fact_pos_b(istep-1)
              source_vector(k+1:k+npassing+1,2)                              &
                   =source_vector(k+1:k+npassing+1,2)                        &
                   +asource(m,2)/1.5d0*q_rip(1:npassing+1,istep-1,2)         &
                   *fact_pos_b(istep-1)
              source_vector(k_prev+npassing_prev+2,2)                        &
                     =source_vector(k_prev+npassing_prev+2,2)                &
                     +asource(m,2)/1.5d0*q_rip(npassing_prev+1,istep-1,2)    &
                     *fact_pos_b(istep-1)
              source_vector(k+1:k+npassing+1,3)                              &
                   =source_vector(k+1:k+npassing+1,3)                        &
                   +asource(m,3)/1.5d0*q_rip(1:npassing+1,istep-1,1)         &
                   *fact_pos_b(istep-1)
              source_vector(k_prev+npassing_prev+2,3)                        &
                     =source_vector(k_prev+npassing_prev+2,3)                &
                     +asource(m,3)/1.5d0*q_rip(npassing_prev+1,istep-1,1)    &
                     *fact_pos_b(istep-1)
            ENDIF
            IF(npassing_next.LE.npassing) THEN
              source_vector(k+1:k+npassing_next+1,1)                         &
                   =source_vector(k+1:k+npassing_next+1,1)                   &
                   -asource(m,1)/12d0*q_rip(1:npassing_next+1,istep-2,1)     &
                   *fact_pos_b(istep-2)
              source_vector(k+1:k+npassing_next+1,2)                         &
                   =source_vector(k+1:k+npassing_next+1,2)                   &
                   -asource(m,2)/12d0*q_rip(1:npassing_next+1,istep-2,2)     &
                   *fact_pos_b(istep-2)
              source_vector(k+1:k+npassing_next+1,3)                         &
                   =source_vector(k+1:k+npassing_next+1,3)                   &
                   -asource(m,3)/12d0*q_rip(1:npassing_next+1,istep-2,1)     &
                   *fact_pos_b(istep-2)
            ELSE
              source_vector(k+1:k+npassing+1,1)                              &
                   =source_vector(k+1:k+npassing+1,1)                        &
                   -asource(m,1)/12d0*q_rip(1:npassing+1,istep-2,1)          &
                   *fact_pos_b(istep-2)
              source_vector(k_prev+npassing_prev+2,1)                        &
                     =source_vector(k_prev+npassing_prev+2,1)                &
                     -asource(m,1)/12d0*q_rip(npassing_next+1,istep-2,1)     &
                     *fact_pos_b(istep-2)
              source_vector(k+1:k+npassing+1,2)                              &
                   =source_vector(k+1:k+npassing+1,2)                        &
                   -asource(m,2)/12d0*q_rip(1:npassing+1,istep-2,2)          &
                   *fact_pos_b(istep-2)
              source_vector(k_prev+npassing_prev+2,2)                        &
                     =source_vector(k_prev+npassing_prev+2,2)                &
                     -asource(m,2)/12d0*q_rip(npassing_next+1,istep-2,2)     &
                     *fact_pos_b(istep-2)
              source_vector(k+1:k+npassing+1,3)                              &
                   =source_vector(k+1:k+npassing+1,3)                        &
                   -asource(m,3)/12d0*q_rip(1:npassing+1,istep-2,1)          &
                   *fact_pos_b(istep-2)
              source_vector(k_prev+npassing_prev+2,3)                        &
                     =source_vector(k_prev+npassing_prev+2,3)                &
                     -asource(m,3)/12d0*q_rip(npassing_next+1,istep-2,1)     &
                     *fact_pos_b(istep-2)
            ENDIF
          ENDIF
        ELSEIF(iplot.EQ.1.AND.isw_axisymm.NE.1) THEN
          source_vector(k+1:k+npassing+1,:)                                &
                         =flux_pl(npass_l*m+1:npass_l*(m+1),:)
        ENDIF
!
        IF(istep.LT.iend) THEN
          npassing_prev=npl(istep+1)
          k_prev=ind_start(istep+1)+2*(npassing_prev+1)*m
          IF(ioddeven.EQ.1) THEN
            npassing_next=npl(istep-1)
            source_vector(k+npassing+2:k+2*npassing+2,1)                     &
              =source_vector(k+npassing+2:k+2*npassing+2,1)                  &
              +asource(m,1)/1.5d0*q_rip(npassing+1:1:-1,istep,1)             &
              *fact_neg_e(istep)
            source_vector(k+npassing+2:k+2*npassing+2,2)                     &
              =source_vector(k+npassing+2:k+2*npassing+2,2)                  &
              -asource(m,2)/1.5d0*q_rip(npassing+1:1:-1,istep,2)             &
              *fact_neg_e(istep)
            source_vector(k+npassing+2:k+2*npassing+2,3)                     &
              =source_vector(k+npassing+2:k+2*npassing+2,3)                  &
              +asource(m,3)/1.5d0*q_rip(npassing+1:1:-1,istep,1)             &
              *fact_neg_e(istep)
            source_vector(k+2*npassing+2-npassing_prev:k+2*npassing+2,1)     &
              =source_vector(k+2*npassing+2-npassing_prev:k+2*npassing+2,1)  &
              +asource(m,1)/2.4d0*q_rip(npassing_prev+1:1:-1,istep+1,1)      &
              *fact_neg_b(istep+1)
            source_vector(k+2*npassing+2-npassing_prev:k+2*npassing+2,2)     &
              =source_vector(k+2*npassing+2-npassing_prev:k+2*npassing+2,2)  &
              -asource(m,2)/2.4d0*q_rip(npassing_prev+1:1:-1,istep+1,2)      &
              *fact_neg_b(istep+1)
            source_vector(k+2*npassing+2-npassing_prev:k+2*npassing+2,3)     &
              =source_vector(k+2*npassing+2-npassing_prev:k+2*npassing+2,3)  &
              +asource(m,3)/2.4d0*q_rip(npassing_prev+1:1:-1,istep+1,1)      &
              *fact_neg_b(istep+1)
            source_vector(k+2*npassing+2-npassing_next:k+2*npassing+2,1)     &
              =source_vector(k+2*npassing+2-npassing_next:k+2*npassing+2,1)  &
              -asource(m,1)/12d0*q_rip(npassing_next+1:1:-1,istep-1,1)       &
              *fact_neg_e(istep-1)
            source_vector(k+2*npassing+2-npassing_next:k+2*npassing+2,2)     &
              =source_vector(k+2*npassing+2-npassing_next:k+2*npassing+2,2)  &
              +asource(m,2)/12d0*q_rip(npassing_next+1:1:-1,istep-1,2)       &
              *fact_neg_e(istep-1)
            source_vector(k+2*npassing+2-npassing_next:k+2*npassing+2,3)     &
              =source_vector(k+2*npassing+2-npassing_next:k+2*npassing+2,3)  &
              -asource(m,3)/12d0*q_rip(npassing_next+1:1:-1,istep-1,1)       &
              *fact_neg_e(istep-1)
          ELSE
            npassing_next=npl(istep+2)
            source_vector(k+npassing+2:k+2*npassing+2,1)                     &
              =source_vector(k+npassing+2:k+2*npassing+2,1)                  &
              +asource(m,1)/2.4d0*q_rip(npassing+1:1:-1,istep,1)             &
              *fact_neg_e(istep)
            source_vector(k+npassing+2:k+2*npassing+2,2)                     &
              =source_vector(k+npassing+2:k+2*npassing+2,2)                  &
              -asource(m,2)/2.4d0*q_rip(npassing+1:1:-1,istep,2)             &
              *fact_neg_e(istep)
            source_vector(k+npassing+2:k+2*npassing+2,3)                     &
              =source_vector(k+npassing+2:k+2*npassing+2,3)                  &
              +asource(m,3)/2.4d0*q_rip(npassing+1:1:-1,istep,1)             &
              *fact_neg_e(istep)
            IF(npassing_prev.LE.npassing) THEN
              source_vector(k+2*npassing+2-npassing_prev:k+2*npassing+2,1)   &
                =source_vector(k+2*npassing+2-npassing_prev:k+2*npassing+2,1)&
                +asource(m,1)/1.5d0*q_rip(npassing_prev+1:1:-1,istep+1,1)    &
                *fact_neg_b(istep+1)
              source_vector(k+2*npassing+2-npassing_prev:k+2*npassing+2,2)   &
                =source_vector(k+2*npassing+2-npassing_prev:k+2*npassing+2,2)&
                -asource(m,2)/1.5d0*q_rip(npassing_prev+1:1:-1,istep+1,2)    &
                *fact_neg_b(istep+1)
              source_vector(k+2*npassing+2-npassing_prev:k+2*npassing+2,3)   &
                =source_vector(k+2*npassing+2-npassing_prev:k+2*npassing+2,3)&
                +asource(m,3)/1.5d0*q_rip(npassing_prev+1:1:-1,istep+1,1)    &
                *fact_neg_b(istep+1)
            ELSE
              source_vector(k+npassing+2:k+2*npassing+2,1)                   &
                =source_vector(k+npassing+2:k+2*npassing+2,1)                &
                +asource(m,1)/1.5d0*q_rip(npassing+1:1:-1,istep+1,1)         &
                *fact_neg_b(istep+1)
              source_vector(k_prev+npassing_prev+1,1)                        &
                  =source_vector(k_prev+npassing_prev+1,1)                   &
                  +asource(m,1)/1.5d0*q_rip(npassing_prev+1,istep+1,1)       &
                  *fact_neg_b(istep+1)
              source_vector(k+npassing+2:k+2*npassing+2,2)                   &
                =source_vector(k+npassing+2:k+2*npassing+2,2)                &
                -asource(m,2)/1.5d0*q_rip(npassing+1:1:-1,istep+1,2)         &
                *fact_neg_b(istep+1)
              source_vector(k_prev+npassing_prev+1,2)                        &
                  =source_vector(k_prev+npassing_prev+1,2)                   &
                  -asource(m,2)/1.5d0*q_rip(npassing_prev+1,istep+1,2)       &
                  *fact_neg_b(istep+1)
              source_vector(k+npassing+2:k+2*npassing+2,3)                   &
                =source_vector(k+npassing+2:k+2*npassing+2,3)                &
                +asource(m,3)/1.5d0*q_rip(npassing+1:1:-1,istep+1,1)         &
                *fact_neg_b(istep+1)
              source_vector(k_prev+npassing_prev+1,3)                        &
                  =source_vector(k_prev+npassing_prev+1,3)                   &
                  +asource(m,3)/1.5d0*q_rip(npassing_prev+1,istep+1,1)       &
                  *fact_neg_b(istep+1)
            ENDIF
            IF(npassing_next.LE.npassing) THEN
              source_vector(k+2*npassing+2-npassing_next:k+2*npassing+2,1)   &
                =source_vector(k+2*npassing+2-npassing_next:k+2*npassing+2,1)&
                -asource(m,1)/12d0*q_rip(npassing_next+1:1:-1,istep+2,1)     &
                *fact_neg_b(istep+2)
              source_vector(k+2*npassing+2-npassing_next:k+2*npassing+2,2)   &
                =source_vector(k+2*npassing+2-npassing_next:k+2*npassing+2,2)&
                +asource(m,2)/12d0*q_rip(npassing_next+1:1:-1,istep+2,2)     &
                *fact_neg_b(istep+2)
              source_vector(k+2*npassing+2-npassing_next:k+2*npassing+2,3)   &
                =source_vector(k+2*npassing+2-npassing_next:k+2*npassing+2,3)&
                -asource(m,3)/12d0*q_rip(npassing_next+1:1:-1,istep+2,1)     &
                *fact_neg_b(istep+2)
            ELSE
              source_vector(k+npassing+2:k+2*npassing+2,1)                   &
                =source_vector(k+npassing+2:k+2*npassing+2,1)                &
                -asource(m,1)/12d0*q_rip(npassing+1:1:-1,istep+2,1)          &
                *fact_neg_b(istep+2)
              source_vector(k_prev+npassing_prev+1,1)                        &
                  =source_vector(k_prev+npassing_prev+1,1)                   &
                  -asource(m,1)/12d0*q_rip(npassing_next+1,istep+2,1)        &
                  *fact_neg_b(istep+2)
              source_vector(k+npassing+2:k+2*npassing+2,2)                   &
                =source_vector(k+npassing+2:k+2*npassing+2,2)                &
                +asource(m,2)/12d0*q_rip(npassing+1:1:-1,istep+2,2)          &
                *fact_neg_b(istep+2)
              source_vector(k_prev+npassing_prev+1,2)                        &
                  =source_vector(k_prev+npassing_prev+1,2)                   &
                  +asource(m,2)/12d0*q_rip(npassing_next+1,istep+2,2)        &
                  *fact_neg_b(istep+2)
              source_vector(k+npassing+2:k+2*npassing+2,3)                   &
                =source_vector(k+npassing+2:k+2*npassing+2,3)                &
                -asource(m,3)/12d0*q_rip(npassing+1:1:-1,istep+2,1)          &
                *fact_neg_b(istep+2)
              source_vector(k_prev+npassing_prev+1,3)                        &
                  =source_vector(k_prev+npassing_prev+1,3)                   &
                  -asource(m,3)/12d0*q_rip(npassing_next+1,istep+2,1)        &
                  *fact_neg_b(istep+2)
            ENDIF
          ENDIF
        ELSEIF(iplot.EQ.1.AND.isw_axisymm.NE.1) THEN
          source_vector(k+npassing+2:k+2*npassing+2,:)                     &
                         =flux_mr(npass_r*(m+1):npass_r*m+1:-1,:)
        ENDIF
!
      ENDDO
    ENDDO
!
    END SUBROUTINE source_flux
!
!------------------------------------------------------------------------
!
    SUBROUTINE add_f01_source
!
    INTEGER :: k_next
!
    DO istep=ibeg,iend
!
      ioddeven=MOD(istep-ibeg,2) !0 for full RK step, 1 for half RK step
!
      npassing=npl(istep)
!
      DO m=0,lag
        k=ind_start(istep)+2*(npassing+1)*m
!
        IF(istep.GT.ibeg) THEN
          npassing_prev=npl(istep-1)
          k_prev=ind_start(istep-1)+2*(npassing_prev+1)*m
          IF(ioddeven.EQ.1) THEN
            npassing_next=npl(istep+1)
            k_next=ind_start(istep+1)+2*(npassing_next+1)*m
            source_vector(k+1:k+npassing+1,:)                                &
              =source_vector(k+1:k+npassing+1,:)                             &
              +(bnoverb0(istep)*f0_coll(k+1:k+npassing+1,:)                  &
              + ttmpfact(istep)*f0_ttmp(k+1:k+npassing+1,:))                 &
              *fact_pos_e(istep)/1.5d0
            source_vector(k+1:k+npassing_prev+1,:)                           &
              =source_vector(k+1:k+npassing_prev+1,:)                        &
              +(bnoverb0(istep-1)*f0_coll(k_prev+1:k_prev+npassing_prev+1,:) &
              + ttmpfact(istep-1)*f0_ttmp(k_prev+1:k_prev+npassing_prev+1,:))&
              *fact_pos_e(istep-1)/2.4d0
            source_vector(k+1:k+npassing_next+1,:)                           &
              =source_vector(k+1:k+npassing_next+1,:)                        &
              -(bnoverb0(istep+1)*f0_coll(k_next+1:k_next+npassing_next+1,:) &
              + ttmpfact(istep+1)*f0_ttmp(k_next+1:k_next+npassing_next+1,:))&
              *fact_pos_e(istep+1)/12d0
          ELSE
            npassing_next=npl(istep-2)
            k_next=ind_start(istep-2)+2*(npassing_next+1)*m
            source_vector(k+1:k+npassing+1,:)                                &
              =source_vector(k+1:k+npassing+1,:)                             &
              +(bnoverb0(istep)*f0_coll(k+1:k+npassing+1,:)                  &
              + ttmpfact(istep)*f0_ttmp(k+1:k+npassing+1,:))                 &
              *fact_pos_e(istep)/2.4d0
            IF(npassing_prev.LE.npassing) THEN
              source_vector(k+1:k+npassing_prev+1,:)                           &
                =source_vector(k+1:k+npassing_prev+1,:)                        &
                +(bnoverb0(istep-1)*f0_coll(k_prev+1:k_prev+npassing_prev+1,:) &
                + ttmpfact(istep-1)*f0_ttmp(k_prev+1:k_prev+npassing_prev+1,:))&
                *fact_pos_b(istep-1)/1.5d0
            ELSE
              source_vector(k+1:k+npassing+1,:)                                &
                =source_vector(k+1:k+npassing+1,:)                             &
                +(bnoverb0(istep-1)*f0_coll(k_prev+1:k_prev+npassing+1,:)      &
                + ttmpfact(istep-1)*f0_ttmp(k_prev+1:k_prev+npassing+1,:))     &
                *fact_pos_b(istep-1)/1.5d0
              source_vector(k_prev+npassing_prev+2,:)                          &
                =source_vector(k_prev+npassing_prev+2,:)                       &
                +(bnoverb0(istep-1)*f0_coll(k_prev+npassing_prev+1,:)          &
                + ttmpfact(istep-1)*f0_ttmp(k_prev+npassing_prev+1,:))         &
                *fact_pos_b(istep-1)/1.5d0
            ENDIF
            IF(npassing_next.LE.npassing) THEN
              source_vector(k+1:k+npassing_next+1,:)                           &
                =source_vector(k+1:k+npassing_next+1,:)                        &
                -(bnoverb0(istep-2)*f0_coll(k_next+1:k_next+npassing_next+1,:) &
                + ttmpfact(istep-2)*f0_ttmp(k_next+1:k_next+npassing_next+1,:))&
                *fact_pos_b(istep-2)/12d0
            ELSE
              source_vector(k+1:k+npassing+1,:)                                &
                =source_vector(k+1:k+npassing+1,:)                             &
                -(bnoverb0(istep-2)*f0_coll(k_next+1:k_next+npassing+1,:)      &
                + ttmpfact(istep-2)*f0_ttmp(k_next+1:k_next+npassing+1,:))     &
                *fact_pos_b(istep-2)/12d0
              source_vector(k_prev+npassing_prev+2,:)                          &
                =source_vector(k_prev+npassing_prev+2,:)                       &
                -(bnoverb0(istep-2)*f0_coll(k_next+npassing_next+1,:)          &
                + ttmpfact(istep-2)*f0_ttmp(k_next+npassing_next+1,:))         &
                *fact_pos_b(istep-2)/12d0
            ENDIF
          ENDIF
        ENDIF
!
        IF(istep.LT.iend) THEN
          npassing_prev=npl(istep+1)
          k_prev=ind_start(istep+1)+2*(npassing_prev+1)*m
          IF(ioddeven.EQ.1) THEN
            npassing_next=npl(istep-1)
            k_next=ind_start(istep-1)+2*(npassing_next+1)*m
            source_vector(k+npassing+2:k+2*npassing+2,:)                     &
              =source_vector(k+npassing+2:k+2*npassing+2,:)                  &
              +(bnoverb0(istep)*f0_coll(k+npassing+2:k+2*npassing+2,:)       &
              + ttmpfact(istep)*f0_ttmp(k+npassing+2:k+2*npassing+2,:))      &
              *fact_neg_e(istep)/1.5d0
            source_vector(k+2*npassing+2-npassing_prev:k+2*npassing+2,:)     &
              =source_vector(k+2*npassing+2-npassing_prev:k+2*npassing+2,:)  &
              +(bnoverb0(istep+1)*f0_coll(k_prev+npassing_prev+2:            &
                                          k_prev+2*npassing_prev+2,:)        &
              + ttmpfact(istep+1)*f0_ttmp(k_prev+npassing_prev+2:            &
                                          k_prev+2*npassing_prev+2,:))       &
              *fact_neg_b(istep+1)/2.4d0
            source_vector(k+2*npassing+2-npassing_next:k+2*npassing+2,:)     &
              =source_vector(k+2*npassing+2-npassing_next:k+2*npassing+2,:)  &
              -(bnoverb0(istep-1)*f0_coll(k_next+npassing_next+2:            &
                                          k_next+2*npassing_next+2,:)        &
              + ttmpfact(istep-1)*f0_ttmp(k_next+npassing_next+2:            &
                                          k_next+2*npassing_next+2,:))       &
              *fact_neg_e(istep-1)/12d0
          ELSE
            npassing_next=npl(istep+2)
            k_next=ind_start(istep+2)+2*(npassing_next+1)*m
            source_vector(k+npassing+2:k+2*npassing+2,:)                     &
              =source_vector(k+npassing+2:k+2*npassing+2,:)                  &
              +(bnoverb0(istep)*f0_coll(k+npassing+2:k+2*npassing+2,:)       &
              + ttmpfact(istep)*f0_ttmp(k+npassing+2:k+2*npassing+2,:))      &
              *fact_neg_e(istep)/2.4d0
            IF(npassing_prev.LE.npassing) THEN
              source_vector(k+2*npassing+2-npassing_prev:k+2*npassing+2,:)   &
                =source_vector(k+2*npassing+2-npassing_prev:k+2*npassing+2,:)&
                +(bnoverb0(istep+1)*f0_coll(k_prev+npassing_prev+2:          &
                                            k_prev+2*npassing_prev+2,:)      &
                + ttmpfact(istep+1)*f0_ttmp(k_prev+npassing_prev+2:          &
                                            k_prev+2*npassing_prev+2,:))     &
                *fact_neg_b(istep+1)/1.5d0
            ELSE
              source_vector(k+npassing+2:k+2*npassing+2,:)                   &
                =source_vector(k+npassing+2:k+2*npassing+2,:)                &
                +(bnoverb0(istep+1)*f0_coll(k_prev+2*npassing_prev+2-npassing: &
                                            k_prev+2*npassing_prev+2,:)      &
                + ttmpfact(istep+1)*f0_ttmp(k_prev+2*npassing_prev+2-npassing: &
                                            k_prev+2*npassing_prev+2,:))     &
                *fact_neg_b(istep+1)/1.5d0
              source_vector(k_prev+npassing_prev+1,:)                        &
                  =source_vector(k_prev+npassing_prev+1,:)                   &
                  +(bnoverb0(istep+1)*f0_coll(k_prev+npassing_prev+2,:)      &
                  + ttmpfact(istep+1)*f0_ttmp(k_prev+npassing_prev+2,:))     &
                  *fact_neg_b(istep+1)/1.5d0
            ENDIF
            IF(npassing_next.LE.npassing) THEN
              source_vector(k+2*npassing+2-npassing_next:k+2*npassing+2,:)   &
                =source_vector(k+2*npassing+2-npassing_next:k+2*npassing+2,:)&
                -(bnoverb0(istep+2)*f0_coll(k_next+npassing_next+2:          &
                                            k_next+2*npassing_next+2,:)      &
                + ttmpfact(istep+2)*f0_ttmp(k_next+npassing_next+2:          &
                                            k_next+2*npassing_next+2,:))     &
                *fact_neg_b(istep+2)/12d0
            ELSE
              source_vector(k+npassing+2:k+2*npassing+2,:)                   &
                =source_vector(k+npassing+2:k+2*npassing+2,:)                &
                -(bnoverb0(istep+2)*f0_coll(k_next+2*npassing_next+2-npassing: &
                                            k_next+2*npassing_next+2,:)      &
                + ttmpfact(istep+2)*f0_ttmp(k_next+2*npassing_next+2-npassing: &
                                            k_next+2*npassing_next+2,:))     &
                *fact_neg_b(istep+2)/12d0
              source_vector(k_prev+npassing_prev+1,:)                        &
                  =source_vector(k_prev+npassing_prev+1,:)                   &
                  -(bnoverb0(istep+2)*f0_coll(k_next+npassing_next+2,:)      &
                  + ttmpfact(istep+2)*f0_ttmp(k_next+npassing_next+2,:))     &
                  *fact_neg_b(istep+2)/12d0
            ENDIF
          ENDIF
        ENDIF
!
      ENDDO
    ENDDO
!
    END SUBROUTINE add_f01_source
!
!------------------------------------------------------------------------
!
    SUBROUTINE iterate(vec_inout)
!
    IMPLICIT NONE
!
    DOUBLE COMPLEX, DIMENSION(n_2d_size) :: vec_inout
!
    bvec_prev=vec_inout
!
    DO iter=1,niter
!
      CALL integral_part(bvec_prev,bvec_iter)
!
      CALL sparse_solve(nrow,ncol,nz,irow(1:nz),ipcol,amat_sp(1:nz),         &
                        bvec_iter,iopt)
!
      IF(precond) THEN
!
        CALL preconditioner( .FALSE. ,bvec_iter,bvec_prev)
!
      ELSE
        bvec_prev=bvec_iter
      ENDIF
!
      vec_inout=vec_inout+bvec_prev
      IF(SUM(ABS(bvec_prev)).LT.SUM(ABS(vec_inout))*epserr_iter) THEN
        PRINT *,'number of iterations = ',iter
!
        IF(.NOT.precond) EXIT
!
        CALL integral_part(vec_inout,bvec_iter)
!
        CALL sparse_solve(nrow,ncol,nz,irow(1:nz),ipcol,amat_sp(1:nz),       &
                        bvec_iter,iopt)
!
        CALL preconditioner( .TRUE. ,bvec_iter,bvec_prev)
!
        EXIT
      ENDIF
PRINT *,iter,SUM(ABS(bvec_prev)),SUM(ABS(vec_inout))*epserr_iter
    ENDDO
!
    END SUBROUTINE iterate
!
!------------------------------------------------------------------------
!
    SUBROUTINE preconditioner(nochange,vec_in,vec_out)
!
    IMPLICIT NONE
!
    LOGICAL :: nochange
!
    INTEGER :: m,k,istep,npassing
!
    DOUBLE COMPLEX, DIMENSION(n_2d_size)           :: vec_in,vec_out
    DOUBLE COMPLEX, DIMENSION(:,:), ALLOCATABLE    :: scalpvec
!
    ALLOCATE(scalpvec(0:lag,ibeg:iend))
!
    DO istep=ibeg,iend
!
      npassing=npl(istep)
!
      DO m=0,lag
        k=ind_start(istep)+2*(npassing+1)*m
        scalpvec(m,istep)=SUM(vec_in(k+1:k+npassing+1))
        k=k+2*(npassing+1)
        scalpvec(m,istep)=scalpvec(m,istep)-SUM(vec_in(k-npassing:k))
      ENDDO
!
    ENDDO
!
    DO m=0,lag
      scalprod(m)=SUM((scalpvec(m,ibeg:iend-1)+scalpvec(m,ibeg+1:iend))     &
                 *    (phi_mfl(ibeg+1:iend)-phi_mfl(ibeg:iend-1)))
    ENDDO
!
    DEALLOCATE(scalpvec)
!
    scalprod=0.5d0*scalprod/(phi_mfl(iend)-phi_mfl(ibeg))
!
    IF(nochange) RETURN
!
    DO istep=ibeg,iend
!
      npassing=npl(istep)
!
      DO m=0,lag
        k=ind_start(istep)+2*(npassing+1)*m
        vec_out(k+1:k+npassing+1)=vec_in(k+1:k+npassing+1)-scalprod(m)     &
                                 *(0.5d0*pleg_ket(2,1:npassing+1,istep)    &
                                 +0.25d0*pleg_ket(0,1:npassing+1,istep))   &
                                 *bhat_mfl(istep)
        k=k+2*(npassing+1)
        vec_out(k:k-npassing:-1) =vec_in(k:k-npassing:-1)+scalprod(m)      &
                                 *(0.5d0*pleg_ket(2,1:npassing+1,istep)    &
                                 +0.25d0*pleg_ket(0,1:npassing+1,istep))   &
                                 *bhat_mfl(istep)
      ENDDO
!
    ENDDO
!
    END SUBROUTINE preconditioner
!
!------------------------------------------------------------------------
!
    SUBROUTINE basevec_precond(m,vec_out)
!
    IMPLICIT NONE
!
    INTEGER :: m,k,istep,npassing
!
    DOUBLE COMPLEX, DIMENSION(n_2d_size)           :: vec_out
!
    vec_out=0.d0
!
    DO istep=ibeg,iend
!
      npassing=npl(istep)
      k=ind_start(istep)+2*(npassing+1)*m
      vec_out(k+1:k+npassing+1)= (0.5d0*pleg_ket(2,1:npassing+1,istep)    &
                               + 0.25d0*pleg_ket(0,1:npassing+1,istep))   &
                               * bhat_mfl(istep)
      k=k+2*(npassing+1)
      vec_out(k:k-npassing:-1) =-(0.5d0*pleg_ket(2,1:npassing+1,istep)    &
                               + 0.25d0*pleg_ket(0,1:npassing+1,istep))   &
                               * bhat_mfl(istep)
    ENDDO
!
    END SUBROUTINE basevec_precond
!
!-------------------------------------------------------------------------
!
    SUBROUTINE integral_part(vec_in,vec_out)
!
    IMPLICIT NONE
!
    INTEGER :: l,m,i,k,istep,npassing,kk
!
    DOUBLE COMPLEX, DIMENSION(n_2d_size)                 :: vec_in,vec_out
    DOUBLE COMPLEX, DIMENSION(0:lag,0:leg)               :: scalprod_pleg
    DOUBLE COMPLEX, DIMENSION(:,:,:), ALLOCATABLE        :: vec_tmp
!
    ALLOCATE(vec_tmp(0:lag,2*(npart+1),ibeg:iend))
    vec_tmp=0.d0
!
    DO istep=ibeg,iend
!
      npassing=npl(istep)
!
IF(.FALSE.) THEN
!if(istep.eq.ibeg) then
DO m=0,lag
DO i=1,npassing
k=ind_start(istep)+i+2*(npassing+1)*m
WRITE(2001,*) vec_in(k)/(eta(i)-eta(i-1))
ENDDO
i=npassing+1
k=ind_start(istep)+i+2*(npassing+1)*m
WRITE(2001,*) vec_in(k)/(1.d0/bhat_mfl(istep)-eta(i-1))
k=ind_start(istep)+i+1+2*(npassing+1)*m
WRITE(2001,*) vec_in(k)/(1.d0/bhat_mfl(istep)-eta(i-1))
DO i=npassing,1,-1
k=ind_start(istep)+2*npassing+3-i+2*(npassing+1)*m
WRITE(2001,*) vec_in(k)/(eta(i)-eta(i-1))
ENDDO
ENDDO
CLOSE(2001)
ENDIF
IF(.FALSE.) THEN
!if(istep.eq.(ibeg+iend)/2) then
DO m=0,lag
DO i=1,npassing
k=ind_start(istep)+i+2*(npassing+1)*m
WRITE(2002,*) vec_in(k)/(eta(i)-eta(i-1))
ENDDO
i=npassing+1
k=ind_start(istep)+i+2*(npassing+1)*m
WRITE(2002,*) vec_in(k)/(1.d0/bhat_mfl(istep)-eta(i-1))
k=ind_start(istep)+i+1+2*(npassing+1)*m
WRITE(2002,*) vec_in(k)/(1.d0/bhat_mfl(istep)-eta(i-1))
DO i=npassing,1,-1
k=ind_start(istep)+2*npassing+3-i+2*(npassing+1)*m
WRITE(2002,*) vec_in(k)/(eta(i)-eta(i-1))
ENDDO
ENDDO
CLOSE(2002)
!pause
ENDIF
!
      DO m=0,lag
        k=ind_start(istep)+2*(npassing+1)*m
        DO l=0,leg
          scalprod_pleg(m,l)=                                             &
              SUM(pleg_bra(l,1:npassing+1,istep)*vec_in(k+1:k+npassing+1))
        ENDDO
        k=k+2*(npassing+1)
        DO l=0,leg,2
          scalprod_pleg(m,l)=scalprod_pleg(m,l)                           &
             +SUM(pleg_bra(l,1:npassing+1,istep)*vec_in(k:k-npassing:-1))
        ENDDO
        DO l=1,leg,2
          scalprod_pleg(m,l)=scalprod_pleg(m,l)                           &
             -SUM(pleg_bra(l,1:npassing+1,istep)*vec_in(k:k-npassing:-1))
        ENDDO
      ENDDO
!
      DO l=0,leg
        scalprod_pleg(0:lag,l)=MATMUL(ailmm(0:lag,0:lag,l),               &
                                      scalprod_pleg(0:lag,l))
      ENDDO
!
      DO m=0,lag
!
        DO l=0,leg
          vec_tmp(m,1:npassing+1,istep)=vec_tmp(m,1:npassing+1,istep)   &
                      +scalprod_pleg(m,l)*pleg_ket(l,1:npassing+1,istep)
        ENDDO
!
        k=2*(npassing+1)
!
        DO l=0,leg,2
          vec_tmp(m,k:k-npassing:-1,istep)=vec_tmp(m,k:k-npassing:-1,istep) &
                      +scalprod_pleg(m,l)*pleg_ket(l,1:npassing+1,istep)
        ENDDO
        DO l=1,leg,2
          vec_tmp(m,k:k-npassing:-1,istep)=vec_tmp(m,k:k-npassing:-1,istep) &
                      -scalprod_pleg(m,l)*pleg_ket(l,1:npassing+1,istep)
        ENDDO
!
      ENDDO
!
    ENDDO
!
    vec_out=0.d0
!
    DO istep=ibeg,iend
!
      ioddeven=MOD(istep-ibeg,2) !0 for full RK step, 1 for half RK step
!
      npassing=npl(istep)
!
      DO m=0,lag
        k=ind_start(istep)+2*(npassing+1)*m
!
        IF(istep.GT.ibeg) THEN
          npassing_prev=npl(istep-1)
          k_prev=ind_start(istep-1)+2*(npassing_prev+1)*m
          IF(ioddeven.EQ.1) THEN
            npassing_next=npl(istep+1)
            vec_out(k+1:k+npassing+1)=vec_out(k+1:k+npassing+1)               &
                                     +vec_tmp(m,1:npassing+1,istep)           &
                                     *fact_pos_e(istep)/1.5d0
            vec_out(k+1:k+npassing_prev+1)=vec_out(k+1:k+npassing_prev+1)     &
                                     +vec_tmp(m,1:npassing_prev+1,istep-1)    &
                                     *fact_pos_b(istep-1)/2.4d0
            vec_out(k+1:k+npassing_next+1)=vec_out(k+1:k+npassing_next+1)     &
                                     -vec_tmp(m,1:npassing_next+1,istep+1)    &
                                     *fact_pos_e(istep+1)/12d0
          ELSE
            npassing_next=npl(istep-2)
            vec_out(k+1:k+npassing+1)=vec_out(k+1:k+npassing+1)               &
                                     +vec_tmp(m,1:npassing+1,istep)           &
                                     *fact_pos_e(istep)/2.4d0
            IF(npassing_prev.LE.npassing) THEN
              vec_out(k+1:k+npassing_prev+1)=vec_out(k+1:k+npassing_prev+1)   &
                                     +vec_tmp(m,1:npassing_prev+1,istep-1)    &
                                     *fact_pos_b(istep-1)/1.5d0
            ELSE
              vec_out(k+1:k+npassing+1)=vec_out(k+1:k+npassing+1)             &
                                       +vec_tmp(m,1:npassing+1,istep-1)       &
                                       *fact_pos_b(istep-1)/1.5d0
              vec_out(k_prev+npassing_prev+2)=vec_out(k_prev+npassing_prev+2) &
                                       +vec_tmp(m,npassing_prev+1,istep-1)    &
                                       *fact_pos_b(istep-1)/1.5d0
            ENDIF
            IF(npassing_next.LE.npassing) THEN
              vec_out(k+1:k+npassing_next+1)=vec_out(k+1:k+npassing_next+1)   &
                                       -vec_tmp(m,1:npassing_next+1,istep-2)  &
                                       *fact_pos_b(istep-2)/12d0
            ELSE
              vec_out(k+1:k+npassing+1)=vec_out(k+1:k+npassing+1)             &
                                       -vec_tmp(m,1:npassing+1,istep-2)       &
                                       *fact_pos_b(istep-2)/12d0
              vec_out(k_prev+npassing_prev+2)=vec_out(k_prev+npassing_prev+2) &
                                       -vec_tmp(m,npassing_next+1,istep-2)    &
                                        *fact_pos_b(istep-2)/12d0
            ENDIF
          ENDIF
        ENDIF
!
        IF(istep.LT.iend) THEN
          npassing_prev=npl(istep+1)
          k_prev=ind_start(istep+1)+2*(npassing_prev+1)*m
          IF(ioddeven.EQ.1) THEN
            npassing_next=npl(istep-1)
            vec_out(k+npassing+2:k+2*npassing+2)                              &
                 =vec_out(k+npassing+2:k+2*npassing+2)                        &
                 +vec_tmp(m,npassing+2:2*npassing+2,istep)                    &
                 *fact_neg_e(istep)/1.5d0
            vec_out(k+2*npassing+2-npassing_prev:k+2*npassing+2)              &
                 =vec_out(k+2*npassing+2-npassing_prev:k+2*npassing+2)        &
                 +vec_tmp(m,npassing_prev+2:2*npassing_prev+2,istep+1)        &
                 *fact_neg_b(istep+1)/2.4d0
            vec_out(k+2*npassing+2-npassing_next:k+2*npassing+2)              &
                 =vec_out(k+2*npassing+2-npassing_next:k+2*npassing+2)        &
                 -vec_tmp(m,npassing_next+2:2*npassing_next+2,istep-1)        &
                 *fact_neg_e(istep-1)/12d0
          ELSE
            npassing_next=npl(istep+2)
            vec_out(k+npassing+2:k+2*npassing+2)                              &
                 =vec_out(k+npassing+2:k+2*npassing+2)                        &
                 +vec_tmp(m,npassing+2:2*npassing+2,istep)                    &
                 *fact_neg_e(istep)/2.4d0
            IF(npassing_prev.LE.npassing) THEN
              vec_out(k+2*npassing+2-npassing_prev:k+2*npassing+2)            &
                 =vec_out(k+2*npassing+2-npassing_prev:k+2*npassing+2)        &
                 +vec_tmp(m,npassing_prev+2:2*npassing_prev+2,istep+1)        &
                 *fact_neg_b(istep+1)/1.5d0
            ELSE
              vec_out(k+npassing+2:k+2*npassing+2)                            &
                 =vec_out(k+npassing+2:k+2*npassing+2)                        &
                 +vec_tmp(m,2*npassing_prev+2-npassing:2*npassing_prev+2,     &
                                                             istep+1)         &
                 *fact_neg_b(istep+1)/1.5d0
              vec_out(k_prev+npassing_prev+1)                                 &
                 =vec_out(k_prev+npassing_prev+1)                             &
                 +vec_tmp(m,npassing_prev+2,istep+1)                          &
                 *fact_neg_b(istep+1)/1.5d0
            ENDIF
            IF(npassing_next.LE.npassing) THEN
              vec_out(k+2*npassing+2-npassing_next:k+2*npassing+2)            &
                 =vec_out(k+2*npassing+2-npassing_next:k+2*npassing+2)        &
                 -vec_tmp(m,npassing_next+2:2*npassing_next+2,istep+2)        &
                 *fact_neg_b(istep+2)/12d0
            ELSE
              vec_out(k+npassing+2:k+2*npassing+2)                            &
                 =vec_out(k+npassing+2:k+2*npassing+2)                        &
                 -vec_tmp(m,2*npassing_next+2-npassing:2*npassing_next+2,     &
                                                             istep+2)         &
                 *fact_neg_b(istep+2)/12d0
              vec_out(k_prev+npassing_prev+1)                                 &
                 =vec_out(k_prev+npassing_prev+1)                             &
                 -vec_tmp(m,npassing_next+2,istep+2)                          &
                 *fact_neg_b(istep+2)/12d0
            ENDIF
          ENDIF
        ENDIF
      ENDDO
!
    ENDDO
!
    DEALLOCATE(vec_tmp)
!
    END SUBROUTINE integral_part
!
END SUBROUTINE ripple_solver
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
SUBROUTINE rearrange_phideps(ibeg,iend,npart,ncomp,nreal,subsqmin,phi_divide, &
                             phi_mfl,bhat_mfl,arr_real,arr_comp,eta,          &
                             delt_pos,delt_neg,                               &
                             fact_pos_b,fact_neg_b,fact_pos_e,fact_neg_e)
!
! Mnemonics: 
! fact_pos_b(i) - integration step in positive direction starts at point i
! fact_pos_e(i) - integration step in positive direction ends at point i
! fact_neg_b(i) - integration step in negative direction starts at point i
! fact_neg_e(i) - integration step in negative direction ends at point i
!
  USE plagrange_mod
!
  IMPLICIT NONE
!
!  logical, parameter :: stepmode=.true.
  LOGICAL, PARAMETER :: stepmode=.FALSE.
  INTEGER, PARAMETER :: npoi=6, nder=0, npoihalf=npoi/2, nstepmin=8
  DOUBLE PRECISION, PARAMETER :: bparabmax=0.2d0
!
  INTEGER :: i,ibeg,iend,npart,istep,ibmin,npassing,npassing_prev
  INTEGER :: ncomp,nreal
  INTEGER :: ncross_l,ncross_r,ib,ie,intb,inte,k,imid,isplit
!
  DOUBLE PRECISION :: subsqmin,ht,ht2,bparab,x1,x2,f1,f2
!
  INTEGER, DIMENSION(1)              :: idummy
  INTEGER, DIMENSION(1:iend)         :: phi_divide
  INTEGER, DIMENSION(:), ALLOCATABLE :: icross_l,icross_r
!
  DOUBLE PRECISION, DIMENSION(npoi)             :: tp,weight
  DOUBLE PRECISION, DIMENSION(0:nder,npoi)      :: coeff
  DOUBLE PRECISION, DIMENSION(0:npart)          :: eta
  DOUBLE PRECISION, DIMENSION(ibeg:iend)        :: phi_mfl,bhat_mfl
  DOUBLE PRECISION, DIMENSION(ibeg:iend,nreal)  :: arr_real
  DOUBLE COMPLEX,   DIMENSION(ibeg:iend,ncomp)  :: arr_comp
  DOUBLE PRECISION, DIMENSION(ibeg:iend)        :: delt_pos,delt_neg
  DOUBLE PRECISION, DIMENSION(ibeg:iend)        :: fact_pos_b,fact_neg_b
  DOUBLE PRECISION, DIMENSION(ibeg:iend)        :: fact_pos_e,fact_neg_e
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE   :: phi_new,bhat_new
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: arr_real_new
  DOUBLE COMPLEX,   DIMENSION(:,:), ALLOCATABLE :: arr_comp_new
!
  CALL fix_phiplacement_problem(ibeg,iend,npart,subsqmin,        &
                                phi_mfl,bhat_mfl,eta)
!
  phi_divide=1
!
  delt_pos(ibeg+1:iend)=phi_mfl(ibeg+1:iend)-phi_mfl(ibeg:iend-1)
  fact_pos_b=1.d0
  fact_pos_e=1.d0
!
! determine level crossings:
!
  idummy=MINLOC(bhat_mfl(ibeg:iend))
  ibmin=idummy(1)+ibeg-1
!
  ncross_l=0
  IF(ibmin.GT.ibeg) THEN
    istep=ibmin
    DO i=0,npart
      IF(1.d0-bhat_mfl(istep)*eta(i).GT.subsqmin) THEN
        npassing=i
      ELSE
        EXIT
      ENDIF
    ENDDO
    npassing_prev=npassing
    DO istep=ibmin-1,ibeg,-1
      DO i=0,npart
        IF(1.d0-bhat_mfl(istep)*eta(i).GT.subsqmin) THEN
          npassing=i
        ELSE
          EXIT
        ENDIF
      ENDDO
      IF(npassing.LT.npassing_prev) THEN
        ncross_l=ncross_l+1
        npassing_prev=npassing
      ENDIF
    ENDDO
    IF(ncross_l.GT.0) THEN
      ALLOCATE(icross_l(ncross_l))
      ncross_l=0
      istep=ibmin
      DO i=0,npart
        IF(1.d0-bhat_mfl(istep)*eta(i).GT.subsqmin) THEN
          npassing=i
        ELSE
          EXIT
        ENDIF
      ENDDO
      npassing_prev=npassing
      DO istep=ibmin-1,ibeg,-1
        DO i=0,npart
          IF(1.d0-bhat_mfl(istep)*eta(i).GT.subsqmin) THEN
            npassing=i
          ELSE
            EXIT
          ENDIF
        ENDDO
        IF(npassing.LT.npassing_prev) THEN
          ncross_l=ncross_l+1
          icross_l(ncross_l)=istep
          npassing_prev=npassing
        ENDIF
      ENDDO
    ENDIF
  ENDIF
!
  ncross_r=0
  IF(ibmin.LT.iend) THEN
    istep=ibmin
    DO i=0,npart
      IF(1.d0-bhat_mfl(istep)*eta(i).GT.subsqmin) THEN
        npassing=i
      ELSE
        EXIT
      ENDIF
    ENDDO
    npassing_prev=npassing
    DO istep=ibmin+1,iend
      DO i=0,npart
        IF(1.d0-bhat_mfl(istep)*eta(i).GT.subsqmin) THEN
          npassing=i
        ELSE
          EXIT
        ENDIF
      ENDDO
      IF(npassing.LT.npassing_prev) THEN
        ncross_r=ncross_r+1
        npassing_prev=npassing
      ENDIF
    ENDDO
    IF(ncross_r.GT.0) THEN
      ALLOCATE(icross_r(ncross_r))
      ncross_r=0
      istep=ibmin
      DO i=0,npart
        IF(1.d0-bhat_mfl(istep)*eta(i).GT.subsqmin) THEN
          npassing=i
        ELSE
          EXIT
        ENDIF
      ENDDO
      npassing_prev=npassing
      DO istep=ibmin+1,iend
        DO i=0,npart
          IF(1.d0-bhat_mfl(istep)*eta(i).GT.subsqmin) THEN
            npassing=i
          ELSE
            EXIT
          ENDIF
        ENDDO
        IF(npassing.LT.npassing_prev) THEN
          ncross_r=ncross_r+1
          icross_r(ncross_r)=istep
          npassing_prev=npassing
        ENDIF
      ENDDO
    ENDIF
  ENDIF
!
! place ibmin to an odd point:
!
  IF(MOD(ibmin-ibeg,2).EQ.1) THEN
    IF(ncross_l.GT.0.AND.ncross_r.GT.0) THEN
      IF(icross_r(1)-ibmin.GT.ibmin-icross_l(1)) THEN
        ibmin=ibmin+1
      ELSE
        ibmin=ibmin-1
      ENDIF
    ELSEIF(ncross_l.GT.0) THEN
      ibmin=ibmin+1
    ELSEIF(ncross_r.GT.0) THEN
      ibmin=ibmin-1
    ENDIF
  ENDIF
!
! check the number of steps in sub-intervals for parabolic bhat term:
!
  IF(ncross_l.GT.0) THEN
    ie=icross_l(1)
    DO i=2,ncross_l
      ib=icross_l(i)
      IF(ie-ib.LT.nstepmin) THEN
        imid=(ib+ie)/2
        x1=phi_mfl(imid)-phi_mfl(ib)
        x2=phi_mfl(ie)-phi_mfl(ib)
        f1=bhat_mfl(imid)-bhat_mfl(ib)
        f2=bhat_mfl(ie)-bhat_mfl(ib)
        bparab=ABS((f1*x2-f2*x1)*x2/((x1-x2)*x1*f2))
        IF(bparab.GT.bparabmax) THEN
          isplit=2*MAX(NINT(0.5*float(nstepmin)/float(ie-ib)),1)
          phi_divide(ib+1:ie)=isplit
        ENDIF
      ENDIF
      ie=ib
    ENDDO
    ib=ibeg
    IF(ie-ib.LT.nstepmin) THEN
      isplit=2*MAX(NINT(0.5*float(nstepmin)/float(ie-ib)),1)
      phi_divide(ib+1:ie)=isplit
    ENDIF
  ENDIF
!
  IF(ncross_r.GT.0) THEN
    ib=icross_r(1)
    DO i=2,ncross_r
      ie=icross_r(i)
      IF(ie-ib.LT.nstepmin) THEN
        imid=(ib+ie)/2
        x1=phi_mfl(imid)-phi_mfl(ib)
        x2=phi_mfl(ie)-phi_mfl(ib)
        f1=bhat_mfl(imid)-bhat_mfl(ib)
        f2=bhat_mfl(ie)-bhat_mfl(ib)
        bparab=ABS((f1*x2-f2*x1)*x2/((x1-x2)*x1*f2))
        IF(bparab.GT.bparabmax) THEN
          isplit=2*MAX(NINT(0.5*float(nstepmin)/float(ie-ib)),1)
          phi_divide(ib+1:ie)=isplit
        ENDIF
      ENDIF
      ib=ie
    ENDDO
    ie=iend
    IF(ie-ib.LT.nstepmin) THEN
      isplit=2*MAX(NINT(0.5*float(nstepmin)/float(ie-ib)),1)
      phi_divide(ib+1:ie)=isplit
    ENDIF
  ENDIF
!
  IF(MAXVAL(phi_divide).GT.1) RETURN
!
! change the integration variable phi -> sqrt(phi-phi0):
!
  IF(stepmode) THEN
!
    ALLOCATE(phi_new(ibeg:iend),bhat_new(ibeg:iend))
    ALLOCATE(arr_real_new(ibeg:iend,nreal))
    ALLOCATE(arr_comp_new(ibeg:iend,ncomp))
    phi_new=phi_mfl
    bhat_new=bhat_mfl
    arr_real_new=arr_real
    arr_comp_new=arr_comp
!
    ie=ibmin
    DO i=1,ncross_l
      ib=icross_l(i)
      ht=SQRT(phi_mfl(ie)-phi_mfl(ib))/float(ie-ib)
      ht2=ht**2
      k=ib
      DO istep=ib+1,ie-1
        phi_new(istep)=phi_mfl(ib)+ht2*float(istep-ib)**2
        fact_pos_e(istep)=2.d0*ht*float(istep-ib)
        delt_pos(istep)=ht
        DO WHILE(phi_mfl(k).LT.phi_new(istep))
          k=k+1
        ENDDO
        intb=MAX(ibeg,MIN(iend-npoi+1,k-npoihalf))
        inte=intb+npoi-1
!
        CALL plagrange_coeff(npoi,nder,phi_new(istep),phi_mfl(intb:inte),coeff)
!
        bhat_new(istep)=SUM(coeff(0,:)*bhat_mfl(intb:inte))
        arr_real_new(istep,:)=MATMUL(coeff(0,:),arr_real(intb:inte,:))
        arr_comp_new(istep,:)=MATMUL(coeff(0,:),arr_comp(intb:inte,:))
      ENDDO
      delt_pos(ie)=ht
      fact_pos_e(ie)=2.d0*ht*float(ie-ib)
      fact_pos_b(ib)=0.d0
      fact_pos_b(ib+1:ie-1)=fact_pos_e(ib+1:ie-1)
      ie=ib
    ENDDO
!
    ib=ibmin
    DO i=1,ncross_r
      ie=icross_r(i)
      ht=SQRT(phi_mfl(ie)-phi_mfl(ib))/float(ie-ib)
      ht2=ht**2
      k=ib
      DO istep=ib+1,ie-1
        phi_new(istep)=phi_mfl(ie)-ht2*float(ie-istep)**2
        delt_pos(istep)=ht
        fact_pos_b(istep)=2.d0*ht*float(ie-istep)
        DO WHILE(phi_mfl(k).LT.phi_new(istep))
          k=k+1
        ENDDO
        intb=MAX(ibeg,MIN(iend-npoi+1,k-npoihalf))
        inte=intb+npoi-1
!
        CALL plagrange_coeff(npoi,nder,phi_new(istep),phi_mfl(intb:inte),coeff)
!
        bhat_new(istep)=SUM(coeff(0,:)*bhat_mfl(intb:inte))
        arr_real_new(istep,:)=MATMUL(coeff(0,:),arr_real(intb:inte,:))
        arr_comp_new(istep,:)=MATMUL(coeff(0,:),arr_comp(intb:inte,:))
      ENDDO
      delt_pos(ie)=ht
      fact_pos_b(ib)=2.d0*ht*float(ie-ib)
      fact_pos_e(ie)=0.d0
      fact_pos_e(ib+1:ie-1)=fact_pos_b(ib+1:ie-1)
      ib=ie
    ENDDO
!
    phi_mfl=phi_new
    bhat_mfl=bhat_new
    arr_real=arr_real_new
    arr_comp=arr_comp_new
!
    DEALLOCATE(phi_new,bhat_new,arr_real_new,arr_comp_new)
!
  ENDIF
!
  delt_neg(ibeg:iend-1)=delt_pos(ibeg+1:iend)
  fact_neg_b=fact_pos_e
  fact_neg_e=fact_pos_b
!
  IF(ALLOCATED(icross_l)) DEALLOCATE(icross_l)
  IF(ALLOCATED(icross_r)) DEALLOCATE(icross_r)
!
END SUBROUTINE rearrange_phideps
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
SUBROUTINE fix_phiplacement_problem(ibeg,iend,npart,subsqmin,        &
                                    phi_mfl,bhat_mfl,eta)
!
  USE device_mod
!
  IMPLICIT NONE
!
  INTEGER :: i,ibeg,iend,npart,istep,ibmin,npassing,npassing_prev
  INTEGER :: ncross_l,ncross_r,ib,ie
!
  DOUBLE PRECISION :: subsqmin
!
  INTEGER, DIMENSION(1)              :: idummy
  INTEGER, DIMENSION(:), ALLOCATABLE :: icross_l,icross_r
!
  DOUBLE PRECISION, DIMENSION(0:npart)        :: eta
  DOUBLE PRECISION, DIMENSION(ibeg:iend)      :: phi_mfl,bhat_mfl
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: eta_cross_l,eta_cross_r
!
! determine level crossings:
!
  idummy=MINLOC(bhat_mfl(ibeg:iend))
  ibmin=idummy(1)+ibeg-1
!
  ncross_l=0
  IF(ibmin.GT.ibeg) THEN
    istep=ibmin
    DO i=0,npart
      IF(1.d0-bhat_mfl(istep)*eta(i).GT.subsqmin) THEN
        npassing=i
      ELSE
        EXIT
      ENDIF
    ENDDO
    npassing_prev=npassing
    DO istep=ibmin-1,ibeg,-1
      DO i=0,npart
        IF(1.d0-bhat_mfl(istep)*eta(i).GT.subsqmin) THEN
          npassing=i
        ELSE
          EXIT
        ENDIF
      ENDDO
      IF(npassing.LT.npassing_prev) THEN
        ncross_l=ncross_l+1
        npassing_prev=npassing
      ENDIF
    ENDDO
    IF(ncross_l.GT.0) THEN
      ALLOCATE(icross_l(ncross_l),eta_cross_l(ncross_l))
      ncross_l=0
      istep=ibmin
      DO i=0,npart
        IF(1.d0-bhat_mfl(istep)*eta(i).GT.subsqmin) THEN
          npassing=i
        ELSE
          EXIT
        ENDIF
      ENDDO
      npassing_prev=npassing
      DO istep=ibmin-1,ibeg,-1
        DO i=0,npart
          IF(1.d0-bhat_mfl(istep)*eta(i).GT.subsqmin) THEN
            npassing=i
          ELSE
            EXIT
          ENDIF
        ENDDO
        IF(npassing.LT.npassing_prev) THEN
          ncross_l=ncross_l+1
          icross_l(ncross_l)=istep
          eta_cross_l(ncross_l)=eta(npassing_prev)
          npassing_prev=npassing
        ENDIF
      ENDDO
      DO i=1,ncross_l
        istep=icross_l(i)
        IF(ABS(bhat_mfl(istep-1)*eta_cross_l(i)-1.d0).LT. &
           ABS(bhat_mfl(istep)  *eta_cross_l(i)-1.d0)) THEN
          OPEN(111,file='phi_placement_problem.dat',position='append')
          WRITE(111,*) ' propagator tag = ',fieldpropagator%tag, &
                       ' step number = ',istep-1,                &
                       ' 1 / bhat = ',1.d0/bhat_mfl(istep-1),    &
                       ' eta = ',eta_cross_l(i)
          CLOSE(111)
          bhat_mfl(istep-1)=1/eta_cross_l(i)
        ELSEIF(ABS(bhat_mfl(istep+1)*eta_cross_l(i)-1.d0).LT. &
               ABS(bhat_mfl(istep)  *eta_cross_l(i)-1.d0)) THEN
          OPEN(111,file='phi_placement_problem.dat',position='append')
          WRITE(111,*) ' propagator tag = ',fieldpropagator%tag, &
                       ' step number = ',istep+1,                &
                       ' 1 / bhat = ',1.d0/bhat_mfl(istep+1),    &
                       ' eta = ',eta_cross_l(i)
          bhat_mfl(istep+1)=1/eta_cross_l(i)
          CLOSE(111)
        ENDIF
      ENDDO
      DEALLOCATE(icross_l,eta_cross_l)
    ENDIF
  ENDIF
!
  ncross_r=0
  IF(ibmin.LT.iend) THEN
    istep=ibmin
    DO i=0,npart
      IF(1.d0-bhat_mfl(istep)*eta(i).GT.subsqmin) THEN
        npassing=i
      ELSE
        EXIT
      ENDIF
    ENDDO
    npassing_prev=npassing
    DO istep=ibmin+1,iend
      DO i=0,npart
        IF(1.d0-bhat_mfl(istep)*eta(i).GT.subsqmin) THEN
          npassing=i
        ELSE
          EXIT
        ENDIF
      ENDDO
      IF(npassing.LT.npassing_prev) THEN
        ncross_r=ncross_r+1
        npassing_prev=npassing
      ENDIF
    ENDDO
    IF(ncross_r.GT.0) THEN
      ALLOCATE(icross_r(ncross_r),eta_cross_r(ncross_r))
      ncross_r=0
      istep=ibmin
      DO i=0,npart
        IF(1.d0-bhat_mfl(istep)*eta(i).GT.subsqmin) THEN
          npassing=i
        ELSE
          EXIT
        ENDIF
      ENDDO
      npassing_prev=npassing
      DO istep=ibmin+1,iend
        DO i=0,npart
          IF(1.d0-bhat_mfl(istep)*eta(i).GT.subsqmin) THEN
            npassing=i
          ELSE
            EXIT
          ENDIF
        ENDDO
        IF(npassing.LT.npassing_prev) THEN
          ncross_r=ncross_r+1
          icross_r(ncross_r)=istep
          eta_cross_r(ncross_r)=eta(npassing_prev)
          npassing_prev=npassing
        ENDIF
      ENDDO
      DO i=1,ncross_r
        istep=icross_r(i)
        IF(ABS(bhat_mfl(istep-1)*eta_cross_r(i)-1.d0).LT. &
           ABS(bhat_mfl(istep)  *eta_cross_r(i)-1.d0)) THEN
          OPEN(111,file='phi_placement_problem.dat',position='append')
          WRITE(111,*) ' propagator tag = ',fieldpropagator%tag, &
                       ' step number = ',istep-1,                &
                       ' 1 / bhat = ',1.d0/bhat_mfl(istep-1),    &
                       ' eta = ',eta_cross_r(i)
          CLOSE(111)
          bhat_mfl(istep-1)=1/eta_cross_r(i)
        ELSEIF(ABS(bhat_mfl(istep+1)*eta_cross_r(i)-1.d0).LT. &
               ABS(bhat_mfl(istep)  *eta_cross_r(i)-1.d0)) THEN
          OPEN(111,file='phi_placement_problem.dat',position='append')
          WRITE(111,*) ' propagator tag = ',fieldpropagator%tag, &
                       ' step number = ',istep+1,                &
                       ' 1 / bhat = ',1.d0/bhat_mfl(istep+1),    &
                       ' eta = ',eta_cross_r(i)
          CLOSE(111)
          bhat_mfl(istep+1)=1/eta_cross_r(i)
        ENDIF
      ENDDO
      DEALLOCATE(icross_r,eta_cross_r)
    ENDIF
  ENDIF
!
END SUBROUTINE fix_phiplacement_problem
