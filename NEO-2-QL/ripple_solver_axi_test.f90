!
  MODULE ntv_eqmat_mod
    INTEGER                                         :: nz_symm,nz_asymm,nz_regper
    INTEGER                                         :: nz_per_pos,nz_per_neg
    INTEGER,          DIMENSION(:),     ALLOCATABLE :: irow_symm,icol_symm
    INTEGER,          DIMENSION(:),     ALLOCATABLE :: irow_regper,icol_regper
    INTEGER,          DIMENSION(:),     ALLOCATABLE :: irow_asymm,icol_asymm
    INTEGER,          DIMENSION(:),     ALLOCATABLE :: irow_per_pos,icol_per_pos
    INTEGER,          DIMENSION(:),     ALLOCATABLE :: irow_per_neg,icol_per_neg
    DOUBLE PRECISION, DIMENSION(:),     ALLOCATABLE :: amat_symm
    DOUBLE PRECISION, DIMENSION(:),     ALLOCATABLE :: amat_regper
    complex(kind=kind(1d0)), DIMENSION(:), ALLOCATABLE :: amat_asymm
    DOUBLE PRECISION, DIMENSION(:,:),   ALLOCATABLE :: f0_coll,f0_ttmp
    DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: f0_coll_all,f0_ttmp_all
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
       nvel,vel_array,num_spec

! collpar - the same as $\kappa$ - inverse mean-free path times 4
  USE lapack_band
  USE rkstep_mod
  USE polleg_mod
  ! USE binarysplit_mod, ONLY : bsfunc_reconstruct_levels !<- LOCAL
  ! new switch in propagator_mod
  USE propagator_mod,ONLY : prop_ripple_plot,prop_reconstruct,flux_mr,flux_pl, &
                            eta_modboundary_l,eta_modboundary_r,               &
                            sw_first_prop,sw_last_prop,                        &
                            prop_reconstruct_levels!,                           &
                            !prop_fileformat
  USE sparse_mod, ONLY : sparse_talk,sparse_solve_method,sparse_solve, &
       column_full2pointer,remap_rc,sparse_solver_test
  USE mag_interface_mod, ONLY: average_bhat,average_one_over_bhat,             &
                               surface_boozer_B00,travis_convfac,              &
                               mag_magfield,boozer_s

  USE development

  !*****************************
  ! HDF5
  !*****************************
  USE hdf5_tools
  !! Modification by Andreas F. Martitsch (28.07.2015)
  ! MPI SUPPORT for multi-species part
  ! (run with, e.g.,  mpiexec -np 3 ./neo2.x)
  USE ntv_eqmat_mod, ONLY : nz_symm,irow_symm,icol_symm,amat_symm,             &
                            nz_regper,irow_regper,icol_regper,amat_regper
  USE mpiprovider_module
  USE collop
  !! End Modification by Andreas F. Martitsch (28.07.2015)

  IMPLICIT NONE
  !INTEGER, PARAMETER :: dp = KIND(1.0d0)

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
  REAL(kind=dp), DIMENSION(:,:),  ALLOCATABLE, INTENT(out)   :: qflux
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
  REAL(kind=dp) :: time_start,time_factorization,time_solver, time1, time2, time3, time4, time5, time6, time7


  !------------------------------------------------------------------------
  ! SERGEI
  !------------------------------------------------------------------------
  ! additional definitions from Sergei
  ! you can also use "double precision" isntead of  "REAL(kind=dp)"
  INTEGER :: itotstep,ntotstep,npart_loc,info
  INTEGER :: ndim,nhalf,istep,nstep,nback,npoints,npassing,ioddeven
  INTEGER :: ipart,nhalfstep,idir,inhom
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
  DOUBLE PRECISION :: diflam_flux,coef_cf
  DOUBLE PRECISION :: hxeta,fun_bound_new,amin2ovb
  DOUBLE PRECISION :: a11,a12,a21,a22,determ,deleta_b,dum_phi,dum_a3
  DOUBLE PRECISION :: exl,exprenorm,accurfac,renormfac
!
  DOUBLE PRECISION, DIMENSION(6)           :: alp,bet,gam,del
  DOUBLE PRECISION, DIMENSION(4)           :: old_fun
  DOUBLE PRECISION, DIMENSION(ndim0)       :: ybeg
  DOUBLE PRECISION, DIMENSION(4,4)         :: a_beg,b_beg,a_end,b_end
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: amat,bvec_lapack,deriv_coef
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: prod_p,prod_m,prod
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
  DOUBLE PRECISION, DIMENSION(:,:),   ALLOCATABLE :: alambd,Vg_vp_over_B
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: redmat_forw_l,redmat_back_l
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: redmat_forw_r,redmat_back_r
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: delta_eta
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: enu_coef        !!!term[1]
  INTEGER :: km1,kp1,m,m1,kb,kb1,ke,ke1,mfactorial,nplp1,k_bound,k1
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE ::  source_pow !!!term[2]
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE ::  conv_pow !!!term[3]
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE ::  alampow !!!term[3]
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE ::  vrecurr !!!term[3]
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE ::  dellampow !!!term[3]
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
  INTEGER :: iter,nphiequi,npassing_next
  !INTEGER :: niter ! now defined via rkstep_mod (neo_mod.f90)
  DOUBLE PRECISION :: delphim1,deloneovb,step_factor_p,step_factor_m
  !DOUBLE PRECISION :: epserr_iter
  INTEGER,          DIMENSION(:),   ALLOCATABLE :: ind_start,irow,icol,ipcol
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: amat_sp,funsol_p,funsol_m
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: bvec_sp,bvec_iter,bvec_lor
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: bvec_prev
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: flux_vector,source_vector
  INTEGER :: isw_lor,isw_ene,isw_intp
  INTEGER,          DIMENSION(:),       ALLOCATABLE :: npl
  DOUBLE PRECISION, DIMENSION(:,:,:),   ALLOCATABLE :: rhs_mat_fzero
  DOUBLE PRECISION, DIMENSION(:,:,:),   ALLOCATABLE :: rhs_mat_lorentz,q_rip
  DOUBLE PRECISION, DIMENSION(:,:,:),   ALLOCATABLE :: rhs_mat_energ
  DOUBLE PRECISION, DIMENSION(:,:,:),   ALLOCATABLE :: pleg_bra,pleg_ket
  DOUBLE PRECISION, DIMENSION(:,:),     ALLOCATABLE :: convol_flux,convol_curr
  DOUBLE PRECISION, DIMENSION(:,:),     ALLOCATABLE :: scalprod_pleg
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: phi_mfl
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: bhat_mfl,geodcu_mfl,h_phi_mfl
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: dlogbdphi_mfl
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: delt_pos,delt_neg
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: fact_pos_b,fact_neg_b
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: fact_pos_e,fact_neg_e

  !! End Modifications by Andreas F. Martitsch (13.06.2014)
  !! Modification by Andreas F. Martitsch (28.07.2015)
  !  multi-species part
  INTEGER :: ispec, ispecp, ispecpp ! species indices
  INTEGER :: drive_spec
  INTEGER :: isw_regper, ipart1
  DOUBLE PRECISION :: deleta_factor
  DOUBLE PRECISION,   DIMENSION(:,:,:), ALLOCATABLE :: source_vector_all
  REAL(kind=dp), DIMENSION(:,:,:,:), ALLOCATABLE :: qflux_allspec, qflux_allspec_tmp
  LOGICAL :: problem_type
  DOUBLE PRECISION, DIMENSION(0:num_spec-1) :: break_cond1
  DOUBLE PRECISION, DIMENSION(0:num_spec-1) :: break_cond2
  INTEGER :: prop_fileformat=0
  !! End Modification by Andreas F. Martitsch (28.07.2015)

  !***************************
  ! HDF5
  !***************************
  integer(HID_T) :: h5id_final_spitzer, h5id_phi_mesh, h5id_dentf, h5id_enetf, h5id_spitf, h5id_sizeplot
  integer(HID_T) :: h5id_bhat_mfl

  !**********************************************************
  ! For faster read/write of whole HDF5 dataset
  !**********************************************************
  double precision, dimension(:), allocatable       :: phi_mfl_h5, bhat_mfl_h5, npassing_h5
  double precision, dimension(:,:,:,:), allocatable :: dentf_p_h5, enetf_p_h5, spitf_p_h5
  double precision, dimension(:,:,:,:), allocatable :: dentf_m_h5, enetf_m_h5, spitf_m_h5

  ! integer :: isw_axisymm=0 ! now in collisionality_mod
  niter=1000
  isw_regper=1       !regulariization by periodic boundary condition

  !! Modification by Andreas F. Martitsch (28.07.2015)
  ! multi-species part - MPI rank determines species
  ispec = mpro%getRank()
  PRINT *,"Species: ", ispec
  CALL collop_set_species(ispec)
  !PRINT *,'asource: ',asource(:,1)
  !PRINT *,'anumm:',anumm(1,:)
  !PRINT *,'denmm:',denmm(1,:)
  !STOP
  !! End Modification by Andreas F. Martitsch (28.07.2015)
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

  if (isw_momentum .eq. 1) then ! Grid
     print *, 'isw_momentum = ',isw_momentum,' not implemented in ripple solver!'
     print *, 'I stop here'
     stop
  end if



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
goto 10
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
10 continue
!
  ! place for boundary
  npass_l = npass_l + 1 - ignore_lb
  npass_r = npass_r + 1 - ignore_rb

  ! allocate and copy the magnetic stuff
  ub_mag = UBOUND(fieldpropagator%coords%x2,1)
  ibeg   = 0 - 2*modify_bl
  iend   = ub_mag + 2*modify_br

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

  IF (ALLOCATED(h_phi_mfl)) DEALLOCATE(h_phi_mfl)
  ALLOCATE(h_phi_mfl(ibeg:iend))
  h_phi_mfl(0:ub_mag) = fieldpropagator%mdata%h_phi
  sign_of_bphi= sign(1.d0,h_phi_mfl(0))                                !08.12.08
  h_phi_mfl(0:ub_mag)=h_phi_mfl(0:ub_mag)*sign_of_bphi                 !08.12.08
  geodcu_mfl(0:ub_mag)=geodcu_mfl(0:ub_mag)*sign_of_bphi               !08.12.08

  IF (ALLOCATED(dlogbdphi_mfl)) DEALLOCATE(dlogbdphi_mfl)
  ALLOCATE(dlogbdphi_mfl(ibeg:iend))
  dlogbdphi_mfl(0:ub_mag) = fieldpropagator%mdata%dlogbdphi

  IF(modify_bl.EQ.1) THEN
    phi_mfl(-2:-1)=phi_mfl(0)
    bhat_mfl(-2:-1)=bhat_changed_l
    geodcu_mfl(-2:-1)=geodcu_mfl(0)
    h_phi_mfl(-2:-1)=h_phi_mfl(0)
    dlogbdphi_mfl(-2:-1)=dlogbdphi_mfl(0)
  ENDIF
  IF(modify_br.EQ.1) THEN
    phi_mfl(iend-1:iend)=phi_mfl(ub_mag)
    bhat_mfl(iend-1:iend)=bhat_changed_r
    geodcu_mfl(iend-1:iend)=geodcu_mfl(ub_mag)
    h_phi_mfl(iend-1:iend)=h_phi_mfl(ub_mag)
    dlogbdphi_mfl(iend-1:iend)=dlogbdphi_mfl(ub_mag)
  ENDIF
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
     solver_talk=0
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

     if (prop_fileformat .eq. 0) then
        OPEN(123,file='bhat_mfl.dat')
        DO i=ibeg,iend
           WRITE(123,*) phi_mfl(i),bhat_mfl(i),geodcu_mfl(i),h_phi_mfl(i),dlogbdphi_mfl(i)
        ENDDO
        CLOSE(123)

     elseif (prop_fileformat .eq. 1) then
        !**********************************************************
        ! HDF5 of bhat_mfl
        !**********************************************************
        write(propname,*) fieldpropagator%tag
        call h5_create('bhat_mfl_' // trim(adjustl(propname)) // '.h5', h5id_bhat_mfl)

        call h5_add(h5id_bhat_mfl, 'phi_mfl', phi_mfl(ibeg:iend), lbound(phi_mfl), ubound(phi_mfl))
        call h5_add(h5id_bhat_mfl, 'bhat_mfl', bhat_mfl(ibeg:iend), lbound(bhat_mfl), ubound(bhat_mfl))
        call h5_add(h5id_bhat_mfl, 'geodcu_mfl', geodcu_mfl(ibeg:iend), lbound(geodcu_mfl), ubound(geodcu_mfl))
        call h5_add(h5id_bhat_mfl, 'h_phi_mfl', h_phi_mfl(ibeg:iend), lbound(h_phi_mfl), ubound(h_phi_mfl))
        call h5_add(h5id_bhat_mfl, 'dlogbdphi_mfl', dlogbdphi_mfl(ibeg:iend), lbound(dlogbdphi_mfl), ubound(dlogbdphi_mfl))

        call h5_close(h5id_bhat_mfl)
  end if

     OPEN(123,file='eta.dat')
     DO i=0,ub_eta
        WRITE(123,*) phi_mfl(ibeg),eta(i)
        WRITE(123,*) phi_mfl(iend),eta(i)
        WRITE(123,*) ' '
     ENDDO
     CLOSE(123)

     ! PAUSE 'bmod and eta written' !Warning in gfortran-4.7
  END IF

  !------------------------------------------------------------------------
  ! SERGEI
  !------------------------------------------------------------------------
!
! Check for axisymmetry:
!
  if(isw_axisymm.eq.1.and.npass_l.ne.npass_r) then
     print *, 'npass_l,npass_r ',npass_l,npass_r
    print *,'ripple_solver: cannot run axisymmetric mode, sizes do not fit'
    ierr=1
    return
  endif
!
  iplot=prop_ripple_plot
!
! Preparation of coefficients for the kinetic equation solver
!
  ALLOCATE(deriv_coef(4,0:npart+1))
  ALLOCATE(enu_coef(4,npart+1))                                    !!!term[1]
  ALLOCATE(alambd(0:npart+3,ibeg:iend),Vg_vp_over_B(0:npart,ibeg:iend))
  ALLOCATE(scalprod_pleg(0:lag,0:legmax))                          !!!term[3]
  ALLOCATE(alampow(legmax+1,0:npart+1))                            !!!terms[2,3]
  ALLOCATE(vrecurr(0:legmax,0:3,1:npart+1))                        !!!terms[2,3]
  ALLOCATE(dellampow(4,1:npart+1))                                 !!!terms[1-3]
  ALLOCATE(convol_polpow(0:legmax,1:npart+3))                      !!!terms[2,3]
  ALLOCATE(pleg_bra(0:legmax,1:npart+1,ibeg:iend))                 !!!terms[2,3]
  ALLOCATE(pleg_ket(0:legmax,1:npart+1,ibeg:iend))                 !!!terms[2,3]
  ALLOCATE(npl(ibeg:iend))
  ALLOCATE(rhs_mat_fzero(4,ibeg:iend,0:1))
  ALLOCATE(rhs_mat_lorentz(5,npart+1,ibeg:iend))
  ALLOCATE(rhs_mat_energ(4,npart+1,ibeg:iend))
  ALLOCATE(q_rip(npart+2,ibeg:iend,0:2))
  ALLOCATE(convol_flux(npart+1,ibeg:iend),convol_curr(npart+1,ibeg:iend))
  ALLOCATE(ind_start(ibeg:iend))
  if (iplot.EQ.1) then
    allocate(derivs_plot(0:3,4,npart+1,ibeg:iend))
    derivs_plot = 1.234e5
  end if
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
  allocate(delt_pos(ibeg:iend),delt_neg(ibeg:iend))
  allocate(fact_pos_b(ibeg:iend),fact_neg_b(ibeg:iend))
  allocate(fact_pos_e(ibeg:iend),fact_neg_e(ibeg:iend))
!
  call rearrange_phideps_old(ibeg,iend,npart,subsqmin,phi_divide,    &
                         phi_mfl,bhat_mfl,geodcu_mfl,h_phi_mfl,eta,  &
                         delt_pos,delt_neg,                          &
                         fact_pos_b,fact_neg_b,fact_pos_e,fact_neg_e)
!
  if(maxval(phi_divide).gt.1) then
    ierr=3
    DEALLOCATE(deriv_coef,npl)
    DEALLOCATE(rhs_mat_lorentz,rhs_mat_energ)
    DEALLOCATE(q_rip)
    DEALLOCATE(convol_flux,convol_curr)
    DEALLOCATE(pleg_bra,pleg_ket,scalprod_pleg)
    return
  endif
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
                             *(4.d0*eta0-eta(i))*geodcu_mfl(istep)/3.d0
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
  do istep=ibeg,iend-1
    ind_start(istep+1)=ind_start(istep)+2*(lag+1)*(npl(istep)+1)
  enddo
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
! term[2] (Legendre polynomials) -  ket-vector
! numbering of levels (N=npassing): 1-N - $f_n$, N+1 - $f^b$, N+2 - $f^a$
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
!
      enu_coef(:,i)=MATMUL(dellampow(:,i),bvec_lapack)*coefenu      !!!term[1]
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
    coef_cf=geodcu_mfl(istep)/bhat_mfl(istep)**2/h_phi_mfl(istep)
    convol_flux(1:npassing+1,istep)                                          &
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
    rhs_mat_fzero(:,istep,1)=deriv_coef(:,npassing+1)
!
! Change of the boundary layer width
!
!
!
! begin term[1]:
!
    rhs_mat_energ(:,1:npassing+1,istep)=enu_coef(:,1:npassing+1)
!
! end term[1]
!
    q_rip(1:npassing,istep,1)                                    &
          =Vg_vp_over_B(1:npassing,istep)-Vg_vp_over_B(0:npassing-1,istep)
    q_rip(npassing+1,istep,1)=-Vg_vp_over_B(npassing,istep)
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
!! Modification by Andreas F. Martitsch (16.09.2015)
! NEO-2 can treat now multiple species
! (move collpar from pleg_bra to pleg_ket to avoid mixing up
! of species-dependent parameters)
  pleg_bra=pleg_bra/collpar
  pleg_ket=pleg_ket*collpar
!! End Modification by Andreas F. Martitsch (16.09.2015)
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
  allocate(flux_vector(3,n_2d_size),source_vector(n_2d_size,3))
  flux_vector=0.d0
  source_vector=0.d0
!
  do istep=ibeg,iend
!
    ioddeven=mod(istep-ibeg,2) !0 for full RK step, 1 for half RK step
!
    if(istep.eq.ibeg) then
      step_factor_p=delt_pos(ibeg+1)*fact_pos_b(ibeg)/3.d0
      step_factor_m=delt_neg(ibeg)*fact_neg_e(ibeg)/3.d0
    elseif(istep.eq.iend) then
      step_factor_p=delt_pos(iend)*fact_pos_e(iend)/3.d0
      step_factor_m=delt_neg(iend-1)*fact_neg_b(iend)/3.d0
    elseif(ioddeven.eq.1) then
      step_factor_p=(delt_pos(istep+1)*fact_pos_b(istep)        &
                   + delt_pos(istep)*fact_pos_e(istep))/1.5d0
      step_factor_m=(delt_neg(istep-1)*fact_neg_b(istep)       &
                   + delt_neg(istep)*fact_neg_e(istep))/1.5d0
    else
      step_factor_p=(delt_pos(istep+1)*fact_pos_b(istep)        &
                   + delt_pos(istep)*fact_pos_e(istep))/3.d0
      step_factor_m=(delt_neg(istep-1)*fact_neg_b(istep)       &
                   + delt_neg(istep)*fact_neg_e(istep))/3.d0
    endif
!
    npassing=npl(istep)
!
    do m=0,lag
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
      if(istep.gt.ibeg) then
        npassing_prev=npl(istep-1)
        k_prev=ind_start(istep-1)+2*(npassing_prev+1)*m
        if(ioddeven.eq.1) then
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
        else
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
          if(npassing_prev.le.npassing) then
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
          else
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
          endif
          if(npassing_next.le.npassing) then
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
          else
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
          endif
        endif
      elseif(iplot.eq.1.and.isw_axisymm.ne.1) then
        source_vector(k+1:k+npassing+1,:)                                &
                       =flux_pl(npass_l*m+1:npass_l*(m+1),:)
      endif
!
      if(istep.lt.iend) then
        npassing_prev=npl(istep+1)
        k_prev=ind_start(istep+1)+2*(npassing_prev+1)*m
        if(ioddeven.eq.1) then
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
        else
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
          if(npassing_prev.le.npassing) then
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
          else
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
          endif
          if(npassing_next.le.npassing) then
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
          else
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
          endif
        endif
      elseif(iplot.eq.1.and.isw_axisymm.ne.1) then
        source_vector(k+npassing+2:k+2*npassing+2,:)                     &
                       =flux_mr(npass_r*(m+1):npass_r*m+1:-1,:)
      endif
!
    enddo
  enddo
!
!
! Determine the size of arrays (number of non-zero elements):
!
  nz=0
  nz_regper=0
!
! Co-passing: sigma=1
!
  istep=ibeg
  npassing=npl(istep)
!
! entry:
!
  do m=0,lag
    k=ind_start(istep)+2*(npassing+1)*m
!
    do ipart=1,npassing+1
      nz=nz+1
!      irow(nz)=k+ipart
!      icol(nz)=k+ipart
!      amat_sp(nz)=1.d0
    enddo
!
    if(isw_axisymm.eq.1) then
      k_prev=ind_start(iend)+2*(npassing+1)*m
!
      do ipart=1,npassing+1
        nz=nz+1
!        irow(nz)=k+ipart
!        icol(nz)=k_prev+ipart
!        amat_sp(nz)=-1.d0
      enddo
!
      if(isw_regper.eq.1.and.m.lt.1) then
        do ipart=1,npassing+1
!          if(ipart.le.npassing) then
!            deleta_factor=(eta(ipart)-eta(ipart-1))*bhat_mfl(istep)
!          else
!            deleta_factor=1.d0-eta(ipart-1)*bhat_mfl(istep)
!          endif
!
          do ipart1=1,npassing+1
            nz_regper=nz_regper+1
!            irow_regper(nz_regper)=k+ipart
!            icol_regper(nz_regper)=k_prev+ipart1
!            amat_regper(nz_regper)=deleta_factor
            nz_regper=nz_regper+1
!            irow_regper(nz_regper)=k+ipart
!            icol_regper(nz_regper)=k+npassing+1+ipart1
!            amat_regper(nz_regper)=deleta_factor
          enddo
!
        enddo
      endif
!
    endif
!
  enddo
!
  do istep=ibeg+1,iend
    npassing_prev=npl(istep-1)
    npassing=npl(istep)
!    delphim1=1.d0/delt_pos(istep)
!    deloneovb=0.5d0*delphim1*(1.d0/bhat_mfl(istep-1)-1.d0/bhat_mfl(istep))
!
    do m=0,lag
      k_prev=ind_start(istep-1)+2*(npassing_prev+1)*m
      k=ind_start(istep)+2*(npassing+1)*m
!
! free flight:
!
      do ipart=1,npassing+1
        nz=nz+1
!        irow(nz)=k+ipart
!        icol(nz)=k+ipart
!        amat_sp(nz)=delphim1
      enddo
!
      do ipart=1,npassing
        nz=nz+1
!        irow(nz)=k+ipart
!        icol(nz)=k_prev+ipart
!        amat_sp(nz)=-delphim1
      enddo
!
      if(npassing_prev.ge.npassing) then
        nz=nz+1
!        irow(nz)=k+npassing+1
!        icol(nz)=k_prev+npassing+1
!        amat_sp(nz)=-delphim1
      endif
!
! mirroring:
!
      if(npassing_prev.eq.npassing) then
!
        do kk=1,4
          nz=nz+1
!          irow(nz)=k+npassing+1
!          icol(nz)=k+npassing+kk-1
!          amat_sp(nz)=deloneovb*rhs_mat_fzero(kk,istep,0)
        enddo
!
        do kk=1,4
          nz=nz+1
!          irow(nz)=k+npassing+1
!          icol(nz)=k_prev+npassing_prev+kk-1
!          amat_sp(nz)=deloneovb*rhs_mat_fzero(kk,istep-1,0)
        enddo
!
      elseif(npassing_prev.gt.npassing) then
        nz=nz+1
!        irow(nz)=k_prev+npassing_prev+2
!        icol(nz)=k_prev+npassing_prev+1
!        amat_sp(nz)=-delphim1
      endif
!
! collisions:
!
      if(fact_pos_e(istep).ne.0.d0) then
!
! Lorentz operator:
!
        if(isw_lor.eq.1) then
!
          do ipart=1,npassing+1
            do kk=1,5
              do mm=0,lag
                nz=nz+1
!                irow(nz)=k+ipart
!                icol(nz)=k+max(0,ipart-3)+kk+2*(npassing+1)*(mm-m)
!                amat_sp(nz)=anumm(m,mm)*rhs_mat_lorentz(kk,ipart,istep) &
!                           *fact_pos_e(istep)
              enddo
            enddo
          enddo
!
! matching collisional flux through the boundary in backward Euler scheme:
!
          do kk=1,4
            do mm=0,lag
              nz=nz+1
!              irow(nz)=k+npassing+1
!              icol(nz)=k+npassing-1+kk+2*(npassing+1)*(mm-m)
!              amat_sp(nz)=-0.5d0*anumm(m,mm)*rhs_mat_fzero(kk,istep,1)  &
!                         *fact_pos_e(istep)
            enddo
          enddo
!
          do kk=1,4
            do mm=0,lag
              nz=nz+1
!              irow(nz)=k+npassing+1
!              icol(nz)=k_prev+npassing_prev-1+kk+2*(npassing_prev+1)*(mm-m)
!              amat_sp(nz)=0.5d0*anumm(m,mm)*rhs_mat_fzero(kk,istep-1,1)
!                         *fact_pos_b(istep-1)
            enddo
          enddo
!
        endif
!
!        nz_beg=nz+1
!
! energy diffusion operator:
!
        if(isw_ene.eq.1) then
!
          do ipart=1,npassing
            do kk=1,4
              do mm=0,lag
                nz=nz+1
!                irow(nz)=k+ipart
!                icol(nz)=k+max(0,ipart-2)+kk+2*(npassing+1)*(mm-m)
!                amat_sp(nz)=denmm(m,mm)*rhs_mat_energ(kk,ipart,istep)
              enddo
            enddo
          enddo
!
          do kk=1,4
            do mm=0,lag
              nz=nz+1
!              irow(nz)=k+npassing+1
!              icol(nz)=k+npassing-1+kk+2*(npassing+1)*(mm-m)
!              amat_sp(nz)=0.5d0*denmm(m,mm)                                  &
!                         *rhs_mat_energ(kk,npassing+1,istep)
            enddo
          enddo
!
          do kk=1,4
            do mm=0,lag
              nz=nz+1
!              irow(nz)=k+npassing+1
!              icol(nz)=k_prev+npassing_prev-1+kk+2*(npassing_prev+1)*(mm-m)
!              amat_sp(nz)=0.5d0*denmm(m,mm)                                  &
!                         *rhs_mat_energ(kk,npassing_prev+1,istep-1)
            enddo
          enddo
!
        endif
!
!        amat_sp(nz_beg:nz)=fact_pos_e(istep)*amat_sp(nz_beg:nz)
      endif
!
    enddo
!
  enddo
!
! Counter-passing: sigma=-1
!
  istep=iend
  npassing=npl(istep)
!
! entry:
!
  do m=0,lag
    k=ind_start(istep)+2*(npassing+1)*m+2*npassing+3
!
    do ipart=1,npassing+1
      nz=nz+1
!      irow(nz)=k-ipart
!      icol(nz)=k-ipart
!      amat_sp(nz)=1.d0
    enddo
!
    if(isw_axisymm.eq.1) then
      k_prev=ind_start(ibeg)+2*(npassing+1)*m+2*npassing+3
!
      do ipart=1,npassing+1
        nz=nz+1
!        irow(nz)=k-ipart
!        icol(nz)=k_prev-ipart
!        amat_sp(nz)=-1.d0
      enddo
!
      if(isw_regper.eq.1.and.m.lt.1) then
        do ipart=1,npassing+1
!          if(ipart.le.npassing) then
!            deleta_factor=(eta(ipart)-eta(ipart-1))*bhat_mfl(istep)
!          else
!            deleta_factor=1.d0-eta(ipart-1)*bhat_mfl(istep)
!          endif
!
          do ipart1=1,npassing+1
            nz_regper=nz_regper+1
!            irow_regper(nz_regper)=k-ipart
!            icol_regper(nz_regper)=k_prev-ipart1
!            amat_regper(nz_regper)=deleta_factor
            nz_regper=nz_regper+1
!            irow_regper(nz_regper)=k-ipart
!            icol_regper(nz_regper)=k-npassing-1-ipart1
!            amat_regper(nz_regper)=deleta_factor
          enddo
!
        enddo
      endif
!
    endif
!
  enddo
!
  do istep=ibeg,iend-1
    npassing_prev=npl(istep+1)
    npassing=npl(istep)
!    delphim1=1.d0/delt_neg(istep)
!    deloneovb=0.5d0*delphim1*(1.d0/bhat_mfl(istep+1)-1.d0/bhat_mfl(istep))
!
    do m=0,lag
      k_prev=ind_start(istep+1)+2*(npassing_prev+1)*m+2*npassing_prev+3
      k=ind_start(istep)+2*(npassing+1)*m+2*npassing+3
!
! free flight:
!
      do ipart=1,npassing+1
        nz=nz+1
!        irow(nz)=k-ipart
!        icol(nz)=k-ipart
!        amat_sp(nz)=delphim1
      enddo
!
      do ipart=1,npassing
        nz=nz+1
!        irow(nz)=k-ipart
!        icol(nz)=k_prev-ipart
!        amat_sp(nz)=-delphim1
      enddo
!
      if(npassing_prev.ge.npassing) then
        nz=nz+1
!        irow(nz)=k-npassing-1
!        icol(nz)=k_prev-npassing-1
!        amat_sp(nz)=-delphim1
      endif
!
! mirroring:
!
      if(npassing_prev.eq.npassing) then
!
        do kk=1,4
          nz=nz+1
!          irow(nz)=k-npassing-1
!          icol(nz)=k-npassing-kk+1
!          amat_sp(nz)=deloneovb*rhs_mat_fzero(kk,istep,0)
        enddo
!
        do kk=1,4
          nz=nz+1
!          irow(nz)=k-npassing-1
!          icol(nz)=k_prev-npassing_prev-kk+1
!          amat_sp(nz)=deloneovb*rhs_mat_fzero(kk,istep+1,0)
        enddo
!
      elseif(npassing_prev.gt.npassing) then
        nz=nz+1
!        irow(nz)=k_prev-npassing_prev-2
!        icol(nz)=k_prev-npassing_prev-1
!        amat_sp(nz)=-delphim1
      endif
!
! collisions:
!
      if(fact_neg_e(istep).ne.0.d0) then
!
! Lorentz operator:
!
        if(isw_lor.eq.1) then
!
          do ipart=1,npassing+1
            do kk=1,5
              do mm=0,lag
                nz=nz+1
!                irow(nz)=k-ipart
!                icol(nz)=k-max(0,ipart-3)-kk+2*(npassing+1)*(mm-m)
!                amat_sp(nz)=anumm(m,mm)*rhs_mat_lorentz(kk,ipart,istep) &
!                           *fact_neg_e(istep)
              enddo
            enddo
          enddo
!
! matching collisional flux through the boundary in backward Euler scheme:
!
          do kk=1,4
            do mm=0,lag
              nz=nz+1
!              irow(nz)=k-npassing-1
!              icol(nz)=k-npassing+1-kk+2*(npassing+1)*(mm-m)
!              amat_sp(nz)=-0.5d0*anumm(m,mm)*rhs_mat_fzero(kk,istep,1) &
!                         *fact_neg_e(istep)
            enddo
          enddo
!
          do kk=1,4
            do mm=0,lag
              nz=nz+1
!              irow(nz)=k-npassing-1
!              icol(nz)=k_prev-npassing_prev+1-kk+2*(npassing_prev+1)*(mm-m)
!              amat_sp(nz)=0.5d0*anumm(m,mm)*rhs_mat_fzero(kk,istep+1,1) &
!                         *fact_neg_b(istep+1)
            enddo
          enddo
!
        endif
!
!        nz_beg=nz+1
!
! energy diffusion operator:
!
        if(isw_ene.eq.1) then
!
          do ipart=1,npassing
            do kk=1,4
              do mm=0,lag
                nz=nz+1
!                irow(nz)=k-ipart
!                icol(nz)=k-max(0,ipart-2)-kk+2*(npassing+1)*(mm-m)
!                amat_sp(nz)=denmm(m,mm)*rhs_mat_energ(kk,ipart,istep)
              enddo
            enddo
          enddo
!
          do kk=1,4
            do mm=0,lag
              nz=nz+1
!              irow(nz)=k-npassing-1
!              icol(nz)=k-npassing+1-kk+2*(npassing+1)*(mm-m)
!              amat_sp(nz)=0.5d0*denmm(m,mm)                                  &
!                         *rhs_mat_energ(kk,npassing+1,istep)
            enddo
          enddo
!
          do kk=1,4
            do mm=0,lag
              nz=nz+1
!              irow(nz)=k-npassing-1
!              icol(nz)=k_prev-npassing_prev+1-kk+2*(npassing_prev+1)*(mm-m)
!              amat_sp(nz)=0.5d0*denmm(m,mm)                                  &
!                         *rhs_mat_energ(kk,npassing_prev+1,istep+1)
            enddo
          enddo
!
        endif
!
!        amat_sp(nz_beg:nz)=fact_neg_e(istep)*amat_sp(nz_beg:nz)
      endif
!
    enddo
!
  enddo
!
!
  allocate(irow(nz),icol(nz),amat_sp(nz),ipcol(ncol),bvec_sp(ncol))
  allocate(irow_regper(nz_regper),icol_regper(nz_regper),amat_regper(nz_regper))
  if(isw_intp.eq.1) allocate(bvec_iter(ncol),bvec_lor(ncol),bvec_prev(ncol))
!
! Fill the arrays:
!
  nz=0
  nz_regper=0
!
! Co-passing: sigma=1
!
  istep=ibeg
  npassing=npl(istep)
!
! entry:
!
  do m=0,lag
    k=ind_start(istep)+2*(npassing+1)*m
!
    do ipart=1,npassing+1
      nz=nz+1
      irow(nz)=k+ipart
      icol(nz)=k+ipart
      amat_sp(nz)=1.d0
    enddo
!
    if(isw_axisymm.eq.1) then
      k_prev=ind_start(iend)+2*(npassing+1)*m
!
      do ipart=1,npassing+1
        nz=nz+1
        irow(nz)=k+ipart
        icol(nz)=k_prev+ipart
        amat_sp(nz)=-1.d0
      enddo
!
      if(isw_regper.eq.1.and.m.lt.1) then
        do ipart=1,npassing+1
          if(ipart.le.npassing) then
            deleta_factor=(eta(ipart)-eta(ipart-1))*bhat_mfl(istep)
          else
            deleta_factor=1.d0-eta(ipart-1)*bhat_mfl(istep)
          endif
!
          do ipart1=1,npassing+1
            nz_regper=nz_regper+1
            irow_regper(nz_regper)=k+ipart
            icol_regper(nz_regper)=k_prev+ipart1
            amat_regper(nz_regper)=deleta_factor
            nz_regper=nz_regper+1
            irow_regper(nz_regper)=k+ipart
            icol_regper(nz_regper)=k+npassing+1+ipart1
            amat_regper(nz_regper)=deleta_factor
          enddo
!
        enddo
      endif
!
    endif
!
  enddo
!
  do istep=ibeg+1,iend
    npassing_prev=npl(istep-1)
    npassing=npl(istep)
    delphim1=1.d0/delt_pos(istep)
    deloneovb=0.5d0*delphim1*(1.d0/bhat_mfl(istep-1)-1.d0/bhat_mfl(istep))
!
    do m=0,lag
      k_prev=ind_start(istep-1)+2*(npassing_prev+1)*m
      k=ind_start(istep)+2*(npassing+1)*m
!
! free flight:
!
      do ipart=1,npassing+1
        nz=nz+1
        irow(nz)=k+ipart
        icol(nz)=k+ipart
        amat_sp(nz)=delphim1
      enddo
!
      do ipart=1,npassing
        nz=nz+1
        irow(nz)=k+ipart
        icol(nz)=k_prev+ipart
        amat_sp(nz)=-delphim1
      enddo
!
      if(npassing_prev.ge.npassing) then
        nz=nz+1
        irow(nz)=k+npassing+1
        icol(nz)=k_prev+npassing+1
        amat_sp(nz)=-delphim1
      endif
!
! mirroring:
!
      if(npassing_prev.eq.npassing) then
!
        do kk=1,4
          nz=nz+1
          irow(nz)=k+npassing+1
          icol(nz)=k+npassing+kk-1
          amat_sp(nz)=deloneovb*rhs_mat_fzero(kk,istep,0)
        enddo
!
        do kk=1,4
          nz=nz+1
          irow(nz)=k+npassing+1
          icol(nz)=k_prev+npassing_prev+kk-1
          amat_sp(nz)=deloneovb*rhs_mat_fzero(kk,istep-1,0)
        enddo
!
      elseif(npassing_prev.gt.npassing) then
        nz=nz+1
        irow(nz)=k_prev+npassing_prev+2
        icol(nz)=k_prev+npassing_prev+1
        amat_sp(nz)=-delphim1
      endif
!
! collisions:
!
      if(fact_pos_e(istep).ne.0.d0) then
!
! Lorentz operator:
!
        if(isw_lor.eq.1) then
!
          do ipart=1,npassing+1
            do kk=1,5
              do mm=0,lag
                nz=nz+1
                irow(nz)=k+ipart
                icol(nz)=k+max(0,ipart-3)+kk+2*(npassing+1)*(mm-m)
                amat_sp(nz)=anumm(m,mm)*rhs_mat_lorentz(kk,ipart,istep) &
                           *fact_pos_e(istep)
              enddo
            enddo
          enddo
!
! matching collisional flux through the boundary in backward Euler scheme:
!
          do kk=1,4
            do mm=0,lag
              nz=nz+1
              irow(nz)=k+npassing+1
              icol(nz)=k+npassing-1+kk+2*(npassing+1)*(mm-m)
              amat_sp(nz)=-0.5d0*anumm(m,mm)*rhs_mat_fzero(kk,istep,1)  &
                         *fact_pos_e(istep)
            enddo
          enddo
!
          do kk=1,4
            do mm=0,lag
              nz=nz+1
              irow(nz)=k+npassing+1
              icol(nz)=k_prev+npassing_prev-1+kk+2*(npassing_prev+1)*(mm-m)
              amat_sp(nz)=0.5d0*anumm(m,mm)*rhs_mat_fzero(kk,istep-1,1) &
                         *fact_pos_b(istep-1)
            enddo
          enddo
!
        endif
!
        nz_beg=nz+1
!
! energy diffusion operator:
!
        if(isw_ene.eq.1) then
!
          do ipart=1,npassing
            do kk=1,4
              do mm=0,lag
                nz=nz+1
                irow(nz)=k+ipart
                icol(nz)=k+max(0,ipart-2)+kk+2*(npassing+1)*(mm-m)
                amat_sp(nz)=denmm(m,mm)*rhs_mat_energ(kk,ipart,istep)
              enddo
            enddo
          enddo
!
          do kk=1,4
            do mm=0,lag
              nz=nz+1
              irow(nz)=k+npassing+1
              icol(nz)=k+npassing-1+kk+2*(npassing+1)*(mm-m)
              amat_sp(nz)=0.5d0*denmm(m,mm)                                  &
                         *rhs_mat_energ(kk,npassing+1,istep)
            enddo
          enddo
!
          do kk=1,4
            do mm=0,lag
              nz=nz+1
              irow(nz)=k+npassing+1
              icol(nz)=k_prev+npassing_prev-1+kk+2*(npassing_prev+1)*(mm-m)
              amat_sp(nz)=0.5d0*denmm(m,mm)                                  &
                         *rhs_mat_energ(kk,npassing_prev+1,istep-1)
            enddo
          enddo
!
        endif
!
        amat_sp(nz_beg:nz)=fact_pos_e(istep)*amat_sp(nz_beg:nz)
      endif
!
    enddo
!
  enddo
!
! Counter-passing: sigma=-1
!
  istep=iend
  npassing=npl(istep)
!
! entry:
!
  do m=0,lag
    k=ind_start(istep)+2*(npassing+1)*m+2*npassing+3
!
    do ipart=1,npassing+1
      nz=nz+1
      irow(nz)=k-ipart
      icol(nz)=k-ipart
      amat_sp(nz)=1.d0
    enddo
!
    if(isw_axisymm.eq.1) then
      k_prev=ind_start(ibeg)+2*(npassing+1)*m+2*npassing+3
!
      do ipart=1,npassing+1
        nz=nz+1
        irow(nz)=k-ipart
        icol(nz)=k_prev-ipart
        amat_sp(nz)=-1.d0
      enddo
!
      if(isw_regper.eq.1.and.m.lt.1) then
        do ipart=1,npassing+1
          if(ipart.le.npassing) then
            deleta_factor=(eta(ipart)-eta(ipart-1))*bhat_mfl(istep)
          else
            deleta_factor=1.d0-eta(ipart-1)*bhat_mfl(istep)
          endif
!
          do ipart1=1,npassing+1
            nz_regper=nz_regper+1
            irow_regper(nz_regper)=k-ipart
            icol_regper(nz_regper)=k_prev-ipart1
            amat_regper(nz_regper)=deleta_factor
            nz_regper=nz_regper+1
            irow_regper(nz_regper)=k-ipart
            icol_regper(nz_regper)=k-npassing-1-ipart1
            amat_regper(nz_regper)=deleta_factor
          enddo
!
        enddo
      endif
!
    endif
!
  enddo
!
  do istep=ibeg,iend-1
    npassing_prev=npl(istep+1)
    npassing=npl(istep)
    delphim1=1.d0/delt_neg(istep)
    deloneovb=0.5d0*delphim1*(1.d0/bhat_mfl(istep+1)-1.d0/bhat_mfl(istep))
!
    do m=0,lag
      k_prev=ind_start(istep+1)+2*(npassing_prev+1)*m+2*npassing_prev+3
      k=ind_start(istep)+2*(npassing+1)*m+2*npassing+3
!
! free flight:
!
      do ipart=1,npassing+1
        nz=nz+1
        irow(nz)=k-ipart
        icol(nz)=k-ipart
        amat_sp(nz)=delphim1
      enddo
!
      do ipart=1,npassing
        nz=nz+1
        irow(nz)=k-ipart
        icol(nz)=k_prev-ipart
        amat_sp(nz)=-delphim1
      enddo
!
      if(npassing_prev.ge.npassing) then
        nz=nz+1
        irow(nz)=k-npassing-1
        icol(nz)=k_prev-npassing-1
        amat_sp(nz)=-delphim1
      endif
!
! mirroring:
!
      if(npassing_prev.eq.npassing) then
!
        do kk=1,4
          nz=nz+1
          irow(nz)=k-npassing-1
          icol(nz)=k-npassing-kk+1
          amat_sp(nz)=deloneovb*rhs_mat_fzero(kk,istep,0)
        enddo
!
        do kk=1,4
          nz=nz+1
          irow(nz)=k-npassing-1
          icol(nz)=k_prev-npassing_prev-kk+1
          amat_sp(nz)=deloneovb*rhs_mat_fzero(kk,istep+1,0)
        enddo
!
      elseif(npassing_prev.gt.npassing) then
        nz=nz+1
        irow(nz)=k_prev-npassing_prev-2
        icol(nz)=k_prev-npassing_prev-1
        amat_sp(nz)=-delphim1
      endif
!
! collisions:
!
      if(fact_neg_e(istep).ne.0.d0) then
!
! Lorentz operator:
!
        if(isw_lor.eq.1) then
!
          do ipart=1,npassing+1
            do kk=1,5
              do mm=0,lag
                nz=nz+1
                irow(nz)=k-ipart
                icol(nz)=k-max(0,ipart-3)-kk+2*(npassing+1)*(mm-m)
                amat_sp(nz)=anumm(m,mm)*rhs_mat_lorentz(kk,ipart,istep) &
                           *fact_neg_e(istep)
              enddo
            enddo
          enddo
!
! matching collisional flux through the boundary in backward Euler scheme:
!
          do kk=1,4
            do mm=0,lag
              nz=nz+1
              irow(nz)=k-npassing-1
              icol(nz)=k-npassing+1-kk+2*(npassing+1)*(mm-m)
              amat_sp(nz)=-0.5d0*anumm(m,mm)*rhs_mat_fzero(kk,istep,1)  &
                         *fact_neg_e(istep)
            enddo
          enddo
!
          do kk=1,4
            do mm=0,lag
              nz=nz+1
              irow(nz)=k-npassing-1
              icol(nz)=k_prev-npassing_prev+1-kk+2*(npassing_prev+1)*(mm-m)
              amat_sp(nz)=0.5d0*anumm(m,mm)*rhs_mat_fzero(kk,istep+1,1) &
                         *fact_neg_b(istep+1)
            enddo
          enddo
!
        endif
!
        nz_beg=nz+1
!
! energy diffusion operator:
!
        if(isw_ene.eq.1) then
!
          do ipart=1,npassing
            do kk=1,4
              do mm=0,lag
                nz=nz+1
                irow(nz)=k-ipart
                icol(nz)=k-max(0,ipart-2)-kk+2*(npassing+1)*(mm-m)
                amat_sp(nz)=denmm(m,mm)*rhs_mat_energ(kk,ipart,istep)
              enddo
            enddo
          enddo
!
          do kk=1,4
            do mm=0,lag
              nz=nz+1
              irow(nz)=k-npassing-1
              icol(nz)=k-npassing+1-kk+2*(npassing+1)*(mm-m)
              amat_sp(nz)=0.5d0*denmm(m,mm)                                  &
                         *rhs_mat_energ(kk,npassing+1,istep)
            enddo
          enddo
!
          do kk=1,4
            do mm=0,lag
              nz=nz+1
              irow(nz)=k-npassing-1
              icol(nz)=k_prev-npassing_prev+1-kk+2*(npassing_prev+1)*(mm-m)
              amat_sp(nz)=0.5d0*denmm(m,mm)                                  &
                         *rhs_mat_energ(kk,npassing_prev+1,istep+1)
            enddo
          enddo
!
        endif
!
        amat_sp(nz_beg:nz)=fact_neg_e(istep)*amat_sp(nz_beg:nz)
      endif
!
    enddo
!
  enddo
!
! Save the symmetric matrix:
!
  nz_symm=nz
  ALLOCATE(irow_symm(nz_symm),icol_symm(nz_symm),amat_symm(nz_symm))
  irow_symm=irow
  icol_symm=icol
  amat_symm=amat_sp
  DEALLOCATE(irow,icol,amat_sp)
!
! End save symmetric matrix
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
! Solve the linear equation set:
!
  call  remap_rc(nz,nz_sq,irow,icol,amat_sp)
!
  print *,'system size = ',n_2d_size
  print *,'non-zeros before and after truncation = ',nz,nz_sq
  nz=nz_sq
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
  iopt=1
!
  call cpu_time(time_start)
  CALL sparse_solve(nrow,ncol,nz,irow(1:nz),ipcol,amat_sp(1:nz),       &
                    bvec_sp,iopt)
!
  call cpu_time(time_factorization)
  print *,'factorization completed ',time_factorization - time_start,' sec'
!
  iopt=2
!
! Solution of inhomogeneus equation (account of sources):
!
call cpu_time(time1)

!! Modification by Andreas F. Martitsch (23.08.2015)
!  multi-species part (allocate storage for source_vector)
IF(ALLOCATED(source_vector_all)) DEALLOCATE(source_vector_all)
ALLOCATE(source_vector_all(n_2d_size,1:3,0:num_spec-1))
source_vector_all=0.0d0
! save solution of the differential part for species=ispec
! (diffusion coeff. driven by thermodyn. forces of other
! species are zero -> interaction through integral part)
source_vector_all(:,1:3,ispec)=source_vector(:,1:3)
! solve system
DO ispecp=0,num_spec-1
  IF(ispecp .NE. ispec) CYCLE
  CALL sparse_solve(nrow,ncol,nz,irow(1:nz),ipcol,amat_sp(1:nz),       &
                    source_vector_all(:,1:3,ispecp),iopt)
ENDDO
!! End Modification by Andreas F. Martitsch (23.08.2015)
CALL cpu_TIME(time2)
!write (*,*) "Time in first solver:", time2-time1
!
! integral part:
!
call cpu_time(time1)
  if(isw_intp.eq.1) then
!
    do k=1,3
!
      do ispecp=0,num_spec-1
        print *,'species',ispecp,':'
        bvec_lor=source_vector_all(:,k,ispecp)
!
        do iter=1,niter
          bvec_prev=source_vector_all(:,k,ispecp)
!
          call integral_part(npart,leg,lag,ibeg,iend,n_2d_size,npl,ind_start,   &
                             phi_mfl,pleg_bra(0:leg,:,:),pleg_ket(0:leg,:,:),   &
                             ailmm,source_vector_all(:,k,ispecp),bvec_iter)
!
          CALL sparse_solve(nrow,ncol,nz,irow(1:nz),ipcol,amat_sp(1:nz),        &
                            bvec_iter,iopt)

!
          source_vector_all(:,k,ispecp)=bvec_lor+bvec_iter
          !! Modification by Andreas F. Martitsch (20.08.2015)
          ! MPI Barrier -> Exchange exit conditions between
          ! different processes
          break_cond1(ispec)=SUM(ABS(source_vector_all(:,k,ispecp)-bvec_prev))
          break_cond2(ispec)=SUM(ABS(bvec_prev))*epserr_iter
          PRINT *,iter,break_cond1(ispec),break_cond2(ispec)
          CALL mpro%allgather_inplace(break_cond1)
          CALL mpro%allgather_inplace(break_cond2)
          IF(ALL(break_cond1 .LT. break_cond2)) EXIT
          !! End Modification by Andreas F. Martitsch (20.08.2015)
        enddo
        if (niter .eq. iter) write (*,*) "Maximum number of iterations reached in ripple solver."
!
      enddo
!
    enddo
!
  endif
  call cpu_time(time2)
  !write (*,*) "Time in integral part:", time2 - time1

!
! Plotting:
!
  if(iplot.eq.1) then
    nphiplot=200    !serves as an upper limit for data output
    delphiplot=maxval(phi_mfl(ibeg+1:iend)-phi_mfl(ibeg:iend-1))
    if(delphiplot.gt.epsilon(1.d0)) then
      nphiequi=(phi_mfl(iend)-phi_mfl(ibeg))/delphiplot
      nphiequi=max(1,nphiequi)
    else
      nphiequi=1
    endif
    delphiplot=(phi_mfl(iend)-phi_mfl(ibeg))/min(nphiplot,nphiequi)
    nplp1=npart_loc+1
    ALLOCATE(fun_write(0:lag,0:3,0:nplp1,3))
    icounter=0
    phiplot=phi_mfl(ibeg)-1.d0
!
    write(propname,*) fieldpropagator%tag

    if (prop_fileformat .eq. 1) then

       !call h5_create('spitzer_' // trim(adjustl(propname)) // '.h5', h5id_final_spitzer)

       ! Create unlimited arrays in HDF5 file
       call h5_create('phi_mesh_' // trim(adjustl(propname)) // '.h5', h5id_phi_mesh)
       !call h5_define_group(h5id_final_spitzer, 'phi_mesh', h5id_phi_mesh)

       call h5_create('dentf_' // trim(adjustl(propname)) // '.h5', h5id_dentf)
       !call h5_define_group(h5id_final, 'dentf', h5id_dentf)

       call h5_create('enetf_' // trim(adjustl(propname)) // '.h5', h5id_enetf)
       !call h5_define_group(h5id_final, 'enetf', h5id_enetf)

       call h5_create('spitf_' // trim(adjustl(propname)) // '.h5', h5id_spitf)
       !call h5_define_group(h5id_final_spitzer, 'spitf', h5id_spitf)

    else

       ! Create ASCII files
       open(iunit_phi,file='phi_mesh.'                        &
            //trim(adjustl(propname))//'.dat')
       open(iunit_dt_p,form='unformatted',file='dentf_p.'     &
            //trim(adjustl(propname))//'.dat')
       open(iunit_dt_m,form='unformatted',file='dentf_m.'     &
            //trim(adjustl(propname))//'.dat')
       open(iunit_sp_p,form='unformatted',file='spitf_p.'     &
            //trim(adjustl(propname))//'.dat')
       open(iunit_sp_m,form='unformatted',file='spitf_m.'     &
            //trim(adjustl(propname))//'.dat')
       open(iunit_et_p,form='unformatted',file='enetf_p.'     &
            //trim(adjustl(propname))//'.dat')
       open(iunit_et_m,form='unformatted',file='enetf_m.'     &
            //trim(adjustl(propname))//'.dat')
    end if
    !

    !**********************************************************
    ! Allocate space for datasets written to HDF5 at once
    !**********************************************************
    if (allocated(phi_mfl_h5))  deallocate(phi_mfl_h5)
    if (allocated(bhat_mfl_h5)) deallocate(bhat_mfl_h5)
    if (allocated(npassing_h5)) deallocate(npassing_h5)
    if (allocated(dentf_p_h5)) deallocate(dentf_p_h5)
    if (allocated(spitf_p_h5)) deallocate(spitf_p_h5)
    if (allocated(enetf_p_h5)) deallocate(enetf_p_h5)
    if (allocated(dentf_m_h5)) deallocate(dentf_m_h5)
    if (allocated(spitf_m_h5)) deallocate(spitf_m_h5)
    if (allocated(enetf_m_h5)) deallocate(enetf_m_h5)
    allocate(phi_mfl_h5(iend-ibeg+1))
    allocate(bhat_mfl_h5(iend-ibeg+1))
    allocate(npassing_h5(iend-ibeg+1))
    allocate(dentf_p_h5(0:lag, 0:3, 0:nplp1, iend-ibeg+1))
    allocate(spitf_p_h5(0:lag, 0:3, 0:nplp1, iend-ibeg+1))
    allocate(enetf_p_h5(0:lag, 0:3, 0:nplp1, iend-ibeg+1))
    allocate(dentf_m_h5(0:lag, 0:3, 0:nplp1, iend-ibeg+1))
    allocate(spitf_m_h5(0:lag, 0:3, 0:nplp1, iend-ibeg+1))
    allocate(enetf_m_h5(0:lag, 0:3, 0:nplp1, iend-ibeg+1))


    do istep=ibeg,iend
      if(phi_mfl(istep).lt.phiplot.and.istep.ne.iend) cycle
      icounter=icounter+1
      phiplot=phi_mfl(istep)+delphiplot
      npassing=npl(istep)
      eta0=1.d0/bhat_mfl(istep)
      !
      if (prop_fileformat .eq. 1) then
         phi_mfl_h5(icounter)  = phi_mfl(istep)
         bhat_mfl_h5(icounter) = bhat_mfl(istep)
         npassing_h5(icounter) = npassing
      else
         write (iunit_phi,*) phi_mfl(istep),npassing,bhat_mfl(istep)
      end if
      !
      fun_write=0.d0
      do m=0,lag
        k=ind_start(istep)+2*(npassing+1)*m
        fun_write(m,:,1,:)=MATMUL(derivs_plot(:,:,1,istep),             &
                                  source_vector(k+1:k+4,:))
        do i=2,npassing+1
          fun_write(m,:,i,:)=MATMUL(derivs_plot(:,:,i,istep),           &
                                    source_vector(k+i-1:k+i+2,:))
        enddo
      enddo
      !
      if (prop_fileformat .eq. 1) then
         dentf_p_h5(:,:,:,icounter) = fun_write(:,:,:,1)
         spitf_p_h5(:,:,:,icounter) = fun_write(:,:,:,2)/surface_boozer_B00
         enetf_p_h5(:,:,:,icounter) = fun_write(:,:,:,3)
      else
         write(iunit_dt_p) fun_write(:,:,:,1)
         write(iunit_sp_p) fun_write(:,:,:,2)/surface_boozer_B00
         write(iunit_et_p) fun_write(:,:,:,3)
      end if
      !
      fun_write=0.d0
      do m=0,lag
        k=ind_start(istep)+2*(npassing+1)*(m+1)
        fun_write(m,:,1,:)=MATMUL(derivs_plot(:,:,1,istep),             &
                                  source_vector(k:k-3:-1,:))
        do i=2,npassing+1
          fun_write(m,:,i,:)=MATMUL(derivs_plot(:,:,i,istep),           &
                                    source_vector(k-i+2:k-i-1:-1,:))
        enddo
      enddo
      !
      if (prop_fileformat .eq. 1) then
         dentf_m_h5(:,:,:,icounter) = fun_write(:,:,:,1)
         spitf_m_h5(:,:,:,icounter) = fun_write(:,:,:,2)/surface_boozer_B00
         enetf_m_h5(:,:,:,icounter) = fun_write(:,:,:,3)
     else
         write(iunit_dt_m) fun_write(:,:,:,1)
         write(iunit_sp_m) fun_write(:,:,:,2)/surface_boozer_B00
         write(iunit_et_m) fun_write(:,:,:,3)
      end if
      !
    enddo

    if (prop_fileformat .eq. 1) then
       call h5_add(h5id_phi_mesh, 'phi_mfl',  phi_mfl_h5(1:icounter),  &
            lbound(phi_mfl_h5(1:icounter)),  ubound(phi_mfl_h5(1:icounter)))
       call h5_add(h5id_phi_mesh, 'bhat_mfl', bhat_mfl_h5(1:icounter), &
            lbound(bhat_mfl_h5(1:icounter)), ubound(bhat_mfl_h5(1:icounter)))
       call h5_add(h5id_phi_mesh, 'npassing', npassing_h5(1:icounter), &
            lbound(npassing_h5(1:icounter)), ubound(npassing_h5(1:icounter)))

       call h5_add(h5id_dentf, 'dentf_p', dentf_p_h5(:,:,:,1:icounter), &
            lbound(dentf_p_h5(:,:,:,1:icounter)), ubound(dentf_p_h5(:,:,:,1:icounter)))
       call h5_add(h5id_spitf, 'spitf_p', spitf_p_h5(:,:,:,1:icounter), &
            lbound(spitf_p_h5(:,:,:,1:icounter)), ubound(spitf_p_h5(:,:,:,1:icounter)))
       call h5_add(h5id_enetf, 'enetf_p', enetf_p_h5(:,:,:,1:icounter), &
            lbound(enetf_p_h5(:,:,:,1:icounter)), ubound(enetf_p_h5(:,:,:,1:icounter)))

       call h5_add(h5id_dentf, 'dentf_m', dentf_m_h5(:,:,:,1:icounter), &
            lbound(dentf_m_h5(:,:,:,1:icounter)), ubound(dentf_m_h5(:,:,:,1:icounter)))
       call h5_add(h5id_spitf, 'spitf_m', spitf_m_h5(:,:,:,1:icounter), &
            lbound(spitf_m_h5(:,:,:,1:icounter)), ubound(spitf_m_h5(:,:,:,1:icounter)))
       call h5_add(h5id_enetf, 'enetf_m', enetf_m_h5(:,:,:,1:icounter), &
            lbound(enetf_m_h5(:,:,:,1:icounter)), ubound(enetf_m_h5(:,:,:,1:icounter)))

       call h5_close(h5id_enetf)
       !call h5_close_group(h5id_spitf)
       call h5_close(h5id_phi_mesh)
       call h5_close(h5id_dentf)
       call h5_close(h5id_spitf)

       deallocate(phi_mfl_h5, bhat_mfl_h5, npassing_h5)
       deallocate(dentf_p_h5, spitf_p_h5, enetf_p_h5)
       deallocate(dentf_m_h5, spitf_m_h5, enetf_m_h5)
    else
       close(iunit_phi)
       close(iunit_dt_p)
       close(iunit_dt_m)
       close(iunit_sp_p)
       close(iunit_sp_m)
       close(iunit_et_p)
       close(iunit_et_m)
    end if
    !
    if (prop_fileformat .eq. 1) then
       call h5_create('sizeplot_etalev_' // trim(adjustl(propname)) // '.h5', h5id_sizeplot)
       !call h5_define_group(h5id_final_spitzer, 'sizeplot_etalev', h5id_sizeplot)
       call h5_add(h5id_sizeplot, 'lag', lag)
       call h5_add(h5id_sizeplot, 'nplp1', nplp1)
       call h5_add(h5id_sizeplot, 'icounter', icounter)
       call h5_add(h5id_sizeplot, 'collpar', collpar)
       call h5_add(h5id_sizeplot, 'travis_convfac', travis_convfac )
       call h5_add(h5id_sizeplot, 'eta', eta(0:nplp1), lbound(eta(0:nplp1)), ubound(eta(0:nplp1)))
       !call h5_add(h5id_sizeplot, 'eta_all', eta, lbound(eta), ubound(eta))
       !call h5_add(h5id_sizeplot, 'eta_glob', eta_glob, lbound(eta_glob), ubound(eta_glob))
       !call h5_add(h5id_sizeplot, 'eta_loc', eta_loc, lbound(eta_loc), ubound(eta_loc))
       !call h5_close_group(h5id_sizeplot)
       call h5_add(h5id_sizeplot, 'ripple_tag',fieldpropagator%ch_act%tag)
       call h5_add(h5id_sizeplot, 'ripple_b_max_l',fieldpropagator%ch_act%b_max_l)
       call h5_add(h5id_sizeplot, 'ripple_b_max_r',fieldpropagator%ch_act%b_max_r)
       call h5_close(h5id_sizeplot)
    else
       open(iunit_sizes,file='sizeplot_etalev.'               &
            //trim(adjustl(propname))//'.dat')
       write(iunit_sizes,*) lag,nplp1,icounter,collpar,travis_convfac
       write(iunit_sizes,*) eta(0:nplp1)
       close(iunit_sizes)
    end if
!
    DEALLOCATE(fun_write)
!
  endif
!
call cpu_time(time1)
!! Modification by Andreas F. Martitsch (23.08.2015)
! old behavior (for a single species)
!qflux=MATMUL(flux_vector,source_vector)
!  multi-species part
IF(ALLOCATED(qflux_allspec)) DEALLOCATE(qflux_allspec)
ALLOCATE(qflux_allspec(1:3,1:3,0:num_spec-1,0:num_spec-1))
qflux_allspec=0.0d0
DO ispecp=0,num_spec-1
   qflux=MATMUL(flux_vector,source_vector_all(:,:,ispecp))
   qflux_allspec(:,:,ispecp,ispec)=qflux
ENDDO
! order of species inidices (ispecp,ispec) interchanged
! (-> easier to handle within mpro%allgather)
CALL mpro%allgather_inplace(qflux_allspec)
! go back to the "natural" order of species indices (ispec,ispecp)
IF(ALLOCATED(qflux_allspec_tmp)) DEALLOCATE(qflux_allspec_tmp)
ALLOCATE(qflux_allspec_tmp(1:3,1:3,0:num_spec-1,0:num_spec-1))
qflux_allspec_tmp=0.0d0
DO ispecp=0,num_spec-1
   DO ispecpp=0,num_spec-1
      qflux_allspec_tmp(:,:,ispecp,ispecpp)=qflux_allspec(:,:,ispecpp,ispecp)
   ENDDO
ENDDO
qflux_allspec=qflux_allspec_tmp
IF(ALLOCATED(qflux_allspec_tmp)) DEALLOCATE(qflux_allspec_tmp)
IF(mpro%getrank() .EQ. 0) THEN
  ! D11
  PRINT *,'qflux(1,1,0,0):'
  PRINT *,qflux_allspec(1,1,0,0)
  PRINT *,'qflux(1,1,1,0):'
  PRINT *,qflux_allspec(1,1,1,0)
  PRINT *,'qflux(1,1,0,1):'
  PRINT *,qflux_allspec(1,1,0,1)
  PRINT *,'qflux(1,1,1,1):'
  PRINT *,qflux_allspec(1,1,1,1)
  ! D12
  PRINT *,'qflux(1,3,0,0):'
  PRINT *,qflux_allspec(1,3,0,0)
  PRINT *,'qflux(1,3,1,0):'
  PRINT *,qflux_allspec(1,3,1,0)
  PRINT *,'qflux(1,3,0,1):'
  PRINT *,qflux_allspec(1,3,0,1)
  PRINT *,'qflux(1,3,1,1):'
  PRINT *,qflux_allspec(1,3,1,1)
  OPEN(070915,file='qflux_symm_allspec.dat')
  WRITE(070915,*) boozer_s, collpar, &
        qflux_allspec(1,1,0,0), qflux_allspec(1,1,1,0), &
        qflux_allspec(1,1,0,1), qflux_allspec(1,1,1,1), &
        qflux_allspec(1,3,0,0), qflux_allspec(1,3,1,0), &
        qflux_allspec(1,3,0,1), qflux_allspec(1,3,1,1)
  CLOSE(070915)
  !STOP
END IF
!! End Modification by Andreas F. Martitsch (23.08.2015)
call cpu_time(time2)
RETURN
!write (*,*) "Time in matmul(): ", time2-time1, iplot
!call gemm(qflux, flux_vector, source_vector)
!write (*,*) flux_vector
!write (*,*) source_vector
!
  do m=0,lag
    do kk=1,3
      k=ind_start(ibeg)+2*npass_l*m
      source_m(npass_l*m+1:npass_l*m+npass_l,kk) &
          =source_vector(k+2*npass_l:k+npass_l+1:-1,kk)
      k=ind_start(iend)+2*npass_r*m
      source_p(npass_r*m+1:npass_r*m+npass_r,kk) &
          =source_vector(k+1:k+npass_r,kk)
    enddo
  enddo
!
  if(iplot.eq.1.or.isw_axisymm.eq.1) then
!
    flux_p=0.d0
    flux_m=0.d0
    amat_plus_plus=0.d0
    amat_plus_minus=0.d0
    amat_minus_minus=0.d0
    amat_minus_plus=0.d0
!
    deallocate(flux_vector,source_vector,irow,icol,amat_sp,ipcol,bvec_sp)
    if(isw_intp.eq.1) deallocate(bvec_iter,bvec_lor,bvec_prev)
    DEALLOCATE(deriv_coef,enu_coef,alambd,Vg_vp_over_B,scalprod_pleg)
    DEALLOCATE(alampow,vrecurr,dellampow,convol_polpow,pleg_bra,pleg_ket)
    DEALLOCATE(npl,rhs_mat_fzero,rhs_mat_lorentz,rhs_mat_energ,q_rip)
    DEALLOCATE(convol_flux,convol_curr,ind_start)
    DEALLOCATE(phi_mfl,bhat_mfl,geodcu_mfl,h_phi_mfl,dlogbdphi_mfl,eta)
    DEALLOCATE(delt_pos,delt_neg,fact_pos_b,fact_neg_b,fact_pos_e,fact_neg_e)
    IF (ALLOCATED(eta_prev)) DEALLOCATE(eta_prev)
    IF (ALLOCATED(eta_next)) DEALLOCATE(eta_next)
!
    call cpu_time(time_solver)
    print *,'solving (1) completed       ',time_solver - time_factorization,' sec'
!
    return
!
  endif
!
! Calculation of propagators:
!
time3 = 0
call cpu_time(time1)
  do m=0,lag
    k=ind_start(ibeg)+2*npass_l*m
    do i=1,npass_l
      bvec_sp=0.d0
      bvec_sp(k+i)=1.d0
!
call cpu_time(time4)
      CALL sparse_solve(nrow,ncol,nz,irow(1:nz),ipcol,amat_sp(1:nz),   &
                        bvec_sp,iopt)
call cpu_time(time5)
time3 = time3 + (time5-time4)
!
! integral part:
!
      if(isw_intp.eq.1) then
        bvec_lor=bvec_sp
!
        do iter=1,niter
          bvec_prev=bvec_sp
!
          call integral_part(npart,leg,lag,ibeg,iend,n_2d_size,npl,ind_start, &
                             phi_mfl,pleg_bra(0:leg,:,:),pleg_ket(0:leg,:,:), &
                             ailmm,bvec_sp,bvec_iter)
!
          call cpu_time(time4)
          CALL sparse_solve(nrow,ncol,nz,irow(1:nz),ipcol,amat_sp(1:nz),      &
                            bvec_iter,iopt)
          call cpu_time(time5)
          time3 = time3 + (time5-time4)
!
          bvec_sp=bvec_lor+bvec_iter
          !write (*,*) "2 ", iter, sum(abs(bvec_sp-bvec_prev)), sum(abs(bvec_prev))*epserr_iter

          if(sum(abs(bvec_sp-bvec_prev)) .lt.                                 &
             sum(abs(bvec_prev))*epserr_iter) then
   	     !write (*,*) "Number of iterations: ", iter
            exit
          endif
!
        enddo
        !
        if (niter .eq. iter) write (*,*) "Maximum number of iterations reached in ripple solver."
      endif
!
      do mm=0,lag
        kk=ind_start(iend)+2*npass_r*mm
        amat_plus_plus(npass_r*mm+1:npass_r*mm+npass_r,npass_l*m+i)    &
                     =bvec_sp(kk+1:kk+npass_r)
      enddo
      do mm=0,lag
        kk=ind_start(ibeg)+2*npass_l*mm
        amat_plus_minus(npass_l*mm+1:npass_l*mm+npass_l,npass_l*m+i)   &
                     =bvec_sp(kk+2*npass_l:kk+npass_l+1:-1)
      enddo
      flux_p(:,npass_l*m+i)=matmul(flux_vector,bvec_sp(:))
    enddo
  enddo
!
  do m=0,lag
    k=ind_start(iend)+2*npass_r*m
    do i=1,npass_r
      bvec_sp=0.d0
      bvec_sp(k+2*npass_r+1-i)=1.d0
!
      call cpu_time(time4)
      CALL sparse_solve(nrow,ncol,nz,irow(1:nz),ipcol,amat_sp(1:nz),   &
                        bvec_sp,iopt)
      call cpu_time(time5)
      time3 = time3 + (time5-time4)
!
! integral part:
!
      if(isw_intp.eq.1) then
        bvec_lor=bvec_sp
!
        do iter=1,niter
          bvec_prev=bvec_sp
!
          call integral_part(npart,leg,lag,ibeg,iend,n_2d_size,npl,ind_start, &
                             phi_mfl,pleg_bra(0:leg,:,:),pleg_ket(0:leg,:,:), &
                             ailmm,bvec_sp,bvec_iter)
!
          call cpu_time(time4)
          CALL sparse_solve(nrow,ncol,nz,irow(1:nz),ipcol,amat_sp(1:nz),      &
                            bvec_iter,iopt)
          call cpu_time(time5)
          time3 = time3 + (time5-time4)
!
          bvec_sp=bvec_lor+bvec_iter
          !write (*,*) "3 ", iter, sum(abs(bvec_sp-bvec_prev)), sum(abs(bvec_prev))*epserr_iter

          if(sum(abs(bvec_sp-bvec_prev)) .lt.                                 &
             sum(abs(bvec_prev))*epserr_iter) then
             !write (*,*) "Number of iterations: ", iter
             exit
          endif
!
        enddo
        !
        if (niter .eq. iter) write (*,*) "Maximum number of iterations reached in ripple solver."
      endif
!
      do mm=0,lag
        kk=ind_start(iend)+2*npass_r*mm
        amat_minus_plus(npass_r*mm+1:npass_r*mm+npass_r,npass_r*m+i)   &
                     =bvec_sp(kk+1:kk+npass_r)
      enddo
      do mm=0,lag
        kk=ind_start(ibeg)+2*npass_l*mm
        amat_minus_minus(npass_l*mm+1:npass_l*mm+npass_l,npass_r*m+i)  &
                     =bvec_sp(kk+2*npass_l:kk+npass_l+1:-1)
      enddo
      flux_m(:,npass_r*m+i)=matmul(flux_vector,bvec_sp(:))
    enddo
  enddo
!
  iopt=3
!
call cpu_time(time2)
write (*,*) "Time in solver in calculating propagator: ", time3
write (*,*) "Time in calculating propagator: ", time2 - time1


call cpu_time(time1)
  CALL sparse_solve(nrow,ncol,nz,irow(1:nz),ipcol,amat_sp(1:nz),bvec_sp,iopt)
call cpu_time(time2)
!write (*,*) "Time in solver:", time2 - time1
!
  deallocate(flux_vector,source_vector,irow,icol,amat_sp,ipcol,bvec_sp)
  if(isw_intp.eq.1) deallocate(bvec_iter,bvec_lor,bvec_prev)
!
!
!
  call cpu_time(time_solver)
  print *,'solving completed (2)       ',time_solver - time_factorization,' sec'
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

#if !defined(MPI_SUPPORT)
OPEN(112,file='onsager.dat')
OPEN(113,file='levhist.dat')
OPEN(111,file='qfluxhist.dat')
ELSE
OPEN(112,file='onsager.dat',position='append')
OPEN(113,file='levhist.dat',position='append')
OPEN(111,file='qfluxhist.dat',position='append')
#endif
ENDIF
facnorm_m=1.d0
facnorm_p=1.d0
!WRITE(112,*) fieldpropagator%tag,qflux(1:2,1),qflux(1:2,2),npart_loc &
!             ,b_max_l,b_max_r,facnorm_p,facnorm_m,ignore_boundary_layer_new &
!             ,ignore_boundary_layer,ifilter_l,ifilter_r
#if !defined(MPI_SUPPORT)
WRITE(112,*) fieldpropagator%tag,qflux(1:2,1),qflux(1:2,2),npart_loc &
             ,b_max_l,b_max_r,facnorm_p,facnorm_m,ignore_lb &
             ,ignore_rb,modify_bl,modify_br
!             ,ignore_rb,ifilter_l,ifilter_r
write(111,*) fieldpropagator%tag,qflux
DO i=0,npart_loc
WRITE(113,*) fieldpropagator%tag,eta(i)
ENDDO
CLOSE(112)
CLOSE(113)
CLOSE(111)
#endif

PRINT *,qflux(1:2,1),qflux(1:2,2)
!!PAUSE
!
GOTO 1
!
#if !defined(MPI_SUPPORT)
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
#endif
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
!PAUSE 'written' ! Warning in gfortran-4.7
1 CONTINUE
!
PRINT *,' '
  DEALLOCATE(deriv_coef,enu_coef,alambd,Vg_vp_over_B,scalprod_pleg)
  DEALLOCATE(alampow,vrecurr,dellampow,convol_polpow,pleg_bra,pleg_ket)
  DEALLOCATE(npl,rhs_mat_fzero,rhs_mat_lorentz,rhs_mat_energ,q_rip)
  DEALLOCATE(convol_flux,convol_curr,ind_start)
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
  IF (ALLOCATED(eta)) DEALLOCATE(eta)
  IF (ALLOCATED(eta_prev)) DEALLOCATE(eta_prev)
  IF (ALLOCATED(eta_next)) DEALLOCATE(eta_next)

  !------------------------------------------------------------------------
  RETURN
END SUBROUTINE ripple_solver
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
SUBROUTINE integral_part(npart,leg,lag,ibeg,iend,n_2d_size,npl,ind_start, &
                         phi_mfl,pleg_bra,pleg_ket,ailmm,vec_in,vec_out)
!
  !! Modification by Andreas F. Martitsch (28.07.2015)
  ! MPI SUPPORT for multi-species part
  ! (run with, e.g.,  mpiexec -np 3 ./neo2.x)
  USE mpiprovider_module
  USE collop, ONLY : ailmm_aa, num_spec
  USE hdf5_tools
  !! End Modification by Andreas F. Martitsch (28.07.2015)
!
  IMPLICIT NONE
!
  INTEGER :: npart,leg,lag,ibeg,iend,n_2d_size
  INTEGER :: l,m,i,k,istep,npassing
!
  INTEGER,          DIMENSION(ibeg:iend)                 :: npl,ind_start
  DOUBLE PRECISION, DIMENSION(ibeg:iend)                 :: phi_mfl
  DOUBLE PRECISION, DIMENSION(n_2d_size)                 :: vec_in,vec_out
  DOUBLE PRECISION, DIMENSION(0:leg,1:npart+1,ibeg:iend) :: pleg_bra,pleg_ket
  DOUBLE PRECISION, DIMENSION(0:lag,0:lag,0:leg)         :: ailmm
  !! Modification by Andreas F. Martitsch (20.08.2015)
  ! Array extended by 3rd (phi-steps) and 4th dimension (species)
  !DOUBLE PRECISION, DIMENSION(0:lag,0:leg)               :: scalprod_pleg
  DOUBLE PRECISION, DIMENSION(0:lag,0:leg,ibeg:iend,0:num_spec-1) :: scalprod_pleg
  DOUBLE PRECISION, DIMENSION(0:lag,0:leg,ibeg:iend,0:num_spec-1) :: scalprod_pleg_tmp
  ! Species index
  INTEGER :: ispec, ispecp
  INTEGER :: h5id_scalprod_pleg
  !! End Modification by Andreas F. Martitsch (20.08.2015)
!
  !! Modification by Andreas F. Martitsch (28.07.2015)
  ! multi-species part - MPI rank determines species
  ispec = mpro%getRank()
  !PRINT *,"Species: ", ispec
  !CALL collop_set_species(ispec)
  !PRINT *,'asource: ',asource(:,1)
  !PRINT *,'anumm:',anumm(1,:)
  !PRINT *,'denmm:',denmm(1,:)
  !STOP
  !! End Modification by Andreas F. Martitsch (28.07.2015)
!
  vec_out=0.d0
!
  DO istep=ibeg,iend
!
    npassing=npl(istep)
!
    DO m=0,lag
      k=ind_start(istep)+2*(npassing+1)*m
      DO l=0,leg
        scalprod_pleg(m,l,istep,ispec)=                                 &
            SUM(pleg_bra(l,1:npassing+1,istep)*vec_in(k+1:k+npassing+1))
      ENDDO
      k=k+2*(npassing+1)
      DO l=0,leg,2
        scalprod_pleg(m,l,istep,ispec)=scalprod_pleg(m,l,istep,ispec)   &
           +SUM(pleg_bra(l,1:npassing+1,istep)*vec_in(k:k-npassing:-1))
      ENDDO
      DO l=1,leg,2
        scalprod_pleg(m,l,istep,ispec)=scalprod_pleg(m,l,istep,ispec)   &
           -SUM(pleg_bra(l,1:npassing+1,istep)*vec_in(k:k-npassing:-1))
      ENDDO
    ENDDO
!
  ENDDO
!
  !! Modification by Andreas F. Martitsch (20.08.2015)
  ! MPI Barrier -> collect scalprod (4D - leg,lag,phi,species)
  ! (mpro%allgather supports 3D and 4D matrices)
  !PRINT *,'mpro%getrank() before:', mpro%getrank()
  CALL mpro%allgather_inplace(scalprod_pleg)
  !PRINT *,'mpro%getrank() after:', mpro%getrank()
  !PRINT *,'scalprod_pleg, species = ',ispec
  !! End Modification by Andreas F. Martitsch (20.08.2015)
!
  DO istep=ibeg,iend
!
    npassing=npl(istep)
!
    !
    ! ailmm is now 5D object of species (alpha,alphap)
    !
    DO l=0,leg
      scalprod_pleg_tmp(0:lag,l,istep,ispec)=0.0d0
      DO ispecp=0,num_spec-1
        scalprod_pleg_tmp(0:lag,l,istep,ispec)=scalprod_pleg_tmp(0:lag,l,istep,ispec)+&
             MATMUL(ailmm_aa(0:lag,0:lag,l,ispec,ispecp),&
             scalprod_pleg(0:lag,l,istep,ispecp))
      ENDDO
      scalprod_pleg(0:lag,l,istep,ispec)=scalprod_pleg_tmp(0:lag,l,istep,ispec)
    ENDDO
    ! old behavior (for a single species)
    !DO l=0,leg
    !  scalprod_pleg(0:lag,l)=MATMUL(ailmm(0:lag,0:lag,l),scalprod_pleg(0:lag,l))
    !ENDDO
    !
    ! end of interaction with rest processors
    !

!
    DO m=0,lag
      k=ind_start(istep)+2*(npassing+1)*m
!
      IF(istep.GT.ibeg) THEN
        DO l=0,leg
          vec_out(k+1:k+npassing+1)=vec_out(k+1:k+npassing+1)            &
                      +scalprod_pleg(m,l,istep,ispec)*pleg_ket(l,1:npassing+1,istep)
        ENDDO
      ENDIF
!
      k=k+2*(npassing+1)
!
      IF(istep.LT.iend) THEN
        DO l=0,leg,2
          vec_out(k:k-npassing:-1)=vec_out(k:k-npassing:-1)              &
                      +scalprod_pleg(m,l,istep,ispec)*pleg_ket(l,1:npassing+1,istep)
        ENDDO
        DO l=1,leg,2
          vec_out(k:k-npassing:-1)=vec_out(k:k-npassing:-1)              &
                      -scalprod_pleg(m,l,istep,ispec)*pleg_ket(l,1:npassing+1,istep)
        ENDDO
      ENDIF
!
    ENDDO
!
  ENDDO
!
END SUBROUTINE integral_part
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine rearrange_phideps_old(ibeg,iend,npart,subsqmin,phi_divide,        &
                             phi_mfl,bhat_mfl,geodcu_mfl,h_phi_mfl,eta,  &
                             delt_pos,delt_neg,                          &
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
  implicit none
!
!  logical, parameter :: stepmode=.true.
  logical, parameter :: stepmode=.false.
  INTEGER, PARAMETER :: npoi=6, nder=0, npoihalf=npoi/2, nstepmin=8
  double precision, parameter :: bparabmax=0.2d0
!
  integer :: i,ibeg,iend,npart,istep,ibmin,npassing,npassing_prev
  integer :: ncross_l,ncross_r,ib,ie,intb,inte,k,imid,isplit
!
  double precision :: subsqmin,ht,ht2,bparab,x1,x2,f1,f2
!
  integer, dimension(1)              :: idummy
  integer, dimension(1:iend)         :: phi_divide
  integer, dimension(:), allocatable :: icross_l,icross_r
!
  DOUBLE PRECISION, DIMENSION(npoi)           :: tp,weight
  DOUBLE PRECISION, DIMENSION(0:nder,npoi)    :: coeff
  double precision, dimension(0:npart)        :: eta
  double precision, dimension(ibeg:iend)      :: phi_mfl,bhat_mfl
  double precision, dimension(ibeg:iend)      :: geodcu_mfl,h_phi_mfl
  double precision, dimension(ibeg:iend)      :: delt_pos,delt_neg
  double precision, dimension(ibeg:iend)      :: fact_pos_b,fact_neg_b
  double precision, dimension(ibeg:iend)      :: fact_pos_e,fact_neg_e
  double precision, dimension(:), allocatable :: phi_new,bhat_new
  double precision, dimension(:), allocatable :: geodcu_new,h_phi_new

  npassing = -1

  call fix_phiplacement_problem_old(ibeg,iend,npart,subsqmin,        &
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
  idummy=minloc(bhat_mfl(ibeg:iend))
  ibmin=idummy(1)+ibeg-1
!
  ncross_l=0
  if(ibmin.gt.ibeg) then
    istep=ibmin
    do i=0,npart
      if(1.d0-bhat_mfl(istep)*eta(i).GT.subsqmin) THEN
        npassing=i
      else
        exit
      endif
    enddo
    npassing_prev=npassing
    do istep=ibmin-1,ibeg,-1
      do i=0,npart
        if(1.d0-bhat_mfl(istep)*eta(i).GT.subsqmin) THEN
          npassing=i
        else
          exit
        endif
      enddo
      if(npassing.lt.npassing_prev) then
        ncross_l=ncross_l+1
        npassing_prev=npassing
      endif
    enddo
    if(ncross_l.gt.0) then
      allocate(icross_l(ncross_l))
      ncross_l=0
      istep=ibmin
      do i=0,npart
        if(1.d0-bhat_mfl(istep)*eta(i).GT.subsqmin) THEN
          npassing=i
        else
          exit
        endif
      enddo
      npassing_prev=npassing
      do istep=ibmin-1,ibeg,-1
        do i=0,npart
          if(1.d0-bhat_mfl(istep)*eta(i).GT.subsqmin) THEN
            npassing=i
          else
            exit
          endif
        enddo
        if(npassing.lt.npassing_prev) then
          ncross_l=ncross_l+1
          icross_l(ncross_l)=istep
          npassing_prev=npassing
        endif
      enddo
    endif
  endif
!
  ncross_r=0
  if(ibmin.lt.iend) then
    istep=ibmin
    do i=0,npart
      if(1.d0-bhat_mfl(istep)*eta(i).GT.subsqmin) THEN
        npassing=i
      else
        exit
      endif
    enddo
    npassing_prev=npassing
    do istep=ibmin+1,iend
      do i=0,npart
        if(1.d0-bhat_mfl(istep)*eta(i).GT.subsqmin) THEN
          npassing=i
        else
          exit
        endif
      enddo
      if(npassing.lt.npassing_prev) then
        ncross_r=ncross_r+1
        npassing_prev=npassing
      endif
    enddo
    if(ncross_r.gt.0) then
      allocate(icross_r(ncross_r))
      ncross_r=0
      istep=ibmin
      do i=0,npart
        if(1.d0-bhat_mfl(istep)*eta(i).GT.subsqmin) THEN
          npassing=i
        else
          exit
        endif
      enddo
      npassing_prev=npassing
      do istep=ibmin+1,iend
        do i=0,npart
          if(1.d0-bhat_mfl(istep)*eta(i).GT.subsqmin) THEN
            npassing=i
          else
            exit
          endif
        enddo
        if(npassing.lt.npassing_prev) then
          ncross_r=ncross_r+1
          icross_r(ncross_r)=istep
          npassing_prev=npassing
        endif
      enddo
    endif
  endif
!
! place ibmin to an odd point:
!
  if(mod(ibmin-ibeg,2).eq.1) then
    if(ncross_l.gt.0.and.ncross_r.gt.0) then
      if(icross_r(1)-ibmin.gt.ibmin-icross_l(1)) then
        ibmin=ibmin+1
      else
        ibmin=ibmin-1
      endif
    elseif(ncross_l.gt.0) then
      ibmin=ibmin+1
    elseif(ncross_r.gt.0) then
      ibmin=ibmin-1
    endif
  endif
!
! check the number of steps in sub-intervals for parabolic bhat term:
!
  if(ncross_l.gt.0) then
    ie=icross_l(1)
    do i=2,ncross_l
      ib=icross_l(i)
      if(ie-ib.lt.nstepmin) then
        imid=(ib+ie)/2
        x1=phi_mfl(imid)-phi_mfl(ib)
        x2=phi_mfl(ie)-phi_mfl(ib)
        f1=bhat_mfl(imid)-bhat_mfl(ib)
        f2=bhat_mfl(ie)-bhat_mfl(ib)
        bparab=abs((f1*x2-f2*x1)*x2/((x1-x2)*x1*f2))
        if(bparab.gt.bparabmax) then
          isplit=2*max(nint(0.5*float(nstepmin)/float(ie-ib)),1)
          phi_divide(ib+1:ie)=isplit
        endif
      endif
      ie=ib
    enddo
    ib=ibeg
    if(ie-ib.lt.nstepmin) then
      isplit=2*max(nint(0.5*float(nstepmin)/float(ie-ib)),1)
      phi_divide(ib+1:ie)=isplit
    endif
  endif
!
  if(ncross_r.gt.0) then
    ib=icross_r(1)
    do i=2,ncross_r
      ie=icross_r(i)
      if(ie-ib.lt.nstepmin) then
        imid=(ib+ie)/2
        x1=phi_mfl(imid)-phi_mfl(ib)
        x2=phi_mfl(ie)-phi_mfl(ib)
        f1=bhat_mfl(imid)-bhat_mfl(ib)
        f2=bhat_mfl(ie)-bhat_mfl(ib)
        bparab=abs((f1*x2-f2*x1)*x2/((x1-x2)*x1*f2))
        if(bparab.gt.bparabmax) then
          isplit=2*max(nint(0.5*float(nstepmin)/float(ie-ib)),1)
          phi_divide(ib+1:ie)=isplit
        endif
      endif
      ib=ie
    enddo
    ie=iend
    if(ie-ib.lt.nstepmin) then
      isplit=2*max(nint(0.5*float(nstepmin)/float(ie-ib)),1)
      phi_divide(ib+1:ie)=isplit
    endif
  endif
!
  if(maxval(phi_divide).gt.1) return
!
! change the integration variable phi -> sqrt(phi-phi0):
!
  if(stepmode) then
!
    allocate(phi_new(ibeg:iend),bhat_new(ibeg:iend))
    allocate(geodcu_new(ibeg:iend),h_phi_new(ibeg:iend))
    phi_new=phi_mfl
    bhat_new=bhat_mfl
    geodcu_new=geodcu_mfl
    h_phi_new=h_phi_mfl
!
    ie=ibmin
    do i=1,ncross_l
      ib=icross_l(i)
      ht=sqrt(phi_mfl(ie)-phi_mfl(ib))/float(ie-ib)
      ht2=ht**2
      k=ib
      do istep=ib+1,ie-1
        phi_new(istep)=phi_mfl(ib)+ht2*float(istep-ib)**2
        fact_pos_e(istep)=2.d0*ht*float(istep-ib)
        delt_pos(istep)=ht
        do while(phi_mfl(k).lt.phi_new(istep))
          k=k+1
        enddo
        intb=max(ibeg,min(iend-npoi+1,k-npoihalf))
        inte=intb+npoi-1
!
        CALL plagrange_coeff(npoi,nder,phi_new(istep),phi_mfl(intb:inte),coeff)
!
        bhat_new(istep)=sum(coeff(0,:)*bhat_mfl(intb:inte))
        geodcu_new(istep)=sum(coeff(0,:)*geodcu_mfl(intb:inte))
        h_phi_new(istep)=sum(coeff(0,:)*h_phi_mfl(intb:inte))
      enddo
      delt_pos(ie)=ht
      fact_pos_e(ie)=2.d0*ht*float(ie-ib)
      fact_pos_b(ib)=0.d0
      fact_pos_b(ib+1:ie-1)=fact_pos_e(ib+1:ie-1)
      ie=ib
    enddo
!
    ib=ibmin
    do i=1,ncross_r
      ie=icross_r(i)
      ht=sqrt(phi_mfl(ie)-phi_mfl(ib))/float(ie-ib)
      ht2=ht**2
      k=ib
      do istep=ib+1,ie-1
        phi_new(istep)=phi_mfl(ie)-ht2*float(ie-istep)**2
        delt_pos(istep)=ht
        fact_pos_b(istep)=2.d0*ht*float(ie-istep)
        do while(phi_mfl(k).lt.phi_new(istep))
          k=k+1
        enddo
        intb=max(ibeg,min(iend-npoi+1,k-npoihalf))
        inte=intb+npoi-1
!
        CALL plagrange_coeff(npoi,nder,phi_new(istep),phi_mfl(intb:inte),coeff)
!
        bhat_new(istep)=sum(coeff(0,:)*bhat_mfl(intb:inte))
        geodcu_new(istep)=sum(coeff(0,:)*geodcu_mfl(intb:inte))
        h_phi_new(istep)=sum(coeff(0,:)*h_phi_mfl(intb:inte))
      enddo
      delt_pos(ie)=ht
      fact_pos_b(ib)=2.d0*ht*float(ie-ib)
      fact_pos_e(ie)=0.d0
      fact_pos_e(ib+1:ie-1)=fact_pos_b(ib+1:ie-1)
      ib=ie
    enddo
!
    phi_mfl=phi_new
    bhat_mfl=bhat_new
    geodcu_mfl=geodcu_new
    h_phi_mfl=h_phi_new
!
    deallocate(phi_new,bhat_new,geodcu_new,h_phi_new)
!
  endif
!
  delt_neg(ibeg:iend-1)=delt_pos(ibeg+1:iend)
  fact_neg_b=fact_pos_e
  fact_neg_e=fact_pos_b
!
  if(allocated(icross_l)) deallocate(icross_l)
  if(allocated(icross_r)) deallocate(icross_r)
!
end subroutine rearrange_phideps_old
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine fix_phiplacement_problem_old(ibeg,iend,npart,subsqmin,        &
                                    phi_mfl,bhat_mfl,eta)
!
  use device_mod
!
  implicit none
!
  integer :: i,ibeg,iend,npart,istep,ibmin,npassing,npassing_prev
  integer :: ncross_l,ncross_r,ib,ie
!
  double precision :: subsqmin
!
  integer, dimension(1)              :: idummy
  integer, dimension(:), allocatable :: icross_l,icross_r
!
  double precision, dimension(0:npart)        :: eta
  double precision, dimension(ibeg:iend)      :: phi_mfl,bhat_mfl
  double precision, dimension(:), allocatable :: eta_cross_l,eta_cross_r
!
! determine level crossings:
!
  idummy=minloc(bhat_mfl(ibeg:iend))
  ibmin=idummy(1)+ibeg-1

  npassing = -1

  ncross_l=0
  if(ibmin.gt.ibeg) then
    istep=ibmin
    do i=0,npart
      if(1.d0-bhat_mfl(istep)*eta(i).GT.subsqmin) THEN
        npassing=i
      else
        exit
      endif
    enddo
    npassing_prev=npassing
    do istep=ibmin-1,ibeg,-1
      do i=0,npart
        if(1.d0-bhat_mfl(istep)*eta(i).GT.subsqmin) THEN
          npassing=i
        else
          exit
        endif
      enddo
      if(npassing.lt.npassing_prev) then
        ncross_l=ncross_l+1
        npassing_prev=npassing
      endif
    enddo
    if(ncross_l.gt.0) then
      allocate(icross_l(ncross_l),eta_cross_l(ncross_l))
      ncross_l=0
      istep=ibmin
      do i=0,npart
        if(1.d0-bhat_mfl(istep)*eta(i).GT.subsqmin) THEN
          npassing=i
        else
          exit
        endif
      enddo
      npassing_prev=npassing
      do istep=ibmin-1,ibeg,-1
        do i=0,npart
          if(1.d0-bhat_mfl(istep)*eta(i).GT.subsqmin) THEN
            npassing=i
          else
            exit
          endif
        enddo
        if(npassing.lt.npassing_prev) then
          ncross_l=ncross_l+1
          icross_l(ncross_l)=istep
          eta_cross_l(ncross_l)=eta(npassing_prev)
          npassing_prev=npassing
        endif
      enddo
      do i=1,ncross_l
        istep=icross_l(i)
        if(abs(bhat_mfl(istep-1)*eta_cross_l(i)-1.d0).lt. &
           abs(bhat_mfl(istep)  *eta_cross_l(i)-1.d0)) then
          open(111,file='phi_placement_problem.dat',position='append')
          write(111,*) ' propagator tag = ',fieldpropagator%tag, &
                       ' step number = ',istep-1,                &
                       ' 1 / bhat = ',1.d0/bhat_mfl(istep-1),    &
                       ' eta = ',eta_cross_l(i)
          close(111)
          bhat_mfl(istep-1)=1/eta_cross_l(i)
        elseif(abs(bhat_mfl(istep+1)*eta_cross_l(i)-1.d0).lt. &
               abs(bhat_mfl(istep)  *eta_cross_l(i)-1.d0)) then
          open(111,file='phi_placement_problem.dat',position='append')
          write(111,*) ' propagator tag = ',fieldpropagator%tag, &
                       ' step number = ',istep+1,                &
                       ' 1 / bhat = ',1.d0/bhat_mfl(istep+1),    &
                       ' eta = ',eta_cross_l(i)
          bhat_mfl(istep+1)=1/eta_cross_l(i)
          close(111)
        endif
      enddo
      deallocate(icross_l,eta_cross_l)
    endif
  endif
!
  ncross_r=0
  if(ibmin.lt.iend) then
    istep=ibmin
    do i=0,npart
      if(1.d0-bhat_mfl(istep)*eta(i).GT.subsqmin) THEN
        npassing=i
      else
        exit
      endif
    enddo
    npassing_prev=npassing
    do istep=ibmin+1,iend
      do i=0,npart
        if(1.d0-bhat_mfl(istep)*eta(i).GT.subsqmin) THEN
          npassing=i
        else
          exit
        endif
      enddo
      if(npassing.lt.npassing_prev) then
        ncross_r=ncross_r+1
        npassing_prev=npassing
      endif
    enddo
    if(ncross_r.gt.0) then
      allocate(icross_r(ncross_r),eta_cross_r(ncross_r))
      ncross_r=0
      istep=ibmin
      do i=0,npart
        if(1.d0-bhat_mfl(istep)*eta(i).GT.subsqmin) THEN
          npassing=i
        else
          exit
        endif
      enddo
      npassing_prev=npassing
      do istep=ibmin+1,iend
        do i=0,npart
          if(1.d0-bhat_mfl(istep)*eta(i).GT.subsqmin) THEN
            npassing=i
          else
            exit
          endif
        enddo
        if(npassing.lt.npassing_prev) then
          ncross_r=ncross_r+1
          icross_r(ncross_r)=istep
          eta_cross_r(ncross_r)=eta(npassing_prev)
          npassing_prev=npassing
        endif
      enddo
      do i=1,ncross_r
        istep=icross_r(i)
        if(abs(bhat_mfl(istep-1)*eta_cross_r(i)-1.d0).lt. &
           abs(bhat_mfl(istep)  *eta_cross_r(i)-1.d0)) then
          open(111,file='phi_placement_problem.dat',position='append')
          write(111,*) ' propagator tag = ',fieldpropagator%tag, &
                       ' step number = ',istep-1,                &
                       ' 1 / bhat = ',1.d0/bhat_mfl(istep-1),    &
                       ' eta = ',eta_cross_r(i)
          close(111)
          bhat_mfl(istep-1)=1/eta_cross_r(i)
        elseif(abs(bhat_mfl(istep+1)*eta_cross_r(i)-1.d0).lt. &
               abs(bhat_mfl(istep)  *eta_cross_r(i)-1.d0)) then
          open(111,file='phi_placement_problem.dat',position='append')
          write(111,*) ' propagator tag = ',fieldpropagator%tag, &
                       ' step number = ',istep+1,                &
                       ' 1 / bhat = ',1.d0/bhat_mfl(istep+1),    &
                       ' eta = ',eta_cross_r(i)
          close(111)
          bhat_mfl(istep+1)=1/eta_cross_r(i)
        endif
      enddo
      deallocate(icross_r,eta_cross_r)
    endif
  endif
!
end subroutine fix_phiplacement_problem_old
