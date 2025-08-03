!Sergei 20.07.2006 : modification of boundary layer is done now locally,
!                    filtering of the magnetic field maxima has been removed,
!                    old unused lines which were commented have been removed

subroutine ripple_solver_ArnoldiO2( &
    npass_l, npass_r, nvelocity, &
    amat_plus_plus, amat_minus_minus, &
    amat_plus_minus, amat_minus_plus, &
    source_p, source_m, &
    flux_p, flux_m, &
    qflux, &
    ierr &
    )

    use device_mod
    use flint_mod, only: phi_divide                      !<-in Winny
    use collisionality_mod, only: collpar, conl_over_mfp, isw_lorentz, &
                                  isw_energy, isw_integral, isw_axisymm, & !<-in Winny
                                  isw_momentum, nvel, num_spec, lsw_multispecies
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
    use lapack_band
    use rkstep_mod
    use polleg_mod
    use propagator_mod, only: prop_ripple_plot, prop_reconstruct, flux_mr, flux_pl, &
                              eta_modboundary_l, eta_modboundary_r, &
                              prop_reconstruct_levels
    use sparse_mod, only: sparse_solve_method, sparse_solve, &
                          column_full2pointer, sparse_solver_test, remap_rc
    use mag_interface_mod, only: surface_boozer_B00, travis_convfac, boozer_s, mag_magfield
    use ntv_eqmat_mod, only: nz_symm, nz_asymm, nz_per_pos, nz_per_neg, &
                             irow_symm, icol_symm, amat_symm, &
                             irow_per_pos, icol_per_pos, &
                             irow_per_neg, icol_per_neg, &
                             irow_asymm, icol_asymm, amat_asymm, &
                             f0_coll, f0_ttmp, f0_coll_all, f0_ttmp_all, &
                             nz_regper, irow_regper, icol_regper, amat_regper
    use partpa_mod, only: bmod0
    use development
  !! Modifications by Andreas F. Martitsch (12.03.2014)
    ! boozer_curr_tor, boozer_curr_pol, boozer_psi_pr,
    ! boozer_sqrtg11 and boozer_isqrg are now converted
    ! to cgs-units within neo_magfie.
    ! This step requires changes within rhs_kin.f90 and
    ! ripple_solver.f90!
    use neo_magfie, only: boozer_iota, boozer_curr_pol_hat, &
                          boozer_curr_tor_hat, boozer_psi_pr_hat, boozer_curr_pol_hat_s, &
                          boozer_curr_tor_hat_s, boozer_iota_s, boozer_isqrg
  !! End Modifications by Andreas F. Martitsch (12.03.2014)
  !! Modifications by Andreas F. Martitsch (12.06.2014)
    ! quantities of the perturbation field extracted
    ! from the Boozer file (e.g., toroidal mode number m_phi)
    use neo_magfie_perturbation, only: calc_bnoverb0_arr, calc_ntv_output
  !! End Modifications by Andreas F. Martitsch (12.06.2014)
  !! Modification by Andreas F. Martitsch (14.07.2015)
    ! Extra input for NTV computations
    use ntv_mod, only: isw_qflux_NA, MtOvR, B_rho_L_loc, &
                       m_phi, qflux_symm, eps_M_2_val, av_gphph_val, av_inv_bhat_val, &
                       qflux_symm_allspec, qflux_ntv_allspec, &
                       MtOvR_spec, isw_calc_Er, B_rho_L_loc_spec, isw_calc_MagDrift
    use ntv_mod, only: get_Er, get_B_rho_L_loc
  !! End Modification by Andreas F. Martitsch (14.07.2015)
    ! MPI SUPPORT for multi-species part
    ! (run with, e.g.,  mpiexec -np 3 ./neo2.x)
    use mpiprovider_module, only: mpro
    ! Load x1mm and x2mm (=energy dependence of drift frequencies)
    ! from collision operator module. This step allows for
    ! support of different basis functions and replaces routine "lagxmm".
    use collop
    use arnoldi_mod, only: iterator, f_init_arnoldi, &
      & lsw_write_flux_surface_distribution, write_flux_surface_distribution

    implicit none
    complex(dp), parameter :: imun = (0.d0, 1.d0)
    real(dp), parameter :: PI = 3.141592653589793238462643383279502884197_dp

    ! parameter list
    integer, intent(out)   :: npass_l
    integer, intent(out)   :: npass_r
    integer, intent(out)   :: nvelocity
    real(dp), dimension(:, :), allocatable, intent(inout) :: amat_plus_plus
    real(dp), dimension(:, :), allocatable, intent(inout) :: amat_minus_minus
    real(dp), dimension(:, :), allocatable, intent(inout) :: amat_plus_minus
    real(dp), dimension(:, :), allocatable, intent(inout) :: amat_minus_plus
    real(dp), dimension(:, :), allocatable, intent(inout) :: source_p
    real(dp), dimension(:, :), allocatable, intent(inout) :: source_m
    real(dp), dimension(:, :), allocatable, intent(inout) :: flux_p
    real(dp), dimension(:, :), allocatable, intent(inout) :: flux_m
    real(dp), dimension(:, :), allocatable, intent(out)   :: qflux
    integer, intent(out)   :: ierr

    ! local stuff

    integer :: add_global_eta
    integer :: ub_eta, ub_mag, ub_eta_loc
    integer :: ub_eta_prev, ub_eta_next, npass_l_out, npass_r_out
    integer :: npart
    integer :: ibeg, iend
    integer :: bin_split_mode
    integer :: i, i_loc
    integer, allocatable :: eta_loc_ind(:), eta_glob_ind(:)
    real(dp) :: phibeg, phiend
    real(dp) :: xetami, xetama
    real(dp) :: rt0
    real(dp) :: b_prop_l, b_prop_r, b_prop_min
    real(dp) :: b_max_l, b_max_r, b_min, width
    real(dp) :: b_l, b_r, eta_l, eta_r
    real(dp), allocatable :: eta(:), eta_prev(:), eta_next(:)
    real(dp), allocatable :: eta_loc(:), eta_glob(:) !<- LOCAL

    ! Timing
    real(dp) :: time_start, time_factorization, time_solver

    !------------------------------------------------------------------------
    ! SERGEI
    !------------------------------------------------------------------------
    ! additional definitions from Sergei

    integer :: npart_loc, info
    integer :: ndim, istep, npassing, ioddeven
    integer :: ipart, ipart1
    integer :: i1
    integer :: k, i1min
    integer :: kmax
    integer, dimension(:), allocatable :: ipivot

    real(dp)                         :: aiota
    real(dp)                         :: eta0
    real(dp)                         :: subsq, subsqmin
    real(dp)                         :: diflam, diflampow, coefdir
    real(dp)                         :: coefenu, coefenu_averb   !!!term[1]
    real(dp) :: amin2ovb
    complex(dp) :: coef_cf

    real(dp), dimension(6)           :: alp, bet, gam, del
    real(dp), dimension(:, :), allocatable :: amat, bvec_lapack, deriv_coef
    complex(dp), dimension(:, :), allocatable :: amat_z, bvec_lapack_z
    real(dp), dimension(:, :), allocatable :: fun_coef
    real(dp), dimension(:, :), allocatable :: alambd
    complex(dp), dimension(:, :), allocatable :: Vg_vp_over_B
    real(dp), dimension(:), allocatable :: delta_eta
    real(dp), dimension(:, :), allocatable :: enu_coef        !!!term[1]
    real(dp), dimension(:, :), allocatable :: enu_coef2       !!!NTV
    integer :: km1, kp1, m, mfactorial, nplp1
    real(dp), dimension(:, :), allocatable ::  alampow !!!term[3]
    real(dp), dimension(:, :, :), allocatable ::  vrecurr !!!term[3]
    real(dp), dimension(:, :), allocatable ::  dellampow !!!term[3]
    real(dp), dimension(:, :), allocatable ::  dellampow2 !!!NTV
    real(dp), dimension(:, :), allocatable ::  convol_polpow !!!term[3]
    real(dp), dimension(:, :), allocatable ::  coefleg      !!!terms[2,3]

    integer :: ntotsize, nts_r, nts_l, kk

    real(dp), dimension(:, :, :, :), allocatable :: derivs_plot, fun_write
    integer :: iplot, nphiplot, iunit_phi, iunit_sizes
    integer :: iunit_dt_p, iunit_dt_m
    integer :: iunit_sp_p, iunit_sp_m
    integer :: iunit_et_p, iunit_et_m
    real(dp) :: phiplot, delphiplot, facnorm_p, facnorm_m
    real(dp) :: boundlayer_ignore
    integer :: ignore_lb, ignore_rb, ignore_lb_out, ignore_rb_out, modify_bl, modify_br
    real(dp) :: bhat_changed_l, bhat_changed_r
    real(dp) :: bhat_changed_l_out, bhat_changed_r_out
    real(dp) :: sign_of_bphi
    integer :: icounter

    CHARACTER(len=100) :: propname
    integer :: n_2d_size, nz_sq, nz_beg, npassing_prev, k_prev, mm
    integer :: iter, nphiequi, npassing_next, n_arnoldi, mode_iter
    integer :: isw_regper, nz_coll, nz_ttmp, nz_coll_beg
    integer :: nrow, ncol, nz, iopt
    real(dp) :: delphim1, deloneovb, step_factor_p, step_factor_m
    real(dp) :: deleta_factor
    integer, dimension(:), allocatable :: ind_start
    integer, dimension(:), allocatable :: irow, icol, ipcol
    integer, dimension(:), allocatable :: irow_coll, icol_coll
    integer, dimension(:), allocatable :: irow_ttmp, icol_ttmp
    real(dp), dimension(:), allocatable :: amat_coll, amat_ttmp
    complex(dp), dimension(:), allocatable :: amat_sp, bvec_sp
    complex(dp), dimension(:), allocatable :: bvec_iter, bvec_prev
    complex(dp), dimension(:), allocatable :: bvec_parflow
    ! Use pre-conditioned iterations:
    ! -> remove null-space of axisymmetric solution (energy conservation)
    complex(dp) :: denom_energ, coef_energ, denom_dens, coef_dens
    complex(dp), dimension(:), allocatable :: energvec_bra, energvec_ket
    complex(dp), dimension(:), allocatable :: densvec_bra, densvec_ket
    ! End Use pre-conditioned iterations
    complex(dp), dimension(:, :), allocatable :: flux_vector, source_vector
    complex(dp), dimension(:, :), allocatable :: flux_vector_plot
    complex(dp), dimension(:, :), allocatable :: basevec_p
    integer :: isw_lor, isw_ene, isw_intp
    integer, dimension(:), allocatable :: npl
    real(dp), dimension(:, :, :), allocatable :: rhs_mat_fzero
    real(dp), dimension(:, :, :), allocatable :: rhs_mat_lorentz
    real(dp), dimension(:, :, :), allocatable :: ttmp_mat
    complex(dp), dimension(:, :, :), allocatable :: q_rip
    complex(dp), dimension(:, :), allocatable :: q_rip_1
    complex(dp), dimension(:, :), allocatable :: q_rip_incompress
    complex(dp), dimension(:, :), allocatable :: q_rip_parflow
    real(dp), dimension(:, :, :), allocatable :: rhs_mat_energ
    real(dp), dimension(:, :, :), allocatable :: rhs_mat_energ2     !NTV
    real(dp), dimension(:, :, :), allocatable :: pleg_bra, pleg_ket
    complex(dp), dimension(:, :), allocatable :: convol_flux, convol_curr
    complex(dp), dimension(:, :), allocatable :: convol_flux_0
    real(dp), dimension(:, :), allocatable :: scalprod_pleg
    complex(dp), dimension(:), allocatable :: scalprod
    real(dp), dimension(:), allocatable :: phi_mfl
    real(dp), dimension(:), allocatable :: bhat_mfl, h_phi_mfl
    complex(dp), dimension(:), allocatable :: geodcu_mfl
    complex(dp), dimension(:), allocatable :: geodcu_forw, geodcu_back
    real(dp), dimension(:), allocatable :: dlogbdphi_mfl
    real(dp), dimension(:), allocatable :: delt_pos, delt_neg
    real(dp), dimension(:), allocatable :: fact_pos_b, fact_neg_b
    real(dp), dimension(:), allocatable :: fact_pos_e, fact_neg_e
    integer          :: nreal, ncomp
    complex(dp) :: expforw, expbackw, perbou_pos, perbou_neg, rotfactor
    real(dp) :: Er, avEparB_ov_avb2 ! radial and inductive electric field
    real(dp) :: a1b, a2b, hatOmegaE, hatOmegaB, denomjac
  !! Modifications by Andreas F. Martitsch (14.03.2014)
    ! Subsequent quantities are given now in cgs-units and they are
    ! renormalized using bmod0 within neo_magfie:
    !real(dp) :: bcovar_theta,bcovar_phi,dbcovar_theta_ds,dbcovar_phi_ds
    ! For this reason these variables are renamed:
    real(dp) :: bcovar_theta_hat, bcovar_phi_hat
    real(dp) :: dbcovar_theta_hat_ds, dbcovar_phi_hat_ds
  !! End Modifications by Andreas F. Martitsch (14.03.2014)
    real(dp) :: scalefac_kG
    real(dp), dimension(:, :), allocatable :: arr_real
    complex(dp), dimension(:, :), allocatable :: arr_comp
  !! Modifications by Andreas F. Martitsch (13.06.2014)
    ! Subsequent quantities (given now in cgs-units) are computed by
    ! magdata_for_particles and stored within the fieldpropagator-structure.
    ! This step required changes within neo_magfie, magfie, mag,
    ! magdata_for_particles, mag_interface_mod, plagrange_mod,
    ! modify_propagator and magnetics_mod.
    real(dp), dimension(:), allocatable :: dlogbds_mfl
    real(dp), dimension(:), allocatable :: bcovar_s_hat_mfl
    real(dp), dimension(:), allocatable :: dbcovar_s_hat_dphi_mfl
    complex(dp), dimension(:), allocatable :: bnoverb0, dbnoverb0_dphi_mfl
    ! For testing you can specify here an artificial perturbation field
    complex(dp), dimension(:), allocatable :: bnoverb0_test, dbnoverb0_dphi_mfl_test
    ! amplitude
    complex(dp) :: bnoverb0_test_val = (1.0d-3, 0.0d-0)
    ! poloidal mode number
    integer :: m_theta = 0
  !! End Modifications by Andreas F. Martitsch (13.06.2014)
  !! Modification by Andreas F. Martitsch (28.07.2015)
    !  multi-species part
    integer :: ispec, ispecp, ispecpp ! species indices
    integer :: drive_spec
    complex(dp), dimension(:, :, :), allocatable :: source_vector_all
    real(dp), dimension(:, :, :, :), allocatable :: qflux_allspec
    logical :: problem_type
    real(dp), dimension(:, :), allocatable :: source_vector_real
  !! End Modification by Andreas F. Martitsch (28.07.2015)
    complex(dp), dimension(:), allocatable :: ttmpfact
    logical :: colltest = .FALSE.
    logical :: ttmptest = .FALSE.
    logical :: nobounceaver = .TRUE.

! DEBUGGING
    integer :: i_ctr = 0, uw, uw_new
    logical :: lsw_debug_distfun = .FALSE.
    logical :: addboucol = .false.

    ! multi-species part - MPI rank determines species
    ispec = mpro%getRank()
    print *, "Species: ", ispec
    call collop_set_species(ispec)

    niter = 100       !maximum number of integral part iterations
    n_arnoldi = 500     !maximum number of Arnoldi iterations
    isw_regper = 1       !regulariization by periodic boundary condition
    epserr_iter = 1.d-7  !relative error of integral part iterations
    sparse_solve_method = 3 !2,3 - with and without iterative refinement, resp.

    addboucol = .true.
    !------------------------------------------------------------------------
    ! END SERGEI
    !------------------------------------------------------------------------

    !------------------------------------------------------------------------
    ! Winny: sparse
    !------------------------------------------------------------------------
    ! ifsparse_solve should talk
    !  sparse_talk = .TRUE. ! default .FALSE. - neo2.in - settings
    !
    ! sparse method - only 1 (SuperLU) implemented
    !  sparse_solve_method = 1 ! default 0 - neo2.in - settings
    ! ifsparse_solve is called with (sparse_solve_method .eq. 0)
    !  program stops
    !  one should call a normal solver in this case
    !
    ! These are three possibilities to call sparse_solve
    !  call sparse_solve(nrow,ncol,nz,irow,icol,val,b)  ! full column index
    !  call sparse_solve(nrow,ncol,nz,irow,pcol,val,b)  ! column pointer
    !  call sparse_solve(A,b)                           ! full matrix
    ! results are returned in b
    !
    ! Input
    !  integer :: nrow,ncol,nz,nrhs
    !  integer, dimension(:), allocatable :: irow,icol,pcol
    !  real(dp), dimension(:), allocatable :: val
    !  real(dp), dimension(:,:), allocatable :: A
    !
    ! In/output
    !  real(dp), dimension(:), allocatable :: b
    ! or
    !  real(dp), dimension(:,:), allocatable :: b
    !
    ! Matrix - square
    !  ncol = nrow
    !  allocate( irow(nz) )     ! row index
    !  allocate( icol(nz) )     ! column index (column ordering)
    !  allocate( pcol(ncol+1) ) ! column pointer (other possibility)
    !  allocate( val(nz) )      ! values
    !  allocate( A(nrow,ncol) ) ! full matrix will be converted to sparse
    !
    ! rhs
    !  allocate( b(nrow) )
    ! or
    !  allocate( b(nrow,nrhs) )
    !
    !------------------------------------------------------------------------

    !------------------------------------------------------------------------
    ! reconstruction work by Winny
    !------------------------------------------------------------------------
    ! values of reconstructed fluxes are put into fun_mr and fun_pl
    ! allocation in the 2nd-dimension is from 0 on
    ! (as requested by Sergei)
    if (prop_reconstruct .EQ. 2) then
        prop_ripple_plot = 1
    end if
    !------------------------------------------------------------------------
    ! end reconstruction work by Winny
    !------------------------------------------------------------------------

    if (isw_momentum .EQ. 1) then ! Grid
        print *, 'ERROR: isw_momentum = ', isw_momentum, ' not implemented in ripple solver!'
        print *, 'I stop here'
        stop
    end if

    ierr = 0
    boundlayer_ignore = 0.01d0

    iunit_phi = 200
    iunit_sizes = 201
    iunit_dt_p = 300
    iunit_dt_m = 301
    iunit_sp_p = 302
    iunit_sp_m = 303
    iunit_et_p = 304
    iunit_et_m = 305

    if (isw_lorentz .EQ. 1) then
        lag = 0
        nvel = 0
        legmax = 2            !Caution: legmax cannot be less than 2 (for convol_flux)
        isw_lor = 1
        isw_ene = 0
        isw_intp = 0
        anumm(0, 0) = 1.d0
        asource(0, 1:3) = 1.d0
        weightlag(1:3, 0) = 1.d0
    else
        legmax = MAX(leg, 2)!Caution: legmax cannot be less than 2 (for convol_flux)
        isw_lor = 1
        isw_ene = isw_energy
        isw_intp = isw_integral
    end if

    nvelocity = nvel

    iprintflag = 1
    !------------------------------------------------------------------------
    ! here i get everything for you what you might need

    ! eta related stuff
    ub_eta = UBOUND(fieldpropagator%ch_act%eta, 1)
    npart = ub_eta
    add_global_eta = 0
    if (prop_reconstruct_levels .EQ. 0) then
        ! This is the old stuff
        ! here eta = eta_glob = eta_loc
        ub_eta_loc = ub_eta

        ! eta
        if (allocated(eta)) deallocate (eta)
        allocate (eta(0:ub_eta))
        eta = fieldpropagator%ch_act%eta
        ! eta_glob = eta
        if (allocated(eta_glob)) deallocate (eta_glob)
        allocate (eta_glob(0:ub_eta))
        eta_glob = fieldpropagator%ch_act%eta

        ! eta_loc = eta
        if (allocated(eta_loc)) deallocate (eta_loc)
        allocate (eta_loc(0:ub_eta_loc))
        eta_loc = fieldpropagator%ch_act%eta

        ! index
        if (allocated(eta_glob_ind)) deallocate (eta_glob_ind)
        allocate (eta_glob_ind(0:ub_eta))
        if (allocated(eta_loc_ind)) deallocate (eta_loc_ind)
        allocate (eta_loc_ind(0:ub_eta_loc))
        do i = 0, ub_eta
            eta_loc_ind(i) = i
            eta_glob_ind(i) = i
        end do

    else
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
        !    1 - ifadditional levels exist; 0 - otherwise
        !
        ! eta_loc_ind(0:ub_eta_loc) - integer:
        !    index where the value of eta_loc can be found in eta_glob
        !
        ! eta_glob_ind(0:ub_eta) - integer:
        !    index where the value of eta_glob can be found in eta_loc
        !    or -1 ifit does not exist in eta_loc
        !
        ! reconstruction has to be done,
        !    if(add_global_eta .gt. 0)
        !
        ! So the usage is now (input file):
        !
        !  prop_reconstruct_levels = 0 : no reconstruction
        !                            1 : reconstruction (ifproper bsfunc_local_solver value)
        !
        !  bsfunc_local_solver = 0 : original stuff (no reconstruction possible)
        !                        1 : only local (no reconstruction possible)
        !                        2 : local + "absolute maximum" (reconstruction possible)
        !                        3 : local + "absolute maximum" + rest (reconstruction possible)
        !                        4 : not recommended (no reconstruction possible)
        !     ifreconstruction is not possible, prop_reconstruct_levels is automatically set to 0
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

        ub_eta_loc = UBOUND(fieldpropagator%ch_act%eta_loc, 1)
        ! eta
        if (allocated(eta)) deallocate (eta)
        allocate (eta(0:ub_eta))
        eta = fieldpropagator%ch_act%eta

        ! eta_glob = eta
        if (allocated(eta_glob)) deallocate (eta_glob)
        allocate (eta_glob(0:ub_eta))
        eta_glob = fieldpropagator%ch_act%eta

        ! eta_loc
        if (allocated(eta_loc)) deallocate (eta_loc)
        allocate (eta_loc(0:ub_eta_loc))
        eta_loc = fieldpropagator%ch_act%eta_loc

        ! index
        if (allocated(eta_glob_ind)) deallocate (eta_glob_ind)
        allocate (eta_glob_ind(0:ub_eta))
        if (allocated(eta_loc_ind)) deallocate (eta_loc_ind)
        allocate (eta_loc_ind(0:ub_eta_loc))
        i_loc = 0
        do i = 0, ub_eta
            if (eta_glob(i) .EQ. eta_loc(i_loc)) then
                eta_loc_ind(i_loc) = i
                eta_glob_ind(i) = i_loc
                i_loc = i_loc + 1
            else
                eta_glob_ind = -1
                add_global_eta = 1
            end if
        end do

    end if
    xetami = eta(0)
    xetama = eta(ub_eta)

    ! previous eta
    if (associated(fieldpropagator%prev)) then
        fieldripple => fieldpropagator%prev%ch_act
    else
        fieldripple => fieldpropagator%parent%ch_las%ch_act
    end if
    ub_eta_prev = UBOUND(fieldripple%eta, 1)
    if (allocated(eta_prev)) deallocate (eta_prev)
    allocate (eta_prev(0:ub_eta_prev))
    eta_prev = fieldripple%eta
    ! and next eta
    if (associated(fieldpropagator%next)) then
        fieldripple => fieldpropagator%next%ch_act
    else
        fieldripple => fieldpropagator%parent%ch_fir%ch_act
    end if
    ub_eta_next = UBOUND(fieldripple%eta, 1)
    if (allocated(eta_next)) deallocate (eta_next)
    allocate (eta_next(0:ub_eta_next))
    eta_next = fieldripple%eta
    ! fieldripple back to original
    fieldripple => fieldpropagator%ch_act

    ! bhat (eta) on the left and right side of propagator
    b_l = fieldpropagator%b_l
    b_r = fieldpropagator%b_r
    eta_l = 1.0_dp/b_l
    eta_r = 1.0_dp/b_r

    ! additional stuff
    ! device
    rt0 = fieldpropagator%parent%parent%parent%parent%r0
    ! fieldpropagator
    b_prop_l = fieldpropagator%b_l
    b_prop_r = fieldpropagator%b_r
    b_prop_min = fieldpropagator%b_min

    ! fieldripple
    bin_split_mode = fieldpropagator%ch_act%bin_split_mode
    b_max_l = fieldpropagator%ch_act%b_max_l
    b_max_r = fieldpropagator%ch_act%b_max_r
    b_min = fieldpropagator%ch_act%b_min
    width = fieldpropagator%ch_act%width

    do i = 1, ub_eta
        if ((1.d0 - b_l*eta(i)) + 10.d0*EPSILON(1.d0) .GT. 0.d0) npass_l = i
        if ((1.d0 - b_r*eta(i)) + 10.d0*EPSILON(1.d0) .GT. 0.d0) npass_r = i
    end do
    do i = 1, ub_eta_prev
        if ((1.d0 - b_l*eta_prev(i)) + 10.d0*EPSILON(1.d0) .GT. 0.d0) npass_l_out = i
    end do
    do i = 1, ub_eta_next
        if ((1.d0 - b_r*eta_next(i)) + 10.d0*EPSILON(1.d0) .GT. 0.d0) npass_r_out = i
    end do

    ! Ignore the boundary layer ifit is too narrow
    ignore_lb = 0
    bhat_changed_l = 0.d0
    ignore_lb_out = 0
    bhat_changed_l_out = 0.d0
    modify_bl = 0
    ignore_rb = 0
    bhat_changed_r = 0.d0
    ignore_rb_out = 0
    bhat_changed_r_out = 0.d0
    modify_br = 0

    if (.false.) then

        ! Left boundary:

        ! check own eta-levels
        if (eta_l - eta(npass_l) .LT. &
            (eta(npass_l) - eta(npass_l - 1))*boundlayer_ignore) then
            ignore_lb = 1
            bhat_changed_l = 1.d0/eta(npass_l) + 100.d0*EPSILON(1.d0)
        else
            ignore_lb = 0
            bhat_changed_l = 0.d0
        end if

        ! check outer eta-levels
        if (eta_l - eta_prev(npass_l_out) .LT. &
            (eta_prev(npass_l_out) - eta_prev(npass_l_out - 1))*boundlayer_ignore) then
            ignore_lb_out = 1
            bhat_changed_l_out = 1.d0/eta_prev(npass_l_out) + 100.d0*EPSILON(1.d0)
        else
            ignore_lb_out = 0
            bhat_changed_l_out = 0.d0
        end if

        ! forbid bhat modification ifa regular band is eliminated in addition to b.l.
        if (1.d0 - bhat_changed_l*eta_prev(npass_l_out - 1) + 10.d0*EPSILON(1.d0) &
            .LE. 0.d0) then
            ignore_lb = 0
            bhat_changed_l = 0.d0
            print *, 'cannot ignore left boundary layer: jump over normal band'
        end if
        if (1.d0 - bhat_changed_l_out*eta(npass_l - 1) + 10.d0*EPSILON(1.d0) &
            .LE. 0.d0) then
            ignore_lb_out = 0
            bhat_changed_l_out = 0.d0
            print *, 'cannot ignore right boundary layer: jump over normal band'
        end if

        ! final value of modified bhat
        if (ignore_lb .EQ. 1 .OR. ignore_lb_out .EQ. 1) then
            bhat_changed_l = MAX(bhat_changed_l, bhat_changed_l_out)
            modify_bl = 1
            print *, 'field at the left boundary modified'
        else
            modify_bl = 0
        end if

        ! final decision on the boundary layer
        if (modify_bl .EQ. 1) then
            if (1.d0 - bhat_changed_l*eta(npass_l) + 10.d0*EPSILON(1.d0) .LE. 0.d0) then
                ignore_lb = 1
                print *, 'left boundary layer ignored'
            else
                ignore_lb = 0
            end if
        end if

        ! Right boundary:

        ! check own eta-levels
        if (eta_r - eta(npass_r) .LT. &
            (eta(npass_r) - eta(npass_r - 1))*boundlayer_ignore) then
            ignore_rb = 1
            bhat_changed_r = 1.d0/eta(npass_r) + 100.d0*EPSILON(1.d0)
        else
            ignore_rb = 0
            bhat_changed_r = 0.d0
        end if

        ! check outer eta-levels
        if (eta_r - eta_next(npass_r_out) .LT. &
            (eta_next(npass_r_out) - eta_next(npass_r_out - 1))*boundlayer_ignore) then
            ignore_rb_out = 1
            bhat_changed_r_out = 1.d0/eta_next(npass_r_out) + 100.d0*EPSILON(1.d0)
        else
            ignore_rb_out = 0
            bhat_changed_r_out = 0.d0
        end if

        ! forbid bhat modification ifa regular band is eliminated in addition to b.l.
        if (1.d0 - bhat_changed_r*eta_next(npass_r_out - 1) + 10.d0*EPSILON(1.d0) &
            .LE. 0.d0) then
            ignore_rb = 0
            bhat_changed_r = 0.d0
            print *, 'cannot ignore right boundary layer: jump over normal band'
        end if
        if (1.d0 - bhat_changed_r_out*eta(npass_r - 1) + 10.d0*EPSILON(1.d0) &
            .LE. 0.d0) then
            ignore_rb_out = 0
            bhat_changed_r_out = 0.d0
            print *, 'cannot ignore left boundary layer: jump over normal band'
        end if

        ! final value of modified bhat
        if (ignore_rb .EQ. 1 .OR. ignore_rb_out .EQ. 1) then
            bhat_changed_r = MAX(bhat_changed_r, bhat_changed_r_out)
            modify_br = 1
            print *, 'field at the right boundary modified'
        else
            modify_br = 0
        end if

        ! final decision on the boundary layer
        if (modify_br .EQ. 1) then
            if (1.d0 - bhat_changed_r*eta(npass_r) + 10.d0*EPSILON(1.d0) .LE. 0.d0) then
                ignore_rb = 1
                print *, 'right boundary layer ignored'
            else
                ignore_rb = 0
            end if
        end if

    end if

    ! place for boundary
    npass_l = npass_l + 1 - ignore_lb
    npass_r = npass_r + 1 - ignore_rb

    ! allocate and copy the magnetic stuff
    ub_mag = UBOUND(fieldpropagator%coords%x2, 1)
    ibeg = 0 - 2*modify_bl
    iend = ub_mag + 2*modify_br
    print *, ub_mag, ibeg, iend

    if (allocated(phi_divide)) deallocate (phi_divide)
    allocate (phi_divide(1:ub_mag))
    phi_divide = 2

    if (allocated(phi_mfl)) deallocate (phi_mfl)
    allocate (phi_mfl(ibeg:iend))
    phi_mfl(0:ub_mag) = fieldpropagator%coords%x2
    phibeg = phi_mfl(0)
    phiend = phi_mfl(ub_mag)

    if (allocated(bhat_mfl)) deallocate (bhat_mfl)
    allocate (bhat_mfl(ibeg:iend))
    bhat_mfl(0:ub_mag) = fieldpropagator%mdata%bhat

    if (allocated(geodcu_mfl)) deallocate (geodcu_mfl)
    allocate (geodcu_mfl(ibeg:iend))
    geodcu_mfl(0:ub_mag) = fieldpropagator%mdata%geodcu

    if (allocated(geodcu_forw)) deallocate (geodcu_forw)
    allocate (geodcu_forw(ibeg:iend))
    if (allocated(geodcu_back)) deallocate (geodcu_back)
    allocate (geodcu_back(ibeg:iend))

    if (allocated(h_phi_mfl)) deallocate (h_phi_mfl)
    allocate (h_phi_mfl(ibeg:iend))
    h_phi_mfl(0:ub_mag) = fieldpropagator%mdata%h_phi
    sign_of_bphi = SIGN(1.d0, h_phi_mfl(0))
    h_phi_mfl(0:ub_mag) = h_phi_mfl(0:ub_mag)*sign_of_bphi

    if (allocated(dlogbdphi_mfl)) deallocate (dlogbdphi_mfl)
    allocate (dlogbdphi_mfl(ibeg:iend))
    dlogbdphi_mfl(0:ub_mag) = fieldpropagator%mdata%dlogbdphi

  !! Modifications by Andreas F. Martitsch (14.03.2014)
    ! Optional output (necessary for modeling the magnetic rotation)
    if (allocated(dlogbds_mfl)) deallocate (dlogbds_mfl)
    allocate (dlogbds_mfl(ibeg:iend))
    dlogbds_mfl(0:ub_mag) = fieldpropagator%mdata%dlogbds

    if (allocated(bcovar_s_hat_mfl)) deallocate (bcovar_s_hat_mfl)
    allocate (bcovar_s_hat_mfl(ibeg:iend))
    bcovar_s_hat_mfl(0:ub_mag) = fieldpropagator%mdata%bcovar_s_hat

    if (allocated(dbcovar_s_hat_dphi_mfl)) deallocate (dbcovar_s_hat_dphi_mfl)
    allocate (dbcovar_s_hat_dphi_mfl(ibeg:iend))
    dbcovar_s_hat_dphi_mfl(0:ub_mag) = fieldpropagator%mdata%dbcovar_s_hat_dphi
  !! End Modifications by Andreas F. Martitsch (14.03.2014)

    if (modify_bl .EQ. 1) then
        phi_mfl(-2:-1) = phi_mfl(0)
        bhat_mfl(-2:-1) = bhat_changed_l
        geodcu_mfl(-2:-1) = geodcu_mfl(0)
        h_phi_mfl(-2:-1) = h_phi_mfl(0)
        dlogbdphi_mfl(-2:-1) = dlogbdphi_mfl(0)
    !! Modifications by Andreas F. Martitsch (14.03.2014)
        ! Optional output (necessary for modeling the magnetic rotation)
        dlogbds_mfl(-2:-1) = dlogbds_mfl(0)
        bcovar_s_hat_mfl(-2:-1) = bcovar_s_hat_mfl(0)
        dbcovar_s_hat_dphi_mfl(-2:-1) = dbcovar_s_hat_dphi_mfl(0)
    !! End Modifications by Andreas F. Martitsch (14.03.2014)
    end if
    if (modify_br .EQ. 1) then
        phi_mfl(iend - 1:iend) = phi_mfl(ub_mag)
        bhat_mfl(iend - 1:iend) = bhat_changed_r
        geodcu_mfl(iend - 1:iend) = geodcu_mfl(ub_mag)
        h_phi_mfl(iend - 1:iend) = h_phi_mfl(ub_mag)
        dlogbdphi_mfl(iend - 1:iend) = dlogbdphi_mfl(ub_mag)
    !! Modifications by Andreas F. Martitsch (14.03.2014)
        ! Optional output (necessary for modeling the magnetic rotation)
        dlogbds_mfl(iend - 1:iend) = dlogbds_mfl(ub_mag)
        bcovar_s_hat_mfl(iend - 1:iend) = bcovar_s_hat_mfl(ub_mag)
        dbcovar_s_hat_dphi_mfl(iend - 1:iend) = dbcovar_s_hat_dphi_mfl(ub_mag)
    !! End Modifications by Andreas F. Martitsch (14.03.2014)
    end if

    ! Computation of the perturbed quantities without
    ! usage of interfaces (fieldpropagator-structure,...)
    if (isw_qflux_na .ne. 0) then

        if (allocated(bnoverb0)) deallocate (bnoverb0)
        allocate (bnoverb0(ibeg:iend))

        if (allocated(dbnoverb0_dphi_mfl)) deallocate (dbnoverb0_dphi_mfl)
        allocate (dbnoverb0_dphi_mfl(ibeg:iend))

        call calc_bnoverb0_arr(phi_mfl, ibeg, iend, bnoverb0, dbnoverb0_dphi_mfl)
        call calc_ntv_output(phi_mfl, bhat_mfl, bnoverb0, ibeg, iend, &
                             eps_M_2_val, av_inv_bhat_val, av_gphph_val)

    end if

    eta_modboundary_l = 1.d0/bhat_mfl(ibeg)
    eta_modboundary_r = 1.d0/bhat_mfl(iend)
    ! allocation
    !  at the moment everything from 1:npart
    !  attention eta is from 0:npart-1
    !
    ! 2-D quantities
    nts_l = (lag + 1)*npass_l
    nts_r = (lag + 1)*npass_r
    if (allocated(amat_plus_plus)) deallocate (amat_plus_plus)
    allocate (amat_plus_plus(nts_r, nts_l))
    if (allocated(amat_minus_minus)) deallocate (amat_minus_minus)
    allocate (amat_minus_minus(nts_l, nts_r))
    if (allocated(amat_plus_minus)) deallocate (amat_plus_minus)
    allocate (amat_plus_minus(nts_l, nts_l))
    if (allocated(amat_minus_plus)) deallocate (amat_minus_plus)
    allocate (amat_minus_plus(nts_r, nts_r))

    if (allocated(source_p)) deallocate (source_p)
    allocate (source_p(nts_r, 3))
    if (allocated(source_m)) deallocate (source_m)
    allocate (source_m(nts_l, 3))

    if (allocated(flux_p)) deallocate (flux_p)
    allocate (flux_p(3, nts_l))
    if (allocated(flux_m)) deallocate (flux_m)
    allocate (flux_m(3, nts_r))

    if (allocated(qflux)) deallocate (qflux)
    allocate (qflux(3, 3))

    print *, 'propagator tag         ', fieldpropagator%tag

    if (solver_talk .EQ. 1 .and. ispec .eq. 0) then
        print *, ' '
        print *, 'I am in ripple_solver'
        print *, ' '
        print *, 'fieldpropagator tag    ', fieldpropagator%tag
        print *, 'fieldripple tag        ', fieldpropagator%ch_act%tag
        print *, ' fieldprop first last  ', fieldpropagator%ch_act%pa_fir%tag, fieldpropagator%ch_act%pa_las%tag
        print *, ' b_prop left,right,min ', b_prop_l, b_prop_r, b_prop_min
        print *, '                       ', fieldpropagator%mdata%bhat(0), fieldpropagator%mdata%bhat(ub_mag)
        print *, ' b_max left,right,min  ', b_max_l, b_max_r, b_min
        print *, ' width                 ', width
        print *, 'bin_split_mode         ', bin_split_mode
        print *, 'phibeg,phiend          ', phibeg, phiend
        print *, '                       ', fieldpropagator%coords%x2(0), fieldpropagator%coords%x2(ub_mag)
        print *, 'ub_mag                 ', ub_mag
        print *, 'ignore_lb,ignore_rb    ', ignore_lb, ignore_rb
        print *, 'ibeg,iend              ', ibeg, iend
        print *, 'ub_eta                 ', ub_eta
        print *, 'npart                  ', npart
        print *, 'npass_l,npass_r        ', npass_l, npass_r
        print *, 'eta_l,eta_r            ', eta_l, eta_r
        print *, 'xetami,xetama          ', xetami, xetama
        print *, 'rt0                    ', rt0
        print *, 'collpar,conl_over_mfp  ', collpar, conl_over_mfp

        ! Optional output (necessary for modeling the magnetic rotation)
        open (123, file='bhat_mfl.dat')
        do i = ibeg, iend
            write (123, *) phi_mfl(i), bhat_mfl(i), REAL(geodcu_mfl(i)), h_phi_mfl(i), dlogbdphi_mfl(i), &
                dlogbds_mfl(i), bcovar_s_hat_mfl(i), dbcovar_s_hat_dphi_mfl(i), &
                REAL(bnoverb0(i)), AIMAG(bnoverb0(i)), REAL(dbnoverb0_dphi_mfl(i)), AIMAG(dbnoverb0_dphi_mfl(i))
        end do
        close (123)

        open (123, file='eta.dat')
        do i = 0, ub_eta
            write (123, *) phi_mfl(ibeg), eta(i)
            write (123, *) phi_mfl(iend), eta(i)
            write (123, *) ' '
        end do
        close (123)

    end if

    !------------------------------------------------------------------------
    ! SERGEI
    !------------------------------------------------------------------------

! Check for axisymmetry:

    if ((isw_axisymm .EQ. 1) .AND. (npass_l .NE. npass_r)) then
        print *, 'ERROR in ripple_solver: cannot run axisymmetric mode, sizes do not fit'
        print *, npass_l, ' =/= ', npass_r
        print *, '  return to calling function.'
        ierr = 1
        return
    end if

    iplot = prop_ripple_plot

! Preparation of coefficients for the kinetic equation solver

    allocate (deriv_coef(4, 0:npart + 1))
    allocate (fun_coef(4, 0:npart + 1))
    allocate (enu_coef(4, npart + 1))                                    !!!term[1]
    allocate (enu_coef2(4, npart + 1))                                   !!!NTV
    allocate (alambd(0:npart + 3, ibeg:iend), Vg_vp_over_B(0:npart, ibeg:iend))
    allocate (scalprod_pleg(0:lag, 0:legmax))                          !!!term[3]
    allocate (alampow(legmax + 1, 0:npart + 1))                            !!!terms[2,3]
    allocate (vrecurr(0:legmax, 0:3, 1:npart + 1))                        !!!terms[2,3]
    allocate (dellampow(4, 1:npart + 1))                                 !!!terms[1-3]
    allocate (dellampow2(4, 1:npart + 1))                                !!!NTV
    allocate (convol_polpow(0:legmax, 1:npart + 3))                      !!!terms[2,3]
    allocate (pleg_bra(0:legmax, 1:npart + 1, ibeg:iend))                 !!!terms[2,3]
    allocate (pleg_ket(0:legmax, 1:npart + 1, ibeg:iend))                 !!!terms[2,3]
    allocate (npl(ibeg:iend))
    allocate (rhs_mat_fzero(4, ibeg:iend, 0:1))
    allocate (rhs_mat_lorentz(5, npart + 1, ibeg:iend))
    allocate (ttmp_mat(5, npart + 1, ibeg:iend))
    allocate (rhs_mat_energ(4, npart + 1, ibeg:iend))
    allocate (rhs_mat_energ2(4, npart + 1, ibeg:iend))              !NTV
    allocate (q_rip(npart + 2, ibeg:iend, 0:2))
    allocate (q_rip_1(npart + 2, ibeg:iend))
    allocate (q_rip_incompress(npart + 2, ibeg:iend))
    allocate (q_rip_parflow(npart + 2, ibeg:iend))
    allocate (convol_flux(npart + 1, ibeg:iend), convol_curr(npart + 1, ibeg:iend))
    allocate (convol_flux_0(npart + 1, ibeg:iend))
    allocate (ind_start(ibeg:iend))
    if (iplot .EQ. 1) allocate (derivs_plot(0:3, 4, npart + 1, ibeg:iend))

!tmp:
    aiota = boozer_iota
! boozer_curr_tor, boozer_curr_pol, boozer_psi_pr,
! boozer_sqrtg11 and boozer_isqrg are now converted
! to cgs-units within neo_magfie.
! This step requires changes within rhs_kin.f90 and
! ripple_solver.f90!
    bcovar_theta_hat = boozer_curr_tor_hat
    bcovar_phi_hat = boozer_curr_pol_hat
    dbcovar_theta_hat_ds = boozer_curr_tor_hat_s
    dbcovar_phi_hat_ds = boozer_curr_pol_hat_s

    allocate (bnoverb0_test(ibeg:iend))
    allocate (dbnoverb0_dphi_mfl_test(ibeg:iend))
    bnoverb0_test = bnoverb0_test_val*EXP(imun*(m_theta*aiota + m_phi)*phi_mfl)
    dbnoverb0_dphi_mfl_test = imun*(m_theta*aiota + m_phi)*bnoverb0_test
    deallocate (bnoverb0_test)
    deallocate (dbnoverb0_dphi_mfl_test)

! boozer_curr_tor, boozer_curr_pol, boozer_psi_pr,
! boozer_sqrtg11 and boozer_isqrg are now converted
! to cgs-units within neo_magfie.
! This step requires changes within rhs_kin.f90 and
! ripple_solver.f90!
!             scalefac_kG=$B_{ref}/(\iota \psi_{tor}^a)$
!This was the old version (conversion to cgs is done here):
!scalefac_kG=1d-4*bmod0/(aiota*boozer_psi_pr)
!Now the quantities are already converted within neo_magfie:

    scalefac_kG = 1.0d0/(aiota*boozer_psi_pr_hat)

    scalefac_kG = scalefac_kG*sign(1.d0, boozer_psi_pr_hat*bcovar_phi_hat*boozer_isqrg)

    rotfactor = imun*m_phi

!end tmp

! Compute coefficients of Legendre polynomials of the order 0,...,legmax:
    call polleg(legmax, coefleg)
! coefleg(l,k) - coefficient of $x^k$ of polynomial $P_l(x)$

    q_rip(:, :, 0) = 0.d0
    do i = 1, npart
        q_rip(i, :, 2) = eta(i) - eta(i - 1)
    end do

    ndim = 4
    allocate (amat(ndim, ndim), bvec_lapack(ndim, ndim), ipivot(ndim))

    npart_loc = 0
    subsqmin = 1.d5*EPSILON(1.d0)

    allocate (delt_pos(ibeg:iend), delt_neg(ibeg:iend))
    allocate (fact_pos_b(ibeg:iend), fact_neg_b(ibeg:iend))
    allocate (fact_pos_e(ibeg:iend), fact_neg_e(ibeg:iend))

    nreal = 1
    ncomp = 1
    allocate (arr_real(ibeg:iend, nreal), arr_comp(ibeg:iend, ncomp))
    arr_real(:, 1) = h_phi_mfl
    arr_comp(:, 1) = geodcu_mfl

    call rearrange_phideps(ibeg, iend, npart, ncomp, nreal, subsqmin, phi_divide, &
                           phi_mfl, bhat_mfl, arr_real, arr_comp, eta, &
                           delt_pos, delt_neg, &
                           fact_pos_b, fact_neg_b, fact_pos_e, fact_neg_e)

    h_phi_mfl = arr_real(:, 1)
    geodcu_mfl = arr_comp(:, 1)
    deallocate (arr_real, arr_comp)

    ! Consistency check of phi_divide, ifnot sucessful, return from this
    ! subroutine.
    if (maxval(phi_divide) > 1) then
        ierr = 3
        deallocate (deriv_coef, npl)
        deallocate (rhs_mat_lorentz, rhs_mat_energ)
        deallocate (fun_coef, ttmp_mat)
        deallocate (rhs_mat_energ2)        !NTV
        deallocate (q_rip, q_rip_1, q_rip_incompress, q_rip_parflow)
        deallocate (convol_flux, convol_curr, convol_flux_0)
        deallocate (pleg_bra, pleg_ket, scalprod_pleg)
        write (*, *) 'ERROR: maxval(phi_divide) > 1, returning.'
        return
    end if

    do istep = ibeg, iend

! semi-levels

        eta0 = 1.d0/bhat_mfl(istep)

        do i = 0, npart
            subsq = 1.d0 - bhat_mfl(istep)*eta(i)
            if (subsq .GT. subsqmin) then
                npassing = i
                alambd(i, istep) = SQRT(subsq)
                Vg_vp_over_B(i, istep) = alambd(i, istep)*eta0/h_phi_mfl(istep) &
                                         *(4.d0*eta0 - eta(i))/3.d0
            else
                alambd(i, istep) = 0.d0
                Vg_vp_over_B(i, istep) = 0.d0
            end if
        end do

        alambd(npart + 1:npart + 3, istep) = 0.d0
        alambd(npassing + 1, istep) = 0.d0
        alambd(npassing + 2, istep) = -alambd(npassing, istep)
        alambd(npassing + 3, istep) = -alambd(npassing - 1, istep)

        npl(istep) = npassing

        npart_loc = MAX(npart_loc, npassing)

        if (istep .EQ. ibeg) then
            npass_l = npassing + 1
        elseif (istep .EQ. iend) then
            npass_r = npassing + 1
        end if

    end do

! compute starting index for 2D vectors

    ind_start(ibeg) = 0

    do istep = ibeg, iend - 1
        ind_start(istep + 1) = ind_start(istep) + 2*(lag + 1)*(npl(istep) + 1)
    end do

    n_2d_size = ind_start(iend) + 2*(lag + 1)*(npl(iend) + 1)

    if (allocated(f_init_arnoldi)) deallocate (f_init_arnoldi)
    allocate (f_init_arnoldi(n_2d_size))

    if (allocated(alam_l)) deallocate (alam_l)
    if (allocated(alam_r)) deallocate (alam_r)
    if (allocated(delta_eta_l)) deallocate (delta_eta_l)
    if (allocated(delta_eta_r)) deallocate (delta_eta_r)
    allocate (delta_eta(npart_loc), alam_l(npass_l), alam_r(npass_r))
    allocate (delta_eta_l(npass_l), delta_eta_r(npass_r))
    delta_eta = eta(1:npart_loc) - eta(0:npart_loc - 1)
    do i = 1, npass_l - 1
        alam_l(i) = SQRT(1.d0 - 0.5d0*(eta(i - 1) + eta(i))*bhat_mfl(ibeg))
        delta_eta_l(i) = eta(i) - eta(i - 1)
    end do
    i = npass_l
    alam_l(i) = SQRT(1.d0 - 0.5d0*(eta(i - 1)*bhat_mfl(ibeg) + 1.d0))
    delta_eta_l(i) = 1.d0/bhat_mfl(ibeg) - eta(i - 1)
    do i = 1, npass_r - 1
        alam_r(i) = SQRT(1.d0 - 0.5d0*(eta(i - 1) + eta(i))*bhat_mfl(iend))
        delta_eta_r(i) = eta(i) - eta(i - 1)
    end do
    i = npass_r
    alam_r(i) = SQRT(1.d0 - 0.5d0*(eta(i - 1)*bhat_mfl(iend) + 1.d0))
    delta_eta_r(i) = 1.d0/bhat_mfl(iend) - eta(i - 1)

! Calculation of the ODE coefficients

    do istep = ibeg, iend

        npassing = npl(istep)
        eta0 = 1.d0/bhat_mfl(istep)
        amin2ovb = -2.d0/bhat_mfl(istep)
        coefdir = 0.5*collpar/h_phi_mfl(istep)
        coefenu_averb = 0.5d0*collpar/h_phi_mfl(istep)           !!!term[1]
        coefenu = -coefenu_averb*2.d0/bhat_mfl(istep)             !!!term[1]

!-----------------------------
! begin terms[2,3]

! here dellampow(m,n)=$(\lambda_{n-1}-\lambda_{n})^m$
        dellampow(1, 1:npassing + 1) &
            = alambd(0:npassing, istep) - alambd(1:npassing + 1, istep)
        do k = 2, 4
            km1 = k - 1
            dellampow(k, 1:npassing + 1) = dellampow(km1, 1:npassing + 1) &
                                           *dellampow(1, 1:npassing + 1)
        end do

        alampow(1, 0:npassing + 1) = alambd(0:npassing + 1, istep)
        do k = 2, legmax + 1
            km1 = k - 1
            alampow(k, 0:npassing + 1) = alambd(0:npassing + 1, istep) &
                                         *alampow(km1, 0:npassing + 1)
        end do

        do k = 1, legmax + 1
! Caution:
! here power index is shifted (instead of term (k) -> term (k-1) is computed)
            km1 = k - 1
            vrecurr(km1, 0, 1:npassing + 1) &
                = (alampow(k, 0:npassing) - alampow(k, 1:npassing + 1))/DBLE(k)
            do m = 1, 3
                vrecurr(km1, m, 1:npassing + 1) &
                    = (alampow(k, 0:npassing)*dellampow(m, 1:npassing + 1) &
                       - DBLE(m)*alambd(1:npassing + 1, istep)*vrecurr(km1, m - 1, 1:npassing + 1)) &
                      /DBLE(k + m)
            end do
! divide by factorial
            mfactorial = 1
            do m = 1, 3
                mfactorial = mfactorial*m
                vrecurr(km1, m, 1:npassing + 1) = vrecurr(km1, m, 1:npassing + 1) &
                                                  /DBLE(mfactorial)
            end do
        end do

! re-definition: here dellampow(m,n)=$(\lambda_{n-1}-\lambda_{n})^m/m!$
! (divided by factorial)

        mfactorial = 1
        do m = 2, 4
            mfactorial = mfactorial*m
            dellampow(m, 1:npassing + 1) = dellampow(m, 1:npassing + 1)/DBLE(mfactorial)
        end do

! new stuff: NTV

        do m = 1, 4
            dellampow2(m, 1:npassing + 1) = dellampow(m, 1:npassing + 1) &
                                            *(alambd(1:npassing + 1, istep)**2 &
                                              + 2.d0*alambd(1:npassing + 1, istep)*dellampow(1, 1:npassing + 1) &
                                              *real(m, kind=kind(0d0))/real(m + 1, kind=kind(0d0)) &
                                              + dellampow(1, 1:npassing + 1)**2*real(m, kind=kind(0d0))/real(m + 2, kind=kind(0d0)))
        end do

! end new stuff: NTV

! term[2] (Legendre polynomials) -  ket-vector
! Caution: ket-vector cooresponds do discretization of
! P_l(lambda)/|lambda|, not P_l(lambda). Scalar prduct with bra-vector
! does not mean the product of Legendre polynomials!
! even powers:
        do m = 0, legmax, 2
            pleg_ket(m, 1:npassing + 1, istep) = amin2ovb*coefleg(m, 0) &
                                                 *(alampow(1, 1:npassing + 1) - alampow(1, 0:npassing))
            do k = 2, m, 2
                kp1 = k + 1
                pleg_ket(m, 1:npassing + 1, istep) &
                    = pleg_ket(m, 1:npassing + 1, istep) + amin2ovb*coefleg(m, k) &
                      *(alampow(kp1, 1:npassing + 1) - alampow(kp1, 0:npassing))/DBLE(kp1)
            end do
        end do
! odd powers:
        do m = 1, legmax, 2
            pleg_ket(m, 1:npassing + 1, istep) = amin2ovb*coefleg(m, 1) &
                                                 *(alampow(2, 1:npassing + 1) - alampow(2, 0:npassing))/2.d0
            do k = 3, m, 2
                kp1 = k + 1
                pleg_ket(m, 1:npassing + 1, istep) &
                    = pleg_ket(m, 1:npassing + 1, istep) + amin2ovb*coefleg(m, k) &
                      *(alampow(kp1, 1:npassing + 1) - alampow(kp1, 0:npassing))/DBLE(kp1)
            end do
        end do

        convol_polpow = 0.d0

! end terms[2,3]
!---------------------------------------

        do i = 1, npassing + 1

            i1min = MAX(0, i - 2)

            kmax = 5

            do k = 1, kmax
                i1 = k - 1 + i1min
                diflam = alambd(i1, istep) - alambd(i, istep)
                diflampow = diflam
                alp(k) = (alambd(i, istep) + diflam/2.d0)*diflampow
                diflampow = diflam*diflampow
                bet(k) = (alambd(i, istep)/2.d0 + diflam/3.d0)*diflampow
                diflampow = diflam*diflampow
                gam(k) = (alambd(i, istep)/6.d0 + diflam/8.d0)*diflampow
                diflampow = diflam*diflampow
                del(k) = (alambd(i, istep)/24.d0 + diflam/30.d0)*diflampow
            end do

            do k = 1, 4
                amat(k, 1) = (alp(k + 1) - alp(k))*amin2ovb
                amat(k, 2) = (bet(k + 1) - bet(k))*amin2ovb
                amat(k, 3) = (gam(k + 1) - gam(k))*amin2ovb
                amat(k, 4) = (del(k + 1) - del(k))*amin2ovb
            end do

            if (i .EQ. npassing) then
                amat(4, :) = -amat(4, :)
            elseif (i .EQ. npassing + 1) then
                amat(3, :) = -amat(3, :)
                amat(4, :) = -amat(4, :)
            end if

            bvec_lapack = 0.d0
            do k = 1, ndim
                bvec_lapack(k, k) = 1.d0
            end do

            call gbsv(ndim, ndim, amat, ipivot, bvec_lapack, info)

!> bvec_lapack(j,k) - contribution to the derivative of the distribution
!> function $\hat f^\sigma$ of the order j-1=0,1,2,3 at the boundary
!> $\lambda=\lambda_i$ (at the level $\eta=\eta_i$) from the band i+k-2,
!> where k=1,2,3,4. ifi=1 contributing bands are i+k-1=1,2,3,4 (shift up by 1).
!> ifi=npassing, sigma=-1 fluxes start contribute:
!> contributions for k=1,2,3,4 come from fun(N-1),fun(N),fun(N+1),fun_b(N*1)
!> ifi=npassing+1
!> contributions for k=1,2,3,4 come from fun(N),fun(N+1),fun_b(N+1),fun_b(N)
!> Actual derivative can be obtained by summation of corresponding
!> band-integrated fluxes, $f_{i+k-2}$, multiplied with these contributions

            if (iplot .EQ. 1) derivs_plot(0:3, 1:4, i, istep) = bvec_lapack(1:4, 1:4)
            deriv_coef(:, i) = bvec_lapack(2, :)*coefdir*MIN(eta(i), eta0)
            fun_coef(:, i) = bvec_lapack(1, :)*MIN(eta(i), eta0)

            enu_coef(:, i) = MATMUL(dellampow(:, i), bvec_lapack)*coefenu      !!!term[1]
            enu_coef2(:, i) = MATMUL(dellampow2(:, i), bvec_lapack)*coefenu    !!!NTV

            if (i .EQ. 1) then
                convol_polpow(0:legmax, i:i + 3) = convol_polpow(0:legmax, i:i + 3) &
                                                   + MATMUL(vrecurr(0:legmax, 0:3, i), bvec_lapack)       !!!term[3]
            else
                convol_polpow(0:legmax, i - 1:i + 2) = convol_polpow(0:legmax, i - 1:i + 2) &
                                                       + MATMUL(vrecurr(0:legmax, 0:3, i), bvec_lapack)       !!!term[3]
            end if
            if (i .EQ. npassing + 1) then
! distribution function at the trapped-passing boundary:
! f0=sum(fun(npassing:npassing+3)*rhs_mat_fzero(:,istep,0))
                rhs_mat_fzero(:, istep, 0) = bvec_lapack(1, :)
            end if

        end do

! Eliminate stepping over the boundary:

        do k = 0, legmax, 2
            convol_polpow(k, npassing) = convol_polpow(k, npassing) &
                                         + convol_polpow(k, npassing + 3)
            convol_polpow(k, npassing + 1) = convol_polpow(k, npassing + 1) &
                                             + convol_polpow(k, npassing + 2)
        end do

        do k = 1, legmax, 2
            convol_polpow(k, npassing) = convol_polpow(k, npassing) &
                                         - convol_polpow(k, npassing + 3)
            convol_polpow(k, npassing + 1) = convol_polpow(k, npassing + 1) &
                                             - convol_polpow(k, npassing + 2)
        end do

! term[3] (Legendre polynomials) -  bra-vector
! numbering of levels (N=npassing): 1-N - $f_n$, N+1 - $f^b$, N+2 - $f^a$
! even powers:
        do m = 0, legmax, 2
            pleg_bra(m, 1:npassing + 1, istep) = coefleg(m, 0) &
                                                 *convol_polpow(0, 1:npassing + 1)
            do k = 2, m, 2
                pleg_bra(m, 1:npassing + 1, istep) &
                    = pleg_bra(m, 1:npassing + 1, istep) + coefleg(m, k) &
                      *convol_polpow(k, 1:npassing + 1)
            end do
        end do
! odd powers:
        do m = 1, legmax, 2
            pleg_bra(m, 1:npassing + 1, istep) = coefleg(m, 1) &
                                                 *convol_polpow(1, 1:npassing + 1)
            do k = 3, m, 2
                pleg_bra(m, 1:npassing + 1, istep) &
                    = pleg_bra(m, 1:npassing + 1, istep) + coefleg(m, k) &
                      *convol_polpow(k, 1:npassing + 1)
            end do
        end do

        pleg_bra(0:legmax, 1:npassing + 1, istep) &
            = pleg_bra(0:legmax, 1:npassing + 1, istep)*coefenu_averb

        coef_cf = (1.d0, 0.d0)/bhat_mfl(istep)**2/h_phi_mfl(istep)
        convol_flux_0(1:npassing + 1, istep) &
            = (convol_polpow(0, 1:npassing + 1) + convol_polpow(2, 1:npassing + 1)) &
              *coef_cf

! levels

        rhs_mat_lorentz(5, 1:npassing + 1, istep) = 0.d0

        rhs_mat_lorentz(1:4, 1, istep) = deriv_coef(:, 1)
        rhs_mat_lorentz(1:4, 2, istep) = deriv_coef(:, 2) - deriv_coef(:, 1)
        rhs_mat_lorentz(1:4, 3:npassing + 1, istep) &
            = -deriv_coef(:, 2:npassing)
        rhs_mat_lorentz(2:5, 3:npassing + 1, istep) &
            = rhs_mat_lorentz(2:5, 3:npassing + 1, istep) &
              + deriv_coef(:, 3:npassing + 1)

        ttmp_mat(5, 1:npassing + 1, istep) = 0.d0

        ttmp_mat(1:4, 1, istep) = fun_coef(:, 1)
        ttmp_mat(1:4, 2, istep) = fun_coef(:, 2) - fun_coef(:, 1)
        ttmp_mat(1:4, 3:npassing + 1, istep) = -fun_coef(:, 2:npassing)
        ttmp_mat(2:5, 3:npassing + 1, istep) &
            = ttmp_mat(2:5, 3:npassing + 1, istep) &
              + fun_coef(:, 3:npassing + 1)

        rhs_mat_fzero(:, istep, 1) = deriv_coef(:, npassing + 1)

! Change of the boundary layer width

! begin term[1]:

        rhs_mat_energ(:, 1:npassing + 1, istep) = enu_coef(:, 1:npassing + 1)
        rhs_mat_energ2(:, 1:npassing + 1, istep) = enu_coef2(:, 1:npassing + 1) !NTV

! end term[1]

        q_rip_1(1:npassing, istep) &
            = Vg_vp_over_B(1:npassing, istep) - Vg_vp_over_B(0:npassing - 1, istep)
        q_rip_1(npassing + 1, istep) = -Vg_vp_over_B(npassing, istep)
        q_rip(npassing + 1, istep, 2) = eta0 - eta(npassing)
        q_rip(1:npassing + 1, istep, 2) &
            = q_rip(1:npassing + 1, istep, 2) &
              *bhat_mfl(istep)/h_phi_mfl(istep)*sign_of_bphi
        q_rip_incompress(1:npassing, istep) &
            = (alambd(0:npassing - 1, istep)**3 - alambd(0:npassing - 1, istep) &
               - alambd(1:npassing, istep)**3 + alambd(1:npassing, istep)) &
              *dlogbdphi_mfl(istep)
        q_rip_incompress(npassing + 1, istep) &
            = (alambd(npassing, istep)**3 - alambd(npassing, istep)) &
              *dlogbdphi_mfl(istep)
        q_rip_parflow(1:npassing, istep) = 2.d0/3.d0 &
                                           *(alambd(0:npassing - 1, istep)**3 - alambd(1:npassing, istep)**3)
        q_rip_parflow(npassing + 1, istep) = 2.d0/3.d0*alambd(npassing, istep)**3

        convol_curr(1:npassing + 1, istep) = bhat_mfl(istep)/h_phi_mfl(istep)*sign_of_bphi

    end do

    deallocate (amat, bvec_lapack, ipivot)

! NEO-2 can treat now multiple species
! (move collpar from pleg_bra to pleg_ket to avoid mixing up
! of species-dependent parameters)
    pleg_bra = pleg_bra/collpar
    pleg_ket = pleg_ket*collpar

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

! Matrix size:

    nrow = n_2d_size
    ncol = n_2d_size

! Compute vectors for convolution of fluxes and source vectors:

    allocate (flux_vector(3, n_2d_size), source_vector(n_2d_size, 4), bvec_parflow(n_2d_size))
    allocate (flux_vector_plot(3, n_2d_size))

    ! Use pre-conditioned iterations:
    ! -> remove null-space of axisymmetric solution (energy conservation)
    allocate (energvec_ket(n_2d_size), energvec_bra(n_2d_size))
    allocate (densvec_bra(n_2d_size), densvec_ket(n_2d_size))
    energvec_ket = 0.d0
    energvec_bra = 0.d0
    denom_energ = 0.d0
    ! End Use pre-conditioned iterations

    if (isw_lorentz .EQ. 1) then
        x1mm(0, 0) = 1.d0
        x2mm(0, 0) = 1.d0
    end if

! Determine the size of arrays (number of non-zero elements):

    nz = 0
    nz_coll = 0
    nz_ttmp = 0
    nz_regper = 0

! Co-passing: sigma=1

    istep = ibeg
    npassing = npl(istep)

! entry:

    do m = 0, lag
        k = ind_start(istep) + 2*(npassing + 1)*m

        do ipart = 1, npassing + 1
            nz = nz + 1
!      irow(nz)=k+ipart
!      icol(nz)=k+ipart
!      amat_sp(nz)=(1.d0,0.d0)
        end do

        if (isw_axisymm .EQ. 1) then
! periodicity:
            k_prev = ind_start(iend) + 2*(npassing + 1)*m

            do ipart = 1, npassing + 1
                nz = nz + 1
!        irow(nz)=k+ipart
!        icol(nz)=k_prev+ipart
!        amat_sp(nz)=(-1.d0,0.d0)
            end do

!! Modifications by Andreas F. Martitsch (12.12.2016)
! -> Removal of null-space of axisymmetric equation with
! full linearized collision operator (wrong)
!      IF(isw_regper.EQ.1.AND.m.LE.1) then
! -> Removal of null-space of axisymmetric equation without
! integral part of the collision operator
! (solution to the problem without iterations --> correct)
            if (isw_regper .EQ. 1 .AND. m .LT. 1) then
                do ipart = 1, 1
!          if(ipart.le.npassing) then
!            deleta_factor=(eta(ipart)-eta(ipart-1))*bhat_mfl(istep)
!          else
!            deleta_factor=1.d0-eta(ipart-1)*bhat_mfl(istep)
!          end if

                    do ipart1 = 1, npassing + 1
                        nz_regper = nz_regper + 1
!            irow_regper(nz_regper)=k+ipart
!            icol_regper(nz_regper)=k_prev+ipart1
!            amat_regper(nz_regper)=deleta_factor
                        nz_regper = nz_regper + 1
!            irow_regper(nz_regper)=k+ipart
!            icol_regper(nz_regper)=k+npassing+1+ipart1
!            amat_regper(nz_regper)=deleta_factor
                    end do
!
                end do
            end if

!! End Modifications by Andreas F. Martitsch (12.12.2016)

        end if

    end do

    do istep = ibeg + 1, iend
        npassing_prev = npl(istep - 1)
        npassing = npl(istep)
!    delphim1=1.d0/delt_pos(istep)
!    deloneovb=0.5d0*delphim1*(1.d0/bhat_mfl(istep-1)-1.d0/bhat_mfl(istep))

        do m = 0, lag
            k_prev = ind_start(istep - 1) + 2*(npassing_prev + 1)*m
            k = ind_start(istep) + 2*(npassing + 1)*m

! free flight:

            do ipart = 1, npassing + 1
                nz = nz + 1
!        irow(nz)=k+ipart
!        icol(nz)=k+ipart
!        amat_sp(nz)=delphim1
                if (colltest) nz_ttmp = nz_ttmp + 1
            end do

            do ipart = 1, npassing
                nz = nz + 1
!        irow(nz)=k+ipart
!        icol(nz)=k_prev+ipart
!        amat_sp(nz)=-delphim1
                if (colltest) nz_ttmp = nz_ttmp + 1
            end do

            if (npassing_prev .GE. npassing) then
                nz = nz + 1
!        irow(nz)=k+npassing+1
!        icol(nz)=k_prev+npassing+1
!        amat_sp(nz)=-delphim1
                if (colltest) nz_ttmp = nz_ttmp + 1
            end if

! mirroring:

            if (npassing_prev .EQ. npassing) then

                do kk = 1, 4
                    nz = nz + 1
!          irow(nz)=k+npassing+1
!          icol(nz)=k+npassing+kk-1
!          amat_sp(nz)=deloneovb*rhs_mat_fzero(kk,istep,0)
                    if (colltest) nz_ttmp = nz_ttmp + 1
                end do

                do kk = 1, 4
                    nz = nz + 1
!          irow(nz)=k+npassing+1
!          icol(nz)=k_prev+npassing_prev+kk-1
!          amat_sp(nz)=deloneovb*rhs_mat_fzero(kk,istep-1,0)
                    if (colltest) nz_ttmp = nz_ttmp + 1
                end do

            elseif (npassing_prev .GT. npassing) then
                nz = nz + 1
!        irow(nz)=k_prev+npassing_prev+2
!        icol(nz)=k_prev+npassing_prev+1
!        amat_sp(nz)=-delphim1
                if (colltest) nz_ttmp = nz_ttmp + 1
            end if

! collisions:

            if (fact_pos_e(istep) .NE. 0.d0) then
! Lorentz operator:

                if (isw_lor .EQ. 1) then

                    do ipart = 1, npassing + 1
                        do kk = 1, 5
                            do mm = 0, lag
                                nz = nz + 1
                                nz_coll = nz_coll + 1
!                irow(nz)=k+ipart
!                icol(nz)=k+max(0,ipart-3)+kk+2*(npassing+1)*(mm-m)
!                amat_sp(nz)=anumm(m,mm)*rhs_mat_lorentz(kk,ipart,istep)   &
!                           *fact_pos_e(istep)*0.5d0

                                if (.NOT. colltest .AND. mm .EQ. m) then
                                    nz_ttmp = nz_ttmp + 1
                                end if

                                if (ipart .LE. npassing_prev + 1) then
                                    nz = nz + 1
!                  irow(nz)=k+ipart
!                  icol(nz)=k_prev+max(0,ipart-3)+kk+2*(npassing_prev+1)*(mm-m)
!                  amat_sp(nz)=anumm(m,mm)*rhs_mat_lorentz(kk,ipart,istep-1) &
!                             *fact_pos_e(istep)*0.5d0

                                    if (.NOT. colltest .AND. mm .EQ. m) then
                                        nz_ttmp = nz_ttmp + 1
                                    end if
                                elseif (ipart .gt. npassing_prev + 1 .and. addboucol) then
                                    nz = nz + 1
!                  irow(nz) = k + ipart
!                  icol(nz) = k + ipart + 4 - kk + 2*(npassing+1)*(mm-m)
!                  amat_sp(nz) = anumm(m,mm)*rhs_mat_lorentz(kk,ipart,istep) &
!                             *fact_pos_e(istep)*0.5d0
                                end if
                            end do
                        end do
                    end do

                end if

!        nz_beg=nz+1

! energy diffusion operator:

                if (isw_ene .EQ. 1) then

                    do ipart = 1, npassing + 1
                        do kk = 1, 4
                            do mm = 0, lag
                                nz = nz + 1
                                nz_coll = nz_coll + 1
!                irow(nz)=k+ipart
!                icol(nz)=k+max(0,ipart-2)+kk+2*(npassing+1)*(mm-m)
!                amat_sp(nz)=denmm(m,mm)*rhs_mat_energ(kk,ipart,istep)*0.5d0
                            end do
                        end do

                        if (ipart .LE. npassing_prev + 1) then
                            do kk = 1, 4
                                do mm = 0, lag
                                    nz = nz + 1
!                  irow(nz)=k+ipart
!                  icol(nz)=k_prev+max(0,ipart-2)+kk+2*(npassing_prev+1)*(mm-m)
!                  amat_sp(nz)=denmm(m,mm)*rhs_mat_energ(kk,ipart,istep-1)*0.5d0
                                end do
                            end do
                        elseif (ipart .gt. npassing_prev + 1 .and. addboucol) then
                            do kk = 1, 4
                                do mm = 0, lag
                                    nz = nz + 1
!                  irow(nz) = k + ipart
!                  icol(nz) = k + ipart + 3 - kk + 2*(npassing+1)*(mm-m)
!                  amat_sp(nz) = denmm(m,mm)*rhs_mat_energ(kk,ipart,istep)*0.5d0
                                end do
                            end do
                        end if

                    end do

                end if

!        amat_sp(nz_beg:nz)=fact_pos_e(istep)*amat_sp(nz_beg:nz)

            end if

        end do

    end do

! Counter-passing: sigma=-1

    istep = iend
    npassing = npl(istep)

! entry:

    do m = 0, lag
        k = ind_start(istep) + 2*(npassing + 1)*m + 2*npassing + 3

        do ipart = 1, npassing + 1
            nz = nz + 1
!      irow(nz)=k-ipart
!      icol(nz)=k-ipart
!      amat_sp(nz)=(1.d0,0.d0)
        end do

        if (isw_axisymm .EQ. 1) then
! periodicity:
            k_prev = ind_start(ibeg) + 2*(npassing + 1)*m + 2*npassing + 3

            do ipart = 1, npassing + 1
                nz = nz + 1
!        irow(nz)=k-ipart
!        icol(nz)=k_prev-ipart
!        amat_sp(nz)=(-1.d0,0.d0)
            end do

!! Modifications by Andreas F. Martitsch (12.12.2016)
! -> Removal of null-space of axisymmetric equation with
! full linearized collision operator (wrong)
!      IF(isw_regper.EQ.1.AND.m.LE.1) then
! -> Removal of null-space of axisymmetric equation without
! integral part of the collision operator
! (solution to the problem without iterations --> correct)
            if (isw_regper .EQ. 1 .AND. m .LT. 1) then
                do ipart = 1, 1
!          if(ipart.le.npassing) then
!            deleta_factor=(eta(ipart)-eta(ipart-1))*bhat_mfl(istep)
!          else
!            deleta_factor=1.d0-eta(ipart-1)*bhat_mfl(istep)
!          end if

                    do ipart1 = 1, npassing + 1
                        nz_regper = nz_regper + 1
!            irow_regper(nz_regper)=k-ipart
!            icol_regper(nz_regper)=k_prev-ipart1
!            amat_regper(nz_regper)=deleta_factor
                        nz_regper = nz_regper + 1
!            irow_regper(nz_regper)=k-ipart
!            icol_regper(nz_regper)=k-npassing-1-ipart1
!            amat_regper(nz_regper)=deleta_factor
                    end do

                end do
            end if
!! End Modifications by Andreas F. Martitsch (12.12.2016)

        end if

    end do

    do istep = ibeg, iend - 1
        npassing_prev = npl(istep + 1)
        npassing = npl(istep)
!    delphim1=1.d0/delt_neg(istep)
!    deloneovb=0.5d0*delphim1*(1.d0/bhat_mfl(istep+1)-1.d0/bhat_mfl(istep))

        do m = 0, lag
            k_prev = ind_start(istep + 1) + 2*(npassing_prev + 1)*m + 2*npassing_prev + 3
            k = ind_start(istep) + 2*(npassing + 1)*m + 2*npassing + 3

! free flight:

            do ipart = 1, npassing + 1
                nz = nz + 1
!        irow(nz)=k-ipart
!        icol(nz)=k-ipart
!        amat_sp(nz)=delphim1
                if (colltest) nz_ttmp = nz_ttmp + 1
            end do
!
            do ipart = 1, npassing
                nz = nz + 1
!        irow(nz)=k-ipart
!        icol(nz)=k_prev-ipart
!        amat_sp(nz)=-delphim1
                if (colltest) nz_ttmp = nz_ttmp + 1
            end do

            if (npassing_prev .GE. npassing) then
                nz = nz + 1
!        irow(nz)=k-npassing-1
!        icol(nz)=k_prev-npassing-1
!        amat_sp(nz)=-delphim1
                if (colltest) nz_ttmp = nz_ttmp + 1
            end if

! mirroring:

            if (npassing_prev .EQ. npassing) then

                do kk = 1, 4
                    nz = nz + 1
!          irow(nz)=k-npassing-1
!          icol(nz)=k-npassing-kk+1
!          amat_sp(nz)=deloneovb*rhs_mat_fzero(kk,istep,0)
                    if (colltest) nz_ttmp = nz_ttmp + 1
                end do

                do kk = 1, 4
                    nz = nz + 1
!          irow(nz)=k-npassing-1
!          icol(nz)=k_prev-npassing_prev-kk+1
!          amat_sp(nz)=deloneovb*rhs_mat_fzero(kk,istep+1,0)
                    if (colltest) nz_ttmp = nz_ttmp + 1
                end do
!
            elseif (npassing_prev .GT. npassing) then
                nz = nz + 1
!        irow(nz)=k_prev-npassing_prev-2
!        icol(nz)=k_prev-npassing_prev-1
!        amat_sp(nz)=-delphim1
                if (colltest) nz_ttmp = nz_ttmp + 1
            end if

! collisions:

            if (fact_neg_e(istep) .NE. 0.d0) then
! Lorentz operator:

                if (isw_lor .EQ. 1) then

                    do ipart = 1, npassing + 1
                        do kk = 1, 5
                            do mm = 0, lag
                                nz = nz + 1
                                nz_coll = nz_coll + 1
!                irow(nz)=k-ipart
!                icol(nz)=k-max(0,ipart-3)-kk+2*(npassing+1)*(mm-m)
!                amat_sp(nz)=anumm(m,mm)*rhs_mat_lorentz(kk,ipart,istep) &
!                           *fact_neg_e(istep)*0.5d0

                                if (.NOT. colltest .AND. mm .EQ. m) then
                                    nz_ttmp = nz_ttmp + 1
                                end if

                                if (ipart .LE. npassing_prev + 1) then
                                    nz = nz + 1
!                  irow(nz)=k-ipart
!                  icol(nz)=k_prev-max(0,ipart-3)-kk+2*(npassing_prev+1)*(mm-m)
!                  amat_sp(nz)=anumm(m,mm)*rhs_mat_lorentz(kk,ipart,istep+1) &
!                             *fact_neg_e(istep)*0.5d0

                                    if (.NOT. colltest .AND. mm .EQ. m) then
                                        nz_ttmp = nz_ttmp + 1
                                    end if
                                elseif (ipart .gt. npassing_prev + 1 .and. addboucol) then
                                    nz = nz + 1
!                  irow(nz) = k - ipart
!                  icol(nz) = k - ipart - 4 + kk + 2*(npassing+1)*(mm-m)
!                  amat_sp(nz) = anumm(m,mm)*rhs_mat_lorentz(kk,ipart,istep) &
!                             *fact_neg_e(istep)*0.5d0
                                end if
                            end do
                        end do
                    end do

                end if

!        nz_beg=nz+1

! energy diffusion operator:

                if (isw_ene .EQ. 1) then

                    do ipart = 1, npassing + 1
                        do kk = 1, 4
                            do mm = 0, lag
                                nz = nz + 1
                                nz_coll = nz_coll + 1
!                irow(nz)=k-ipart
!                icol(nz)=k-max(0,ipart-2)-kk+2*(npassing+1)*(mm-m)
!                amat_sp(nz)=denmm(m,mm)*rhs_mat_energ(kk,ipart,istep)*0.5d0
                            end do
                        end do

                        if (ipart .LE. npassing_prev + 1) then
                            do kk = 1, 4
                                do mm = 0, lag
                                    nz = nz + 1
!                  irow(nz)=k-ipart
!                  icol(nz)=k_prev-max(0,ipart-2)-kk+2*(npassing_prev+1)*(mm-m)
!                  amat_sp(nz)=denmm(m,mm)*rhs_mat_energ(kk,ipart,istep+1)*0.5d0
                                end do
                            end do
                        elseif (ipart .gt. npassing_prev + 1 .and. addboucol) then
                            do kk = 1, 4
                                do mm = 0, lag
                                    nz = nz + 1
!                  irow(nz) = k - ipart
!                  icol(nz) = k - ipart - 3 + kk + 2*(npassing+1)*(mm-m)
!                  amat_sp(nz) = denmm(m,mm)*rhs_mat_energ(kk,ipart,istep)*0.5d0
                                end do
                            end do
                        end if

                    end do

                end if

!        amat_sp(nz_beg:nz)=fact_neg_e(istep)*amat_sp(nz_beg:nz)

            end if

        end do

    end do

    allocate (irow(nz), icol(nz), amat_sp(nz))
    allocate (irow_coll(nz_coll), icol_coll(nz_coll), amat_coll(nz_coll))
    allocate (irow_ttmp(nz_ttmp), icol_ttmp(nz_ttmp), amat_ttmp(nz_ttmp))
    allocate (irow_regper(nz_regper), icol_regper(nz_regper), amat_regper(nz_regper))

! Fill the arrays:

    nz = 0
    nz_coll = 0
    nz_ttmp = 0
    nz_regper = 0

! Co-passing: sigma=1

    istep = ibeg
    npassing = npl(istep)

! entry:

    do m = 0, lag
        k = ind_start(istep) + 2*(npassing + 1)*m

        do ipart = 1, npassing + 1
            nz = nz + 1
            irow(nz) = k + ipart
            icol(nz) = k + ipart
            amat_sp(nz) = (1.d0, 0.d0)
        end do

        if (isw_axisymm .EQ. 1) then
! periodicity:
            k_prev = ind_start(iend) + 2*(npassing + 1)*m

            do ipart = 1, npassing + 1
                nz = nz + 1
                irow(nz) = k + ipart
                icol(nz) = k_prev + ipart
                amat_sp(nz) = (-1.d0, 0.d0)
            end do

!! Modifications by Andreas F. Martitsch (12.12.2016)
! -> Removal of null-space of axisymmetric equation with
! full linearized collision operator (wrong)
!      IF(isw_regper.EQ.1.AND.m.LE.1) then
! -> Removal of null-space of axisymmetric equation without
! integral part of the collision operator
! (solution to the problem without iterations --> correct)
            if (isw_regper .EQ. 1 .AND. m .LT. 1) then
                do ipart = 1, 1
                    if (ipart .LE. npassing) then
                        deleta_factor = (eta(ipart) - eta(ipart - 1))*bhat_mfl(istep)
                    else
                        deleta_factor = 1.d0 - eta(ipart - 1)*bhat_mfl(istep)
                    end if

                    do ipart1 = 1, npassing + 1
                        nz_regper = nz_regper + 1
                        irow_regper(nz_regper) = k + ipart
                        icol_regper(nz_regper) = k_prev + ipart1
                        amat_regper(nz_regper) = deleta_factor
                        nz_regper = nz_regper + 1
                        irow_regper(nz_regper) = k + ipart
                        icol_regper(nz_regper) = k + npassing + 1 + ipart1
                        amat_regper(nz_regper) = deleta_factor
                    end do

                end do
            end if
!! End Modifications by Andreas F. Martitsch (12.12.2016)

        end if

    end do

    do istep = ibeg + 1, iend
        npassing_prev = npl(istep - 1)
        npassing = npl(istep)
        delphim1 = 1.d0/delt_pos(istep)
        deloneovb = 0.5d0*delphim1*(1.d0/bhat_mfl(istep - 1) - 1.d0/bhat_mfl(istep))

        do m = 0, lag
            k_prev = ind_start(istep - 1) + 2*(npassing_prev + 1)*m
            k = ind_start(istep) + 2*(npassing + 1)*m

! free flight:

            do ipart = 1, npassing + 1
                nz = nz + 1
                irow(nz) = k + ipart
                icol(nz) = k + ipart
                amat_sp(nz) = delphim1
                if (colltest) then
                    nz_ttmp = nz_ttmp + 1
                    irow_ttmp(nz_ttmp) = irow(nz)
                    icol_ttmp(nz_ttmp) = icol(nz)
                    amat_ttmp(nz_ttmp) = REAL(amat_sp(nz), dp)
                end if
            end do

            do ipart = 1, npassing
                nz = nz + 1
                irow(nz) = k + ipart
                icol(nz) = k_prev + ipart
                amat_sp(nz) = -delphim1
                if (colltest) then
                    nz_ttmp = nz_ttmp + 1
                    irow_ttmp(nz_ttmp) = irow(nz)
                    icol_ttmp(nz_ttmp) = icol(nz)
                    amat_ttmp(nz_ttmp) = REAL(amat_sp(nz), dp)
                end if
            end do

            if (npassing_prev .GE. npassing) then
                nz = nz + 1
                irow(nz) = k + npassing + 1
                icol(nz) = k_prev + npassing + 1
                amat_sp(nz) = -delphim1
                if (colltest) then
                    nz_ttmp = nz_ttmp + 1
                    irow_ttmp(nz_ttmp) = irow(nz)
                    icol_ttmp(nz_ttmp) = icol(nz)
                    amat_ttmp(nz_ttmp) = REAL(amat_sp(nz), dp)
                end if
            end if

! mirroring:

            if (npassing_prev .EQ. npassing) then

                do kk = 1, 4
                    nz = nz + 1
                    irow(nz) = k + npassing + 1
                    icol(nz) = k + npassing + kk - 1
                    amat_sp(nz) = deloneovb*rhs_mat_fzero(kk, istep, 0)
                    if (colltest) then
                        nz_ttmp = nz_ttmp + 1
                        irow_ttmp(nz_ttmp) = irow(nz)
                        icol_ttmp(nz_ttmp) = icol(nz)
                        amat_ttmp(nz_ttmp) = REAL(amat_sp(nz), dp)
                    end if
                end do

                do kk = 1, 4
                    nz = nz + 1
                    irow(nz) = k + npassing + 1
                    icol(nz) = k_prev + npassing_prev + kk - 1
                    amat_sp(nz) = deloneovb*rhs_mat_fzero(kk, istep - 1, 0)
                    if (colltest) then
                        nz_ttmp = nz_ttmp + 1
                        irow_ttmp(nz_ttmp) = irow(nz)
                        icol_ttmp(nz_ttmp) = icol(nz)
                        amat_ttmp(nz_ttmp) = REAL(amat_sp(nz), dp)
                    end if
                end do

            elseif (npassing_prev .GT. npassing) then
                nz = nz + 1
                irow(nz) = k_prev + npassing_prev + 2
                icol(nz) = k_prev + npassing_prev + 1
                amat_sp(nz) = -delphim1
                if (colltest) then
                    nz_ttmp = nz_ttmp + 1
                    irow_ttmp(nz_ttmp) = irow(nz)
                    icol_ttmp(nz_ttmp) = icol(nz)
                    amat_ttmp(nz_ttmp) = REAL(amat_sp(nz), dp)
                end if
            end if

! collisions:

            if (fact_pos_e(istep) .NE. 0.d0) then

! Lorentz operator:

                if (isw_lor .EQ. 1) then

                    do ipart = 1, npassing + 1
                        do kk = 1, 5
                            do mm = 0, lag
                                nz = nz + 1
                                irow(nz) = k + ipart
                                icol(nz) = k + MAX(0, ipart - 3) + kk + 2*(npassing + 1)*(mm - m)
                                amat_sp(nz) = anumm(m, mm)*rhs_mat_lorentz(kk, ipart, istep) &
                                              *fact_pos_e(istep)*0.5d0
                                nz_coll = nz_coll + 1
                                irow_coll(nz_coll) = irow(nz)
                                icol_coll(nz_coll) = icol(nz)
                                amat_coll(nz_coll) = anumm(m, mm)*rhs_mat_lorentz(kk, ipart, istep)

                                if (.NOT. colltest .AND. mm .EQ. m) then
                                    nz_ttmp = nz_ttmp + 1
                                    irow_ttmp(nz_ttmp) = irow(nz)
                                    icol_ttmp(nz_ttmp) = icol(nz)
                                    amat_ttmp(nz_ttmp) = -ttmp_mat(kk, ipart, istep)*0.5d0
                                    if (irow(nz) .EQ. icol(nz)) then
                                        amat_ttmp(nz_ttmp) = amat_ttmp(nz_ttmp) + 1.d0
                                    end if
                                end if

                                if (ipart .LE. npassing_prev + 1) then
                                    nz = nz + 1
                                    irow(nz) = k + ipart
                                    icol(nz) = k_prev + MAX(0, ipart - 3) + kk + 2*(npassing_prev + 1)*(mm - m)
                                    amat_sp(nz) = anumm(m, mm)*rhs_mat_lorentz(kk, ipart, istep - 1) &
                                                  *fact_pos_e(istep)*0.5d0
                                    if (.NOT. colltest .AND. mm .EQ. m) then
                                        nz_ttmp = nz_ttmp + 1
                                        irow_ttmp(nz_ttmp) = irow(nz)
                                        icol_ttmp(nz_ttmp) = icol(nz)
                                        amat_ttmp(nz_ttmp) = -ttmp_mat(kk, ipart, istep - 1)*0.5d0
                                    end if
                                elseif (ipart .gt. npassing_prev + 1 .and. addboucol) then
                                    nz = nz + 1
                                    irow(nz) = k + ipart
                                    icol(nz) = k + ipart + 4 - kk + 2*(npassing + 1)*(mm - m)
                                    amat_sp(nz) = anumm(m, mm)*rhs_mat_lorentz(kk, ipart, istep) &
                                                  *fact_pos_e(istep)*0.5d0
                                end if
                            end do
                        end do
                    end do

                end if

                nz_beg = nz + 1
                nz_coll_beg = nz_coll + 1

! energy diffusion operator:

                if (isw_ene .EQ. 1) then

                    do ipart = 1, npassing + 1
                        do kk = 1, 4
                            do mm = 0, lag
                                nz = nz + 1
                                irow(nz) = k + ipart
                                icol(nz) = k + MAX(0, ipart - 2) + kk + 2*(npassing + 1)*(mm - m)
                                amat_sp(nz) = denmm(m, mm)*rhs_mat_energ(kk, ipart, istep)*0.5d0
                                nz_coll = nz_coll + 1
                                irow_coll(nz_coll) = irow(nz)
                                icol_coll(nz_coll) = icol(nz)
                                ! fixed warning: Possible change of value in conversion
                                ! from COMPLEX(8) to REAL(8)
                                amat_coll(nz_coll) = REAL(amat_sp(nz), dp)*2.d0
                            end do
                        end do

                        if (ipart .LE. npassing_prev + 1) then
                            do kk = 1, 4
                                do mm = 0, lag
                                    nz = nz + 1
                                    irow(nz) = k + ipart
                                    icol(nz) = k_prev + MAX(0, ipart - 2) + kk + 2*(npassing_prev + 1)*(mm - m)
                                    amat_sp(nz) = denmm(m, mm)*rhs_mat_energ(kk, ipart, istep - 1)*0.5d0
                                end do
                            end do
                        elseif (ipart .gt. npassing_prev + 1 .and. addboucol) then
                            do kk = 1, 4
                                do mm = 0, lag
                                    nz = nz + 1
                                    irow(nz) = k + ipart
                                    icol(nz) = k + ipart + 3 - kk + 2*(npassing + 1)*(mm - m)
                                    amat_sp(nz) = denmm(m, mm)*rhs_mat_energ(kk, ipart, istep)*0.5d0
                                end do
                            end do
                        end if

                    end do

                end if

                amat_sp(nz_beg:nz) = fact_pos_e(istep)*amat_sp(nz_beg:nz)

            end if

        end do

    end do

! Counter-passing: sigma=-1

    istep = iend
    npassing = npl(istep)

! entry:

    do m = 0, lag
        k = ind_start(istep) + 2*(npassing + 1)*m + 2*npassing + 3

        do ipart = 1, npassing + 1
            nz = nz + 1
            irow(nz) = k - ipart
            icol(nz) = k - ipart
            amat_sp(nz) = (1.d0, 0.d0)
        end do

        if (isw_axisymm .EQ. 1) then
! periodicity:
            k_prev = ind_start(ibeg) + 2*(npassing + 1)*m + 2*npassing + 3

            do ipart = 1, npassing + 1
                nz = nz + 1
                irow(nz) = k - ipart
                icol(nz) = k_prev - ipart
                amat_sp(nz) = (-1.d0, 0.d0)
            end do

!! Modifications by Andreas F. Martitsch (12.12.2016)
! -> Removal of null-space of axisymmetric equation with
! full linearized collision operator (wrong)
!      IF(isw_regper.EQ.1.AND.m.LE.1) then
! -> Removal of null-space of axisymmetric equation without
! integral part of the collision operator
! (solution to the problem without iterations --> correct)
            if (isw_regper .EQ. 1 .AND. m .LT. 1) then
                do ipart = 1, 1
                    if (ipart .LE. npassing) then
                        deleta_factor = (eta(ipart) - eta(ipart - 1))*bhat_mfl(istep)
                    else
                        deleta_factor = 1.d0 - eta(ipart - 1)*bhat_mfl(istep)
                    end if

                    do ipart1 = 1, npassing + 1
                        nz_regper = nz_regper + 1
                        irow_regper(nz_regper) = k - ipart
                        icol_regper(nz_regper) = k_prev - ipart1
                        amat_regper(nz_regper) = deleta_factor
                        nz_regper = nz_regper + 1
                        irow_regper(nz_regper) = k - ipart
                        icol_regper(nz_regper) = k - npassing - 1 - ipart1
                        amat_regper(nz_regper) = deleta_factor
                    end do

                end do
            end if
!! End Modifications by Andreas F. Martitsch (12.12.2016)

        end if

    end do

    do istep = ibeg, iend - 1
        npassing_prev = npl(istep + 1)
        npassing = npl(istep)
        delphim1 = 1.d0/delt_neg(istep)
        deloneovb = 0.5d0*delphim1*(1.d0/bhat_mfl(istep + 1) - 1.d0/bhat_mfl(istep))

        do m = 0, lag
            k_prev = ind_start(istep + 1) + 2*(npassing_prev + 1)*m + 2*npassing_prev + 3
            k = ind_start(istep) + 2*(npassing + 1)*m + 2*npassing + 3

! free flight:

            do ipart = 1, npassing + 1
                nz = nz + 1
                irow(nz) = k - ipart
                icol(nz) = k - ipart
                amat_sp(nz) = delphim1
                if (colltest) then
                    nz_ttmp = nz_ttmp + 1
                    irow_ttmp(nz_ttmp) = irow(nz)
                    icol_ttmp(nz_ttmp) = icol(nz)
                    amat_ttmp(nz_ttmp) = REAL(amat_sp(nz), dp)
                end if
            end do

            do ipart = 1, npassing
                nz = nz + 1
                irow(nz) = k - ipart
                icol(nz) = k_prev - ipart
                amat_sp(nz) = -delphim1
                if (colltest) then
                    nz_ttmp = nz_ttmp + 1
                    irow_ttmp(nz_ttmp) = irow(nz)
                    icol_ttmp(nz_ttmp) = icol(nz)
                    amat_ttmp(nz_ttmp) = REAL(amat_sp(nz), dp)
                end if
            end do

            if (npassing_prev .GE. npassing) then
                nz = nz + 1
                irow(nz) = k - npassing - 1
                icol(nz) = k_prev - npassing - 1
                amat_sp(nz) = -delphim1
                if (colltest) then
                    nz_ttmp = nz_ttmp + 1
                    irow_ttmp(nz_ttmp) = irow(nz)
                    icol_ttmp(nz_ttmp) = icol(nz)
                    amat_ttmp(nz_ttmp) = REAL(amat_sp(nz), dp)
                end if
            end if

! mirroring:

            if (npassing_prev .EQ. npassing) then

                do kk = 1, 4
                    nz = nz + 1
                    irow(nz) = k - npassing - 1
                    icol(nz) = k - npassing - kk + 1
                    amat_sp(nz) = deloneovb*rhs_mat_fzero(kk, istep, 0)
                    if (colltest) then
                        nz_ttmp = nz_ttmp + 1
                        irow_ttmp(nz_ttmp) = irow(nz)
                        icol_ttmp(nz_ttmp) = icol(nz)
                        amat_ttmp(nz_ttmp) = REAL(amat_sp(nz), dp)
                    end if
                end do

                do kk = 1, 4
                    nz = nz + 1
                    irow(nz) = k - npassing - 1
                    icol(nz) = k_prev - npassing_prev - kk + 1
                    amat_sp(nz) = deloneovb*rhs_mat_fzero(kk, istep + 1, 0)
                    if (colltest) then
                        nz_ttmp = nz_ttmp + 1
                        irow_ttmp(nz_ttmp) = irow(nz)
                        icol_ttmp(nz_ttmp) = icol(nz)
                        amat_ttmp(nz_ttmp) = REAL(amat_sp(nz), dp)
                    end if
                end do

            elseif (npassing_prev .GT. npassing) then
                nz = nz + 1
                irow(nz) = k_prev - npassing_prev - 2
                icol(nz) = k_prev - npassing_prev - 1
                amat_sp(nz) = -delphim1
                if (colltest) then
                    nz_ttmp = nz_ttmp + 1
                    irow_ttmp(nz_ttmp) = irow(nz)
                    icol_ttmp(nz_ttmp) = icol(nz)
                    amat_ttmp(nz_ttmp) = REAL(amat_sp(nz), dp)
                end if
            end if

! collisions:

            if (fact_neg_e(istep) .NE. 0.d0) then
! Lorentz operator:

                if (isw_lor .EQ. 1) then

                    do ipart = 1, npassing + 1
                        do kk = 1, 5
                            do mm = 0, lag
                                nz = nz + 1
                                irow(nz) = k - ipart
                                icol(nz) = k - MAX(0, ipart - 3) - kk + 2*(npassing + 1)*(mm - m)
                                amat_sp(nz) = anumm(m, mm)*rhs_mat_lorentz(kk, ipart, istep) &
                                              *fact_neg_e(istep)*0.5d0
                                nz_coll = nz_coll + 1
                                irow_coll(nz_coll) = irow(nz)
                                icol_coll(nz_coll) = icol(nz)
                                amat_coll(nz_coll) = anumm(m, mm)*rhs_mat_lorentz(kk, ipart, istep)

                                if (.NOT. colltest .AND. mm .EQ. m) then
                                    nz_ttmp = nz_ttmp + 1
                                    irow_ttmp(nz_ttmp) = irow(nz)
                                    icol_ttmp(nz_ttmp) = icol(nz)
                                    amat_ttmp(nz_ttmp) = ttmp_mat(kk, ipart, istep)*0.5d0
                                    if (irow(nz) .EQ. icol(nz)) then
                                        amat_ttmp(nz_ttmp) = amat_ttmp(nz_ttmp) - 1.d0
                                    end if
                                end if

                                if (ipart .LE. npassing_prev + 1) then
                                    nz = nz + 1
                                    irow(nz) = k - ipart
                                    icol(nz) = k_prev - MAX(0, ipart - 3) - kk + 2*(npassing_prev + 1)*(mm - m)
                                    amat_sp(nz) = anumm(m, mm)*rhs_mat_lorentz(kk, ipart, istep + 1) &
                                                  *fact_neg_e(istep)*0.5d0
                                    if (.NOT. colltest .AND. mm .EQ. m) then
                                        nz_ttmp = nz_ttmp + 1
                                        irow_ttmp(nz_ttmp) = irow(nz)
                                        icol_ttmp(nz_ttmp) = icol(nz)
                                        amat_ttmp(nz_ttmp) = ttmp_mat(kk, ipart, istep + 1)*0.5d0
                                    end if
                                elseif (ipart .gt. npassing_prev + 1 .and. addboucol) then
                                    nz = nz + 1
                                    irow(nz) = k - ipart
                                    icol(nz) = k - ipart - 4 + kk + 2*(npassing + 1)*(mm - m)
                                    amat_sp(nz) = anumm(m, mm)*rhs_mat_lorentz(kk, ipart, istep) &
                                                  *fact_neg_e(istep)*0.5d0
                                end if

                            end do
                        end do
                    end do

                end if

                nz_beg = nz + 1
                nz_coll_beg = nz_coll + 1

! energy diffusion operator:

                if (isw_ene .EQ. 1) then

                    do ipart = 1, npassing + 1
                        do kk = 1, 4
                            do mm = 0, lag
                                nz = nz + 1
                                irow(nz) = k - ipart
                                icol(nz) = k - MAX(0, ipart - 2) - kk + 2*(npassing + 1)*(mm - m)
                                amat_sp(nz) = denmm(m, mm)*rhs_mat_energ(kk, ipart, istep)*0.5d0
                                nz_coll = nz_coll + 1
                                irow_coll(nz_coll) = irow(nz)
                                icol_coll(nz_coll) = icol(nz)
                                ! fixed warning: Possible change of value in conversion
                                ! from COMPLEX(8) to REAL(8)
                                amat_coll(nz_coll) = REAL(amat_sp(nz), dp)*2.d0
                            end do
                        end do

                        if (ipart .LE. npassing_prev + 1) then
                            do kk = 1, 4
                                do mm = 0, lag
                                    nz = nz + 1
                                    irow(nz) = k - ipart
                                    icol(nz) = k_prev - MAX(0, ipart - 2) - kk + 2*(npassing_prev + 1)*(mm - m)
                                    amat_sp(nz) = denmm(m, mm)*rhs_mat_energ(kk, ipart, istep + 1)*0.5d0
                                end do
                            end do
                        elseif (ipart .gt. npassing_prev + 1 .and. addboucol) then
                            do kk = 1, 4
                                do mm = 0, lag
                                    nz = nz + 1
                                    irow(nz) = k - ipart
                                    icol(nz) = k - ipart - 3 + kk + 2*(npassing + 1)*(mm - m)
                                    amat_sp(nz) = denmm(m, mm)*rhs_mat_energ(kk, ipart, istep)*0.5d0
                                end do
                            end do
                        end if

                    end do

                end if

                amat_sp(nz_beg:nz) = fact_neg_e(istep)*amat_sp(nz_beg:nz)

            end if

        end do

    end do

    if (.false.) then
        ! test particle conservation by the equation matrix (stops the execution after plotting).
        call test_conservation(nz, irow, icol, amat_sp)
    end if

! Save the symmetric matrix:

    call remap_rc(nz, nz_symm, irow, icol, amat_sp)

    allocate (irow_symm(nz_symm), icol_symm(nz_symm), amat_symm(nz_symm))
    irow_symm = irow(1:nz_symm)
    icol_symm = icol(1:nz_symm)
    ! fixed warning: Possible change of value in conversion
    ! from COMPLEX(8) to REAL(8)
    amat_symm = real(amat_sp(1:nz_symm), dp)

! End save symmetric matrix

    deallocate (irow, icol, amat_sp)

    !------------------------------------------------------------------------
    ! Solve axi-symmetric equation set
    !------------------------------------------------------------------------
! For the computation of hatOmegaE from the profile
! one must know the solution of the axisymmetric
! equation set.
    allocate (ipcol(ncol), bvec_sp(ncol))

! Solve the axisymmetric equation set:

    if (allocated(qflux_symm)) deallocate (qflux_symm)
    !  multi-species part (allocate storage for source_vector)
    if (allocated(source_vector_all)) deallocate (source_vector_all)
    allocate (source_vector_all(n_2d_size, 1:4, 0:num_spec - 1))
    source_vector_all = (0.0d0, 0.0d0)

  !! Modification by Andreas F. Martitsch (23.08.2015)
    ! NEO-2 can treat now multiple species -> qflux is now a 4D array
    ! (at the moment these arrays cannot be handled correctly using the
    ! propagator structure -> global variables used):
    if (allocated(qflux_symm_allspec)) deallocate (qflux_symm_allspec)
    allocate (qflux_symm_allspec(1:3, 1:3, 0:num_spec - 1, 0:num_spec - 1))
    qflux_symm_allspec = 0.0d0
  !! End Modification by Andreas F. Martitsch (23.08.2015)
    if (nobounceaver) then

        nz = nz_symm + nz_regper

        allocate (irow(nz), icol(nz), amat_sp(nz))

        irow(1:nz_symm) = irow_symm
        icol(1:nz_symm) = icol_symm
        amat_sp(1:nz_symm) = amat_symm
        if (nz_regper .GT. 0) then
            irow(nz_symm + 1:nz) = irow_regper
            icol(nz_symm + 1:nz) = icol_regper
            amat_sp(nz_symm + 1:nz) = amat_regper
        end if

    !! Modifications by Andreas F. Martitsch (28.08.2014)
        !> geodesic curvature for the axisymmetric field computed by
        !> external routines
        !> geodcu_forw=geodcu_mfl
        !> geodcu_back=geodcu_mfl
        !> computation of geodesic curvature according to
        !> \f[ \|{\nabla}s\| k_{G0} = - \frac{B_\phi}{\iota B_{\theta}+B_{\phi}}
        !> \frac{B_{\rm ref}}{\psi_{\rm tor}^{a}} \frac{\partial B_0}{\partial \theta} \f]
        denomjac = -scalefac_kG*bcovar_phi_hat/(aiota*bcovar_theta_hat + bcovar_phi_hat)
        !> geodcu_forw used for computation of q_rip(1:npassing+1,istep,1),
        !> which in turn enters the source_vector

        if (mag_magfield .ne. 3) then
            geodcu_forw = denomjac*dlogbdphi_mfl*bhat_mfl
        else
            ! Overwrite geodcu_forw in the case of EFIT input
            geodcu_forw = geodcu_mfl
        end if

        ! geodcu_back used for the computation of convol_flux, which enters q_flux
        ! via flux_vector(1,:) and flux_vector(3,:)
        !--> Computation of D31/D32 not affected ( flux_vector(2,:) determined by convol_curr )
        geodcu_back = geodcu_forw
    !! End Modifications by Andreas F. Martitsch (28.08.2014)

        call source_flux
    !! Modification by Andreas F. Martitsch (23.08.2015)
        ! save solution of the differential part for species=ispec
        ! (diffusion coeff. driven by thermodyn. forces of other
        ! species are zero -> interaction through integral part)
        source_vector_all(:, 1:4, ispec) = source_vector(:, 1:4)
    !! End Modification by Andreas F. Martitsch (23.08.2015)

        problem_type = .TRUE.
        call solve_eqs(.TRUE.)

        ! Debugging - plot distribution function (axisymmetric problem)
        if (lsw_debug_distfun) then
            do ispecp = 0, num_spec - 1
                uw = 10000*(num_spec*ispec + ispecp + 1)
                istep = (ibeg + iend)/2
                uw_new = uw
                call plotsource(uw_new, REAL(source_vector_all(:, :, ispec)))
                uw_new = uw + 1000
                call plotsource(uw_new, aimag(source_vector_all(:, :, ispec)))
                istep = ibeg
                uw_new = uw + 10
                call plotsource(uw_new, REAL(source_vector_all(:, :, ispec)))
                uw_new = uw + 1010
                call plotsource(uw_new, aimag(source_vector_all(:, :, ispec)))
                istep = iend
                uw_new = uw + 20
                call plotsource(uw_new, REAL(source_vector_all(:, :, ispec)))
                uw_new = uw + 1020
                call plotsource(uw_new, aimag(source_vector_all(:, :, ispec)))
                istep = ibeg + 1
                uw_new = uw + 30
                call plotsource(uw_new, REAL(source_vector_all(:, :, ispec)))
                uw_new = uw + 1030
                call plotsource(uw_new, aimag(source_vector_all(:, :, ispec)))
            end do
        end if

        deallocate (irow, icol, amat_sp)

! Caution!!! factor 2 is not needed!!!
        qflux = 2.d0*qflux
        if (.NOT. lsw_multispecies) then ! single-species output
            open (1234, file='qflux_symm.dat')
            write (1234, *) conl_over_mfp/boozer_iota
            write (1234, *) - qflux(1, 1:3)
            write (1234, *) - qflux(2, 1:3)
            write (1234, *) - qflux(3, 1:3)
            close (1234)
        end if

    !! Modifications by Andreas F. Martitsch (28.07.2014)
        ! Save here the solution of the axisymmetric case (ifnon-axisymmetric
        ! solution is also computed). Otherwise exit here ripple_solver and
        ! return to calling routine propagator_solver.
        if (isw_qflux_NA .EQ. 0) then
            ! compute qflux only for the axisymmetric equation set
       !! Modification by Andreas F. Martitsch (23.08.2015)
            ! old behavior:
            ! --> qflux already available from prop_a%p%qflux
            !RETURN
            ! NEO-2 can treat now multiple species -> qflux is now a 4D array
            ! (at the moment these arrays cannot be handled correctly using the
            ! propagator structure -> global variables used):
            if (.not. allocated(qflux_allspec)) stop "ERROR: Axisymm. solution does not exist!"
            qflux_allspec = 2.0d0*qflux_allspec ! Caution!!! factor 2 is not needed!!!
            qflux_symm_allspec = qflux_allspec
            if (allocated(qflux_allspec)) deallocate (qflux_allspec)
            call save_qflux_symm_allspec()
            RETURN
       !! End Modification by Andreas F. Martitsch (23.08.2015)
        else if (isw_qflux_NA .EQ. 1) then
            ! save qflux for the axisymmetric equation set
            ! and proceed with the solution of the non-axisymmetric
            ! equation set (stored within prop_a%p%qflux)
            allocate (qflux_symm(3, 3))
            qflux_symm = qflux
       !! Modification by Andreas F. Martitsch (23.08.2015)
            ! NEO-2 can treat now multiple species -> qflux is now a 4D array
            ! (at the moment these arrays cannot be handled correctly using the
            ! propagator structure -> global variables used):
            if (.not. allocated(qflux_allspec)) stop "ERROR: Axisymm. solution does not exist!"
            qflux_allspec = 2.0d0*qflux_allspec ! Caution!!! factor 2 is not needed!!!
            qflux_symm_allspec = qflux_allspec
            if (allocated(qflux_allspec)) deallocate (qflux_allspec)
            call save_qflux_symm_allspec()
       !! End Modification by Andreas F. Martitsch (23.08.2015)
        else
            stop "ERROR: Invalid input for isw_qflux_symm (0/1)!"
        end if
    !! End Modifications by Andreas F. Martitsch (28.07.2014)

    end if

    if (lsw_multispecies .AND. isw_calc_Er .EQ. 1) then
        print *, 'Compute radial electric field ...'
        if (num_spec .EQ. 1) then
            call get_Er(qflux_symm_allspec, Er)
            print *, 'Er: ', Er
        else
            call get_Er(qflux_symm_allspec, Er, avEparB_ov_avb2)
            print *, 'Er, avEparB_ov_avb2: ', Er, avEparB_ov_avb2
        end if
        MtOvR = MtOvR_spec(ispec)
    end if
    if (lsw_multispecies .AND. isw_calc_MagDrift .EQ. 1) then
        print *, 'Compute hatOmegaB ...'
        call get_B_rho_L_loc()
        B_rho_L_loc = B_rho_L_loc_spec(ispec)
    end if

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
    hatOmegaE = MtOvR/collpar
    ! normalized magnetic rotation frequency ($\hat{\Omega}_{tB}^{\rm ref}$)
    ! specified via Larmor radius associated with $B_{00}^{Booz}$ (rho_L_loc)
    !definition with conl_over_mfp suffers again from wrong conversion to
    !collpar
    !hatOmegaB=((device%r0*PI)/(2.0d0*boozer_psi_pr_hat))*&
    !     (B_rho_L_loc/(avbhat*conl_over_mfp))
    !correct definition with collpar=4/$l_c$ ($l_c=2 v_{Ta} \tau_{aa}$)

    hatOmegaB = B_rho_L_loc/(2.0d0*boozer_psi_pr_hat*(bmod0*1.0d4)*collpar)

    hatOmegaB = hatOmegaB*sign(1.d0, boozer_psi_pr_hat*bcovar_phi_hat*boozer_isqrg)

    print *, 'hatOmegaB,hatOmegaE: ', hatOmegaB, hatOmegaE

  !! End Modifications by Andreas F. Martitsch (14.07.2015)

    !------------------------------------------------------------------------
    ! End Solve axi-symmetric equation set
    !------------------------------------------------------------------------

! Nox-axisymmetric matrices:

! Periodicity for Co-passing, sigma=1

! Determine the size of arrays:

    nz = nz_symm

    istep = ibeg
    npassing = npl(istep)

    do m = 0, lag
        k = ind_start(istep) + 2*(npassing + 1)*m
        k_prev = ind_start(iend) + 2*(npassing + 1)*m

        do ipart = 1, npassing + 1
            nz = nz + 1
!      irow_per_pos(nz)=k+ipart
!      icol_per_pos(nz)=k_prev+ipart
        end do

    end do

    nz_per_pos = nz

    allocate (irow_per_pos(nz_symm + 1:nz_per_pos))
    allocate (icol_per_pos(nz_symm + 1:nz_per_pos))

! Fill the arrays:

    nz = nz_symm

    istep = ibeg
    npassing = npl(istep)

    do m = 0, lag
        k = ind_start(istep) + 2*(npassing + 1)*m
        k_prev = ind_start(iend) + 2*(npassing + 1)*m

        do ipart = 1, npassing + 1
            nz = nz + 1
            irow_per_pos(nz) = k + ipart
            icol_per_pos(nz) = k_prev + ipart
        end do

    end do

! End periodicity for Co-passing, sigma=1

! Periodicity for Counter-passing, sigma=-1

! Determine the size of arrays:

    nz = nz_per_pos

    istep = iend
    npassing = npl(istep)

    do m = 0, lag
        k = ind_start(istep) + 2*(npassing + 1)*m + 2*npassing + 3
        k_prev = ind_start(ibeg) + 2*(npassing + 1)*m + 2*npassing + 3

        do ipart = 1, npassing + 1
            nz = nz + 1
!      irow_per_neg(nz)=k-ipart
!      icol_per_neg(nz)=k_prev-ipart
        end do

    end do

    nz_per_neg = nz

    allocate (irow_per_neg(nz_per_pos + 1:nz_per_neg))
    allocate (icol_per_neg(nz_per_pos + 1:nz_per_neg))

! Fill the arrays:

    nz = nz_per_pos

    istep = iend
    npassing = npl(istep)

    do m = 0, lag
        k = ind_start(istep) + 2*(npassing + 1)*m + 2*npassing + 3
        k_prev = ind_start(ibeg) + 2*(npassing + 1)*m + 2*npassing + 3

        do ipart = 1, npassing + 1
            nz = nz + 1
            irow_per_neg(nz) = k - ipart
            icol_per_neg(nz) = k_prev - ipart
        end do

    end do

! End periodicity for Counter-passing, sigma=-1

! Rotation matrix:

! Determine the size of arrays:

    nz = nz_per_neg

! Co-passing: sigma=1

    do istep = ibeg + 1, iend
        npassing_prev = npl(istep - 1)
        npassing = npl(istep)

        do m = 0, lag
            k_prev = ind_start(istep - 1) + 2*(npassing_prev + 1)*m
            k = ind_start(istep) + 2*(npassing + 1)*m

            if (fact_pos_e(istep) .NE. 0.d0) then
!        nz_beg=nz+1

! Toroidal rotation:

                do ipart = 1, npassing
                    do kk = 1, 4
                        do mm = 0, lag
                            nz = nz + 1
!              irow_asymm(nz)=k+ipart
!              icol_asymm(nz)=k+max(0,ipart-2)+kk+2*(npassing+1)*(mm-m)
!              amat_asymm(nz)=x1mm(m,mm)*rhs_mat_energ(kk,ipart,istep)
                        end do
                    end do
                end do

                do kk = 1, 4
                    do mm = 0, lag
                        nz = nz + 1
!            irow_asymm(nz)=k+npassing+1
!            icol_asymm(nz)=k+npassing-1+kk+2*(npassing+1)*(mm-m)
!            amat_asymm(nz)=0.5d0*x1mm(m,mm)                                  &
!                       *rhs_mat_energ(kk,npassing+1,istep)
                    end do
                end do

                do kk = 1, 4
                    do mm = 0, lag
                        nz = nz + 1
!            irow_asymm(nz)=k+npassing+1
!            icol_asymm(nz)=k_prev+npassing_prev-1+kk+2*(npassing_prev+1)*(mm-m)
!            amat_asymm(nz)=0.5d0*x1mm(m,mm)                                  &
!                       *rhs_mat_energ(kk,npassing_prev+1,istep-1)
                    end do
                end do

!        amat_asymm(nz_beg:nz)=-fact_pos_e(istep)*amat_asymm(nz_beg:nz)

            end if

        end do

    end do

! Counter-passing: sigma=-1

    do istep = ibeg, iend - 1
        npassing_prev = npl(istep + 1)
        npassing = npl(istep)

        do m = 0, lag
            k_prev = ind_start(istep + 1) + 2*(npassing_prev + 1)*m + 2*npassing_prev + 3
            k = ind_start(istep) + 2*(npassing + 1)*m + 2*npassing + 3

            if (fact_neg_e(istep) .NE. 0.d0) then
!        nz_beg=nz+1

! Toroidal rotation:

                do ipart = 1, npassing
                    do kk = 1, 4
                        do mm = 0, lag
                            nz = nz + 1
!              irow_asymm(nz)=k-ipart
!              icol_asymm(nz)=k-max(0,ipart-2)-kk+2*(npassing+1)*(mm-m)
!              amat_asymm(nz)=x1mm(m,mm)*rhs_mat_energ(kk,ipart,istep)
                        end do
                    end do
                end do

                do kk = 1, 4
                    do mm = 0, lag
                        nz = nz + 1
!            irow_asymm(nz)=k-npassing-1
!            icol_asymm(nz)=k-npassing+1-kk+2*(npassing+1)*(mm-m)
!            amat_asymm(nz)=0.5d0*x1mm(m,mm)                                  &
!                       *rhs_mat_energ(kk,npassing+1,istep)
                    end do
                end do

                do kk = 1, 4
                    do mm = 0, lag
                        nz = nz + 1
!            irow_asymm(nz)=k-npassing-1
!            icol_asymm(nz)=k_prev-npassing_prev+1-kk+2*(npassing_prev+1)*(mm-m)
!            amat_asymm(nz)=0.5d0*x1mm(m,mm)                                  &
!                       *rhs_mat_energ(kk,npassing_prev+1,istep+1)
                    end do
                end do

!        amat_asymm(nz_beg:nz)=-fact_neg_e(istep)*amat_asymm(nz_beg:nz)

            end if

        end do

    end do

    nz_asymm = nz

    allocate (irow_asymm(nz_per_neg + 1:nz_asymm))
    allocate (icol_asymm(nz_per_neg + 1:nz_asymm))
    allocate (amat_asymm(nz_per_neg + 1:nz_asymm))

! Fill the arrays:
!
! Notes on normalization:
! matrix rhs_mat_energ corresponds to discretization over $\eta$
! of the following function, $-\kappa f/(2 h^\varphi |\lambda|)$
! matrix rhs_mat_energ2 corresponds to discretization over $\eta$
! of the following function, $-\kappa f |\lambda|/(2 h^\varphi)$
! where $\kappa$=collpar

  !! Modifications by Andreas F. Martitsch (01.04.2015)
    ! Definition of $\kappa$ specified in the comment above originates
    ! from the ntv_booz (Version Nov 2013), which differs by a factor 2
    ! from the quantity collpar that is used here. Therefore, rhs_mat_energ
    ! and rhs_mat_energ2 are corrected by a factor two (no 2 in the denominator)
  !! End Modifications by Andreas F. Martitsch (01.04.2015)

    denomjac = aiota*bcovar_theta_hat + bcovar_phi_hat

    nz = nz_per_neg

! Co-passing: sigma=1

    do istep = ibeg + 1, iend

    !! Modifications by Andreas F. Martitsch (14.03.2014)
        ! Optional output (necessary for modeling the magnetic rotation)
        a1b = (bcovar_s_hat_mfl(istep)*dlogbdphi_mfl(istep)/denomjac &
               - dlogbds_mfl(istep))/aiota
    !! Modifications by Andreas F. Martitsch (17.03.2016)
        ! derivative of iota for non-local NTV computations
        ! (with magnetic shear)
        !-> old:
        !a2b=a1b+2.d0*(dbcovar_theta_hat_ds+dbcovar_phi_hat_ds/aiota                  &
        !     -          dbcovar_s_hat_dphi_mfl(istep)/aiota)/denomjac
        !-> new [include radial derivative of iota if
        !-> isw_mag_shear .eq. 0; otherwise set to zero]
        a2b = a1b + 2.d0*(dbcovar_theta_hat_ds + dbcovar_phi_hat_ds/aiota &
                          - bcovar_phi_hat*boozer_iota_s/(aiota**2) &
                          - dbcovar_s_hat_dphi_mfl(istep)/aiota)/denomjac
    !! End Modifications by Andreas F. Martitsch (17.03.2016)
    !! End Modifications by Andreas F. Martitsch (14.03.2014)

        npassing_prev = npl(istep - 1)
        npassing = npl(istep)

        do m = 0, lag
            k_prev = ind_start(istep - 1) + 2*(npassing_prev + 1)*m
            k = ind_start(istep) + 2*(npassing + 1)*m

            if (fact_pos_e(istep) .NE. 0.d0) then
                nz_beg = nz + 1

! Toroidal rotation:

                do ipart = 1, npassing
                    do kk = 1, 4
                        do mm = 0, lag
                            nz = nz + 1
                            irow_asymm(nz) = k + ipart
                            icol_asymm(nz) = k + MAX(0, ipart - 2) + kk + 2*(npassing + 1)*(mm - m)
              !! Modifications by Andreas F. Martitsch (01.04.2015)
                            ! rhs_mat_energ and rhs_mat_energ2 are corrected by a factor two
                            ! (see comments related to "DEFINITION OF INPUT QUANTITIES" and
                            ! "Fill the arrays")
                            amat_asymm(nz) = (hatOmegaE*x1mm(m, mm) + hatOmegaB*a1b*x2mm(m, mm)) &
                                             *(rhs_mat_energ(kk, ipart, istep)*2.0d0) &
                                             + hatOmegaB*a2b*x2mm(m, mm) &
                                             *(rhs_mat_energ2(kk, ipart, istep)*2.0d0)
              !! End Modifications by Andreas F. Martitsch (01.04.2015)
                        end do
                    end do
                end do

                do kk = 1, 4
                    do mm = 0, lag
                        nz = nz + 1
                        irow_asymm(nz) = k + npassing + 1
                        icol_asymm(nz) = k + npassing - 1 + kk + 2*(npassing + 1)*(mm - m)
            !! Modifications by Andreas F. Martitsch (01.04.2015)
                        ! rhs_mat_energ and rhs_mat_energ2 are corrected by a factor two
                        ! (see comments related to "DEFINITION OF INPUT QUANTITIES" and
                        ! "Fill the arrays")
                        amat_asymm(nz) = 0.5d0*( &
                                         (hatOmegaE*x1mm(m, mm) + hatOmegaB*a1b*x2mm(m, mm)) &
                                         *(rhs_mat_energ(kk, npassing + 1, istep)*2.0d0) &
                                         + hatOmegaB*a2b*x2mm(m, mm) &
                                         *(rhs_mat_energ2(kk, npassing + 1, istep)*2.0d0) &
                                         )
            !! End Modifications by Andreas F. Martitsch (01.04.2015)
                    end do
                end do

                do kk = 1, 4
                    do mm = 0, lag
                        nz = nz + 1
                        irow_asymm(nz) = k + npassing + 1
                        icol_asymm(nz) = k_prev + npassing_prev - 1 + kk + 2*(npassing_prev + 1)*(mm - m)
            !! Modifications by Andreas F. Martitsch (01.04.2015)
                        ! rhs_mat_energ and rhs_mat_energ2 are corrected by a factor two
                        ! (see comments related to "DEFINITION OF INPUT QUANTITIES" and
                        ! "Fill the arrays")
                        amat_asymm(nz) = 0.5d0*( &
                                         (hatOmegaE*x1mm(m, mm) + hatOmegaB*a1b*x2mm(m, mm)) &
                                         *(rhs_mat_energ(kk, npassing_prev + 1, istep - 1)*2.0d0) &
                                         + hatOmegaB*a2b*x2mm(m, mm) &
                                         *(rhs_mat_energ2(kk, npassing_prev + 1, istep - 1)*2.0d0) &
                                         )
            !! End Modifications by Andreas F. Martitsch (01.04.2015)
                    end do
                end do

                amat_asymm(nz_beg:nz) = -fact_pos_e(istep)*amat_asymm(nz_beg:nz)

            end if

        end do

    end do

! Counter-passing: sigma=-1

    do istep = ibeg, iend - 1

     !! Modifications by Andreas F. Martitsch (14.03.2014)
        ! Optional output (necessary for modeling the magnetic rotation)
        a1b = (bcovar_s_hat_mfl(istep)*dlogbdphi_mfl(istep)/denomjac &
               - dlogbds_mfl(istep))/aiota
     !! Modifications by Andreas F. Martitsch (17.03.2016)
        ! derivative of iota for non-local NTV computations
        ! (with magnetic shear)
        !-> old:
        !a2b=a1b+2.d0*(dbcovar_theta_hat_ds+dbcovar_phi_hat_ds/aiota                  &
        !     -          dbcovar_s_hat_dphi_mfl(istep)/aiota)/denomjac
        !-> new [include radial derivative of iota if
        !-> isw_mag_shear .eq. 0; otherwise set to zero]
        a2b = a1b + 2.d0*(dbcovar_theta_hat_ds + dbcovar_phi_hat_ds/aiota &
                          - bcovar_phi_hat*boozer_iota_s/(aiota**2) &
                          - dbcovar_s_hat_dphi_mfl(istep)/aiota)/denomjac
     !! End Modifications by Andreas F. Martitsch (17.03.2016)
     !! End Modifications by Andreas F. Martitsch (14.03.2014)

        npassing_prev = npl(istep + 1)
        npassing = npl(istep)

        do m = 0, lag
            k_prev = ind_start(istep + 1) + 2*(npassing_prev + 1)*m + 2*npassing_prev + 3
            k = ind_start(istep) + 2*(npassing + 1)*m + 2*npassing + 3

            if (fact_neg_e(istep) .NE. 0.d0) then
                nz_beg = nz + 1

! Toroidal rotation:

                do ipart = 1, npassing
                    do kk = 1, 4
                        do mm = 0, lag
                            nz = nz + 1
                            irow_asymm(nz) = k - ipart
                            icol_asymm(nz) = k - MAX(0, ipart - 2) - kk + 2*(npassing + 1)*(mm - m)
              !! Modifications by Andreas F. Martitsch (01.04.2015)
                            ! rhs_mat_energ and rhs_mat_energ2 are corrected by a factor two
                            ! (see comments related to "DEFINITION OF INPUT QUANTITIES" and
                            ! "Fill the arrays")
                            amat_asymm(nz) = (hatOmegaE*x1mm(m, mm) + hatOmegaB*a1b*x2mm(m, mm)) &
                                             *(rhs_mat_energ(kk, ipart, istep)*2.0d0) &
                                             + hatOmegaB*a2b*x2mm(m, mm) &
                                             *(rhs_mat_energ2(kk, ipart, istep)*2.0d0)
              !! End Modifications by Andreas F. Martitsch (01.04.2015)
                        end do
                    end do
                end do

                do kk = 1, 4
                    do mm = 0, lag
                        nz = nz + 1
                        irow_asymm(nz) = k - npassing - 1
                        icol_asymm(nz) = k - npassing + 1 - kk + 2*(npassing + 1)*(mm - m)
            !! Modifications by Andreas F. Martitsch (01.04.2015)
                        ! rhs_mat_energ and rhs_mat_energ2 are corrected by a factor two
                        ! (see comments related to "DEFINITION OF INPUT QUANTITIES" and
                        ! "Fill the arrays")
                        amat_asymm(nz) = 0.5d0*( &
                                         (hatOmegaE*x1mm(m, mm) + hatOmegaB*a1b*x2mm(m, mm)) &
                                         *(rhs_mat_energ(kk, npassing + 1, istep)*2.0d0) &
                                         + hatOmegaB*a2b*x2mm(m, mm) &
                                         *(rhs_mat_energ2(kk, npassing + 1, istep)*2.0d0) &
                                         )
            !! End Modifications by Andreas F. Martitsch (01.04.2015)
                    end do
                end do
!
                do kk = 1, 4
                    do mm = 0, lag
                        nz = nz + 1
                        irow_asymm(nz) = k - npassing - 1
                        icol_asymm(nz) = k_prev - npassing_prev + 1 - kk + 2*(npassing_prev + 1)*(mm - m)
            !! Modifications by Andreas F. Martitsch (01.04.2015)
                        ! rhs_mat_energ and rhs_mat_energ2 are corrected by a factor two
                        ! (see comments related to "DEFINITION OF INPUT QUANTITIES" and
                        ! "Fill the arrays")
                        amat_asymm(nz) = 0.5d0*( &
                                         (hatOmegaE*x1mm(m, mm) + hatOmegaB*a1b*x2mm(m, mm)) &
                                         *(rhs_mat_energ(kk, npassing_prev + 1, istep + 1)*2.0d0) &
                                         + hatOmegaB*a2b*x2mm(m, mm) &
                                         *(rhs_mat_energ2(kk, npassing_prev + 1, istep + 1)*2.0d0) &
                                         )
            !! End Modifications by Andreas F. Martitsch (01.04.2015)
                    end do
                end do

                amat_asymm(nz_beg:nz) = -fact_neg_e(istep)*amat_asymm(nz_beg:nz)

            end if

        end do

    end do

! Use solution of axisymmetric equation set:

    if (nobounceaver) then

        allocate (f0_coll(n_2d_size, 3), f0_ttmp(n_2d_size, 3))
        allocate (f0_coll_all(n_2d_size, 3, 0:num_spec - 1), f0_ttmp_all(n_2d_size, 3, 0:num_spec - 1))

        do ispecp = 0, num_spec - 1
            f0_coll = 0.d0
            f0_ttmp = 0.d0

! Here f0_coll=$-\frac{1}{h^\varphi}\sum_{m^\prime}\hat L_{mm^\prime}^c
!                bar f_{m^\prime}^{\sigma (k)}$ :

            do nz = 1, nz_coll
                ! fixed warning: Possible change of value in conversion
                ! from COMPLEX(8) to REAL(8)
                f0_coll(irow_coll(nz), :) = f0_coll(irow_coll(nz), :) + amat_coll(nz) &
                                            *REAL(source_vector_all(icol_coll(nz), 1:3, ispecp), dp)
            end do

            if (isw_intp .EQ. 1) then
                allocate (bvec_iter(ncol), bvec_prev(ncol))

                do i = 1, 3
                    bvec_prev = source_vector_all(:, i, ispecp)

                    call integral_part(bvec_prev, bvec_iter)

                    ! fixed warning: Possible change of value in conversion
                    ! from COMPLEX(8) to REAL(8)
                    f0_coll(:, i) = f0_coll(:, i) - REAL(bvec_iter, dp)
                end do

                deallocate (bvec_iter, bvec_prev)
            end if

! Here f0_ttmp=$-\sigma \eta \difp{}{\eta} f_{m^\prime}^{\sigma (k)}$ :

            do nz = 1, nz_ttmp
                ! fixed warning: Possible change of value in conversion
                ! from COMPLEX(8) to REAL(8)
                f0_ttmp(irow_ttmp(nz), :) = f0_ttmp(irow_ttmp(nz), :) + amat_ttmp(nz) &
                                            *REAL(source_vector_all(icol_ttmp(nz), 1:3, ispecp), dp)
            end do

            f0_ttmp_all(:, :, ispecp) = f0_ttmp(:, :)
            f0_coll_all(:, :, ispecp) = f0_coll(:, :)

        end do

        if (colltest) then

            call source_flux

            istep = (ibeg + iend)/3
            npassing = npl(istep)

            do m = 0, lag
                k = ind_start(istep) + 2*(npassing + 1)*m
                do i = k + 1, k + 2*(npassing + 1)
                    write (7000 + m, *) sngl(f0_coll(i, :) + f0_ttmp(i, :))
                    write (8000 + m, *) sngl(REAL(source_vector(i, 1:3)))
                end do
            end do
        end if

        if (ttmptest) then

! Plot the mirroring force $-\lambda \eta \difp{}{\eta} f_{m^\prime}^{\sigma (k)}$
! as function of $\lambda$ :

            istep = (ibeg + iend)/3

            call plotsource(9000, f0_ttmp)

        end if

    end if

! End Use solution axisymmetric equation set

! Solve the non-axisymmetric equation set:

  !! Modifications by Andreas F. Martitsch (13.06.2014)
    ! quantities of the perturbation field are extracted
    ! from the Boozer file
    ! (-> no further need to specify them manually +
    ! tests using an artificial perturbation field can
    ! be done inside tmp section (see above))
    !bnoverb0=bnoverb0*EXP(imun*m_phi*phi_mfl)
    !dbnoverb0_dtheta=dbnoverb0_dtheta*EXP(imun*m_phi*phi_mfl)
  !! End Modifications by Andreas F. Martitsch (13.06.2014)
    denomjac = aiota*bcovar_theta_hat + bcovar_phi_hat

    geodcu_back = scalefac_kG*imun*m_phi*bnoverb0*bhat_mfl
    geodcu_forw = geodcu_back

    ! Ware pinch effect
    do istep = ibeg, iend
        q_rip(:, istep, 2) = -2.d0*q_rip(:, istep, 2)*bnoverb0(istep)
    end do

    if (nobounceaver) then
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
        geodcu_forw = geodcu_forw - scalefac_kG*(bcovar_phi_hat/denomjac) &
                      *(dbnoverb0_dphi_mfl - bnoverb0*dlogbdphi_mfl)*bhat_mfl
    !! End Modifications by Andreas F. Martitsch (28.08.2014)
    end if

    do ispecp = 0, num_spec - 1
        if (nobounceaver) then
            f0_ttmp(:, :) = f0_ttmp_all(:, :, ispecp)
            f0_coll(:, :) = f0_coll_all(:, :, ispecp)
        end if

        if (ispecp .EQ. ispec) then
            call source_flux
        else
            source_vector = 0.0d0
        end if

        if (nobounceaver) then
            allocate (ttmpfact(ibeg:iend))
        !! Modifications by Andreas F. Martitsch (13.06.2014)
            ! derivative along the periodic Boozer angle theta has
            ! been redefined to a derivative along the field line (phi_mfl)
            ! (changes concerning dbnoverb0_dtheta are required!)
            !ttmpfact=aiota*dbnoverb0_dtheta+imun*m_phi*bnoverb0
            ttmpfact = dbnoverb0_dphi_mfl
        !! End Modifications by Andreas F. Martitsch (13.06.2014)

            call add_f01_source

            deallocate (ttmpfact)
        end if
        source_vector_all(:, :, ispecp) = source_vector(:, 1:4)
    end do

    expforw = EXP(imun*m_phi*(phi_mfl(iend) - phi_mfl(ibeg)))
    expbackw = (1.d0, 0.d0)/expforw
    perbou_pos = (1.d0, 0.d0) - expbackw
    perbou_neg = (1.d0, 0.d0) - expforw

    nz = nz_asymm

    allocate (irow(nz), icol(nz), amat_sp(nz))

    irow(1:nz_symm) = irow_symm
    icol(1:nz_symm) = icol_symm
    amat_sp(1:nz_symm) = amat_symm

    irow(nz_symm + 1:nz_per_pos) = irow_per_pos(nz_symm + 1:nz_per_pos)
    icol(nz_symm + 1:nz_per_pos) = icol_per_pos(nz_symm + 1:nz_per_pos)
    amat_sp(nz_symm + 1:nz_per_pos) = perbou_pos

    irow(nz_per_pos + 1:nz_per_neg) = irow_per_neg(nz_per_pos + 1:nz_per_neg)
    icol(nz_per_pos + 1:nz_per_neg) = icol_per_neg(nz_per_pos + 1:nz_per_neg)
    amat_sp(nz_per_pos + 1:nz_per_neg) = perbou_neg

    irow(nz_per_neg + 1:nz_asymm) = irow_asymm(nz_per_neg + 1:nz_asymm)
    icol(nz_per_neg + 1:nz_asymm) = icol_asymm(nz_per_neg + 1:nz_asymm)
    amat_sp(nz_per_neg + 1:nz_asymm) = amat_asymm(nz_per_neg + 1:nz_asymm)*rotfactor

    problem_type = .FALSE.
    call solve_eqs(.TRUE.)
  !! Modification by Andreas F. Martitsch (23.08.2015)
    ! NEO-2 can treat now multiple species -> qflux is now a 4D array
    ! (at the moment these arrays cannot be handled correctly using the
    ! propagator structure -> global variables used):
    if (allocated(qflux_ntv_allspec)) deallocate (qflux_ntv_allspec)
    allocate (qflux_ntv_allspec(1:3, 1:3, 0:num_spec - 1, 0:num_spec - 1))
    if (.NOT. allocated(qflux_allspec)) STOP "Non-Axisymm. solution does not exist!"
    qflux_ntv_allspec = qflux_allspec
  !! End Modification by Andreas F. Martitsch (23.08.2015)
  !! Modification by Andreas F. Martitsch (23.08.2015)
    !  multi-species part (ifclean is true, deallocate memory)
    if (allocated(source_vector_all)) deallocate (source_vector_all)
    if (allocated(qflux_allspec)) deallocate (qflux_allspec)
  !! End Modification by Andreas F. Martitsch (23.08.2015)

    if (.NOT. lsw_multispecies) then ! single-species output
        open (1234, file='qflux_ntv.dat')
        write (1234, *) conl_over_mfp/boozer_iota
        write (1234, *) - qflux(1, 1:3)
        write (1234, *) - qflux(2, 1:3)
        write (1234, *) - qflux(3, 1:3)
        close (1234)
    end if

    if (lsw_multispecies .AND. mpro%getrank() .EQ. 0) then ! multi-species output
        open (070915, file='qflux_ntv_allspec.dat')
        write (070915, *) '% boozer_s, collpar'
        write (070915, *) boozer_s, collpar
        write (070915, *) '% qflux(1,1,a,b}'
        do ispecp = 0, num_spec - 1
            write (070915, *) (qflux_ntv_allspec(1, 1, ispecp, ispecpp), &
                               ispecpp=0, num_spec - 1)
        end do
        write (070915, *) '% qflux(1,3,a,b}'
        do ispecp = 0, num_spec - 1
            write (070915, *) (qflux_ntv_allspec(1, 3, ispecp, ispecpp), &
                               ispecpp=0, num_spec - 1)
        end do
        close (070915)

    end if

    iopt = 3

    call sparse_solve(nrow, ncol, nz, irow(1:nz), ipcol, amat_sp(1:nz), bvec_sp, iopt)

    deallocate (flux_vector, source_vector, irow, icol, amat_sp, ipcol, bvec_sp, bvec_parflow)
    if (allocated(flux_vector_plot)) deallocate (flux_vector_plot)

    deallocate (energvec_ket, energvec_bra)
    deallocate (densvec_ket, densvec_bra)

    call CPU_TIME(time_solver)
    print *, 'solving completed       ', time_solver - time_factorization, ' sec'

    deallocate (deriv_coef, enu_coef, alambd, Vg_vp_over_B, scalprod_pleg)
    deallocate (alampow, vrecurr, dellampow, convol_polpow, pleg_bra, pleg_ket)
    deallocate (enu_coef2, dellampow2, rhs_mat_energ2)
    deallocate (npl, rhs_mat_fzero, rhs_mat_lorentz, rhs_mat_energ, q_rip, q_rip_1)
    deallocate (q_rip_incompress, q_rip_parflow)
    deallocate (fun_coef, ttmp_mat)
    deallocate (convol_flux, convol_curr, ind_start, convol_flux_0)
    deallocate (delt_pos, delt_neg, fact_pos_b, fact_neg_b, fact_pos_e, fact_neg_e)

    !------------------------------------------------------------------------
    ! END SERGEI
    !------------------------------------------------------------------------

    !------------------------------------------------------------------------
    ! final deallocation of initial things provided by winny
    if (allocated(phi_mfl)) deallocate (phi_mfl)
    if (allocated(bhat_mfl)) deallocate (bhat_mfl)
    if (allocated(geodcu_mfl)) deallocate (geodcu_mfl)
    if (allocated(h_phi_mfl)) deallocate (h_phi_mfl)
    if (allocated(dlogbdphi_mfl)) deallocate (dlogbdphi_mfl)
    ! Optional output (necessary for modeling the magnetic rotation)
    ! --> deallocate the arrays
    if (allocated(dlogbds_mfl)) deallocate (dlogbds_mfl)
    if (allocated(bcovar_s_hat_mfl)) deallocate (bcovar_s_hat_mfl)
    if (allocated(dbcovar_s_hat_dphi_mfl)) deallocate (dbcovar_s_hat_dphi_mfl)

    if (allocated(eta)) deallocate (eta)
    if (allocated(eta_prev)) deallocate (eta_prev)
    if (allocated(eta_next)) deallocate (eta_next)

    !------------------------------------------------------------------------

CONTAINS

!------------------------------------------------------------------------
    subroutine plotsource(iunit_base, sourcevec_tmp)

        integer, intent(in) :: iunit_base
        double precision, dimension(n_2d_size, 3), intent(in) :: sourcevec_tmp

        double precision :: alambd_save1

        npassing = npl(istep)
        delta_eta = eta(1:npassing) - eta(0:npassing - 1)
        eta0 = 1.d0/bhat_mfl(istep)
        do m = 0, lag
            k = ind_start(istep) + 2*(npassing + 1)*m
            do i = 1, npassing + 1
                if (i .LE. npassing) then
                    alambd_save1 = 0.5d0*(alambd(i, istep) + alambd(i - 1, istep))
                    write (iunit_base + m, *) - alambd_save1, alambd_save1 &
                        *sngl(sourcevec_tmp(k + 2*(npassing + 1) - i + 1, :))/delta_eta(i)
                else
                    alambd_save1 = 0.5d0*alambd(i - 1, istep)
                    write (iunit_base + m, *) - alambd_save1, alambd_save1 &
                        *sngl(sourcevec_tmp(k + 2*(npassing + 1) - i + 1, :))/(eta0 - eta(i - 1))
                end if
            end do
            do i = npassing + 1, 1, -1
                if (i .LE. npassing) then
                    alambd_save1 = 0.5d0*(alambd(i, istep) + alambd(i - 1, istep))
                    write (iunit_base + m, *) alambd_save1, alambd_save1 &
                        *sngl(sourcevec_tmp(k + i, :))/delta_eta(i)
                else
                    alambd_save1 = 0.5d0*alambd(i - 1, istep)
                    write (iunit_base + m, *) alambd_save1, alambd_save1 &
                        *sngl(sourcevec_tmp(k + i, :))/(eta0 - eta(i - 1))
                end if
            end do

            flush (iunit_base + m)

        end do

    end subroutine plotsource

    subroutine fluxdenplot(sourcevec_tmp, asymplot, prefix, ind_spec)

        implicit none

        complex(dp), dimension(n_2d_size), intent(in) :: sourcevec_tmp
        logical, intent(in) :: asymplot
        character(len=2), intent(in) :: prefix
        integer, intent(in) :: ind_spec

        integer, parameter :: nx = 100

        integer :: iunit_base, nmax, m, i, k
        complex(dp), dimension(3) :: ab, abexp
        real(kind=kind(1d0)), dimension(:, :), allocatable :: phi_mat, alam_mat, fun_mat
        real(kind=kind(1d0)), dimension(:), allocatable :: exp2inp
        real(kind=kind(1d0)), dimension(:, :), allocatable :: dummy2d
        real(kind=kind(1d0)) :: hx
        logical :: dir_exisits

        inquire (file=char(48 + ispec), exist=dir_exisits)
        if (.not. dir_exisits) return

        iunit_base = 12345

        if (.not. asymplot) then
            open (iunit_base, file=char(48 + ispec)//'/fluxden'//trim(prefix) &
              &                  //char(48 + ind_spec)//'.dat')
            do istep = ibeg, iend
                k = ind_start(istep)
                npassing = npl(istep)
                write (iunit_base, *) aiota*phi_mfl(istep), &
                    dble(matmul(conjg(flux_vector_plot(1:3, k + 1:k + 2*(npassing + 1)*(lag + 1))), &
                                sourcevec_tmp(k + 1:k + 2*(npassing + 1)*(lag + 1))))
            end do
            close (iunit_base)

        else
            allocate (exp2inp(0:nx))
            hx = 8.d0*atan(1.d0)/real(nx, kind=kind(1d0))
            open (iunit_base, file=char(48 + ispec)//'/phi1d.dat')
            do i = 0, nx
                exp2inp(i) = exp(cmplx(0.d0, 2.d0*hx*real(i*m_phi, kind=kind(1d0))))
                write (iunit_base, *) hx*real(i, kind=kind(1d0))
            end do
            close (iunit_base)
            open (iunit_base, file=char(48 + ispec)//'/fluxden'//trim(prefix) &
              &                 //char(48 + ind_spec)//'.dat')
            open (iunit_base + 1, file=char(48 + ispec)//'/theta1d.dat')
            do istep = ibeg, iend
                k = ind_start(istep)
                npassing = npl(istep)
                ab = matmul(conjg(flux_vector_plot(1:3, k + 1:k + 2*(npassing + 1)*(lag + 1))), &
                   & sourcevec_tmp(k + 1:k + 2*(npassing + 1)*(lag + 1)))
                abexp = matmul(flux_vector_plot(1:3, k + 1:k + 2*(npassing + 1)*(lag + 1)), &
                   &  sourcevec_tmp(k + 1:k + 2*(npassing + 1)*(lag + 1)))               &
                   & *exp(cmplx(0.d0, -real(2*m_phi, kind=kind(1d0))*phi_mfl(istep)))
                write (iunit_base, *) real(ab(1) + abexp(1)*exp2inp, kind=kind(1d0))
                write (iunit_base + 1, *) aiota*phi_mfl(istep)
            end do
            close (iunit_base)
            close (iunit_base + 1)
            deallocate (exp2inp)

        end if

    end subroutine fluxdenplot

    subroutine solve_eqs(clean)
! Solve the linear equation set:

        use arnoldi_mod, only: iterator, f_init_arnoldi

        logical :: clean
        real(dp), dimension(:), allocatable :: bvec_sp_real
        !  multi-species part
        integer :: ispecpp ! species indices (loop over sources)
        real(dp), dimension(:, :, :, :), allocatable :: qflux_allspec_tmp

        integer :: i_src ! own source index for plotting the distribution function

        if (isw_intp .EQ. 1) allocate (bvec_iter(ncol), bvec_prev(ncol))

        call remap_rc(nz, nz_sq, irow, icol, amat_sp)

        print *, 'information for species = ', ispec
        print *, 'system size = ', n_2d_size
        print *, 'non-zeros before and after truncation = ', nz, nz_sq
        print *, 'maximum value = ', maxval(abs(amat_sp))
        nz = nz_sq

        call column_full2pointer(icol(1:nz), ipcol)

! There are now three different calls sparse_solve
!   iopt = 1 ; factorization
!   iopt = 2 ; solve
!   iopt = 3 ; free memory
! without iopt or with iopt = 0: old behaviour (factorize,solve,free)

! factorization:

        bvec_sp = 0.d0
        iopt = 1

        call CPU_TIME(time_start)
        if (problem_type) then
            allocate (bvec_sp_real(ncol))
            bvec_sp_real = 0.0d0
            call sparse_solve(nrow, ncol, nz, irow(1:nz), ipcol, DBLE(amat_sp(1:nz)), &
                              bvec_sp_real, iopt)
        else
            call sparse_solve(nrow, ncol, nz, irow(1:nz), ipcol, amat_sp(1:nz), &
                              bvec_sp, iopt)
        end if

        call CPU_TIME(time_factorization)
        print *, 'factorization completed ', time_factorization - time_start, ' sec'

        iopt = 2

! Solution of inhomogeneus equation (account of sources):

        denom_energ = sum(energvec_bra*energvec_ket)
        denom_dens = sum(densvec_bra*densvec_ket)

        if (isw_intp .EQ. 0) then
            ! no integral part:

            mode_iter = 1

            do k = 1, 3
                f_init_arnoldi = source_vector_all(:, 1, ispec)
                call iterator(mode_iter, n_2d_size, n_arnoldi, epserr_iter, niter,&
                            & source_vector_all(:, k, ispec), ispec, next_iteration)
            end do

        elseif (isw_intp .EQ. 1) then
            ! with integral part:

            mode_iter = 1

            do k = 1, 3
                !  multi-species part:
                do ispecp = 0, num_spec - 1
                    f_init_arnoldi = source_vector_all(:, 1, ispec)
                    if (ispec .eq. 0) print *, 'source ', k, ' driving species', ispecp, ':'
                    drive_spec = ispecp
                    call iterator(mode_iter, n_2d_size, n_arnoldi, epserr_iter, niter,&
                                & source_vector_all(:, k, ispecp), ispec, next_iteration)
                end do

            end do

        end if

        if (lsw_write_flux_surface_distribution) then
            call write_flux_surface_distribution(flux_vector_plot, source_vector_all, &
              & ind_start, npl, delt_pos, phi_mfl, iend - ibeg + 1)
        end if

        i_src = 1
        if (problem_type) then
            call matlabplot_allm(real(source_vector_all(:, i_src, ispec)), .true., .true., '  ')
            do ispecp = 0, num_spec - 1
                call fluxdenplot(source_vector_all(:, i_src, ispec), .false., 'ax', ispecp)
            end do
        else
            call matlabplot_allm(real(source_vector_all(:, i_src, ispec)), .true., .false., 're')
            call matlabplot_allm(dimag(source_vector_all(:, i_src, ispec)), .true., .false., 'im')
            do ispecp = 0, num_spec - 1
                call fluxdenplot(source_vector_all(:, i_src, ispec), .true., 'na', ispecp)
            end do
        end if

        if (clean) then
            mode_iter = 3

            call iterator(mode_iter, n_2d_size, n_arnoldi, epserr_iter, niter, &
              & source_vector(:, k), ispec, next_iteration)

        end if

! Plotting:

        if (iplot .EQ. 1) then
            nphiplot = 200    !serves as an upper limit for data output
            delphiplot = MAXVAL(phi_mfl(ibeg + 1:iend) - phi_mfl(ibeg:iend - 1))
            if (delphiplot .GT. EPSILON(1.d0)) then
                ! fixed warning: Possible change of value in conversion
                ! from REAL(8) to integer(4)
                nphiequi = INT((phi_mfl(iend) - phi_mfl(ibeg))/delphiplot)

                nphiequi = MAX(1, nphiequi)
            else
                nphiequi = 1
            end if
            delphiplot = (phi_mfl(iend) - phi_mfl(ibeg))/MIN(nphiplot, nphiequi)
            nplp1 = npart_loc + 1
            allocate (fun_write(0:lag, 0:3, 0:nplp1, 3))
            icounter = 0
            phiplot = phi_mfl(ibeg) - 1.d0

            write (propname, *) fieldpropagator%tag
            open (iunit_phi, file='phi_mesh.' &
                  //TRIM(ADJUSTL(propname))//'.dat')
            open (iunit_dt_p, form='unformatted', file='dentf_p.' &
                  //TRIM(ADJUSTL(propname))//'.dat')
            open (iunit_dt_m, form='unformatted', file='dentf_m.' &
                  //TRIM(ADJUSTL(propname))//'.dat')
            open (iunit_sp_p, form='unformatted', file='spitf_p.' &
                  //TRIM(ADJUSTL(propname))//'.dat')
            open (iunit_sp_m, form='unformatted', file='spitf_m.' &
                  //TRIM(ADJUSTL(propname))//'.dat')
            open (iunit_et_p, form='unformatted', file='enetf_p.' &
                  //TRIM(ADJUSTL(propname))//'.dat')
            open (iunit_et_m, form='unformatted', file='enetf_m.' &
                  //TRIM(ADJUSTL(propname))//'.dat')

            do istep = ibeg, iend
                if (phi_mfl(istep) .LT. phiplot .AND. istep .NE. iend) CYCLE
                icounter = icounter + 1
                phiplot = phi_mfl(istep) + delphiplot
                npassing = npl(istep)
                eta0 = 1.d0/bhat_mfl(istep)
                write (iunit_phi, *) phi_mfl(istep), npassing, bhat_mfl(istep)

                fun_write = 0.d0
                do m = 0, lag
                    k = ind_start(istep) + 2*(npassing + 1)*m
                    ! fixed warning: Possible change of value in conversion
                    ! from COMPLEX(8) to REAL(8)
                    fun_write(m, :, 1, :) = REAL(MATMUL(derivs_plot(:, :, 1, istep), &
                                                        source_vector(k + 1:k + 4, 1:3)), dp)
                    do i = 2, npassing + 1
                        ! fixed warning: Possible change of value in conversion
                        ! from COMPLEX(8) to REAL(8)
                        fun_write(m, :, i, :) = REAL(MATMUL(derivs_plot(:, :, i, istep), &
                                                            source_vector(k + i - 1:k + i + 2, 1:3)), dp)
                    end do
                end do
                write (iunit_dt_p) fun_write(:, :, :, 1)
                write (iunit_sp_p) fun_write(:, :, :, 2)/surface_boozer_B00
                write (iunit_et_p) fun_write(:, :, :, 3)

                fun_write = 0.d0
                do m = 0, lag
                    k = ind_start(istep) + 2*(npassing + 1)*(m + 1)
                    ! fixed warning: Possible change of value in conversion
                    ! from COMPLEX(8) to REAL(8)
                    fun_write(m, :, 1, :) = REAL(MATMUL(derivs_plot(:, :, 1, istep), &
                                                        source_vector(k:k - 3:-1, 1:3)), dp)
                    do i = 2, npassing + 1
                        ! fixed warning: Possible change of value in conversion
                        ! from COMPLEX(8) to REAL(8)
                        fun_write(m, :, i, :) = REAL(MATMUL(derivs_plot(:, :, i, istep), &
                                                            source_vector(k - i + 2:k - i - 1:-1, 1:3)), dp)
                    end do
                end do
                write (iunit_dt_m) fun_write(:, :, :, 1)
                write (iunit_sp_m) fun_write(:, :, :, 2)/surface_boozer_B00
                write (iunit_et_m) fun_write(:, :, :, 3)

            end do

            close (iunit_phi)
            close (iunit_dt_p)
            close (iunit_dt_m)
            close (iunit_sp_p)
            close (iunit_sp_m)
            close (iunit_et_p)
            close (iunit_et_m)
            open (iunit_sizes, file='sizeplot_etalev.' &
                  //TRIM(ADJUSTL(propname))//'.dat')
            write (iunit_sizes, *) lag, nplp1, icounter, collpar, travis_convfac
            write (iunit_sizes, *) eta(0:nplp1)
            close (iunit_sizes)

            deallocate (fun_write)

        end if

    !! Modification by Andreas F. Martitsch (23.08.2015)
        ! old behavior (for a single species)
        !qflux=0.5d0*REAL(MATMUL(CONJG(flux_vector),source_vector(:,1:3)))
        !  multi-species part
        if (allocated(qflux_allspec)) deallocate (qflux_allspec)
        allocate (qflux_allspec(1:3, 1:3, 0:num_spec - 1, 0:num_spec - 1))
        qflux_allspec = 0.0d0
        do ispecp = 0, num_spec - 1
            qflux = 0.5d0*REAL(MATMUL(CONJG(flux_vector), source_vector_all(:, 1:3, ispecp)), dp)
            qflux_allspec(:, :, ispecp, ispec) = qflux
        end do
        ! order of species inidices (ispecp,ispec) interchanged
        ! (-> easier to handle within mpro%allgather)
        call mpro%allgather_inplace(qflux_allspec)
        ! go back to the "natural" order of species indices (ispec,ispecp)
        if (allocated(qflux_allspec_tmp)) deallocate (qflux_allspec_tmp)
        allocate (qflux_allspec_tmp(1:3, 1:3, 0:num_spec - 1, 0:num_spec - 1))
        qflux_allspec_tmp = 0.0d0
        do ispecp = 0, num_spec - 1
            do ispecpp = 0, num_spec - 1
                qflux_allspec_tmp(:, :, ispecp, ispecpp) = qflux_allspec(:, :, ispecpp, ispecp)
            end do
        end do
        qflux_allspec = qflux_allspec_tmp
        if (allocated(qflux_allspec_tmp)) deallocate (qflux_allspec_tmp)
    !! End Modification by Andreas F. Martitsch (23.08.2015)

        if (clean) then

            iopt = 3

            if (problem_type) then
                call sparse_solve(nrow, ncol, nz, irow(1:nz), ipcol, DBLE(amat_sp(1:nz)), &
                                  bvec_sp_real, iopt)
                deallocate (bvec_sp_real)
            else
                call sparse_solve(nrow, ncol, nz, irow(1:nz), ipcol, amat_sp(1:nz), bvec_sp, iopt)
            end if

            if (isw_intp .EQ. 1) then
                deallocate (bvec_iter, bvec_prev)
            end if

        end if

    end subroutine solve_eqs

!------------------------------------------------------------------------
    subroutine source_flux

        do istep = ibeg, iend
            npassing = npl(istep)
            q_rip(1:npassing + 1, istep, 1) = q_rip_1(1:npassing + 1, istep) &
                                              *geodcu_forw(istep)
            convol_flux(1:npassing + 1, istep) = convol_flux_0(1:npassing + 1, istep) &
                                                 *geodcu_back(istep)
        end do

        flux_vector = 0.d0
        source_vector = 0.d0
        bvec_parflow = 0.d0
        energvec_ket = 0.d0
        energvec_bra = 0.d0
        densvec_bra = 0.d0

        call generate_maxwellian(densvec_ket)

        do istep = ibeg, iend

            ioddeven = MOD(istep - ibeg, 2) !0 for full RK step, 1 for half RK step

            if (istep .EQ. ibeg) then
                step_factor_p = delt_pos(ibeg + 1)*fact_pos_b(ibeg)/3.d0
                step_factor_m = delt_neg(ibeg)*fact_neg_e(ibeg)/3.d0
            elseif (istep .EQ. iend) then
                step_factor_p = delt_pos(iend)*fact_pos_e(iend)/3.d0
                step_factor_m = delt_neg(iend - 1)*fact_neg_b(iend)/3.d0
            elseif (ioddeven .EQ. 1) then
                step_factor_p = (delt_pos(istep + 1)*fact_pos_b(istep) &
                                 + delt_pos(istep)*fact_pos_e(istep))/1.5d0
                step_factor_m = (delt_neg(istep - 1)*fact_neg_b(istep) &
                                 + delt_neg(istep)*fact_neg_e(istep))/1.5d0
            else
                step_factor_p = (delt_pos(istep + 1)*fact_pos_b(istep) &
                                 + delt_pos(istep)*fact_pos_e(istep))/3.d0
                step_factor_m = (delt_neg(istep - 1)*fact_neg_b(istep) &
                                 + delt_neg(istep)*fact_neg_e(istep))/3.d0
            end if

            npassing = npl(istep)

            do m = 0, lag
                k = ind_start(istep) + 2*(npassing + 1)*m

                flux_vector(1, k + 1:k + npassing + 1) = &
                    step_factor_p*weightlag(1, m)*convol_flux(1:npassing + 1, istep)
                flux_vector(1, k + npassing + 2:k + 2*npassing + 2) = &
                    step_factor_m*weightlag(1, m)*convol_flux(npassing + 1:1:-1, istep)

                flux_vector(2, k + 1:k + npassing + 1) = &
                    step_factor_p*weightlag(2, m)*convol_curr(1:npassing + 1, istep)
                flux_vector(2, k + npassing + 2:k + 2*npassing + 2) = &
                    -step_factor_m*weightlag(2, m)*convol_curr(npassing + 1:1:-1, istep)

                flux_vector(3, k + 1:k + npassing + 1) = &
                    step_factor_p*weightlag(3, m)*convol_flux(1:npassing + 1, istep)
                flux_vector(3, k + npassing + 2:k + 2*npassing + 2) = &
                    step_factor_m*weightlag(3, m)*convol_flux(npassing + 1:1:-1, istep)

                flux_vector_plot(1, k + 1:k + npassing + 1) = &
                    weightlag(1, m)*convol_flux(1:npassing + 1, istep)
                flux_vector_plot(1, k + npassing + 2:k + 2*npassing + 2) = &
                    weightlag(1, m)*convol_flux(npassing + 1:1:-1, istep)
                flux_vector_plot(2, k + 1:k + npassing + 1) = &
                    weightlag(2, m)*convol_curr(1:npassing + 1, istep)
                flux_vector_plot(2, k + npassing + 2:k + 2*npassing + 2) = &
                    -weightlag(2, m)*convol_curr(npassing + 1:1:-1, istep)
                flux_vector_plot(3, k + 1:k + npassing + 1) = &
                    weightlag(3, m)*convol_flux(1:npassing + 1, istep)
                flux_vector_plot(3, k + npassing + 2:k + 2*npassing + 2) = &
                    weightlag(3, m)*convol_flux(npassing + 1:1:-1, istep)

                ! Use pre-conditioned iterations (not necessary/depricated):
                ! -> remove parallel flow from solution
                ! Computation of bvec_parflow generalized to non-orthogonal polynomials
                !-> old version (orthogonal test functions):
                !bvec_parflow(k+1:k+npassing+1)           = asource(m,1)*q_rip_parflow(1:npassing+1,istep)
                !bvec_parflow(k+npassing+2:k+2*npassing+2)=-asource(m,1)*q_rip_parflow(npassing+1:1:-1,istep)
                !-> new version (general test functions):
                bvec_parflow(k + 1:k + npassing + 1) = weightparflow(m)*q_rip_parflow(1:npassing + 1, istep)
                bvec_parflow(k + npassing + 2:k + 2*npassing + 2) = -weightparflow(m)*q_rip_parflow(npassing + 1:1:-1, istep)
                ! End Use pre-conditioned iterations (not necessary/depricated)

                ! Use pre-conditioned iterations:
                ! -> remove null-space of axisymmetric solution (energy conservation)
                energvec_ket(k + 1:k + npassing) = &
                    weightenerg(m)*(eta(1:npassing) - eta(0:npassing - 1))
                energvec_ket(k + 2*npassing + 2:k + npassing + 3:-1) = &
                    weightenerg(m)*(eta(1:npassing) - eta(0:npassing - 1))

                energvec_ket(k + npassing + 1) = &
                    weightenerg(m)*((1.d0/bhat_mfl(istep)) - eta(npassing))
                energvec_ket(k + npassing + 2) = &
                    weightenerg(m)*((1.d0/bhat_mfl(istep)) - eta(npassing))

                energvec_bra(k + 1:k + npassing + 1) = &
                    step_factor_p*(weightlag(1, m) - 1.5d0*weightden(m))*pleg_bra(0, 1:npassing + 1, istep)
                densvec_bra(k + 1:k + npassing + 1) =                                &
                  & step_factor_p*weightden(m)*pleg_bra(0, 1:npassing + 1, istep)
                energvec_bra(k + npassing + 2:k + 2*npassing + 2) = &
                    step_factor_m*(weightlag(1, m) - 1.5d0*weightden(m))*pleg_bra(0, npassing + 1:1:-1, istep)
                densvec_bra(k + npassing + 2:k + 2*npassing + 2) =                     &
                  & step_factor_m*weightden(m)*pleg_bra(0, npassing + 1:1:-1, istep)

                energvec_bra(k + 1:k + 2*npassing + 2) = &
                    energvec_bra(k + 1:k + 2*npassing + 2)/(bhat_mfl(istep))
                densvec_bra(k + 1:k + 2*npassing + 2) =                              &
                  & densvec_bra(k + 1:k + 2*npassing + 2)/(bhat_mfl(istep))
                ! End Use pre-conditioned iterations

                if (istep .GT. ibeg) then
                    npassing_prev = npl(istep - 1)
                    k_prev = ind_start(istep - 1) + 2*(npassing_prev + 1)*m

                    if (ioddeven .EQ. 1) then
                        npassing_next = npl(istep + 1)
                        source_vector(k + 1:k + npassing + 1, 1) &
                            = source_vector(k + 1:k + npassing + 1, 1) &
                              + asource(m, 1)/1.5d0*q_rip(1:npassing + 1, istep, 1) &
                              *fact_pos_e(istep)
                        source_vector(k + 1:k + npassing + 1, 2) &
                            = source_vector(k + 1:k + npassing + 1, 2) &
                              + asource(m, 2)/1.5d0*q_rip(1:npassing + 1, istep, 2) &
                              *fact_pos_e(istep)
                        source_vector(k + 1:k + npassing + 1, 3) &
                            = source_vector(k + 1:k + npassing + 1, 3) &
                              + asource(m, 3)/1.5d0*q_rip(1:npassing + 1, istep, 1) &
                              *fact_pos_e(istep)
                        source_vector(k + 1:k + npassing + 1, 4) &
                            = source_vector(k + 1:k + npassing + 1, 4) &
                              + asource(m, 1)/1.5d0*q_rip_incompress(1:npassing + 1, istep) &
                              *fact_pos_e(istep)

                        source_vector(k + 1:k + npassing_prev + 1, 1) &
                            = source_vector(k + 1:k + npassing_prev + 1, 1) &
                              + asource(m, 1)/2.4d0*q_rip(1:npassing_prev + 1, istep - 1, 1) &
                              *fact_pos_b(istep - 1)
                        source_vector(k + 1:k + npassing_prev + 1, 2) &
                            = source_vector(k + 1:k + npassing_prev + 1, 2) &
                              + asource(m, 2)/2.4d0*q_rip(1:npassing_prev + 1, istep - 1, 2) &
                              *fact_pos_b(istep - 1)
                        source_vector(k + 1:k + npassing_prev + 1, 3) &
                            = source_vector(k + 1:k + npassing_prev + 1, 3) &
                              + asource(m, 3)/2.4d0*q_rip(1:npassing_prev + 1, istep - 1, 1) &
                              *fact_pos_b(istep - 1)
                        source_vector(k + 1:k + npassing_prev + 1, 4) &
                            = source_vector(k + 1:k + npassing_prev + 1, 4) &
                              + asource(m, 1)/2.4d0*q_rip_incompress(1:npassing_prev + 1, istep - 1) &
                              *fact_pos_b(istep - 1)

                        source_vector(k + 1:k + npassing_next + 1, 1) &
                            = source_vector(k + 1:k + npassing_next + 1, 1) &
                              - asource(m, 1)/12d0*q_rip(1:npassing_next + 1, istep + 1, 1) &
                              *fact_pos_e(istep + 1)
                        source_vector(k + 1:k + npassing_next + 1, 2) &
                            = source_vector(k + 1:k + npassing_next + 1, 2) &
                              - asource(m, 2)/12d0*q_rip(1:npassing_next + 1, istep + 1, 2) &
                              *fact_pos_e(istep + 1)
                        source_vector(k + 1:k + npassing_next + 1, 3) &
                            = source_vector(k + 1:k + npassing_next + 1, 3) &
                              - asource(m, 3)/12d0*q_rip(1:npassing_next + 1, istep + 1, 1) &
                              *fact_pos_e(istep + 1)
                        source_vector(k + 1:k + npassing_next + 1, 4) &
                            = source_vector(k + 1:k + npassing_next + 1, 4) &
                              - asource(m, 1)/12d0*q_rip_incompress(1:npassing_next + 1, istep + 1) &
                              *fact_pos_e(istep + 1)
                    else
                        npassing_next = npl(istep - 2)
                        source_vector(k + 1:k + npassing + 1, 1) &
                            = source_vector(k + 1:k + npassing + 1, 1) &
                              + asource(m, 1)/2.4d0*q_rip(1:npassing + 1, istep, 1) &
                              *fact_pos_e(istep)
                        source_vector(k + 1:k + npassing + 1, 2) &
                            = source_vector(k + 1:k + npassing + 1, 2) &
                              + asource(m, 2)/2.4d0*q_rip(1:npassing + 1, istep, 2) &
                              *fact_pos_e(istep)
                        source_vector(k + 1:k + npassing + 1, 3) &
                            = source_vector(k + 1:k + npassing + 1, 3) &
                              + asource(m, 3)/2.4d0*q_rip(1:npassing + 1, istep, 1) &
                              *fact_pos_e(istep)
                        source_vector(k + 1:k + npassing + 1, 4) &
                            = source_vector(k + 1:k + npassing + 1, 4) &
                              + asource(m, 1)/2.4d0*q_rip_incompress(1:npassing + 1, istep) &
                              *fact_pos_e(istep)
                        if (npassing_prev .LE. npassing) then
                            source_vector(k + 1:k + npassing_prev + 1, 1) &
                                = source_vector(k + 1:k + npassing_prev + 1, 1) &
                                  + asource(m, 1)/1.5d0*q_rip(1:npassing_prev + 1, istep - 1, 1) &
                                  *fact_pos_b(istep - 1)
                            source_vector(k + 1:k + npassing_prev + 1, 2) &
                                = source_vector(k + 1:k + npassing_prev + 1, 2) &
                                  + asource(m, 2)/1.5d0*q_rip(1:npassing_prev + 1, istep - 1, 2) &
                                  *fact_pos_b(istep - 1)
                            source_vector(k + 1:k + npassing_prev + 1, 3) &
                                = source_vector(k + 1:k + npassing_prev + 1, 3) &
                                  + asource(m, 3)/1.5d0*q_rip(1:npassing_prev + 1, istep - 1, 1) &
                                  *fact_pos_b(istep - 1)
                            source_vector(k + 1:k + npassing_prev + 1, 4) &
                                = source_vector(k + 1:k + npassing_prev + 1, 4) &
                                  + asource(m, 1)/1.5d0*q_rip_incompress(1:npassing_prev + 1, istep - 1) &
                                  *fact_pos_b(istep - 1)
                        else
                            source_vector(k + 1:k + npassing + 1, 1) &
                                = source_vector(k + 1:k + npassing + 1, 1) &
                                  + asource(m, 1)/1.5d0*q_rip(1:npassing + 1, istep - 1, 1) &
                                  *fact_pos_b(istep - 1)
                            source_vector(k_prev + npassing_prev + 2, 1) &
                                = source_vector(k_prev + npassing_prev + 2, 1) &
                                  + asource(m, 1)/1.5d0*q_rip(npassing_prev + 1, istep - 1, 1) &
                                  *fact_pos_b(istep - 1)
                            source_vector(k + 1:k + npassing + 1, 2) &
                                = source_vector(k + 1:k + npassing + 1, 2) &
                                  + asource(m, 2)/1.5d0*q_rip(1:npassing + 1, istep - 1, 2) &
                                  *fact_pos_b(istep - 1)
                            source_vector(k_prev + npassing_prev + 2, 2) &
                                = source_vector(k_prev + npassing_prev + 2, 2) &
                                  + asource(m, 2)/1.5d0*q_rip(npassing_prev + 1, istep - 1, 2) &
                                  *fact_pos_b(istep - 1)
                            source_vector(k + 1:k + npassing + 1, 3) &
                                = source_vector(k + 1:k + npassing + 1, 3) &
                                  + asource(m, 3)/1.5d0*q_rip(1:npassing + 1, istep - 1, 1) &
                                  *fact_pos_b(istep - 1)
                            source_vector(k_prev + npassing_prev + 2, 3) &
                                = source_vector(k_prev + npassing_prev + 2, 3) &
                                  + asource(m, 3)/1.5d0*q_rip(npassing_prev + 1, istep - 1, 1) &
                                  *fact_pos_b(istep - 1)
                            source_vector(k + 1:k + npassing + 1, 4) &
                                = source_vector(k + 1:k + npassing + 1, 4) &
                                  + asource(m, 1)/1.5d0*q_rip_incompress(1:npassing + 1, istep - 1) &
                                  *fact_pos_b(istep - 1)
                            source_vector(k_prev + npassing_prev + 2, 4) &
                                = source_vector(k_prev + npassing_prev + 2, 4) &
                                  + asource(m, 1)/1.5d0*q_rip_incompress(npassing_prev + 1, istep - 1) &
                                  *fact_pos_b(istep - 1)
                        end if
                        if (npassing_next .LE. npassing) then
                            source_vector(k + 1:k + npassing_next + 1, 1) &
                                = source_vector(k + 1:k + npassing_next + 1, 1) &
                                  - asource(m, 1)/12d0*q_rip(1:npassing_next + 1, istep - 2, 1) &
                                  *fact_pos_b(istep - 2)
                            source_vector(k + 1:k + npassing_next + 1, 2) &
                                = source_vector(k + 1:k + npassing_next + 1, 2) &
                                  - asource(m, 2)/12d0*q_rip(1:npassing_next + 1, istep - 2, 2) &
                                  *fact_pos_b(istep - 2)
                            source_vector(k + 1:k + npassing_next + 1, 3) &
                                = source_vector(k + 1:k + npassing_next + 1, 3) &
                                  - asource(m, 3)/12d0*q_rip(1:npassing_next + 1, istep - 2, 1) &
                                  *fact_pos_b(istep - 2)
                            source_vector(k + 1:k + npassing_next + 1, 4) &
                                = source_vector(k + 1:k + npassing_next + 1, 4) &
                                  - asource(m, 1)/12d0*q_rip_incompress(1:npassing_next + 1, istep - 2) &
                                  *fact_pos_b(istep - 2)
                        else
                            source_vector(k + 1:k + npassing + 1, 1) &
                                = source_vector(k + 1:k + npassing + 1, 1) &
                                  - asource(m, 1)/12d0*q_rip(1:npassing + 1, istep - 2, 1) &
                                  *fact_pos_b(istep - 2)
                            source_vector(k_prev + npassing_prev + 2, 1) &
                                = source_vector(k_prev + npassing_prev + 2, 1) &
                                  - asource(m, 1)/12d0*q_rip(npassing_next + 1, istep - 2, 1) &
                                  *fact_pos_b(istep - 2)
                            source_vector(k + 1:k + npassing + 1, 2) &
                                = source_vector(k + 1:k + npassing + 1, 2) &
                                  - asource(m, 2)/12d0*q_rip(1:npassing + 1, istep - 2, 2) &
                                  *fact_pos_b(istep - 2)
                            source_vector(k_prev + npassing_prev + 2, 2) &
                                = source_vector(k_prev + npassing_prev + 2, 2) &
                                  - asource(m, 2)/12d0*q_rip(npassing_next + 1, istep - 2, 2) &
                                  *fact_pos_b(istep - 2)
                            source_vector(k + 1:k + npassing + 1, 3) &
                                = source_vector(k + 1:k + npassing + 1, 3) &
                                  - asource(m, 3)/12d0*q_rip(1:npassing + 1, istep - 2, 1) &
                                  *fact_pos_b(istep - 2)
                            source_vector(k_prev + npassing_prev + 2, 3) &
                                = source_vector(k_prev + npassing_prev + 2, 3) &
                                  - asource(m, 3)/12d0*q_rip(npassing_next + 1, istep - 2, 1) &
                                  *fact_pos_b(istep - 2)
                            source_vector(k + 1:k + npassing + 1, 4) &
                                = source_vector(k + 1:k + npassing + 1, 4) &
                                  - asource(m, 1)/12d0*q_rip_incompress(1:npassing + 1, istep - 2) &
                                  *fact_pos_b(istep - 2)
                            source_vector(k_prev + npassing_prev + 2, 4) &
                                = source_vector(k_prev + npassing_prev + 2, 4) &
                                  - asource(m, 1)/12d0*q_rip_incompress(npassing_next + 1, istep - 2) &
                                  *fact_pos_b(istep - 2)
                        end if
                    end if
                elseif (iplot .EQ. 1 .AND. isw_axisymm .NE. 1) then
                    source_vector(k + 1:k + npassing + 1, 1:3) &
                        = flux_pl(npass_l*m + 1:npass_l*(m + 1), :)
                end if

                if (istep .LT. iend) then
                    npassing_prev = npl(istep + 1)
                    k_prev = ind_start(istep + 1) + 2*(npassing_prev + 1)*m
                    if (ioddeven .EQ. 1) then
                        npassing_next = npl(istep - 1)
                        source_vector(k + npassing + 2:k + 2*npassing + 2, 1) &
                            = source_vector(k + npassing + 2:k + 2*npassing + 2, 1) &
                              + asource(m, 1)/1.5d0*q_rip(npassing + 1:1:-1, istep, 1) &
                              *fact_neg_e(istep)
                        source_vector(k + npassing + 2:k + 2*npassing + 2, 2) &
                            = source_vector(k + npassing + 2:k + 2*npassing + 2, 2) &
                              - asource(m, 2)/1.5d0*q_rip(npassing + 1:1:-1, istep, 2) &
                              *fact_neg_e(istep)
                        source_vector(k + npassing + 2:k + 2*npassing + 2, 3) &
                            = source_vector(k + npassing + 2:k + 2*npassing + 2, 3) &
                              + asource(m, 3)/1.5d0*q_rip(npassing + 1:1:-1, istep, 1) &
                              *fact_neg_e(istep)
                        source_vector(k + npassing + 2:k + 2*npassing + 2, 4) &
                            = source_vector(k + npassing + 2:k + 2*npassing + 2, 4) &
                              + asource(m, 1)/1.5d0*q_rip_incompress(npassing + 1:1:-1, istep) &
                              *fact_neg_e(istep)

                        source_vector(k + 2*npassing + 2 - npassing_prev:k + 2*npassing + 2, 1) &
                            = source_vector(k + 2*npassing + 2 - npassing_prev:k + 2*npassing + 2, 1) &
                              + asource(m, 1)/2.4d0*q_rip(npassing_prev + 1:1:-1, istep + 1, 1) &
                              *fact_neg_b(istep + 1)
                        source_vector(k + 2*npassing + 2 - npassing_prev:k + 2*npassing + 2, 2) &
                            = source_vector(k + 2*npassing + 2 - npassing_prev:k + 2*npassing + 2, 2) &
                              - asource(m, 2)/2.4d0*q_rip(npassing_prev + 1:1:-1, istep + 1, 2) &
                              *fact_neg_b(istep + 1)
                        source_vector(k + 2*npassing + 2 - npassing_prev:k + 2*npassing + 2, 3) &
                            = source_vector(k + 2*npassing + 2 - npassing_prev:k + 2*npassing + 2, 3) &
                              + asource(m, 3)/2.4d0*q_rip(npassing_prev + 1:1:-1, istep + 1, 1) &
                              *fact_neg_b(istep + 1)
                        source_vector(k + 2*npassing + 2 - npassing_prev:k + 2*npassing + 2, 4) &
                            = source_vector(k + 2*npassing + 2 - npassing_prev:k + 2*npassing + 2, 4) &
                              + asource(m, 1)/2.4d0*q_rip_incompress(npassing_prev + 1:1:-1, istep + 1) &
                              *fact_neg_b(istep + 1)

                        source_vector(k + 2*npassing + 2 - npassing_next:k + 2*npassing + 2, 1) &
                            = source_vector(k + 2*npassing + 2 - npassing_next:k + 2*npassing + 2, 1) &
                              - asource(m, 1)/12d0*q_rip(npassing_next + 1:1:-1, istep - 1, 1) &
                              *fact_neg_e(istep - 1)
                        source_vector(k + 2*npassing + 2 - npassing_next:k + 2*npassing + 2, 2) &
                            = source_vector(k + 2*npassing + 2 - npassing_next:k + 2*npassing + 2, 2) &
                              + asource(m, 2)/12d0*q_rip(npassing_next + 1:1:-1, istep - 1, 2) &
                              *fact_neg_e(istep - 1)
                        source_vector(k + 2*npassing + 2 - npassing_next:k + 2*npassing + 2, 3) &
                            = source_vector(k + 2*npassing + 2 - npassing_next:k + 2*npassing + 2, 3) &
                              - asource(m, 3)/12d0*q_rip(npassing_next + 1:1:-1, istep - 1, 1) &
                              *fact_neg_e(istep - 1)
                        source_vector(k + 2*npassing + 2 - npassing_next:k + 2*npassing + 2, 4) &
                            = source_vector(k + 2*npassing + 2 - npassing_next:k + 2*npassing + 2, 4) &
                              - asource(m, 1)/12d0*q_rip_incompress(npassing_next + 1:1:-1, istep - 1) &
                              *fact_neg_e(istep - 1)
                    else
                        npassing_next = npl(istep + 2)
                        source_vector(k + npassing + 2:k + 2*npassing + 2, 1) &
                            = source_vector(k + npassing + 2:k + 2*npassing + 2, 1) &
                              + asource(m, 1)/2.4d0*q_rip(npassing + 1:1:-1, istep, 1) &
                              *fact_neg_e(istep)
                        source_vector(k + npassing + 2:k + 2*npassing + 2, 2) &
                            = source_vector(k + npassing + 2:k + 2*npassing + 2, 2) &
                              - asource(m, 2)/2.4d0*q_rip(npassing + 1:1:-1, istep, 2) &
                              *fact_neg_e(istep)
                        source_vector(k + npassing + 2:k + 2*npassing + 2, 3) &
                            = source_vector(k + npassing + 2:k + 2*npassing + 2, 3) &
                              + asource(m, 3)/2.4d0*q_rip(npassing + 1:1:-1, istep, 1) &
                              *fact_neg_e(istep)
                        source_vector(k + npassing + 2:k + 2*npassing + 2, 4) &
                            = source_vector(k + npassing + 2:k + 2*npassing + 2, 4) &
                              + asource(m, 1)/2.4d0*q_rip_incompress(npassing + 1:1:-1, istep) &
                              *fact_neg_e(istep)
                        if (npassing_prev .LE. npassing) then
                            source_vector(k + 2*npassing + 2 - npassing_prev:k + 2*npassing + 2, 1) &
                                = source_vector(k + 2*npassing + 2 - npassing_prev:k + 2*npassing + 2, 1) &
                                  + asource(m, 1)/1.5d0*q_rip(npassing_prev + 1:1:-1, istep + 1, 1) &
                                  *fact_neg_b(istep + 1)
                            source_vector(k + 2*npassing + 2 - npassing_prev:k + 2*npassing + 2, 2) &
                                = source_vector(k + 2*npassing + 2 - npassing_prev:k + 2*npassing + 2, 2) &
                                  - asource(m, 2)/1.5d0*q_rip(npassing_prev + 1:1:-1, istep + 1, 2) &
                                  *fact_neg_b(istep + 1)
                            source_vector(k + 2*npassing + 2 - npassing_prev:k + 2*npassing + 2, 3) &
                                = source_vector(k + 2*npassing + 2 - npassing_prev:k + 2*npassing + 2, 3) &
                                  + asource(m, 3)/1.5d0*q_rip(npassing_prev + 1:1:-1, istep + 1, 1) &
                                  *fact_neg_b(istep + 1)
                            source_vector(k + 2*npassing + 2 - npassing_prev:k + 2*npassing + 2, 4) &
                                = source_vector(k + 2*npassing + 2 - npassing_prev:k + 2*npassing + 2, 4) &
                                  + asource(m, 1)/1.5d0*q_rip_incompress(npassing_prev + 1:1:-1, istep + 1) &
                                  *fact_neg_b(istep + 1)
                        else
                            source_vector(k + npassing + 2:k + 2*npassing + 2, 1) &
                                = source_vector(k + npassing + 2:k + 2*npassing + 2, 1) &
                                  + asource(m, 1)/1.5d0*q_rip(npassing + 1:1:-1, istep + 1, 1) &
                                  *fact_neg_b(istep + 1)
                            source_vector(k_prev + npassing_prev + 1, 1) &
                                = source_vector(k_prev + npassing_prev + 1, 1) &
                                  + asource(m, 1)/1.5d0*q_rip(npassing_prev + 1, istep + 1, 1) &
                                  *fact_neg_b(istep + 1)
                            source_vector(k + npassing + 2:k + 2*npassing + 2, 2) &
                                = source_vector(k + npassing + 2:k + 2*npassing + 2, 2) &
                                  - asource(m, 2)/1.5d0*q_rip(npassing + 1:1:-1, istep + 1, 2) &
                                  *fact_neg_b(istep + 1)
                            source_vector(k_prev + npassing_prev + 1, 2) &
                                = source_vector(k_prev + npassing_prev + 1, 2) &
                                  - asource(m, 2)/1.5d0*q_rip(npassing_prev + 1, istep + 1, 2) &
                                  *fact_neg_b(istep + 1)
                            source_vector(k + npassing + 2:k + 2*npassing + 2, 3) &
                                = source_vector(k + npassing + 2:k + 2*npassing + 2, 3) &
                                  + asource(m, 3)/1.5d0*q_rip(npassing + 1:1:-1, istep + 1, 1) &
                                  *fact_neg_b(istep + 1)
                            source_vector(k_prev + npassing_prev + 1, 3) &
                                = source_vector(k_prev + npassing_prev + 1, 3) &
                                  + asource(m, 3)/1.5d0*q_rip(npassing_prev + 1, istep + 1, 1) &
                                  *fact_neg_b(istep + 1)
                            source_vector(k + npassing + 2:k + 2*npassing + 2, 4) &
                                = source_vector(k + npassing + 2:k + 2*npassing + 2, 4) &
                                  + asource(m, 1)/1.5d0*q_rip_incompress(npassing + 1:1:-1, istep + 1) &
                                  *fact_neg_b(istep + 1)
                            source_vector(k_prev + npassing_prev + 1, 4) &
                                = source_vector(k_prev + npassing_prev + 1, 4) &
                                  + asource(m, 1)/1.5d0*q_rip_incompress(npassing_prev + 1, istep + 1) &
                                  *fact_neg_b(istep + 1)
                        end if
                        if (npassing_next .LE. npassing) then
                            source_vector(k + 2*npassing + 2 - npassing_next:k + 2*npassing + 2, 1) &
                                = source_vector(k + 2*npassing + 2 - npassing_next:k + 2*npassing + 2, 1) &
                                  - asource(m, 1)/12d0*q_rip(npassing_next + 1:1:-1, istep + 2, 1) &
                                  *fact_neg_b(istep + 2)
                            source_vector(k + 2*npassing + 2 - npassing_next:k + 2*npassing + 2, 2) &
                                = source_vector(k + 2*npassing + 2 - npassing_next:k + 2*npassing + 2, 2) &
                                  + asource(m, 2)/12d0*q_rip(npassing_next + 1:1:-1, istep + 2, 2) &
                                  *fact_neg_b(istep + 2)
                            source_vector(k + 2*npassing + 2 - npassing_next:k + 2*npassing + 2, 3) &
                                = source_vector(k + 2*npassing + 2 - npassing_next:k + 2*npassing + 2, 3) &
                                  - asource(m, 3)/12d0*q_rip(npassing_next + 1:1:-1, istep + 2, 1) &
                                  *fact_neg_b(istep + 2)
                            source_vector(k + 2*npassing + 2 - npassing_next:k + 2*npassing + 2, 4) &
                                = source_vector(k + 2*npassing + 2 - npassing_next:k + 2*npassing + 2, 4) &
                                  - asource(m, 1)/12d0*q_rip_incompress(npassing_next + 1:1:-1, istep + 2) &
                                  *fact_neg_b(istep + 2)
                        else
                            source_vector(k + npassing + 2:k + 2*npassing + 2, 1) &
                                = source_vector(k + npassing + 2:k + 2*npassing + 2, 1) &
                                  - asource(m, 1)/12d0*q_rip(npassing + 1:1:-1, istep + 2, 1) &
                                  *fact_neg_b(istep + 2)
                            source_vector(k_prev + npassing_prev + 1, 1) &
                                = source_vector(k_prev + npassing_prev + 1, 1) &
                                  - asource(m, 1)/12d0*q_rip(npassing_next + 1, istep + 2, 1) &
                                  *fact_neg_b(istep + 2)
                            source_vector(k + npassing + 2:k + 2*npassing + 2, 2) &
                                = source_vector(k + npassing + 2:k + 2*npassing + 2, 2) &
                                  + asource(m, 2)/12d0*q_rip(npassing + 1:1:-1, istep + 2, 2) &
                                  *fact_neg_b(istep + 2)
                            source_vector(k_prev + npassing_prev + 1, 2) &
                                = source_vector(k_prev + npassing_prev + 1, 2) &
                                  + asource(m, 2)/12d0*q_rip(npassing_next + 1, istep + 2, 2) &
                                  *fact_neg_b(istep + 2)
                            source_vector(k + npassing + 2:k + 2*npassing + 2, 3) &
                                = source_vector(k + npassing + 2:k + 2*npassing + 2, 3) &
                                  - asource(m, 3)/12d0*q_rip(npassing + 1:1:-1, istep + 2, 1) &
                                  *fact_neg_b(istep + 2)
                            source_vector(k_prev + npassing_prev + 1, 3) &
                                = source_vector(k_prev + npassing_prev + 1, 3) &
                                  - asource(m, 3)/12d0*q_rip(npassing_next + 1, istep + 2, 1) &
                                  *fact_neg_b(istep + 2)
                            source_vector(k + npassing + 2:k + 2*npassing + 2, 4) &
                                = source_vector(k + npassing + 2:k + 2*npassing + 2, 4) &
                                  - asource(m, 1)/12d0*q_rip_incompress(npassing + 1:1:-1, istep + 2) &
                                  *fact_neg_b(istep + 2)
                            source_vector(k_prev + npassing_prev + 1, 4) &
                                = source_vector(k_prev + npassing_prev + 1, 4) &
                                  - asource(m, 1)/12d0*q_rip_incompress(npassing_next + 1, istep + 2) &
                                  *fact_neg_b(istep + 2)
                        end if
                    end if
                elseif (iplot .EQ. 1 .AND. isw_axisymm .NE. 1) then
                    source_vector(k + npassing + 2:k + 2*npassing + 2, 1:3) &
                        = flux_mr(npass_r*(m + 1):npass_r*m + 1:-1, :)
                end if

            end do
        end do

    end subroutine source_flux

!------------------------------------------------------------------------
    subroutine add_f01_source

        integer :: k_next

        do istep = ibeg, iend

            ioddeven = MOD(istep - ibeg, 2) !0 for full RK step, 1 for half RK step

            npassing = npl(istep)

            do m = 0, lag
                k = ind_start(istep) + 2*(npassing + 1)*m

                if (istep .GT. ibeg) then
                    npassing_prev = npl(istep - 1)
                    k_prev = ind_start(istep - 1) + 2*(npassing_prev + 1)*m
                    if (ioddeven .EQ. 1) then
                        npassing_next = npl(istep + 1)
                        k_next = ind_start(istep + 1) + 2*(npassing_next + 1)*m
                        source_vector(k + 1:k + npassing + 1, 1:3) &
                            = source_vector(k + 1:k + npassing + 1, 1:3) &
                              + (bnoverb0(istep)*f0_coll(k + 1:k + npassing + 1, :) &
                                 + ttmpfact(istep)*f0_ttmp(k + 1:k + npassing + 1, :)) &
                              *fact_pos_e(istep)/1.5d0
                        source_vector(k + 1:k + npassing_prev + 1, 1:3) &
                            = source_vector(k + 1:k + npassing_prev + 1, 1:3) &
                              + (bnoverb0(istep - 1)*f0_coll(k_prev + 1:k_prev + npassing_prev + 1, :) &
                                 + ttmpfact(istep - 1)*f0_ttmp(k_prev + 1:k_prev + npassing_prev + 1, :)) &
                              *fact_pos_e(istep - 1)/2.4d0
                        source_vector(k + 1:k + npassing_next + 1, 1:3) &
                            = source_vector(k + 1:k + npassing_next + 1, 1:3) &
                              - (bnoverb0(istep + 1)*f0_coll(k_next + 1:k_next + npassing_next + 1, :) &
                                 + ttmpfact(istep + 1)*f0_ttmp(k_next + 1:k_next + npassing_next + 1, :)) &
                              *fact_pos_e(istep + 1)/12d0
                    else
                        npassing_next = npl(istep - 2)
                        k_next = ind_start(istep - 2) + 2*(npassing_next + 1)*m
                        source_vector(k + 1:k + npassing + 1, 1:3) &
                            = source_vector(k + 1:k + npassing + 1, 1:3) &
                              + (bnoverb0(istep)*f0_coll(k + 1:k + npassing + 1, :) &
                                 + ttmpfact(istep)*f0_ttmp(k + 1:k + npassing + 1, :)) &
                              *fact_pos_e(istep)/2.4d0
                        if (npassing_prev .LE. npassing) then
                            source_vector(k + 1:k + npassing_prev + 1, 1:3) &
                                = source_vector(k + 1:k + npassing_prev + 1, 1:3) &
                                  + (bnoverb0(istep - 1)*f0_coll(k_prev + 1:k_prev + npassing_prev + 1, :) &
                                     + ttmpfact(istep - 1)*f0_ttmp(k_prev + 1:k_prev + npassing_prev + 1, :)) &
                                  *fact_pos_b(istep - 1)/1.5d0
                        else
                            source_vector(k + 1:k + npassing + 1, 1:3) &
                                = source_vector(k + 1:k + npassing + 1, 1:3) &
                                  + (bnoverb0(istep - 1)*f0_coll(k_prev + 1:k_prev + npassing + 1, :) &
                                     + ttmpfact(istep - 1)*f0_ttmp(k_prev + 1:k_prev + npassing + 1, :)) &
                                  *fact_pos_b(istep - 1)/1.5d0
                            source_vector(k_prev + npassing_prev + 2, 1:3) &
                                = source_vector(k_prev + npassing_prev + 2, 1:3) &
                                  + (bnoverb0(istep - 1)*f0_coll(k_prev + npassing_prev + 1, :) &
                                     + ttmpfact(istep - 1)*f0_ttmp(k_prev + npassing_prev + 1, :)) &
                                  *fact_pos_b(istep - 1)/1.5d0
                        end if
                        if (npassing_next .LE. npassing) then
                            source_vector(k + 1:k + npassing_next + 1, 1:3) &
                                = source_vector(k + 1:k + npassing_next + 1, 1:3) &
                                  - (bnoverb0(istep - 2)*f0_coll(k_next + 1:k_next + npassing_next + 1, :) &
                                     + ttmpfact(istep - 2)*f0_ttmp(k_next + 1:k_next + npassing_next + 1, :)) &
                                  *fact_pos_b(istep - 2)/12d0
                        else
                            source_vector(k + 1:k + npassing + 1, 1:3) &
                                = source_vector(k + 1:k + npassing + 1, 1:3) &
                                  - (bnoverb0(istep - 2)*f0_coll(k_next + 1:k_next + npassing + 1, :) &
                                     + ttmpfact(istep - 2)*f0_ttmp(k_next + 1:k_next + npassing + 1, :)) &
                                  *fact_pos_b(istep - 2)/12d0
                            source_vector(k_prev + npassing_prev + 2, 1:3) &
                                = source_vector(k_prev + npassing_prev + 2, 1:3) &
                                  - (bnoverb0(istep - 2)*f0_coll(k_next + npassing_next + 1, :) &
                                     + ttmpfact(istep - 2)*f0_ttmp(k_next + npassing_next + 1, :)) &
                                  *fact_pos_b(istep - 2)/12d0
                        end if
                    end if
                end if

                if (istep .LT. iend) then
                    npassing_prev = npl(istep + 1)
                    k_prev = ind_start(istep + 1) + 2*(npassing_prev + 1)*m
                    if (ioddeven .EQ. 1) then
                        npassing_next = npl(istep - 1)
                        k_next = ind_start(istep - 1) + 2*(npassing_next + 1)*m
                        source_vector(k + npassing + 2:k + 2*npassing + 2, 1:3) &
                            = source_vector(k + npassing + 2:k + 2*npassing + 2, 1:3) &
                              + (bnoverb0(istep)*f0_coll(k + npassing + 2:k + 2*npassing + 2, :) &
                                 + ttmpfact(istep)*f0_ttmp(k + npassing + 2:k + 2*npassing + 2, :)) &
                              *fact_neg_e(istep)/1.5d0
                        source_vector(k + 2*npassing + 2 - npassing_prev:k + 2*npassing + 2, 1:3) &
                            = source_vector(k + 2*npassing + 2 - npassing_prev:k + 2*npassing + 2, 1:3) &
                              + (bnoverb0(istep + 1)*f0_coll(k_prev + npassing_prev + 2: &
                                                             k_prev + 2*npassing_prev + 2, :) &
                                 + ttmpfact(istep + 1)*f0_ttmp(k_prev + npassing_prev + 2: &
                                                               k_prev + 2*npassing_prev + 2, :)) &
                              *fact_neg_b(istep + 1)/2.4d0
                        source_vector(k + 2*npassing + 2 - npassing_next:k + 2*npassing + 2, 1:3) &
                            = source_vector(k + 2*npassing + 2 - npassing_next:k + 2*npassing + 2, 1:3) &
                              - (bnoverb0(istep - 1)*f0_coll(k_next + npassing_next + 2: &
                                                             k_next + 2*npassing_next + 2, :) &
                                 + ttmpfact(istep - 1)*f0_ttmp(k_next + npassing_next + 2: &
                                                               k_next + 2*npassing_next + 2, :)) &
                              *fact_neg_e(istep - 1)/12d0
                    else
                        npassing_next = npl(istep + 2)
                        k_next = ind_start(istep + 2) + 2*(npassing_next + 1)*m
                        source_vector(k + npassing + 2:k + 2*npassing + 2, 1:3) &
                            = source_vector(k + npassing + 2:k + 2*npassing + 2, 1:3) &
                              + (bnoverb0(istep)*f0_coll(k + npassing + 2:k + 2*npassing + 2, :) &
                                 + ttmpfact(istep)*f0_ttmp(k + npassing + 2:k + 2*npassing + 2, :)) &
                              *fact_neg_e(istep)/2.4d0
                        if (npassing_prev .LE. npassing) then
                            source_vector(k + 2*npassing + 2 - npassing_prev:k + 2*npassing + 2, 1:3) &
                                = source_vector(k + 2*npassing + 2 - npassing_prev:k + 2*npassing + 2, 1:3) &
                                  + (bnoverb0(istep + 1)*f0_coll(k_prev + npassing_prev + 2: &
                                                                 k_prev + 2*npassing_prev + 2, :) &
                                     + ttmpfact(istep + 1)*f0_ttmp(k_prev + npassing_prev + 2: &
                                                                   k_prev + 2*npassing_prev + 2, :)) &
                                  *fact_neg_b(istep + 1)/1.5d0
                        else
                            source_vector(k + npassing + 2:k + 2*npassing + 2, 1:3) &
                                = source_vector(k + npassing + 2:k + 2*npassing + 2, 1:3) &
                                  + (bnoverb0(istep + 1)*f0_coll(k_prev + 2*npassing_prev + 2 - npassing: &
                                                                 k_prev + 2*npassing_prev + 2, :) &
                                     + ttmpfact(istep + 1)*f0_ttmp(k_prev + 2*npassing_prev + 2 - npassing: &
                                                                   k_prev + 2*npassing_prev + 2, :)) &
                                  *fact_neg_b(istep + 1)/1.5d0
                            source_vector(k_prev + npassing_prev + 1, 1:3) &
                                = source_vector(k_prev + npassing_prev + 1, 1:3) &
                                  + (bnoverb0(istep + 1)*f0_coll(k_prev + npassing_prev + 2, :) &
                                     + ttmpfact(istep + 1)*f0_ttmp(k_prev + npassing_prev + 2, :)) &
                                  *fact_neg_b(istep + 1)/1.5d0
                        end if
                        if (npassing_next .LE. npassing) then
                            source_vector(k + 2*npassing + 2 - npassing_next:k + 2*npassing + 2, 1:3) &
                                = source_vector(k + 2*npassing + 2 - npassing_next:k + 2*npassing + 2, 1:3) &
                                  - (bnoverb0(istep + 2)*f0_coll(k_next + npassing_next + 2: &
                                                                 k_next + 2*npassing_next + 2, :) &
                                     + ttmpfact(istep + 2)*f0_ttmp(k_next + npassing_next + 2: &
                                                                   k_next + 2*npassing_next + 2, :)) &
                                  *fact_neg_b(istep + 2)/12d0
                        else
                            source_vector(k + npassing + 2:k + 2*npassing + 2, 1:3) &
                                = source_vector(k + npassing + 2:k + 2*npassing + 2, 1:3) &
                                  - (bnoverb0(istep + 2)*f0_coll(k_next + 2*npassing_next + 2 - npassing: &
                                                                 k_next + 2*npassing_next + 2, :) &
                                     + ttmpfact(istep + 2)*f0_ttmp(k_next + 2*npassing_next + 2 - npassing: &
                                                                   k_next + 2*npassing_next + 2, :)) &
                                  *fact_neg_b(istep + 2)/12d0
                            source_vector(k_prev + npassing_prev + 1, 1:3) &
                                = source_vector(k_prev + npassing_prev + 1, 1:3) &
                                  - (bnoverb0(istep + 2)*f0_coll(k_next + npassing_next + 2, :) &
                                     + ttmpfact(istep + 2)*f0_ttmp(k_next + npassing_next + 2, :)) &
                                  *fact_neg_b(istep + 2)/12d0
                        end if
                    end if
                end if

            end do
        end do

    end subroutine add_f01_source

!------------------------------------------------------------------------
    subroutine integral_part(vec_in, vec_out)

        implicit none

        integer :: l, m, i, k, istep, npassing, k_prev

        complex(dp), dimension(n_2d_size) :: vec_in, vec_out
    !! Modification by Andreas F. Martitsch (20.08.2015)
        ! Array extended by 3rd (phi-steps) and 4th dimension (species)
        complex(dp), dimension(0:lag, 0:leg, ibeg:iend, 0:num_spec - 1) :: scalprod_pleg
        complex(dp), dimension(0:lag, 0:leg, ibeg:iend, 0:num_spec - 1) :: scalprod_pleg_tmp
        ! Species index
        integer :: ispecp
    !! End Modification by Andreas F. Martitsch (20.08.2015)
        complex(dp), dimension(:, :, :), allocatable :: vec_tmp

        allocate (vec_tmp(0:lag, 2*(npart + 1), ibeg:iend))
        vec_tmp = 0.d0

!! Modification by Andreas F. Martitsch (20.08.2015)
! Array scalprod_pleg extended by 3rd (phi-steps) and
! 4th dimension (species)
        do istep = ibeg, iend

            npassing = npl(istep)

            do m = 0, lag
                k = ind_start(istep) + 2*(npassing + 1)*m
                do l = 0, leg
                    scalprod_pleg(m, l, istep, ispec) = &
                        SUM(pleg_bra(l, 1:npassing + 1, istep)*vec_in(k + 1:k + npassing + 1))
                end do
                k = k + 2*(npassing + 1)
                do l = 0, leg, 2
                    scalprod_pleg(m, l, istep, ispec) = scalprod_pleg(m, l, istep, ispec) &
                                                        + SUM(pleg_bra(l, 1:npassing + 1, istep)*vec_in(k:k - npassing:-1))
                end do
                do l = 1, leg, 2
                    scalprod_pleg(m, l, istep, ispec) = scalprod_pleg(m, l, istep, ispec) &
                                                        - SUM(pleg_bra(l, 1:npassing + 1, istep)*vec_in(k:k - npassing:-1))
                end do
            end do
        end do

! Finish filling-up array scalprod_pleg
!! End Modification by Andreas F. Martitsch (20.08.2015)

!! Modification by Andreas F. Martitsch (20.08.2015)
! MPI Barrier -> collect scalprod (4D - leg,lag,phi,species)
! (mpro%allgather supports 3D and 4D matrices)
        call mpro%allgather_inplace(scalprod_pleg)

        do istep = ibeg, iend

            npassing = npl(istep)

! ailmm is now 5D object of species (alpha,alphap)

            do l = 0, leg
                scalprod_pleg_tmp(0:lag, l, istep, ispec) = 0.0d0
                do ispecp = 0, num_spec - 1
                    scalprod_pleg_tmp(0:lag, l, istep, ispec) = scalprod_pleg_tmp(0:lag, l, istep, ispec) + &
                                                                MATMUL(ailmm_aa(0:lag, 0:lag, l, ispec, ispecp), &
                                                                       scalprod_pleg(0:lag, l, istep, ispecp))
                end do
                scalprod_pleg(0:lag, l, istep, ispec) = scalprod_pleg_tmp(0:lag, l, istep, ispec)
                ! old behavior (for a single species)
                !scalprod_pleg(0:lag,l)=MATMUL(ailmm(0:lag,0:lag,l),  &
                !                              scalprod_pleg(0:lag,l))
            end do

! end of interaction with rest processors

            do m = 0, lag

                do l = 0, leg
                    vec_tmp(m, 1:npassing + 1, istep) = vec_tmp(m, 1:npassing + 1, istep) &
                                                        + scalprod_pleg(m, l, istep, ispec)*pleg_ket(l, 1:npassing + 1, istep)
                end do

                k = 2*(npassing + 1)

                do l = 0, leg, 2
                    vec_tmp(m, k:k - npassing:-1, istep) = vec_tmp(m, k:k - npassing:-1, istep) &
                                                           + scalprod_pleg(m, l, istep, ispec)*pleg_ket(l, 1:npassing + 1, istep)
                end do
                do l = 1, leg, 2
                    vec_tmp(m, k:k - npassing:-1, istep) = vec_tmp(m, k:k - npassing:-1, istep) &
                                                           - scalprod_pleg(m, l, istep, ispec)*pleg_ket(l, 1:npassing + 1, istep)
                end do

            end do

        end do

! Finish computations with scalprod_pleg
!! End Modification by Andreas F. Martitsch (20.08.2015)

        vec_tmp = 0.5d0*vec_tmp

        vec_out = 0.d0

! forwards:
        do istep = ibeg + 1, iend

            do m = 0, lag
                npassing = npl(istep)
                k = ind_start(istep) + 2*(npassing + 1)*m

                vec_out(k + 1:k + npassing + 1) = vec_out(k + 1:k + npassing + 1) &
                                                  + vec_tmp(m, 1:npassing + 1, istep) &
                                                  *fact_pos_e(istep)

                npassing = MIN(npassing, npl(istep - 1))
                vec_out(k + 1:k + npassing + 1) = vec_out(k + 1:k + npassing + 1) &
                                                  + vec_tmp(m, 1:npassing + 1, istep - 1) &
                                                  *fact_pos_e(istep)
                npassing = npl(istep)
                if (npassing .gt. npl(istep - 1) .and. addboucol) then
                    vec_out(k + npassing + 1) = vec_out(k + npassing + 1)       &
                                         & + vec_tmp(m, npassing + 2, istep) &
                                         & *fact_pos_e(istep)
                end if
            end do

        end do

! backwards:
        do istep = ibeg, iend - 1

            do m = 0, lag
                npassing = npl(istep)
                k_prev = 2*(npassing + 1)
                k = ind_start(istep) + 2*(npassing + 1)*(m + 1)

                vec_out(k - npassing:k) = vec_out(k - npassing:k) &
                                          + vec_tmp(m, k_prev - npassing:k_prev, istep) &
                                          *fact_neg_e(istep)

                npassing = MIN(npassing, npl(istep + 1))
                k_prev = 2*(npl(istep + 1) + 1)
                vec_out(k - npassing:k) = vec_out(k - npassing:k) &
                                          + vec_tmp(m, k_prev - npassing:k_prev, istep + 1) &
                                          *fact_neg_e(istep)
                npassing = npl(istep)
                if ((npassing .gt. npl(istep + 1)) .and. addboucol) then
                    vec_out(k - npassing) = vec_out(k - npassing)         &
                                        & + vec_tmp(m, npassing + 1, istep) &
                                        & *fact_neg_e(istep)
                end if
            end do

        end do

        deallocate (vec_tmp)

    end subroutine integral_part

    !---------------------------------------------------------------------------------
    !> Single step of a preconditioned Richardson iteration which for equation set
    !>
    !>     L f = q,
    !>
    !> where L = L_V - L_D - L_I is evolution operator with L_V being Vlasov operator
    !> L_D - differential part of collision operator and L_I - integral part of collision operator.
    !>
    !> defines next iteration f_new via previous iteration f_old as follows:
    !>
    !>     f_new = f_old + Lbar (q - L f_old)
    !>
    !> Here Lbar (preconditioner) is the approximate inverse to L_V-L_D computed by sparse solver
    !> In case mode=2 skips q in the iteration (used by Arnoldi iterations)
    subroutine next_iteration(n, fold, fnew)

        use arnoldi_mod, only: fzero, mode, tol

        implicit none

        integer :: n
        complex(dp), dimension(n) :: fold, fnew
        real(dp), dimension(:), allocatable :: fnew_real, fnew_imag

        if (isw_intp .eq. 1) then
            call integral_part(fold, fnew)
        else
            fnew = (0.0d0, 0.0d0)
        end if

        if (mode .ne. 2) then
            fnew = fnew + fzero
        end if

        if (problem_type) then
            do i = 1, nz_symm
                fnew(irow_symm(i)) = fnew(irow_symm(i)) - amat_symm(i)*fold(icol_symm(i))
            end do
            allocate (fnew_real(n), fnew_imag(n))
            fnew_real = DBLE(fnew)
            fnew_imag = AIMAG(fnew)
            call sparse_solve(nrow, ncol, nz, irow(1:nz), ipcol, DBLE(amat_sp(1:nz)), &
                              fnew_real, iopt)
            call sparse_solve(nrow, ncol, nz, irow(1:nz), ipcol, DBLE(amat_sp(1:nz)), &
                              fnew_imag, iopt)
            fnew = fnew_real + (0.d0, 1.d0)*fnew_imag
            deallocate (fnew_real, fnew_imag)

            fnew = fnew + fold

            ! -> remove null-space of axisymmetric
            ! solution (particle and energy conservation)
            coef_dens = sum(densvec_bra*fnew)/denom_dens
            fnew = fnew - (1.d0 - tol)*coef_dens*densvec_ket

            if (isw_intp .eq. 1) then
                coef_energ = sum(energvec_bra*fnew)/denom_energ
                fnew = fnew - (1.d0 - tol)*coef_energ*energvec_ket
            end if
        else

            do i = 1, nz
                fnew(irow(i)) = fnew(irow(i)) - amat_sp(i)*fold(icol(i))
            end do
            call sparse_solve(nrow, ncol, nz, irow(1:nz), ipcol, amat_sp(1:nz), &
                              fnew, iopt)

            fnew = fnew + fold
        end if

    end subroutine next_iteration

    subroutine save_qflux_symm_allspec()
        implicit none

        if (lsw_multispecies .AND. mpro%getrank() .EQ. 0) then
            open (070915, file='qflux_symm_allspec.dat')
            write (070915, *) '% boozer_s, collpar'
            write (070915, *) boozer_s, collpar
            write (070915, *) '% qflux(1,1,a,b}'
            do ispecp = 0, num_spec - 1
                write (070915, *) (qflux_symm_allspec(1, 1, ispecp, ispecpp),&
                  & ispecpp=0, num_spec - 1)
            end do
            write (070915, *) '% qflux(1,3,a,b}'
            do ispecp = 0, num_spec - 1
                write (070915, *) (qflux_symm_allspec(1, 3, ispecp, ispecpp),&
                  & ispecpp=0, num_spec - 1)
            end do
            write (070915, *) '% qflux(1,2,a,b}'
            do ispecp = 0, num_spec - 1
                write (070915, *) (qflux_symm_allspec(1, 2, ispecp, ispecpp),&
                  & ispecpp=0, num_spec - 1)
            end do
            write (070915, *) '% qflux(2,2,a,b}'
            do ispecp = 0, num_spec - 1
                write (070915, *) (qflux_symm_allspec(2, 2, ispecp, ispecpp),&
                  & ispecpp=0, num_spec - 1)
            end do
            write (070915, *) '% qflux(3,3,a,b}'
            do ispecp = 0, num_spec - 1
                write (070915, *) (qflux_symm_allspec(3, 3, ispecp, ispecpp),&
                  & ispecpp=0, num_spec - 1)
            end do
            close (070915)

        end if

    end subroutine save_qflux_symm_allspec

    !> ifwrite_coords=.true. writes phase space coordinates (phi, lambda, base functions)
    !> ifwrite_coords=.false. writes only the distribution function
    !
    !> prefix is added to file names of the distribution function
    !> ifnormalize_output=.true.  normalizes output as distribution function (divides f_k by delta eta)
    !> ifnormalize_output=.false. no normalization (plots f_k as is)
    subroutine matlabplot_allm(sourcevec_tmp, normalize_output, write_coords, prefix)

        use collisionality_mod, only: phi_x_max
        use collop_bspline, only: init_phi_bspline, phi_bspline

        implicit none

        logical, intent(in) :: normalize_output, write_coords
        double precision, dimension(n_2d_size), intent(in) :: sourcevec_tmp
        character(len=2), intent(in) :: prefix

        integer, parameter :: nx = 100

        logical :: dir_exisits

        integer :: iunit_base, nmax, m, i, k
        double precision, dimension(:, :), allocatable :: phi_mat, alam_mat, fun_mat
        double precision   :: hx
        double precision, dimension(:), allocatable :: xi
        double precision, dimension(:, :), allocatable :: dummy2d

        ! Check iffolder with species number exists, ifnot return witout
        ! writing the output.
        inquire (file=char(48 + ispec), exist=dir_exisits)
        if (.not. dir_exisits) return

        nmax = maxval(npl) + 1
        allocate (phi_mat(ibeg:iend, -nmax:nmax))
        allocate (alam_mat(ibeg:iend, -nmax:nmax))
        allocate (fun_mat(ibeg:iend, -nmax:nmax))

        iunit_base = 12345
        alam_mat = 10.d0

        do m = 0, lag
            fun_mat = 0.d0

            do istep = ibeg, iend

                phi_mat(istep, :) = phi_mfl(istep)

                npassing = npl(istep)
                delta_eta = eta(1:npassing) - eta(0:npassing - 1)
                eta0 = 1.0d0/bhat_mfl(istep)
                k = ind_start(istep) + 2*(npassing + 1)*m
                do i = 1, npassing + 1
                    if (i .le. npassing) then
                        alam_mat(istep, i - npassing - 1) = 0.5d0*(alambd(i, istep) + alambd(i - 1, istep))
                        if (normalize_output) then
                            fun_mat(istep, i - npassing - 1) = sourcevec_tmp(k + 2*(npassing + 1) - i + 1)/delta_eta(i)
                        else
                            fun_mat(istep, i - npassing - 1) = sourcevec_tmp(k + 2*(npassing + 1) - i + 1)
                        end if
                    else
                        alam_mat(istep, i - npassing - 1) = 0.5d0*alambd(i - 1, istep)
                        if (normalize_output) then
                            fun_mat(istep, i - npassing - 1) = sourcevec_tmp(k + 2*(npassing + 1) - i + 1)/(eta0 - eta(i - 1))
                        else
                            fun_mat(istep, i - npassing - 1) = sourcevec_tmp(k + 2*(npassing + 1) - i + 1)
                        end if
                    end if
                end do
                do i = npassing + 1, 1, -1
                    if (i .le. npassing) then
                        alam_mat(istep, npassing + 2 - i) = -0.5d0*(alambd(i, istep) + alambd(i - 1, istep))
                        if (normalize_output) then
                            fun_mat(istep, npassing + 2 - i) = sourcevec_tmp(k + i)/delta_eta(i)
                        else
                            fun_mat(istep, npassing + 2 - i) = sourcevec_tmp(k + i)
                        end if
                    else
                        alam_mat(istep, npassing + 2 - i) = -0.5d0*alambd(i - 1, istep)
                        if (normalize_output) then
                            fun_mat(istep, npassing + 2 - i) = sourcevec_tmp(k + i)/(eta0 - eta(i - 1))
                        else
                            fun_mat(istep, npassing + 2 - i) = sourcevec_tmp(k + i)
                        end if
                    end if
                end do
            end do

            open (iunit_base, file=char(48 + ispec)//'/fun_matlab'//trim(prefix)//char(48 + m)//'.dat')
            do istep = ibeg, iend
                write (iunit_base, *) fun_mat(istep, :)
            end do
            close (iunit_base)
        end do

        if (write_coords) then
            open (iunit_base, file=char(48 + ispec)//'/phi_matlab.dat')
            do istep = ibeg, iend
                write (iunit_base, *) phi_mat(istep, :)
            end do
            close (iunit_base)

            open (iunit_base, file=char(48 + ispec)//'/lambda_matlab.dat')
            do istep = ibeg, iend
                write (iunit_base, *) alam_mat(istep, :)
            end do
            close (iunit_base)

            call init_phi_bspline(lag, 0)

            allocate (xi(0:nx), dummy2d(0:nx, 0:lag))
            hx = phi_x_max/dble(nx)
            do i = 0, nx
                xi(i) = dble(i)*hx
                do m = 0, lag
                    dummy2d(i, m) = phi_bspline(m, xi(i))
                end do
            end do
            open (iunit_base, file=char(48 + ispec)//'/exp_basis.dat')
            do i = 0, nx
                write (iunit_base, *) xi(i), dummy2d(i, :)
            end do
            close (iunit_base)
            deallocate (xi, dummy2d)
        end if

    end subroutine matlabplot_allm

    !> Generates Maxwellian distribution in discretized form. Needed for filtering out null-space and for testing.
    subroutine generate_maxwellian(fmaxw)

        implicit none

        complex(dp), dimension(n_2d_size), intent(out) :: fmaxw

        do istep = ibeg, iend

            npassing = npl(istep)

            do m = 0, lag
                k = ind_start(istep) + 2*(npassing + 1)*m
                k_prev = k + 2*npassing + 3
                do ipart = 1, npassing
                    fmaxw(k + ipart) = (eta(ipart) - eta(ipart - 1))*asource(m, 2)
                    fmaxw(k_prev - ipart) = fmaxw(k + ipart)
                end do
                fmaxw(k + npassing + 1) = (1.d0/bhat_mfl(istep) - eta(npassing))*asource(m, 2)
                fmaxw(k + npassing + 2) = fmaxw(k + npassing + 1)
            end do
        end do

    end subroutine generate_maxwellian

    !> Test of particle conservation by the finite difference scheme
    !> Computes integrals over phase space (velocity space and flux tube) of the RHS
    !> of kinetic equation, (L_V-L_D) f, where L_V - Vlasov operator and L_D is the
    !> differential part of kinetic equation for arbitrary distribution functions f.
    !> Integral is over phi, eta and v with Jacobian
    !> $\pi B^\varphi \sqrt{g} v^2 /(h^\varphi |\lambda|) = C v^2 /(h^\varphi |\lambda|)$
    !> where C=const for the flux tube. Arbitrary distribution functions f are
    !> represented by f_ikm = 1 for i=i0, k=k0, m=m0 and f_ikm = 0 for all other ikm
    !> where i enumerates phi knots, k enumerates eta levels and m enumerates base functions.
    !> Due to (independent) particle conservation by L_V and L_D all such integrals should
    !> be zero. In absence of energy diffusion operator the following relations should be held:
    !> $\sum_I \delta \varphi_I A_{IJ}=0$
    !> for any J where $A_{IJ}$ is the RHS matrix for L_V-L_D and $\delta \varphi_I$ is integration
    !> step over $\varphi$ for the generalized index I. Generalized indices I and J are 3D vector
    !> indices, i.e. I=ikm. Namely, $\delta \varphi_I = \delta \varphi_i^\sigma$ where $\sigma$ is
    !> parallel velocity sign, $\delta \varphi_i^+=\varphi_i-\varphi_{i-1}$ and
    !> $\delta \varphi_i^-=\varphi_i-\varphi_{i+1}$.
    !> In presence of energy diffusion operator the above relation is generalized to
    !> $\sum_I \delta \varphi_I w_I A_{IJ}=0$
    !> where $w_I = w_m$ is velocity module weight for computation of particle density from the source term
    !> corresponding to the definition of the scalar product. E.g., for the default scalar product,
    !> $<a|b> = C \int \rd x x^4 exp(-x^2) a(x) b(x)$ weight $w_m$ corresponds to expansion coefficients
    !> of function 1/x in basis functions, $1/x = \sum_m w_m \phi_m(x)$.
    subroutine test_conservation(nz, irow, icol, amat)

        use collop, only: scalprod_alpha

        implicit none

        integer, intent(in) :: nz
        integer, dimension(nz), intent(in) :: irow, icol
        double complex, dimension(nz), intent(in) :: amat

        integer :: i, k, m, npassing
        double complex, dimension(:), allocatable :: fun
        double precision, dimension(:), allocatable :: step_of_irow
        double precision, dimension(0:lag) :: densint

        if (scalprod_alpha .eq. 0.0d0) then
            densint = weightlag(2, :)  !<= valid only for default scalar product
        else
            densint = 1  !<= valid only for alpha=-1 (?)
        end if

        allocate (step_of_irow(n_2d_size))
        step_of_irow = 0.0d0

        do istep = ibeg, iend
            npassing = npl(istep)
            do m = 0, lag
                k = ind_start(istep) + 2*(npassing + 1)*m
                step_of_irow(k + 1:k + npassing + 1) = delt_pos(istep)*densint(m)
                step_of_irow(k + npassing + 2:k + 2*(npassing + 1)) = delt_neg(istep)*densint(m)
            end do
        end do

        allocate (fun(n_2d_size))
        fun = (0.d0, 0.d0)

        do i = 1, nz
            fun(icol(i)) = fun(icol(i)) + amat(i)*step_of_irow(irow(i))
        end do

        call matlabplot_allm(real(fun), .false., .true., '  ')

        deallocate (fun, step_of_irow)

    end subroutine test_conservation

end subroutine ripple_solver_ArnoldiO2

!> Simple wrapper to make Arnoldi's sparse solving testable
!! This extracts just the core sparse solve logic used by Arnoldi
!! without requiring all the NEO-2 specific setup
subroutine arnoldi_sparse_solve_test(nrow, ncol, nz, irow, pcol, val, b)
    use sparse_mod, only: sparse_solve_method, sparse_solve
    implicit none
    
    ! Arguments
    integer, intent(in) :: nrow, ncol, nz
    integer, dimension(:), intent(in) :: irow, pcol
    real(kind=8), dimension(:), intent(in) :: val
    real(kind=8), dimension(:), intent(inout) :: b
    
    ! Local variables
    integer :: old_method, iopt
    
    ! Save current method
    old_method = sparse_solve_method
    
    ! Use Arnoldi's preferred sparse method (UMFPACK without iterative refinement)
    sparse_solve_method = 3
    
    ! Use the same solve approach as in Arnoldi
    iopt = 0  ! Full solve (factorize + solve + free)
    
    call sparse_solve(nrow, ncol, nz, irow, pcol, val, b, iopt)
    
    ! Restore original method
    sparse_solve_method = old_method
    
end subroutine arnoldi_sparse_solve_test

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!! Modifications by Andreas F. Martitsch (27.07.2015)
! Multiple definitions avoided
subroutine rearrange_phideps(ibeg, iend, npart, ncomp, nreal, subsqmin, phi_divide, &
                             phi_mfl, bhat_mfl, arr_real, arr_comp, eta, &
                             delt_pos, delt_neg, &
                             fact_pos_b, fact_neg_b, fact_pos_e, fact_neg_e)

! Mnemonics:
! fact_pos_b(i) - integration step in positive direction starts at point i
! fact_pos_e(i) - integration step in positive direction ends at point i
! fact_neg_b(i) - integration step in negative direction starts at point i
! fact_neg_e(i) - integration step in negative direction ends at point i

    use plagrange_mod

    implicit none

    integer, parameter :: dp = kind(1.0d0)

    logical, parameter :: stepmode = .FALSE.
    integer, parameter :: npoi = 6, nder = 0, npoihalf = npoi/2, nstepmin = 8
    real(dp), parameter :: bparabmax = 0.2d0

    integer :: i, ibeg, iend, npart, istep, ibmin, npassing, npassing_prev
    integer :: ncomp, nreal
    integer :: ncross_l, ncross_r, ib, ie, intb, inte, k, imid, isplit

    real(dp) :: subsqmin, ht, ht2, bparab, x1, x2, f1, f2

    integer, dimension(1)              :: idummy
    integer, dimension(1:iend)         :: phi_divide
    integer, dimension(:), allocatable :: icross_l, icross_r

    real(dp), dimension(0:nder, npoi)      :: coeff
    real(dp), dimension(0:npart)          :: eta
    real(dp), dimension(ibeg:iend)        :: phi_mfl, bhat_mfl
    real(dp), dimension(ibeg:iend, nreal)  :: arr_real
    complex(dp), dimension(ibeg:iend, ncomp)  :: arr_comp
    real(dp), dimension(ibeg:iend)        :: delt_pos, delt_neg
    real(dp), dimension(ibeg:iend)        :: fact_pos_b, fact_neg_b
    real(dp), dimension(ibeg:iend)        :: fact_pos_e, fact_neg_e
    real(dp), dimension(:), allocatable   :: phi_new, bhat_new
    real(dp), dimension(:, :), allocatable :: arr_real_new
    complex(dp), dimension(:, :), allocatable :: arr_comp_new

    npassing = -1

    call fix_phiplacement_problem(ibeg, iend, npart, subsqmin, &
                                  phi_mfl, bhat_mfl, eta)

    phi_divide = 1

    delt_pos(ibeg + 1:iend) = phi_mfl(ibeg + 1:iend) - phi_mfl(ibeg:iend - 1)
    fact_pos_b = 1.d0
    fact_pos_e = 1.d0

! determine level crossings:

    idummy = MINLOC(bhat_mfl(ibeg:iend))
    ibmin = idummy(1) + ibeg - 1

    ncross_l = 0
    if (ibmin .GT. ibeg) then
        istep = ibmin
        do i = 0, npart
            if (1.d0 - bhat_mfl(istep)*eta(i) .GT. subsqmin) then
                npassing = i
            else
                EXIT
            end if
        end do
        npassing_prev = npassing
        do istep = ibmin - 1, ibeg, -1
            do i = 0, npart
                if (1.d0 - bhat_mfl(istep)*eta(i) .GT. subsqmin) then
                    npassing = i
                else
                    EXIT
                end if
            end do
            if (npassing .LT. npassing_prev) then
                ncross_l = ncross_l + 1
                npassing_prev = npassing
            end if
        end do
        if (ncross_l .GT. 0) then
            allocate (icross_l(ncross_l))
            ncross_l = 0
            istep = ibmin
            do i = 0, npart
                if (1.d0 - bhat_mfl(istep)*eta(i) .GT. subsqmin) then
                    npassing = i
                else
                    EXIT
                end if
            end do
            npassing_prev = npassing
            do istep = ibmin - 1, ibeg, -1
                do i = 0, npart
                    if (1.d0 - bhat_mfl(istep)*eta(i) .GT. subsqmin) then
                        npassing = i
                    else
                        EXIT
                    end if
                end do
                if (npassing .LT. npassing_prev) then
                    ncross_l = ncross_l + 1
                    icross_l(ncross_l) = istep
                    npassing_prev = npassing
                end if
            end do
        end if
    end if

    ncross_r = 0
    if (ibmin .LT. iend) then
        istep = ibmin
        do i = 0, npart
            if (1.d0 - bhat_mfl(istep)*eta(i) .GT. subsqmin) then
                npassing = i
            else
                EXIT
            end if
        end do
        npassing_prev = npassing
        do istep = ibmin + 1, iend
            do i = 0, npart
                if (1.d0 - bhat_mfl(istep)*eta(i) .GT. subsqmin) then
                    npassing = i
                else
                    EXIT
                end if
            end do
            if (npassing .LT. npassing_prev) then
                ncross_r = ncross_r + 1
                npassing_prev = npassing
            end if
        end do
        if (ncross_r .GT. 0) then
            allocate (icross_r(ncross_r))
            ncross_r = 0
            istep = ibmin
            do i = 0, npart
                if (1.d0 - bhat_mfl(istep)*eta(i) .GT. subsqmin) then
                    npassing = i
                else
                    EXIT
                end if
            end do
            npassing_prev = npassing
            do istep = ibmin + 1, iend
                do i = 0, npart
                    if (1.d0 - bhat_mfl(istep)*eta(i) .GT. subsqmin) then
                        npassing = i
                    else
                        EXIT
                    end if
                end do
                if (npassing .LT. npassing_prev) then
                    ncross_r = ncross_r + 1
                    icross_r(ncross_r) = istep
                    npassing_prev = npassing
                end if
            end do
        end if
    end if

! place ibmin to an odd point:

    if (MOD(ibmin - ibeg, 2) .EQ. 1) then
        if (ncross_l .GT. 0 .AND. ncross_r .GT. 0) then
            if (icross_r(1) - ibmin .GT. ibmin - icross_l(1)) then
                ibmin = ibmin + 1
            else
                ibmin = ibmin - 1
            end if
        elseif (ncross_l .GT. 0) then
            ibmin = ibmin + 1
        elseif (ncross_r .GT. 0) then
            ibmin = ibmin - 1
        end if
    end if

! check the number of steps in sub-intervals for parabolic bhat term:

    if (ncross_l .GT. 0) then
        ie = icross_l(1)
        do i = 2, ncross_l
            ib = icross_l(i)
            if (ie - ib .LT. nstepmin) then
                imid = (ib + ie)/2
                x1 = phi_mfl(imid) - phi_mfl(ib)
                x2 = phi_mfl(ie) - phi_mfl(ib)
                f1 = bhat_mfl(imid) - bhat_mfl(ib)
                f2 = bhat_mfl(ie) - bhat_mfl(ib)
                bparab = ABS((f1*x2 - f2*x1)*x2/((x1 - x2)*x1*f2))
                if (bparab .GT. bparabmax) then
                    isplit = 2*MAX(NINT(0.5*float(nstepmin)/float(ie - ib)), 1)
                    phi_divide(ib + 1:ie) = isplit
                end if
            end if
            ie = ib
        end do
        ib = ibeg
        if (ie - ib .LT. nstepmin) then
            isplit = 2*MAX(NINT(0.5*float(nstepmin)/float(ie - ib)), 1)
            phi_divide(ib + 1:ie) = isplit
        end if
    end if

    if (ncross_r .GT. 0) then
        ib = icross_r(1)
        do i = 2, ncross_r
            ie = icross_r(i)
            if (ie - ib .LT. nstepmin) then
                imid = (ib + ie)/2
                x1 = phi_mfl(imid) - phi_mfl(ib)
                x2 = phi_mfl(ie) - phi_mfl(ib)
                f1 = bhat_mfl(imid) - bhat_mfl(ib)
                f2 = bhat_mfl(ie) - bhat_mfl(ib)
                bparab = ABS((f1*x2 - f2*x1)*x2/((x1 - x2)*x1*f2))
                if (bparab .GT. bparabmax) then
                    isplit = 2*MAX(NINT(0.5*float(nstepmin)/float(ie - ib)), 1)
                    phi_divide(ib + 1:ie) = isplit
                end if
            end if
            ib = ie
        end do
        ie = iend
        if (ie - ib .LT. nstepmin) then
            isplit = 2*MAX(NINT(0.5*float(nstepmin)/float(ie - ib)), 1)
            phi_divide(ib + 1:ie) = isplit
        end if
    end if

    if (MAXVAL(phi_divide) .GT. 1) RETURN

! change the integration variable phi -> sqrt(phi-phi0):

    if (stepmode) then

        allocate (phi_new(ibeg:iend), bhat_new(ibeg:iend))
        allocate (arr_real_new(ibeg:iend, nreal))
        allocate (arr_comp_new(ibeg:iend, ncomp))
        phi_new = phi_mfl
        bhat_new = bhat_mfl
        arr_real_new = arr_real
        arr_comp_new = arr_comp

        ie = ibmin
        do i = 1, ncross_l
            ib = icross_l(i)
            ht = SQRT(phi_mfl(ie) - phi_mfl(ib))/float(ie - ib)
            ht2 = ht**2
            k = ib
            do istep = ib + 1, ie - 1
                phi_new(istep) = phi_mfl(ib) + ht2*float(istep - ib)**2
                fact_pos_e(istep) = 2.d0*ht*float(istep - ib)
                delt_pos(istep) = ht
                do WHILE (phi_mfl(k) .LT. phi_new(istep))
                    k = k + 1
                end do
                intb = MAX(ibeg, MIN(iend - npoi + 1, k - npoihalf))
                inte = intb + npoi - 1

                call plagrange_coeff(npoi, nder, phi_new(istep), phi_mfl(intb:inte), coeff)

                bhat_new(istep) = SUM(coeff(0, :)*bhat_mfl(intb:inte))
                arr_real_new(istep, :) = MATMUL(coeff(0, :), arr_real(intb:inte, :))
                arr_comp_new(istep, :) = MATMUL(coeff(0, :), arr_comp(intb:inte, :))
            end do
            delt_pos(ie) = ht
            fact_pos_e(ie) = 2.d0*ht*float(ie - ib)
            fact_pos_b(ib) = 0.d0
            fact_pos_b(ib + 1:ie - 1) = fact_pos_e(ib + 1:ie - 1)
            ie = ib
        end do

        ib = ibmin
        do i = 1, ncross_r
            ie = icross_r(i)
            ht = SQRT(phi_mfl(ie) - phi_mfl(ib))/float(ie - ib)
            ht2 = ht**2
            k = ib
            do istep = ib + 1, ie - 1
                phi_new(istep) = phi_mfl(ie) - ht2*float(ie - istep)**2
                delt_pos(istep) = ht
                fact_pos_b(istep) = 2.d0*ht*float(ie - istep)
                do WHILE (phi_mfl(k) .LT. phi_new(istep))
                    k = k + 1
                end do
                intb = MAX(ibeg, MIN(iend - npoi + 1, k - npoihalf))
                inte = intb + npoi - 1

                call plagrange_coeff(npoi, nder, phi_new(istep), phi_mfl(intb:inte), coeff)

                bhat_new(istep) = SUM(coeff(0, :)*bhat_mfl(intb:inte))
                arr_real_new(istep, :) = MATMUL(coeff(0, :), arr_real(intb:inte, :))
                arr_comp_new(istep, :) = MATMUL(coeff(0, :), arr_comp(intb:inte, :))
            end do
            delt_pos(ie) = ht
            fact_pos_b(ib) = 2.d0*ht*float(ie - ib)
            fact_pos_e(ie) = 0.d0
            fact_pos_e(ib + 1:ie - 1) = fact_pos_b(ib + 1:ie - 1)
            ib = ie
        end do

        phi_mfl = phi_new
        bhat_mfl = bhat_new
        arr_real = arr_real_new
        arr_comp = arr_comp_new

        deallocate (phi_new, bhat_new, arr_real_new, arr_comp_new)

    end if

    delt_neg(ibeg:iend - 1) = delt_pos(ibeg + 1:iend)
    fact_neg_b = fact_pos_e
    fact_neg_e = fact_pos_b

    if (allocated(icross_l)) deallocate (icross_l)
    if (allocated(icross_r)) deallocate (icross_r)

end subroutine rearrange_phideps

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine fix_phiplacement_problem(ibeg, iend, npart, subsqmin, &
                                    phi_mfl, bhat_mfl, eta)

    use device_mod

    implicit none

    integer, parameter :: dp = kind(1.0d0)

    integer :: i, ibeg, iend, npart, istep, ibmin, npassing, npassing_prev
    integer :: ncross_l, ncross_r

    real(dp) :: subsqmin

    integer, dimension(1)              :: idummy
    integer, dimension(:), allocatable :: icross_l, icross_r

    real(dp), dimension(0:npart)        :: eta
    real(dp), dimension(ibeg:iend)      :: phi_mfl, bhat_mfl
    real(dp), dimension(:), allocatable :: eta_cross_l, eta_cross_r

    npassing = -1

! determine level crossings:

    idummy = MINLOC(bhat_mfl(ibeg:iend))
    ibmin = idummy(1) + ibeg - 1

    ncross_l = 0
    if (ibmin .GT. ibeg) then
        istep = ibmin
        do i = 0, npart
            if (1.d0 - bhat_mfl(istep)*eta(i) .GT. subsqmin) then
                npassing = i
            else
                EXIT
            end if
        end do
        npassing_prev = npassing
        do istep = ibmin - 1, ibeg, -1
            do i = 0, npart
                if (1.d0 - bhat_mfl(istep)*eta(i) .GT. subsqmin) then
                    npassing = i
                else
                    EXIT
                end if
            end do
            if (npassing .LT. npassing_prev) then
                ncross_l = ncross_l + 1
                npassing_prev = npassing
            end if
        end do
        if (ncross_l .GT. 0) then
            allocate (icross_l(ncross_l), eta_cross_l(ncross_l))
            ncross_l = 0
            istep = ibmin
            do i = 0, npart
                if (1.d0 - bhat_mfl(istep)*eta(i) .GT. subsqmin) then
                    npassing = i
                else
                    EXIT
                end if
            end do
            npassing_prev = npassing
            do istep = ibmin - 1, ibeg, -1
                do i = 0, npart
                    if (1.d0 - bhat_mfl(istep)*eta(i) .GT. subsqmin) then
                        npassing = i
                    else
                        EXIT
                    end if
                end do
                if (npassing .LT. npassing_prev) then
                    ncross_l = ncross_l + 1
                    icross_l(ncross_l) = istep
                    eta_cross_l(ncross_l) = eta(npassing_prev)
                    npassing_prev = npassing
                end if
            end do
            do i = 1, ncross_l
                istep = icross_l(i)
                if (ABS(bhat_mfl(istep - 1)*eta_cross_l(i) - 1.d0) .LT. &
                    ABS(bhat_mfl(istep)*eta_cross_l(i) - 1.d0)) then
                    open (111, file='phi_placement_problem.dat', position='append')
                    write (111, *) ' propagator tag = ', fieldpropagator%tag, &
                        ' step number = ', istep - 1, &
                        ' 1 / bhat = ', 1.d0/bhat_mfl(istep - 1), &
                        ' eta = ', eta_cross_l(i)
                    close (111)
                    bhat_mfl(istep - 1) = 1/eta_cross_l(i)
                elseif (ABS(bhat_mfl(istep + 1)*eta_cross_l(i) - 1.d0) .LT. &
                        ABS(bhat_mfl(istep)*eta_cross_l(i) - 1.d0)) then
                    open (111, file='phi_placement_problem.dat', position='append')
                    write (111, *) ' propagator tag = ', fieldpropagator%tag, &
                        ' step number = ', istep + 1, &
                        ' 1 / bhat = ', 1.d0/bhat_mfl(istep + 1), &
                        ' eta = ', eta_cross_l(i)
                    bhat_mfl(istep + 1) = 1/eta_cross_l(i)
                    close (111)
                end if
            end do
            deallocate (icross_l, eta_cross_l)
        end if
    end if

    ncross_r = 0
    if (ibmin .LT. iend) then
        istep = ibmin
        do i = 0, npart
            if (1.d0 - bhat_mfl(istep)*eta(i) .GT. subsqmin) then
                npassing = i
            else
                EXIT
            end if
        end do
        npassing_prev = npassing
        do istep = ibmin + 1, iend
            do i = 0, npart
                if (1.d0 - bhat_mfl(istep)*eta(i) .GT. subsqmin) then
                    npassing = i
                else
                    EXIT
                end if
            end do
            if (npassing .LT. npassing_prev) then
                ncross_r = ncross_r + 1
                npassing_prev = npassing
            end if
        end do
        if (ncross_r .GT. 0) then
            allocate (icross_r(ncross_r), eta_cross_r(ncross_r))
            ncross_r = 0
            istep = ibmin
            do i = 0, npart
                if (1.d0 - bhat_mfl(istep)*eta(i) .GT. subsqmin) then
                    npassing = i
                else
                    EXIT
                end if
            end do
            npassing_prev = npassing
            do istep = ibmin + 1, iend
                do i = 0, npart
                    if (1.d0 - bhat_mfl(istep)*eta(i) .GT. subsqmin) then
                        npassing = i
                    else
                        EXIT
                    end if
                end do
                if (npassing .LT. npassing_prev) then
                    ncross_r = ncross_r + 1
                    icross_r(ncross_r) = istep
                    eta_cross_r(ncross_r) = eta(npassing_prev)
                    npassing_prev = npassing
                end if
            end do
            do i = 1, ncross_r
                istep = icross_r(i)
                if (ABS(bhat_mfl(istep - 1)*eta_cross_r(i) - 1.d0) .LT. &
                    ABS(bhat_mfl(istep)*eta_cross_r(i) - 1.d0)) then
                    open (111, file='phi_placement_problem.dat', position='append')
                    write (111, *) ' propagator tag = ', fieldpropagator%tag, &
                        ' step number = ', istep - 1, &
                        ' 1 / bhat = ', 1.d0/bhat_mfl(istep - 1), &
                        ' eta = ', eta_cross_r(i)
                    close (111)
                    bhat_mfl(istep - 1) = 1/eta_cross_r(i)
                elseif (ABS(bhat_mfl(istep + 1)*eta_cross_r(i) - 1.d0) .LT. &
                        ABS(bhat_mfl(istep)*eta_cross_r(i) - 1.d0)) then
                    open (111, file='phi_placement_problem.dat', position='append')
                    write (111, *) ' propagator tag = ', fieldpropagator%tag, &
                        ' step number = ', istep + 1, &
                        ' 1 / bhat = ', 1.d0/bhat_mfl(istep + 1), &
                        ' eta = ', eta_cross_r(i)
                    close (111)
                    bhat_mfl(istep + 1) = 1/eta_cross_r(i)
                end if
            end do
            deallocate (icross_r, eta_cross_r)
        end if
    end if

end subroutine fix_phiplacement_problem
!! End Modifications by Andreas F. Martitsch (27.07.2015)
