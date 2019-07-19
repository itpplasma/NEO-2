MODULE neo_spline_data
  ! Splines along s
  USE neo_precision
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: a_rmnc,b_rmnc,c_rmnc,d_rmnc 
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: a_zmnc,b_zmnc,c_zmnc,d_zmnc 
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: a_lmnc,b_lmnc,c_lmnc,d_lmnc 
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: a_bmnc,b_bmnc,c_bmnc,d_bmnc
  !! Modifications by Andreas F. Martitsch (06.08.2014)
  ! Additional data from Boozer files without Stellarator symmetry
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: a_rmns,b_rmns,c_rmns,d_rmns 
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: a_zmns,b_zmns,c_zmns,d_zmns 
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: a_lmns,b_lmns,c_lmns,d_lmns 
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: a_bmns,b_bmns,c_bmns,d_bmns
  !! End Modifications by Andreas F. Martitsch (06.08.2014)

  REAL(kind=dp), DIMENSION(:),   ALLOCATABLE :: a_iota,b_iota
  REAL(kind=dp), DIMENSION(:),   ALLOCATABLE :: c_iota,d_iota
  REAL(kind=dp), DIMENSION(:),   ALLOCATABLE :: a_pprime,b_pprime
  REAL(kind=dp), DIMENSION(:),   ALLOCATABLE :: c_pprime,d_pprime
  REAL(kind=dp), DIMENSION(:),   ALLOCATABLE :: a_sqrtg00,b_sqrtg00
  REAL(kind=dp), DIMENSION(:),   ALLOCATABLE :: c_sqrtg00,d_sqrtg00
  REAL(kind=dp), DIMENSION(:),   ALLOCATABLE :: a_curr_tor,b_curr_tor
  REAL(kind=dp), DIMENSION(:),   ALLOCATABLE :: c_curr_tor,d_curr_tor
  REAL(kind=dp), DIMENSION(:),   ALLOCATABLE :: a_curr_pol,b_curr_pol
  REAL(kind=dp), DIMENSION(:),   ALLOCATABLE :: c_curr_pol,d_curr_pol
 
  REAL(kind=dp), DIMENSION(:),   ALLOCATABLE :: r_m, r_mhalf
  INTEGER(I4B),  DIMENSION(:),   ALLOCATABLE :: sp_index

END MODULE neo_spline_data

MODULE neo_spline_b00
  ! spline for b00
  USE neo_precision
  REAL(kind=dp), DIMENSION(:),   ALLOCATABLE :: a_b00, b_b00, c_b00, d_b00
  
END MODULE neo_spline_b00

MODULE neo_work
! Working parameters
  USE neo_precision

  REAL(kind=dp) ::   theta_start
  REAL(kind=dp) ::   theta_end
  REAL(kind=dp) ::   theta_int

  REAL(kind=dp) ::   phi_start
  REAL(kind=dp) ::   phi_end
  REAL(kind=dp) ::   phi_int

! REAL(kind=dp),    DIMENSION(:,:,:),   ALLOCATABLE :: cosval,sinval
  REAL(kind=dp),    DIMENSION(:,:),     ALLOCATABLE :: cosmth,sinmth
  REAL(kind=dp),    DIMENSION(:,:),     ALLOCATABLE :: cosnph,sinnph
  REAL(kind=dp),    DIMENSION(:),       ALLOCATABLE :: theta_arr,phi_arr
  REAL(kind=dp),    DIMENSION(:,:),     ALLOCATABLE :: r,z,l,b
  REAL(kind=dp),    DIMENSION(:,:),     ALLOCATABLE :: r_tb,z_tb,p_tb,b_tb
  REAL(kind=dp),    DIMENSION(:,:),     ALLOCATABLE :: r_pb,z_pb,p_pb,b_pb
  REAL(kind=dp),    DIMENSION(:,:),     ALLOCATABLE :: gtbtb,gpbpb,gtbpb
  REAL(kind=dp),    DIMENSION(:,:),     ALLOCATABLE :: isqrg,sqrg11,kg,pard
  REAL(kind=dp),    DIMENSION(:,:),     ALLOCATABLE :: bqtphi
  REAL(kind=dp),    DIMENSION(:,:),     ALLOCATABLE :: r_nabpsi
  REAL(kind=dp),    DIMENSION(:,:),     ALLOCATABLE :: psi_r, psi_z
END MODULE neo_work

MODULE neo_exchange
! Working parameters
  USE neo_precision

  REAL(kind=dp)              :: b_min, b_max
  INTEGER                    :: nper
  REAL(kind=dp)              :: rt0, rt0_g 
  REAL(kind=dp)              :: bmref, bmref_g, bmref_a
  INTEGER                    :: nstep_per, nstep_min, nstep_max
  INTEGER                    :: write_integrate
  INTEGER                    :: write_diagnostic
  INTEGER                    :: write_cur_inte
  INTEGER                    :: write_pla_inte
  REAL(kind=dp)              :: acc_req
  INTEGER                    :: no_bins
  INTEGER                    :: psi_ind
  INTEGER                    :: calc_nstep_max
  REAL(kind=dp)              :: theta_bmin, phi_bmin
  REAL(kind=dp)              :: theta_bmax, phi_bmax
  REAL(kind=dp),    DIMENSION(:),     ALLOCATABLE :: iota
  REAL(kind=dp),    DIMENSION(:),     ALLOCATABLE :: curr_pol
  REAL(kind=dp),    DIMENSION(:),     ALLOCATABLE :: curr_tor
  REAL(kind=dp)              :: fac
  INTEGER                    :: calc_cur, calc_eps, calc_pla, calc_van
  INTEGER                    :: hit_rat, nfp_rat, nfl_rat
  REAL(kind=dp)              :: delta_theta_rat 
  REAL(kind=dp)              :: delta_cur_fac  
  INTEGER                    :: cutoff_cur_int
  INTEGER                    :: write_cur_disp 
END MODULE neo_exchange

MODULE neo_output
  USE neo_precision
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE ::  epspar
  REAL(kind=dp)                            ::  epstot,ctrone,ctrtot
  REAL(kind=dp)                            ::  epstothat
  REAL(kind=dp)                            ::  bareph,barept,drdpsi
  REAL(kind=dp)                            ::  yps
  INTEGER                                  ::  nintfp
  INTEGER                                  ::  ierr  
  REAL(kind=dp)                            ::  lambda_b
  REAL(kind=dp)                            ::  lambda_b1, lambda_b2
  REAL(kind=dp)                            ::  lambda_ps1, lambda_ps2
  REAL(kind=dp)                            ::  lambda_del
  REAL(kind=dp)                            ::  avnabpsi,rfint
  REAL(kind=dp)                            ::  avb2, f_c, f_p
  REAL(kind=dp)                            ::  lambda_pla  
  REAL(kind=dp)                            ::  delta_cur_max, typ_cur_len
END MODULE neo_output

MODULE neo_eval_switch
  ! Evaluation
  USE neo_precision
  INTEGER, DIMENSION(6)    ::   eval_switch
END MODULE neo_eval_switch

MODULE neo_spline
! Spline arrays
  USE neo_precision
  INTEGER, PARAMETER       ::   mt = 1
  INTEGER, PARAMETER       ::   mp = 1
  INTEGER                  ::   theta_ind, phi_ind
  INTEGER                  ::   ierr

  REAL(kind=dp) ::  theta_d, phi_d 

! Spline array for modb
  REAL(kind=dp),    DIMENSION(:,:,:,:), ALLOCATABLE :: b_spl 
! Spline array for geodesic curviture
  REAL(kind=dp),    DIMENSION(:,:,:,:), ALLOCATABLE :: k_spl
! Spline array for sqrg11
  REAL(kind=dp),    DIMENSION(:,:,:,:), ALLOCATABLE :: g_spl
! Spline array for parallel derivative
  REAL(kind=dp),    DIMENSION(:,:,:,:), ALLOCATABLE :: p_spl
! Spline array for quasi-toroidal phi component of b
  REAL(kind=dp),    DIMENSION(:,:,:,:), ALLOCATABLE :: q_spl
! Spline array for r_nabpsi
  REAL(kind=dp),    DIMENSION(:,:,:,:), ALLOCATABLE :: r_spl
END MODULE neo_spline

MODULE neo_control
! Control parameters from input file
  USE neo_precision
  CHARACTER(20)                      :: in_file
  INTEGER                            :: theta_n
  INTEGER                            :: phi_n
  INTEGER                            :: s_ind_in
  INTEGER                            :: write_progress
  INTEGER                            :: write_output_files
  INTEGER                            :: calc_fourier
  INTEGER                            :: spline_test
  INTEGER                            :: max_m_mode, max_n_mode
  INTEGER                            :: lab_swi, inp_swi, ref_swi, eout_swi
  INTEGER                            :: chk_swi
  INTEGER                            :: fluxs_interp
  INTEGER                            :: s_num
  REAL(kind=dp)                      :: s_start, s_end
  INTEGER                            :: g11_swi
  INTEGER                            :: eval_mode
  INTEGER                            :: no_fluxs, no_fluxs_s
  INTEGER, DIMENSION(:), ALLOCATABLE :: fluxs_arr
END MODULE neo_control

MODULE neo_units
! Units and Formats
  USE neo_precision
  INTEGER, PARAMETER ::   r_u1   = 3
  INTEGER, PARAMETER ::   r_u2   = 4
  INTEGER, PARAMETER ::   r_us   = 5
  INTEGER, PARAMETER ::   r_u23  = 23
  INTEGER, PARAMETER ::   r_ua   = 21
  INTEGER, PARAMETER ::   w_us   = 6
  INTEGER, PARAMETER ::   w_u1   = 7
  INTEGER, PARAMETER ::   w_u2   = 8
  INTEGER, PARAMETER ::   w_u3   = 9
  INTEGER, PARAMETER ::   w_u4   = 10
  INTEGER, PARAMETER ::   w_u5   = 11
  INTEGER, PARAMETER ::   w_u6   = 12
  INTEGER, PARAMETER ::   w_u7   = 13
  INTEGER, PARAMETER ::   w_u8   = 14
  INTEGER, PARAMETER ::   w_u9   = 15
  INTEGER, PARAMETER ::   w_u10  = 16
  INTEGER, PARAMETER ::   w_u11  = 17
  INTEGER, PARAMETER ::   w_u12  = 18
  INTEGER, PARAMETER ::   w_u13  = 19
  INTEGER, PARAMETER ::   w_u14  = 20
  INTEGER, PARAMETER ::   w_u15  = 21
  INTEGER, PARAMETER ::   w_u16  = 22
  INTEGER, PARAMETER ::   w_u17  = 23

  INTEGER            ::   w_u6_open
  CHARACTER(20),PARAMETER :: format220="(500d18.5)"

  CHARACTER(30)                      :: base_file
  CHARACTER(30)                      :: out_file
  CHARACTER(30)                      :: chk_file
  CHARACTER(30)                      :: epslog_file
  CHARACTER(30)                      :: epscon_file
  CHARACTER(30)                      :: epsdia_file
  CHARACTER(30)                      :: curcon_file
  CHARACTER(30)                      :: curint_file
  CHARACTER(30)                      :: curdis_file
  CHARACTER(30)                      :: epsadd_file
  CHARACTER(30)                      :: cur_file
  CHARACTER(30)                      :: pla_file
  CHARACTER(30)                      :: sbc_file

END MODULE neo_units

MODULE sizey_bo
  USE neo_precision
! Definition for rk4d_bo also used in main routine neo
  INTEGER            ::  npart
  INTEGER            ::  multra
  INTEGER            ::  ndim
  INTEGER, PARAMETER ::  npq = 7
END MODULE sizey_bo

MODULE partpa_bo
  USE neo_precision
! Exchange between flint_bo and rhs_bo
  USE sizey_bo
  INTEGER                                  :: ipmax
  INTEGER,       DIMENSION(:), ALLOCATABLE :: isw,ipa,icount
  REAL(kind=dp)                            :: pard0,bmod0
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: eta
END MODULE partpa_bo

MODULE sizey_cur
  USE neo_precision
! Definition for rk4d_bo also used in main routine neo
  INTEGER            ::  npart_cur
  INTEGER            ::  ndim_cur
  INTEGER, PARAMETER ::  npq_cur = 11
  INTEGER            ::  alpha_cur
END MODULE sizey_cur

MODULE partpa_cur
  USE neo_precision
! Exchange between flint_cur and rhs_cur
  USE sizey_cur
  REAL(kind=dp)                            :: bmod0
  REAL(kind=dp)                            :: gamma_cur
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: y_part, yfac, sqyfac
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: yfac_xin, yfac_xid
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: delta_cur
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: k_fac1, k_fac2
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: contrif
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: fcontrif
  INTEGER                                  :: write_curint
END MODULE partpa_cur

MODULE sizey_pla
  USE neo_precision
  ! Definition for rk4d_pla 
  INTEGER            ::  npart_pla
  INTEGER            ::  ndim_pla
  INTEGER, PARAMETER ::  npq_pla = 3
  REAL(kind=dp)      ::  lamup_pla
  REAL(kind=dp)      ::  lambda_alpha
  REAL(kind=dp)      ::  nufac_pla  
END MODULE sizey_pla

MODULE partpa_pla
  USE neo_precision
  ! Exchange between flint_pla and rhs_pla
  USE sizey_pla
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: p_fac
END MODULE partpa_pla

MODULE neo_van_exchange
  USE neo_precision
  REAL(kind=dp)                            :: rho_fac, jperp_fac
  REAL(kind=dp)                            :: theta_fac, phi_fac
END MODULE neo_van_exchange

MODULE neo_conversion
  USE neo_precision
  REAL(kind=dp), PARAMETER                 :: mc_o_e        = 5.6856793d-8 ! cgs
  REAL(kind=dp), PARAMETER                 :: mc_o_2e       = 2.8428397d-8 ! cgs
  REAL(kind=dp), PARAMETER                 :: e_o_mc        = 1.7588048d7  ! cgs
  REAL(kind=dp), PARAMETER                 :: b_convfac     = 1.0d4  ! SI to cgs
  REAL(kind=dp), PARAMETER                 :: i_convfac     = 1.0d6  ! SI to cgs
  REAL(kind=dp), PARAMETER                 :: sqg11_convfac = 1.0d6  ! SI to cgs
  REAL(kind=dp), PARAMETER                 :: len_convfac   = 1.0d2  ! SI to cgs
  REAL(kind=dp), PARAMETER                 :: te_to_vte     = 4.19d7 ! ev to cm/s
  
END MODULE neo_conversion
