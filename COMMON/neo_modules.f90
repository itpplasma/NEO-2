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
