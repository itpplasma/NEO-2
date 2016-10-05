SUBROUTINE compute_resonance_conditions(sw_comp_bounce_time)
  !
  ! interface to PSPLINE library (also for comparison with plagrange_mod)
  USE pspline_routines_mod 
  ! user-defined interface to GSL root-finding routines 
  USE gsl_roots_routines_mod
  ! user-defined interface to GSL integration routines 
  USE gsl_integration_routines_mod
  ! get parameter "Lc/lc"
  USE collisionality_mod, ONLY : conl_over_mfp
  ! get numerical constants (pi,euler,...)
  use nrtype
  ! contains type definitions (e.g., fieldpropagator_struct)
  USE magnetics_mod 
  ! contains global structures defining the device
  USE device_mod 
  ! normalization for "hat"-quantities
  USE partpa_mod,  ONLY : bmod0 
  ! flux surface label,...
  USE mag_interface_mod, ONLY: boozer_s,boozer_theta_beg,boozer_phi_beg 
  USE neo_magfie_mod, ONLY: boozer_iota,boozer_curr_pol_hat,&
       boozer_curr_tor_hat,boozer_psi_pr_hat,boozer_curr_pol_hat_s,&
       boozer_curr_tor_hat_s
  ! direct computation of B-field components
  USE magfie_mod, ONLY : magfie 
  ! or Lagrange interpolation of the B-field components
  USE plagrange_mod
  !
  IMPLICIT NONE
  !
  ! switch on/off computation
  LOGICAL, INTENT(in) ::  sw_comp_bounce_time
  !
  ! information about the B-field (from fieldpropagator)
  ! local handle for global variable
  TYPE(fieldperiod_struct), POINTER     :: fieldperiod_loc
  TYPE(fieldpropagator_struct), POINTER :: fieldprop_loc 
  ! number of discrete phi_mfl values
  INTEGER :: ind_phi_mfl, num_phi_mfl 
  ! variable along the field line
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: phi_mfl 
  ! normalized value of |B|
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: bhat_mfl 
  ! minimum / maximum value of |B| from fieldripple
  DOUBLE PRECISION :: b_max_l, b_max_r, b_min 
  ! find phi-value corresponding to minimum / maximum value of |B|
  DOUBLE PRECISION :: phi_mfl_max_l_tmp, phi_mfl_max_r_tmp, phi_mfl_min_tmp
  DOUBLE PRECISION :: phi_mfl_l, phi_mfl_r, phi_mfl_min
  ! major radius R0 (from device structure)
  DOUBLE PRECISION :: R0
  !
  ! grid for the velocity space coordinates ($\eta$,$z$)
  INTEGER :: ind_z, num_z, ind_eta, num_eta
  DOUBLE PRECISION :: z_val, z_delta, z_max, z_min
  DOUBLE PRECISION :: eta_val, eta_delta, eta_max, eta_min
  DOUBLE PRECISION :: tau_bounce_val
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: tau_bounce_arr,&
       tau_bounce_arr_adapt
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: res_cond_arr,&
       res_cond_arr_adapt
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: res_eta_vec, res_z_vec
  !
  ! print information for debugging
  LOGICAL :: comp_talk
  !
  ! computation of the bounce-average / superbanana plateau resonances
  DOUBLE PRECISION :: phi_max, phi_min
  DOUBLE PRECISION :: phi_bound_acc_fac
  DOUBLE PRECISION :: phi_step_im1, phi_step_i
  DOUBLE PRECISION :: bhat_step_im1, bhat_step_i
  DOUBLE PRECISION :: term_step_im1, term_step_i
  DOUBLE PRECISION :: tau_bounce_phi_step
  DOUBLE PRECISION :: res_cond_phi_step,res_cond_phi_step2
  DOUBLE PRECISION :: denomjac, prefac
  DOUBLE PRECISION :: C_factor
  !
  ! initialize plagrange_interp routine
  INTEGER, PARAMETER :: nlagrange = 5
  !
  ! initialize Q(uadrature) A(daptive) G(eneral integrand)
  ! integration PROCEDUR with known singular points
  DOUBLE PRECISION, PARAMETER :: epsabs=1.0d-10, epsrel=1.0d-10
  INTEGER, PARAMETER :: ilong=SELECTED_INT_KIND(15)
  ! 2 singular points (boundaries of the integration interval)
  INTEGER(ilong) :: siz_pts=2
  ! storage array for the result and the singular points
  DOUBLE PRECISION, DIMENSION(2) :: res_fint1d, singular_pts
  ! initialize CQUAD doubly-adaptive integration
  DOUBLE PRECISION, DIMENSION(3) :: res2_fint1d
  !
!!$  ! Example for using routine magfie
!!$  DOUBLE PRECISION :: phi_val,bmoda,sqrtg
!!$  DOUBLE PRECISION, DIMENSION(3) :: x,bder,hcovar,hctrvr,hcurl
!!$  x(1)=boozer_s
!!$  x(2)=phi_val
!!$  x(3)=boozer_theta_beg+(x(2)-boozer_phi_beg)*boozer_iota
!!$  CALL magfie( x, bmoda, sqrtg, bder, hcovar, hctrvr, hcurl )
  !
  comp_talk=.TRUE.
  !comp_talk=.FALSE.
  !
  PRINT *,'--------------------------------------------'
  PRINT *,'Compute resonance conditions'
  !
  ! store fieldpropagator/fieldperiod locally
  fieldperiod_loc => fieldpropagator%parent
  fieldprop_loc => fieldpropagator
  num_phi_mfl=SIZE(fieldprop_loc%coords%x2,1)
  !
  ! major radius R0 (from device-structure)
  R0=device%r0
  !PRINT *,'R0: ',R0
  !
  ! factor for the magnetic rotation
  C_factor=3.0d0
  !
  ! pre-factor for bounce time, denominator of Jacobian
  denomjac=boozer_iota*boozer_curr_tor_hat+boozer_curr_pol_hat
  prefac=(conl_over_mfp/(PI*R0))*denomjac
  !
  ! mesh along the field line (refined after ripple_solver)
  ! --> simple Simpson Integrator should be sufficient
  IF (ALLOCATED(phi_mfl)) DEALLOCATE(phi_mfl)
  ALLOCATE(phi_mfl(num_phi_mfl))
  phi_mfl=fieldprop_loc%coords%x2
  IF (ALLOCATED(bhat_mfl)) DEALLOCATE(bhat_mfl)
  ALLOCATE(bhat_mfl(num_phi_mfl))
  bhat_mfl=fieldprop_loc%mdata%bhat
  !
  ! detect minimum / maximum value (from fieldripple)
  b_max_l = fieldprop_loc%ch_act%b_max_l
  b_max_r = fieldprop_loc%ch_act%b_max_r
  b_min   = fieldprop_loc%ch_act%b_min
  PRINT *,'b_max_l, b_min, b_max_r:       ',b_max_l,b_min,b_max_r
  phi_mfl_l = fieldprop_loc%phi_l
  phi_mfl_r = fieldprop_loc%phi_r
  phi_mfl_min = fieldprop_loc%phi_min
  PRINT *,'phi_max_l, phi_min, phi_max_r: ',phi_mfl_l,phi_mfl_min,phi_mfl_r
!!$  ! check the maximum / minimum value (using numerical routines)
!!$  ! --> maximum values are identically the same as the 
!!$  ! --> left/right phi-value of the propagator
!!$  PRINT *,'phi_min_tmp,b_min_tmp:          ',phi_mfl_min,fieldprop_loc%b_min
!!$  phi_mfl_max_l_tmp = fzero1d_bisec(find_phi_at_Bval,&
!!$       b_max_l,phi_mfl_l,phi_mfl_min)
!!$  PRINT *,'phi_max_l_tmp,b_max_l_tmp:     ',phi_mfl_max_l_tmp,b_max_l
!!$  phi_mfl_max_r_tmp = fzero1d_bisec(find_phi_at_Bval,&
!!$       b_max_r,phi_mfl_min,phi_mfl_r)
!!$  PRINT *,'phi_max_r_tmp,b_max_r_tmp:     ',phi_mfl_max_r_tmp,b_max_r
  !
  ! generate mesh for ($\eta$,$z$)-space
  num_z=101
  z_max=10.0d0
  z_min=1.0d-1
  z_delta=(z_max-z_min)/DBLE(num_z-1)
  !
  num_eta=101
  eta_max=(1.0d0-1.0d-8)/MINVAL(bhat_mfl)
  !eta_max=(1.0d0+1.0d-7)/MAXVAL(bhat_mfl)
  eta_min=(1.0d0+1.0d-8)/MAXVAL(bhat_mfl)
  eta_delta=(eta_max-eta_min)/DBLE(num_eta-1)
  !
  ! compute the bounce time and the resonance function
  ! for the superbanana plateau as a function of ($z$,$\eta$)
  ! using a simple Simpson integration formula
  PRINT *,'Compute bounce time using a Simpsons formula'
  IF (ALLOCATED(tau_bounce_arr)) DEALLOCATE(tau_bounce_arr)
  ALLOCATE(tau_bounce_arr(num_z,num_eta))
  tau_bounce_arr=0.0d0
  IF (ALLOCATED(res_cond_arr)) DEALLOCATE(res_cond_arr)
  ALLOCATE(res_cond_arr(num_z,num_eta))
  res_cond_arr=0.0d0
  eta_val=eta_min
  all_eta: DO ind_eta=1,num_eta
     phi_max=fzero1d_bisec(find_phiMaxMin_for_etaVal,&
          eta_val,phi_mfl_min,phi_mfl_r)
     !PRINT *,'eta_val,phi_max: ',eta_val,phi_max
     phi_min=fzero1d_bisec(find_phiMaxMin_for_etaVal,&
          eta_val,phi_mfl_l,phi_mfl_min)
     !PRINT *,'eta_val,phi_min: ',eta_val,phi_min
     tau_bounce_phi_step=0.0d0
     res_cond_phi_step=0.0d0
     all_phi: DO ind_phi_mfl=1,num_phi_mfl-1
        IF (phi_mfl(ind_phi_mfl) .GT. phi_min) THEN
           IF (phi_mfl(ind_phi_mfl+1) .LT. phi_max) THEN
              phi_step_im1=phi_mfl(ind_phi_mfl)
              bhat_step_im1=bhat_mfl(ind_phi_mfl)
              phi_step_i=phi_mfl(ind_phi_mfl+1)
              bhat_step_i=bhat_mfl(ind_phi_mfl+1)
              term_step_im1=1.0d0/(bhat_step_im1*SQRT(1-eta_val*bhat_step_im1))
              !PRINT *,phi_step_im1,1.0d0-eta_val*bhat_step_im1,&
              !     1/bhat_step_im1,term_step_im1
              term_step_i=1.0d0/(bhat_step_i*SQRT(1-eta_val*bhat_step_i))
              tau_bounce_phi_step=tau_bounce_phi_step+&
                   0.5d0*(phi_step_i-phi_step_im1)*(term_step_i+term_step_im1)
              term_step_im1=&
                   calc_aB1_bounce_av(phi_step_im1,eta_val)+&
                   calc_aB2_bounce_av(phi_step_im1,eta_val)
              term_step_i=&
                   calc_aB1_bounce_av(phi_step_i,eta_val)+&
                   calc_aB2_bounce_av(phi_step_i,eta_val)
              res_cond_phi_step=res_cond_phi_step+&
                   0.5d0*(phi_step_i-phi_step_im1)*(term_step_i+term_step_im1)
           END IF
        END IF
     END DO all_phi
     !PRINT *,phi_step_i,1.0d0-eta_val*bhat_step_i,1/bhat_step_i,term_step_i
     !PRINT *,'tau_bounce_phi_step        : ',tau_bounce_phi_step
     z_val=z_min
     all_z: DO ind_z=1,num_z
        tau_bounce_arr(ind_z,ind_eta)=&
             (prefac/SQRT(z_val))*tau_bounce_phi_step
        res_cond_arr(ind_z,ind_eta)=&
             (res_cond_phi_step/tau_bounce_phi_step)+1.0d0/(C_factor*z_val)
        z_val=z_val+z_delta
     END DO all_z
     eta_val=eta_val+eta_delta
  END DO all_eta
  !
  ! write file containing the computation of the bounce time
  OPEN(unit=28041443,file='bounce_time_simpson_formula.dat')
  WRITE(unit=28041443,fmt='(A,1X,D16.8,1X,D16.8)') '%',z_min,z_max
  WRITE(unit=28041443,fmt='(A,1X,D16.8,1X,D16.8)') '%',eta_min,eta_max
  DO ind_z=1,num_z
     WRITE(unit=28041443,fmt=*) &
          (tau_bounce_arr(ind_z,ind_eta),ind_eta=1,num_eta)
  END DO
  CLOSE(unit=28041443)
  !
  ! write file containing the computation of the resonance function
  OPEN(unit=28041710,file='res_fun_simpson_formula.dat')
  WRITE(unit=28041710,fmt='(A,1X,D16.8,1X,D16.8)') '%',z_min,z_max
  WRITE(unit=28041710,fmt='(A,1X,D16.8,1X,D16.8)') '%',eta_min,eta_max
  DO ind_z=1,num_z
     WRITE(unit=28041710,fmt=*) &
          (res_cond_arr(ind_z,ind_eta),ind_eta=1,num_eta)
  END DO
  CLOSE(unit=28041710)
  !
  ! compute the bounce time and the resonance function
  ! for the superbanana plateau as a function of ($z$,$\eta$)
  ! using an advanced numerical integration scheme
  PRINT *,'Compute bounce time using an advanced numerical integration scheme'
  IF (ALLOCATED(tau_bounce_arr_adapt)) DEALLOCATE(tau_bounce_arr_adapt)
  ALLOCATE(tau_bounce_arr_adapt(num_z,num_eta))
  tau_bounce_arr_adapt=0.0d0
  IF (ALLOCATED(res_cond_arr_adapt)) DEALLOCATE(res_cond_arr_adapt)
  ALLOCATE(res_cond_arr_adapt(num_z,num_eta))
  res_cond_arr_adapt=0.0d0
  eta_val=eta_min
  all_eta_2: DO ind_eta=1,num_eta
     phi_max=fzero1d_bisec(find_phiMaxMin_for_etaVal,&
          eta_val,phi_mfl_min,phi_mfl_r)
     !PRINT *,'eta_val,phi_max: ',eta_val,phi_max
     phi_min=fzero1d_bisec(find_phiMaxMin_for_etaVal,&
          eta_val,phi_mfl_l,phi_mfl_min)
     !PRINT *,'eta_val,phi_min: ',eta_val,phi_min
     ! use qagp-integrator
     tau_bounce_phi_step=0.0d0
     singular_pts=(/phi_min,phi_max/)
     res_fint1d = fint1d_qagp(bounce_time_fun,eta_val,&
          singular_pts,siz_pts,epsabs,epsrel)
     tau_bounce_phi_step=res_fint1d(1)
     res_fint1d = fint1d_qagp(calc_aB1_bounce_av,eta_val,&
          singular_pts,siz_pts,epsabs,epsrel)
     res_cond_phi_step=res_fint1d(1)
     res_fint1d = fint1d_qagp(calc_aB2_bounce_av,eta_val,&
          singular_pts,siz_pts,epsabs,epsrel)
     res_cond_phi_step2=res_fint1d(1)
     !PRINT *,'tau_bounce_phi_step (qaqp) : ',tau_bounce_phi_step
     ! use cquad-integrator (oscillating behavior near singularity)
     !res2_fint1d = fint1d_cquad(bounce_time_fun,eta_val,&
     !  phi_min,phi_max,epsabs,epsrel)
     !tau_bounce_phi_step=res2_fint1d(1)
     !PRINT *,'tau_bounce_phi_step (cquad): ',tau_bounce_phi_step
     z_val=z_min
     all_z_2: DO ind_z=1,num_z
        tau_bounce_arr_adapt(ind_z,ind_eta)=&
             (prefac/SQRT(z_val))*tau_bounce_phi_step
        res_cond_arr_adapt(ind_z,ind_eta)=&
             ((res_cond_phi_step+res_cond_phi_step2)/tau_bounce_phi_step)+&
             1.0d0/(C_factor*z_val)
        z_val=z_val+z_delta
     END DO all_z_2
     eta_val=eta_val+eta_delta
  END DO all_eta_2
  !
  ! write file containing the computation of the bounce time
  OPEN(unit=28041459,file='bounce_time_adv_integration.dat')
  WRITE(unit=28041459,fmt='(A,1X,D16.8,1X,D16.8)') '%',z_min,z_max
  WRITE(unit=28041459,fmt='(A,1X,D16.8,1X,D16.8)') '%',eta_min,eta_max
  DO ind_z=1,num_z
     WRITE(unit=28041459,fmt=*) &
          (tau_bounce_arr_adapt(ind_z,ind_eta),ind_eta=1,num_eta)
  END DO
  CLOSE(unit=28041459)
  !
  ! write file containing the computation of the function
  OPEN(unit=28041722,file='res_fun_adv_integration.dat')
  WRITE(unit=28041722,fmt='(A,1X,D16.8,1X,D16.8)') '%',z_min,z_max
  WRITE(unit=28041722,fmt='(A,1X,D16.8,1X,D16.8)') '%',eta_min,eta_max
  DO ind_z=1,num_z
     WRITE(unit=28041722,fmt=*) &
          (res_cond_arr_adapt(ind_z,ind_eta),ind_eta=1,num_eta)
  END DO
  CLOSE(unit=28041722)
  !
!!$  ! computation of the resonance condition
!!$  IF (ALLOCATED(res_eta_vec)) DEALLOCATE(res_eta_vec)
!!$  ALLOCATE(res_eta_vec(num_z))
!!$  IF (ALLOCATED(res_z_vec)) DEALLOCATE(res_z_vec)
!!$  ALLOCATE(res_z_vec(num_z))
!!$  z_val=z_min
!!$  all_z_3: DO ind_z=1,num_z
!!$     eta_val=fzero1d_bisec(resonance_fun_adapt,z_val,eta_min,eta_max)
!!$     res_eta_vec(ind_z)=eta_val
!!$     res_z_vec(ind_z)=z_val
!!$     z_val=z_val+z_delta
!!$  END DO all_z_3
!!$  !
!!$  ! write file containing the computation of the resonance condition
!!$  OPEN(unit=28041812,file='res_cond_adv_integration.dat')
!!$  DO ind_z=1,num_z
!!$     WRITE(unit=28041812,fmt=*) res_z_vec(ind_z),res_eta_vec(ind_z)
!!$  END DO
!!$  CLOSE(unit=28041812)
  !
  ! deallocate storage arrays / nullify pointer
  IF (ALLOCATED(tau_bounce_arr)) DEALLOCATE(tau_bounce_arr)
  IF (ALLOCATED(tau_bounce_arr_adapt)) DEALLOCATE(tau_bounce_arr_adapt)
  IF (ALLOCATED(res_cond_arr)) DEALLOCATE(res_cond_arr)
  IF (ALLOCATED(res_cond_arr_adapt)) DEALLOCATE(res_cond_arr_adapt)
  IF (ALLOCATED(res_eta_vec)) DEALLOCATE(res_eta_vec)
  IF (ALLOCATED(res_z_vec)) DEALLOCATE(res_z_vec)
  IF (ALLOCATED(phi_mfl)) DEALLOCATE(phi_mfl)
  IF (ALLOCATED(bhat_mfl)) DEALLOCATE(bhat_mfl)
  NULLIFY(fieldperiod_loc)
  NULLIFY(fieldprop_loc)
  !
  ! stop / leave routine
  PRINT *,'--------------------------------------------'
  STOP "Exit compute_resonance_conditions"
  !RETURN
  !
CONTAINS
  !
  FUNCTION find_phi_at_Bval(phi,Bval)
    ! use splines to interpolate data to ensure a consistent 
    ! computation of B-field values
    ! -> using magfie, the results of the computations differ 
    ! -> from the splined values (difference at the last digits)
    DOUBLE PRECISION :: find_phi_at_Bval
    DOUBLE PRECISION :: phi, Bval
    find_phi_at_Bval=cub1d_pspline_0(phi_mfl,bhat_mfl,phi,1)-Bval
  END FUNCTION find_phi_at_Bval
  !
  FUNCTION find_phiMaxMin_for_etaVal(phi,etaVal)
    ! use splines to interpolate data to ensure a consistent 
    ! computation of B-field values
    ! -> using magfie, the results of the computations differ 
    ! -> from the splined values (difference at the last digits)
    DOUBLE PRECISION :: find_phiMaxMin_for_etaVal
    DOUBLE PRECISION :: phi, etaVal
    DOUBLE PRECISION :: geodcu_loc,h_phi_loc,dlogbdphi_loc,&
         bcovar_s_hat_loc,dlogbds_loc,dbcovar_s_hat_dphi_loc,&
         bhat_loc,x1_loc,x3_loc
    ! compute the result either using the PSPLINE-library
    !find_phiMaxMin_for_etaVal=&
    !     1.0d0-cub1d_pspline_0(phi_mfl,bhat_mfl,phi,1)*etaVal
    ! or "our" Lagrange interpolation
    CALL plagrange_interp(fieldperiod_loc,phi,nlagrange,&
         x1_loc,x3_loc,bhat_loc,geodcu_loc,h_phi_loc,&
         dlogbdphi_loc,dbcovar_s_hat_dphi_loc,&
         bcovar_s_hat_loc,dlogbds_loc)
    find_phiMaxMin_for_etaVal=1.0d0-bhat_loc*etaVal
  END FUNCTION find_phiMaxMin_for_etaVal
  !
  FUNCTION bounce_time_fun(phi,etaVal)
    DOUBLE PRECISION :: bounce_time_fun
    DOUBLE PRECISION :: phi, etaVal
    DOUBLE PRECISION :: geodcu_loc,h_phi_loc,dlogbdphi_loc,&
         bcovar_s_hat_loc,dlogbds_loc,dbcovar_s_hat_dphi_loc,&
         bhat_loc,x1_loc,x3_loc
    ! compute the result either using the PSPLINE-library
    !find_phiMaxMin_for_etaVal=&
    !     1.0d0-cub1d_pspline_0(phi_mfl,bhat_mfl,phi,1)*etaVal
    ! or "our" Lagrange interpolation
    CALL plagrange_interp(fieldperiod_loc,phi,nlagrange,&
         x1_loc,x3_loc,bhat_loc,geodcu_loc,h_phi_loc,&
         dlogbdphi_loc,dbcovar_s_hat_dphi_loc,&
         bcovar_s_hat_loc,dlogbds_loc)
    bounce_time_fun=1.0d0/(bhat_loc*SQRT(1.0d0-bhat_loc*etaVal))
  END FUNCTION bounce_time_fun
  !
  FUNCTION calc_aB1_bounce_av(phi,etaVal)
    DOUBLE PRECISION :: calc_aB1_bounce_av
    DOUBLE PRECISION :: phi, etaVal
    DOUBLE PRECISION :: geodcu_loc,h_phi_loc,dlogbdphi_loc,&
         bcovar_s_hat_loc,dlogbds_loc,dbcovar_s_hat_dphi_loc,&
         bhat_loc,x1_loc,x3_loc
    ! compute the result either using the PSPLINE-library
    !find_phiMaxMin_for_etaVal=&
    !     1.0d0-cub1d_pspline_0(phi_mfl,bhat_mfl,phi,1)*etaVal
    ! or "our" Lagrange interpolation
    CALL plagrange_interp(fieldperiod_loc,phi,nlagrange,&
         x1_loc,x3_loc,bhat_loc,geodcu_loc,h_phi_loc,&
         dlogbdphi_loc,dbcovar_s_hat_dphi_loc,&
         bcovar_s_hat_loc,dlogbds_loc)
    calc_aB1_bounce_av=&
         (bcovar_s_hat_loc*dlogbdphi_loc/denomjac-dlogbds_loc)/boozer_iota
    calc_aB1_bounce_av=calc_aB1_bounce_av*bounce_time_fun(phi,etaVal)
  END FUNCTION calc_aB1_bounce_av
  !
  FUNCTION calc_aB2_bounce_av(phi,etaVal)
    DOUBLE PRECISION :: calc_aB2_bounce_av
    DOUBLE PRECISION :: phi, etaVal
    DOUBLE PRECISION :: geodcu_loc,h_phi_loc,dlogbdphi_loc,&
         bcovar_s_hat_loc,dlogbds_loc,dbcovar_s_hat_dphi_loc,&
         bhat_loc,x1_loc,x3_loc
    ! compute the result either using the PSPLINE-library
    !find_phiMaxMin_for_etaVal=&
    !     1.0d0-cub1d_pspline_0(phi_mfl,bhat_mfl,phi,1)*etaVal
    ! or "our" Lagrange interpolation
    CALL plagrange_interp(fieldperiod_loc,phi,nlagrange,&
         x1_loc,x3_loc,bhat_loc,geodcu_loc,h_phi_loc,&
         dlogbdphi_loc,dbcovar_s_hat_dphi_loc,&
         bcovar_s_hat_loc,dlogbds_loc)
    calc_aB2_bounce_av=&
         (2.0d0/denomjac)*(boozer_curr_tor_hat_s+boozer_curr_pol_hat_s/&
         boozer_iota-dbcovar_s_hat_dphi_loc/boozer_iota)
    calc_aB2_bounce_av=&
         calc_aB2_bounce_av*SQRT(1.0d0-bhat_loc*etaVal)/bhat_loc
    calc_aB2_bounce_av=calc_aB2_bounce_av+&
         calc_aB1_bounce_av(phi,etaVal)*(1.0d0-bhat_loc*etaVal)
  END FUNCTION calc_aB2_bounce_av
  !
  FUNCTION resonance_fun_adapt(etaVal,z)
    DOUBLE PRECISION :: resonance_fun_adapt
    DOUBLE PRECISION :: etaVal,z
    DOUBLE PRECISION :: phi_min, phi_max
    DOUBLE PRECISION :: tau_bounce, aB1, aB2
    DOUBLE PRECISION, DIMENSION(2) :: result_tmp, sing_pts
    ! find the boundaries of the bounce period
    phi_max=fzero1d_bisec(find_phiMaxMin_for_etaVal,&
         etaVal,phi_mfl_min,phi_mfl_r)
    phi_min=fzero1d_bisec(find_phiMaxMin_for_etaVal,&
         etaVal,phi_mfl_l,phi_mfl_min)
    ! compute the bounce average
    sing_pts=(/phi_min,phi_max/)
    !
    result_tmp = fint1d_qagp(bounce_time_fun,etaVal,&
         sing_pts,siz_pts,epsabs,epsrel)
    tau_bounce=result_tmp(1)
    !
    result_tmp = fint1d_qagp(calc_aB1_bounce_av,etaVal,&
         sing_pts,siz_pts,epsabs,epsrel)
    aB1=result_tmp(1)
    !
    result_tmp = fint1d_qagp(calc_aB2_bounce_av,etaVal,&
         sing_pts,siz_pts,epsabs,epsrel)
    aB2=result_tmp(1)
    !
    resonance_fun_adapt=((aB1+aB2)/tau_bounce)+1.0d0/(C_factor*z)
    !
  END FUNCTION resonance_fun_adapt
  !
END SUBROUTINE compute_resonance_conditions
