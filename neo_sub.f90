MODULE neo_sub_mod

CONTAINS

SUBROUTINE neo_init(npsi)
! Initialization Routine
! **********************************************************************
! Modules
! **********************************************************************
  USE neo_precision
  USE neo_input
  USE neo_work
  USE neo_exchange
  USE neo_units
  USE neo_parameters
  USE neo_control
  USE neo_spline
! **********************************************************************
! Local Definitions
! **********************************************************************
  IMPLICIT NONE
  INTEGER, INTENT(out)       :: npsi
  INTEGER                    :: imn
! **********************************************************************
! Read input from data file and allocate necessary arrays
! **********************************************************************
  IF (write_progress .NE. 0) WRITE (w_us,*) 'before neo_filenames'
  CALL neo_filenames
  IF (write_progress .NE. 0) WRITE (w_us,*) 'after  neo_filenames'
! **********************************************************************
! Read input from data file and allocate necessary arrays
! **********************************************************************
  IF (write_progress .NE. 0) WRITE (w_us,*) 'before neo_read'
  CALL neo_read
  IF (write_progress .NE. 0) WRITE (w_us,*) 'after  neo_read'
! **********************************************************************
  npsi = ns
! **********************************************************************
! Prepare b00 for radius correction
! **********************************************************************
  IF (write_progress .NE. 0) WRITE (w_us,*) 'before neo_prep_b00'
  CALL neo_prep_b00
  IF (write_progress .NE. 0) WRITE (w_us,*) 'after  neo_prep_b00'
! **********************************************************************
! Allocate and prepare necessary arrays
! **********************************************************************
  IF (write_progress .NE. 0) WRITE (w_us,*) 'before neo_prep'
  CALL neo_prep
  IF (write_progress .NE. 0) WRITE (w_us,*) 'after  neo_prep'
! **********************************************************************
! Allocate and prepare spline along s
! **********************************************************************
  IF (fluxs_interp .NE. 0) THEN
     IF (write_progress .NE. 0) WRITE (w_us,*) 'before neo_init_spline'
     CALL neo_init_spline()
     IF (write_progress .NE. 0) WRITE (w_us,*) 'after  neo_init_spline'
  END IF
! **********************************************************************
! Calculation of rt0 and bmref (innermost flux surface)
! might be changed later
! **********************************************************************
  rt0=0.0_dp
  bmref=0.0_dp
  DO imn=1,mnmax
     IF(ixm(imn).EQ.0 .AND. ixn(imn).EQ.0) THEN
        rt0 = rmnc(1,imn)
        bmref = bmnc(1,imn)

        rt0_g = rt0
        bmref_g = bmref
     ENDIF
  ENDDO
  IF(rt0.EQ.0.0_dp .OR. bmref.EQ.0.0_dp) THEN
    WRITE (w_us,*) ' NEO_INIT: Fatal problem setting rt0 or bmref'
    STOP
  ENDIF
!
  nper = nfp
! **********************************************************************
  w_u6_open = 0
! **********************************************************************
  RETURN
END SUBROUTINE neo_init
! **********************************************************************

! **********************************************************************
SUBROUTINE neo_init_spline()
! Initialization for splines along s
! **********************************************************************
! Modules
! **********************************************************************
  USE neo_precision
  USE neo_spline_data
  USE neo_input
  USE neo_exchange
  USE inter_interfaces, ONLY: splinecof3_hi_driv, splinecof3, tf
!  Test
!  USE inter_interfaces, ONLY: splinecof3_hi_driv, splinecof3, tf,      &
!       splint_horner3, tfp, tfpp, tfppp 
!  Test End
  IMPLICIT NONE
  INTEGER      :: i
  INTEGER(I4B) :: sw1, sw2
  REAL(dp)     :: m0, c1, cn
  REAL(dp), DIMENSION(:), ALLOCATABLE :: lambda
  INTEGER(I4B) :: m
  INTEGER,  PARAMETER                 :: m_max_sp = 12
!
! Testing
!
!!$  INTEGER(I4B) :: swd
!!$  INTEGER      :: k
!!$  REAL(dp)     :: f_es, f_bmnc, dummy
!
! Testing End
!
  ALLOCATE ( a_rmnc(ns,mnmax), b_rmnc(ns,mnmax) )
  ALLOCATE ( c_rmnc(ns,mnmax), d_rmnc(ns,mnmax) )
  ALLOCATE ( a_zmnc(ns,mnmax), b_zmnc(ns,mnmax) )
  ALLOCATE ( c_zmnc(ns,mnmax), d_zmnc(ns,mnmax) )
  ALLOCATE ( a_lmnc(ns,mnmax), b_lmnc(ns,mnmax) )
  ALLOCATE ( c_lmnc(ns,mnmax), d_lmnc(ns,mnmax) )
  ALLOCATE ( a_bmnc(ns,mnmax), b_bmnc(ns,mnmax) ) 
  ALLOCATE ( c_bmnc(ns,mnmax), d_bmnc(ns,mnmax) ) 

  ALLOCATE ( a_iota(ns), b_iota(ns) )
  ALLOCATE ( c_iota(ns), d_iota(ns) )
  ALLOCATE ( a_pprime(ns), b_pprime(ns) )
  ALLOCATE ( c_pprime(ns), d_pprime(ns) )
  ALLOCATE ( a_sqrtg00(ns), b_sqrtg00(ns) )
  ALLOCATE ( c_sqrtg00(ns), d_sqrtg00(ns) )
  ALLOCATE ( a_curr_tor(ns), b_curr_tor(ns) )
  ALLOCATE ( c_curr_tor(ns), d_curr_tor(ns) )
  ALLOCATE ( a_curr_pol(ns), b_curr_pol(ns) )
  ALLOCATE ( c_curr_pol(ns), d_curr_pol(ns) )
 
  ALLOCATE ( r_m(mnmax), r_mhalf(mnmax) )
  ALLOCATE ( sp_index(ns) )

  DO i = 1,mnmax
     m = ixm(i)
     IF (m .LE. m_max_sp) THEN
        r_m(i)     = DBLE(m)
     ELSE
        IF (MODULO(m,2) .EQ. 1) THEN
           r_m(i)     = DBLE(m_max_sp+1)
        ELSE
           r_m(i)     = DBLE(m_max_sp)
        END IF
     END IF
     r_mhalf(i) = r_m(i) / 2._dp
  END DO
  sp_index = (/ (i, i=1,ns) /) 

  ! 1-d splines of 2-d arrays
  CALL splinecof3_hi_driv(es, rmnc, r_mhalf,                         &
       a_rmnc, b_rmnc, c_rmnc, d_rmnc, sp_index, tf)
  CALL splinecof3_hi_driv(es, zmnc, r_mhalf,                         &
       a_zmnc, b_zmnc, c_zmnc, d_zmnc, sp_index, tf)
  CALL splinecof3_hi_driv(es, lmnc, r_mhalf,                         &
       a_lmnc, b_lmnc, c_lmnc, d_lmnc, sp_index, tf)
  CALL splinecof3_hi_driv(es, bmnc, r_mhalf,                         &
       a_bmnc, b_bmnc, c_bmnc, d_bmnc, sp_index, tf)
  !
  ! Testing
  !
!!$  swd = 0 ! no derivatives
!!$  k = 3
!!$  f_es = es(k)
!!$  DO i=1,mnmax
!!$     m0 = r_mhalf(i)
!!$     CALL splint_horner3(es,a_bmnc(:,i),b_bmnc(:,i),            &
!!$          c_bmnc(:,i),d_bmnc(:,i),swd,m0,                       &
!!$          f_es,tf,tfp,tfpp,tfppp,                               &
!!$          f_bmnc,dummy,dummy,dummy)
!!$     PRINT *, ixm(i),ixn(i),m0,bmnc(k,i),f_bmnc 
!!$  END DO  
  !
  ! End Testing
  !

  ! boundary types (natural spline)
  sw1 = 2
  sw2 = 4
  ! input for test function for spline
  m0  = 0.0_dp
! boundary condition for spline
  c1 = 0.0_dp  
  cn = 0.0_dp
! we use no smoothing for spline
  ALLOCATE ( lambda(ns) )
  lambda = 1.0D0
! 1-d splines of 1-d arrays
  CALL splinecof3(es, iota, c1, cn, lambda, sp_index, sw1, sw2, &
       a_iota, b_iota, c_iota, d_iota, m0, tf)
  CALL splinecof3(es, pprime, c1, cn, lambda, sp_index, sw1, sw2, &
       a_pprime, b_pprime, c_pprime, d_pprime, m0, tf)
  CALL splinecof3(es, sqrtg00, c1, cn, lambda, sp_index, sw1, sw2, &
       a_sqrtg00, b_sqrtg00, c_sqrtg00, d_sqrtg00, m0, tf)
  CALL splinecof3(es, curr_tor, c1, cn, lambda, sp_index, sw1, sw2, &
       a_curr_tor, b_curr_tor, c_curr_tor, d_curr_tor, m0, tf)
  CALL splinecof3(es, curr_pol, c1, cn, lambda, sp_index, sw1, sw2, &
       a_curr_pol, b_curr_pol, c_curr_pol, d_curr_pol, m0, tf)
!
  DEALLOCATE( lambda )

END SUBROUTINE neo_init_spline

! **********************************************************************
SUBROUTINE neo_init_s(psi,dpsi)
! Initialization for Specific Magnetic Surface
! **********************************************************************
! Modules
! **********************************************************************
  USE neo_precision
  USE neo_input
  USE neo_work
  USE neo_exchange
  USE neo_units
  USE neo_parameters
  USE neo_control
  USE neo_spline
  USE neo_eval_switch
  USE neo_actual_fluxs
  USE neo_actual_spectra
! **********************************************************************
! Local Definitions
! **********************************************************************
  IMPLICIT NONE

  REAL(kind=dp),                INTENT(out)     :: psi, dpsi 
  INTEGER,       DIMENSION(2)                   :: b_minpos, b_maxpos 

  REAL(kind=dp), PARAMETER                      :: eps_newt = 1.0d-10
  INTEGER                                       :: iter, error
  INTEGER                                       :: imn
  INTEGER,       PARAMETER                      :: iterma_newt = 100
  INTEGER                                       :: ff_theta_n, ff_phi_n
  REAL(kind=dp)                                 :: gval_bmin
  REAL(kind=dp)                                 :: kval_bmin,pval_bmin
  REAL(kind=dp)                                 :: gval_bmax
  REAL(kind=dp)                                 :: kval_bmax,pval_bmax
  REAL(kind=dp)                                 :: qval_bmin,qval_bmax
  REAL(kind=dp)                                 :: rval_bmin,rval_bmax
  LOGICAL                                       :: ff_exist
  CHARACTER(30)                                 :: flux_file
  CHARACTER(3)                                  :: flux_numc
!  REAL(kind=dp)  :: tht, pht, f, g, dfdx, dfdy, dgdx, dgdy
! **********************************************************************
  eval_switch = (/1,0,0,0,0,0/)
!
! **********************************************************************
! Artificial Cutting of modes (wrong place has to be changed)
! **********************************************************************
!  DO imn=1,mnmax
!     IF(ixm(imn).GE.5 .OR. ixn(imn).LE. -11*nper                       &
!        .OR. ixn(imn).GE. 12*nper                                      &
!        .OR. abs(bmnc(psi_ind,imn)) .LE. 1.d-4*bmref_a ) THEN
!        bmnc(psi_ind,imn) = 0.0_dp
!     ENDIF
!  ENDDO
!
! **********************************************************************
! Reading or writing of Data
! **********************************************************************
  ff_exist = .FALSE.
  IF (calc_fourier .EQ. 0 .AND. fluxs_interp .EQ. 0) THEN
     WRITE(flux_numc,'(i3)') psi_ind
     flux_file = TRIM(in_file)//'.'//TRIM(ADJUSTL(flux_numc))//'.fsn'
     INQUIRE(FILE=flux_file,EXIST=ff_exist)
  END IF
  
  ! Fourier should be done
  IF ((calc_fourier .EQ. 1) .OR.                                       &
       (calc_fourier .EQ. 0 .AND. .NOT. ff_exist) .OR.                 &
       (fluxs_interp .NE. 0) .OR.                                      &
       (eval_mode .NE. 0) ) THEN
     through_fourier = 1
     ! get spectra and other quantities on flux surface
     IF (write_progress .NE. 0) WRITE (w_us,*) 'before neo_init_fluxsurface'
     CALL neo_init_fluxsurface()
     IF (write_progress .NE. 0) WRITE (w_us,*) 'after  neo_init_fluxsurface'
     ! **************************************************************
     ! Calculate Fourier sums and derived quantities
     ! **************************************************************
     IF (write_progress .NE. 0) WRITE (w_us,*) 'before neo_fourier'
     CALL neo_fourier
     IF (write_progress .NE. 0) WRITE (w_us,*) 'after  neo_fourier'
     ! *****************************************************************
     ! Calculation of bmref_a (average on flux surface)
     ! *****************************************************************
     bmref_a=0.0_dp
     DO imn=1,mnmax
        IF(ixm(imn).EQ.0 .AND. ixn(imn).EQ.0) THEN
           bmref_a = s_bmnc(imn)
        ENDIF
     ENDDO
     IF(bmref_a .EQ. 0.0_dp) THEN
        WRITE (w_us,*) ' NEO_INIT: Fatal problem setting bmref_a'
        STOP
     ENDIF
     ! writing of files to avoid neo_fourier
     IF (calc_fourier .EQ. 0) THEN
        IF (write_progress .NE. 0) WRITE (w_us,*) '  Writing: ',flux_file
        OPEN(unit=w_u10,file=flux_file,status='replace',               &
             form='unformatted',action='write')
        WRITE (w_u10) theta_n, phi_n
        WRITE (w_u10) b,sqrg11,kg,pard,bqtphi,r_nabpsi
        WRITE (w_u10) s_iota,s_pprime,s_sqrtg00,s_curr_tor,s_curr_pol
        WRITE (w_u10) s_es,s_b00,s_b00_s
        WRITE (w_u10) bmref_a
        CLOSE(unit=w_u10)
     END IF
  ELSE ! read data from pre-computed files
     through_fourier = 0
     IF (write_progress .NE. 0) WRITE (w_us,*) '  Reading: ',flux_file
     OPEN(unit=w_u10,file=flux_file,status='old',                   &
             form='unformatted',action='read')
     READ  (w_u10) ff_theta_n, ff_phi_n
     IF (ff_theta_n .NE. theta_n .OR. ff_phi_n .NE. phi_n) THEN
        WRITE (w_us,*) 'FATAL: Wrong contents in File ',flux_file
        WRITE (w_us,*) 'Please delete this file'
        STOP
     END IF
     READ  (w_u10) b,sqrg11,kg,pard,bqtphi,r_nabpsi
     READ  (w_u10) s_iota,s_pprime,s_sqrtg00,s_curr_tor,s_curr_pol
     READ  (w_u10) s_es,s_b00,s_b00_s
     READ  (w_u10) bmref_a
     CLOSE(unit=w_u10)
  END IF
! **********************************************************************
! Initilaze spline arrays (2d-periodic)
! **********************************************************************
  ! IF (write_progress .NE. 0) WRITE (w_us,*) 'before neo_spline2d'
  CALL neo_spline2d
  ! IF (write_progress .NE. 0) WRITE (w_us,*) 'after  neo_spline2d'
! **********************************************************************
! Calculate absolute minimum and maximum of b and its location (theta, phi)
! **********************************************************************
  b_minpos   = MINLOC(b)
  b_min      = b(b_minpos(1),b_minpos(2))
  theta_bmin = theta_arr(b_minpos(1))
  phi_bmin   = phi_arr(b_minpos(2))
  
  b_maxpos   = MAXLOC(b)
  b_max      = b(b_maxpos(1),b_maxpos(2))
  theta_bmax = theta_arr(b_maxpos(1))
  phi_bmax   = phi_arr(b_maxpos(2))

  ! IF (write_progress .NE. 0) WRITE (w_us,*) 'before neo_zeros2d (min)'
  CALL neo_zeros2d(theta_bmin, phi_bmin, eps_newt, iterma_newt, iter, error)
  CALL neo_eval(theta_bmin,phi_bmin,b_min,gval_bmin,kval_bmin,&
       pval_bmin,qval_bmin,rval_bmin)
  ! IF (write_progress .NE. 0) WRITE (w_us,*) 'after  neo_zeros2d'

  ! IF (write_progress .NE. 0) WRITE (w_us,*) 'before neo_zeros2d (max)'
  CALL neo_zeros2d(theta_bmax, phi_bmax, eps_newt, iterma_newt, iter, error)
  CALL neo_eval(theta_bmax,phi_bmax,b_max,gval_bmax,kval_bmax,&
       pval_bmax,qval_bmax,rval_bmax)
  ! IF (write_progress .NE. 0) WRITE (w_us,*) 'after  neo_zeros2d'
! **********************************************************************
! Calculation of dpsi (only on flussurface
! **********************************************************************
  psi = s_es * psi_pr
  IF (psi_ind < 1000) THEN
     IF (psi_ind .EQ. 1) THEN
        dpsi = psi
     ELSE
        dpsi = psi - es(psi_ind-1)*psi_pr
     ENDIF
  ELSE
     dpsi = 0.0_dp
  END IF
! **********************************************************************
! Set bmref to the absolute maximum of b on flux surface
! This is absolutely necessary for internal routines
! Rescaling is done at the end of the main program neo.f90
! **********************************************************************
  bmref = b_max
!
  IF (write_output_files .NE. 0) THEN
     IF (write_progress .NE. 0) WRITE (w_us,*) 'write plot files'
     CALL neo_plot_files
  ENDIF
!
! **********************************************************************
!  Deallocation of unnecessary arrays
! **********************************************************************
  IF (through_fourier .EQ. 1) THEN
     DEALLOCATE (r,z,l)
     DEALLOCATE (r_tb,z_tb,p_tb,b_tb)
     DEALLOCATE (r_pb,z_pb,p_pb,b_pb)
     DEALLOCATE (gtbtb,gpbpb,gtbpb,isqrg)
     DEALLOCATE (psi_r,psi_z)
  END IF
! **********************************************************************
  RETURN
END SUBROUTINE neo_init_s
! **********************************************************************

! **********************************************************************
! Write optional output for Plotting
! **********************************************************************
SUBROUTINE neo_plot_files
  ! **********************************************************************
  ! Modules
  ! **********************************************************************
  USE neo_precision
  USE neo_support
  USE neo_units
  USE neo_input
  USE neo_exchange
  USE neo_control
  USE neo_work
  USE neo_actual_fluxs

  IMPLICIT NONE
  CHARACTER(len=50)   :: plot_file
  CHARACTER(len=50)   :: plot_file1
  CHARACTER(len=10)   :: c_ns,c_phi_n,c_theta_n
  CHARACTER(len=20)   :: f_ns,f_phi_n,f_theta_n,f,f_f
  INTEGER             :: uw = 101
  INTEGER             :: irecl

  CALL add_extension(base_file,psi_ind,plot_file1)
  CALL add_extension(plot_file1,'pdat',plot_file)
  CALL unit_check(uw)

  WRITE(c_ns,*) ns
  WRITE(c_phi_n,*) phi_n
  WRITE(c_theta_n,*) theta_n

  WRITE(f_f,*) 'e15.5E3'
  WRITE(f,*) '(1x,',TRIM(ADJUSTL(f_f)),')'
  WRITE(f_ns,*) '(1x,',TRIM(ADJUSTL(c_ns)),TRIM(ADJUSTL(f_f)),')'
  WRITE(f_phi_n,*) '(1x,',TRIM(ADJUSTL(c_phi_n)),TRIM(ADJUSTL(f_f)),')'
  WRITE(f_theta_n,*) '(1x,',TRIM(ADJUSTL(c_theta_n)),TRIM(ADJUSTL(f_f)),')'
  

  INQUIRE(IOLENGTH=irecl ) b
  OPEN(unit=uw, file=plot_file, status='replace', action='write', recl=irecl)

  WRITE (uw,*) '#v version = 3.3'
  WRITE (uw,*) '#e equilibrium = ',base_file
  WRITE (uw,*) '#s psi_ind = ',psi_ind
  WRITE (uw,*) '#s ns = ',ns
  WRITE (uw,*) '#s mnmax = ',mnmax
  WRITE (uw,*) '#s nfp = ',nfp
  WRITE (uw,*) '#s m_max = ',m_max
  WRITE (uw,*) '#s n_max = ',n_max
  WRITE (uw,*) '#s mnmax = ',mnmax
  WRITE (uw,*) '#s flux = ',flux
  WRITE (uw,*) '#s psi_pr = ',psi_pr
  WRITE (uw,*) '#s s = ',s_es
  WRITE (uw,*) '#s psi = ',s_es*psi_pr
  WRITE (uw,*) '#s b_min = ',b_min
  WRITE (uw,*) '#s theta_bmin = ',theta_bmin
  WRITE (uw,*) '#s phi_bmin = ',phi_bmin
  WRITE (uw,*) '#s b_max = ',b_max
  WRITE (uw,*) '#s theta_bmax = ',theta_bmax
  WRITE (uw,*) '#s phi_bmax = ',phi_bmax

  WRITE (uw,*) '#s phi_n = ',phi_n
  WRITE (uw,*) '#s theta_n = ',theta_n

  WRITE (uw,*) '#g es'
  WRITE (uw,f_ns) es
  WRITE (uw,*) '#g iota'
  WRITE (uw,f_ns) iota
  WRITE (uw,*) '#g curr_pol'
  WRITE (uw,f_ns) curr_pol
  WRITE (uw,*) '#g curr_tor'
  WRITE (uw,f_ns) curr_tor
  WRITE (uw,*) '#g pprime'
  WRITE (uw,f_ns) pprime
  WRITE (uw,*) '#g sqrtg00'
  WRITE (uw,f_ns) sqrtg00

  WRITE (uw,*) '#g phi'
  WRITE (uw,f_phi_n) phi_arr
  WRITE (uw,*) '#g theta'
  WRITE (uw,f_theta_n) theta_arr

  WRITE (uw,*) '#a bfield'
  WRITE (uw,*) '#c (theta_n, phi_n)'
  WRITE (uw,*) '#k modb'
  WRITE (uw,f_phi_n) TRANSPOSE(b)
 
  WRITE (uw,*) '#a geom'
  WRITE (uw,*) '#c (theta_n, phi_n)'
  WRITE (uw,*) '#c r [m], z [m], d_phi = (phib-phi)*nper/twopi'
  WRITE (uw,*) '#k r'
  WRITE (uw,f_phi_n) TRANSPOSE(r)
  WRITE (uw,*) '#k z'
  WRITE (uw,f_phi_n) TRANSPOSE(z)
  WRITE (uw,*) '#k d_phi'
  WRITE (uw,f_phi_n) TRANSPOSE(l)

  WRITE (uw,*) '#a derived'
  WRITE (uw,*) '#c (theta_n, phi_n)'
  WRITE (uw,*) '#c derived quantities'
  WRITE (uw,*) '#k isqrg'
  WRITE (uw,f_phi_n) TRANSPOSE(isqrg)
  WRITE (uw,*) '#k sqrg11'
  WRITE (uw,f_phi_n) TRANSPOSE(sqrg11)
  WRITE (uw,*) '#k kg'
  WRITE (uw,f_phi_n) TRANSPOSE(kg)
  WRITE (uw,*) '#k pard'
  WRITE (uw,f_phi_n) TRANSPOSE(pard)


!!$  WRITE(w_u1,*)ixm(j),ixn(j)
!!$  WRITE(w_u1,*) rmnc(i,j)
!!$  WRITE(w_u1,*) zmnc(i,j)
!!$  WRITE(w_u1,*) lmnc(i,j)
!!$  WRITE(w_u1,*) bmnc(i,j)
!!$  WRITE(w_u1,*) theta_arr(j)
!!$  WRITE(w_u1,*) phi_arr(k)
!!$ 
!!$
!!$  DO i=1,theta_n
!!$     DO j=1,phi_n
!!$        WRITE(w_u1,*)  r(i,j)
!!$        WRITE(w_u1,*)  z(i,j)
!!$        WRITE(w_u1,*)  l(i,j)
!!$        WRITE(w_u1,*)  b(i,j)
!!$        WRITE(w_u1,*)  isqrg(i,j)
!!$        WRITE(w_u1,*)  sqrg11(i,j)
!!$        WRITE(w_u1,*)  kg(i,j)
!!$        WRITE(w_u1,*)  pard(i,j)
!!$     END DO

  CLOSE(unit=uw)

END SUBROUTINE neo_plot_files

! **********************************************************************
! Initialization for Specific Magnetic Surface
! **********************************************************************
SUBROUTINE neo_init_fluxsurface
  ! **********************************************************************
  ! Modules
  ! **********************************************************************
  USE neo_precision
  USE neo_input
  USE neo_control
  USE neo_exchange
  USE neo_actual_fluxs
  USE neo_actual_spectra
  USE neo_spline_data
  USE inter_interfaces, ONLY: splint_horner3, tf, tfp, tfpp, tfppp 
! **********************************************************************
! Local Definitions
! **********************************************************************
  INTEGER       :: i
  INTEGER(I4B)  :: swd
  REAL(dp)      :: m0
  REAL(dp)      :: yp, ypp, yppp   ! dummies for derivatives
!
! direct from data on fluxsurface
  IF (fluxs_interp .EQ. 0) THEN
     s_rmnc(:)  = rmnc(psi_ind,:)
     s_zmnc(:)  = zmnc(psi_ind,:)
     s_lmnc(:)  = lmnc(psi_ind,:)
     s_bmnc(:)  = bmnc(psi_ind,:)
     s_curr_pol = curr_pol(psi_ind)
     s_curr_tor = curr_tor(psi_ind)
     s_pprime   = pprime(psi_ind)
     s_sqrtg00  = sqrtg00(psi_ind)
     s_iota     = iota(psi_ind)
     CALL neo_get_b00
     !s_b00      = b00(psi_ind)
  ELSE ! from splines
     ! avaluation of 2-d arrays
     swd = 0 ! no derivatives
     DO i=1,mnmax
        m0 = r_mhalf(i)
        CALL splint_horner3(es,a_rmnc(:,i),b_rmnc(:,i),            &
             c_rmnc(:,i),d_rmnc(:,i),swd,m0,                       &
             s_es,tf,tfp,tfpp,tfppp,                               &
             s_rmnc(i),yp,ypp,yppp)
        CALL splint_horner3(es,a_zmnc(:,i),b_zmnc(:,i),            &
             c_zmnc(:,i),d_zmnc(:,i),swd,m0,                       &
             s_es,tf,tfp,tfpp,tfppp,                               &
             s_zmnc(i),yp,ypp,yppp)
        CALL splint_horner3(es,a_lmnc(:,i),b_lmnc(:,i),            &
             c_lmnc(:,i),d_lmnc(:,i),swd,m0,                       &
             s_es,tf,tfp,tfpp,tfppp,                               &
             s_lmnc(i),yp,ypp,yppp)
        CALL splint_horner3(es,a_bmnc(:,i),b_bmnc(:,i),            &
             c_bmnc(:,i),d_bmnc(:,i),swd,m0,                       &
             s_es,tf,tfp,tfpp,tfppp,                               &
             s_bmnc(i),yp,ypp,yppp)
     END DO
     ! evaluation of 1-d arrays
     swd = 0
     m0  = 0.0_dp
     CALL splint_horner3(es,a_iota,b_iota,c_iota,d_iota,swd,m0, &
          s_es,tf,tfp,tfpp,tfppp,                               &
          s_iota,yp,ypp,yppp)
     CALL splint_horner3(es,a_pprime,b_pprime,c_pprime,d_pprime,swd,m0, &
          s_es,tf,tfp,tfpp,tfppp,                               &
          s_pprime,yp,ypp,yppp)
     CALL splint_horner3(es,a_sqrtg00,b_sqrtg00,c_sqrtg00,d_sqrtg00,swd,m0, &
          s_es,tf,tfp,tfpp,tfppp,                               &
          s_sqrtg00,yp,ypp,yppp)
     CALL splint_horner3(es,a_curr_tor,b_curr_tor,c_curr_tor,d_curr_tor,swd,m0, &
          s_es,tf,tfp,tfpp,tfppp,                               &
          s_curr_tor,yp,ypp,yppp)
     CALL splint_horner3(es,a_curr_pol,b_curr_pol,c_curr_pol,d_curr_pol,swd,m0, &
          s_es,tf,tfp,tfpp,tfppp,                               &
          s_curr_pol,yp,ypp,yppp)    
     CALL neo_get_b00
  END IF
! 
END SUBROUTINE neo_init_fluxsurface

SUBROUTINE neo_read_control
! Read Control File
!***********************************************************************
! Modules
!***********************************************************************
  USE neo_units
  USE neo_control
  USE neo_input
  USE neo_exchange
  USE sizey_bo
  USE sizey_cur
  USE sizey_pla
  USE neo_van
!***********************************************************************
! Local definitions
!***********************************************************************
  IMPLICIT NONE
  CHARACTER(1)          :: dummy
  INTEGER               :: i,n,stat
  INTEGER, DIMENSION(3) :: iarr
!***********************************************************************
! Open input-unit and read data
!***********************************************************************
  OPEN(unit=r_u1,file='neo.in',status='old',form='formatted')
  READ (r_u1,*) dummy
  READ (r_u1,*) dummy
  READ (r_u1,*) dummy
  READ (r_u1,*) in_file

  READ (r_u1,*) fluxs_interp  
  READ (r_u1,*) no_fluxs
  IF (fluxs_interp .EQ. 0) THEN
     IF (no_fluxs .EQ. 0) THEN
        READ (r_u1,*) dummy
     ELSE IF (no_fluxs .LT. 0) THEN
        READ (r_u1,*) iarr
        no_fluxs = (iarr(2)-iarr(1))/iarr(3) + 1
        ALLOCATE ( fluxs_arr(no_fluxs) )
        n = 0
        DO i = iarr(1),iarr(2),iarr(3)
           n = n + 1
           fluxs_arr(n) = i
        END DO
     ELSE
        ALLOCATE ( fluxs_arr(no_fluxs) )
        READ (r_u1,*,iostat=stat) fluxs_arr
        IF (stat .NE. 0) THEN
           WRITE(w_us,*) 'FATAL: Not enough flux surfaces in the input file!'
           STOP
        END IF
     END IF
     no_fluxs_s = no_fluxs
     READ (r_u1,*) dummy
     READ (r_u1,*) dummy
     READ (r_u1,*) dummy
  ELSE
     READ (r_u1,*) dummy
     READ (r_u1,*) s_start
     READ (r_u1,*) s_end
     READ (r_u1,*) s_num
     ALLOCATE ( fluxs_arr(s_num) )
     n = 1000
     DO i = 1,s_num
        n = n + 1
        fluxs_arr(i) = n 
     END DO
     no_fluxs = s_num
  END IF
  no_fluxs_s = no_fluxs

  READ (r_u1,*) theta_n
  READ (r_u1,*) phi_n
  READ (r_u1,*) max_m_mode
  READ (r_u1,*) max_n_mode
  READ (r_u1,*) calc_eps
  READ (r_u1,*) npart 
  READ (r_u1,*) multra
  READ (r_u1,*) acc_req
  READ (r_u1,*) no_bins
  READ (r_u1,*) nstep_per
  READ (r_u1,*) nstep_min
  READ (r_u1,*) nstep_max
  READ (r_u1,*) calc_nstep_max
  READ (r_u1,*) eout_swi
  READ (r_u1,*) lab_swi
  READ (r_u1,*) inp_swi
  READ (r_u1,*) g11_swi
  READ (r_u1,*) ref_swi
  READ (r_u1,*) eval_mode
  READ (r_u1,*) write_progress
  READ (r_u1,*) write_output_files
  READ (r_u1,*) spline_test
  READ (r_u1,*) write_integrate
  READ (r_u1,*) write_diagnostic
  READ (r_u1,*) calc_fourier
  READ (r_u1,*) chk_swi
  READ (r_u1,*) dummy
  READ (r_u1,*) dummy
  READ (r_u1,*) dummy
  READ (r_u1,*) calc_cur
  READ (r_u1,*) npart_cur
  READ (r_u1,*) delta_cur_fac
  READ (r_u1,*) cutoff_cur_int
  READ (r_u1,*) alpha_cur
  READ (r_u1,*) write_cur_inte
  READ (r_u1,*) write_cur_disp
  READ (r_u1,*) dummy
  READ (r_u1,*) dummy
  READ (r_u1,*) dummy
  READ (r_u1,*) calc_pla
  READ (r_u1,*) npart_pla
  READ (r_u1,*) lamup_pla
  READ (r_u1,*) lambda_alpha
  READ (r_u1,*) nufac_pla
  READ (r_u1,*) write_pla_inte
  READ (r_u1,*) dummy
  READ (r_u1,*) dummy
  READ (r_u1,*) dummy
  READ (r_u1,*) calc_van

  READ (r_u1,*) no_minima
  IF (no_minima .EQ. 0) THEN
     READ (r_u1,*) dummy
  ELSE IF (no_minima .LT. 0) THEN
     READ (r_u1,*) iarr
     no_minima = (iarr(2)-iarr(1))/iarr(3) + 1
     ALLOCATE ( li_minima(no_minima) )
     n = 0
     DO i = iarr(1),iarr(2),iarr(3)
        n = n + 1
        li_minima(n) = i
     END DO
  ELSE
     ALLOCATE ( li_minima(no_minima) )
     READ (r_u1,*,iostat=stat) li_minima
     IF (stat .NE. 0) THEN
        WRITE(w_us,*) 'FATAL: You did not specify enough minima!'
        STOP
     END IF
  END IF

  READ (r_u1,*) v_phi0
  READ (r_u1,*) v_theta0
  READ (r_u1,*) v_nper
  READ (r_u1,*) v_steps
  READ (r_u1,*) bmin_tol
  READ (r_u1,*) v_num_mm
  READ (r_u1,*) no_gamma
  READ (r_u1,*) lambda_fac
  READ (r_u1,*) temp_e  
  READ (r_u1,*) tau_num
  READ (r_u1,*) gamma_eps
  READ (r_u1,*) phi_eps
  READ (r_u1,*) tau_max_iter

!
  CLOSE (unit=r_u1)
! **********************************************************************
  RETURN

END SUBROUTINE neo_read_control
! **********************************************************************

SUBROUTINE neo_read
! Read Boozer Files
!***********************************************************************
! Modules
!***********************************************************************
  USE neo_precision
  USE neo_input
  USE neo_units
  USE neo_control
  USE neo_work
  USE neo_exchange
!***********************************************************************
! Local definitions
!***********************************************************************
  IMPLICIT NONE

  INTEGER :: i,j,j_m,j_n
  INTEGER :: m,n,num_m,num_n,m_found,n_found
  INTEGER :: mm,nn
  INTEGER :: i_alloc
  INTEGER :: id1,id2,id3,id4,id5,id6,id7
  INTEGER :: extra_count
  LOGICAL :: extra_zero
  CHARACTER(5) :: dummy
  CHARACTER(45) :: cdum
  REAL(kind=dp) :: xra, xrm
  REAL(kind=dp) :: r_small, r_big
!***********************************************************************
! Open input-unit and read first quantities
!***********************************************************************
  OPEN(unit=r_u1,file=in_file,status='old',form='formatted')
!***********************************************************************
  IF (inp_swi .EQ. 1) THEN        !Princeton Boozer file        
     READ (r_u1,*) dummy
     READ (r_u1,*) m0b,n0b,ns,nfp,flux
     m_max = m0b+1
     n_max = 2*n0b+1
     mnmax = m_max*n_max
! **********************************************************************
! Allocate storage arrays
! **********************************************************************
     ALLOCATE(ixm(mnmax), ixn(mnmax), stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for integer arrays failed!'

     ALLOCATE(pixm(mnmax), pixn(mnmax), stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for integer arrays pointers failed!'

     ALLOCATE(i_m(m_max), i_n(n_max), stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for integer arrays failed!'

     ALLOCATE(es(ns), iota(ns), curr_pol(ns), curr_tor(ns),               & 
          pprime(ns), sqrtg00(ns), b00(ns), stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for real arrays failed!'

     ALLOCATE(rmnc(ns,mnmax), zmnc(ns,mnmax), lmnc(ns,mnmax),             &
          bmnc(ns,mnmax),                                                 &
          stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for fourier arrays (1) failed!'
!***********************************************************************
! Read input arrays
!***********************************************************************
     DO i =1, ns
        READ(r_u1,*) dummy
        READ(r_u1,*) es(i),iota(i),curr_pol(i),curr_tor(i),               &
             pprime(i),sqrtg00(i)
        READ(r_u1,*) dummy
        DO j=1,mnmax
           READ(r_u1,*) ixm(j),ixn(j),                                    &
                rmnc(i,j),zmnc(i,j),lmnc(i,j),                            &
                bmnc(i,j)
        END DO
     END DO

  ELSE IF(inp_swi .EQ. 2) THEN        !ORNL Boozer file
     READ(r_u1,'(a)') cdum
     READ(r_u1,'(a)') cdum
     READ(r_u1,'(a45,10i5)') cdum,ns,id1,id2,id3,id4,m0b,n0b,id5,nfp,mnmax
     READ(r_u1,'(a44,10i5)') cdum,id1,id2,id3,id4,id5,id6,id7
     DO i=1,4
        READ(r_u1,'(a)')cdum
     END DO
     ns = ns - 1
     m_max = m0b+1
     n_max = 2*n0b+1
!     mnmax = m_max*n_max - n0b
! **********************************************************************
! Allocate storage arrays
! **********************************************************************
     ALLOCATE(ixm(mnmax), ixn(mnmax), stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for integer arrays failed!'

     ALLOCATE(pixm(mnmax), pixn(mnmax), stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for integer arrays pointers failed!'

     ALLOCATE(i_m(m_max), i_n(n_max), stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for integer arrays failed!'

     ALLOCATE(es(ns), iota(ns), curr_pol(ns), curr_tor(ns),               & 
          pprime(ns), sqrtg00(ns), b00(ns), stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for real arrays failed!'

     ALLOCATE(rmnc(ns,mnmax), zmnc(ns,mnmax), lmnc(ns,mnmax),             &
          bmnc(ns,mnmax),                                                 &
          stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for fourier arrays (1) failed!'
!***********************************************************************
! Read input arrays
!***********************************************************************
     DO  i = 1,ns
        READ(r_u1,'(a)') cdum
        READ(r_u1,'(5e12.4)') es(i),iota(i),curr_pol(i),curr_tor(i),flux
        pprime(i) = 0.; sqrtg00(i) = 0.
        READ(r_u1,'(a)') cdum
        READ(r_u1,"(2i5,1p,4e16.8)") (ixm(j),ixn(j),rmnc(i,j),zmnc(i,j),  &
             lmnc(i,j),bmnc(i,j),j=1,mnmax)
        READ(r_u1,'(a)') cdum
     END DO

  ELSE IF(inp_swi .EQ. 3) THEN        !LHD Boozer file
     READ (r_u1,*) mnmax,ns,n0b,nfp,xra,xrm
     n_max = 2*n0b+1
     m_max = (mnmax-n0b-1)/n_max + 1
     m0b = m_max-1
     mnmax = mnmax + n0b ! negativ n for m=0
! **********************************************************************
! Allocate storage arrays
! **********************************************************************
     ALLOCATE(ixm(mnmax), ixn(mnmax), stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for integer arrays failed!'

     ALLOCATE(pixm(mnmax), pixn(mnmax), stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for integer arrays pointers failed!'

     ALLOCATE(i_m(m_max), i_n(n_max), stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for integer arrays failed!'

     ALLOCATE(es(ns), iota(ns), curr_pol(ns), curr_tor(ns),               & 
          pprime(ns), sqrtg00(ns), b00(ns), stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for real arrays failed!'

     ALLOCATE(rmnc(ns,mnmax), zmnc(ns,mnmax), lmnc(ns,mnmax),             &
          bmnc(ns,mnmax),                                                 &
          stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for fourier arrays (1) failed!'
!***********************************************************************
! Read input arrays
!***********************************************************************
     DO i =1, ns
        READ(r_u1,*) es(i),iota(i),curr_tor(i),curr_pol(i)
     END DO
     pprime  = 0.0_dp ! not given
     sqrtg00 = 0.0_dp
     flux = es(ns)
     DO i =1, ns      ! normalize
        es(i) = es(i) / flux
     END DO
     
     DO j = 1,n0b
        ixm(j) = 0
        ixn(j) = (-n0b + j - 1)*nfp
        bmnc(:,j) = 0.0_dp
        rmnc(:,j) = 0.0_dp
        zmnc(:,j) = 0.0_dp
        lmnc(:,j) = 0.0_dp
     END DO
     DO j = n0b+1,mnmax
        READ(r_u1,'(2i4)') mm,nn
        ixm(j) = mm
        ixn(j) = nn
        DO i=1,ns
           READ(r_u1,*) bmnc(i,j),rmnc(i,j),zmnc(i,j),lmnc(i,j)
        END DO
     END DO

  ELSEIF (inp_swi .EQ. 4) THEN        !W7-X File        
     READ (r_u1,*) dummy
     READ (r_u1,*) m0b,n0b,ns,nfp,flux,r_small,r_big
     m_max = m0b+1
     n_max = 2*n0b+1
     mnmax = m_max*n_max

!!$     PRINT *, 'm0b     ',m0b
!!$     PRINT *, 'n0b     ',n0b
!!$     PRINT *, 'ns      ',ns
!!$     PRINT *, 'nfp     ',nfp
!!$     PRINT *, 'flux    ',flux
!!$     PRINT *, 'r_small ',r_small
!!$     PRINT *, 'r_big   ',r_big
!!$     PRINT *, 'm_max   ',m_max
!!$     PRINT *, 'n_max   ',n_max
!!$     PRINT *, 'mnmax   ',mnmax
!!$     PAUSE

! **********************************************************************
! Allocate storage arrays
! **********************************************************************
     ALLOCATE(ixm(mnmax), ixn(mnmax), stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for integer arrays failed!'

     ALLOCATE(pixm(mnmax), pixn(mnmax), stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for integer arrays pointers failed!'

     ALLOCATE(i_m(m_max), i_n(n_max), stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for integer arrays failed!'

     ALLOCATE(es(ns), iota(ns), curr_pol(ns), curr_tor(ns),               & 
          pprime(ns), sqrtg00(ns), b00(ns), stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for real arrays failed!'

     ALLOCATE(rmnc(ns,mnmax), zmnc(ns,mnmax), lmnc(ns,mnmax),             &
          bmnc(ns,mnmax),                                                 &
          stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for fourier arrays (1) failed!'
!***********************************************************************
! Read input arrays
!***********************************************************************
     DO i =1, ns
        READ(r_u1,*) dummy
        READ(r_u1,*) es(i),iota(i),curr_pol(i),curr_tor(i),               &
             pprime(i),sqrtg00(i)
        READ(r_u1,*) dummy
        DO j=1,mnmax
           READ(r_u1,*) ixm(j),ixn(j),                                    &
                rmnc(i,j),zmnc(i,j),lmnc(i,j),                            &
                bmnc(i,j)
        END DO
     END DO

  ELSE IF(inp_swi .EQ. 5) THEN        !CHS Boozer file
     READ (r_u1,*) mnmax,ns,n0b,nfp
     n_max = 2*n0b+1
     m_max = (mnmax-n0b-1)/n_max + 1
     m0b = m_max-1
     mnmax = mnmax + n0b ! negativ n for m=0
! **********************************************************************
! Allocate storage arrays
! **********************************************************************
     ALLOCATE(ixm(mnmax), ixn(mnmax), stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for integer arrays failed!'

     ALLOCATE(pixm(mnmax), pixn(mnmax), stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for integer arrays pointers failed!'

     ALLOCATE(i_m(m_max), i_n(n_max), stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for integer arrays failed!'

     ALLOCATE(es(ns), iota(ns), curr_pol(ns), curr_tor(ns),               & 
          pprime(ns), sqrtg00(ns), b00(ns), stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for real arrays failed!'

     ALLOCATE(rmnc(ns,mnmax), zmnc(ns,mnmax), lmnc(ns,mnmax),             &
          bmnc(ns,mnmax),                                                 &
          stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for fourier arrays (1) failed!'
!***********************************************************************
! Read input arrays
!***********************************************************************
     DO i =1, ns
        READ(r_u1,*) es(i),iota(i),curr_tor(i),curr_pol(i)
     END DO
     pprime  = 0.0_dp ! not given
     sqrtg00 = 0.0_dp
     flux = es(ns)
     DO i =1, ns      ! normalize
        es(i) = es(i) / flux
     END DO
     
     DO j = 1,n0b
        ixm(j) = 0
        ixn(j) = (-n0b + j - 1)*nfp
        bmnc(:,j) = 0.0_dp
        rmnc(:,j) = 0.0_dp
        zmnc(:,j) = 0.0_dp
        lmnc(:,j) = 0.0_dp
     END DO
     DO j = n0b+1,mnmax
        READ(r_u1,'(2i4)') mm,nn
        ixm(j) = mm
        ixn(j) = nn
        DO i=1,ns
           READ(r_u1,*) bmnc(i,j),rmnc(i,j),zmnc(i,j),lmnc(i,j)
        END DO
     END DO

  ELSEIF (inp_swi .EQ. 6) THEN        ! NEW IPP, HSX        
     READ (r_u1,*) dummy
     READ (r_u1,*) dummy
     READ (r_u1,*) dummy
     READ (r_u1,*) dummy
     READ (r_u1,*) dummy
     READ (r_u1,*) m0b,n0b,ns,nfp,flux,r_small,r_big
     m_max = m0b+1
     n_max = 2*n0b+1
     mnmax = m_max*n_max
     ! m = 0 , n only >= 0
     ! mnmax = m0b * n_max + n0b + 1

! **********************************************************************
! Allocate storage arrays
! **********************************************************************
     ALLOCATE(ixm(mnmax), ixn(mnmax), stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for integer arrays failed!'

     ALLOCATE(pixm(mnmax), pixn(mnmax), stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for integer arrays pointers failed!'

     ALLOCATE(i_m(m_max), i_n(n_max), stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for integer arrays failed!'

     ALLOCATE(es(ns), iota(ns), curr_pol(ns), curr_tor(ns),               & 
          pprime(ns), sqrtg00(ns), b00(ns), stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for real arrays failed!'

     ALLOCATE(rmnc(ns,mnmax), zmnc(ns,mnmax), lmnc(ns,mnmax),             &
          bmnc(ns,mnmax),                                                 &
          stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for fourier arrays (1) failed!'
!***********************************************************************
! Read input arrays
!***********************************************************************
     DO i =1, ns
        READ(r_u1,*) dummy
        READ(r_u1,*) dummy
        READ(r_u1,*) es(i),iota(i),curr_pol(i),curr_tor(i),               &
             pprime(i),sqrtg00(i)
        READ(r_u1,*) dummy
        
        extra_zero = .FALSE.
        extra_count = 0
        DO j=1,mnmax
           IF (j .GT. 1 .AND. ixm(j-1) .EQ. 0 .AND. ixn(j-1) .EQ. 0) THEN
               extra_zero = .TRUE.
           ENDIF
           IF (extra_zero) THEN
               extra_count =  extra_count + 1
               IF (extra_count .EQ. n0b) extra_zero = .FALSE.
               ixm(j) = 0
               ixn(j) = -extra_count
               rmnc(i,j) = 0.0d0 
               zmnc(i,j) = 0.0d0 
               lmnc(i,j) = 0.0d0 
               bmnc(i,j) = 0.0d0 
           ELSE
               READ(r_u1,*) ixm(j),ixn(j),                                    &
                    rmnc(i,j),zmnc(i,j),lmnc(i,j),                            &
                    bmnc(i,j)
           ENDIF
        END DO
     END DO

  ELSEIF (inp_swi .EQ. 7) THEN        !QPS File        
     READ (r_u1,*) dummy
     READ (r_u1,*) m0b,n0b,ns,nfp,flux,r_small,r_big
     m_max = m0b+1
     n_max = 2*n0b+1
     mnmax = m_max*n_max

!!$     PRINT *, 'm0b     ',m0b
!!$     PRINT *, 'n0b     ',n0b
!!$     PRINT *, 'ns      ',ns
!!$     PRINT *, 'nfp     ',nfp
!!$     PRINT *, 'flux    ',flux
!!$     PRINT *, 'r_small ',r_small
!!$     PRINT *, 'r_big   ',r_big
!!$     PRINT *, 'm_max   ',m_max
!!$     PRINT *, 'n_max   ',n_max
!!$     PRINT *, 'mnmax   ',mnmax
!!$     PAUSE

! **********************************************************************
! Allocate storage arrays
! **********************************************************************
     ALLOCATE(ixm(mnmax), ixn(mnmax), stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for integer arrays failed!'

     ALLOCATE(pixm(mnmax), pixn(mnmax), stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for integer arrays pointers failed!'

     ALLOCATE(i_m(m_max), i_n(n_max), stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for integer arrays failed!'

     ALLOCATE(es(ns), iota(ns), curr_pol(ns), curr_tor(ns),               & 
          pprime(ns), sqrtg00(ns), b00(ns), stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for real arrays failed!'

     ALLOCATE(rmnc(ns,mnmax), zmnc(ns,mnmax), lmnc(ns,mnmax),             &
          bmnc(ns,mnmax),                                                 &
          stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for fourier arrays (1) failed!'
!***********************************************************************
! Read input arrays
!***********************************************************************
     DO i =1, ns
        READ(r_u1,*) dummy
        READ(r_u1,*) es(i),iota(i),curr_pol(i),curr_tor(i),               &
             pprime(i)
        sqrtg00(i) = 0.0d0
        READ(r_u1,*) dummy
        do j=1,n0b
           ixm(j) = 0
           ixn(j) = -nfp*(n0b-j+1)
           rmnc(i,j) = 0.0d0
           zmnc(i,j) = 0.0d0
           lmnc(i,j) = 0.0d0
           bmnc(i,j) = 0.0d0
        end do
        DO j=n0b+1,mnmax
           READ(r_u1,*) ixm(j),ixn(j),                                    &
                rmnc(i,j),zmnc(i,j),lmnc(i,j),                            &
                bmnc(i,j)
        END DO
     END DO
  ELSEIF (inp_swi .EQ. 8) THEN        ! NEW IPP TOKAMAK        
     READ (r_u1,*) dummy
     READ (r_u1,*) dummy
     READ (r_u1,*) dummy
     READ (r_u1,*) dummy
     READ (r_u1,*) dummy
     READ (r_u1,*) m0b,n0b,ns,nfp,flux,r_small,r_big
     m_max = m0b+1
     n_max = 2*n0b+1
     mnmax = m_max*n_max
     ! print *, 'm_max,n_max,mnmax: ',m_max,n_max,mnmax
     ! m = 0 , n only >= 0
     ! mnmax = m0b * n_max + n0b + 1

! **********************************************************************
! Allocate storage arrays
! **********************************************************************
     ALLOCATE(ixm(mnmax), ixn(mnmax), stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for integer arrays failed!'

     ALLOCATE(pixm(mnmax), pixn(mnmax), stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for integer arrays pointers failed!'

     ALLOCATE(i_m(m_max), i_n(n_max), stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for integer arrays failed!'

     ALLOCATE(es(ns), iota(ns), curr_pol(ns), curr_tor(ns),               & 
          pprime(ns), sqrtg00(ns), b00(ns), stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for real arrays failed!'

     ALLOCATE(rmnc(ns,mnmax), zmnc(ns,mnmax), lmnc(ns,mnmax),             &
          bmnc(ns,mnmax),                                                 &
          stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for fourier arrays (1) failed!'
!***********************************************************************
! Read input arrays
!***********************************************************************
     DO i =1, ns
        READ(r_u1,*) dummy
        READ(r_u1,*) dummy
        READ(r_u1,*) es(i),iota(i),curr_pol(i),curr_tor(i),               &
             pprime(i),sqrtg00(i)
        READ(r_u1,*) dummy
        
        DO j=1,mnmax
           !print *, 'j: ',j
           READ(r_u1,*) ixm(j),ixn(j),                                    &
                rmnc(i,j),zmnc(i,j),lmnc(i,j),                            &
                bmnc(i,j)
           !print *, 'ixm,ixn: ',ixm(j),ixn(j)
        END DO
     END DO

  ELSE
     WRITE (w_us,*) 'FATAL: There is yet no other input type defined'
     STOP
  END IF
!
! Filling of i_m and i_n 
! and pointers pixm from ixm to i_m, and pixn from ixn to i_n
  DO j = 1,mnmax
     m = ixm(j)
     n = ixn(j)
     IF (j .EQ. 1) THEN
        num_m = 1
        i_m(num_m) = m
        pixm(j) = num_m
        num_n = 1
        i_n(num_n) = n
        pixn(j) = num_n
     ELSE
        m_found = 0
        DO j_m = 1, num_m
           IF (m .EQ. i_m(j_m)) THEN
              pixm(j) = j_m
              m_found = 1
           END IF
        END DO
        IF (m_found .EQ. 0) THEN
           num_m = num_m + 1
           i_m(num_m) = m
           pixm(j) = num_m
        END IF
        n_found = 0
        DO j_n = 1, num_n
           IF (n .EQ. i_n(j_n)) THEN
              pixn(j) = j_n
              n_found = 1
           END IF
        END DO
        IF (n_found .EQ. 0) THEN
           num_n = num_n + 1
           i_n(num_n) = n
           pixn(j) = num_n
        END IF
     END IF
  END DO

!!$  PRINT *, 'i_m'
!!$  PRINT *, i_m
!!$  PRINT *, 'i_n'
!!$  PRINT *, i_n
!!$  PRINT *, 'ixm(1:mnmax:2*n0b+1)'
!!$  PRINT *, ixm(1:mnmax:2*n0b+1)
!!$  PRINT *, 'ixn(1:2*n0b+1)'
!!$  PRINT *, ixn(1:2*n0b+1)
!!$  PAUSE
!!$  PRINT *, 'pixm'
!!$  PRINT *, pixm
!!$  PRINT *, 'pixn'
!!$  PRINT *, pixn
!!$  PAUSE

  IF (lab_swi .EQ. 1) THEN              ! NCSX Boozer file
! 
! ATTENTION: Switch n TO -n
!            Toroidal mode numbers have to multiplied by number of field periods
!            Change iota to iota*nfp (PRINCETON)
!            it is named iota but actually it is iotabar
     ixn =  - ixn * nfp
     i_n =  - i_n * nfp
     max_n_mode = max_n_mode * nfp
     iota = iota*nfp
     psi_pr = ABS(flux) / twopi 
     curr_pol = curr_pol / twopi * nfp
     curr_tor = curr_tor / twopi
  ELSE IF  (lab_swi .EQ. 2) THEN         !ORNL Boozer file
     DO i=1,ns
        DO j=1,mnmax
           lmnc(i,j) = REAL(nfp)*lmnc(i,j)/TWOPI
        END DO
     END DO
     psi_pr = ABS(flux)
  ELSE IF  (lab_swi .EQ. 3) THEN         !LHD Boozer file
     psi_pr = ABS(flux)
     ixn =  - ixn ! change n to -n
     i_n =  - i_n
     ixm =  - ixm ! change m to -m (rhs coord. syst)
     i_m =  - i_m
  ELSE IF  (lab_swi .EQ. 4) THEN         ! W7X file
     curr_pol = - curr_pol * 2.d-7 * nfp   ! ? -
     curr_tor = curr_tor * 2.d-7 
     max_n_mode = max_n_mode * nfp
     ixn =  ixn * nfp
     i_n =  i_n * nfp
     ixm =  ixm 
     i_m =  i_m
     psi_pr = ABS(flux) / twopi
  ELSE IF  (lab_swi .EQ. 5) THEN         !CHS Boozer file
     psi_pr = ABS(flux)
     ixn =  - ixn ! change n to -n
     i_n =  - i_n
     ixm =  - ixm ! change m to -m
     i_m =  - i_m
     curr_tor = curr_tor
     curr_pol = curr_pol
     iota     = iota
  ELSE IF  (lab_swi .EQ. 6) THEN         ! NEW IPP, HSX
     curr_pol = - curr_pol * 2.d-7 * nfp   ! ? -
     curr_tor = curr_tor * 2.d-7 * nfp ! Henning  
     max_n_mode = max_n_mode * nfp
     ixn =  ixn * nfp
     i_n =  i_n * nfp
     ixm =  ixm 
     i_m =  i_m
     psi_pr = ABS(flux) / twopi
  ELSE IF  (lab_swi .EQ. 7) THEN         !QPS Boozer file
     DO i=1,ns
        DO j=1,mnmax
           lmnc(i,j) = REAL(nfp)*lmnc(i,j)/TWOPI
        END DO
     END DO
     psi_pr = ABS(flux) / twopi * nfp
  ELSE IF  (lab_swi .EQ. 8) THEN         ! NEW IPP TOKAMAKK
     curr_pol = - curr_pol * 2.d-7 * nfp   ! ? -
     curr_tor = curr_tor * 2.d-7 * nfp ! Henning  
     max_n_mode = max_n_mode * nfp
     ixn =  ixn * nfp
     i_n =  i_n * nfp
     ixm =  ixm 
     i_m =  i_m
     psi_pr = ABS(flux) / twopi
  ELSE
     WRITE(w_us,*) 'FATAL: There is yet no other Laboratory defined!'
     STOP
  END IF
!
! ATTENTION THIS IS JUST FOR TESTING
!
! For scaling of B change the following three with the same factor
! eps_eff should then stay unchanged if the reference for B and R is the same!
!
! bmnc = bmnc * 2.0_dp
! curr_pol = curr_pol * 2.0_dp
! curr_tor = curr_tor * 2.0_dp
!
! For scaling of R change the following four with the same factor
! eps_eff should then stay unchanged if the reference for B and R is the same!
!
! rmnc = rmnc * 2.0_dp
! zmnc = zmnc * 2.0_dp
! curr_pol = curr_pol * 2.0_dp
! curr_tor = curr_tor * 2.0_dp
! 
  CLOSE (unit=r_u1)
! **********************************************************************
! Write optional output for Plotting
! **********************************************************************
!!$  IF (write_output_files .NE. 0) THEN
!!$     IF (write_progress .NE. 0) WRITE (w_us,*) 'write dimension.dat'
!!$     OPEN(unit=w_u1,file='dimension.dat',status='unknown',form='formatted')
!!$     WRITE (w_u1,*) ns
!!$     WRITE (w_u1,*) mnmax
!!$     WRITE (w_u1,*) nfp
!!$     WRITE (w_u1,*) theta_n
!!$     WRITE (w_u1,*) phi_n
!!$     WRITE (w_u1,*) s_ind_in
!!$     CLOSE(unit=w_u1)
!!$  ENDIF
!!$  IF (write_output_files .NE. 0) THEN
!!$     IF (write_progress .NE. 0) WRITE (w_us,*) 'write es_arr.dat'
!!$     OPEN(unit=w_u1,file='es_arr.dat',status='unknown',form='formatted')
!!$     DO j=1,ns
!!$        WRITE(w_u1,format220) es(j),iota(j),curr_pol(j),curr_tor(j),       &
!!$             pprime(j),sqrtg00(j)
!!$     ENDDO
!!$     CLOSE(unit=w_u1)
!!$
!!$     IF (write_progress .NE. 0) WRITE (w_us,*) 'write mn_arr.dat'
!!$     OPEN(unit=w_u1,file='mn_arr.dat',status='unknown',form='formatted')
!!$     DO j = 1,mnmax
!!$        WRITE(w_u1,*)ixm(j),ixn(j)
!!$     ENDDO
!!$     CLOSE(unit=w_u1)
!!$

!!$  ! Write rmnc, zmnc, lmnc, bmnc
!!$     IF (write_progress .NE. 0) WRITE (w_us,*) 'write rmnc_arr.dat'
!!$     OPEN(unit=w_u1,file='rmnc_arr.dat',status='unknown',form='formatted')
!!$     WRITE(w_u1,*) ns
!!$     WRITE(w_u1,*) mnmax
!!$     DO i=1,ns
!!$        DO j=1,mnmax
!!$           WRITE(w_u1,*) rmnc(i,j)
!!$        END DO
!!$     END DO
!!$     CLOSE(unit=w_u1)
!!$
!!$     IF (write_progress .NE. 0) WRITE (w_us,*) 'write zmnc_arr.dat'
!!$     OPEN(unit=w_u1,file='zmnc_arr.dat',status='unknown',form='formatted')
!!$     WRITE(w_u1,*) ns
!!$     WRITE(w_u1,*) mnmax
!!$     DO i=1,ns
!!$        DO j=1,mnmax
!!$           WRITE(w_u1,*) zmnc(i,j)
!!$        END DO
!!$     END DO
!!$     CLOSE(unit=w_u1)
!!$
!!$     IF (write_progress .NE. 0) WRITE (w_us,*) 'write lmnc_arr.dat'
!!$     OPEN(unit=w_u1,file='lmnc_arr.dat',status='unknown',form='formatted')
!!$     WRITE(w_u1,*) ns
!!$     WRITE(w_u1,*) mnmax
!!$     DO i=1,ns
!!$        DO j=1,mnmax
!!$           WRITE(w_u1,*) lmnc(i,j)
!!$        END DO
!!$     END DO
!!$     CLOSE(unit=w_u1)
!!$
!!$     IF (write_progress .NE. 0) WRITE (w_us,*) 'write bmnc_arr.dat'
!!$     OPEN(unit=w_u1,file='bmnc_arr.dat',status='unknown',form='formatted')
!!$     WRITE(w_u1,*) ns
!!$     WRITE(w_u1,*) mnmax
!!$     DO i=1,ns
!!$        DO j=1,mnmax
!!$           WRITE(w_u1,*) bmnc(i,j)
!!$        END DO
!!$     END DO
!!$     CLOSE(unit=w_u1)
!!$
!!$     STOP

!!$  ENDIF
! **********************************************************************



  RETURN
END SUBROUTINE neo_read
! **********************************************************************

SUBROUTINE neo_prep_b00
  ! Preparation of b00 with splines
  !***********************************************************************
  ! Modules
  USE neo_precision
  USE neo_input
  USE neo_spline_b00
  USE inter_interfaces, ONLY: splinecof3, tf 
  INTEGER                             :: i, j
  INTEGER(I4B) :: sw1, sw2
  REAL(dp)     :: m0, c1, cn

  REAL(dp), DIMENSION(:), ALLOCATABLE :: lambda
  INTEGER,  DIMENSION(:), ALLOCATABLE :: index_i
  
  ! collect the b00 in a vector
  ALLOCATE ( lambda(ns) )
  ALLOCATE ( index_i(ns) )
  ALLOCATE ( a_b00(ns), b_b00(ns) )
  ALLOCATE ( c_b00(ns), d_b00(ns) )

  DO i = 1, ns
     DO j= 1, mnmax
        IF (ixm(j) .EQ. 0 .AND. ixn(j) .EQ. 0) THEN
           b00(i) = bmnc(i,j)
           EXIT
        END IF
     END DO
  END DO
  ! spline b00 along s
  ! boundary types (natural spline)
  sw1 = 2
  sw2 = 4
  ! input for test function for spline
  m0  = 0.0_dp
! boundary condition for spline
  c1 = 0.0_dp  
  cn = 0.0_dp
! we use no smoothing for spline
  lambda = 1.0D0
  index_i = (/ (i, i=1,ns) /) 
  CALL splinecof3(es, b00, c1, cn, lambda, index_i, sw1, sw2, &
       a_b00, b_b00, c_b00, d_b00, m0, tf)

  DEALLOCATE( lambda )
  DEALLOCATE( index_i )

END SUBROUTINE neo_prep_b00

SUBROUTINE neo_get_b00

  USE neo_precision
  USE neo_input
  USE neo_actual_fluxs
  USE neo_spline_b00
  USE inter_interfaces, ONLY: splint_horner3, tf, tfp, tfpp, tfppp 

  IMPLICIT NONE
  INTEGER(I4B)  :: swd = 1
  REAL(dp)      :: m0  = 0.0_dp
  REAL(dp)      :: ypp, yppp   ! dummies for derivatives

  CALL splint_horner3(es,a_b00,b_b00,c_b00,d_b00,swd,m0,     &
       s_es,tf,tfp,tfpp,tfppp,                               &
       s_b00,s_b00_s,ypp,yppp)

END SUBROUTINE neo_get_b00

SUBROUTINE neo_prep
! Preparation of Arrays
!***********************************************************************
! Modules
!***********************************************************************
  USE neo_input
  USE neo_work
  USE neo_exchange
  USE neo_parameters
  USE neo_control
  USE neo_units
  USE neo_spline
  USE neo_actual_spectra
!***********************************************************************
! Local definitions
!***********************************************************************
  IMPLICIT NONE

  INTEGER :: i_alloc
  INTEGER :: ip, it
! INTEGER :: k, j
  INTEGER :: im, in
  INTEGER :: m, n
! **********************************************************************
! Allocate Storage Arrays
! **********************************************************************
! ALLOCATE(cosval(phi_n,theta_n,mnmax),                                &
!          sinval(phi_n,theta_n,mnmax),                                &
!          stat = i_alloc)
  ALLOCATE(cosmth(theta_n,m_max),                                      &
           sinmth(theta_n,m_max),                                      &
           cosnph(phi_n,  n_max),                                      &
           sinnph(phi_n,  n_max),                                      &
           stat = i_alloc)
  IF(i_alloc /= 0) STOP 'Allocation for cos/sin-arrays failed!'
  ALLOCATE(theta_arr(theta_n),                                         &
           phi_arr(phi_n),                                             &
           stat = i_alloc)
  IF(i_alloc /= 0) STOP 'Allocation for theta/phi-arrays failed!'
! **********************************************************************
! Allocation for arrays for output quantities 
! **********************************************************************
  ALLOCATE(b(theta_n,phi_n),stat = i_alloc)
  IF(i_alloc /= 0) STOP 'Allocation for b-array failed!'
  ALLOCATE(sqrg11(theta_n,phi_n),stat = i_alloc)
  IF(i_alloc /= 0) STOP 'Allocation for sqrg11-array failed!'
  ALLOCATE(kg(theta_n,phi_n),stat = i_alloc)
  IF(i_alloc /= 0) STOP 'Allocation for kg-array failed!'
  ALLOCATE(pard(theta_n,phi_n),stat = i_alloc)
  IF(i_alloc /= 0) STOP 'Allocation for pard-array failed!'
  IF (calc_cur .EQ. 1 .OR. calc_fourier .EQ. 0) THEN
     ALLOCATE(bqtphi(theta_n,phi_n),stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for bqtphi-array failed!'
  END IF
  ALLOCATE(r_nabpsi(theta_n,phi_n),stat = i_alloc)
  IF(i_alloc /= 0) STOP 'Allocation for r_nabpsi-array failed!'
! **********************************************************************
! Allocation for spline arrays
! **********************************************************************
  ALLOCATE(b_spl(4,4,theta_n,phi_n),                                   &
           k_spl(4,4,theta_n,phi_n),                                   &
           g_spl(4,4,theta_n,phi_n),                                   &
           p_spl(4,4,theta_n,phi_n),                                   &
           r_spl(4,4,theta_n,phi_n),                                   &
           stat = i_alloc)
  IF(i_alloc /= 0) STOP 'Allocation for spline-arrays failed!'
  IF (calc_cur .EQ. 1) THEN
     ALLOCATE(q_spl(4,4,theta_n,phi_n),                                 &
              stat = i_alloc)
     IF(i_alloc /= 0) STOP 'Allocation for q_spl--array failed!'
  END IF
! **********************************************************************
! Allocation of spectra on flux surfaces
! **********************************************************************
  ALLOCATE(s_rmnc(mnmax),                                              &
           s_zmnc(mnmax),                                              &
           s_lmnc(mnmax),                                              &
           s_bmnc(mnmax),                                              &
          stat = i_alloc)
  IF(i_alloc /= 0) STOP 'Allocation for spectra on flux surfaces!'
! **********************************************************************
! Some initial work
! **********************************************************************
  theta_start = 0.0
  theta_end   = twopi
  theta_int   = (theta_end-theta_start)/(theta_n-1)
  phi_start   = 0.0
  phi_end     = twopi / nfp
  phi_int     = (phi_end-phi_start)/(phi_n-1)
! **********************************************************************
! Preparation of arrays
! **********************************************************************
  DO it=1,theta_n
    theta_arr(it) = theta_start + theta_int*(it-1)
  ENDDO

  DO ip=1,phi_n
    phi_arr(ip) = phi_start + phi_int*(ip-1)
  ENDDO

! DO imn=1,mnmax
!   ixm_i = ixm(imn)
!   ixn_i = ixn(imn)
!   DO it=1,theta_n
!     DO ip=1,phi_n
!       sinval(ip,it,imn) = SIN(ixm_i*theta_arr(it) - ixn_i*phi_arr(ip))
!       cosval(ip,it,imn) = COS(ixm_i*theta_arr(it) - ixn_i*phi_arr(ip))
!     ENDDO
!   ENDDO
! ENDDO
  DO im = 1,m_max
     m = i_m(im)
     IF (ABS(m) .LE. max_m_mode) THEN
        DO it=1,theta_n
           sinmth(it,im) = SIN( m * theta_arr(it) )
           cosmth(it,im) = COS( m * theta_arr(it) )
        ENDDO
     END IF
  ENDDO
  DO in = 1,n_max
     n = i_n(in)
     IF (ABS(n) .LE. max_n_mode) THEN
        DO ip=1,phi_n
           sinnph(ip,in) = SIN( n * phi_arr(ip) )
           cosnph(ip,in) = COS( n * phi_arr(ip) )
        ENDDO
     END IF
  ENDDO
! **********************************************************************
! Write optional output
! **********************************************************************
!!$  IF (write_output_files .NE. 0) THEN
!!$     IF (write_progress .NE. 0) WRITE (w_us,*) 'write theta_arr.dat'
!!$     OPEN(unit=w_u1,file='theta_arr.dat',status='unknown',form='formatted')
!!$     DO j=1,theta_n
!!$        WRITE(w_u1,*) theta_arr(j)
!!$     ENDDO
!!$     CLOSE(unit=w_u1)
!!$
!!$     IF (write_progress .NE. 0) WRITE (w_us,*) 'write phi_arr.dat'
!!$     OPEN(unit=w_u1,file='phi_arr.dat',status='unknown',form='formatted')
!!$     DO k=1,phi_n
!!$        WRITE(w_u1,*) phi_arr(k)
!!$     ENDDO
!!$     CLOSE(unit=w_u1)
!!$  ENDIF
! **********************************************************************
  RETURN
END SUBROUTINE neo_prep
! **********************************************************************

SUBROUTINE neo_fourier
! Summation of Fourier Sums and Computation of Derived Quantities
!***********************************************************************
! Modules
!***********************************************************************
  USE neo_precision
  USE neo_input
  USE neo_work
  USE neo_exchange
  USE neo_parameters
  USE neo_control
  USE neo_units
  USE neo_support
  USE neo_actual_fluxs
  USE neo_actual_spectra
!***********************************************************************
! Local definitions
!***********************************************************************
  IMPLICIT NONE

  INTEGER       :: i_alloc
  INTEGER       :: im, in, m, n
! INTEGER       :: i, j
  INTEGER       :: it, ip, imn
  INTEGER       :: uw = 100
  REAL(kind=dp) :: ri, zi, li, bi
  REAL(kind=dp) :: cosv, sinv
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: hh1, hh2
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: sqrga,sqrgb
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: currpol_b,currtor_b
  REAL(kind=dp) :: shh1, shh2, mhh, mhh2, shh
  REAL(kind=dp) :: scurrpol_b, scurrpol_b2, scpol
  REAL(kind=dp) :: scurrtor_b, scurrtor_b2, sctor
  REAL(kind=dp) :: mcurrpol_b, mcurrpol_b2
  REAL(kind=dp) :: mcurrtor_b, mcurrtor_b2
!***********************************************************************
! Allocation of arrays
!***********************************************************************
  ALLOCATE(r(theta_n,phi_n),                                           &
           z(theta_n,phi_n),                                           &
           l(theta_n,phi_n),                                           &
           stat = i_alloc)
  IF(i_alloc /= 0) STOP 'Allocation for sum-arrays failed!'
  ALLOCATE(r_tb(theta_n,phi_n),                                        &
           z_tb(theta_n,phi_n),                                        &
           p_tb(theta_n,phi_n),                                        &
           b_tb(theta_n,phi_n),                                        &
           stat = i_alloc)
  IF(i_alloc /= 0) STOP 'Allocation for sum-arrays failed!'
  ALLOCATE(r_pb(theta_n,phi_n),                                        &
           z_pb(theta_n,phi_n),                                        &
           p_pb(theta_n,phi_n),                                        &
           b_pb(theta_n,phi_n),                                        &
           stat = i_alloc)
  IF(i_alloc /= 0) STOP 'Allocation for sum-arrays failed!'
  ALLOCATE(gtbtb(theta_n,phi_n),                                       &
           gpbpb(theta_n,phi_n),                                       &
           gtbpb(theta_n,phi_n),                                       &
           isqrg(theta_n,phi_n),                                       &
           stat = i_alloc)
  IF(i_alloc /= 0) STOP 'Allocation for sum-arrays failed!'
  ALLOCATE(psi_r(theta_n,phi_n),                                       &
           psi_z(theta_n,phi_n),                                       &
           stat = i_alloc)
  IF(i_alloc /= 0) STOP 'Allocation for sum-arrays failed!'



!***********************************************************************
! Summation of Fourier components
!***********************************************************************
  r = 0.0d0
  z = 0.0d0
  l = 0.0d0
  b = 0.0d0

  r_tb = 0.0d0
  z_tb = 0.0d0
  p_tb = 0.0d0
  b_tb = 0.0d0
  r_pb = 0.0d0
  z_pb = 0.0d0
  p_pb = 0.0d0
  b_pb = 0.0d0
  
  DO imn=1,mnmax
     ri = s_rmnc(imn)
     zi = s_zmnc(imn)
     li = s_lmnc(imn)
     bi = s_bmnc(imn)
     m = ixm(imn)
     n = ixn(imn)
     im = pixm(imn)
     in = pixn(imn)
     ! PRINT *, m,n,im,in
     IF (ABS(m) .LE. max_m_mode .AND. ABS(n) .LE. max_n_mode) THEN
        DO ip=1,phi_n
           DO it=1,theta_n
!             cosv = cosval(ip,it,imn)
!             sinv = sinval(ip,it,imn)
              cosv = cosmth(it,im) * cosnph(ip,in) + sinmth(it,im) * sinnph(ip,in)
              sinv = sinmth(it,im) * cosnph(ip,in) - cosmth(it,im) * sinnph(ip,in)

              r(it,ip) = r(it,ip) + ri*cosv
              z(it,ip) = z(it,ip) + zi*sinv
              l(it,ip) = l(it,ip) + li*sinv
              b(it,ip) = b(it,ip) + bi*cosv

              r_tb(it,ip) = r_tb(it,ip) - m*ri*sinv
              r_pb(it,ip) = r_pb(it,ip) + n*ri*sinv
              z_tb(it,ip) = z_tb(it,ip) + m*zi*cosv
              z_pb(it,ip) = z_pb(it,ip) - n*zi*cosv
              p_tb(it,ip) = p_tb(it,ip) - m*li*cosv ! -l_tb
              p_pb(it,ip) = p_pb(it,ip) + n*li*cosv ! -l_pb
              b_tb(it,ip) = b_tb(it,ip) - m*bi*sinv
              b_pb(it,ip) = b_pb(it,ip) + n*bi*sinv
           END DO
        END DO
     END IF
  END DO

  IF (lab_swi .EQ. 5 .OR. lab_swi .EQ. 3) THEN ! CHS, LHD
      p_tb = - p_tb
      p_pb = 1 - p_pb
  ELSE
     p_tb = p_tb * twopi / nfp
     p_pb = 1.0_dp + p_pb * twopi / nfp
  END IF

! **********************************************************************
! Ensure periodicity boundaries to be the same
! **********************************************************************
  r(theta_n,:) = r(1,:)
  r(:,phi_n)   = r(:,1)
  z(theta_n,:) = z(1,:)
  z(:,phi_n)   = z(:,1)
  l(theta_n,:) = l(1,:)
  l(:,phi_n)   = l(:,1)
  b(theta_n,:) = b(1,:)
  b(:,phi_n)   = b(:,1)
  r_tb(theta_n,:) = r_tb(1,:)
  r_tb(:,phi_n)   = r_tb(:,1)
  r_pb(theta_n,:) = r_pb(1,:)
  r_pb(:,phi_n)   = r_pb(:,1)
  z_tb(theta_n,:) = z_tb(1,:)
  z_tb(:,phi_n)   = z_tb(:,1)
  z_pb(theta_n,:) = z_pb(1,:)
  z_pb(:,phi_n)   = z_pb(:,1)
  p_tb(theta_n,:) = p_tb(1,:)
  p_tb(:,phi_n)   = p_tb(:,1)
  p_pb(theta_n,:) = p_pb(1,:)
  p_pb(:,phi_n)   = p_pb(:,1)
  b_tb(theta_n,:) = b_tb(1,:)
  b_tb(:,phi_n)   = b_tb(:,1)
  b_pb(theta_n,:) = b_pb(1,:)
  b_pb(:,phi_n)   = b_pb(:,1)

  !PRINT *, theta_n
  !DO it = 1,theta_n
  !   PRINT *, it,r(it,10),z(it,10),l(it,10)
  !END DO
! **********************************************************************
! Derived quantities
! **********************************************************************
! metric tensor
  gtbtb = r_tb*r_tb + z_tb*z_tb + r*r*p_tb*p_tb  
  gpbpb = r_pb*r_pb + z_pb*z_pb + r*r*p_pb*p_pb  
  gtbpb = r_tb*r_pb + z_tb*z_pb + r*r*p_tb*p_pb  
! $1/sqrt(g)$
  fac = s_curr_pol + s_iota * s_curr_tor  ! (J + iota I)
  isqrg  = b*b / fac 
! $sqrt(g^{11})$
  IF (g11_swi .EQ. 0) THEN
     sqrg11 = 1.0_dp
  ELSE
     sqrg11 = SQRT( (gtbtb*gpbpb - gtbpb*gtbpb ) * isqrg**2)
  END IF
!  PRINT *, 'fac: ',fac,s_curr_pol,s_iota,s_curr_tor
! geodesic curvature term $k_G |\nabla \psi|$ 
  kg = (s_curr_tor * b_pb - s_curr_pol * b_tb) / fac
! parallel derivative of mod-B
  pard = b_pb + s_iota * b_tb
! quasi-toroidal phi component of b (only for parallel current)
  IF (calc_cur .EQ. 1 .OR. calc_fourier .EQ. 0) THEN
     bqtphi = isqrg * (p_pb + s_iota * p_tb)
  END IF

! Effective Radius Calculations
  psi_r = p_tb * z_pb - z_tb * p_pb
  psi_r = psi_r * r * isqrg
  psi_z = r_tb * p_pb - p_tb * r_pb
  psi_z = psi_z * r * isqrg
  r_nabpsi = r * psi_r + z * psi_z

! Testing of consistency
  IF (chk_swi .EQ. 1) THEN
     ALLOCATE (hh1(theta_n,phi_n),hh2(theta_n,phi_n))
     ALLOCATE (sqrga(theta_n,phi_n))
     ALLOCATE (sqrgb(theta_n,phi_n))
     ALLOCATE (currpol_b(theta_n,phi_n))
     ALLOCATE (currtor_b(theta_n,phi_n))

     hh1 = (fac/b)**2
     hh2 = s_iota**2 * gtbtb + 2.0_dp * s_iota * gtbpb + gpbpb
  
     sqrga = 1 / isqrg
     sqrgb = SQRT(hh2) / b
     currpol_b = (gpbpb + s_iota * gtbpb) / sqrgb
     currtor_b = (gtbpb + s_iota * gtbtb) / sqrgb

     shh1 = 0.0_dp
     shh2 = 0.0_dp
  
     scurrpol_b  = 0.0_dp
     scurrpol_b2 = 0.0_dp
     scurrtor_b  = 0.0_dp
     scurrtor_b2 = 0.0_dp

     DO ip=1,phi_n
        DO it=1,theta_n
           shh1  = shh1 + hh1(it,ip)/hh2(it,ip)
           shh2  = shh2 + (hh1(it,ip)/hh2(it,ip))**2
           scurrpol_b  = scurrpol_b  + currpol_b(it,ip)
           scurrpol_b2 = scurrpol_b2 + currpol_b(it,ip)**2
           scurrtor_b  = scurrtor_b  + currtor_b(it,ip)
           scurrtor_b2 = scurrtor_b2 + currtor_b(it,ip)**2
        END DO
     END DO
     mhh   = shh1 / phi_n / theta_n
     mhh2  = shh2 / phi_n / theta_n
     shh   = SQRT(ABS(mhh2 - mhh**2))
     
     mcurrpol_b  = scurrpol_b  / phi_n / theta_n
     mcurrpol_b2 = scurrpol_b2 / phi_n / theta_n
     mcurrtor_b  = scurrtor_b  / phi_n / theta_n
     mcurrtor_b2 = scurrtor_b2 / phi_n / theta_n
     
     scpol = SQRT(ABS(mcurrpol_b2 - mcurrpol_b**2))
     sctor = SQRT(ABS(mcurrtor_b2 - mcurrtor_b**2))

     CALL unit_check(uw)
     OPEN(unit=uw, file=chk_file, status='old',                     &
          position='append', action='write')
     WRITE(uw,'(18(1x,e17.10))')                                    &
          s_es, SQRT(mhh), shh,                                     &
          s_curr_pol, mcurrpol_b, scpol,                            &
          s_curr_tor, mcurrtor_b, sctor
     CLOSE(unit=uw)

     DEALLOCATE (hh1)
     DEALLOCATE (hh2)
     DEALLOCATE (sqrga)
     DEALLOCATE (sqrgb)
     DEALLOCATE (currpol_b)
     DEALLOCATE (currtor_b)
  END IF

! **********************************************************************
! Optional Output
! **********************************************************************
!!$  IF (write_output_files .NE. 0) THEN
!!$     IF (write_progress .NE. 0) WRITE (w_us,*) 'write b_s_arr.dat'
!!$     OPEN(unit=w_u1,file='b_s_arr.dat',status='unknown',form='formatted')
!!$     DO i=1,theta_n
!!$        DO j=1,phi_n
!!$           WRITE(w_u1,*)  b(i,j)
!!$        END DO
!!$     ENDDO
!!$     CLOSE(unit=w_u1)
!!$
!!$     IF (write_progress .NE. 0) WRITE (w_us,*) 'write r_s_arr.dat'
!!$     OPEN(unit=w_u1,file='r_s_arr.dat',status='unknown',form='formatted')
!!$     DO i=1,theta_n
!!$        DO j=1,phi_n
!!$           WRITE(w_u1,*)  r(i,j)
!!$        END DO
!!$     ENDDO
!!$     CLOSE(unit=w_u1)
!!$
!!$     IF (write_progress .NE. 0) WRITE (w_us,*) 'write z_s_arr.dat'
!!$     OPEN(unit=w_u1,file='z_s_arr.dat',status='unknown',form='formatted')
!!$     DO i=1,theta_n
!!$        DO j=1,phi_n
!!$           WRITE(w_u1,*)  z(i,j)
!!$        END DO
!!$     ENDDO
!!$     CLOSE(unit=w_u1)
!!$
!!$     IF (write_progress .NE. 0) WRITE (w_us,*) 'write l_s_arr.dat'
!!$     OPEN(unit=w_u1,file='l_s_arr.dat',status='unknown',form='formatted')
!!$     DO i=1,theta_n
!!$        DO j=1,phi_n
!!$           WRITE(w_u1,*)  l(i,j)
!!$        END DO
!!$     ENDDO
!!$     CLOSE(unit=w_u1)
!!$
!!$     IF (write_progress .NE. 0) WRITE (w_us,*) 'write isqrg_arr.dat'
!!$     OPEN(unit=w_u1,file='isqrg_arr.dat',status='unknown',form='formatted')
!!$     DO i=1,theta_n
!!$        DO j=1,phi_n
!!$           WRITE(w_u1,*)  isqrg(i,j)
!!$        END DO
!!$     ENDDO
!!$     CLOSE(unit=w_u1)
!!$
!!$     IF (write_progress .NE. 0) WRITE (w_us,*) 'write sqrg11_arr.dat'
!!$     OPEN(unit=w_u1,file='sqrg11_arr.dat',status='unknown',form='formatted')
!!$     DO i=1,theta_n
!!$        DO j=1,phi_n
!!$           WRITE(w_u1,*)  sqrg11(i,j)
!!$        END DO
!!$     ENDDO
!!$     CLOSE(unit=w_u1)
!!$
!!$     IF (write_progress .NE. 0) WRITE (w_us,*) 'write kg_arr.dat'
!!$     OPEN(unit=w_u1,file='kg_arr.dat',status='unknown',form='formatted')
!!$     DO i=1,theta_n
!!$        DO j=1,phi_n
!!$           WRITE(w_u1,*)  kg(i,j)
!!$        END DO
!!$     ENDDO
!!$     CLOSE(unit=w_u1)
!!$
!!$     IF (write_progress .NE. 0) WRITE (w_us,*) 'write pard_arr.dat'
!!$     OPEN(unit=w_u1,file='pard_arr.dat',status='unknown',form='formatted')
!!$     DO i=1,theta_n
!!$        DO j=1,phi_n
!!$           WRITE(w_u1,*)  pard(i,j)
!!$        END DO
!!$     ENDDO
!!$     CLOSE(unit=w_u1)
!!$  ENDIF
! **********************************************************************
  RETURN
END SUBROUTINE neo_fourier
! **********************************************************************

SUBROUTINE neo_spline2d
! Creation of Spline Arrays
! **********************************************************************
! Modules
! **********************************************************************
  USE neo_precision
  USE neo_input
  USE neo_work
  USE neo_exchange
  USE neo_units
  USE neo_parameters
  USE neo_control
  USE neo_spline
  USE spline_mod, ONLY: spl2d
! **********************************************************************
! Local Definitions
! **********************************************************************
  IMPLICIT NONE
! **********************************************************************
! Allocation of spline arrays is done in neo_prep
! **********************************************************************
! **********************************************************************
! Double periodic splines (parameter mt=1 and mp=1) 
! **********************************************************************
! Spline for mod b
! IF (write_progress .NE. 0) WRITE (w_us,*) 'before spl2d'
  CALL spl2d(theta_n,phi_n,theta_int,phi_int,mt,mp,b,b_spl)
! IF (write_progress .NE. 0) WRITE (w_us,*) 'after spl2d'
! Spline for sqrg11
! IF (write_progress .NE. 0) WRITE (w_us,*) 'before spl2d'
  CALL spl2d(theta_n,phi_n,theta_int,phi_int,mt,mp,sqrg11,g_spl)
! IF (write_progress .NE. 0) WRITE (w_us,*) 'after spl2d'
! Spline for geodesic curviture
! IF (write_progress .NE. 0) WRITE (w_us,*) 'before spl2d'
  CALL spl2d(theta_n,phi_n,theta_int,phi_int,mt,mp,kg,k_spl)
! IF (write_progress .NE. 0) WRITE (w_us,*) 'after spl2d'
! Spline for parallel derivative
! IF (write_progress .NE. 0) WRITE (w_us,*) 'before spl2d'
  CALL spl2d(theta_n,phi_n,theta_int,phi_int,mt,mp,pard,p_spl)
! IF (write_progress .NE. 0) WRITE (w_us,*) 'after spl2d'
! Spline for quasi-toroidal phi component of b
  IF (calc_cur .EQ. 1) THEN
!    IF (write_progress .NE. 0) WRITE (w_us,*) 'before spl2d'
     CALL spl2d(theta_n,phi_n,theta_int,phi_int,mt,mp,bqtphi,q_spl)
!    IF (write_progress .NE. 0) WRITE (w_us,*) 'after spl2d'
  END IF
! Spline for r_nabpsi
! IF (write_progress .NE. 0) WRITE (w_us,*) 'before spl2d'
  CALL spl2d(theta_n,phi_n,theta_int,phi_int,mt,mp,r_nabpsi,r_spl)
! IF (write_progress .NE. 0) WRITE (w_us,*) 'after spl2d'
! **********************************************************************
! Spline test
! **********************************************************************
  IF (spline_test .GT. 0) THEN
     IF (write_progress .NE. 0) WRITE (w_us,*) 'before neo_spline_test'
     CALL neo_spline_test
     IF (write_progress .NE. 0) WRITE (w_us,*) 'after neo_spline_test'
  ENDIF
! **********************************************************************
  RETURN
END SUBROUTINE neo_spline2d
! **********************************************************************

SUBROUTINE neo_eval(theta,phi,bval,gval,kval,pval,qval,rval)
! Evaluation of 2-D Splines
! **********************************************************************
! Modules
! **********************************************************************
  USE neo_precision
  USE neo_input
  USE neo_work
  USE neo_exchange
  USE neo_units
  USE neo_parameters
  USE neo_control
  USE neo_spline
  USE spline_mod, ONLY: eva2d, poi2d
  USE neo_eval_switch
  USE neo_actual_fluxs
  USE neo_actual_spectra
! **********************************************************************
! Local Definitions
! **********************************************************************
  IMPLICIT NONE

  REAL(kind=dp), INTENT(in)  ::   theta, phi 
  REAL(kind=dp), INTENT(out) ::   bval, gval, kval, pval, qval, rval
!
! New definitions for direct evaluation
!
  INTEGER       :: i_alloc
  INTEGER       :: imn, m, n, im, in
  REAL(kind=dp) :: ri, zi, li, bi
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: cosval,sinval
  REAL(kind=dp) :: cosv, sinv
  REAL(kind=dp) :: rr, zz, ll, bb
  REAL(kind=dp) :: rr_tb, zz_tb, bb_tb, pp_tb
  REAL(kind=dp) :: rr_pb, zz_pb, bb_pb, pp_pb
  REAL(kind=dp) :: ggtbtb, ggpbpb, ggtbpb, ffac, iisqrg
  REAL(kind=dp) :: ppsi_r, ppsi_z
!
  ALLOCATE(cosval(mnmax), sinval(mnmax), stat = i_alloc)
  IF(i_alloc /= 0) STOP 'Allocation for integer arrays failed!'
!
! **********************************************************************
! Using splines for evaluation
! **********************************************************************
  IF (eval_mode .EQ. 0) THEN
! **********************************************************************
! Evaluation of pointer
! **********************************************************************
     CALL poi2d(theta_int,phi_int,mt,mp,                               &
          theta_start,theta_end,phi_start,phi_end,                     &
          theta,phi,theta_ind,phi_ind,theta_d,phi_d,ierr)
! **********************************************************************
! Evaluation of 2d-splines
! **********************************************************************
     IF (eval_switch(1) .EQ. 1) THEN
        CALL eva2d(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,         &
             b_spl,bval)
     END IF
     IF (eval_switch(2) .EQ. 1) THEN
     CALL eva2d(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,         &
          g_spl,gval)
     END IF
     IF (eval_switch(3) .EQ. 1) THEN
     CALL eva2d(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,         &
          k_spl,kval)
     END IF
     IF (eval_switch(4) .EQ. 1) THEN
     CALL eva2d(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,         &
          p_spl,pval)
     END IF
     IF (eval_switch(5) .EQ. 1) THEN
        CALL eva2d(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,      &
             q_spl,qval)
     END IF
     IF (eval_switch(6) .EQ. 1) THEN
     CALL eva2d(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,         &
          r_spl,rval)
     END IF
! **********************************************************************
! Using direct evaluation
! **********************************************************************
  ELSEIF (eval_mode .EQ. 1) THEN

     sinval = SIN(ixm*theta - ixn*phi)
     cosval = COS(ixm*theta - ixn*phi)

!***********************************************************************
! Summation of Fourier components
!***********************************************************************
     rr = 0.0d0
     zz = 0.0d0
     ll = 0.0d0
     bb = 0.0d0
     
     rr_tb = 0.0d0
     zz_tb = 0.0d0
     pp_tb = 0.0d0
     bb_tb = 0.0d0
     rr_pb = 0.0d0
     zz_pb = 0.0d0
     pp_pb = 0.0d0
     bb_pb = 0.0d0
  
     DO imn = 1,mnmax
        ri = s_rmnc(imn)
        zi = s_zmnc(imn)
        li = s_lmnc(imn)
        bi = s_bmnc(imn)
        m = ixm(imn)
        n = ixn(imn)
        im = pixm(imn)
        in = pixn(imn)

        IF (ABS(m) .LE. max_m_mode .AND. ABS(n) .LE. max_n_mode) THEN
           cosv = cosval(imn)
           sinv = sinval(imn)

           rr = rr + ri*cosv
           zz = zz + zi*sinv
           ll = ll + li*sinv
           bb = bb + bi*cosv

           rr_tb = rr_tb - m*ri*sinv
           rr_pb = rr_pb + n*ri*sinv
           zz_tb = zz_tb + m*zi*cosv
           zz_pb = zz_pb - n*zi*cosv
           pp_tb = pp_tb - m*li*cosv
           pp_pb = pp_pb + n*li*cosv
           bb_tb = bb_tb - m*bi*sinv
           bb_pb = bb_pb + n*bi*sinv
        END IF
     END DO
     pp_tb = pp_tb * twopi / nfp
     pp_pb = 1.0_dp + pp_pb * twopi / nfp
! **********************************************************************
! Derived quantities
! **********************************************************************
     bval = bb
! metric tensor
     ggtbtb = rr_tb*rr_tb + zz_tb*zz_tb + rr*rr*pp_tb*pp_tb  
     ggpbpb = rr_pb*rr_pb + zz_pb*zz_pb + rr*rr*pp_pb*pp_pb  
     ggtbpb = rr_tb*rr_pb + zz_tb*zz_pb + rr*rr*pp_tb*pp_pb  
! $1/sqrt(g)$
     ffac = s_curr_pol + s_iota * s_curr_tor
     iisqrg  = bb*bb / ffac 
! $sqrt(g^{11})$
     IF (g11_swi .EQ. 0) THEN
        gval = 1.0_dp          ! sqrg11
     ELSE
        gval = SQRT( (ggtbtb*ggpbpb - ggtbpb*ggtbpb ) * iisqrg**2)
     END IF
! geodesic curvature term $k_G |\nabla \psi|$ 
     kval = (s_curr_tor * bb_pb - s_curr_pol*bb_tb) / ffac ! kg
! parallel derivative of mod-B
     pval = bb_pb + s_iota * bb_tb  ! pard
! quasi-toroidal phi component of b (only for parallel current)
     IF (calc_cur .EQ. 1) THEN
        qval = iisqrg * (pp_pb + s_iota * pp_tb) ! bqtphi
     END IF
! Effective Radius Calculations
     ppsi_r = (pp_tb * zz_pb - zz_tb * pp_pb) * rr * iisqrg
     ppsi_z = (rr_tb * pp_pb - pp_tb * rr_pb) * rr * iisqrg
     rval   = rr * ppsi_r + zz * ppsi_z    
!
  END IF
!
  RETURN
END SUBROUTINE neo_eval

! **********************************************************************
SUBROUTINE neo_bderiv(theta, phi, f, g, dfdx, dfdy, dgdx, dgdy)
!
! Calculates first and second derivatives of b using the 2d-splines
!
! Input:  theta, phi
! Output: f      db/dt              (t = theta)
!         g      db/dp              (p = phi)
!         dfdx   d^2b/dt^2
!         dfdy   d^2b/(dt dp)
!         dgdx   d^2b/(dt dp)
!         dgdy   d^2b/dp^2
! 
! Input/output consistent for neo_zeros2d
! **********************************************************************
! Modules
! **********************************************************************
  USE neo_precision
  USE neo_input
  USE neo_work
  USE neo_units
  USE neo_parameters
  USE neo_control
  USE neo_spline
  USE spline_mod, ONLY: poi2d, eva2d_sd, eva2d_fd
! **********************************************************************
! Local Definitions
! **********************************************************************
  IMPLICIT NONE

  REAL(kind=dp), INTENT(in)    ::   theta, phi 
  REAL(kind=dp), INTENT(out)   ::   f, g, dfdx, dfdy, dgdx, dgdy

  REAL(kind=dp), DIMENSION(2)  ::   fderiv 
  REAL(kind=dp), DIMENSION(3)  ::   sderiv 
! **********************************************************************
! Evaluation of pointer
! **********************************************************************
  CALL poi2d(theta_int,phi_int,mt,mp,                                  &
             theta_start,theta_end,phi_start,phi_end,                  &
             theta,phi,theta_ind,phi_ind,theta_d,phi_d,ierr)
! **********************************************************************
! Evaluation of 2d-splines (first and second derivatives)
! **********************************************************************
!  print *, 'before eva2d_fd'
  CALL eva2d_fd(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,         &
             b_spl,fderiv)
!  print *, 'before eva2d_sd'
  CALL eva2d_sd(theta_n,phi_n,theta_ind,phi_ind,theta_d,phi_d,         &
             b_spl,sderiv)  
!  print *, 'after eva2d_sd'
! **********************************************************************
! Outputvalues (for neo_zeros2d)
! **********************************************************************
  f    = fderiv(1)
  g    = fderiv(2)
  dfdx = sderiv(1)
  dfdy = sderiv(2)
  dgdx = sderiv(2)
  dgdy = sderiv(3)
END SUBROUTINE neo_bderiv
! **********************************************************************
SUBROUTINE neo_dealloc
! **********************************************************************
! Modules
! **********************************************************************
  USE neo_precision
  USE neo_input
  USE neo_work
  USE neo_exchange
  USE neo_units
  USE neo_parameters
  USE neo_control
  USE neo_spline
  USE neo_van
  USE neo_actual_spectra
  USE neo_spline_data
! **********************************************************************
! Local Definitions
! **********************************************************************
  IMPLICIT NONE
! *******************************************************************
! DeAllocate Storage Arrays
! *******************************************************************
  DEALLOCATE (es,iota,curr_pol,curr_tor,pprime,sqrtg00)
  DEALLOCATE (theta_arr,phi_arr)
  DEALLOCATE (rmnc,zmnc,lmnc,bmnc)
! DEALLOCATE (cosval,sinval)  
  DEALLOCATE (ixm, ixn)
  DEALLOCATE (pixm, pixn)
  DEALLOCATE (i_m, i_n)
  DEALLOCATE (b_spl,k_spl,g_spl,p_spl,r_spl)
  DEALLOCATE (b,sqrg11,kg,pard)
  IF (calc_cur .EQ. 1 .OR. calc_fourier .EQ. 0) THEN
     DEALLOCATE (bqtphi)
  END IF
  IF (calc_cur .EQ. 1) THEN
     DEALLOCATE (q_spl)
  END IF
  DEALLOCATE (fluxs_arr)
  DEALLOCATE (li_minima)
  DEALLOCATE (s_rmnc, s_zmnc, s_lmnc, s_bmnc)

  IF (fluxs_interp .NE. 0) THEN
     DEALLOCATE ( a_rmnc, b_rmnc, c_rmnc, d_rmnc )
     DEALLOCATE ( a_zmnc, b_zmnc, c_zmnc, d_zmnc )
     DEALLOCATE ( a_lmnc, b_lmnc, c_lmnc, d_lmnc )
     DEALLOCATE ( a_bmnc, b_bmnc, c_bmnc, d_bmnc ) 
     
     DEALLOCATE ( a_iota, b_iota )
     DEALLOCATE ( c_iota, d_iota )
     DEALLOCATE ( a_pprime, b_pprime )
     DEALLOCATE ( c_pprime, d_pprime )
     DEALLOCATE ( a_sqrtg00, b_sqrtg00 )
     DEALLOCATE ( c_sqrtg00, d_sqrtg00 )
     DEALLOCATE ( a_curr_tor, b_curr_tor )
     DEALLOCATE ( c_curr_tor, d_curr_tor )
     DEALLOCATE ( a_curr_pol, b_curr_pol )
     DEALLOCATE ( c_curr_pol, d_curr_pol )
 
     DEALLOCATE ( r_m, r_mhalf )
     DEALLOCATE ( sp_index )
  END IF
! *******************************************************************
  RETURN
END SUBROUTINE neo_dealloc

! *******************************************************************
SUBROUTINE neo_spline_test
! Test of Spline Routine
! **********************************************************************
! Modules
! **********************************************************************
  USE neo_precision
  USE neo_input
  USE neo_work
  USE neo_units
  USE neo_parameters
  USE neo_control
  USE neo_eval_switch  
! **********************************************************************
! Local Definitions
! **********************************************************************
  IMPLICIT NONE
! **********************************************************************
  INTEGER            :: i, j, n
  REAL(kind=dp)      :: theta, phi, bval, gval, kval, pval, qval, rval
  INTEGER, PARAMETER :: div = 5, ip = 43, it = 88
  REAL(kind=dp)      :: td, pd 
! **********************************************************************
  eval_switch = (/1,1,1,1,0,0/)
  n = MIN(theta_n,phi_n)

  IF (spline_test .EQ. 1) THEN ! along given phi
     OPEN(unit=w_u1,file='sptest1.dat',status='unknown',form='formatted')
     OPEN(unit=w_u2,file='sptest2.dat',status='unknown',form='formatted')
     td = theta_end - theta_start
     phi = phi_arr(ip)
     theta = theta_start-td 
     DO WHILE (theta < theta_end+td)
        CALL neo_eval(theta,phi,bval,gval,kval,pval,qval,rval)
        WRITE(w_u1,*) theta,phi,bval,gval,kval
        theta = theta + theta_int/div
     END DO
     DO j = -1,1
        DO i = 1, theta_n-1
           theta = theta_arr(i) + j*td
           bval = b(i,ip)
           gval = sqrg11(i,ip)
           kval = kg(i,ip)
           WRITE(w_u2,*) theta,phi,bval,gval,kval
        END DO
     END DO
     CLOSE(unit=w_u2)
     CLOSE(unit=w_u1)
  ELSEIF (spline_test .EQ. 2) THEN ! along given theta
     OPEN(unit=w_u1,file='sptest1.dat',status='unknown',form='formatted')
     OPEN(unit=w_u2,file='sptest2.dat',status='unknown',form='formatted')
     pd = phi_end - phi_start
     theta = theta_arr(it)
     phi = phi_start-pd
     DO WHILE (phi < phi_end+pd)
        CALL neo_eval(theta,phi,bval,gval,kval,pval,qval,rval)
        WRITE(w_u1,*) theta,phi,bval,gval,kval
        phi = phi + phi_int/div
     END DO
     DO j = -1,1
        DO i = 1, phi_n-1
           phi = phi_arr(i) + j*pd
           bval = b(it,i)
           gval = sqrg11(it,i)
           kval = kg(it,i)
           WRITE(w_u2,*) theta,phi,bval,gval,kval
        END DO
     END DO
     CLOSE(unit=w_u2)
     CLOSE(unit=w_u1)
  ELSEIF (spline_test .EQ. 3) THEN ! diagonal
     OPEN(unit=w_u1,file='sptest1.dat',status='unknown',form='formatted')
     DO i = 1, n*div
        theta = theta_start + i*theta_int/div
        phi   = phi_start + i*phi_int/div
        CALL neo_eval(theta,phi,bval,gval,kval,pval,qval,rval)
        WRITE(w_u1,*) theta,phi,bval,gval,kval
     END DO
     CLOSE(unit=w_u1)
     OPEN(unit=w_u1,file='sptest2.dat',status='unknown',form='formatted')
     DO i = 1, n
        theta = theta_arr(i)
        phi   = phi_arr(i)
        bval = b(i,i)
        gval = sqrg11(i,i)
        kval = kg(i,i)
        WRITE(w_u1,*) theta,phi,bval,gval,kval
     END DO
     CLOSE(unit=w_u1)
  ENDIF
! **********************************************************************
  RETURN
END SUBROUTINE neo_spline_test
! **********************************************************************
SUBROUTINE neo_zeros2d(x, y, eps, iter_ma, iter, error)
!
  USE neo_precision

  IMPLICIT NONE

  INTEGER,        INTENT(out)   :: error, iter
  INTEGER,        INTENT(in)    :: iter_ma
  REAL (kind=dp), INTENT(in)    :: eps
  REAL (kind=dp), INTENT(inout) :: x, y

  REAL (kind=dp) :: x_n, y_n
  REAL (kind=dp) :: f, dfdx, dfdy, g, dgdx, dgdy
  REAL (kind=dp) :: f_n,g_n
  REAL (kind=dp) :: det
  REAL (kind=dp) :: x_err, y_err

  error = 0

! compute f(x,y), g(x,y) and all first derivatives
  CALL neo_bderiv(x, y, f, g, dfdx, dfdy, dgdx, dgdy)
  DO iter = 1, iter_ma

     det = dfdx * dgdy - dfdy * dgdx
!!$     PRINT *, 'f,    g    ', f, g
!!$     PRINT *, 'dfdx, dgdx ', dfdx, dgdx
!!$     PRINT *, 'dfdy, dgdy ', dfdy, dgdy
!!$     PRINT *, 'det        ', det
     IF (det .NE. 0.0_dp) THEN
        x_n = x + ( dfdy *  g   -  f   * dgdy ) / det
        y_n = y + (  f   * dgdx - dfdx *  g   ) / det
     ELSE
        x_n = x - f / dfdx
        y_n = y
     END IF
!!$     PRINT *, 'x,    y    ', x, y
!!$     PRINT *, 'x_n,  y_n  ', x_n, y_n

!    compute f(x,y), g(x,y) and all first derivatives
!    at the new positions
!     PRINT *, x_n, y_n
     CALL neo_bderiv(x_n, y_n, f_n, g_n, dfdx, dfdy, dgdx, dgdy)

!    remaining relatve errors
!     IF (x_n .NE. 0.0d0) THEN
     IF (ABS(x_n) .GT. eps) THEN
       x_err = ABS ( (x_n - x) / x_n )
     ELSE
       x_err = ABS ( x_n - x )
     END IF
!     IF (y_n .NE. 0.0d0) THEN
     IF (ABS(y_n) .GT. eps) THEN
       y_err = ABS ( (y_n - y) / y_n )
     ELSE
       y_err = ABS ( y_n - y )
     END IF

!    new values
     f = f_n
     g = g_n
     x = x_n
     y = y_n

!    exit if error is small enough
     IF ( MAX ( x_err, y_err ) < eps ) RETURN

  END DO

  error = 1

  RETURN
END SUBROUTINE neo_zeros2d
! **********************************************************************

SUBROUTINE neo_filenames

  USE neo_precision
  USE neo_units
  USE neo_support
  USE neo_control
  USE neo_exchange

  IMPLICIT NONE
  CHARACTER(len=50)   :: o_file1, o_file2
! **********************************************************************
! File handling
! **********************************************************************
  CALL strip_extension(in_file,'bc',o_file1)
  IF (fluxs_interp .NE. 0) THEN
     CALL add_extension(o_file1,'sp',o_file2)
  ELSE
     o_file2 = o_file1
  END IF

  base_file = o_file2

  CALL add_extension(o_file2,'chk.dat',chk_file)
  IF (chk_swi .EQ. 1) THEN
     OPEN(unit=w_u3, file=chk_file, status='replace', action='write')
     CLOSE(unit=w_u3)
  END IF
  CALL add_extension(o_file2,'dat',out_file)
  CALL add_extension(o_file2,'log.dat',epslog_file)
  CALL add_extension(o_file2,'con.dat',epscon_file)
  CALL add_extension(o_file2,'dia.dat',epsdia_file)
  CALL add_extension(o_file2,'add.dat',epsadd_file)
  CALL add_extension(o_file2,'cur.con.dat',curcon_file)
  CALL add_extension(o_file2,'cur.int.dat',curint_file)
  CALL add_extension(o_file2,'cur.dis.dat',curdis_file)
  IF (calc_eps .EQ. 1) THEN
     OPEN(unit=w_u3, file=out_file, status='replace', action='write')
     CLOSE(unit=w_u3)
     OPEN(unit=w_u3, file=epslog_file, status='replace', action='write')
     CLOSE(unit=w_u3)
  END IF

  CALL add_extension(o_file2,'cur.dat',cur_file)
  IF (calc_cur .EQ. 1) THEN
     OPEN(unit=w_u9, file=cur_file, status='replace', action='write')
     CLOSE(unit=w_u9)
     IF (write_cur_inte .EQ. 1) THEN
        OPEN(unit=w_u9, file=curcon_file, status='replace', action='write')
        CLOSE(unit=w_u9)
        OPEN(unit=w_u9, file=curint_file, status='replace', action='write')
        CLOSE(unit=w_u9)
     END IF
     IF (write_cur_disp .EQ. 1) THEN
        OPEN(unit=w_u9, file=curdis_file, status='replace', action='write')
        CLOSE(unit=w_u9)
     END IF
  END IF

  CALL add_extension(o_file2,'pla.dat',pla_file)
  IF (calc_pla .EQ. 1) THEN
     OPEN(unit=w_u9, file=pla_file, status='replace', action='write')
     CLOSE(unit=w_u9)
  END IF

  CALL add_extension(o_file1,'s.bc',sbc_file)

END SUBROUTINE neo_filenames


SUBROUTINE neo_write_bc
! Write Standard Boozer Files (W7-x)
!***********************************************************************
! Modules
!***********************************************************************
  USE neo_precision
  USE neo_input
  USE neo_units
  USE neo_exchange
  USE neo_control
!***********************************************************************
! Local definitions
!***********************************************************************
  IMPLICIT NONE
 
  REAL(kind=dp)              :: r_small,r_big
  INTEGER                    :: m0b_out,n0b_out
  INTEGER                    :: ii,i,j

  m0b_out = m_max -1
  n0b_out = (n_max-1) / 2
  mnmax = m_max*n_max

  r_small = 0.0_dp
  r_big   = 0.0_dp
  OPEN(unit=w_u1, file=sbc_file, status='replace', action='write')

  WRITE (w_u1,*) ' m0b  n0b nsurf nper flux/[Tm^2]     a/[m]     R/[m]'
  WRITE (w_u1,'(4(i5),3(e13.4))')                                      &
       m0b_out,n0b_out,no_fluxs,nfp,flux,r_small,r_big
  DO ii = 1, no_fluxs
     i = fluxs_arr(ii)
     WRITE(w_u1,*) '     s          iota      curr_pol    curr_tor      pprime  sqrt g(0,0)'
     WRITE(w_u1,'(6(e13.4))')                                          &
          es(i),                                                       &
          iota(i),                                                     &
          curr_pol(i) / 2.d-7 / nfp,                                   &
          curr_tor(i) / 2.d-7,                                         &
          pprime(i),                                                   &
          sqrtg00(i)
     WRITE(w_u1,*) '    m    n        r/[m]           z/[m] (phib-phi)*nper/twopi     bmn/[T]'
     DO j=1,mnmax
        WRITE(w_u1,'(2(i5),4(e17.8))')                                 &
             ixm(j), ixn(j)/nfp,                                       &
             rmnc(i,j), zmnc(i,j), lmnc(i,j),                          &
             bmnc(i,j)
     END DO
  END DO
  CLOSE(w_u1)


END SUBROUTINE neo_write_bc

END MODULE neo_sub_mod







