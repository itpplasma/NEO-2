MODULE neo_magfie_perturbation
  !
  ! module containing numerical constants
  USE neo_precision
  ! module containing switches from the input file (neo.in)
  USE neo_control, ONLY: lab_swi, inp_swi
  ! interface to spline routines
  USE inter_interfaces, ONLY: splinecof3_hi_driv,&
       tf, tfp, tfpp, tfppp, splint_horner3
  ! get boozer_iota (from the axiymmetric file - is the same)
  USE neo_magfie_mod, ONLY: boozer_iota, compute_RZ
  ! magfie
  USE magfie_mod, ONLY : magfie
  ! used for normalization (hat-quantities)
  USE partpa_mod,  ONLY : bmod0
  !! Modification by Andreas F. Martitsch (17.07.2014)
  ! Extra input for NTV computations
  USE ntv_mod, ONLY : in_file_pert, m_phi
  !! End Modification by Andreas F. Martitsch (17.07.2014)
  ! coordinates of starting point
  USE mag_interface_mod, ONLY : boozer_s, boozer_theta_beg,&
       boozer_phi_beg
  ! routine mag for the computation of bmod
  USE mag_sub, ONLY: mag
  USE mpiprovider_module
  !
  IMPLICIT NONE
  !
  ! define kind of double complex, quad, quad complex
  PRIVATE dcp, qp, qcp
  INTEGER, PARAMETER :: dcp=KIND((1.0_dp,1.0_dp))
  INTEGER, PARAMETER :: qp=SELECTED_REAL_KIND(33, 4931)
  INTEGER, PARAMETER :: qcp=KIND((1.0_qp,1.0_qp))
  !
  ! internal (private) storage arrays for the
  ! Fourier spectra from the input file
  PRIVATE mnmax_pert, ns_pert
  INTEGER :: mnmax_pert, ns_pert
  PRIVATE ixm_pert, ixn_pert
  INTEGER, DIMENSION(:), ALLOCATABLE :: ixm_pert, ixn_pert
  ! arrays for s-values and Fourier modes
  PRIVATE es_pert
  PRIVATE rmnc_pert, zmnc_pert, lmnc_pert, bmnc_pert
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: es_pert
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: rmnc_pert
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: zmnc_pert
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: lmnc_pert
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: bmnc_pert
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: rmns_pert
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: zmns_pert
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: lmns_pert
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: bmns_pert
  !
  ! storage arrays for the splines of the 2d
  ! Fourier arrays
  !PRIVATE a_rmnc_pert, b_rmnc_pert, c_rmnc_pert, d_rmnc_pert
  !PRIVATE a_zmnc_pert, b_zmnc_pert, c_zmnc_pert, d_zmnc_pert
  !PRIVATE a_lmnc_pert, b_lmnc_pert, c_lmnc_pert, d_lmnc_pert
  PRIVATE a_bmnc_pert, b_bmnc_pert, c_bmnc_pert, d_bmnc_pert
  !PRIVATE a_rmns_pert, b_rmns_pert, c_rmns_pert, d_rmns_pert
  !PRIVATE a_zmns_pert, b_zmns_pert, c_zmns_pert, d_zmns_pert
  !PRIVATE a_lmns_pert, b_lmns_pert, c_lmns_pert, d_lmns_pert
  PRIVATE a_bmns_pert, b_bmns_pert, c_bmns_pert, d_bmns_pert
  !REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: a_rmnc_pert, b_rmnc_pert
  !REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: c_rmnc_pert, d_rmnc_pert
  !REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: a_zmnc_pert, b_zmnc_pert
  !REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: c_zmnc_pert, d_zmnc_pert
  !REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: a_lmnc_pert, b_lmnc_pert
  !REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: c_lmnc_pert, d_lmnc_pert
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: a_bmnc_pert, b_bmnc_pert
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: c_bmnc_pert, d_bmnc_pert
  !REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: a_rmns_pert, b_rmns_pert
  !REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: c_rmns_pert, d_rmns_pert
  !REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: a_zmns_pert, b_zmns_pert
  !REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: c_zmns_pert, d_zmns_pert
  !REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: a_lmns_pert, b_lmns_pert
  !REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: c_lmns_pert, d_lmns_pert
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: a_bmns_pert, b_bmns_pert
  REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: c_bmns_pert, d_bmns_pert
  !
  ! storage arrays for r_mhalf and sp_index
  ! requested by the spline routines
  PRIVATE r_m_pert, r_mhalf_pert, sp_index_pert
  REAL(kind=dp), DIMENSION(:),   ALLOCATABLE :: r_m_pert, r_mhalf_pert
  INTEGER(I4B),  DIMENSION(:),   ALLOCATABLE :: sp_index_pert
  !
  PUBLIC neo_magfie_pert
  PRIVATE neo_magfie_pert_a, neo_magfie_pert_b
  INTERFACE neo_magfie_pert
     MODULE PROCEDURE neo_magfie_pert_a, neo_magfie_pert_b
  END INTERFACE neo_magfie_pert
  !
CONTAINS
  !
  ! Read Boozer files for the perturbation field
  SUBROUTINE neo_read_pert()
    !
    ! local definitions
    !
    ! specify mode numbers and number of flux surfaces
    INTEGER :: m0b_pert, n0b_pert
    INTEGER :: m_max_pert, n_max_pert
    INTEGER :: nfp_pert
    ! loop counters, status variables
    INTEGER :: i,j
    INTEGER :: i_alloc, r_un
    ! dummy variables
    INTEGER :: int_dummy
    REAL(kind=dp) :: real_dummy
    CHARACTER(5) :: char_dummy
    !
    ! open Boozer file and read first quantities
    r_un=27052014
    in_file_pert=TRIM(ADJUSTL(in_file_pert))
    OPEN(unit=r_un,file=in_file_pert,status='old',form='formatted')
    !
    IF (inp_swi .EQ. 8) THEN ! NEW IPP TOKAMAK        
       READ (r_un,*) char_dummy
       READ (r_un,*) char_dummy
       READ (r_un,*) char_dummy
       READ (r_un,*) char_dummy
       READ (r_un,*) char_dummy
       READ (r_un,*) m0b_pert,n0b_pert,ns_pert,&
            nfp_pert,real_dummy,real_dummy,real_dummy
       m_max_pert = m0b_pert+1
       n_max_pert = n0b_pert+1
       mnmax_pert = m_max_pert*n_max_pert
       !
       ! allocate storage arrays
       ALLOCATE(ixm_pert(mnmax_pert), ixn_pert(mnmax_pert), stat = i_alloc)
       IF(i_alloc /= 0) THEN
          STOP "Allocation for the arrays containing the mode numbers&
               & of the perturbation field failed!"
       END IF
       !
       ALLOCATE(es_pert(ns_pert), stat = i_alloc)
       IF(i_alloc /= 0) THEN
          STOP "Allocation for the real array containing&
               & the s-values of the perturbation field failed!"
       END IF
       !
       ALLOCATE(rmnc_pert(ns_pert,mnmax_pert), zmnc_pert(ns_pert,mnmax_pert),&
            lmnc_pert(ns_pert,mnmax_pert), bmnc_pert(ns_pert,mnmax_pert), stat = i_alloc)
       IF(i_alloc /= 0) THEN
          STOP "Allocation for the Fourier arrays for the perturbation field failed!"
       END IF
       !
       ! read input arrays
       DO i =1, ns_pert
          READ(r_un,*) char_dummy
          READ(r_un,*) char_dummy
          READ(r_un,*) es_pert(i),real_dummy,real_dummy,real_dummy,real_dummy,real_dummy
          READ(r_un,*) char_dummy

          DO j=1,mnmax_pert
             !print *, 'j: ',j
             READ(r_un,*) ixm_pert(j),ixn_pert(j),&
                  rmnc_pert(i,j),zmnc_pert(i,j),lmnc_pert(i,j),&
                  bmnc_pert(i,j)
             !PRINT *, 'ixm,ixn,bmnc: ',ixm_pert(j),ixn_pert(j),bmnc_pert(i,j)
          END DO
       END DO
    ELSEIF (inp_swi .EQ. 9) THEN ! ASDEX-U (E. Strumberger)     
       READ (r_un,*) char_dummy
       READ (r_un,*) char_dummy
       READ (r_un,*) char_dummy
       READ (r_un,*) char_dummy
       READ (r_un,*) char_dummy
       READ (r_un,*) m0b_pert,n0b_pert,ns_pert,&
            nfp_pert,real_dummy,real_dummy,real_dummy
       m_max_pert = m0b_pert+1
       n_max_pert = n0b_pert+1
       mnmax_pert = m_max_pert*n_max_pert
       !
       ! allocate storage arrays
       ALLOCATE(ixm_pert(mnmax_pert), ixn_pert(mnmax_pert), stat = i_alloc)
       IF(i_alloc /= 0) THEN
          STOP "Allocation for the arrays containing the mode numbers&
               & of the perturbation field failed!"
       END IF
       !
       ALLOCATE(es_pert(ns_pert), stat = i_alloc)
       IF(i_alloc /= 0) THEN
          STOP "Allocation for the real array containing&
               & the s-values of the perturbation field failed!"
       END IF
       !
       ALLOCATE(rmnc_pert(ns_pert,mnmax_pert), zmnc_pert(ns_pert,mnmax_pert),&
            lmnc_pert(ns_pert,mnmax_pert), bmnc_pert(ns_pert,mnmax_pert), stat = i_alloc)
       IF(i_alloc /= 0) THEN
          STOP "Allocation for the Fourier arrays for the perturbation field failed!"
       END IF
       ALLOCATE(rmns_pert(ns_pert,mnmax_pert), zmns_pert(ns_pert,mnmax_pert),&
            lmns_pert(ns_pert,mnmax_pert), bmns_pert(ns_pert,mnmax_pert), stat = i_alloc)
       IF(i_alloc /= 0) THEN
          STOP "Allocation for the Fourier arrays for the perturbation field failed!"
       END IF
       !
       ! read input arrays
       DO i =1, ns_pert
          READ(r_un,*) char_dummy
          READ(r_un,*) char_dummy
          READ(r_un,*) es_pert(i),real_dummy,real_dummy,real_dummy,real_dummy,real_dummy
          READ(r_un,*) char_dummy

          DO j=1,mnmax_pert
             !print *, 'j: ',j
             READ(r_un,*) ixm_pert(j),ixn_pert(j),&
                  rmnc_pert(i,j),rmns_pert(i,j),zmnc_pert(i,j),zmns_pert(i,j),&
                  lmnc_pert(i,j),lmns_pert(i,j),bmnc_pert(i,j),bmns_pert(i,j)
             !PRINT *, 'ixm,ixn,bmnc: ',ixm_pert(j),ixn_pert(j),bmnc_pert(i,j),bmns_pert(i,j),es_pert(i)
          END DO
       END DO
    ELSE
       PRINT *,'FATAL: There is yet no other input type for the perturbed field defined'
       STOP
    END IF
    !
    IF (lab_swi .EQ. 8) THEN ! NEW IPP TOKAMAK
       ixn_pert =  ixn_pert * nfp_pert
       m_phi = ixn_pert(1)
       ixm_pert =  ixm_pert
    ELSEIF ((lab_swi .EQ. 9) .OR. (lab_swi .EQ. 10)) THEN ! ASDEX-U (E. Strumberger)
          ixn_pert =  ixn_pert * nfp_pert
          m_phi = ixn_pert(1)
          ixm_pert =  ixm_pert
    ELSE
       PRINT *,'FATAL: There is yet no other Laboratory for the perturbed field defined!'
       STOP
    END IF
    ! close Boozer file
    CLOSE (unit=r_un)
    RETURN
  END SUBROUTINE neo_read_pert
  !
  ! Initialization for splines along s
  SUBROUTINE neo_init_spline_pert()
    !
    ! local definitions
    !
    ! loop variables, poloidal mode number
    INTEGER :: i
    INTEGER(I4B) :: m
    ! some maximum value of poloidal mode numbers
    ! (used for the computation of r_mhalf)
    INTEGER, PARAMETER :: m_max_sp = 12
    !
    ! testing
    !
!!$    ! define switch and select a flux surface
!!$    INTEGER(I4B) :: swd
!!$    INTEGER      :: k
!!$    REAL(dp)     :: m0, f_es_pert, f_bmnc_pert, dummy
    !
    ! allocate 2d arrays for the splines of the Fourier coefficients
    !at the moment unused stuff
    !ALLOCATE ( a_rmnc_pert(ns_pert,mnmax_pert), b_rmnc_pert(ns_pert,mnmax_pert) )
    !ALLOCATE ( c_rmnc_pert(ns_pert,mnmax_pert), d_rmnc_pert(ns_pert,mnmax_pert) )
    !ALLOCATE ( a_zmnc_pert(ns_pert,mnmax_pert), b_zmnc_pert(ns_pert,mnmax_pert) )
    !ALLOCATE ( c_zmnc_pert(ns_pert,mnmax_pert), d_zmnc_pert(ns_pert,mnmax_pert) )
    !ALLOCATE ( a_lmnc_pert(ns_pert,mnmax_pert), b_lmnc_pert(ns_pert,mnmax_pert) )
    !ALLOCATE ( c_lmnc_pert(ns_pert,mnmax_pert), d_lmnc_pert(ns_pert,mnmax_pert) )
    ALLOCATE ( a_bmnc_pert(ns_pert,mnmax_pert), b_bmnc_pert(ns_pert,mnmax_pert) ) 
    ALLOCATE ( c_bmnc_pert(ns_pert,mnmax_pert), d_bmnc_pert(ns_pert,mnmax_pert) )
    ! Additional data from Boozer files without Stellarator symmetry
    IF (inp_swi .EQ. 9) THEN        ! ASDEX-U (E. Strumberger)
       !at the moment unused stuff
       !ALLOCATE ( a_rmns_pert(ns_pert,mnmax_pert), b_rmns_pert(ns_pert,mnmax_pert) )
       !ALLOCATE ( c_rmns_pert(ns_pert,mnmax_pert), d_rmns_pert(ns_pert,mnmax_pert) )
       !ALLOCATE ( a_zmns_pert(ns_pert,mnmax_pert), b_zmns_pert(ns_pert,mnmax_pert) )
       !ALLOCATE ( c_zmns_pert(ns_pert,mnmax_pert), d_zmns_pert(ns_pert,mnmax_pert) )
       !ALLOCATE ( a_lmns_pert(ns_pert,mnmax_pert), b_lmns_pert(ns_pert,mnmax_pert) )
       !ALLOCATE ( c_lmns_pert(ns_pert,mnmax_pert), d_lmns_pert(ns_pert,mnmax_pert) )
       ALLOCATE ( a_bmns_pert(ns_pert,mnmax_pert), b_bmns_pert(ns_pert,mnmax_pert) ) 
       ALLOCATE ( c_bmns_pert(ns_pert,mnmax_pert), d_bmns_pert(ns_pert,mnmax_pert) )
    END IF
    ! allocate arrays requested by the spline routines
    ALLOCATE ( r_m_pert(mnmax_pert), r_mhalf_pert(mnmax_pert) )
    ALLOCATE ( sp_index_pert(ns_pert) )
    !
    ! compute r_mhalf_pert and sp_index_pert requested by the spline routines
    ! (r_mhalf_pert enters the test function tf, which is used within the spline
    ! routines - tf(x) = x^r_mhalf_pert, x is here boozer_s. Since mode number m
    ! has to be symmetric around 0 for the representation of the perturbation
    ! field, r_mhalf_pert should be also symmetric -> m = abs(ixm_pert(i)).
    ! Without abs() the value of the test function grows unbounded for
    ! boozer_s << 1 and r_mhalf_pert < 0)
    DO i = 1,mnmax_pert
       ! abs() - because of negative m-values
       ! (not necessary for Boozer files with m>=0, see neo_init_spline/neo_sub.f90)
       m = ABS(ixm_pert(i))
       IF (m .LE. m_max_sp) THEN
          r_m_pert(i)     = DBLE(m)
       ELSE
          IF (MODULO(m,2) .EQ. 1) THEN
             r_m_pert(i)     = DBLE(m_max_sp+1)
          ELSE
             r_m_pert(i)     = DBLE(m_max_sp)
          END IF
       END IF
       r_mhalf_pert(i) = r_m_pert(i) / 2._dp
    END DO
    sp_index_pert = (/ (i, i=1,ns_pert) /) 
    !
    ! 1-d splines of 2-d arrays
    !at the moment unused stuff
    !CALL splinecof3_hi_driv(es_pert, rmnc_pert, r_mhalf_pert,&
    !     a_rmnc_pert, b_rmnc_pert, c_rmnc_pert, d_rmnc_pert, sp_index_pert, tf)
    !CALL splinecof3_hi_driv(es_pert, zmnc_pert, r_mhalf_pert,&
    !     a_zmnc_pert, b_zmnc_pert, c_zmnc_pert, d_zmnc_pert, sp_index_pert, tf)
    !CALL splinecof3_hi_driv(es_pert, lmnc_pert, r_mhalf_pert,&
    !     a_lmnc_pert, b_lmnc_pert, c_lmnc_pert, d_lmnc_pert, sp_index_pert, tf)
    CALL splinecof3_hi_driv(es_pert, bmnc_pert, r_mhalf_pert,&
         a_bmnc_pert, b_bmnc_pert, c_bmnc_pert, d_bmnc_pert, sp_index_pert, tf)
    ! Additional data from Boozer files without Stellarator symmetry
    IF (inp_swi .EQ. 9) THEN        ! ASDEX-U (E. Strumberger)
       !at the moment unused stuff
       !CALL splinecof3_hi_driv(es_pert, rmns_pert, r_mhalf_pert,&
       !     a_rmns_pert, b_rmns_pert, c_rmns_pert, d_rmns_pert, sp_index_pert, tf)
       !CALL splinecof3_hi_driv(es_pert, zmns_pert, r_mhalf_pert,&
       !     a_zmns_pert, b_zmns_pert, c_zmns_pert, d_zmns_pert, sp_index_pert, tf)
       !CALL splinecof3_hi_driv(es_pert, lmns_pert, r_mhalf_pert,&
       !     a_lmns_pert, b_lmns_pert, c_lmns_pert, d_lmns_pert, sp_index_pert, tf)
       CALL splinecof3_hi_driv(es_pert, bmns_pert, r_mhalf_pert,&
            a_bmns_pert, b_bmns_pert, c_bmns_pert, d_bmns_pert, sp_index_pert, tf)
    END IF
    !
    ! Testing
    !
!!$    swd = 0 ! no derivatives
!!$    k = 1
!!$    f_es_pert = es_pert(k)
!!$    DO i=1,mnmax_pert
!!$       m0 = r_mhalf_pert(i)
!!$       CALL splint_horner3(es_pert,a_bmnc_pert(:,i), b_bmnc_pert(:,i), c_bmnc_pert(:,i),&
!!$            d_bmnc_pert(:,i), swd, m0, f_es_pert, tf, tfp, tfpp, tfppp, f_bmnc_pert,&
!!$            dummy,dummy,dummy)
!!$       PRINT *, ixm_pert(i),ixn_pert(i),bmnc_pert(k,i),f_bmnc_pert  
!!$    END DO
    !
    RETURN
  END SUBROUTINE neo_init_spline_pert
  !
  ! Compute the perturbation field for a certain x-value
  SUBROUTINE neo_magfie_pert_a( x, bn_hat_pert )
    !
    ! input / output
    !
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: x
    COMPLEX(kind=dcp), INTENT(out) :: bn_hat_pert
    !
    ! local definitions
    !
    COMPLEX(kind=dcp), PARAMETER :: imun=(0.0_dp,1.0_dp)
    INTEGER(i4b) :: swd
    INTEGER :: i, m, n
    REAL(kind=dp) :: yp, ypp, yppp
    REAL(kind=dp) :: bmnc_pert_val, bmns_pert_val
    COMPLEX(kind=dcp) :: expv
    !
    ! read Boozer file and prepare spline routines (1st call)
    !
    IF (.NOT. ALLOCATED(es_pert)) THEN
       CALL neo_read_pert()
       CALL neo_init_spline_pert()
    END IF
    !
    ! direct summation of Fourier components
    !
    bn_hat_pert = (0.0_dp,0.0_dp)
    !PRINT *,'mnmax_pert: ', mnmax_pert
    DO i = 1, mnmax_pert
       swd = 1
       CALL splint_horner3(es_pert, a_bmnc_pert(:,i), b_bmnc_pert(:,i),&
            c_bmnc_pert(:,i), d_bmnc_pert(:,i), swd, r_mhalf_pert(i),&
            x(1), tf, tfp, tfpp, tfppp, bmnc_pert_val, yp, ypp, yppp)
       ! Additional data from Boozer files without Stellarator symmetry
       IF (inp_swi .EQ. 9) THEN ! ASDEX-U (E. Strumberger)
          CALL splint_horner3(es_pert, a_bmns_pert(:,i), b_bmns_pert(:,i),&
               c_bmns_pert(:,i), d_bmns_pert(:,i), swd, r_mhalf_pert(i),&
               x(1), tf, tfp, tfpp, tfppp, bmns_pert_val, yp, ypp, yppp)
       END IF

       IF (inp_swi .EQ. 8) THEN ! NEW IPP TOKAMAK
          ! requested representation of the perturbation field
          ! $B_n = \sum_{m>-\infty} \tilde{b}_{mn} \exp{i(m\vartheta_B+n\varphi_B)}$
          n = ixn_pert(i)
          !PRINT *,'n: ',n
          ! minus sign due to spectrum of the form
          ! $B = \sum_{mn} b_{mn} cos(m\vartheta_B-n\varphi_B)$
          m = (-1)*ixm_pert(i)
          !PRINT *,'m: ',m

          ! perturbation field is represented in terms of an expansion
          ! over the periodic toroidal Boozer angle times a phase
          ! (Note that x(3) is here a periodic angle - see mag_interface.f90 -
          ! while x(3) is the coordinate along the field line! Make sure that
          ! the different periodicities are treated correctly)
          expv = EXP(imun*(m*x(3)+n*x(2)))
          !PRINT *,'exponential (expv): ',expv

          ! part corresponding to the expansion over the periodic angle
          bn_hat_pert = bn_hat_pert + bmnc_pert_val * expv

          ! alternative definition (for large aspect ratios, there is a factor 2 difference)
          !m = ixm_pert(i)
          !bn_hat_pert = bn_hat_pert + bmnc_pert_val * COS((m*boozer_iota-n)*x(2)) - &
          !     imun * bmnc_pert_val * SIN((m*boozer_iota-n)*x(2))
       ELSEIF (inp_swi .EQ. 9) THEN ! ASDEX-U (E. Strumberger)
          ! requested representation of the perturbation field
          ! $B_n = \sum_{m>-\infty} \tilde{b}_{mn} \exp{i(m\vartheta_B+n\varphi_B)}$
          n = ixn_pert(i)
          m = ixm_pert(i)
          !PRINT *,'m: ',m

          ! perturbation field is represented in terms of an expansion
          ! over the periodic toroidal Boozer angle times a phase
          ! (Note that x(3) is here a periodic angle - see mag_interface.f90 -
          ! while x(3) is the coordinate along the field line! Make sure that
          ! the different periodicities are treated correctly)
          expv = EXP(imun*(m*x(3)+n*x(2)))
          !PRINT *,'exponential (expv): ',expv

          ! part corresponding to the expansion over the periodic angle
          bn_hat_pert = bn_hat_pert + (bmnc_pert_val - imun * bmns_pert_val) * expv

          ! alternative definition
          !bn_hat_pert = bn_hat_pert + bmnc_pert_val * COS((m*boozer_iota+n)*x(2)) + &
          !     bmns_pert_val * SIN((m*boozer_iota+n)*x(2)) + imun * ( bmnc_pert_val * &
          !     SIN((m*boozer_iota+n)*x(2)) - bmns_pert_val * COS((m*boozer_iota+n)*x(2)) )
       END IF

    END DO
    bn_hat_pert = bn_hat_pert / bmod0
    !PRINT *,'bhat, bn_hat_pert: ',(bn_hat_pert/bhat),(1.0d-3*EXP(imun*m_phi*x(2)))
    !PRINT *,'bn_hat_pert: ',bn_hat_pert
    !PRINT *,'bmod0: ',bmod0
  END SUBROUTINE neo_magfie_pert_a
  !
  ! Compute the perturbation field for a certain x-value and its 
  ! derivative over theta
  SUBROUTINE neo_magfie_pert_b( x, bn_hat_pert, dbn_hat_pert_dtheta )
    !
    ! input / output
    !
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: x
    COMPLEX(kind=dcp), INTENT(out) :: bn_hat_pert
    COMPLEX(kind=dcp), INTENT(out) :: dbn_hat_pert_dtheta
    !
    ! local definitions
    !
    COMPLEX(kind=dcp), PARAMETER :: imun=(0.0_dp,1.0_dp)
    INTEGER(i4b) :: swd
    INTEGER :: i, m, n
    REAL(kind=dp) :: yp, ypp, yppp
    REAL(kind=dp) :: bmnc_pert_val, bmns_pert_val
    COMPLEX(kind=dcp) :: expv
    !
    ! read Boozer file and prepare spline routines (1st call)
    !
    IF (.NOT. ALLOCATED(es_pert)) THEN
       CALL neo_read_pert()
       CALL neo_init_spline_pert()
    END IF
    !
    ! direct summation of Fourier components
    !
    bn_hat_pert = (0.0_dp,0.0_dp)
    dbn_hat_pert_dtheta = (0.0_dp,0.0_dp)
    !PRINT *,'mnmax_pert: ', mnmax_pert
    DO i = 1, mnmax_pert
       swd = 1
       CALL splint_horner3(es_pert, a_bmnc_pert(:,i), b_bmnc_pert(:,i),&
            c_bmnc_pert(:,i), d_bmnc_pert(:,i), swd, r_mhalf_pert(i),&
            x(1), tf, tfp, tfpp, tfppp, bmnc_pert_val, yp, ypp, yppp)
       ! Additional data from Boozer files without Stellarator symmetry
       IF (inp_swi .EQ. 9) THEN ! ASDEX-U (E. Strumberger)
          CALL splint_horner3(es_pert, a_bmns_pert(:,i), b_bmns_pert(:,i),&
               c_bmns_pert(:,i), d_bmns_pert(:,i), swd, r_mhalf_pert(i),&
               x(1), tf, tfp, tfpp, tfppp, bmns_pert_val, yp, ypp, yppp)
       END IF

       IF (inp_swi .EQ. 8) THEN ! NEW IPP TOKAMAK
          ! requested representation of the perturbation field
          ! $B_n = \sum_{m>-\infty} \tilde{b}_{mn} \exp{i(m\vartheta_B+n\varphi_B)}$
          n = ixn_pert(i)
          !PRINT *,'n: ',n
          ! minus sign due to spectrum of the form
          ! $B = \sum_{mn} b_{mn} cos(m\vartheta_B-n\varphi_B)$
          m = (-1)*ixm_pert(i)
          !PRINT *,'m: ',m

          ! perturbation field is represented in terms of an expansion
          ! over the periodic toroidal Boozer angle times a phase
          ! (Note that x(3) is here a periodic angle - see mag_interface.f90 -
          ! while x(3) is the coordinate along the field line! Make sure that
          ! the different periodicities are treated correctly)
          expv = EXP(imun*(m*x(3)+n*x(2)))
          !PRINT *,'exponential (expv): ',expv

          ! part corresponding to the expansion over the periodic angle
          bn_hat_pert = bn_hat_pert + bmnc_pert_val * expv
          dbn_hat_pert_dtheta = dbn_hat_pert_dtheta + imun * m * bmnc_pert_val * expv
       ELSEIF (inp_swi .EQ. 9) THEN ! ASDEX-U (E. Strumberger)
          ! requested representation of the perturbation field
          ! $B_n = \sum_{m>-\infty} \tilde{b}_{mn} \exp{i(m\vartheta_B+n\varphi_B)}$
          n = ixn_pert(i)
          m = ixm_pert(i)
          !PRINT *,'m: ',m

          ! perturbation field is represented in terms of an expansion
          ! over the periodic toroidal Boozer angle times a phase
          ! (Note that x(3) is here a periodic angle - see mag_interface.f90 -
          ! while x(3) is the coordinate along the field line! Make sure that
          ! the different periodicities are treated correctly)
          expv = EXP(imun*(m*x(3)+n*x(2)))
          !PRINT *,'exponential (expv): ',expv

          ! part corresponding to the expansion over the periodic angle
          bn_hat_pert = bn_hat_pert + (bmnc_pert_val - imun * bmns_pert_val) * expv
          !PRINT *,'m, n, expv, bn_hat_pert, bmnc, bmns: ',m,n,expv,bn_hat_pert,bmnc_pert_val,bmns_pert_val
          dbn_hat_pert_dtheta = dbn_hat_pert_dtheta +&
               imun * m * (bmnc_pert_val - imun * bmns_pert_val) * expv
       END IF

    END DO
    bn_hat_pert = bn_hat_pert / bmod0
    dbn_hat_pert_dtheta = dbn_hat_pert_dtheta / bmod0
    !PRINT *,'bhat, bn_hat_pert: ',(bn_hat_pert/bhat),(1.0d-3*EXP(imun*m_phi*x(2)))
    !PRINT *,'bn_hat_pert: ',bn_hat_pert
    !PRINT *,'bmod0: ',bmod0
  END SUBROUTINE neo_magfie_pert_b
  !
  SUBROUTINE calc_bnoverb0_arr(phi_arr,ibeg,iend,bnoverb0_arr,dbnoverb0_dphi_arr)
    !
    ! input / output
    !
    INTEGER, INTENT(in) :: ibeg, iend
    REAL(kind=dp), DIMENSION(ibeg:iend), INTENT(in) :: phi_arr
    COMPLEX(kind=dcp), DIMENSION(ibeg:iend), INTENT(inout) :: bnoverb0_arr
    COMPLEX(kind=dcp), DIMENSION(ibeg:iend), INTENT(inout) :: dbnoverb0_dphi_arr
    !
    ! local definitions
    !
    COMPLEX(kind=dcp), PARAMETER :: imun=(0.0_dp,1.0_dp)
    COMPLEX(kind=dcp) :: bn_hat_pert, dbn_hat_pert_dtheta
    REAL(kind=dp) :: phi_val, theta_val, bmod, sqrtg
    REAL(kind=dp), DIMENSION(3) :: x, bder, hcovar, hctrvr
    REAL(kind=dp), DIMENSION(3,3) :: hcoder,hctder
    INTEGER :: ind_arr
    !
    DO ind_arr = ibeg,iend
       !
       phi_val=phi_arr(ind_arr)
       theta_val=boozer_theta_beg+(phi_val-boozer_phi_beg)*boozer_iota
       IF (theta_val .GT. 2.0_dp*PI) theta_val = theta_val - 2.0_dp*PI
       x = (/boozer_s,phi_val,theta_val/)
       !
       CALL mag(x,bmod,sqrtg,bder,hcovar,hctrvr,hcoder,hctder)
       CALL neo_magfie_pert(x,bn_hat_pert,dbn_hat_pert_dtheta)
       !
       bnoverb0_arr(ind_arr) = (bn_hat_pert * bmod0) / bmod
       dbnoverb0_dphi_arr(ind_arr) = &
            (dbn_hat_pert_dtheta * bmod0 / bmod) * boozer_iota + &
            bnoverb0_arr(ind_arr) * (imun * m_phi - bder(3) * boozer_iota)
       !
    END DO
    !
  END SUBROUTINE calc_bnoverb0_arr
  !
  SUBROUTINE calc_ntv_output(phi_arr,bhat_arr,bnoverb0_arr,ibeg,iend,eps_M_2,av_inv_bhat,av_gphph)
    !
    ! input / output
    !
    INTEGER, INTENT(in) :: ibeg, iend
    REAL(kind=dp), DIMENSION(ibeg:iend), INTENT(in) :: phi_arr, bhat_arr
    COMPLEX(kind=dcp), DIMENSION(ibeg:iend), INTENT(in) :: bnoverb0_arr
    REAL(kind=dp), INTENT(out) :: eps_M_2, av_gphph, av_inv_bhat
    !
    ! local definitions
    !
    INTEGER :: k
    REAL(kind=dp) :: bnoverb0_k, bnoverb0_kp1
    REAL(kind=dp), DIMENSION(3) :: x
    COMPLEX(kind=dcp), PARAMETER :: imun=(0.0_dp,1.0_dp)
    !
    ! compute $\bar{1/\hat{B}}$
    av_inv_bhat=0.0_dp
    DO k = ibeg,iend-1
       av_inv_bhat = av_inv_bhat + (phi_arr(k+1)-phi_arr(k)) * &
            ( (1.0_dp/(bhat_arr(k+1)**2.0_dp)) + (1.0_dp/(bhat_arr(k)**2.0_dp)) )
    END DO
    av_inv_bhat = 0.5_dp * av_inv_bhat / (phi_arr(iend)-phi_arr(ibeg))
    !
    ! compute $\bar{ (({\delta}B/B)^2) * (1/\hat{B}) }$
    eps_M_2 = 0.0_dp
    DO k = ibeg,iend-1
       bnoverb0_k = (REAL( bnoverb0_arr(k)*EXP(-imun*m_phi*phi_arr(k)) ))**2.0_dp +&
            (AIMAG( bnoverb0_arr(k)*EXP(-imun*m_phi*phi_arr(k)) ))**2.0_dp
       bnoverb0_kp1 = (REAL( bnoverb0_arr(k+1)*EXP(-imun*m_phi*phi_arr(k+1)) ))**2.0_dp +&
            (AIMAG( bnoverb0_arr(k+1)*EXP(-imun*m_phi*phi_arr(k+1)) ))**2.0_dp
       ! simplified
       !eps_M_2 = eps_M_2 + (phi_arr(k+1)-phi_arr(k)) * &
       !     ( bnoverb0_kp1 + bnoverb0_k )
       ! full geometry
       eps_M_2 = eps_M_2 + (phi_arr(k+1)-phi_arr(k)) * &
            ( (bnoverb0_kp1/(bhat_arr(k+1)**3.0_dp)) + (bnoverb0_k/(bhat_arr(k)**3.0_dp)) )
    END DO
    eps_M_2 = 0.5_dp * eps_M_2 / (phi_arr(iend)-phi_arr(ibeg))
    eps_M_2 = eps_M_2 / av_inv_bhat
    !
    ! compute $\langle g_{\varphi\varphi} \rangle$
    x = (/boozer_s,boozer_phi_beg,boozer_theta_beg/)
    CALL calc_av_gphph(x,phi_arr,ibeg,iend,av_gphph)
    !
  END SUBROUTINE calc_ntv_output
  !
  ! compute the flux surface average of $g_{\varphi\varphi}$
  ! for symmetry flux coordinates
  SUBROUTINE calc_av_gphph(x_start,phi_arr,ibeg,iend,av_gphph)
    !
    ! input / output
    !
    INTEGER, INTENT(in) :: ibeg, iend
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: x_start
    REAL(kind=dp), DIMENSION(ibeg:iend), INTENT(in) :: phi_arr
    REAL(kind=dp), INTENT(out) :: av_gphph
    !
    ! local definitions:
    !
    INTEGER :: k, w_un
    REAL(kind=dp) :: bmoda, sqrtg, R_val, Z_val, fac1, fac2
    REAL(kind=dp), DIMENSION(3) :: x, bder, hcovar, hctrvr
    REAL(kind=dp), DIMENSION(3) :: hcurl, bcovar_s_hat_der
    REAL(kind=dp), DIMENSION(ibeg:iend) :: R_arr, bmod_arr
    !
    ! compute R, bmod as a function of phi
    w_un=14112014
    IF ( mpro%isMaster() ) &
      OPEN(unit=w_un,file='fluxsurface.dat',status='replace',form='formatted')
    DO k = ibeg,iend
       x(1)=x_start(1)
       x(2)=phi_arr(k)
       x(3)=x_start(3)+(x(2)-x_start(2))*boozer_iota
       CALL magfie( x, bmoda, sqrtg, bder, hcovar, hctrvr,&
            hcurl, bcovar_s_hat_der, R_val, Z_val )
       R_arr(k) = R_val
       bmod_arr(k) = bmoda
       IF ( mpro%isMaster() ) WRITE(w_un,*) x(2), R_val, Z_val, bmoda
    END DO
    IF ( mpro%isMaster() ) CLOSE(unit=w_un)
    !
    ! evaluate flux surface average
    fac1 = 0.0_dp 
    fac2 = 0.0_dp 
    DO k = ibeg,iend-1
       fac1 = fac1 + (phi_arr(k+1)-phi_arr(k)) * &
            ( (R_arr(k+1)/bmod_arr(k+1))**2.0_dp + (R_arr(k)/bmod_arr(k))**2.0_dp )
       fac2 = fac2 + (phi_arr(k+1)-phi_arr(k)) * &
            ( (1.0_dp/bmod_arr(k+1))**2.0_dp + (1.0_dp/bmod_arr(k))**2.0_dp )
       
    END DO
    fac1 = 0.5_dp * fac1
    fac2 = 0.5_dp * fac2
    av_gphph = fac1 / fac2
    !PRINT *,'av_gphph: ',av_gphph
    !
  END SUBROUTINE calc_av_gphph
!
END MODULE neo_magfie_perturbation
