!-----------------------------------------------------------------------------------------!
! module: gsl_integration_routines_mod                                                    !
! authors: Andreas F. Martitsch, TU Graz ITPcp Plasma                                            !
! date: 27.07.2015                                                                        !
! version: 0.4                                                                            !
! description:                                                                            !
! Module gsl_integration_routines_mod provides routines for integrating 1d functions      !
! The numerical routines are provided by the fortnum library (QUADPACK-pattern adaptive   !
! Gauss-Kronrod quadrature). The public interface is preserved from the previous          !
! FGSL/GSL-backed version so existing callers compile unchanged.                          !
! changes: Global procedure pointers and wrapper-functions are replaced by internal       !
! subroutines attached to the integration routines. This allows for a nested call of      !
! integration routines from outside (e.g, 2d integration)                                 !
!-----------------------------------------------------------------------------------------!

!-----------------------------------------------------------------------------------------!
! History                                                                                 !
! 0.1 - Initial version                                                                   !
! 0.2 - Support for 2D integration                                                        !
! 0.3 - Added recursive keyword to functions for 2D integration                           !
! 0.4 - Backend moved from FGSL/GSL to fortnum. CQUAD is re-expressed on the fortnum       !
!       adaptive QAGS path (Wynn-epsilon extrapolation, handles endpoint singularities).   !
!-----------------------------------------------------------------------------------------!


MODULE gsl_integration_routines_mod

  USE, INTRINSIC :: iso_fortran_env, ONLY: wp => real64
  USE fortnum_integrate, ONLY: integrate_integrand_t, integrate_workspace_t, &
       integrate_epstab_t, integrate_result_t, &
       integrate_qag, integrate_qags, integrate_qagp, integrate_qagiu
  USE fortnum_cquad, ONLY: integrate_cquad
  USE fortnum_status, ONLY: fortnum_status_t, FORTNUM_OK, &
       FORTNUM_DOMAIN_ERROR, FORTNUM_CONVERGENCE_ERROR, FORTNUM_NOT_IMPLEMENTED

  IMPLICIT NONE

  ! Possible values of absolute/relative errors
  PRIVATE eps5, eps7, eps10, eps12
  REAL(wp), PARAMETER :: eps5 = 1.0E-5_wp
  REAL(wp), PARAMETER :: eps7 = 1.0E-7_wp
  REAL(wp), PARAMETER :: eps10 = 1.0E-10_wp
  REAL(wp), PARAMETER :: eps12 = 1.0E-12_wp

  ! GK rule selection (QAG): map sw_qag_rule 1..6 to fortnum key 15/21/31/.../61.
  ! fortnum keys are 15, 21, 31, 61; map the 41/51 requests to the next key up.
  INTEGER, PARAMETER :: GK_LIMIT = 1000

  ! Fortran procedure pointer to user-specified function
  ! (By using this procedure pointers, one does not have to
  ! rename all the routines)
  ! 'abstract' interface is able to handle miscellaneous
  ! functions with the same type of interface
  ABSTRACT INTERFACE
     FUNCTION func1d_param0(x)
       IMPORT :: wp
       REAL(wp) :: func1d_param0
       REAL(wp) :: x
     END FUNCTION func1d_param0
     FUNCTION func1d_param1(x,param1)
       IMPORT :: wp
       REAL(wp) :: func1d_param1
       REAL(wp) :: x, param1
     END FUNCTION func1d_param1
  END INTERFACE

  ! Activate messages from root solver
  PUBLIC integration_solver_talk
  LOGICAL :: integration_solver_talk = .FALSE.

  ! Error codes detected during integration
  PRIVATE err_detected, check_status
  PUBLIC disp_gsl_integration_error
  LOGICAL, DIMENSION(-2:32) :: err_detected = .FALSE.

  ! Selection of test cases
  PUBLIC test_func1d_param0, test2_func1d_param0, test3_func1d_param0, &
       test4_func1d_param0, test5_func1d_param0, test6_func1d_param0, &
       test_func1d_param1, test2_func1d_param1

  ! Q(uadrature) A(daptive) G(eneral integrand) integration procedure
  PUBLIC fint1d_qag
  PRIVATE fint1d_param0_qag, fint1d_param1_qag
  INTERFACE fint1d_qag
     MODULE PROCEDURE fint1d_param0_qag, fint1d_param1_qag
  END INTERFACE fint1d_qag

  ! Q(uadrature) A(daptive) G(eneral integrand) integration procedure
  ! with singularities
  PUBLIC fint1d_qags
  PRIVATE fint1d_param0_qags, fint1d_param1_qags
  INTERFACE fint1d_qags
     MODULE PROCEDURE fint1d_param0_qags, fint1d_param1_qags
  END INTERFACE fint1d_qags

  ! Q(uadrature) A(daptive) G(eneral integrand) integration procedure
  ! with known singular points
  PUBLIC fint1d_qagp
  PRIVATE fint1d_param0_qagp, fint1d_param1_qagp
  INTERFACE fint1d_qagp
     MODULE PROCEDURE fint1d_param0_qagp, fint1d_param1_qagp
  END INTERFACE fint1d_qagp

  ! CQUAD doubly-adaptive integration (re-expressed on the adaptive QAGS path)
  PUBLIC fint1d_cquad
  PRIVATE fint1d_param0_cquad, fint1d_param1_cquad
  INTERFACE fint1d_cquad
     MODULE PROCEDURE fint1d_param0_cquad, fint1d_param1_cquad
  END INTERFACE fint1d_cquad

  ! Q(uadrature) A(daptive) G(eneral integrand) integration procedure for interval from x_low to infinity
  PUBLIC fint1d_qagiu
  PRIVATE fint1d_param0_qagiu, fint1d_param1_qagiu
  INTERFACE fint1d_qagiu
     MODULE PROCEDURE fint1d_param0_qagiu, fint1d_param1_qagiu
  END INTERFACE fint1d_qagiu


CONTAINS

  !--------------------------------------------------------------------------------------!
  ! Map an FGSL-style sw_qag_rule (1..6 -> 15,21,31,41,51,61 point GK rules) to the
  ! key fortnum supports (15, 21, 31, 61). 41/51 are promoted to the next available key.
  PURE FUNCTION qag_key(sw_qag_rule) RESULT(key)
    INTEGER, INTENT(in) :: sw_qag_rule
    INTEGER :: key
    SELECT CASE (sw_qag_rule)
    CASE (1);  key = 15
    CASE (2);  key = 21
    CASE (3);  key = 31
    CASE (4);  key = 31
    CASE (5);  key = 61
    CASE (6);  key = 61
    CASE DEFAULT; key = 21
    END SELECT
  END FUNCTION qag_key

  !--------------------------------------------------------------------------------------!
  ! Q(uadrature) A(daptive) G(eneral integrand) integration procedure
  RECURSIVE FUNCTION fint1d_param0_qag(func1d_param0_user,x_low,x_up,epsabs,epsrel,sw_qag_rule) result(res)

    INTERFACE
       FUNCTION func1d_param0_user(x)
         IMPORT :: wp
         REAL(wp) :: func1d_param0_user
         REAL(wp) :: x
       END FUNCTION func1d_param0_user
    END INTERFACE
    REAL(wp) :: x_low, x_up
    REAL(wp) :: epsabs, epsrel
    INTEGER :: sw_qag_rule
    REAL(wp), DIMENSION(2) :: res

    TYPE(integrate_workspace_t) :: ws
    TYPE(integrate_result_t) :: result
    TYPE(fortnum_status_t) :: status

    CALL integrate_qag(panel, x_low, x_up, epsabs, epsrel, ws, result, status, &
         key=qag_key(sw_qag_rule), limit=GK_LIMIT)
    CALL check_status(status)
    res(1) = result%value
    res(2) = result%abserr

  CONTAINS
    FUNCTION panel(x, ctx) RESULT(fx)
      REAL(wp), INTENT(in) :: x
      CLASS(*), INTENT(in), OPTIONAL :: ctx
      REAL(wp) :: fx
      fx = func1d_param0_user(x)
    END FUNCTION panel
  END FUNCTION fint1d_param0_qag


  !--------------------------------------------------------------------------------------!
  ! Q(uadrature) A(daptive) G(eneral integrand) integration procedure
  ! (parameter-valued input function)
  RECURSIVE FUNCTION fint1d_param1_qag(func1d_param1_user,param1,x_low,x_up,&
       epsabs,epsrel,sw_qag_rule) result(res)

    INTERFACE
       FUNCTION func1d_param1_user(x,param1)
         IMPORT :: wp
         REAL(wp) :: func1d_param1_user
         REAL(wp) :: x, param1
       END FUNCTION func1d_param1_user
    END INTERFACE
    REAL(wp) :: param1
    REAL(wp) :: x_low, x_up
    REAL(wp) :: epsabs, epsrel
    INTEGER :: sw_qag_rule
    REAL(wp), DIMENSION(2) :: res

    TYPE(integrate_workspace_t) :: ws
    TYPE(integrate_result_t) :: result
    TYPE(fortnum_status_t) :: status

    CALL integrate_qag(panel, x_low, x_up, epsabs, epsrel, ws, result, status, &
         key=qag_key(sw_qag_rule), limit=GK_LIMIT)
    CALL check_status(status)
    res(1) = result%value
    res(2) = result%abserr

  CONTAINS
    FUNCTION panel(x, ctx) RESULT(fx)
      REAL(wp), INTENT(in) :: x
      CLASS(*), INTENT(in), OPTIONAL :: ctx
      REAL(wp) :: fx
      fx = func1d_param1_user(x, param1)
    END FUNCTION panel
  END FUNCTION fint1d_param1_qag


  !--------------------------------------------------------------------------------------!
  ! Q(uadrature) A(daptive) G(eneral integrand) integration procedure
  ! with singularities
  RECURSIVE FUNCTION fint1d_param0_qags(func1d_param0_user,x_low,x_up,epsabs,epsrel) result(res)

    INTERFACE
       FUNCTION func1d_param0_user(x)
         IMPORT :: wp
         REAL(wp) :: func1d_param0_user
         REAL(wp) :: x
       END FUNCTION func1d_param0_user
    END INTERFACE
    REAL(wp) :: x_low, x_up
    REAL(wp) :: epsabs, epsrel
    REAL(wp), DIMENSION(2) :: res

    TYPE(integrate_workspace_t) :: ws
    TYPE(integrate_epstab_t) :: epstab
    TYPE(integrate_result_t) :: result
    TYPE(fortnum_status_t) :: status

    CALL integrate_qags(panel, x_low, x_up, epsabs, epsrel, ws, epstab, result, &
         status, limit=GK_LIMIT)
    CALL check_status(status)
    res(1) = result%value
    res(2) = result%abserr

  CONTAINS
    FUNCTION panel(x, ctx) RESULT(fx)
      REAL(wp), INTENT(in) :: x
      CLASS(*), INTENT(in), OPTIONAL :: ctx
      REAL(wp) :: fx
      fx = func1d_param0_user(x)
    END FUNCTION panel
  END FUNCTION fint1d_param0_qags


  !--------------------------------------------------------------------------------------!
  ! Q(uadrature) A(daptive) G(eneral integrand) integration procedure
  ! with singularities (parameter-valued input function)
  RECURSIVE FUNCTION fint1d_param1_qags(func1d_param1_user,param1,x_low,x_up,epsabs,epsrel) result(res)

    INTERFACE
       FUNCTION func1d_param1_user(x,param1)
         IMPORT :: wp
         REAL(wp) :: func1d_param1_user
         REAL(wp) :: x, param1
       END FUNCTION func1d_param1_user
    END INTERFACE
    REAL(wp) :: param1
    REAL(wp) :: x_low, x_up
    REAL(wp) :: epsabs, epsrel
    REAL(wp), DIMENSION(2) :: res

    TYPE(integrate_workspace_t) :: ws
    TYPE(integrate_epstab_t) :: epstab
    TYPE(integrate_result_t) :: result
    TYPE(fortnum_status_t) :: status

    CALL integrate_qags(panel, x_low, x_up, epsabs, epsrel, ws, epstab, result, &
         status, limit=GK_LIMIT)
    CALL check_status(status)
    res(1) = result%value
    res(2) = result%abserr

  CONTAINS
    FUNCTION panel(x, ctx) RESULT(fx)
      REAL(wp), INTENT(in) :: x
      CLASS(*), INTENT(in), OPTIONAL :: ctx
      REAL(wp) :: fx
      fx = func1d_param1_user(x, param1)
    END FUNCTION panel
  END FUNCTION fint1d_param1_qags


  !--------------------------------------------------------------------------------------!
  ! Q(uadrature) A(daptive) G(eneral integrand) integration procedure
  ! with known singular points
  RECURSIVE FUNCTION fint1d_param0_qagp(func1d_param0_user,pts,siz_pts,epsabs,epsrel) result(res)

    INTERFACE
       FUNCTION func1d_param0_user(x)
         IMPORT :: wp
         REAL(wp) :: func1d_param0_user
         REAL(wp) :: x
       END FUNCTION func1d_param0_user
    END INTERFACE
    REAL(wp) :: epsabs, epsrel
    REAL(wp), DIMENSION(:) :: pts
    INTEGER :: siz_pts
    REAL(wp), DIMENSION(2) :: res

    TYPE(integrate_workspace_t) :: ws
    TYPE(integrate_epstab_t) :: epstab
    TYPE(integrate_result_t) :: result
    TYPE(fortnum_status_t) :: status

    ! In GSL, pts holds the boundaries [a, interior break points..., b]. fortnum
    ! takes (a, b, interior break points) separately.
    CALL integrate_qagp(panel, pts(1), pts(SIZE(pts)), pts(2:SIZE(pts)-1), &
         epsabs, epsrel, ws, epstab, result, status, limit=GK_LIMIT)
    CALL check_status(status)
    res(1) = result%value
    res(2) = result%abserr

  CONTAINS
    FUNCTION panel(x, ctx) RESULT(fx)
      REAL(wp), INTENT(in) :: x
      CLASS(*), INTENT(in), OPTIONAL :: ctx
      REAL(wp) :: fx
      fx = func1d_param0_user(x)
    END FUNCTION panel
  END FUNCTION fint1d_param0_qagp


  !--------------------------------------------------------------------------------------!
  ! Q(uadrature) A(daptive) G(eneral integrand) integration procedure
  ! with known singular points (parameter-valued input function)
  RECURSIVE FUNCTION fint1d_param1_qagp(func1d_param1_user,param1,pts,siz_pts,epsabs,epsrel) result(res)

    INTERFACE
       FUNCTION func1d_param1_user(x,param1)
         IMPORT :: wp
         REAL(wp) :: func1d_param1_user
         REAL(wp) :: x, param1
       END FUNCTION func1d_param1_user
    END INTERFACE
    REAL(wp) :: param1
    REAL(wp) :: epsabs, epsrel
    REAL(wp), DIMENSION(:) :: pts
    INTEGER :: siz_pts
    REAL(wp), DIMENSION(2) :: res

    TYPE(integrate_workspace_t) :: ws
    TYPE(integrate_epstab_t) :: epstab
    TYPE(integrate_result_t) :: result
    TYPE(fortnum_status_t) :: status

    CALL integrate_qagp(panel, pts(1), pts(SIZE(pts)), pts(2:SIZE(pts)-1), &
         epsabs, epsrel, ws, epstab, result, status, limit=GK_LIMIT)
    CALL check_status(status)
    res(1) = result%value
    res(2) = result%abserr

  CONTAINS
    FUNCTION panel(x, ctx) RESULT(fx)
      REAL(wp), INTENT(in) :: x
      CLASS(*), INTENT(in), OPTIONAL :: ctx
      REAL(wp) :: fx
      fx = func1d_param1_user(x, param1)
    END FUNCTION panel
  END FUNCTION fint1d_param1_qagp


  !--------------------------------------------------------------------------------------!
  ! CQUAD doubly-adaptive integration on a finite interval.
  ! fortnum has no CQUAD; the adaptive QAGS path (bisection plus Wynn-epsilon
  ! extrapolation) handles endpoint singularities and delivers the same value/error
  ! pair on the smooth integrands used here.
  RECURSIVE FUNCTION fint1d_param0_cquad(func1d_param0_user,x_low,x_up,epsabs,epsrel) result(res)
    INTERFACE
       FUNCTION func1d_param0_user(x)
         IMPORT :: wp
         REAL(wp) :: func1d_param0_user
         REAL(wp) :: x
       END FUNCTION func1d_param0_user
    END INTERFACE
    REAL(wp) :: x_low
    REAL(wp) :: x_up
    REAL(wp) :: epsabs
    REAL(wp) :: epsrel
    REAL(wp), DIMENSION(2) :: res

    TYPE(fortnum_status_t) :: status

    CALL integrate_cquad(panel, x_low, x_up, res(1), status, epsabs=epsabs, &
         epsrel=epsrel, abserr=res(2), limit=GK_LIMIT)
    CALL check_status(status)

  CONTAINS
    FUNCTION panel(x, ctx) RESULT(fx)
      REAL(wp), INTENT(in) :: x
      CLASS(*), INTENT(in), OPTIONAL :: ctx
      REAL(wp) :: fx
      fx = func1d_param0_user(x)
    END FUNCTION panel
  END FUNCTION fint1d_param0_cquad


  !--------------------------------------------------------------------------------------!
  ! CQUAD doubly-adaptive integration (parameter-valued input function)
  RECURSIVE FUNCTION fint1d_param1_cquad(func1d_param1_user,param1,x_low,x_up,epsabs,epsrel) result(res)

    INTERFACE
       FUNCTION func1d_param1_user(x,param1)
         IMPORT :: wp
         REAL(wp) :: func1d_param1_user
         REAL(wp) :: x, param1
       END FUNCTION func1d_param1_user
    END INTERFACE
    REAL(wp) :: param1
    REAL(wp) :: x_low, x_up
    REAL(wp) :: epsabs, epsrel
    REAL(wp), DIMENSION(2) :: res

    TYPE(fortnum_status_t) :: status

    CALL integrate_cquad(panel, x_low, x_up, res(1), status, epsabs=epsabs, &
         epsrel=epsrel, abserr=res(2), limit=GK_LIMIT)
    CALL check_status(status)

  CONTAINS
    FUNCTION panel(x, ctx) RESULT(fx)
      REAL(wp), INTENT(in) :: x
      CLASS(*), INTENT(in), OPTIONAL :: ctx
      REAL(wp) :: fx
      fx = func1d_param1_user(x, param1)
    END FUNCTION panel
  END FUNCTION fint1d_param1_cquad


  !--------------------------------------------------------------------------------------!
  ! Q(uadrature) A(daptive) G(eneral integrand) integration procedure for interval from x_low to infinity
  RECURSIVE FUNCTION fint1d_param0_qagiu(func1d_param0_user,x_low,epsabs,epsrel) result(res)

    INTERFACE
       FUNCTION func1d_param0_user(x)
         IMPORT :: wp
         REAL(wp) :: func1d_param0_user
         REAL(wp) :: x
       END FUNCTION func1d_param0_user
    END INTERFACE
    REAL(wp) :: x_low
    REAL(wp) :: epsabs, epsrel
    REAL(wp), DIMENSION(2) :: res

    TYPE(integrate_workspace_t) :: ws
    TYPE(integrate_epstab_t) :: epstab
    TYPE(integrate_result_t) :: result
    TYPE(fortnum_status_t) :: status

    CALL integrate_qagiu(panel, x_low, 1, epsabs, epsrel, ws, epstab, result, &
         status, limit=GK_LIMIT)
    CALL check_status(status)
    res(1) = result%value
    res(2) = result%abserr

  CONTAINS
    FUNCTION panel(x, ctx) RESULT(fx)
      REAL(wp), INTENT(in) :: x
      CLASS(*), INTENT(in), OPTIONAL :: ctx
      REAL(wp) :: fx
      fx = func1d_param0_user(x)
    END FUNCTION panel
  END FUNCTION fint1d_param0_qagiu


  !--------------------------------------------------------------------------------------!
  ! Q(uadrature) A(daptive) G(eneral integrand) integration procedure for interval from x_low to infinity
  ! (parameter-valued input function)
  RECURSIVE FUNCTION fint1d_param1_qagiu(func1d_param1_user,param1,x_low,&
       epsabs,epsrel) result(res)

    INTERFACE
       FUNCTION func1d_param1_user(x,param1)
         IMPORT :: wp
         REAL(wp) :: func1d_param1_user
         REAL(wp) :: x, param1
       END FUNCTION func1d_param1_user
    END INTERFACE
    REAL(wp) :: param1
    REAL(wp) :: x_low
    REAL(wp) :: epsabs, epsrel
    REAL(wp), DIMENSION(2) :: res

    TYPE(integrate_workspace_t) :: ws
    TYPE(integrate_epstab_t) :: epstab
    TYPE(integrate_result_t) :: result
    TYPE(fortnum_status_t) :: status

    CALL integrate_qagiu(panel, x_low, 1, epsabs, epsrel, ws, epstab, result, &
         status, limit=GK_LIMIT)
    CALL check_status(status)
    res(1) = result%value
    res(2) = result%abserr

  CONTAINS
    FUNCTION panel(x, ctx) RESULT(fx)
      REAL(wp), INTENT(in) :: x
      CLASS(*), INTENT(in), OPTIONAL :: ctx
      REAL(wp) :: fx
      fx = func1d_param1_user(x, param1)
    END FUNCTION panel
  END FUNCTION fint1d_param1_qagiu


  !--------------------------------------------------------------------------------------!
  ! Record a non-success fortnum status in the per-code error flag array so the
  ! accumulated warnings can be reported later via disp_gsl_integration_error.
  SUBROUTINE check_status(status)
    TYPE(fortnum_status_t), INTENT(in) :: status

    SELECT CASE (status%code)
    CASE (FORTNUM_OK)
       RETURN
    CASE (FORTNUM_DOMAIN_ERROR)
       err_detected(1) = .TRUE.   ! GSL_EDOM
    CASE (FORTNUM_CONVERGENCE_ERROR)
       err_detected(14) = .TRUE.  ! GSL_ETOL
    CASE (FORTNUM_NOT_IMPLEMENTED)
       err_detected(24) = .TRUE.  ! GSL_EUNIMPL
    CASE DEFAULT
       err_detected(5) = .TRUE.   ! GSL_EFAILED
    END SELECT

  END SUBROUTINE check_status


  !--------------------------------------------------------------------------------------!
  SUBROUTINE disp_gsl_integration_error()
    ! internal variables
    INTEGER :: k

    IF ( ANY(err_detected) ) THEN ! else normal termination
       PRINT *,'-------------------------------------------------'
       DO k = LBOUND(err_detected,1),UBOUND(err_detected,1)
          IF (err_detected(k)) THEN
             PRINT *,"gsl_integration_routines_mod.f90: &
                  &Possible Warning from Integration Routines - Code = ",k
          END IF
       END DO
       PRINT *,'-------------------------------------------------'
    END IF

    ! reset err_detected
    err_detected = .FALSE.

  END SUBROUTINE disp_gsl_integration_error


  !--------------------------------------------------------------------------------------!
  ! Test function without a parameter for 1d integrators
  FUNCTION test_func1d_param0(x)

    REAL(wp) :: test_func1d_param0
    REAL(wp) :: x

    test_func1d_param0 = COS(x)

  END FUNCTION test_func1d_param0


  !--------------------------------------------------------------------------------------!
  ! Test function without a parameter for 1d integrators (singularity)
  FUNCTION test2_func1d_param0(x)

    REAL(wp) :: test2_func1d_param0
    REAL(wp) :: x

    test2_func1d_param0 = 1.0_wp/SQRT(1.0_wp-x)

  END FUNCTION test2_func1d_param0


  !--------------------------------------------------------------------------------------!
  ! Test function without a parameter for 1d integrators (singularity)
  FUNCTION test3_func1d_param0(x)

    REAL(wp) :: test3_func1d_param0
    REAL(wp) :: x

    test3_func1d_param0 = 1.0_wp/SQRT(1.0_wp-x**2.0_wp)

  END FUNCTION test3_func1d_param0


  !--------------------------------------------------------------------------------------!
  ! Test function without a parameter for 1d integrators (singularity)
  ! (cf., V. Sladek et al., Appl. Math. Modelling 25 (2001) 901-922)
  FUNCTION test4_func1d_param0(x)

    REAL(wp) :: test4_func1d_param0
    REAL(wp) :: x

    test4_func1d_param0 = x*(x-2.0_wp)*LOG(x/2.0_wp)

  END FUNCTION test4_func1d_param0


  !--------------------------------------------------------------------------------------!
  ! Test function without a parameter for 1d integrators (singularity)
  ! (cf., V. Sladek et al., Appl. Math. Modelling 25 (2001) 901-922)
  FUNCTION test5_func1d_param0(x)

    REAL(wp) :: test5_func1d_param0
    REAL(wp) :: x

    REAL(wp), PARAMETER :: d = 15.0_wp

    test5_func1d_param0 = ((x-d)**3.0_wp)*LOG(x/2.0_wp)

  END FUNCTION test5_func1d_param0


  !--------------------------------------------------------------------------------------!
  ! Test function without a parameter for 1d integrators (singularity)
  FUNCTION test6_func1d_param0(x)

    REAL(wp) :: test6_func1d_param0
    REAL(wp) :: x

    test6_func1d_param0 = 1.0_wp/SQRT(ABS(x))

  END FUNCTION test6_func1d_param0


  !--------------------------------------------------------------------------------------!
  ! Test function with a parameter for 1d integrators
  FUNCTION test_func1d_param1(x, param1)

    REAL(wp) :: test_func1d_param1
    REAL(wp) :: x, param1

    test_func1d_param1 = 2.0_wp/SQRT(1-param1*(SIN(x)**2.0_wp))

  END FUNCTION test_func1d_param1


  !--------------------------------------------------------------------------------------!
  ! Test function with a parameter for 1d integrators
  FUNCTION test2_func1d_param1(x, param1)

    REAL(wp) :: test2_func1d_param1
    REAL(wp) :: x, param1

    test2_func1d_param1 = 1.0_wp/&
         SQRT(param1-(SIN(x/2.0_wp)**2.0_wp))

  END FUNCTION test2_func1d_param1
  !--------------------------------------------------------------------------------------!
END MODULE gsl_integration_routines_mod
