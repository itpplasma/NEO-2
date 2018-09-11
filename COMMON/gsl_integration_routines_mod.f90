!-----------------------------------------------------------------------------------------!
! module: gsl_integration_routines_mod                                                    !
! authors: Andreas F. Martitsch, Gernot Kapper                                            !
! date: 27.07.2015                                                                        !
! version: 0.3                                                                            !
! description:                                                                            !
! Module gsl_integration_routines_mod provides routines for integrating 1d functions      !
! The numerical routines used here rely on the GNU Scientific Library (GSL).              !
! The C-interface between the GSL and the Fortran code is provided by FGSL                !
! (See http://www.lrz.de/services/software/mathematik/gsl/fortran/                        !
! for further information)                                                                !
! changes: Global procedure pointers and wrapper-functions are replaced by internal       !
! subroutines attached to the integration routines. This allows for a nested call of      !
! integration routines from outside (e.g, 2d integration)                                 !
!-----------------------------------------------------------------------------------------!

!-----------------------------------------------------------------------------------------!
! History                                                                                 !
! 0.1 - Initial version                                                                   !
! 0.2 - Support for 2D integration                                                        !
! 0.3 - Added recursive keyword to functions for 2D integration                           !
!-----------------------------------------------------------------------------------------!


MODULE gsl_integration_routines_mod
  !
  USE fgsl ! Fortran interface of the GSL Library
  USE, INTRINSIC :: iso_c_binding
  !
  IMPLICIT NONE
  !
  ! Possible values of absolute/relative errors (FGSL)
  ! (Note: fgsl_double (etc.) originate from c-binding.
  ! This is only to make sure that C and Fortran routines
  ! work with the same definition of float numbers (etc.).
  ! For practical applications, one can simply replace 
  ! "fgsl_double" by "double precision" (etc.). 
  PRIVATE eps5, eps7, eps10, eps12
  REAL(fgsl_double), PARAMETER :: eps5 = 1.0E-5_fgsl_double
  REAL(fgsl_double), PARAMETER :: eps7 = 1.0E-7_fgsl_double
  REAL(fgsl_double), PARAMETER :: eps10 = 1.0E-10_fgsl_double
  REAL(fgsl_double), PARAMETER :: eps12 = 1.0E-12_fgsl_double
  !
  ! Fortran procedure pointer to user-specified function
  ! (By using this procedure pointers, one does not have to
  ! rename all the routines)
  ! 'abstract' interface is able to handle miscellaneous
  ! functions with the same type of interface
  ABSTRACT INTERFACE  
     FUNCTION func1d_param0(x)
       USE fgsl
       REAL(fgsl_double) :: func1d_param0
       REAL(fgsl_double) :: x
     END FUNCTION func1d_param0
     FUNCTION func1d_param1(x,param1)
       USE fgsl
       REAL(fgsl_double) :: func1d_param1
       REAL(fgsl_double) :: x, param1
     END FUNCTION func1d_param1
  END INTERFACE
  ! This part is replaced by internal subroutines (01.04.2015)
!!$  ! declare global procedure pointers (maybe there is a better solution)
!!$  PRIVATE func1d_param0_ptr
!!$  PROCEDURE(func1d_param0), POINTER :: func1d_param0_ptr => NULL()
!!$  PRIVATE func1d_param1_ptr
!!$  PROCEDURE(func1d_param1), POINTER :: func1d_param1_ptr => NULL()
  !
  ! Activate messages from root solver
  PUBLIC integration_solver_talk
  LOGICAL :: integration_solver_talk = .FALSE.
  !
  ! GSL error codes
  PRIVATE gsl_err_detected, check_error
  PUBLIC disp_gsl_integration_error
  LOGICAL, DIMENSION(-2:32) :: gsl_err_detected = .FALSE.
  !
  ! This part is replaced by internal subroutines (01.04.2015)
!!$  ! Internal Fortran-to-C wrapper functions (public, because of c-binding)
!!$  PUBLIC f2c_wrapper_func1d_param0, f2c_wrapper_func1d_param1
  !
  ! Selection of test cases
  PUBLIC test_func1d_param0, test2_func1d_param0, test3_func1d_param0, &
       test4_func1d_param0, test5_func1d_param0, test6_func1d_param0, &
       test_func1d_param1, test2_func1d_param1
  !
  ! Q(uadrature) A(daptive) G(eneral integrand) integration procedure
  PUBLIC fint1d_qag
  PRIVATE fint1d_param0_qag, fint1d_param1_qag
  INTERFACE fint1d_qag
     MODULE PROCEDURE fint1d_param0_qag, fint1d_param1_qag
  END INTERFACE fint1d_qag
  !
  ! Q(uadrature) A(daptive) G(eneral integrand) integration procedure
  ! with singularities
  PUBLIC fint1d_qags
  PRIVATE fint1d_param0_qags, fint1d_param1_qags
  INTERFACE fint1d_qags
     MODULE PROCEDURE fint1d_param0_qags, fint1d_param1_qags
  END INTERFACE fint1d_qags
  !
  ! Q(uadrature) A(daptive) G(eneral integrand) integration procedure
  ! with known singular points
  PUBLIC fint1d_qagp
  PRIVATE fint1d_param0_qagp, fint1d_param1_qagp
  INTERFACE fint1d_qagp
     MODULE PROCEDURE fint1d_param0_qagp, fint1d_param1_qagp
  END INTERFACE fint1d_qagp
  !
  ! CQUAD doubly-adaptive integration
  PUBLIC fint1d_cquad
  PRIVATE fint1d_param0_cquad, fint1d_param1_cquad
  INTERFACE fint1d_cquad
     MODULE PROCEDURE fint1d_param0_cquad, fint1d_param1_cquad
  END INTERFACE fint1d_cquad
  !
  ! Q(uadrature) A(daptive) G(eneral integrand) integration procedure for interval from x_low to infinity
  PUBLIC fint1d_qagiu
  PRIVATE fint1d_param0_qagiu, fint1d_param1_qagiu
  INTERFACE fint1d_qagiu
     MODULE PROCEDURE fint1d_param0_qagiu, fint1d_param1_qagiu
  END INTERFACE fint1d_qagiu
  !
  
CONTAINS
  !--------------------------------------------------------------------------------------!
  ! This part is replaced by internal subroutines (01.04.2015)
!!$  ! Fortran-to-C wrapper function in order to guarantee interoperability between
!!$  ! Fortran and C
!!$  FUNCTION f2c_wrapper_func1d_param0(x, params) BIND(c)
!!$    !
!!$    REAL(c_double), VALUE :: x
!!$    TYPE(c_ptr), VALUE :: params
!!$    REAL(c_double) :: f2c_wrapper_func1d_param0
!!$    !
!!$    ! Check, whether C-pointer 'params' is not associated (no parameter passed)
!!$    IF(C_ASSOCIATED(params)) STOP '***Error*** in f2c_wrapper_func1d_param0'
!!$    !
!!$    ! Wrap user-specified function to a C-interoperable function
!!$    f2c_wrapper_func1d_param0 = func1d_param0_ptr(x)
!!$    !
!!$  END FUNCTION f2c_wrapper_func1d_param0
  !--------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------!
  ! This part is replaced by internal subroutines (01.04.2015)
!!$  ! Fortran-to-C wrapper function in order to guarantee interoperability between
!!$  ! Fortran and C
!!$  FUNCTION f2c_wrapper_func1d_param1(x, params) BIND(c)
!!$    !
!!$    REAL(c_double), VALUE :: x
!!$    TYPE(c_ptr), VALUE :: params
!!$    REAL(c_double) :: f2c_wrapper_func1d_param1
!!$    !
!!$    REAL(fgsl_double), POINTER :: p ! This is the type of C-pointer to be casted
!!$    !
!!$    ! Check, whether C-pointer 'params' is associated
!!$    IF(.NOT. C_ASSOCIATED(params)) STOP '***Error*** in f2c_wrapper_func1d_param1'
!!$    !
!!$    ! Cast C-pointer to the above-defined Fortran pointer
!!$    CALL C_F_POINTER(params, p)
!!$    ! Wrap user-specified function to a C-interoperable function
!!$    f2c_wrapper_func1d_param1 = func1d_param1_ptr(x,p)
!!$    !
!!$  END FUNCTION f2c_wrapper_func1d_param1
  !--------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------!
  ! Q(uadrature) A(daptive) G(eneral integrand) integration procedure
  RECURSIVE FUNCTION fint1d_param0_qag(func1d_param0_user,x_low,x_up,epsabs,epsrel,sw_qag_rule) result(res)
    !
    INTERFACE  
       FUNCTION func1d_param0_user(x)
         USE fgsl
         REAL(fgsl_double) :: func1d_param0_user
         REAL(fgsl_double) :: x
       END FUNCTION func1d_param0_user
    END INTERFACE
    ! x_low: lower boundary, x_up: upper boundary
    REAL(fgsl_double) :: x_low, x_up
    ! epsabs: absolute error, epsrel: relative error
    REAL(fgsl_double) :: epsabs, epsrel
    ! sw_qag_rule: 1=15, 2=21, 3=31, 4=41, 5=51 and 6=61 point Gauss-Kronrod rules
    INTEGER(fgsl_int) :: sw_qag_rule 
    REAL(fgsl_double), DIMENSION(2) :: res
    !
    INTEGER(fgsl_size_t), PARAMETER :: limit = 1000_fgsl_size_t
    REAL(fgsl_double) :: ra, rda ! result and absolute error
    TYPE(fgsl_error_handler_t) :: std
    INTEGER(fgsl_int) :: status
    TYPE(fgsl_function) :: stdfunc
    TYPE(fgsl_integration_workspace) :: integ_wk
    TYPE(c_ptr) :: param0_ptr ! This pointer holds the C-location of user-specified
    !                         ! parameter (here is no parameter specified -> c_null_ptr)
    !
    ! Turn off error handler
    std = fgsl_set_error_handler_off()
    !
    ! This part is replaced by internal subroutines (01.04.2015)
!!$    ! Associate global procedure pointer with user-specified function
!!$    func1d_param0_ptr => func1d_param0_user
!!$    IF(.NOT. ASSOCIATED(func1d_param0_ptr)) &
!!$         STOP '***Error*** in fint1d_param0_qag'
    !
    ! Nullify param0_ptr
    param0_ptr = c_null_ptr
    IF(C_ASSOCIATED(param0_ptr)) &
         STOP '***Error*** in fint1d_param0_qag'
    !
    ! Initialize fgsl_function and workspace
    stdfunc = fgsl_function_init(f2c_wrapper_func1d_param0, param0_ptr)
    integ_wk = fgsl_integration_workspace_alloc(limit)
    !
    ! Initialize solver 'fgsl_integration_qag' to use the function 'stdfunc' and
    ! the user-specified parameters
    status = fgsl_integration_qag(stdfunc, x_low, x_up, &
       epsabs, epsrel, limit, sw_qag_rule, integ_wk, ra, rda)
    CALL check_error(status)
    !
    ! Return the results (ra,rda)
    res(1) = ra
    res(2) = rda
    !
    ! Free memory
    CALL fgsl_function_free(stdfunc)
    CALL fgsl_integration_workspace_free(integ_wk)
    !
  CONTAINS
    ! Fortran-to-C wrapper function in order to guarantee interoperability between
    ! Fortran and C
    RECURSIVE FUNCTION f2c_wrapper_func1d_param0(x, params) BIND(c)
      !
      REAL(c_double), VALUE :: x
      TYPE(c_ptr), VALUE :: params
      REAL(c_double) :: f2c_wrapper_func1d_param0
      !
      ! Check, whether C-pointer 'params' is not associated (no parameter passed)
      IF(C_ASSOCIATED(params)) STOP '***Error*** in f2c_wrapper_func1d_param0'
      !
      ! Wrap user-specified function to a C-interoperable function
      f2c_wrapper_func1d_param0 = func1d_param0_user(x)
      !
    END FUNCTION f2c_wrapper_func1d_param0
  END FUNCTION fint1d_param0_qag
  !--------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------!
  ! Q(uadrature) A(daptive) G(eneral integrand) integration procedure
  ! (parameter-valued input function)
  RECURSIVE FUNCTION fint1d_param1_qag(func1d_param1_user,param1,x_low,x_up,&
       epsabs,epsrel,sw_qag_rule) result(res)
    !
    INTERFACE  
       FUNCTION func1d_param1_user(x,param1)
         USE fgsl
         REAL(fgsl_double) :: func1d_param1_user
         REAL(fgsl_double) :: x, param1
       END FUNCTION func1d_param1_user
    END INTERFACE
    REAL(fgsl_double), TARGET :: param1
    ! x_low: lower boundary, x_up: upper boundary
    REAL(fgsl_double) :: x_low, x_up
    ! epsabs: absolute error, epsrel: relative error
    REAL(fgsl_double) :: epsabs, epsrel
    ! sw_qag_rule: 1=15, 2=21, 3=31, 4=41, 5=51 and 6=61 point Gauss-Kronrod rules
    INTEGER(fgsl_int) :: sw_qag_rule 
    REAL(fgsl_double), DIMENSION(2) :: res
    !
    INTEGER(fgsl_size_t), PARAMETER :: limit = 1000_fgsl_size_t
    REAL(fgsl_double) :: ra, rda ! result and absolute error
    TYPE(fgsl_error_handler_t) :: std
    INTEGER(fgsl_int) :: status
    TYPE(fgsl_function) :: stdfunc
    TYPE(fgsl_integration_workspace) :: integ_wk
    TYPE(c_ptr) :: param1_ptr ! This pointer holds the C-location of user-specified
    !                         ! parameter
    !
    ! Turn off error handler
    std = fgsl_set_error_handler_off()
    !
    ! This part is replaced by internal subroutines (01.04.2015)
!!$    ! Associate global procedure pointer with user-specified function
!!$    func1d_param1_ptr => func1d_param1_user
!!$    IF(.NOT. ASSOCIATED(func1d_param1_ptr)) &
!!$         STOP '***Error*** in fint1d_param1_qag'
    !
    ! Get C-location of 'param1'
    param1_ptr = c_LOC(param1)
    IF(.NOT. C_ASSOCIATED(param1_ptr)) &
         STOP '***Error*** in fint1d_param1_qag'
    !
    ! Initialize fgsl_function and workspace
    stdfunc = fgsl_function_init(f2c_wrapper_func1d_param1, param1_ptr)
    integ_wk = fgsl_integration_workspace_alloc(limit)
    !
    ! Initialize solver 'fgsl_integration_qag' to use the function 'stdfunc' and
    ! the user-specified parameters
    status = fgsl_integration_qag(stdfunc, x_low, x_up, &
       epsabs, epsrel, limit, sw_qag_rule, integ_wk, ra, rda)
    CALL check_error(status)
    !
    ! Return the results (ra,rda)
    res(1) = ra
    res(2) = rda
    !
    ! Free memory
    CALL fgsl_function_free(stdfunc)
    CALL fgsl_integration_workspace_free(integ_wk)
    !
  CONTAINS
    ! Fortran-to-C wrapper function in order to guarantee interoperability between
    ! Fortran and C
    RECURSIVE FUNCTION f2c_wrapper_func1d_param1(x, params) BIND(c)
      !
      REAL(c_double), VALUE :: x
      TYPE(c_ptr), VALUE :: params
      REAL(c_double) :: f2c_wrapper_func1d_param1
      !
      REAL(fgsl_double), POINTER :: p ! This is the type of C-pointer to be casted
      !
      ! Check, whether C-pointer 'params' is associated
      IF(.NOT. C_ASSOCIATED(params)) STOP '***Error*** in f2c_wrapper_func1d_param1'
      !
      ! Cast C-pointer to the above-defined Fortran pointer
      CALL C_F_POINTER(params, p)
      ! Wrap user-specified function to a C-interoperable function
      f2c_wrapper_func1d_param1 = func1d_param1_user(x,p)
      !
    END FUNCTION f2c_wrapper_func1d_param1
  END FUNCTION fint1d_param1_qag
  !--------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------!
  ! Q(uadrature) A(daptive) G(eneral integrand) integration procedure
  ! with singularities
  RECURSIVE FUNCTION fint1d_param0_qags(func1d_param0_user,x_low,x_up,epsabs,epsrel) result(res)
    !
    INTERFACE  
       FUNCTION func1d_param0_user(x)
         USE fgsl
         REAL(fgsl_double) :: func1d_param0_user
         REAL(fgsl_double) :: x
       END FUNCTION func1d_param0_user
    END INTERFACE
    ! x_low: lower boundary, x_up: upper boundary
    REAL(fgsl_double) :: x_low, x_up
    ! epsabs: absolute error, epsrel: relative error
    REAL(fgsl_double) :: epsabs, epsrel
    REAL(fgsl_double), DIMENSION(2) :: res
    !
    INTEGER(fgsl_size_t), PARAMETER :: limit = 1000_fgsl_size_t
    REAL(fgsl_double) :: ra, rda ! result and absolute error
    TYPE(fgsl_error_handler_t) :: std
    INTEGER(fgsl_int) :: status
    TYPE(fgsl_function) :: stdfunc
    TYPE(fgsl_integration_workspace) :: integ_wk
    TYPE(c_ptr) :: param0_ptr ! This pointer holds the C-location of user-specified
    !                         ! parameter (here is no parameter specified -> c_null_ptr)
    !
    ! Turn off error handler
    std = fgsl_set_error_handler_off()
    !
    ! This part is replaced by internal subroutines (01.04.2015)
!!$    ! Associate global procedure pointer with user-specified function
!!$    func1d_param0_ptr => func1d_param0_user
!!$    IF(.NOT. ASSOCIATED(func1d_param0_ptr)) &
!!$         STOP '***Error*** in fint1d_param0_qags'
    !
    ! Nullify param0_ptr
    param0_ptr = c_null_ptr
    IF(C_ASSOCIATED(param0_ptr)) &
         STOP '***Error*** in fint1d_param0_qags'
    !
    ! Initialize fgsl_function and workspace
    stdfunc = fgsl_function_init(f2c_wrapper_func1d_param0, param0_ptr)
    integ_wk = fgsl_integration_workspace_alloc(limit)
    !
    ! Initialize solver 'fgsl_integration_qags' to use the function 'stdfunc' and
    ! the user-specified parameters
    status = fgsl_integration_qags(stdfunc, x_low, x_up, &
       epsabs, epsrel, limit, integ_wk, ra, rda)
    CALL check_error(status)
    !
    ! Return the results (ra,rda)
    res(1) = ra
    res(2) = rda
    !
    ! Free memory
    CALL fgsl_function_free(stdfunc)
    CALL fgsl_integration_workspace_free(integ_wk)
    !
  CONTAINS
    ! Fortran-to-C wrapper function in order to guarantee interoperability between
    ! Fortran and C
    RECURSIVE FUNCTION f2c_wrapper_func1d_param0(x, params) BIND(c)
      !
      REAL(c_double), VALUE :: x
      TYPE(c_ptr), VALUE :: params
      REAL(c_double) :: f2c_wrapper_func1d_param0
      !
      ! Check, whether C-pointer 'params' is not associated (no parameter passed)
      IF(C_ASSOCIATED(params)) STOP '***Error*** in f2c_wrapper_func1d_param0'
      !
      ! Wrap user-specified function to a C-interoperable function
      f2c_wrapper_func1d_param0 = func1d_param0_user(x)
      !
    END FUNCTION f2c_wrapper_func1d_param0
  END FUNCTION fint1d_param0_qags
  !--------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------!
  ! Q(uadrature) A(daptive) G(eneral integrand) integration procedure
  ! with singularities (parameter-valued input function)
  RECURSIVE FUNCTION fint1d_param1_qags(func1d_param1_user,param1,x_low,x_up,epsabs,epsrel) result(res)
    !
    INTERFACE  
       FUNCTION func1d_param1_user(x,param1)
         USE fgsl
         REAL(fgsl_double) :: func1d_param1_user
         REAL(fgsl_double) :: x, param1
       END FUNCTION func1d_param1_user
    END INTERFACE
    REAL(fgsl_double), TARGET :: param1
    ! x_low: lower boundary, x_up: upper boundary
    REAL(fgsl_double) :: x_low, x_up
    ! epsabs: absolute error, epsrel: relative error
    REAL(fgsl_double) :: epsabs, epsrel
    REAL(fgsl_double), DIMENSION(2) :: res
    !
    INTEGER(fgsl_size_t), PARAMETER :: limit = 1000_fgsl_size_t
    REAL(fgsl_double) :: ra, rda ! result and absolute error
    TYPE(fgsl_error_handler_t) :: std
    INTEGER(fgsl_int) :: status
    TYPE(fgsl_function) :: stdfunc
    TYPE(fgsl_integration_workspace) :: integ_wk
    TYPE(c_ptr) :: param1_ptr ! This pointer holds the C-location of user-specified
    !                         ! parameter
    !
    ! Turn off error handler
    std = fgsl_set_error_handler_off()
    !
    ! This part is replaced by internal subroutines (01.04.2015)
!!$    ! Associate global procedure pointer with user-specified function
!!$    func1d_param1_ptr => func1d_param1_user
!!$    IF(.NOT. ASSOCIATED(func1d_param1_ptr)) &
!!$         STOP '***Error*** in fint1d_param1_qags'
    !
    ! Get C-location of 'param1'
    param1_ptr = c_LOC(param1)
    IF(.NOT. C_ASSOCIATED(param1_ptr)) &
         STOP '***Error*** in fint1d_param1_qags'
    !
    ! Initialize fgsl_function and workspace
    stdfunc = fgsl_function_init(f2c_wrapper_func1d_param1, param1_ptr)
    integ_wk = fgsl_integration_workspace_alloc(limit)
    !
    ! Initialize solver 'fgsl_integration_qags' to use the function 'stdfunc' and
    ! the user-specified parameters
    status = fgsl_integration_qags(stdfunc, x_low, x_up, &
       epsabs, epsrel, limit, integ_wk, ra, rda)
    CALL check_error(status)
    !
    ! Return the results (ra,rda)
    res(1) = ra
    res(2) = rda
    !
    ! Free memory
    CALL fgsl_function_free(stdfunc)
    CALL fgsl_integration_workspace_free(integ_wk)
    !
    CONTAINS
    ! Fortran-to-C wrapper function in order to guarantee interoperability between
    ! Fortran and C
    RECURSIVE FUNCTION f2c_wrapper_func1d_param1(x, params) BIND(c)
      !
      REAL(c_double), VALUE :: x
      TYPE(c_ptr), VALUE :: params
      REAL(c_double) :: f2c_wrapper_func1d_param1
      !
      REAL(fgsl_double), POINTER :: p ! This is the type of C-pointer to be casted
      !
      ! Check, whether C-pointer 'params' is associated
      IF(.NOT. C_ASSOCIATED(params)) STOP '***Error*** in f2c_wrapper_func1d_param1'
      !
      ! Cast C-pointer to the above-defined Fortran pointer
      CALL C_F_POINTER(params, p)
      ! Wrap user-specified function to a C-interoperable function
      f2c_wrapper_func1d_param1 = func1d_param1_user(x,p)
      !
    END FUNCTION f2c_wrapper_func1d_param1
  END FUNCTION fint1d_param1_qags
  !--------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------!
  ! Q(uadrature) A(daptive) G(eneral integrand) integration procedure
  ! with known singular points
  RECURSIVE FUNCTION fint1d_param0_qagp(func1d_param0_user,pts,siz_pts,epsabs,epsrel) result(res)
    !
    INTERFACE  
       FUNCTION func1d_param0_user(x)
         USE fgsl
         REAL(fgsl_double) :: func1d_param0_user
         REAL(fgsl_double) :: x
       END FUNCTION func1d_param0_user
    END INTERFACE
    ! epsabs: absolute error, epsrel: relative error
    REAL(fgsl_double) :: epsabs, epsrel
    ! specify singular points
    REAL(fgsl_double), DIMENSION(:) :: pts
    INTEGER(fgsl_size_t) :: siz_pts
    REAL(fgsl_double), DIMENSION(2) :: res
    !
    INTEGER(fgsl_size_t), PARAMETER :: limit = 1000_fgsl_size_t
    REAL(fgsl_double) :: ra, rda ! result and absolute error
    TYPE(fgsl_error_handler_t) :: std
    INTEGER(fgsl_int) :: status
    TYPE(fgsl_function) :: stdfunc
    TYPE(fgsl_integration_workspace) :: integ_wk
    TYPE(c_ptr) :: param0_ptr ! This pointer holds the C-location of user-specified
    !                         ! parameter (here is no parameter specified -> c_null_ptr)
    !
    ! Turn off error handler
    std = fgsl_set_error_handler_off()
    !
    ! This part is replaced by internal subroutines (01.04.2015)
!!$    ! Associate global procedure pointer with user-specified function
!!$    func1d_param0_ptr => func1d_param0_user
!!$    IF(.NOT. ASSOCIATED(func1d_param0_ptr)) &
!!$         STOP '***Error*** in fint1d_param0_qagp'
    !
    ! Nullify param0_ptr
    param0_ptr = c_null_ptr
    IF(C_ASSOCIATED(param0_ptr)) &
         STOP '***Error*** in fint1d_param0_qagp'
    !
    ! Initialize fgsl_function and workspace
    stdfunc = fgsl_function_init(f2c_wrapper_func1d_param0, param0_ptr)
    integ_wk = fgsl_integration_workspace_alloc(limit)
    !
    ! Initialize solver 'fgsl_integration_qagp' to use the function 'stdfunc' and
    ! the user-specified parameters
    status = fgsl_integration_qagp(stdfunc, pts, & !siz_pts, &
         epsabs, epsrel, limit, integ_wk, ra, rda)
    CALL check_error(status)
    !
    ! Return the results (ra,rda)
    res(1) = ra
    res(2) = rda
    !
    ! Free memory
    CALL fgsl_function_free(stdfunc)
    CALL fgsl_integration_workspace_free(integ_wk)
    !
  CONTAINS
    ! Fortran-to-C wrapper function in order to guarantee interoperability between
    ! Fortran and C
    RECURSIVE FUNCTION f2c_wrapper_func1d_param0(x, params) BIND(c)
      !
      REAL(c_double), VALUE :: x
      TYPE(c_ptr), VALUE :: params
      REAL(c_double) :: f2c_wrapper_func1d_param0
      !
      ! Check, whether C-pointer 'params' is not associated (no parameter passed)
      IF(C_ASSOCIATED(params)) STOP '***Error*** in f2c_wrapper_func1d_param0'
      !
      ! Wrap user-specified function to a C-interoperable function
      f2c_wrapper_func1d_param0 = func1d_param0_user(x)
      !
    END FUNCTION f2c_wrapper_func1d_param0
  END FUNCTION fint1d_param0_qagp
  !--------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------!
  ! Q(uadrature) A(daptive) G(eneral integrand) integration procedure
  ! with known singular points (parameter-valued input function)
  RECURSIVE FUNCTION fint1d_param1_qagp(func1d_param1_user,param1,pts,siz_pts,epsabs,epsrel) result(res)
    !
    INTERFACE  
       FUNCTION func1d_param1_user(x,param1)
         USE fgsl
         REAL(fgsl_double) :: func1d_param1_user
         REAL(fgsl_double) :: x, param1
       END FUNCTION func1d_param1_user
    END INTERFACE
    REAL(fgsl_double), TARGET :: param1
    ! epsabs: absolute error, epsrel: relative error
    REAL(fgsl_double) :: epsabs, epsrel
    ! specify singular points
    REAL(fgsl_double), DIMENSION(:) :: pts
    INTEGER(fgsl_size_t) :: siz_pts
    REAL(fgsl_double), DIMENSION(2) :: res
    !
    INTEGER(fgsl_size_t), PARAMETER :: limit = 1000_fgsl_size_t
    REAL(fgsl_double) :: ra, rda ! result and absolute error
    TYPE(fgsl_error_handler_t) :: std
    INTEGER(fgsl_int) :: status
    TYPE(fgsl_function) :: stdfunc
    TYPE(fgsl_integration_workspace) :: integ_wk
    TYPE(c_ptr) :: param1_ptr ! This pointer holds the C-location of user-specified
    !                         ! parameter
    !
    ! Turn off error handler
    std = fgsl_set_error_handler_off()
    !
    ! This part is replaced by internal subroutines (01.04.2015)
!!$    ! Associate global procedure pointer with user-specified function
!!$    func1d_param1_ptr => func1d_param1_user
!!$    IF(.NOT. ASSOCIATED(func1d_param1_ptr)) &
!!$         STOP '***Error*** in fint1d_param1_qagp'
    !
    ! Get C-location of 'param1'
    param1_ptr = c_LOC(param1)
    IF(.NOT. C_ASSOCIATED(param1_ptr)) &
         STOP '***Error*** in fint1d_param1_qagp'
    !
    ! Initialize fgsl_function and workspace
    stdfunc = fgsl_function_init(f2c_wrapper_func1d_param1, param1_ptr)
    integ_wk = fgsl_integration_workspace_alloc(limit)
    !
    ! Initialize solver 'fgsl_integration_qagp' to use the function 'stdfunc' and
    ! the user-specified parameters
    status = fgsl_integration_qagp(stdfunc, pts, & !siz_pts, &
         epsabs, epsrel, limit, integ_wk, ra, rda)
    CALL check_error(status)
    !
    ! Return the results (ra,rda)
    res(1) = ra
    res(2) = rda
    !
    ! Free memory
    CALL fgsl_function_free(stdfunc)
    CALL fgsl_integration_workspace_free(integ_wk)
    !
  CONTAINS
    ! Fortran-to-C wrapper function in order to guarantee interoperability between
    ! Fortran and C
    RECURSIVE FUNCTION f2c_wrapper_func1d_param1(x, params) BIND(c)
      !
      REAL(c_double), VALUE :: x
      TYPE(c_ptr), VALUE :: params
      REAL(c_double) :: f2c_wrapper_func1d_param1
      !
      REAL(fgsl_double), POINTER :: p ! This is the type of C-pointer to be casted
      !
      ! Check, whether C-pointer 'params' is associated
      IF(.NOT. C_ASSOCIATED(params)) STOP '***Error*** in f2c_wrapper_func1d_param1'
      !
      ! Cast C-pointer to the above-defined Fortran pointer
      CALL C_F_POINTER(params, p)
      ! Wrap user-specified function to a C-interoperable function
      f2c_wrapper_func1d_param1 = func1d_param1_user(x,p)
      !
    END FUNCTION f2c_wrapper_func1d_param1
  END FUNCTION fint1d_param1_qagp
  !--------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------!
  ! CQUAD doubly-adaptive integration
  RECURSIVE FUNCTION fint1d_param0_cquad(func1d_param0_user,x_low,x_up,epsabs,epsrel) result(res)
    !
    INTERFACE  
       FUNCTION func1d_param0_user(x)
         USE fgsl
         REAL(fgsl_double) :: func1d_param0_user
         REAL(fgsl_double) :: x
       END FUNCTION func1d_param0_user
    END INTERFACE
    ! x_low: lower boundary, x_up: upper boundary
    REAL(fgsl_double) :: x_low, x_up
    ! epsabs: absolute error, epsrel: relative error
    REAL(fgsl_double) :: epsabs, epsrel
    REAL(fgsl_double), DIMENSION(2) :: res
    !
    INTEGER(fgsl_size_t), PARAMETER :: limit_cq = 100_fgsl_size_t
    INTEGER(fgsl_size_t) :: neval ! nmber of function evaluations
    REAL(fgsl_double) :: ra, rda ! result and absolute error
    TYPE(fgsl_error_handler_t) :: std
    INTEGER(fgsl_int) :: status
    TYPE(fgsl_function) :: stdfunc
    TYPE(fgsl_integration_cquad_workspace) :: integ_cq
    TYPE(c_ptr) :: param0_ptr ! This pointer holds the C-location of user-specified
    !                         ! parameter (here is no parameter specified -> c_null_ptr)
    !
    ! Turn off error handler
    std = fgsl_set_error_handler_off()
    !
    ! This part is replaced by internal subroutines (01.04.2015)
!!$    ! Associate global procedure pointer with user-specified function
!!$    func1d_param0_ptr => func1d_param0_user
!!$    IF(.NOT. ASSOCIATED(func1d_param0_ptr)) &
!!$         STOP '***Error*** in fint1d_param0_cquad'
    !
    ! Nullify param0_ptr
    param0_ptr = c_null_ptr
    IF(C_ASSOCIATED(param0_ptr)) &
         STOP '***Error*** in fint1d_param0_cquad'
    !
    ! Initialize fgsl_function and workspace
    stdfunc = fgsl_function_init(f2c_wrapper_func1d_param0, param0_ptr)
    integ_cq = fgsl_integration_cquad_workspace_alloc(limit_cq)
    !
    ! Initialize solver 'fgsl_integration_cquad' to use the function 'stdfunc' and
    ! the user-specified parameters
    status = fgsl_integration_cquad(stdfunc, x_low, x_up, &
       epsabs, epsrel, integ_cq, ra, rda, neval)
    CALL check_error(status)
    !
    ! Return the results (ra,rda,neval)
    res(1) = ra
    res(2) = rda
    !fint1d_param0_cquad(3) = REAL(neval,fgsl_double)
    !
    ! Free memory
    CALL fgsl_function_free(stdfunc)
    CALL fgsl_integration_cquad_workspace_free(integ_cq)
    !
  CONTAINS
    ! Fortran-to-C wrapper function in order to guarantee interoperability between
    ! Fortran and C
    RECURSIVE FUNCTION f2c_wrapper_func1d_param0(x, params) BIND(c)
      !
      REAL(c_double), VALUE :: x
      TYPE(c_ptr), VALUE :: params
      REAL(c_double) :: f2c_wrapper_func1d_param0
      !
      ! Check, whether C-pointer 'params' is not associated (no parameter passed)
      IF(C_ASSOCIATED(params)) STOP '***Error*** in f2c_wrapper_func1d_param0'
      !
      ! Wrap user-specified function to a C-interoperable function
      f2c_wrapper_func1d_param0 = func1d_param0_user(x)
      !
    END FUNCTION f2c_wrapper_func1d_param0
  END FUNCTION fint1d_param0_cquad
  !--------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------!
  ! CQUAD doubly-adaptive integration (parameter-valued input function)
  RECURSIVE FUNCTION fint1d_param1_cquad(func1d_param1_user,param1,x_low,x_up,epsabs,epsrel) result(res)
    !
    INTERFACE  
       FUNCTION func1d_param1_user(x,param1)
         USE fgsl
         REAL(fgsl_double) :: func1d_param1_user
         REAL(fgsl_double) :: x, param1
       END FUNCTION func1d_param1_user
    END INTERFACE
    REAL(fgsl_double), TARGET :: param1
    ! x_low: lower boundary, x_up: upper boundary
    REAL(fgsl_double) :: x_low, x_up
    ! epsabs: absolute error, epsrel: relative error
    REAL(fgsl_double) :: epsabs, epsrel
    REAL(fgsl_double), DIMENSION(2) :: res
    !
    INTEGER(fgsl_size_t), PARAMETER :: limit_cq = 1000_fgsl_size_t
    INTEGER(fgsl_size_t) :: neval ! nmber of function evaluations
    REAL(fgsl_double) :: ra, rda ! result and absolute error
    TYPE(fgsl_error_handler_t) :: std
    INTEGER(fgsl_int) :: status
    TYPE(fgsl_function) :: stdfunc
    TYPE(fgsl_integration_cquad_workspace) :: integ_cq
    TYPE(c_ptr) :: param1_ptr ! This pointer holds the C-location of user-specified
    !                         ! parameter 
    !
    ! Turn off error handler
    std = fgsl_set_error_handler_off()
    !
    ! This part is replaced by internal subroutines (01.04.2015)
!!$    ! Associate global procedure pointer with user-specified function
!!$    func1d_param1_ptr => func1d_param1_user
!!$    IF(.NOT. ASSOCIATED(func1d_param1_ptr)) &
!!$         STOP '***Error*** in fint1d_param1_cquad'
    !
    ! Get C-location of 'param1'
    param1_ptr = c_LOC(param1)
    IF(.NOT. C_ASSOCIATED(param1_ptr)) &
         STOP '***Error*** in fint1d_param1_cquad'
    !
    ! Initialize fgsl_function and workspace
    stdfunc = fgsl_function_init(f2c_wrapper_func1d_param1, param1_ptr)
    integ_cq = fgsl_integration_cquad_workspace_alloc(limit_cq)
    !
    ! Initialize solver 'fgsl_integration_cquad' to use the function 'stdfunc' and
    ! the user-specified parameters
    status = fgsl_integration_cquad(stdfunc, x_low, x_up, &
       epsabs, epsrel, integ_cq, ra, rda, neval)
    CALL check_error(status)
    !
    ! Return the results (ra,rda,neval)
    res(1) = ra
    res(2) = rda
    !fint1d_param1_cquad(3) = REAL(neval,fgsl_double)
    !
    ! Free memory
    CALL fgsl_function_free(stdfunc)
    CALL fgsl_integration_cquad_workspace_free(integ_cq)
    !
  CONTAINS
    ! Fortran-to-C wrapper function in order to guarantee interoperability between
    ! Fortran and C
    RECURSIVE FUNCTION f2c_wrapper_func1d_param1(x, params) BIND(c)
      !
      REAL(c_double), VALUE :: x
      TYPE(c_ptr), VALUE :: params
      REAL(c_double) :: f2c_wrapper_func1d_param1
      !
      REAL(fgsl_double), POINTER :: p ! This is the type of C-pointer to be casted
      !
      ! Check, whether C-pointer 'params' is associated
      IF(.NOT. C_ASSOCIATED(params)) STOP '***Error*** in f2c_wrapper_func1d_param1'
      !
      ! Cast C-pointer to the above-defined Fortran pointer
      CALL C_F_POINTER(params, p)
      ! Wrap user-specified function to a C-interoperable function
      f2c_wrapper_func1d_param1 = func1d_param1_user(x,p)
      !
    END FUNCTION f2c_wrapper_func1d_param1
  END FUNCTION fint1d_param1_cquad
  !--------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------!
  ! Q(uadrature) A(daptive) G(eneral integrand) integration procedure for interval from x_low to infinity
  RECURSIVE FUNCTION fint1d_param0_qagiu(func1d_param0_user,x_low,epsabs,epsrel) result(res)
    !
    INTERFACE  
       FUNCTION func1d_param0_user(x)
         USE fgsl
         REAL(fgsl_double) :: func1d_param0_user
         REAL(fgsl_double) :: x
       END FUNCTION func1d_param0_user
    END INTERFACE
    ! x_low: lower boundary, x_up: infinity
    REAL(fgsl_double) :: x_low
    ! epsabs: absolute error, epsrel: relative error
    REAL(fgsl_double) :: epsabs, epsrel
    REAL(fgsl_double), DIMENSION(2) :: res 
    !
    INTEGER(fgsl_size_t), PARAMETER :: limit = 1000_fgsl_size_t
    REAL(fgsl_double) :: ra, rda ! result and absolute error
    TYPE(fgsl_error_handler_t) :: std
    INTEGER(fgsl_int) :: status
    TYPE(fgsl_function) :: stdfunc
    TYPE(fgsl_integration_workspace) :: integ_wk
    TYPE(c_ptr) :: param0_ptr ! This pointer holds the C-location of user-specified
    !                         ! parameter (here is no parameter specified -> c_null_ptr)
    !
    ! Turn off error handler
    std = fgsl_set_error_handler_off()
    !
    ! This part is replaced by internal subroutines (01.04.2015)
!!$    ! Associate global procedure pointer with user-specified function
!!$    func1d_param0_ptr => func1d_param0_user
!!$    IF(.NOT. ASSOCIATED(func1d_param0_ptr)) &
!!$         STOP '***Error*** in fint1d_param0_qag'
    !
    ! Nullify param0_ptr
    param0_ptr = c_null_ptr
    IF(C_ASSOCIATED(param0_ptr)) &
         STOP '***Error*** in fint1d_param0_qagiu'
    !
    ! Initialize fgsl_function and workspace
    stdfunc = fgsl_function_init(f2c_wrapper_func1d_param0, param0_ptr)
    integ_wk = fgsl_integration_workspace_alloc(limit)
    !
    ! Initialize solver 'fgsl_integration_qagiu' to use the function 'stdfunc' and
    ! the user-specified parameters
    status = fgsl_integration_qagiu(stdfunc, x_low, &
       epsabs, epsrel, limit, integ_wk, ra, rda)
    CALL check_error(status)
    !
    ! Return the results (ra,rda)
    res(1) = ra
    res(2) = rda
    !
    ! Free memory
    CALL fgsl_function_free(stdfunc)
    CALL fgsl_integration_workspace_free(integ_wk)
    !
  CONTAINS
    ! Fortran-to-C wrapper function in order to guarantee interoperability between
    ! Fortran and C
    RECURSIVE FUNCTION f2c_wrapper_func1d_param0(x, params) BIND(c)
      !
      REAL(c_double), VALUE :: x
      TYPE(c_ptr), VALUE :: params
      REAL(c_double) :: f2c_wrapper_func1d_param0
      !
      ! Check, whether C-pointer 'params' is not associated (no parameter passed)
      IF(C_ASSOCIATED(params)) STOP '***Error*** in f2c_wrapper_func1d_param0'
      !
      ! Wrap user-specified function to a C-interoperable function
      f2c_wrapper_func1d_param0 = func1d_param0_user(x)
      !
    END FUNCTION f2c_wrapper_func1d_param0
  END FUNCTION fint1d_param0_qagiu
  !--------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------!
  ! Q(uadrature) A(daptive) G(eneral integrand) integration procedure for interval from x_low to infinity
  ! (parameter-valued input function)
  RECURSIVE FUNCTION fint1d_param1_qagiu(func1d_param1_user,param1,x_low,&
       epsabs,epsrel) result(res)
    !
    INTERFACE  
       FUNCTION func1d_param1_user(x,param1)
         USE fgsl
         REAL(fgsl_double) :: func1d_param1_user
         REAL(fgsl_double) :: x, param1
       END FUNCTION func1d_param1_user
    END INTERFACE
    REAL(fgsl_double), TARGET :: param1
    ! x_low: lower boundary, x_up: infinity
    REAL(fgsl_double) :: x_low
    ! epsabs: absolute error, epsrel: relative error
    REAL(fgsl_double) :: epsabs, epsrel
    REAL(fgsl_double), DIMENSION(2) :: res
    !
    INTEGER(fgsl_size_t), PARAMETER :: limit = 1000_fgsl_size_t
    REAL(fgsl_double) :: ra, rda ! result and absolute error
    TYPE(fgsl_error_handler_t) :: std
    INTEGER(fgsl_int) :: status
    TYPE(fgsl_function) :: stdfunc
    TYPE(fgsl_integration_workspace) :: integ_wk
    TYPE(c_ptr) :: param1_ptr ! This pointer holds the C-location of user-specified
    !                         ! parameter
    !
    ! Turn off error handler
    std = fgsl_set_error_handler_off()
    !
    ! This part is replaced by internal subroutines (01.04.2015)
!!$    ! Associate global procedure pointer with user-specified function
!!$    func1d_param1_ptr => func1d_param1_user
!!$    IF(.NOT. ASSOCIATED(func1d_param1_ptr)) &
!!$         STOP '***Error*** in fint1d_param1_qag'
    !
    ! Get C-location of 'param1'
    param1_ptr = c_LOC(param1)
    IF(.NOT. C_ASSOCIATED(param1_ptr)) &
         STOP '***Error*** in fint1d_param1_qag'
    !
    ! Initialize fgsl_function and workspace
    stdfunc = fgsl_function_init(f2c_wrapper_func1d_param1, param1_ptr)
    integ_wk = fgsl_integration_workspace_alloc(limit)
    !
    ! Initialize solver 'fgsl_integration_qag' to use the function 'stdfunc' and
    ! the user-specified parameters
    status = fgsl_integration_qagiu(stdfunc, x_low, &
       epsabs, epsrel, limit, integ_wk, ra, rda)
    CALL check_error(status)
    !
    ! Return the results (ra,rda)
    res(1) = ra
    res(2) = rda
    !
    ! Free memory
    CALL fgsl_function_free(stdfunc)
    CALL fgsl_integration_workspace_free(integ_wk)
    !
  CONTAINS
    ! Fortran-to-C wrapper function in order to guarantee interoperability between
    ! Fortran and C
    RECURSIVE FUNCTION f2c_wrapper_func1d_param1(x, params) BIND(c)
      !
      REAL(c_double), VALUE :: x
      TYPE(c_ptr), VALUE :: params
      REAL(c_double) :: f2c_wrapper_func1d_param1
      !
      REAL(fgsl_double), POINTER :: p ! This is the type of C-pointer to be casted
      !
      ! Check, whether C-pointer 'params' is associated
      IF(.NOT. C_ASSOCIATED(params)) STOP '***Error*** in f2c_wrapper_func1d_param1'
      !
      ! Cast C-pointer to the above-defined Fortran pointer
      CALL C_F_POINTER(params, p)
      ! Wrap user-specified function to a C-interoperable function
      f2c_wrapper_func1d_param1 = func1d_param1_user(x,p)
      !
    END FUNCTION f2c_wrapper_func1d_param1
  END FUNCTION fint1d_param1_qagiu
  !--------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------!
  ! Check error-codes from GSL:
!!$  GSL_SUCCESS  = 0, 
!!$  GSL_FAILURE  = -1,
!!$  GSL_CONTINUE = -2,  /* iteration has not converged */
!!$  GSL_EDOM     = 1,   /* input domain error, e.g sqrt(-1) */
!!$  GSL_ERANGE   = 2,   /* output range error, e.g. exp(1e100) */
!!$  GSL_EFAULT   = 3,   /* invalid pointer */
!!$  GSL_EINVAL   = 4,   /* invalid argument supplied by user */
!!$  GSL_EFAILED  = 5,   /* generic failure */
!!$  GSL_EFACTOR  = 6,   /* factorization failed */
!!$  GSL_ESANITY  = 7,   /* sanity check failed - shouldn't happen */
!!$  GSL_ENOMEM   = 8,   /* malloc failed */
!!$  GSL_EBADFUNC = 9,   /* problem with user-supplied function */
!!$  GSL_ERUNAWAY = 10,  /* iterative process is out of control */
!!$  GSL_EMAXITER = 11,  /* exceeded max number of iterations */
!!$  GSL_EZERODIV = 12,  /* tried to divide by zero */
!!$  GSL_EBADTOL  = 13,  /* user specified an invalid tolerance */
!!$  GSL_ETOL     = 14,  /* failed to reach the specified tolerance */
!!$  GSL_EUNDRFLW = 15,  /* underflow */
!!$  GSL_EOVRFLW  = 16,  /* overflow  */
!!$  GSL_ELOSS    = 17,  /* loss of accuracy */
!!$  GSL_EROUND   = 18,  /* failed because of roundoff error */
!!$  GSL_EBADLEN  = 19,  /* matrix, vector lengths are not conformant */
!!$  GSL_ENOTSQR  = 20,  /* matrix not square */
!!$  GSL_ESING    = 21,  /* apparent singularity detected */
!!$  GSL_EDIVERGE = 22,  /* integral or series is divergent */
!!$  GSL_EUNSUP   = 23,  /* requested feature is not supported by the hardware */
!!$  GSL_EUNIMPL  = 24,  /* requested feature not (yet) implemented */
!!$  GSL_ECACHE   = 25,  /* cache limit exceeded */
!!$  GSL_ETABLE   = 26,  /* table limit exceeded */
!!$  GSL_ENOPROG  = 27,  /* iteration is not making progress towards solution */
!!$  GSL_ENOPROGJ = 28,  /* jacobian evaluations are not improving the solution */
!!$  GSL_ETOLF    = 29,  /* cannot reach the specified tolerance in F */
!!$  GSL_ETOLX    = 30,  /* cannot reach the specified tolerance in X */
!!$  GSL_ETOLG    = 31,  /* cannot reach the specified tolerance in gradient */
!!$  GSL_EOF      = 32   /* end of file */
  !--------------------------------------------------------------------------------------!
  SUBROUTINE check_error(err_gsl)
    ! input
    INTEGER , INTENT(in) :: err_gsl
    !
    IF (err_gsl .EQ. 0) RETURN ! GSL_SUCCESS
    !
    IF ( (err_gsl .LT. -2) .OR.  (err_gsl .GT. 32)) THEN
       PRINT *, 'gsl_integration_routines_mod.f90: Unknown error code from GSL!'
       STOP
    END IF
    !
    gsl_err_detected(err_gsl) = .TRUE.
    !
  END SUBROUTINE check_error
  !--------------------------------------------------------------------------------------!
  SUBROUTINE disp_gsl_integration_error()
    ! internal variables
    INTEGER :: k
    !
    IF ( ANY(gsl_err_detected) ) THEN ! else normal termination
       PRINT *,'-------------------------------------------------'
       DO k = LBOUND(gsl_err_detected,1),UBOUND(gsl_err_detected,1)
          IF (gsl_err_detected(k)) THEN
             PRINT *,"gsl_integration_routines_mod.f90: &
                  &Possible Warning from GSL Integration Routines - Code = ",k
          END IF
       END DO
       PRINT *,'-------------------------------------------------'
    END IF
    !
    ! reset gsl_err_detected
    gsl_err_detected = .FALSE.
    !
  END SUBROUTINE disp_gsl_integration_error
  !--------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------!
  ! Test function without a parameter for 1d integrators
  FUNCTION test_func1d_param0(x)
    !
    REAL(fgsl_double) :: test_func1d_param0
    REAL(fgsl_double) :: x
    !
    test_func1d_param0 = COS(x)
    !
  END FUNCTION test_func1d_param0
  !--------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------!
  ! Test function without a parameter for 1d integrators (singularity)
  FUNCTION test2_func1d_param0(x)
    !
    REAL(fgsl_double) :: test2_func1d_param0
    REAL(fgsl_double) :: x
    !
    test2_func1d_param0 = 1.0_fgsl_double/SQRT(1.0_fgsl_double-x)
    !
  END FUNCTION test2_func1d_param0
  !--------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------!
  ! Test function without a parameter for 1d integrators (singularity)
  FUNCTION test3_func1d_param0(x)
    !
    REAL(fgsl_double) :: test3_func1d_param0
    REAL(fgsl_double) :: x
    !
    test3_func1d_param0 = 1.0_fgsl_double/SQRT(1.0_fgsl_double-x**2.0_fgsl_double)
    !
  END FUNCTION test3_func1d_param0
  !--------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------!
  ! Test function without a parameter for 1d integrators (singularity)
  ! (cf., V. Sladek et al., Appl. Math. Modelling 25 (2001) 901-922)
  FUNCTION test4_func1d_param0(x)
    !
    REAL(fgsl_double) :: test4_func1d_param0
    REAL(fgsl_double) :: x
    !
    test4_func1d_param0 = x*(x-2.0_fgsl_double)*LOG(x/2.0_fgsl_double)
    !
  END FUNCTION test4_func1d_param0
  !--------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------!
  ! Test function without a parameter for 1d integrators (singularity)
  ! (cf., V. Sladek et al., Appl. Math. Modelling 25 (2001) 901-922)
  FUNCTION test5_func1d_param0(x)
    !
    REAL(fgsl_double) :: test5_func1d_param0
    REAL(fgsl_double) :: x
    !
    REAL(fgsl_double), PARAMETER :: d = 15.0_fgsl_double
    !
    test5_func1d_param0 = ((x-d)**3.0_fgsl_double)*LOG(x/2.0_fgsl_double)
    !
  END FUNCTION test5_func1d_param0
  !--------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------!
  ! Test function without a parameter for 1d integrators (singularity)
  FUNCTION test6_func1d_param0(x)
    !
    REAL(fgsl_double) :: test6_func1d_param0
    REAL(fgsl_double) :: x
    !
    test6_func1d_param0 = 1.0_fgsl_double/SQRT(ABS(x))
    !
  END FUNCTION test6_func1d_param0
  !--------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------!
  ! Test function with a parameter for 1d integrators
  FUNCTION test_func1d_param1(x, param1)
    !
    REAL(fgsl_double) :: test_func1d_param1
    REAL(fgsl_double) :: x, param1
    !
    test_func1d_param1 = 2.0_fgsl_double/SQRT(1-param1*(SIN(x)**2.0_fgsl_double))
    !
  END FUNCTION test_func1d_param1
  !--------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------!
  ! Test function with a parameter for 1d integrators
  FUNCTION test2_func1d_param1(x, param1)
    !
    REAL(fgsl_double) :: test2_func1d_param1
    REAL(fgsl_double) :: x, param1
    !
    test2_func1d_param1 = 1.0_fgsl_double/&
         SQRT(param1-(SIN(x/2.0_fgsl_double)**2.0_fgsl_double))
    !
  END FUNCTION test2_func1d_param1
  !--------------------------------------------------------------------------------------!
END MODULE gsl_integration_routines_mod
