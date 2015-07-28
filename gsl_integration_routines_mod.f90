!-----------------------------------------------------------------------------------------!
! module: gsl_integration_routines_mod                                                    !
! author: Andreas F. Martitsch                                                            !
! date: 23.01.2014                                                                        !
! version: 0.1                                                                            !
! description:                                                                            !
! Module gsl_integration_routines_mod provides routines for integrating 1d functions      !
! The numerical routines used here rely on the GNU Scientific Library (GSL).              !
! The C-interface between the GSL and the Fortran code is provided by FGSL                !
! (See http://www.lrz.de/services/software/mathematik/gsl/fortran/                        !
! for further information)                                                                !
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
  ! declare global procedure pointers (maybe there is a better solution)
  PRIVATE func1d_param0_ptr_int_rout
  PROCEDURE(func1d_param0), POINTER :: func1d_param0_ptr_int_rout => NULL()
  PRIVATE func1d_param1_ptr_int_rout
  PROCEDURE(func1d_param1), POINTER :: func1d_param1_ptr_int_rout => NULL()
  !
  ! Activate messages from root solver
  PUBLIC integration_solver_talk
  LOGICAL :: integration_solver_talk = .FALSE.
  !
  ! Internal Fortran-to-C wrapper functions (public, because of c-binding)
  PUBLIC f2c_wrapper_func1d_param0_int_rout, f2c_wrapper_func1d_param1_int_rout
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
CONTAINS
  !--------------------------------------------------------------------------------------!
  ! Fortran-to-C wrapper function in order to guarantee interoperability between
  ! Fortran and C
  FUNCTION f2c_wrapper_func1d_param0_int_rout(x, params) BIND(c)
    !
    REAL(c_double), VALUE :: x
    TYPE(c_ptr), VALUE :: params
    REAL(c_double) :: f2c_wrapper_func1d_param0_int_rout
    !
    ! Check, whether C-pointer 'params' is not associated (no parameter passed)
    IF(C_ASSOCIATED(params)) STOP '***Error*** in f2c_wrapper_func1d_param0_int_rout'
    !
    ! Wrap user-specified function to a C-interoperable function
    f2c_wrapper_func1d_param0_int_rout = func1d_param0_ptr_int_rout(x)
    !
  END FUNCTION f2c_wrapper_func1d_param0_int_rout
  !--------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------!
  ! Fortran-to-C wrapper function in order to guarantee interoperability between
  ! Fortran and C
  FUNCTION f2c_wrapper_func1d_param1_int_rout(x, params) BIND(c)
    !
    REAL(c_double), VALUE :: x
    TYPE(c_ptr), VALUE :: params
    REAL(c_double) :: f2c_wrapper_func1d_param1_int_rout
    !
    REAL(fgsl_double), POINTER :: p ! This is the type of C-pointer to be casted
    !
    ! Check, whether C-pointer 'params' is associated
    IF(.NOT. C_ASSOCIATED(params)) STOP '***Error*** in f2c_wrapper_func1d_param1_int_rout'
    !
    ! Cast C-pointer to the above-defined Fortran pointer
    CALL C_F_POINTER(params, p)
    ! Wrap user-specified function to a C-interoperable function
    f2c_wrapper_func1d_param1_int_rout = func1d_param1_ptr_int_rout(x,p)
    !
  END FUNCTION f2c_wrapper_func1d_param1_int_rout
  !--------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------!
  ! Q(uadrature) A(daptive) G(eneral integrand) integration procedure
  FUNCTION fint1d_param0_qag(func1d_param0_user,x_low,x_up,epsabs,epsrel,sw_qag_rule)
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
    REAL(fgsl_double), DIMENSION(2) :: fint1d_param0_qag
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
    ! Associate global procedure pointer with user-specified function
    func1d_param0_ptr_int_rout => func1d_param0_user
    IF(.NOT. ASSOCIATED(func1d_param0_ptr_int_rout)) &
         STOP '***Error*** in fint1d_param0_qag'
    !
    ! Nullify param0_ptr
    param0_ptr = c_null_ptr
    IF(C_ASSOCIATED(param0_ptr)) &
         STOP '***Error*** in fint1d_param0_qag'
    !
    ! Initialize fgsl_function and workspace
    stdfunc = fgsl_function_init(f2c_wrapper_func1d_param0_int_rout, param0_ptr)
    integ_wk = fgsl_integration_workspace_alloc(limit)
    !
    ! Initialize solver 'fgsl_integration_qag' to use the function 'stdfunc' and
    ! the user-specified parameters
    status = fgsl_integration_qag(stdfunc, x_low, x_up, &
       epsabs, epsrel, limit, sw_qag_rule, integ_wk, ra, rda)
    !
    ! Return the results (ra,rda)
    fint1d_param0_qag(1) = ra
    fint1d_param0_qag(2) = rda
    !
    ! Free memory
    CALL fgsl_function_free(stdfunc)
    CALL fgsl_integration_workspace_free(integ_wk)
    !
  END FUNCTION fint1d_param0_qag
  !--------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------!
  ! Q(uadrature) A(daptive) G(eneral integrand) integration procedure
  ! (parameter-valued input function)
  FUNCTION fint1d_param1_qag(func1d_param1_user,param1,x_low,x_up,&
       epsabs,epsrel,sw_qag_rule)
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
    REAL(fgsl_double), DIMENSION(2) :: fint1d_param1_qag
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
    ! Associate global procedure pointer with user-specified function
    func1d_param1_ptr_int_rout => func1d_param1_user
    IF(.NOT. ASSOCIATED(func1d_param1_ptr_int_rout)) &
         STOP '***Error*** in fint1d_param1_qag'
    !
    ! Get C-location of 'param1'
    param1_ptr = c_loc(param1)
    IF(.NOT. C_ASSOCIATED(param1_ptr)) &
         STOP '***Error*** in fint1d_param1_qag'
    !
    ! Initialize fgsl_function and workspace
    stdfunc = fgsl_function_init(f2c_wrapper_func1d_param1_int_rout, param1_ptr)
    integ_wk = fgsl_integration_workspace_alloc(limit)
    !
    ! Initialize solver 'fgsl_integration_qag' to use the function 'stdfunc' and
    ! the user-specified parameters
    status = fgsl_integration_qag(stdfunc, x_low, x_up, &
       epsabs, epsrel, limit, sw_qag_rule, integ_wk, ra, rda)
    !
    ! Return the results (ra,rda)
    fint1d_param1_qag(1) = ra
    fint1d_param1_qag(2) = rda
    !
    ! Free memory
    CALL fgsl_function_free(stdfunc)
    CALL fgsl_integration_workspace_free(integ_wk)
    !
  END FUNCTION fint1d_param1_qag
  !--------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------!
  ! Q(uadrature) A(daptive) G(eneral integrand) integration procedure
  ! with singularities
  FUNCTION fint1d_param0_qags(func1d_param0_user,x_low,x_up,epsabs,epsrel)
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
    REAL(fgsl_double), DIMENSION(2) :: fint1d_param0_qags
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
    ! Associate global procedure pointer with user-specified function
    func1d_param0_ptr_int_rout => func1d_param0_user
    IF(.NOT. ASSOCIATED(func1d_param0_ptr_int_rout)) &
         STOP '***Error*** in fint1d_param0_qags'
    !
    ! Nullify param0_ptr
    param0_ptr = c_null_ptr
    IF(C_ASSOCIATED(param0_ptr)) &
         STOP '***Error*** in fint1d_param0_qags'
    !
    ! Initialize fgsl_function and workspace
    stdfunc = fgsl_function_init(f2c_wrapper_func1d_param0_int_rout, param0_ptr)
    integ_wk = fgsl_integration_workspace_alloc(limit)
    !
    ! Initialize solver 'fgsl_integration_qags' to use the function 'stdfunc' and
    ! the user-specified parameters
    status = fgsl_integration_qags(stdfunc, x_low, x_up, &
       epsabs, epsrel, limit, integ_wk, ra, rda)
    !
    ! Return the results (ra,rda)
    fint1d_param0_qags(1) = ra
    fint1d_param0_qags(2) = rda
    !
    ! Free memory
    CALL fgsl_function_free(stdfunc)
    CALL fgsl_integration_workspace_free(integ_wk)
    !
  END FUNCTION fint1d_param0_qags
  !--------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------!
  ! Q(uadrature) A(daptive) G(eneral integrand) integration procedure
  ! with singularities (parameter-valued input function)
  FUNCTION fint1d_param1_qags(func1d_param1_user,param1,x_low,x_up,epsabs,epsrel)
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
    REAL(fgsl_double), DIMENSION(2) :: fint1d_param1_qags
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
    ! Associate global procedure pointer with user-specified function
    func1d_param1_ptr_int_rout => func1d_param1_user
    IF(.NOT. ASSOCIATED(func1d_param1_ptr_int_rout)) &
         STOP '***Error*** in fint1d_param1_qags'
    !
    ! Get C-location of 'param1'
    param1_ptr = c_loc(param1)
    IF(.NOT. C_ASSOCIATED(param1_ptr)) &
         STOP '***Error*** in fint1d_param1_qags'
    !
    ! Initialize fgsl_function and workspace
    stdfunc = fgsl_function_init(f2c_wrapper_func1d_param1_int_rout, param1_ptr)
    integ_wk = fgsl_integration_workspace_alloc(limit)
    !
    ! Initialize solver 'fgsl_integration_qags' to use the function 'stdfunc' and
    ! the user-specified parameters
    status = fgsl_integration_qags(stdfunc, x_low, x_up, &
       epsabs, epsrel, limit, integ_wk, ra, rda)
    !
    ! Return the results (ra,rda)
    fint1d_param1_qags(1) = ra
    fint1d_param1_qags(2) = rda
    !
    ! Free memory
    CALL fgsl_function_free(stdfunc)
    CALL fgsl_integration_workspace_free(integ_wk)
    !
  END FUNCTION fint1d_param1_qags
  !--------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------!
  ! Q(uadrature) A(daptive) G(eneral integrand) integration procedure
  ! with known singular points
  FUNCTION fint1d_param0_qagp(func1d_param0_user,pts,siz_pts,epsabs,epsrel)
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
    REAL(fgsl_double), DIMENSION(2) :: fint1d_param0_qagp
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
    ! Associate global procedure pointer with user-specified function
    func1d_param0_ptr_int_rout => func1d_param0_user
    IF(.NOT. ASSOCIATED(func1d_param0_ptr_int_rout)) &
         STOP '***Error*** in fint1d_param0_qagp'
    !
    ! Nullify param0_ptr
    param0_ptr = c_null_ptr
    IF(C_ASSOCIATED(param0_ptr)) &
         STOP '***Error*** in fint1d_param0_qagp'
    !
    ! Initialize fgsl_function and workspace
    stdfunc = fgsl_function_init(f2c_wrapper_func1d_param0_int_rout, param0_ptr)
    integ_wk = fgsl_integration_workspace_alloc(limit)
    !
    ! Initialize solver 'fgsl_integration_qagp' to use the function 'stdfunc' and
    ! the user-specified parameters
    status = fgsl_integration_qagp(stdfunc, pts, siz_pts, &
         epsabs, epsrel, limit, integ_wk, ra, rda)
    !
    ! Return the results (ra,rda)
    fint1d_param0_qagp(1) = ra
    fint1d_param0_qagp(2) = rda
    !
    ! Free memory
    CALL fgsl_function_free(stdfunc)
    CALL fgsl_integration_workspace_free(integ_wk)
    !
  END FUNCTION fint1d_param0_qagp
  !--------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------!
  ! Q(uadrature) A(daptive) G(eneral integrand) integration procedure
  ! with known singular points (parameter-valued input function)
  FUNCTION fint1d_param1_qagp(func1d_param1_user,param1,pts,siz_pts,epsabs,epsrel)
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
    REAL(fgsl_double), DIMENSION(2) :: fint1d_param1_qagp
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
    ! Associate global procedure pointer with user-specified function
    func1d_param1_ptr_int_rout => func1d_param1_user
    IF(.NOT. ASSOCIATED(func1d_param1_ptr_int_rout)) &
         STOP '***Error*** in fint1d_param1_qagp'
    !
    ! Get C-location of 'param1'
    param1_ptr = c_loc(param1)
    IF(.NOT. C_ASSOCIATED(param1_ptr)) &
         STOP '***Error*** in fint1d_param1_qagp'
    !
    ! Initialize fgsl_function and workspace
    stdfunc = fgsl_function_init(f2c_wrapper_func1d_param1_int_rout, param1_ptr)
    integ_wk = fgsl_integration_workspace_alloc(limit)
    !
    ! Initialize solver 'fgsl_integration_qagp' to use the function 'stdfunc' and
    ! the user-specified parameters
    status = fgsl_integration_qagp(stdfunc, pts, siz_pts, &
         epsabs, epsrel, limit, integ_wk, ra, rda)
    !
    ! Return the results (ra,rda)
    fint1d_param1_qagp(1) = ra
    fint1d_param1_qagp(2) = rda
    !
    ! Free memory
    CALL fgsl_function_free(stdfunc)
    CALL fgsl_integration_workspace_free(integ_wk)
    !
  END FUNCTION fint1d_param1_qagp
  !--------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------!
  ! CQUAD doubly-adaptive integration
  FUNCTION fint1d_param0_cquad(func1d_param0_user,x_low,x_up,epsabs,epsrel)
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
    REAL(fgsl_double), DIMENSION(3) :: fint1d_param0_cquad
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
    ! Associate global procedure pointer with user-specified function
    func1d_param0_ptr_int_rout => func1d_param0_user
    IF(.NOT. ASSOCIATED(func1d_param0_ptr_int_rout)) &
         STOP '***Error*** in fint1d_param0_cquad'
    !
    ! Nullify param0_ptr
    param0_ptr = c_null_ptr
    IF(C_ASSOCIATED(param0_ptr)) &
         STOP '***Error*** in fint1d_param0_cquad'
    !
    ! Initialize fgsl_function and workspace
    stdfunc = fgsl_function_init(f2c_wrapper_func1d_param0_int_rout, param0_ptr)
    integ_cq = fgsl_integration_cquad_workspace_alloc(limit_cq)
    !
    ! Initialize solver 'fgsl_integration_cquad' to use the function 'stdfunc' and
    ! the user-specified parameters
    status = fgsl_integration_cquad(stdfunc, x_low, x_up, &
       epsabs, epsrel, integ_cq, ra, rda, neval)
    !
    ! Return the results (ra,rda,neval)
    fint1d_param0_cquad(1) = ra
    fint1d_param0_cquad(2) = rda
    fint1d_param0_cquad(3) = REAL(neval,fgsl_double)
    !
    ! Free memory
    CALL fgsl_function_free(stdfunc)
    CALL fgsl_integration_cquad_workspace_free(integ_cq)
    !
  END FUNCTION fint1d_param0_cquad
  !--------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------!
  ! CQUAD doubly-adaptive integration (parameter-valued input function)
  FUNCTION fint1d_param1_cquad(func1d_param1_user,param1,x_low,x_up,epsabs,epsrel)
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
    REAL(fgsl_double), DIMENSION(3) :: fint1d_param1_cquad
    !
    INTEGER(fgsl_size_t), PARAMETER :: limit_cq = 100_fgsl_size_t
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
    ! Associate global procedure pointer with user-specified function
    func1d_param1_ptr_int_rout => func1d_param1_user
    IF(.NOT. ASSOCIATED(func1d_param1_ptr_int_rout)) &
         STOP '***Error*** in fint1d_param1_cquad'
    !
    ! Get C-location of 'param1'
    param1_ptr = c_loc(param1)
    IF(.NOT. C_ASSOCIATED(param1_ptr)) &
         STOP '***Error*** in fint1d_param1_cquad'
    !
    ! Initialize fgsl_function and workspace
    stdfunc = fgsl_function_init(f2c_wrapper_func1d_param1_int_rout, param1_ptr)
    integ_cq = fgsl_integration_cquad_workspace_alloc(limit_cq)
    !
    ! Initialize solver 'fgsl_integration_cquad' to use the function 'stdfunc' and
    ! the user-specified parameters
    status = fgsl_integration_cquad(stdfunc, x_low, x_up, &
       epsabs, epsrel, integ_cq, ra, rda, neval)
    !
    ! Return the results (ra,rda,neval)
    fint1d_param1_cquad(1) = ra
    fint1d_param1_cquad(2) = rda
    fint1d_param1_cquad(3) = REAL(neval,fgsl_double)
    !
    ! Free memory
    CALL fgsl_function_free(stdfunc)
    CALL fgsl_integration_cquad_workspace_free(integ_cq)
    !
  END FUNCTION fint1d_param1_cquad
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
