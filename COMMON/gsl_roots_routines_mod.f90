!-----------------------------------------------------------------------------------------!
! module: gsl_roots_routines_mod                                                          !
! author: Andreas F. Martitsch                                                            !
! date: 23.01.2014                                                                        !
! version: 0.1                                                                            !
! description:                                                                            !
! Module gsl_roots_routines_mod provides routines for finding roots of arbitrary 1d       !
! functions. The numerical routines used here rely on the GNU Scientific Library (GSL).   !
! The C-interface between the GSL and the Fortran code is provided by FGSL                !
! (See http://www.lrz.de/services/software/mathematik/gsl/fortran/                        !
! for further information)                                                                !
!-----------------------------------------------------------------------------------------!

MODULE gsl_roots_routines_mod
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
     FUNCTION func1d_param1(x,param1) ! used for bisection method
       USE fgsl
       REAL(fgsl_double) :: func1d_param1
       REAL(fgsl_double) :: x, param1
     END FUNCTION func1d_param1
     FUNCTION fdfunc1d_param1(x,param1) ! used for Newton method
       USE fgsl
       REAL(fgsl_double), DIMENSION(2) :: fdfunc1d_param1 ! 1st entry: f, 2nd entry: df
       REAL(fgsl_double) :: x, param1
     END FUNCTION fdfunc1d_param1
  END INTERFACE
  ! declare global procedure pointers (maybe there is a better solution)
  PRIVATE func1d_param1_ptr
  PROCEDURE(func1d_param1), POINTER :: func1d_param1_ptr => NULL()
  PRIVATE fdfunc1d_param1_ptr
  PROCEDURE(fdfunc1d_param1), POINTER :: fdfunc1d_param1_ptr => NULL()
  !
  ! Activate messages from root solver
  PUBLIC root_solver_talk
  LOGICAL :: root_solver_talk = .FALSE.
  !
  ! Internal Fortran-to-C wrapper functions (public, because of c-binding)
  PUBLIC f2c_wrapper_func1d_param1, f2c_wrapper_newton_func1d_param1, &
       f2c_wrapper_newton_dfunc1d_param1, f2c_wrapper_newton_fdfunc1d_param1
  !
  ! Selection of test cases
  PUBLIC test_func1d_param1, test_fdfunc1d_param1
  !
  ! Find root of a 1d function by a bisection algorithm
  PUBLIC fzero1d_bisec
  PRIVATE fzero1d_param1_bisec
  INTERFACE fzero1d_bisec
     module procedure fzero1d_param1_bisec
  END INTERFACE fzero1d_bisec
  !
  ! Find root of a 1d function by a Newton algorithm
  PUBLIC fzero1d_newton
  PRIVATE fzero1d_param1_newton
  INTERFACE fzero1d_newton
     module procedure fzero1d_param1_newton
  END INTERFACE fzero1d_newton
  !
CONTAINS
  !--------------------------------------------------------------------------------------!
  ! Fortran-to-C wrapper function in order to guarantee interoperability between
  ! Fortran and C
  FUNCTION f2c_wrapper_func1d_param1(x, params) BIND(c)
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
    f2c_wrapper_func1d_param1 = func1d_param1_ptr(x,p)
    !
  END FUNCTION f2c_wrapper_func1d_param1
  !--------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------!
  ! Fortran-to-C wrapper function in order to guarantee interoperability between
  ! Fortran and C
  FUNCTION f2c_wrapper_newton_func1d_param1(x, params) BIND(c)
    !
    REAL(c_double), VALUE :: x
    TYPE(c_ptr), VALUE :: params
    REAL(c_double) :: f2c_wrapper_newton_func1d_param1
    !
    REAL(fgsl_double), POINTER :: p ! This is the type of C-pointer to be casted
    REAL(fgsl_double), DIMENSION(2) :: fdf_val_temp
    !
    ! Check, whether C-pointer 'params' is associated
    IF(.NOT. C_ASSOCIATED(params)) STOP '***Error*** in f2c_wrapper_dfunc1d_param1'
    !
    ! Cast C-pointer to the above-defined Fortran pointer
    CALL C_F_POINTER(params, p)
    ! Wrap user-specified function to a C-interoperable function
    fdf_val_temp = fdfunc1d_param1_ptr(x,p)
    f2c_wrapper_newton_func1d_param1 = fdf_val_temp(1)
    !
  END FUNCTION f2c_wrapper_newton_func1d_param1
  !--------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------!
  ! Fortran-to-C wrapper function in order to guarantee interoperability between
  ! Fortran and C
  FUNCTION f2c_wrapper_newton_dfunc1d_param1(x, params) BIND(c)
    !
    REAL(c_double), VALUE :: x
    TYPE(c_ptr), VALUE :: params
    REAL(c_double) :: f2c_wrapper_newton_dfunc1d_param1
    !
    REAL(fgsl_double), POINTER :: p ! This is the type of C-pointer to be casted
    REAL(fgsl_double), DIMENSION(2) :: fdf_val_temp
    !
    ! Check, whether C-pointer 'params' is associated
    IF(.NOT. C_ASSOCIATED(params)) STOP '***Error*** in f2c_wrapper_dfunc1d_param1'
    !
    ! Cast C-pointer to the above-defined Fortran pointer
    CALL C_F_POINTER(params, p)
    ! Wrap user-specified function to a C-interoperable function
    fdf_val_temp = fdfunc1d_param1_ptr(x,p)
    f2c_wrapper_newton_dfunc1d_param1 = fdf_val_temp(2)
    !
  END FUNCTION f2c_wrapper_newton_dfunc1d_param1
  !--------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------!
  ! Fortran-to-C wrapper function in order to guarantee interoperability between
  ! Fortran and C
  SUBROUTINE f2c_wrapper_newton_fdfunc1d_param1(x, params, y, dy) BIND(c)
    !
    REAL(c_double), VALUE :: x
    TYPE(c_ptr), VALUE :: params
    REAL(c_double), INTENT(out) :: y, dy
    !
    REAL(fgsl_double), POINTER :: p ! This is the type of C-pointer to be casted
    REAL(fgsl_double), DIMENSION(2) :: fdf_val_temp
    !
    ! Check, whether C-pointer 'params' is associated
    IF(.NOT. C_ASSOCIATED(params)) STOP '***Error*** in f2c_wrapper_fdfunc1d_param1'
    !
    ! Cast C-pointer to the above-defined Fortran pointer
    CALL C_F_POINTER(params, p)
    ! Wrap user-specified function to a C-interoperable function
    fdf_val_temp = fdfunc1d_param1_ptr(x,p)
    y = fdf_val_temp(1)
    dy = fdf_val_temp(2)
    !
  END SUBROUTINE f2c_wrapper_newton_fdfunc1d_param1
  !--------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------!
  ! Find root of a 1d function by a bisection algorithm
  FUNCTION fzero1d_param1_bisec(func1d_param1_user,param1,x_low,x_up)
    !
    INTERFACE
       FUNCTION func1d_param1_user(x,param1)
         USE fgsl
         REAL(fgsl_double) :: func1d_param1_user
         REAL(fgsl_double) :: x, param1
       END FUNCTION func1d_param1_user
    END INTERFACE
    REAL(fgsl_double), TARGET :: param1
    REAL(fgsl_double) :: x_low, x_up
    REAL(fgsl_double) :: fzero1d_param1_bisec
    !
    INTEGER(fgsl_int), PARAMETER :: itmax_root = 50
    INTEGER(fgsl_int) :: status, i
    REAL(fgsl_double) :: xlo, xup
    TYPE(fgsl_function) :: stdfunc
    TYPE(fgsl_root_fsolver) :: root_fslv
    TYPE(c_ptr) :: param1_ptr ! This pointer holds the C-location of 'param1'
    !
    ! Associate global procedure pointer with user-specified function
    func1d_param1_ptr => func1d_param1_user
    IF(.NOT. ASSOCIATED(func1d_param1_ptr)) &
         STOP '***Error*** in fzero1d_bisec'
    !
    ! Get C-location of 'param1'
    param1_ptr = c_loc(param1)
    IF(.NOT. C_ASSOCIATED(param1_ptr)) STOP '***Error*** in fzero1d_bisec'
    !
    ! Initialize fgsl_function and solver routine
    stdfunc = fgsl_function_init(f2c_wrapper_func1d_param1, param1_ptr)
    root_fslv = fgsl_root_fsolver_alloc(fgsl_root_fsolver_bisection)
    ! Initialize solver 'root_fslv' to use the function 'stdfunc' and
    ! the initial search interval [x_low,x_up]
    status = fgsl_root_fsolver_set(root_fslv, stdfunc, x_low, x_up)
    ! Start iteration
    i = 0
    DO
       i = i + 1
       ! Perform iteration
       status = fgsl_root_fsolver_iterate(root_fslv)
       IF (status /= fgsl_success .OR. i > itmax_root) THEN
          PRINT *, '***Error*** in fzero1d_bisec'
          PRINT *, 'An error occurred during iteration or'
          PRINT *, 'maximum number of iterations exceeded'
          EXIT ! exit, if an error happened or
               ! if max. number of iterations exceeded
       END IF
       ! Return current root
       fzero1d_param1_bisec = fgsl_root_fsolver_root(root_fslv)
       ! Return current search interval
       xlo = fgsl_root_fsolver_x_lower(root_fslv)
       xup = fgsl_root_fsolver_x_upper(root_fslv)
       ! Check for convergence
       status = fgsl_root_test_interval (xlo, xup, 0.0_fgsl_double, eps5)
       IF (status == fgsl_success) EXIT ! in case of convergence, exit
    END DO
    ! Print number of iterations
    IF (root_solver_talk) THEN
       PRINT *, 'number of iterations: ',i
    END IF
    !
    ! Free memory
    CALL fgsl_root_fsolver_free(root_fslv)
    CALL fgsl_function_free(stdfunc)
    !
  END FUNCTION fzero1d_param1_bisec
  !--------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------!
  ! Find root of a 1d function by a Newton algorithm
  FUNCTION fzero1d_param1_newton(fdfunc1d_param1_user,param1,x0)
    !
    INTERFACE
       FUNCTION fdfunc1d_param1_user(x,param1)
         USE fgsl
         REAL(fgsl_double), DIMENSION(2) :: fdfunc1d_param1_user
         REAL(fgsl_double) :: x, param1
       END FUNCTION fdfunc1d_param1_user
    END INTERFACE
    REAL(fgsl_double), TARGET :: param1
    REAL(fgsl_double) :: x0
    REAL(fgsl_double) :: fzero1d_param1_newton
    !
    INTEGER(fgsl_int), PARAMETER :: itmax_root = 50
    INTEGER(fgsl_int) :: status, i
    REAL(fgsl_double) :: ra, ri
    TYPE(fgsl_function_fdf) :: stdfunc_fdf
    TYPE(fgsl_root_fdfsolver) :: root_fdfslv
    TYPE(c_ptr) :: param1_ptr ! This pointer holds the C-location of 'param1'
    !
    ! Associate global procedure pointer with user-specified function
    fdfunc1d_param1_ptr => fdfunc1d_param1_user
    IF(.NOT. ASSOCIATED(fdfunc1d_param1_ptr)) &
         STOP '***Error*** in fzero1d_newton'
    !
    ! Get C-location of 'param1'
    param1_ptr = c_loc(param1)
    IF(.NOT. C_ASSOCIATED(param1_ptr)) STOP '***Error*** in fzero1d_newton'
    !
    ! Initialize fgsl_function and solver routine
    stdfunc_fdf = fgsl_function_fdf_init(f2c_wrapper_newton_func1d_param1, &
         f2c_wrapper_newton_dfunc1d_param1, f2c_wrapper_newton_fdfunc1d_param1,&
         param1_ptr)
    root_fdfslv = fgsl_root_fdfsolver_alloc(fgsl_root_fdfsolver_newton)
    ! Initialize solver 'root_fdfslv' to use the function 'stdfunc_fdf' and
    ! the initial guess 'x0'
    status = fgsl_root_fdfsolver_set(root_fdfslv, stdfunc_fdf, x0)
    ! Start iteration
    i = 0
    ra = 1.0E100_fgsl_double
    DO
       i = i + 1
       ! Perform iteration
       status = fgsl_root_fdfsolver_iterate(root_fdfslv)
       IF (status /= fgsl_success .OR. i > itmax_root) THEN
          PRINT *, '***Error*** in fzero1d_newton'
          PRINT *, 'An error occurred during iteration or'
          PRINT *, 'maximum number of iterations exceeded'
          EXIT ! exit, if an error happened or
               ! if max. number of iterations exceeded
       END IF
       ! Return current root
       ri = ra
       ra = fgsl_root_fdfsolver_root(root_fdfslv)
       ! Check for convergence
       status = fgsl_root_test_delta (ra, ri, 0.0_fgsl_double, eps5)
       IF (status == fgsl_success) EXIT ! in case of convergence, exit
    END DO
    !
    ! Return the root
    fzero1d_param1_newton = ra
    ! Print number of iterations
    IF (root_solver_talk) THEN
       PRINT *, 'number of iterations: ',i
    END IF
    !
    ! Free memory
    CALL fgsl_root_fdfsolver_free(root_fdfslv)
    CALL fgsl_function_fdf_free(stdfunc_fdf)
    !
  END FUNCTION fzero1d_param1_newton
  !--------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------!
  ! 1d test function with a parameter for bisection method
  FUNCTION test_func1d_param1(x, param1)
    !
    REAL(fgsl_double) :: test_func1d_param1
    REAL(fgsl_double) :: x, param1
    !
    test_func1d_param1 = COS(param1*x)
    !
  END FUNCTION test_func1d_param1
  !--------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------!
  ! 1d test function with a parameter for Newton method
  FUNCTION test_fdfunc1d_param1(x, param1)
    !
    REAL(fgsl_double), DIMENSION(2) :: test_fdfunc1d_param1
    REAL(fgsl_double) :: x, param1
    !
    test_fdfunc1d_param1(1) = COS(param1*x)
    test_fdfunc1d_param1(2) = -SIN(param1*x)
    !
  END FUNCTION test_fdfunc1d_param1
  !--------------------------------------------------------------------------------------!
END MODULE gsl_roots_routines_mod
