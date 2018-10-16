!-----------------------------------------------------------------------------------------!
! module: pspline_routines_mod                                                            !
! author: Andreas F. Martitsch                                                            !
! date: 18.01.2014                                                                        !
! version: 0.1                                                                            !
! description: pspline_routines_mod provides a user-friendly interface to spline routines !
! of PSPLINE (cf., http://w3.pppl.gov/ntcc/PSPLINE/). At the moment only routines for 1d  !
! are supported by this module, which can be readily extended to 2d and 3d.               !
!-----------------------------------------------------------------------------------------!

MODULE pspline_routines_mod
  !
  USE EZspline_obj
  USE EZspline
  use nrtype, only : dp
  !
  IMPLICIT NONE

  public cub1d_pspline_0
  private cub1d_pspline_0_sca, cub1d_pspline_0_vec
  interface cub1d_pspline_0
     module procedure cub1d_pspline_0_sca, cub1d_pspline_0_vec
  end interface cub1d_pspline_0
  !
  public cub1d_pspline_1
  private cub1d_pspline_1_sca, cub1d_pspline_1_vec
  interface cub1d_pspline_1
     module procedure cub1d_pspline_1_sca, cub1d_pspline_1_vec
  end interface cub1d_pspline_1
  !
  public cub1d_pspline_2
  private cub1d_pspline_2_sca, cub1d_pspline_2_vec
  interface cub1d_pspline_2
     module procedure cub1d_pspline_2_sca, cub1d_pspline_2_vec
  end interface cub1d_pspline_2
  !
  public cub1d_pspline_allder1
  private cub1d_pspline_allder1_sca, cub1d_pspline_allder1_vec
  interface cub1d_pspline_allder1
     module procedure cub1d_pspline_allder1_sca, cub1d_pspline_allder1_vec
  end interface cub1d_pspline_allder1
  !
  public cub1d_pspline_allder2
  private cub1d_pspline_allder2_sca, cub1d_pspline_allder2_vec
  interface cub1d_pspline_allder2
     module procedure cub1d_pspline_allder2_sca, cub1d_pspline_allder2_vec
  end interface cub1d_pspline_allder2
  !
CONTAINS
  !--------------------------------------------------------------------------------------!
  ! Evaluate 1d cubic spline (scalar)
  FUNCTION cub1d_pspline_0_sca(xd,yd,x,isw)
    !
    ! Function cub1d_pspline_0_sca interpolates a given set of data (xd,yd) by
    ! a cubic spline at position x (scalar)
    !
    REAL(dp), DIMENSION(:) :: xd,yd
    REAL(dp) :: x
    REAL(dp) :: cub1d_pspline_0_sca
    INTEGER :: isw
    !
    TYPE(EZspline1_dp) :: spl
    INTEGER :: nmax, ier, bcs1(2)
    !
    nmax=SIZE(xd,1)
    !
    ! Select boundary condition type
    SELECT CASE(isw)
       CASE(1)
          bcs1 = (/ 0, 0 /) ! not-a-knot boundary condition
          ! Initialize the spline routine
          ier = 0 ! no error at the beginning
          CALL EZspline_init(spl, nmax, bcs1, ier)
          CALL EZspline_error(ier)
          IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_0_sca'
       CASE(2)
          bcs1 = (/ -1, -1 /) ! periodic boundary condition
          ! Initialize the spline routine
          ier = 0 ! no error at the beginning
          CALL EZspline_init(spl, nmax, bcs1, ier)
          CALL EZspline_error(ier)
          IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_0_sca'
       CASE(3)
          bcs1 = (/ 1, 1/) ! fix 1st derivative
          ! Initialize the spline routine
          ier = 0 ! no error at the beginning
          CALL EZspline_init(spl, nmax, bcs1, ier)
          CALL EZspline_error(ier)
          IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_0_sca'
          spl%bcval1min = 1.0_dp ! value of the 1st derivative
          spl%bcval1max = 1.0_dp ! at the boundary
       CASE(4)
          bcs1 = (/ 2, 2/) ! fix 2nd derivative
          ! Initialize the spline routine
          ier = 0 ! no error at the beginning
          CALL EZspline_init(spl, nmax, bcs1, ier)
          CALL EZspline_error(ier)
          IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_0_sca'
          spl%bcval1min = 0.0_dp ! value of the 2nd derivative at the boundary
          spl%bcval1max = 0.0_dp ! (if zero, natural boundary condition)
       CASE default
          STOP '**ERROR** in cub1d_pspline_0_sca'
    END SELECT
    !
    ! Set data for x1-coordinate
    spl%x1 = xd
    ! Check the grid
    CALL EZspline_isGridRegular(spl, ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_0_sca'
    ! Set data for y1-coordinate
    CALL EZspline_setup(spl, yd, ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_0_sca'
    ! Check, whether point to interpolate (x) is in the domain
    CALL EZspline_isInDomain(spl, x, ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_0_sca'
    ! Evaluate the spline polynomial
    CALL EZspline_interp(spl, x, cub1d_pspline_0_sca, ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_0_sca'
    ! Free memory
    CALL EZspline_free(spl, ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_0_sca'
    !
  END FUNCTION cub1d_pspline_0_sca
  !--------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------!
  ! Evaluate 1d cubic spline (vector)
  FUNCTION cub1d_pspline_0_vec(xd,yd,x,isw)
    !
    ! Function cub1d_pspline_0_vec interpolates a given set of data (xd,yd) by
    ! a cubic spline at position x (vector)
    !
    REAL(dp), DIMENSION(:) :: xd,yd
    REAL(dp), DIMENSION(:) :: x
    REAL(dp), DIMENSION(SIZE(x,1)) :: cub1d_pspline_0_vec
    INTEGER :: isw
    !
    TYPE(EZspline1_dp) :: spl
    INTEGER :: nmax, num_x, ier, bcs1(2)
    !
    nmax=SIZE(xd,1)
    num_x=SIZE(x,1)
    !
    ! Select boundary condition type
    SELECT CASE(isw)
       CASE(1)
          bcs1 = (/ 0, 0 /) ! not-a-knot boundary condition
          ! Initialize the spline routine
          ier = 0 ! no error at the beginning
          CALL EZspline_init(spl, nmax, bcs1, ier)
          CALL EZspline_error(ier)
          IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_0_vec'
       CASE(2)
          bcs1 = (/ -1, -1 /) ! periodic boundary condition
          ! Initialize the spline routine
          ier = 0 ! no error at the beginning
          CALL EZspline_init(spl, nmax, bcs1, ier)
          CALL EZspline_error(ier)
          IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_0_vec'
       CASE(3)
          bcs1 = (/ 1, 1/) ! fix 1st derivative
          ! Initialize the spline routine
          ier = 0 ! no error at the beginning
          CALL EZspline_init(spl, nmax, bcs1, ier)
          CALL EZspline_error(ier)
          IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline__vec'
          spl%bcval1min = 1.0_dp ! value of the 1st derivative
          spl%bcval1max = 1.0_dp ! at the boundary
       CASE(4)
          bcs1 = (/ 2, 2/) ! fix 2nd derivative
          ! Initialize the spline routine
          ier = 0 ! no error at the beginning
          CALL EZspline_init(spl, nmax, bcs1, ier)
          CALL EZspline_error(ier)
          IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_0_vec'
          spl%bcval1min = 0.0_dp ! value of the 2nd derivative at the boundary
          spl%bcval1max = 0.0_dp ! (if zero, natural boundary condition)
       CASE default
          STOP '**ERROR** in cub1d_pspline_0_vec'
    END SELECT
    !
    ! Set data for x1-coordinate
    spl%x1 = xd
    ! Check the grid
    CALL EZspline_isGridRegular(spl, ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_0_vec'
    ! Set data for y1-coordinate
    CALL EZspline_setup(spl, yd, ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_0_vec'
    ! Check, whether point to interpolate (x) is in the domain
    CALL EZspline_isInDomain(spl, num_x, x, ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_0_vec'
    ! Evaluate the spline polynomial
    CALL EZspline_interp(spl, num_x, x, cub1d_pspline_0_vec, ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_0_vec'
    ! Free memory
    CALL EZspline_free(spl, ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_0_vec'
    !
  END FUNCTION cub1d_pspline_0_vec
  !--------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------!
  ! Evaluate derivative from 1d cubic spline (scalar)
  FUNCTION cub1d_pspline_1_sca(xd,yd,x,isw)
    !
    ! Function cub1d_pspline_1_sca interpolates a given set of data (xd,yd) by
    ! a cubic spline and evaluates the derivative at position x (scalar)
    !
    REAL(dp), DIMENSION(:) :: xd,yd
    REAL(dp) :: x
    REAL(dp) :: cub1d_pspline_1_sca
    INTEGER :: isw
    !
    TYPE(EZspline1_dp) :: spl
    INTEGER :: nmax, ier, bcs1(2)
    !
    nmax=SIZE(xd,1)
    !
    ! Select boundary condition type
    SELECT CASE(isw)
       CASE(1)
          bcs1 = (/ 0, 0 /) ! not-a-knot boundary condition
          ! Initialize the spline routine
          ier = 0 ! no error at the beginning
          CALL EZspline_init(spl, nmax, bcs1, ier)
          CALL EZspline_error(ier)
          IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_1_sca'
       CASE(2)
          bcs1 = (/ -1, -1 /) ! periodic boundary condition
          ! Initialize the spline routine
          ier = 0 ! no error at the beginning
          CALL EZspline_init(spl, nmax, bcs1, ier)
          CALL EZspline_error(ier)
          IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_1_sca'
       CASE(3)
          bcs1 = (/ 1, 1/) ! fix 1st derivative
          ! Initialize the spline routine
          ier = 0 ! no error at the beginning
          CALL EZspline_init(spl, nmax, bcs1, ier)
          CALL EZspline_error(ier)
          IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_1_sca'
          spl%bcval1min = 1.0_dp ! value of the 1st derivative
          spl%bcval1max = 1.0_dp ! at the boundary
       CASE(4)
          bcs1 = (/ 2, 2/) ! fix 2nd derivative
          ! Initialize the spline routine
          ier = 0 ! no error at the beginning
          CALL EZspline_init(spl, nmax, bcs1, ier)
          CALL EZspline_error(ier)
          IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_1_sca'
          spl%bcval1min = 0.0_dp ! value of the 2nd derivative at the boundary
          spl%bcval1max = 0.0_dp ! (if zero, natural boundary condition)
       CASE default
          STOP '**ERROR** in cub1d_pspline_1_sca'
    END SELECT
    !
    ! Set data for x1-coordinate
    spl%x1 = xd
    ! Check the grid
    CALL EZspline_isGridRegular(spl, ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_1_sca'
    ! Set data for y1-coordinate
    CALL EZspline_setup(spl, yd, ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_1_sca'
    ! Check, whether point to interpolate (x) is in the domain
    CALL EZspline_isInDomain(spl, x, ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_1_sca'
    ! Evaluate the derivative from the spline polynomial
    CALL EZspline_derivative(spl, 1, x, cub1d_pspline_1_sca, ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_1_sca'
    ! Free memory
    CALL EZspline_free(spl, ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_1_sca'
    !
  END FUNCTION cub1d_pspline_1_sca
  !--------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------!
  ! Evaluate derivative from 1d cubic spline (vector)
  FUNCTION cub1d_pspline_1_vec(xd,yd,x,isw)
    !
    ! Function cub1d_pspline_1_vec interpolates a given set of data (xd,yd) by
    ! a cubic spline and evaluates the derivative at position x (vector)
    !
    REAL(dp), DIMENSION(:) :: xd,yd
    REAL(dp), DIMENSION(:) :: x
    REAL(dp), DIMENSION(SIZE(x,1)) :: cub1d_pspline_1_vec
    INTEGER :: isw
    !
    TYPE(EZspline1_dp) :: spl
    INTEGER :: nmax, num_x, ier, bcs1(2)
    !
    nmax=SIZE(xd,1)
    num_x=SIZE(x,1)
    !
    ! Select boundary condition type
    SELECT CASE(isw)
       CASE(1)
          bcs1 = (/ 0, 0 /) ! not-a-knot boundary condition
          ! Initialize the spline routine
          ier = 0 ! no error at the beginning
          CALL EZspline_init(spl, nmax, bcs1, ier)
          CALL EZspline_error(ier)
          IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_1_vec'
       CASE(2)
          bcs1 = (/ -1, -1 /) ! periodic boundary condition
          ! Initialize the spline routine
          ier = 0 ! no error at the beginning
          CALL EZspline_init(spl, nmax, bcs1, ier)
          CALL EZspline_error(ier)
          IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_1_vec'
       CASE(3)
          bcs1 = (/ 1, 1/) ! fix 1st derivative
          ! Initialize the spline routine
          ier = 0 ! no error at the beginning
          CALL EZspline_init(spl, nmax, bcs1, ier)
          CALL EZspline_error(ier)
          IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_1_vec'
          spl%bcval1min = 1.0_dp ! value of the 1st derivative
          spl%bcval1max = 1.0_dp ! at the boundary
       CASE(4)
          bcs1 = (/ 2, 2/) ! fix 2nd derivative
          ! Initialize the spline routine
          ier = 0 ! no error at the beginning
          CALL EZspline_init(spl, nmax, bcs1, ier)
          CALL EZspline_error(ier)
          IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_1_vec'
          spl%bcval1min = 0.0_dp ! value of the 2nd derivative at the boundary
          spl%bcval1max = 0.0_dp ! (if zero, natural boundary condition)
       CASE default
          STOP '**ERROR** in cub1d_pspline_1_vec'
    END SELECT
    !
    ! Set data for x1-coordinate
    spl%x1 = xd
    ! Check the grid
    CALL EZspline_isGridRegular(spl, ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_1_vec'
    ! Set data for y1-coordinate
    CALL EZspline_setup(spl, yd, ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_1_vec'
    ! Check, whether point to interpolate (x) is in the domain
    CALL EZspline_isInDomain(spl, num_x, x, ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_1_vec'
    ! Evaluate the derivative from the spline polynomial
    CALL EZspline_derivative(spl, 1, num_x, x, cub1d_pspline_1_vec, ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_1_vec'
    ! Free memory
    CALL EZspline_free(spl, ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_1_vec'
    !
  END FUNCTION cub1d_pspline_1_vec
  !--------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------!
  ! Evaluate 2nd derivative from 1d cubic spline (scalar)
  FUNCTION cub1d_pspline_2_sca(xd,yd,x,isw)
    !
    ! Function cub1d_pspline_2_sca interpolates a given set of data (xd,yd) by
    ! a cubic spline and evaluates the 2nd derivative at position x (scalar)
    !
    REAL(dp), DIMENSION(:) :: xd,yd
    REAL(dp) :: x
    REAL(dp) :: cub1d_pspline_2_sca
    INTEGER :: isw
    !
    TYPE(EZspline1_dp) :: spl
    INTEGER :: nmax, ier, bcs1(2)
    !
    nmax=SIZE(xd,1)
    !
    ! Select boundary condition type
    SELECT CASE(isw)
       CASE(1)
          bcs1 = (/ 0, 0 /) ! not-a-knot boundary condition
          ! Initialize the spline routine
          ier = 0 ! no error at the beginning
          CALL EZspline_init(spl, nmax, bcs1, ier)
          CALL EZspline_error(ier)
          IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_2_sca'
       CASE(2)
          bcs1 = (/ -1, -1 /) ! periodic boundary condition
          ! Initialize the spline routine
          ier = 0 ! no error at the beginning
          CALL EZspline_init(spl, nmax, bcs1, ier)
          CALL EZspline_error(ier)
          IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_2_sca'
       CASE(3)
          bcs1 = (/ 1, 1/) ! fix 1st derivative
          ! Initialize the spline routine
          ier = 0 ! no error at the beginning
          CALL EZspline_init(spl, nmax, bcs1, ier)
          CALL EZspline_error(ier)
          IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_2_sca'
          spl%bcval1min = 1.0_dp ! value of the 1st derivative
          spl%bcval1max = 1.0_dp ! at the boundary
       CASE(4)
          bcs1 = (/ 2, 2/) ! fix 2nd derivative
          ! Initialize the spline routine
          ier = 0 ! no error at the beginning
          CALL EZspline_init(spl, nmax, bcs1, ier)
          CALL EZspline_error(ier)
          IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_2_sca'
          spl%bcval1min = 0.0_dp ! value of the 2nd derivative at the boundary
          spl%bcval1max = 0.0_dp ! (if zero, natural boundary condition)
       CASE default
          STOP '**ERROR** in cub1d_pspline_2_sca'
    END SELECT
    !
    ! Set data for x1-coordinate
    spl%x1 = xd
    ! Check the grid
    CALL EZspline_isGridRegular(spl, ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_2_sca'
    ! Set data for y1-coordinate
    CALL EZspline_setup(spl, yd, ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_2_sca'
    ! Check, whether point to interpolate (x) is in the domain
    CALL EZspline_isInDomain(spl, x, ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_2_sca'
    ! Evaluate the derivative from the spline polynomial
    CALL EZspline_derivative(spl, 2, x, cub1d_pspline_2_sca, ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_2_sca'
    ! Free memory
    CALL EZspline_free(spl, ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_2_sca'
    !
  END FUNCTION cub1d_pspline_2_sca
  !--------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------!
  ! Evaluate 2nd derivative from 1d cubic spline (vector)
  FUNCTION cub1d_pspline_2_vec(xd,yd,x,isw)
    !
    ! Function cub1d_pspline_2_vec interpolates a given set of data (xd,yd) by
    ! a cubic spline and evaluates the 2nd derivative at position x (vector)
    !
    REAL(dp), DIMENSION(:) :: xd,yd
    REAL(dp), DIMENSION(:) :: x
    REAL(dp), DIMENSION(SIZE(x,1)) :: cub1d_pspline_2_vec
    INTEGER :: isw
    !
    TYPE(EZspline1_dp) :: spl
    INTEGER :: nmax, num_x, ier, bcs1(2)
    !
    nmax=SIZE(xd,1)
    num_x=SIZE(x,1)
    !
    ! Select boundary condition type
    SELECT CASE(isw)
       CASE(1)
          bcs1 = (/ 0, 0 /) ! not-a-knot boundary condition
          ! Initialize the spline routine
          ier = 0 ! no error at the beginning
          CALL EZspline_init(spl, nmax, bcs1, ier)
          CALL EZspline_error(ier)
          IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_2_vec'
       CASE(2)
          bcs1 = (/ -1, -1 /) ! periodic boundary condition
          ! Initialize the spline routine
          ier = 0 ! no error at the beginning
          CALL EZspline_init(spl, nmax, bcs1, ier)
          CALL EZspline_error(ier)
          IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_2_vec'
       CASE(3)
          bcs1 = (/ 1, 1/) ! fix 1st derivative
          ! Initialize the spline routine
          ier = 0 ! no error at the beginning
          CALL EZspline_init(spl, nmax, bcs1, ier)
          CALL EZspline_error(ier)
          IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_2_vec'
          spl%bcval1min = 1.0_dp ! value of the 1st derivative
          spl%bcval1max = 1.0_dp ! at the boundary
       CASE(4)
          bcs1 = (/ 2, 2/) ! fix 2nd derivative
          ! Initialize the spline routine
          ier = 0 ! no error at the beginning
          CALL EZspline_init(spl, nmax, bcs1, ier)
          CALL EZspline_error(ier)
          IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_2_vec'
          spl%bcval1min = 0.0_dp ! value of the 2nd derivative at the boundary
          spl%bcval1max = 0.0_dp ! (if zero, natural boundary condition)
       CASE default
          STOP '**ERROR** in cub1d_pspline_2_vec'
    END SELECT
    !
    ! Set data for x1-coordinate
    spl%x1 = xd
    ! Check the grid
    CALL EZspline_isGridRegular(spl, ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_2_vec'
    ! Set data for y1-coordinate
    CALL EZspline_setup(spl, yd, ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_2_vec'
    ! Check, whether point to interpolate (x) is in the domain
    CALL EZspline_isInDomain(spl, num_x, x, ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_2_vec'
    ! Evaluate the derivative from the spline polynomial
    CALL EZspline_derivative(spl, 2, num_x, x, cub1d_pspline_2_vec, ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_2_vec'
    ! Free memory
    CALL EZspline_free(spl, ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_2_vec'
    !
  END FUNCTION cub1d_pspline_2_vec
  !--------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------!
  ! Evaluate 1d cubic spline and its derivative (scalar)
  FUNCTION cub1d_pspline_allder1_sca(xd,yd,x,isw)
    !
    ! Function cub1d_pspline_allder1_sca interpolates a given set of data (xd,yd) by
    ! a cubic spline, and computes the value and its derivative at position x (scalar)
    !
    REAL(dp), DIMENSION(:) :: xd,yd
    REAL(dp) :: x
    REAL(dp), DIMENSION(2) :: cub1d_pspline_allder1_sca
    INTEGER :: isw
    !
    TYPE(EZspline1_dp) :: spl
    INTEGER :: nmax, ier, bcs1(2)
    !
    nmax=SIZE(xd,1)
    !
    ! Select boundary condition type
    SELECT CASE(isw)
       CASE(1)
          bcs1 = (/ 0, 0 /) ! not-a-knot boundary condition
          ! Initialize the spline routine
          ier = 0 ! no error at the beginning
          CALL EZspline_init(spl, nmax, bcs1, ier)
          CALL EZspline_error(ier)
          IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_allder1_sca'
       CASE(2)
          bcs1 = (/ -1, -1 /) ! periodic boundary condition
          ! Initialize the spline routine
          ier = 0 ! no error at the beginning
          CALL EZspline_init(spl, nmax, bcs1, ier)
          CALL EZspline_error(ier)
          IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_allder1_sca'
       CASE(3)
          bcs1 = (/ 1, 1/) ! fix 1st derivative
          ! Initialize the spline routine
          ier = 0 ! no error at the beginning
          CALL EZspline_init(spl, nmax, bcs1, ier)
          CALL EZspline_error(ier)
          IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_allder1_sca'
          spl%bcval1min = 1.0_dp ! value of the 1st derivative
          spl%bcval1max = 1.0_dp ! at the boundary
       CASE(4)
          bcs1 = (/ 2, 2/) ! fix 2nd derivative
          ! Initialize the spline routine
          ier = 0 ! no error at the beginning
          CALL EZspline_init(spl, nmax, bcs1, ier)
          CALL EZspline_error(ier)
          IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_allder1_sca'
          spl%bcval1min = 0.0_dp ! value of the 2nd derivative at the boundary
          spl%bcval1max = 0.0_dp ! (if zero, natural boundary condition)
       CASE default
          STOP '**ERROR** in cub1d_pspline_allder1_sca'
    END SELECT
    !
    ! Set data for x1-coordinate
    spl%x1 = xd
    ! Check the grid
    CALL EZspline_isGridRegular(spl, ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_allder1_sca'
    ! Set data for y1-coordinate
    CALL EZspline_setup(spl, yd, ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_allder1_sca'
    ! Check, whether point to interpolate (x) is in the domain
    CALL EZspline_isInDomain(spl, x, ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_allder1_sca'
    ! Evaluate the spline polynomial
    CALL EZspline_interp(spl, x, cub1d_pspline_allder1_sca(1), ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_allder1_sca'
    ! Evaluate the derivative from the spline polynomial
    CALL EZspline_derivative(spl, 1, x, cub1d_pspline_allder1_sca(2), ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_allder1_sca'
    ! Free memory
    CALL EZspline_free(spl, ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_allder1_sca'
    !
  END FUNCTION cub1d_pspline_allder1_sca
  !--------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------!
  ! Evaluate 1d cubic spline and its derivative (vector)
  FUNCTION cub1d_pspline_allder1_vec(xd,yd,x,isw)
    !
    ! Function cub1d_pspline_allder1_vec interpolates a given set of data (xd,yd) by
    ! a cubic spline, and computes the value and its derivative at position x (vector)
    !
    REAL(dp), DIMENSION(:) :: xd,yd
    REAL(dp), DIMENSION(:) :: x
    REAL(dp), DIMENSION(2,SIZE(x,1)) :: cub1d_pspline_allder1_vec
    INTEGER :: isw
    !
    TYPE(EZspline1_dp) :: spl
    INTEGER :: nmax, num_x, ier, bcs1(2)
    !
    nmax=SIZE(xd,1)
    num_x=SIZE(x,1)
    !
    ! Select boundary condition type
    SELECT CASE(isw)
       CASE(1)
          bcs1 = (/ 0, 0 /) ! not-a-knot boundary condition
          ! Initialize the spline routine
          ier = 0 ! no error at the beginning
          CALL EZspline_init(spl, nmax, bcs1, ier)
          CALL EZspline_error(ier)
          IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_allder1_vec'
       CASE(2)
          bcs1 = (/ -1, -1 /) ! periodic boundary condition
          ! Initialize the spline routine
          ier = 0 ! no error at the beginning
          CALL EZspline_init(spl, nmax, bcs1, ier)
          CALL EZspline_error(ier)
          IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_allder1_vec'
       CASE(3)
          bcs1 = (/ 1, 1/) ! fix 1st derivative
          ! Initialize the spline routine
          ier = 0 ! no error at the beginning
          CALL EZspline_init(spl, nmax, bcs1, ier)
          CALL EZspline_error(ier)
          IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_allder1_vec'
          spl%bcval1min = 1.0_dp ! value of the 1st derivative
          spl%bcval1max = 1.0_dp ! at the boundary
       CASE(4)
          bcs1 = (/ 2, 2/) ! fix 2nd derivative
          ! Initialize the spline routine
          ier = 0 ! no error at the beginning
          CALL EZspline_init(spl, nmax, bcs1, ier)
          CALL EZspline_error(ier)
          IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_allder1_vec'
          spl%bcval1min = 0.0_dp ! value of the 2nd derivative at the boundary
          spl%bcval1max = 0.0_dp ! (if zero, natural boundary condition)
       CASE default
          STOP '**ERROR** in cub1d_pspline_allder1_vec'
    END SELECT
    !
    ! Set data for x1-coordinate
    spl%x1 = xd
    ! Check the grid
    CALL EZspline_isGridRegular(spl, ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_allder1_vec'
    ! Set data for y1-coordinate
    CALL EZspline_setup(spl, yd, ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_allde1_vec'
    ! Check, whether point to interpolate (x) is in the domain
    CALL EZspline_isInDomain(spl, num_x, x, ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_allder1_vec'
    ! Evaluate the spline polynomial
    CALL EZspline_interp(spl, num_x, x, cub1d_pspline_allder1_vec(1,:), ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_allder1_vec'
    ! Evaluate the derivative from the spline polynomial
    CALL EZspline_derivative(spl, 1, num_x, x, cub1d_pspline_allder1_vec(2,:), ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_allder1_vec'
    ! Free memory
    CALL EZspline_free(spl, ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_allder1_vec'
    !
  END FUNCTION cub1d_pspline_allder1_vec
  !--------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------!
  ! Evaluate 1d cubic spline, the derivative and the second derivative (scalar)
  FUNCTION cub1d_pspline_allder2_sca(xd,yd,x,isw)
    !
    ! Function cub1d_pspline_allder2_sca interpolates a given set of data (xd,yd) by
    ! a cubic spline, and computes the value, the derivative and the
    ! second derivative at position x (scalar)
    !
    REAL(dp), DIMENSION(:) :: xd,yd
    REAL(dp) :: x
    REAL(dp), DIMENSION(3) :: cub1d_pspline_allder2_sca
    INTEGER :: isw
    !
    TYPE(EZspline1_dp) :: spl
    INTEGER :: nmax, ier, bcs1(2)
    !
    nmax=SIZE(xd,1)
    !
    ! Select boundary condition type
    SELECT CASE(isw)
       CASE(1)
          bcs1 = (/ 0, 0 /) ! not-a-knot boundary condition
          ! Initialize the spline routine
          ier = 0 ! no error at the beginning
          CALL EZspline_init(spl, nmax, bcs1, ier)
          CALL EZspline_error(ier)
          IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_allder2_sca'
       CASE(2)
          bcs1 = (/ -1, -1 /) ! periodic boundary condition
          ! Initialize the spline routine
          ier = 0 ! no error at the beginning
          CALL EZspline_init(spl, nmax, bcs1, ier)
          CALL EZspline_error(ier)
          IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_allder2_sca'
       CASE(3)
          bcs1 = (/ 1, 1/) ! fix 1st derivative
          ! Initialize the spline routine
          ier = 0 ! no error at the beginning
          CALL EZspline_init(spl, nmax, bcs1, ier)
          CALL EZspline_error(ier)
          IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_allder2_sca'
          spl%bcval1min = 1.0_dp ! value of the 1st derivative
          spl%bcval1max = 1.0_dp ! at the boundary
       CASE(4)
          bcs1 = (/ 2, 2/) ! fix 2nd derivative
          ! Initialize the spline routine
          ier = 0 ! no error at the beginning
          CALL EZspline_init(spl, nmax, bcs1, ier)
          CALL EZspline_error(ier)
          IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_allder2_sca'
          spl%bcval1min = 0.0_dp ! value of the 2nd derivative at the boundary
          spl%bcval1max = 0.0_dp ! (if zero, natural boundary condition)
       CASE default
          STOP '**ERROR** in cub1d_pspline_allder2_sca'
    END SELECT
    !
    ! Set data for x1-coordinate
    spl%x1 = xd
    ! Check the grid
    CALL EZspline_isGridRegular(spl, ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_allder2_sca'
    ! Set data for y1-coordinate
    CALL EZspline_setup(spl, yd, ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_allder2_sca'
    ! Check, whether point to interpolate (x) is in the domain
    CALL EZspline_isInDomain(spl, x, ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_allder2_sca'
    ! Evaluate the spline polynomial
    CALL EZspline_interp(spl, x, cub1d_pspline_allder2_sca(1), ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_allder2_sca'
    ! Evaluate the derivative from the spline polynomial
    CALL EZspline_derivative(spl, 1, x, cub1d_pspline_allder2_sca(2), ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_allder2_sca'
    ! Evaluate the derivative from the spline polynomial
    CALL EZspline_derivative(spl, 2, x, cub1d_pspline_allder2_sca(3), ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_allder2_sca'
    ! Free memory
    CALL EZspline_free(spl, ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_allder2_sca'
    !
  END FUNCTION cub1d_pspline_allder2_sca
  !--------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------!
  ! Evaluate 1d cubic spline, the derivative and the second derivative (vector)
  FUNCTION cub1d_pspline_allder2_vec(xd,yd,x,isw)
    !
    ! Function cub1d_pspline_allder2_vec interpolates a given set of data (xd,yd) by
    ! a cubic spline, and computes the value, the derivative and the
    ! second derivative at position x (vector)
    !
    REAL(dp), DIMENSION(:) :: xd,yd
    REAL(dp), DIMENSION(:) :: x
    REAL(dp), DIMENSION(3,SIZE(x,1)) :: cub1d_pspline_allder2_vec
    INTEGER :: isw
    !
    TYPE(EZspline1_dp) :: spl
    INTEGER :: nmax, num_x, ier, bcs1(2)
    !
    nmax=SIZE(xd,1)
    num_x=SIZE(x,1)
    !
    ! Select boundary condition type
    SELECT CASE(isw)
       CASE(1)
          bcs1 = (/ 0, 0 /) ! not-a-knot boundary condition
          ! Initialize the spline routine
          ier = 0 ! no error at the beginning
          CALL EZspline_init(spl, nmax, bcs1, ier)
          CALL EZspline_error(ier)
          IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_allder2_vec'
       CASE(2)
          bcs1 = (/ -1, -1 /) ! periodic boundary condition
          ! Initialize the spline routine
          ier = 0 ! no error at the beginning
          CALL EZspline_init(spl, nmax, bcs1, ier)
          CALL EZspline_error(ier)
          IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_allder2_vec'
       CASE(3)
          bcs1 = (/ 1, 1/) ! fix 1st derivative
          ! Initialize the spline routine
          ier = 0 ! no error at the beginning
          CALL EZspline_init(spl, nmax, bcs1, ier)
          CALL EZspline_error(ier)
          IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_allder2_vec'
          spl%bcval1min = 1.0_dp ! value of the 1st derivative
          spl%bcval1max = 1.0_dp ! at the boundary
       CASE(4)
          bcs1 = (/ 2, 2/) ! fix 2nd derivative
          ! Initialize the spline routine
          ier = 0 ! no error at the beginning
          CALL EZspline_init(spl, nmax, bcs1, ier)
          CALL EZspline_error(ier)
          IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_allder2_vec'
          spl%bcval1min = 0.0_dp ! value of the 2nd derivative at the boundary
          spl%bcval1max = 0.0_dp ! (if zero, natural boundary condition)
       CASE default
          STOP '**ERROR** in cub1d_pspline_allder2_vec'
    END SELECT
    !
    ! Set data for x1-coordinate
    spl%x1 = xd
    ! Check the grid
    CALL EZspline_isGridRegular(spl, ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_allder2_vec'
    ! Set data for y1-coordinate
    CALL EZspline_setup(spl, yd, ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_allder2_vec'
    ! Check, whether point to interpolate (x) is in the domain
    CALL EZspline_isInDomain(spl, num_x, x, ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_allder2_vec'
    ! Evaluate the spline polynomial
    CALL EZspline_interp(spl, num_x, x, cub1d_pspline_allder2_vec(1,:), ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_allder2_vec'
    ! Evaluate the derivative from the spline polynomial
    CALL EZspline_derivative(spl, 1, num_x, x, cub1d_pspline_allder2_vec(2,:), ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_allder2_vec'
    ! Evaluate the derivative from the spline polynomial
    CALL EZspline_derivative(spl, 2, num_x, x, cub1d_pspline_allder2_vec(3,:), ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_allder2_vec'
    ! Free memory
    CALL EZspline_free(spl, ier)
    CALL EZspline_error(ier)
    IF(ier /=0 ) STOP '**ERROR** in cub1d_pspline_allder2_vec'
    !
  END FUNCTION cub1d_pspline_allder2_vec
  !--------------------------------------------------------------------------------------!
END MODULE pspline_routines_mod
