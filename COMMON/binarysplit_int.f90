MODULE binarysplit_int

  USE hdf5_tools
  USE hdf5_tools_f2003
  IMPLICIT NONE

  PRIVATE dp
  INTEGER, PARAMETER    :: dp = KIND(1.0d0)

  PUBLIC bsfunc_message
  INTEGER               :: bsfunc_message = 0

  INTEGER, PUBLIC       :: bsfunc_modelfunc = 1


  PRIVATE bsfunc_evaldegree
  INTEGER               :: bsfunc_evaldegree = 2 ! 1   ! 3 or 4

  PRIVATE longint
  INCLUDE 'longint.f90'

  ! parameters for evaluation functions
  ! have to be set with (.e.g. for gauss)
  !   CALL construct_bsfunc(x0,s)
  PRIVATE bsfunc_p1,bsfunc_p2
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: bsfunc_p1,bsfunc_p2

  ! parameters for limit functions
  ! have to be set with
  !   CALL create_bslimit( ..... )
  ! optional arguments
  !   total_err,local_err,min_distance,max_index,max_splitlevel
  REAL(kind=dp), PUBLIC                    :: bsfunc_total_err
  REAL(kind=dp), PUBLIC                    :: bsfunc_local_err
  REAL(kind=dp), PUBLIC                    :: bsfunc_min_distance
  INTEGER,       PUBLIC                    :: bsfunc_max_index
  INTEGER,       PUBLIC                    :: bsfunc_max_splitlevel

  REAL(kind=dp), PUBLIC                    :: bsfunc_base_distance = 0.05d0
  REAL(kind=dp), PUBLIC                    :: bsfunc_mult_constant = 0.1d0
  REAL(kind=dp), PUBLIC                    :: bsfunc_sigma_multiplier = 1.147202690439877 !sqrt(3.0d0)




  ! additional quantities which belon logically to splitting of levels
  ! but they are not used internally
  REAL(kind=dp), PUBLIC                    :: bsfunc_sigma_mult
  REAL(kind=dp), PUBLIC                    :: bsfunc_sigma_min
  INTEGER,       PUBLIC                    :: bsfunc_local_solver

  PUBLIC construct_bsfunc
  PRIVATE construct_bsfunc_2
  INTERFACE construct_bsfunc
    MODULE PROCEDURE construct_bsfunc_2
  END INTERFACE

  PUBLIC destruct_bsfunc
  PRIVATE destruct_bsfunc_2
  INTERFACE destruct_bsfunc
    MODULE PROCEDURE destruct_bsfunc_2
  END INTERFACE

  PUBLIC eval_bsfunc
  PRIVATE eval_bsfunc_d1, eval_bsfunc_d2
  INTERFACE eval_bsfunc
    MODULE PROCEDURE eval_bsfunc_d1, eval_bsfunc_d2
  END INTERFACE

  PUBLIC eval_bsfitsplit
  PRIVATE eval_bsfitsplit_4
  INTERFACE eval_bsfitsplit
    MODULE PROCEDURE eval_bsfitsplit_4
  END INTERFACE

  PUBLIC eval_bsfitsplit_external
  PRIVATE eval_bsfitsplit_ext
  INTERFACE eval_bsfitsplit_external
    MODULE PROCEDURE eval_bsfitsplit_ext
  END INTERFACE

  PUBLIC eval_bsfitjoin_external
  PRIVATE eval_bsfitjoin_ext
  INTERFACE eval_bsfitjoin_external
    MODULE PROCEDURE eval_bsfitjoin_ext
  END INTERFACE

  PUBLIC linspace
  PRIVATE linspace_d
  INTERFACE linspace
    MODULE PROCEDURE linspace_d
  END INTERFACE

  PUBLIC eval_bsinterr
  PRIVATE eval_bsinterr_d1
  INTERFACE eval_bsinterr
    MODULE PROCEDURE eval_bsinterr_d1,eval_bsinterr_test
  END INTERFACE

  PUBLIC eval_bslimit
  PRIVATE eval_bslimit_1, eval_bslimit_3
  INTERFACE eval_bslimit
    MODULE PROCEDURE eval_bslimit_1, eval_bslimit_3
  END INTERFACE

  PUBLIC create_bslimit
  PRIVATE create_bslimit_opt
  INTERFACE create_bslimit
    MODULE PROCEDURE create_bslimit_opt
  END INTERFACE

  PUBLIC  disp
  PRIVATE disp_i8, disp_i81, disp_c
  INTERFACE disp
    MODULE PROCEDURE disp_i8, disp_i81, disp_c
  END INTERFACE

  PRIVATE gauss
  PRIVATE gauss_d0, gauss_d1
  INTERFACE gauss
    MODULE PROCEDURE gauss_d0, gauss_d1
  END INTERFACE

  PRIVATE plagrange_coeff
  PRIVATE plag_coeff
  INTERFACE plagrange_coeff
    MODULE PROCEDURE plag_coeff
  END INTERFACE

  PRIVATE plagrange_stencil
  PRIVATE plag_stencil
  INTERFACE plagrange_stencil
    MODULE PROCEDURE plag_coeff
  END INTERFACE

  PUBLIC plagrange_test
  PRIVATE plag_test
  INTERFACE plagrange_test
    MODULE PROCEDURE plag_test
  END INTERFACE

CONTAINS

  ! construct
  ! puts the parameters p1, p2
  ! into the private variables bsfunc_p1, .....
  ! which are used internally
  SUBROUTINE construct_bsfunc_2(p1,p2)
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(in) :: p1,p2

    IF (ALLOCATED(bsfunc_p1)) DEALLOCATE(bsfunc_p1)
    ALLOCATE(bsfunc_p1(SIZE(p1,1)))
    bsfunc_p1 = p1

    IF (ALLOCATED(bsfunc_p2)) DEALLOCATE(bsfunc_p2)
    ALLOCATE(bsfunc_p2(SIZE(p2,1)))
    bsfunc_p2 = p2
  END SUBROUTINE construct_bsfunc_2

  ! destruct
  SUBROUTINE destruct_bsfunc_2()
    IF (ALLOCATED(bsfunc_p1)) DEALLOCATE(bsfunc_p1)
    IF (ALLOCATED(bsfunc_p2)) DEALLOCATE(bsfunc_p2)
  END SUBROUTINE destruct_bsfunc_2

  ! evaluation
  ! here the gauss function is called for the 2 parameters
  !  x0 (bsfunc_p1) and s (bsfunc_p1)
  !  has to be changed if another function has to be used
  SUBROUTINE eval_bsfunc_d1(x,y)
    REAL(kind=dp), INTENT(in)    :: x
    REAL(kind=dp), INTENT(inout) :: y
    CALL gauss(x,bsfunc_p1,bsfunc_p2,y)
  END SUBROUTINE eval_bsfunc_d1

  SUBROUTINE eval_bsfunc_d2(x,y)
    REAL(kind=dp), DIMENSION(:), INTENT(in)    :: x ! 0:
    REAL(kind=dp), DIMENSION(:), INTENT(inout) :: y ! 0:
    INTEGER                                    :: k
    DO k = LBOUND(x,1), UBOUND(x,1)
      CALL eval_bsfunc(x(k),y(k))
    END DO
  END SUBROUTINE eval_bsfunc_d2

  ! procedure when you split
  ! handles all tasks when a level is split into two
  !
  !  internal stuff: function value and integeral
  !  with polynomial of degree 3
  !
  ! should be adopted to the interr
  SUBROUTINE eval_bsfitsplit_4(x,y,i,e,loc)
    REAL(kind=dp), DIMENSION(0:), INTENT(in)    :: x ! 0:
    REAL(kind=dp), DIMENSION(0:), INTENT(inout) :: y ! 0:
    REAL(kind=dp), DIMENSION(0:), INTENT(inout) :: i ! 0:
    REAL(kind=dp), DIMENSION(0:), INTENT(inout) :: e ! 0:
    INTEGER,                      INTENT(in)   :: loc

    REAL(kind=dp) :: dx,dx2,dx3,dx4,dy,dydx1,dydx2,b,c

    ! internal stuff (parabola of degree 3)
    IF (loc .GT. 2) THEN
      dydx1 = (y(loc-1)-y(loc-2)) / (x(loc-1)-x(loc-2))
    ELSE
      dydx1 = (y(loc)-y(loc-1)) / (x(loc)-x(loc-1))
    END IF
    dx    = x(loc+1) - x(loc-1)
    dy    = y(loc+1) - y(loc-1)
    dydx2 = dy / dx
    dx2   = dx * dx
    dx3   = dx2 * dx
    b = -(dydx2-dydx1)/dx  +3.0_dp*(dy-dydx1*dx)/dx2
    c = +(dydx2-dydx1)/dx2 -2.0_dp*(dy-dydx1*dx)/dx3

    dx  = x(loc) - x(loc-1)
    dx2 = dx * dx
    dx3 = dx2 * dx
    dx4 = dx2 * dx2
    y(loc) = y(loc-1) + dydx1*dx + b*dx2 + c*dx3
    i(loc) = y(loc-1)*dx + dydx1*dx2 + b*dx3 + c*dx4
    i(loc+1) = i(loc+1) - i(loc)
    e(loc) = e(loc+1)/2.0_dp
    e(loc+1) = e(loc)
    CALL eval_bsfitsplit_external(x,loc)
  END SUBROUTINE eval_bsfitsplit_4

  ! external procedure when you split
  ! has to be adopted for your problem
  SUBROUTINE eval_bsfitsplit_ext(x,loc)
    REAL(kind=dp), DIMENSION(0:), INTENT(in)    :: x !0:
    INTEGER,                     INTENT(in)    :: loc
    INTERFACE join_ripples_bsfitsplit
      SUBROUTINE join_ripples_bsfitsplit(x,loc)
        INTEGER, PARAMETER :: dp = KIND(1.0d0)
        REAL(kind=dp), DIMENSION(0:), INTENT(in)    :: x !0:
        INTEGER,                     INTENT(in)    :: loc
      END SUBROUTINE join_ripples_bsfitsplit
    END INTERFACE
    ! add here what is necessary for external stuff
    CALL join_ripples_bsfitsplit(x,loc) !EXTERNAL
  END SUBROUTINE eval_bsfitsplit_ext

  ! external procedure when you join
  ! has to be adopted for your problem
  SUBROUTINE eval_bsfitjoin_ext(x,loc)
    REAL(kind=dp), DIMENSION(0:), INTENT(in)    :: x !0:
    INTEGER,                      INTENT(in)   :: loc
    INTERFACE join_ripples_bsfitjoin
      SUBROUTINE join_ripples_bsfitjoin(x,loc)
        INTEGER, PARAMETER :: dp = KIND(1.0d0)
        REAL(kind=dp), DIMENSION(0:), INTENT(in)    :: x !0:
        INTEGER,                      INTENT(in)    :: loc
      END SUBROUTINE join_ripples_bsfitjoin
    END INTERFACE
    ! add here what is necessary for external stuff
    CALL join_ripples_bsfitjoin(x,loc) !EXTERNAL
  END SUBROUTINE eval_bsfitjoin_ext

  ! limits
  SUBROUTINE create_bslimit_opt(total_err,local_err,min_distance, &
       max_index,max_splitlevel)
    REAL(kind=dp), INTENT(in), OPTIONAL :: total_err
    REAL(kind=dp), INTENT(in), OPTIONAL :: local_err
    REAL(kind=dp), INTENT(in), OPTIONAL :: min_distance
    INTEGER,       INTENT(in), OPTIONAL :: max_index
    INTEGER,       INTENT(in), OPTIONAL :: max_splitlevel

    IF (PRESENT(total_err)) THEN
      bsfunc_total_err = total_err
    ELSE
      bsfunc_total_err = 0.0_dp
    END IF
    IF (PRESENT(local_err)) THEN
      bsfunc_local_err = local_err
    ELSE
      bsfunc_local_err = 0.0_dp
    END IF
    IF (PRESENT(min_distance)) THEN
      bsfunc_min_distance = min_distance
    ELSE
      bsfunc_min_distance = 0.0_dp
    END IF
    IF (PRESENT(max_index)) THEN
      bsfunc_max_index = max_index
    ELSE
      bsfunc_max_index = 1000
    END IF
    IF (PRESENT(max_splitlevel)) THEN
      bsfunc_max_splitlevel = max_splitlevel
    ELSE
      bsfunc_max_splitlevel = 32
    END IF
  END SUBROUTINE create_bslimit_opt

  FUNCTION eval_bslimit_1(err) RESULT(l)
    REAL(kind=dp), DIMENSION(0:), INTENT(in) :: err  !0:
    LOGICAL                                 :: l, l1, l2
    !l1 = SUM(err,1) .GT. bsfunc_total_err
    l2 = MAXVAL(err,1) .GT. bsfunc_local_err
    !PRINT *, 'Local Error ',MAXVAL(err,1),bsfunc_local_err
    !l  = l1 .OR. l2
    l = l2
!!$    IF (.NOT. l1) THEN
!!$       IF (bsfunc_message .GT. 0) THEN
!!$          PRINT *, 'Message: Total Error Limits reached!'
!!$       END IF
!!$    END IF
    IF (.NOT. l2) THEN
      IF (bsfunc_message .GT. 0) THEN
        PRINT *, 'Message: Local Error Limits reached!'
      END IF
    END IF
  END FUNCTION eval_bslimit_1

  FUNCTION eval_bslimit_3(dist,split,maxind) RESULT(l)
    REAL(kind=dp),                INTENT(in) :: dist
    INTEGER,                      INTENT(in) :: split
    INTEGER,                      INTENT(in) :: maxind
    LOGICAL                                  :: l, l1, l2, l3
    l1 = dist   .GE. bsfunc_min_distance
    l2 = split  .LE. bsfunc_max_splitlevel
    l3 = maxind .LE. bsfunc_max_index
    l  = l1 .AND. l2 .AND. l3
    IF (.NOT. l1) THEN
      IF (bsfunc_message .GT. 0) PRINT *, 'Message: Minimum distance reached!'
    END IF
    IF (.NOT. l2) THEN
      IF (bsfunc_message .GT. 0) PRINT *, 'Message: Maximum split level reached!'
    END IF
    IF (.NOT. l3) THEN
      PRINT *, 'Message: Maximum index reached!'
      STOP
    END IF
  END FUNCTION eval_bslimit_3

  !> Equivalent to matlab/octave linspace, with output as last parameter.
  !>
  !> This subroutine is an equivalent to the matlab/octave function
  !> linspace, just the return value is here the last parameter.
  SUBROUTINE linspace_d(x_beg,x_end,x_num,x)
    REAL(kind=dp),                            INTENT(in)     :: x_beg
    REAL(kind=dp),                            INTENT(in)     :: x_end
    INTEGER,                                  INTENT(in)     :: x_num
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(inout)  :: x !0:

    INTEGER                                                  :: i_beg_l
    INTEGER                                                  :: i_end_l
    INTEGER                                                  :: k,i
    REAL(kind=dp)                                            :: x_del

    i_beg_l = 0
    i_end_l = x_num + i_beg_l -1
    IF (x_num .LE. 1) i_end_l = i_beg_l
    IF (ALLOCATED(x)) DEALLOCATE(x)
    ALLOCATE (x(i_beg_l:i_end_l))

    IF (x_num .LE. 1) THEN
      x(i_beg_l) = x_beg
    ELSE
      x_del = (x_end-x_beg) / DBLE(x_num-1)

      i = i_beg_l
      x(i) = x_beg
      DO k = 2, x_num-1
        i = i + 1
        x(i) = x(i-1) + x_del
      END DO
      x(i_end_l) = x_end
    END IF
  END SUBROUTINE linspace_d


  ! gauss (helper routine)
  SUBROUTINE gauss_d0(x,x0,s,g)
    REAL(kind=dp),               INTENT(in)   :: x
    REAL(kind=dp), DIMENSION(:), INTENT(in)   :: x0
    REAL(kind=dp), DIMENSION(:), INTENT(in)   :: s
    REAL(kind=dp),               INTENT(out)  :: g

    REAL(kind=dp),               PARAMETER    :: pi=3.14159265358979_dp
    REAL(kind=dp),               PARAMETER    :: os2pi=0.39894228040143_dp
    REAL(kind=dp),               PARAMETER    :: os8=0.35355339059327_dp
    REAL(kind=dp),               PARAMETER    :: sq2=1.41421356237310_dp
    INTEGER                                   :: k

    g = 0.0_dp
    IF (bsfunc_modelfunc .EQ. 1) THEN
      ! Winny - This is now the return of the normalized Gauss-Function
      ! changed to -> NO NORMALIZATION
      DO k = 1, SIZE(x0,1)
        g = g + EXP(- (x-x0(k))**2 / s(k)**2 / 2.0_dp) ! / s(k)
      END DO
      ! g = g * os2pi

!!$       IF (bsfunc_evaldegree .NE. 2) THEN
!!$          DO k = 1, SIZE(x0,1)
!!$             g = g + EXP(- ((x-x0(k)) / s(k))**2 / 2) / s(k)
!!$          END DO
!!$          g = g * os2pi
!!$       ELSE
!!$          DO k = 1, SIZE(x0,1)
!!$             g = g + EXP(- ((x-x0(k)) / s(k))**2 / 2)
!!$          END DO
!!$       END IF
    ELSEIF (bsfunc_modelfunc .EQ. 2) THEN
      DO k = 1, SIZE(x0,1)
        g = g + SQRT( s(k)**2/(2.0_dp*(x-x0(k))**2+s(k)**2) )
      END DO
       !g = g/size(x0,1)
!       IF (bsfunc_evaldegree .NE. 2) THEN
!          DO k = 1, SIZE(x0,1)
!             !->out          g = g + EXP(- ABS(x-x0(k)) / s(k) / sq2) / s(k)
!             g=g+s(k)/((x-x0(k))**2+s(k)**2/0.12732200375004d0)              !<-in
!          END DO
!          !->out       g = g * os8
!          g = g / (pi*0.35682208977310d0)
!       ELSE
!          DO k = 1, SIZE(x0,1)
!             g=g+s(k)/((x-x0(k))**2+s(k)**2/0.12732200375004d0) / 0.12732200375004d0 * s(k)
!          END DO
!       END IF
    ELSEIF (bsfunc_modelfunc .EQ. 3) THEN
      IF (bsfunc_evaldegree .NE. 2) THEN
        DO k = 1, SIZE(x0,1)
          !->out          g = g + s(k) / ( (x-x0(k))**2 + s(k)**2 )
          g=g+s(k)/((x-x0(k))**2+s(k)**2/0.87267799624997d0)              !<-in
        END DO
        !->out       g = g / pi
        g = g / (pi*0.93417235896272)                                      !<-in
      ELSE
        DO k = 1, SIZE(x0,1)
          g=g+s(k)/((x-x0(k))**2+s(k)**2/0.87267799624997d0) / 0.87267799624997d0 * s(k)
        END DO
      END IF
    ELSEIF (bsfunc_modelfunc .EQ. 4) THEN
      g = 1.0_dp
      DO k = 1, SIZE(x0,1)
        g = g * (1.0_dp + EXP(- ((x-x0(k)) / s(k))**2 / 2) / s(k))
      END DO
    ELSEIF (bsfunc_modelfunc .EQ. 5) THEN
      g = 1.0_dp
      DO k = 1, SIZE(x0,1)
        !g = g * (1.0_dp + EXP(- ((x-x0(k)) / s(k))**2 / 2) / s(k))
         g = g * SIGN(x-x0(k),1.0d0) * (s(k) * ABS(x-x0(k)))**(8.0d0/8.0d0)
      END DO
      g = (s(k) * ABS(x0(k)))**(8.0d0/8.0d0) - g
    ELSE
      PRINT *,'Error from binarysplit: bsfunc_modelfunc wrong: ', &
            bsfunc_modelfunc
      STOP
    END IF
  END SUBROUTINE gauss_d0

  SUBROUTINE gauss_d1(x,x0,s,g)
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(in)    :: x !0:
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE,  INTENT(in)    :: x0
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE,  INTENT(in)    :: s
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(inout) :: g !0:

    INTEGER                                     :: k

    DO k = LBOUND(x,1), UBOUND(x,1)
      CALL gauss(x(k),x0,s,g(k))
    END DO
  END SUBROUTINE gauss_d1

  ! evaluation of error for splitting
  SUBROUTINE eval_bsinterr_d1(x,y,int0,err,splitloc)
    !USE plagrange_mod
    REAL(kind=dp), DIMENSION(0:), INTENT(in)    :: x !0:
    REAL(kind=dp), DIMENSION(0:), INTENT(in)    :: y !0:
    REAL(kind=dp), DIMENSION(0:), INTENT(inout) :: int0 !0:
    REAL(kind=dp), DIMENSION(0:), INTENT(inout) :: err !0:
    INTEGER,                      INTENT(out)   :: splitloc

    INTEGER                                   :: lb, ub, k, i
    INTEGER                                   :: in, in1, in2
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE  :: yp, ypp
    REAL(kind=dp)                             :: dx, dx2, dx3, dx4, dx5
    REAL(kind=dp)                             :: dy, a, b, c, d, e, aa, bb
    REAL(kind=dp)                             :: xloc,yloc,sint0,fac

    INTEGER :: k1,k2,i1,i2,ii
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE  :: xlag, ylag
    REAL(kind=dp) :: err_int

    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: coeff
    REAL(kind=dp) :: clag,dlag
    REAL(kind=dp) :: dloc,d1,kd
    INTEGER, PARAMETER :: npoi = 6
    INTEGER, PARAMETER :: nder = 3



    lb = LBOUND(x,1)
    ub = UBOUND(x,1)
    int0(lb) = 0.0_dp
    err(lb)  = 0.0_dp
    ALLOCATE(yp(lb:ub))
    ALLOCATE(ypp(lb:ub))

    !PRINT *, 'eval_bsinterr_d1'
    IF (bsfunc_evaldegree .EQ. 1) THEN
      sint0 = 0.0_dp
      in = 8
      in1 = 8
      in2 = 64 ! around max
      fac = 2.0_dp
      !DO WHILE (ABS(sint0) .LT. 1.0d-4)
        !PRINT *, 'ub, in ',ub,in
      DO k = 1, ub
        dx = x(k)-x(k-1)
        IF (x(k) .LT. bsfunc_p1(1) - fac*bsfunc_p2(1) .OR. &
          & x(k-1) .GT. bsfunc_p1(1) + fac*bsfunc_p2(1) ) THEN
          int0(k) = 0.0_dp
          DO i = 0, in1
            xloc = x(k-1) + DBLE(i) * dx / DBLE(in1)
            CALL eval_bsfunc(xloc,yloc)
            IF (i .EQ. 0 .OR. i .EQ. in1) yloc = yloc / 2.0_dp
            int0(k) = int0(k) + yloc
          END DO
          int0(k) = int0(k) * dx / DBLE(in1)
        ELSEIF (x(k-1) .LE. bsfunc_p1(1) - fac*bsfunc_p2(1) .AND. &
          x(k) .GE. bsfunc_p1(1) + fac*bsfunc_p2(1) ) THEN
          int0(k) = 0.0_dp
          DO i = 0, in2
            xloc = x(k-1) + DBLE(i) * dx / DBLE(in2)
            CALL eval_bsfunc(xloc,yloc)
            IF (i .EQ. 0 .OR. i .EQ. in2) yloc = yloc / 2.0_dp
            int0(k) = int0(k) + yloc
          END DO
          int0(k) = int0(k) * dx / DBLE(in2)
          IF (int0(k) .LT. 1.d-4) THEN
            int0(k) = 1.0_dp
            !PRINT *, 'set to 1.0'
          END IF
        ELSE
          int0(k) = 0.0_dp
          DO i = 0, in
            xloc = x(k-1) + DBLE(i) * dx / DBLE(in)
            CALL eval_bsfunc(xloc,yloc)
            IF (i .EQ. 0 .OR. i .EQ. in) yloc = yloc / 2.0_dp
            int0(k) = int0(k) + yloc
          END DO
          int0(k) = int0(k) * dx / DBLE(in)
        END IF
        err(k)  = (y(k-1)+y(k))*dx/2.0_dp
        err(k)  = ABS(err(k)-int0(k))
      END DO
      sint0 = SUM(int0)
      !PRINT *, 'int0'
      !PRINT *, int0
      !PRINT *, 'err'
      !PRINT *, err
      !PRINT *, 'sint0 ',sint0,in
      !PRINT *, 'max err ',MAXVAL(err)
      !PAUSE
      !   in = in * 10
      !END DO

    ELSEIF (bsfunc_evaldegree .EQ. 2) THEN
      ! new stuff with third derivative
      err = 0.0_dp
      !print *,'lb,ub ',lbound(err,1),ubound(err,1)
      int0 = 0.0_dp
      ALLOCATE( coeff(0:nder,npoi) )
      DO k = 1, ub

        CALL plag_stencil(ub,npoi,k,k1,k2,i1,i2)
        !print *, ub,npoi,k,k1,k2,i1,i2
        ALLOCATE(xlag(i1:i2))
        ALLOCATE(ylag(i1:i2))
        xlag = x(k1:k2)
        ylag = y(k1:k2)

        xloc = (xlag(0) + xlag(1)) / 2.0_dp
        CALL plagrange_coeff(npoi,nder,xloc,xlag,coeff)

        dloc = ABS(xloc-xlag(0))
        d1   = ABS(xlag(1)-xlag(0))

        dlag = SUM(coeff(3,:)*ylag) / 6.0_dp
        ! ori
        clag = ( SUM(coeff(2,:)*ylag) - dlag*6.0_dp*dloc ) / 2.0_dp
        ! modified
        ! clag = ( SUM(coeff(2,:)*ylag) ) / 2.0_dp
        ! ori
        err(k) = ( ABS(clag) * d1**2 + ABS(dlag) * d1**3 ) / 2.0_dp
        ! modified
        ! err(k) = ( ABS(clag) * d1**2 + ABS(dlag) * d1**3 )
        ! err(k) = ABS( clag*d1**2 + dlag*d1**3 )
        ! err(k) = max( ABS(clag) * d1**2 , ABS(dlag) * d1**3 )

        ! relative error for second model_func
        IF (bsfunc_modelfunc .EQ. 2) THEN
          err(k) = err(k) /  ( ABS( SUM(coeff(0,:)*ylag) ) )
        END IF
        !print *, 'k,cla,dlag,err ',k,ABS(clag),ABS(dlag),d1,err(k)
        !print *, 'xlag ',xlag
        !print *, 'ylag ',ylag
        DEALLOCATE(xlag)
        DEALLOCATE(ylag)
      END DO
      !PRINT *, ''
      !PRINT *,'x'
      !PRINT *,x
      !PRINT *, 'y'
      !PRINT *,y
      !PRINT *,'err'
      !PRINT *,err
      !print *, 'sumerr ',sum(err)
      !print *, 'ub splitloc maxerr', ubound(err,1),maxloc(err)-1,maxval(err)
      !PRINT *, ''
      !pause
      DEALLOCATE(coeff)
      ! end of new stuff with third derivative

    ELSEIF (bsfunc_evaldegree .EQ. 3) THEN
      DO k = 1, ub
        yp(k) = (y(k)-y(k-1))/(x(k)-x(k-1))
      END DO
      yp(0) = yp(1)

      DO k = 1, ub
        dx = x(k)-x(k-1)
        dx2 = dx * dx
        dx3 = dx2 * dx
        dx4 = dx2 * dx2

        a = y(k-1)
        b = yp(k-1)

        aa = y(k) - b*dx - a
        bb = (yp(k) - b)*dx

        c = ( 3.0_dp*aa - bb) / dx2
        d = (-2.0_dp*aa + bb) / dx3

        int0(k) = d*dx4/4.0_dp + c*dx3/3.0_dp + b*dx2/2.0_dp + a*dx
        err(k)  = (y(k-1)+y(k))*dx/2.0_dp
        err(k)  = ABS(err(k)-int0(k))
      END DO
    ELSE IF (bsfunc_evaldegree .EQ. 4) THEN
      DO k = 1, ub
        yp(k) = (y(k)-y(k-1))/(x(k)-x(k-1))
      END DO
      yp(0) = yp(1)
      DO k = 1, ub
        ypp(k) = (yp(k)-yp(k-1))/(x(k)-x(k-1))
      END DO
      ypp(0) = ypp(1)


      DO k = 1, ub
        dx = x(k)-x(k-1)
        dx2 = dx * dx
        dx3 = dx2 * dx
        dx4 = dx2 * dx2
        dx5 = dx4 * dx

        a = y(k-1)
        b = yp(k-1)
        c = ypp(k-1) / 2.0_dp

        aa = y(k) - c*dx2 - b*dx - a
        bb = (yp(k) - 2.0_dp*c*dx - b)*dx

        d = ( 4.0_dp*aa - bb) / dx3
        e = (-3.0_dp*aa + bb) / dx4

        int0(k) = e*dx5/5.0_dp + d*dx4/4.0_dp + c*dx3/3.0_dp + b*dx2/2.0_dp + a*dx
        err(k)  = (y(k-1)+y(k))*dx/2.0_dp
        err(k)  = ABS(err(k)-int0(k))
      END DO
    END IF
    splitloc = MAXLOC(err,1)-1
    !PRINT *, 'loc ',splitloc, err(splitloc), SUM(err),x(splitloc)
    !PRINT *, ' '
    !PAUSE
    DEALLOCATE(yp)
    DEALLOCATE(ypp)
  END SUBROUTINE eval_bsinterr_d1

  ! just for some test, puts back interpolated values
  SUBROUTINE eval_bsinterr_test(x,y,int0,err)
    REAL(kind=dp), DIMENSION(0:), INTENT(inout) :: x !0:
    REAL(kind=dp), DIMENSION(0:), INTENT(inout) :: y !0:
    REAL(kind=dp), DIMENSION(0:), INTENT(inout) :: int0 !0:
    REAL(kind=dp), DIMENSION(0:), INTENT(inout) :: err !0:

    INTEGER                                   :: lb, ub, k
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE  :: yp, ypp
    REAL(kind=dp)                             :: dx, dx2, dx3, dx4, dx5
    REAL(kind=dp)                             :: dy, a, b, c, d, e, aa, bb

    lb = LBOUND(x,1)
    ub = UBOUND(x,1)
    int0(lb) = 0.0_dp
    err(lb)  = 0.0_dp
    ALLOCATE(yp(lb:ub))
    ALLOCATE(ypp(lb:ub))

    IF (bsfunc_evaldegree .EQ. 3) THEN
      DO k = 1, ub
        yp(k) = (y(k)-y(k-1))/(x(k)-x(k-1))
      END DO
      yp(0) = yp(1)

      DO k = 1, ub
        dx = x(k)-x(k-1)
        dx2 = dx * dx
        dx3 = dx2 * dx
        dx4 = dx2 * dx2

        a = y(k-1)
        b = yp(k-1)

        aa = y(k) - b*dx - a
        bb = (yp(k) - b)*dx

        c = ( 3.0_dp*aa - bb) / dx2
        d = (-2.0_dp*aa + bb) / dx3

        int0(k) = d*dx4/4.0_dp + c*dx3/3.0_dp + b*dx2/2.0_dp + a*dx
        err(k)  = (y(k-1)+y(k))*dx/2.0_dp
        err(k)  = ABS(err(k)-int0(k))

        dx = dx / 2.0_dp
        dx2 = dx * dx
        dx3 = dx2 * dx
        x(k-1) = x(k-1) + dx
        y(k-1) = d*dx3 + c*dx2 + b*dx + a
      END DO
    ELSE IF (bsfunc_evaldegree .EQ. 4 .OR. bsfunc_evaldegree .EQ. 1) THEN
      DO k = 1, ub
        yp(k) = (y(k)-y(k-1))/(x(k)-x(k-1))
      END DO
      yp(0) = yp(1)
      DO k = 1, ub
        ypp(k) = (yp(k)-yp(k-1))/(x(k)-x(k-1))
      END DO
      ypp(0) = ypp(1)


      DO k = 1, ub
        dx = x(k)-x(k-1)
        dx2 = dx * dx
        dx3 = dx2 * dx
        dx4 = dx2 * dx2
        dx5 = dx4 * dx

        a = y(k-1)
        b = yp(k-1)
        c = ypp(k-1) / 2.0_dp

        aa = y(k) - c*dx2 - b*dx - a
        bb = (yp(k) - 2.0_dp*c*dx - b)*dx

        d = ( 4.0_dp*aa - bb) / dx3
        e = (-3.0_dp*aa + bb) / dx4

        int0(k) = e*dx5/5.0_dp + d*dx4/4.0_dp + c*dx3/3.0_dp + b*dx2/2.0_dp + a*dx
        err(k)  = (y(k-1)+y(k))*dx/2.0_dp
        err(k)  = ABS(err(k)-int0(k))

        dx = dx / 2.0_dp
        dx2 = dx * dx
        dx3 = dx2 * dx
        dx4 = dx2 * dx2
        x(k-1) = x(k-1) + dx
        y(k-1) = e*dx4 + d*dx3 + c*dx2 + b*dx + a
      END DO
    END IF
    !splitloc = MAXLOC(err,1)-1
    DEALLOCATE(yp)
    DEALLOCATE(ypp)
  END SUBROUTINE eval_bsinterr_test

  ! display (helper routine)
  SUBROUTINE disp_i8(i)
    INTEGER(kind=longint), INTENT(in)  :: i
    INTEGER                      :: k
    INTEGER                      :: is
    CHARACTER(len=BIT_SIZE(i))   :: c
    CHARACTER(len=20)            :: form
    is = BIT_SIZE(i)

    WRITE(form,'(a,i2,a)') '(',is,'i1)'
    WRITE(c,TRIM(form)) (IBITS(i,k,1), k=is-1,0,-1)
    CALL disp(c)
  END SUBROUTINE disp_i8

  SUBROUTINE disp_i81(i)
    INTEGER(kind=longint), DIMENSION(:), INTENT(in)  :: i

    INTEGER                                    :: k,ks
    INTEGER                                    :: is,s1
    CHARACTER(len=SIZE(i,1)*BIT_SIZE(i(1)))    :: c
    CHARACTER(len=20)                          :: form
    is = BIT_SIZE(i(1))

    s1 = SIZE(i,1)
    WRITE(form,'(a,i2,a)') '(',s1*is,'i1)'
    WRITE(c,TRIM(form)) ( (IBITS(i(ks),k,1),k=is-1,0,-1),ks=s1,1,-1 )
    CALL disp(c)
  END SUBROUTINE disp_i81

  SUBROUTINE disp_c(c)
    CHARACTER(len=*), INTENT(in) :: c
    PRINT *, c
  END SUBROUTINE disp_c


  SUBROUTINE plag_coeff(npoi,nder,x,xp,coef)
    !
    ! npoi - number of points (determines the order of Lagrange polynomial
    ! which is equal npoi-1)
    ! nder - number of derivatives computed 0 - function only, 1 - first
    ! derivative
    ! x - actual point where function and derivatives are evaluated
    ! xp(npoi) - array of points where function is known
    ! coef(0:nder,npoi) - weights for computation of function and derivatives,
    ! f=sum(fun(1:npoi)*coef(0,1:npoi) gives the function value
    ! df=sum(fun(1:npoi)*coef(1,1:npoi) gives the derivative value value
    !
    !
    INTEGER, INTENT(in)                                :: npoi,nder
    REAL(kind=dp), INTENT(in)                          :: x
    REAL(kind=dp), DIMENSION(npoi), INTENT(in)         :: xp
    REAL(kind=dp), DIMENSION(0:nder,npoi), INTENT(out) :: coef
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE           :: dummy
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE           :: fak_i
    !
    INTEGER                                            :: i,k,j,l,m
    REAL(kind=dp)                                      :: fac
    REAL(kind=dp)                                      :: j_sum,l_sum,m_sum,k_prod
    !
    DO i=1,npoi
      coef(0,i)=1.d0
      DO k=1,npoi
        IF(k.EQ.i) CYCLE
        coef(0,i)=coef(0,i)*(x-xp(k))/(xp(i)-xp(k))
      ENDDO
    ENDDO
    !
    IF(nder.EQ.0) RETURN
    !
    ALLOCATE(dummy(npoi))
    !
    DO i=1,npoi
      dummy=1.d0
      dummy(i)=0.d0
      DO k=1,npoi
        IF(k.EQ.i) CYCLE
        fac=(x-xp(k))/(xp(i)-xp(k))
        DO j=1,npoi
          IF(j.EQ.k) THEN
            dummy(j)=dummy(j)/(xp(i)-xp(k))
          ELSE
            dummy(j)=dummy(j)*fac
          ENDIF
        ENDDO
      ENDDO
      coef(1,i)=SUM(dummy)
    ENDDO
    !
    DEALLOCATE(dummy)
    !
    ! second derivative
    IF(nder.LE.1) RETURN
    !
    ALLOCATE(fak_i(npoi))
    do_i: DO i = 1,npoi
      fak_i = 0.0d0
      do_prep: DO k = 1,npoi
        IF (k .EQ. i) CYCLE
        fak_i(k) = (x-xp(k)) / (xp(i)-xp(k))
      END DO do_prep
      j_sum = 0.0d0
      do_j: DO j =1,npoi
        IF (j .EQ. i) CYCLE
        l_sum = 0.0d0
        do_l: DO l = 1,npoi
          IF (l .EQ. i .OR. l .EQ. j) CYCLE
          k_prod = 1.0d0
          do_k: DO k =1,npoi
            IF (k .EQ. i .OR. k .EQ. j .OR. k .EQ. l) CYCLE
            k_prod = k_prod * fak_i(k)
          END DO do_k
          l_sum = l_sum + k_prod / (xp(i)-xp(l))
        END DO do_l
        j_sum = j_sum + l_sum / (xp(i)-xp(j))
      END DO do_j
      coef(2,i)=j_sum
    END DO do_i
    DEALLOCATE(fak_i)

    ! third derivative
    IF(nder.LE.2) RETURN
    !
    ALLOCATE(fak_i(npoi))
    do_i3: DO i = 1,npoi
      fak_i = 0.0d0
      do_prep3: DO k = 1,npoi
        IF (k .EQ. i) CYCLE
        fak_i(k) = (x-xp(k)) / (xp(i)-xp(k))
      END DO do_prep3
      j_sum = 0.0d0
      do_j3: DO j =1,npoi
        IF (j .EQ. i) CYCLE
        l_sum = 0.0d0
        do_l3: DO l = 1,npoi
           IF (l .EQ. i .OR. l .EQ. j) CYCLE
           m_sum = 0.0d0
           do_m3: DO m = 1,npoi
             IF (m .EQ. i .OR. m .EQ. j .OR. m .EQ. l) CYCLE
             k_prod = 1.0d0
             do_k3: DO k =1,npoi
               IF (k .EQ. i .OR. k .EQ. j .OR. k .EQ. l .OR. k .EQ. m) CYCLE
               k_prod = k_prod * fak_i(k)
             END DO do_k3
             m_sum = m_sum + k_prod / (xp(i)-xp(m))
           END DO do_m3
           l_sum = l_sum + m_sum / (xp(i)-xp(l))
        END DO do_l3
        j_sum = j_sum + l_sum / (xp(i)-xp(j))
      END DO do_j3
      coef(3,i)=j_sum
    END DO do_i3
    DEALLOCATE(fak_i)

    RETURN
  END SUBROUTINE plag_coeff

  SUBROUTINE plag_stencil(ub,npoi,k,k1,k2,i1,i2)
    INTEGER, INTENT(in)  :: ub,npoi,k
    INTEGER, INTENT(out) :: k1,k2,i1,i2
    INTEGER :: kd

    k1 = k - npoi/2
    i1 = 1 - npoi/2
    IF (k1 .LT. 0) THEN
      kd = -k1
      k1 = k1 + kd
      i1 = i1 + kd
    ELSEIF (k1 .GT. ub-npoi+1) THEN
      kd = k1 - (ub-npoi+1)
      k1 = k1 - kd
      i1 = i1 - kd
    END IF
    k2 = k1 + npoi - 1
    i2 = i1 + npoi - 1

    RETURN
  END SUBROUTINE plag_stencil

  SUBROUTINE plag_test
    INTEGER, PARAMETER :: unitno = 9999
    INTEGER :: i
    REAL(kind=dp), PARAMETER :: pi = 3.141592653589793d0
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: x_ori,y0_ori,y1_ori,y2_ori,y3_ori
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: x,y0,y1,y2,y3
    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: coeff
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE  :: xlag, ylag
    REAL(kind=dp) :: xloc
    INTEGER :: lb,ub
    INTEGER :: k,k1,k2,i1,i2
    INTEGER, PARAMETER :: npoi = 6
    INTEGER, PARAMETER :: nder = 3
    INTEGER, PARAMETER :: ndata = 25

    ! basic data
    !CALL linspace(-pi,pi,ndata,x_ori)

    ALLOCATE(x_ori(0:ndata))
    x_ori(0) = -pi
    DO k = 1,ndata
      x_ori(k) = x_ori(0) + 2*pi*(DBLE(k)/DBLE(ndata))**1
    END DO

    lb = LBOUND(x_ori,1)
    ub = UBOUND(x_ori,1)
    ALLOCATE(y0_ori(lb:ub))
    ALLOCATE(y1_ori(lb:ub))
    ALLOCATE(y2_ori(lb:ub))
    ALLOCATE(y3_ori(lb:ub))
    y0_ori =  SIN(x_ori)
    y1_ori =  COS(x_ori)
    y2_ori = -SIN(x_ori)
    y3_ori = -COS(x_ori)
    ! plot basic data
    OPEN(file='plag_ori.dat',unit=unitno)
    DO i = lb,ub
      WRITE(unitno,*) x_ori(i),y0_ori(i),y1_ori(i),y2_ori(i),y3_ori(i)
    END DO
    CLOSE(unitno)

    ! lagrange coefficent
    ALLOCATE( coeff(0:nder,npoi) )


    ! interpolation data
    ALLOCATE(x(1:ub))
    ALLOCATE(y0(1:ub))
    ALLOCATE(y1(1:ub))
    ALLOCATE(y2(1:ub))
    ALLOCATE(y3(1:ub))

    DO k = 1,ub

      CALL plag_stencil(ub,npoi,k,k1,k2,i1,i2)
      ALLOCATE(xlag(i1:i2))
      ALLOCATE(ylag(i1:i2))
      xlag = x_ori(k1:k2)
      ylag = y0_ori(k1:k2)

      xloc = (xlag(0) + xlag(1)) / 2.0_dp
      CALL plagrange_coeff(npoi,nder,xloc,xlag,coeff)

      x(k) = xloc
      y0(k) = SUM(coeff(0,:)*ylag)
      y1(k) = SUM(coeff(1,:)*ylag)
      y2(k) = SUM(coeff(2,:)*ylag)
      y3(k) = SUM(coeff(3,:)*ylag)
      DEALLOCATE(xlag,ylag)
    END DO

    ! plot basic data
    OPEN(file='plag_int.dat',unit=unitno)
    DO i = 1,ub
      WRITE(unitno,*) x(i),y0(i),y1(i),y2(i),y3(i)
    END DO
    CLOSE(unitno)

    DEALLOCATE(x_ori)
    DEALLOCATE(y0_ori,y1_ori,y2_ori,y3_ori)
    DEALLOCATE(x)
    DEALLOCATE(y0,y1,y2,y3)
    DEALLOCATE(coeff)
    RETURN
  END SUBROUTINE plag_test

END MODULE binarysplit_int
