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

  ! linspace (helper routine)
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
                  x(k-1) .GT. bsfunc_p1(1) + fac*bsfunc_p2(1) ) THEN
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

MODULE binarysplit_mod
  ! This module handles all the splitting and joining stuff
  ! No changes should be necessary here
  USE binarysplit_int

  USE sparsevec_mod
  
  IMPLICIT NONE

  PRIVATE dp
  INTEGER, PARAMETER    :: dp = KIND(1.0d0)

  PRIVATE binarysplit_message
  INTEGER :: binarysplit_message = 0

  ! for forced splitting
  PUBLIC binarysplit_fsplitdepth
  INTEGER :: binarysplit_fsplitdepth = 9
  PRIVATE binarysplit_y0limfac
  REAL(kind=dp) :: binarysplit_y0limfac = 1.d-1

  PRIVATE longint
  INCLUDE 'longint.f90'

  PRIVATE binarysplit_limit,binarysplit_checklimit
  INTEGER               :: binarysplit_limit,binarysplit_checklimit
  PUBLIC binarysplit
  TYPE binarysplit
     INTEGER                                        :: n_ori
     INTEGER                                        :: n_split
     !INTEGER(kind=longint), dimension(:,:), allocatable   :: x_ori_bin
     TYPE(sparsevec), DIMENSION(:), ALLOCATABLE     :: x_ori_bin_sparse
     INTEGER,         DIMENSION(:),   ALLOCATABLE   :: x_ori_poi
     INTEGER,         DIMENSION(:),   ALLOCATABLE   :: x_poi
     INTEGER,         DIMENSION(:),   ALLOCATABLE   :: x_split
     INTEGER,         DIMENSION(:),   ALLOCATABLE   :: x_pos     
     REAL(kind=dp),   DIMENSION(:),   ALLOCATABLE   :: x
     REAL(kind=dp),   DIMENSION(:),   ALLOCATABLE   :: y
     REAL(kind=dp),   DIMENSION(:),   ALLOCATABLE   :: int
     REAL(kind=dp),   DIMENSION(:),   ALLOCATABLE   :: err     
  END TYPE binarysplit

  PUBLIC construct_binarysplit
  PRIVATE construct_binsplit, construct_binsplit_d1
  INTERFACE construct_binarysplit
     MODULE PROCEDURE construct_binsplit, construct_binsplit_d1
  END INTERFACE

  PUBLIC ASSIGNMENT(=)
  PRIVATE assign_binsplit
  INTERFACE ASSIGNMENT(=)
     MODULE PROCEDURE assign_binsplit
  END INTERFACE

  PUBLIC deconstruct_binarysplit
  PRIVATE deconstruct_binsplit
  INTERFACE deconstruct_binarysplit
     MODULE PROCEDURE deconstruct_binsplit
  END INTERFACE

  PUBLIC reallocate_binarysplit
  PRIVATE reallocate_binsplit
  INTERFACE reallocate_binarysplit
     MODULE PROCEDURE reallocate_binsplit
  END INTERFACE

  PUBLIC reposition_binarysplit
  PRIVATE reposition_binsplit, reposition_binsplit_bin
  INTERFACE reposition_binarysplit
     MODULE PROCEDURE reposition_binsplit, reposition_binsplit_bin
  END INTERFACE

  PUBLIC get_binarysplit
  PRIVATE get_binsplit_i1, get_binsplit_i18, get_binsplit_d1
  INTERFACE get_binarysplit
     MODULE PROCEDURE get_binsplit_i1, get_binsplit_i18, &
          get_binsplit_d1
  END INTERFACE

  PUBLIC split_binarysplit
  PRIVATE split_binsplit
  INTERFACE split_binarysplit
     MODULE PROCEDURE split_binsplit
  END INTERFACE

  PUBLIC multiple_binarysplit
  PRIVATE  multiple_binsplit
  INTERFACE multiple_binarysplit
     MODULE PROCEDURE multiple_binsplit
  END INTERFACE

  PUBLIC find_binarysplit
  PRIVATE find_binsplit
  INTERFACE find_binarysplit
     MODULE PROCEDURE find_binsplit,find_binsplit_check
  END INTERFACE

  PUBLIC dosplit_binarysplit
  PRIVATE dosplit_binsplit
  INTERFACE dosplit_binarysplit
     MODULE PROCEDURE dosplit_binsplit
  END INTERFACE

  PUBLIC join_binarysplit
  PRIVATE join_binsplit_v, join_binsplit_va
  INTERFACE join_binarysplit
     MODULE PROCEDURE join_binsplit_v, join_binsplit_va
  END INTERFACE

  PUBLIC compare_binarysplit
  PRIVATE compare_binsplit
  INTERFACE compare_binarysplit
     MODULE PROCEDURE compare_binsplit
  END INTERFACE

  PUBLIC printsummary_binarysplit
  PRIVATE printsummary_binsplit
  INTERFACE printsummary_binarysplit
     MODULE PROCEDURE printsummary_binsplit
  END INTERFACE

  PUBLIC printbin_binarysplit
  PRIVATE printbin_binsplit, printbin_binsplit_bin
  INTERFACE printbin_binarysplit
     MODULE PROCEDURE printbin_binsplit, printbin_binsplit_bin
  END INTERFACE

  PUBLIC plotfile_binarysplit
  PRIVATE plotfile_binsplit
  INTERFACE plotfile_binarysplit
     MODULE PROCEDURE plotfile_binsplit, plotfile_binsplit_ori, plotfile_binsplit_test
  END INTERFACE

  PUBLIC test_binarysplit
  PRIVATE test_binsplit
  INTERFACE test_binarysplit
     MODULE PROCEDURE test_binsplit
  END INTERFACE

CONTAINS

  ! construct
  SUBROUTINE construct_binsplit(x_beg,x_end,n_ori,xbs,n_split_max_o)
    REAL(kind=dp),               INTENT(in)       :: x_beg
    REAL(kind=dp),               INTENT(in)       :: x_end
    INTEGER,                     INTENT(in)       :: n_ori
    TYPE(binarysplit),           INTENT(inout)    :: xbs
    INTEGER,           OPTIONAL, INTENT(in)       :: n_split_max_o
    
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE      :: x
    INTEGER                                       :: n_split_max = 100

    IF (PRESENT(n_split_max_o)) n_split_max = n_split_max_o
    CALL linspace(x_beg,x_end,n_ori+1,x)
    CALL construct_binarysplit(x,xbs,n_split_max) 
    DEALLOCATE(x)
  END SUBROUTINE construct_binsplit

  SUBROUTINE construct_binsplit_d1(x,xbs,n_split_max_o)
    REAL(kind=dp),     DIMENSION(:), ALLOCATABLE, INTENT(in)        :: x !0:
    TYPE(binarysplit),               INTENT(inout)     :: xbs
    INTEGER,           OPTIONAL,     INTENT(in)        :: n_split_max_o
     
    INTEGER                                            :: k, n_ori
    INTEGER                                            :: n_split_max = 100

    IF (PRESENT(n_split_max_o)) n_split_max = n_split_max_o
    
    n_ori = SIZE(x,1)-1
    xbs%n_ori   = n_ori
    xbs%n_split = 0
    
    !IF (ALLOCATED(xbs%x_ori_bin)) DEALLOCATE(xbs%x_ori_bin)
    !ALLOCATE(xbs%x_ori_bin(0:0,0:n_ori))
    !xbs%x_ori_bin = 0
    !DO k = 0, n_ori
    !   xbs%x_ori_bin(0,k) = ibset(xbs%x_ori_bin(0,k),0)
    !END DO

    !**********************************************************
    ! Sparse x_ori_bin
    !**********************************************************
    IF (ALLOCATED(xbs%x_ori_bin_sparse)) DEALLOCATE(xbs%x_ori_bin_sparse)
    ALLOCATE(xbs%x_ori_bin_sparse(0:n_ori))
    DO k = 0, n_ori
       CALL xbs%x_ori_bin_sparse(k)%modify(0, IBSET(xbs%x_ori_bin_sparse(k)%get(0),0))
       xbs%x_ori_bin_sparse(k)%len = 0
    END DO
        
    IF (ALLOCATED(xbs%x_ori_poi)) DEALLOCATE(xbs%x_ori_poi)
    ALLOCATE(xbs%x_ori_poi(0:n_ori))
    xbs%x_ori_poi(0:n_ori) = (/ (k, k=0,n_ori) /)
    
    IF (ALLOCATED(xbs%x_poi)) DEALLOCATE(xbs%x_poi)
    ALLOCATE (xbs%x_poi(0:n_ori+n_split_max))
    xbs%x_poi(0:n_ori) = xbs%x_ori_poi(0:n_ori)
    xbs%x_poi(n_ori+1:n_ori+n_split_max) = 0
    
    IF (ALLOCATED(xbs%x_split)) DEALLOCATE(xbs%x_split)
    ALLOCATE (xbs%x_split(0:n_ori+n_split_max))
    xbs%x_split(0:n_ori) = 0
    
    IF (ALLOCATED(xbs%x_pos)) DEALLOCATE(xbs%x_pos)
    ALLOCATE (xbs%x_pos(0:n_ori+n_split_max))
    xbs%x_pos(0:n_ori) = 0
    
    IF (ALLOCATED(xbs%x)) DEALLOCATE(xbs%x)
    ALLOCATE (xbs%x(0:n_ori+n_split_max))
    xbs%x(0:n_ori) = x(0:n_ori)
    xbs%x(n_ori+1:n_ori+n_split_max) = 0.0_dp

    IF (ALLOCATED(xbs%y)) DEALLOCATE(xbs%y)
    ALLOCATE (xbs%y(0:n_ori+n_split_max))
    xbs%y = 0.0_dp
  
    IF (ALLOCATED(xbs%int)) DEALLOCATE(xbs%int)
    ALLOCATE (xbs%INT(0:n_ori+n_split_max))
    xbs%int = 0.0_dp
    
    IF (ALLOCATED(xbs%err)) DEALLOCATE(xbs%err)
    ALLOCATE (xbs%err(0:n_ori+n_split_max))
    xbs%err = 0.0_dp
    
  END SUBROUTINE construct_binsplit_d1

  ! assign
  SUBROUTINE assign_binsplit(xbs1,xbs2)
    TYPE(binarysplit), INTENT(inout)  :: xbs1
    TYPE(binarysplit), INTENT(in)     :: xbs2

    INTEGER                           :: n_ori,n_split,n_tot,s1, k

    ! This is not the full check but it works
    IF (ALLOCATED(xbs2%x) .AND. ALLOCATED(xbs2%x_ori_bin_sparse)) THEN
       n_ori = xbs2%n_ori
       n_split = xbs2%n_split
       n_tot = UBOUND(xbs2%x,1)
       !s1 = UBOUND(xbs2%x_ori_bin,1)
       s1 = xbs2%x_ori_bin_sparse(0)%len
       xbs1%n_ori = n_ori
       xbs1%n_split = n_split
       
       !IF (ALLOCATED(xbs1%x_ori_bin)) DEALLOCATE(xbs1%x_ori_bin)
       !ALLOCATE(xbs1%x_ori_bin(0:s1,0:n_ori))
       !xbs1%x_ori_bin = xbs2%x_ori_bin 

       IF (ALLOCATED(xbs1%x_ori_bin_sparse)) DEALLOCATE(xbs1%x_ori_bin_sparse)
       ALLOCATE(xbs1%x_ori_bin_sparse(0:n_ori))
       DO k = 0, n_ori
          CALL xbs1%x_ori_bin_sparse(k)%ASSIGN(xbs2%x_ori_bin_sparse(k))
       END DO
       !xbs1%x_ori_bin_sparse = xbs2%x_ori_bin_sparse
       
       IF (ALLOCATED(xbs1%x_ori_poi)) DEALLOCATE(xbs1%x_ori_poi)
       ALLOCATE(xbs1%x_ori_poi(0:n_ori))
       xbs1%x_ori_poi = xbs2%x_ori_poi 

       IF (ALLOCATED(xbs1%x_poi)) DEALLOCATE(xbs1%x_poi)
       ALLOCATE (xbs1%x_poi(0:n_tot))
       xbs1%x_poi = xbs2%x_poi
       
       IF (ALLOCATED(xbs1%x_split)) DEALLOCATE(xbs1%x_split)
       ALLOCATE (xbs1%x_split(0:n_tot))
       xbs1%x_split = xbs2%x_split
       
       IF (ALLOCATED(xbs1%x_pos)) DEALLOCATE(xbs1%x_pos)
       ALLOCATE (xbs1%x_pos(0:n_tot))
       xbs1%x_pos = xbs2%x_pos
       
       IF (ALLOCATED(xbs1%x)) DEALLOCATE(xbs1%x)
       ALLOCATE (xbs1%x(0:n_tot))
       xbs1%x =  xbs2%x
       
       IF (ALLOCATED(xbs1%y)) DEALLOCATE(xbs1%y)
       ALLOCATE (xbs1%y(0:n_tot))
       xbs1%y =  xbs2%y

       IF (ALLOCATED(xbs1%int)) DEALLOCATE(xbs1%int)
       ALLOCATE (xbs1%INT(0:n_tot))
       xbs1%int =  xbs2%int

       IF (ALLOCATED(xbs1%err)) DEALLOCATE(xbs1%err)
       ALLOCATE (xbs1%err(0:n_tot))
       xbs1%err =  xbs2%err    
    END IF
    
  END SUBROUTINE assign_binsplit



  ! deconstruct
  SUBROUTINE deconstruct_binsplit(xbs)
    TYPE(binarysplit), INTENT(inout)    :: xbs
    
    xbs%n_ori   = 0
    xbs%n_split = 0
    
    !IF (ALLOCATED(xbs%x_ori_bin)) DEALLOCATE(xbs%x_ori_bin)
    IF (ALLOCATED(xbs%x_ori_poi)) DEALLOCATE(xbs%x_ori_poi)
    IF (ALLOCATED(xbs%x_poi)) DEALLOCATE(xbs%x_poi)
    IF (ALLOCATED(xbs%x_split)) DEALLOCATE(xbs%x_split)
    IF (ALLOCATED(xbs%x_pos)) DEALLOCATE(xbs%x_pos)
    IF (ALLOCATED(xbs%x)) DEALLOCATE(xbs%x)
    IF (ALLOCATED(xbs%y)) DEALLOCATE(xbs%y)    
    IF (ALLOCATED(xbs%int)) DEALLOCATE(xbs%int)    
    IF (ALLOCATED(xbs%err)) DEALLOCATE(xbs%err)

    !**********************************************************
    ! Sparse
    !**********************************************************
    IF (ALLOCATED(xbs%x_ori_bin_sparse)) DEALLOCATE(xbs%x_ori_bin_sparse)
  END SUBROUTINE deconstruct_binsplit

  ! reallocate
  SUBROUTINE reallocate_binsplit(xbs,endall)
    TYPE(binarysplit),           INTENT(inout) :: xbs
    INTEGER,                     INTENT(in)    :: endall
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE   :: dvec
    INTEGER,       DIMENSION(:), ALLOCATABLE   :: ivec
    INTEGER                                    :: n_ori, n_split
    INTEGER                                    :: n_split_max
        
    n_ori   = xbs%n_ori
    n_split = xbs%n_split
    IF (endall .EQ. 1) THEN
       n_split_max = n_split
    ELSE
       n_split_max = n_split + 100
    END IF

    ALLOCATE(ivec(0:n_ori+n_split))
    ivec(0:n_ori+n_split) = xbs%x_poi(0:n_ori+n_split)
    DEALLOCATE(xbs%x_poi)
    ALLOCATE (xbs%x_poi(0:n_ori+n_split_max))
    xbs%x_poi(0:n_ori+n_split) = ivec(0:n_ori+n_split)

    ivec(0:n_ori+n_split) = xbs%x_split(0:n_ori+n_split)
    DEALLOCATE(xbs%x_split)
    ALLOCATE (xbs%x_split(0:n_ori+n_split_max))
    xbs%x_split(0:n_ori+n_split) = ivec(0:n_ori+n_split)
    
    ivec(0:n_ori+n_split) = xbs%x_pos(0:n_ori+n_split)
    DEALLOCATE(xbs%x_pos)
    ALLOCATE (xbs%x_pos(0:n_ori+n_split_max))
    xbs%x_pos(0:n_ori+n_split) = ivec(0:n_ori+n_split)
    DEALLOCATE(ivec)
     
    ALLOCATE(dvec(0:n_ori+n_split))
    dvec(0:n_ori+n_split) = xbs%x(0:n_ori+n_split)
    DEALLOCATE(xbs%x)
    ALLOCATE (xbs%x(0:n_ori+n_split_max))
    xbs%x(0:n_ori+n_split) = dvec(0:n_ori+n_split)

    dvec(0:n_ori+n_split) = xbs%y(0:n_ori+n_split)
    DEALLOCATE(xbs%y)
    ALLOCATE (xbs%y(0:n_ori+n_split_max))
    xbs%y(0:n_ori+n_split) = dvec(0:n_ori+n_split)
     
    dvec(0:n_ori+n_split) = xbs%INT(0:n_ori+n_split)
    DEALLOCATE(xbs%int)
    ALLOCATE (xbs%INT(0:n_ori+n_split_max))
    xbs%INT(0:n_ori+n_split) = dvec(0:n_ori+n_split)
     
    dvec(0:n_ori+n_split) = xbs%err(0:n_ori+n_split)
    DEALLOCATE(xbs%err)
    ALLOCATE (xbs%err(0:n_ori+n_split_max))
    xbs%err(0:n_ori+n_split) = dvec(0:n_ori+n_split)
    DEALLOCATE(dvec)
     
  END SUBROUTINE reallocate_binsplit
    
  ! reposition
  SUBROUTINE reposition_binsplit(xbs)   
    TYPE(binarysplit),               INTENT(inout) :: xbs
    
    INTEGER                                        :: n_ori,maxdim,s1,is
    INTEGER                                        :: k,n,i
    INTEGER                                        :: pos_n,k_n,s1_n
    INTEGER(kind=longint), DIMENSION(:,:), ALLOCATABLE   :: v

    INTEGER(HID_T) :: h5id
    CHARACTER(len=1024) :: h5_ds_name

    INTEGER(kind=longint) :: zero
    INTEGER :: sparselen
    TYPE(sparsevec), DIMENSION(:), ALLOCATABLE :: v_sparse

!!$    type(sparsevec) :: testvec
!!$    zero = 0
!!$    call testvec%modify(0,zero)
!!$    
!!$    zero = 10
!!$    call testvec%modify(1,zero)
!!$    
!!$    zero = 50
!!$    call testvec%modify(5,zero)
!!$    
!!$    zero = 20
!!$    call testvec%modify(2,zero)
!!$    
!!$    zero = 11
!!$    call testvec%modify(1,zero)
!!$    
!!$    zero = 55
!!$    call testvec%modify(5,zero)
!!$    
!!$    zero = -1
!!$    call testvec%modify(0,zero)
!!$    
!!$    write (*,*) "*** Testvec ***"
!!$    write (*,*) testvec%idxvec
!!$    write (*,*) testvec%values
!!$    write (*,*) "*** *** ***"
!!$    !stop

    
    n_ori  = xbs%n_ori
    maxdim = n_ori + xbs%n_split
    xbs%x_pos(0:maxdim) = 2*xbs%x_pos(0:maxdim)
    
    !s1 = UBOUND(xbs%x_ori_bin,1)
    s1 = xbs%x_ori_bin_sparse(0)%len
    s1_n = 2**(NINT(LOG(DBLE(s1+1))/LOG(2.0_dp))+1)-1

    !ALLOCATE(v(0:s1,0:n_ori))
    !v(0:s1,0:n_ori) = xbs%x_ori_bin(0:s1,0:n_ori)
    !DEALLOCATE(xbs%x_ori_bin)
    !ALLOCATE(xbs%x_ori_bin(0:s1_n,0:n_ori))
    !xbs%x_ori_bin = 0
    !is = BIT_size(v(0,0))
    is = bit_SIZE(xbs%x_ori_bin_sparse(0)%get(0))

    !**********************************************************
    ! Original version
    !**********************************************************
    !do n = 0, n_ori
    !   do k = s1,0,-1
    !      do i = is-1,0,-1
    !         if (btest(v(k,n),i)) then
    !            pos_n = 2*(k*is+i)
    !            k_n   = floor(dble(pos_n)/dble(is))
    !            pos_n = mod(pos_n,is)
    !            xbs%x_ori_bin(k_n,n) = ibset(xbs%x_ori_bin(k_n,n),pos_n)
    !         end if
    !      end do
    !   end do
    !end do
    !deallocate(v)

    !**********************************************************
    ! Sparse version
    !**********************************************************
    ALLOCATE(v_sparse(0:n_ori))
    v_sparse%len = s1
    xbs%x_ori_bin_sparse%len = s1_n
    DO n = 0, n_ori
       CALL v_sparse(n)%ASSIGN(xbs%x_ori_bin_sparse(n))
       CALL xbs%x_ori_bin_sparse(n)%clear()
       
       !do k = s1,0,-1
       DO k = 1, v_sparse(n)%len_sparse
          DO i = is-1,0,-1
             !if (btest(v_sparse(n)%get(k),i)) then
             IF (BTEST(v_sparse(n)%get(v_sparse(n)%idxvec(k)),i)) THEN
                !pos_n = 2*(k*is+i)
                pos_n = 2*(v_sparse(n)%idxvec(k)*is+i)
                k_n   = FLOOR(DBLE(pos_n)/DBLE(is))
                pos_n = MOD(pos_n,is)
                
                !call xbs%x_ori_bin_sparse(n)%modify(k_n, ibset(xbs%x_ori_bin_sparse(n)%get(k_n), pos_n))
                CALL xbs%x_ori_bin_sparse(n)%IBSET(k_n, pos_n)
             END IF
          END DO
       END DO
    END DO
    DEALLOCATE(v_sparse)
    
!!$    if (.false.) then
!!$       call h5_open_rw('binarysplit_reposition.h5', h5id)
!!$       write (h5_ds_name, '(A,I0)') "x_ori_bin_", INT(log(s1_n+1d0)/log(2d0))
!!$       call h5_add(h5id, h5_ds_name, xbs%x_ori_bin, lbound(xbs%x_ori_bin), ubound(xbs%x_ori_bin))
!!$       do k=0,n_ori
!!$          write (h5_ds_name, '(A,I0,A,I0)') "x_ori_bin_sparse_idxvec", int(log(s1_n+1d0)/log(2d0)),"_",k
!!$          call h5_add(h5id, h5_ds_name, xbs%x_ori_bin_sparse(k)%idxvec, lbound(xbs%x_ori_bin_sparse(k)%idxvec),&
!!$               ubound(xbs%x_ori_bin_sparse(k)%idxvec))
!!$          write (h5_ds_name, '(A,I0,A,I0)') "x_ori_bin_sparse_values", int(log(s1_n+1d0)/log(2d0)),"_",k
!!$          call h5_add(h5id, h5_ds_name, xbs%x_ori_bin_sparse(k)%values, lbound(xbs%x_ori_bin_sparse(k)%values),&
!!$               ubound(xbs%x_ori_bin_sparse(k)%values))
!!$
!!$       end do
!!$       call h5_close(h5id)
!!$    end if
    
  END SUBROUTINE reposition_binsplit

  SUBROUTINE reposition_binsplit_bin(bin_sparse)
    !INTEGER(kind=longint), DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: bin !0:
    TYPE(sparsevec), DIMENSION(:), ALLOCATABLE, INTENT(inout) :: bin_sparse
    
    INTEGER                                        :: n_ori,s1,is
    INTEGER                                        :: k,n,i
    INTEGER                                        :: pos_n,k_n,s1_n
    !INTEGER(kind=longint), DIMENSION(:,:), ALLOCATABLE   :: v
    TYPE(sparsevec), DIMENSION(:), ALLOCATABLE :: v_sparse
    
    !n_ori  = UBOUND(bin,2)
    !s1 = UBOUND(bin,1)
    n_ori = UBOUND(bin_sparse, 1)
    s1 = bin_sparse(0)%len
    s1_n = 2**(NINT(LOG(DBLE(s1+1))/LOG(2.0_dp))+1)-1
    
    !ALLOCATE(v(0:s1,0:n_ori))
    !v(0:s1,0:n_ori) = bin(0:s1,0:n_ori)
    !DEALLOCATE(bin)
    !ALLOCATE(bin(0:s1_n,0:n_ori))
    !bin = 0
    !is = BIT_SIZE(v(0,0))

    ALLOCATE(v_sparse(0:n_ori))
    v_sparse%len = s1
    is = bit_SIZE(v_sparse(0)%get(0))

    bin_sparse%len = s1_n
    DO n = 0, n_ori
       CALL v_sparse(n)%ASSIGN(bin_sparse(n))
       CALL bin_sparse(n)%clear()
       !do k = s1,0,-1
       DO k = 1, bin_sparse(n)%len_sparse
          DO i = is-1,0,-1
             !IF (BTEST(v(k,n),i)) THEN
             !if (btest(v_sparse(n)%get(k),i)) then
             IF (BTEST(v_sparse(n)%get(v_sparse(n)%idxvec(k)),i)) THEN
                !pos_n = 2*(k*is+i)
                pos_n = 2*(v_sparse(n)%idxvec(k)*is+i)
                k_n   = FLOOR(DBLE(pos_n)/DBLE(is))
                pos_n = MOD(pos_n,is)
                !bin(k_n,n) = ibset(bin(k_n,n),pos_n)
                !call bin_sparse(n)%modify(k_n, ibset(bin_sparse(n)%get(k_n), pos_n))
                CALL bin_sparse(n)%IBSET(k_n, pos_n)
             END IF
          END DO
       END DO
    END DO
    !deallocate(v)
    DEALLOCATE(v_sparse)

  END SUBROUTINE reposition_binsplit_bin

  ! extract
  SUBROUTINE get_binsplit_d1(xbs,v,c)
    TYPE(binarysplit),                        INTENT(in)    :: xbs
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(inout) :: v
    CHARACTER(len=*),                         INTENT(in)    :: c
  
    IF (ALLOCATED(v)) DEALLOCATE(v)
    ALLOCATE(v(0:xbs%n_ori+xbs%n_split))
  
    IF (c .EQ. 'x') THEN
       v(0:xbs%n_ori+xbs%n_split) = xbs%x(0:xbs%n_ori+xbs%n_split)
    ELSEIF (c .EQ. 'y') THEN
       v(0:xbs%n_ori+xbs%n_split) = xbs%y(0:xbs%n_ori+xbs%n_split)
    ELSEIF (c .EQ. 'int') THEN
       v(0:xbs%n_ori+xbs%n_split) = xbs%INT(0:xbs%n_ori+xbs%n_split)
    ELSEIF (c .EQ. 'err') THEN
       v(0:xbs%n_ori+xbs%n_split) = xbs%err(0:xbs%n_ori+xbs%n_split)
    ELSE
       IF (binarysplit_message .GT. 0) PRINT *, c,' can not be extracted'
    END IF
  END SUBROUTINE get_binsplit_d1

  SUBROUTINE get_binsplit_i1(xbs,v,c)
    TYPE(binarysplit),                        INTENT(in)    :: xbs
    INTEGER,       DIMENSION(:), ALLOCATABLE, INTENT(inout) :: v
    CHARACTER(len=*),                         INTENT(in)    :: c
  
    IF (ALLOCATED(v)) DEALLOCATE(v)
  
    IF (c .EQ. 'x_ori_poi') THEN
       ALLOCATE(v(0:xbs%n_ori))
       v(0:xbs%n_ori) = xbs%x_ori_poi(0:xbs%n_ori)
    ELSEIF (c .EQ. 'x_poi') THEN
       ALLOCATE(v(0:xbs%n_ori+xbs%n_split))
       v(0:xbs%n_ori+xbs%n_split) = xbs%x_poi(0:xbs%n_ori+xbs%n_split)
    ELSEIF (c .EQ. 'x_split') THEN
       ALLOCATE(v(0:xbs%n_ori+xbs%n_split))
       v(0:xbs%n_ori+xbs%n_split) = xbs%x_split(0:xbs%n_ori+xbs%n_split)
    ELSEIF (c .EQ. 'x_pos') THEN
       ALLOCATE(v(0:xbs%n_ori+xbs%n_split))
       v(0:xbs%n_ori+xbs%n_split) = xbs%x_pos(0:xbs%n_ori+xbs%n_split)
    ELSE
       IF (binarysplit_message .GT. 0) PRINT *, c,' can not be extracted'
    END IF
  END SUBROUTINE get_binsplit_i1
  
  SUBROUTINE get_binsplit_i18(xbs,v_sparse,c)
    TYPE(binarysplit),                            INTENT(in)    :: xbs
    !INTEGER(kind=longint), DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: v
    TYPE(sparsevec), DIMENSION(:), ALLOCATABLE, INTENT(inout)   :: v_sparse
    CHARACTER(len=*),                             INTENT(in)    :: c

    INTEGER                                                     :: s1
    INTEGER :: k
    
    !s1 = UBOUND(xbs%x_ori_bin,1)
    s1 = xbs%x_ori_bin_sparse(0)%len
    !IF (ALLOCATED(v)) DEALLOCATE(v)
    IF (ALLOCATED(v_sparse)) DEALLOCATE(v_sparse)
    IF (c .EQ. 'x_ori_bin') THEN
       !ALLOCATE(v(0:s1,0:xbs%n_ori))
       !v(:,0:xbs%n_ori) = xbs%x_ori_bin(:,0:xbs%n_ori)

       ALLOCATE(v_sparse(0:xbs%n_ori))
       DO k = 0, xbs%n_ori
          CALL v_sparse(k)%ASSIGN(xbs%x_ori_bin_sparse(k))
       END DO
       
    ELSE
       IF (binarysplit_message .GT. 0) PRINT *, c,' can not be extracted'
    END IF
  END SUBROUTINE get_binsplit_i18

  ! split
  SUBROUTINE split_binsplit(xbs,splitloc)
    TYPE(binarysplit), INTENT(inout)              :: xbs
    INTEGER,           INTENT(in)                 :: splitloc

    INTEGER                                       :: maxind,k,k_n
    INTEGER                                       :: n_ori,n_split,s1
    INTEGER                                       :: nb2,nb,is
    INTEGER                                       :: ival,split,pos,pos_n
    !INTEGER(kind=longint),   DIMENSION(:), ALLOCATABLE  :: bin
    TYPE(sparsevec), DIMENSION(:), ALLOCATABLE    :: bin_sparse
    INTEGER                                       :: nb1 
    REAL(kind=dp)                                 :: dval,dist

    n_ori       = xbs%n_ori
    n_split     = xbs%n_split
    maxind      = n_ori + n_split
    !s1  = UBOUND(xbs%x_ori_bin,1)
    !is  = BIT_SIZE(xbs%x_ori_bin(0,0))
    s1   = xbs%x_ori_bin_sparse(0)%len
    is   = bit_SIZE(xbs%x_ori_bin_sparse(0)%get(0))
    
    nb1 = NINT(LOG(DBLE(is))/LOG(2.0_dp))
    nb2 = NINT(LOG(DBLE(s1+1))/LOG(2.0_dp)) 
    nb  = nb1 + nb2

    ! new value and new split level
    dval = (xbs%x(splitloc-1)+xbs%x(splitloc))/2.0_dp
    dist = ABS(dval - xbs%x(splitloc-1))
    split = xbs%x_split(splitloc)+1
    ! limits
    IF (binarysplit_checklimit .EQ. 1) THEN
       IF (.NOT. eval_bslimit(dist,split,maxind+1)) THEN
          binarysplit_limit = 1
          RETURN
       END IF
    END IF

    IF (maxind+1 .GT. UBOUND(xbs%x,1)) THEN
       CALL reallocate_binarysplit(xbs,0)
       IF (binarysplit_message .GT. 0) PRINT *, 'Message: reallocate binarysplit'
    END IF
    xbs%n_split = n_split + 1
    
    IF (split .GT. nb) THEN
       CALL reposition_binarysplit(xbs)
       !s1  = UBOUND(xbs%x_ori_bin,1)
       s1   = xbs%x_ori_bin_sparse(0)%len
       nb2 = NINT(LOG(DBLE(s1+1))/LOG(2.0_dp)) 
       nb  = nb1 + nb2
       IF (binarysplit_message .GT. 0) &
            PRINT *, 'Message: double x_ori_bin, nb: ', nb
    END IF
    
    ! x_ori_poi
    DO k = xbs%x_poi(splitloc), n_ori
       xbs%x_ori_poi(k) = xbs%x_ori_poi(k) + 1
    END DO
    ! x_poi
    xbs%x_poi(splitloc+1:maxind+1) = xbs%x_poi(splitloc:maxind)
    ! x_split
    xbs%x_split(splitloc+1:maxind+1) = xbs%x_split(splitloc:maxind)
    xbs%x_split(splitloc:splitloc+1) = split
    ! x_ori_bin, x_pos
    !ALLOCATE(bin(0:s1))
    xbs%x_pos(splitloc+1:maxind+1) = xbs%x_pos(splitloc:maxind)
    pos = xbs%x_pos(splitloc)
    !bin(:) = xbs%x_ori_bin(:,xbs%x_poi(splitloc))

    pos_n = pos+2**(nb-split)
    xbs%x_pos(splitloc) = pos_n
    k_n   = FLOOR(DBLE(pos_n)/DBLE(is))
    pos_n = MOD(pos_n,is)
       
    !IF (BTEST(bin(k_n),pos_n)) THEN
    !   IF (binarysplit_message .GT. 0) PRINT *, 'Warning x_ori_bin already set'
    !END IF
    !bin(k_n) = IBSET(bin(k_n),pos_n)
    !xbs%x_ori_bin(:,xbs%x_poi(splitloc)) = bin(:)
    !DEALLOCATE(bin)
    !write (*,*) "split_binsplit: E"

    !**********************************************************
    ! Sparse version
    !**********************************************************
    !call xbs%x_ori_bin_sparse(xbs%x_poi(splitloc))%modify(k_n, &
    !     ibset(xbs%x_ori_bin_sparse(xbs%x_poi(splitloc))%get(k_n), pos_n))
    IF (BTEST(xbs%x_ori_bin_sparse(xbs%x_poi(splitloc))%get(k_n),pos_n)) THEN
       IF (binarysplit_message .GT. 0) PRINT *, 'Warning x_ori_bin already set'
    END IF
    CALL xbs%x_ori_bin_sparse(xbs%x_poi(splitloc))%IBSET(k_n, pos_n)
    
    ! values
    xbs%x(splitloc+1:maxind+1) = xbs%x(splitloc:maxind)
    xbs%x(splitloc) = dval
    ! func (value has to be computed externally)
    xbs%y(splitloc+1:maxind+1) = xbs%y(splitloc:maxind)
    xbs%y(splitloc) = 0.0_dp
    ! int and err (value has to be computed externally)
    xbs%INT(splitloc+1:maxind+1) = xbs%INT(splitloc:maxind)
    xbs%INT(splitloc) = 0.0_dp
    xbs%err(splitloc+1:maxind+1) = xbs%err(splitloc:maxind)
    xbs%err(splitloc) = 0.0_dp

  END SUBROUTINE split_binsplit

!!$  SUBROUTINE split_binsplit_v(v,splitloc,maxind)
!!$    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(inout) :: v
!!$    INTEGER,                                   INTENT(in)    :: splitloc
!!$    INTEGER,                                   INTENT(in)    :: maxind
!!$
!!$    INTEGER                                                  :: ub
!!$    REAL(kind=dp), DIMENSION(:),  ALLOCATABLE                :: vh
!!$    
!!$    IF (maxind+1 .GT. UBOUND(v,1)) THEN
!!$       ub = UBOUND(v,1)
!!$       ALLOCATE(vh(0:ub+100))
!!$       vh(0:ub) = v
!!$       DEALLOCATE(v)
!!$       ALLOCATE(v(0:ub+100))
!!$       v = vh
!!$       DEALLOCATE(vh)
!!$    END IF
!!$    v(splitloc+1:maxind+1) = v(splitloc:maxind)
!!$    v(splitloc) = 0.0_dp
!!$  END SUBROUTINE split_binsplit_v

  ! do the splits
  SUBROUTINE find_binsplit(xbs)
    TYPE(binarysplit),                         INTENT(inout)  :: xbs

    INTEGER :: splitloc,maxind

    maxind = xbs%n_ori + xbs%n_split
    CALL eval_bsfunc(xbs%x(0:maxind),xbs%y(0:maxind))
    CALL eval_bsinterr(xbs%x(0:maxind),xbs%y(0:maxind),   &
         xbs%INT(0:maxind),xbs%err(0:maxind),splitloc)
    binarysplit_checklimit = 1
    binarysplit_limit = 0
    DO WHILE (eval_bslimit(xbs%err(0:maxind)))
       CALL split_binarysplit(xbs,splitloc)
       !write (*,*) "Iteration in find_binsplit() after split_binarysplit()"
       IF (binarysplit_limit .EQ. 1) THEN
          binarysplit_limit = 0
          EXIT
       END IF
       maxind = xbs%n_ori + xbs%n_split
       CALL eval_bsfunc(xbs%x(splitloc),xbs%y(splitloc))
       !write (*,*) "Iteration in find_binsplit() after eval_bsfunc()"

       CALL eval_bsinterr(xbs%x(0:maxind),xbs%y(0:maxind),   &
            xbs%INT(0:maxind),xbs%err(0:maxind),splitloc)
    END DO
    CALL reallocate_binarysplit(xbs,1)
    
  END SUBROUTINE find_binsplit

  ! do the splits and make checks
  SUBROUTINE find_binsplit_check(xbs,x0)
    TYPE(binarysplit),                         INTENT(inout)  :: xbs
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE,  INTENT(in)     :: x0
    
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE                  :: x,y
    REAL(kind=dp)                                             :: x0_a,y0_a
    REAL(kind=dp)                                             :: y0_lo,y0_up
    INTEGER                                                   :: splitloc,maxind
    INTEGER                                                   :: int_ref, int_act
    INTEGER                                                   :: ready
    INTEGER                                                   :: counter
    INTEGER                                                   :: ix0,i_up,i_lo,i_mi
    INTEGER                                                   :: split_c
    
    int_ref = SIZE(x0)
    ready = 0
    counter = 0
    outer: DO WHILE (ready .EQ. 0) 
       counter = counter + 1
       CALL find_binarysplit(xbs)
       maxind = xbs%n_ori + xbs%n_split
       int_act = NINT(SUM(xbs%int,1))
       IF (int_act .EQ. int_ref) THEN
          ready = 1
       ELSE IF (counter .GT. binarysplit_fsplitdepth) THEN
          ready = 1
       ELSE ! improve
          IF (binarysplit_message .GT. 0) &
               PRINT *, 'trying to improve:', counter
          split_c = 0
          DO ix0 = 1, int_ref
             CALL get_binarysplit(xbs,x,'x')
             CALL get_binarysplit(xbs,y,'y')
             x0_a = x0(ix0)
             i_up = xbs%n_ori + xbs%n_split
             i_lo = 0
             DO WHILE ((i_up-i_lo) .GT. 1)
                i_mi = (i_up+i_lo)/2
                IF (x0_a .LT. x(i_mi)) THEN
                   i_up = i_mi
                ELSE IF (x0_a .GT. x(i_mi)) THEN
                   i_lo = i_mi
                ELSE
                   i_up = i_mi
                   i_lo = i_mi -1
                END IF
             END DO
             CALL eval_bsfunc(x0_a,y0_a)
             y0_lo = y(i_lo)
             y0_up = y(i_up)
             
             y0_a = y0_a * binarysplit_y0limfac
             !PRINT *, i_lo,i_up,x(i_lo),x0_a,x(i_up)
             !PRINT *, y0_lo,y0_a,y0_up
             IF (y0_lo .LT. y0_a .OR. y0_up .LT. y0_a) THEN
                ! do a forced split here
                IF (binarysplit_message .GT. 0) PRINT *, 'forced split'
                binarysplit_checklimit = 1 ! communication with the splitter
                binarysplit_limit = 0
                CALL split_binarysplit(xbs,i_up)
                IF (binarysplit_limit .EQ. 1) THEN
                   binarysplit_limit = 0
                   EXIT outer
                END IF
                maxind = xbs%n_ori + xbs%n_split                
                CALL eval_bsfunc(xbs%x(i_up),xbs%y(i_up))
                split_c = split_c + 1
             END IF
             !PRINT *, 'split_counter ',split_c
             !PRINT *, 'maxind        ',maxind
             !PAUSE
          END DO
          IF (split_c .EQ. 0) ready = 1
       END IF
    END DO outer
    IF (ALLOCATED(x)) DEALLOCATE(x)
    IF (ALLOCATED(y)) DEALLOCATE(y)
    CALL reallocate_binarysplit(xbs,1)
  END SUBROUTINE find_binsplit_check

  ! do multiple splits with n (optional) 
  SUBROUTINE multiple_binsplit(xbs,x0,s,n_opt)
    TYPE(binarysplit),                        INTENT(inout)  :: xbs
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(in)     :: x0
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(in)     :: s
    INTEGER, OPTIONAL, INTENT(in) :: n_opt
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: s1

    !REAL(kind=dp) :: e1, e2, lambda
    !REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: d_hlp
    !REAL(kind=dp) :: d

    INTEGER :: i,n
    INTEGER :: save_binarysplit_fsplitdepth
    REAL(kind=dp) :: save_bsfunc_local_err,save_bsfunc_base_distance
    !integer :: maxind

    save_bsfunc_base_distance = bsfunc_base_distance
    save_binarysplit_fsplitdepth = binarysplit_fsplitdepth
    save_bsfunc_local_err = bsfunc_local_err

    !maxind = xbs%n_ori + xbs%n_split
    !print *, xbs%x(0:maxind)
    
    IF (PRESENT(n_opt)) THEN
       n = n_opt
    ELSE
       n = MAX(2,INT( LOG( 1.0d0 * bsfunc_base_distance/(s(1)*SQRT(-2*LOG(bsfunc_mult_constant))) ) / &
            LOG(bsfunc_sigma_multiplier) ) + 2) 
    END IF
    !print *, 'n = ',n
    !print *, 'bsfunc_base_distance ',bsfunc_base_distance
    !print *, 's ',s(1)
    !print *, 'bsfunc_mult_constant ',bsfunc_mult_constant
    !print *, 'bsfunc_sigma_multiplier ',bsfunc_sigma_multiplier

    bsfunc_modelfunc = 1
    binarysplit_fsplitdepth = AINT( LOG( bsfunc_base_distance/(s(1) ) / LOG(2.0d0) ) ) + 1
    !print *, 'binarysplit_fsplitdepth ',binarysplit_fsplitdepth
    CALL construct_bsfunc(x0,s)
    CALL find_binarysplit(xbs,x0)
    binarysplit_fsplitdepth = 0
    !CALL printsummary_binarysplit(xbs)

    !e2 = bsfunc_local_err
    !e1 = e2 / 4.0d0
    !lambda = 0.25d0
    
    !bsfunc_modelfunc = 2
    !CALL construct_bsfunc(x0,s)
    !CALL find_binarysplit(xbs,x0)
    !CALL printsummary_binarysplit(xbs)

    ALLOCATE(s1(LBOUND(s,1):UBOUND(s,1)))
    s1 = s * bsfunc_sigma_multiplier**n
    DO i = 1,n
       s1 = s1 / bsfunc_sigma_multiplier
       !bsfunc_local_err = e1 + (e2-e1)*(1.0d0 - exp(-lambda*dble(i-1)) + exp(-lambda*dble(n-1)))
       !print *, 'n,i,s1 = ',n,i,s1,bsfunc_local_err
       CALL construct_bsfunc(x0,s1)
       CALL find_binarysplit(xbs,x0)
    END DO
    !CALL printsummary_binarysplit(xbs)

    bsfunc_base_distance = save_bsfunc_base_distance
    binarysplit_fsplitdepth = save_binarysplit_fsplitdepth
    bsfunc_local_err = save_bsfunc_local_err
    DEALLOCATE(s1)
    RETURN
  END SUBROUTINE multiple_binsplit



 ! split at forced positions
  SUBROUTINE dosplit_binsplit(xbs1,xbs2,bin_sparse)
    TYPE(binarysplit),                           INTENT(inout) :: xbs1
    TYPE(binarysplit),                           INTENT(in)    :: xbs2
    !INTEGER(kind=longint),DIMENSION(:,:),ALLOCATABLE,INTENT(inout)   :: bin !0:
    TYPE(sparsevec), DIMENSION(:), ALLOCATABLE :: bin_sparse
    
    INTEGER :: n_ori,s1,is,nb1,nb2,nb,n_split
    INTEGER :: splitloc,x_ori_poi
    INTEGER :: n,sk,ik,k
    INTEGER :: pos_s,pos_d,pos,pos_n,k_n
    INTEGER :: count,p_c,p_nc,k_nc
    INTEGER :: ori_sw
    INTEGER :: counter,maxind

    REAL    :: start, finish
    INTEGER :: n_last, k_n_last, k_nc_last
    INTEGER(kind=longint) :: bin_sparse_last, xbs1_x_ori_bin_sparse

    binarysplit_checklimit = 0
    xbs1 = xbs2
    DO k = 0, xbs1%n_ori
       CALL xbs1%x_ori_bin_sparse(k)%ASSIGN(xbs2%x_ori_bin_sparse(k))
    END DO
    n_ori = xbs1%n_ori
    counter = 0

    CALL cpu_TIME(start)
    forever: DO
       !s1  = UBOUND(bin,1)
       !is  = BIT_SIZE(bin(0,0))

       s1 = bin_sparse(0)%len
       is = bit_SIZE(bin_sparse(0)%get(0))
       nb1 = NINT(LOG(DBLE(is))/LOG(2.0_dp))
       nb2 = NINT(LOG(DBLE(s1+1))/LOG(2.0_dp)) 
       nb  = nb1 + nb2

       n_last = -1
       k_n    = -1
       ori: DO n = 0, n_ori
          ori_sw = 0
          DO sk = nb-1, 0, -1
             pos_s = 2**sk
             pos_d = pos_s * 2
             ik = -1
             DO
                ik = ik + 1
                pos = pos_s + ik * pos_d
                IF (pos .GT. is*(s1+1)-1) EXIT
                pos_n = MOD(pos,is)
                k_n   = FLOOR(DBLE(pos)/DBLE(is))

                IF ((n .NE. n_last) .OR. (k_n .NE. k_n_last)) THEN
                   n_last = n
                   k_n_last = k_n
                   bin_sparse_last = bin_sparse(n)%get(k_n)
                ELSE
                   !write (*,*) "Reading binsplit from cache."
                END IF
                  
                IF (BTEST(bin_sparse_last,pos_n)) THEN
                !IF (BTEST(bin_sparse(n)%get(k_n),pos_n)) THEN
                   counter = counter + 1
                   !bin(k_n,n) = IBCLR(bin(k_n,n),pos_n)
                   !call bin_sparse(n)%modify(k_n, ibclr(bin_sparse(n)%get(k_n),pos_n))
                   CALL bin_sparse(n)%IBCLR(k_n, pos_n)

                   x_ori_poi = xbs1%x_ori_poi(n)
                   count = 0
                   k_nc_last = -1
                   DO p_c = pos-1,1,-1
                      p_nc = MOD(p_c,is)
                      k_nc   = FLOOR(DBLE(p_c)/DBLE(is))
                      !if (btest(xbs1%x_ori_bin(k_nc,n),p_nc)) count = count + 1

                      !**********************************************************
                      ! Sparse version
                      !**********************************************************
                      IF (k_nc .NE. k_nc_last) THEN
                         xbs1_x_ori_bin_sparse = xbs1%x_ori_bin_sparse(n)%get(k_nc)
                         k_nc_last = k_nc
                      END IF
                      
                      IF (BTEST(xbs1_x_ori_bin_sparse,p_nc)) count = count + 1 
                      !write (*,*) "Loop 2 in do split", k_nc
                   END DO
                   splitloc = x_ori_poi - count 
                   CALL split_binarysplit(xbs1,splitloc)
                   maxind = xbs1%n_ori + xbs1%n_split
                   CALL eval_bsfitsplit(                              &
                        xbs1%x(0:maxind),xbs1%y(0:maxind),            &
                        xbs1%INT(0:maxind),xbs1%err(0:maxind),        &
                        splitloc)
                   !IF (UBOUND(xbs1%x_ori_bin,1) .GT. s1) THEN 
                   IF (xbs1%x_ori_bin_sparse(0)%len .GT. s1) THEN 
                      IF (binarysplit_message .GT. 0) &
                           PRINT *, 'Message: I leave loop'
                      ori_sw = 1
                      EXIT ori
                   END IF
                END IF
             END DO
          END DO
       END DO ori
       IF (ori_sw .EQ. 1) THEN
          CALL reposition_binarysplit(bin_sparse)
       ELSE
          EXIT forever
       END IF
    END DO forever

    CALL reallocate_binarysplit(xbs1,1)
    CALL cpu_TIME(finish)
    PRINT '("Time in dosplit_binsplit = ",f12.3," seconds.")',finish-start
    
  END SUBROUTINE dosplit_binsplit

  ! join 
  SUBROUTINE join_binsplit_va(xbs1,xbs2)
    TYPE(binarysplit),                           INTENT(inout) :: xbs1
    TYPE(binarysplit),                           INTENT(in)    :: xbs2
     
    !INTEGER(kind=longint),DIMENSION(:,:), ALLOCATABLE          :: bin
    TYPE(sparsevec), DIMENSION(:), ALLOCATABLE :: bin_sparse
    INTEGER                                                    :: n

    !CALL get_binarysplit(xbs2,bin,'x_ori_bin')
    CALL get_binarysplit(xbs2,bin_sparse,'x_ori_bin')
    
    DO n = 0, xbs2%n_ori
       !bin(0,n) = IBCLR(bin(0,n),0)
       !call bin_sparse(n)%modify(0, ibclr(bin_sparse(n)%get(0),0))
       CALL bin_sparse(n)%IBCLR(0,0)
    END DO
    !CALL join_binarysplit(xbs1,xbs2,bin)
    CALL join_binarysplit(xbs1,xbs2,bin_sparse)
  END SUBROUTINE join_binsplit_va

  SUBROUTINE join_binsplit_v(xbs1,xbs2,bin_sparse)
    TYPE(binarysplit),                           INTENT(inout) :: xbs1
    TYPE(binarysplit),                           INTENT(in)    :: xbs2
    !INTEGER(kind=longint),DIMENSION(:,:), ALLOCATABLE,INTENT(in)    :: bin !0:
    TYPE(sparsevec), DIMENSION(:), ALLOCATABLE :: bin_sparse
    
    INTEGER                                                    :: n_ori,n_split
    INTEGER                                                    :: maxind
    INTEGER                                                    :: s1,is
    INTEGER                                                    :: nb1,nb2,nb
    INTEGER                                                    :: n,ik,sk
    INTEGER                                                    :: pos,pos_n
    INTEGER                                                    :: pos_s,pos_d
    INTEGER                                                    :: k_n, k
    INTEGER                                                    :: p_c,p_nc,k_nc
    INTEGER                                                    :: count
    INTEGER                                                    :: x_ori_poi
    INTEGER                                                    :: joinloc
    INTEGER                                                    :: s1bs
    INTEGER(kind=longint)                                      :: bin_sparse_last, xbs1_x_ori_bin_sparse
    INTEGER :: k_n_last, n_last, k_nc_last
    REAL    :: start, finish
    
    CALL cpu_TIME(start)
    
    xbs1 = xbs2
    n_ori = xbs2%n_ori
    n_split = xbs2%n_split
    maxind  = n_ori + n_split
    !s1  = UBOUND(bin,1)
    s1 = bin_sparse(0)%len
    !is  = BIT_size(bin(lbound(bin,1),lbound(bin,1)))
    is = bit_SIZE(bin_sparse(0)%get(0))

    nb1 = NINT(LOG(DBLE(is))/LOG(2.0_dp))
    nb2 = NINT(LOG(DBLE(s1+1))/LOG(2.0_dp)) 
    nb  = nb1 + nb2
    
    !s1bs = ubound(xbs1%x_ori_bin,1)
    s1bs = xbs1%x_ori_bin_sparse(0)%len
    IF (s1bs .GT. s1) THEN
       IF (binarysplit_message .GT. 0) &
            PRINT *, 'Message join_binarysplit: size not ok'
    ELSEIF (s1bs .LT. s1) THEN
       DO
          CALL reposition_binsplit(xbs1)
          !s1bs = UBOUND(xbs1%x_ori_bin,1)
          s1bs = xbs1%x_ori_bin_sparse(0)%len
          IF (s1bs .GE. s1) EXIT
       END DO
    END IF
    
    n_last = -1
    k_n    = -1
    DO n = 0,n_ori
       DO sk = 0,nb-1
          pos_s = 2**sk
          pos_d = pos_s * 2
          ik = -1
          DO 
             ik = ik + 1
             pos = pos_s + ik * pos_d
             IF (pos .GT. is*(s1+1)-1) EXIT
             pos_n = MOD(pos,is)
             k_n   = FLOOR(DBLE(pos)/DBLE(is))

             IF ((n .NE. n_last) .OR. (k_n .NE. k_n_last)) THEN
                n_last = n
                k_n_last = k_n
                bin_sparse_last = bin_sparse(n)%get(k_n)
             ELSE
                !write (*,*) "Reading binsplit from cache."
             END IF
             
             IF (BTEST(bin_sparse_last, pos_n)) THEN             
             !if (btest(bin_sparse(n)%get(k_n),pos_n)) then             
!             IF (BTEST(bin(k_n,n),pos_n)) THEN
                x_ori_poi = xbs1%x_ori_poi(n)
                n_split = n_split - 1
                maxind  = maxind - 1
                xbs1%n_split = n_split
                
                !xbs1%x_ori_bin(k_n,n) = ibclr(xbs1%x_ori_bin(k_n,n),pos_n)
                !**********************************************************
                ! Sparse
                !**********************************************************
                !call xbs1%x_ori_bin_sparse(n)%modify(k_n, &
                !     ibclr(xbs1%x_ori_bin_sparse(n)%get(k_n),pos_n))
                CALL xbs1%x_ori_bin_sparse(n)%IBCLR(k_n,pos_n)
                
                DO k = n, n_ori
                   xbs1%x_ori_poi(k) = xbs1%x_ori_poi(k) - 1
                END DO
                count = 0
                k_nc_last = -1
                DO p_c = pos-1,1,-1
                   p_nc = MOD(p_c,is)
                   k_nc   = FLOOR(DBLE(p_c)/DBLE(is))
                   !if (btest(xbs1%x_ori_bin(k_nc,n),p_nc)) count = count + 1
                   !**********************************************************
                   ! Sparse
                   !**********************************************************
                   IF (k_nc .NE. k_nc_last) THEN
                      xbs1_x_ori_bin_sparse = xbs1%x_ori_bin_sparse(n)%get(k_nc)
                      k_nc_last = k_nc
                   END IF
                   
                   IF (BTEST(xbs1_x_ori_bin_sparse, p_nc)) count = count + 1
                   !write (*,*) "Loop 2 in join_binsplit_v", k_nc
                END DO
                joinloc = x_ori_poi - count - 1
                xbs1%x_poi(joinloc:maxind) = xbs1%x_poi(joinloc+1:maxind+1)
                xbs1%x_poi(maxind+1) = 0.0_dp
                xbs1%x_split(joinloc:maxind) = xbs1%x_split(joinloc+1:maxind+1)
                xbs1%x_split(maxind+1) = 0.0_dp
                xbs1%x_split(joinloc) = xbs1%x_split(joinloc) - 1
                xbs1%x_pos(joinloc:maxind) = xbs1%x_pos(joinloc+1:maxind+1)
                xbs1%x_pos(maxind+1) = 0.0_dp
                xbs1%x(joinloc:maxind) = xbs1%x(joinloc+1:maxind+1)
                xbs1%x(maxind+1) = 0.0_dp
                xbs1%y(joinloc:maxind) = xbs1%y(joinloc+1:maxind+1)
                xbs1%y(maxind+1) = 0.0_dp
              
                xbs1%INT(joinloc) = xbs1%INT(joinloc) + xbs1%INT(joinloc+1)
                xbs1%INT(joinloc+1:maxind) = xbs1%INT(joinloc+2:maxind+1) 
                xbs1%err(joinloc) = xbs1%err(joinloc) + xbs1%err(joinloc+1)
                xbs1%err(joinloc+1:maxind) = xbs1%err(joinloc+2:maxind+1) 
                
                CALL eval_bsfitjoin_external(xbs1%x(0:maxind),joinloc)
            END IF
          END DO
       END DO
    END DO
    CALL reallocate_binarysplit(xbs1,1)

    CALL cpu_TIME(finish)
    PRINT '("Time in join_binsplit_v = ",f12.3," seconds.")',finish-start
  END SUBROUTINE join_binsplit_v

  ! difference (what splits are in 1 and not in 2)
  SUBROUTINE compare_binsplit(xbs1,xbs2,bin_sparse,c)
    TYPE(binarysplit),                           INTENT(inout) :: xbs1
    TYPE(binarysplit),                           INTENT(inout) :: xbs2
    !integer(kind=longint),dimension(:,:),allocatable,  intent(inout) :: bin !0:
    TYPE(sparsevec), DIMENSION(:), ALLOCATABLE, INTENT(inout)  :: bin_sparse
    CHARACTER(len=*),                            INTENT(in)    :: c

    INTEGER                                                    :: s1_1,s1_2,s1
    INTEGER                                                    :: n,n_ori
    INTEGER                                                    :: is,k,i
    INTEGER(kind=longint)                                      :: xbs1_temp, xbs2_temp

    REAL :: start, finish
    
    !s1_1 = UBOUND(xbs1%x_ori_bin,1)
    !s1_2 = ubound(xbs2%x_ori_bin,1)
    s1_1 = xbs1%x_ori_bin_sparse(0)%len
    s1_2 = xbs2%x_ori_bin_sparse(0)%len

    CALL cpu_TIME(start)
    
    DO
       IF (s1_1 .GE. s1_2) EXIT
       CALL reposition_binsplit(xbs1)
       !s1_1 = ubound(xbs1%x_ori_bin,1)
       s1_1 = xbs1%x_ori_bin_sparse(0)%len
    END DO
    DO
       IF (s1_2 .GE. s1_1) EXIT
       CALL reposition_binsplit(xbs2)
       !s1_2 = UBOUND(xbs2%x_ori_bin,1)
       s1_2 = xbs2%x_ori_bin_sparse(0)%len
    END DO
    s1 = s1_1
    n_ori = xbs1%n_ori

    !IF (ALLOCATED(bin)) DEALLOCATE(bin)
    !ALLOCATE(bin(0:s1,0:n_ori))
    !bin = 0
    !is  = BIT_SIZE(bin(0,0))
    IF (ALLOCATED(bin_sparse)) DEALLOCATE(bin_sparse)
    ALLOCATE(bin_sparse(0:n_ori))
    bin_sparse%len = s1
    is   = bit_SIZE(bin_sparse(0)%get(0))    

    DO n = 0,n_ori
       IF (c .EQ. 'diffdir') THEN
          DO k = s1, 0, -1
             !do k = lbound(xbs1%x_ori_bin_sparse(n)%idxvec,1), ubound(xbs1%x_ori_bin_sparse(n)%idxvec,1)
             DO i = 0, is-1, 1
                !IF ( BTEST(xbs1%x_ori_bin(k,n),i) .AND. &
                !     .NOT. BTEST(xbs2%x_ori_bin(k,n),i) ) THEN
                !   bin(k,n) = IBSET(bin(k,n),i)
                !END IF
                IF ( BTEST(xbs1%x_ori_bin_sparse(n)%get(k),i) .AND. &
                     .NOT. BTEST(xbs2%x_ori_bin_sparse(n)%get(k),i) ) THEN
                   !call bin_sparse(n)%modify(xbs1%x_ori_bin_sparse(n)%idxvec(k), &
                   !     ibset(bin_sparse(n)%get(xbs1%x_ori_bin_sparse(n)%idxvec(k)),i))
                   CALL bin_sparse(n)%IBSET(k, i)
                END IF
             END DO
          END DO
       ELSEIF (c .EQ. 'diff') THEN
          !do k = s1, 0, -1
          DO k = 1, xbs1%x_ori_bin_sparse(n)%len_sparse
             !bin(k,n) = IEOR(xbs1%x_ori_bin(k,n), &
             !     IAND(xbs1%x_ori_bin(k,n),xbs2%x_ori_bin(k,n)))
             xbs1_temp = xbs1%x_ori_bin_sparse(n)%get(xbs1%x_ori_bin_sparse(n)%idxvec(k))
             xbs2_temp = xbs2%x_ori_bin_sparse(n)%get(xbs1%x_ori_bin_sparse(n)%idxvec(k))
             CALL bin_sparse(n)%modify(xbs1%x_ori_bin_sparse(n)%idxvec(k), &
                  IEOR(xbs1_temp, IAND(xbs1_temp, xbs2_temp)))
          END DO
       END IF
    END DO

    CALL cpu_TIME(finish)
    PRINT '("Time in compare_binsplit_v = ",f12.3," seconds.")',finish-start
    
  END SUBROUTINE compare_binsplit

  ! print summary
  SUBROUTINE printsummary_binsplit(xbs)
    TYPE(binarysplit),                           INTENT(in) :: xbs

    INTEGER :: maxind
    maxind = xbs%n_ori + xbs%n_split
    PRINT *, 'Sum Error: ',SUM(xbs%err(0:maxind),1)
    PRINT *, 'Max Error: ',MAXVAL(xbs%err,1)
    PRINT *, 'Integral:  ',SUM(xbs%int,1)
    PRINT *, 'Size:      ',SIZE(xbs%x,1)
    PRINT *, 'Datasize:  ',maxind+1
    PRINT *, ' '
  END SUBROUTINE printsummary_binsplit

  ! print bin
  SUBROUTINE printbin_binsplit(xbs)
    TYPE(binarysplit),                           INTENT(in) :: xbs
    INTEGER                                                 :: k
    !DO k = LBOUND(xbs%x_ori_bin,2), UBOUND(xbs%x_ori_bin,2)
    !   CALL disp(xbs%x_ori_bin(:,k))
    !END DO

    !**********************************************************
    ! Sparse
    !**********************************************************
    DO k = LBOUND(xbs%x_ori_bin_sparse,1), UBOUND(xbs%x_ori_bin_sparse,1)
       !call disp(xbs%x_ori_bin(:,k))
    END DO

    WRITE (*,*) "Print sparse binsplit not yet implemented"
    STOP
    
    PRINT *, ' '
  END SUBROUTINE printbin_binsplit

  SUBROUTINE printbin_binsplit_bin(bin_sparse)
    !integer(kind=longint), dimension(:,:),             intent(in) :: bin
    TYPE(sparsevec), DIMENSION(:),             INTENT(in) :: bin_sparse 

    INTEGER                                                 :: k
    !DO k = LBOUND(bin,2), UBOUND(bin,2)
    !   CALL disp(bin(:,k))
    !END DO

    WRITE (*,*) "Print sparse binsplit not yet implemented"
    STOP
    PRINT *, ' '
  END SUBROUTINE printbin_binsplit_bin

  ! create plotfile
  SUBROUTINE plotfile_binsplit(xbs,filename,unitno)
    TYPE(binarysplit),                           INTENT(in) :: xbs
    CHARACTER(len=*),                            INTENT(in) :: filename
    INTEGER,                                     INTENT(in) :: unitno
    INTEGER                                                 :: k
    OPEN(file=filename,unit=unitno)
    DO k = 0, xbs%n_ori+xbs%n_split
       WRITE(unitno,*) xbs%x(k),xbs%y(k),xbs%INT(k),xbs%err(k)
    END DO
    CLOSE(unit=unitno)
  END SUBROUTINE plotfile_binsplit

  SUBROUTINE plotfile_binsplit_ori(xbs,filename,unitno,n_points)
    TYPE(binarysplit),                           INTENT(in) :: xbs
    CHARACTER(len=*),                            INTENT(in) :: filename
    INTEGER,                                     INTENT(in) :: unitno
    INTEGER,                                     INTENT(in) :: n_points
    INTEGER                                                 :: k,maxdim
    
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE                :: x,y
    
    maxdim = xbs%n_ori + xbs%n_split
    CALL linspace(xbs%x(0),xbs%x(maxdim),n_points,x)
    ALLOCATE(y(0:UBOUND(x,1)))
    CALL eval_bsfunc(x(:),y(:))
    OPEN(file=filename,unit=unitno)
    DO k = 0, UBOUND(x,1)
       WRITE(unitno,*) x(k),y(k)
    END DO
    CLOSE(unit=unitno)
    DEALLOCATE(x)
    DEALLOCATE(y)
  END SUBROUTINE plotfile_binsplit_ori

  ! tests interpolation
  SUBROUTINE plotfile_binsplit_test(xbs,filename,unitno,c)
    TYPE(binarysplit),                           INTENT(in) :: xbs
    CHARACTER(len=*),                            INTENT(in) :: filename
    INTEGER,                                     INTENT(in) :: unitno
    CHARACTER(len=*),                            INTENT(in) :: c
    INTEGER                                                 :: k,maxdim
    
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE                :: x,y,int0,err
    
    maxdim = xbs%n_ori + xbs%n_split
    CALL get_binarysplit(xbs,x,'x')
    CALL get_binarysplit(xbs,y,'y')
    CALL get_binarysplit(xbs,int0,'int')
    CALL get_binarysplit(xbs,err,'err')
    CALL eval_bsinterr_test(x(:),y(:),int0(:),err(:))
    OPEN(file=filename,unit=unitno)
    DO k = 0, UBOUND(x,1)
       WRITE(unitno,*) x(k),y(k)
    END DO
    CLOSE(unit=unitno)
    DEALLOCATE(x)
    DEALLOCATE(y)
  END SUBROUTINE plotfile_binsplit_test

  ! testing
  SUBROUTINE test_binsplit(test_join,test_printbin)
    INTEGER, INTENT(in) :: test_join 
    INTEGER, INTENT(in) :: test_printbin

    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: x0_1
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: s_1
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: x0_2
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: s_2

    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: x_ori
    INTEGER(kind=longint),DIMENSION(:,:), ALLOCATABLE :: bin_1,bin_2,bin_3
    INTEGER :: n,i
    REAL(kind=dp) :: d,c,g
    TYPE(binarysplit)                        :: cbs_1,cbs_2,cbs_3
    TYPE(binarysplit)                        :: cbs_1a,cbs_2a

    ! create linspace for tests
    !CALL linspace(0.0_dp,1.0_dp,20,x_ori)
    CALL linspace(0.0_dp,1.054104448087870_dp,16,x_ori)


    
    ALLOCATE( x0_1(5), s_1(5) )
    x0_1=(/0.2_dp,0.35_dp,0.8_dp,0.6_dp,0.5_dp/)
    s_1 =(/0.001_dp,0.07_dp,0.01_dp,0.02_dp,0.015_dp/)
!!$    ALLOCATE( x0_2(4), s_2(4) )
!!$    x0_2=(/0.2_dp,0.35_dp,0.8_dp,0.6_dp/)
!!$    s_2 =(/0.001_dp,0.07_dp,0.01_dp,0.02_dp/)
    ALLOCATE( x0_2(1), s_2(1) )
    !x0_2=(/0.468_dp/)
    !x0_2=(/0.473684210526316/)!-5.0d-2
    !s_2 =(/1.d-5/)

    !x0_2=(/0.96196138353943894d0/)-2.0d-1
    x0_2=(/0.01/)
    x0_2 = x_ori(LBOUND(x_ori,1)+1)+0.02
    x0_2 = x_ori(UBOUND(x_ori,1)-7)+0.1
    !x0_2 = ( x_ori(lbound(x_ori,1)) + x_ori(lbound(x_ori,1)+1) ) / 2.0d0
    !x0_2 = ( x_ori(ubound(x_ori,1)) + x_ori(ubound(x_ori,1)-1) ) / 2.0d0
    !x0_2 = x_ori(ubound(x_ori,1)) - 0.01
    s_2 =(/7.88038877731418762d-6/)
    s_2 = (/1.0d-2/)
    PRINT *, 'x0,s ',x0_2,s_2


    ! settings
    bsfunc_modelfunc = 1

    ! create limits 
    CALL create_bslimit(                      &
         total_err = 0.0d0,                   &
         local_err = 1.0d-1,                  &
         min_distance = 0.0d0,                &
         max_index = 1000,                    &
         max_splitlevel = 32                  &
         )
    !binarysplit_fsplitdepth = 0
    

!!$    ! do it for a first type
!!$    CALL construct_bsfunc(x0_1,s_1)
!!$    CALL construct_binarysplit(x_ori,cbs_1)
!!$    CALL find_binarysplit(cbs_1,x0_1)

    
    ! do it for cbs_2
    !binarysplit_fsplitdepth = 0
    bsfunc_base_distance = MAXVAL(x_ori(1:UBOUND(x_ori,1))-x_ori(0:UBOUND(x_ori,1)-1))
    PRINT *, 'bsfunc_base_distance ',bsfunc_base_distance
    CALL construct_binarysplit(x_ori,cbs_2)
    bsfunc_modelfunc = 1
    CALL multiple_binarysplit(cbs_2,x0_2,s_2)
    CALL printsummary_binarysplit(cbs_2)
    CALL plotfile_binarysplit(cbs_2,'gauss_1.dat',9)

!!$    bsfunc_modelfunc = 2
!!$    CALL construct_bsfunc(x0_2,s_2)
!!$    CALL find_binarysplit(cbs_2,x0_2)
!!$    bsfunc_modelfunc = 3
!!$    CALL construct_bsfunc(x0_2,s_2)
!!$    CALL find_binarysplit(cbs_2,x0_2)
!!$    bsfunc_modelfunc = 1

!!$    ! compare
!!$    IF (test_printbin .EQ. 1) THEN
!!$       CALL compare_binarysplit(cbs_1,cbs_2,bin_3,'diff')
!!$       CALL printbin_binarysplit(cbs_1)
!!$       CALL printbin_binarysplit(cbs_2)
!!$       CALL printbin_binarysplit(bin_1)
!!$    END IF
!!$
!!$    
!!$    ! Joining or splitting
!!$    CALL compare_binarysplit(cbs_1,cbs_2,bin_1,'diff')
!!$    CALL compare_binarysplit(cbs_2,cbs_1,bin_2,'diff')
!!$    IF (test_join .EQ. 1) THEN
!!$       CALL join_binarysplit(cbs_1a,cbs_1,bin_1)
!!$       CALL join_binarysplit(cbs_2a,cbs_2,bin_2)
!!$    ELSE
!!$       !cbs_1a = cbs_1
!!$       CALL dosplit_binarysplit(cbs_1a,cbs_1,bin_2)
!!$       CALL dosplit_binarysplit(cbs_2a,cbs_2,bin_1)
!!$    END IF
!!$    
!!$    CALL printsummary_binarysplit(cbs_1)
!!$    CALL printsummary_binarysplit(cbs_1a)
!!$    CALL printsummary_binarysplit(cbs_2)
!!$    CALL printsummary_binarysplit(cbs_2a)
!!$    
!!$    CALL plotfile_binarysplit(cbs_1,'gaussb_1.dat',9)
!!$    CALL plotfile_binarysplit(cbs_2,'gaussb_2.dat',9)
    
    CALL destruct_bsfunc

!!$    CALL deconstruct_binarysplit(cbs_1)
    CALL deconstruct_binarysplit(cbs_2)
!!$    CALL deconstruct_binarysplit(cbs_1a)
!!$    CALL deconstruct_binarysplit(cbs_2a)
!!$    CALL deconstruct_binarysplit(cbs_3)
    
  END SUBROUTINE test_binsplit


END MODULE binarysplit_mod
