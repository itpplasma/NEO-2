module binarysplit_int

  use hdf5_tools
  use hdf5_tools_f2003
  use nrtype, only : dp

  implicit none

  public bsfunc_message
  integer         :: bsfunc_message = 0

  integer, public :: bsfunc_modelfunc = 1


  private bsfunc_evaldegree
  integer         :: bsfunc_evaldegree = 2 ! 1   ! 3 or 4

  ! parameters for evaluation functions
  ! have to be set with (.e.g. for gauss)
  !   call construct_bsfunc(x0,s)
  private bsfunc_p1,bsfunc_p2
  real(kind=dp), dimension(:), allocatable :: bsfunc_p1,bsfunc_p2

  ! parameters for limit functions
  ! have to be set with
  !   call create_bslimit( ..... )
  ! optional arguments
  !   total_err,local_err,min_distance,max_index,max_splitlevel
  real(kind=dp), public :: bsfunc_total_err
  real(kind=dp), public :: bsfunc_local_err
  real(kind=dp), public :: bsfunc_min_distance
  integer,       public :: bsfunc_max_index
  integer,       public :: bsfunc_max_splitlevel

  real(kind=dp), public :: bsfunc_base_distance = 0.05d0
  real(kind=dp), public :: bsfunc_mult_constant = 0.1d0
  real(kind=dp), public :: bsfunc_sigma_multiplier = 1.147202690439877 !sqrt(3.0d0)




  ! additional quantities which belon logically to splitting of levels
  ! but they are not used internally
  real(kind=dp), public                    :: bsfunc_sigma_mult
  real(kind=dp), public                    :: bsfunc_sigma_min
  integer,       public                    :: bsfunc_local_solver

  public construct_bsfunc
  private construct_bsfunc_2
  interface construct_bsfunc
    module procedure construct_bsfunc_2
  end interface

  public destruct_bsfunc
  private destruct_bsfunc_2
  interface destruct_bsfunc
    module procedure destruct_bsfunc_2
  end interface

  public eval_bsfunc
  private eval_bsfunc_d1, eval_bsfunc_d2
  interface eval_bsfunc
    module procedure eval_bsfunc_d1, eval_bsfunc_d2
  end interface

  public eval_bsfitsplit
  private eval_bsfitsplit_4
  interface eval_bsfitsplit
    module procedure eval_bsfitsplit_4
  end interface

  public eval_bsfitsplit_external
  private eval_bsfitsplit_ext
  interface eval_bsfitsplit_external
    module procedure eval_bsfitsplit_ext
  end interface

  public eval_bsfitjoin_external
  private eval_bsfitjoin_ext
  interface eval_bsfitjoin_external
    module procedure eval_bsfitjoin_ext
  end interface

  public linspace
  private linspace_d
  interface linspace
    module procedure linspace_d
  end interface

  public eval_bsinterr
  private eval_bsinterr_d1
  interface eval_bsinterr
    module procedure eval_bsinterr_d1,eval_bsinterr_test
  end interface

  public eval_bslimit
  private eval_bslimit_1, eval_bslimit_3
  interface eval_bslimit
    module procedure eval_bslimit_1, eval_bslimit_3
  end interface

  public create_bslimit
  private create_bslimit_opt
  interface create_bslimit
    module procedure create_bslimit_opt
  end interface

  public  disp
  private disp_i8, disp_i81, disp_c
  interface disp
    module procedure disp_i8, disp_i81, disp_c
  end interface

  private gauss
  private gauss_d0, gauss_d1
  interface gauss
    module procedure gauss_d0, gauss_d1
  end interface

  private plagrange_coeff
  private plag_coeff
  interface plagrange_coeff
    module procedure plag_coeff
  end interface

  private plagrange_stencil
  private plag_stencil
  interface plagrange_stencil
    module procedure plag_coeff
  end interface

  public plagrange_test
  private plag_test
  interface plagrange_test
    module procedure plag_test
  end interface

contains

  ! construct
  ! puts the parameters p1, p2
  ! into the private variables bsfunc_p1, .....
  ! which are used internally
  subroutine construct_bsfunc_2(p1,p2)
    real(kind=dp), dimension(:), allocatable, intent(in) :: p1,p2

    if (allocateD(bsfunc_p1)) deallocate(bsfunc_p1)
    allocate(bsfunc_p1(size(p1,1)))
    bsfunc_p1 = p1

    if (allocateD(bsfunc_p2)) deallocate(bsfunc_p2)
    allocate(bsfunc_p2(size(p2,1)))
    bsfunc_p2 = p2
  end subroutine construct_bsfunc_2

  ! destruct
  subroutine destruct_bsfunc_2()
    if (allocateD(bsfunc_p1)) deallocate(bsfunc_p1)
    if (allocateD(bsfunc_p2)) deallocate(bsfunc_p2)
  end subroutine destruct_bsfunc_2

  ! evaluation
  ! here the gauss function is called for the 2 parameters
  !  x0 (bsfunc_p1) and s (bsfunc_p1)
  !  has to be changed if another function has to be used
  subroutine eval_bsfunc_d1(x,y)
    real(kind=dp), intent(in)    :: x
    real(kind=dp), intent(inout) :: y
    call gauss(x,bsfunc_p1,bsfunc_p2,y)
  end subroutine eval_bsfunc_d1

  subroutine eval_bsfunc_d2(x,y)
    real(kind=dp), dimension(:), intent(in)    :: x ! 0:
    real(kind=dp), dimension(:), intent(inout) :: y ! 0:
    integer                                    :: k
    do k = LBOUND(x,1), UBOUND(x,1)
      call eval_bsfunc(x(k),y(k))
    end do
  end subroutine eval_bsfunc_d2

  ! procedure when you split
  ! handles all tasks when a level is split into two
  !
  !  internal stuff: function value and integeral
  !  with polynomial of degree 3
  !
  ! should be adopted to the interr
  subroutine eval_bsfitsplit_4(x,y,i,e,loc)
    real(kind=dp), dimension(0:), intent(in)    :: x ! 0:
    real(kind=dp), dimension(0:), intent(inout) :: y ! 0:
    real(kind=dp), dimension(0:), intent(inout) :: i ! 0:
    real(kind=dp), dimension(0:), intent(inout) :: e ! 0:
    integer,                      intent(in)   :: loc

    real(kind=dp) :: dx,dx2,dx3,dx4,dy,dydx1,dydx2,b,c

    ! internal stuff (parabola of degree 3)
    if (loc .GT. 2) then
      dydx1 = (y(loc-1)-y(loc-2)) / (x(loc-1)-x(loc-2))
    else
      dydx1 = (y(loc)-y(loc-1)) / (x(loc)-x(loc-1))
    end if
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
    call eval_bsfitsplit_external(x,loc)
  end subroutine eval_bsfitsplit_4

  ! external procedure when you split
  ! has to be adopted for your problem
  subroutine eval_bsfitsplit_ext(x,loc)
    real(kind=dp), dimension(0:), intent(in)    :: x !0:
    integer,                     intent(in)    :: loc
    interface join_ripples_bsfitsplit
      subroutine join_ripples_bsfitsplit(x,loc)
        use nrtype, only : dp
        real(kind=dp), dimension(0:), intent(in)    :: x !0:
        integer,                     intent(in)    :: loc
      end subroutine join_ripples_bsfitsplit
    end interface
    ! add here what is necessary for external stuff
    call join_ripples_bsfitsplit(x,loc) !EXTERNAL
  end subroutine eval_bsfitsplit_ext

  ! external procedure when you join
  ! has to be adopted for your problem
  subroutine eval_bsfitjoin_ext(x,loc)
    real(kind=dp), dimension(0:), intent(in)    :: x !0:
    integer,                      intent(in)   :: loc
    interface join_ripples_bsfitjoin
      subroutine join_ripples_bsfitjoin(x,loc)
        use nrtype, only : dp
        real(kind=dp), dimension(0:), intent(in)    :: x !0:
        integer,                      intent(in)    :: loc
      end subroutine join_ripples_bsfitjoin
    end interface
    ! add here what is necessary for external stuff
    call join_ripples_bsfitjoin(x,loc) !EXTERNAL
  end subroutine eval_bsfitjoin_ext

  ! limits
  subroutine create_bslimit_opt(total_err,local_err,min_distance, &
       max_index,max_splitlevel)
    real(kind=dp), intent(in), OPTIONAL :: total_err
    real(kind=dp), intent(in), OPTIONAL :: local_err
    real(kind=dp), intent(in), OPTIONAL :: min_distance
    integer,       intent(in), OPTIONAL :: max_index
    integer,       intent(in), OPTIONAL :: max_splitlevel

    if (present(total_err)) then
      bsfunc_total_err = total_err
    else
       bsfunc_total_err = 0.0_dp
    end if
    if (present(local_err)) then
      bsfunc_local_err = local_err
    else
      bsfunc_local_err = 0.0_dp
    end if
    if (present(min_distance)) then
      bsfunc_min_distance = min_distance
    else
      bsfunc_min_distance = 0.0_dp
    end if
    if (present(max_index)) then
      bsfunc_max_index = max_index
    else
      bsfunc_max_index = 1000
    end if
    if (present(max_splitlevel)) then
      bsfunc_max_splitlevel = max_splitlevel
    else
      bsfunc_max_splitlevel = 32
    end if
  end subroutine create_bslimit_opt

  function eval_bslimit_1(err) result(l)
    real(kind=dp), dimension(0:), intent(in) :: err  !0:
    logical                                 :: l, l1, l2
    !l1 = SUM(err,1) .GT. bsfunc_total_err
    l2 = MAXVAL(err,1) .GT. bsfunc_local_err
    !print *, 'Local Error ',MAXVAL(err,1),bsfunc_local_err
    !l  = l1 .OR. l2
    l = l2
!!$    if (.NOT. l1) then
!!$       if (bsfunc_message .GT. 0) then
!!$          print *, 'Message: Total Error Limits reached!'
!!$       end if
!!$    end if
    if (.NOT. l2) then
      if (bsfunc_message .GT. 0) then
        print *, 'Message: Local Error Limits reached!'
      end if
    end if
  end function eval_bslimit_1

  function eval_bslimit_3(dist,split,maxind) result(l)
    real(kind=dp),                intent(in) :: dist
    integer,                      intent(in) :: split
    integer,                      intent(in) :: maxind
    logical                                  :: l, l1, l2, l3
    l1 = dist   .GE. bsfunc_min_distance
    l2 = split  .LE. bsfunc_max_splitlevel
    l3 = maxind .LE. bsfunc_max_index
    l  = l1 .AND. l2 .AND. l3
    if (.NOT. l1) then
      if (bsfunc_message .GT. 0) print *, 'Message: Minimum distance reached!'
    end if
    if (.NOT. l2) then
      if (bsfunc_message .GT. 0) print *, 'Message: Maximum split level reached!'
    end if
    if (.NOT. l3) then
      print *, 'Message: Maximum index reached!'
      stop
    end if
  end function eval_bslimit_3

  ! linspace (helper routine)
  subroutine linspace_d(x_beg,x_end,x_num,x)
    real(kind=dp),                            intent(in)     :: x_beg
    real(kind=dp),                            intent(in)     :: x_end
    integer,                                  intent(in)     :: x_num
    real(kind=dp), dimension(:), allocatable, intent(inout)  :: x !0:

    integer                                                  :: i_beg_l
    integer                                                  :: i_end_l
    integer                                                  :: k,i
    real(kind=dp)                                            :: x_del

    i_beg_l = 0
    i_end_l = x_num + i_beg_l -1
    if (x_num .LE. 1) i_end_l = i_beg_l
    if (allocateD(x)) deallocate(x)
    allocate (x(i_beg_l:i_end_l))

    if (x_num .LE. 1) then
      x(i_beg_l) = x_beg
    else
      x_del = (x_end-x_beg) / DBLE(x_num-1)

      i = i_beg_l
      x(i) = x_beg
      do k = 2, x_num-1
        i = i + 1
        x(i) = x(i-1) + x_del
      end do
      x(i_end_l) = x_end
    end if
  end subroutine linspace_d


  ! gauss (helper routine)
  subroutine gauss_d0(x,x0,s,g)
    use math_constants, only : pi, sqrt2

    real(kind=dp),               intent(in)   :: x
    real(kind=dp), dimension(:), intent(in)   :: x0
    real(kind=dp), dimension(:), intent(in)   :: s
    real(kind=dp),               intent(out)  :: g


    real(kind=dp),               parameter    :: os2pi=0.39894228040143_dp
    real(kind=dp),               parameter    :: os8=0.35355339059327_dp

    integer                                   :: k

    g = 0.0_dp
    if (bsfunc_modelfunc .EQ. 1) then
       ! Winny - This is now the return of the normalized Gauss-function
       ! changed to -> NO NORMALIZATION
      do k = 1, size(x0,1)
        g = g + EXP(- (x-x0(k))**2 / s(k)**2 / 2.0_dp) ! / s(k)
      end do
       ! g = g * os2pi

!!$       if (bsfunc_evaldegree .NE. 2) then
!!$          do k = 1, size(x0,1)
!!$             g = g + EXP(- ((x-x0(k)) / s(k))**2 / 2) / s(k)
!!$          end do
!!$          g = g * os2pi
!!$       else
!!$          do k = 1, size(x0,1)
!!$             g = g + EXP(- ((x-x0(k)) / s(k))**2 / 2)
!!$          end do
!!$       end if
    elseif (bsfunc_modelfunc .EQ. 2) then
      do k = 1, size(x0,1)
        g = g + SQRT( s(k)**2/(2.0_dp*(x-x0(k))**2+s(k)**2) )
      end do
       !g = g/size(x0,1)
!       if (bsfunc_evaldegree .NE. 2) then
!          do k = 1, size(x0,1)
!             !->out          g = g + EXP(- ABS(x-x0(k)) / s(k) / sq2) / s(k)
!             g=g+s(k)/((x-x0(k))**2+s(k)**2/0.12732200375004d0)              !<-in
!          end do
!          !->out       g = g * os8
!          g = g / (pi*0.35682208977310d0)
!       else
!          do k = 1, size(x0,1)
!             g=g+s(k)/((x-x0(k))**2+s(k)**2/0.12732200375004d0) / 0.12732200375004d0 * s(k)
!          end do
!       end if
    elseif (bsfunc_modelfunc .EQ. 3) then
      if (bsfunc_evaldegree .NE. 2) then
        do k = 1, size(x0,1)
          !->out          g = g + s(k) / ( (x-x0(k))**2 + s(k)**2 )
          g=g+s(k)/((x-x0(k))**2+s(k)**2/0.87267799624997d0)              !<-in
        end do
          !->out       g = g / pi
          g = g / (pi*0.93417235896272)                                      !<-in
      else
        do k = 1, size(x0,1)
          g=g+s(k)/((x-x0(k))**2+s(k)**2/0.87267799624997d0) / 0.87267799624997d0 * s(k)
        end do
      end if
    elseif (bsfunc_modelfunc .EQ. 4) then
      g = 1.0_dp
      do k = 1, size(x0,1)
        g = g * (1.0_dp + EXP(- ((x-x0(k)) / s(k))**2 / 2) / s(k))
      end do
    elseif (bsfunc_modelfunc .EQ. 5) then
      g = 1.0_dp
      do k = 1, size(x0,1)
        !g = g * (1.0_dp + EXP(- ((x-x0(k)) / s(k))**2 / 2) / s(k))
        g = g * sign(x-x0(k),1.0d0) * (s(k) * ABS(x-x0(k)))**(8.0d0/8.0d0)
      end do
      g = (s(k) * ABS(x0(k)))**(8.0d0/8.0d0) - g
    else
      print *,'Error from binarysplit: bsfunc_modelfunc wrong: ', &
          bsfunc_modelfunc
      stop
    end if
  end subroutine gauss_d0

  subroutine gauss_d1(x,x0,s,g)
    real(kind=dp), dimension(:), allocatable, intent(in)    :: x !0:
    real(kind=dp), dimension(:), allocatable,  intent(in)    :: x0
    real(kind=dp), dimension(:), allocatable,  intent(in)    :: s
    real(kind=dp), dimension(:), allocatable, intent(inout) :: g !0:

    integer                                     :: k

    do k = LBOUND(x,1), UBOUND(x,1)
      call gauss(x(k),x0,s,g(k))
    end do
  end subroutine gauss_d1

  ! evaluation of error for splitting
  subroutine eval_bsinterr_d1(x,y,int0,err,splitloc)
    !USE plagrange_mod
    real(kind=dp), dimension(0:), intent(in)    :: x !0:
    real(kind=dp), dimension(0:), intent(in)    :: y !0:
    real(kind=dp), dimension(0:), intent(inout) :: int0 !0:
    real(kind=dp), dimension(0:), intent(inout) :: err !0:
    integer,                      intent(out)   :: splitloc

    integer                                   :: lb, ub, k, i
    integer                                   :: in, in1, in2
    real(kind=dp), dimension(:), allocatable  :: yp, ypp
    real(kind=dp)                             :: dx, dx2, dx3, dx4, dx5
    real(kind=dp)                             :: dy, a, b, c, d, e, aa, bb
    real(kind=dp)                             :: xloc,yloc,sint0,fac

    integer :: k1,k2,i1,i2,ii
    real(kind=dp), dimension(:), allocatable  :: xlag, ylag
    real(kind=dp) :: err_int

    real(kind=dp), dimension(:,:), allocatable :: coeff
    real(kind=dp) :: clag,dlag
    real(kind=dp) :: dloc,d1,kd
    integer, parameter :: npoi = 6
    integer, parameter :: nder = 3


    lb = LBOUND(x,1)
    ub = UBOUND(x,1)
    int0(lb) = 0.0_dp
    err(lb)  = 0.0_dp
    allocate(yp(lb:ub))
    allocate(ypp(lb:ub))

    !print *, 'eval_bsinterr_d1'
    if (bsfunc_evaldegree .EQ. 1) then
      sint0 = 0.0_dp
      in = 8
      in1 = 8
      in2 = 64 ! around max
      fac = 2.0_dp
      !do WHILE (ABS(sint0) .LT. 1.0d-4)
         !print *, 'ub, in ',ub,in
      do k = 1, ub
        dx = x(k)-x(k-1)
        if (x(k) .LT. bsfunc_p1(1) - fac*bsfunc_p2(1) .OR. &
            x(k-1) .GT. bsfunc_p1(1) + fac*bsfunc_p2(1) ) then
          int0(k) = 0.0_dp
          do i = 0, in1
            xloc = x(k-1) + DBLE(i) * dx / DBLE(in1)
            call eval_bsfunc(xloc,yloc)
            if (i .EQ. 0 .OR. i .EQ. in1) yloc = yloc / 2.0_dp
            int0(k) = int0(k) + yloc
          end do
          int0(k) = int0(k) * dx / DBLE(in1)
        elseif (x(k-1) .LE. bsfunc_p1(1) - fac*bsfunc_p2(1) .AND. &
                x(k) .GE. bsfunc_p1(1) + fac*bsfunc_p2(1) ) then
          int0(k) = 0.0_dp
          do i = 0, in2
            xloc = x(k-1) + DBLE(i) * dx / DBLE(in2)
            call eval_bsfunc(xloc,yloc)
            if (i .EQ. 0 .OR. i .EQ. in2) yloc = yloc / 2.0_dp
            int0(k) = int0(k) + yloc
          end do
          int0(k) = int0(k) * dx / DBLE(in2)
          if (int0(k) .LT. 1.d-4) then
            int0(k) = 1.0_dp
            !print *, 'set to 1.0'
          end if
        else
          int0(k) = 0.0_dp
          do i = 0, in
            xloc = x(k-1) + DBLE(i) * dx / DBLE(in)
            call eval_bsfunc(xloc,yloc)
            if (i .EQ. 0 .OR. i .EQ. in) yloc = yloc / 2.0_dp
            int0(k) = int0(k) + yloc
          end do
          int0(k) = int0(k) * dx / DBLE(in)
        end if
        err(k)  = (y(k-1)+y(k))*dx/2.0_dp
        err(k)  = ABS(err(k)-int0(k))
      end do
      sint0 = SUM(int0)
      !print *, 'int0'
      !print *, int0
      !print *, 'err'
      !print *, err
      !print *, 'sint0 ',sint0,in
      !print *, 'max err ',MAXVAL(err)
      !PAUSE
      !   in = in * 10
      !end do

    elseif (bsfunc_evaldegree .EQ. 2) then
      ! new stuff with third derivative
      err = 0.0_dp
      !print *,'lb,ub ',lbound(err,1),ubound(err,1)
      int0 = 0.0_dp
      allocate( coeff(0:nder,npoi) )
      do k = 1, ub

        call plag_stencil(ub,npoi,k,k1,k2,i1,i2)
        !print *, ub,npoi,k,k1,k2,i1,i2
        allocate(xlag(i1:i2))
        allocate(ylag(i1:i2))
        xlag = x(k1:k2)
        ylag = y(k1:k2)

        xloc = (xlag(0) + xlag(1)) / 2.0_dp
        call plagrange_coeff(npoi,nder,xloc,xlag,coeff)

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
        if (bsfunc_modelfunc .EQ. 2) then
          err(k) = err(k) /  ( ABS( SUM(coeff(0,:)*ylag) ) )
        end if
        !print *, 'k,cla,dlag,err ',k,ABS(clag),ABS(dlag),d1,err(k)
        !print *, 'xlag ',xlag
        !print *, 'ylag ',ylag
        deallocate(xlag)
        deallocate(ylag)
      end do
      !print *, ''
      !print *,'x'
      !print *,x
      !print *, 'y'
      !print *,y
      !print *,'err'
      !print *,err
      !print *, 'sumerr ',sum(err)
      !print *, 'ub splitloc maxerr', ubound(err,1),maxloc(err)-1,maxval(err)
      !print *, ''
      !pause
      deallocate(coeff)
      ! end of new stuff with third derivative

    elseif (bsfunc_evaldegree .EQ. 3) then
      do k = 1, ub
        yp(k) = (y(k)-y(k-1))/(x(k)-x(k-1))
      end do
      yp(0) = yp(1)

      do k = 1, ub
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
      end do
    else if (bsfunc_evaldegree .EQ. 4) then
      do k = 1, ub
        yp(k) = (y(k)-y(k-1))/(x(k)-x(k-1))
      end do
      yp(0) = yp(1)
      do k = 1, ub
        ypp(k) = (yp(k)-yp(k-1))/(x(k)-x(k-1))
      end do
      ypp(0) = ypp(1)


      do k = 1, ub
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
      end do
    end if
    splitloc = MAXLOC(err,1)-1
    !print *, 'loc ',splitloc, err(splitloc), SUM(err),x(splitloc)
    !print *, ' '
    !PAUSE
    deallocate(yp)
    deallocate(ypp)
  end subroutine eval_bsinterr_d1

  ! just for some test, puts back interpolated values
  subroutine eval_bsinterr_test(x,y,int0,err)
    real(kind=dp), dimension(0:), intent(inout) :: x !0:
    real(kind=dp), dimension(0:), intent(inout) :: y !0:
    real(kind=dp), dimension(0:), intent(inout) :: int0 !0:
    real(kind=dp), dimension(0:), intent(inout) :: err !0:

    integer                                   :: lb, ub, k
    real(kind=dp), dimension(:), allocatable  :: yp, ypp
    real(kind=dp)                             :: dx, dx2, dx3, dx4, dx5
    real(kind=dp)                             :: dy, a, b, c, d, e, aa, bb

    lb = LBOUND(x,1)
    ub = UBOUND(x,1)
    int0(lb) = 0.0_dp
    err(lb)  = 0.0_dp
    allocate(yp(lb:ub))
    allocate(ypp(lb:ub))

    if (bsfunc_evaldegree .EQ. 3) then
      do k = 1, ub
        yp(k) = (y(k)-y(k-1))/(x(k)-x(k-1))
      end do
      yp(0) = yp(1)

      do k = 1, ub
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
      end do
    else if (bsfunc_evaldegree .EQ. 4 .OR. bsfunc_evaldegree .EQ. 1) then
       do k = 1, ub
         yp(k) = (y(k)-y(k-1))/(x(k)-x(k-1))
       end do
       yp(0) = yp(1)
       do k = 1, ub
         ypp(k) = (yp(k)-yp(k-1))/(x(k)-x(k-1))
       end do
       ypp(0) = ypp(1)


       do k = 1, ub
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
      end do
    end if
    !splitloc = MAXLOC(err,1)-1
    deallocate(yp)
    deallocate(ypp)
  end subroutine eval_bsinterr_test




  ! display (helper routine)
  subroutine disp_i8(i)
    use nrtype, only : longint

    integer(kind=longint), intent(in)  :: i
    integer                      :: k
    integer                      :: is
    character(len=BIT_size(i))   :: c
    character(len=20)            :: form
    is = BIT_size(i)

    write(form,'(a,i2,a)') '(',is,'i1)'
    write(c,TRIM(form)) (IBITS(i,k,1), k=is-1,0,-1)
    call disp(c)
  end subroutine disp_i8

  subroutine disp_i81(i)
    use nrtype, only : longint

    integer(kind=longint), dimension(:), intent(in)  :: i

    integer                                    :: k,ks
    integer                                    :: is,s1
    character(len=size(i,1)*BIT_size(i(1)))    :: c
    character(len=20)                          :: form
    is = BIT_size(i(1))

    s1 = size(i,1)
    write(form,'(a,i2,a)') '(',s1*is,'i1)'
    write(c,TRIM(form)) ( (IBITS(i(ks),k,1),k=is-1,0,-1),ks=s1,1,-1 )
    call disp(c)
  end subroutine disp_i81

  subroutine disp_c(c)
    character(len=*), intent(in) :: c
    print *, c
  end subroutine disp_c


  subroutine plag_coeff(npoi,nder,x,xp,coef)
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

    integer, intent(in)                                :: npoi,nder
    real(kind=dp), intent(in)                          :: x
    real(kind=dp), dimension(npoi), intent(in)         :: xp
    real(kind=dp), dimension(0:nder,npoi), intent(out) :: coef
    real(kind=dp), dimension(:), allocatable           :: dummy
    real(kind=dp), dimension(:), allocatable           :: fak_i

    integer                                            :: i,k,j,l,m
    real(kind=dp)                                      :: fac
    real(kind=dp)                                      :: j_sum,l_sum,m_sum,k_prod

    do i=1,npoi
      coef(0,i)=1.d0
      do k=1,npoi
        if(k.EQ.i) cycle
        coef(0,i)=coef(0,i)*(x-xp(k))/(xp(i)-xp(k))
      end do
    end do

    if(nder.EQ.0) return

    allocate(dummy(npoi))

    do i=1,npoi
      dummy=1.d0
      dummy(i)=0.d0
      do k=1,npoi
        if(k.EQ.i) cycle
        fac=(x-xp(k))/(xp(i)-xp(k))
        do j=1,npoi
          if (j.EQ.k) then
            dummy(j)=dummy(j)/(xp(i)-xp(k))
          else
            dummy(j)=dummy(j)*fac
          endif
        end do
      end do
      coef(1,i)=SUM(dummy)
    end do

    deallocate(dummy)

    ! second derivative
    if(nder.LE.1) return

    allocate(fak_i(npoi))
    do_i: do i = 1,npoi
      fak_i = 0.0d0
      do_prep: do k = 1,npoi
        if (k .EQ. i) cycle
        fak_i(k) = (x-xp(k)) / (xp(i)-xp(k))
      end do do_prep
      j_sum = 0.0d0
      do_j: do j =1,npoi
        if (j .EQ. i) cycle
        l_sum = 0.0d0
        do_l: do l = 1,npoi
          if (l .EQ. i .OR. l .EQ. j) cycle
          k_prod = 1.0d0
          do_k: do k =1,npoi
            if (k .EQ. i .OR. k .EQ. j .OR. k .EQ. l) cycle
            k_prod = k_prod * fak_i(k)
          end do do_k
          l_sum = l_sum + k_prod / (xp(i)-xp(l))
        end do do_l
        j_sum = j_sum + l_sum / (xp(i)-xp(j))
      end do do_j
      coef(2,i)=j_sum
    end do do_i
    deallocate(fak_i)

    ! third derivative
    if(nder.LE.2) return

    allocate(fak_i(npoi))
    do_i3: do i = 1,npoi
      fak_i = 0.0d0
      do_prep3: do k = 1,npoi
        if (k .EQ. i) cycle
        fak_i(k) = (x-xp(k)) / (xp(i)-xp(k))
      end do do_prep3
      j_sum = 0.0d0
      do_j3: do j =1,npoi
        if (j .EQ. i) cycle
        l_sum = 0.0d0
        do_l3: do l = 1,npoi
          if (l .EQ. i .OR. l .EQ. j) cycle
          m_sum = 0.0d0
          do_m3: do m = 1,npoi
            if (m .EQ. i .OR. m .EQ. j .OR. m .EQ. l) cycle
            k_prod = 1.0d0
            do_k3: do k =1,npoi
              if (k .EQ. i .OR. k .EQ. j .OR. k .EQ. l .OR. k .EQ. m) cycle
              k_prod = k_prod * fak_i(k)
            end do do_k3
            m_sum = m_sum + k_prod / (xp(i)-xp(m))
          end do do_m3
          l_sum = l_sum + m_sum / (xp(i)-xp(l))
        end do do_l3
        j_sum = j_sum + l_sum / (xp(i)-xp(j))
      end do do_j3
      coef(3,i)=j_sum
    end do do_i3
    deallocate(fak_i)

  end subroutine plag_coeff

  subroutine plag_stencil(ub,npoi,k,k1,k2,i1,i2)
    integer, intent(in)  :: ub,npoi,k
    integer, intent(out) :: k1,k2,i1,i2
    integer :: kd

    k1 = k - npoi/2
    i1 = 1 - npoi/2
    if (k1 .LT. 0) then
      kd = -k1
      k1 = k1 + kd
      i1 = i1 + kd
    elseif (k1 .GT. ub-npoi+1) then
      kd = k1 - (ub-npoi+1)
      k1 = k1 - kd
      i1 = i1 - kd
    end if
    k2 = k1 + npoi - 1
    i2 = i1 + npoi - 1

  end subroutine plag_stencil

  subroutine plag_test
    use math_constants, only : pi

    integer, parameter :: unitno = 9999
    integer :: i
    real(kind=dp), dimension(:), allocatable :: x_ori,y0_ori,y1_ori,y2_ori,y3_ori
    real(kind=dp), dimension(:), allocatable :: x,y0,y1,y2,y3
    real(kind=dp), dimension(:,:), allocatable :: coeff
    real(kind=dp), dimension(:), allocatable  :: xlag, ylag
    real(kind=dp) :: xloc
    integer :: lb,ub
    integer :: k,k1,k2,i1,i2
    integer, parameter :: npoi = 6
    integer, parameter :: nder = 3
    integer, parameter :: ndata = 25

    ! basic data
    !call linspace(-pi,pi,ndata,x_ori)

    allocate(x_ori(0:ndata))
    x_ori(0) = -pi
    do k = 1,ndata
      x_ori(k) = x_ori(0) + 2*pi*(DBLE(k)/DBLE(ndata))**1
    end do

    lb = LBOUND(x_ori,1)
    ub = UBOUND(x_ori,1)
    allocate(y0_ori(lb:ub))
    allocate(y1_ori(lb:ub))
    allocate(y2_ori(lb:ub))
    allocate(y3_ori(lb:ub))
    y0_ori =  SIN(x_ori)
    y1_ori =  COS(x_ori)
    y2_ori = -SIN(x_ori)
    y3_ori = -COS(x_ori)
    ! plot basic data
    open(file='plag_ori.dat',unit=unitno)
    do i = lb,ub
      write(unitno,*) x_ori(i),y0_ori(i),y1_ori(i),y2_ori(i),y3_ori(i)
    end do
    close(unitno)

    ! lagrange coefficent
    allocate( coeff(0:nder,npoi) )


    ! interpolation data
    allocate(x(1:ub))
    allocate(y0(1:ub))
    allocate(y1(1:ub))
    allocate(y2(1:ub))
    allocate(y3(1:ub))

    do k = 1,ub

      call plag_stencil(ub,npoi,k,k1,k2,i1,i2)
      allocate(xlag(i1:i2))
      allocate(ylag(i1:i2))
      xlag = x_ori(k1:k2)
      ylag = y0_ori(k1:k2)

      xloc = (xlag(0) + xlag(1)) / 2.0_dp
      call plagrange_coeff(npoi,nder,xloc,xlag,coeff)

      x(k) = xloc
      y0(k) = SUM(coeff(0,:)*ylag)
      y1(k) = SUM(coeff(1,:)*ylag)
      y2(k) = SUM(coeff(2,:)*ylag)
      y3(k) = SUM(coeff(3,:)*ylag)
      deallocate(xlag,ylag)
    end do

    ! plot basic data
    open(file='plag_int.dat',unit=unitno)
    do i = 1,ub
      write(unitno,*) x(i),y0(i),y1(i),y2(i),y3(i)
    end do
    close(unitno)

    deallocate(x_ori)
    deallocate(y0_ori,y1_ori,y2_ori,y3_ori)
    deallocate(x)
    deallocate(y0,y1,y2,y3)
    deallocate(coeff)

  end subroutine plag_test



end module binarysplit_int
