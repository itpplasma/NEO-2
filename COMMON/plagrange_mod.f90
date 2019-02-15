module plagrange_mod

  use magnetics_mod

  implicit none

  public plagrange_coeff
  private plag_coeff
  interface plagrange_coeff
    module procedure plag_coeff
  end interface

  public plagrange_test
  private plag_test,plag_testperiod_2
  interface plagrange_test
    module procedure plag_test,plag_testperiod,plag_testperiod_2
  end interface

  public plagrange_stencel
  private plag_stencel
  interface plagrange_stencel
    module procedure plag_stencel
  end interface

  public plagrange_value
  private plag_value_1,plag_value_2,plag_value_all,plag_value_all2
  interface plagrange_value
    module procedure plag_value_1,plag_value_2,plag_value_all,&
          plag_value_all2
  end interface

  public plagrange_interp
  private plag_interp_2,plag_interp_3,plag_interp_all,plag_interp_all2
  interface plagrange_interp
    module procedure plag_interp_2,plag_interp_3,plag_interp_all,&
          plag_interp_all2
  end interface

contains

  !> \brief Determine lagrange coefficients?
  !>
  !> \param npoi [in] - number of points (determines the order of Lagrange polynomial
  !>   which is equal npoi-1)
  !> \param nder [in] - number of derivatives computed 0 - function only, 1 - first
  !>   derivative
  !> \param  x [in] - actual point where function and derivatives are evaluated
  !> \param  xp(npoi) [in] - array of points where function is known
  !> \param  coef(0:nder,npoi) [out] - weights for computation of
  !>   function and derivatives,
  !>   f=sum(fun(1:npoi)*coef(0,1:npoi) gives the function value
  !>   df=sum(fun(1:npoi)*coef(1,1:npoi) gives the value of the derivative
  subroutine plag_coeff(npoi,nder,x,xp,coef)
    use nrtype, only : dp

    integer, intent(in)                                :: npoi,nder
    real(kind=dp), intent(in)                          :: x
    real(kind=dp), dimension(npoi), intent(in)         :: xp
    real(kind=dp), dimension(0:nder,npoi), intent(out) :: coef
    real(kind=dp), dimension(:), allocatable           :: dummy
    real(kind=dp), dimension(:), allocatable           :: fak_i

    integer                                            :: i,k,j,l
    real(kind=dp)                                      :: fac
    real(kind=dp)                                      :: j_sum,l_sum,k_prod

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
          if(j.EQ.k) then
            dummy(j)=dummy(j)/(xp(i)-xp(k))
          else
            dummy(j)=dummy(j)*fac
          end if
        end do
      end do
      coef(1,i)=SUM(dummy)
    end do

    deallocate(dummy)

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

  end subroutine plag_coeff

  !--------------------------------------------------------------------
  subroutine plag_test(npoi,imax)
    use math_constants, only : pi
    use nrtype, only : dp

    integer, intent(in) :: npoi
    integer, intent(in) :: imax

    integer, parameter                         :: nder=1
    integer                                    :: i
    real(kind=dp)                              :: u,umax
    real(kind=dp), dimension(:), allocatable   :: up,fun
    real(kind=dp), dimension(:,:), allocatable :: coeff

    umax=pi

    allocate( up(npoi) )
    allocate( fun(npoi) )
    allocate( coeff(0:nder,npoi) )

    open(1,file='input_plag.dat')
    do i=1,npoi
      up(i)=umax*(float(i-1)/float(npoi-1))
      fun(i)=SIN(up(i))
      write (1,*) up(i),fun(i)
    end do
    close(1)

    open(1,file='sinus.dat')
    do i=1,imax
      u=umax*(float(i-1)/float(imax-1))
      call plagrange_coeff(npoi,nder,u,up(:),coeff(:,:))
      write (1,*) u,SIN(u),SUM(coeff(0,:)*fun),COS(u),SUM(coeff(1,:)*fun)
    end do
    close(1)

  end subroutine plag_test

  !---------------------------------------------------------------------
  subroutine plag_testperiod(fieldperiod,nlagrange,ndata)
    use nrtype, only : dp

    type(fieldperiod_struct), pointer :: fieldperiod
    integer, intent(in)      :: nlagrange,ndata

    integer, parameter :: nder=1
    integer :: nstep,i
    integer, save :: u1 = 117

    real(kind=dp) :: phi_start,phi_end,phi_span,phi
    real(kind=dp) :: bhat,bhatder
    real(kind=dp) :: x1,x3,geodcu,h_phi,dlogbdphi
    real(kind=dp) :: dbcovar_s_hat_dphi,bcovar_s_hat,dlogbds

    nstep = fieldperiod%parent%parent%nstep
    phi_start = fieldperiod%coords%x2(0)
    phi_end   = fieldperiod%coords%x2(nstep)
    phi_span  = phi_end - phi_start

    open(unit=u1,file='testperiod.dat')
    do i = 0, nstep
      !! Modifications by Andreas F. Martitsch (11.06.2014)
      ! Optional output (necessary for modeling the magnetic rotation)
      if ( allocated(fieldperiod%mdata%dbcovar_s_hat_dphi) .AND. &
           allocated(fieldperiod%mdata%bcovar_s_hat)       .AND. &
           allocated(fieldperiod%mdata%dlogbds) ) then
        write(u1,*) &
             fieldperiod%coords%x1(i),fieldperiod%coords%x2(i), &
             fieldperiod%coords%x3(i), &
             fieldperiod%mdata%bhat(i),fieldperiod%mdata%geodcu(i), &
             fieldperiod%mdata%h_phi(i),fieldperiod%mdata%dlogbdphi(i), &
             fieldperiod%mdata%dbcovar_s_hat_dphi(i), &
             fieldperiod%mdata%bcovar_s_hat(i),fieldperiod%mdata%dlogbds(i)
      else ! This is the old version:
        write(u1,*) &
             fieldperiod%coords%x1(i),fieldperiod%coords%x2(i),&
             fieldperiod%coords%x3(i), &
             fieldperiod%mdata%bhat(i),fieldperiod%mdata%geodcu(i), &
             fieldperiod%mdata%h_phi(i),fieldperiod%mdata%dlogbdphi(i)
      end if
      !! end Modifications by Andreas F. Martitsch (11.06.2014)
    end do
    close(unit=u1)

    open(unit=u1,file='testlagrange.dat')
    do i = 0, ndata
      phi = phi_start + phi_span * DBLE(i) / DBLE(ndata)
      call plagrange_interp(fieldperiod,phi,nlagrange,bhat,bhatder)
      write(u1,*) phi,bhat,bhatder
    end do
    close(unit=u1)
    !print *, 'Derivative'
    !PAUSE

    open(unit=u1,file='testlagall.dat')
    do i = 0, ndata
      phi = phi_start + phi_span * DBLE(i) / DBLE(ndata)
      !! Modifications by Andreas F. Martitsch (11.06.2014)
      ! Optional output (necessary for modeling the magnetic rotation)
      if ( allocated(fieldperiod%mdata%dbcovar_s_hat_dphi) .AND. &
           allocated(fieldperiod%mdata%bcovar_s_hat)       .AND. &
           allocated(fieldperiod%mdata%dlogbds) ) then
        call plagrange_interp(fieldperiod,phi,nlagrange,x1,x3,bhat,geodcu,h_phi,&
             dlogbdphi,dbcovar_s_hat_dphi,bcovar_s_hat,dlogbds)
        write(u1,*) x1,phi,x3,bhat,geodcu,h_phi,dlogbdphi,&
             dbcovar_s_hat_dphi,bcovar_s_hat,dlogbds
      else ! This is the old version:
        call plagrange_interp(fieldperiod,phi,nlagrange,x1,x3,bhat,geodcu,h_phi,dlogbdphi)
        write(u1,*) x1,phi,x3,bhat,geodcu,h_phi,dlogbdphi
      end if
      !! end Modifications by Andreas F. Martitsch (11.06.2014)
    end do
    close(unit=u1)

    print *, '------ plag_testperiod ------'
    print *, 'fieldperiod%tag', fieldperiod%tag
    print *, 'nstep ',nstep

  end subroutine plag_testperiod

  !---------------------------------------------------------------------
  subroutine plag_testperiod_2(fieldperiod)
    use nrtype, only : dp

    type(fieldperiod_struct), pointer :: fieldperiod
    type(fieldpropagator_struct), pointer :: fieldpropagator
    type(fieldripple_struct), pointer :: fieldripple

    integer :: u1 = 117
    integer :: u2 = 118
    integer :: u3 = 119
    integer :: i,count
    integer :: pend_tag
    integer :: pstart_tag

    real(kind=dp) :: r0

    pend_tag = 15
    pstart_tag = -1

    fieldperiod => fieldperiod%parent%ch_fir
    fieldpropagator => fieldperiod%ch_ext
    fieldripple => fieldpropagator%ch_act

    r0 = fieldperiod%parent%parent%parent%r0


    open(unit=u1,file='periods.dat')
    open(unit=u2,file='limits.dat')
    count = 0
    do
      if (fieldperiod%tag .LT. pstart_tag) then
        if (ASSOCIATED(fieldperiod%next)) then
          fieldperiod => fieldperiod%next
        else
          exit
        end if
        cycle
      end if
      do i = count, UBOUND(fieldperiod%coords%x2,1)
        write(u1,'(1000f25.16)') &
             fieldperiod%coords%x1(i),fieldperiod%coords%x2(i),fieldperiod%coords%x3(i), &
             fieldperiod%mdata%bhat(i),fieldperiod%mdata%geodcu(i), &
             fieldperiod%mdata%h_phi(i),fieldperiod%mdata%dlogbdphi(i)
      end do
      count = 1
      write(u2,'(1000f25.16)') fieldperiod%phi_l
      write(u2,'(1000f25.16)') fieldperiod%phi_r
      if (ASSOCIATED(fieldperiod%next)) then
        fieldperiod => fieldperiod%next
      else
        exit
      end if
      if (fieldperiod%tag .GE. pend_tag) exit
    end do
    close(unit=u1)
    close(unit=u2)

    open(unit=u1,file='props.dat')
    open(unit=u2,file='extrema.dat')
    open(unit=u3,file='proplim.dat')
    count = 0
    do
      if (fieldpropagator%parent%tag .LT. pstart_tag) then
        if (ASSOCIATED(fieldpropagator%next)) then
          fieldpropagator => fieldpropagator%next
        else
          exit
        end if
        cycle
      end if
      if (fieldpropagator%parent%tag .GT. pend_tag) exit
      do i = count, UBOUND(fieldpropagator%coords%x2,1)
        write(u1,'(1000f25.16)') &
             fieldpropagator%coords%x1(i),fieldpropagator%coords%x2(i),fieldpropagator%coords%x3(i), &
             fieldpropagator%mdata%bhat(i),fieldpropagator%mdata%geodcu(i), &
             fieldpropagator%mdata%h_phi(i),fieldpropagator%mdata%dlogbdphi(i)
      end do
      count = 1
      write(u2,'(1000f25.16)') fieldpropagator%phi_l,fieldpropagator%b_l
      write(u2,'(1000f25.16)') fieldpropagator%phi_min,fieldpropagator%b_min
      write(u2,'(1000f25.16)') fieldpropagator%phi_r,fieldpropagator%b_r
      write(u3,'(1000f25.16)') fieldpropagator%phi_l
      write(u3,'(1000f25.16)') fieldpropagator%phi_r
      if (ASSOCIATED(fieldpropagator%next)) then
        fieldpropagator => fieldpropagator%next
      else
        exit
      end if
    end do
    close(unit=u1)
    close(unit=u2)
    close(unit=u3)

    open(unit=u1,file='ripple.dat')
    do
      if (fieldripple%pa_fir%parent%tag .LT. pstart_tag) then
        if (ASSOCIATED(fieldripple%next)) then
          fieldripple => fieldripple%next
        else
          exit
        end if
        cycle
      end if
      if (fieldripple%pa_fir%parent%tag .GT. pend_tag) exit
      write(u1,'(1000f25.16)') fieldripple%pa_fir%phi_l,fieldripple%b_max_l
      write(u1,'(1000f25.16)') fieldripple%pa_fir%phi_l + fieldripple%width_l / r0 ,fieldripple%b_min
      write(u1,'(1000f25.16)') fieldripple%pa_las%phi_r,fieldripple%b_max_r
      if (ASSOCIATED(fieldripple%next)) then
        fieldripple => fieldripple%next
      else
        exit
      end if
    end do
    close(unit=u1)

  end subroutine plag_testperiod_2

  !---------------------------------------------------------------------
  subroutine plag_stencel(fieldperiod,phi,nlagrange,stencel)
    use nrtype, only : dp

    type(fieldperiod_struct), pointer                 :: fieldperiod, fp
    real(kind=dp), intent(in)                         :: phi
    integer, intent(in)                               :: nlagrange
    integer, dimension(:), allocatable, intent(inout) :: stencel

    real(kind=dp) :: phi_start,phi_end,phi_span
    integer       :: nstep,mid,i


    !integer :: nstep
    real(kind=dp) ::  phi_l
    real(kind=dp) ::  phi_last,phi_first

    fp => fieldperiod

    nstep = UBOUND(fp%coords%x2,1)
    phi_l = phi
    ! fix phi_l locally to be within range
    phi_first = fp%parent%ch_fir%coords%x2(0)
    phi_last  = fp%parent%ch_las%coords%x2(nstep)
    phi_span  = phi_last - phi_first
    do while (phi_l .GT. phi_last)
      phi_l = phi_l - phi_span
    end do
    do while (phi_l .LT. phi_first)
       phi_l = phi_l + phi_span
    end do
    ! fix period if necessary
    phi_last  = fp%coords%x2(nstep)
    do while (phi_l .GT. phi_last)
      fp => fp%next
      phi_last  = fp%coords%x2(nstep)
    end do
    phi_first = fp%coords%x2(0)
    do while (phi_l .LT. phi_first)
      fp => fp%prev
      phi_first  = fp%coords%x2(0)
    end do

    nstep = fp%parent%parent%nstep
    phi_start = fp%coords%x2(0)
    phi_end   = fp%coords%x2(nstep)
    phi_span  = phi_end - phi_start

    if (allocated(stencel)) deallocate(stencel)
    allocate(stencel(nlagrange+1))

    mid = INT( (phi_l - phi_start) / phi_span * DBLE(nstep) )
    do i = 1,nlagrange+1
      stencel(i) = mid - (nlagrange+1)/2 + i
    end do

  end subroutine plag_stencel

  !---------------------------------------------------------------------
  subroutine plag_value_1(fieldperiod,stencel,phi_arr)
    use nrtype, only : dp

    type(fieldperiod_struct), pointer                       :: fieldperiod
    integer, dimension(:), allocatable, intent(inout)       :: stencel
    real(kind=dp), dimension(:), allocatable, intent(inout) :: phi_arr

    real(kind=dp), dimension(:), allocatable :: hlp

    integer :: nstep

    if (allocated(phi_arr)) deallocate(phi_arr)
    allocate(phi_arr(LBOUND(stencel,1):UBOUND(stencel,1)))

    nstep = fieldperiod%parent%parent%nstep

    if (MINVAL(stencel) .LT. 0) then
      allocate(hlp(-nstep:nstep))
      if (ASSOCIATED(fieldperiod%prev)) then
        hlp(-nstep:0) = fieldperiod%prev%coords%x2(0:nstep)
        hlp(0:nstep)  = fieldperiod%coords%x2(0:nstep)
      else
        hlp(-nstep:0) = fieldperiod%parent%ch_las%coords%x2(0:nstep) &
             - fieldperiod%parent%ch_las%coords%x2(nstep) &
             + fieldperiod%coords%x2(0)
        hlp(0:nstep)  = fieldperiod%coords%x2(0:nstep)
      end if
    elseif (MAXVAL(stencel) .GT. nstep) then
      allocate(hlp(0:2*nstep))
      if (ASSOCIATED(fieldperiod%next)) then
        hlp(0:nstep)  = fieldperiod%coords%x2(0:nstep)
        hlp(nstep:2*nstep) = fieldperiod%next%coords%x2(0:nstep)
      else
        hlp(0:nstep)  = fieldperiod%coords%x2(0:nstep)
        hlp(nstep:2*nstep) = fieldperiod%parent%ch_fir%coords%x2(0:nstep) &
             - fieldperiod%parent%ch_fir%coords%x2(0) &
             + fieldperiod%coords%x2(nstep)
      end if
    else
      allocate(hlp(0:nstep))
      hlp(0:nstep) = fieldperiod%coords%x2(0:nstep)
    end if
    phi_arr = hlp(stencel)
    deallocate(hlp)

  end subroutine plag_value_1

  !---------------------------------------------------------------------
  subroutine plag_value_2(fieldperiod,stencel,phi_arr,bhat_arr)
    use nrtype, only : dp

    type(fieldperiod_struct), pointer                       :: fieldperiod
    integer, dimension(:), allocatable, intent(inout)       :: stencel
    real(kind=dp), dimension(:), allocatable, intent(inout) :: phi_arr
    real(kind=dp), dimension(:), allocatable, intent(inout) :: bhat_arr

    real(kind=dp), dimension(:), allocatable :: hlp

    integer :: nstep

    nstep = fieldperiod%parent%parent%nstep

    ! phi_arr
    if (allocated(phi_arr)) deallocate(phi_arr)
    allocate(phi_arr(LBOUND(stencel,1):UBOUND(stencel,1)))
    if (MINVAL(stencel) .LT. 0) then
      allocate(hlp(-nstep:nstep))
      if (ASSOCIATED(fieldperiod%prev)) then
        hlp(-nstep:0) = fieldperiod%prev%coords%x2(0:nstep)
        hlp(0:nstep)  = fieldperiod%coords%x2(0:nstep)
      else
        hlp(-nstep:0) = fieldperiod%parent%ch_las%coords%x2(0:nstep) &
             - fieldperiod%parent%ch_las%coords%x2(nstep) &
             + fieldperiod%coords%x2(0)
        hlp(0:nstep)  = fieldperiod%coords%x2(0:nstep)
      end if
    elseif (MAXVAL(stencel) .GT. nstep) then
      allocate(hlp(0:2*nstep))
      if (ASSOCIATED(fieldperiod%next)) then
        hlp(0:nstep)  = fieldperiod%coords%x2(0:nstep)
        hlp(nstep:2*nstep) = fieldperiod%next%coords%x2(0:nstep)
      else
        hlp(0:nstep)  = fieldperiod%coords%x2(0:nstep)
        hlp(nstep:2*nstep) = fieldperiod%parent%ch_fir%coords%x2(0:nstep) &
             - fieldperiod%parent%ch_fir%coords%x2(0) &
             + fieldperiod%coords%x2(nstep)
      end if
    else
      allocate(hlp(0:nstep))
      hlp(0:nstep) = fieldperiod%coords%x2(0:nstep)
    end if
    phi_arr = hlp(stencel)
    deallocate(hlp)

    ! bhat_arr
    if (allocated(bhat_arr)) deallocate(bhat_arr)
    allocate(bhat_arr(LBOUND(stencel,1):UBOUND(stencel,1)))
    if (MINVAL(stencel) .LT. 0) then
      allocate(hlp(-nstep:nstep))
      if (ASSOCIATED(fieldperiod%prev)) then
        hlp(-nstep:0) = fieldperiod%prev%mdata%bhat(0:nstep)
        hlp(0:nstep)  = fieldperiod%mdata%bhat(0:nstep)
      else
        hlp(-nstep:0) = fieldperiod%parent%ch_las%mdata%bhat(0:nstep)
        hlp(0:nstep)  = fieldperiod%mdata%bhat(0:nstep)
      end if
    elseif (MAXVAL(stencel) .GT. nstep) then
      allocate(hlp(0:2*nstep))
      if (ASSOCIATED(fieldperiod%next)) then
        hlp(0:nstep)  = fieldperiod%mdata%bhat(0:nstep)
        hlp(nstep:2*nstep) = fieldperiod%next%mdata%bhat(0:nstep)
      else
        hlp(0:nstep)  = fieldperiod%mdata%bhat(0:nstep)
        hlp(nstep:2*nstep) = fieldperiod%parent%ch_fir%mdata%bhat(0:nstep)
      end if
    else
      allocate(hlp(0:nstep))
      hlp(0:nstep) = fieldperiod%mdata%bhat(0:nstep)
    end if
    bhat_arr = hlp(stencel)
    deallocate(hlp)
    ! phi_arr, coords, x2
  end subroutine plag_value_2

  !---------------------------------------------------------------------
  subroutine plag_value_all(fieldperiod,stencel,phi_arr, &
       x1_arr,x3_arr,bhat_arr,geodcu_arr,h_phi_arr,dlogbdphi_arr)
    use nrtype, only : dp

    type(fieldperiod_struct), pointer                       :: fieldperiod
    integer, dimension(:), allocatable, intent(inout)       :: stencel
    real(kind=dp), dimension(:), allocatable, intent(inout) :: phi_arr
    real(kind=dp), dimension(:), allocatable, intent(inout) :: x1_arr,x3_arr,bhat_arr
    real(kind=dp), dimension(:), allocatable, intent(inout) :: geodcu_arr,h_phi_arr,dlogbdphi_arr

    real(kind=dp), dimension(:), allocatable :: hlp

    integer :: nstep

    nstep = fieldperiod%parent%parent%nstep

    ! phi_arr
    if (allocated(phi_arr)) deallocate(phi_arr)
    allocate(phi_arr(LBOUND(stencel,1):UBOUND(stencel,1)))
    if (MINVAL(stencel) .LT. 0) then
      allocate(hlp(-nstep:nstep))
      if (ASSOCIATED(fieldperiod%prev)) then
        hlp(-nstep:0) = fieldperiod%prev%coords%x2(0:nstep)
        hlp(0:nstep)  = fieldperiod%coords%x2(0:nstep)
      else
        hlp(-nstep:0) = fieldperiod%parent%ch_las%coords%x2(0:nstep) &
             - fieldperiod%parent%ch_las%coords%x2(nstep) &
             + fieldperiod%coords%x2(0)
        hlp(0:nstep)  = fieldperiod%coords%x2(0:nstep)
      end if
    elseif (MAXVAL(stencel) .GT. nstep) then
      allocate(hlp(0:2*nstep))
      if (ASSOCIATED(fieldperiod%next)) then
        hlp(0:nstep)  = fieldperiod%coords%x2(0:nstep)
        hlp(nstep:2*nstep) = fieldperiod%next%coords%x2(0:nstep)
      else
        hlp(0:nstep)  = fieldperiod%coords%x2(0:nstep)
        hlp(nstep:2*nstep) = fieldperiod%parent%ch_fir%coords%x2(0:nstep) &
             - fieldperiod%parent%ch_fir%coords%x2(0) &
             + fieldperiod%coords%x2(nstep)
      end if
    else
      allocate(hlp(0:nstep))
      hlp(0:nstep) = fieldperiod%coords%x2(0:nstep)
    end if
    phi_arr = hlp(stencel)
    deallocate(hlp)

    ! x1_arr
    if (allocated(x1_arr)) deallocate(x1_arr)
    allocate(x1_arr(LBOUND(stencel,1):UBOUND(stencel,1)))
    if (MINVAL(stencel) .LT. 0) then
      allocate(hlp(-nstep:nstep))
      if (ASSOCIATED(fieldperiod%prev)) then
        hlp(-nstep:0) = fieldperiod%prev%coords%x1(0:nstep)
        hlp(0:nstep)  = fieldperiod%coords%x1(0:nstep)
      else
        hlp(-nstep:0) = fieldperiod%parent%ch_las%coords%x1(0:nstep)
        hlp(0:nstep)  = fieldperiod%coords%x1(0:nstep)
      end if
    elseif (MAXVAL(stencel) .GT. nstep) then
      allocate(hlp(0:2*nstep))
      if (ASSOCIATED(fieldperiod%next)) then
        hlp(0:nstep)  = fieldperiod%coords%x1(0:nstep)
        hlp(nstep:2*nstep) = fieldperiod%next%coords%x1(0:nstep)
      else
        hlp(0:nstep)  = fieldperiod%coords%x1(0:nstep)
        hlp(nstep:2*nstep) = fieldperiod%parent%ch_fir%coords%x1(0:nstep)
      end if
    else
      allocate(hlp(0:nstep))
      hlp(0:nstep) = fieldperiod%coords%x1(0:nstep)
    end if
    x1_arr = hlp(stencel)
    deallocate(hlp)
    ! x1

    ! x3_arr
    if (allocated(x3_arr)) deallocate(x3_arr)
    allocate(x3_arr(LBOUND(stencel,1):UBOUND(stencel,1)))
    if (MINVAL(stencel) .LT. 0) then
      allocate(hlp(-nstep:nstep))
      if (ASSOCIATED(fieldperiod%prev)) then
        hlp(-nstep:0) = fieldperiod%prev%coords%x3(0:nstep)
        hlp(0:nstep)  = fieldperiod%coords%x3(0:nstep)
      else
        hlp(-nstep:0) = fieldperiod%parent%ch_las%coords%x3(0:nstep)
        hlp(0:nstep)  = fieldperiod%coords%x3(0:nstep)
      end if
    elseif (MAXVAL(stencel) .GT. nstep) then
      allocate(hlp(0:2*nstep))
      if (ASSOCIATED(fieldperiod%next)) then
        hlp(0:nstep)  = fieldperiod%coords%x3(0:nstep)
        hlp(nstep:2*nstep) = fieldperiod%next%coords%x3(0:nstep)
      else
        hlp(0:nstep)  = fieldperiod%coords%x3(0:nstep)
        hlp(nstep:2*nstep) = fieldperiod%parent%ch_fir%coords%x3(0:nstep)
      end if
    else
      allocate(hlp(0:nstep))
      hlp(0:nstep) = fieldperiod%coords%x3(0:nstep)
    end if
    x3_arr = hlp(stencel)
    deallocate(hlp)
    ! x3

    ! bhat_arr
    if (allocated(bhat_arr)) deallocate(bhat_arr)
    allocate(bhat_arr(LBOUND(stencel,1):UBOUND(stencel,1)))
    if (MINVAL(stencel) .LT. 0) then
      allocate(hlp(-nstep:nstep))
      if (ASSOCIATED(fieldperiod%prev)) then
        hlp(-nstep:0) = fieldperiod%prev%mdata%bhat(0:nstep)
        hlp(0:nstep)  = fieldperiod%mdata%bhat(0:nstep)
      else
        hlp(-nstep:0) = fieldperiod%parent%ch_las%mdata%bhat(0:nstep)
        hlp(0:nstep)  = fieldperiod%mdata%bhat(0:nstep)
      end if
    elseif (MAXVAL(stencel) .GT. nstep) then
      allocate(hlp(0:2*nstep))
      if (ASSOCIATED(fieldperiod%next)) then
        hlp(0:nstep)  = fieldperiod%mdata%bhat(0:nstep)
        hlp(nstep:2*nstep) = fieldperiod%next%mdata%bhat(0:nstep)
      else
        hlp(0:nstep)  = fieldperiod%mdata%bhat(0:nstep)
        hlp(nstep:2*nstep) = fieldperiod%parent%ch_fir%mdata%bhat(0:nstep)
      end if
    else
      allocate(hlp(0:nstep))
      hlp(0:nstep) = fieldperiod%mdata%bhat(0:nstep)
    end if
    bhat_arr = hlp(stencel)
    deallocate(hlp)
    ! bhat_arr

    ! geodcu_arr
    if (allocated(geodcu_arr)) deallocate(geodcu_arr)
    allocate(geodcu_arr(LBOUND(stencel,1):UBOUND(stencel,1)))
    if (MINVAL(stencel) .LT. 0) then
      allocate(hlp(-nstep:nstep))
      if (ASSOCIATED(fieldperiod%prev)) then
        hlp(-nstep:0) = fieldperiod%prev%mdata%geodcu(0:nstep)
        hlp(0:nstep)  = fieldperiod%mdata%geodcu(0:nstep)
      else
        hlp(-nstep:0) = fieldperiod%parent%ch_las%mdata%geodcu(0:nstep)
        hlp(0:nstep)  = fieldperiod%mdata%geodcu(0:nstep)
      end if
    elseif (MAXVAL(stencel) .GT. nstep) then
      allocate(hlp(0:2*nstep))
      if (ASSOCIATED(fieldperiod%next)) then
        hlp(0:nstep)  = fieldperiod%mdata%geodcu(0:nstep)
        hlp(nstep:2*nstep) = fieldperiod%next%mdata%geodcu(0:nstep)
      else
        hlp(0:nstep)  = fieldperiod%mdata%geodcu(0:nstep)
        hlp(nstep:2*nstep) = fieldperiod%parent%ch_fir%mdata%geodcu(0:nstep)
      end if
    else
      allocate(hlp(0:nstep))
      hlp(0:nstep) = fieldperiod%mdata%geodcu(0:nstep)
    end if
    geodcu_arr = hlp(stencel)
    deallocate(hlp)
    ! geodcu_arr

    ! h_phi_arr
    if (allocated(h_phi_arr)) deallocate(h_phi_arr)
    allocate(h_phi_arr(LBOUND(stencel,1):UBOUND(stencel,1)))
    if (MINVAL(stencel) .LT. 0) then
      allocate(hlp(-nstep:nstep))
      if (ASSOCIATED(fieldperiod%prev)) then
        hlp(-nstep:0) = fieldperiod%prev%mdata%h_phi(0:nstep)
        hlp(0:nstep)  = fieldperiod%mdata%h_phi(0:nstep)
      else
        hlp(-nstep:0) = fieldperiod%parent%ch_las%mdata%h_phi(0:nstep)
        hlp(0:nstep)  = fieldperiod%mdata%h_phi(0:nstep)
      end if
    elseif (MAXVAL(stencel) .GT. nstep) then
      allocate(hlp(0:2*nstep))
      if (ASSOCIATED(fieldperiod%next)) then
        hlp(0:nstep)  = fieldperiod%mdata%h_phi(0:nstep)
        hlp(nstep:2*nstep) = fieldperiod%next%mdata%h_phi(0:nstep)
      else
        hlp(0:nstep)  = fieldperiod%mdata%h_phi(0:nstep)
        hlp(nstep:2*nstep) = fieldperiod%parent%ch_fir%mdata%h_phi(0:nstep)
      end if
    else
      allocate(hlp(0:nstep))
      hlp(0:nstep) = fieldperiod%mdata%h_phi(0:nstep)
    end if
    h_phi_arr = hlp(stencel)
    deallocate(hlp)
    ! h_phi_arr

    ! dlogbdphi_arr
    if (allocated(dlogbdphi_arr)) deallocate(dlogbdphi_arr)
    allocate(dlogbdphi_arr(LBOUND(stencel,1):UBOUND(stencel,1)))
    if (MINVAL(stencel) .LT. 0) then
      allocate(hlp(-nstep:nstep))
      if (ASSOCIATED(fieldperiod%prev)) then
        hlp(-nstep:0) = fieldperiod%prev%mdata%dlogbdphi(0:nstep)
        hlp(0:nstep)  = fieldperiod%mdata%dlogbdphi(0:nstep)
      else
        hlp(-nstep:0) = fieldperiod%parent%ch_las%mdata%dlogbdphi(0:nstep)
        hlp(0:nstep)  = fieldperiod%mdata%dlogbdphi(0:nstep)
      end if
    elseif (MAXVAL(stencel) .GT. nstep) then
      allocate(hlp(0:2*nstep))
      if (ASSOCIATED(fieldperiod%next)) then
        hlp(0:nstep)  = fieldperiod%mdata%dlogbdphi(0:nstep)
        hlp(nstep:2*nstep) = fieldperiod%next%mdata%dlogbdphi(0:nstep)
      else
        hlp(0:nstep)  = fieldperiod%mdata%dlogbdphi(0:nstep)
        hlp(nstep:2*nstep) = fieldperiod%parent%ch_fir%mdata%dlogbdphi(0:nstep)
      end if
    else
      allocate(hlp(0:nstep))
      hlp(0:nstep) = fieldperiod%mdata%dlogbdphi(0:nstep)
    end if
    dlogbdphi_arr = hlp(stencel)
    deallocate(hlp)
    ! dlogbdphi_arr

  end subroutine plag_value_all

  !---------------------------------------------------------------------
  !! Modifications by Andreas F. Martitsch (13.03.2014)
  ! Optional output (necessary for modeling the magnetic rotation)
  subroutine plag_value_all2(fieldperiod,stencel,phi_arr, &
       x1_arr,x3_arr,bhat_arr,geodcu_arr,h_phi_arr,dlogbdphi_arr,&
       dbcovar_s_hat_dphi_arr,bcovar_s_hat_arr,dlogbds_arr)
    use nrtype, only : dp

    type(fieldperiod_struct), pointer                       :: fieldperiod
    integer, dimension(:), allocatable, intent(inout)       :: stencel
    real(kind=dp), dimension(:), allocatable, intent(inout) :: phi_arr
    real(kind=dp), dimension(:), allocatable, intent(inout) :: x1_arr,x3_arr,bhat_arr
    real(kind=dp), dimension(:), allocatable, intent(inout) :: geodcu_arr,h_phi_arr
    real(kind=dp), dimension(:), allocatable, intent(inout) :: dlogbdphi_arr
    real(kind=dp), dimension(:), allocatable, intent(inout) :: dbcovar_s_hat_dphi_arr
    real(kind=dp), dimension(:), allocatable, intent(inout) :: bcovar_s_hat_arr
    real(kind=dp), dimension(:), allocatable, intent(inout) :: dlogbds_arr

    real(kind=dp), dimension(:), allocatable :: hlp

    integer :: nstep

    nstep = fieldperiod%parent%parent%nstep

    ! Compute phi_arr, x1_arr, x3_arr, bhat_arr, geodcu_arr, h_phi_arr, dlogbdphi_arr,
    ! dbcovar_s_hat_dphi_arr, bcovar_s_hat_arr, dlogbds_arr
    call plag_value_all(fieldperiod,stencel,phi_arr, &
         x1_arr,x3_arr,bhat_arr,geodcu_arr,h_phi_arr,dlogbdphi_arr)

    ! dbcovar_s_hat_dphi_arr
    if (allocated(dbcovar_s_hat_dphi_arr)) deallocate(dbcovar_s_hat_dphi_arr)
    allocate(dbcovar_s_hat_dphi_arr(LBOUND(stencel,1):UBOUND(stencel,1)))
    if (MINVAL(stencel) .LT. 0) then
      allocate(hlp(-nstep:nstep))
      if (ASSOCIATED(fieldperiod%prev)) then
        hlp(-nstep:0) = fieldperiod%prev%mdata%dbcovar_s_hat_dphi(0:nstep)
        hlp(0:nstep)  = fieldperiod%mdata%dbcovar_s_hat_dphi(0:nstep)
      else
        hlp(-nstep:0) = fieldperiod%parent%ch_las%mdata%dbcovar_s_hat_dphi(0:nstep)
        hlp(0:nstep)  = fieldperiod%mdata%dbcovar_s_hat_dphi(0:nstep)
      end if
    elseif (MAXVAL(stencel) .GT. nstep) then
      allocate(hlp(0:2*nstep))
      if (ASSOCIATED(fieldperiod%next)) then
        hlp(0:nstep)  = fieldperiod%mdata%dbcovar_s_hat_dphi(0:nstep)
        hlp(nstep:2*nstep) = fieldperiod%next%mdata%dbcovar_s_hat_dphi(0:nstep)
      else
        hlp(0:nstep)  = fieldperiod%mdata%dbcovar_s_hat_dphi(0:nstep)
        hlp(nstep:2*nstep) = fieldperiod%parent%ch_fir%mdata%dbcovar_s_hat_dphi(0:nstep)
      end if
    else
      allocate(hlp(0:nstep))
      hlp(0:nstep) = fieldperiod%mdata%dbcovar_s_hat_dphi(0:nstep)
    end if
    dbcovar_s_hat_dphi_arr = hlp(stencel)
    deallocate(hlp)
    ! dbcovar_s_hat_dphi_arr

    ! bcovar_s_hat_arr
    if (allocated(bcovar_s_hat_arr)) deallocate(bcovar_s_hat_arr)
    allocate(bcovar_s_hat_arr(LBOUND(stencel,1):UBOUND(stencel,1)))
    if (MINVAL(stencel) .LT. 0) then
      allocate(hlp(-nstep:nstep))
      if (ASSOCIATED(fieldperiod%prev)) then
        hlp(-nstep:0) = fieldperiod%prev%mdata%bcovar_s_hat(0:nstep)
        hlp(0:nstep)  = fieldperiod%mdata%bcovar_s_hat(0:nstep)
      else
        hlp(-nstep:0) = fieldperiod%parent%ch_las%mdata%bcovar_s_hat(0:nstep)
        hlp(0:nstep)  = fieldperiod%mdata%bcovar_s_hat(0:nstep)
      end if
    elseif (MAXVAL(stencel) .GT. nstep) then
      allocate(hlp(0:2*nstep))
      if (ASSOCIATED(fieldperiod%next)) then
        hlp(0:nstep)  = fieldperiod%mdata%bcovar_s_hat(0:nstep)
        hlp(nstep:2*nstep) = fieldperiod%next%mdata%bcovar_s_hat(0:nstep)
      else
        hlp(0:nstep)  = fieldperiod%mdata%bcovar_s_hat(0:nstep)
        hlp(nstep:2*nstep) = fieldperiod%parent%ch_fir%mdata%bcovar_s_hat(0:nstep)
      end if
    else
      allocate(hlp(0:nstep))
      hlp(0:nstep) = fieldperiod%mdata%bcovar_s_hat(0:nstep)
    end if
    bcovar_s_hat_arr = hlp(stencel)
    deallocate(hlp)
    ! bcovar_s_hat_arr

    ! dlogbds_arr
    if (allocated(dlogbds_arr)) deallocate(dlogbds_arr)
    allocate(dlogbds_arr(LBOUND(stencel,1):UBOUND(stencel,1)))
    if (MINVAL(stencel) .LT. 0) then
      allocate(hlp(-nstep:nstep))
      if (ASSOCIATED(fieldperiod%prev)) then
        hlp(-nstep:0) = fieldperiod%prev%mdata%dlogbds(0:nstep)
        hlp(0:nstep)  = fieldperiod%mdata%dlogbds(0:nstep)
      else
        hlp(-nstep:0) = fieldperiod%parent%ch_las%mdata%dlogbds(0:nstep)
        hlp(0:nstep)  = fieldperiod%mdata%dlogbds(0:nstep)
      end if
    elseif (MAXVAL(stencel) .GT. nstep) then
      allocate(hlp(0:2*nstep))
      if (ASSOCIATED(fieldperiod%next)) then
        hlp(0:nstep)  = fieldperiod%mdata%dlogbds(0:nstep)
        hlp(nstep:2*nstep) = fieldperiod%next%mdata%dlogbds(0:nstep)
      else
        hlp(0:nstep)  = fieldperiod%mdata%dlogbds(0:nstep)
        hlp(nstep:2*nstep) = fieldperiod%parent%ch_fir%mdata%dlogbds(0:nstep)
      end if
    else
      allocate(hlp(0:nstep))
      hlp(0:nstep) = fieldperiod%mdata%dlogbds(0:nstep)
    end if
    dlogbds_arr = hlp(stencel)
    deallocate(hlp)
    ! dlogbds_arr

  end subroutine plag_value_all2
  !! end Modifications by Andreas F. Martitsch (13.03.2014)

  !---------------------------------------------------------------------
  subroutine plag_interp_2(fieldperiod,phi,nlagrange,bhat,bhatder)
    use nrtype, only : dp

    type(fieldperiod_struct), pointer                 :: fieldperiod, fp
    real(kind=dp), intent(in)                         :: phi
    integer, intent(in)                               :: nlagrange
    real(kind=dp), intent(out)                        :: bhat,bhatder

    integer, parameter :: nder = 1
    integer, dimension(:), allocatable :: stencel
    real(kind=dp), dimension(:), allocatable :: phi_arr,bhat_arr
    real(kind=dp), dimension(:,:), allocatable :: coeff


    integer :: nstep
    real(kind=dp) ::  phi_l
    real(kind=dp) ::  phi_span,phi_last,phi_first

    fp => fieldperiod

    nstep = UBOUND(fp%coords%x2,1)
    phi_l = phi
    ! fix phi_l locally to be within range
    phi_first = fp%parent%ch_fir%coords%x2(0)
    phi_last  = fp%parent%ch_las%coords%x2(nstep)
    phi_span  = phi_last - phi_first
    do while (phi_l .GT. phi_last)
      phi_l = phi_l - phi_span
    end do
    do while (phi_l .LT. phi_first)
      phi_l = phi_l + phi_span
    end do
    ! fix period if necessary
    phi_last  = fp%coords%x2(nstep)
    do while (phi_l .GT. phi_last)
      fp => fp%next
      phi_last  = fp%coords%x2(nstep)
    end do
    phi_first = fp%coords%x2(0)
    do while (phi_l .LT. phi_first)
      fp => fp%prev
      phi_first  = fp%coords%x2(0)
    end do

    call plagrange_stencel(fp,phi_l,nlagrange,stencel)
    call plagrange_value(fp,stencel,phi_arr,bhat_arr)
    allocate( coeff(0:nder,nlagrange+1) )
    call plagrange_coeff(nlagrange+1,nder,phi_l,phi_arr(:),coeff(:,:))
    bhat = SUM(coeff(0,:)*bhat_arr)
    bhatder = SUM(coeff(1,:)*bhat_arr)

    if (allocated(stencel)) deallocate(stencel)
    if (allocated(phi_arr)) deallocate(phi_arr)
    if (allocated(bhat_arr)) deallocate(bhat_arr)
    if (allocated(coeff)) deallocate(coeff)

  end subroutine plag_interp_2

  !---------------------------------------------------------------------
  subroutine plag_interp_3(fieldperiod,phi,nlagrange,bhat,bhatder,bhatdder)
    use nrtype, only : dp

    type(fieldperiod_struct), pointer                 :: fieldperiod, fp
    real(kind=dp), intent(in)                         :: phi
    integer, intent(in)                               :: nlagrange
    real(kind=dp), intent(out)                        :: bhat,bhatder,bhatdder

    integer, parameter :: nder = 2
    integer, dimension(:), allocatable :: stencel
    real(kind=dp), dimension(:), allocatable :: phi_arr,bhat_arr
    real(kind=dp), dimension(:,:), allocatable :: coeff


    integer :: nstep
    real(kind=dp) ::  phi_l
    real(kind=dp) ::  phi_span,phi_last,phi_first

    fp => fieldperiod

    nstep = UBOUND(fp%coords%x2,1)
    phi_l = phi
    ! fix phi_l locally to be within range
    phi_first = fp%parent%ch_fir%coords%x2(0)
    phi_last  = fp%parent%ch_las%coords%x2(nstep)
    phi_span  = phi_last - phi_first
    do while (phi_l .GT. phi_last)
      phi_l = phi_l - phi_span
    end do
    do while (phi_l .LT. phi_first)
      phi_l = phi_l + phi_span
    end do
    ! fix period if necessary
    phi_last  = fp%coords%x2(nstep)
    do while (phi_l .GT. phi_last)
      fp => fp%next
      phi_last  = fp%coords%x2(nstep)
    end do
    phi_first = fp%coords%x2(0)
    do while (phi_l .LT. phi_first)
      fp => fp%prev
      phi_first  = fp%coords%x2(0)
    end do

    call plagrange_stencel(fp,phi_l,nlagrange,stencel)
    call plagrange_value(fp,stencel,phi_arr,bhat_arr)
    allocate( coeff(0:nder,nlagrange+1) )
    call plagrange_coeff(nlagrange+1,nder,phi_l,phi_arr(:),coeff(:,:))
    bhat = SUM(coeff(0,:)*bhat_arr)
    bhatder = SUM(coeff(1,:)*bhat_arr)
    bhatdder = SUM(coeff(2,:)*bhat_arr)

    if (allocated(stencel)) deallocate(stencel)
    if (allocated(phi_arr)) deallocate(phi_arr)
    if (allocated(bhat_arr)) deallocate(bhat_arr)
    if (allocated(coeff)) deallocate(coeff)

  end subroutine plag_interp_3

  !---------------------------------------------------------------------
  subroutine plag_interp_all(fieldperiod,phi,nlagrange,x1,x3,bhat,geodcu,h_phi,dlogbdphi)
    use nrtype, only : dp

    type(fieldperiod_struct), pointer                 :: fieldperiod, fp
    real(kind=dp), intent(in)                         :: phi
    integer, intent(in)                               :: nlagrange
    real(kind=dp), intent(out)                        :: x1,x3,bhat,geodcu,h_phi,dlogbdphi

    integer, parameter :: nder = 0
    integer, dimension(:), allocatable :: stencel
    real(kind=dp), dimension(:), allocatable :: phi_arr
    real(kind=dp), dimension(:), allocatable :: x1_arr,x3_arr,bhat_arr
    real(kind=dp), dimension(:), allocatable :: geodcu_arr,h_phi_arr,dlogbdphi_arr
    real(kind=dp), dimension(:,:), allocatable :: coeff

    integer :: nstep
    real(kind=dp) ::  phi_l
    real(kind=dp) ::  phi_span,phi_last,phi_first
    fp => fieldperiod
    nstep = UBOUND(fp%coords%x2,1)
    phi_l = phi
    ! fix phi_l locally to be within range
    phi_first = fp%parent%ch_fir%coords%x2(0)
    phi_last  = fp%parent%ch_las%coords%x2(nstep)
    phi_span  = phi_last - phi_first
    do while (phi_l .GT. phi_last)
      phi_l = phi_l - phi_span
    end do
    do while (phi_l .LT. phi_first)
      phi_l = phi_l + phi_span
    end do
    ! fix period if necessary
    phi_last  = fp%coords%x2(nstep)
    do while (phi_l .GT. phi_last)
      fp => fp%next
      phi_last  = fp%coords%x2(nstep)
    end do
    phi_first = fp%coords%x2(0)
    do while (phi_l .LT. phi_first)
      fp => fp%prev
      phi_first  = fp%coords%x2(0)
    end do


    call plagrange_stencel(fp,phi_l,nlagrange,stencel)
    call plagrange_value(fp,stencel,phi_arr, &
         x1_arr,x3_arr,bhat_arr,geodcu_arr,h_phi_arr,dlogbdphi_arr)
    allocate( coeff(0:nder,nlagrange+1) )
    call plagrange_coeff(nlagrange+1,nder,phi_l,phi_arr(:),coeff(:,:))
    x1   = SUM(coeff(0,:)*x1_arr)
    x3   = SUM(coeff(0,:)*x3_arr)
    bhat = SUM(coeff(0,:)*bhat_arr)
    geodcu = SUM(coeff(0,:)*geodcu_arr)
    h_phi = SUM(coeff(0,:)*h_phi_arr)
    dlogbdphi = SUM(coeff(0,:)*dlogbdphi_arr)

    if (allocated(stencel)) deallocate(stencel)
    if (allocated(phi_arr)) deallocate(phi_arr)
    if (allocated(x1_arr)) deallocate(x1_arr)
    if (allocated(x3_arr)) deallocate(x3_arr)
    if (allocated(bhat_arr)) deallocate(bhat_arr)
    if (allocated(geodcu_arr)) deallocate(geodcu_arr)
    if (allocated(h_phi_arr)) deallocate(h_phi_arr)
    if (allocated(dlogbdphi_arr)) deallocate(dlogbdphi_arr)
    if (allocated(coeff)) deallocate(coeff)
  end subroutine plag_interp_all

  !---------------------------------------------------------------------
  !! Modifications by Andreas F. Martitsch (13.03.2014)
  ! Optional output (necessary for modeling the magnetic rotation)
  subroutine plag_interp_all2(fieldperiod,phi,nlagrange,x1,x3,bhat,geodcu,&
       h_phi,dlogbdphi,dbcovar_s_hat_dphi,bcovar_s_hat,dlogbds)
    use nrtype, only : dp

    type(fieldperiod_struct), pointer                 :: fieldperiod, fp
    real(kind=dp), intent(in)                         :: phi
    integer, intent(in)                               :: nlagrange
    real(kind=dp), intent(out)                        :: x1,x3,bhat,geodcu,h_phi,dlogbdphi
    real(kind=dp), intent(out)                        :: dbcovar_s_hat_dphi,bcovar_s_hat
    real(kind=dp), intent(out)                        :: dlogbds

    integer, parameter :: nder = 0
    integer, dimension(:), allocatable :: stencel
    real(kind=dp), dimension(:), allocatable :: phi_arr
    real(kind=dp), dimension(:), allocatable :: x1_arr,x3_arr,bhat_arr
    real(kind=dp), dimension(:), allocatable :: geodcu_arr,h_phi_arr,dlogbdphi_arr
    real(kind=dp), dimension(:), allocatable :: dbcovar_s_hat_dphi_arr,bcovar_s_hat_arr
    real(kind=dp), dimension(:), allocatable :: dlogbds_arr
    real(kind=dp), dimension(:,:), allocatable :: coeff

    integer :: nstep
    real(kind=dp) ::  phi_l
    real(kind=dp) ::  phi_span,phi_last,phi_first
    fp => fieldperiod
    nstep = UBOUND(fp%coords%x2,1)
    phi_l = phi
    ! fix phi_l locally to be within range
    phi_first = fp%parent%ch_fir%coords%x2(0)
    phi_last  = fp%parent%ch_las%coords%x2(nstep)
    phi_span  = phi_last - phi_first
    do while (phi_l .GT. phi_last)
      phi_l = phi_l - phi_span
    end do
    do while (phi_l .LT. phi_first)
      phi_l = phi_l + phi_span
    end do
    ! fix period if necessary
    phi_last  = fp%coords%x2(nstep)
    do while (phi_l .GT. phi_last)
      fp => fp%next
      phi_last  = fp%coords%x2(nstep)
    end do
    phi_first = fp%coords%x2(0)
    do while (phi_l .LT. phi_first)
      fp => fp%prev
      phi_first  = fp%coords%x2(0)
    end do


    call plagrange_stencel(fp,phi_l,nlagrange,stencel)
    call plagrange_value(fp,stencel,phi_arr, &
         x1_arr,x3_arr,bhat_arr,geodcu_arr,h_phi_arr,dlogbdphi_arr,&
         dbcovar_s_hat_dphi_arr,bcovar_s_hat_arr,dlogbds_arr)
    allocate( coeff(0:nder,nlagrange+1) )
    call plagrange_coeff(nlagrange+1,nder,phi_l,phi_arr(:),coeff(:,:))
    x1                 = SUM(coeff(0,:)*x1_arr)
    x3                 = SUM(coeff(0,:)*x3_arr)
    bhat               = SUM(coeff(0,:)*bhat_arr)
    geodcu             = SUM(coeff(0,:)*geodcu_arr)
    h_phi              = SUM(coeff(0,:)*h_phi_arr)
    dlogbdphi          = SUM(coeff(0,:)*dlogbdphi_arr)
    dbcovar_s_hat_dphi = SUM(coeff(0,:)*dbcovar_s_hat_dphi_arr)
    bcovar_s_hat       = SUM(coeff(0,:)*bcovar_s_hat_arr)
    dlogbds            = SUM(coeff(0,:)*dlogbds_arr)

    if (allocated(stencel)) deallocate(stencel)
    if (allocated(phi_arr)) deallocate(phi_arr)
    if (allocated(x1_arr)) deallocate(x1_arr)
    if (allocated(x3_arr)) deallocate(x3_arr)
    if (allocated(bhat_arr)) deallocate(bhat_arr)
    if (allocated(geodcu_arr)) deallocate(geodcu_arr)
    if (allocated(h_phi_arr)) deallocate(h_phi_arr)
    if (allocated(dlogbdphi_arr)) deallocate(dlogbdphi_arr)
    if (allocated(dbcovar_s_hat_dphi_arr)) deallocate(dbcovar_s_hat_dphi_arr)
    if (allocated(bcovar_s_hat_arr)) deallocate(bcovar_s_hat_arr)
    if (allocated(dlogbds_arr)) deallocate(dlogbds_arr)
    if (allocated(coeff)) deallocate(coeff)
  end subroutine plag_interp_all2
  !! end Modifications by Andreas F. Martitsch (13.03.2014)
  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
end module plagrange_mod
