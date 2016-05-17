module collop_laguerre
  use nrtype, only : dp, pi

  implicit none
  
  real(kind=dp), dimension(:),     allocatable :: phi_hm
  real(kind=dp), dimension(:,:,:), allocatable :: gencoeflag, phitrans

contains

  !**********************************************************
  ! Laguerre specific - BEGIN
  !**********************************************************
  subroutine init_laguerre(lagmax, legmax)
    integer :: lagmax, legmax
    integer :: i, j, k, l
    real(kind=dp) :: bincoef

    if (allocated(gencoeflag)) return
    write (*,*) "Initializing Laguerre coefficients..."
    allocate(gencoeflag(0:lagmax, 0:lagmax, 0:legmax))
    gencoeflag=0.d0
    gencoeflag(0,0,:)=1.d0

    do l=1, legmax
       do i=1, lagmax
          do k=0, i
             bincoef=1.d0
             do j=1,i-k
                bincoef=bincoef*(1.d0+(k+l+0.5d0)/j)
             enddo
             do j=1,k
                bincoef=-bincoef/j
             enddo
             gencoeflag(i,k,l)=bincoef
          enddo
       enddo
    enddo
  end subroutine init_laguerre

  subroutine init_phi_laguerre(lagmax, legmax)
    integer :: lagmax, legmax
    integer :: m

    call init_laguerre(lagmax, legmax)

    if (allocated(phi_hm)) deallocate(phi_hm)
    allocate(phi_hm(0:lagmax))

    do m = 0, lagmax
       phi_hm(m) = 1d0/sqrt(gamma(m + 2.5d0) / (2d0 * gamma(m+1d0)))
    end do
    
  end subroutine init_phi_laguerre

  function phi_laguerre(m, x)
    integer       :: m, k
    real(kind=dp) :: phi_laguerre
    real(kind=dp) :: x, plag, xpow, add
    real(kind=dp) :: x2

    x2 = x**2
    plag=gencoeflag(m,0,1)
    xpow=1.d0

    do k=1,m
       add=gencoeflag(m,k,1)*xpow
       plag=plag+add*x2
       xpow=xpow*x2
    enddo

    phi_laguerre = plag * pi**(3d0/4d0) * phi_hm(m)
  end function phi_laguerre

  function d_phi_laguerre(m, x)
    integer       :: m, k
    real(kind=dp) :: d_phi_laguerre
    real(kind=dp) :: x, dplag, xpow, add
    real(kind=dp) :: x2

    x2 = x**2
    dplag=0.d0
    xpow=1.d0

    do k=1,m
       add=gencoeflag(m,k,1)*xpow
       dplag=dplag+k*add
       xpow=xpow*x2
    enddo

    d_phi_laguerre = 2d0 * x * dplag * pi**(3d0/4d0) * phi_hm(m)

    !write (*,*) x2, dplag, xpow, add, phi_hm(m)
  end function d_phi_laguerre

  function dd_phi_laguerre(m, x)
    integer :: m,k
    real(kind=dp) :: dd_phi_laguerre
    real(kind=dp) :: x, dplag, ddplag, xpow,add
    real(kind=dp) :: x2

    x2 = x**2
    dplag=0.d0
    ddplag=0d0
    xpow=1.d0

    do k=1,m
       add=gencoeflag(m,k,1)*xpow
       dplag=dplag+k*add
       if (x .ne. 0d0) ddplag = ddplag + k*(k-1)*add/x2
       xpow=xpow*x2
    enddo

    dd_phi_laguerre = (2d0 * dplag + 4d0*x**2*ddplag) * pi**(3d0/4d0) * phi_hm(m)
  end function dd_phi_laguerre

  !**********************************************************
  ! Laguerre specific - END
  !**********************************************************

end module collop_laguerre
