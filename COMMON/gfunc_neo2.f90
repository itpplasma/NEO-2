!
  module spitzer_fun_mod
    integer :: lag,nplp1,nphi
    INTEGER,          DIMENSION(:),        ALLOCATABLE :: npass
    DOUBLE PRECISION, DIMENSION(:),        ALLOCATABLE :: phi_mfl,bhat_mfl,eta
    DOUBLE PRECISION, DIMENSION(:,:),      ALLOCATABLE :: coeflag
    DOUBLE PRECISION, DIMENSION(:,:,:,:),  ALLOCATABLE :: fun_p,fun_m
    DOUBLE PRECISION, DIMENSION(:),        ALLOCATABLE :: spfun,dspfun
  end module spitzer_fun_mod
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine load_spitzer_fun(iphi0)
!
  use spitzer_fun_mod
!
  IMPLICIT NONE
  INTEGER, PARAMETER :: dp = KIND(1.0d0)
  DOUBLE PRECISION, PARAMETER  :: pi=3.14159265358979d0
!
  INTEGER :: iunit,iunit_p,iunit_m,iphi,iphi0,m
  DOUBLE PRECISION :: phi,alam,collpar
!
!
  iunit=721
  iunit_p=722
  iunit_m=723
!
  open(iunit,file='NEO2_DATA/sizeplot_etalev.dat')
  read(iunit,*) lag,nplp1,nphi,collpar
  allocate(eta(0:nplp1))
  read(iunit,*) eta(0:nplp1)
  close(iunit)
!
  allocate(phi_mfl(nphi),npass(nphi),bhat_mfl(nphi))
  allocate(fun_p(0:lag,0:3,0:nplp1,nphi),fun_m(0:lag,0:3,0:nplp1,nphi))
  open(iunit,file='NEO2_DATA/phi_mesh.dat')
  open(iunit_p,form='unformatted',file='NEO2_DATA/spitf_p.dat')
  open(iunit_m,form='unformatted',file='NEO2_DATA/spitf_m.dat')
  do iphi=1,nphi
  read(iunit,*) phi_mfl(iphi),npass(iphi),bhat_mfl(iphi)
  read(iunit_p) fun_p(0:lag,0:3,0:nplp1,iphi)
  read(iunit_m) fun_m(0:lag,0:3,0:nplp1,iphi)
  enddo
  close(iunit)
  close(iunit_p)
  close(iunit_m)
  collpar=collpar*3.d0*sqrt(pi)/8.d0
  fun_p=fun_p*collpar
  fun_m=fun_m*collpar
!
  allocate(coeflag(0:lag,0:lag))
!
  call pollag(lag,coeflag)
!
  allocate(spfun(0:lag),dspfun(0:lag))
!
  print *,phi_mfl(1),' < phi < ',phi_mfl(nphi)
  print *,'phi?'
  read *,phi
!
  if(phi.lt.phi_mfl(1).or.phi.gt.phi_mfl(nphi)) then
    print *,'phi out of range'
    stop
  endif
  do iphi=2,nphi
    if(phi_mfl(iphi).gt.phi) then
      if(abs(phi_mfl(iphi)-phi).lt.abs(phi_mfl(iphi-1)-phi)) then
        iphi0=iphi
      else
        iphi0=iphi-1
      endif
      exit
    endif
  enddo
!
  end subroutine load_spitzer_fun
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  SUBROUTINE spfun_lambda(iphi,alam)
!
  use spitzer_fun_mod
!
  IMPLICIT NONE
!
  INTEGER :: m,i,iphi
  DOUBLE PRECISION :: alam
  DOUBLE PRECISION :: etafix,dellam
!
  i=nplp1
  etafix=(1.d0-alam**2)/bhat_mfl(iphi)
  if(etafix.lt.0.d0) then
    spfun=0.d0
    dspfun=0.d0
    return
  endif
  DO WHILE(eta(i).GT.etafix)
    IF(i.LE.1) EXIT
    i=i-1
  ENDDO
  dellam=abs(alam)-SQRT(1.d0-eta(i)*bhat_mfl(iphi))
  IF(alam.gt.0) THEN
    DO m=0,lag
      spfun(m) =             fun_p(m,0,i,iphi)                  &
               +dellam*     (fun_p(m,1,i,iphi)                  &
               +dellam/2.d0*(fun_p(m,2,i,iphi)                  &
               +dellam/3.d0*(fun_p(m,3,i,iphi))))
      dspfun(m)=(            fun_p(m,1,i,iphi)                  &
               +dellam*     (fun_p(m,2,i,iphi)                  &
               +dellam/2.d0*(fun_p(m,3,i,iphi))))
    ENDDO
  ELSE
    DO m=0,lag
      spfun(m) =             fun_m(m,0,i,iphi)                  &
               +dellam*     (fun_m(m,1,i,iphi)                  &
               +dellam/2.d0*(fun_m(m,2,i,iphi)                  &
               +dellam/3.d0*(fun_m(m,3,i,iphi))))
      dspfun(m)=-(           fun_m(m,1,i,iphi)                  &
               +dellam*     (fun_m(m,2,i,iphi)                  &
               +dellam/2.d0*(fun_m(m,3,i,iphi))))
    ENDDO
  ENDIF
!
  RETURN
  END SUBROUTINE spfun_lambda
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  SUBROUTINE pollag(n,coeflag)
!
! Computes coefficients of Laguerre polynomials of orders from 0 to n

!
! Input parameters:
!           Formal: n            - maximum order of Leguerre polynomials
! Output parameters:
!           Formal: coeflag(i,j) - j-th coefficient of re-normalized Leguerre 
!                                  polynomial (function $\varphi$)
!                                  of the order i
  IMPLICIT NONE
!
  INTEGER :: n,i,j,k
!
  DOUBLE PRECISION, PARAMETER  :: pi=3.14159265358979d0
  DOUBLE PRECISION :: bincoef,hm
  DOUBLE PRECISION, DIMENSION(0:n,0:n) :: coeflag
!
  coeflag=0.d0
!
  hm=3.d0/(8.d0*pi)
  coeflag(0,0)=1.d0/sqrt(hm)
!
  DO i=1,n
    hm=hm*(1.d0+1.5d0/i)
    DO k=0,i
      bincoef=1.d0
      do j=1,i-k
        bincoef=bincoef*(1.d0+(k+1.5d0)/j)
      enddo
      do j=1,k
        bincoef=-bincoef/j
      enddo
      coeflag(i,k)=bincoef
    ENDDO
    coeflag(i,0:i)=coeflag(i,0:i)/sqrt(hm)
  ENDDO
!
  RETURN
END SUBROUTINE pollag
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine laguerre(m,x,plag,dplag)
!
  use spitzer_fun_mod
!
  IMPLICIT NONE
!
  integer :: m,k
  double precision :: x,plag,dplag,xpow,add
!
  if(m.gt.lag) then
    print *,'laguerre : order of polynomial out of range'
    stop
  endif
!
  plag=coeflag(m,0)
  dplag=0.d0
  xpow=1.d0
!
  do k=1,m
    add=coeflag(m,k)*xpow
    dplag=dplag+k*add
    plag=plag+add*x
    xpow=xpow*x
  enddo
!
  return
  end subroutine laguerre
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine gfunc_neo2(xtr,xpa,g,gtr,gpa)
!
  use spitzer_fun_mod
!
  IMPLICIT NONE
!
  integer, parameter :: npmod=100, nlam=99
!
  integer :: i,j,m
  DOUBLE PRECISION :: alam,x,xx,plag,g,plag0,xperp,xpar,gg,h,dplag
  DOUBLE PRECISION :: xtr,xpa,gtr,gpa,dg_dlam,dg_dx,epst
  DOUBLE PRECISION, dimension(3,-nlam:nlam) :: arr
  logical, save :: prop=.true.
  integer, save :: iphi
  double precision, save :: transfac
!
  if(prop) then
    prop=.false.
!
    call load_spitzer_fun(iphi)
!
    epst=0.25d0
    transfac=1.d0/sqrt(1.d0-epst**2)
  endif
!
  x=sqrt(xtr**2+xpa**2)
  if(x.ne.0.d0) then
    alam=xpa/x
  else
    g=0.d0
    gtr=0.d0
    gpa=0.d0
    return
  endif
  xx=x**2
  call spfun_lambda(iphi,alam)
  g=0.d0
  dg_dx=0.d0
  dg_dlam=0.d0
  do m=0,lag
    call laguerre(m,xx,plag,dplag)
    g=g+plag*spfun(m)
    dg_dx=dg_dx+2.d0*x*dplag*spfun(m)
    dg_dlam=dg_dlam+plag*dspfun(m)
  enddo
  gtr=sqrt(abs(1.d0-alam**2))*(dg_dx-dg_dlam*alam/x)
  gpa=alam*dg_dx+dg_dlam*(1.d0-alam**2)/x
!
  g=g*transfac
  gtr=gtr*transfac
  gpa=gpa*transfac
!
  end subroutine gfunc_neo2
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  implicit none
!
  integer, parameter :: npmod=100, nlam=99
  integer :: i,j
  DOUBLE PRECISION :: hx,halam,xtr,xpa,x,alam,g,gtr,gpa
  DOUBLE PRECISION, dimension(3,-nlam:nlam) :: arr
!
  hx=5.d0/npmod
  halam=1.d0/nlam
!
  do i=0,npmod
    x=i*hx
    do j=-nlam,nlam
      alam=j*halam
      xtr=x*sqrt(1.d0-alam**2)
      xpa=x*alam
      call gfunc_neo2(xtr,xpa,g,gtr,gpa)
      arr(1,j)=g
      arr(2,j)=gtr
      arr(3,j)=gpa
      write (1,*) x,alam,g
    enddo
    write (21,*) arr(1,:)
    write (22,*) arr(2,:)
    write (23,*) arr(3,:)
    write (31,*) x,arr(:,nlam)
  enddo
!
  end

