!
  MODULE rel_kernels_mod
!
  integer :: legmax=-1
!
  integer,          dimension(:),   allocatable :: nodes
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: coefleg,coefdleg,weights
!
  CONTAINS
!
  SUBROUTINE polleg_coefs(n)
!
! Computes coefficients of Legendre polynomials of orders from 0 to n

!
! Input parameters:
!           Formal: n             - maximum order of Legendre polynomials
! Output parameters:
!           module: coefleg(i,j)  - j-th coefficient of Legendre polynomial
!                   coefdleg(i,j) - j-th coefficient of derivative of Legendre 
!                                   polynomial of the order i
  IMPLICIT NONE
!
  INTEGER :: n,i,j
!
  DOUBLE PRECISION :: frontfac,rearfac
!
  IF(ALLOCATED(coefleg)) DEALLOCATE(coefleg)
  ALLOCATE(coefleg(0:n,0:n))
  IF(ALLOCATED(coefdleg)) DEALLOCATE(coefdleg)
  ALLOCATE(coefdleg(0:n,0:n))
!
  coefleg=0.d0
  coefdleg=0.d0
  coefleg(0,0)=1.d0
  coefleg(1,1)=1.d0
  frontfac=1.d0
!
  DO i=2,n
    frontfac=frontfac*(2.d0-1.d0/DBLE(i))
    rearfac=frontfac
    coefleg(i,i)=rearfac
    DO j=i-2,0,-2
      rearfac=-rearfac*DBLE(j+1)*DBLE(j+2)/DBLE(i-j)/DBLE(i+j+1)
      coefleg(i,j)=rearfac
    ENDDO
  ENDDO
!
  DO j=1,n
    coefdleg(:,j-1)=coefleg(:,j)*DBLE(j)
  ENDDO
!
  RETURN
!
  END SUBROUTINE polleg_coefs
!
  SUBROUTINE load_rel_kernels(leg)
!
  IMPLICIT NONE
!
  INTEGER, parameter :: nodemin=5
  INTEGER            :: leg,n,i,j,nodemax,ndim,info,ndh
  INTEGER,          DIMENSION(:),   ALLOCATABLE :: ipivot
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: xnode
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: amat,bvec,wint
!
  if(leg.le.legmax) return
!
  legmax=max(1,leg)
!
  call polleg_coefs(legmax)  
!
  IF(ALLOCATED(nodes)) DEALLOCATE(nodes)
  ALLOCATE(nodes(0:legmax))
!
  nodes=100
  nodemax=maxval(nodes)
!
  IF(ALLOCATED(weights)) DEALLOCATE(weights)
  ALLOCATE(weights(0:legmax,0:nodemax))
!
  weights=1.d0
  weights(:,0)=0.5d0
  weights(:,nodemax)=0.5d0
!
  ndh=7
  ndim=2*ndh
  allocate(amat(ndim,ndim),bvec(ndim,ndim),ipivot(ndim),xnode(ndim),wint(ndim-1,ndim))
!
  do j=1,ndim-1
    xnode(1)=1.d0-dble(j)
    do i=2,ndim
      xnode(i)=xnode(i-1)+1.d0
    enddo
    amat(:,1)=1.d0
    do i=2,ndim
      amat(:,i)=amat(:,i-1)*xnode
    enddo
    bvec=0.d0
    do i=1,ndim
      bvec(i,i)=1.d0
    enddo
!
    call dgesv(ndim,ndim,amat,ndim,ipivot,bvec,ndim,info)
!
    wint(j,:)=bvec(1,:)
    do i=2,ndim
      wint(j,:)=wint(j,:)+bvec(i,:)/dfloat(i)
    enddo
  enddo
!
  weights=0.d0
  do i=1,ndh
    weights(0,0:ndim-1)=weights(0,0:ndim-1)+wint(i,:)
  enddo
  do i=ndh,nodemax-ndh
    weights(0,i-ndh+1:i+ndh)=weights(0,i-ndh+1:i+ndh)+wint(ndh,:)
  enddo
  do i=ndh+1,ndim-1
    weights(0,nodemax-ndim+1:nodemax)=weights(0,nodemax-ndim+1:nodemax)+wint(i,:)
  enddo
!
  do i=1,legmax
     weights(i,:)=weights(0,:)
  enddo
!
  deallocate(amat,bvec,ipivot,xnode,wint)
!
  END SUBROUTINE load_rel_kernels
!
  END MODULE rel_kernels_mod
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine rel_kernels(leg,z_in,zp_in,R_11,R_10,R_01,R_00)
!
  use rel_kernels_mod, only : load_rel_kernels,coefleg,coefdleg,nodes,weights
!
  implicit none
!
  integer :: i,j,n,leg
  double precision :: z,zp,R_11,R_10,R_01,R_00,z_in,zp_in
  double precision :: z2,zp2,gam,gamp,pleg,dpleg,ggpm1,dyovsq,plegllp1,gz_10,gz_01
  double precision :: ymin,ymax,hy,y,y2,r,xi,cxi0,cxi1,sqrp1,onemxi2,comfac,dfllp1
!
  z=max(z_in,1d-33)
  zp=max(zp_in,1d-33)
!
  call load_rel_kernels(leg)
!
  n=nodes(leg)
!
  z2=z*z
  zp2=zp*zp
  gam=sqrt(1.d0+z2)
  gamp=sqrt(1.d0+zp2)
  ggpm1=(z2+zp2+z2*zp2)/(gam*gamp+1.d0)
  ymin=sqrt(ggpm1-z*zp)
  ymax=sqrt(ggpm1+z*zp)
  hy=(ymax-ymin)/dfloat(n)
  cxi1=1.d0/(z*zp)
  cxi0=ggpm1*cxi1
  
!
  R_11=0.d0
  R_10=0.d0
  R_01=0.d0
  R_00=0.d0
!
  do i=0,n
    y=ymin+hy*dfloat(i)
    y2=y*y
    r=1.d0+y2
    dyovsq=weights(leg,i)/sqrt(r+1.d0)
    xi=cxi0-cxi1*y2
    onemxi2=1.d0-xi*xi
!
    pleg=coefleg(leg,leg)
    dpleg=0.d0
!
    do j=leg-1,0,-1
      pleg=pleg*xi+coefleg(leg,j)
      dpleg=dpleg*xi+coefdleg(leg,j)
    enddo
!
    dfllp1=dfloat(leg*(leg+1))
    plegllp1=pleg*dfllp1
    gz_10=gam*zp/(gamp*z)
    gz_01=1.d0/gz_10
    gz_10=gz_10-xi
    gz_01=gz_01-xi
!
    R_11=R_11+dyovsq*((2.d0*r*xi+z*zp*onemxi2)*pleg-r*onemxi2*dpleg)
    R_10=R_10+dyovsq*(r*gz_10*plegllp1+(r+z*zp*gz_10)*onemxi2*dpleg)
    R_01=R_01+dyovsq*(r*gz_01*plegllp1+(r+z*zp*gz_01)*onemxi2*dpleg)
    R_00=R_00+dyovsq*onemxi2*r*(r*dpleg*dfllp1+2.d0*z*zp*(xi*dpleg-plegllp1))
  enddo
!
  comfac=0.5d0*hy*cxi1
  R_11=R_11*comfac
  R_10=R_10*comfac/zp
  R_01=R_01*comfac/z
  R_00=R_00*comfac*cxi1/(gam*gamp)
!
  end subroutine rel_kernels
