SUBROUTINE plot_distrf(source_p,source_m,eta_l,eta_r,                         &
                       eta_boundary_l,eta_boundary_r)
!
  USE lapack_band
  !USE rkstep_mod, ONLY : lag
  USe collisionality_mod, ONLY : nvel
  use nrtype, only : dp


  REAL(kind=dp), DIMENSION(:,:),   ALLOCATABLE, INTENT(inout) :: source_p
  REAL(kind=dp), DIMENSION(:,:),   ALLOCATABLE, INTENT(inout) :: source_m
  REAL(kind=dp), DIMENSION(:),     ALLOCATABLE, INTENT(inout) :: eta_l ! 0:
  REAL(kind=dp), DIMENSION(:),     ALLOCATABLE, INTENT(inout) :: eta_r ! 0:
  REAL(kind=dp)                               , INTENT(inout) :: eta_boundary_l
  REAL(kind=dp)                               , INTENT(inout) :: eta_boundary_r
!
  INTEGER, PARAMETER :: ndim=4
  INTEGER            :: nts_l,nts_r,nl,nr,i,i1,k,kk,kmax,i1min,info,m,inhom,nlam
  INTEGER            :: isig,iunit_eta_p,iunit_eta_m
  REAL(kind=dp)      :: diflam,diflampow,amin2ovb
  INTEGER,       DIMENSION(ndim)      :: ipivot
  INTEGER,       DIMENSION(1:3)       :: iunitp,iunitm,iunit_lam
  REAL(kind=dp), DIMENSION(ndim,ndim) :: amat,bvec_lapack
  REAL(kind=dp), DIMENSION(6)         :: alp,bet,gam,del
  REAL(kind=dp), DIMENSION(:),       ALLOCATABLE :: alambd
  REAL(kind=dp), DIMENSION(:,:),     ALLOCATABLE :: fun
  REAL(kind=dp), DIMENSION(:,:,:),   ALLOCATABLE :: derivs_mat
  REAL(kind=dp), DIMENSION(:,:,:),   ALLOCATABLE :: fun_write
  REAL(kind=dp), DIMENSION(:,:,:,:), ALLOCATABLE :: fun_lambda
!
  nlam=200
!
  DO k=1,3
    iunitp(k)=100+10*k
    iunitm(k)=100+10*k+1
    iunit_lam(k)=100+10*k+2
  ENDDO
  iunit_eta_p=200
  iunit_eta_m=201
!
  OPEN(iunitp(1),file='qfplu1.dat')
  OPEN(iunitm(1),file='qfmin1.dat')
  OPEN(iunit_lam(1),file='qflam1.dat')
  OPEN(iunitp(2),file='qfplu2.dat')
  OPEN(iunitm(2),file='qfmin2.dat')
  OPEN(iunit_lam(2),file='qflam2.dat')
  OPEN(iunitp(3),file='qfplu3.dat')
  OPEN(iunitm(3),file='qfmin3.dat')
  OPEN(iunit_lam(3),file='qflam3.dat')
  OPEN(iunit_eta_p,file='eta_p.dat')
  OPEN(iunit_eta_m,file='eta_m.dat')
!
  nts_l=SIZE(source_p,1)
  nts_r=SIZE(source_m,1)
  nl=nts_l/(nvel+1)
  nr=nts_r/(nvel+1)
!
  WRITE(iunit_eta_p,*) eta_l(0:nl-1),eta_boundary_l
  WRITE(iunit_eta_m,*) eta_r(0:nr-1),eta_boundary_r
  CLOSE(iunit_eta_p)
  CLOSE(iunit_eta_m)
!
  ALLOCATE(fun_lambda(0:nvel,0:3,-nlam:nlam,3))
!
!
! Left boundary:
!
!
  ALLOCATE(alambd(0:nl+2),derivs_mat(0:3,4,0:nl),fun(0:nvel,nl+2))
  ALLOCATE(fun_write(0:nvel,0:3,0:nl))
!
  DO i=0,nl-1
    alambd(i)=SQRT(1.d0-eta_l(i)/eta_boundary_l+10.d0*EPSILON(1.d0))
  ENDDO
  alambd(nl)=0.d0
  alambd(nl+1)=-SQRT(1.d0-eta_r(nr-1)/eta_boundary_r+10.d0*EPSILON(1.d0))
  alambd(nl+2)=-SQRT(1.d0-eta_r(nr-2)/eta_boundary_r+10.d0*EPSILON(1.d0))
!
  amin2ovb=-2.d0*eta_boundary_l
!
  DO i=0,nl
!
    i1min=MAX(0,i-2)
    kmax=5
! 
    DO k=1,kmax
      i1=k-1+i1min
      diflam=alambd(i1)-alambd(i)
      diflampow=diflam
      alp(k)=(alambd(i)+diflam/2.d0)*diflampow
      diflampow=diflam*diflampow
      bet(k)=(alambd(i)/2.d0+diflam/3.d0)*diflampow
      diflampow=diflam*diflampow
      gam(k)=(alambd(i)/6.d0+diflam/8.d0)*diflampow
      diflampow=diflam*diflampow
      del(k)=(alambd(i)/24.d0+diflam/30.d0)*diflampow
    ENDDO
!
    DO k=1,4
      amat(k,1)=(alp(k+1)-alp(k))*amin2ovb
      amat(k,2)=(bet(k+1)-bet(k))*amin2ovb
      amat(k,3)=(gam(k+1)-gam(k))*amin2ovb
      amat(k,4)=(del(k+1)-del(k))*amin2ovb
    ENDDO
!
    bvec_lapack=0.d0
    DO k=1,ndim
      bvec_lapack(k,k)=1.d0
    ENDDO
!
    CALL gbsv(ndim,ndim,amat,ipivot,bvec_lapack,info)
!
    derivs_mat(0:3,1:4,i)=bvec_lapack(1:4,1:4)
!
  ENDDO
!
  DO inhom=1,3
!
    DO m=0,nvel
      k=m*nl
      kk=m*nr
      fun(m,1:nl)=source_p(k+1:k+nl,inhom)
      fun(m,nl+1)=-source_m(kk+nr,inhom)
      fun(m,nl+2)=-source_m(kk+nr-1,inhom)
!
      fun_write(m,0:3,0)=MATMUL(derivs_mat(0:3,1:4,0),fun(m,1:4))
      fun_write(m,0:3,1)=MATMUL(derivs_mat(0:3,1:4,1),fun(m,1:4))
!
      DO i=2,nl
        fun_write(m,0:3,i)=MATMUL(derivs_mat(0:3,1:4,i),fun(m,i-1:i+2))
      ENDDO
!
    ENDDO
!
    WRITE(iunitp(inhom),*) fun_write(0,0,:)
!
    isig=1
!
!    CALL make_fun_lambda(nvel,nl,nlam,isig,eta_l(0:nl),eta_boundary_l,    &
!                         fun_write(0:nvel,0:3,0:nl),                      &
!                         fun_lambda(0:nvel,0:3,0:nlam,inhom))
!
  ENDDO
!
  DEALLOCATE(alambd,derivs_mat,fun,fun_write)
!
!
! Right boundary:
!
!
  ALLOCATE(alambd(0:nr+2),derivs_mat(0:3,4,0:nr),fun(0:nvel,nr+2))
  ALLOCATE(fun_write(0:nvel,0:3,0:nr))
!
  DO i=0,nr-1
    alambd(i)=SQRT(1.d0-eta_r(i)/eta_boundary_r+10.d0*EPSILON(1.d0))
  ENDDO
  alambd(nr)=0.d0
  alambd(nr+1)=-SQRT(1.d0-eta_l(nl-1)/eta_boundary_l+10.d0*EPSILON(1.d0))
  alambd(nr+2)=-SQRT(1.d0-eta_l(nl-2)/eta_boundary_l+10.d0*EPSILON(1.d0))
!
  amin2ovb=-2.d0*eta_boundary_r
!
  DO i=0,nr
!
    i1min=MAX(0,i-2)
    kmax=5
! 
    DO k=1,kmax
      i1=k-1+i1min
      diflam=alambd(i1)-alambd(i)
      diflampow=diflam
      alp(k)=(alambd(i)+diflam/2.d0)*diflampow
      diflampow=diflam*diflampow
      bet(k)=(alambd(i)/2.d0+diflam/3.d0)*diflampow
      diflampow=diflam*diflampow
      gam(k)=(alambd(i)/6.d0+diflam/8.d0)*diflampow
      diflampow=diflam*diflampow
      del(k)=(alambd(i)/24.d0+diflam/30.d0)*diflampow
    ENDDO
!
    DO k=1,4
      amat(k,1)=(alp(k+1)-alp(k))*amin2ovb
      amat(k,2)=(bet(k+1)-bet(k))*amin2ovb
      amat(k,3)=(gam(k+1)-gam(k))*amin2ovb
      amat(k,4)=(del(k+1)-del(k))*amin2ovb
    ENDDO
!
    bvec_lapack=0.d0
    DO k=1,ndim
      bvec_lapack(k,k)=1.d0
    ENDDO
!
    CALL gbsv(ndim,ndim,amat,ipivot,bvec_lapack,info)
!
    derivs_mat(0:3,1:4,i)=bvec_lapack(1:4,1:4)
!
  ENDDO
!
  DO inhom=1,3
!
    DO m=0,nvel
      k=m*nr
      kk=m*nl
      fun(m,1:nr)=source_m(k+1:k+nr,inhom)
      fun(m,nr+1)=-source_p(kk+nl,inhom)
      fun(m,nr+2)=-source_p(kk+nl-1,inhom)
!
      fun_write(m,0:3,0)=MATMUL(derivs_mat(0:3,1:4,0),fun(m,1:4))
      fun_write(m,0:3,1)=MATMUL(derivs_mat(0:3,1:4,1),fun(m,1:4))
!
      DO i=2,nr
        fun_write(m,0:3,i)=MATMUL(derivs_mat(0:3,1:4,i),fun(m,i-1:i+2))
      ENDDO
!
    ENDDO
!
    WRITE(iunitm(inhom),*) fun_write(0,0,:)
!
    isig=-1
!
!    CALL make_fun_lambda(nvel,nr,nlam,isig,eta_r(0:nr),eta_boundary_r,    &
!                         fun_write(0:nvel,0:3,0:nr),                      &
!                         fun_lambda(0:nvel,0:3,-nlam:0,inhom))
!
  ENDDO
!
  DEALLOCATE(alambd,derivs_mat,fun,fun_write)
!
  DO inhom=1,3
!    WRITE(iunit_lam(inhom),*) fun_lambda(0,0,:,inhom)
  ENDDO
!
  DEALLOCATE(fun_lambda)
!
  DO k=1,3
    CLOSE(iunitp(k))
    CLOSE(iunitm(k))
    CLOSE(iunit_lam(k))
  ENDDO
!
END SUBROUTINE plot_distrf
