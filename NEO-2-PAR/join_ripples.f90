! THIS IS A NOBAND VERSION - MODIFIED CALL
!  modified for _g and _e drive
!  comments added by Winny
!  checked with write-up
!
SUBROUTINE join_ripples(ierr)

!
! Joins the results from the "old ripple", n, (naming o%...)
! and the "new ripple", n+1, (naming n%...).
! Puts the result in "old ripple" variables.
!
! WINNY:
!
! c_forward and c_backward should be allocatable arrays
!  they will be passed most likely through the parameter list
!  do you really want to make matmul in cases when binarysplit is not done?
!
!  c_forward  will have dim (n%p%npass_l , o%p%npass_r)
!  c_backward will have dim (o%p%npass_r , n%p%npass_l)
!
! 

  USE propagator_mod
  USE lapack_band
  USE collisionality_mod, ONLY : isw_lorentz
USE development
!
  IMPLICIT NONE
  INTEGER, PARAMETER :: dp = KIND(1.0d0)
!
  INTEGER, INTENT(out) :: ierr
!
  TYPE(propagator), POINTER            :: o
  TYPE(propagator), POINTER            :: n

!
  INTEGER :: ndim,ndim1,k,k1,i,i1
  INTEGER :: nvel,m
  INTEGER :: info
  INTEGER,          DIMENSION(:),   ALLOCATABLE :: ipivot
!
  integer, parameter :: joinripples_write=0
!
  ! WINNY
  ! made it compatible with propagator_mod
  !
  ! One can finally remove this c_forward and c_backward if we just
  ! use o%p%cmat (for forward) and n%p%cmat (for backward)
  !
  ! DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: c_forward,c_backward
  DOUBLE PRECISION, DIMENSION(:,:), POINTER :: c_forward,c_backward

  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: bvec,bvec_min
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: amat,bvec_lapack,prod
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: bvec_min_lapack
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: dummat1,dummat2
  DOUBLE PRECISION :: facnorm

  integer :: icount = 0
  !
  ! WINNY for SERGEI
  !  The following lines were forgotten in hurry in the evening so 
  !  o (old) and n (new) pointed nowhere. My mistake.
  !
  ! initialize
  ierr = 0
  !o => prop_c%prev
  !n => prop_c
  o => prop_c_old
  n => prop_c_new

  c_forward  => o%p%cmat
  c_backward => n%p%cmat
  !
  i=SIZE(c_forward,1)
  i1=SIZE(c_forward,2)
  IF(SIZE(c_backward,1).NE.i1.OR.SIZE(c_backward,2).NE.i) THEN
    PRINT *,'join_ripples : sizes of c_forward and c_backward do not fit'
    ierr=1
    RETURN
  ELSEIF(n%p%npass_l.NE.i) THEN
    PRINT *,'join_ripples : size 1 of c_forward differs from new npass_l'
    ierr=1
    RETURN
  ELSEIF(o%p%npass_r.NE.i1) THEN
    PRINT *,'join_ripples : size 2 of c_forward differs from old npass_r'
    ierr=1
    RETURN
  ELSE
    ierr=0
  ENDIF
!
  nvel=n%p%nvelocity
!
  ndim=n%p%npass_l*(nvel+1)
  ndim1=o%p%npass_r*(nvel+1)
  ALLOCATE(amat(ndim,ndim),bvec_lapack(ndim,3),bvec_min_lapack(ndim,3))
  ALLOCATE(ipivot(ndim))
  !
  ! coefficient matrix A 
  amat=0.d0
  DO i=1,ndim
    amat(i,i)=1.d0
  ENDDO
!  amat=amat-MATMUL(c_forward,MATMUL(o%p%amat_m_p,                   &
!            MATMUL(c_backward,n%p%amat_p_m)))
!
  ALLOCATE(dummat1(ndim,ndim1),dummat2(ndim1,ndim))
!
  DO m=0,nvel
    k=m*n%p%npass_l
    k1=m*o%p%npass_r
    dummat1(k+1:k+n%p%npass_l,:)                                    &
           =MATMUL(c_forward,o%p%amat_m_p(k1+1:k1+o%p%npass_r,:))
    dummat2(k1+1:k1+o%p%npass_r,:)                                  &
           =MATMUL(c_backward,n%p%amat_p_m(k+1:k+n%p%npass_l,:))
  ENDDO
!
  amat=amat-MATMUL(dummat1,dummat2)
!
  DEALLOCATE(dummat1,dummat2)
!
  !
  ! source terms
  ! integrals of particle and heat fluxes and current over fieldline
  ! 
  ! gradient drive
  ! rhs of algebraic equation
!
!   bvec=MATMUL(c_forward,o%p%source_p_g+MATMUL(o%p%amat_m_p,        &
!        MATMUL(c_backward,n%p%source_m_g)))
!
  ALLOCATE(dummat1(ndim1,3),dummat2(ndim1,3))
!
  DO m=0,nvel
    k=m*n%p%npass_l
    k1=m*o%p%npass_r
    dummat1(k1+1:k1+o%p%npass_r,:)                                  &
           =MATMUL(c_backward,n%p%source_m(k+1:k+n%p%npass_l,:))
  ENDDO
!
  dummat2=o%p%source_p+MATMUL(o%p%amat_m_p,dummat1)
!
  DO m=0,nvel
    k=m*n%p%npass_l
    k1=m*o%p%npass_r
    bvec_lapack(k+1:k+n%p%npass_l,:)                                &
           =MATMUL(c_forward,dummat2(k1+1:k1+o%p%npass_r,:))
  ENDDO
!
  DEALLOCATE(dummat2)
!
  !
  ! solution bvec -> $\ifour{\fvec}{+}{}{o+1}{g,l}$
!
  CALL gbsv(ndim,ndim,amat,ipivot,bvec_lapack,info)
!
!  CALL gbsv(ndim,ndim,amat,ipivot,bvec,info)
  IF(info.NE.0) THEN
    ierr=2
    RETURN
  ENDIF
  !
  ! bvec_min -> $\ifour{\fvec}{-}{}{o}{g,r}$
!  bvec_min=n%p%source_m_g+MATMUL(n%p%amat_p_m,bvec)
  bvec_min_lapack=n%p%source_m+MATMUL(n%p%amat_p_m,bvec_lapack)
  !
  ! sources
!  DEALLOCATE(o%p%source_p_g)
!  ALLOCATE(o%p%source_p_g(n%p%npass_r))
  DEALLOCATE(o%p%source_p)
  ALLOCATE(o%p%source_p(n%p%npass_r*(nvel+1),3))
!  o%p%source_p_g=n%p%source_p_g+MATMUL(n%p%amat_p_p,bvec)
  o%p%source_p=n%p%source_p+MATMUL(n%p%amat_p_p,bvec_lapack)
  !
!  o%p%source_m_g=o%p%source_m_g+MATMUL(o%p%amat_m_m,MATMUL(c_backward,bvec_min))
!
  DO m=0,nvel
    k=m*n%p%npass_l
    k1=m*o%p%npass_r
    dummat1(k1+1:k1+o%p%npass_r,:)                                  &
           =MATMUL(c_backward,bvec_min_lapack(k+1:k+n%p%npass_l,:))
  ENDDO
!
  o%p%source_m=o%p%source_m+MATMUL(o%p%amat_m_m,dummat1)
!
  !
  ! integrals of particle flux and current over fieldline 
!  o%p%qflux_g=o%p%qflux_g+n%p%qflux_g                               &
!             +SUM(o%p%flux_m*MATMUL(c_backward,bvec_min))           &
!             +SUM(n%p%flux_p*bvec)
!  o%p%qcurr_g=o%p%qcurr_g + n%p%qcurr_g                             &
!             +SUM(o%p%curr_m*MATMUL(c_backward,bvec_min))           &
!             +SUM(n%p%curr_p*bvec)
  o%p%qflux=o%p%qflux+n%p%qflux                                     &
             +MATMUL(o%p%flux_m,dummat1)                            &
             +MATMUL(n%p%flux_p,bvec_lapack)
  DEALLOCATE(dummat1)
  !open(1003,file='testonsager.dat',position='append')
  !icount = icount + 1
  !write(1003,*) icount,o%p%qflux(1:2,1:2)
  !close(1003)
  !
  ! source terms
  ! integrals of particle flux and current over fieldline
  ! 
!  ! parallel electric field drive
!  ! rhs of algebraic equation
!  bvec=MATMUL(c_forward,o%p%source_p_e+MATMUL(o%p%amat_m_p,         &
!       MATMUL(c_backward,n%p%source_m_e)))
!  !
!  ! solution bvec -> $\ifour{\fvec}{+}{}{o+1}{e,l}$
!  CALL gbsv(ndim,ndim,amat,ipivot,bvec,info)
!  IF(info.NE.0) THEN
!    ierr=2
!    RETURN
!  ENDIF
!  !
!  ! bvec_min -> $\ifour{\fvec}{-}{}{o}{e,r}$
!  bvec_min=n%p%source_m_e+MATMUL(n%p%amat_p_m,bvec)
!  !
!  ! sources
!  DEALLOCATE(o%p%source_p_e)
!  ALLOCATE(o%p%source_p_e(n%p%npass_r))
!  o%p%source_p_e=n%p%source_p_e+MATMUL(n%p%amat_p_p,bvec)
!  !
!  o%p%source_m_e=o%p%source_m_e+MATMUL(o%p%amat_m_m,MATMUL(c_backward,bvec_min))
!  !
!  ! integrals of particle flux and current over fieldline 
!  o%p%qflux_e=o%p%qflux_e+n%p%qflux_e+SUM(o%p%flux_m                &
!             *MATMUL(c_backward,bvec_min))+SUM(n%p%flux_p*bvec)
!  o%p%qcurr_e=o%p%qcurr_e+n%p%qcurr_e+SUM(o%p%curr_m                &
!             *MATMUL(c_backward,bvec_min))+SUM(n%p%curr_p*bvec)

  !
  ! drive: f^+ from period n (from the left)
  !
  ! rhs of equation
!  ALLOCATE(bvec_lapack(ndim,o%p%npass_l),bvec_min_lapack(ndim,o%p%npass_l))
!  bvec_lapack=MATMUL(c_forward,o%p%amat_p_p)
  DEALLOCATE(bvec_lapack,bvec_min_lapack)
  ALLOCATE(bvec_lapack(ndim,o%p%npass_l*(nvel+1)))
  ALLOCATE(bvec_min_lapack(ndim,o%p%npass_l*(nvel+1)))
!
  DO m=0,nvel
    k=m*n%p%npass_l
    k1=m*o%p%npass_r
    bvec_lapack(k+1:k+n%p%npass_l,:)                                 &
        =MATMUL(c_forward,o%p%amat_p_p(k1+1:k1+o%p%npass_r,:))
  ENDDO
  !
  ! solution bvec_lapack -> $\ifour{\Fmat}{+}{}{o+1}{l,l}$
  CALL gbsv(ndim,ndim,amat,ipivot,bvec_lapack,info)
  IF(info.NE.0) THEN
    ierr=3
    RETURN
  ENDIF
  !
  ! bvec_lapack = f^+ from period n+1 driven by f^+ from period n
  ! bvec_min_lapack = f^- from period n driven by f^+ from period n
  !
  bvec_min_lapack=MATMUL(n%p%amat_p_m,bvec_lapack)
  !
!  o%p%amat_p_m=o%p%amat_p_m+MATMUL(o%p%amat_m_m,                     &
!               MATMUL(c_backward,bvec_min_lapack))
  !
  ALLOCATE(dummat1(ndim1,o%p%npass_l*(nvel+1)))
!
  DO m=0,nvel
    k=m*n%p%npass_l
    k1=m*o%p%npass_r
    dummat1(k1+1:k1+o%p%npass_r,:)                                   &
          =MATMUL(c_backward,bvec_min_lapack(k+1:k+n%p%npass_l,:))
  ENDDO
  !
  o%p%amat_p_m=o%p%amat_p_m+MATMUL(o%p%amat_m_m,dummat1) 
!
  DEALLOCATE(o%p%amat_p_p)
!  ALLOCATE(o%p%amat_p_p(n%p%npass_r,o%p%npass_l))
  ALLOCATE(o%p%amat_p_p(n%p%npass_r*(nvel+1),o%p%npass_l*(nvel+1)))
  o%p%amat_p_p=MATMUL(n%p%amat_p_p,bvec_lapack)
  !
!  o%p%flux_p=o%p%flux_p+MATMUL(n%p%flux_p,bvec_lapack)               &
!            +MATMUL(o%p%flux_m,MATMUL(c_backward,bvec_min_lapack))
!  o%p%curr_p=o%p%curr_p+MATMUL(n%p%curr_p,bvec_lapack)               &
!            +MATMUL(o%p%curr_m,MATMUL(c_backward,bvec_min_lapack))
  o%p%flux_p=o%p%flux_p+MATMUL(n%p%flux_p,bvec_lapack)               &
            +MATMUL(o%p%flux_m,dummat1)
!  DEALLOCATE(bvec_lapack,bvec_min_lapack)
  DEALLOCATE(bvec_lapack,bvec_min_lapack,dummat1)

  !
  ! drive: f^- from period n+1 (from right)
  !
!  ALLOCATE(bvec_lapack(ndim,n%p%npass_r))
!  ALLOCATE(bvec_min_lapack(ndim,n%p%npass_r))
  ALLOCATE(bvec_lapack(ndim,n%p%npass_r*(nvel+1)))
  ALLOCATE(bvec_min_lapack(ndim,n%p%npass_r*(nvel+1)))
  ALLOCATE(dummat1(ndim,ndim1),dummat2(ndim1,n%p%npass_r*(nvel+1)))
  !
  ! rhs of equation
!  bvec_lapack=MATMUL(c_forward,MATMUL(o%p%amat_m_p,                  &
!              MATMUL(c_backward,n%p%amat_m_m)))
  DO m=0,nvel
    k=m*n%p%npass_l
    k1=m*o%p%npass_r
    dummat1(k+1:k+n%p%npass_l,:)                                    &
           =MATMUL(c_forward,o%p%amat_m_p(k1+1:k1+o%p%npass_r,:))
    dummat2(k1+1:k1+o%p%npass_r,:)                                  &
           =MATMUL(c_backward,n%p%amat_m_m(k+1:k+n%p%npass_l,:))
  ENDDO
  !
  bvec_lapack=MATMUL(dummat1,dummat2)
  !
  ! solution bvec_lapack -> $\ifour{\Fmat}{+}{}{o+1}{r,l}$
  CALL gbsv(ndim,ndim,amat,ipivot,bvec_lapack,info)
  IF(info.NE.0) THEN
    ierr=4
    RETURN
  ENDIF
  !
  bvec_min_lapack=n%p%amat_m_m+MATMUL(n%p%amat_p_m,bvec_lapack)
  !
  DEALLOCATE(amat)
!  ALLOCATE(amat(o%p%npass_l,n%p%npass_r))
  ALLOCATE(amat(o%p%npass_l*(nvel+1),n%p%npass_r*(nvel+1)))
!  amat=MATMUL(o%p%amat_m_m,MATMUL(c_backward,bvec_min_lapack))
!
  DO m=0,nvel
    k=m*n%p%npass_l
    k1=m*o%p%npass_r
    dummat2(k1+1:k1+o%p%npass_r,:)                                  &
           =MATMUL(c_backward,bvec_min_lapack(k+1:k+n%p%npass_l,:))
  ENDDO
!
  amat=MATMUL(o%p%amat_m_m,dummat2)
!
  DEALLOCATE(o%p%amat_m_m)
  ALLOCATE(o%p%amat_m_m(o%p%npass_l*(nvel+1),n%p%npass_r*(nvel+1)))
  o%p%amat_m_m=amat
  !
  DEALLOCATE(o%p%amat_m_p)
!  ALLOCATE(o%p%amat_m_p(n%p%npass_r,n%p%npass_r))
  ALLOCATE(o%p%amat_m_p(n%p%npass_r*(nvel+1),n%p%npass_r*(nvel+1)))
  o%p%amat_m_p=n%p%amat_m_p+MATMUL(n%p%amat_p_p,bvec_lapack) 
  !
!  n%p%flux_m=n%p%flux_m+MATMUL(n%p%flux_p,bvec_lapack)             &
!            +MATMUL(o%p%flux_m,MATMUL(c_backward,bvec_min_lapack))
!  n%p%curr_m=n%p%curr_m+MATMUL(n%p%curr_p,bvec_lapack)             &
!            +MATMUL(o%p%curr_m,MATMUL(c_backward,bvec_min_lapack))
  n%p%flux_m=n%p%flux_m+MATMUL(n%p%flux_p,bvec_lapack)             &
            +MATMUL(o%p%flux_m,dummat2)
  !
!  DEALLOCATE(o%p%flux_m,o%p%curr_m)
  DEALLOCATE(o%p%flux_m)
!  ALLOCATE(o%p%flux_m(n%p%npass_r),o%p%curr_m(n%p%npass_r))
  ALLOCATE(o%p%flux_m(3,n%p%npass_r*(nvel+1)))
  DEALLOCATE(dummat1,dummat2)
  !
  o%p%flux_m=n%p%flux_m
!  o%p%curr_m=n%p%curr_m
  !
  o%p%npass_r=n%p%npass_r
  !
!  DEALLOCATE(amat,bvec,bvec_min,bvec_lapack,bvec_min_lapack)
  DEALLOCATE(amat,bvec_lapack,bvec_min_lapack)
  DEALLOCATE(ipivot)
  ! WINNY
  ! DEALLOCATE(c_forward,c_backward)
  NULLIFY(c_forward)
  NULLIFY(c_backward)
!
!
! Correction of flux conservation (for Lorentz model only)
!
  IF(isw_lorentz.EQ.1) THEN
    DO i=1,o%p%npass_l
      facnorm=SUM(o%p%amat_p_p(:,i))+SUM(o%p%amat_p_m(:,i))
      o%p%amat_p_p(:,i)=o%p%amat_p_p(:,i)/facnorm
      o%p%amat_p_m(:,i)=o%p%amat_p_m(:,i)/facnorm
    ENDDO
    DO i=1,o%p%npass_r
      facnorm=SUM(o%p%amat_m_p(:,i))+SUM(o%p%amat_m_m(:,i))
      o%p%amat_m_p(:,i)=o%p%amat_m_p(:,i)/facnorm
      o%p%amat_m_m(:,i)=o%p%amat_m_m(:,i)/facnorm
    ENDDO
  ENDIF
!
if(joinripples_write.eq.1) then
OPEN(111,file='tot_flux_p.dat')
WRITE(111,'(3(1x,e12.5))') (o%p%flux_p(:,i) ,i=1,o%p%npass_l)
CLOSE(111)
OPEN(111,file='tot_flux_m.dat')
WRITE(111,'(3(1x,e12.5))') (o%p%flux_m(:,i) ,i=1,o%p%npass_r)
CLOSE(111)
OPEN(111,file='tot_source_p.dat')
WRITE(111,'(3(1x,e12.5))') (o%p%source_p(i,:) ,i=1,o%p%npass_r)
CLOSE(111)
OPEN(111,file='tot_source_m.dat')
WRITE(111,'(3(1x,e12.5))') (o%p%source_m(i,:) ,i=1,o%p%npass_l)
CLOSE(111)
OPEN(111,file='tot_qflux.dat',position='append')
WRITE(111,'(4(1x,e12.5))') o%p%qflux(1:2,1:2)
CLOSE(111)
OPEN(111,file='tot_amat_p_p.dat')
DO i=1,o%p%npass_r
WRITE(111,*) o%p%amat_p_p(i,1:o%p%npass_l)
ENDDO
CLOSE(111)
OPEN(111,file='tot_amat_p_m.dat')
DO i=1,o%p%npass_l
WRITE(111,*) o%p%amat_p_m(i,1:o%p%npass_l)
ENDDO
CLOSE(111)
OPEN(111,file='tot_amat_m_p.dat')
DO i=1,o%p%npass_r
WRITE(111,*) o%p%amat_m_p(i,1:o%p%npass_r)
ENDDO
CLOSE(111)
OPEN(111,file='tot_amat_m_m.dat')
DO i=1,o%p%npass_l
WRITE(111,*) o%p%amat_m_m(i,1:o%p%npass_r)
ENDDO
CLOSE(111)
!pause 'join ripples written' !Warning in gfortran-4.7
endif

!
  !
  RETURN
END SUBROUTINE join_ripples
!
!
!
SUBROUTINE duplicate_ripple(nvel,ndim,ndim1,npass_l,npass_r,                   &
                            c_forward,c_backward,                             &
                            amat_p_p,amat_p_m,amat_m_m,amat_m_p,              &
                            source_p,source_m,flux_p,flux_m,qflux,            &
                            damat_p_p,damat_p_m,damat_m_m,damat_m_p,          &
                            dsource_p,dsource_m,dflux_p,dflux_m,dqflux,       &
                            ierr)
!
! Joins the periodic propagator to itself
! Puts the result into a "double" propagator.
!
  USE lapack_band
  USE collisionality_mod, ONLY : isw_lorentz
!
  IMPLICIT NONE
  INTEGER, PARAMETER :: dp = KIND(1.0d0)
!
  INTEGER, INTENT(out) :: ierr
!
  INTEGER :: ndim,ndim1,k,k1,i,i1,npass_l,npass_r
  INTEGER :: nvel,m
  INTEGER :: info
  INTEGER,          DIMENSION(:),   ALLOCATABLE :: ipivot
!
  DOUBLE PRECISION, DIMENSION(npass_l,npass_r) :: c_forward
  DOUBLE PRECISION, DIMENSION(npass_r,npass_l) :: c_backward
!
!  ndim=npass_l*(nvel+1)
!  ndim1=npass_r*(nvel+1)
!
  DOUBLE PRECISION, DIMENSION(3,3)         :: qflux,dqflux
  DOUBLE PRECISION, DIMENSION(ndim,3)      :: source_m,dsource_m
  DOUBLE PRECISION, DIMENSION(ndim1,3)     :: source_p,dsource_p
  DOUBLE PRECISION, DIMENSION(3,ndim)      :: flux_p,dflux_p
  DOUBLE PRECISION, DIMENSION(3,ndim1)     :: flux_m,dflux_m
  DOUBLE PRECISION, DIMENSION(ndim1,ndim)  :: amat_p_p,damat_p_p
  DOUBLE PRECISION, DIMENSION(ndim,ndim)   :: amat_p_m,damat_p_m
  DOUBLE PRECISION, DIMENSION(ndim,ndim1)  :: amat_m_m,damat_m_m
  DOUBLE PRECISION, DIMENSION(ndim1,ndim1) :: amat_m_p,damat_m_p
!
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: bvec,bvec_min
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: amat,bvec_lapack,prod
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: bvec_min_lapack
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: dummat1,dummat2
  DOUBLE PRECISION :: facnorm
!
! initialize
  ierr = 0
!
!
  ALLOCATE(amat(ndim,ndim),bvec_lapack(ndim,3),bvec_min_lapack(ndim,3))
  ALLOCATE(ipivot(ndim))
!
! coefficient matrix A 
  amat=0.d0
  DO i=1,ndim
    amat(i,i)=1.d0
  ENDDO
!
  ALLOCATE(dummat1(ndim,ndim1),dummat2(ndim1,ndim))
!
  DO m=0,nvel
    k=m*npass_l
    k1=m*npass_r
    dummat1(k+1:k+npass_l,:)                                    &
           =MATMUL(c_forward,amat_m_p(k1+1:k1+npass_r,:))
    dummat2(k1+1:k1+npass_r,:)                                  &
           =MATMUL(c_backward,amat_p_m(k+1:k+npass_l,:))
  ENDDO
!
  amat=amat-MATMUL(dummat1,dummat2)
!
  DEALLOCATE(dummat1,dummat2)
!
!
! source terms
! integrals of particle and heat fluxes and current over fieldline
! 
! gradient drive
! rhs of algebraic equation
!
  ALLOCATE(dummat1(ndim1,3),dummat2(ndim1,3))
!
  DO m=0,nvel
    k=m*npass_l
    k1=m*npass_r
    dummat1(k1+1:k1+npass_r,:)                                  &
           =MATMUL(c_backward,source_m(k+1:k+npass_l,:))
  ENDDO
!
  dummat2=source_p+MATMUL(amat_m_p,dummat1)
!
  DO m=0,nvel
    k=m*npass_l
    k1=m*npass_r
    bvec_lapack(k+1:k+npass_l,:)                                &
           =MATMUL(c_forward,dummat2(k1+1:k1+npass_r,:))
  ENDDO
!
  DEALLOCATE(dummat2)
!
!
! solution bvec -> $\ifour{\fvec}{+}{}{o+1}{g,l}$
!
  CALL gbsv(ndim,ndim,amat,ipivot,bvec_lapack,info)
!
  IF(info.NE.0) THEN
    ierr=2
    RETURN
  ENDIF
!
! bvec_min -> $\ifour{\fvec}{-}{}{o}{g,r}$
  bvec_min_lapack=source_m+MATMUL(amat_p_m,bvec_lapack)
!
! sources
  dsource_p=source_p+MATMUL(amat_p_p,bvec_lapack)
!
!
  DO m=0,nvel
    k=m*npass_l
    k1=m*npass_r
    dummat1(k1+1:k1+npass_r,:)                                  &
           =MATMUL(c_backward,bvec_min_lapack(k+1:k+npass_l,:))
  ENDDO
!
  dsource_m=source_m+MATMUL(amat_m_m,dummat1)
!
!
! integrals of particle flux and current over fieldline 
  dqflux=2.d0*qflux+MATMUL(flux_m,dummat1)+MATMUL(flux_p,bvec_lapack)
  DEALLOCATE(dummat1)
!
! source terms
! integrals of particle flux and current over fieldline
! 
!
! drive: f^+ from period n (from the left)
!
! rhs of equation
  DEALLOCATE(bvec_lapack,bvec_min_lapack)
  ALLOCATE(bvec_lapack(ndim,npass_l*(nvel+1)))
  ALLOCATE(bvec_min_lapack(ndim,npass_l*(nvel+1)))
!
  DO m=0,nvel
    k=m*npass_l
    k1=m*npass_r
    bvec_lapack(k+1:k+npass_l,:)                                 &
        =MATMUL(c_forward,amat_p_p(k1+1:k1+npass_r,:))
  ENDDO
!
! solution bvec_lapack -> $\ifour{\Fmat}{+}{}{o+1}{l,l}$
  CALL gbsv(ndim,ndim,amat,ipivot,bvec_lapack,info)
  IF(info.NE.0) THEN
    ierr=3
    RETURN
  ENDIF
!
! bvec_lapack = f^+ from period n+1 driven by f^+ from period n
! bvec_min_lapack = f^- from period n driven by f^+ from period n
!
  bvec_min_lapack=MATMUL(amat_p_m,bvec_lapack)
!
!
  ALLOCATE(dummat1(ndim1,ndim))
!
  DO m=0,nvel
    k=m*npass_l
    k1=m*npass_r
    dummat1(k1+1:k1+npass_r,:)                                   &
          =MATMUL(c_backward,bvec_min_lapack(k+1:k+npass_l,:))
  ENDDO
!
  damat_p_m=amat_p_m+MATMUL(amat_m_m,dummat1) 
!
  damat_p_p=MATMUL(amat_p_p,bvec_lapack)
!
  dflux_p=flux_p+MATMUL(flux_p,bvec_lapack)+MATMUL(flux_m,dummat1)
  DEALLOCATE(bvec_lapack,bvec_min_lapack,dummat1)

!
! drive: f^- from period n+1 (from right)
!
  ALLOCATE(bvec_lapack(ndim,ndim1))
  ALLOCATE(bvec_min_lapack(ndim,ndim1))
  ALLOCATE(dummat1(ndim,ndim1),dummat2(ndim1,ndim1))
!
! rhs of equation
  DO m=0,nvel
    k=m*npass_l
    k1=m*npass_r
    dummat1(k+1:k+npass_l,:)                                    &
           =MATMUL(c_forward,amat_m_p(k1+1:k1+npass_r,:))
    dummat2(k1+1:k1+npass_r,:)                                  &
           =MATMUL(c_backward,amat_m_m(k+1:k+npass_l,:))
  ENDDO
!
  bvec_lapack=MATMUL(dummat1,dummat2)
!
! solution bvec_lapack -> $\ifour{\Fmat}{+}{}{o+1}{r,l}$
  CALL gbsv(ndim,ndim,amat,ipivot,bvec_lapack,info)
  IF(info.NE.0) THEN
    ierr=4
    RETURN
  ENDIF
!
  bvec_min_lapack=amat_m_m+MATMUL(amat_p_m,bvec_lapack)
!
!
  DO m=0,nvel
    k=m*npass_l
    k1=m*npass_r
    dummat2(k1+1:k1+npass_r,:)                                  &
           =MATMUL(c_backward,bvec_min_lapack(k+1:k+npass_l,:))
  ENDDO
!
  damat_m_m=MATMUL(amat_m_m,dummat2)
!
!
  damat_m_p=amat_m_p+MATMUL(amat_p_p,bvec_lapack) 
!
  dflux_m=flux_m+MATMUL(flux_p,bvec_lapack)+MATMUL(flux_m,dummat2)
!
  DEALLOCATE(dummat1,dummat2)
  DEALLOCATE(amat,bvec_lapack,bvec_min_lapack)
  DEALLOCATE(ipivot)
!
! Correction of flux conservation (for Lorentz model only)
!
  IF(isw_lorentz.EQ.1) THEN
    DO i=1,npass_l
      facnorm=SUM(damat_p_p(:,i))+SUM(damat_p_m(:,i))
      damat_p_p(:,i)=damat_p_p(:,i)/facnorm
      damat_p_m(:,i)=damat_p_m(:,i)/facnorm
    ENDDO
    DO i=1,npass_r
      facnorm=SUM(damat_m_p(:,i))+SUM(damat_m_m(:,i))
      damat_m_p(:,i)=damat_m_p(:,i)/facnorm
      damat_m_m(:,i)=damat_m_m(:,i)/facnorm
    ENDDO
  ENDIF
!
  RETURN
END SUBROUTINE duplicate_ripple
!
! WINNY - THIS IS NEVER USED
!
!!$SUBROUTINE join_ripples_simple(ierr)
!!$
!!$!
!!$! Joins the results from the "old ripple", n, (naming o%...)
!!$! and the "new ripple", n+1, (naming n%...).
!!$! Puts the result in "old ripple" variables.
!!$!
!!$! WINNY:
!!$!
!!$! c_forward and c_backward should be allocatable arrays
!!$!  they will be passed most likely through the parameter list
!!$!  do you really want to make matmul in cases when binarysplit is not done?
!!$!
!!$!  c_forward  will have dim (n%p%npass_l , o%p%npass_r)
!!$!  c_backward will have dim (o%p%npass_r , n%p%npass_l)
!!$!
!!$! 
!!$  USE lapack_band
!!$  USE join_ripples_simple_mod
!!$  USE collisionality_mod, ONLY : isw_lorentz
!!$!
!!$  IMPLICIT NONE
!!$  INTEGER, PARAMETER :: dp = KIND(1.0d0)
!!$!
!!$  INTEGER, INTENT(out) :: ierr
!!$!
!!$  INTEGER :: ndim,ndim1,k,k1,i,i1
!!$  INTEGER :: lag,m
!!$  INTEGER :: info
!!$  INTEGER,          DIMENSION(:),   ALLOCATABLE :: ipivot
!!$!
!!$  integer, parameter :: joinripplessimple_write=0
!!$!
!!$  ! WINNY
!!$  ! made it compatible with propagator_mod
!!$  !
!!$  ! One can finally remove this c_forward and c_backward if we just
!!$  ! use o%p%cmat (for forward) and n%p%cmat (for backward)
!!$  !
!!$  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: c_forward,c_backward
!!$
!!$  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: bvec,bvec_min
!!$  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: amat,bvec_lapack,prod
!!$  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: bvec_min_lapack
!!$  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: dummat1,dummat2
!!$  DOUBLE PRECISION :: facnorm
!!$  !
!!$if(joinripplessimple_write.eq.1) then
!!$OPEN(111,file='o_flux_p.dat')
!!$WRITE(111,'(3(1x,e12.5))') (o_flux_p(:,i) ,i=1,o_l)
!!$CLOSE(111)
!!$OPEN(111,file='o_flux_m.dat')
!!$WRITE(111,'(3(1x,e12.5))') (o_flux_m(:,i) ,i=1,o_r)
!!$CLOSE(111)
!!$OPEN(111,file='o_source_p.dat')
!!$WRITE(111,'(3(1x,e12.5))') (o_source_p(i,:) ,i=1,o_r)
!!$CLOSE(111)
!!$OPEN(111,file='o_source_m.dat')
!!$WRITE(111,'(3(1x,e12.5))') (o_source_m(i,:) ,i=1,o_l)
!!$CLOSE(111)
!!$OPEN(111,file='o_qflux.dat')
!!$WRITE(111,'(4(1x,e12.5))') o_qflux(1:2,1:2)
!!$CLOSE(111)
!!$OPEN(111,file='o_amat_p_p.dat')
!!$DO i=1,o_r
!!$WRITE(111,*) o_amat_p_p(i,1:o_l)
!!$ENDDO
!!$CLOSE(111)
!!$OPEN(111,file='o_amat_p_m.dat')
!!$DO i=1,o_l
!!$WRITE(111,*) o_amat_p_m(i,1:o_l)
!!$ENDDO
!!$CLOSE(111)
!!$OPEN(111,file='o_amat_m_p.dat')
!!$DO i=1,o_r
!!$WRITE(111,*) o_amat_m_p(i,1:o_r)
!!$ENDDO
!!$CLOSE(111)
!!$OPEN(111,file='o_amat_m_m.dat')
!!$DO i=1,o_l
!!$WRITE(111,*) o_amat_m_m(i,1:o_r)
!!$ENDDO
!!$CLOSE(111)
!!$!
!!$OPEN(111,file='n_flux_p.dat')
!!$WRITE(111,'(3(1x,e12.5))') (n_flux_p(:,i) ,i=1,n_l)
!!$CLOSE(111)
!!$OPEN(111,file='n_flux_m.dat')
!!$WRITE(111,'(3(1x,e12.5))') (n_flux_m(:,i) ,i=1,n_r)
!!$CLOSE(111)
!!$OPEN(111,file='n_source_p.dat')
!!$WRITE(111,'(3(1x,e12.5))') (n_source_p(i,:) ,i=1,n_r)
!!$CLOSE(111)
!!$OPEN(111,file='n_source_m.dat')
!!$WRITE(111,'(3(1x,e12.5))') (n_source_m(i,:) ,i=1,n_l)
!!$CLOSE(111)
!!$OPEN(111,file='n_qflux.dat')
!!$WRITE(111,'(4(1x,e12.5))') n_qflux(1:2,1:2)
!!$CLOSE(111)
!!$OPEN(111,file='n_amat_p_p.dat')
!!$DO i=1,n_r
!!$WRITE(111,*) n_amat_p_p(i,1:n_l)
!!$ENDDO
!!$CLOSE(111)
!!$OPEN(111,file='n_amat_p_m.dat')
!!$DO i=1,n_l
!!$WRITE(111,*) n_amat_p_m(i,1:n_l)
!!$ENDDO
!!$CLOSE(111)
!!$OPEN(111,file='n_amat_m_p.dat')
!!$DO i=1,n_r
!!$WRITE(111,*) n_amat_m_p(i,1:n_r)
!!$ENDDO
!!$CLOSE(111)
!!$OPEN(111,file='n_amat_m_m.dat')
!!$DO i=1,n_l
!!$WRITE(111,*) n_amat_m_m(i,1:n_r)
!!$ENDDO
!!$CLOSE(111)
!!$endif
!!$!
!!$  ! initialize
!!$  ierr = 0
!!$
!!$  if(o_r.ne.n_l) then
!!$    ierr = 1
!!$    return
!!$  endif
!!$  allocate(c_forward(o_r,o_r),c_backward(o_r,o_r))
!!$  c_forward=0.d0
!!$  c_backward=0.d0
!!$  do i=1,o_r
!!$    c_forward(i,i)=1.d0
!!$    c_backward(i,i)=1.d0
!!$  enddo
!!$  !
!!$!
!!$  lag=n%p%nvelocity
!!$!
!!$  ndim=n_l*(lag+1)
!!$  ndim1=o_r*(lag+1)
!!$  ALLOCATE(amat(ndim,ndim),bvec_lapack(ndim,3),bvec_min_lapack(ndim,3))
!!$  ALLOCATE(ipivot(ndim))
!!$  !
!!$  ! coefficient matrix A 
!!$  amat=0.d0
!!$  DO i=1,ndim
!!$    amat(i,i)=1.d0
!!$  ENDDO
!!$!
!!$  ALLOCATE(dummat1(ndim,ndim1),dummat2(ndim1,ndim))
!!$!
!!$  DO m=0,lag
!!$    k=m*n_l
!!$    k1=m*o_r
!!$    dummat1(k+1:k+n_l,:)                                    &
!!$           =MATMUL(c_forward,o_amat_m_p(k1+1:k1+o_r,:))
!!$    dummat2(k1+1:k1+o_r,:)                                  &
!!$           =MATMUL(c_backward,n_amat_p_m(k+1:k+n_l,:))
!!$  ENDDO
!!$!
!!$  amat=amat-MATMUL(dummat1,dummat2)
!!$!
!!$  DEALLOCATE(dummat1,dummat2)
!!$!
!!$  !
!!$  ! source terms
!!$  ! integrals of particle and heat fluxes and current over fieldline
!!$  ! 
!!$  ! gradient drive
!!$  ! rhs of algebraic equation
!!$!
!!$!
!!$  ALLOCATE(dummat1(ndim1,3),dummat2(ndim1,3))
!!$!
!!$  DO m=0,lag
!!$    k=m*n_l
!!$    k1=m*o_r
!!$    dummat1(k1+1:k1+o_r,:)                                  &
!!$           =MATMUL(c_backward,n_source_m(k+1:k+n_l,:))
!!$  ENDDO
!!$!
!!$  dummat2=o_source_p+MATMUL(o_amat_m_p,dummat1)
!!$!
!!$  DO m=0,lag
!!$    k=m*n_l
!!$    k1=m*o_r
!!$    bvec_lapack(k+1:k+n_l,:)                                &
!!$           =MATMUL(c_forward,dummat2(k1+1:k1+o_r,:))
!!$  ENDDO
!!$!
!!$  DEALLOCATE(dummat2)
!!$!
!!$  !
!!$  ! solution bvec -> $\ifour{\fvec}{+}{}{o+1}{g,l}$
!!$!
!!$  CALL gbsv(ndim,ndim,amat,ipivot,bvec_lapack,info)
!!$!
!!$  IF(info.NE.0) THEN
!!$    ierr=2
!!$    RETURN
!!$  ENDIF
!!$  !
!!$  ! bvec_min -> $\ifour{\fvec}{-}{}{o}{g,r}$
!!$  bvec_min_lapack=n_source_m+MATMUL(n_amat_p_m,bvec_lapack)
!!$  !
!!$  ! sources
!!$  DEALLOCATE(o_source_p)
!!$  ALLOCATE(o_source_p(n_r*(lag+1),3))
!!$  o_source_p=n_source_p+MATMUL(n_amat_p_p,bvec_lapack)
!!$  !
!!$!
!!$  DO m=0,lag
!!$    k=m*n_l
!!$    k1=m*o_r
!!$    dummat1(k1+1:k1+o_r,:)                                  &
!!$           =MATMUL(c_backward,bvec_min_lapack(k+1:k+n_l,:))
!!$  ENDDO
!!$!
!!$  o_source_m=o_source_m+MATMUL(o_amat_m_m,dummat1)
!!$!
!!$  !
!!$  ! integrals of particle flux and current over fieldline 
!!$  o_qflux=o_qflux+n_qflux                                     &
!!$             +MATMUL(o_flux_m,dummat1)                            &
!!$             +MATMUL(n_flux_p,bvec_lapack)
!!$  DEALLOCATE(dummat1)
!!$  !
!!$  ! source terms
!!$  ! integrals of particle flux and current over fieldline
!!$  ! 
!!$  !
!!$  ! drive: f^+ from period n (from the left)
!!$  !
!!$  ! rhs of equation
!!$  DEALLOCATE(bvec_lapack,bvec_min_lapack)
!!$  ALLOCATE(bvec_lapack(ndim,o_l*(lag+1)))
!!$  ALLOCATE(bvec_min_lapack(ndim,o_l*(lag+1)))
!!$!
!!$  DO m=0,lag
!!$    k=m*n_l
!!$    k1=m*o_r
!!$    bvec_lapack(k+1:k+n_l,:)                                 &
!!$        =MATMUL(c_forward,o_amat_p_p(k1+1:k1+o_r,:))
!!$  ENDDO
!!$  !
!!$  ! solution bvec_lapack -> $\ifour{\Fmat}{+}{}{o+1}{l,l}$
!!$  CALL gbsv(ndim,ndim,amat,ipivot,bvec_lapack,info)
!!$  IF(info.NE.0) THEN
!!$    ierr=3
!!$    RETURN
!!$  ENDIF
!!$  !
!!$  ! bvec_lapack = f^+ from period n+1 driven by f^+ from period n
!!$  ! bvec_min_lapack = f^- from period n driven by f^+ from period n
!!$  !
!!$  bvec_min_lapack=MATMUL(n_amat_p_m,bvec_lapack)
!!$  !
!!$  !
!!$  ALLOCATE(dummat1(ndim1,o_l*(lag+1)))
!!$!
!!$  DO m=0,lag
!!$    k=m*n_l
!!$    k1=m*o_r
!!$    dummat1(k1+1:k1+o_r,:)                                   &
!!$          =MATMUL(c_backward,bvec_min_lapack(k+1:k+n_l,:))
!!$  ENDDO
!!$  !
!!$  o_amat_p_m=o_amat_p_m+MATMUL(o_amat_m_m,dummat1) 
!!$!
!!$  DEALLOCATE(o_amat_p_p)
!!$  ALLOCATE(o_amat_p_p(n_r*(lag+1),o_l*(lag+1)))
!!$  o_amat_p_p=MATMUL(n_amat_p_p,bvec_lapack)
!!$  !
!!$  o_flux_p=o_flux_p+MATMUL(n_flux_p,bvec_lapack)               &
!!$            +MATMUL(o_flux_m,dummat1)
!!$  DEALLOCATE(bvec_lapack,bvec_min_lapack,dummat1)
!!$  !
!!$  ! drive: f^- from period n+1 (from right)
!!$  !
!!$  ALLOCATE(bvec_lapack(ndim,n_r*(lag+1)))
!!$  ALLOCATE(bvec_min_lapack(ndim,n_r*(lag+1)))
!!$  ALLOCATE(dummat1(ndim,ndim1),dummat2(ndim1,n_r*(lag+1)))
!!$  !
!!$  ! rhs of equation
!!$  DO m=0,lag
!!$    k=m*n_l
!!$    k1=m*o_r
!!$    dummat1(k+1:k+n_l,:)                                    &
!!$           =MATMUL(c_forward,o_amat_m_p(k1+1:k1+o_r,:))
!!$    dummat2(k1+1:k1+o_r,:)                                  &
!!$           =MATMUL(c_backward,n_amat_m_m(k+1:k+n_l,:))
!!$  ENDDO
!!$  !
!!$  bvec_lapack=MATMUL(dummat1,dummat2)
!!$  !
!!$  ! solution bvec_lapack -> $\ifour{\Fmat}{+}{}{o+1}{r,l}$
!!$  CALL gbsv(ndim,ndim,amat,ipivot,bvec_lapack,info)
!!$  IF(info.NE.0) THEN
!!$    ierr=4
!!$    RETURN
!!$  ENDIF
!!$  !
!!$  bvec_min_lapack=n_amat_m_m+MATMUL(n_amat_p_m,bvec_lapack)
!!$  !
!!$  DEALLOCATE(amat)
!!$  ALLOCATE(amat(o_l*(lag+1),n_r*(lag+1)))
!!$!
!!$  DO m=0,lag
!!$    k=m*n_l
!!$    k1=m*o_r
!!$    dummat2(k1+1:k1+o_r,:)                                  &
!!$           =MATMUL(c_backward,bvec_min_lapack(k+1:k+n_l,:))
!!$  ENDDO
!!$!
!!$  amat=MATMUL(o_amat_m_m,dummat2)
!!$!
!!$  DEALLOCATE(o_amat_m_m)
!!$  ALLOCATE(o_amat_m_m(o_l*(lag+1),n_r*(lag+1)))
!!$  o_amat_m_m=amat
!!$  !
!!$  DEALLOCATE(o_amat_m_p)
!!$  ALLOCATE(o_amat_m_p(n_r*(lag+1),n_r*(lag+1)))
!!$  o_amat_m_p=n_amat_m_p+MATMUL(n_amat_p_p,bvec_lapack) 
!!$  !
!!$  n_flux_m=n_flux_m+MATMUL(n_flux_p,bvec_lapack)             &
!!$            +MATMUL(o_flux_m,dummat2)
!!$  !
!!$  DEALLOCATE(o_flux_m)
!!$  ALLOCATE(o_flux_m(3,n_r*(lag+1)))
!!$  DEALLOCATE(dummat1,dummat2)
!!$  !
!!$  o_flux_m=n_flux_m
!!$  !
!!$  o_r=n_r
!!$  !
!!$  DEALLOCATE(amat,bvec_lapack,bvec_min_lapack)
!!$  DEALLOCATE(ipivot)
!!$  ! WINNY
!!$  ! DEALLOCATE(c_forward,c_backward)
!!$  DEALLOCATE(c_forward)
!!$  DEALLOCATE(c_backward)
!!$!
!!$!
!!$! Correction of flux conservation (for Lorentz model only)
!!$!
!!$  IF(isw_lorentz.EQ.1) THEN
!!$    DO i=1,o_l
!!$      facnorm=SUM(o_amat_p_p(:,i))+SUM(o_amat_p_m(:,i))
!!$      o_amat_p_p(:,i)=o_amat_p_p(:,i)/facnorm
!!$      o_amat_p_m(:,i)=o_amat_p_m(:,i)/facnorm
!!$    ENDDO
!!$    DO i=1,o_r
!!$      facnorm=SUM(o_amat_m_p(:,i))+SUM(o_amat_m_m(:,i))
!!$      o_amat_m_p(:,i)=o_amat_m_p(:,i)/facnorm
!!$      o_amat_m_m(:,i)=o_amat_m_m(:,i)/facnorm
!!$    ENDDO
!!$  ENDIF
!!$!
!!$!
!!$if(joinripplessimple_write.eq.1) then
!!$OPEN(111,file='j_flux_p.dat')
!!$WRITE(111,'(3(1x,e12.5))') (o_flux_p(:,i) ,i=1,o_l)
!!$CLOSE(111)
!!$OPEN(111,file='j_flux_m.dat')
!!$WRITE(111,'(3(1x,e12.5))') (o_flux_m(:,i) ,i=1,o_r)
!!$CLOSE(111)
!!$OPEN(111,file='j_source_p.dat')
!!$WRITE(111,'(3(1x,e12.5))') (o_source_p(i,:) ,i=1,o_r)
!!$CLOSE(111)
!!$OPEN(111,file='j_source_m.dat')
!!$WRITE(111,'(3(1x,e12.5))') (o_source_m(i,:) ,i=1,o_l)
!!$CLOSE(111)
!!$OPEN(111,file='j_qflux.dat')
!!$WRITE(111,'(4(1x,e12.5))') o_qflux(1:2,1:2)
!!$CLOSE(111)
!!$OPEN(111,file='j_amat_p_p.dat')
!!$DO i=1,o_r
!!$WRITE(111,*) o_amat_p_p(i,1:o_l)
!!$ENDDO
!!$CLOSE(111)
!!$OPEN(111,file='j_amat_p_m.dat')
!!$DO i=1,o_l
!!$WRITE(111,*) o_amat_p_m(i,1:o_l)
!!$ENDDO
!!$CLOSE(111)
!!$OPEN(111,file='j_amat_m_p.dat')
!!$DO i=1,o_r
!!$WRITE(111,*) o_amat_m_p(i,1:o_r)
!!$ENDDO
!!$CLOSE(111)
!!$OPEN(111,file='j_amat_m_m.dat')
!!$DO i=1,o_l
!!$WRITE(111,*) o_amat_m_m(i,1:o_r)
!!$ENDDO
!!$CLOSE(111)
!!$!pause 'join_ripples_simple: written'
!!$endif
!!$!
!!$  !
!!$  RETURN
!!$END SUBROUTINE join_ripples_simple
!!$!
