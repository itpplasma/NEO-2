! THIS IS A NOBAND VERSION - MODIFIED CALL
!  modified for _g and _e drive
!  comments added by Winny
!  checked with write-up
!
SUBROUTINE join_ends(ierr)

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
! things with npart will not work properly
!  there is now o%p%npart and n%p%npart
!  but matrix sizes are anyway as they should be
!
!

  USE propagator_mod
  USE lapack_band
  use collisionality_mod, only : isw_lorentz, v_max_resolution, isw_energy, isw_integral
  !**********************************************************
  ! Definition of base functions
  !**********************************************************
  use rkstep_mod, only : lag
  use collop, only : collop_construct, collop_load
  use collop_compute, only : phi_exp, phi_prj
!
  IMPLICIT NONE
  INTEGER, PARAMETER :: dp = KIND(1.0d0)
!
  INTEGER, INTENT(out) :: ierr
!
  TYPE(propagator), POINTER            :: o
  TYPE(propagator), POINTER            :: n

!
  INTEGER :: ndim,i,i1,ntranseq,info,nl,nr,nts_r,nts_l,kl,kr,m,kr1,nvel
  INTEGER,          DIMENSION(:),   ALLOCATABLE :: ipivot,iminvec,imaxvec

  ! WINNY
  ! made it compatible with propagator_mod
  !
  ! One can finally remove this c_forward and c_backward if we just
  ! use o%p%cmat (for forward) and n%p%cmat (for backward)
  !
  DOUBLE PRECISION, DIMENSION(:,:), POINTER :: c_forward,c_backward
  !
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: totfun,totlev
  !
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: delta_eta_r,delta_eta_l
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: alam_l,alam_r
  !
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: amat,bvec_lapack,transmat
!
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: qflux,dqflux
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: source_m,dsource_m
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: source_p,dsource_p
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: flux_p,dflux_p
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: flux_m,dflux_m
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: amat_p_p,damat_p_p
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: amat_p_m,damat_p_m
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: amat_m_m,damat_m_m
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: amat_m_p,damat_m_p
  INTEGER :: ndupl,idupl
  DOUBLE PRECISION :: facnorm

  real(kind=dp) :: x, r

  !**********************************************************
  ! Automatically choose prop_finaljoin_mode
  !**********************************************************
  if (isw_energy .eq. 0 .and. isw_integral .eq. 0) then
     if (prop_finaljoin_mode .ne. 3) then
        write (*,*) "WARNING: Overwriting prop_finaljoin_mode to 3"
        prop_finaljoin_mode = 3
     end if
  else
     if (prop_finaljoin_mode .ne. 1) then
        !write (*,*) "WARNING: Overwriting prop_finaljoin_mode to 1"
        !prop_finaljoin_mode = 1
     end if
  end if

  if (prop_finaljoin_mode .eq. 2) then
     call collop_construct
     call collop_load
     nvel = lag
  end if
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
    PRINT *,'join_ends : sizes of c_forward and c_backward do not fit'
    ierr=1
    RETURN
  ELSEIF(n%p%npass_l.NE.i) THEN
    PRINT *,'join_ends : size 1 of c_forward differs from new npass_l'
    ierr=1
    RETURN
  ELSEIF(o%p%npass_r.NE.i1) THEN
    PRINT *,'join_ends : size 2 of c_forward differs from old npass_r'
    ierr=1
    RETURN
  ELSE
    ierr=0
  ENDIF
!
  nl=o%p%npass_l
  nr=o%p%npass_r
PRINT *,'nl,nr = ',nl,nr
  nvel=o%p%nvelocity
  nts_l=nl*(nvel+1)
  nts_r=nr*(nvel+1)
!
!
! Correction of sorce and flux symmetry due to Pfirsch-Schlueter current closure
! condition and asymmetry of electric drive (for Lorentz model only)
!
  IF(isw_lorentz.EQ.1) THEN
!
    ALLOCATE(delta_eta_r(nr),delta_eta_l(nl))
!
    delta_eta_l(1:nl-1)=o%p%eta_l(1:nl-1)-o%p%eta_l(0:nl-2)
    delta_eta_l(nl) = o%p%eta_boundary_l - o%p%eta_l(nl-1)
    delta_eta_r(1:nr-1)=o%p%eta_r(1:nr-1)-o%p%eta_r(0:nr-2)
    delta_eta_r(nr) = o%p%eta_boundary_r - o%p%eta_r(nr-1)
!
    DO i=1,3
!
      facnorm=(SUM(o%p%source_p(:,i))+SUM(o%p%source_m(:,i)))          &
             *0.5d0/o%p%eta_boundary_r
      o%p%source_p(:,i)=o%p%source_p(:,i)-delta_eta_r*facnorm
      o%p%source_m(:,i)=o%p%source_m(:,i)-delta_eta_l*facnorm
!
      facnorm=(SUM(o%p%flux_p(i,:)*delta_eta_l)                        &
             +SUM(o%p%flux_m(i,:)*delta_eta_r))                        &
             *0.5d0/o%p%eta_boundary_r
      o%p%flux_p(i,:)=o%p%flux_p(i,:)-facnorm
      o%p%flux_m(i,:)=o%p%flux_m(i,:)-facnorm
!
    ENDDO
!
  ENDIF
!
  ndupl=0
!
  IF(ndupl.GT.0) THEN
    ALLOCATE(qflux(3,3),dqflux(3,3))
    ALLOCATE(source_m(nts_l,3),dsource_m(nts_l,3))
    ALLOCATE(source_p(nts_r,3),dsource_p(nts_r,3))
    ALLOCATE(flux_p(3,nts_l),dflux_p(3,nts_l))
    ALLOCATE(flux_m(3,nts_r),dflux_m(3,nts_r))
    ALLOCATE(amat_p_p(nts_r,nts_l),damat_p_p(nts_r,nts_l))
    ALLOCATE(amat_p_m(nts_l,nts_l),damat_p_m(nts_l,nts_l))
    ALLOCATE(amat_m_m(nts_l,nts_r),damat_m_m(nts_l,nts_r))
    ALLOCATE(amat_m_p(nts_r,nts_r),damat_m_p(nts_r,nts_r))
!
    qflux=o%p%qflux
    source_m=o%p%source_m
    source_p=o%p%source_p
    flux_p=o%p%flux_p
    flux_m=o%p%flux_m
    amat_p_p=o%p%amat_p_p
    amat_p_m=o%p%amat_p_m
    amat_m_m=o%p%amat_m_m
    amat_m_p=o%p%amat_m_p
!
    DO idupl=1,ndupl
!
      CALL  duplicate_ripple(nvel,nts_l,nts_r,nl,nr,                          &
                             c_forward,c_backward,                           &
                             amat_p_p,amat_p_m,amat_m_m,amat_m_p,            &
                             source_p,source_m,flux_p,flux_m,qflux,          &
                             damat_p_p,damat_p_m,damat_m_m,damat_m_p,        &
                             dsource_p,dsource_m,dflux_p,dflux_m,dqflux,     &
                             ierr)
!
      PRINT *,'Duplication ',idupl,' :'
      PRINT *,'before : ',qflux(1:2,1:2)
      PRINT *,'after  : ',dqflux(1:2,1:2)
!
! GOTO 1
! #if !defined(MPI_SUPPORT)
OPEN(111,file='delta_eta_l.dat')
DO i=1,nl
WRITE(111,*) delta_eta_l(i)
ENDDO
CLOSE(111)
OPEN(111,file='delta_eta_r.dat')
DO i=1,nr
WRITE(111,*) delta_eta_r(i)
ENDDO
CLOSE(111)
!
OPEN(111,file='flux_p.dat')
DO i=1,nl
WRITE(111,*) flux_p(:,i)
ENDDO
CLOSE(111)
OPEN(111,file='flux_m.dat')
DO i=1,nr
WRITE(111,*) flux_m(:,i)
ENDDO
CLOSE(111)
OPEN(111,file='source_p.dat')
DO i=1,nr
WRITE(111,*) source_p(i,:)
ENDDO
CLOSE(111)
OPEN(111,file='source_m.dat')
DO i=1,nl
WRITE(111,*) source_m(i,:)
ENDDO
CLOSE(111)
!
OPEN(111,file='amat_p_p.dat')
DO i=1,nr
WRITE(111,*) amat_p_p(i,1:nl)
ENDDO
CLOSE(111)
OPEN(111,file='amat_p_m.dat')
DO i=1,nl
WRITE(111,*) amat_p_m(i,1:nl)
ENDDO
CLOSE(111)
OPEN(111,file='amat_m_p.dat')
DO i=1,nr
WRITE(111,*) amat_m_p(i,1:nr)
ENDDO
CLOSE(111)
OPEN(111,file='amat_m_m.dat')
DO i=1,nl
WRITE(111,*) amat_m_m(i,1:nr)
ENDDO
CLOSE(111)
! #endif

!PAUSE 'written'
1 CONTINUE

      PRINT *,' '
!
      qflux=dqflux
      source_m=dsource_m
      source_p=dsource_p
      flux_p=dflux_p
      flux_m=dflux_m
      amat_p_p=damat_p_p
      amat_p_m=damat_p_m
      amat_m_m=damat_m_m
      amat_m_p=damat_m_p
!
    ENDDO
!
    o%p%qflux=dqflux
    o%p%source_m=dsource_m
    o%p%source_p=dsource_p
    o%p%flux_p=dflux_p
    o%p%flux_m=dflux_m
    o%p%amat_p_p=damat_p_p
    o%p%amat_p_m=damat_p_m
    o%p%amat_m_m=damat_m_m
    o%p%amat_m_p=damat_m_p
!
  ENDIF
!
!
!
  ndim=nts_l+nts_r
  ntranseq=3
  ALLOCATE(amat(ndim,ndim),bvec_lapack(ndim,ntranseq),ipivot(ndim),totlev(ndim))
  ALLOCATE(transmat(ntranseq,ntranseq),totfun(ntranseq))
  !
  ! coefficient matrix A
!
  amat=0.d0
!
  DO i=1,ndim
    amat(i,i)=1.d0
  ENDDO
!
!  amat(1:nl,1:nl)          =amat(1:nl,1:nl)                                 &
!                           -MATMUL(c_forward,o%p%amat_p_p)
!  amat(1:nl,nl+1:ndim)     =amat(1:nl,nl+1:ndim)                            &
!                           -MATMUL(c_forward,o%p%amat_m_p)
!  amat(nl+1:ndim,1:nl)     =amat(nl+1:ndim,1:nl)                            &
!                           -MATMUL(c_backward,o%p%amat_p_m)
!  amat(nl+1:ndim,nl+1:ndim)=amat(nl+1:ndim,nl+1:ndim)                       &
!                           -MATMUL(c_backward,o%p%amat_m_m)
!
  !
  ! source terms
  ! integrals of particle flux and current over fieldline
  !
  ! drive
  ! rhs of algebraic equation
!   bvec_lapack(1:nl,1:3)     =MATMUL(c_forward,o%p%source_p)
!   bvec_lapack(nl+1:ndim,1:3)=MATMUL(c_backward,o%p%source_m)
  !
  DO m=0,nvel
    kl=m*nl
    kr=m*nr
    kr1=m*nr+nts_l
    amat(kl+1:kl+nl,1:nts_l)=amat(kl+1:kl+nl,1:nts_l)                          &
                     -MATMUL(c_forward,o%p%amat_p_p(kr+1:kr+nr,1:nts_l))
    amat(kl+1:kl+nl,nts_l+1:ndim)=amat(kl+1:kl+nl,nts_l+1:ndim)                &
                     -MATMUL(c_forward,o%p%amat_m_p(kr+1:kr+nr,1:nts_r))
    amat(kr1+1:kr1+nr,1:nts_l)=amat(kr1+1:kr1+nr,1:nts_l)                      &
                     -MATMUL(c_backward,o%p%amat_p_m(kl+1:kl+nl,1:nts_l))
    amat(kr1+1:kr1+nr,nts_l+1:ndim)=amat(kr1+1:kr1+nr,nts_l+1:ndim)            &
                     -MATMUL(c_backward,o%p%amat_m_m(kl+1:kl+nl,1:nts_r))
    bvec_lapack(kl+1:kl+nl,1:3)                                                &
                     =MATMUL(c_forward,o%p%source_p(kr+1:kr+nr,1:3))
    bvec_lapack(kr1+1:kr1+nr,1:3)                                              &
                     =MATMUL(c_backward,o%p%source_m(kl+1:kl+nl,1:3))
  ENDDO
  !
! #if !defined(MPI_SUPPORT)
OPEN(751,file='amat_before.dat')
OPEN(752,file='bvec_before.dat')
DO i=1,ndim
WRITE(751,*) amat(i,:)
WRITE(752,*) bvec_lapack(i,:)
ENDDO
CLOSE(751)
CLOSE(752)
! #endif
! System is overdefined, one equation is redundant. Replace the redundant
! equation with zero perturbed density condition at the joining point.
!  amat(nl,:)=0.d0
!  amat(nl,nl)=1.d0
!  bvec_lapack(nl,:)=0.d0
!!  amat(ndim,:)=0.d0
!!  amat(ndim,ndim)=1.d0
!!  bvec_lapack(ndim,:)=0.d0
!!!  amat(1,:)=0.d0
!!!  amat(1,1)=1.d0
!!!  bvec_lapack(1,:)=0.d0
!!!!  amat(1,:)=0.d0
!!!!  amat(1,1:nl)=1.d0
!!!!  amat(1,nts_l+1:nts_l+nr)=1.d0
!!!!  bvec_lapack(1,:)=0.d0
if (prop_finaljoin_mode .eq. 0) then
   amat(nl,:)=amat(nl,:)-amat(nts_l+nr,:)
   bvec_lapack(nl,:)=bvec_lapack(nl,:)-bvec_lapack(nts_l+nr,:)
   !  amat(nts_l+nr,1:nl)=0.d0
   amat(nts_l+nr,1:nl)=1.d0
   amat(nts_l+nr,nts_l+1:nts_l+nr)=1.d0
   bvec_lapack(nts_l+nr,:)=0.d0
   if(nvel.gt.0 .and. isw_integral .eq. 1) then
      amat(2*nl,:)=amat(2*nl,:)-amat(nts_l+2*nr,:)
      bvec_lapack(2*nl,:)=bvec_lapack(2*nl,:)-bvec_lapack(nts_l+2*nr,:)
      !  amat(nts_l+2*nr,1:nl)=0.d0
      amat(nts_l+2*nr,nl+1:2*nl)=1.d0
      amat(nts_l+2*nr,nts_l+nr+1:nts_l+2*nr)=1.d0
      bvec_lapack(nts_l+2*nr,:)=0.d0
   endif
elseif (prop_finaljoin_mode .eq. 1) then
   amat(nl,:)=amat(nl,:)-amat(nts_l+nr,:)
   bvec_lapack(nl,:)=bvec_lapack(nl,:)-bvec_lapack(nts_l+nr,:)
   amat(nts_l+nr,:)=0.d0
   amat(nts_l+nr,1:nl)=1.d0
   amat(nts_l+nr,nts_l+1:nts_l+nr)=1.d0
   bvec_lapack(nts_l+nr,:)=0.d0
   if(nvel.gt.0 .and. isw_integral .eq. 1) then
      amat(2*nl,:)=amat(2*nl,:)-amat(nts_l+2*nr,:)
      bvec_lapack(2*nl,:)=bvec_lapack(2*nl,:)-bvec_lapack(nts_l+2*nr,:)
      amat(nts_l+2*nr,:)=0.d0
      amat(nts_l+2*nr,nl+1:2*nl)=1.d0
      amat(nts_l+2*nr,nts_l+nr+1:nts_l+2*nr)=1.d0
      bvec_lapack(nts_l+2*nr,:)=0.d0
   endif
elseif (prop_finaljoin_mode .eq. 2) then
   !**********************************************************
   ! Additional conditions:
   ! Should be deactivated if particle sink in ripple_solver is active
   !**********************************************************
   amat(nl,:)=amat(nl,:)-amat(nts_l+nr,:)
   bvec_lapack(nl,:)=bvec_lapack(nl,:)-bvec_lapack(nts_l+nr,:)
   amat(nts_l+nr,:)=0.d0
   !
   x  = 0d0                                       !<=new
   do m=0,nvel                                    !<=new
      r  = phi_exp(m, x)                           !<=new
      amat(nts_l+nr,nl*m+1:nl*m+nl)=r              !<=new
      amat(nts_l+nr,nts_l+nr*m+1:nts_l+nr*m+nr)=r  !<=new
   enddo                                          !<=new
   !
   bvec_lapack(nts_l+nr,:)=0.d0
   if(nvel.gt.0 .and. isw_integral .eq. 1) then
      amat(2*nl,:)=amat(2*nl,:)-amat(nts_l+2*nr,:)
      bvec_lapack(2*nl,:)=bvec_lapack(2*nl,:)-bvec_lapack(nts_l+2*nr,:)
      amat(nts_l+2*nr,:)=0.d0
      !
      x  = v_max_resolution                            !<=new
      do m=0,nvel                                      !<=new
         r  = phi_exp(m, x)                             !<=new
         amat(nts_l+2*nr,nl*m+1:nl*m+nl)=r              !<=new
         amat(nts_l+2*nr,nts_l+nr*m+1:nts_l+nr*m+nr)=r  !<=new
      enddo                                            !<=new
      !
      bvec_lapack(nts_l+2*nr,:)=0.d0
   endif
   !*****************************************************************
elseif (prop_finaljoin_mode .eq. 3) then                                           !<=new 05.04.18
   do m=0,nvel                                                                     !<=new 05.04.18
     amat(nl*(m+1),:)=amat(nl*(m+1),:)-amat(nts_l+nr*(m+1),:)                      !<=new 05.04.18
     bvec_lapack(nl*(m+1),:)=bvec_lapack(nl*(m+1),:)-bvec_lapack(nts_l+nr*(m+1),:) !<=new 05.04.18
     amat(nts_l+nr*(m+1),:)=0.d0                                                   !<=new 05.04.18
     amat(nts_l+nr*(m+1),nl*m+1:nl*(m+1))=1.d0                                     !<=new 05.04.18
     amat(nts_l+nr*(m+1),nts_l+nr*m+1:nts_l+nr*(m+1))=1.d0                         !<=new 05.04.18
     bvec_lapack(nts_l+nr*(m+1),:)=0.d0                                            !<=new 05.04.18
   enddo                                                                           !<=new 05.04.18
end if
! #if !defined(MPI_SUPPORT)
OPEN(751,file='amat_after.dat')
OPEN(752,file='bvec_after.dat')
DO i=1,ndim
WRITE(751,*) amat(i,:)
WRITE(752,*) bvec_lapack(i,:)
ENDDO
CLOSE(751)
CLOSE(752)
!#endif
  !
  CALL gbsv(ndim,ndim,amat,ipivot,bvec_lapack,info)
  IF(info.NE.0) THEN
    ierr=2
    RETURN
  ENDIF
  !
!  ! solution with zero average
!  ALLOCATE(delta_eta_r(nr),alam_r(nr))
!  DO i = 1, nr - 1
!     delta_eta_r(i) = o%p%eta_r(i) - o%p%eta_r(i-1)
!     alam_r(i)=sqrt(1.d0-0.5d0*(o%p%eta_r(i)+o%p%eta_r(i-1))/o%p%eta_boundary_r)
!  END DO
!  delta_eta_r(nr) = o%p%eta_boundary_r - o%p%eta_r(nr-1)
!  alam_r(nr)=sqrt(1.d0-0.5d0*(o%p%eta_boundary_r                              &
!            +o%p%eta_r(nr-1))/o%p%eta_boundary_r)
!
!  ALLOCATE(delta_eta_l(nl),alam_l(nl))
!  DO i = 1, nl - 1
!     delta_eta_l(i) = o%p%eta_l(i) - o%p%eta_l(i-1)
!     alam_l(i)=sqrt(1.d0-0.5d0*(o%p%eta_l(i)+o%p%eta_l(i-1))/o%p%eta_boundary_l)
!  END DO
!  delta_eta_l(nl) = o%p%eta_boundary_l - o%p%eta_l(nl-1)
!  alam_l(nl)=sqrt(1.d0-0.5d0*(o%p%eta_boundary_l                              &
!            +o%p%eta_l(nl-1))/o%p%eta_boundary_l)
!  !
!  totlev(1:nl)=delta_eta_l(1:nl)
!!  totlev(nl+1:ndim)=delta_eta_r(1:nr)
!  totlev(nl+1:nl+nr)=delta_eta_r(1:nr)
!  totfun=sum(bvec_lapack,1)/sum(totlev)
!  do i=1,ntranseq
!    bvec_lapack(:,i)=bvec_lapack(:,i)-totlev(:)*totfun(i)
!  enddo
! #if !defined(MPI_SUPPORT)
OPEN(752,file='bvec_solution.dat')
DO i=1,ndim
WRITE(752,*) bvec_lapack(i,:)
ENDDO
CLOSE(752)
! #endif
  !
  ! Now sources from outgoing become incoming, dimensions switch places !!!
  !
  DEALLOCATE(o%p%source_p,o%p%source_m)
!  allocate(o%p%source_p(nl,3),o%p%source_m(nr,3))
  ALLOCATE(o%p%source_p(nts_l,3),o%p%source_m(nts_r,3))
  !
  o%p%source_p=bvec_lapack(1:nts_l,1:3)
  o%p%source_m=bvec_lapack(nts_l+1:ndim,1:3)
  !
  ! integrals of particle flux and current over fieldline
  transmat=o%p%qflux+MATMUL(o%p%flux_p,o%p%source_p)           &
                    +MATMUL(o%p%flux_m,o%p%source_m)

! #if !defined(MPI_SUPPORT)
OPEN(751,file='fin_source_p.dat')
OPEN(752,file='fin_flux_p.dat')
DO i=1,nts_l
WRITE(751,*) o%p%source_p(i,:)
WRITE(752,*) o%p%flux_p(:,i)
ENDDO
CLOSE(751)
CLOSE(752)
OPEN(751,file='fin_source_m.dat')
OPEN(752,file='fin_flux_m.dat')
DO i=1,nts_r
WRITE(751,*) o%p%source_m(i,:)
WRITE(752,*) o%p%flux_m(:,i)
ENDDO
CLOSE(751)
CLOSE(752)
OPEN(751,file='fin_dims.dat')
WRITE(751,*) nts_l,nts_r
CLOSE(751)
! #endif
  !
  o%p%qflux=transmat
  !
  DEALLOCATE(amat,bvec_lapack,ipivot)
  NULLIFY(c_forward)
  NULLIFY(c_backward)
  !
  !
!  DEALLOCATE (delta_eta_r)
!  DEALLOCATE (delta_eta_l)
!  DEALLOCATE (alam_l,alam_r)

  RETURN
END SUBROUTINE join_ends
