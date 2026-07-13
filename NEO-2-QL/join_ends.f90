! THIS IS A NOBAND VERSION - MODIFIED CALL
!  modified for _g and _e drive
!  comments added by Winny
!  checked with write-up
SUBROUTINE join_ends(ierr)
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

  USE propagator_mod
  USE lapack_band
  USE collisionality_mod, ONLY : isw_lorentz
  USE join_diagnostics_mod, ONLY : join_end_trace_enabled, &
                                   record_join_end_compatibility

  IMPLICIT NONE
  INTEGER, PARAMETER :: dp = KIND(1.0d0)

  INTEGER, INTENT(out) :: ierr

  TYPE(propagator), POINTER            :: o
  TYPE(propagator), POINTER            :: n


  INTEGER :: ndim,i,i1,ntranseq,info,nl,nr,nts_r,nts_l,kl,kr,m,kr1,nvel
  INTEGER :: idiag,jdiag,nconstraints,trace_info
  INTEGER,          DIMENSION(:),   ALLOCATABLE :: ipivot,iminvec,imaxvec
  INTEGER,          DIMENSION(:),   ALLOCATABLE :: trace_pivot

  ! WINNY
  ! made it compatible with propagator_mod
  !
  ! One can finally remove this c_forward and c_backward if we just
  ! use o%p%cmat (for forward) and n%p%cmat (for backward)
  DOUBLE PRECISION, DIMENSION(:,:), POINTER :: c_forward,c_backward

  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: totfun,totlev

  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: delta_eta_r,delta_eta_l
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: alam_l,alam_r

  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: amat,bvec_lapack,transmat
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: compatibility,dropped_p
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: dropped_m,left_null
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: compatibility_scale
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: dropped_p_scale
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: dropped_m_scale
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: right_null,periodic_matrix
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: periodic_rhs,trace_matrix
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: left_residual,right_residual
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: left_scale,right_scale

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
  DOUBLE PRECISION :: measure_sum,source_after(3),source_before(3)
  DOUBLE PRECISION :: source_factor(3),source_scale(3),transfer_error(2)
  LOGICAL :: trace_join_end

  ! initialize
  ierr = 0
  trace_join_end = join_end_trace_enabled(ierr)
  IF(ierr.NE.0) RETURN
  o => prop_c%prev
  n => prop_c

  c_forward  => o%p%cmat
  c_backward => n%p%cmat

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

  nl=o%p%npass_l
  nr=o%p%npass_r
PRINT *,'nl,nr = ',nl,nr
  nvel=o%p%nvelocity
  nts_l=nl*(nvel+1)
  nts_r=nr*(nvel+1)

  IF(trace_join_end) THEN
    source_factor=0.d0
    source_before=SUM(o%p%source_p,DIM=1)+SUM(o%p%source_m,DIM=1)
    source_scale=SUM(ABS(o%p%source_p),DIM=1) &
                +SUM(ABS(o%p%source_m),DIM=1)
    source_scale=MAX(source_scale,TINY(1.d0))
    source_after=source_before
    measure_sum=0.d0
    transfer_error=0.d0
    DO idiag=1,SIZE(c_forward,1)
      DO jdiag=1,SIZE(c_forward,2)
        IF(idiag.EQ.jdiag) THEN
          transfer_error(1)=MAX(transfer_error(1), &
                                ABS(c_forward(idiag,jdiag)-1.d0))
        ELSE
          transfer_error(1)=MAX(transfer_error(1), &
                                ABS(c_forward(idiag,jdiag)))
        ENDIF
      ENDDO
    ENDDO
    DO idiag=1,SIZE(c_backward,1)
      DO jdiag=1,SIZE(c_backward,2)
        IF(idiag.EQ.jdiag) THEN
          transfer_error(2)=MAX(transfer_error(2), &
                                ABS(c_backward(idiag,jdiag)-1.d0))
        ELSE
          transfer_error(2)=MAX(transfer_error(2), &
                                ABS(c_backward(idiag,jdiag)))
        ENDIF
      ENDDO
    ENDDO
  ENDIF


! Correction of sorce and flux symmetry due to Pfirsch-Schlueter current closure
! condition and asymmetry of electric drive (for Lorentz model only)

  IF(isw_lorentz.EQ.1) THEN

    ALLOCATE(delta_eta_r(nr),delta_eta_l(nl))

    delta_eta_l(1:nl-1)=o%p%eta_l(1:nl-1)-o%p%eta_l(0:nl-2)
    delta_eta_l(nl) = o%p%eta_boundary_l - o%p%eta_l(nl-1)
    delta_eta_r(1:nr-1)=o%p%eta_r(1:nr-1)-o%p%eta_r(0:nr-2)
    delta_eta_r(nr) = o%p%eta_boundary_r - o%p%eta_r(nr-1)

    IF(trace_join_end) measure_sum=SUM(delta_eta_l)+SUM(delta_eta_r)

    DO i=1,3

      facnorm=(SUM(o%p%source_p(:,i))+SUM(o%p%source_m(:,i)))          &
             *0.5d0/o%p%eta_boundary_r
      IF(trace_join_end) source_factor(i)=facnorm
      o%p%source_p(:,i)=o%p%source_p(:,i)-delta_eta_r*facnorm
      o%p%source_m(:,i)=o%p%source_m(:,i)-delta_eta_l*facnorm

      facnorm=(SUM(o%p%flux_p(i,:)*delta_eta_l)                        &
             +SUM(o%p%flux_m(i,:)*delta_eta_r))                        &
             *0.5d0/o%p%eta_boundary_r
      o%p%flux_p(i,:)=o%p%flux_p(i,:)-facnorm
      o%p%flux_m(i,:)=o%p%flux_m(i,:)-facnorm

    ENDDO

    IF(trace_join_end) &
      source_after=SUM(o%p%source_p,DIM=1)+SUM(o%p%source_m,DIM=1)

  ENDIF

  ndupl=0

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

    qflux=o%p%qflux
    source_m=o%p%source_m
    source_p=o%p%source_p
    flux_p=o%p%flux_p
    flux_m=o%p%flux_m
    amat_p_p=o%p%amat_p_p
    amat_p_m=o%p%amat_p_m
    amat_m_m=o%p%amat_m_m
    amat_m_p=o%p%amat_m_p

    DO idupl=1,ndupl

      CALL  duplicate_ripple(nvel,nts_l,nts_r,nl,nr,                          &
                             c_forward,c_backward,                           &
                             amat_p_p,amat_p_m,amat_m_m,amat_m_p,            &
                             source_p,source_m,flux_p,flux_m,qflux,          &
                             damat_p_p,damat_p_m,damat_m_m,damat_m_p,        &
                             dsource_p,dsource_m,dflux_p,dflux_m,dqflux,     &
                             ierr)

      PRINT *,'Duplication ',idupl,' :'
      PRINT *,'before : ',qflux(1:2,1:2)
      PRINT *,'after  : ',dqflux(1:2,1:2)

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

1 CONTINUE

      PRINT *,' '

      qflux=dqflux
      source_m=dsource_m
      source_p=dsource_p
      flux_p=dflux_p
      flux_m=dflux_m
      amat_p_p=damat_p_p
      amat_p_m=damat_p_m
      amat_m_m=damat_m_m
      amat_m_p=damat_m_p

    ENDDO

    o%p%qflux=dqflux
    o%p%source_m=dsource_m
    o%p%source_p=dsource_p
    o%p%flux_p=dflux_p
    o%p%flux_m=dflux_m
    o%p%amat_p_p=damat_p_p
    o%p%amat_p_m=damat_p_m
    o%p%amat_m_m=damat_m_m
    o%p%amat_m_p=damat_m_p

  ENDIF


  ndim=nts_l+nts_r
  ntranseq=3
  ALLOCATE(amat(ndim,ndim),bvec_lapack(ndim,ntranseq),ipivot(ndim),totlev(ndim))
  ALLOCATE(transmat(ntranseq,ntranseq),totfun(ntranseq))

  ! coefficient matrix A 

  amat=0.d0

  DO i=1,ndim
    amat(i,i)=1.d0
  ENDDO

  ! source terms
  ! integrals of particle flux and current over fieldline
  ! 
  ! drive
  ! rhs of algebraic equation
  !   bvec_lapack(1:nl,1:3)     =MATMUL(c_forward,o%p%source_p)
  !   bvec_lapack(nl+1:ndim,1:3)=MATMUL(c_backward,o%p%source_m)
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

  IF(trace_join_end) THEN
    nconstraints=1
    IF(nvel.GT.0) nconstraints=2
    ALLOCATE(periodic_matrix(ndim,ndim),periodic_rhs(ndim,ntranseq))
    ALLOCATE(left_null(ndim,nconstraints),right_null(ndim,nconstraints))
    ALLOCATE(compatibility(nconstraints,ntranseq))
    ALLOCATE(dropped_p(nconstraints,ntranseq))
    ALLOCATE(dropped_m(nconstraints,ntranseq))
    ALLOCATE(compatibility_scale(nconstraints,ntranseq))
    ALLOCATE(dropped_p_scale(nconstraints,ntranseq))
    ALLOCATE(dropped_m_scale(nconstraints,ntranseq))
    ALLOCATE(left_residual(nconstraints),right_residual(nconstraints))
    ALLOCATE(left_scale(nconstraints),right_scale(nconstraints))
    ALLOCATE(trace_matrix(ndim,ndim),trace_pivot(ndim))
    periodic_matrix=amat
    periodic_rhs=bvec_lapack
    left_null=0.d0
    DO m=0,nconstraints-1
      left_null(m*nl+1:(m+1)*nl,m+1)=1.d0
      left_null(nts_l+m*nr+1:nts_l+(m+1)*nr,m+1)=1.d0
    ENDDO
    compatibility=MATMUL(TRANSPOSE(left_null),periodic_rhs)
    DO m=1,nconstraints
      DO i=1,ntranseq
        compatibility_scale(m,i)=MAX(SUM(ABS(periodic_rhs(:,i)) &
                                      *ABS(left_null(:,m))),TINY(1.d0))
      ENDDO
      left_residual(m)=MAXVAL(ABS(MATMUL(TRANSPOSE(periodic_matrix), &
                                       left_null(:,m))))
      left_scale(m)=MAX(MAXVAL(MATMUL(TRANSPOSE(ABS(periodic_matrix)), &
                                     ABS(left_null(:,m)))),TINY(1.d0))
    ENDDO
  ENDIF

OPEN(751,file='amat_before.dat')
OPEN(752,file='bvec_before.dat')
DO i=1,ndim
WRITE(751,*) amat(i,:)
WRITE(752,*) bvec_lapack(i,:)
ENDDO
CLOSE(751)
CLOSE(752)

  amat(nl,:)=amat(nl,:)-amat(nts_l+nr,:)
  bvec_lapack(nl,:)=bvec_lapack(nl,:)-bvec_lapack(nts_l+nr,:)

  amat(nts_l+nr,1:nl)=1.d0
  amat(nts_l+nr,nts_l+1:nts_l+nr)=1.d0
  bvec_lapack(nts_l+nr,:)=0.d0
if(nvel.gt.0) then
  amat(2*nl,:)=amat(2*nl,:)-amat(nts_l+2*nr,:)
  bvec_lapack(2*nl,:)=bvec_lapack(2*nl,:)-bvec_lapack(nts_l+2*nr,:)

  amat(nts_l+2*nr,nl+1:2*nl)=1.d0
  amat(nts_l+2*nr,nts_l+nr+1:nts_l+2*nr)=1.d0
  bvec_lapack(nts_l+2*nr,:)=0.d0
endif

  IF(trace_join_end) THEN
    trace_matrix=periodic_matrix
    right_null=0.d0
    trace_matrix(nts_l+nr,:)=0.d0
    trace_matrix(nts_l+nr,1:nl)=1.d0
    trace_matrix(nts_l+nr,nts_l+1:nts_l+nr)=1.d0
    right_null(nts_l+nr,1)=1.d0
    IF(nconstraints.GT.1) THEN
      trace_matrix(nts_l+2*nr,:)=0.d0
      trace_matrix(nts_l+2*nr,nl+1:2*nl)=1.d0
      trace_matrix(nts_l+2*nr,nts_l+nr+1:nts_l+2*nr)=1.d0
      right_null(nts_l+2*nr,2)=1.d0
    ENDIF
    CALL gbsv(ndim,ndim,trace_matrix,trace_pivot,right_null,trace_info)
    IF(trace_info.NE.0) THEN
      ierr=9
      RETURN
    ENDIF
    DO m=1,nconstraints
      right_residual(m)=MAXVAL(ABS(MATMUL(periodic_matrix, &
                                        right_null(:,m))))
      right_scale(m)=MAX(MAXVAL(MATMUL(ABS(periodic_matrix), &
                                      ABS(right_null(:,m)))),TINY(1.d0))
    ENDDO
  ENDIF
OPEN(751,file='amat_after.dat')
OPEN(752,file='bvec_after.dat')
DO i=1,ndim
WRITE(751,*) amat(i,:)
WRITE(752,*) bvec_lapack(i,:)
ENDDO
CLOSE(751)
CLOSE(752)

  CALL gbsv(ndim,ndim,amat,ipivot,bvec_lapack,info)
  IF(info.NE.0) THEN
    ierr=2
    RETURN
  ENDIF

  IF(trace_join_end) THEN
    DO m=0,nconstraints-1
      dropped_p(m+1,:)=MATMUL(periodic_matrix((m+1)*nl,:),bvec_lapack) &
                       -periodic_rhs((m+1)*nl,:)
      dropped_m(m+1,:)=MATMUL(periodic_matrix(nts_l+(m+1)*nr,:), &
                              bvec_lapack) &
                       -periodic_rhs(nts_l+(m+1)*nr,:)
      DO i=1,ntranseq
        dropped_p_scale(m+1,i)=MAX( &
          DOT_PRODUCT(ABS(periodic_matrix((m+1)*nl,:)), &
                      ABS(bvec_lapack(:,i))) &
          +ABS(periodic_rhs((m+1)*nl,i)),TINY(1.d0))
        dropped_m_scale(m+1,i)=MAX( &
          DOT_PRODUCT(ABS(periodic_matrix(nts_l+(m+1)*nr,:)), &
                      ABS(bvec_lapack(:,i))) &
          +ABS(periodic_rhs(nts_l+(m+1)*nr,i)),TINY(1.d0))
      ENDDO
    ENDDO
    CALL record_join_end_compatibility(source_factor,source_before, &
         source_after,source_scale,measure_sum,compatibility,dropped_p, &
         dropped_m,compatibility_scale,dropped_p_scale,dropped_m_scale, &
         left_null,right_null,left_residual,right_residual,left_scale, &
         right_scale,transfer_error,ierr)
    IF(ierr.NE.0) RETURN
  ENDIF

OPEN(752,file='bvec_solution.dat')
DO i=1,ndim
WRITE(752,*) bvec_lapack(i,:)
ENDDO
CLOSE(752)

  ! Now sources from outgoing become incoming, dimensions switch places !!!

  DEALLOCATE(o%p%source_p,o%p%source_m)
  ALLOCATE(o%p%source_p(nts_l,3),o%p%source_m(nts_r,3))

  o%p%source_p=bvec_lapack(1:nts_l,1:3)
  o%p%source_m=bvec_lapack(nts_l+1:ndim,1:3)

  ! integrals of particle flux and current over fieldline 
  transmat=o%p%qflux+MATMUL(o%p%flux_p,o%p%source_p)           &
                    +MATMUL(o%p%flux_m,o%p%source_m)
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

  o%p%qflux=transmat

  DEALLOCATE(amat,bvec_lapack,ipivot)
  NULLIFY(c_forward)
  NULLIFY(c_backward)

END SUBROUTINE join_ends
