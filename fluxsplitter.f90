MODULE fluxsplit_mod

  USE lapack_band
!
  IMPLICIT NONE

  ! private variables
  INTEGER, PARAMETER, PRIVATE :: dp = KIND(1.0D0)

  PUBLIC fluxsplitter
  INTERFACE fluxsplitter
     MODULE PROCEDURE fluxsplitter
  END INTERFACE

CONTAINS

!
  SUBROUTINE fluxsplitter(mode,ishift,wsplit,eta_in,splitmatrix,ierr)
    !
    ! Splits the flux in the band in two parts. Uses the flux values
    ! in neigbouring bands for determination of the split ratio.
    ! The original bands are numbered from -2 to 2, their boundaries
    ! (levels), eta, are numbered from -3 to 2 so that band #0 is in between
    ! levels -1 and 0. The band to be split is indicated by its relative
    ! position with respect to band #0, ishift (non-zero ishift is needed 
    ! for splitting few first and few last bands). The band is split into
    ! two with relative widths wsplit and 1-wsplit for the left and right
    ! band, respectively. Standard case - wsplit = 0.5d0. 
    ! Allows 3 interpolation modes with different accuracy.
    ! The result is put into splitmatrix so that in order to split flux
    ! f(n+ishift) into fm and fp one has to multiply it like follows.
    ! For mode=1 and mode=0:
    !    fm = sum(splitmatrix(1,-1:1)*f(n-1:n+1))
    !    fp = sum(splitmatrix(2,-1:1)*f(n-1:n+1))
    ! For mode=2 and mode=0:
    !    fm = sum(splitmatrix(1,-2:2)*f(n-2:n+2))
    !    fp = sum(splitmatrix(2,-2:2)*f(n-2:n+2))
    ! Here mode=0 is added for checks
    ! In order to pass to the routine necessaly eta-levels from the array
    ! eta, put at the place of argument eta_in 
    ! in mode=1 : eta(n-2,n+1);  in mode=2 :  eta(n-3,n+2)
    !
    ! Input arguments:
    !          Formal: mode    - split mode: 0 - 0-order (no neigbours)
    !                                        1 - 2-order (two neigbours)
    !                                        2 - 4-order (four neigbours)
    !                  ishift  - index of the band to be split 
    !                  wsplit  - weight of the first split interval
    !                  eta_in  - array of levels (band bondaries)
    ! Output arguments:
    !          Formal: splitmatrix - splitting matrix
    !                  ierr        - error code (0 - normal work)
    !
    !
    INTEGER, INTENT(in)                     :: mode,ishift
    INTEGER, INTENT(out)                    :: ierr
    REAL(kind=dp), INTENT(in)               :: wsplit
    REAL(kind=dp), DIMENSION(*), INTENT(in) :: eta_in
    ! WINNY
    ! make it allocatable
    ! DOUBLE PRECISION, DIMENSION(2,-2:2) :: splitmatrix
    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: splitmatrix
    !
    INTEGER                                       :: i,info,ndim,k
    INTEGER, DIMENSION(:), ALLOCATABLE            :: ipivot
    REAL(kind=dp)                                 :: deleta,deletapow
    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE    :: amat,bvec,etapow
!
    ierr=0
    
    ! WINNY: allocate
    IF (ALLOCATED(splitmatrix)) DEALLOCATE(splitmatrix)
    ALLOCATE(splitmatrix(2,-mode:mode))
    
    !
    IF(mode.EQ.0) THEN
       IF(ABS(ishift).GT.2) THEN
          PRINT *,'fluxsplitter : shift is too large, ishift = ',ishift
          ierr=2
          RETURN
       ENDIF
       splitmatrix=0.d0
       splitmatrix(1,ishift)=wsplit
       splitmatrix(2,ishift)=1.d0-wsplit
    ELSEIF(mode.EQ.1) THEN
       IF(ABS(ishift).GT.1) THEN
          PRINT *,'fluxsplitter : shift is too large, ishift = ',ishift
          ierr=2
          RETURN
       ENDIF
       splitmatrix=0.d0
       ndim=3
       ALLOCATE(amat(ndim,ndim),bvec(ndim,ndim),ipivot(ndim),etapow(-2:1,ndim))
       etapow(-2:1,1)=eta_in(1:4)-eta_in(3)
       DO i=2,ndim
          etapow(:,i)=etapow(:,i-1)*etapow(:,1)
       ENDDO
       bvec=0.d0
       DO i=1,ndim
          bvec(i,i)=1.d0
       ENDDO
       amat=etapow(-1:1,:)-etapow(-2:0,:)
       CALL gbsv(ndim,ndim,amat,ipivot,bvec,info)
       IF(info.NE.0) THEN
          PRINT *,'fluxsplitter : gbsv error = ',info
          ierr=3
          DEALLOCATE(amat,bvec,etapow)
          RETURN
       ENDIF
       etapow(ishift,1)=etapow(ishift-1,1)*(1.d0-wsplit)+etapow(ishift,1)*wsplit
       splitmatrix(1,-1:1)=(etapow(ishift,1)-etapow(ishift-1,1))*bvec(1,:)
       DO i=2,ndim
          etapow(ishift,i)=etapow(ishift,i-1)*etapow(ishift,1)
          splitmatrix(1,-1:1)=splitmatrix(1,-1:1)                                  &
               +(etapow(ishift,i)-etapow(ishift-1,i))*bvec(i,:)
       ENDDO
       splitmatrix(2,-1:1)=-splitmatrix(1,-1:1)
       splitmatrix(2,ishift)=splitmatrix(2,ishift)+1.d0
       DEALLOCATE(amat,bvec,etapow)
    ELSEIF(mode.EQ.2) THEN
       IF(ABS(ishift).GT.2) THEN
          PRINT *,'fluxsplitter : shift is too large, ishift = ',ishift
          ierr=2
          RETURN
       ENDIF
       splitmatrix=0.d0
       ndim=5
       ALLOCATE(amat(ndim,ndim),bvec(ndim,ndim),ipivot(ndim),etapow(-3:2,ndim))
       etapow(-3:2,1)=eta_in(1:6)-eta_in(4)
       DO i=2,ndim
          etapow(:,i)=etapow(:,i-1)*etapow(:,1)
       ENDDO
       bvec=0.d0
       DO i=1,ndim
          bvec(i,i)=1.d0
       ENDDO
       amat=etapow(-2:2,:)-etapow(-3:1,:)
       CALL gbsv(ndim,ndim,amat,ipivot,bvec,info)
       IF(info.NE.0) THEN
          PRINT *,'fluxsplitter : gbsv error = ',info
          ierr=3
          DEALLOCATE(amat,bvec,etapow)
          RETURN
       ENDIF
       etapow(ishift,1)=etapow(ishift-1,1)*(1.d0-wsplit)+etapow(ishift,1)*wsplit
       splitmatrix(1,-2:2)=(etapow(ishift,1)-etapow(ishift-1,1))*bvec(1,:)
       DO i=2,ndim
          etapow(ishift,i)=etapow(ishift,i-1)*etapow(ishift,1)
          splitmatrix(1,-2:2)=splitmatrix(1,-2:2)                                  &
               +(etapow(ishift,i)-etapow(ishift-1,i))*bvec(i,:)
       ENDDO
       splitmatrix(2,-2:2)=-splitmatrix(1,-2:2)
       splitmatrix(2,ishift)=splitmatrix(2,ishift)+1.d0
       DEALLOCATE(amat,bvec,etapow)
    ELSE
       PRINT *,'fluxsplitter : unknown mode, mode = ',mode
       ierr=1
    ENDIF
    !
    RETURN
  END SUBROUTINE fluxsplitter
  
END MODULE fluxsplit_mod
