MODULE arnoldi_mod
  INTEGER :: ngrow,ierr
  INTEGER :: ntol,mode=0
  DOUBLE PRECISION :: tol
  DOUBLE COMPLEX, DIMENSION(:),   ALLOCATABLE :: fzero,ritznum
  DOUBLE COMPLEX, DIMENSION(:,:), ALLOCATABLE :: eigvecs

contains
  !---------------------------------------------------------------------------------
  !
  !> Computes eigenvalues, ritznum, of the upper Hessenberg matrix hmat
  !> of the dimension (m,m), orders eigenvelues into the decreasing by module
  !> sequence and computes the eigenvectors, eigh, for eigenvalues exceeding
  !> the tolerance tol (number of these eigenvalues is ngrowing)
  !>
  !> Input arguments:
  !>          Formal: m        - matrix size
  !>                  tol      - tolerance
  !>                  hmat     - upper Hessenberg matrix
  !> Output arguments:
  !>          Formal: ngrowing - number of exceeding the tolerance
  !>                  ritznum  - eigenvalues
  !>                  eigh     - eigenvectors
  !>                  ierr     - error code (0 - normal work)
  subroutine try_eigvecvals(m,tol,hmat,ngrowing,ritznum,eigh,ierr)

    implicit none

    integer :: m,ngrowing,ierr,k,j,lwork,info

    double precision, parameter :: tiny_diff = 1.d-12

    double precision :: tol
    double complex   :: tmp

    double complex, dimension(m)   :: ritznum
    !! Modification by Andreas F. Martitsch (19.10.2016)
    ! old:
    !DOUBLE COMPLEX, DIMENSION(m,m) :: hmat,eigh
    ! new:
    DOUBLE COMPLEX, DIMENSION(:,:) :: hmat
    DOUBLE COMPLEX, DIMENSION(m,m) :: eigh
    !! End Modification by Andreas F. Martitsch (19.10.2016)

    LOGICAL,          DIMENSION(:),   ALLOCATABLE :: selec
    INTEGER,          DIMENSION(:),   ALLOCATABLE :: ifailr
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: rwork
    DOUBLE COMPLEX,   DIMENSION(:),   ALLOCATABLE :: work,rnum
    DOUBLE COMPLEX,   DIMENSION(:,:), ALLOCATABLE :: hmat_work

    ierr=0

    ALLOCATE(hmat_work(m,m))

    hmat_work=hmat

    ALLOCATE(work(1))
    lwork=-1

    CALL zhseqr('E','N',m,1,m,hmat_work,m,ritznum,hmat_work,m,work,lwork,info)

    IF(info.NE.0) THEN
      IF(info.GT.0) THEN
         PRINT *,'arnoldi: zhseqr failed to compute all eigenvalues'
         PRINT *,'info: ',info
      ELSE
        PRINT *,'arnoldi: argument ',-info,' has illigal value in zhseqr'
      ENDIF
      DEALLOCATE(hmat_work,work)
      ierr=1
      RETURN
    ENDIF

    lwork=work(1)
    DEALLOCATE(work)
    ALLOCATE(work(lwork))

    CALL zhseqr('E','N',m,1,m,hmat_work,m,ritznum,hmat_work,m,work,lwork,info)

    IF(info.NE.0) THEN
      IF(info.GT.0) THEN
        PRINT *,'arnoldi: zhseqr failed to compute all eigenvalues'
      ELSE
        PRINT *,'arnoldi: argument ',-info,' has illigal value in zhseqr'
      ENDIF
      DEALLOCATE(hmat_work,work)
      ierr=1
      RETURN
    ENDIF

    ! Sort eigenvalues acording to magnitude.
    DO k=1,m
      info=0
      DO j=2,m
        IF(ABS(ritznum(j)).GT.ABS(ritznum(j-1))) THEN
          info=1
          tmp=ritznum(j-1)
          ritznum(j-1)=ritznum(j)
          ritznum(j)=tmp
        ENDIF
      ENDDO
      IF(info.EQ.0) EXIT
    ENDDO

    ! Sort eigenvalues so that for conjugate pairs, the one with imaginary
    ! part > 0 comes first.
    ! If this is not done, they can be in arbitrary order in subsequent
    ! iterations, which will cause problems with checking of convergence.
    do j=2,m
      if(abs(ritznum(j)-conjg(ritznum(j-1))) .lt. tiny_diff) then
        if(dimag(ritznum(j)).gt.0.d0) then
          tmp = ritznum(j-1)
          ritznum(j-1) = ritznum(j)
          ritznum(j) = tmp
        end if
      end if
    end do

    DEALLOCATE(work)

    ! compute how many eigenvalues exceed the tolerance (TOL):
    ALLOCATE(selec(m),rnum(m))
    selec=.FALSE.
    ngrowing=0
    DO j=1,m
      IF(ABS(ritznum(j)).LT.tol) EXIT
      ngrow=ngrow+1
      selec(j)=.TRUE.
    ENDDO
    rnum=ritznum
    hmat_work=hmat
    ALLOCATE(work(m*m),rwork(m),ifailr(m))
    eigh=(0.d0,0.d0)

    CALL zhsein('R','Q','N',selec,m,hmat_work,m,rnum,rnum,1,eigh(:,1:ngrowing),m,  &
              ngrowing,ngrowing,work,rwork,ifailr,ifailr,info)

    IF(info.NE.0) THEN
      IF(info.GT.0) THEN
        PRINT *,'arnoldi: ',info,' eigenvectors not converged in zhsein'
      ELSE
        PRINT *,'arnoldi: argument ',-info,' has illigal value in zhsein'
      ENDIF
      ierr=1
    ENDIF

    DEALLOCATE(hmat_work,work,rwork,selec,rnum,ifailr)

  END SUBROUTINE try_eigvecvals

END MODULE arnoldi_mod
