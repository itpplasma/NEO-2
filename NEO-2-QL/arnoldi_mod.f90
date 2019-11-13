MODULE arnoldi_mod
  INTEGER :: ngrow,ierr
  INTEGER :: ntol,mode=0
  DOUBLE PRECISION :: tol
  complex(kind=kind(1d0)), DIMENSION(:),   ALLOCATABLE :: fzero,ritznum
  complex(kind=kind(1d0)), DIMENSION(:,:), ALLOCATABLE :: eigvecs

contains
  !---------------------------------------------------------------------------------
  !> \brief Solves the equation f=Af+q where A is a matrix and q is a given vector.
  !>
  !> Iterates the system which may be unstable at direct iterations
  !> using subtraction of unstable eigenvectors. Iterations are terminated
  !> when relative error, defined as $\sum(|f_n-f_{n-1}|)/\sum(|f_n|)$
  !> is below the input value or maximum number of combined iterations
  !> is reached.
  !>
  !> Method is a brand of preconditioned iteration:
  !>
  !>         f_new = f_old + P (A f_old + q - f_old)
  !>
  !> where P is a preconditioner which transforms the unstable Krylov subspace
  !> to make it stable. Namely, preconditioner acts on unstbale eigenvecotrs
  !> phi_m as follows
  !>
  !>         P phi_m = phi_m / (1 - lambda_m)
  !>
  !> where lambda_m is the corresponding eigenvalue, and it acts on stable
  !> eigenvectors as follows
  !>
  !>         P phi_m = phi_m
  !>
  !> i.e. leaves them untouched.
  !> Here this method is applied to preconditioned Richardson iteration given
  !> by the routine "next_iteration". The resulting combination corresponds,
  !> again, to Richardson iteration with a combined preconditioner:
  !>
  !>         f_new = f_old + P Lbar (q - L f_old)
  !>
  !> where the rest notation is described in comments to "next_iteration" routine.
  !>
  !> Input  parameters:
  !>            Formal: mode_in        - iteration mode (0 - direct iterations,
  !>                                                         "next_iteration" provides Af
  !>                                                         and q is provided as an input via
  !>                                                         "result"
  !>                                                     1 - "next_iteration" provides Af+q,
  !>                                                         preconditioner stays allocated.
  !>                                                         If the routine is re-entered
  !>                                                         with this mode, old preconditioner
  !>                                                         is used, otherwise it is updated
  !>                                                     3 - just dealocates preconditioner
  !>                    n              - system size
  !>                    narn           - maximum number of Arnoldi iterations
  !>                    relerr         - relative error
  !>                    itermax        - maximum number of combined iterations
  !>                    result_        - source vector q for kinetic equation L f = q, see
  !>                                     comments to "next_iteration" routine.
  !>                    next_iteration - routine computing next iteration, "fnew",
  !>                                     of the solution from the previous, "fold",
  !>                                          fnew = A fold + q
  !>                                     in case of mode_in=2
  !>                                          fnew = A fold
  !>                                     call next_iteration(n,fold,fnew)
  !>                                     where "n" is a vector size
  !>
  !> Output parameters:
  !>            Formal: result         - solution vector
  SUBROUTINE iterator(mode_in,n,narn,relerr,itermax,RESULT_, ispec, &
    & next_iteration)

    use mpiprovider_module, only : mpro
    use collisionality_mod, only : num_spec

    IMPLICIT NONE

    ! tol0 - largest eigenvalue tolerated in combined iterations:
    INTEGER,          PARAMETER :: ntol0=10
    DOUBLE PRECISION, PARAMETER :: tol0=0.5d0

    interface
      subroutine next_iteration(n,fold,fnew)
        integer :: n
        complex(kind=kind(1d0)), dimension(n) :: fold,fnew
      end subroutine next_iteration
    end interface

    integer :: ispec
    INTEGER :: mode_in,n,narn,itermax,i,j,iter,nsize,info,iarnflag
    DOUBLE PRECISION :: relerr
    COMPLEX(kind=kind(1d0)), DIMENSION(n), intent(inout) :: RESULT_
    INTEGER,        DIMENSION(:),   ALLOCATABLE :: ipiv
    complex(kind=kind(1d0)), DIMENSION(:),   ALLOCATABLE :: fold,fnew
    complex(kind=kind(1d0)), DIMENSION(:),   ALLOCATABLE :: coefren
    complex(kind=kind(1d0)), DIMENSION(:,:), ALLOCATABLE :: amat,bvec
    DOUBLE PRECISION, DIMENSION(0:num_spec-1) :: break_cond1
    DOUBLE PRECISION, DIMENSION(0:num_spec-1) :: break_cond2
    complex(kind=kind(1d0)), DIMENSION(:), ALLOCATABLE :: coefren_spec
    complex(kind=kind(1d0)), DIMENSION(:), ALLOCATABLE :: amat_spec

    IF(mode_in.EQ.3) THEN
      mode=mode_in
      IF(ALLOCATED(ritznum)) DEALLOCATE(eigvecs,ritznum)
      RETURN
    ELSEIF(mode_in.EQ.1) THEN
      IF(mode.EQ.1) THEN
        iarnflag=0
      ELSE
        iarnflag=1
        IF(ALLOCATED(ritznum)) DEALLOCATE(eigvecs,ritznum)
        ALLOCATE(ritznum(narn))
      ENDIF
    ELSEIF(mode_in.EQ.0) THEN
      iarnflag=0
      ngrow=0
    ELSE
      PRINT *,'unknown mode'
      RETURN
    ENDIF

    ALLOCATE(fzero(n))
    fzero = RESULT_

    IF(iarnflag.EQ.1) THEN

      ! estimate NARN largest eigenvalues:
      tol=tol0
      ntol=ntol0

      CALL arnoldi(n,narn, ispec, next_iteration)

      IF(ngrow .GT. 0) PRINT *,'ritznum = ',ritznum(1:ngrow)

      IF(ierr.NE.0) THEN
        PRINT *,'iterator: error in arnoldi'
        DEALLOCATE(fzero,eigvecs,ritznum)
        RETURN
      ENDIF

    ENDIF

    if (ispec .eq. 0) print *,'iterator: number of bad modes = ', ngrow
    nsize=ngrow

    ALLOCATE(fold(n),fnew(n))

    mode = mode_in

    ! there are no bad eigenvalues, use direct iterations:
    IF(ngrow.EQ.0) THEN

      if (mode .eq. 2) then
        fold = fzero
      else
        fold = (0.0d0, 0.0d0)
      end if

      DO iter=1,itermax
        CALL next_iteration(n,fold,fnew)
        if (mode .EQ. 2) fnew = fnew + fzero
        break_cond1(ispec)=SUM(ABS(fnew-fold))
        break_cond2(ispec)=relerr*SUM(ABS(fnew))
        PRINT *,iter,ispec,break_cond1(ispec),break_cond2(ispec)
        CALL mpro%allgather_inplace(break_cond1)
        CALL mpro%allgather_inplace(break_cond2)
        IF(ALL(break_cond1 .LE. break_cond2)) EXIT
        !! End Modification by Andreas F. Martitsch (20.08.2015)
        fold=fnew
        IF(iter.EQ.itermax) PRINT *, &
                'iterator: maximum number of iterations reached'
      ENDDO
  !
      PRINT *,'iterator: number of direct iterations = ',iter-1

      RESULT_=fnew
      DEALLOCATE(fold,fnew,fzero)
      RETURN

    ENDIF

    ! compute subtraction matrix:
    ALLOCATE(amat(nsize,nsize),bvec(nsize,nsize),ipiv(nsize),coefren(nsize))
    bvec=(0.d0,0.d0)
    IF(ALLOCATED(coefren_spec)) DEALLOCATE(coefren_spec)
    ALLOCATE(coefren_spec(0:num_spec-1))
    IF(ALLOCATED(amat_spec)) DEALLOCATE(amat_spec)
    ALLOCATE(amat_spec(0:num_spec-1))

    DO i=1,nsize
      bvec(i,i)=(1.d0,0.d0)
      DO j=1,nsize
        amat_spec(ispec)=SUM(CONJG(eigvecs(:,i))*eigvecs(:,j))
        CALL mpro%allgather_inplace(amat_spec)
        amat(i,j)=SUM(amat_spec)*(ritznum(j)-(1.d0,0.d0))
      ENDDO
    ENDDO

    CALL zgesv(nsize,nsize,amat,nsize,ipiv,bvec,nsize,info)

    IF(info.NE.0) THEN
      IF(info.GT.0) THEN
        PRINT *,'iterator: singular matrix in zgesv'
      ELSE
        PRINT *,'iterator: argument ',-info,' has illigal value in zgesv'
      ENDIF
      DEALLOCATE(ritznum,coefren,amat,bvec,ipiv)
      DEALLOCATE(fzero,fold,fnew,eigvecs)
      RETURN
    ENDIF

    ! iterate the solution:
    fold=fzero

    DO iter=1,itermax
      CALL next_iteration(n,fold,fnew)
      IF(mode.EQ.2) fnew=fnew+fzero
      DO j=1,nsize
        coefren_spec(ispec)=SUM(bvec(j,:)                           &
                  *MATMUL(TRANSPOSE(CONJG(eigvecs(:,1:nsize))),fnew-fold))
        CALL mpro%allgather_inplace(coefren_spec)
        coefren(j)=ritznum(j)*SUM(coefren_spec)
      ENDDO
      fnew=fnew-MATMUL(eigvecs(:,1:nsize),coefren)
      break_cond1(ispec)=SUM(ABS(fnew-fold))
      break_cond2(ispec)=relerr*SUM(ABS(fnew))
      print *,iter,break_cond1(ispec),break_cond2(ispec),ispec
      CALL mpro%allgather_inplace(break_cond1)
      CALL mpro%allgather_inplace(break_cond2)
      IF(ALL(break_cond1 .LE. break_cond2)) EXIT
      fold = fnew
      IF(iter.EQ.itermax) PRINT *,'iterator: maximum number of iterations reached'
    ENDDO

    if (ispec .eq. 0) print *,'iterator: number of stabilized iterations = ', iter-1

    RESULT_=fnew

    DEALLOCATE(coefren,amat,bvec,ipiv,fold,fnew,fzero)

  END SUBROUTINE iterator

  !-----------------------------------------------------------------------------
  !> Computes m Ritz eigenvalues (approximations to extreme eigenvalues)
  !> of the iteration procedure of the vector with dimension n.
  !> Eigenvalues are computed by means of Arnoldi iterations.
  !> Optionally computes Ritz vectors (approximation to eigenvectors).
  !>
  !> Input  parameters:
  !> Formal:             n              - system dimension
  !>                     mmax           - maximum number of Ritz eigenvalues
  !>                     ispec
  !>                     next_iteration - routine computing next iteration
  !>                                      of the solution from the previous
  !> Module arnoldi_mod: tol            - eigenvectors are not computed for
  !>                                      eigenvalues smaller than this number
  !> Output parameters:
  !> Module arnoldi_mod: ngrow          - number of eigenvalues larger or equal
  !>                                      to TOL
  !>                     ritznum        - Ritz eigenvalues
  !>                     eigvecs        - array of eigenvectors, size - (m,ngrow)
  !>                     ierr           - error code (0 - normal work, 1 - error)
  SUBROUTINE arnoldi(n, mmax, ispec, next_iteration)

    use mpiprovider_module, only : mpro

    use collisionality_mod, only : num_spec

    IMPLICIT NONE

    ! set independent accuracy-level for eigenvalue computation
    DOUBLE PRECISION, PARAMETER :: epserr_ritznum=1.0d-3

    interface
      subroutine next_iteration(n,fold,fnew)
        integer :: n
        complex(kind=kind(1d0)), dimension(n) :: fold,fnew
      end subroutine next_iteration
    end interface

    integer :: ispec
    INTEGER                                       :: n,m,k,j,mmax,mbeg,ncount
    INTEGER :: driv_spec
    complex(kind=kind(1d0)),   DIMENSION(:),   ALLOCATABLE :: fold,fnew,ritznum_prev
    complex(kind=kind(1d0)),   DIMENSION(:,:), ALLOCATABLE :: qvecs,hmat,eigh,qvecs_prev
    complex(kind=kind(1d0)),   DIMENSION(:), ALLOCATABLE :: q_spec, h_spec

    INTEGER :: m_tol, m_ind
    complex(kind=kind(1d0)), DIMENSION(500) :: ritzum_write

    ALLOCATE(fold(n),fnew(n))
    ALLOCATE(qvecs_prev(n,1),ritznum_prev(mmax))
    ALLOCATE(hmat(mmax,mmax))

    fold=(0.d0,0.d0)

    hmat=(0.d0,0.d0)

    mode = 1

    CALL next_iteration(n,fold,fnew)

    mode = 2

    IF(ALLOCATED(q_spec)) DEALLOCATE(q_spec)
    ALLOCATE(q_spec(0:num_spec-1))
    q_spec=0.0d0
    IF(ALLOCATED(h_spec)) DEALLOCATE(h_spec)
    ALLOCATE(h_spec(0:num_spec-1))
    h_spec=0.0d0
    q_spec(ispec)=SUM(CONJG(fnew)*fnew)
    CALL mpro%allgather_inplace(q_spec)
    qvecs_prev(:,1)=fnew/SQRT(SUM(q_spec))
    ierr=0
    mbeg=2
    ncount=0

    DO m=2,mmax

      ALLOCATE(qvecs(n,m))
      qvecs(:,1:m-1)=qvecs_prev(:,1:m-1)

      ALLOCATE(eigh(m,m))

      DO k=mbeg,m
        fold=qvecs(:,k-1)
        call next_iteration(n, fold, fnew)
        qvecs(:,k)=fnew
        DO j=1,k-1
          h_spec=0.0d0
          h_spec(ispec)=SUM(CONJG(qvecs(:,j))*qvecs(:,k))
          CALL mpro%allgather_inplace(h_spec)
          hmat(j,k-1)=SUM(h_spec)
          qvecs(:,k)=qvecs(:,k)-hmat(j,k-1)*qvecs(:,j)
        ENDDO
        h_spec=0.0d0
        h_spec(ispec)=SUM(CONJG(qvecs(:,k))*qvecs(:,k))
        CALL mpro%allgather_inplace(h_spec)
        hmat(k,k-1)=SQRT(SUM(h_spec))
        qvecs(:,k)=qvecs(:,k)/hmat(k,k-1)
      ENDDO

      CALL try_eigvecvals(m,tol,hmat(1:m,1:m),ngrow,ritznum(1:m),eigh,ierr)

      IF(m.GT.2) THEN

        ! check for convergence of ritznum exceeding tol
        m_tol=m-1
        DO m_ind = 1,m-1
          IF(ABS(ritznum(m_ind)) .LT. tol) THEN
            m_tol = m_ind
            EXIT
          ENDIF
        ENDDO
        IF(SUM(ABS(ritznum(1:m_tol)-ritznum_prev(1:m_tol))).LT.(epserr_ritznum*m_tol)) THEN
          ncount=ncount+1
        ELSE
          ncount=0
        ENDIF
      ENDIF
      ritznum_prev(1:m)=ritznum(1:m)

      IF(ncount.GE.ntol.OR.m.EQ.mmax) THEN
        IF(ALLOCATED(eigvecs)) DEALLOCATE(eigvecs)
        ALLOCATE(eigvecs(n,ngrow))

        eigvecs=MATMUL(qvecs(:,1:m),eigh(1:m,1:ngrow))

        PRINT *,'arnoldi: number of iterations = ',m
        EXIT
      ENDIF

      DEALLOCATE(qvecs_prev)
      ALLOCATE(qvecs_prev(n,m))
      qvecs_prev=qvecs
      DEALLOCATE(qvecs)

      DEALLOCATE(eigh)
      mbeg=m+1

    ENDDO

    DEALLOCATE(fold,fnew,qvecs,qvecs_prev,hmat)
    DEALLOCATE(q_spec,h_spec)


  END SUBROUTINE arnoldi


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

    integer, intent(in) :: m
    integer, intent(out) :: ngrowing,ierr
    double precision, intent(in) :: tol
    complex(kind=kind(1d0)), dimension(:,:), intent(in) :: hmat
    complex(kind=kind(1d0)), dimension(m,m), intent(out) :: eigh
    complex(kind=kind(1d0)), dimension(m) :: ritznum

    integer :: k,j,lwork,info

    double precision, parameter :: tiny_diff = 1.d-12

    complex(kind=kind(1d0)) :: tmp

    logical,          dimension(:),   allocatable :: selec
    integer,          dimension(:),   allocatable :: ifailr
    double precision, dimension(:),   allocatable :: rwork
    complex(kind=kind(1d0)), dimension(:),   allocatable :: work,rnum
    complex(kind=kind(1d0)), dimension(:,:), allocatable :: hmat_work

    ierr=0

    allocate(hmat_work(m,m))

    hmat_work=hmat

    allocate(work(1))
    lwork=-1

    call zhseqr('E','N',m,1,m,hmat_work,m,ritznum,hmat_work,m,work,lwork,info)

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

    lwork = int(work(1))
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
        if(aimag(ritznum(j)).gt.0.d0) then
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
      ngrowing=ngrowing+1
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
