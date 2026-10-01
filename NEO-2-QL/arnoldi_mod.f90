MODULE arnoldi_mod
  INTEGER :: ngrow,ierr
  INTEGER :: ntol,mode=0
  DOUBLE PRECISION :: tol
  complex(kind=kind(1d0)), DIMENSION(:),   ALLOCATABLE :: fzero,ritznum
  complex(kind=kind(1d0)), DIMENSION(:),   ALLOCATABLE :: f_init_arnoldi
  complex(kind=kind(1d0)), DIMENSION(:,:), ALLOCATABLE :: eigvecs

  logical :: lsw_write_flux_surface_distribution = .true.

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
  !> to make it stable. Namely, preconditioner acts on unstable eigenvecotrs
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
        PRINT *,'ERROR in iterator: error in arnoldi'
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

      fold = (0.0d0, 0.0d0)

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
        IF(iter.EQ.itermax) CALL warn_not_converged('iterator (direct)', &
                relerr*MAXVAL(break_cond1/MAX(break_cond2,TINY(1d0))), relerr)
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
      IF(iter.EQ.itermax) CALL warn_not_converged('iterator (deflated)', &
              relerr*MAXVAL(break_cond1/MAX(break_cond2,TINY(1d0))), relerr)
    ENDDO

    if (ispec .eq. 0) print *,'iterator: number of stabilized iterations = ', iter-1

    RESULT_=fnew

    DEALLOCATE(coefren,amat,bvec,ipiv,fold,fnew,fzero)

  END SUBROUTINE iterator

  !---------------------------------------------------------------------------------
  !> \brief Solves the fixed-point problem f = A f + c by restarted GMRES.
  !>
  !> "next_iteration" is the preconditioned Richardson step of the ripple
  !> solver, fnew = A fold + c with c = Lbar q (module mode /= 2) and
  !> fnew = A fold (mode = 2). The fixed point solves (I - A) f = c with
  !> I - A = Lbar L, i.e. the kinetic equation L f = q left-preconditioned
  !> by the sparse LU of L_V - L_D. GMRES needs no estimate of the
  !> unstable eigenvalues of A, so it converges where the deflated
  !> Richardson iteration of "iterator" stalls or diverges (eigenvalues of A
  !> on or outside the unit circle, inaccurate Ritz values close to 1).
  !> One GMRES iteration costs one call of next_iteration, the same as one
  !> Richardson step, and no Arnoldi eigenvalue estimate is needed.
  !>
  !> In the multi-species case each MPI rank holds the part of the vector
  !> for its species ispec; all inner products are reduced over species, so
  !> all ranks take identical decisions and call next_iteration equally often.
  !>
  !> Residual: ||N(f) - f||_2 / ||f||_2 (global over species), the 2-norm
  !> analogue of the criterion of "iterator", checked with the true residual
  !> after each restart cycle (the recursive GMRES estimate only ends a
  !> cycle). I - A can be nearly singular (Ritz values of A within 1e-5 of
  !> 1, |f|/|c| ~ 1e3), so a residual of relerr does not yet give an error
  !> of relerr in f. The iteration therefore continues to reltarget < relerr
  !> or until a cycle no longer halves the residual (noise floor of the
  !> sparse solves), and warns only if the residual stays above relerr.
  !> The returned vector is one more step N(f).
  !>
  !> Input:  n         - local system size
  !>         mrestart  - Krylov subspace size between restarts
  !>         relerr    - relative residual that must be reached
  !>         reltarget - relative residual aimed at (<= relerr)
  !>         itermax   - maximum number of GMRES iterations
  !>         result_   - source vector q (see next_iteration)
  !> Output: result_   - solution f
  subroutine gmres_iterator(n, mrestart, relerr, reltarget, itermax, result_, ispec, &
      next_iteration)

    integer, intent(in) :: n, mrestart, itermax, ispec
    double precision, intent(in) :: relerr, reltarget
    complex(kind=kind(1d0)), dimension(n), intent(inout) :: result_

    interface
      subroutine next_iteration(n,fold,fnew)
        integer :: n
        complex(kind=kind(1d0)), dimension(n) :: fold,fnew
      end subroutine next_iteration
    end interface

    complex(kind=kind(1d0)), dimension(:,:), allocatable :: v, h
    complex(kind=kind(1d0)), dimension(:), allocatable :: x, w, g, sn, y, hcol
    double precision, dimension(:), allocatable :: cs
    double precision :: cnorm, beta, relres, xnorm, relres_prev
    complex(kind=kind(1d0)) :: temp
    integer :: m, j, k, iter, ncycle
    logical :: converged

    m = max(1, min(mrestart, itermax))
    allocate(fzero(n), x(n), w(n), v(n, m+1), h(m+1, m), g(m+1), sn(m), &
      cs(m), y(m), hcol(m+1))
    fzero = result_

    ! c = N(0): right-hand side and initial residual for x = 0.
    x = (0.d0, 0.d0)
    mode = 1
    call next_iteration(n, x, w)
    cnorm = global_norm(w, ispec)
    if (cnorm .eq. 0.d0) then
      result_ = (0.d0, 0.d0)
      call finish()
      return
    end if

    iter = 0
    ncycle = 0
    xnorm = 0.d0
    relres_prev = huge(1.d0)
    converged = .false.
    beta = cnorm
    do
      ncycle = ncycle + 1
      v(:, 1) = w/beta
      g = (0.d0, 0.d0)
      g(1) = beta
      mode = 2
      do j = 1, m
        iter = iter + 1
        call next_iteration(n, v(:, j), w)
        w = v(:, j) - w
        ! classical Gram-Schmidt, twice (one reduction per pass)
        call global_dot_block(v(:, 1:j), w, hcol(1:j), ispec)
        w = w - matmul(v(:, 1:j), hcol(1:j))
        h(1:j, j) = hcol(1:j)
        call global_dot_block(v(:, 1:j), w, hcol(1:j), ispec)
        w = w - matmul(v(:, 1:j), hcol(1:j))
        h(1:j, j) = h(1:j, j) + hcol(1:j)
        h(j+1, j) = global_norm(w, ispec)
        if (abs(h(j+1, j)) .gt. 0.d0) v(:, j+1) = w/h(j+1, j)
        ! apply previous Givens rotations, then eliminate h(j+1,j)
        do k = 1, j - 1
          temp = cs(k)*h(k, j) + sn(k)*h(k+1, j)
          h(k+1, j) = -conjg(sn(k))*h(k, j) + cs(k)*h(k+1, j)
          h(k, j) = temp
        end do
        call zlartg(h(j, j), h(j+1, j), cs(j), sn(j), temp)
        h(j, j) = temp
        h(j+1, j) = (0.d0, 0.d0)
        g(j+1) = -conjg(sn(j))*g(j)
        g(j) = cs(j)*g(j)
        ! y = R^{-1} g (cheap, j <= m) to estimate ||x + V y||
        do k = j, 1, -1
          y(k) = (g(k) - sum(h(k, k+1:j)*y(k+1:j)))/h(k, k)
        end do
        relres = abs(g(j+1))/max(xnorm, sqrt(sum(abs(y(1:j))**2)))
        if (ispec .eq. 0) print *, 'gmres_iterator: ', iter, relres
        if (relres .le. reltarget .or. iter .ge. itermax) exit
      end do
      j = min(j, m)
      x = x + matmul(v(:, 1:j), y(1:j))
      xnorm = global_norm(x, ispec)
      ! True residual r = N(x) - x decides convergence (the recursive GMRES
      ! estimate can be too optimistic for this ill-conditioned operator)
      ! and is the start vector of the next cycle.
      mode = 1
      call next_iteration(n, x, w)
      w = w - x
      beta = global_norm(w, ispec)
      relres = beta/xnorm
      converged = relres .le. reltarget
      if (converged .or. iter .ge. itermax) exit
      ! A cycle that does not halve the true residual has reached the noise
      ! floor of next_iteration (rounding in the sparse solves); stop.
      if (relres .gt. 0.5d0*relres_prev) then
        if (ispec .eq. 0) print *, 'gmres_iterator: stagnation at the noise floor'
        exit
      end if
      relres_prev = relres
    end do

    result_ = x + w
    if (ispec .eq. 0) print *, 'gmres_iterator: iterations = ', iter, &
      ' restarts = ', ncycle - 1, ' true relative residual = ', relres, &
      ' |f|/|c| = ', xnorm/cnorm
    if (relres .gt. relerr) call warn_not_converged('gmres_iterator', &
      relres, relerr)
    call finish()

  contains

    subroutine finish()
      deallocate(fzero, x, w, v, h, g, sn, cs, y, hcol)
      ! Do not let a later call of iterator reuse a stale preconditioner.
      mode = 0
    end subroutine finish

  end subroutine gmres_iterator

  !> Global 2-norm of a vector distributed over species (one rank per species).
  function global_norm(vec, ispec) result(nrm)
    use mpiprovider_module, only : mpro
    use collisionality_mod, only : num_spec
    complex(kind=kind(1d0)), dimension(:), intent(in) :: vec
    integer, intent(in) :: ispec
    double precision :: nrm
    double precision, dimension(0:num_spec-1) :: buf

    buf = 0.d0
    buf(ispec) = sum(dble(vec)**2 + aimag(vec)**2)
    call mpro%allgather_inplace(buf)
    nrm = sqrt(sum(buf))
  end function global_norm

  !> Inner products dots(i) = <basis(:,i), vec>, reduced over species.
  subroutine global_dot_block(basis, vec, dots, ispec)
    use mpiprovider_module, only : mpro
    use collisionality_mod, only : num_spec
    complex(kind=kind(1d0)), dimension(:,:), intent(in) :: basis
    complex(kind=kind(1d0)), dimension(:), intent(in) :: vec
    complex(kind=kind(1d0)), dimension(:), intent(out) :: dots
    integer, intent(in) :: ispec
    complex(kind=kind(1d0)), dimension(size(dots), 0:num_spec-1) :: buf

    buf = (0.d0, 0.d0)
    buf(:, ispec) = conjg(matmul(conjg(vec), basis))
    call mpro%allgather_inplace(buf)
    dots = sum(buf, dim=2)
  end subroutine global_dot_block

  !> Loud warning for an iteration that stopped before reaching relerr.
  subroutine warn_not_converged(name, achieved, relerr)
    use, intrinsic :: iso_fortran_env, only : error_unit
    character(len=*), intent(in) :: name
    double precision, intent(in) :: achieved, relerr

    write (*, '(a,a,a,es10.3,a,es10.3)') ' WARNING: ', name, &
      ': NOT CONVERGED, relative error ', achieved, ' > tolerance ', relerr
    write (error_unit, '(a,a,a,es10.3,a,es10.3)') 'WARNING: ', name, &
      ': NOT CONVERGED, relative error ', achieved, ' > tolerance ', relerr
  end subroutine warn_not_converged

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

    integer, intent(in) :: ispec, n, mmax
    INTEGER                                       :: m,k,j,mbeg,ncount
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

    fnew = f_init_arnoldi

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

    ! Standard function(?) compute the eigenvalues of a Hessenberg matrix
    call zhseqr('E','N',m,1,m,hmat_work,m,ritznum,hmat_work,m,work,lwork,info)

    IF(info.NE.0) THEN
      IF(info.GT.0) THEN
         write(*,*)'try_eigvecvals: zhseqr failed to determine optimal value for lwork'
      ELSE
        PRINT *,'try_eigvecvals: argument ',-info,' has illegal value in zhseqr'
      ENDIF
      DEALLOCATE(hmat_work,work)
      ierr=1
      RETURN
    ENDIF

    lwork = int(work(1))
    DEALLOCATE(work)
    ALLOCATE(work(lwork))

    ! Standard function(?) compute the eigenvalues of a Hessenberg matrix
    CALL zhseqr('E','N',m,1,m,hmat_work,m,ritznum,hmat_work,m,work,lwork,info)

    IF(info.NE.0) THEN
      IF(info.GT.0) THEN
        PRINT *,'try_eigvecvals: zhseqr failed to compute all eigenvalues'
      ELSE
        PRINT *,'try_eigvecvals: argument ',-info,' has illegal value in zhseqr'
      ENDIF
      DEALLOCATE(hmat_work,work)
      ierr=1
      RETURN
    ENDIF

    ! Sort eigenvalues acording to magnitude, and so that for conjugate
    ! pairs, the one with imaginary part > 0 comes first.
    ! If this is not done, they can be in arbitrary order in subsequent
    ! iterations, which will cause problems with checking of convergence.
    call sort_eigenvalues(ritznum)

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

  !---------------------------------------------------------------------
  !> Write the distribution for the flux surface to a file. So far the
  !> settings are all hardcoded (e.g. filename, filenumber, if to write).
  !> If this subroutine is called should depend on the lsw_* variable.
  !>
  !> Side effects:
  !> Creates one file for each species with output data.
  !>
  !> input:
  !> ------
  !> geometrical_factor: factors from integration for each point.
  !> distribution_function:
  !> ind_start: array with starting indices
  !> number_passing_particles:
  !> delphi: array with Delta phi values for each point.
  !> phi_grid: array with phi values for each point.
  !> number_elements: number of phi points (and size of delphi, phi_grid).
  subroutine write_flux_surface_distribution(geometrical_factor, distribution_function, &
      & ind_start, number_passing_particles, delphi, phi_grid, number_elements)

    use mpiprovider_module, only : mpro

    use neo_magfie, only: boozer_iota
    use collisionality_mod, only : num_spec
    use rkstep_mod, only : lag

    implicit none

    complex(kind=kind(1d0)), dimension(:,:), intent(in) :: geometrical_factor
    complex(kind=kind(1d0)), dimension(:,:,:), intent(in) :: distribution_function
    real(kind=kind(1d0)), dimension(:), intent(in) :: delphi, phi_grid

    integer, dimension(:), intent(in) :: ind_start, number_passing_particles
    integer, intent(in) :: number_elements

    integer :: range_start, range_end
    integer :: npassing, istep, ispec

    integer, save :: number_of_call = 1

    character(len=2) :: species_number_as_string
    character(len=2) :: call_number_as_string

    if (any(lbound(phi_grid) /= 1)) then
      stop 'ERROR write_flux_surface_distribution: lbound has not expected value.'
    end if
    if (any(ubound(phi_grid) /= number_elements)) then
      stop 'ERROR write_flux_surface_distribution: ubound has not expected value.'
    end if

    ispec = mpro%getRank()+1

    write(species_number_as_string,'(I2.2)') ispec
    write(call_number_as_string,'(I2.2)') number_of_call
    open(543+ispec,file="flux_surface_distribution_spec"//species_number_as_string//"_"//call_number_as_string//".dat")

    do istep = 1,number_elements
      npassing = number_passing_particles(istep)

      range_start = ind_start(istep) + 1
      range_end = range_start + 2*(lag + 1)*(npassing+1) - 1

      if (istep < number_elements) then
        if (range_end /= ind_start(istep+1)-1+1) then
          write(*,*) 'ERROR: sanity check of range_end failed in subroutine'
          write(*,*) '       write_flux_surface_distribution.'
          write(*,*) 'istep: ', istep
          write(*,*) 'range_start: ', range_start
          write(*,*) 'range_end: ', range_end
          write(*,*) 'ind_start(istep+1)...: ', ind_start(istep+1)-1+1
          stop
        end if
      end if
      write(543+ispec,*) boozer_iota*phi_grid(istep), &
        & real(matmul(conjg(geometrical_factor(1:3, range_start:range_end)), &
        &     real(distribution_function(range_start:range_end, 1:3, ispec))))
    end do

    close(543+ispec)

    number_of_call = number_of_call + 1

  end subroutine write_flux_surface_distribution

  !> \brief Simple sorting routine for small 1d arrays (decreasing magnitude).
  subroutine sort_complex_array(array)

    implicit none

    complex(kind=kind(1d0)), dimension(:), intent(inout) :: array

    logical :: entries_swapped
    integer :: m, k, j
    complex(kind=kind(1d0)) :: temp

    m = ubound(array, 1)

    do k=1,m
      entries_swapped = .false.
      do j=2,m
        if (abs(array(j)) .GT. abs(array(j-1))) then
          entries_swapped = .true.
          temp = array(j-1)
          array(j-1) = array(j)
          array(j) = temp
        end if
      end do
      if (.not. entries_swapped) exit
    end do

  end subroutine sort_complex_array

  !> \brief Sort 1d array of complex eigenvalues.
  !>
  !> This subroutine sorts a 1d array of complex eigenvalues. First
  !> entries are sorted according to magnitude (largest first), then is
  !> made sure, that conjugate pairs of eigenvalues always have the same
  !> order, with one with positive imaginary coming first.
  subroutine sort_eigenvalues(eiv)

    implicit none

    complex(kind=kind(1d0)), dimension(:), intent(inout) :: eiv

    ! Absolute accuracy for determining conjugate pairs.
    double precision, parameter :: tiny_diff = 1.d-12

    integer :: j, m
    complex(kind=kind(1d0)) :: temp

    call sort_complex_array(eiv)

    m = ubound(eiv, 1)

    ! Additional sorting: make sure complex conjugate pairs do not have
    ! random order.
    do j=2,m
      if (abs(eiv(j) - conjg(eiv(j-1))) .lt. tiny_diff) then
        if (aimag(eiv(j)) .gt. 0.d0) then
          temp = eiv(j-1)
          eiv(j-1) = eiv(j)
          eiv(j) = temp
        end if
      end if
    end do

  end subroutine sort_eigenvalues

END MODULE arnoldi_mod
