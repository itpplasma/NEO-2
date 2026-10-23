MODULE arnoldi_mod
  INTEGER :: ngrow,ierr
  INTEGER :: ntol,mode=0
  DOUBLE PRECISION :: tol
  complex(kind=kind(1d0)), DIMENSION(:),   ALLOCATABLE :: fzero,ritznum
  complex(kind=kind(1d0)), DIMENSION(:),   ALLOCATABLE :: f_init_arnoldi
  complex(kind=kind(1d0)), DIMENSION(:,:), ALLOCATABLE :: eigvecs

  logical :: stabilization_probe = .false., stabilization_stage_trace = .false.
  logical :: lsw_write_flux_surface_distribution = .true.

contains
  !---------------------------------------------------------------------------------
  !> \brief Solves the equation f=Af+q where A is a matrix and q is a given vector.
  !>
  !> Iterates the system which may be unstable at direct iterations
  !> using subtraction of unstable eigenvectors. Iterations are terminated
  !> when the supplied original-equation residual is below relerr.
  !> Without a residual callback, uses the unaccelerated fixed-point defect.
  !> Iteration exhaustion returns a nonzero status; it is not convergence.
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
    & next_iteration, status, iterations, equation_residual, evaluate_residual, physical_action, &
    & early_original_accept, evaluate_constraints)

    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
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

    interface
      subroutine evaluate_residual(n, solution, residual)
        integer, intent(in) :: n
        complex(kind=kind(1d0)), intent(in) :: solution(n)
        double precision, intent(out) :: residual
      end subroutine evaluate_residual
    end interface
    optional :: evaluate_residual
    logical, optional, intent(in) :: early_original_accept
    logical :: early_accept
    interface
      subroutine evaluate_constraints(n, candidate, constraint)
        integer, intent(in) :: n
        complex(kind=kind(1d0)), intent(in) :: candidate(n)
        double precision, intent(out) :: constraint
      end subroutine evaluate_constraints
    end interface
    optional :: evaluate_constraints
    interface
      subroutine physical_action(v, lv, entry_sum)
        complex(kind=kind(1d0)), intent(in) :: v(:)
        complex(kind=kind(1d0)), intent(out) :: lv(:)
        double precision, intent(out) :: entry_sum
      end subroutine physical_action
    end interface
    optional :: physical_action
    integer, optional, intent(out) :: status, iterations
    double precision, optional, intent(out) :: equation_residual
    integer :: solve_status, completed_iterations
    double precision :: residual_value
    logical :: converged
    integer :: saved_mode
    double precision, allocatable :: gsv(:), ssv(:)
    double precision :: residual_bound, inverse_defect
    double precision :: global_norms(0:num_spec-1), full_defect_sq
    complex(kind=kind(1d0)), allocatable :: gram_probe(:,:), s_probe(:,:)
    complex(kind=kind(1d0)), allocatable :: before_probe(:), action_probe(:)
    double precision :: entry_bound, defect_norm, vector_norm

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

    if (stabilization_probe) then
      if (.not. present(evaluate_residual)) error stop 'probe needs residual'
      if (.not. present(physical_action)) error stop 'probe needs physical L'
      if (mode_in /= 1) error stop 'probe requires first stabilized AX solve'
    end if
    early_accept=.false.
    if (present(early_original_accept)) early_accept=early_original_accept
    solve_status = 1
    completed_iterations = 0
    residual_value = huge(1.0d0)
    if (present(status)) status = solve_status
    if (present(iterations)) iterations = completed_iterations
    if (present(equation_residual)) equation_residual = residual_value
    if (early_accept .and. mode_in /= 3) then
      if (.not. present(evaluate_residual) .or. .not. present(evaluate_constraints)) then
        if (present(status)) status=4
        return
      end if
    end if
    if (mode_in /= 3) then
      if (itermax < 1 .or. (mode_in == 1 .and. narn < 2) .or. .not. ieee_is_finite(relerr) .or. relerr <= 0.0d0) then
        if (present(status)) status = 4
        return
      end if
    end if

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
      if (present(status)) status = 4
      RETURN
    ENDIF

    ALLOCATE(fzero(n))
    fzero = RESULT_

    IF(iarnflag.EQ.1) THEN

      ! estimate NARN largest eigenvalues:
      tol=tol0
      ntol=ntol0

      CALL arnoldi(n,narn, ispec, next_iteration)

      IF(ierr.NE.0) THEN
        PRINT *,'ERROR in iterator: error in arnoldi'
        if (present(status)) status = 2
        DEALLOCATE(fzero,ritznum)
        IF(ALLOCATED(eigvecs)) DEALLOCATE(eigvecs)
        RETURN
      ENDIF

      IF(ngrow .GT. 0) PRINT *,'ritznum = ',ritznum(1:ngrow)

    ENDIF

    if (ispec .eq. 0) print *,'iterator: number of bad modes = ', ngrow
    nsize=ngrow

    ALLOCATE(fold(n),fnew(n))

    mode = mode_in

    if (stabilization_probe .and. ngrow /= 4) &
      error stop 'diagnostic expected measured four-mode space'
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
        call check_equation(fnew, converged)
        if (converged .or. solve_status /= 1) exit
        !! End Modification by Andreas F. Martitsch (20.08.2015)
        fold=fnew
        IF(iter.EQ.itermax) PRINT *, &
                'iterator: maximum number of iterations reached'
      ENDDO
  !
      PRINT *,'iterator: number of direct iterations = ',completed_iterations
      call publish_status()

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

    if (stabilization_probe) then
      allocate(gsv(nsize),ssv(nsize))
      allocate(gram_probe(nsize,nsize), s_probe(nsize,nsize))
      allocate(before_probe(n), action_probe(n))
    end if
    DO i=1,nsize
      bvec(i,i)=(1.d0,0.d0)
      DO j=1,nsize
        amat_spec(ispec)=SUM(CONJG(eigvecs(:,i))*eigvecs(:,j))
        CALL mpro%allgather_inplace(amat_spec)
        amat(i,j)=SUM(amat_spec)*(ritznum(j)-(1.d0,0.d0))
        if (stabilization_probe) gram_probe(i,j)=sum(amat_spec)
      ENDDO
    ENDDO

    if (stabilization_probe) s_probe=amat
    CALL zgesv(nsize,nsize,amat,nsize,ipiv,bvec,nsize,info)
    if (stabilization_probe) then
      call probe_singular_values(gram_probe,gsv(1:nsize))
      call probe_singular_values(s_probe,ssv(1:nsize))
      if (ispec == 0) then
        print *, 'STAB_PROBE G singular_values ',gsv(1:nsize)
        print *, 'STAB_PROBE S singular_values ',ssv(1:nsize)
        print *, 'STAB_PROBE G/S reciprocal_condition2 ', &
          gsv(nsize)/max(gsv(1),tiny(1.0d0)), &
          ssv(nsize)/max(ssv(1),tiny(1.0d0))
        print *, 'STAB_PROBE zgesv info ',info
        do i=1,nsize
          print *, 'STAB_PROBE matrix row/G/S/inverse ',i, &
            gram_probe(i,:),s_probe(i,:),bvec(i,:)
        end do
      end if
      if (info /= 0) error stop 'diagnostic projected inverse failed'
      inverse_defect=probe_inverse_defect(s_probe,bvec)
      if (ispec == 0) then
        print *, 'STAB_PROBE inverse I-SB Frobenius ',inverse_defect
        print *, 'STAB_PROBE scaled inverse residual ',inverse_defect/ &
          max(sqrt(sum(abs(s_probe)**2))*sqrt(sum(abs(bvec)**2))+ &
          sqrt(dble(nsize)),tiny(1.0d0))
      end if
      saved_mode=mode
      mode=2 ! Exact native homogeneous map, source omitted only here.
      full_defect_sq=0.0d0
      do j=1,nsize
        before_probe=eigvecs(:,j) ! Callback fold has native inout semantics.
        call next_iteration(n,before_probe,fnew)
        call probe_action_metrics(eigvecs(:,j),fnew,ritznum(j), &
          defect_norm,vector_norm)
        global_norms(ispec)=defect_norm**2
        call mpro%allgather_inplace(global_norms)
        defect_norm=sqrt(sum(global_norms))
        full_defect_sq=full_defect_sq+defect_norm**2
        global_norms(ispec)=vector_norm**2
        call mpro%allgather_inplace(global_norms)
        vector_norm=sqrt(sum(global_norms))
        call physical_action(eigvecs(:,j),action_probe,entry_bound)
        print *, 'STAB_PROBE Ritz species/mode/L1/L2/B_M ',ispec,j, &
          sum(abs(action_probe)),sqrt(sum(abs(action_probe)**2)),entry_bound
        if (ispec == 0) print *, 'STAB_PROBE full_T defect/vnorm/gap/relative/gap_ratio ',j, &
          defect_norm,vector_norm,abs(ritznum(j)-1.0d0), &
          defect_norm/max(vector_norm,tiny(1.0d0)), &
          defect_norm/max(vector_norm*abs(ritznum(j)-1.0d0),tiny(1.0d0))
      end do
      mode=saved_mode
      ! A posteriori exact-arithmetic bound on ||(T V-V Lambda)V^+||2.
      ! Rounded matvecs/SVD make this an observed indicator, not certification.
      residual_bound=huge(1.0d0)
      if (gsv(nsize)>0.0d0) residual_bound=sqrt(full_defect_sq/gsv(nsize))
      if (ispec == 0) print *, 'STAB_PROBE invariance_bound_indicator ',residual_bound
    end if

    IF(info.NE.0) THEN
      IF(info.GT.0) THEN
        PRINT *,'iterator: singular matrix in zgesv'
      ELSE
        PRINT *,'iterator: argument ',-info,' has illigal value in zgesv'
      ENDIF
      if (present(status)) status = 3
      DEALLOCATE(ritznum,coefren,amat,bvec,ipiv)
      DEALLOCATE(fzero,fold,fnew,eigvecs)
      RETURN
    ENDIF

    ! iterate the solution:
    fold=fzero

    DO iter=1,itermax
      if (stabilization_probe) then
        call probe_stage('fold',n,fold,fold,ispec,evaluate_residual,physical_action)
        stabilization_stage_trace=.true.
      end if
      CALL next_iteration(n,fold,fnew)
      stabilization_stage_trace=.false.
      if (stabilization_probe) then
        before_probe=fnew
        call probe_stage('before_deflation',n,fnew,fold,ispec, &
          evaluate_residual,physical_action)
      end if
      IF(mode.EQ.2) fnew=fnew+fzero
      if (early_accept) then
        ! A valid original-equation candidate needs no approximate Ritz repair.
        ! Both original physical residual and native constraints must pass.
        call check_equation(fnew, converged)
        if (converged .or. solve_status /= 1) then
          print *, 'native original predeflation decision: ', &
            iter,ispec,residual_value,converged,solve_status
          exit
        end if
      end if
      DO j=1,nsize
        coefren_spec(ispec)=SUM(bvec(j,:)                           &
                  *MATMUL(TRANSPOSE(CONJG(eigvecs(:,1:nsize))),fnew-fold))
        CALL mpro%allgather_inplace(coefren_spec)
        coefren(j)=ritznum(j)*SUM(coefren_spec)
      ENDDO
      fnew=fnew-MATMUL(eigvecs(:,1:nsize),coefren)
      if (stabilization_probe) then
        call probe_stage('after_deflation',n,fnew,before_probe,ispec, &
          evaluate_residual,physical_action)
        call probe_stage('stabilized_update',n,fnew,fold,ispec, &
          evaluate_residual,physical_action)
        if (ispec == 0) print *, 'NEO2_DIAGNOSTIC_ONLY stabilization NO_TORQUE'
        call mpro%deinit()
        stop ! No check_equation/report_solve/transport admission on this path.
      end if
      break_cond1(ispec)=SUM(ABS(fnew-fold))
      break_cond2(ispec)=relerr*SUM(ABS(fnew))
      print *,iter,break_cond1(ispec),break_cond2(ispec),ispec
      CALL mpro%allgather_inplace(break_cond1)
      CALL mpro%allgather_inplace(break_cond2)
      call check_equation(fnew, converged)
      if (converged .or. solve_status /= 1) exit
      fold = fnew
      IF(iter.EQ.itermax) PRINT *,'iterator: maximum number of iterations reached'
    ENDDO

    if (ispec .eq. 0) print *,'iterator: number of stabilized iterations = ', completed_iterations
    call publish_status()

    RESULT_=fnew

    DEALLOCATE(coefren,amat,bvec,ipiv,fold,fnew,fzero)

  contains

    subroutine check_equation(candidate, accepted)
      complex(kind=kind(1d0)), intent(in) :: candidate(n)
      logical, intent(out) :: accepted
      double precision :: local_residual, constraints(0:num_spec-1)
      completed_iterations = iter
      if (present(evaluate_residual)) then
        call evaluate_residual(n, candidate, local_residual)
      else
        ! Generic iterator contract: the unaccelerated fixed-point equation.
        ! The kinetic caller supplies the original assembled-equation oracle.
        call next_iteration(n, candidate, fold)
        if (mode == 2) fold = fold + fzero
        local_residual = sum(abs(fold-candidate)) / &
          max(sum(abs(fold))+sum(abs(candidate)), tiny(1.0d0))
      end if
      if (.not. ieee_is_finite(local_residual) .or. local_residual < 0.0d0) then
        local_residual = huge(1.0d0)
      end if
      break_cond1(ispec) = local_residual
      call mpro%allgather_inplace(break_cond1)
      residual_value = maxval(break_cond1)
      if (residual_value == huge(1.0d0)) solve_status = 5
      accepted = residual_value <= relerr
      if (early_accept) then
        ! The callback follows the residual call, so its native moment readback
        ! is for this candidate. Never relabel a moment as an equation residual.
        call evaluate_constraints(n,candidate,constraints(ispec))
        if (.not. ieee_is_finite(constraints(ispec)) .or. constraints(ispec)<0.0d0) &
          constraints(ispec)=huge(1.0d0)
        call mpro%allgather_inplace(constraints)
        if (maxval(constraints)==huge(1.0d0)) solve_status=5
        accepted=accepted .and. maxval(constraints)<=relerr
        print *, 'native original acceptance constraints: ',ispec,maxval(constraints)
      end if
      if (accepted) solve_status = 0
    end subroutine check_equation

    subroutine publish_status()
      if (present(status)) status = solve_status
      if (present(iterations)) iterations = completed_iterations
      if (present(equation_residual)) equation_residual = residual_value
    end subroutine publish_status

  END SUBROUTINE iterator

  subroutine probe_singular_values(a, singular)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    complex(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: singular(:)
    complex(dp) :: copy(size(a,1),size(a,2)), dummy_u(1,1), dummy_v(1,1)
    complex(dp) :: work(8*size(a,1))
    real(dp) :: rwork(5*size(a,1))
    integer :: k, info
    k=size(a,1)
    copy=a
    call zgesvd('N','N',k,k,copy,k,singular,dummy_u,1,dummy_v,1, &
      work,size(work),rwork,info)
    if (info /= 0) error stop 'diagnostic SVD failed'
  end subroutine probe_singular_values

  double precision function probe_inverse_defect(a,b) result(defect)
    implicit none
    complex(kind=kind(1d0)), intent(in) :: a(:,:),b(:,:)
    complex(kind=kind(1d0)) :: error(size(a,1),size(a,1))
    integer :: i
    error=-matmul(a,b)
    do i=1,size(a,1)
      error(i,i)=error(i,i)+1.0d0
    end do
    defect=sqrt(sum(abs(error)**2))
  end function probe_inverse_defect

  subroutine probe_action_metrics(v,tv,lambda,defect,vnorm)
    implicit none
    complex(kind=kind(1d0)), intent(in) :: v(:),tv(:),lambda
    double precision, intent(out) :: defect,vnorm
    defect=sqrt(sum(abs(tv-lambda*v)**2))
    vnorm=sqrt(sum(abs(v)**2))
  end subroutine probe_action_metrics

  subroutine probe_stage(label,n,v,reference,ispec,residual,action)
    implicit none
    character(len=*), intent(in) :: label
    integer, intent(in) :: n,ispec
    complex(kind=kind(1d0)), intent(in) :: v(n),reference(n)
    interface
      subroutine residual(n,v,rho)
        integer, intent(in) :: n
        complex(kind=kind(1d0)), intent(in) :: v(n)
        double precision, intent(out) :: rho
      end subroutine residual
      subroutine action(v,lv,bound)
        complex(kind=kind(1d0)), intent(in) :: v(:)
        complex(kind=kind(1d0)), intent(out) :: lv(:)
        double precision, intent(out) :: bound
      end subroutine action
    end interface
    complex(kind=kind(1d0)) :: lv(n),delta(n),ld(n)
    double precision :: rho,bound,delta_bound
    call residual(n,v,rho) ! Native q, original BC and max-species denominator.
    call action(v,lv,bound)
    delta=v-reference
    call action(delta,ld,delta_bound)
    print *, 'STAB_PROBE stage species/label/rho ',ispec,trim(label),rho
    print *, 'STAB_PROBE vector L1/L2/Lv_L1/Lv_L2/B_M ',sum(abs(v)), &
      sqrt(sum(abs(v)**2)),sum(abs(lv)),sqrt(sum(abs(lv)**2)),bound
    print *, 'STAB_PROBE delta L1/L2/Ldelta_L1/Ldelta_L2/B_M ',sum(abs(delta)), &
      sqrt(sum(abs(delta)**2)),sum(abs(ld)),sqrt(sum(abs(ld)**2)),delta_bound
  end subroutine probe_stage

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

    USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
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
    INTEGER                                       :: m,j,ncount,m_ritz
    INTEGER :: driv_spec
    complex(kind=kind(1d0)),   DIMENSION(:),   ALLOCATABLE :: fold,fnew,ritznum_prev
    complex(kind=kind(1d0)),   DIMENSION(:,:), ALLOCATABLE :: qvecs,hmat,eigh,qvecs_prev
    complex(kind=kind(1d0)),   DIMENSION(:), ALLOCATABLE :: q_spec, h_spec
    ! Norm of the candidate vector before orthogonalization.
    DOUBLE PRECISION :: q_norm
    LOGICAL :: subspace_closed, trace_projection
    character(len=32) :: trace_value
    integer :: trace_status
    double precision :: projection_error

    INTEGER :: m_tol, m_ind
    complex(kind=kind(1d0)), DIMENSION(500) :: ritzum_write

    ALLOCATE(fold(n),fnew(n))
    ALLOCATE(qvecs_prev(n,1),ritznum_prev(mmax))
    ALLOCATE(hmat(mmax+1,mmax))

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
    ierr=0
    q_norm=SQRT(SUM(q_spec))
    IF(.NOT.ieee_is_finite(q_norm)) THEN
      ierr=1
      ngrow=0
      IF(ALLOCATED(eigvecs)) DEALLOCATE(eigvecs)
      ALLOCATE(eigvecs(n,0))
      PRINT *,'arnoldi: initial vector contains non-finite values'
      DEALLOCATE(fold,fnew,qvecs_prev,ritznum_prev,hmat)
      DEALLOCATE(q_spec,h_spec)
      RETURN
    ENDIF
    IF(q_norm.LE.TINY(1.0d0)) THEN
      ngrow=0
      IF(ALLOCATED(eigvecs)) DEALLOCATE(eigvecs)
      ALLOCATE(eigvecs(n,0))
      PRINT *,'arnoldi: initial vector is zero; no Krylov iterations needed'
      DEALLOCATE(fold,fnew,qvecs_prev,ritznum_prev,hmat)
      DEALLOCATE(q_spec,h_spec)
      RETURN
    ENDIF
    qvecs_prev(:,1)=fnew/q_norm
    ncount=0
    trace_projection = .false.
    call get_environment_variable('NEO2_AX_TRACE', trace_value, status=trace_status)
    if (trace_status == 0) trace_projection = trim(trace_value) == '1'

    DO m=1,mmax

      subspace_closed=.FALSE.
      m_ritz=m

      ! Retain one extra vector for the Arnoldi residual.  Every column of
      ! the m-dimensional projection requires its own operator application,
      ! including the final column when the subspace has not closed.
      ALLOCATE(qvecs(n,m+1))
      qvecs(:,1:m)=qvecs_prev(:,1:m)

      ALLOCATE(eigh(m,m))

      fold=qvecs(:,m)
      call next_iteration(n, fold, fnew)
      qvecs(:,m+1)=fnew
      h_spec=0.0d0
      h_spec(ispec)=SUM(CONJG(fnew)*fnew)
      CALL mpro%allgather_inplace(h_spec)
      q_norm=SQRT(SUM(h_spec))
      DO j=1,m
        h_spec=0.0d0
        h_spec(ispec)=SUM(CONJG(qvecs(:,j))*qvecs(:,m+1))
        CALL mpro%allgather_inplace(h_spec)
        hmat(j,m)=SUM(h_spec)
        qvecs(:,m+1)=qvecs(:,m+1)-hmat(j,m)*qvecs(:,j)
      ENDDO
      h_spec=0.0d0
      h_spec(ispec)=SUM(CONJG(qvecs(:,m+1))*qvecs(:,m+1))
      CALL mpro%allgather_inplace(h_spec)
      hmat(m+1,m)=SQRT(SUM(h_spec))
      IF(ABS(hmat(m+1,m)).LE.SQRT(EPSILON(1.0d0))*MAX(1.0d0,q_norm)) THEN
        subspace_closed=.TRUE.
      ELSE
        qvecs(:,m+1)=qvecs(:,m+1)/hmat(m+1,m)
      ENDIF

      CALL try_eigvecvals(m_ritz,tol,hmat(1:m_ritz,1:m_ritz),ngrow, &
        & ritznum(1:m_ritz),eigh,ierr)

      IF(ierr.NE.0) THEN
        PRINT *,'arnoldi: Ritz solve failed at dimension ',m_ritz
        DEALLOCATE(fold,fnew,qvecs,qvecs_prev,ritznum_prev,hmat,eigh)
        DEALLOCATE(q_spec,h_spec)
        RETURN
      ENDIF

      IF(m_ritz.GT.2) THEN

        ! check for convergence of ritznum exceeding tol
        m_tol=m_ritz-1
        DO m_ind = 1,m_ritz-1
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
      ritznum_prev(1:m_ritz)=ritznum(1:m_ritz)

      IF(subspace_closed.OR.ncount.GE.ntol.OR.m.EQ.mmax) THEN
        IF(ALLOCATED(eigvecs)) DEALLOCATE(eigvecs)
        ALLOCATE(eigvecs(n,ngrow))

        eigvecs=MATMUL(qvecs(:,1:m_ritz),eigh(:,1:ngrow))

        if (trace_projection) then
          print *, 'native Arnoldi: dimension ', m_ritz, ' closed ', subspace_closed, &
            ' stable_count ', ncount, ' stability_target ', ntol, &
            ' subdiagonal ', abs(hmat(m_ritz+1,m_ritz)), ' modes ', ngrow
          do j=1,ngrow
            projection_error = sum(abs(matmul(hmat(1:m_ritz,1:m_ritz),eigh(:,j)) &
              -ritznum(j)*eigh(:,j)))
            print *, 'native Ritz: mode ', j, ' projected_defect_l1 ', projection_error, &
              ' Arnoldi_tail ', abs(hmat(m_ritz+1,m_ritz)*eigh(m_ritz,j))
          end do
        end if
        PRINT *,'arnoldi: number of iterations = ',m_ritz
        EXIT
      ENDIF

      DEALLOCATE(qvecs_prev)
      ALLOCATE(qvecs_prev(n,m+1))
      qvecs_prev=qvecs
      DEALLOCATE(qvecs)

      DEALLOCATE(eigh)

    ENDDO

    DEALLOCATE(fold,fnew,qvecs,qvecs_prev,ritznum_prev,hmat)
    IF(ALLOCATED(eigh)) DEALLOCATE(eigh)
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
        PRINT *,'try_eigvecvals: zhseqr failed to compute all eigenvalues, info = ',info
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
