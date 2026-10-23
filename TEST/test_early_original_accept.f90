! Independent physical M is singular; native homogeneous T is nonnormal.
! A deliberately imperfect reused Ritz pair reproduces destructive subtraction.
program test_early_original_accept
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan,ieee_value
    use arnoldi_mod, only: iterator,mode,ngrow,eigvecs,ritznum,fzero
    use collisionality_mod, only: num_spec
    use mpiprovider_module, only: mpro
    implicit none
    complex(kind=kind(1d0)) :: q(2),solution(2),legacy(2),scale
    double precision :: rho,legacy_rho
    integer :: rank,status,iterations,legacy_status
    logical :: bad_moment,bad_residual,finite_bad_moment
    call mpro%init()
    num_spec=mpro%getNumProcs();rank=mpro%getRank()
    scale=cmplx(1.0d0,0.5d0,kind(1d0))
    q=[scale,2.0d0*scale]
    bad_moment=.false.;bad_residual=.false.;finite_bad_moment=.false.
    ! M=[[0,1],[0,2]], K=diag(1,2); q=(1,2)*scale.
    ! T=I-K^-1 M=[[1,-1],[0,0]], nonnormal with exact unit eigenvalue.
    ! Raw h=q+K^-1(q-Mq)=(0,scale) is an actual M solution and has zero x moment.
    call retained_pair()
    legacy=q
    call iterator(1,2,2,1.0d-5,1,legacy,rank,next_step,status,iterations,rho,residual)
    legacy_status=status;legacy_rho=rho
    if (status==0 .or. rho<0.9d0) error stop 'old destructive correction not reproduced'
    call retained_pair()
    solution=q
    call iterator(1,2,2,1.0d-5,1,solution,rank,next_step,status,iterations,rho,residual, &
        early_original_accept=.false.)
    if (status/=legacy_status .or. rho/=legacy_rho .or. any(solution/=legacy)) &
        error stop 'disabled/default changed'
    call retained_pair()
    solution=q
    call iterator(1,2,2,1.0d-5,1,solution,rank,next_step,status,iterations,rho,residual, &
        early_original_accept=.true.,evaluate_constraints=moment)
    if (status/=0 .or. iterations/=1 .or. rho>1.0d-5) error stop 'valid M candidate rejected'
    if (abs(solution(1))>1.0d-13 .or. abs(solution(2)-scale)>1.0d-13) &
        error stop 'wrong physical solution'
    ! A finite failing constraint on one species must prevent early acceptance.
    finite_bad_moment=.true.;call retained_pair();solution=q
    call iterator(1,2,2,1.0d-5,1,solution,rank,next_step,status,iterations,rho,residual, &
        early_original_accept=.true.,evaluate_constraints=moment)
    if (status==0) error stop 'one species moment failed but accepted'
    finite_bad_moment=.false.
    ! Invalid physical residual OR moment on one rank rejects every rank.
    bad_residual=.true.;call retained_pair();solution=q
    call iterator(1,2,2,1.0d-5,1,solution,rank,next_step,status,iterations,rho,residual, &
        early_original_accept=.true.,evaluate_constraints=moment)
    if (status/=5) error stop 'nonfinite residual accepted'
    bad_residual=.false.;bad_moment=.true.;call retained_pair();solution=q
    call iterator(1,2,2,1.0d-5,1,solution,rank,next_step,status,iterations,rho,residual, &
        early_original_accept=.true.,evaluate_constraints=moment)
    if (status/=5) error stop 'nonfinite moment accepted'
    bad_moment=.false.;call retained_pair();solution=q
    call iterator(1,2,2,1.0d-5,1,solution,rank,next_step,status,iterations,rho,residual, &
        early_original_accept=.true.)
    if (status/=4) error stop 'missing constraint authority accepted'
    call iterator(3,2,2,1.0d-5,1,solution,rank,next_step)
    call mpro%deinit()
    print *, 'All tests passed!'
contains
    subroutine retained_pair()
        mode=1;ngrow=1
        if (allocated(eigvecs)) deallocate(eigvecs)
        if (allocated(ritznum)) deallocate(ritznum)
        allocate(eigvecs(2,1),ritznum(1))
        eigvecs(:,1)=[cmplx(1.0d0,0.0d0,kind(1d0)),cmplx(1.0d-6,0.0d0,kind(1d0))]
        ritznum(1)=cmplx(1.0d0+1.0d-12,0.0d0,kind(1d0))
    end subroutine retained_pair
    subroutine next_step(n,input,output)
        integer :: n
        complex(kind=kind(1d0)) :: input(n),output(n)
        output=[input(1)-input(2),cmplx(0.0d0,0.0d0,kind(1d0))]
        if (mode/=2) output=output+[fzero(1),0.5d0*fzero(2)]
    end subroutine next_step
    subroutine residual(n,x,value)
        integer, intent(in) :: n
        complex(kind=kind(1d0)), intent(in) :: x(n)
        double precision, intent(out) :: value
        complex(kind=kind(1d0)) :: mx(2)
        mx=[x(2),2.0d0*x(2)]
        value=sum(abs(q-mx))/max(sum(abs(q))+sum(abs(mx)),tiny(1.0d0))
        if (bad_residual .and. rank==num_spec-1) &
            value=ieee_value(0.0d0,ieee_quiet_nan)
    end subroutine residual
    subroutine moment(n,x,value)
        integer, intent(in) :: n
        complex(kind=kind(1d0)), intent(in) :: x(n)
        double precision, intent(out) :: value
        value=abs(x(1))/max(maxval(abs(x)),tiny(1.0d0))
        if (finite_bad_moment .and. rank==num_spec-1) value=0.1d0
        if (bad_moment .and. rank==num_spec-1) value=ieee_value(0.0d0,ieee_quiet_nan)
    end subroutine moment
end program test_early_original_accept
