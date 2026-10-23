program test_stabilization_probe
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use arnoldi_mod, only: probe_singular_values, probe_inverse_defect, &
        probe_action_metrics, probe_stage
    implicit none
    complex(dp) :: a(2,2),b(2,2),v(2),tv(2)
    real(dp) :: singular(2),defect,vnorm
    integer :: residual_calls,action_calls
    ! Analytic singular values 3 and 1, complex diagonal phase retained.
    a=(0.0_dp,0.0_dp); b=a
    a(1,1)=(0.0_dp,3.0_dp); a(2,2)=1.0_dp
    b(1,1)=(0.0_dp,-1.0_dp/3.0_dp); b(2,2)=1.0_dp
    call probe_singular_values(a,singular)
    if (maxval(abs(singular-[3.0_dp,1.0_dp]))>1.0e-13_dp) error stop 'SVD'
    if (probe_inverse_defect(a,b)>1.0e-13_dp) error stop 'inverse'
    b(2,2)=1.5_dp
    if (abs(probe_inverse_defect(a,b)-0.5_dp)>1.0e-13_dp) error stop 'bad inverse'
    ! Nonnormal T=[[1,2],[0,2]]: (2,1) is a true eigenvector at lambda2;
    ! (1,1) has true full defect1 although a scalar Rayleigh value can mislead.
    v=[2.0_dp,1.0_dp]; tv=[4.0_dp,2.0_dp]
    call probe_action_metrics(v,tv,cmplx(2.0_dp,0.0_dp,dp),defect,vnorm)
    if (defect>1.0e-13_dp .or. abs(vnorm-sqrt(5.0_dp))>1.0e-13_dp) error stop 'exact mode'
    v=[1.0_dp,1.0_dp]; tv=[3.0_dp,2.0_dp]
    call probe_action_metrics(v,tv,cmplx(2.0_dp,0.0_dp,dp),defect,vnorm)
    if (abs(defect-1.0_dp)>1.0e-13_dp) error stop 'false mode'
    ! Native diagnostic stage must call original residual and both L actions,
    ! including the non-null filter delta; callbacks enforce known vectors.
    residual_calls=0; action_calls=0
    call probe_stage('analytic_half_filter',2,[cmplx(0.5_dp,0.0_dp,dp), &
        cmplx(1.0_dp,0.0_dp,dp)],v,0,known_residual,known_action)
    if (residual_calls/=1 .or. action_calls/=2) error stop 'stage coverage'
    print *, 'All tests passed!'
contains
    subroutine known_residual(n,x,rho)
        integer, intent(in) :: n
        complex(dp), intent(in) :: x(n)
        real(dp), intent(out) :: rho
        if (maxval(abs(x-[0.5_dp,1.0_dp]))>1.0e-13_dp) error stop 'stage vector'
        ! L=diag(2,3), q=(2,3): residual=1/(5+4)=1/9 after density half-filter.
        rho=1.0_dp/9.0_dp
        residual_calls=residual_calls+1
    end subroutine known_residual
    subroutine known_action(x,lx,bound)
        complex(dp), intent(in) :: x(:)
        complex(dp), intent(out) :: lx(:)
        real(dp), intent(out) :: bound
        action_calls=action_calls+1
        if (action_calls==2) then
            if (maxval(abs(x-[-0.5_dp,0.0_dp]))>1.0e-13_dp) error stop 'filter delta'
        end if
        lx=[2.0_dp,3.0_dp]*x
        bound=sum(abs(lx))
    end subroutine known_action
end program test_stabilization_probe
