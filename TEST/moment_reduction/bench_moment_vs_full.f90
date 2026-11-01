! Local synthetic benchmark of the actual collision U,C,V implementation.
! Both formulations use the SAME P^-1 and original L = P - UCV.
module moment_vs_full_bench_state
    use, intrinsic :: iso_fortran_env, only: real64
    use collision_moment_map
    use moment_schur_solver
    implicit none
    integer, parameter :: dp = real64
    type(collision_moment_map_t) :: map
    integer :: ns, nphi, npass, nenergy, nleg, local_n, nm, n, r
    complex(dp), allocatable :: workspace_m(:, :), workspace_c(:, :)
    real(dp), allocatable :: coeff(:, :, :, :, :), diag(:)
contains
    subroutine build_case(nphi_in, npass_in, ne_in, nl_in, ns_in, gain)
        integer, intent(in) :: nphi_in, npass_in, ne_in, nl_in, ns_in
        real(dp), intent(in) :: gain
        integer, allocatable :: passing(:)
        real(dp), allocatable :: bra(:, :, :), ket(:, :, :), fwd(:), bwd(:)
        integer :: s,t,i,j,l,p,info
        ns=ns_in; nphi=nphi_in; npass=npass_in; nenergy=ne_in; nleg=nl_in
        allocate(passing(nphi),bra(nleg,npass+1,nphi),ket(nleg,npass+1,nphi),fwd(nphi),bwd(nphi))
        passing=npass
        fwd=0.85_dp; bwd=0.82_dp
        do i=1,nphi
            do l=1,nleg
                do p=1,npass+1
                    bra(l,p,i) = cos(0.091_dp*p+0.31_dp*l+0.04_dp*i)/sqrt(real(npass+1,dp))
                    ket(l,p,i) = sin(0.087_dp*p+0.37_dp*l+0.06_dp*i)/sqrt(real(npass+1,dp))
                end do
            end do
        end do
        call init_collision_moment_map(map,nenergy,passing,bra,ket,fwd,bwd,.true.,info)
        if(info /= 0) error stop 'bad map'
        local_n=map%field_size(); nm=map%moment_size(); n=local_n*ns; r=nm*ns
        allocate(coeff(nenergy,nenergy,nleg,ns,ns),workspace_m(nm,ns),workspace_c(nm,ns),diag(n))
        coeff=0
        do s=1,ns
            do t=1,ns
                do l=1,nleg
                    do j=1,nenergy
                        do i=1,nenergy
                            if(s == t .and. i==j) then
                                coeff(i,j,l,s,t)=gain*(0.72_dp + 0.04_dp*l)
                            else
                                coeff(i,j,l,s,t)=gain*0.028_dp*cos(real(3*i+5*j+l+2*s+t,dp))
                            end if
                        end do
                    end do
                end do
            end do
        end do
        do i=1,n
            diag(i)=2.0_dp+0.1_dp*sin(0.09_dp*i)
        end do
    end subroutine
    subroutine vin(x,y,info)
        complex(dp),intent(in):: x(:)
        complex(dp),intent(out):: y(:)
        integer,intent(out):: info
        integer :: i
        do i=1,ns
            call map%project(x((i-1)*local_n+1:i*local_n), &
                  y((i-1)*nm+1:i*nm), info)
            if(info/=0) return
        end do
    end subroutine
    subroutine ucin(x,y,info)
        complex(dp),intent(in):: x(:)
        complex(dp),intent(out):: y(:)
        integer,intent(out):: info
        integer :: i
        workspace_m=reshape(x,shape(workspace_m))
        call mix_collision_moments(coeff,workspace_m,workspace_c,info)
        if(info/=0) return
        do i=1,ns
            call map%lift(workspace_c(:,i),y((i-1)*local_n+1:i*local_n),info)
            if(info/=0) return
        end do
    end subroutine
    subroutine collision(x,y,info)
        complex(dp),intent(in):: x(:)
        complex(dp),intent(out):: y(:)
        integer,intent(out):: info
        complex(dp) :: moments(r)
        call vin(x,moments,info)
        if(info/=0) return
        call ucin(moments,y,info)
    end subroutine
    subroutine identity(x,y,info)
        complex(dp),intent(in):: x(:)
        complex(dp),intent(out):: y(:)
        integer,intent(out):: info
        y=x;info=0
    end subroutine
    subroutine pinv(x,y,info)
        complex(dp),intent(in):: x(:)
        complex(dp),intent(out):: y(:)
        integer,intent(out):: info
        integer::i,t
        do t=0,ns-1
            y(t*local_n+1)=x(t*local_n+1)/diag(t*local_n+1)
            do i=t*local_n+2,(t+1)*local_n
                y(i)=(x(i)+0.3_dp*y(i-1))/diag(i)
            end do
        end do
        info=0
    end subroutine
    subroutine original(x,y,info)
        complex(dp),intent(in):: x(:)
        complex(dp),intent(out):: y(:)
        integer,intent(out):: info
        integer:: i,t
        call collision(x,y,info)
        if(info/=0) return
        y=-y
        do t=0,ns-1
            y(t*local_n+1)=y(t*local_n+1)+diag(t*local_n+1)*x(t*local_n+1)
            do i=t*local_n+2,(t+1)*local_n
                y(i)=y(i)+diag(i)*x(i)-0.3_dp*x(i-1)
            end do
        end do
    end subroutine
end module

program benchmark
    use moment_vs_full_bench_state
    implicit none
    type(moment_schur_operator_t) :: op
    type(moment_schur_options_t):: opts
    type(moment_schur_report_t):: rep
    complex(dp), allocatable:: truth(:), q(:), full(:), reduced(:)
    complex(dp):: gempty(0), lempty(0)
    real(dp):: t1,t2,tf(12),tr(12),err
    integer:: mode,iter,round,caseid,numruns
    integer :: NPs(3), Phis(3)
    real(dp) :: Gains(3)
    NPs=[16,50,90]; Phis=[36,64,100]; Gains=[0.7_dp,1.15_dp,1.2_dp]
    numruns=8
    print '(a)', 'case,nfull,nmoment,method,runs,median_cpu_s,min_cpu_s,max_cpu_s,steps,applications,inv_calls,residual_norm,field_difference'
    do caseid=1,3
        call build_case(Phis(caseid),NPs(caseid),4,4,2,Gains(caseid))
        allocate(truth(n),q(n),full(n),reduced(n))
        do iter=1,n
            truth(iter)=cmplx(sin(real(iter,dp)*0.017_dp),cos(real(iter,dp)*0.031_dp),dp)
        end do
        call original(truth,q,iter)
        if(iter/=0) error stop 'original failed'
        opts%rtol=2.0e-10_dp; opts%atol=0; opts%restart=20; opts%maxiter=120
        do mode=1,2
            op%nfield=n
            if(mode==1) then
                op%nmoment=n
                op%moments=>identity
                op%feedback=>collision
            else
                op%nmoment=r
                op%moments=>vin
                op%feedback=>ucin
            end if
            op%inverse_p=>pinv
            op%original=>original
            ! Warmup is not counted. Both methods use the same original RHS.
            if(mode==1) then
                call solve_moment_schur(op,q,gempty,full,lempty,opts,rep)
            else
                call solve_moment_schur(op,q,gempty,reduced,lempty,opts,rep)
            end if
            if(rep%status/=schur_ok) error stop 'warmup failed'
            do round=1,numruns
                call cpu_time(t1)
                if(mode==1) then
                    call solve_moment_schur(op,q,gempty,full,lempty,opts,rep)
                else
                    call solve_moment_schur(op,q,gempty,reduced,lempty,opts,rep)
                end if
                call cpu_time(t2)
                if(rep%status/=schur_ok) then
                    print *, 'FAIL: case=',caseid,' mode=',mode,' status=',rep%status,' residual=',rep%original_residual
                    error stop 'no fair comparison'
                end if
                if(mode==1) tf(round)=t2-t1
                if(mode==2) tr(round)=t2-t1
            end do
            err=maxval(abs((merge(full,reduced,mode==1)-truth)))
            if(err>1e-8_dp) error stop 'bad recovery'
            if(mode==1) then
                call report_times('full',tf(1:numruns),err)
            else
                call report_times('moment',tr(1:numruns),err)
            end if
        end do
        print '(a,i0,a,f9.3,a,es10.3)', 'ratio case ', caseid, ': full/moment=', &
           median_time(tf(1:numruns))/median_time(tr(1:numruns)), ' max_field_diff=',maxval(abs(full-reduced))
        deallocate(truth,q,full,reduced,diag,coeff,workspace_m,workspace_c)
    end do
contains
    real(dp) function median_time(times) result(median)
        real(dp), intent(in):: times(:)
        real(dp)::tmp(size(times)),value
        integer :: i,j
        tmp=times
        do i=2,size(tmp)
            value=tmp(i);j=i-1
            do while(j>=1)
                if(tmp(j)<=value) exit
                tmp(j+1)=tmp(j);j=j-1
            end do
            tmp(j+1)=value
        end do
        median=0.5_dp*(tmp(size(tmp)/2)+tmp(size(tmp)/2+1))
    end function
    subroutine report_times(label,times,error)
        character(*),intent(in)::label
        real(dp),intent(in)::times(:),error
        real(dp):: tmp(size(times)),h
        integer::i,j
        tmp=times
        do i=2,size(tmp)
            h=tmp(i);j=i-1
            do while(j>=1)
                if(tmp(j)<=h) exit
                tmp(j+1)=tmp(j);j=j-1
            end do
            tmp(j+1)=h
        end do
        print '(i0,a,i0,a,i0,a,a,a,i0,a,f10.6,a,f10.6,a,f10.6,a,i0,a,i0,a,i0,a,es10.3,a,es10.3)', &
           caseid,',',n,',',r,',',label,',',size(times),',',0.5_dp*(tmp(size(tmp)/2)+tmp(size(tmp)/2+1)),',', &
           minval(tmp),',',maxval(tmp),',',rep%iterations,',',rep%applications,',',rep%inverse_calls,',', &
           rep%original_residual,',',error
    end subroutine
end program benchmark
