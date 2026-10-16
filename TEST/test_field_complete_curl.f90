program test_field_complete_curl
    ! Fixture coefficients and gates frozen before candidate execution.
    ! See DOC/complete-covariant-curl.md; no field-equation oracle is copied here.
    use iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use neo_field_provider, only: neo_field_load, neo_field_prepare, neo_field_free
    use magfie_mod, only: magfie
    use partpa_mod, only: bmod0
    use nrtype, only: pi
    implicit none
    integer, parameter :: grids(3) = [65,129,257]
    real(dp), parameter :: points(2,3) = reshape( &
        [0.37_dp,0.79_dp,1.21_dp,2.11_dp,4.33_dp,5.27_dp], [2,3])
    real(dp), parameter :: relative_gate=1.0e-8_dp, floor=1.0e-8_dp
    real(dp) :: first(3,6), actual(3,6), final_error, material, errors(3)
    integer :: orientation, iota_sign, grid

    do orientation=-1,1,2
        do iota_sign=-1,1,2
            call write_chart('complete_curl.bc',orientation,iota_sign,2.0_dp)
            do grid=1,size(grids)
                call check_grid('complete_curl.bc',grids(grid),2.0_dp,actual,final_error,material)
                errors(grid)=final_error
                print '(a,3i5,a,es12.4,a,es12.4)', 'curl grid/signs ', &
                    grids(grid),orientation,iota_sign,' scaled error=',final_error, &
                    ' omitted-C/gate=',material
            end do
            if (material <= 100.0_dp) error stop 'ineffective radial covariant fixture'
            if (errors(3)>relative_gate) &
                error stop 'full covariant curl interpolation/derivative consistency'
            first=actual
            call write_chart('complete_curl_second.bc',orientation,iota_sign,3.0_dp)
            call check_grid('complete_curl_second.bc',129,3.0_dp,actual,final_error,material)
            call check_grid('complete_curl.bc',257,2.0_dp,actual,final_error,material)
            if (maxval(abs(actual-first))>1.0e-12_dp) &
                error stop 'stale full-curl profile/grid cache'
        end do
    end do
    call neo_field_free()
    print *, 'All complete provider curl tests passed!'
contains
    subroutine check_grid(path,n,field_t,central_curl,max_error,material)
        character(*), intent(in) :: path
        integer, intent(in) :: n
        real(dp), intent(in) :: field_t
        real(dp), intent(out) :: central_curl(3,6),max_error,material
        real(dp) :: location(3,6), jac, ignored(3), field(3,6,5)
        real(dp) :: dfield(3,3,6), d_half(3,3,6), plus(3),minus(3),pp(3),mm(3)
        real(dp) :: x(3), expected(3), expected_half(3), h
        real(dp) :: saved_jac(6),c_partial(3),scale(3),magnitude,expected_b,hc(3)
        real(dp) :: point_error,largest_step_error,largest_actual,largest_expected
        integer :: worst_point,worst_component,component
        integer :: point,axis,shift,step
        real(dp), parameter :: shifts(5)=[0.0_dp,1.0_dp,-1.0_dp,0.5_dp,-0.5_dp]
        call neo_field_load(path,9,theta_points=n,phi_points=n)
        bmod0=1.0e-4_dp
        location(1,:)=0.5_dp
        location(2:3,1:3)=points
        location(2:3,4:6)=real(nint(points*(n-1)/(2*pi)),dp)*2*pi/(n-1)
        max_error=0; material=0; largest_step_error=0
        do step=1,2
            h=1.0e-3_dp/2**(step-1)
            ! Rebuild a radial cache only five times per step, not for each angle.
            do shift=1,5
                call neo_field_prepare(0.5_dp+h*shifts(shift))
                if (shift==1) then
                    x=[0.5_dp,0.0_dp,0.0_dp]
                    call evaluate(x,hc,jac,ignored,magnitude)
                    expected_b=field_t*sqrt(1+(0.5_dp*sqrt(0.5_dp) &
                        *0.2_dp*(1+0.1_dp*0.5_dp)/5)**2) &
                        *(1-0.12_dp*sqrt(0.5_dp))
                    if (abs(magnitude-expected_b)>1.0e-12_dp*expected_b) &
                        error stop 'profile reload or field unit normalization'
                end if
                do point=1,6
                    x=location(:,point);x(1)=x(1)+h*shifts(shift)
                    call evaluate(x,field(:,point,shift),jac,ignored)
                    if (shift==1) then
                        central_curl(:,point)=ignored
                        saved_jac(point)=jac
                    end if
                end do
            end do
            do point=1,6
                dfield(:,1,point)=(4*(field(:,point,4)-field(:,point,5))/h &
                    -(field(:,point,2)-field(:,point,3))/(2*h))/3
            end do
            call neo_field_prepare(0.5_dp)
            do point=1,6
                do axis=2,3
                    x=location(:,point);x(axis)=x(axis)+h
                    call evaluate(x,plus,jac,ignored)
                    x(axis)=x(axis)-2*h
                    call evaluate(x,minus,jac,ignored)
                    x(axis)=x(axis)+1.5_dp*h
                    call evaluate(x,pp,jac,ignored)
                    x(axis)=x(axis)-h
                    call evaluate(x,mm,jac,ignored)
                    dfield(:,axis,point)=(4*(pp-mm)/h-(plus-minus)/(2*h))/3
                end do
                expected=[dfield(3,2,point)-dfield(2,3,point), &
                    dfield(1,3,point)-dfield(3,1,point), &
                    dfield(2,1,point)-dfield(1,2,point)]/saved_jac(point)
                scale=max(abs(expected),floor)
                component=maxloc(abs(central_curl(:,point)-expected)/scale,dim=1)
                point_error=abs(central_curl(component,point)-expected(component))/scale(component)
                if (point_error>=max_error) then
                    max_error=point_error
                    worst_point=point
                    worst_component=component
                    largest_actual=central_curl(component,point)
                    largest_expected=expected(component)
                end if
                c_partial=[0.0_dp,dfield(1,3,point),-dfield(1,2,point)]/saved_jac(point)
                material=max(material,maxval(abs(c_partial)/(relative_gate*scale)))
                if (step==2) then
                    expected_half=[d_half(3,2,point)-d_half(2,3,point), &
                        d_half(1,3,point)-d_half(3,1,point), &
                        d_half(2,1,point)-d_half(1,2,point)]/saved_jac(point)
                    largest_step_error=max(largest_step_error, &
                        maxval(abs(expected-expected_half)/scale))
                end if
            end do
            d_half=dfield
        end do
        print *, 'FD step scaled difference:',largest_step_error
        if (max_error>relative_gate) then
            print *, 'worst grid/point/component:',n,worst_point,worst_component
            print *, 'computed and independent FD:',largest_actual,largest_expected
        end if
        if (n==grids(size(grids))) then
            if (largest_step_error>relative_gate) &
                error stop 'full-curl independent FD step not converged'
        end if
        call neo_field_free()
    end subroutine check_grid
    subroutine evaluate(x,hc,jac,curl,magnitude)
        real(dp),intent(in)::x(3)
        real(dp),intent(out)::hc(3),jac,curl(3)
        real(dp),intent(out),optional::magnitude
        real(dp)::b,bd(3),ht(3)
        call magfie(x,b,jac,bd,hc,ht,curl)
        if (.not. all(ieee_is_finite([b,jac,hc,curl]))) &
            error stop 'nonfinite provider oracle'
        if (present(magnitude)) magnitude=b
    end subroutine evaluate
    subroutine write_chart(path,orientation,iota_sign,b0)
        character(*),intent(in)::path
        integer,intent(in)::orientation,iota_sign
        real(dp),intent(in)::b0
        integer::unit,j,k,handedness
        real(dp)::s,r,iota,g,i,flux,b
        handedness=-orientation
        flux=-real(handedness,dp)*pi*0.5_dp**2*b0
        open(newunit=unit,file=path,status='replace')
        do k=1,5
            write(unit,'(a)') 'CC frozen manufactured 3D full-curl oracle'
        end do
        write(unit,*) 1,1,7,1,flux,0.5_dp,5.0_dp
        do j=0,6
            s=0.2_dp+0.1_dp*j;r=0.5_dp*sqrt(s)
            iota=iota_sign*0.2_dp*(1+0.1_dp*s)
            g=b0*5*(1+0.05_dp*s)
            i=b0*0.5_dp**2*iota*s/5
            b=b0*sqrt(1+(r*iota/5)**2)
            write(unit,'(a)') 'CC s iota Jpol Itor pprime sqrtg00'
            write(unit,'(a)') 'CC units A A Pa m3'
            write(unit,*) s,iota,g*1.0e6_dp/(0.2_dp*handedness), &
                i*1.0e6_dp/(0.2_dp*handedness),0.0_dp,0.0_dp
            write(unit,'(a)') 'CC m n rmnc rmns zmnc zmns vmnc vmns bmnc bmns'
            write(unit,*) 0,-1,(0.0_dp,k=1,8)
            write(unit,*) 0,0,5.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                0.0_dp,0.0_dp,b,0.0_dp
            write(unit,*) 0,1,(0.0_dp,k=1,8)
            ! Format9 evaluates m*theta+n*phi; n=-1 realizes frozen theta-phi.
            write(unit,*) 1,-1,0.18_dp*r,0.0_dp,0.0_dp,orientation*0.11_dp*r, &
                0.0_dp,0.0_dp,0.0_dp,0.06_dp*sqrt(s)*b
            write(unit,*) 1,0,r,0.0_dp,0.0_dp,orientation*r, &
                0.0_dp,0.0_dp,-0.12_dp*sqrt(s)*b,0.0_dp
            write(unit,*) 1,1,(0.0_dp,k=1,8)
        end do
        close(unit)
    end subroutine write_chart
end program test_field_complete_curl
