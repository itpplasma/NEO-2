program test_electromagnetic_work_geometry
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use nrtype, only: dp
    use mag_interface_mod, only: mag_coordinates, boozer_s, boozer_bmod0, &
        make_magnetics
    use neo_control, only: lab_swi, inp_swi
    use neo_input, only: psi_pr, n0b
    use neo_magfie, only: neo_magfie_calc, bmod0, boozer_psi_pr_hat, &
        boozer_isqrg, boozer_curr_pol_hat, boozer_curr_tor_hat, boozer_iota, &
        curr_pol_array, curr_tor_array, curr_pol_s_array, curr_tor_s_array
    use electromagnetic_work_mod, only: work_geometry_supported, &
        work_lab10_geometry_supported, work_geometry_scale
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), a = 0.62_dp, rmajor = 6.2_dp
    real(dp), parameter :: bzero = 5.3_dp, surface = 0.09_dp
    real(dp), parameter :: angles(4) = [0.3_dp, 0.9_dp, 2.2_dp, 4.8_dp]
    real(dp) :: x(3), b, jac, bd(3), hc(3), ht(3), curl(3)
    real(dp) :: original_curl(3), original_h(3), original_jac
    real(dp) :: scale, native_source, physical_source, cboozer
    real(dp) :: es(3), et(3), ep(3), gs(3), gt(3), h(3), jd, radius, theta
    real(dp) :: polarity, zsign
    integer :: expected_lab, k
    character(len=8) :: argument

    call get_command_argument(1, argument)
    read(argument, *) expected_lab
    polarity = 1.0_dp
    if (expected_lab == 10) polarity = -1.0_dp
    zsign = -polarity
    call write_axis()
    mag_coordinates = 1
    boozer_s = surface
    call make_magnetics('magnetic-work signed geometry regression')
    bmod0 = boozer_bmod0
    if (lab_swi /= expected_lab .or. inp_swi /= 9 .or. n0b /= 0) &
        error stop 'FAIL: native input branch changed'

    do k = 1, size(angles)
        ! The positive-field partner reverses theta as well as the field.
        ! Physical cylindrical phi has the same orientation in both cases.
        theta = -polarity*angles(k)
        x = [surface, 0.0_dp, theta]
        call neo_magfie_calc(x, b, jac, bd, hc, ht, curl)
        if (expected_lab == 10) then
            if (.not. lab10_supported()) error stop 'FAIL: signed lab10 rejected'
            call check_rejections()
        else
            if (.not. work_geometry_supported(1, lab_swi, jac, ht(2), &
                    boozer_curr_pol_hat, boozer_iota, boozer_psi_pr_hat, &
                    0.01_dp, 100.0_dp)) error stop 'FAIL: positive branch rejected'
        end if
        scale = work_geometry_scale(boozer_iota, boozer_psi_pr_hat, &
            boozer_curr_pol_hat, boozer_isqrg)
        if (scale*polarity <= 0.0_dp) error stop 'FAIL: native scale sign'
        cboozer = boozer_curr_pol_hat + boozer_iota*boozer_curr_tor_hat
        ! The actual axisymmetric Arnoldi source, lines 2495/2500, uses
        ! dlogB/dphi=iota*dlogB/dtheta and the full signed scale.
        native_source = -scale*boozer_curr_pol_hat/cboozer &
            *(boozer_iota*bd(3))*(b/bmod0)

        ! Independent Cartesian map and physical signed-flux B direction.
        ! The manufactured circular scalar B is not force balanced: compare
        ! source signs, not amplitudes or a converged TC24 response.
        radius = a*sqrt(surface)
        es = [a/(2.0_dp*sqrt(surface))*cos(theta), 0.0_dp, &
            zsign*a/(2.0_dp*sqrt(surface))*sin(theta)]
        et = [-radius*sin(theta), 0.0_dp, zsign*radius*cos(theta)]
        ep = [0.0_dp, rmajor + radius*cos(theta), 0.0_dp]
        jd = dot_product(es, cross(et, ep))
        gs = cross(et, ep)/jd
        gt = cross(ep, es)/jd
        h = psi_pr/abs(jd)*(ep + boozer_iota*et)
        h = h/sqrt(sum(h**2))
        physical_source = dot_product(gs, cross(h, bd(3)*gt))*1.0e-4_dp
        if (native_source*physical_source <= 0.0_dp) &
            error stop 'FAIL: physical signed source disagrees'
        print '(a,i2,a,f5.2,a,3es16.7)', 'lab ', lab_swi, ', theta=', theta, &
            ', scale/source/Cartesian sign oracle:', &
            scale, native_source, physical_source

        if (expected_lab == 10) then
            ! Fixed Cartesian map, physical field reversal: all signed
            ! covariant B components and signed flux reverse together.
            ! Exercise native curl parity, without asserting its absolute
            ! orientation is the work-source convention.
            original_curl = curl
            original_h = ht
            original_jac = jac
            call reverse_signed_field()
            call neo_magfie_calc(x, b, jac, bd, hc, ht, curl)
            call close('curl reversal', curl, -original_curl)
            call close('field reversal', ht, -original_h)
            call close('fixed-map Jacobian', [jac], [original_jac])
            call reverse_signed_field()
        end if
    end do
    print *, 'All tests passed!'

contains

    logical function lab10_supported() result(supported)
        supported = work_lab10_geometry_supported(1, lab_swi, inp_swi, n0b, &
            jac, ht(2), boozer_curr_pol_hat, boozer_curr_tor_hat, &
            boozer_iota, boozer_psi_pr_hat, boozer_isqrg, b/bmod0, &
            0.01_dp, 100.0_dp)
    end function lab10_supported

    subroutine check_rejections()
        real(dp) :: saved

        saved = jac
        jac = -jac
        if (lab10_supported()) error stop 'FAIL: wrong Jacobian admitted'
        jac = saved*1.01_dp
        if (lab10_supported()) error stop 'FAIL: inconsistent Jacobian admitted'
        jac = saved
        saved = ht(2)
        ht(2) = saved*1.01_dp
        if (lab10_supported()) error stop 'FAIL: inconsistent hphi admitted'
        ht(2) = saved
        saved = boozer_psi_pr_hat
        boozer_psi_pr_hat = -saved
        if (lab10_supported()) error stop 'FAIL: absolute flux admitted'
        boozer_psi_pr_hat = saved
        saved = boozer_iota
        boozer_iota = 0.0_dp
        if (lab10_supported()) error stop 'FAIL: zero iota admitted'
        boozer_iota = saved
        saved = boozer_isqrg
        boozer_isqrg = ieee_value(0.0_dp, ieee_quiet_nan)
        if (lab10_supported()) error stop 'FAIL: nonfinite geometry admitted'
        boozer_isqrg = saved
        n0b = 1
        if (lab10_supported()) error stop 'FAIL: 3D background admitted'
        n0b = 0
        if (work_lab10_geometry_supported(1, 10, 8, 0, jac, ht(2), &
                boozer_curr_pol_hat, boozer_curr_tor_hat, boozer_iota, &
                boozer_psi_pr_hat, boozer_isqrg, b/bmod0, 0.01_dp, 100.0_dp)) &
            error stop 'FAIL: untested input format admitted'
    end subroutine check_rejections

    subroutine reverse_signed_field()
        psi_pr = -psi_pr
        curr_pol_array = -curr_pol_array
        curr_tor_array = -curr_tor_array
        curr_pol_s_array = -curr_pol_s_array
        curr_tor_s_array = -curr_tor_s_array
    end subroutine reverse_signed_field

    pure function cross(u, v) result(w)
        real(dp), intent(in) :: u(3), v(3)
        real(dp) :: w(3)
        w = [u(2)*v(3) - u(3)*v(2), u(3)*v(1) - u(1)*v(3), &
            u(1)*v(2) - u(2)*v(1)]
    end function cross

    subroutine close(label, actual, expected)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: actual(:), expected(:)
        if (maxval(abs(actual - expected)) > &
                2.0e-12_dp*max(1.0e-30_dp, maxval(abs(expected)))) then
            print *, 'FAIL: ', label, actual, expected
            error stop 'Native field reversal parity failed'
        end if
    end subroutine close

    subroutine write_axis()
        real(dp), parameter :: surfaces(8) = &
            [0.01_dp, 0.04_dp, 0.09_dp, 0.16_dp, 0.25_dp, 0.49_dp, 0.81_dp, 1.0_dp]
        real(dp) :: s, r, gboozer, iboozer
        integer :: unit, j

        open(newunit=unit, file='axis.bc', status='replace', action='write')
        do j = 1, 5
            write(unit, '(a)') 'CC manufactured orientation; no force-balance claim'
        end do
        write(unit, '(4i4,3es24.16)') 1, 0, 8, 1, &
            polarity*pi*bzero*a**2, a, rmajor
        do j = 1, size(surfaces)
            s = surfaces(j)
            r = a*sqrt(s)
            gboozer = polarity*bzero*rmajor
            iboozer = 0.5_dp*bzero*r**2/rmajor
            write(unit, '(a)') 'CC s iota Jpol Itor pprime sqrtg00'
            write(unit, '(a)') 'CC native input 9'
            write(unit, '(6es24.16)') s, polarity*0.5_dp, &
                -gboozer/2.0e-7_dp, -iboozer/2.0e-7_dp, 0.0_dp, 0.0_dp
            write(unit, '(a)') 'CC m n Rcos Rsin Zcos Zsin Lcos Lsin Bcos Bsin'
            write(unit, '(2i4,8es24.16)') 0, 0, rmajor, 0.0_dp, 0.0_dp, &
                0.0_dp, 0.0_dp, 0.0_dp, bzero, 0.0_dp
            write(unit, '(2i4,8es24.16)') 1, 0, r, 0.0_dp, 0.0_dp, &
                zsign*r, 0.0_dp, 0.0_dp, -bzero*r/rmajor, 0.0_dp
        end do
        close(unit)
    end subroutine write_axis
end program test_electromagnetic_work_geometry
