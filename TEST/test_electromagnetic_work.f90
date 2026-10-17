program test_electromagnetic_work
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use electromagnetic_work_mod, only: work_back_geometry, &
        work_pitch_convolution, compute_electromagnetic_work, &
        work_geometry_supported
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), c = 2.9979e10_dp
    real(dp), parameter :: elementary_charge = 4.8032e-10_dp
    real(dp), parameter :: bref = 2.3e4_dp, g = 0.012_dp
    real(dp), parameter :: eps = 0.2_dp, kappa = 0.021_dp
    real(dp) :: scale
    real(dp), parameter :: mass(0:2) = &
        [3.343583719e-24_dp, 5.0073567e-24_dp, 9.109382150e-28_dp]
    real(dp), parameter :: charge(0:2) = [1.0_dp, 1.0_dp, -1.0_dp]
    real(dp), parameter :: density(0:2) = [4.0e12_dp, 6.0e12_dp, 1.0e13_dp]
    real(dp), parameter :: temperature(0:2) = &
        [1500.0_dp, 2400.0_dp, 1800.0_dp]*1.602176634e-12_dp
    ! Integral_0^infinity exp(-u^2)*u^(4+2p) du / sqrt(pi), p=0,1,2.
    ! Together with the independent 2/sqrt(pi) Maxwellian velocity measure
    ! below, these moments check the absolute factor two in the observer.
    real(dp), parameter :: gaussian(3) = &
        [3.0_dp/8.0_dp, 15.0_dp/16.0_dp, 105.0_dp/32.0_dp]
    real(dp) :: qflux(3,3,0:2,0:2), torque(0:2), expected(0:2)
    real(dp) :: source_weight(3,0:2), physical_weight(3,0:2)
    real(dp) :: dn_ds(0:2), dT_ds(0:2), er, epar, dl_over_b, max_error
    integer :: drive, harmonic, m_phi, orientation

    max_error = 0.0_dp
    ! hphi=kappa*Bhat is the Boozer-consistent choice for
    ! Bhat=1+eps*cos(theta), whose uniform-angle mean is exactly one.
    dl_over_b = 2.0_dp*pi/(kappa*(1.0_dp - eps**2)**1.5_dp)
    call check_native_orientation()
    do orientation = -1, 1, 2
        scale = orientation*0.003_dp
        call check_original_operation_order()
        do drive = 1, 3
            call set_physical_drive(drive)
            do harmonic = 1, 2
                m_phi = 2
                if (harmonic == 2) m_phi = -3
                call build_native_response(m_phi)
                call compute_electromagnetic_work(qflux, mass, charge, density, &
                    temperature, dn_ds, dT_ds, er, epar, bref, g, dl_over_b, &
                    scale, torque, source_weight)
                call four_phase_physical_torque(m_phi, expected)
                call check_close('physical source amplitudes', &
                    reshape(source_weight, [9]), reshape(physical_weight, [9]))
                call check_close('four-phase CGS torque', torque, expected)
                print '(a,i2,a,i1,a,i3,a,3es16.7)', 'scale sign ', orientation, &
                    ', drive ', drive, ', n=', m_phi, ', torque [erg/cm^3]:', torque
            end do
        end do
    end do
    call build_native_response(0)
    call compute_electromagnetic_work(qflux, mass, charge, density, &
        temperature, dn_ds, dT_ds, er, epar, bref, g, dl_over_b, scale, &
        torque, source_weight)
    if (any(abs(torque) > 0.0_dp)) error stop 'Axisymmetric work is nonzero'
    print '(a,es12.4)', 'Maximum relative oracle error: ', max_error
    print *, 'All tests passed!'

contains

    subroutine check_native_orientation()
        real(dp), parameter :: p = 0.62_dp**2*5.3_dp/2.0_dp
        real(dp), parameter :: cboozer = 5.3_dp*6.2_dp &
            + 0.5_dp*(5.3_dp*0.62_dp**2*0.09_dp/(2.0_dp*6.2_dp))
        real(dp), parameter :: b = 5.141_dp, iota = 0.5_dp
        real(dp) :: jac, hphi, psi_hat, bphi_hat
        integer :: lab

        ! Independent Boozer identities in native (s,phi,theta), SI to CGS.
        ! P=abs(Phi_edge)/(2*pi), C=G+iota*I, hphi=B/C.
        jac = -cboozer*p/b**2*1.0e6_dp
        hphi = b/cboozer*0.01_dp
        psi_hat = p/5.3_dp*1.0e4_dp
        bphi_hat = 620.0_dp
        do lab = 8, 9
            if (.not. work_geometry_supported(1, lab, jac, hphi, &
                bphi_hat, iota, psi_hat, g, dl_over_b)) &
                error stop 'Native Boozer branch rejected'
            jac = -jac
            if (work_geometry_supported(1, lab, jac, hphi, bphi_hat, &
                iota, psi_hat, g, dl_over_b)) &
                error stop 'Opposite Jacobian was admitted'
            jac = -jac
            hphi = -hphi
            if (work_geometry_supported(1, lab, jac, hphi, bphi_hat, &
                iota, psi_hat, g, dl_over_b)) &
                error stop 'Reversed field line was admitted'
            hphi = -hphi
            bphi_hat = -bphi_hat
            if (work_geometry_supported(1, lab, jac, hphi, bphi_hat, &
                iota, psi_hat, g, dl_over_b)) &
                error stop 'Reversed toroidal field was admitted'
            bphi_hat = -bphi_hat
            psi_hat = -psi_hat
            if (work_geometry_supported(1, lab, jac, hphi, bphi_hat, &
                iota, psi_hat, g, dl_over_b)) &
                error stop 'Signed lab10 flux was admitted'
            psi_hat = -psi_hat
        end do
        if (work_geometry_supported(1, 10, jac, hphi, bphi_hat, iota, &
                psi_hat, g, dl_over_b)) error stop 'Lab10 was admitted'
        if (work_geometry_supported(0, 8, jac, hphi, bphi_hat, iota, &
                psi_hat, g, dl_over_b)) error stop 'Cylindrical case admitted'
        if (work_geometry_supported(1, 8, jac, hphi, bphi_hat, 0.0_dp, &
                psi_hat, g, dl_over_b)) error stop 'Zero iota was admitted'
        if (.not. work_geometry_supported(1, 8, jac, hphi, bphi_hat, -iota, &
                psi_hat, g, dl_over_b)) error stop 'Negative iota was rejected'
    end subroutine check_native_orientation

    subroutine set_physical_drive(which)
        integer, intent(in) :: which
        real(dp) :: gyro_length, q, log_density, log_temperature
        integer :: b

        dn_ds = density*[-0.17_dp, 0.23_dp, 0.41_dp]
        dT_ds = temperature*[0.13_dp, -0.21_dp, 0.31_dp]
        er = 4.0e-4_dp
        epar = 2.1e-5_dp
        if (which == 2) then
            dn_ds = 0.0_dp
            dT_ds = 0.0_dp
            er = 0.0_dp
        else if (which == 3) then
            epar = 0.0_dp
        end if
        do b = 0, 2
            q = charge(b)*elementary_charge
            gyro_length = c*sqrt(2.0_dp*mass(b)*temperature(b))/(q*bref)
            ! Differentiate the physical Maxwellian at fixed energy:
            ! d_s log(F_M) = d_s log(n) + (u^2-3/2)*d_s log(T).
            ! The radial electric term is a physical derivative per cm;
            ! division by g converts it to the same s derivative.
            log_density = dn_ds(b)/density(b)
            log_temperature = dT_ds(b)/temperature(b)
            physical_weight(1,b) = gyro_length*(log_density &
                - 1.5_dp*log_temperature - q*er/(g*temperature(b)))
            physical_weight(2,b) = q*bref*epar/temperature(b)
            physical_weight(3,b) = gyro_length*log_temperature
        end do
        if (which == 2) then
            if (physical_weight(2,0)*physical_weight(2,2) >= 0.0_dp) &
                error stop 'Inductive source lost the electron charge sign'
        end if
    end subroutine set_physical_drive

    subroutine build_native_response(n)
        integer, intent(in) :: n
        integer, parameter :: ntheta = 64
        integer :: j, a, b, channel
        real(dp) :: theta, bhat, hphi, dtheta, pol0, pol2
        complex(dp) :: back, pitch, amplitude, bn

        ! Non-radial rows are deliberate sentinels: the observer must only
        ! contract row one and must preserve the two distinct species axes.
        qflux = 1.0e12_dp
        qflux(1,:,:,:) = 0.0_dp
        dtheta = 2.0_dp*pi/ntheta
        do j = 0, ntheta - 1
            theta = j*dtheta
            bhat = 1.0_dp + eps*cos(theta)
            hphi = kappa*bhat
            bn = field_amplitude(theta)
            back = work_back_geometry(scale, n, bn, bhat)
            do b = 0, 2
                do a = 0, 2
                    do channel = 1, 3
                        ! Integral of p(lambda) and lambda^2*p(lambda).
                        ! The odd term cancels on the two velocity signs.
                        pol0 = 2.0_dp + 0.2_dp*channel
                        pol2 = 2.0_dp/3.0_dp + 0.12_dp*channel
                        pitch = work_pitch_convolution(pol0, pol2, bhat, hphi)
                        amplitude = response(a, b, channel, theta)
                        qflux(1,channel,a,b) = qflux(1,channel,a,b) &
                            + 0.5_dp*real(conjg(back*pitch)*amplitude, dp) &
                            *gaussian(channel)*dtheta
                    end do
                end do
            end do
        end do
    end subroutine build_native_response

    subroutine four_phase_physical_torque(n, result)
        integer, intent(in) :: n
        real(dp), intent(out) :: result(0:2)
        integer, parameter :: ntheta = 64
        real(dp), parameter :: lambda(3) = &
            [-sqrt(0.6_dp), 0.0_dp, sqrt(0.6_dp)]
        real(dp), parameter :: weight(3) = &
            [5.0_dp/9.0_dp, 8.0_dp/9.0_dp, 5.0_dp/9.0_dp]
        real(dp) :: theta, bhat, chi, radial_measure, pitch_shape
        real(dp) :: physical_distribution, hamiltonian_derivative, contribution
        complex(dp) :: phase, bn, amplitude
        integer :: j, l, phase_index, a, b, channel

        result = 0.0_dp
        do j = 0, ntheta - 1
            theta = j*2.0_dp*pi/ntheta
            bhat = 1.0_dp + eps*cos(theta)
            radial_measure = 2.0_dp*pi/(ntheta*kappa*bhat**2*dl_over_b)
            bn = field_amplitude(theta)
            do a = 0, 2
                do l = 1, 3
                    do phase_index = 0, 3
                        chi = phase_index*pi/2.0_dp
                        phase = cmplx(cos(chi), sin(chi), dp)
                        ! H_1/T_a = u^2*(1+lambda^2)*Re(bn*exp(i*chi)).
                        ! Gaussian integration supplies u^2 from H_1 and
                        ! u^2 from d^3v; the remaining speed powers belong
                        ! to each manufactured SOL response channel.
                        hamiltonian_derivative = temperature(a) &
                            *(1.0_dp + lambda(l)**2) &
                            *real(cmplx(0.0_dp, real(n, dp), dp)*bn*phase, dp)
                        contribution = 0.0_dp
                        do b = 0, 2
                            do channel = 1, 3
                                pitch_shape = 1.0_dp + 0.3_dp*channel*lambda(l)**2 &
                                    + (0.07_dp*b - 0.03_dp*a)*lambda(l)
                                amplitude = response(a, b, channel, theta)
                                physical_distribution = density(a) &
                                    *physical_weight(channel,b)*pitch_shape &
                                    *real(amplitude*phase, dp)
                                contribution = contribution - 2.0_dp &
                                    *gaussian(channel)*physical_distribution &
                                    *hamiltonian_derivative
                            end do
                        end do
                        result(a) = result(a) &
                            + radial_measure*weight(l)*contribution/4.0_dp
                    end do
                end do
            end do
        end do
    end subroutine four_phase_physical_torque

    pure function field_amplitude(theta) result(bn)
        real(dp), intent(in) :: theta
        complex(dp) :: bn
        bn = cmplx(0.017_dp*(1.0_dp + 0.1_dp*cos(theta)), &
            -0.009_dp*(1.0_dp + 0.2_dp*sin(theta)), dp)
    end function field_amplitude

    pure function response(a, b, channel, theta) result(amplitude)
        integer, intent(in) :: a, b, channel
        real(dp), intent(in) :: theta
        complex(dp) :: amplitude
        amplitude = cmplx(0.013_dp*(1 + 2*a - b + channel), &
            0.019_dp*(2 - a + 3*b - channel), dp) &
            *(1.0_dp + 0.15_dp*sin(theta))
    end function response

    subroutine check_original_operation_order()
        complex(dp) :: bn, coefficient, original, actual
        real(dp) :: bhat, hphi, pol0, pol2
        integer :: j

        do j = 1, 17
            bhat = 0.61_dp + 0.071_dp*j
            hphi = -0.13_dp + 0.0031_dp*j
            bn = cmplx(0.017_dp*j, -0.031_dp/(j + 1), dp)
            pol0 = 1.0_dp/j
            pol2 = -0.017_dp*j
            original = scale*cmplx(0.0_dp, 1.0_dp, dp)*(-3)*bn*bhat
            actual = work_back_geometry(scale, -3, bn, bhat)
            call check_bits('original back geometry', actual, original)
            coefficient = cmplx(1.0_dp, 0.0_dp, dp)/bhat**2/hphi
            original = (pol0 + pol2)*coefficient
            actual = work_pitch_convolution(pol0, pol2, bhat, hphi)
            call check_bits('original pitch convolution', actual, original)
        end do
    end subroutine check_original_operation_order

    subroutine check_bits(label, actual, expected_value)
        character(len=*), intent(in) :: label
        complex(dp), intent(in) :: actual, expected_value
        if (any(transfer(actual, [0_int64, 0_int64]) &
                /= transfer(expected_value, [0_int64, 0_int64]))) then
            print *, 'FAIL: ', label, actual, expected_value
            error stop 'Default arithmetic changed'
        end if
    end subroutine check_bits

    subroutine check_close(label, actual, expected_value)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: actual(:), expected_value(:)
        real(dp) :: difference

        difference = maxval(abs(actual - expected_value)) &
            /max(1.0e-25_dp, maxval(abs(expected_value)))
        max_error = max(max_error, difference)
        if (difference > 3.0e-12_dp) then
            print *, 'FAIL: ', label, difference
            error stop 'Physical work oracle mismatch'
        end if
    end subroutine check_close

end program test_electromagnetic_work
