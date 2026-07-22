program probe_neo2_perturbation_fourier_kernel
    use nrtype, only: dp
    use neo_control, only: inp_swi, lab_swi
    use neo_magfie_perturbation, only: neo_magfie_pert
    use neo_spline_data, only: lsw_linear_boozer
    use ntv_mod, only: in_file_pert
    use partpa_mod, only: bmod0

    implicit none

    character(len=*), parameter :: fixture = 'manufactured_neo2_perturbation.bc'
    real(dp), parameter :: theta = 0.6_dp, phi = -0.4_dp
    real(dp), parameter :: cosine_amplitude = 2.0e-4_dp
    real(dp), parameter :: sine_amplitude = 0.7e-4_dp
    complex(dp) :: actual
    real(dp) :: expected_real, expected_imag, wrong_sign_real, x(3)

    call write_fixture(fixture)
    inp_swi = 9
    lab_swi = 9
    in_file_pert = fixture
    bmod0 = 1.0_dp
    lsw_linear_boozer = .true.

    x = [0.5_dp, phi, theta]
    call neo_magfie_pert(x, actual)
    expected_real = real(expected_value(theta + 3.0_dp*phi))
    expected_imag = aimag(expected_value(theta + 3.0_dp*phi))
    wrong_sign_real = real(expected_value(theta - 3.0_dp*phi))

    write (*, '(*(g0,1x))') 'NEO2_PERT_KERNEL', 'theta', theta, &
        'phi_ccw', phi, 'actual_real', real(actual), 'actual_imag', aimag(actual), &
        'expected_plus_n_real', expected_real, &
        'expected_plus_n_imag', expected_imag, 'wrong_minus_n_real', wrong_sign_real

contains

    complex(dp) function expected_value(phase)
        real(dp), intent(in) :: phase

        expected_value = cmplx(cosine_amplitude, -sine_amplitude, kind=dp) &
            * exp(cmplx(0.0_dp, phase, kind=dp))
    end function expected_value

    subroutine write_fixture(path)
        character(len=*), intent(in) :: path

        integer :: unit, surface
        real(dp), parameter :: surfaces(3) = [0.2_dp, 0.5_dp, 0.8_dp]

        open (newunit=unit, file=path, status='replace', action='write')
        write (unit, '(a)') 'CC manufactured perturbation kernel'
        write (unit, '(a)') 'CC physical phi is CCW'
        write (unit, '(a)') 'CC kernel m*theta+n*phi'
        write (unit, '(a)') 'CC complex coefficient bmnc-i*bmns'
        write (unit, '(a)') 'm0b n0b nsurf nper flux a R'
        write (unit, '(*(g0,1x))') 2, 0, 3, 1, 1.0_dp, 1.0_dp, 1.0_dp
        do surface = 1, size(surfaces)
            call write_surface(unit, surfaces(surface))
        end do
        close (unit)
    end subroutine write_fixture

    subroutine write_surface(unit, radial_coordinate)
        integer, intent(in) :: unit
        real(dp), intent(in) :: radial_coordinate

        write (unit, '(a)') 's iota Jpol Itor pprime sqrtg'
        write (unit, '(a)') 'units A A Pa m3'
        write (unit, '(*(es24.16e3,1x))') radial_coordinate, 1.0_dp, &
            0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp
        write (unit, '(a)') 'm n rmnc rmns zmnc zmns vmnc vmns bmnc bmns'
        call write_mode(unit, -1, 0.0_dp, 0.0_dp)
        call write_mode(unit, 0, 0.0_dp, 0.0_dp)
        call write_mode(unit, 1, cosine_amplitude, sine_amplitude)
    end subroutine write_surface

    subroutine write_mode(unit, poloidal_mode, bmnc, bmns)
        integer, intent(in) :: unit, poloidal_mode
        real(dp), intent(in) :: bmnc, bmns

        write (unit, '(*(g0,1x))') poloidal_mode, 3, 0.0_dp, 0.0_dp, &
            0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, bmnc, bmns
    end subroutine write_mode

end program probe_neo2_perturbation_fourier_kernel
