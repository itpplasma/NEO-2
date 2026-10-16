program test_field_provider
    ! Independent large-aspect-ratio circular oracle:
    ! R = R0 + a*sqrt(s)*cos(theta), Z = sigma_theta*a*sqrt(s)*sin(theta).
    ! The metric gives B_phi = sigma_phi*B0*R0,
    ! B_theta = r**2*iota*sigma_phi*B0/R0, and
    ! sqrt(g) = handedness*R0*a**2/2.
    use iso_fortran_env, only: dp => real64
    use neo_field_provider, only: neo_field_load, neo_field_prepare, neo_field_free
    use magfie_mod, only: magfie
    use neo_magfie, only: Bphcov => boozer_curr_pol_hat, Bthcov => boozer_curr_tor_hat
    use neo_magfie_perturbation, only: neo_load_pert, neo_free_pert, &
        neo_magfie_pert, neo_magfie_pert_amplitude
    use partpa_mod, only: bmod0
    use nrtype, only: pi

    implicit none

    character(32) :: mode_arg
    real(dp) :: field_t = 2.0_dp, mirror_fraction = 0.0_dp
    integer :: orientation, iota_sign
    real(dp), parameter :: major_radius_m = 5.0_dp
    real(dp), parameter :: minor_radius_m = 0.5_dp
    real(dp), parameter :: iota_abs = 0.2_dp
    real(dp), parameter :: evaluation_s = 0.5_dp
    real(dp), parameter :: current_to_covar = 0.2_dp
    real(dp), parameter :: meter_to_cm = 100.0_dp
    real(dp), parameter :: tesla_to_gauss = 1.0e4_dp
    real(dp), parameter :: tolerance = 1.0e-12_dp

    call get_command_argument(1, mode_arg)
    if (mode_arg == 'bad-format') call neo_field_load('unused', 0)
    call write_chart("manufactured_left.bc", 1, -1, -1)
    call write_chart("manufactured_left_positive_iota.bc", 1, 1, -1)
    call write_chart("manufactured_right.bc", -1, 1, 1)
    if (mode_arg == 'outside') then
        call neo_field_load('manufactured_left.bc', 9)
        call neo_field_prepare(2.0_dp)
        error stop 'outside surface accepted'
    end if
    call check_chart("manufactured_left.bc", 1, -1, -1)
    call check_chart("manufactured_left_positive_iota.bc", 1, 1, -1)
    call check_chart("manufactured_right.bc", -1, 1, 1)
    field_t = 3.0_dp
    call write_chart("manufactured_second_field.bc", 1, 1, -1)
    call check_chart("manufactured_second_field.bc", 1, 1, -1)
    field_t = 2.0_dp
    call check_chart("manufactured_left.bc", 1, -1, -1)
    call check_perturbation(8, 1.0_dp)
    call check_perturbation(9, 2.0_dp)
    call check_perturbation(8, 3.0_dp)
    mirror_fraction = 0.15_dp
    do orientation = -1, 1, 2
        do iota_sign = -1, 1, 2
            call write_chart('manufactured_curl.bc', orientation, iota_sign, -orientation)
            call check_curl_orientation('manufactured_curl.bc')
        end do
    end do
    call check_curl_orientation('manufactured_curl.bc', 8)
    call neo_free_pert()
    call neo_free_pert()
    call neo_field_free()
    call neo_field_free()

contains

    subroutine check_curl_orientation(path, laboratory)
        use neo_control, only: lab_swi
        character(*), intent(in) :: path
        integer, intent(in), optional :: laboratory
        real(dp), parameter :: angles(5) = [0.0_dp, 0.7_dp, 1.9_dp, 3.4_dp, 5.1_dp]
        real(dp), parameter :: h = 1.0e-3_dp
        real(dp) :: x(3), jacobian, hc(3), curl(3)
        real(dp) :: plus(3), minus(3), plus_half(3), minus_half(3)
        real(dp) :: derivative, expected, measured, saved_curl(3), saved_jacobian
        integer :: angle, direction
        call neo_field_load(path, 9)
        if (present(laboratory)) lab_swi = laboratory
        bmod0 = 1.0e-4_dp
        do angle = 1, size(angles)
            x = [evaluation_s, 0.0_dp, angles(angle)]
            call curl_evaluate(x, hc, jacobian, curl)
            saved_curl = curl
            saved_jacobian = jacobian
            do direction = 1, 3, 2
                x(direction) = x(direction) + h
                call curl_evaluate(x, plus, jacobian, curl)
                x(direction) = x(direction) - 2*h
                call curl_evaluate(x, minus, jacobian, curl)
                x(direction) = x(direction) + 1.5_dp*h
                call curl_evaluate(x, plus_half, jacobian, curl)
                x(direction) = x(direction) - h
                call curl_evaluate(x, minus_half, jacobian, curl)
                x(direction) = x(direction) + 0.5_dp*h
                derivative = (4*(plus_half(2)-minus_half(2))/h &
                    - (plus(2)-minus(2))/(2*h))/3
                expected = derivative/saved_jacobian
                measured = saved_curl(3)
                if (direction == 3) then
                    expected = -expected
                    measured = saved_curl(1)
                end if
                if (abs(measured-expected) > 1.0e-8_dp*max(abs(expected),1.0e-8_dp)) &
                    error stop 'curl orientation or independent covariant derivative'
            end do
        end do
    end subroutine check_curl_orientation

    subroutine curl_evaluate(position, covariant, jacobian, curl)
        real(dp), intent(in) :: position(3)
        real(dp), intent(out) :: covariant(3), jacobian, curl(3)
        real(dp) :: bmod, bd(3), ht(3)
        call neo_field_prepare(position(1))
        call magfie(position, bmod, jacobian, bd, covariant, ht, curl)
    end subroutine curl_evaluate

    subroutine check_perturbation(input_switch, scale)
        integer, intent(in) :: input_switch
        real(dp), intent(in) :: scale
        real(dp), parameter :: radial_grid(4) = [0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp]
        real(dp), parameter :: cosine(3) = [0.1_dp, 0.2_dp, 0.3_dp]
        real(dp), parameter :: sine(3) = [0.4_dp, -0.2_dp, 0.15_dp]
        complex(dp), parameter :: imaginary = (0.0_dp, 1.0_dp)
        real(dp) :: x(3), radial_factor, c, d
        complex(dp) :: actual, expected, full_field
        integer :: unit, j, k, m, n
        open(newunit=unit, file="provider_pert.bc", status="replace")
        do j = 1, 5
            write(unit, '(a)') 'CC manufactured signed-m perturbation'
        end do
        write(unit, *) 2, 0, 4, 1, 1.0_dp, 0.5_dp, 5.0_dp
        do j = 1, 4
            write(unit, '(a)') 'CC surface'
            write(unit, '(a)') 'CC s iota currents'
            write(unit, *) radial_grid(j), 0.2_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp
            write(unit, '(a)') 'CC m n R Z v B'
            do k = 1, 3
                m = k - 2
                n = 2
                if (mode_arg == 'multi-n') n = k
                radial_factor = radial_grid(j)**(0.5_dp*abs(m))
                c = scale*cosine(k)*radial_factor
                d = scale*sine(k)*radial_factor
                if (input_switch == 8) then
                    write(unit, *) m, n, 0.0_dp, 0.0_dp, 0.0_dp, c
                else
                    write(unit, *) m, n, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                        0.0_dp, 0.0_dp, c, d
                end if
            end do
        end do
        close(unit)
        call neo_load_pert('provider_pert.bc', input_switch)
        x = [0.5_dp, 0.37_dp, 0.63_dp]
        expected = (0.0_dp, 0.0_dp)
        do k = 1, 3
            m = k - 2
            radial_factor = x(1)**(0.5_dp*abs(m))
            c = scale*cosine(k)*radial_factor
            d = scale*sine(k)*radial_factor
            if (input_switch == 8) then
                expected = expected + c*exp(-imaginary*m*x(3))
            else
                expected = expected + (c - imaginary*d)*exp(imaginary*m*x(3))
            end if
        end do
        expected = expected/bmod0
        call neo_magfie_pert_amplitude(x, actual)
        if (abs(actual - expected) > 1.0e-10_dp*max(1.0_dp, abs(expected))) &
            error stop 'signed-m amplitude or unit normalization'
        call neo_magfie_pert(x, full_field)
        if (abs(full_field - expected*exp(imaginary*2*x(2))) > &
            1.0e-10_dp*max(1.0_dp, abs(expected))) error stop 'toroidal phase'
    end subroutine check_perturbation

    subroutine check_chart(path, orientation, iota_sign, handedness)
        character(len=*), intent(in) :: path
        integer, intent(in) :: orientation, iota_sign, handedness

        real(dp) :: bder(3), bmod, hcovar(3), hctrvr(3), hcurl(3)
        real(dp) :: expected_Bphcov, expected_Bthcov, expected_sqrtg
        real(dp) :: iota_chart, radius_cm, sigma_phi, sqrtg, x(3)

        sigma_phi = -real(handedness, dp)/real(orientation, dp)
        iota_chart = real(iota_sign, dp)*iota_abs
        radius_cm = minor_radius_m*meter_to_cm*sqrt(evaluation_s)
        expected_Bphcov = sigma_phi*field_t*tesla_to_gauss*major_radius_m*meter_to_cm
        expected_Bthcov = field_t*tesla_to_gauss*radius_cm**2*iota_chart*sigma_phi &
            /(major_radius_m*meter_to_cm)
        expected_sqrtg = real(handedness, dp)*major_radius_m*meter_to_cm &
            *(minor_radius_m*meter_to_cm)**2/2.0_dp

        call neo_field_load(path, 9)
        bmod0 = 1.0e-4_dp
        call neo_field_prepare(evaluation_s)
        x = [evaluation_s, 0.0_dp, 0.0_dp]
        call magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)

        call assert_close(path//" B_theta", Bthcov, expected_Bthcov)
        call assert_close(path//" B_phi", Bphcov, expected_Bphcov)
        call assert_close(path//" h_theta", hcovar(3)*bmod*1.0e4_dp, expected_Bthcov)
        call assert_close(path//" h_phi", hcovar(2)*bmod*1.0e4_dp, expected_Bphcov)
        call assert_close(path//" sqrt(g)", sqrtg, expected_sqrtg)
    end subroutine check_chart

    subroutine write_chart(path, orientation, iota_sign, handedness)
        character(len=*), intent(in) :: path
        integer, intent(in) :: orientation, iota_sign, handedness

        integer :: file_unit, surface_index
        real(dp) :: Bph_si, Bth_si, flux, iota_chart, radius_m
        real(dp) :: raw_Jpol, raw_Itor, s_surface, sigma_phi
        real(dp), parameter :: surfaces(3) = [0.25_dp, 0.5_dp, 0.75_dp]

        sigma_phi = -real(handedness, dp)/real(orientation, dp)
        iota_chart = real(iota_sign, dp)*iota_abs
        flux = -real(handedness, dp)*sigma_phi*pi*minor_radius_m**2*field_t
        Bph_si = sigma_phi*field_t*major_radius_m

        open (newunit=file_unit, file=path, status="replace", action="write")
        call write_header(file_unit, flux)
        do surface_index = 1, size(surfaces)
            s_surface = surfaces(surface_index)
            radius_m = minor_radius_m*sqrt(s_surface)
            Bth_si = field_t*radius_m**2*iota_chart*sigma_phi/major_radius_m
            raw_Jpol = Bph_si*tesla_to_gauss*meter_to_cm &
                /(current_to_covar*real(handedness, dp))
            raw_Itor = Bth_si*tesla_to_gauss*meter_to_cm &
                /(current_to_covar*real(handedness, dp))
            call write_surface(file_unit, s_surface, iota_chart, raw_Jpol, &
                raw_Itor, radius_m, orientation)
        end do
        close (file_unit)
    end subroutine write_chart

    subroutine write_header(file_unit, flux)
        integer, intent(in) :: file_unit
        real(dp), intent(in) :: flux

        write (file_unit, "(a)") "CC manufactured large-aspect-ratio circular chart"
        write (file_unit, "(a)") "CC independent analytic field and geometry"
        write (file_unit, "(a)") "CC m0b n0b nsurf nper flux a R"
        write (file_unit, "(a)") "CC"
        write (file_unit, "(a)") "CC"
        write (file_unit, *) 1, 0, 3, 1, flux, minor_radius_m, major_radius_m
    end subroutine write_header

    subroutine write_surface(file_unit, s_surface, iota_chart, Jpol, Itor, &
            radius_m, orientation)
        integer, intent(in) :: file_unit
        integer, intent(in) :: orientation
        real(dp), intent(in) :: s_surface, iota_chart, Jpol, Itor, radius_m

        real(dp) :: field_magnitude

        field_magnitude = field_t*sqrt(1.0_dp + (radius_m*iota_abs/major_radius_m)**2)
        write (file_unit, "(a)") "CC s iota Jpol Itor pprime sqrtg00"
        write (file_unit, "(a)") "CC units A A Pa m3"
        write (file_unit, *) s_surface, iota_chart, Jpol, Itor, 0.0_dp, 0.0_dp
        write (file_unit, "(a)") "CC m n rmnc rmns zmnc zmns vmnc vmns bmnc bmns"
        write (file_unit, *) 0, 0, major_radius_m, 0.0_dp, 0.0_dp, 0.0_dp, &
            0.0_dp, 0.0_dp, field_magnitude, 0.0_dp
        write (file_unit, *) 1, 0, radius_m, 0.0_dp, 0.0_dp, &
            real(orientation, dp)*radius_m, &
            0.0_dp, 0.0_dp, -mirror_fraction*field_magnitude*sqrt(s_surface), 0.0_dp
    end subroutine write_surface

    subroutine assert_close(label, actual, expected)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: actual, expected

        real(dp) :: relative_error

        relative_error = abs(actual - expected)/max(abs(expected), tiny(expected))
        if (relative_error > tolerance) then
            print *, trim(label), " expected", expected, "got", actual
            error stop "manufactured Boozer oracle failed"
        end if
    end subroutine assert_close

end program test_field_provider
