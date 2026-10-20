program test_perturbation_repeated_surfaces
    use nrtype, only: dp
    use neo_control, only: lab_swi
    use partpa_mod, only: bmod0
    use neo_magfie_perturbation, only: neo_load_pert, neo_magfie_pert_amplitude
    implicit none
    real(dp), parameter :: surfaces(5) = [0.2_dp, 0.5_dp, 0.8_dp, 0.8_dp, 0.8_dp]
    complex(dp), parameter :: coefficient = (1.0e-4_dp, -0.3e-4_dp)
    complex(dp) :: amplitude, baseline(2), expected
    real(dp) :: x(3), bmnc
    integer :: unit, j, k, count, variant
    character(16) :: mode
    character(64) :: path

    call get_command_argument(1, mode)
    path = "repeated_perturbation"//trim(mode)//".bc"
    lab_swi = 10
    bmod0 = 1.0e-4_dp
    do variant = 1, 2
        count = 3
        if (variant == 2) count = 5
        open(newunit=unit, file=path, status='replace')
        do j = 1, 5
            write(unit, '(a)') 'CC manufactured repeated perturbation surfaces'
        end do
        write(unit, *) 0, 0, count, 1, 0.0_dp, 0.0_dp, 0.0_dp
        do j = 1, count
            write(unit, '(a)') 'CC s iota Jpol Itor pprime sqrtg'
            write(unit, '(a)') 'CC units'
            write(unit, *) surfaces(j), 0.2_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp
            write(unit, '(a)') 'CC m n Rcos Rsin Zcos Zsin vcos vsin Bcos Bsin'
            bmnc = real(coefficient, dp)
            if (mode == 'conflict' .and. j == 5) bmnc = 2.0_dp*bmnc
            write(unit, *) 0, 3, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                0.0_dp, 0.0_dp, bmnc, -aimag(coefficient)
        end do
        close(unit)
        call neo_load_pert(trim(path), 9)
        expected = coefficient/bmod0
        do k = 1, 2
            x = [0.5_dp + 0.2_dp*(k-1), 0.4_dp, 0.6_dp]
            call neo_magfie_pert_amplitude(x, amplitude)
            if (abs(amplitude-expected) > 1.0e-12_dp) error stop "manufactured amplitude"
            if (variant == 1) then
                baseline(k) = amplitude
            else
                if (abs(amplitude-baseline(k)) > 1.0e-12_dp) error stop "repeated knot changes field"
            end if
        end do
    end do
    open(newunit=unit, file=path, status='old')
    close(unit, status='delete')
    print *, "repeated perturbation surfaces: PASS"
end program test_perturbation_repeated_surfaces
