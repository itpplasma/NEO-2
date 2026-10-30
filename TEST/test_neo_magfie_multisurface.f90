! Real native Strumberger.bc bootstrap. Separate-process single surfaces are
! the independent lifecycle oracle for simultaneous initialization.
program test_neo_magfie_multisurface
    use nrtype, only: dp
    use neo_magfie, only: neo_magfie_calc, magfie_spline, magfie_result, &
        magfie_sarray
    use partpa_mod, only: bmod0
    use neo_exchange, only: b_min, b_max
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    character(32) :: mode
    integer :: count, first, j, k, fo, u1, u2, um, ios
    real(dp) :: expected(17), actual(17), relerr, sample_min, sample_max
    character(16) :: filename
    real(dp) :: surfaces(2), x(3), B, Jlegacy, bd(3), hc(3), ht(3), curl(3)
    real(dp), parameter :: pi = acos(-1.0_dp)
    call get_command_argument(1, mode)
    if (trim(mode) == 'compare') then
        open (newunit=u1, file='one_first.dat', status='old')
        open (newunit=u2, file='one_second.dat', status='old')
        open (newunit=um, file='both.dat', status='old')
        do j = 1, 2
            do k = 0, 32
                if (j == 1) read (u1, *) expected
                if (j == 2) read (u2, *) expected
                read (um, *) actual
                if (.not. all(ieee_is_finite([actual, expected]))) then
                    error stop 'FAIL nonfinite'
                end if
                relerr = maxval(abs(actual - expected)/max(1.0_dp, abs(expected)))
                if (relerr > 1.0e-12_dp) error stop 'FAIL single/multi disagreement'
            end do
        end do
        read (um, *, iostat=ios) actual
        if (ios >= 0) error stop 'FAIL unexpected extra multisurface rows'
        close (u1); close (u2); close (um)
        print *, 'All tests passed!'
        stop
    end if
    surfaces = [0.00475_dp, 0.00525_dp]
    select case (trim(mode))
    case ('one_first')
        count = 1; first = 1
    case ('one_second')
        count = 1; first = 2
    case ('both')
        count = 2; first = 1
    case default
        error stop 'EXPECTED one_first/one_second/both'
    end select
    allocate (magfie_sarray(count))
    magfie_sarray = surfaces(first:first + count - 1)
    magfie_spline = 1; magfie_result = 0; bmod0 = 1.0e-4_dp
    filename = trim(mode)//'.dat'
    open (newunit=fo, file=trim(filename), status='replace')
    do j = first, first + count - 1
        do k = 0, 32
            x = [surfaces(j), 0.0_dp, 2*pi*k/33]
            call neo_magfie_calc(x, B, Jlegacy, bd, hc, ht, curl)
            if (.not. all(ieee_is_finite([B, Jlegacy, bd, hc, ht, curl]))) &
                error stop 'FAIL nonfinite native field'
            write (fo, '(17(1x,es25.16e3))') x, B, Jlegacy, &
                bd, hc, ht, curl
        end do
    end do
    close (fo)
    if (count == 1) then
        ! Existing axisymmetric fixture: sample every angular initialization node.
        sample_min = huge(1.0_dp); sample_max = -huge(1.0_dp)
        do k = 0, 128
            x = [surfaces(first), 0.0_dp, 2*pi*k/128]
            call neo_magfie_calc(x, B, Jlegacy, bd, hc, ht, curl)
            sample_min = min(sample_min, B); sample_max = max(sample_max, B)
        end do
        if (abs(b_min/sample_min - 1.0_dp) > 1.0e-12_dp) &
            error stop 'FAIL single-surface minimum'
        if (abs(b_max/sample_max - 1.0_dp) > 1.0e-12_dp) &
            error stop 'FAIL single-surface maximum'
    end if
end program test_neo_magfie_multisurface
