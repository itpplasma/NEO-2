module native_moment_fixture
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use collision_moment_map
    use moment_schur_solver
    use native_integral_oracle
    use moment_test_support
    implicit none
    type(collision_moment_map_t) :: map, invalid_map
    type(moment_schur_operator_t) :: op
    type(moment_schur_options_t) :: options
    type(moment_schur_report_t) :: report
    complex(dp), allocatable :: qfull(:), ffull(:), truth(:)
    complex(dp) :: no_constraints(0)
    complex(dp), allocatable :: fields(:, :), moments(:, :), mixed(:, :), expected(:), actual(:), zerof(:)
    real(dp) :: error, worst, before, changed
    real(dp), allocatable :: saved_coeff(:, :, :, :, :)
    integer :: v, l, p, phi, a, b, j, k, info, ib, trial, trials

contains
    subroutine run_native_moment_tests()
        lag = 2; leg = 3; ibeg = -2; iend = 2; npart = 4; num_spec = 3
        allocate (npl(ibeg:iend), ind_start(ibeg:iend))
        npl = [0, 3, 1, 4, 2]
        ind_start(ibeg) = 0
        do phi = ibeg + 1, iend
            ind_start(phi) = ind_start(phi - 1) + 2*(npl(phi - 1) + 1)*(lag + 1)
        end do
        n_2d_size = ind_start(iend) + 2*(npl(iend) + 1)*(lag + 1)
        allocate (pleg_bra(0:leg, npart + 1, ibeg:iend), pleg_ket(0:leg, npart + 1, ibeg:iend))
        allocate (ailmm_aa(0:lag, 0:lag, 0:leg, 0:num_spec - 1, 0:num_spec - 1))
        allocate (fact_pos_e(ibeg:iend), fact_neg_e(ibeg:iend))
        allocate (gathered(0:lag, 0:leg, ibeg:iend, 0:num_spec - 1))
        do phi = ibeg, iend
            fact_pos_e(phi) = 0.4_dp + 0.07_dp*(phi - ibeg)
            fact_neg_e(phi) = -0.9_dp + 0.03_dp*(phi - ibeg)
            do p = 1, npart + 1
                do l = 0, leg
                    ! Deliberately nonorthogonal, phi-dependent bras/kets: no VU=I assumption.
                    pleg_bra(l, p, phi) = sin(real(1 + l + 2*p + phi, dp))/3.0_dp
                    pleg_ket(l, p, phi) = cos(real(2 + 2*l + p - phi, dp))/4.0_dp
                end do
            end do
        end do
        do a = 0, num_spec - 1
            do b = 0, num_spec - 1
                do l = 0, leg
                    do v = 0, lag
                        do j = 0, lag
                            ailmm_aa(v, j, l, a, b) = sin(real(1 + v + 3*j + 7*l + 11*a + 13*b, dp))/8.0_dp
                        end do
                    end do
                end do
            end do
        end do
        k = (lag + 1)*(leg + 1)*(iend - ibeg + 1)
        allocate (fields(n_2d_size, num_spec), moments(k, num_spec), mixed(k, num_spec))
        allocate (expected(n_2d_size), actual(n_2d_size), zerof(n_2d_size))
        worst = 0; trials = 0
        do ib = 0, 1
            addboucol = ib == 1
            call init_collision_moment_map(map, lag + 1, npl, pleg_bra, pleg_ket, &
                                          fact_pos_e, fact_neg_e, addboucol, info)
            call check(info == 0, 'valid ragged map')
            call check(map%field_size() == n_2d_size .and. map%moment_size() == k, 'map dimensions')
            do trial = 1, 7
                do a = 1, num_spec
                    do j = 1, n_2d_size
                        fields(j, a) = cmplx(sin(real(j + 17*a + trial, dp)), cos(real(3*j - a + trial, dp)), dp)
                    end do
                end do
                if (trial == 1) fields = 0
                if (trial == 2) fields(:, 2:) = 0 ! drive one species, retain interspecies response
                if (trial == 3) then
                    fields = 0
                    fields(1, 1) = (1.0_dp, -0.5_dp)
                    fields(n_2d_size, 3) = (-2.0_dp, 0.7_dp)
                end if
                do a = 1, num_spec
                    call map%project(fields(:, a), moments(:, a), info)
                    call check(info == 0, 'project')
                end do
                gathered = reshape(moments, shape(gathered))
                call mix_collision_moments(ailmm_aa, moments, mixed, info)
                call check(info == 0, 'mix all species')
                do ispec = 0, num_spec - 1
                    call integral_part(fields(:, ispec + 1), expected)
                    call map%lift(mixed(:, ispec + 1), actual, info)
                    call check(info == 0, 'lift')
                    error = maxval(abs(actual - expected))/max(1.0_dp, maxval(abs(expected)))
                    worst = max(worst, error)
                    call check(error < 2.0e-13_dp, 'native UCV equality including bounce and endpoints')
                    trials = trials + 1
                end do
            end do
        end do
        call check(mpro%calls == trials, 'each emulated species collective was exercised')

        ! Negative controls: the comparisons must detect dropped interspecies coupling and bounce terms.
        call map%lift(mixed(:, 1), expected, info)
        before = maxval(abs(expected))
        call init_collision_moment_map(invalid_map, lag + 1, npl, pleg_bra, pleg_ket, &
                                      fact_pos_e, fact_neg_e, .false., info)
        call invalid_map%lift(mixed(:, 1), actual, info)
        call check(maxval(abs(actual - expected)) > 1.0e-6_dp*before, 'fixture detects missing bounce terms')
        saved_coeff = ailmm_aa
        ailmm_aa(:, :, :, 0, 1:) = 0
        call mix_collision_moments(ailmm_aa, moments, mixed, info)
        call map%lift(mixed(:, 1), actual, info)
        changed = maxval(abs(actual - expected))
        call check(changed > 1.0e-6_dp*before, 'fixture detects missing interspecies coupling')
        ailmm_aa = saved_coeff

        call map%project(fields(1:n_2d_size - 1, 1), moments(:, 1), info)
        call check(info == moment_map_invalid, 'reject wrong field size')
        call map%lift(mixed(1:k - 1, 1), actual, info)
        call check(info == moment_map_invalid, 'reject wrong moment size')
        fields(1, 1) = cmplx(ieee_value(0.0_dp, ieee_quiet_nan), 0.0_dp, dp)
        call map%project(fields(:, 1), moments(:, 1), info)
        call check(info == moment_map_invalid, 'reject nonfinite fields')
        call init_collision_moment_map(invalid_map, 0, npl, pleg_bra, pleg_ket, &
                                      fact_pos_e, fact_neg_e, .true., info)
        call check(info == moment_map_invalid, 'reject invalid velocity count')
        call invalid_map%project(fields(:, 2), moments(:, 2), info)
        call check(info == moment_map_invalid, 'reject uninitialized map')
        zerof = 0
        call map%project(zerof, moments(:, 1), info)
        call map%lift(moments(:, 1), actual, info)
        call check(maxval(abs(actual)) < tiny(1.0_dp), 'exact zero source')
        print '(a,i0,a,es14.6)', 'native comparisons=', trials, ' max_scaled_error=', worst

        ! End-to-end manufactured kinetic problem using the REAL collision kernel.
        ! P is a known nonsingular diagonal surrogate, NOT a production streaming matrix.
        allocate (qfull(n_2d_size*num_spec), ffull(n_2d_size*num_spec), truth(n_2d_size*num_spec))
        do j = 1, size(truth)
            truth(j) = cmplx(sin(real(j, dp)), cos(real(3*j, dp)), dp)
        end do
        op%nfield = size(truth); op%nmoment = k*num_spec
        op%inverse_p => inverse_diagonal
        op%moments => project_all
        op%feedback => lift_all
        op%original => original_all
        call original_all(truth, qfull, info)
        call check(info == 0, 'native manufactured RHS')
        options%rtol = 1.0e-12_dp
        call solve_moment_schur(op, qfull, no_constraints, ffull, no_constraints, options, report)
        call check(report%status == schur_ok, 'end-to-end native collision reduced solve')
        error = maxval(abs(ffull - truth))
        call check(error < 2.0e-11_dp, 'native collision full-distribution recovery')
        print '(a,i0,a,i0,a,i0,a,es14.6,a,es14.6)', &
            'native field dimension=', size(truth), ' moment dimension=', op%nmoment, &
            ' GMRES steps=', report%iterations, ' field_error=', error, ' raw_residual=', report%original_residual
        call finished('QL collision moment map')
    end subroutine run_native_moment_tests
    subroutine project_all(x, y, ierr)
        complex(dp), intent(in) :: x(:)
        complex(dp), intent(out) :: y(:)
        integer, intent(out) :: ierr
        integer :: a
        do a = 1, num_spec
            call map%project(x((a - 1)*n_2d_size + 1:a*n_2d_size), y((a - 1)*k + 1:a*k), ierr)
            if (ierr /= 0) return
        end do
    end subroutine
    subroutine lift_all(x, y, ierr)
        complex(dp), intent(in) :: x(:)
        complex(dp), intent(out) :: y(:)
        integer, intent(out) :: ierr
        complex(dp) :: mm(k, num_spec), cc(k, num_spec)
        integer :: a
        mm = reshape(x, shape(mm))
        call mix_collision_moments(ailmm_aa, mm, cc, ierr)
        if (ierr /= 0) return
        do a = 1, num_spec
            call map%lift(cc(:, a), y((a - 1)*n_2d_size + 1:a*n_2d_size), ierr)
            if (ierr /= 0) return
        end do
    end subroutine
    subroutine inverse_diagonal(x, y, ierr)
        complex(dp), intent(in) :: x(:)
        complex(dp), intent(out) :: y(:)
        integer, intent(out) :: ierr
        integer :: j
        do j = 1, size(x)
            y(j) = x(j)/(2.0_dp + real(j, dp)/size(x))
        end do
        ierr = 0
    end subroutine
    subroutine original_all(x, y, ierr)
        complex(dp), intent(in) :: x(:)
        complex(dp), intent(out) :: y(:)
        integer, intent(out) :: ierr
        complex(dp) :: mx(k*num_spec), native_out(n_2d_size)
        integer :: a, j
        call project_all(x, mx, ierr)
        if (ierr /= 0) return
        gathered = reshape(mx, shape(gathered))
        do a = 1, num_spec
            ispec = a - 1
            call integral_part(x((a - 1)*n_2d_size + 1:a*n_2d_size), native_out)
            y((a - 1)*n_2d_size + 1:a*n_2d_size) = -native_out
        end do
        do j = 1, size(x)
            y(j) = y(j) + (2.0_dp + real(j, dp)/size(x))*x(j)
        end do
    end subroutine
end module native_moment_fixture

program test_collision_moment_map
    use native_moment_fixture, only: run_native_moment_tests
    implicit none
    call run_native_moment_tests()
end program test_collision_moment_map
