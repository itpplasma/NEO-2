! Application-specific SAME-L_N diagnostic. No new collision model or solver.
module work_only_response_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use hdf5_tools, only: hid_t, h5_add
    use response_witness_mod, only: response_witness_complex
    implicit none
    private
    public :: lsw_work_only_response, work_only_scope_supported
    public :: work_only_constant_b, evaluate_work_only_response

    logical :: lsw_work_only_response = .false.

    abstract interface
        subroutine solve_work_rhs(value)
            import dp
            complex(dp), intent(inout) :: value(:)
        end subroutine solve_work_rhs
        subroutine apply_native_ln(value, action)
            import dp
            complex(dp), intent(in) :: value(:)
            complex(dp), intent(out) :: action(:)
        end subroutine apply_native_ln
        function sum_receiving_ranks(value) result(total)
            import dp
            real(dp), intent(in) :: value
            real(dp) :: total
        end function sum_receiving_ranks
    end interface
contains
    pure logical function work_only_scope_supported(capture, multispecies, &
            solver, qflux_na, er_mode, coordinates, laboratory, shear, axisym, &
            equilibrium_n) result(supported)
        logical, intent(in) :: capture, multispecies
        integer, intent(in) :: solver, qflux_na, er_mode, coordinates
        integer, intent(in) :: laboratory, shear, axisym, equilibrium_n

        supported = capture .and. multispecies .and. solver == 3 .and. &
            qflux_na == 1 .and. er_mode == 2 .and. coordinates /= 0 .and. &
            (laboratory == 8 .or. laboratory == 9) .and. shear == 0 .and. &
            axisym == 1 .and. equilibrium_n == 0
    end function work_only_scope_supported

    pure logical function work_only_constant_b(b, db, phi_span) result(supported)
        complex(dp), intent(in) :: b(:), db(:)
        real(dp), intent(in) :: phi_span
        real(dp) :: threshold

        supported = .false.
        if (size(b) == 0 .or. size(db) /= size(b)) return
        if (.not. all(ieee_is_finite(real(b, dp)))) return
        if (.not. all(ieee_is_finite(aimag(b)))) return
        if (.not. all(ieee_is_finite(real(db, dp)))) return
        if (.not. all(ieee_is_finite(aimag(db)))) return
        if (.not. ieee_is_finite(phi_span)) return
        if (abs(b(1)) <= tiny(1.0_dp)) return
        threshold = 128.0_dp*epsilon(1.0_dp)*abs(b(1))
        if (maxval(abs(b-b(1))) > threshold) return
        if (maxval(abs(db))*max(1.0_dp, abs(phi_span)) > threshold) return
        supported = .true.
    end function work_only_constant_b

    ! Callbacks retain the native coupled solver and its unpreconditioned
    ! operator. All ranks enter every callback, including zero local RHS ranks.
    ! Vector algebra here is assembly in native coefficient coordinates.
    subroutine evaluate_work_only_response(file_id, bwork, bcomplete, y1, &
            shift, work_bra, solve_rhs, apply_ln, global_sum)
        integer(hid_t), intent(in) :: file_id
        complex(dp), intent(in) :: bwork(:), bcomplete(:), y1(:), shift(:)
        complex(dp), intent(in) :: work_bra(:)
        procedure(solve_work_rhs) :: solve_rhs
        procedure(apply_native_ln) :: apply_ln
        procedure(sum_receiving_ranks) :: global_sum
        complex(dp), allocatable :: ywork(:), action(:), remainder(:), z(:)
        real(dp) :: norms(4), moments(5), shift_scale, cpu_start, cpu_end
        integer(int64) :: wall_start, wall_end, wall_rate

        call cpu_time(cpu_start)
        call system_clock(wall_start, wall_rate)
        allocate(ywork(size(bwork)), action(size(bwork)), &
            remainder(size(bwork)), z(size(bwork)))
        call h5_add(file_id, 'work_only_schema_version', 1)
        call h5_add(file_id, 'work_only_native_solve_qualified', 0)
        call h5_add(file_id, 'work_only_nr_continuum_identification', 0)
        call h5_add(file_id, 'work_only_operator_applications', 4)
        call h5_add(file_id, 'work_only_additional_solve_count', 1)
        call h5_add(file_id, 'work_only_operator', &
            'L_N*x=A_n*x-integral_part(x); full nonaxisymmetric native rows')
        call h5_add(file_id, 'work_only_cost_scope', &
            'one cached solve, four L_N applications and diagnostic vector I/O')
        call response_witness_complex(file_id, 'work_only_Bwork', bwork)
        call response_witness_complex(file_id, 'work_only_Bcomplete', bcomplete)
        call response_witness_complex(file_id, 'work_only_Y1physical', y1)
        call response_witness_complex(file_id, 'work_only_M', shift)

        ywork = bwork
        call solve_rhs(ywork)
        call response_witness_complex(file_id, 'work_only_Ywork', ywork)
        call apply_ln(ywork, action)
        call response_witness_complex(file_id, 'work_only_LnYwork', action)
        action = action-bwork
        call response_witness_complex(file_id, 'work_only_residual', action)
        norms(1) = sqrt(global_sum(sum(abs(action)**2)))
        norms(2) = sqrt(global_sum(sum(abs(bwork)**2)))

        call apply_ln(shift, action)
        call response_witness_complex(file_id, 'work_only_LnM', action)
        remainder = bcomplete+action-bwork
        call response_witness_complex(file_id, 'work_only_RN', remainder)
        z = y1+shift-ywork
        call response_witness_complex(file_id, 'work_only_Z', z)
        call apply_ln(z, action)
        call response_witness_complex(file_id, 'work_only_LnZ', action)
        action = action-remainder
        call response_witness_complex(file_id, 'work_only_identity_residual', action)
        norms(3) = sqrt(global_sum(sum(abs(action)**2)))
        call apply_ln(y1, action)
        action = action-bcomplete
        call response_witness_complex(file_id, 'work_only_forward_residual', action)
        norms(4) = sqrt(global_sum(sum(abs(action)**2)))
        call h5_add(file_id, 'work_only_global_l2_norms', norms, [1], [4])

        ! Same raw receiving moment as native qflux(1,...), including 1/2.
        ! DOT_PRODUCT conjugates its first complex argument exactly once.
        moments(1) = 0.5_dp*real(dot_product(work_bra, y1), dp)
        moments(2) = 0.5_dp*real(dot_product(work_bra, ywork), dp)
        moments(3) = 0.5_dp*real(dot_product(work_bra, shift), dp)
        moments(4) = 0.5_dp*real(dot_product(work_bra, z), dp)
        moments(5) = moments(4)-moments(1)+moments(2)
        shift_scale = 0.5_dp*sum(abs(work_bra)*abs(shift))
        call h5_add(file_id, 'work_only_receiving_raw_moments', &
            moments, [1], [5])
        call h5_add(file_id, 'work_only_shift_work_scale', shift_scale)
        call h5_add(file_id, 'work_only_shift_work_gate_pass', &
            merge(1, 0, abs(moments(3)) <= &
                128.0_dp*epsilon(1.0_dp)*max(shift_scale, tiny(1.0_dp))))
        call cpu_time(cpu_end)
        call system_clock(wall_end)
        call h5_add(file_id, 'work_only_cpu_seconds', cpu_end-cpu_start)
        call h5_add(file_id, 'work_only_wall_seconds', &
            real(wall_end-wall_start, dp)/real(wall_rate, dp))
        call h5_add(file_id, 'work_only_diagnostic_complete', 1)
    end subroutine evaluate_work_only_response
end module work_only_response_mod
