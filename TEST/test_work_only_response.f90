program test_work_only_response
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use work_only_response_mod
    use response_witness_mod
    use hdf5_tools, only: hid_t, h5_open, h5_get, h5_close
    implicit none
    complex(dp) :: a(3,3), integral(3,3), ln(3,3), bw(3), bc(3), y1(3)
    complex(dp) :: shift(3), bra(3), exact(3), got(3), phase, error_vector(3)
    complex(dp) :: bnodes(3), derivative(3), phases(4)
    real(dp) :: saved(3), norms(4), moments(5), reference_work, saved_cpu
    real(dp) :: expected_norm, nan
    integer(hid_t) :: file_id
    integer :: j, unit, qualified, gate, solves, applications
    character(len=160) :: filename
    logical :: inaccurate

    if (lsw_work_only_response) error stop 'FAIL enabled by default'
    if (.not. work_only_scope_supported(.true., .true., 3, 1, 2, 1, 8, 0, 1, 0)) &
        error stop 'FAIL supported scope rejected'
    if (work_only_scope_supported(.false., .true., 3, 1, 2, 1, 8, 0, 1, 0)) &
        error stop 'FAIL missing raw capture accepted'
    if (work_only_scope_supported(.true., .false., 3, 1, 2, 1, 8, 0, 1, 0)) &
        error stop 'FAIL single species scope accepted'
    if (work_only_scope_supported(.true., .true., 2, 1, 2, 1, 8, 0, 1, 0)) &
        error stop 'FAIL different solver accepted'
    if (work_only_scope_supported(.true., .true., 3, 0, 2, 1, 8, 0, 1, 0)) &
        error stop 'FAIL missing NA response accepted'
    if (work_only_scope_supported(.true., .true., 3, 1, 1, 1, 8, 0, 1, 0)) &
        error stop 'FAIL self-consistent Er accepted'
    if (work_only_scope_supported(.true., .true., 3, 1, 2, 0, 8, 0, 1, 0)) &
        error stop 'FAIL unsupported coordinates accepted'
    if (work_only_scope_supported(.true., .true., 3, 1, 2, 1, 10, 0, 1, 0)) &
        error stop 'FAIL lab10 accepted'
    if (work_only_scope_supported(.true., .true., 3, 1, 2, 1, 8, 1, 1, 0)) &
        error stop 'FAIL magnetic shear accepted'
    if (work_only_scope_supported(.true., .true., 3, 1, 2, 1, 8, 0, 0, 0)) &
        error stop 'FAIL nonaxisymmetric background accepted'
    if (work_only_scope_supported(.true., .true., 3, 1, 2, 1, 8, 0, 1, 1)) &
        error stop 'FAIL nonaxisymmetric harmonic accepted'
    nan = ieee_value(0.0_dp, ieee_quiet_nan)
    derivative = 0.0_dp
    bnodes = cmplx(0.0_dp, 1.0e-3_dp, dp)
    if (.not. work_only_constant_b(bnodes, derivative, 10.0_dp)) &
        error stop 'FAIL constant complex phase rejected'
    bnodes(3) = -bnodes(3)
    if (work_only_constant_b(bnodes, derivative, 10.0_dp)) &
        error stop 'FAIL amplitude-only gate missed varying phase'
    bnodes = cmplx(0.0_dp, 1.0e-3_dp, dp)
    derivative(1) = 1.0e-6_dp
    if (work_only_constant_b(bnodes, derivative, 10.0_dp)) &
        error stop 'FAIL nonconstant derivative accepted'
    derivative = 0.0_dp
    bnodes(3) = cmplx(nan, 0.0_dp, dp)
    if (work_only_constant_b(bnodes, derivative, 10.0_dp)) &
        error stop 'FAIL nonfinite perturbation accepted'
    bnodes = 0.0_dp
    if (work_only_constant_b(bnodes, derivative, 10.0_dp)) &
        error stop 'FAIL zero perturbation accepted'

    ! Differential/streaming A and all-pair integral I are deliberately separate.
    ! Row1,col3 is a complex periodic seam; the remaining off-diagonal entries
    ! couple receiving blocks. Dropping I or any row changes the oracle.
    a = 0.0_dp
    a(1,1) = cmplx(2.0_dp, 1.0_dp, dp)
    a(2,2) = cmplx(3.0_dp, -1.0_dp, dp)
    a(3,3) = cmplx(4.0_dp, 0.5_dp, dp)
    a(1,3) = cmplx(0.0_dp, 0.3_dp, dp)
    integral = 0.0_dp
    integral(2,1) = 0.2_dp
    integral(3,2) = cmplx(0.0_dp, -0.4_dp, dp)
    ln = a-integral
    phases = [cmplx(1.0_dp,0.0_dp,dp), cmplx(0.0_dp,1.0_dp,dp), &
              cmplx(-1.0_dp,0.0_dp,dp), cmplx(0.0_dp,-1.0_dp,dp)]
    lsw_response_witness = .true.
    do j = 1, 5
        phase = phases(mod(j-1,4)+1)
        bw = phase*[cmplx(1.0_dp,0.0_dp,dp), &
                    cmplx(0.0_dp,2.0_dp,dp), cmplx(-1.0_dp,0.5_dp,dp)]
        ! Closed elimination oracle, independent of the callback's pivot solve.
        exact(1) = (bw(1)-ln(1,3)*bw(3)/ln(3,3) + &
            ln(1,3)*ln(3,2)*bw(2)/(ln(2,2)*ln(3,3))) / &
            (ln(1,1)+ln(1,3)*ln(3,2)*ln(2,1)/(ln(2,2)*ln(3,3)))
        exact(2) = (bw(2)-ln(2,1)*exact(1))/ln(2,2)
        exact(3) = (bw(3)-ln(3,2)*exact(2))/ln(3,3)
        y1 = phase*[cmplx(0.2_dp,0.7_dp,dp), &
                    cmplx(-0.3_dp,0.4_dp,dp), cmplx(0.1_dp,-0.6_dp,dp)]
        bc = matmul(ln,y1)
        shift = phase*[0.2_dp, -0.5_dp, 0.3_dp]
        bra = cmplx(0.0_dp,1.0_dp,dp)*phase*[0.7_dp, 1.2_dp, -0.8_dp]
        inaccurate = j == 5
        error_vector = 0.0_dp
        if (inaccurate) error_vector(3) = cmplx(0.0_dp,0.125_dp,dp)
        solves = 0
        applications = 0
        call response_witness_begin(j, 0, file_id, filename)
        call evaluate_work_only_response(file_id, bw, bc, y1, shift, bra, &
            solve, apply_ln, sum_two_ranks)
        call response_witness_end(file_id)
        if (solves /= 1 .or. applications /= 4) &
            error stop 'FAIL diagnostic solve/operator count'
        call h5_open(trim(filename), file_id)
        call h5_get(file_id, 'work_only_Ywork_real', saved)
        got = saved
        call h5_get(file_id, 'work_only_Ywork_imag', saved)
        got = got+cmplx(0.0_dp,1.0_dp,dp)*saved
        if (maxval(abs(got-exact-error_vector)) > 2.0e-14_dp) &
            error stop 'FAIL coupled solution differs from analytic oracle'
        call h5_get(file_id, 'work_only_global_l2_norms', norms)
        expected_norm = sqrt(2.0_dp*sum(abs(matmul(ln,error_vector))**2))
        if (abs(norms(1)-expected_norm) > 2.0e-14_dp) &
            error stop 'FAIL true global unpreconditioned residual'
        if (abs(norms(2)-sqrt(2.0_dp*sum(abs(bw)**2))) > 2.0e-14_dp) &
            error stop 'FAIL all receiving contributions missing'
        if (abs(norms(3)-expected_norm) > 2.0e-14_dp .or. norms(4) > 2.0e-14_dp) &
            error stop 'FAIL discrete remainder identity'
        call h5_get(file_id, 'work_only_receiving_raw_moments', moments)
        if (abs(moments(2)-0.5_dp*real(dot_product(bra,exact+error_vector),dp)) &
                > 2.0e-14_dp) error stop 'FAIL work phase or conjugation'
        if (abs(moments(3)) > 2.0e-14_dp .or. abs(moments(5)) > 2.0e-14_dp) &
            error stop 'FAIL native zero shift work'
        if (j == 1) reference_work = moments(2)
        if (j <= 4 .and. abs(moments(2)-reference_work) > 2.0e-14_dp) &
            error stop 'FAIL four-phase invariance'
        call h5_get(file_id, 'work_only_shift_work_gate_pass', gate)
        if (gate /= 1) error stop 'FAIL shift-work gate'
        call h5_get(file_id, 'work_only_native_solve_qualified', qualified)
        if (qualified /= 0) error stop 'FAIL premature native qualification'
        call h5_get(file_id, 'work_only_cpu_seconds', saved_cpu)
        if (saved_cpu < 0.0_dp) error stop 'FAIL missing cost'
        call h5_close(file_id)
        open(newunit=unit, file=trim(filename), status='old')
        close(unit, status='delete')
    end do
    print *, 'All tests passed!'
contains
    subroutine apply_ln(value, action)
        complex(dp), intent(in) :: value(:)
        complex(dp), intent(out) :: action(:)
        applications = applications+1
        action = matmul(a,value)-matmul(integral,value)
    end subroutine apply_ln

    subroutine solve(value)
        complex(dp), intent(inout) :: value(:)
        complex(dp) :: matrix(3,3), factor
        integer :: row, column

        solves = solves+1
        matrix = ln
        do column = 1, 2
            do row = column+1, 3
                factor = matrix(row,column)/matrix(column,column)
                matrix(row,column:3) = &
                    matrix(row,column:3)-factor*matrix(column,column:3)
                value(row) = value(row)-factor*value(column)
            end do
        end do
        do row = 3, 1, -1
            value(row) = (value(row)-sum(matrix(row,row+1:3)*value(row+1:3))) &
                /matrix(row,row)
        end do
        value = value+error_vector
    end subroutine solve

    function sum_two_ranks(value) result(total)
        real(dp), intent(in) :: value
        real(dp) :: total
        ! Test the receiving-rank reduction contract; this is not an MPI test.
        total = 2.0_dp*value
    end function sum_two_ranks
end program test_work_only_response
