program test_ntv_work_export
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan, &
                                            ieee_is_finite
    use ntv_work_export, only: ntv_work_export_t, work_begin, work_basis, &
                               work_vector, work_end, work_node_deltas
    implicit none
    integer, parameter :: n = 8201
    type(ntv_work_export_t) :: export
    complex(dp) :: flux(3, n), rhs(n), solution(n), solution3(n)
    real(dp) :: coefficients(4, 4), x, q1, q3
    real(dp) :: delta_positive(0:2), delta_negative(0:2), active_deltas(2)
    integer :: i, j, unit
    logical :: exists
    character(128) :: prefix

    delta_positive = [ieee_value(0.0_dp, ieee_quiet_nan), 1.0_dp, 2.0_dp]
    delta_negative = [3.0_dp, 4.0_dp, ieee_value(0.0_dp, ieee_quiet_nan)]
    call work_node_deltas(0, 0, 2, delta_positive, delta_negative, active_deltas)
    call check_deltas(active_deltas, [0.0_dp, 3.0_dp])
    call work_node_deltas(2, 0, 2, delta_positive, delta_negative, active_deltas)
    call check_deltas(active_deltas, [2.0_dp, 0.0_dp])
    call work_node_deltas(1, 0, 2, delta_positive, delta_negative, active_deltas)
    call check_deltas(active_deltas, [1.0_dp, 4.0_dp])
    delta_positive(1) = ieee_value(0.0_dp, ieee_quiet_nan)
    call work_node_deltas(1, 0, 2, delta_positive, delta_negative, active_deltas)
    if (.not. ieee_is_nan(active_deltas(1))) error stop 'Active NaN masked'
    if (.not. ieee_is_finite(active_deltas(2))) error stop 'Finite delta lost'
    if (abs(active_deltas(2) - 4.0_dp) > 0.0_dp) error stop 'Active delta changed'

    prefix = 'NEO2_WORK_s0_j1'
    call remove('NEO2_WORK.REQUEST')
    call remove(trim(prefix)//'.txt')
    call remove(trim(prefix)//'_basis.bin')
    call remove(trim(prefix)//'_FLUX_k1_b0.bin')
    call remove(trim(prefix)//'_RHS_k1_b0.bin')
    call remove(trim(prefix)//'_RHS_k3_b0.bin')
    call remove(trim(prefix)//'_SOL_k1_b0.bin')
    call remove(trim(prefix)//'_SOL_k3_b0.bin')
    call work_begin(export, 0.09_dp, 0, n, 2)
    if (export%enabled) error stop 'Default export must be disabled'
    inquire(file=trim(prefix)//'.txt', exist=exists)
    if (exists) error stop 'Disabled export wrote a file'
    open(newunit=unit, file='NEO2_WORK.REQUEST', status='new')
    write(unit, *) 0.09_dp, 0
    close(unit)
    call work_begin(export, 0.09_dp, 0, n, 2)
    if (.not. export%enabled) error stop 'Selected work export not enabled'
    flux = cmplx(ieee_value(0.0_dp, ieee_quiet_nan), 0.0_dp, dp)
    do i = 1, n
        x = real(i, dp)
        flux(1, i) = cmplx(x, 2*x - 3, dp)
        rhs(i) = cmplx(3*x + 7, x*x - 1, dp)
        solution(i) = cmplx(-x + 2, 4*x + 1, dp)
        solution3(i) = cmplx(2*x + 3, -x + 5, dp)
    end do
    call work_vector(export, 'RHS', 1, 0, rhs)
    call work_vector(export, 'RHS', 3, 0, rhs)
    rhs = solution
    ! Actual rank-two, strided flux row; inactive rows deliberately contain NaNs.
    call work_vector(export, 'FLUX', 1, 0, flux(1, :))
    call work_vector(export, 'SOL', 1, 0, rhs)
    call work_vector(export, 'SOL', 3, 0, solution3)
    do j = 1, 4
        do i = 1, 4
            coefficients(i, j) = real(i + 10*j, dp)
        end do
    end do
    call work_basis(export, 0, 1, 1, 0.5_dp, 1.0_dp, 2.0_dp, coefficients)
    call work_basis(export, 1, 2, 1, 0.25_dp, 1.0_dp, 2.0_dp, coefficients)
    q1 = 0.5_dp*real(sum(conjg(flux(1, :))*solution), dp)
    q3 = 0.5_dp*real(sum(conjg(flux(1, :))*solution3), dp)
    write(export%text_unit, *) 'QFLUX ', 0, q1, q3
    call work_end(export)
    print *, 'All tests passed!'

contains
    subroutine check_deltas(values, expected)
        real(dp), intent(in) :: values(2), expected(2)

        if (.not. all(ieee_is_finite(values))) error stop 'Inactive delta read'
        if (any(abs(values - expected) > 0.0_dp)) error stop 'Delta changed'
    end subroutine

    subroutine remove(path)
        character(*), intent(in) :: path
        integer :: file_unit, ios

        open(newunit=file_unit, file=path, status='old', iostat=ios)
        if (ios /= 0) return
        close(file_unit, status='delete')
    end subroutine
end program
