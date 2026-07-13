program test_cquad_ieee_guard
    use, intrinsic :: ieee_arithmetic, only: ieee_get_flag, ieee_invalid, &
        ieee_set_flag, ieee_set_halting_mode
    use gsl_integration_routines_mod, only: fint1d_cquad
    implicit none

    character(len=1024) :: argument, command, executable
    integer :: command_status, exit_status
    logical :: invalid_raised
    real(kind(1.0d0)) :: result(2)

    call get_command_argument(1, argument)
    if (trim(argument) == 'invalid-child') then
        result = fint1d_cquad(invalid_integrand, 0.0d0, 1.0d0, 1.0d-13, 1.0d-13)
        error stop 'FAIL: CQUAD accepted an invalid integrand'
    end if

    call ieee_set_flag(ieee_invalid, .false.)
    call ieee_set_halting_mode(ieee_invalid, .true.)
    result = fint1d_cquad(zero_integrand, 0.0d0, 1.0d0, 1.0d-13, 1.0d-13)
    call ieee_get_flag(ieee_invalid, invalid_raised)
    if (invalid_raised) error stop 'FAIL: CQUAD leaked an internal invalid flag'
    if (any(result /= 0.0d0)) error stop 'FAIL: zero integrand gave a nonzero result'

    call get_command_argument(0, executable)
    command = '"' // trim(executable) // '" invalid-child'
    call execute_command_line(trim(command), wait=.true., exitstat=exit_status, &
        cmdstat=command_status)
    if (command_status /= 0) error stop 'FAIL: invalid-integrand child did not launch'
    if (exit_status == 0) error stop 'FAIL: invalid integrand was not rejected'
    print *, 'All tests passed!'

contains

    function zero_integrand(x) result(value)
        real(kind(1.0d0)) :: x
        real(kind(1.0d0)) :: value

        value = 0.0d0*x
    end function zero_integrand

    function invalid_integrand(x) result(value)
        real(kind(1.0d0)) :: x
        real(kind(1.0d0)), volatile :: denominator
        real(kind(1.0d0)) :: value

        denominator = x - x
        value = denominator/denominator
    end function invalid_integrand
end program test_cquad_ieee_guard
