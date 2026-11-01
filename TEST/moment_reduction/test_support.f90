module moment_test_support
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
    integer :: checks = 0
contains
    subroutine check(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label
        checks = checks + 1
        if (.not. condition) then
            print *, 'FAIL: ', label
            error stop 1
        end if
    end subroutine
    subroutine finished(label)
        character(*), intent(in) :: label
        print '(a,a,i0,a)', trim(label), ': All tests passed! (checks=', checks, ')'
    end subroutine
end module
