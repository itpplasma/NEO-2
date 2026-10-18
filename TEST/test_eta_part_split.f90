program test_eta_part_split
    use flint_mod, only: split_eta_part_trapped
    implicit none
    character(len=32) :: mode
    integer :: n_passing, n_trapped

    mode = ''
    if (command_argument_count() > 0) call get_command_argument(1, mode)

    select case (trim(mode))
    case ('')
        ! eta_part is the total; the trapped request moves levels from passing.
        call split_eta_part_trapped(150, 70, n_passing, n_trapped)
        call check(n_passing == 80 .and. n_trapped == 70, '150/70 split')
        ! Holding the passing count while refining trapped levels.
        call split_eta_part_trapped(360, 280, n_passing, n_trapped)
        call check(n_passing == 80 .and. n_trapped == 280, '360/280 split')
        call split_eta_part_trapped(6, 3, n_passing, n_trapped)
        call check(n_passing == 3 .and. n_trapped == 3, 'minimal split')
        print *, 'All tests passed!'
    case ('starved-passing')
        ! Former silent clamp: 3 passing + 210 trapped, total 213 /= 150.
        call split_eta_part_trapped(150, 210, n_passing, n_trapped)
    case ('two-passing')
        call split_eta_part_trapped(150, 148, n_passing, n_trapped)
    case ('few-trapped')
        call split_eta_part_trapped(150, 2, n_passing, n_trapped)
    case default
        error stop 'unknown mode'
    end select
contains
    subroutine check(condition, label)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: label
        if (.not. condition) then
            print *, 'FAIL: ', label
            error stop 1
        end if
    end subroutine check
end program test_eta_part_split
