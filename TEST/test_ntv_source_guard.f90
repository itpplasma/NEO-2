program test_ntv_source_guard
    use neo2_ql, only: init, check
    use ntv_mod, only: isw_qflux_NA, in_file_pert, isw_hel_drive, &
        hel_brad_re, hel_brad_im, hel_phim_re, hel_phim_im
    implicit none

    character(len=16) :: test_case

    call get_command_argument(1, test_case)
    call init()

    isw_qflux_NA = 1
    in_file_pert = 'none'
    isw_hel_drive = 0
    hel_brad_re = 0.0d0
    hel_brad_im = 0.0d0
    hel_phim_re = 0.0d0
    hel_phim_im = 0.0d0

    select case (trim(test_case))
    case ('zero')
        call check()
        error stop 'returned from rejected source'
    case ('driven')
        isw_hel_drive = 1
        hel_phim_re = 1.0d0
        call check()
        print *, 'Driven source accepted!'
    case default
        error stop 'FAIL: unknown test case'
    end select
end program test_ntv_source_guard
