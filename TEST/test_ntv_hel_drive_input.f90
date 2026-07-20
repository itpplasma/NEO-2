program test_ntv_hel_drive_input
    use nrtype, only: dp
    use ntv_mod, only: isw_hel_drive, hel_brad_re, hel_brad_im, &
        hel_phim_re, hel_phim_im, has_helical_drive_source
    implicit none

    integer :: failures

    failures = 0

    call set_drive(0, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp)
    call expect_source(.false., failures)

    call set_drive(1, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp)
    call expect_source(.false., failures)

    call set_drive(1, 1.0e-3_dp, 0.0_dp, 0.0_dp, 0.0_dp)
    call expect_source(.true., failures)

    call set_drive(1, 0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp)
    call expect_source(.true., failures)

    if (failures == 0) then
        print *, 'All tests passed!'
    else
        error stop 'helical-drive source detection failed'
    end if

contains

    subroutine set_drive(switch, brad_re, brad_im, phim_re, phim_im)
        integer, intent(in) :: switch
        real(dp), intent(in) :: brad_re, brad_im, phim_re, phim_im

        isw_hel_drive = switch
        hel_brad_re = brad_re
        hel_brad_im = brad_im
        hel_phim_re = phim_re
        hel_phim_im = phim_im
    end subroutine set_drive

    subroutine expect_source(expected, status)
        logical, intent(in) :: expected
        integer, intent(inout) :: status

        if (has_helical_drive_source() .neqv. expected) then
            print *, 'FAIL: unexpected helical-drive source state'
            status = status + 1
        end if
    end subroutine expect_source

end program test_ntv_hel_drive_input
