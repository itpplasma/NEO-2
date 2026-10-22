! Round trip through the production reader and harmonic evaluator. The caller
! repeats the single-surface spectrum at four radii for the native cubic spline.
program test_phi1_read
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neo_control, only: inp_swi, lab_swi
    use ntv_mod, only: in_file_pert
    use neo_magfie_perturbation, only: neo_magfie_pert
    use neo_spline_data, only: lsw_linear_boozer
    use partpa_mod, only: bmod0
    implicit none
    integer :: j
    complex(dp) :: value
    call get_command_argument(1,in_file_pert)
    inp_swi = 9
    lab_swi = 10
    lsw_linear_boozer = .false.
    bmod0 = 1.0_dp
    do j = 1, 3
        call neo_magfie_pert([.5_dp,.79_dp*j,.31_dp*j],value)
        write(*,'(2(es25.17,1x))') real(value,dp), aimag(value)
    end do
end program
