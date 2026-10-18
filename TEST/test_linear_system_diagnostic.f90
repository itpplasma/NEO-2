program test_linear_system_diagnostic
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use linear_system_diagnostic_mod, only: diagnostic_system, diagnostic_array
    implicit none
    integer, parameter :: rows(4) = [1,2,1,2], cols(4) = [1,1,2,2]
    complex(dp) :: values(4), bra(3,2), rhs(2,3), solution(2,3)
    integer :: channel
    character(1) :: label

    ! Nonsymmetric complex triangular system with a closed-form solution.
    ! The Python verifier is a separate dense LAPACK solve of the exported A,b.
    values = [cmplx(2.0_dp,1.0_dp,dp), cmplx(0.0_dp,0.0_dp,dp), &
              cmplx(3.0_dp,-1.0_dp,dp), cmplx(4.0_dp,2.0_dp,dp)]
    bra(1,:) = [cmplx(1.0_dp,2.0_dp,dp), cmplx(3.0_dp,-1.0_dp,dp)]
    bra(2,:) = [cmplx(-2.0_dp,0.0_dp,dp), cmplx(1.0_dp,4.0_dp,dp)]
    bra(3,:) = [cmplx(0.0_dp,1.0_dp,dp), cmplx(-1.0_dp,2.0_dp,dp)]
    rhs(:,1) = [cmplx(1.0_dp,1.0_dp,dp), cmplx(2.0_dp,-1.0_dp,dp)]
    rhs(:,2) = 0.0_dp
    rhs(:,3) = [cmplx(-1.0_dp,2.0_dp,dp), cmplx(3.0_dp,1.0_dp,dp)]
    do channel = 1, 3
        solution(2,channel) = rhs(2,channel)/values(4)
        solution(1,channel) = &
            (rhs(1,channel)-values(3)*solution(2,channel))/values(1)
    end do
    call diagnostic_system(.true., 0, 2, rows, cols, values, bra)
    do channel = 1, 3
        write(label, '(i1)') channel
        call diagnostic_array(.true., 0, 'source'//label//'_rhs', &
                              rhs(:,channel:channel))
        call diagnostic_array(.true., 0, 'source'//label//'_solution', &
                              solution(:,channel:channel))
    end do
    call diagnostic_system(.false., 0, 2, rows, cols, values, bra)
    call diagnostic_array(.false., 0, 'drive0_geometric_source', 0.3_dp*rhs)
    call diagnostic_array(.false., 0, 'drive0_total_source', rhs)
    do channel = 1, 3
        write(label, '(i1)') channel
        call diagnostic_array(.false., 0, 'source'//label//'_rhs', &
                              rhs(:,channel:channel))
        call diagnostic_array(.false., 0, 'source'//label//'_solution', &
                              solution(:,channel:channel))
    end do
end program test_linear_system_diagnostic
