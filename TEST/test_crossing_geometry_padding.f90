program test_crossing_geometry_padding
    ! Regression for the modify_br = 1 crossing path. rearrange_phideps
    ! receives the padded array extent iend = ub_mag + 2 while phi_divide
    ! is allocated (1:ub_mag). The crossing scan must ignore the duplicated
    ! padding nodes, alignment must match the unpadded result bit for bit,
    ! and every refinement flag must land inside phi_divide(1:ub_mag) so
    ! the caller's ierr = 3 gate keeps firing.
    !
    ! Profile: phi(i) = 0.05 i on 0:40, bhat = 1 + (phi - 1)^2/2, one pitch
    ! band at bhat_cross = 1.1. The two simple crossings sit at
    ! phi = 1 -/+ sqrt(0.2), aligned onto nodes 11 and 29, both more than
    ! nstepmin nodes away from either edge so no legacy sub-interval split
    ! can trigger.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use device_mod, only: fieldpropagator

    implicit none

    integer, parameter :: ub_mag = 40, npart = 1
    real(dp), parameter :: bhat_cross = 1.1_dp
    real(dp), parameter :: root_tolerance = 1.0e-9_dp
    real(dp), parameter :: band_tolerance = 1.0e-12_dp

    real(dp) :: phi_ref(0:ub_mag), bhat_ref(0:ub_mag)
    real(dp) :: phi_out(0:ub_mag), bhat_out(0:ub_mag)
    real(dp) :: root_left, root_right
    integer :: divide_ref(1:ub_mag), divide_out(1:ub_mag)

    allocate (fieldpropagator)

    root_left = 1.0_dp - sqrt(0.2_dp)
    root_right = 1.0_dp + sqrt(0.2_dp)

    ! Control: no padding (modify_br = 0 layout).
    call run_case(.false., 0.0_dp, .false., phi_ref, bhat_ref, divide_ref)
    if (maxval(divide_ref) /= 1) &
        error stop 'FAIL: unpadded case requested refinement'
    if (abs(bhat_ref(11)/bhat_cross - 1.0_dp) > band_tolerance) &
        error stop 'FAIL: unpadded left crossing was not aligned'
    if (abs(bhat_ref(29)/bhat_cross - 1.0_dp) > band_tolerance) &
        error stop 'FAIL: unpadded right crossing was not aligned'
    if (abs(phi_ref(11) - root_left) > root_tolerance) &
        error stop 'FAIL: unpadded left crossing node misses the root'
    if (abs(phi_ref(29) - root_right) > root_tolerance) &
        error stop 'FAIL: unpadded right crossing node misses the root'

    ! Padded, boundary field raised above every physical value.
    call run_case(.true., 1.6_dp, .false., phi_out, bhat_out, divide_out)
    if (maxval(divide_out) /= 1) &
        error stop 'FAIL: benign padding triggered spurious refinement'
    if (any(phi_out /= phi_ref)) &
        error stop 'FAIL: benign padding changed the aligned coordinates'
    if (any(bhat_out /= bhat_ref)) &
        error stop 'FAIL: benign padding changed the aligned field'

    ! Padded, boundary field below the crossing band: the padding nodes
    ! look like an extra pitch crossing and must still be ignored.
    call run_case(.true., 1.05_dp, .false., phi_out, bhat_out, divide_out)
    if (maxval(divide_out) /= 1) &
        error stop 'FAIL: low padding triggered spurious refinement'
    if (any(phi_out /= phi_ref)) &
        error stop 'FAIL: low padding changed the aligned coordinates'
    if (any(bhat_out /= bhat_ref)) &
        error stop 'FAIL: low padding changed the aligned field'

    ! Padded, non-monotone physical field: the rejection flag must land
    ! inside phi_divide(1:ub_mag) so maxval(phi_divide) > 1 in the caller.
    call run_case(.true., 1.6_dp, .true., phi_out, bhat_out, divide_out)
    if (divide_out(3) /= 2) &
        error stop 'FAIL: physical rejection flag missed interval 3'
    if (sum(divide_out) /= ub_mag + 1) &
        error stop 'FAIL: physical rejection flagged extra intervals'

    print '(a)', 'All tests passed!'

contains

    subroutine run_case(padded, bhat_pad, bump, phi_physical, bhat_physical, &
            divide)
        logical, intent(in) :: padded, bump
        real(dp), intent(in) :: bhat_pad
        real(dp), intent(out) :: phi_physical(0:ub_mag)
        real(dp), intent(out) :: bhat_physical(0:ub_mag)
        integer, intent(out) :: divide(1:ub_mag)

        interface
            subroutine rearrange_phideps(ibeg, iend, ub_mag, npart, ncomp, &
                    nreal, bcovar_column, subsqmin, phi_divide, phi_mfl, &
                    bhat_mfl, dlogbdphi_mfl, dbcovar_s_hat_dphi_mfl, &
                    arr_real, arr_comp, eta, delt_pos, delt_neg, &
                    fact_pos_b, fact_neg_b, fact_pos_e, fact_neg_e)
                import :: dp
                integer :: ibeg, iend, npart, ncomp, nreal, bcovar_column
                integer, intent(in) :: ub_mag
                real(dp) :: subsqmin
                integer, dimension(1:ub_mag) :: phi_divide
                real(dp), dimension(ibeg:iend) :: phi_mfl, bhat_mfl
                real(dp), dimension(ibeg:iend) :: dlogbdphi_mfl
                real(dp), dimension(ibeg:iend) :: dbcovar_s_hat_dphi_mfl
                real(dp), dimension(ibeg:iend, nreal) :: arr_real
                complex(dp), dimension(ibeg:iend, ncomp) :: arr_comp
                real(dp), dimension(0:npart) :: eta
                real(dp), dimension(ibeg:iend) :: delt_pos, delt_neg
                real(dp), dimension(ibeg:iend) :: fact_pos_b, fact_neg_b
                real(dp), dimension(ibeg:iend) :: fact_pos_e, fact_neg_e
            end subroutine rearrange_phideps
        end interface

        real(dp), allocatable :: phi(:), bhat(:), dlogbdphi(:), dbcovar(:)
        real(dp), allocatable :: delt_pos(:), delt_neg(:)
        real(dp), allocatable :: fact_pos_b(:), fact_neg_b(:)
        real(dp), allocatable :: fact_pos_e(:), fact_neg_e(:)
        real(dp), allocatable :: arr_real(:, :)
        complex(dp), allocatable :: arr_comp(:, :)
        real(dp) :: eta(0:npart), subsqmin
        integer :: i, iend

        iend = ub_mag
        if (padded) iend = ub_mag + 2

        allocate (phi(0:iend), bhat(0:iend), dlogbdphi(0:iend))
        allocate (dbcovar(0:iend), delt_pos(0:iend), delt_neg(0:iend))
        allocate (fact_pos_b(0:iend), fact_neg_b(0:iend))
        allocate (fact_pos_e(0:iend), fact_neg_e(0:iend))
        allocate (arr_real(0:iend, 1), arr_comp(0:iend, 1))

        do i = 0, ub_mag
            phi(i) = 0.05_dp*real(i, dp)
            bhat(i) = 1.0_dp + 0.5_dp*(phi(i) - 1.0_dp)**2
        end do
        if (bump) bhat(2) = 1.05_dp
        if (padded) then
            phi(ub_mag + 1:iend) = phi(ub_mag)
            bhat(ub_mag + 1:iend) = bhat_pad
        end if
        dlogbdphi = (phi - 1.0_dp)/bhat
        dbcovar = 0.0_dp
        arr_real(:, 1) = 1.0_dp
        arr_comp(:, 1) = cmplx(phi, 0.0_dp, dp)

        eta(0) = 0.0_dp
        eta(1) = 1.0_dp/bhat_cross
        subsqmin = 1.0e5_dp*epsilon(1.0_dp)
        divide = 2

        call rearrange_phideps(0, iend, ub_mag, npart, 1, 1, 0, subsqmin, &
            divide, phi, bhat, dlogbdphi, dbcovar, arr_real, arr_comp, eta, &
            delt_pos, delt_neg, fact_pos_b, fact_neg_b, fact_pos_e, &
            fact_neg_e)

        if (padded) then
            if (any(phi(ub_mag + 1:iend) /= phi(ub_mag))) &
                error stop 'FAIL: alignment moved a padding coordinate'
            if (any(bhat(ub_mag + 1:iend) /= bhat_pad)) &
                error stop 'FAIL: alignment changed a padding field value'
        end if
        if (maxval(divide) == 1) then
            if (any(delt_pos(1:iend) /= phi(1:iend) - phi(0:iend - 1))) &
                error stop 'FAIL: integration steps mismatch moved nodes'
            if (any(fact_pos_b(0:iend) /= 1.0_dp)) &
                error stop 'FAIL: step factors were modified'
        end if

        phi_physical = phi(0:ub_mag)
        bhat_physical = bhat(0:ub_mag)
    end subroutine run_case

end program test_crossing_geometry_padding
