program test_axisym_solver_controls
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use axisym_solver_controls_mod, only: read_axisym_solver_controls, &
        regper_speed_selected, apply_ax_filters, configure_ax_feasibility, &
        record_ax_feasibility, ax_feasibility_complete, ax_feasibility_enabled
    use sparse_mod, only: sparse_solve, sparse_solve_method, sparse_coo_action
    implicit none
    character(len=32) :: mode
    complex(dp) :: f(2), old(2), d(2), e(2), one, dd, ed, q(2), expected(2)
    complex(dp) :: dbra(2), dket(2), ebra(2), eket(2), original_action(8)
    real(dp) :: a(2,2), residual, m(8,8), r(8,8), k(8,8)
    real(dp) :: exact(8), rhs(8)
    logical :: no_filters, all_speeds
    integer :: speed, i, j
    call get_command_argument(1,mode)
    if (index(mode,'feas_') == 1) then
        call feasibility_guard_oracle(trim(mode))
        print *, 'All tests passed!'
        stop
    end if
    select case (trim(mode))
    case ('reject_fp')
        call read_axisym_solver_controls(0,1,1)
    case ('reject_test_particle')
        call read_axisym_solver_controls(0,1,0)
    case ('reject_lorentz')
        call read_axisym_solver_controls(1,0,0)
    case ('none','energy_full','energy_invalid','energy_conflict')
        call read_axisym_solver_controls(0,1,1)
    case ('energy_reject_pitch')
        call read_axisym_solver_controls(0,0,0)
    case ('energy_reject_test_particle')
        call read_axisym_solver_controls(0,1,0)
    case ('energy_reject_lorentz')
        call read_axisym_solver_controls(1,0,0)
    case default
        call read_axisym_solver_controls(0,0,0)
    end select
    ! WILL_FAIL must reject success if an invalid/misused option was accepted.
    if (trim(mode) == 'invalid' .or. index(mode,'reject_') == 1) stop 0
    if (index(mode,'energy_') == 1 .and. trim(mode) /= 'energy_full') stop 0
    if (trim(mode) == 'energy_full') then
        call full_energy_oracles()
        print *, 'All tests passed!'
        stop
    end if
    no_filters = trim(mode) == 'none'
    all_speeds = trim(mode) == 'all'
    one = (1.0_dp,0.0_dp)
    d = [one,(0.0_dp,0.0_dp)]; e = [d(2),d(1)]
    ! Actual helper versus legacy sequential density/energy arithmetic, complex input.
    dbra = d; dket = [(1.0_dp,0.0_dp),(0.25_dp,0.0_dp)]
    ebra = e
    eket = [(0.5_dp,0.0_dp),(1.0_dp,0.0_dp)]
    dd = sum(dbra*dket); ed = sum(ebra*eket)
    f = [(1.0_dp,0.5_dp),(-0.5_dp,0.25_dp)]; old = f; expected = f
    if (.not. no_filters) then
        expected = expected - 0.5_dp*(sum(dbra*expected)/dd)*dket
        expected = expected - 0.5_dp*(sum(ebra*expected)/ed)*eket
    end if
    call apply_ax_filters(f,dbra,dket,ebra,eket,dd,ed,0.5_dp,.true.,.true.)
    if (any(f /= expected)) error stop 'legacy/default filter parity'
    f = old
    call apply_ax_filters(f,dbra,dket,ebra,eket,dd,ed,0.5_dp,.true.,.false.)
    if (any(f /= old)) error stop 'AX filter option changed NA'
    ! Analytic nonnormal iteration: T=I-A is nilpotent, but F*T has eigenvalue1.
    a = reshape([3.0_dp,2.0_dp,-2.0_dp,-1.0_dp],[2,2])
    q = [(-2.0_dp,0.0_dp),(-1.0_dp,0.0_dp)]; f = (0.0_dp,0.0_dp)
    do i = 1, 3
        f = f+q-matmul(a,f)
        call apply_ax_filters(f,d,d,e,e,one,one,0.5_dp,.false.,.true.)
    end do
    residual = sum(abs(q-matmul(a,f)))/(sum(abs(q))+sum(abs(matmul(a,f))))
    if (no_filters) then
        if (any(f /= e) .or. residual /= 0.0_dp .or. abs(f(1)) /= 0.0_dp) &
            error stop 'unfiltered original equation and moment'
    else
        if (any(f /= (-1.0_dp,0.0_dp)) .or. abs(residual-0.2_dp) > 1e-14_dp) &
            error stop 'analytic filtered fixed point bias'
    end if
    ! One non-null plus seven speed blocks; the actual guard selects R support.
    m = 0.0_dp; m(1,1) = 1.0_dp; r = 0.0_dp
    do speed = 0, 6
        if (regper_speed_selected(speed,.true.)) r(speed+2,speed+2) = 1.0_dp
        if (regper_speed_selected(speed,.true.) .neqv. &
            (speed == 0 .or. all_speeds)) error stop 'AX speed support'
        if (regper_speed_selected(speed,.false.) .neqv. (speed == 0)) &
            error stop 'AX regularizer option changed NA'
    end do
    k = m+r
    if (all_speeds) then
        do j = 1, 8
            do i = 1, 8
                if (k(i,j) /= merge(1.0_dp,0.0_dp,i == j)) &
                    error stop 'all-speed analytic inverse matrix'
            end do
        end do
        exact = [(real(i,dp),i=1,8)]; rhs = matmul(k,exact)
        sparse_solve_method = 2
        call sparse_solve(k,rhs)
        if (maxval(abs(rhs-exact)) > 1e-12_dp) error stop 'all-speed inverse'
        if (maxval(abs(matmul(k,rhs)-exact)) > 1e-12_dp) &
            error stop 'independent inverse action'
        call sparse_coo_action([1],[1],[1.0_dp],cmplx(rhs,0.0_dp,dp), &
            original_action)
        if (any(original_action /= cmplx(matmul(m,exact),0.0_dp,dp))) &
            error stop 'unregularized operator action after inverse'
    else
        if (k(2,2) /= 1.0_dp .or. any(k(3:8,3:8) /= 0.0_dp)) &
            error stop 'legacy speed regularizer parity'
    end if
    print *, 'All tests passed!'
contains
    subroutine feasibility_guard_oracle(selection)
        character(len=*), intent(in) :: selection
        integer :: source, driver
        call read_axisym_solver_controls(0,1,1)
        if (selection == 'feas_default' .or. selection == 'feas_off') then
            ! Unset control must leave the default route untouched, even when
            ! the deliberately incompatible diagnostic-only arguments are given.
            call configure_ax_feasibility(1,0,0,1,.false.,1.0_dp,1,1,.true.)
            if (ax_feasibility_enabled .or. ax_feasibility_complete()) &
                error stop 'Default AX feasibility changed normal execution'
            call record_ax_feasibility(-1,-1,1)
            return
        end if
        call configure_ax_feasibility(0,1,1,3,.true.,1.0e-5_dp,100,500, &
            selection == 'feas_incompatible')
        if (selection == 'feas_invalid' .or. &
            selection == 'feas_incompatible') stop 0
        if (.not. ax_feasibility_enabled) error stop 'AX feasibility not enabled'
        do source = 1, 3
            do driver = 0, 2
                if (selection == 'feas_partial' .and. &
                    source == 3 .and. driver == 2) cycle
                if (selection == 'feas_failed') then
                    call record_ax_feasibility(source,driver,1)
                    stop 0
                end if
                call record_ax_feasibility(source,driver,0)
                if (selection == 'feas_duplicate') then
                    call record_ax_feasibility(source,driver,0)
                    stop 0
                end if
            end do
        end do
        if (selection == 'feas_partial') then
            if (ax_feasibility_complete()) stop 0
            stop 0
        end if
        if (.not. ax_feasibility_complete()) error stop 'Nine AX reports incomplete'
    end subroutine feasibility_guard_oracle

    subroutine full_energy_oracles()
        complex(dp) :: v(2), initial(2), bd(2), kd(2), be(2), ke(2)
        complex(dp) :: dn, en
        ! Independent normalized Maxwellian moments: <1>=1,
        ! <x^2-3/2>=0 and <(x^2-3/2)^2>=3/2.
        bd = [(1.0_dp,0.0_dp),(0.0_dp,0.0_dp)]
        kd = bd
        be = [(0.0_dp,0.0_dp),(1.5_dp,0.0_dp)]
        ke = [(0.0_dp,0.0_dp),(1.0_dp,0.0_dp)]
        dn = (1.0_dp,0.0_dp)
        en = (1.5_dp,0.0_dp)
        initial = [(2.0_dp,1.0_dp),(-3.0_dp,4.0_dp)]
        v = initial
        call apply_ax_filters(v,bd,kd,be,ke,dn,en,0.5_dp,.true.,.true.)
        if (any(v /= [(1.0_dp,0.5_dp),(0.0_dp,0.0_dp)])) &
            error stop 'full energy complex Maxwellian oracle'
        ! Actual helper NA branch must preserve the complex vector exactly.
        v = initial
        call apply_ax_filters(v,bd,kd,be,ke,dn,en,0.5_dp,.true.,.false.)
        if (any(v /= initial)) error stop 'full energy policy changed NA'
        ! Nonorthogonal rows: full energy removal changes density too. Literal
        ! rational result checks ordering and complex phase without reusing code.
        bd = [(1.0_dp,0.0_dp),(0.0_dp,0.0_dp)]
        kd = [(1.0_dp,0.0_dp),(0.25_dp,0.0_dp)]
        be = [(0.0_dp,0.0_dp),(1.0_dp,0.0_dp)]
        ke = [(0.5_dp,0.0_dp),(1.0_dp,0.0_dp)]
        initial = [(1.0_dp,0.5_dp),(-0.5_dp,0.25_dp)]
        v = initial
        call apply_ax_filters(v,bd,kd,be,ke,dn,dn,0.5_dp,.true.,.true.)
        if (any(v /= [(0.8125_dp,0.15625_dp),(0.0_dp,0.0_dp)])) &
            error stop 'full energy nonorthogonal literal oracle'
        if (sum(be*v) /= (0.0_dp,0.0_dp)) error stop 'full energy moment'
        if (sum(bd*v) == (0.0_dp,0.0_dp)) &
            error stop 'density must remain a separate enabled constraint'
    end subroutine full_energy_oracles
end program test_axisym_solver_controls
