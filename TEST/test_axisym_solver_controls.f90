program test_axisym_solver_controls
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use axisym_solver_controls_mod, only: read_axisym_solver_controls, &
        regper_speed_selected, apply_ax_filters
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
    select case (trim(mode))
    case ('reject_fp')
        call read_axisym_solver_controls(0,1,1)
    case ('reject_test_particle')
        call read_axisym_solver_controls(0,1,0)
    case ('reject_lorentz')
        call read_axisym_solver_controls(1,0,0)
    case ('none')
        call read_axisym_solver_controls(0,1,1)
    case default
        call read_axisym_solver_controls(0,0,0)
    end select
    ! WILL_FAIL must reject success if an invalid/misused option was accepted.
    if (trim(mode) == 'invalid' .or. index(mode,'reject_') == 1) stop 0
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
end program test_axisym_solver_controls
