program test_fullfp_stage_probe
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fullfp_stage_probe_mod
    use axisym_solver_controls_mod, only: read_axisym_solver_controls, apply_ax_filters
    use mpiprovider_module, only: mpro
    use collisionality_mod, only: num_spec
    implicit none
    complex(dp) :: q(2), before(2), candidate(2), d(2), e(2), expected(2)
    complex(dp) :: basis(2,1), ritz(1), observed(2)
    character(len=1024) :: directory, path, line
    character(len=32) :: argument, rank_text, label, keyword
    real(dp) :: stats(8), scale, total_scale
    real(dp), allocatable :: physical(:,:)
    integer :: rank, unit, status, i, source, driver, iteration

    call get_command_argument(1, argument)
    call mpro%init()
    num_spec = mpro%getNumProcs()
    rank = mpro%getRank()
    if (trim(argument) == 'reject_pitch') then
        call fullfp_probe_configure(0, 0, 0, .true.)
        stop 0
    end if
    call fullfp_probe_configure(0, 1, 1, .true.)
    if (trim(argument) == 'off') then
        if (fullfp_probe_enabled) error stop 'unset observer enabled'
        call fullfp_probe_begin(2, 3, 1, 30)
        if (fullfp_probe_active) error stop 'disabled observer active'
        call mpro%deinit()
        stop
    end if
    if (.not. fullfp_probe_enabled) error stop 'observer not enabled'
    call get_environment_variable('NEO2_FULLFP_STAGE_PROBE_DIR', directory)
    write(rank_text, '(i0)') rank
    path = trim(directory)//'/rank'//trim(rank_text)
    call clean_outputs()
    call read_axisym_solver_controls(0, 1, 1)
    if (fullfp_probe_at_target(1, 1)) error stop 'old target guard retained'
    if (fullfp_probe_at_target(3, 0)) error stop 'wrong driver guard'
    if (.not. fullfp_probe_at_target(3, 1)) error stop 'missing target guard'
    call fullfp_probe_begin(2, 1, 1, 30)
    if (fullfp_probe_active) error stop 'wrong source observed'
    call fullfp_probe_begin(2, 3, 0, 30)
    if (fullfp_probe_active) error stop 'wrong driver observed'
    call fullfp_probe_begin(2, 3, 1, 29)
    if (fullfp_probe_active) error stop 'wrong iteration observed'
    call fullfp_probe_begin(2, 3, 1, 30)
    if (.not. fullfp_probe_active) error stop 'target observer not active'
    scale = real(rank + 1, dp)
    total_scale = real(num_spec*(num_spec + 1), dp)/2.0_dp
    q = [2.0_dp, 3.0_dp]
    before = scale
    fullfp_probe_inverse = scale*[0.25_dp, -0.125_dp]
    fullfp_probe_seed = scale*[1.0_dp, 2.0_dp]
    fullfp_probe_raw = before + fullfp_probe_inverse
    d = [1.0_dp, 0.0_dp]
    e = [0.0_dp, 1.0_dp]
    candidate = fullfp_probe_raw
    call apply_ax_filters(candidate, d, d, e, e, &
        cmplx(1.0_dp, 0.0_dp, dp), cmplx(1.0_dp, 0.0_dp, dp), &
        0.5_dp, .true., .true., fullfp_probe_density, fullfp_probe_energy)
    if (any(fullfp_probe_density /= scale*[0.625_dp, 0.875_dp])) &
        error stop 'actual density stage'
    if (any(fullfp_probe_energy /= scale*[0.625_dp, 0.4375_dp])) &
        error stop 'actual energy stage'
    candidate = candidate + scale*[0.125_dp, -0.125_dp]
    basis(:,1) = [1.0_dp, 0.0_dp]
    ritz = 0.5_dp
    call fullfp_probe_finish(q, before, candidate, ritz, basis, d, e, coupled_action)
    allocate(physical(6,0:num_spec-1))
    physical = 1.0_dp
    do i = 0, num_spec-1
        physical(2,i) = real(i + 2, dp)
    end do
    call fullfp_probe_moments(q + cmplx(0.0_dp,1.0_dp,dp)*[1.0_dp,-2.0_dp], &
        before, candidate, d, d, e, cmplx(1.0_dp,1.0_dp,dp)*e, &
        cmplx(1.0_dp,0.0_dp,dp), cmplx(1.0_dp,0.0_dp,dp), d, e, &
        physical, split_action)
    open(newunit=unit, file=trim(path)//'-energy-species0-M.bin', &
        access='stream', form='unformatted', status='old')
    read(unit) observed
    close(unit)
    expected = (0.0_dp,0.0_dp)
    if (rank == 0) expected = [(1.0_dp,1.0_dp),(3.0_dp,3.0_dp)]
    if (any(observed /= expected)) error stop 'individual energy M action'
    open(newunit=unit, file=trim(path)//'-energy-species0-C.bin', &
        access='stream', form='unformatted', status='old')
    read(unit) observed
    close(unit)
    if (any(observed /= [(0.25_dp,0.25_dp),(0.0_dp,0.0_dp)])) &
        error stop 'coupled energy C action'
    open(newunit=unit, file=trim(path)//'-energy-common-dT-R.bin', &
        access='stream', form='unformatted', status='old')
    read(unit) observed
    close(unit)
    expected = [(0.0_dp,0.0_dp),(0.125_dp,0.125_dp)]* &
        physical(2,min(1,num_spec-1))/physical(2,rank)
    if (maxval(abs(observed-expected)) > 1.0e-14_dp) &
        error stop 'physical common-dT scale and separate regularizer'
    open(newunit=unit, file=trim(path)//'-moments.txt', status='old')
    read(unit, '(a)') line
    read(unit, *) keyword, label, argument, source, driver, iteration
    if (source /= 3 .or. driver /= 1 .or. iteration /= 30) &
        error stop 'wrong moment target metadata'
    do
        read(unit, '(a)', iostat=status) line
        if (status /= 0) error stop 'missing complex source moments'
        if (index(line, 'moment source ') /= 1) cycle
        read(line, *) keyword, label, stats(1:6)
        if (any(abs(stats(1:6)-[2.0_dp,1.0_dp,3.0_dp,-2.0_dp,6.0_dp,-0.5_dp]) &
            > 1.0e-14_dp)) error stop 'complex native moment phase'
        exit
    end do
    close(unit)
    open(newunit=unit, file=trim(path)//'-energy-residual.bin', &
        access='stream', form='unformatted', status='old')
    read(unit) observed
    close(unit)
    ! Independent analytic coupled C across MPI ranks, M=[[2,1],[0,3]].
    expected = q - scale*[1.6875_dp, 1.3125_dp] + &
        total_scale*[0.109375_dp, -0.3125_dp]
    if (maxval(abs(observed - expected)) > 1.0e-14_dp) &
        error stop 'stored original coupled residual'
    open(newunit=unit, file=trim(path)//'-summary.txt', status='old')
    read(unit, '(a)') line
    read(unit, *) keyword, label, argument, source, driver, iteration
    if (source /= 3 .or. driver /= 1 .or. iteration /= 30) &
        error stop 'wrong stage target metadata'
    i = 0
    do
        read(unit, '(a)', iostat=status) line
        if (status /= 0) exit
        if (index(line, 'stage ') /= 1) cycle
        read(line, *) keyword, label, stats
        i = i + 1
        if (stats(8) > 32.0_dp) error stop 'linear residual identity'
        if (label == 'energy') then
            if (abs(stats(1) - sum(abs(expected))) > 1.0e-14_dp) &
                error stop 'residual numerator'
            if (abs(stats(4) - 1.0_dp) > 1.0e-14_dp) error stop 'density moment'
            if (abs(stats(5) - 0.7_dp) > 1.0e-14_dp) error stop 'energy moment'
        end if
    end do
    close(unit)
    if (i /= 4) error stop 'missing observed stages'
    call test_weak_rows()
    call clean_outputs()
    call mpro%deinit()
    print *, 'All tests passed!'
contains
    subroutine clean_outputs()
        character(len=12), parameter :: stages(4) = &
            [character(len=12) :: 'raw', 'density', 'energy', 'deflation']
        character(len=12), parameter :: fixed(6) = &
            [character(len=12) :: 'source', 'fold', 'inverse', 'seed', 'ritz', 'basis1']
        character(len=12), parameter :: suffix(4) = &
            [character(len=12) :: '', '-residual', '-delta', '-Ldelta']
        character(len=32) :: name
        integer :: a, b
        call delete_output('-moments.txt')
        call delete_output('-equation-number-bra.bin')
        call delete_output('-equation-energy-bra.bin')
        do a = 1, 4
            select case (a)
            case (1); name = 'density-bra'
            case (2); name = 'density-ket'
            case (3); name = 'energy-bra'
            case (4); name = 'energy-ket'
            end select
            call delete_output('-'//trim(name)//'.bin')
        end do
        do a = 0, num_spec
            if (a < num_spec) then
                write(name,'(a,i0)') 'energy-species', a
            else
                name = 'energy-common-dT'
            end if
            call delete_output('-'//trim(name)//'-M.bin')
            call delete_output('-'//trim(name)//'-C.bin')
            call delete_output('-'//trim(name)//'-R.bin')
        end do
        do a = 1, 3
            select case (a)
            case (1); name = 'density'
            case (2); name = 'energy'
            case (3); name = 'energy-delta'
            end select
            call delete_output('-'//trim(name)//'-M.bin')
            call delete_output('-'//trim(name)//'-C.bin')
        end do
        call delete_output('-summary.txt')
        do a = 1, size(fixed)
            call delete_output('-'//trim(fixed(a))//'.bin')
        end do
        do a = 1, size(stages)
            do b = 1, size(suffix)
                call delete_output('-'//trim(stages(a))//trim(suffix(b))//'.bin')
            end do
        end do
    end subroutine clean_outputs
    subroutine delete_output(suffix)
        character(len=*), intent(in) :: suffix
        integer :: file_unit
        logical :: exists
        inquire(file=trim(path)//suffix, exist=exists)
        if (.not. exists) return
        open(newunit=file_unit, file=trim(path)//suffix, status='old')
        close(file_unit, status='delete')
    end subroutine delete_output
    subroutine test_weak_rows()
        use collop_compute, only: collision_weak_moment_rows, phi_exp, lagmax, &
            make_ortho, isw_relativistic, integral_cutoff, lsw_split_interval, &
            M_transform_inv
        real(dp) :: rows(2,0:1)
        lagmax = 1
        phi_exp => monomial
        make_ortho = .true.
        isw_relativistic = 0
        integral_cutoff = .false.
        lsw_split_interval = .false.
        call collision_weak_moment_rows(rows)
        ! Independent Gaussian integrals: integral x^(2k+1) exp(-x*x)=k!/2.
        if (maxval(abs(rows(:,0)-[0.5_dp,1.0_dp])) > 1.0e-11_dp .or. &
            maxval(abs(rows(:,1)-[1.0_dp,3.0_dp])) > 1.0e-11_dp) &
            error stop 'native collision weak rows'
        allocate(M_transform_inv(0:1,0:1))
        M_transform_inv = 0.0_dp
        M_transform_inv(0,0) = 2.0_dp
        M_transform_inv(1,1) = 3.0_dp
        make_ortho = .false.
        call collision_weak_moment_rows(rows)
        if (maxval(abs(rows(:,0)-[1.0_dp,2.0_dp])) > 1.0e-11_dp .or. &
            maxval(abs(rows(:,1)-[3.0_dp,9.0_dp])) > 1.0e-11_dp) &
            error stop 'projected collision weak rows'
    end subroutine test_weak_rows
    function monomial(mode, x) result(value)
        integer :: mode
        real(dp) :: x, value
        value = x**(2*mode)
    end function monomial
    subroutine split_action(vector, matrix, collision, regularizer)
        complex(dp), intent(in) :: vector(:)
        complex(dp), intent(out) :: matrix(:), collision(:), regularizer(:)
        call coupled_action(vector, matrix, collision)
        regularizer = [0.5_dp*vector(1), 0.125_dp*vector(2)]
    end subroutine split_action
    subroutine coupled_action(vector, matrix, collision)
        complex(dp), intent(in) :: vector(:)
        complex(dp), intent(out) :: matrix(:), collision(:)
        complex(dp) :: gathered(2,0:num_spec-1)
        gathered(:,rank) = vector
        call mpro%allgather_inplace(gathered)
        matrix = [2.0_dp*vector(1) + vector(2), 3.0_dp*vector(2)]
        collision = [0.25_dp*sum(gathered(2,:)), -0.5_dp*sum(gathered(1,:))]
    end subroutine coupled_action
end program test_fullfp_stage_probe
