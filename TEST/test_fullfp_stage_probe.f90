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
    integer :: rank, unit, status, i

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
        call fullfp_probe_begin(2, 1, 1, 10)
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
    call fullfp_probe_begin(2, 1, 0, 10)
    if (fullfp_probe_active) error stop 'wrong driver observed'
    call fullfp_probe_begin(2, 1, 1, 9)
    if (fullfp_probe_active) error stop 'wrong iteration observed'
    call fullfp_probe_begin(2, 1, 1, 10)
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
        integer :: a, b
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
