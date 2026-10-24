module fullfp_stage_probe_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use mpiprovider_module, only: mpro
    use collisionality_mod, only: num_spec
    implicit none
    private
    public :: fullfp_probe_configure, fullfp_probe_begin, fullfp_probe_finish
    public :: fullfp_probe_at_target
    public :: fullfp_probe_enabled, fullfp_probe_active, fullfp_probe_seed
    public :: fullfp_probe_raw, fullfp_probe_density, fullfp_probe_energy
    public :: fullfp_probe_inverse, stage_identity_units, fullfp_probe_moments
    logical :: fullfp_probe_enabled = .false., fullfp_probe_active = .false.
    complex(dp), allocatable :: fullfp_probe_seed(:), fullfp_probe_raw(:)
    complex(dp), allocatable :: fullfp_probe_density(:), fullfp_probe_energy(:)
    complex(dp), allocatable :: fullfp_probe_inverse(:)
    character(len=1024) :: output_directory = ''
contains
    subroutine fullfp_probe_configure(lorentz, energy, integral, early_accept)
        integer, intent(in) :: lorentz, energy, integral
        logical, intent(in) :: early_accept
        character(len=32) :: policy
        integer :: status, i
        character(len=24), parameter :: controls(2) = &
            [character(len=24) :: 'NEO2_AX_FILTER_POLICY', 'NEO2_AX_REGPER_SPEED']

        fullfp_probe_enabled = .false.
        fullfp_probe_active = .false.
        call get_environment_variable('NEO2_FULLFP_STAGE_PROBE_DIR', &
            output_directory, status=status)
        if (status == 1) return
        if (status /= 0) error stop 'Full-FP stage directory read failed'
        if (len_trim(output_directory) == 0) error stop 'Empty full-FP stage directory'
        if (output_directory(1:1) /= '/') error stop 'Stage directory must be absolute'
        if (any([lorentz, energy, integral] /= [0, 1, 1])) &
            error stop 'Full-FP stage observer requires full_fp'
        if (.not. early_accept) error stop 'Stage observer requires original acceptance'
        do i = 1, size(controls)
            call get_environment_variable(trim(controls(i)), policy, status=status)
            if (status == 1) cycle
            if (status /= 0) error stop 'Full-FP policy read failed'
            if (trim(policy) /= 'native') &
                error stop 'Stage observer requires native policy'
        end do
        fullfp_probe_enabled = .true.
    end subroutine fullfp_probe_configure

    pure logical function fullfp_probe_at_target(source, driver) result(at_target)
        integer, intent(in) :: source, driver
        at_target = source == 3 .and. driver == 1
    end function fullfp_probe_at_target

    subroutine fullfp_probe_begin(n, source, driver, iteration)
        integer, intent(in) :: n, source, driver, iteration
        fullfp_probe_active = .false.
        if (.not. fullfp_probe_enabled) return
        if (.not. fullfp_probe_at_target(source, driver)) return
        if (iteration /= 30) return
        if (allocated(fullfp_probe_raw)) error stop 'Stage observer repeated'
        allocate(fullfp_probe_raw(n), fullfp_probe_density(n), &
            fullfp_probe_energy(n), fullfp_probe_inverse(n))
        fullfp_probe_active = .true.
    end subroutine fullfp_probe_begin

    ! Called only after the actual fixed candidate and original stage observation.
    ! Bras are bilinear native rows: sum(bra*vector), without conjugation.
    subroutine fullfp_probe_moments(q, before, candidate, dbra, dket, ebra, eket, &
            dd, ed, equation_number, equation_energy, physical, action)
        complex(dp), intent(in) :: q(:), before(:), candidate(:)
        complex(dp), intent(in) :: dbra(:), dket(:), ebra(:), eket(:), dd, ed
        complex(dp), intent(in) :: equation_number(:), equation_energy(:)
        real(dp), intent(in) :: physical(6,0:num_spec-1)
        interface
            subroutine action(vector, matrix, collision, regularizer)
                import :: dp
                complex(dp), intent(in) :: vector(:)
                complex(dp), intent(out) :: matrix(:), collision(:), regularizer(:)
            end subroutine action
        end interface
        complex(dp), allocatable :: vector(:), matrix(:), collision(:), regularizer(:)
        character(len=1200) :: prefix
        character(len=32) :: rank_text, label
        integer :: rank, unit, ios, direction, stage

        if (.not. fullfp_probe_active) error stop 'Inactive full-FP moment observation'
        rank = mpro%getRank()
        write(rank_text, '(i0)') rank
        prefix = trim(output_directory)//'/rank'//trim(rank_text)
        open(newunit=unit, file=trim(prefix)//'-moments.txt', status='new', &
            action='write', iostat=ios)
        if (ios /= 0) error stop 'Cannot create full-FP moment summary'
        write(unit, '(a)') 'NEO2_FULLFP_MOMENT_V1 DIAGNOSTIC_ONLY NO_TORQUE'
        write(unit, '(a,3(1x,i0))') 'source driving_species iteration', 3, 1, 30
        write(unit, '(a,i0)') 'rank ', rank
        write(unit, '(a)') 'physical rows: n_cm-3 T_erg m_g Z kappa_cm-1 vT_cm_s'
        do direction = 0, num_spec-1
            write(unit, '(a,i0,6(1x,es24.16e3))') &
                'physical ', direction, physical(:,direction)
        end do
        write(unit, '(a,4(1x,es24.16e3))') 'filter_denominators ', &
            real(dd,dp), aimag(dd), real(ed,dp), aimag(ed)
        write(unit, '(a)') 'moment rows: label density_re density_im '// &
            'centered_energy_re centered_energy_im '// &
            'uncentered_energy_re uncentered_energy_im'
        call vector_file('density-bra', dbra)
        call vector_file('density-ket', dket)
        call vector_file('energy-bra', ebra)
        call vector_file('energy-ket', eket)
        call vector_file('equation-number-bra', equation_number)
        call vector_file('equation-energy-bra', equation_energy)
        call moments('source', q)
        call moments('fold', before)
        call moments('raw', fullfp_probe_raw)
        call moments('density', fullfp_probe_density)
        call moments('energy', fullfp_probe_energy)
        call moments('deflation', candidate)
        call moments('density-ket', dket)
        call moments('energy-ket', eket)
        allocate(vector(size(q)), matrix(size(q)), collision(size(q)), &
            regularizer(size(q)))
        ! Separately saved terms permit the native denominator for fixed beta readback.
        do stage = 1, 3
            select case (stage)
            case (1)
                label = 'density'
                vector = fullfp_probe_density
            case (2)
                label = 'energy'
                vector = fullfp_probe_energy
            case (3)
                label = 'energy-delta'
                vector = fullfp_probe_energy - fullfp_probe_density
            end select
            call action(vector, matrix, collision, regularizer)
            call vector_file(trim(label)//'-M', matrix)
            call vector_file(trim(label)//'-C', collision)
            call action_moments(trim(label))
        end do
        ! Individual consumed energy directions followed by a prescribed common-dT
        ! direction, scaled by T(reference rank 1)/T(a). No nullspace is assumed.
        do direction = 0, num_spec
            vector = (0.0_dp, 0.0_dp)
            if (direction < num_spec) then
                if (rank == direction) vector = eket
                write(label, '(a,i0)') 'energy-species', direction
            else
                vector = eket*physical(2,min(1,num_spec-1))/physical(2,rank)
                label = 'energy-common-dT'
            end if
            call moments(trim(label)//'-input', vector)
            call action(vector, matrix, collision, regularizer)
            call vector_file(trim(label)//'-M', matrix)
            call vector_file(trim(label)//'-C', collision)
            call vector_file(trim(label)//'-R', regularizer)
            call action_moments(trim(label))
        end do
        close(unit)
    contains
        subroutine moments(name, v)
            character(len=*), intent(in) :: name
            complex(dp), intent(in) :: v(:)
            complex(dp) :: density, energy
            density = sum(dbra*v)
            energy = sum(ebra*v)
            write(unit, '(a,a,6(1x,es24.16e3))') 'moment ', name, &
                real(density,dp), aimag(density), real(energy,dp), aimag(energy), &
                real(energy + 1.5_dp*density,dp), aimag(energy + 1.5_dp*density)
            density = sum(equation_number*v)
            energy = sum(equation_energy*v)
            write(unit, '(a,a,4(1x,es24.16e3))') 'equation_moment ', name, &
                real(density,dp), aimag(density), real(energy,dp), aimag(energy)
        end subroutine moments
        subroutine action_moments(name)
            character(len=*), intent(in) :: name
            ! These are the actual filter functionals on equation vectors. They
            ! are deliberately not labelled physical collision conservation rows.
            call moments(name//'-M', matrix)
            call moments(name//'-C', collision)
            call moments(name//'-L', matrix - collision)
            call moments(name//'-R', regularizer)
        end subroutine action_moments
        subroutine vector_file(name, v)
            character(len=*), intent(in) :: name
            complex(dp), intent(in) :: v(:)
            integer :: output, status
            open(newunit=output, file=trim(prefix)//'-'//name//'.bin', &
                access='stream', form='unformatted', status='new', &
                action='write', iostat=status)
            if (status /= 0) error stop 'Cannot create full-FP moment vector'
            write(output, iostat=status) v
            if (status /= 0) error stop 'Cannot write full-FP moment vector'
            close(output)
        end subroutine vector_file
    end subroutine fullfp_probe_moments

    real(dp) function stage_identity_units(residual, previous, delta_m, delta_c, &
            scale) result(units)
        complex(dp), intent(in) :: residual(:), previous(:), delta_m(:), delta_c(:)
        real(dp), intent(in) :: scale
        units = sum(abs(residual - previous + delta_m - delta_c)) / &
            max(epsilon(1.0_dp)*scale, tiny(1.0_dp))
    end function stage_identity_units

    subroutine fullfp_probe_finish(q, before, candidate, ritz, basis, dbra, ebra, &
            action)
        complex(dp), intent(in) :: q(:), before(:), candidate(:), ritz(:), basis(:,:)
        complex(dp), intent(in) :: dbra(:), ebra(:)
        interface
            subroutine action(vector, matrix, collision)
                import :: dp
                complex(dp), intent(in) :: vector(:)
                complex(dp), intent(out) :: matrix(:), collision(:)
            end subroutine action
        end interface
        complex(dp), allocatable :: stage(:), previous(:), residual(:), old_residual(:)
        complex(dp), allocatable :: matrix(:), collision(:), delta(:), dm(:), dc(:)
        real(dp) :: stats(8,0:num_spec-1), old_scale, scale, qnorm
        character(len=12), parameter :: names(4) = &
            [character(len=12) :: 'raw', 'density', 'energy', 'deflation']
        character(len=32) :: rank_text
        character(len=1200) :: prefix
        integer :: n, rank, i, unit, ios
        integer(int64) :: dump_bytes

        if (.not. fullfp_probe_active) error stop 'Inactive full-FP observation'
        n = size(q)
        dump_bytes = 16_int64*n*num_spec* &
            (35_int64 + 3_int64*num_spec + size(ritz)) + 2_int64*1024**2
        if (dump_bytes > 2_int64*1024**3) error stop 'Stage dump exceeds 2 GiB'
        if (.not. allocated(fullfp_probe_seed)) &
            error stop 'Missing original Arnoldi seed'
        if (size(fullfp_probe_seed) /= n) error stop 'Stage seed size mismatch'
        rank = mpro%getRank()
        write(rank_text, '(i0)') rank
        prefix = trim(output_directory)//'/rank'//trim(rank_text)
        open(newunit=unit, file=trim(prefix)//'-summary.txt', status='new', &
            action='write', iostat=ios)
        if (ios /= 0) error stop 'Cannot create full-FP stage summary'
        write(unit, '(a)') 'NEO2_FULLFP_STAGE_V1 DIAGNOSTIC_ONLY NO_TORQUE'
        write(unit, '(a,3(1x,i0))') 'source driving_species iteration', 3, 1, 30
        write(unit, '(a,3(1x,i0))') 'vector_size ranks modes', n, num_spec, size(ritz)
        write(unit, '(a)') 'vectors: native complex(real64), stream, Fortran order'
        write(unit, '(a)') 'identity_units is a floating-point indicator, not a bound'
        write(unit, '(a)') 'metrics: residual_L1 denominator normalized density energy '// &
            'delta_L1 Ldelta_L1 identity_units'
        call write_vector('source', q)
        call write_vector('fold', before)
        call write_vector('inverse', fullfp_probe_inverse)
        call write_vector('seed', fullfp_probe_seed)
        call write_vector('ritz', ritz)
        do i = 1, size(ritz)
            write(rank_text, '(i0)') i
            call write_vector('basis'//trim(rank_text), basis(:,i))
        end do
        allocate(stage(n), previous(n), residual(n), old_residual(n), &
            matrix(n), collision(n), delta(n), dm(n), dc(n))
        previous = before
        call action(previous, matrix, collision)
        old_residual = q + collision - matrix
        qnorm = sum(abs(q))
        old_scale = qnorm + sum(abs(matrix)) + sum(abs(collision))
        do i = 1, size(names)
            select case (i)
            case (1); stage = fullfp_probe_raw
            case (2); stage = fullfp_probe_density
            case (3); stage = fullfp_probe_energy
            case (4); stage = candidate
            end select
            delta = stage - previous
            call action(stage, matrix, collision)
            residual = q + collision - matrix
            scale = qnorm + sum(abs(matrix)) + sum(abs(collision))
            call action(delta, dm, dc)
            stats(:,rank) = [sum(abs(residual)), scale, &
                sum(abs(residual))/max(scale, tiny(1.0_dp)), &
                abs(sum(dbra*stage))/max(sum(abs(dbra))*maxval(abs(stage)), &
                    tiny(1.0_dp)), &
                abs(sum(ebra*stage))/max(sum(abs(ebra))*maxval(abs(stage)), &
                    tiny(1.0_dp)), &
                sum(abs(delta)), sum(abs(dm - dc)), &
                stage_identity_units(residual, old_residual, dm, dc, &
                    scale + old_scale + sum(abs(dm)) + sum(abs(dc)))]
            call mpro%allgather_inplace(stats)
            if (.not. all(ieee_is_finite(stats))) error stop 'Nonfinite stage metrics'
            write(unit, '(a,1x,a,8(1x,es24.16e3))') &
                'stage', trim(names(i)), stats(:,rank)
            if (rank == 0) then
                print *, 'FULLFP_STAGE ', trim(names(i)), &
                    ' species residual ', stats(3,:)
                print *, 'FULLFP_STAGE ', trim(names(i)), ' identity_units ', stats(8,:)
            end if
            call write_vector(trim(names(i)), stage)
            call write_vector(trim(names(i))//'-residual', residual)
            call write_vector(trim(names(i))//'-delta', delta)
            call write_vector(trim(names(i))//'-Ldelta', dm - dc)
            previous = stage
            old_residual = residual
            old_scale = scale
        end do
        close(unit)
    contains
        subroutine write_vector(label, vector)
            character(len=*), intent(in) :: label
            complex(dp), intent(in) :: vector(:)
            integer :: output, status
            open(newunit=output, file=trim(prefix)//'-'//label//'.bin', &
                access='stream', form='unformatted', status='new', &
                action='write', iostat=status)
            if (status /= 0) error stop 'Cannot create full-FP stage vector'
            write(output, iostat=status) vector
            if (status /= 0) error stop 'Cannot write full-FP stage vector'
            close(output)
        end subroutine write_vector
    end subroutine fullfp_probe_finish
end module fullfp_stage_probe_mod
