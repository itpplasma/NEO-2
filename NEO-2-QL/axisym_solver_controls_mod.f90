! Opt-in AX iteration/preconditioner experiments, never physical operator changes.
module axisym_solver_controls_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private
    public :: read_axisym_solver_controls, regper_speed_selected, apply_ax_filters
    logical :: omit_filters = .false., all_speeds = .false.
    logical :: full_energy_projection = .false.
contains
    logical function alternative_selected(name, alternative) result(selected)
        character(len=*), intent(in) :: name, alternative
        character(len=32) :: value
        integer :: status
        call get_environment_variable(name, value, status=status)
        if (status == 1) value = 'native'
        if (status /= 0 .and. status /= 1) error stop 'AX policy read failed'
        selected = .false.
        select case (trim(value))
        case ('native')
        case default
            if (trim(value) /= alternative) error stop 'Invalid AX solver policy'
            selected = .true.
        end select
    end function alternative_selected

    subroutine read_axisym_solver_controls(lorentz, energy, integral)
        integer, intent(in) :: lorentz, energy, integral
        omit_filters = alternative_selected('NEO2_AX_FILTER_POLICY', 'none')
        all_speeds = alternative_selected('NEO2_AX_REGPER_SPEED', 'all')
        full_energy_projection = &
            alternative_selected('NEO2_AX_ENERGY_FILTER_POLICY', 'full')
        if (full_energy_projection) then
            if (any([lorentz, energy, integral] /= [0, 1, 1])) &
                error stop 'AX full energy projection requires full_fp'
            if (omit_filters) &
                error stop 'AX full energy projection conflicts with omitted filters'
            print *, 'native AX energy filter policy: full'
        end if
        if (all_speeds) then
            if (any([lorentz,energy,integral] /= 0)) &
                error stop 'AX all-speed regularizer requires effective pitch model'
        end if
        if (omit_filters) then
            print *, 'native AX filter policy: none'
        else
            print *, 'native AX filter policy: native'
        end if
        if (all_speeds) then
            print *, 'native AX regper speed: all'
        else
            print *, 'native AX regper speed: native'
        end if
    end subroutine read_axisym_solver_controls

    logical function regper_speed_selected(m, axisymmetric) result(selected)
        integer, intent(in) :: m
        logical, intent(in) :: axisymmetric
        selected = m < 1
        if (axisymmetric .and. all_speeds) selected = .true.
    end function regper_speed_selected

    subroutine apply_ax_filters(f, dbra, dket, ebra, eket, dd, ed, tol, &
                                energy_enabled, axisymmetric, density_stage, energy_stage)
        complex(dp), intent(inout) :: f(:)
        complex(dp), intent(in) :: dbra(:), dket(:), ebra(:), eket(:), dd, ed
        real(dp), intent(in) :: tol
        logical, intent(in) :: energy_enabled, axisymmetric
        complex(dp), optional, intent(out) :: density_stage(:), energy_stage(:)
        complex(dp) :: coef_dens, coef_energ
        if (present(density_stage)) density_stage = f
        if (present(energy_stage)) energy_stage = f
        ! NA has no such filter; an AX option must not change or abort its path.
        if (.not. axisymmetric) return
        if (omit_filters) return
        coef_dens = sum(dbra*f)/dd
        f = f - (1.0_dp - tol)*coef_dens*dket
        if (present(density_stage)) density_stage = f
        if (energy_enabled) then
            coef_energ = sum(ebra*f)/ed
            if (full_energy_projection) then
                f = f - coef_energ*eket
            else
                f = f - (1.0_dp - tol)*coef_energ*eket
            end if
        end if
        if (present(energy_stage)) energy_stage = f
    end subroutine apply_ax_filters
end module axisym_solver_controls_mod
