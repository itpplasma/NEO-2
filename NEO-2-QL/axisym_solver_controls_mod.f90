! Opt-in AX iteration/preconditioner experiments, never physical operator changes.
module axisym_solver_controls_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private
    public :: read_axisym_solver_controls, regper_speed_selected, apply_ax_filters
    public :: configure_ax_feasibility, record_ax_feasibility, ax_feasibility_complete
    public :: ax_feasibility_enabled
    logical, protected :: ax_feasibility_enabled = .false.
    logical :: accepted_ax(3,0:2) = .false.
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

    subroutine configure_ax_feasibility(lorentz, energy, integral, species, &
                                        early_accept, tolerance, iterations, &
                                        arnoldi, probe)
        integer, intent(in) :: lorentz, energy, integral, species, iterations, arnoldi
        logical, intent(in) :: early_accept, probe
        real(dp), intent(in) :: tolerance
        character(len=8) :: value
        integer :: status
        ax_feasibility_enabled = .false.
        accepted_ax = .false.
        call get_environment_variable('NEO2_AX_ONLY_FEASIBILITY', value, status=status)
        if (status == 1) return
        if (status /= 0) error stop 'AX feasibility control read failed'
        select case (trim(value))
        case ('0')
            return
        case ('1')
        case default
            error stop 'AX feasibility control must be0or1'
        end select
        if (any([lorentz, energy, integral] /= [0, 1, 1]) .or. species /= 3) &
            error stop 'AX feasibility requires three-species full_fp'
        if (.not. full_energy_projection .or. omit_filters .or. all_speeds) &
            error stop 'AX feasibility requires full energy and native other policies'
        if (.not. early_accept .or. probe) &
            error stop 'AX feasibility requires original acceptance without probes'
        if (tolerance /= 1.0e-5_dp .or. iterations /= 100 .or. arnoldi /= 500) &
            error stop 'AX feasibility requires unchanged1e-5/100/500 controls'
        ax_feasibility_enabled = .true.
    end subroutine configure_ax_feasibility

    subroutine record_ax_feasibility(source, driver, status)
        integer, intent(in) :: source, driver, status
        if (.not. ax_feasibility_enabled) return
        if (source < 1 .or. source > 3 .or. driver < 0 .or. driver > 2) &
            error stop 'Invalid AX feasibility source/driver'
        if (accepted_ax(source,driver)) error stop 'Duplicate AX feasibility solve'
        ! Status0 is the existing original-equation/all-enabled-moment decision.
        if (status /= 0) error stop 'Unconverged AX feasibility solve'
        accepted_ax(source,driver) = .true.
    end subroutine record_ax_feasibility

    logical function ax_feasibility_complete() result(complete)
        complete = .false.
        if (.not. ax_feasibility_enabled) return
        if (.not. all(accepted_ax)) error stop 'Partial AX feasibility is not success'
        complete = .true.
    end function ax_feasibility_complete

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
