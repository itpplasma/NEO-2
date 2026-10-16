! Explicit, single-instance library lifecycle for the historic Boozer provider.
! Evaluations use NEO-2's existing splines; no field equations are duplicated.
module neo_field_provider
    use nrtype, only: dp
    implicit none
    private
    public :: neo_field_load, neo_field_prepare, neo_field_free
    public :: field_minor_radius_m, field_major_radius_m, field_s_min, field_s_max
    real(dp), protected :: field_minor_radius_m = 0.0_dp, field_major_radius_m = 0.0_dp
    real(dp), protected :: field_s_min = 0.0_dp, field_s_max = 0.0_dp
contains
    subroutine neo_field_load(path, input_switch)
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        use neo_control
        use neo_input, only: es
        use neo_sub_mod, only: neo_init
        use neo_spline_data, only: lsw_linear_boozer
        use neo_exchange, only: calc_cur
        use neo_magfie, only: magfie_spline, magfie_result, isw_mag_shear
        character(*), intent(in) :: path
        integer, intent(in) :: input_switch
        integer :: unit, npsi, j, i1, i2, i3, i4
        real(dp) :: flux
        character(1024) :: header

        if (input_switch /= 8 .and. input_switch /= 9) &
            error stop "neo_field_load supports Boozer formats 8 and 9"
        if (len_trim(path) > len(in_file)) error stop "Boozer input path too long"
        call neo_field_free()
        open(newunit=unit, file=path, status='old', action='read')
        do j = 1, 5
            read(unit, '(a)') header
        end do
        read(unit, *) i1, i2, i3, i4, flux, field_minor_radius_m, field_major_radius_m
        close(unit)
        if (.not. ieee_is_finite(field_minor_radius_m + field_major_radius_m)) &
            error stop "nonfinite Boozer file-header radii"
        if (field_minor_radius_m <= 0.0_dp .or. field_major_radius_m <= 0.0_dp) &
            error stop "neo_field_load requires positive file-header radii in metres"
        in_file = path
        inp_swi = input_switch
        ! Laboratory 10 determines the Jacobian/current sign from geometry
        ! (neo_magfie.f90), for either supported .bc column format.
        lab_swi = 10
        ! Existing cubic angular interpolation; 1025 points retain the
        ! consumer's independent 1e-8 curl-consistency gate for tested charts.
        theta_n = 1025
        phi_n = 9
        max_m_mode = huge(1)
        max_n_mode = huge(1)
        write_progress = 0
        write_output_files = 0
        calc_fourier = 0
        calc_cur = 0
        fluxs_interp = 1
        no_fluxs = 0
        ref_swi = 1
        lsw_linear_boozer = .false.
        magfie_spline = 1
        magfie_result = 0
        isw_mag_shear = 1
        call neo_init(npsi)
        field_s_min = es(1)
        field_s_max = es(size(es))
    end subroutine neo_field_load

    subroutine neo_field_prepare(s)
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        use neo_input, only: es
        use neo_magfie, only: magfie_sarray
        use magfie_mod, only: magfie_deallocate
        real(dp), intent(in) :: s
        if (.not. allocated(es)) error stop "neo_field_prepare before load"
        if (.not. ieee_is_finite(s)) error stop "nonfinite Boozer surface"
        if (s < es(1) .or. s > es(size(es))) error stop "Boozer surface outside data"
        call magfie_deallocate()
        if (allocated(magfie_sarray)) deallocate(magfie_sarray)
        allocate(magfie_sarray(1))
        magfie_sarray = s
    end subroutine neo_field_prepare

    subroutine neo_field_free()
        use neo_sub_mod, only: neo_dealloc
        use magfie_mod, only: magfie_deallocate
        use neo_magfie, only: magfie_sarray
        call magfie_deallocate()
        if (allocated(magfie_sarray)) deallocate(magfie_sarray)
        call neo_dealloc()
        field_minor_radius_m = 0.0_dp
        field_major_radius_m = 0.0_dp
        field_s_min = 0.0_dp
        field_s_max = 0.0_dp
    end subroutine neo_field_free
end module neo_field_provider
