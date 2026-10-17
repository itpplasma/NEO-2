module collision_model_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use collisionality_mod, only: collision_model, krook_nu, num_spec, &
        isw_lorentz, isw_energy, isw_integral, isw_momentum, isw_relativistic, &
        lsw_nbi
    implicit none
    private
    public :: resolve_collision_model, uses_krook, krook_carrier_scale
    public :: krook_weight_leaf, krook_pitch_mass
contains
    logical function uses_krook()
        uses_krook = trim(collision_model) == 'krook_constant' .or. &
            trim(collision_model) == 'krook_inverse_cube'
    end function uses_krook

    subroutine resolve_collision_model(read_precom, write_precom)
        logical, intent(in) :: read_precom, write_precom

        select case (trim(collision_model))
        case ('native')
            ! Preserve every existing switch and the legacy default exactly.
        case ('pitch_angle')
            isw_lorentz = 0
            isw_energy = 0
            isw_integral = 0
        case ('test_particle')
            isw_lorentz = 0
            isw_energy = 1
            isw_integral = 0
        case ('full_fp')
            isw_lorentz = 0
            isw_energy = 1
            isw_integral = 1
        case ('monoenergetic_lorentz')
            isw_lorentz = 1
            isw_energy = 0
            isw_integral = 0
        case ('krook_constant', 'krook_inverse_cube')
            if (isw_lorentz /= 0 .or. isw_energy /= 0 .or. isw_integral /= 0) &
                error stop 'Krook requires isw_lorentz=isw_energy=isw_integral=0'
            if (isw_momentum /= 0 .or. isw_relativistic /= 0) &
                error stop 'Krook requires nonrelativistic energy-basis representation'
            if (read_precom .or. write_precom) &
                error stop 'FP collision caches cannot be used with Krook'
            if (lsw_nbi) error stop 'Krook with NBI is not supported'
            if (.not. allocated(krook_nu)) error stop 'Missing Krook species rates'
            if (num_spec < 1 .or. num_spec > size(krook_nu)) &
                error stop 'Invalid Krook species count'
            if (.not. all(ieee_is_finite(krook_nu(1:num_spec)))) &
                error stop 'Krook rates must be finite'
            if (any(krook_nu(1:num_spec) <= 0.0_dp)) &
                error stop 'Krook requires a positive physical rate for every species'
        case default
            error stop 'Unknown collision_model'
        end select

        if (.not. uses_krook()) then
            if (allocated(krook_nu)) then
                if (any(krook_nu /= 0.0_dp)) &
                    error stop 'krook_nu is only valid for a Krook model'
            end if
        end if
    end subroutine resolve_collision_model

    function krook_carrier_scale(nu, kappa, temperature, mass) result(scale)
        real(dp), intent(in) :: nu, kappa, temperature, mass
        real(dp) :: scale

        if (.not. ieee_is_finite(nu) .or. .not. ieee_is_finite(kappa) .or. &
            .not. ieee_is_finite(temperature) .or. .not. ieee_is_finite(mass)) &
            error stop 'Nonfinite Krook normalization input'
        if (min(nu, kappa, temperature, mass) <= 0.0_dp) &
            error stop 'Nonpositive Krook normalization input'
        call krook_scale_leaf(nu, kappa, temperature, mass, scale)
        if (.not. ieee_is_finite(scale)) error stop 'Nonfinite Krook carrier scale'
        if (scale <= 0.0_dp) error stop 'Unrepresentable Krook carrier scale'
    end function krook_carrier_scale

    function krook_pitch_mass(kappa, h_phi, bhat, moments, reconstruction) &
            result(row)
        real(dp), intent(in) :: kappa, h_phi, bhat
        real(dp), intent(in) :: moments(:), reconstruction(:, :)
        real(dp) :: row(size(reconstruction, 2)), coefficient
        integer :: j

        ! Same signed weak pitch mass as the native energy carrier. This
        ! assembly carries a reaction, with no velocity derivative.
        call krook_pitch_leaf(kappa, h_phi, bhat, coefficient)
        do j = 1, size(row)
            row(j) = dot_product(moments,reconstruction(:,j))*coefficient
        end do
    end function krook_pitch_mass

    include 'generated/krook_leaf.inc'
end module collision_model_mod
