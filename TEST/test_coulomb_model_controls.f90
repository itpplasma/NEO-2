program test_coulomb_model_controls
    use collisionality_mod, only: collision_model, isw_lorentz, &
        isw_energy, isw_integral
    use rkstep_mod, only: lag
    use collision_model_mod, only: resolve_collision_model
    implicit none
    character(len=32) :: mode
    integer :: expected(3)

    call get_command_argument(1, mode)
    lag = 6
    isw_lorentz = 7
    isw_energy = 8
    isw_integral = 9
    collision_model = mode
    call resolve_collision_model()
    select case (trim(mode))
    case ('native')
        expected = [7, 8, 9]
    case ('pitch_angle')
        expected = [0, 0, 0]
    case ('test_particle')
        expected = [0, 1, 0]
    case ('full_fp')
        expected = [0, 1, 1]
    case ('monoenergetic_lorentz')
        expected = [1, 0, 0]
    case default
        error stop 'invalid selector was accepted'
    end select
    if (any([isw_lorentz, isw_energy, isw_integral] /= expected)) &
        error stop 'Coulomb switch mismatch'
    if (lag /= 6) error stop 'selector changed requested energy resolution'
end program test_coulomb_model_controls
