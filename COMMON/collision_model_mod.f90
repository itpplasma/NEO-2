! Coulomb selector subset of owning PR214 commit55f183f3. No Krook assembly.
module collision_model_mod
    use collisionality_mod, only: collision_model, isw_lorentz, &
        isw_energy, isw_integral
    implicit none
    private
    public :: resolve_collision_model
contains
    subroutine resolve_collision_model()
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
        case default
            error stop 'Unknown collision_model'
        end select
    end subroutine resolve_collision_model
end module collision_model_mod
