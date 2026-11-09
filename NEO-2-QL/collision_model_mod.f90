!> One switch for the collision operator of NEO-2-QL, collision_model in
!> &collision, in place of the three legacy flags isw_lorentz, isw_energy and
!> isw_integral (still accepted alone, for existing decks):
!>   full_fp                full linearised Fokker-Planck operator      (0, 1, 1)
!>   test_particle          pitch-angle and energy scattering, no
!>                          field-particle part                         (0, 1, 0)
!>   pitch_angle            pitch-angle scattering only                 (0, 0, 0)
!>   monoenergetic_lorentz  the monoenergetic mode of isw_lorentz = 1   (1, 0, 0)
!> as (isw_lorentz, isw_energy, isw_integral). Note that isw_lorentz = 1 is
!> not pitch-angle scattering of a multispecies plasma but the monoenergetic
!> mode. Giving collision_model together with any legacy flag is an error.
module collision_model_mod
   implicit none
   private

   public :: resolve_collision_model, COLLISION_FLAG_UNSET

   integer, parameter :: COLLISION_FLAG_UNSET = -1

contains

   !> The three legacy flags from collision_model, or their legacy defaults
   !> (1, 0, 0) where unset when no model is given. message is empty on
   !> success.
   subroutine resolve_collision_model(model, isw_lorentz, isw_energy, isw_integral, message)
      character(*), intent(in) :: model
      integer, intent(inout) :: isw_lorentz, isw_energy, isw_integral
      character(len=*), intent(out) :: message
      logical :: legacy
      message = ''
      legacy = isw_lorentz /= COLLISION_FLAG_UNSET .or. isw_energy /= COLLISION_FLAG_UNSET &
         .or. isw_integral /= COLLISION_FLAG_UNSET
      if (len_trim(model) == 0) then
         if (isw_lorentz == COLLISION_FLAG_UNSET) isw_lorentz = 1
         if (isw_energy == COLLISION_FLAG_UNSET) isw_energy = 0
         if (isw_integral == COLLISION_FLAG_UNSET) isw_integral = 0
         return
      end if
      if (legacy) then
         message = 'collision_model replaces isw_lorentz, isw_energy and isw_integral; give one or the other'
         return
      end if
      select case (trim(model))
      case ('full_fp')
         call set(0, 1, 1)
      case ('test_particle')
         call set(0, 1, 0)
      case ('pitch_angle')
         call set(0, 0, 0)
      case ('monoenergetic_lorentz')
         call set(1, 0, 0)
      case default
         message = 'unknown collision_model "'//trim(model)// &
            '"; use full_fp, test_particle, pitch_angle or monoenergetic_lorentz'
      end select
   contains
      subroutine set(l, e, i)
         integer, intent(in) :: l, e, i
         isw_lorentz = l
         isw_energy = e
         isw_integral = i
      end subroutine set
   end subroutine resolve_collision_model

end module collision_model_mod
