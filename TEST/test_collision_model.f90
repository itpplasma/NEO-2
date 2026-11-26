!> collision_model of NEO-2-QL: each model sets the legacy flags, the legacy
!> flags alone keep their meaning and defaults, and both together are refused.
program test_collision_model
   use collision_model_mod, only: resolve_collision_model, COLLISION_FLAG_UNSET
   implicit none
   integer :: failures = 0

   call expect('full_fp', COLLISION_FLAG_UNSET, COLLISION_FLAG_UNSET, COLLISION_FLAG_UNSET, 0, 1, 1)
   call expect('test_particle', COLLISION_FLAG_UNSET, COLLISION_FLAG_UNSET, COLLISION_FLAG_UNSET, 0, 1, 0)
   call expect('pitch_angle', COLLISION_FLAG_UNSET, COLLISION_FLAG_UNSET, COLLISION_FLAG_UNSET, 0, 0, 0)
   call expect('monoenergetic_lorentz', COLLISION_FLAG_UNSET, COLLISION_FLAG_UNSET, &
      COLLISION_FLAG_UNSET, 1, 0, 0)
   ! Legacy decks: defaults (1, 0, 0) where unset, given flags unchanged.
   call expect('', COLLISION_FLAG_UNSET, COLLISION_FLAG_UNSET, COLLISION_FLAG_UNSET, 1, 0, 0)
   call expect('', 0, 1, 1, 0, 1, 1)
   call expect('', 0, COLLISION_FLAG_UNSET, COLLISION_FLAG_UNSET, 0, 0, 0)
   call refuse('pitch_angle', 0, COLLISION_FLAG_UNSET, COLLISION_FLAG_UNSET)
   call refuse('krook', COLLISION_FLAG_UNSET, COLLISION_FLAG_UNSET, COLLISION_FLAG_UNSET)
   if (failures > 0) then
      print '(a)', 'FAIL'
      error stop 1
   end if
   print '(a)', 'All tests passed!'

contains

   subroutine expect(model, l0, e0, i0, l, e, i)
      character(*), intent(in) :: model
      integer, intent(in) :: l0, e0, i0, l, e, i
      integer :: ll, ee, ii
      character(len=200) :: message
      ll = l0; ee = e0; ii = i0
      call resolve_collision_model(model, ll, ee, ii, message)
      if (len_trim(message) /= 0 .or. ll /= l .or. ee /= e .or. ii /= i) then
         failures = failures + 1
         print '(a,a,3i3,1x,a)', 'FAIL ', model, ll, ee, ii, trim(message)
      end if
   end subroutine expect

   subroutine refuse(model, l0, e0, i0)
      character(*), intent(in) :: model
      integer, intent(in) :: l0, e0, i0
      integer :: ll, ee, ii
      character(len=200) :: message
      ll = l0; ee = e0; ii = i0
      call resolve_collision_model(model, ll, ee, ii, message)
      if (len_trim(message) == 0) then
         failures = failures + 1
         print '(a,a)', 'FAIL not refused: ', model
      end if
   end subroutine refuse
end program test_collision_model
