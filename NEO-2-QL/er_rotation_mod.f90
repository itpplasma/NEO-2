module er_rotation_mod
   use nrtype, only: dp
   implicit none
   private

   public :: Om_tE_to_MtOvR_spec, MtOvR_spec_to_Om_tE, check_Om_tE_consistency

contains

   pure function Om_tE_to_MtOvR_spec(Om_tE, T_spec, m_spec) result(MtOvR_spec)
      real(dp), intent(in) :: Om_tE
      real(dp), intent(in) :: T_spec(:)
      real(dp), intent(in) :: m_spec(:)
      real(dp) :: MtOvR_spec(size(T_spec))

      MtOvR_spec = Om_tE / sqrt(2.0_dp * T_spec / m_spec)
   end function Om_tE_to_MtOvR_spec

   pure function MtOvR_spec_to_Om_tE(MtOvR_spec, T_spec, m_spec) result(Om_tE)
      real(dp), intent(in) :: MtOvR_spec(:)
      real(dp), intent(in) :: T_spec(:)
      real(dp), intent(in) :: m_spec(:)
      real(dp) :: Om_tE

      if (size(MtOvR_spec) < 1) then
         Om_tE = 0.0_dp
         return
      end if
      Om_tE = MtOvR_spec(1) * sqrt(2.0_dp * T_spec(1) / m_spec(1))
   end function MtOvR_spec_to_Om_tE

   pure function check_Om_tE_consistency(MtOvR_spec, T_spec, m_spec, tol) &
      result(is_consistent)
      real(dp), intent(in) :: MtOvR_spec(:)
      real(dp), intent(in) :: T_spec(:)
      real(dp), intent(in) :: m_spec(:)
      real(dp), intent(in) :: tol
      logical :: is_consistent

      real(dp) :: Om_tE_ref, Om_tE_i
      integer :: i

      is_consistent = .true.
      if (size(MtOvR_spec) < 2) return
      if (size(T_spec) /= size(MtOvR_spec) .or. &
          size(m_spec) /= size(MtOvR_spec)) then
         is_consistent = .false.
         return
      end if

      Om_tE_ref = MtOvR_spec(1) * sqrt(2.0_dp * T_spec(1) / m_spec(1))

      do i = 2, size(MtOvR_spec)
         Om_tE_i = MtOvR_spec(i) * sqrt(2.0_dp * T_spec(i) / m_spec(i))
         if (abs(Om_tE_ref) > 0.0_dp) then
            if (abs(Om_tE_i - Om_tE_ref) / abs(Om_tE_ref) > tol) then
               is_consistent = .false.
               return
            end if
         else
            if (abs(Om_tE_i) > tol) then
               is_consistent = .false.
               return
            end if
         end if
      end do
   end function check_Om_tE_consistency

end module er_rotation_mod
