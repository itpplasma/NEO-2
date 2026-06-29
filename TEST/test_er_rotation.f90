program test_er_rotation
   use nrtype, only: dp
   use er_rotation_mod, only: Om_tE_to_MtOvR_spec, MtOvR_spec_to_Om_tE, &
                              check_Om_tE_consistency
   implicit none

   integer :: test_status

   test_status = 0

   call test_round_trip(test_status)
   call test_known_values(test_status)
   call test_multispecies_consistency(test_status)
   call test_consistency_check(test_status)
   call test_mode1_to_mode2_roundtrip(test_status)
   call test_half_omte_gives_half_mach(test_status)

   if (test_status == 0) then
      print *, "All tests passed!"
   else
      print *, "Some tests failed. Status:", test_status
      error stop
   end if

contains

   subroutine test_round_trip(status)
      integer, intent(inout) :: status
      real(dp) :: Om_tE_in, Om_tE_out
      real(dp) :: T_spec(2), m_spec(2), MtOvR_spec(2)
      real(dp) :: rel_err

      print *, "Testing round-trip conversion..."

      Om_tE_in = 1.0d3

      T_spec = [9.1495920740775243d-9, 7.1751185399075415d-9]
      m_spec = [9.1094d-28, 3.3436d-24]

      MtOvR_spec = Om_tE_to_MtOvR_spec(Om_tE_in, T_spec, m_spec)
      Om_tE_out = MtOvR_spec_to_Om_tE(MtOvR_spec, T_spec, m_spec)

      rel_err = abs(Om_tE_out - Om_tE_in) / abs(Om_tE_in)
      if (rel_err > epsilon(1.0_dp) * 10.0_dp) then
         print *, "FAIL: round-trip test"
         print *, "  Input Om_tE:", Om_tE_in
         print *, "  Output Om_tE:", Om_tE_out
         print *, "  Relative error:", rel_err
         status = status + 1
      else
         print *, "PASS: round-trip test"
      end if
   end subroutine test_round_trip

   subroutine test_known_values(status)
      integer, intent(inout) :: status
      real(dp) :: Om_tE, T_spec(1), m_spec(1), MtOvR_spec(1)
      real(dp) :: v_th, expected_MtOvR, rel_err

      print *, "Testing known values..."

      Om_tE = 1.0d3
      T_spec = [1.6d-9]
      m_spec = [1.672621637d-24]

      v_th = sqrt(2.0_dp * T_spec(1) / m_spec(1))
      expected_MtOvR = Om_tE / v_th

      MtOvR_spec = Om_tE_to_MtOvR_spec(Om_tE, T_spec, m_spec)

      rel_err = abs(MtOvR_spec(1) - expected_MtOvR) / abs(expected_MtOvR)
      if (rel_err > epsilon(1.0_dp) * 10.0_dp) then
         print *, "FAIL: known values test"
         print *, "  Expected MtOvR:", expected_MtOvR
         print *, "  Got MtOvR:", MtOvR_spec(1)
         print *, "  Relative error:", rel_err
         status = status + 1
      else
         print *, "PASS: known values test"
      end if

      if (abs(v_th - 4.3739733630d7) / 4.3739733630d7 > 1.0d-8) then
         print *, "FAIL: thermal velocity sanity check"
         print *, "  v_th:", v_th, " expected 4.3739733630e7 cm/s"
         status = status + 1
      else
         print *, "PASS: thermal velocity sanity check (v_th = 4.374e7 cm/s)"
      end if
   end subroutine test_known_values

   subroutine test_multispecies_consistency(status)
      integer, intent(inout) :: status
      real(dp) :: Om_tE_in
      real(dp) :: T_spec(2), m_spec(2), MtOvR_spec(2)
      real(dp) :: Om_tE_from_electron, Om_tE_from_proton, rel_err

      print *, "Testing multi-species consistency..."

      Om_tE_in = 5.0d4

      T_spec = [9.1495920740775243d-9, 7.1751185399075415d-9]
      m_spec = [9.1094d-28, 3.3436d-24]

      MtOvR_spec = Om_tE_to_MtOvR_spec(Om_tE_in, T_spec, m_spec)

      if (abs(MtOvR_spec(1) - MtOvR_spec(2)) < epsilon(1.0_dp)) then
         print *, "FAIL: MtOvR_spec should differ between electron and proton"
         status = status + 1
      else
         print *, "PASS: MtOvR_spec differs between species (electron/proton)"
      end if

      Om_tE_from_electron = MtOvR_spec(1) * sqrt(2.0_dp * T_spec(1) / m_spec(1))
      Om_tE_from_proton = MtOvR_spec(2) * sqrt(2.0_dp * T_spec(2) / m_spec(2))

      rel_err = abs(Om_tE_from_electron - Om_tE_from_proton) / abs(Om_tE_in)
      if (rel_err > epsilon(1.0_dp) * 100.0_dp) then
         print *, "FAIL: Om_tE should be species-independent"
         print *, "  From electron:", Om_tE_from_electron
         print *, "  From proton:", Om_tE_from_proton
         print *, "  Relative difference:", rel_err
         status = status + 1
      else
         print *, "PASS: Om_tE is species-independent"
      end if
   end subroutine test_multispecies_consistency

   subroutine test_consistency_check(status)
      integer, intent(inout) :: status
      real(dp) :: Om_tE
      real(dp) :: T_spec(2), m_spec(2), MtOvR_spec(2)
      logical :: is_consistent

      print *, "Testing consistency check..."

      Om_tE = 1.0d3
      T_spec = [9.1495920740775243d-9, 7.1751185399075415d-9]
      m_spec = [9.1094d-28, 3.3436d-24]

      MtOvR_spec = Om_tE_to_MtOvR_spec(Om_tE, T_spec, m_spec)
      is_consistent = check_Om_tE_consistency(MtOvR_spec, T_spec, m_spec, 1.0d-12)

      if (.not. is_consistent) then
         print *, "FAIL: consistent data should pass check"
         status = status + 1
      else
         print *, "PASS: consistent data passes check"
      end if

      MtOvR_spec(2) = MtOvR_spec(2) * 1.1_dp
      is_consistent = check_Om_tE_consistency(MtOvR_spec, T_spec, m_spec, 1.0d-12)

      if (is_consistent) then
         print *, "FAIL: inconsistent data should fail check"
         status = status + 1
      else
         print *, "PASS: inconsistent data fails check"
      end if
   end subroutine test_consistency_check

   subroutine test_mode1_to_mode2_roundtrip(status)
      ! Simulates the isw_calc_Er=1 -> isw_calc_Er=2 roundtrip:
      ! mode 1 computes Om_tE, derives MtOvR_spec;
      ! mode 2 takes the same Om_tE and must produce identical MtOvR_spec.
      integer, intent(inout) :: status
      real(dp) :: Om_tE_from_mode1
      real(dp) :: T_spec(2), m_spec(2)
      real(dp) :: MtOvR_mode1(2), MtOvR_mode2(2)
      real(dp) :: rel_err
      integer :: i

      print *, "Testing mode 1 -> mode 2 roundtrip..."

      T_spec = [9.1495920740775243d-9, 7.1751185399075415d-9]
      m_spec = [9.1094d-28, 3.3436d-24]

      Om_tE_from_mode1 = 4.237d4

      MtOvR_mode1 = Om_tE_to_MtOvR_spec(Om_tE_from_mode1, T_spec, m_spec)
      MtOvR_mode2 = Om_tE_to_MtOvR_spec(Om_tE_from_mode1, T_spec, m_spec)

      do i = 1, 2
         if (abs(MtOvR_mode1(i)) > 0.0_dp) then
            rel_err = abs(MtOvR_mode2(i) - MtOvR_mode1(i))/abs(MtOvR_mode1(i))
         else
            rel_err = abs(MtOvR_mode2(i))
         end if
         if (rel_err > epsilon(1.0_dp)) then
            print *, "FAIL: mode1->mode2 roundtrip, species", i
            print *, "  Mode 1 MtOvR:", MtOvR_mode1(i)
            print *, "  Mode 2 MtOvR:", MtOvR_mode2(i)
            status = status + 1
            return
         end if
      end do
      print *, "PASS: mode1->mode2 roundtrip (MtOvR_spec identical)"
   end subroutine test_mode1_to_mode2_roundtrip

   subroutine test_half_omte_gives_half_mach(status)
      ! Verify that halving Om_tE halves each species MtOvR_spec
      ! (MtOvR_spec is linear in Om_tE).
      integer, intent(inout) :: status
      real(dp) :: Om_tE_full, Om_tE_half
      real(dp) :: T_spec(2), m_spec(2)
      real(dp) :: MtOvR_full(2), MtOvR_half(2)
      real(dp) :: ratio, rel_err
      integer :: i

      print *, "Testing half Om_tE gives half MtOvR..."

      T_spec = [9.1495920740775243d-9, 7.1751185399075415d-9]
      m_spec = [9.1094d-28, 3.3436d-24]

      Om_tE_full = 4.237d4
      Om_tE_half = Om_tE_full*0.5_dp

      MtOvR_full = Om_tE_to_MtOvR_spec(Om_tE_full, T_spec, m_spec)
      MtOvR_half = Om_tE_to_MtOvR_spec(Om_tE_half, T_spec, m_spec)

      do i = 1, 2
         ratio = MtOvR_half(i)/MtOvR_full(i)
         rel_err = abs(ratio - 0.5_dp)/0.5_dp
         if (rel_err > epsilon(1.0_dp)*10.0_dp) then
            print *, "FAIL: half Om_tE test, species", i
            print *, "  Full MtOvR:", MtOvR_full(i)
            print *, "  Half MtOvR:", MtOvR_half(i)
            print *, "  Ratio:", ratio, " expected 0.5"
            status = status + 1
            return
         end if
      end do

      if (abs(MtOvR_full(1) - MtOvR_half(1)) < epsilon(1.0_dp)) then
         print *, "FAIL: full and half Om_tE should produce different MtOvR"
         status = status + 1
         return
      end if

      print *, "PASS: half Om_tE gives half MtOvR (linear scaling)"
   end subroutine test_half_omte_gives_half_mach

end program test_er_rotation
