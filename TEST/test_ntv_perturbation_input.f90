program test_ntv_perturbation_input
  use ntv_mod, only: in_file_pert, isw_m_phi_input, m_phi_input, m_phi, &
                     has_perturbation_file
  use neo_magfie_perturbation, only: neo_read_pert
  implicit none

  integer :: test_status

  test_status = 0

  call test_none_file(test_status)
  call test_file_flag(test_status)

  if (test_status == 0) then
     print *, "All tests passed!"
  else
     print *, "Some tests failed. Status:", test_status
     error stop
  end if

contains

  subroutine test_none_file(status)
    integer, intent(inout) :: status

    in_file_pert = 'none'
    isw_m_phi_input = 1
    m_phi_input = -1
    m_phi = 0

    call neo_read_pert()

    if (has_perturbation_file()) then
       print *, "FAIL: none should disable perturbation-file input"
       status = status + 1
    else
       print *, "PASS: none disables perturbation-file input"
    end if

    if (m_phi /= -1) then
       print *, "FAIL: prescribed m_phi was not used"
       print *, "  Expected:", -1
       print *, "  Got:", m_phi
       status = status + 1
    else
       print *, "PASS: prescribed m_phi was used"
    end if
  end subroutine test_none_file

  subroutine test_file_flag(status)
    integer, intent(inout) :: status

    in_file_pert = 'perturbation.bc'

    if (has_perturbation_file()) then
       print *, "PASS: regular file name enables perturbation-file input"
    else
       print *, "FAIL: regular file name should enable perturbation-file input"
       status = status + 1
    end if
  end subroutine test_file_flag

end program test_ntv_perturbation_input
