! Behavioral test for neo_read_boozmn (INP_SWI_BOOZMN = 10).
!
! Reads a boozmn NetCDF and checks:
!   - arrays are allocated and sized consistently
!   - (m=0, n=0) Fourier coefficient of |B| is positive on every surface
!   - iota is non-zero on every surface
!   - psi_pr > 0
!
! The test path is: neo_read_boozmn -> neo_input arrays.
! A cross-check against the .bc path requires the libneo converters
! (itpplasma/libneo#344, #345) which are not yet merged.  This test
! validates the reading path and physical sanity; the FP equality gate
! is deferred pending those PRs.
program test_boozmn_read
  use nrtype
  use neo_input
  use neo_exchange, only : iota, curr_pol, curr_tor
  use neo_sub_mod, only : neo_read_boozmn
  use neo_control, only : in_file, inp_swi, INP_SWI_BOOZMN, &
    & max_m_mode, max_n_mode
  implicit none

  integer :: i, j, imn00
  real(dp) :: b00_test
  integer :: failures

  failures = 0
  inp_swi = INP_SWI_BOOZMN
  max_m_mode = 9999
  max_n_mode = 9999
  in_file = 'boozmn_test.nc'

  call neo_read_boozmn()

  if (.not. allocated(bmnc)) then
    print *, 'FAIL: bmnc not allocated after neo_read_boozmn'
    failures = failures + 1
  end if
  if (.not. allocated(ixm)) then
    print *, 'FAIL: ixm not allocated after neo_read_boozmn'
    failures = failures + 1
  end if
  if (failures > 0) then
    print *, 'FAIL: allocation failures, cannot continue'
    error stop
  end if

  if (ns <= 0) then
    print *, 'FAIL: ns =', ns, '(expected > 0)'
    failures = failures + 1
  end if
  if (mnmax <= 0) then
    print *, 'FAIL: mnmax =', mnmax, '(expected > 0)'
    failures = failures + 1
  end if
  if (nfp <= 0) then
    print *, 'FAIL: nfp =', nfp, '(expected > 0)'
    failures = failures + 1
  end if
  if (psi_pr <= 0.0_dp) then
    print *, 'FAIL: psi_pr =', psi_pr, '(expected > 0)'
    failures = failures + 1
  end if

  ! Find (m=0, n=0) mode index
  imn00 = 0
  do j = 1, mnmax
    if (ixm(j) == 0 .and. ixn(j) == 0) then
      imn00 = j
      exit
    end if
  end do
  if (imn00 == 0) then
    print *, 'FAIL: no (m=0,n=0) mode found in ixm/ixn'
    failures = failures + 1
  end if

  if (imn00 > 0) then
    do i = 1, ns
      b00_test = bmnc(i, imn00)
      if (b00_test <= 0.0_dp) then
        print *, 'FAIL: bmnc(0,0) <= 0 on surface', i, ':', b00_test
        failures = failures + 1
        exit
      end if
    end do
  end if

  do i = 1, ns
    if (abs(iota(i)) < 1.0e-10_dp) then
      print *, 'FAIL: |iota| < 1e-10 on surface', i, ':', iota(i)
      failures = failures + 1
      exit
    end if
  end do

  if (failures == 0) then
    print *, 'All tests passed!'
    print *, '  ns =', ns, ' mnmax =', mnmax, ' nfp =', nfp
    print *, '  psi_pr =', psi_pr
    print *, '  bmnc(0,0) surface 1 =', bmnc(1, imn00)
    print *, '  iota surface 1 =', iota(1)
  else
    print *, 'FAIL'
    error stop
  end if

end program test_boozmn_read
