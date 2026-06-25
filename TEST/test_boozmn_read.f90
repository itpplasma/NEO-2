! Behavioral test for neo_read_boozmn (INP_SWI_BOOZMN = 10).
!
! Reads the committed fixture TEST/fixtures/boozmn_test.nc (generated once
! from NEO-RT/examples/circ.bc via bc_to_booz_xform.py at PR time).
! The fixture path is passed via the BOOZMN_TEST_FILE environment variable
! set by CMake; the test does not access /tmp or other repos at runtime.
!
! Tests:
!   1. Allocation and dimension sanity after a full-tolerance read.
!   2. m_max / n_max set from mboz_b / nboz_b (not 1 as before the fix).
!   3. pixm / pixn are valid indices into i_m / i_n.
!   4. B(0,0) matches the circ.bc reference to TOL_EXACT.
!   5. iota, curr_pol, curr_tor match committed fixture values to TOL_LOOSE.
!   6. Truncation regression: when max_m_mode is set below the data range
!      the reader must honour it and not clip it up to MAXVAL(|ixm|).
!
! circ.bc / boozmn_test.nc reference values (surface index 1, jlist[0]=2):
!   mboz_b = 18  -> m_max = 19
!   nboz_b = 0   -> n_max = 1  (axisymmetric)
!   iota(1)     = 0.9000140542
!   curr_pol(1) = bvco_b[jlist[0]-1] = -3.552357024  (poloidal covariant B)
!   curr_tor(1) = buco_b[jlist[0]-1] = -3.093771076e-3 (toroidal covariant B)
!   bmnc(1, m=0,n=0) = 1.96335464 T  (copied without interpolation)
program test_boozmn_read
  use nrtype
  use neo_input
  use neo_exchange, only : iota, curr_pol, curr_tor
  use neo_sub_mod, only : neo_read_boozmn
  use neo_control, only : in_file, inp_swi, INP_SWI_BOOZMN, &
    & max_m_mode, max_n_mode
  implicit none

  character(len=512) :: boozmn_path
  integer :: i, j, imn00, failures, trunc_m
  real(dp) :: b00_val, rel_err

  real(dp), parameter :: REF_B00   =  1.96335464_dp
  real(dp), parameter :: REF_IOTA  =  0.9000140542_dp
  real(dp), parameter :: REF_CPOL  = -3.552357024_dp
  real(dp), parameter :: REF_CTOR  = -3.093771076e-3_dp
  real(dp), parameter :: TOL_EXACT = 1.0e-6_dp
  real(dp), parameter :: TOL_LOOSE = 1.0e-4_dp

  failures = 0

  call get_environment_variable('BOOZMN_TEST_FILE', boozmn_path)
  if (len_trim(boozmn_path) == 0) then
    boozmn_path = 'boozmn_test.nc'
  end if

  ! -----------------------------------------------------------------------
  ! Pass 1: read with no truncation (max_m_mode well above data range).
  ! -----------------------------------------------------------------------
  inp_swi    = INP_SWI_BOOZMN
  max_m_mode = 9999
  max_n_mode = 9999
  in_file    = trim(boozmn_path)

  call neo_read_boozmn()

  ! allocation checks
  if (.not. allocated(bmnc)) then
    print *, 'FAIL: bmnc not allocated'; failures = failures + 1
  end if
  if (.not. allocated(ixm)) then
    print *, 'FAIL: ixm not allocated'; failures = failures + 1
  end if
  if (.not. allocated(i_m)) then
    print *, 'FAIL: i_m not allocated'; failures = failures + 1
  end if
  if (.not. allocated(i_n)) then
    print *, 'FAIL: i_n not allocated'; failures = failures + 1
  end if
  if (failures > 0) then
    print *, 'FAIL: allocation failures, cannot continue'
    error stop
  end if

  ! dimension sanity
  if (ns <= 0) then
    print *, 'FAIL: ns =', ns; failures = failures + 1
  end if
  if (mnmax <= 0) then
    print *, 'FAIL: mnmax =', mnmax; failures = failures + 1
  end if
  if (nfp <= 0) then
    print *, 'FAIL: nfp =', nfp; failures = failures + 1
  end if
  if (psi_pr <= 0.0_dp) then
    print *, 'FAIL: psi_pr =', psi_pr; failures = failures + 1
  end if

  ! m_max / n_max must come from mboz_b / nboz_b, not default to 1
  if (m_max /= 19) then
    print *, 'FAIL: m_max =', m_max, '(expected 19 from mboz_b=18)'
    failures = failures + 1
  end if
  if (n_max /= 1) then
    print *, 'FAIL: n_max =', n_max, '(expected 1 from nboz_b=0, axisymmetric)'
    failures = failures + 1
  end if

  ! pixm / pixn valid
  do j = 1, mnmax
    if (pixm(j) < 1 .or. pixm(j) > size(i_m)) then
      print *, 'FAIL: pixm(', j, ') =', pixm(j), 'out of range 1..', size(i_m)
      failures = failures + 1
      exit
    end if
    if (pixn(j) < 1 .or. pixn(j) > size(i_n)) then
      print *, 'FAIL: pixn(', j, ') =', pixn(j), 'out of range 1..', size(i_n)
      failures = failures + 1
      exit
    end if
  end do

  ! find (m=0, n=0) mode
  imn00 = 0
  do j = 1, mnmax
    if (ixm(j) == 0 .and. ixn(j) == 0) then
      imn00 = j
      exit
    end if
  end do
  if (imn00 == 0) then
    print *, 'FAIL: no (m=0,n=0) mode found'; failures = failures + 1
  end if

  ! B_00 from fixture (no interpolation in converter)
  if (imn00 > 0) then
    b00_val = bmnc(1, imn00)
    rel_err = abs(b00_val - REF_B00) / abs(REF_B00)
    if (rel_err > TOL_EXACT) then
      print *, 'FAIL: bmnc(1,m=0,n=0) =', b00_val
      print *, '      ref =', REF_B00, '  rel_err =', rel_err
      failures = failures + 1
    end if
  end if

  ! iota
  rel_err = abs(iota(1) - REF_IOTA) / abs(REF_IOTA)
  if (rel_err > TOL_LOOSE) then
    print *, 'FAIL: iota(1) =', iota(1), '  ref =', REF_IOTA
    failures = failures + 1
  end if

  ! curr_pol = bvco_b (poloidal covariant B)
  rel_err = abs(curr_pol(1) - REF_CPOL) / abs(REF_CPOL)
  if (rel_err > TOL_LOOSE) then
    print *, 'FAIL: curr_pol(1) =', curr_pol(1), '  ref =', REF_CPOL
    failures = failures + 1
  end if

  ! curr_tor = buco_b (toroidal covariant B)
  rel_err = abs(curr_tor(1) - REF_CTOR) / abs(REF_CTOR)
  if (rel_err > TOL_LOOSE) then
    print *, 'FAIL: curr_tor(1) =', curr_tor(1), '  ref =', REF_CTOR
    failures = failures + 1
  end if

  ! B(0,0) positive on all surfaces
  if (imn00 > 0) then
    do i = 1, ns
      if (bmnc(i, imn00) <= 0.0_dp) then
        print *, 'FAIL: bmnc(m=0,n=0) <= 0 on surface', i, ':', bmnc(i, imn00)
        failures = failures + 1
        exit
      end if
    end do
  end if

  ! -----------------------------------------------------------------------
  ! Pass 2: truncation regression.
  ! Set max_m_mode = 5, which is below MAXVAL(|ixm|) = 18.
  ! The reader must honour the user limit: after the call max_m_mode == 5.
  ! -----------------------------------------------------------------------
  trunc_m = 5

  ! Deallocate neo_input arrays from pass 1 before re-reading.
  if (allocated(bmnc))    deallocate(bmnc)
  if (allocated(rmnc))    deallocate(rmnc)
  if (allocated(zmnc))    deallocate(zmnc)
  if (allocated(lmnc))    deallocate(lmnc)
  if (allocated(ixm))     deallocate(ixm)
  if (allocated(ixn))     deallocate(ixn)
  if (allocated(pixm))    deallocate(pixm)
  if (allocated(pixn))    deallocate(pixn)
  if (allocated(i_m))     deallocate(i_m)
  if (allocated(i_n))     deallocate(i_n)
  if (allocated(es))      deallocate(es)
  if (allocated(iota))    deallocate(iota)
  if (allocated(curr_pol)) deallocate(curr_pol)
  if (allocated(curr_tor)) deallocate(curr_tor)
  if (allocated(pprime))  deallocate(pprime)
  if (allocated(sqrtg00)) deallocate(sqrtg00)
  if (allocated(b00))     deallocate(b00)

  inp_swi    = INP_SWI_BOOZMN
  max_m_mode = trunc_m
  max_n_mode = 9999
  in_file    = trim(boozmn_path)

  call neo_read_boozmn()

  ! max_m_mode must not have been clipped UP to MAXVAL(|ixm|)=18.
  if (max_m_mode /= trunc_m) then
    print *, 'FAIL: truncation not respected: max_m_mode =', max_m_mode, &
             '(expected', trunc_m, 'after MIN(5,18)=5)'
    failures = failures + 1
  end if

  if (failures == 0) then
    print *, 'All tests passed!'
    print *, '  pass 1: ns =', ns, '  mnmax =', mnmax, '  nfp =', nfp
    print *, '  m_max =', m_max, '  n_max =', n_max
    print *, '  psi_pr =', psi_pr
    print *, '  pass 2: max_m_mode after trunc =', max_m_mode
  else
    print *, 'FAIL'
    error stop
  end if

end program test_boozmn_read
