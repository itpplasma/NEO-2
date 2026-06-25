! Behavioral test for neo_read_boozmn (INP_SWI_BOOZMN = 10).
!
! Reads a boozmn NetCDF generated from circ.bc (NEO-RT/examples/circ.bc)
! by bc_to_booz_xform.py and cross-checks that the reader produces values
! consistent with the file contents.
!
! Reference values are taken from the boozmn file itself (generated from circ.bc).
! The converter directly copies half-grid Fourier coefficients (bmnc_b), so
! bmnc(1, m=0,n=0) must match the .bc B_00 = 1.96335464 T to machine precision.
! Full-grid profiles (iota, buco_b, bvco_b) are interpolated by the converter;
! values read by neo_read_boozmn from the boozmn must equal the boozmn file values.
!
! After the curr_pol/curr_tor swap fix:
!   curr_pol(i) = bvco_b(jlist(i))  [poloidal current, T*m]
!   curr_tor(i) = buco_b(jlist(i))  [toroidal current, T*m]
!
! From circ.bc boozmn (jlist[1]=2, 0-based index 1 of 65-point full grid):
!   iota(1)    = 0.9000140542
!   curr_pol(1) = -3.552357024   (bvco_b[1])
!   curr_tor(1) = -0.003093771076  (buco_b[1])
!   bmnc(1,0)  = 1.96335464   (bmnc_b[0,0], directly from .bc)
!
! .bc reference for cross-check at half-grid surface 0 (s=0.023):
!   iota = 0.89997    -> boozmn interpolated 0.9000 differs < 0.05%
!   B_00 = 1.96335464 -> boozmn matches exactly (no interpolation for Fourier)
program test_boozmn_read
  use nrtype
  use neo_input
  use neo_exchange, only : iota, curr_pol, curr_tor
  use neo_sub_mod, only : neo_read_boozmn
  use neo_control, only : in_file, inp_swi, INP_SWI_BOOZMN, &
    & max_m_mode, max_n_mode
  implicit none

  character(len=200) :: boozmn_path
  integer :: i, j, imn00, failures
  real(dp) :: b00_val, rel_err

  ! Values directly from the circ.bc-derived boozmn at jlist[0]-1 = 1.
  ! bmnc_b is copied without interpolation: exact to machine epsilon.
  real(dp), parameter :: REF_B00    =  1.96335464_dp
  ! Full-grid profile references (interpolated by converter from .bc half-grid).
  ! iota: .bc surface 0 at s=0.023 has iota=0.89997; boozmn jlist[0] at s=0.0156
  ! gives 0.9000 via CubicSpline - within 0.05% of .bc.
  real(dp), parameter :: REF_IOTA   =  0.9000140542_dp
  ! curr_pol = bvco_b = (Jpol/nper)*nfp*mu0/(2*pi); .bc val * nfp * 2e-7 = -3.5508
  real(dp), parameter :: REF_CPOL   = -3.552357024_dp
  ! curr_tor = buco_b = Itor*mu0/(2*pi); .bc val * 2e-7 = -0.0047874
  real(dp), parameter :: REF_CTOR   = -3.093771076e-3_dp
  ! Tight tolerance for exact copies; loose for interpolated profiles.
  real(dp), parameter :: TOL_EXACT  = 1.0e-6_dp
  real(dp), parameter :: TOL_LOOSE  = 1.0e-4_dp

  failures = 0

  call get_environment_variable('BOOZMN_TEST_FILE', boozmn_path)
  if (len_trim(boozmn_path) == 0) then
    boozmn_path = 'boozmn_test.nc'
  end if

  inp_swi    = INP_SWI_BOOZMN
  max_m_mode = 9999
  max_n_mode = 9999
  in_file    = trim(boozmn_path)

  call neo_read_boozmn()

  ! --- allocation checks ---
  if (.not. allocated(bmnc)) then
    print *, 'FAIL: bmnc not allocated'
    failures = failures + 1
  end if
  if (.not. allocated(ixm)) then
    print *, 'FAIL: ixm not allocated'
    failures = failures + 1
  end if
  if (.not. allocated(i_m)) then
    print *, 'FAIL: i_m not allocated'
    failures = failures + 1
  end if
  if (.not. allocated(i_n)) then
    print *, 'FAIL: i_n not allocated'
    failures = failures + 1
  end if
  if (failures > 0) then
    print *, 'FAIL: allocation failures, cannot continue'
    error stop
  end if

  ! --- dimension sanity ---
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

  ! Bug 1: m_max / n_max were never set before fix.
  if (m_max <= 0) then
    print *, 'FAIL: m_max =', m_max, '(should be m0b+1 = 19)'
    failures = failures + 1
  end if
  if (n_max <= 0) then
    print *, 'FAIL: n_max =', n_max, '(should be 2*n0b+1 = 1 for axisymmetric)'
    failures = failures + 1
  end if

  ! --- pixm/pixn are valid indices into i_m/i_n ---
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

  ! --- find (m=0, n=0) mode index ---
  imn00 = 0
  do j = 1, mnmax
    if (ixm(j) == 0 .and. ixn(j) == 0) then
      imn00 = j
      exit
    end if
  end do
  if (imn00 == 0) then
    print *, 'FAIL: no (m=0,n=0) mode found'
    failures = failures + 1
  end if

  ! --- B_00 cross-check: converter copies bmnc_b directly, exact match ---
  if (imn00 > 0) then
    b00_val = bmnc(1, imn00)
    rel_err = abs(b00_val - REF_B00) / abs(REF_B00)
    if (rel_err > TOL_EXACT) then
      print *, 'FAIL: bmnc(1,m=0,n=0) =', b00_val
      print *, '      expected (from circ.bc via boozmn) =', REF_B00
      print *, '      rel. error =', rel_err, '> TOL_EXACT =', TOL_EXACT
      failures = failures + 1
    end if
  end if

  ! Bug 1 consequence: iota cross-check against boozmn value at jlist[0].
  rel_err = abs(iota(1) - REF_IOTA) / abs(REF_IOTA)
  if (rel_err > TOL_LOOSE) then
    print *, 'FAIL: iota(1) =', iota(1)
    print *, '      expected (boozmn jlist[0]) =', REF_IOTA
    print *, '      rel. error =', rel_err, '> TOL_LOOSE =', TOL_LOOSE
    failures = failures + 1
  end if

  ! Bug 2: curr_pol was buco_b, now bvco_b. Cross-check exact boozmn value.
  rel_err = abs(curr_pol(1) - REF_CPOL) / abs(REF_CPOL)
  if (rel_err > TOL_LOOSE) then
    print *, 'FAIL: curr_pol(1) =', curr_pol(1)
    print *, '      expected bvco_b (poloidal covariant B) =', REF_CPOL
    print *, '      rel. error =', rel_err, '> TOL_LOOSE =', TOL_LOOSE
    failures = failures + 1
  end if

  ! Bug 2: curr_tor was bvco_b, now buco_b. Cross-check exact boozmn value.
  rel_err = abs(curr_tor(1) - REF_CTOR) / abs(REF_CTOR)
  if (rel_err > TOL_LOOSE) then
    print *, 'FAIL: curr_tor(1) =', curr_tor(1)
    print *, '      expected buco_b (toroidal covariant B) =', REF_CTOR
    print *, '      rel. error =', rel_err, '> TOL_LOOSE =', TOL_LOOSE
    failures = failures + 1
  end if

  ! Bug 3: max_n_mode must not exceed what's in the data.
  if (max_n_mode > MAXVAL(ABS(ixn))) then
    print *, 'FAIL: max_n_mode =', max_n_mode, '> MAXVAL(|ixn|) =', MAXVAL(ABS(ixn))
    failures = failures + 1
  end if
  if (max_m_mode > MAXVAL(ABS(ixm))) then
    print *, 'FAIL: max_m_mode =', max_m_mode, '> MAXVAL(|ixm|) =', MAXVAL(ABS(ixm))
    failures = failures + 1
  end if

  ! --- B(0,0) positive on all surfaces ---
  if (imn00 > 0) then
    do i = 1, ns
      if (bmnc(i, imn00) <= 0.0_dp) then
        print *, 'FAIL: bmnc(m=0,n=0) <= 0 on surface', i, ':', bmnc(i, imn00)
        failures = failures + 1
        exit
      end if
    end do
  end if

  if (failures == 0) then
    print *, 'All tests passed!'
    print *, '  ns =', ns, '  mnmax =', mnmax, '  nfp =', nfp
    print *, '  m_max =', m_max, '  n_max =', n_max
    print *, '  m0b =', m0b, '  n0b =', n0b
    print *, '  psi_pr =', psi_pr
    print *, '  max_m_mode =', max_m_mode, '  max_n_mode =', max_n_mode
    print *, '  bmnc(1,m=0,n=0) =', bmnc(1, imn00), '  [ref =', REF_B00, ']'
    print *, '  iota(1) =', iota(1), '  [ref =', REF_IOTA, ']'
    print *, '  curr_pol(1) =', curr_pol(1), '  [ref bvco =', REF_CPOL, ']'
    print *, '  curr_tor(1) =', curr_tor(1), '  [ref buco =', REF_CTOR, ']'
  else
    print *, 'FAIL'
    error stop
  end if

end program test_boozmn_read
