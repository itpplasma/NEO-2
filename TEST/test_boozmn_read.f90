! Behavioral test for neo_read_boozmn (INP_SWI_BOOZMN = 10).
!
! Reads the committed fixture TEST/fixtures/boozmn_test.nc (axisymmetric,
! lasym=0, nboz_b=0, nfp=1, mboz_b=18, 63 computed surfaces).
! The fixture path is passed via the BOOZMN_TEST_FILE environment variable
! set by CMake; the test does not access /tmp or other repos at runtime.
!
! Tests:
!   1. Allocation and dimension sanity after a full-tolerance read.
!   2. m_max / n_max set from mboz_b / nboz_b (not 1 as before the fix).
!   3. pixm / pixn are valid indices into i_m / i_n.
!   4. B(0,0) matches the circ.bc reference to TOL_EXACT.
!   5. iota, curr_pol, curr_tor match committed fixture values to TOL_FP.
!   6. Truncation regression: when max_m_mode is set below the data range
!      the reader must honour it and not clip it up to MAXVAL(|ixm|).
!   7. Asymmetric booz_xform modes are imported with NEO-2's INP_SWI_TOK
!      sign convention: file ixn_b=+1 becomes ixn=-1 before evaluation.
!
! circ.bc / boozmn_test.nc reference values (surface index 1, jlist[0]=3):
!   mboz_b = 18  -> m_max = 19
!   nboz_b = 0   -> n_max = 1  (axisymmetric)
!   iota(1)     = 0.89997   (bc.iota[0], written directly at jlist(1)=3)
!   curr_pol(1) = -3.5508   (Jpol/nper * nfp * mu0/(2pi); written at jlist(1))
!   curr_tor(1) = -4.7874e-3 (Itor * mu0/(2pi); written at jlist(1))
!   bmnc(1, m=0,n=0) = 1.96335464 T  (copied without interpolation)
!
! Fixture regenerated with libneo bc_to_booz_xform after PR #347 radial-grid
! fix; full-grid values at jlist(i) are written from .bc arrays directly, so
! the reader recovers them to floating-point precision (TOL_FP = 1e-6).
program test_boozmn_read
  use nrtype
  use neo_input
  use neo_exchange, only : iota, curr_pol, curr_tor
  use neo_sub_mod, only : neo_read_boozmn
  use neo_control, only : in_file, inp_swi, INP_SWI_BOOZMN, &
    & INP_SWI_TOK, max_m_mode, max_n_mode
  use netcdf, only : nf90_clobber, nf90_close, nf90_create, nf90_def_dim, &
    & nf90_def_var, nf90_double, nf90_enddef, nf90_int, nf90_noerr, &
    & nf90_put_var, nf90_strerror
  implicit none

  character(len=512) :: boozmn_path
  integer :: i, j, imn00, failures, trunc_m
  real(dp) :: b00_val, rel_err

  real(dp), parameter :: REF_B00   =  1.96335464_dp
  real(dp), parameter :: REF_IOTA  =  0.89997_dp
  real(dp), parameter :: REF_CPOL  = -3.5508_dp
  real(dp), parameter :: REF_CTOR  = -4.7874e-3_dp
  real(dp), parameter :: TOL_EXACT = 1.0e-6_dp
  real(dp), parameter :: TOL_FP    = 1.0e-6_dp
  character(len=*), parameter :: ASYM_BOOZMN_FILE = 'boozmn_asym_sign.nc'
  real(dp), parameter :: ASYM_THETA = 0.3_dp
  real(dp), parameter :: ASYM_PHI   = 0.7_dp
  real(dp), parameter :: ASYM_B00   = 1.25_dp
  real(dp), parameter :: ASYM_BC11  = 0.4_dp
  real(dp), parameter :: ASYM_BS11  = -0.2_dp

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
  if (rel_err > TOL_FP) then
    print *, 'FAIL: iota(1) =', iota(1), '  ref =', REF_IOTA
    failures = failures + 1
  end if

  ! curr_pol = bvco_b (poloidal covariant B)
  rel_err = abs(curr_pol(1) - REF_CPOL) / abs(REF_CPOL)
  if (rel_err > TOL_FP) then
    print *, 'FAIL: curr_pol(1) =', curr_pol(1), '  ref =', REF_CPOL
    failures = failures + 1
  end if

  ! curr_tor = buco_b (toroidal covariant B)
  rel_err = abs(curr_tor(1) - REF_CTOR) / abs(REF_CTOR)
  if (rel_err > TOL_FP) then
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

  ! -----------------------------------------------------------------------
  ! Pass 3: asymmetric booz_xform sign conversion.
  ! -----------------------------------------------------------------------
  call reset_neo_input()
  call write_asymmetric_boozmn_fixture(ASYM_BOOZMN_FILE)

  inp_swi    = INP_SWI_BOOZMN
  max_m_mode = 9999
  max_n_mode = 9999
  in_file    = ASYM_BOOZMN_FILE

  call neo_read_boozmn()
  call check_asymmetric_sign_conversion(failures)

  if (failures == 0) then
    print *, 'All tests passed!'
    print *, '  pass 1: ns =', ns, '  mnmax =', mnmax, '  nfp =', nfp
    print *, '  m_max =', m_max, '  n_max =', n_max
    print *, '  psi_pr =', psi_pr
    print *, '  pass 2: max_m_mode after trunc =', max_m_mode
    print *, '  pass 3: asymmetric ixn sign conversion verified'
  else
    print *, 'FAIL'
    error stop
  end if

contains

  subroutine reset_neo_input()
    if (allocated(bmnc))     deallocate(bmnc)
    if (allocated(rmnc))     deallocate(rmnc)
    if (allocated(zmnc))     deallocate(zmnc)
    if (allocated(lmnc))     deallocate(lmnc)
    if (allocated(bmns))     deallocate(bmns)
    if (allocated(rmns))     deallocate(rmns)
    if (allocated(zmns))     deallocate(zmns)
    if (allocated(lmns))     deallocate(lmns)
    if (allocated(ixm))      deallocate(ixm)
    if (allocated(ixn))      deallocate(ixn)
    if (allocated(pixm))     deallocate(pixm)
    if (allocated(pixn))     deallocate(pixn)
    if (allocated(i_m))      deallocate(i_m)
    if (allocated(i_n))      deallocate(i_n)
    if (allocated(es))       deallocate(es)
    if (allocated(iota))     deallocate(iota)
    if (allocated(curr_pol)) deallocate(curr_pol)
    if (allocated(curr_tor)) deallocate(curr_tor)
    if (allocated(pprime))   deallocate(pprime)
    if (allocated(sqrtg00))  deallocate(sqrtg00)
    if (allocated(b00))      deallocate(b00)
  end subroutine reset_neo_input

  subroutine check_asymmetric_sign_conversion(failures)
    integer, intent(inout) :: failures

    integer :: j, imn11
    real(dp) :: actual_b, expected_b, angle

    if (inp_swi /= INP_SWI_TOK) then
      print *, 'FAIL: asymmetric boozmn did not switch to INP_SWI_TOK'
      failures = failures + 1
    end if
    if (.not. allocated(bmns)) then
      print *, 'FAIL: asymmetric bmns not allocated'
      failures = failures + 1
      return
    end if

    imn11 = 0
    actual_b = 0.0_dp
    do j = 1, mnmax
      if (ixm(j) == 1 .and. ixn(j) == -1) imn11 = j
      angle = ixm(j)*ASYM_THETA + ixn(j)*ASYM_PHI
      actual_b = actual_b + bmnc(1,j)*cos(angle) + bmns(1,j)*sin(angle)
    end do

    if (imn11 == 0) then
      print *, 'FAIL: no converted (m=1, ixn=-1) asymmetric mode found'
      failures = failures + 1
    end if

    angle = ASYM_THETA - ASYM_PHI
    expected_b = ASYM_B00 + ASYM_BC11*cos(angle) + ASYM_BS11*sin(angle)
    if (abs(actual_b - expected_b) > 1.0e-12_dp) then
      print *, 'FAIL: asymmetric B sign convention mismatch'
      print *, '      actual =', actual_b, ' expected =', expected_b
      failures = failures + 1
    end if
  end subroutine check_asymmetric_sign_conversion

  subroutine write_asymmetric_boozmn_fixture(path)
    character(len=*), intent(in) :: path

    integer, parameter :: ns_full = 3
    integer, parameter :: nsurf_half = 1
    integer, parameter :: mn_modes = 2
    integer :: ncid, dim_radius, dim_mode, dim_surf
    integer :: var_ns, var_nfp, var_mboz, var_nboz, var_mnboz, var_lasym
    integer :: var_jlist, var_ixm, var_ixn
    integer :: var_iota, var_buco, var_bvco, var_phi
    integer :: var_bmnc, var_rmnc, var_zmns, var_pmns
    integer :: var_bmns, var_rmns, var_zmnc, var_pmnc
    integer, dimension(nsurf_half) :: jlist_file
    integer, dimension(mn_modes) :: ixm_file, ixn_file
    real(dp), dimension(ns_full) :: iota_file, buco_file, bvco_file, phi_file
    real(dp), dimension(mn_modes, nsurf_half) :: coeff

    call check_nc(nf90_create(path, nf90_clobber, ncid), 'create boozmn')
    call check_nc(nf90_def_dim(ncid, 'radius', ns_full, dim_radius), &
      & 'define radius dimension')
    call check_nc(nf90_def_dim(ncid, 'mn_mode', mn_modes, dim_mode), &
      & 'define mode dimension')
    call check_nc(nf90_def_dim(ncid, 'comput_surfs', nsurf_half, dim_surf), &
      & 'define surface dimension')

    call check_nc(nf90_def_var(ncid, 'ns_b', nf90_int, varid=var_ns), &
      & 'define ns_b')
    call check_nc(nf90_def_var(ncid, 'nfp_b', nf90_int, varid=var_nfp), &
      & 'define nfp_b')
    call check_nc(nf90_def_var(ncid, 'mboz_b', nf90_int, varid=var_mboz), &
      & 'define mboz_b')
    call check_nc(nf90_def_var(ncid, 'nboz_b', nf90_int, varid=var_nboz), &
      & 'define nboz_b')
    call check_nc(nf90_def_var(ncid, 'mnboz_b', nf90_int, varid=var_mnboz), &
      & 'define mnboz_b')
    call check_nc(nf90_def_var(ncid, 'lasym__logical__', nf90_int, &
      & varid=var_lasym), 'define lasym')
    call check_nc(nf90_def_var(ncid, 'jlist', nf90_int, [dim_surf], &
      & var_jlist), 'define jlist')
    call check_nc(nf90_def_var(ncid, 'ixm_b', nf90_int, [dim_mode], &
      & var_ixm), 'define ixm_b')
    call check_nc(nf90_def_var(ncid, 'ixn_b', nf90_int, [dim_mode], &
      & var_ixn), 'define ixn_b')
    call check_nc(nf90_def_var(ncid, 'iota_b', nf90_double, [dim_radius], &
      & var_iota), 'define iota_b')
    call check_nc(nf90_def_var(ncid, 'buco_b', nf90_double, [dim_radius], &
      & var_buco), 'define buco_b')
    call check_nc(nf90_def_var(ncid, 'bvco_b', nf90_double, [dim_radius], &
      & var_bvco), 'define bvco_b')
    call check_nc(nf90_def_var(ncid, 'phi_b', nf90_double, [dim_radius], &
      & var_phi), 'define phi_b')
    call check_nc(nf90_def_var(ncid, 'bmnc_b', nf90_double, &
      & [dim_mode, dim_surf], var_bmnc), 'define bmnc_b')
    call check_nc(nf90_def_var(ncid, 'rmnc_b', nf90_double, &
      & [dim_mode, dim_surf], var_rmnc), 'define rmnc_b')
    call check_nc(nf90_def_var(ncid, 'zmns_b', nf90_double, &
      & [dim_mode, dim_surf], var_zmns), 'define zmns_b')
    call check_nc(nf90_def_var(ncid, 'pmns_b', nf90_double, &
      & [dim_mode, dim_surf], var_pmns), 'define pmns_b')
    call check_nc(nf90_def_var(ncid, 'bmns_b', nf90_double, &
      & [dim_mode, dim_surf], var_bmns), 'define bmns_b')
    call check_nc(nf90_def_var(ncid, 'rmns_b', nf90_double, &
      & [dim_mode, dim_surf], var_rmns), 'define rmns_b')
    call check_nc(nf90_def_var(ncid, 'zmnc_b', nf90_double, &
      & [dim_mode, dim_surf], var_zmnc), 'define zmnc_b')
    call check_nc(nf90_def_var(ncid, 'pmnc_b', nf90_double, &
      & [dim_mode, dim_surf], var_pmnc), 'define pmnc_b')

    call check_nc(nf90_enddef(ncid), 'end boozmn definitions')

    jlist_file = [2]
    ixm_file = [0, 1]
    ixn_file = [0, 1]
    iota_file = [0.5_dp, 0.55_dp, 0.6_dp]
    buco_file = [0.0_dp, 0.01_dp, 0.02_dp]
    bvco_file = [0.2_dp, 0.25_dp, 0.3_dp]
    phi_file = [0.0_dp, 0.5_dp, 1.0_dp]

    call check_nc(nf90_put_var(ncid, var_ns, ns_full), 'write ns_b')
    call check_nc(nf90_put_var(ncid, var_nfp, 1), 'write nfp_b')
    call check_nc(nf90_put_var(ncid, var_mboz, 1), 'write mboz_b')
    call check_nc(nf90_put_var(ncid, var_nboz, 1), 'write nboz_b')
    call check_nc(nf90_put_var(ncid, var_mnboz, mn_modes), 'write mnboz_b')
    call check_nc(nf90_put_var(ncid, var_lasym, 1), 'write lasym')
    call check_nc(nf90_put_var(ncid, var_jlist, jlist_file), 'write jlist')
    call check_nc(nf90_put_var(ncid, var_ixm, ixm_file), 'write ixm_b')
    call check_nc(nf90_put_var(ncid, var_ixn, ixn_file), 'write ixn_b')
    call check_nc(nf90_put_var(ncid, var_iota, iota_file), 'write iota_b')
    call check_nc(nf90_put_var(ncid, var_buco, buco_file), 'write buco_b')
    call check_nc(nf90_put_var(ncid, var_bvco, bvco_file), 'write bvco_b')
    call check_nc(nf90_put_var(ncid, var_phi, phi_file), 'write phi_b')

    coeff = 0.0_dp
    coeff(1,1) = ASYM_B00
    coeff(2,1) = ASYM_BC11
    call check_nc(nf90_put_var(ncid, var_bmnc, coeff), 'write bmnc_b')
    coeff = 0.0_dp
    coeff(1,1) = 1.0_dp
    coeff(2,1) = 0.1_dp
    call check_nc(nf90_put_var(ncid, var_rmnc, coeff), 'write rmnc_b')
    coeff = 0.0_dp
    coeff(2,1) = 0.05_dp
    call check_nc(nf90_put_var(ncid, var_zmns, coeff), 'write zmns_b')
    coeff = 0.0_dp
    call check_nc(nf90_put_var(ncid, var_pmns, coeff), 'write pmns_b')

    coeff = 0.0_dp
    coeff(2,1) = ASYM_BS11
    call check_nc(nf90_put_var(ncid, var_bmns, coeff), 'write bmns_b')
    coeff = 0.0_dp
    coeff(2,1) = 0.03_dp
    call check_nc(nf90_put_var(ncid, var_rmns, coeff), 'write rmns_b')
    coeff = 0.0_dp
    coeff(2,1) = 0.04_dp
    call check_nc(nf90_put_var(ncid, var_zmnc, coeff), 'write zmnc_b')
    coeff = 0.0_dp
    call check_nc(nf90_put_var(ncid, var_pmnc, coeff), 'write pmnc_b')
    call check_nc(nf90_close(ncid), 'close boozmn')
  end subroutine write_asymmetric_boozmn_fixture

  subroutine check_nc(status, context)
    integer, intent(in) :: status
    character(len=*), intent(in) :: context

    if (status /= nf90_noerr) then
      print *, 'FAIL: NetCDF ', trim(context), ': ', trim(nf90_strerror(status))
      error stop
    end if
  end subroutine check_nc

end program test_boozmn_read
