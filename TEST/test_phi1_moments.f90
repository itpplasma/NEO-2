program test_phi1_moments
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use phi1_mod, only: phi1_exb_band, phi1_density
    use collop_compute, only: compute_phi1_speed, compute_sources, phi_prj, &
        phi_exp, lagmax, alpha, beta, make_ortho, M_transform_inv, &
        weightenerg_offset
    use collop_polynomial, only: phi_polynomial_2
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), tol = 2.0e-12_dp
    real(dp) :: bref(3), psi_prime(3), btheta(3), bphi(3), iota(3), bhat(3)
    real(dp) :: theta(3), phi(3), vthermal(3), hphi(3), rho(3)
    real(dp) :: aE1(0:1), source(2,3), flux(2,3), den(2), par(2), en(2)
    real(dp) :: pitch(5)
    complex(dp) :: expected_us(3), expected_q(3), expected_a1(2,3)
    complex(dp) :: expected_a2(2,3), us, q(1), psi, y(10,2), value
    complex(dp) :: expected_density, expected_energy
    integer :: m(3), n(3), j, unit
    character(1024) :: fixture
    namelist /exb_fixture/ bref, psi_prime, btheta, bphi, iota, bhat, m, n, &
        theta, phi, vthermal, hphi, rho, expected_us, expected_q, &
        expected_a1, expected_a2
    namelist /density_fixture/ pitch, y, expected_density, expected_energy

    call get_command_argument(1,fixture)
    open(newunit=unit,file=trim(fixture),status='old')
    read(unit,nml=exb_fixture)
    close(unit)
    ! This fixture uses (1,x^2), independently of the production default basis.
    lagmax = 1
    alpha = 0.0_dp
    beta = 0.0_dp
    make_ortho = .true.
    phi_prj => phi_polynomial_2
    phi_exp => phi_polynomial_2
    allocate(M_transform_inv(2,2), weightenerg_offset(0:1))
    M_transform_inv = reshape([35.0_dp,-10.0_dp,-10.0_dp,4.0_dp],[2,2]) &
        *(4.0_dp*pi/15.0_dp)
    call compute_phi1_speed(aE1)
    call compute_sources(source,flux,den,par,en)
    ! compute_sources returns pre-swap A1,A2,A3. A1 already projects x.
    do j = 1, 3
        psi = exp(cmplx(0.0_dp,m(j)*theta(j)+n(j)*phi(j),dp))
        call phi1_exb_band(psi,m(j),n(j),btheta(j),bphi(j),iota(j), &
            bref(j),psi_prime(j),vthermal(j),rho(j),bhat(j),hphi(j), &
            [0.8_dp],[0.3_dp],us,q)
        call require(abs(us-expected_us(j)) < tol*abs(expected_us(j)), 'ExB drift')
        call require(abs(q(1)-expected_q(j)) < tol*abs(expected_q(j)), 'band source')
        call require(maxval(abs(aE1*q(1)-expected_a1(:,j))) &
            < tol*maxval(abs(expected_a1(:,j))), 'A1 x^-1 projection')
        call require(maxval(abs(source(:,1)*q(1)-expected_a2(:,j))) &
            < tol*maxval(abs(expected_a2(:,j))), 'A2 x projection')
    end do

    ! Both sigma branches, P0/P1/P2 and an exactly integrated turning cell.
    call get_command_argument(2,fixture)
    open(newunit=unit,file=trim(fixture),status='old')
    read(unit,nml=density_fixture)
    close(unit)
    value = phi1_density(y,pitch,[.25_dp,.375_dp],.005_dp)
    call require(abs(value-expected_density) < tol, 'density moment')
    call require(abs(2.0e12_dp*.02_dp*value &
        -cmplx(5.2e10_dp,4.0e9_dp,dp)) < 1.0_dp, 'physical density')
    value = phi1_density(y,pitch,[.375_dp,.9375_dp],.005_dp)
    call require(abs(value-expected_energy) < tol, 'energy moment')
    print *, 'All tests passed!'
contains
    subroutine require(condition, name)
        logical, intent(in) :: condition
        character(*), intent(in) :: name
        if (.not. condition) then
            print *, 'FAIL: ', name
            error stop 1
        end if
    end subroutine require
end program test_phi1_moments
