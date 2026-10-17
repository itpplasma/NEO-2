program test_collision_krook
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, &
        ieee_positive_inf
    use collision_model_mod, only: resolve_collision_model, uses_krook, &
        krook_carrier_scale, krook_pitch_mass
    use collisionality_mod, only: collision_model, krook_nu, num_spec, &
        krook_geometry_trapped, krook_geometry_passing, &
        isw_lorentz, isw_energy, isw_integral, isw_momentum, isw_relativistic, &
        lsw_nbi, lsw_multispecies, collpar_spec, m_spec, T_spec, n_spec, z_spec, &
        phi_x_max, collop_bspline_order, collop_bspline_dist, &
        collop_bspline_taylor
    use collop, only: collop_load, collop_construct, collop_base_prj, &
        collop_base_exp, scalprod_alpha, scalprod_beta, lsw_read_precom, &
        lsw_write_precom
    use rkstep_mod, only: lag, leg, anumm_a, denmm_a, ailmm_aa
    use hdf5_tools, only: h5_init, h5_deinit
    use mpiprovider_module, only: mpro
    implicit none
    integer :: mode, a, power
    real(dp) :: expected(2,2), norm(2), i0, i2, i4, scale, vth
    real(dp) :: max_error
    character(len=40) :: arg
    interface
        subroutine dgesv(n,nrhs,matrix,lda,ipiv,rhs,ldb,info)
            import dp
            integer, intent(in) :: n, nrhs, lda, ldb
            integer, intent(out) :: ipiv(*), info
            real(dp), intent(inout) :: matrix(lda,*), rhs(ldb,*)
        end subroutine dgesv
        subroutine zgesv(n,nrhs,matrix,lda,ipiv,rhs,ldb,info)
            import dp
            integer, intent(in) :: n, nrhs, lda, ldb
            integer, intent(out) :: ipiv(*), info
            complex(dp), intent(inout) :: matrix(lda,*), rhs(ldb,*)
        end subroutine zgesv
    end interface

    allocate(krook_nu(3))
    krook_nu = [10.0_dp, 20.0_dp, 30.0_dp]
    num_spec = 3
    isw_lorentz = 0
    isw_energy = 0
    isw_integral = 0
    isw_momentum = 0
    isw_relativistic = 0
    lsw_nbi = .false.
    collision_model = 'krook_constant'
    call get_command_argument(1, arg)
    if (len_trim(arg) /= 0) then
        select case (trim(arg))
        case ('unknown-model')
            collision_model = 'krook_harmonic'
        case ('negative-rate')
            krook_nu(2) = -1.0_dp
        case ('missing-rate')
            krook_nu(3) = 0.0_dp
        case ('nan-rate')
            krook_nu(1) = ieee_value(0.0_dp, ieee_quiet_nan)
        case ('zero-trapped')
            krook_geometry_trapped = 0.0_dp
        case ('negative-passing')
            krook_geometry_passing = -1.0_dp
        case ('nan-trapped')
            krook_geometry_trapped = ieee_value(0.0_dp, ieee_quiet_nan)
        case ('nan-passing')
            krook_geometry_passing = ieee_value(0.0_dp, ieee_quiet_nan)
        case ('infinite-trapped')
            krook_geometry_trapped = ieee_value(0.0_dp, ieee_positive_inf)
        case ('infinite-passing')
            krook_geometry_passing = ieee_value(0.0_dp, ieee_positive_inf)
        case ('integral')
            isw_integral = 1
        case ('grid')
            isw_momentum = 1
        case ('cache')
            call resolve_collision_model(.true., .false.)
        case default
            error stop 'Unknown test rejection mode'
        end select
        call resolve_collision_model(.false., .false.)
        print *, 'FAIL: rejection test unexpectedly accepted input'
        stop 0
    end if
    call resolve_collision_model(.false., .false.)
    if (.not. uses_krook()) error stop 'Krook dispatch absent'
    call check_named_modes()
    call mpro%init()
    call h5_init()
    lsw_multispecies = .true.
    allocate(collpar_spec(0:2), m_spec(0:2), T_spec(0:2), &
        n_spec(0:2), z_spec(0:2))
    collpar_spec = [2.0e-8_dp, 7.0e-8_dp, 4.0e-8_dp]
    m_spec = [3.343583719e-24_dp, 5.0073567e-24_dp, 9.109382150e-28_dp]
    T_spec = 2000.0_dp*1.602176634e-12_dp
    n_spec = [4.0e9_dp, 6.0e9_dp, 1.0e10_dp]
    z_spec = [1.0_dp, 1.0_dp, -1.0_dp]
    lag = 1
    leg = 2
    collop_base_prj = 0
    collop_base_exp = 0
    scalprod_alpha = 0.0_dp
    scalprod_beta = 0.0_dp
    lsw_read_precom = .false.
    lsw_write_precom = .false.
    phi_x_max = 5.0_dp
    collop_bspline_order = 4
    collop_bspline_dist = 1.0_dp
    collop_bspline_taylor = .true.
    max_error = 0.0_dp
    do mode = 1, 2
        if (mode == 1) then
            collision_model = 'krook_constant'
            power = 3
        else
            collision_model = 'krook_inverse_cube'
            power = 0
        end if
        call collop_construct()
        call collop_load()
        ! Independent exact moments of the first two normalized Laguerre
        ! polynomials: L0=1, L1=5/2-u^2, with Gaussian weight.
        norm = sqrt(2.0_dp/[gamma(2.5_dp), gamma(3.5_dp)])
        i0 = 0.5_dp*gamma(0.5_dp*(power+1))
        i2 = 0.5_dp*gamma(0.5_dp*(power+3))
        i4 = 0.5_dp*gamma(0.5_dp*(power+5))
        expected(1,1) = norm(1)**2*i0
        expected(1,2) = norm(1)*norm(2)*(2.5_dp*i0-i2)
        expected(2,1) = expected(1,2)
        expected(2,2) = norm(2)**2*(6.25_dp*i0-5.0_dp*i2+i4)
        do a = 0, 2
            vth = sqrt(2.0_dp*T_spec(a)/m_spec(a))
            scale = 2.0_dp*krook_nu(a+1)/(collpar_spec(a)*vth)
            max_error = max(max_error, maxval(abs(denmm_a(:,:,a)/scale+expected)))
            if (maxval(abs(denmm_a(:,:,a)/scale+expected)) > 2.0e-11_dp) &
                error stop 'Native reaction matrix failed analytic Laguerre oracle'
            call check_energy_resolvent(a,expected)
        end do
        if (any(anumm_a /= 0.0_dp)) error stop 'Krook retains pitch diffusion'
        if (any(ailmm_aa /= 0.0_dp)) error stop 'Krook retains interspecies FP integral'
        call check_pitch_action(expected(1,1))
        call check_geometry_projection()
    end do
    call h5_deinit()
    call mpro%deinit()
    print *, 'Maximum native matrix error: ', max_error
    print *, 'All tests passed!'
contains
    subroutine check_energy_resolvent(species, reaction)
        integer, intent(in) :: species
        real(dp), intent(in) :: reaction(2,2)
        real(dp) :: constant(2,2), n0, n1, nu, speed, kappa
        complex(dp) :: matrix(2,2), analytic(2,2), rhs(2), expected(2), forcing(2)
        complex(dp) :: source(2), determinant
        integer :: ipiv(2), info

        n0 = sqrt(2.0_dp/gamma(2.5_dp))
        n1 = sqrt(2.0_dp/gamma(3.5_dp))
        constant(1,1) = n0*n0/2.0_dp
        constant(1,2) = n0*n1/4.0_dp
        constant(2,1) = constant(1,2)
        constant(2,2) = n1*n1*1.125_dp
        source = [cmplx(3.0_dp,1.0_dp,dp), cmplx(-1.0_dp,0.5_dp,dp)]
        forcing(1) = constant(1,1)*source(1)+constant(1,2)*source(2)
        forcing(2) = constant(2,1)*source(1)+constant(2,2)*source(2)
        rhs = forcing
        nu = krook_nu(species+1)
        kappa = collpar_spec(species)
        speed = sqrt(2.0_dp*T_spec(species)/m_spec(species))
        matrix = -0.5_dp*kappa*speed*denmm_a(:,:,species) + &
            cmplx(0.0_dp,-7.0_dp,dp)*constant
        analytic = nu*reaction + cmplx(0.0_dp,-7.0_dp,dp)*constant
        determinant = analytic(1,1)*analytic(2,2)-analytic(1,2)*analytic(2,1)
        expected(1) = (analytic(2,2)*rhs(1)-analytic(1,2)*rhs(2))/determinant
        expected(2) = (analytic(1,1)*rhs(2)-analytic(2,1)*rhs(1))/determinant
        call zgesv(2,1,matrix,2,ipiv,rhs,2,info)
        if (info /= 0) error stop 'Native reaction resolvent is singular'
        if (maxval(abs(rhs-expected)) > 2.0e-11_dp) &
            error stop 'Native two-energy-mode resolvent failed analytic oracle'
        if (real(dot_product(forcing,rhs),dp) <= 0.0_dp) &
            error stop 'Native reaction resolvent has negative dissipation'
    end subroutine check_energy_resolvent

    subroutine check_named_modes()
        character(len=32), parameter :: names(5) = [character(len=32) :: &
            'native', 'pitch_angle', 'test_particle', 'full_fp', &
            'monoenergetic_lorentz']
        integer, parameter :: switches(3,5) = reshape( &
            [1,0,0, 0,0,0, 0,1,0, 0,1,1, 1,0,0], [3,5])
        integer :: j
        krook_nu = 0.0_dp
        do j = 1, size(names)
            collision_model = names(j)
            isw_lorentz = 1
            isw_energy = 0
            isw_integral = 0
            call resolve_collision_model(.false., .false.)
            if (isw_lorentz /= switches(1,j) .or. &
                isw_energy /= switches(2,j) .or. isw_integral /= switches(3,j)) &
                error stop 'Named collision switch resolution is incorrect'
        end do
        collision_model = 'krook_constant'
        isw_lorentz = 0
        krook_nu = [10.0_dp, 20.0_dp, 30.0_dp]
    end subroutine check_named_modes

    subroutine check_geometry_projection()
        real(dp), parameter :: bounds(0:4) = [1.0_dp,0.8_dp,0.6_dp,0.4_dp,0.0_dp]
        real(dp), parameter :: cuts(5) = [0.0_dp,0.3_dp,0.5_dp,0.6_dp,1.0_dp]
        real(dp) :: bands(4,4), inverse(4,4), moments(4), flux(4), row(4)
        real(dp) :: baseline(4), legacy(4), scaled(4), reverse(4), drift(4)
        real(dp) :: lower, upper, cut, bhat, eta_t, factorial, value, oracle
        real(dp) :: trapped_integral, passing_integral, error, matrix_error
        integer :: band, j, n, degree, icut, ipiv(4), info

        error = 0.0_dp
        matrix_error = 0.0_dp
        bhat = 1.7_dp
        do band = 1, 4
            lower = bounds(band)
            upper = bounds(band-1)
            ! Independent reconstruction: flux=2/B integral lambda*f d lambda.
            ! Taylor basis is local to this row, as in the production solvers.
            factorial = 1.0_dp
            do n = 0, 3
                if (n > 0) factorial = factorial*real(n,dp)
                do j = 1, 4
                    bands(j,n+1) = 2.0_dp/bhat/factorial*( &
                        ((bounds(j-1)-lower)**(n+2)-(bounds(j)-lower)**(n+2)) &
                        /real(n+2,dp) + lower*((bounds(j-1)-lower)**(n+1) &
                        -(bounds(j)-lower)**(n+1))/real(n+1,dp))
                end do
                moments(n+1) = (upper-lower)**(n+1)/(factorial*real(n+1,dp))
            end do
            inverse = 0.0_dp
            do j = 1, 4
                inverse(j,j) = 1.0_dp
            end do
            call dgesv(4,4,bands,4,ipiv,inverse,4,info)
            if (info /= 0) error stop 'Geometry reconstruction is singular'
            krook_geometry_trapped = 1.0_dp
            krook_geometry_passing = 1.0_dp
            baseline = krook_pitch_mass(0.37_dp,1.3_dp,bhat,moments,inverse, &
                lower,upper,0.75_dp/bhat)
            do j = 1, 4
                legacy(j) = dot_product(moments,inverse(:,j))* &
                    ((bhat*1.3_dp)**(-1.0_dp)*0.37_dp*(-1.0_dp))
            end do
            if (any(transfer(baseline,[0_int64],4) /= transfer(legacy,[0_int64],4))) &
                error stop 'Default geometry changed legacy pitch arithmetic'
            krook_geometry_trapped = 4.0_dp
            krook_geometry_passing = 4.0_dp
            scaled = krook_pitch_mass(0.37_dp,1.3_dp,bhat,moments,inverse, &
                lower,upper,0.75_dp/bhat)
            if (any(scaled /= 4.0_dp*baseline)) &
                error stop 'Equal geometry factors failed exact matrix scaling'
            ! Tensor with each actual species/energy matrix for BOTH rate laws.
            do j = 1, 4
                matrix_error = max(matrix_error, maxval(abs( &
                    denmm_a*scaled(j)-4.0_dp*(denmm_a*baseline(j)))) &
                    /max(1.0_dp,maxval(abs(4.0_dp*denmm_a*baseline(j)))))
            end do
            if (matrix_error > 32.0_dp*epsilon(1.0_dp)) &
                error stop 'Equal factors failed native energy-pitch tensor scaling'
            drift = krook_pitch_mass(0.37_dp,1.3_dp,bhat,moments,inverse, &
                lower,upper,0.75_dp/bhat,with_geometry=.false.)
            if (any(transfer(drift,[0_int64],4) /= transfer(baseline,[0_int64],4))) &
                error stop 'Common geometry factor changed rotation/drift mass'
            krook_geometry_trapped = 2.5_dp
            krook_geometry_passing = 0.75_dp
            do icut = 1, size(cuts)
                eta_t = (1.0_dp-cuts(icut)**2)/bhat
                row = krook_pitch_mass(0.37_dp,1.3_dp,bhat,moments,inverse, &
                    lower,upper,eta_t)
                drift = krook_pitch_mass(0.37_dp,1.3_dp,bhat,moments,inverse, &
                    lower,upper,eta_t,with_geometry=.false.)
                if (any(transfer(drift,[0_int64],4) /= transfer(baseline,[0_int64],4))) &
                    error stop 'Pitch geometry step changed rotation/drift mass'
                reverse = krook_pitch_mass(0.37_dp,-1.3_dp,bhat,moments,inverse, &
                    lower,upper,eta_t)
                if (maxval(abs(row+reverse)) > 3.0e-14_dp) &
                    error stop 'Geometry projection failed direction reversal'
                cut = max(lower,min(upper,cuts(icut)))
                do degree = 0, 3
                    do j = 1, 4
                        flux(j) = 2.0_dp/bhat*(bounds(j-1)**(degree+2) - &
                            bounds(j)**(degree+2))/real(degree+2,dp)
                    end do
                    ! Direct global-power integrals, independent of the
                    ! generated local Taylor moments (including crossing bands).
                    trapped_integral = (cut**(degree+1)-lower**(degree+1)) &
                        /real(degree+1,dp)
                    passing_integral = (upper**(degree+1)-cut**(degree+1)) &
                        /real(degree+1,dp)
                    oracle = -0.37_dp/(1.3_dp*bhat)*( &
                        2.5_dp*trapped_integral+0.75_dp*passing_integral)
                    value = dot_product(row,flux)
                    error = max(error,abs(value-oracle))
                    if (abs(value-oracle) > 3.0e-12_dp) &
                        error stop 'Geometry projection failed analytic step oracle'
                end do
            end do
        end do
        krook_geometry_trapped = 1.0_dp
        krook_geometry_passing = 1.0_dp
        print *, 'Geometry projection: 80 integrals, max error = ', error
        print *, 'Equal-factor native tensor max relative error = ', matrix_error
    end subroutine check_geometry_projection

    subroutine check_pitch_action(velocity_moment)
        real(dp), intent(in) :: velocity_moment
        real(dp), parameter :: bounds(0:4) = [1.0_dp,0.8_dp,0.6_dp,0.4_dp,0.0_dp]
        real(dp) :: bands(4,4), inverse(4,4), moments(4), flux(4), row(4)
        real(dp) :: reverse_row(4)
        real(dp) :: width, actual, oracle, coefficient, physical_rate, streaming
        real(dp) :: source, omega, pitch_mass, measured_rate
        complex(dp) :: solution, expected_solution
        integer :: j, degree, ipiv(4), info, species

        ! Build the exact native band-integrated-flux reconstruction for
        ! polynomials about each band's lower lambda endpoint. This includes
        ! the final band ending at the turning point lambda=0.
        do j = 1, 4
            do degree = 0, 3
                bands(j,degree+1) = 2.0_dp*(bounds(j-1)**(degree+2) - &
                    bounds(j)**(degree+2))/real(degree+2,dp)
            end do
        end do
        inverse = 0.0_dp
        do j = 1, 4
            inverse(j,j) = 1.0_dp
        end do
        call dgesv(4,4,bands,4,ipiv,inverse,4,info)
        if (info /= 0) error stop 'Independent polynomial reconstruction failed'
        do j = 1, 4
            do degree = 0, 3
                moments(degree+1) = (bounds(j-1)**(degree+1)- &
                    bounds(j)**(degree+1))/real(degree+1,dp)
            end do
            row = krook_pitch_mass(0.37_dp, 1.3_dp, 1.0_dp, moments, inverse, &
                bounds(j), bounds(j-1), 0.75_dp)
            reverse_row = krook_pitch_mass(0.37_dp,-1.3_dp,1.0_dp,moments,inverse, &
                bounds(j), bounds(j-1), 0.75_dp)
            if (maxval(abs(reverse_row+row)) > 3.0e-14_dp) &
                error stop 'Native pitch mass failed field-direction reversal'
            do degree = 0, 3
                do species = 1, 4
                    flux(species) = 2.0_dp*(bounds(species-1)**(degree+2)- &
                        bounds(species)**(degree+2))/real(degree+2,dp)
                end do
                actual = dot_product(row,flux)
                oracle = -0.37_dp/1.3_dp*(bounds(j-1)**(degree+1)- &
                    bounds(j)**(degree+1))/real(degree+1,dp)
                if (abs(actual-oracle) > 3.0e-12_dp) &
                    error stop 'Native pitch mass failed polynomial weak-form oracle'
            end do
            ! Constant source periodic Fourier solution. Normalize native
            ! reaction and streaming by the same independently integrated
            ! pitch mass; species-dependent vth and arbitrary kappa cancel.
            do species = 0, 2
                physical_rate = krook_nu(species+1)
                width = sqrt(2.0_dp*T_spec(species)/m_spec(species))
                coefficient = krook_carrier_scale(physical_rate, &
                    collpar_spec(species),T_spec(species),m_spec(species))
                pitch_mass = 2.0_dp*(bounds(j-1)-bounds(j))
                actual = -coefficient*velocity_moment* &
                    (-collpar_spec(species)*pitch_mass/2.0_dp)
                streaming = (4.0_dp/(3.0_dp*sqrt(acos(-1.0_dp))))*pitch_mass/width
                measured_rate = actual/streaming
                physical_rate = physical_rate*velocity_moment/ &
                    (4.0_dp/(3.0_dp*sqrt(acos(-1.0_dp))))
                if (abs(measured_rate/physical_rate-1.0_dp) > 5.0e-14_dp) &
                    error stop 'Physical damping normalization failed'
                source = 3.0_dp
                omega = -7.0_dp
                solution = source/cmplx(measured_rate,omega,dp)
                expected_solution = source*cmplx(physical_rate,-omega,dp)/ &
                    (physical_rate**2+omega**2)
                if (abs(solution-expected_solution) > 2.0e-14_dp) &
                    error stop 'Native manufactured resolvent failed'
                if (real(solution,dp) <= 0.0_dp) error stop 'Negative damping'
            end do
        end do
    end subroutine check_pitch_action
end program test_collision_krook
