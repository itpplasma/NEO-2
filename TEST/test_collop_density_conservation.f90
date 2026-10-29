program test_collop_density_conservation
    ! Per-species particle conservation of the energy projection of the collision
    ! operator. With the Galerkin weight x^(4+alpha) exp(-x^2) the vector
    ! W_m = int_0^inf x^3 exp(-x^2) phi_m(x) dx is a left null vector of the energy
    ! scattering matrix denmm and of the field-particle matrix ailmm at l = 0. That
    ! requires x^(-1-alpha) in the span of the basis, which the partition of unity of
    ! the B-spline basis gives at scalprod_alpha = -1.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use collop, only: scalprod_alpha, init_collop, compute_collop, compute_source, &
        phi_exp, m_ele, m_d, integral_cutoff, lsw_split_interval, &
        num_sub_intervals, num_sub_intervals_cutoff, x_cutoff, epsabs, epsrel
    use collisionality_mod, only: phi_x_max, collop_bspline_order, collop_bspline_dist, &
        collop_bspline_taylor, isw_relativistic
    use rkstep_mod, only: lag, leg
    implicit none

    integer, parameter :: ncase = 5, ngl = 16, npanel0 = 20
    real(dp), parameter :: x_tail = 8.0d0, pi = 3.14159265358979323846d0
    ! ~10x the measured 1.0751e-10 defect from quadrature (1e-13) and Taylor expansions.
    real(dp), parameter :: tolerance = 1.1d-9
    real(dp), parameter :: w_converged = 1.0d-14

    character(len=3) :: tag_a(ncase), tag_b(ncase)
    real(dp) :: mass_a(ncase), mass_b(ncase), temp_a(ncase), temp_b(ncase)
    real(dp), allocatable :: w(:)
    real(dp), allocatable :: denmm_t(:,:), ailmm_t(:,:,:), anumm_t(:,:)
    real(dp), allocatable :: Amm(:,:), asource(:,:), weightlag(:,:)
    real(dp), allocatable :: weightden(:), weightparflow(:), weightenerg(:)
    real(dp) :: rel_den, rel_ail, worst, w_change
    integer :: c, npanel_used
    logical :: ok

    lag = 6
    leg = 3
    phi_x_max = 5.0d0
    collop_bspline_order = 4
    collop_bspline_dist = 1.0d0
    collop_bspline_taylor = .true.
    isw_relativistic = 0
    ! Same integration settings as collop_load uses in multispecies mode.
    integral_cutoff = .true.
    lsw_split_interval = .true.
    num_sub_intervals = 50
    num_sub_intervals_cutoff = 40
    x_cutoff = 3.0d3
    epsabs = 1.0d-13
    epsrel = 1.0d-13

    ! Equal temperatures, then unequal ones for the e-D pair.
    tag_a = ['e  ', 'e  ', 'D  ', 'D  ', 'e  ']
    tag_b = ['e  ', 'D  ', 'e  ', 'D  ', 'D  ']
    mass_a = [m_ele, m_ele, m_d, m_d, m_ele]
    mass_b = [m_ele, m_d, m_ele, m_d, m_d]
    temp_a = 1.0d0
    temp_b = [1.0d0, 1.0d0, 1.0d0, 1.0d0, 2.0d0]

    print '(A,ES10.3)', 'scalprod_alpha (module default) = ', scalprod_alpha

    call init_collop(11, 11, scalprod_alpha, 0.0d0)

    ! Storage of collop_compute is one element per basis function; its dummy
    ! arguments run from 1 to lag+1 while the matrices of rkstep_mod are 0-based.
    allocate(Amm(1:lag+1, 1:lag+1), asource(1:lag+1, 3), weightlag(3, 1:lag+1))
    allocate(weightden(1:lag+1), weightparflow(1:lag+1), weightenerg(1:lag+1))
    call compute_source(asource, weightlag, weightden, weightparflow, weightenerg, Amm)

    ! W is integrated here, independently of the weights of collop_compute.
    allocate(w(1:lag+1))
    call density_weights(w, npanel_used, w_change)
    print '(A,I0,A,ES12.5)', 'Gauss-Legendre panels: ', npanel_used, &
        ' , |change of W on doubling| = ', w_change

    allocate(denmm_t(1:lag+1, 1:lag+1), anumm_t(1:lag+1, 1:lag+1))
    allocate(ailmm_t(1:lag+1, 1:lag+1, 0:leg))
    ok = .true.
    worst = 0.0d0
    do c = 1, ncase
        print '(A,I0,A,F4.1)', 'case ', c, ', T_b/T_a = ', temp_b(c)/temp_a(c)
        call compute_collop(tag_a(c), tag_b(c), mass_a(c), mass_b(c), temp_a(c), temp_b(c), &
            anumm_t, denmm_t, ailmm_t)
        rel_den = defect(denmm_t, w)
        rel_ail = defect(ailmm_t(:, :, 0), w)
        print '(5A,ES12.5)', 'denmm  ', tag_a(c), '-', tag_b(c), &
            ' rel defect = ', rel_den
        print '(5A,ES12.5)', 'ailmm0 ', tag_a(c), '-', tag_b(c), &
            ' rel defect = ', rel_ail
        if (.not. ieee_is_finite(rel_den) .or. .not. ieee_is_finite(rel_ail)) then
            print '(A)', 'FAIL: non-finite density defect'
            ok = .false.
        else if (max(rel_den, rel_ail) > tolerance) then
            print '(A,I0)', 'FAIL: density conservation in case ', c
            ok = .false.
        end if
        worst = max(worst, rel_den, rel_ail)
    end do

    print '(A,ES12.5)', 'largest relative defect = ', worst
    if (.not. ok) error stop 1
    print '(A)', 'All tests passed!'

contains

    ! Scale by the largest absolute weighted column to avoid tiny-column noise.
    function defect(a, weights) result(rel)
        real(dp), intent(in) :: a(:,:), weights(:)
        real(dp) :: rel, den

        if (.not. all(ieee_is_finite(a))) error stop 'FAIL: non-finite matrix'
        den = maxval(matmul(abs(weights), abs(a)))
        if (den <= tiny(den)) error stop 'FAIL: zero matrix scale'
        rel = maxval(abs(matmul(weights, a)))/den
    end function defect

    ! W_m = int_0^inf x^3 exp(-x^2) phi_m(x) dx by composite Gauss-Legendre over
    ! [0, phi_x_max] and over the Taylor continuation of the basis up to x_tail,
    ! with the panels doubled until W has settled to w_converged.
    subroutine density_weights(weights, npanel_used, change)
        real(dp), intent(out) :: weights(:)
        integer, intent(out) :: npanel_used
        real(dp), intent(out) :: change
        real(dp), allocatable :: nodes(:), wts(:), previous(:)
        real(dp) :: a, b, h, total, xm
        integer :: m, i, j

        allocate(nodes(ngl), wts(ngl), previous(lbound(weights, 1):ubound(weights, 1)))
        call gauss_legendre(ngl, nodes, wts)
        npanel_used = npanel0
        weights = 0.0d0
        do
            previous = weights
            do m = lbound(weights, 1), ubound(weights, 1)
                total = 0.0d0
                h = phi_x_max/real(npanel_used, dp)
                do i = 1, npanel_used
                    a = real(i - 1, dp)*h
                    b = a + h
                    do j = 1, ngl
                        xm = 0.5d0*((b - a)*nodes(j) + a + b)
                        total = total + 0.5d0*(b - a)*wts(j)*xm**3*exp(-xm**2)*phi_exp(m-1, xm)
                    end do
                end do
                h = (x_tail - phi_x_max)/real(npanel_used, dp)
                do i = 1, npanel_used
                    a = phi_x_max + real(i - 1, dp)*h
                    b = a + h
                    do j = 1, ngl
                        xm = 0.5d0*((b - a)*nodes(j) + a + b)
                        total = total + 0.5d0*(b - a)*wts(j)*xm**3*exp(-xm**2)*phi_exp(m-1, xm)
                    end do
                end do
                weights(m) = total
            end do
            if (.not. all(ieee_is_finite(weights))) &
                error stop 'FAIL: non-finite density weights'
            change = maxval(abs(weights - previous))
            if (change < w_converged) exit
            if (npanel_used >= 5120) error stop 'FAIL: density weights did not converge'
            npanel_used = 2*npanel_used
        end do
    end subroutine density_weights

    ! Nodes and weights of the n-point Gauss-Legendre rule, the nodes from Newton's
    ! method on the Legendre polynomial P_n.
    subroutine gauss_legendre(n, nodes, wts)
        integer, intent(in) :: n
        real(dp), intent(out) :: nodes(:), wts(:)
        real(dp) :: p, dpn, guess, t
        integer :: i, iter

        do i = 1, n
            t = cos(pi*(real(i, dp) - 0.25d0)/real(n + 0.5, dp))
            do iter = 1, 100
                call legendre(n, t, p, dpn)
                guess = t - p/dpn
                if (abs(guess - t) <= 4.0d0*epsilon(1.0d0)*abs(t)) exit
                t = guess
            end do
            if (iter > 100) error stop 'FAIL: Gauss-Legendre root did not converge'
            call legendre(n, t, p, dpn)
            nodes(i) = t
            wts(i) = 2.0d0/((1.0d0 - t*t)*dpn*dpn)
        end do
    end subroutine gauss_legendre

    subroutine legendre(n, x, p, dpn)
        integer, intent(in) :: n
        real(dp), intent(in) :: x
        real(dp), intent(out) :: p, dpn
        real(dp) :: pm, pp
        integer :: k

        pm = 1.0d0
        p = x
        do k = 1, n - 1
            pp = pm
            pm = p
            p = ((2*k + 1)*x*pm - k*pp)/(k + 1)
        end do
        dpn = n*(x*p - pm)/(x*x - 1.0d0)
    end subroutine legendre
end program test_collop_density_conservation
