program test_collop_energy_conservation
    ! The linearised Fokker-Planck operator conserves total energy pairwise:
    ! the energy that a perturbation of species a loses through its
    ! test-particle part (denmm_ab) is gained by species b through its
    ! field-particle part (ailmm_ba, l = 0). With the NEO-2-QL row scaling
    ! (collpar_a, taa_ov_tab, kinetic equation divided by v) the balance reads
    !   (v_ta/T_a) E' denmm_ab / 2 + (v_tb/T_b) E' ailmm_ba = 0,
    ! where E is the bra of x^2 - 3/2 in the B-spline basis (exact for
    ! scalprod_alpha = -1). Checked for e-e, D-D, e-D and D-e.
    use rkstep_mod, only: lag, leg
    use collisionality_mod, only: phi_x_max, collop_bspline_order, collop_bspline_dist, &
        collop_bspline_taylor, isw_relativistic
    use collop_compute, only: init_collop, compute_source, compute_collop, inv, &
        M_transform_inv, integral_cutoff, lsw_split_interval, num_sub_intervals, &
        num_sub_intervals_cutoff, x_cutoff, epsabs, epsrel, m_ele, m_d
    implicit none
    integer, parameter :: dp = kind(1.0d0), ns = 2
    real(dp), parameter :: tol = 1.0d-8
    real(dp) :: mass(ns), temp(ns), weight(ns), worst
    real(dp), allocatable :: asource(:,:), weightlag(:,:), bzero(:), weightparflow(:)
    real(dp), allocatable :: weightenerg(:), amm(:,:), mass_matrix(:,:), bra(:)
    real(dp), allocatable :: anumm(:,:,:,:), denmm(:,:,:,:), ailmm(:,:,:,:,:), tp(:), fp(:)
    integer :: a, b, n
    logical :: failed

    lag = 3
    leg = 3
    phi_x_max = 4.0d0
    collop_bspline_order = 3
    collop_bspline_dist = 1.0d0
    collop_bspline_taylor = .true.
    isw_relativistic = 0
    ! Multispecies integration settings of collision_operator_mems.f90
    integral_cutoff = .true.
    lsw_split_interval = .true.
    num_sub_intervals = 50
    num_sub_intervals_cutoff = 40
    x_cutoff = 3.0d3
    epsabs = 1.0d-13
    epsrel = 1.0d-13

    n = lag + 1
    mass = [m_ele, m_d]
    temp = [0.96d0, 1.0d0]
    allocate(asource(n, 3), weightlag(3, n), bzero(n), weightparflow(n), weightenerg(n))
    allocate(amm(n, n), mass_matrix(n, n), bra(n), tp(n), fp(n))
    allocate(anumm(n, n, ns, ns), denmm(n, n, ns, ns), ailmm(n, n, leg + 1, ns, ns))

    call init_collop(11, 11, -1.0d0, 0.0d0)
    call compute_source(asource, weightlag, bzero, weightparflow, weightenerg, amm)
    mass_matrix = M_transform_inv
    call inv(mass_matrix)
    bra = matmul(transpose(mass_matrix), weightenerg)
    do a = 1, ns
        do b = 1, ns
            call compute_collop('a', 'b', mass(a), mass(b), temp(a), temp(b), &
                anumm(:, :, a, b), denmm(:, :, a, b), ailmm(:, :, :, a, b))
        end do
    end do

    weight = sqrt(2.0d0*temp/mass)/temp
    failed = .false.
    do a = 1, ns
        do b = 1, ns
            tp = weight(a)*0.5d0*matmul(bra, denmm(:, :, a, b))
            fp = weight(b)*matmul(bra, ailmm(:, :, 1, b, a))
            worst = maxval(abs(tp + fp)/(abs(tp) + abs(fp)))
            write (*, '(a, 2i2, a, es10.2)') ' perturbed, field species', a, b, &
                ': relative energy defect ', worst
            if (.not. (worst <= tol)) failed = .true.
        end do
    end do

    if (failed) then
        print *, 'FAIL: collision matrices do not conserve total energy'
        stop 1
    end if
    print *, 'All tests passed!'
end program test_collop_energy_conservation
