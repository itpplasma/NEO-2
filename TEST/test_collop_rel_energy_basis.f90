! The relativistic energy-scattering matrix must approach the
! non-relativistic one for T_e -> 0, also when the projection and
! expansion bases differ.
program test_collop_rel_energy_basis
    use collop_compute, only: init_collop, compute_Minv, compute_energyscattering, &
        M_transform, M_transform_inv, isw_relativistic, rmu, norm_maxwell, &
        gamma_ab, m_a, m_b, T_a, T_b, integral_cutoff, lsw_split_interval, &
        x_cutoff, epsabs, epsrel, c, m_ele, eV, DBESK
    use rkstep_mod, only: lag, leg
    use collisionality_mod, only: phi_x_max, collop_bspline_order, collop_bspline_dist
    implicit none

    integer, parameter :: dp = kind(1.0d0)
    real(dp), parameter :: pi = 3.14159265358979d0
    real(dp), parameter :: T_e = 100.0_dp
    integer, parameter :: bases(2, 3) = reshape([11, 11, 0, 1, 1, 0], [2, 3])
    real(dp), allocatable :: den_nr(:,:), den_rel(:,:)
    real(dp) :: err, k2
    integer :: k, ierr
    logical :: ok

    lag = 3
    leg = 3
    collop_bspline_order = 3
    collop_bspline_dist = 1.0_dp
    phi_x_max = 4.0_dp
    integral_cutoff = .true.
    lsw_split_interval = .false.
    x_cutoff = 20.0_dp
    epsabs = 1d-11
    epsrel = 1d-11
    gamma_ab = 1.0_dp
    m_a = 1.0_dp
    m_b = 1.0_dp
    T_a = 1.0_dp
    T_b = 1.0_dp
    allocate(den_nr(lag + 1, lag + 1), den_rel(lag + 1, lag + 1))
    allocate(M_transform(0:lag, 0:lag), M_transform_inv(0:lag, 0:lag))

    ok = .true.
    do k = 1, size(bases, 2)
        call init_collop(bases(1, k), bases(2, k), 0.0_dp, 0.0_dp)

        isw_relativistic = 0
        call compute_Minv(M_transform_inv)
        call compute_energyscattering(den_nr)

        isw_relativistic = 1
        rmu = (c**2*m_ele)/(eV*T_e)
        call DBESK(rmu, 2, k2, ierr)
        norm_maxwell = sqrt(pi/2.0_dp)/k2
        call compute_Minv(M_transform_inv)
        call compute_energyscattering(den_rel)

        err = maxval(abs(den_rel - den_nr))/maxval(abs(den_nr))
        print '(a,2i4,a,es12.4)', 'bases prj/exp', bases(:, k), &
            ': relative difference', err
        if (err > 1d-2) then
            print *, 'FAIL'
            ok = .false.
        end if
    end do

    if (ok) then
        print *, 'All tests passed!'
    else
        error stop 1
    end if
end program test_collop_rel_energy_basis
