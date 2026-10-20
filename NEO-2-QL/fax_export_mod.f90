! Lossless, opt-in export at the axisymmetric solve boundary; no solver changes.
module fax_export_mod
    use nrtype, only: dp, pi
    use hdf5_tools
    implicit none
    private
    public :: export_fax, fax_basis
contains
    subroutine fax_basis(coeff, edges, centers, extend)
        use collop, only: collop_base_exp, phi_exp, d_phi_exp, dd_phi_exp
        use collop_laguerre, only: gencoeflag, phi_hm
        use collop_bspline, only: xknots
        use collisionality_mod, only: isw_lorentz, isw_relativistic, &
            collop_bspline_order, collop_bspline_taylor
        use rkstep_mod, only: lag
        real(dp), allocatable, intent(out) :: coeff(:, :, :), edges(:), centers(:)
        logical, intent(out) :: extend
        integer :: m, j

        if (isw_relativistic /= 0) error stop 'f_AX export: relativistic basis unsupported'
        extend = .true.
        if (isw_lorentz == 1) then
            allocate(coeff(0:0, 0:0, 1), edges(2), centers(1))
            coeff = 1.0_dp
            edges = [0.0_dp, 1.0_dp]
            centers = 0.0_dp
        else if (collop_base_exp == 0) then
            allocate(coeff(0:2*lag, 0:lag, 1), edges(2), centers(1))
            coeff = 0.0_dp
            do m = 0, lag
                coeff(0:2*m:2, m, 1) = gencoeflag(m, 0:m, 1) &
                    * pi**0.75_dp * phi_hm(m)
            end do
            edges = [0.0_dp, 1.0_dp]
            centers = 0.0_dp
        else if (collop_base_exp == 11) then
            if (collop_bspline_order /= 3) &
                error stop 'f_AX export: B-spline order must be 3'
            allocate(coeff(0:2, 0:lag, size(xknots)-1))
            allocate(edges(size(xknots)), centers(size(xknots)-1))
            edges = xknots
            centers = 0.5_dp * (edges(:size(edges)-1) + edges(2:))
            extend = collop_bspline_taylor
            do j = 1, size(centers)
                do m = 0, lag
                    coeff(0, m, j) = phi_exp(m, centers(j))
                    coeff(1, m, j) = d_phi_exp(m, centers(j))
                    coeff(2, m, j) = 0.5_dp * dd_phi_exp(m, centers(j))
                end do
            end do
        else
            error stop 'f_AX export: unsupported expansion basis'
        end if
    end subroutine fax_basis

    subroutine export_fax(ibeg, iend, eta, npl, offsets, phi, bhat, dlogb, &
            hphi, solution, flux, qflux, scalefac, btheta, bphi)
        use collop, only: collop_base_exp, collop_base_prj, scalprod_alpha, &
            scalprod_beta, phi_exp
        use rkstep_mod, only: lag, weightlag
        use collisionality_mod, only: num_spec, species_tag, m_spec, z_spec, &
            n_spec, T_spec, isw_lorentz, collop_bspline_order
        use mpiprovider_module, only: mpro
        use mag_interface_mod, only: boozer_s
        use neo_magfie, only: boozer_psi_pr_hat, boozer_iota
        use partpa_mod, only: bmod0
        use ntv_mod, only: dn_spec_ov_ds, dT_spec_ov_ds, Om_tE, isw_calc_Er, &
            y_ntv_mod
        integer, intent(in) :: ibeg, iend, npl(ibeg:), offsets(ibeg:)
        real(dp), intent(in) :: eta(0:), phi(ibeg:), bhat(ibeg:), dlogb(ibeg:)
        real(dp), intent(in) :: hphi(ibeg:), scalefac, btheta, bphi
        complex(dp), intent(in) :: solution(:, :, 0:), flux(:, :)
        real(dp), intent(in) :: qflux(:, :, 0:)
        real(dp), allocatable :: coeff(:, :, :), edges(:), centers(:), samples(:, :)
        real(dp), parameter :: sample_x(6) = [0.0_dp, 0.17_dp, 0.93_dp, &
            1.7_dp, 3.1_dp, 4.8_dp]
        integer(hid_t) :: file_id
        integer :: rank, m, j, k, b, t, row, index0, nband
        integer :: sample_row(8), sample_step(8), sample_band(8), sample_sign(8)
        real(dp) :: sample_g(3, num_spec, 6, 8)
        real(dp), allocatable :: geometry(:, :)
        logical :: extend
        character(len=128) :: fname

        if (mpro%getNumProcs() /= num_spec) &
            error stop 'f_AX export requires one MPI rank per field species'
        rank = mpro%getRank()
        call fax_basis(coeff, edges, centers, extend)
        allocate(samples(0:lag, size(sample_x)))
        do j = 1, size(sample_x)
            do m = 0, lag
                if (isw_lorentz == 1) then
                    samples(m, j) = 1.0_dp
                else
                    samples(m, j) = phi_exp(m, sample_x(j))
                end if
            end do
        end do
        allocate(geometry(3, size(solution, 1)))
        do t = ibeg, iend
            nband = npl(t)+1
            do m = 0, lag
                do row = 1, 2*nband
                    geometry(:, offsets(t)+2*nband*m+row) = &
                        real(flux(:, offsets(t)+row), dp)/weightlag(:, 0)
                end do
            end do
        end do
        do k = 1, 8
            t = ibeg + ((k-1)/2)*(iend-ibeg)/3
            nband = npl(t)+1
            sample_step(k) = t-ibeg
            sample_band(k) = max(1, nband/2)
            sample_sign(k) = 1-2*mod(k, 2)
            index0 = offsets(t)+sample_band(k)
            if (sample_sign(k) < 0) &
                index0 = offsets(t)+2*nband-sample_band(k)+1
            sample_row(k) = index0
            do j = 1, 6
                do b = 1, num_spec
                    sample_g(:, b, j, k) = 0.0_dp
                    do m = 0, lag
                        sample_g(:, b, j, k) = sample_g(:, b, j, k) &
                            + real(solution(index0+2*nband*m, 1:3, b-1), dp) &
                            * samples(m, j)
                    end do
                end do
            end do
        end do
        write(fname, '(a, i0, a)') 'fax_background_spec', rank, '.h5'
        call h5_create(trim(fname), file_id)
        call h5_add(file_id, 'format_version', 1)
        call h5_add(file_id, 'field_species', rank)
        call h5_add(file_id, 'num_species', num_spec)
        call h5_add(file_id, 'species_tag', species_tag)
        call h5_add(file_id, 'mass_g', m_spec)
        call h5_add(file_id, 'charge_number', z_spec)
        call h5_add(file_id, 'density_cm3', n_spec)
        call h5_add(file_id, 'temperature_erg', T_spec)
        if (allocated(dn_spec_ov_ds)) &
            call h5_add(file_id, 'dn_ds_cm3', dn_spec_ov_ds)
        if (allocated(dT_spec_ov_ds)) &
            call h5_add(file_id, 'dT_ds_erg', dT_spec_ov_ds)
        call h5_add(file_id, 'Om_tE', Om_tE)
        call h5_add(file_id, 'isw_calc_Er', isw_calc_Er)
        call h5_add(file_id, 's_tor', boozer_s)
        call h5_add(file_id, 'B_ref_G', bmod0*1.0e4_dp)
        call h5_add(file_id, 'psi_pr_hat', boozer_psi_pr_hat)
        call h5_add(file_id, 'iota', boozer_iota)
        call h5_add(file_id, 'scalefac_kG', scalefac)
        call h5_add(file_id, 'bcovar_theta_hat', btheta)
        call h5_add(file_id, 'bcovar_phi_hat', bphi)
        call h5_add(file_id, 'y_ntv', y_ntv_mod)
        call h5_add(file_id, 'lag', lag)
        call h5_add(file_id, 'isw_lorentz', isw_lorentz)
        call h5_add(file_id, 'collop_base_exp', collop_base_exp)
        call h5_add(file_id, 'collop_base_prj', collop_base_prj)
        call h5_add(file_id, 'scalprod_alpha', scalprod_alpha)
        call h5_add(file_id, 'scalprod_beta', scalprod_beta)
        call h5_add(file_id, 'bspline_order', collop_bspline_order)
        call h5_add(file_id, 'basis_coeff', coeff, &
            lbound(coeff), ubound(coeff))
        call h5_add(file_id, 'basis_edges', edges)
        call h5_add(file_id, 'basis_centers', centers)
        call h5_add(file_id, 'basis_extend', extend)
        call h5_add(file_id, 'sample_x', sample_x, &
            [1], [size(sample_x)])
        call h5_add(file_id, 'sample_basis', samples, &
            lbound(samples), ubound(samples))
        call h5_add(file_id, 'weightlag', weightlag, &
            lbound(weightlag), ubound(weightlag))
        call h5_add(file_id, 'ibeg', ibeg)
        call h5_add(file_id, 'eta', eta, &
            [0], [size(eta)-1])
        call h5_add(file_id, 'npassing', npl(ibeg:iend), &
            [ibeg], [iend])
        call h5_add(file_id, 'offset', offsets(ibeg:iend), &
            [ibeg], [iend])
        call h5_add(file_id, 'phi', phi(ibeg:iend), &
            [ibeg], [iend])
        call h5_add(file_id, 'bhat', bhat(ibeg:iend), &
            [ibeg], [iend])
        call h5_add(file_id, 'dlogb_dphi', dlogb(ibeg:iend), &
            [ibeg], [iend])
        call h5_add(file_id, 'h_phi', hphi(ibeg:iend), &
            [ibeg], [iend])
        call h5_add(file_id, 'g_real', real(solution(:, 1:3, :), dp), &
            [1,1,0], [size(solution,1),3,num_spec-1])
        call h5_add(file_id, 'g_imag', aimag(solution(:, 1:3, :)), &
            [1,1,0], [size(solution,1),3,num_spec-1])
        call h5_add(file_id, 'flux_real', real(flux, dp), &
            [1,1], shape(flux))
        call h5_add(file_id, 'flux_imag', aimag(flux), &
            [1,1], shape(flux))
        ! Axisymmetric public qflux is twice solve_eqs' half contraction.
        call h5_add(file_id, 'qflux', 2.0_dp*qflux, &
            [1,1,0], [3,3,num_spec-1])
        call h5_add(file_id, 'flux_geometry', geometry, &
            [1,1], shape(geometry))
        call h5_add(file_id, 'sample_g', sample_g, [1,1,1,1], shape(sample_g))
        call h5_add(file_id, 'sample_row', sample_row, [1], [8])
        call h5_add(file_id, 'sample_step', sample_step, [1], [8])
        call h5_add(file_id, 'sample_band', sample_band, [1], [8])
        call h5_add(file_id, 'sample_sign', sample_sign, [1], [8])
        call h5_close(file_id)
    end subroutine export_fax
end module fax_export_mod
