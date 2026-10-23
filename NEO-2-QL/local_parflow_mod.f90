module local_parflow_mod
    ! Local parallel flow on the flux surface for the standard QL drives.
    ! The ripple solver stores, per RK step, the same parallel-flow moment
    ! that solve_eqs sums into qflux(2,:). The flux-surface average of the
    ! written fields reproduces D3j_AX/NA and ParFlow_AX/NA_spec.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use hdf5_tools, only: hid_t, h5_create, h5_add, h5_close
    implicit none
    private

    logical, public :: lsw_local_parflow = .false.
    public :: store_local_parflow, write_local_parflow

    ! NEO-2 source k = 1, 2, 3 drives A1, A3, A2 (ntv_mod ind_map).
    integer, parameter :: source_of_drive(3) = [1, 3, 2]
    real(dp), allocatable :: phi(:), theta(:), bhat(:), h_phi(:), step_weight(:)
    ! (step, source k, drive species, flux species)
    complex(dp), allocatable :: moment_ax(:, :, :, :), moment_na(:, :, :, :)

contains

    subroutine store_local_parflow(axisymmetric, phi_mfl, theta_mfl, bhat_mfl, &
                                   h_phi_mfl, weight, moment)
        logical, intent(in) :: axisymmetric
        real(dp), intent(in) :: phi_mfl(:), theta_mfl(:), bhat_mfl(:), &
                                h_phi_mfl(:), weight(:)
        complex(dp), intent(in) :: moment(:, :, :, :)

        phi = phi_mfl
        theta = theta_mfl
        bhat = bhat_mfl
        h_phi = h_phi_mfl
        step_weight = weight
        if (axisymmetric) then
            moment_ax = moment
        else
            moment_na = moment
        end if
    end subroutine store_local_parflow

    subroutine write_local_parflow(coef_ax, coef_na, forces, have_forces, &
                                   volume, z, density, charge, bref, n)
        ! coef(j, b, a): Dijab(3, j) per unit qflux(2, ind_map(j)) from
        ! compute_Dijab, so the conversion is the existing one.
        real(dp), intent(in) :: coef_ax(:, 0:, 0:), coef_na(:, 0:, 0:)
        real(dp), intent(in) :: forces(:, 0:) ! (A1, A2, A3; drive species)
        logical, intent(in) :: have_forces
        real(dp), intent(in) :: volume, z(0:), density(0:), charge, bref
        integer, intent(in) :: n
        integer(hid_t) :: file_id
        real(dp), allocatable :: weight(:)

        if (.not. allocated(moment_ax)) return
        weight = step_weight/(bhat*h_phi)
        call h5_create('local_parallel_flow.h5', file_id)
        call h5_add(file_id, 'phi_mfl', phi, lbound(phi), ubound(phi), unit='rad', &
                    comment='Boozer phi along the field line')
        call h5_add(file_id, 'theta', theta, lbound(theta), ubound(theta), unit='rad', &
                    comment='Boozer theta=theta_beg+iota*(phi-phi_beg), not reduced mod 2pi')
        call h5_add(file_id, 'bhat', bhat, lbound(bhat), ubound(bhat), unit='1', &
                    comment='B/Bref')
        call h5_add(file_id, 'Bref', bref, unit='G')
        call h5_add(file_id, 'surface_weight', weight, lbound(weight), &
                    ubound(weight), unit='cm', &
                    comment='RK weight dphi/(bhat*h^phi); <x>=sum(surface_weight*x)/surface_volume')
        call h5_add(file_id, 'surface_volume', volume, unit='cm', &
                    comment='y(6) of rhs_kin, the normalisation of D3j and ParFlow')
        call h5_add(file_id, 'n', n, comment='NA field = Re[value*exp(i*n*(phi-phi_mfl))]')
        if (have_forces) call h5_add(file_id, 'A_forces', forces, lbound(forces), &
                                     ubound(forces), comment='(A1, A2, A3; species)')
        call write_problem('AX', moment_ax, coef_ax)
        if (allocated(moment_na)) call write_problem('NA', moment_na, coef_na)
        call h5_close(file_id)

    contains

        subroutine write_problem(suffix, moment, coef)
            character(len=2), intent(in) :: suffix
            complex(dp), intent(in) :: moment(:, :, 0:, 0:)
            real(dp), intent(in) :: coef(:, 0:, 0:)
            complex(dp), allocatable :: d3(:, :, :, :), vparb(:, :), jpar(:)
            integer :: j, a, b

            allocate(d3, mold=moment)
            do a = lbound(moment, 4), ubound(moment, 4)
                do b = lbound(moment, 3), ubound(moment, 3)
                    do j = 1, 3
                        d3(:, j, b, a) = coef(j, b, a)*volume*bhat*h_phi &
                                         *moment(:, source_of_drive(j), b, a)
                    end do
                end do
            end do
            call h5_add(file_id, 'D3j_'//suffix//'_real', real(d3, dp), lbound(d3), &
                        ubound(d3), comment='(step, j, drive species, flux species)')
            call h5_add(file_id, 'D3j_'//suffix//'_imag', aimag(d3), lbound(d3), ubound(d3))
            if (.not. have_forces) return
            ! Same contraction as compute_ParFlow: -sum_b sum_j D3j*Aj.
            allocate(vparb(size(d3, 1), lbound(d3, 4):ubound(d3, 4)))
            vparb = (0.0_dp, 0.0_dp)
            do a = lbound(d3, 4), ubound(d3, 4)
                do b = lbound(d3, 3), ubound(d3, 3)
                    vparb(:, a) = vparb(:, a) - matmul(d3(:, :, b, a), forces(:, b))
                end do
            end do
            jpar = charge*matmul(vparb, z*density)/(bref*bhat)
            call h5_add(file_id, 'VparB_'//suffix//'_real', real(vparb, dp), &
                        lbound(vparb), ubound(vparb), unit='cm*G/s', &
                        comment='(step, species); positive along B')
            call h5_add(file_id, 'VparB_'//suffix//'_imag', aimag(vparb), &
                        lbound(vparb), ubound(vparb), unit='cm*G/s')
            call h5_add(file_id, 'j_par_'//suffix//'_real', real(jpar, dp), &
                        lbound(jpar), ubound(jpar), unit='statA/cm^2', &
                        comment='sum_a e_a n_a V_par,a, positive along B')
            call h5_add(file_id, 'j_par_'//suffix//'_imag', aimag(jpar), &
                        lbound(jpar), ubound(jpar), unit='statA/cm^2')
        end subroutine write_problem
    end subroutine write_local_parflow
end module local_parflow_mod
