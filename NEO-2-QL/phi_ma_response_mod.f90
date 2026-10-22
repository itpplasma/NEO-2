module phi_ma_response_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use hdf5_tools, only: hid_t, h5_create, h5_add, h5_close
    use ntv_mod, only: e
    implicit none
    private

    ! Re[Phi_m exp(i(m theta+n phi))], n is the magnetic perturbation's m_phi.
    complex(dp), public :: phi_m = (0.0_dp, 0.0_dp) ! statvolt
    complex(dp), public :: delta_m = (0.0_dp, 0.0_dp) ! cm
    integer, public :: phi_ma_m = 0
    ! One units boundary, matching derivation/L1/mapping.py. Keep NEO's
    ! internal finite-precision charge e, but use exact SI conversions.
    real(dp), parameter, public :: volt_per_statvolt = 299.792458_dp
    real(dp), parameter, public :: ampere_per_statampere = 1.0_dp/2997924580.0_dp
    real(dp), parameter, public :: ampere_m2_per_statampere_cm2 = &
        1.0e4_dp*ampere_per_statampere
    public :: misaligned_potential, write_phi_ma_current, normalize_phi_ma_j_hat

contains

    pure function misaligned_potential(er) result(phi_ma)
        real(dp), intent(in) :: er ! statvolt/cm, E_r = -dPhi_0/dr
        complex(dp) :: phi_ma

        phi_ma = phi_m + er*delta_m
    end function misaligned_potential

    pure subroutine normalize_phi_ma_j_hat(jpar, z, density, temperature, mass, &
                                           a1, potential, er, jhat, valid)
        complex(dp), intent(in) :: jpar(:), potential
        real(dp), intent(in) :: z, density, temperature, mass, a1, er
        complex(dp), intent(out) :: jhat(:)
        logical, intent(out) :: valid
        real(dp) :: vt_neo, vt_contract

        jhat = (0.0_dp, 0.0_dp)
        valid = potential /= (0.0_dp, 0.0_dp) .and. a1 /= 0.0_dp &
                .and. er /= 0.0_dp .and. z /= 0.0_dp &
                .and. density > 0.0_dp .and. temperature > 0.0_dp .and. mass > 0.0_dp
        if (.not. valid) return
        vt_neo = sqrt(2.0_dp*temperature/mass)
        vt_contract = vt_neo/sqrt(2.0_dp)
        jhat = jpar/(z*e*density*vt_contract*a1*potential/er)
    end subroutine normalize_phi_ma_j_hat

    subroutine write_phi_ma_current(phi, theta, bhat, weight, jpar, row, &
                                    phi_ma, er, n, jhat, jhat_valid)
        real(dp), intent(in) :: phi(:), theta(:), bhat(:), weight(:), er
        complex(dp), intent(in) :: jpar(:), row, phi_ma
        integer, intent(in) :: n
        complex(dp), intent(in) :: jhat(:)
        logical, intent(in) :: jhat_valid
        integer(hid_t) :: file_id
        complex(dp) :: average, b_average

        average = sum(weight*jpar)/sum(weight)
        b_average = sum(weight*bhat*jpar)/sum(weight)
        call h5_create('phi_ma_current.h5', file_id)
        call h5_add(file_id, 'phi_mfl', phi, lbound(phi), ubound(phi), unit='rad')
        call h5_add(file_id, 'theta', theta, lbound(theta), ubound(theta), unit='rad')
        call h5_add(file_id, 'bhat', bhat, lbound(bhat), ubound(bhat), unit='1')
        call h5_add(file_id, 'surface_weight', weight, lbound(weight), &
                    ubound(weight), unit='cm', &
                    comment='dphi/(Bhat*abs(h^phi)); RK weights')
        call h5_add(file_id, 'j_par_real', real(jpar, dp), lbound(jpar), &
                    ubound(jpar), unit='statA/cm^2', &
                    comment='Positive along equilibrium B; real part of complex local current')
        call h5_add(file_id, 'j_par_imag', aimag(jpar), lbound(jpar), &
                    ubound(jpar), unit='statA/cm^2', &
                    comment='Implied real field: Re[j_par(theta,phi)]; includes exp(i*n*phi)')
        call h5_add(file_id, 'j_par_line_average_real', real(average, dp), &
                    unit='statA/cm^2')
        call h5_add(file_id, 'j_par_line_average_imag', aimag(average), &
                    unit='statA/cm^2')
        call h5_add(file_id, 'j_par_bhat_line_average_real', real(b_average, dp), &
                    unit='statA/cm^2')
        call h5_add(file_id, 'j_par_bhat_line_average_imag', aimag(b_average), &
                    unit='statA/cm^2')
        ! On this axisymmetric background the full-surface average of a
        ! single n/=0 physical harmonic vanishes, unlike its line average.
        if (n /= 0) average = (0.0_dp, 0.0_dp)
        call h5_add(file_id, 'j_par_fsa_real', real(average, dp), unit='statA/cm^2', &
                    comment='Full toroidal average of the physical single-n current')
        call h5_add(file_id, 'j_par_fsa_imag', aimag(average), unit='statA/cm^2')
        call h5_add(file_id, 'current_row_real', real(row, dp), unit='statA/cm^2', &
                    comment='Existing current-row contraction divided by Bref and line volume')
        call h5_add(file_id, 'current_row_imag', aimag(row), unit='statA/cm^2')
        call h5_add(file_id, 'm', phi_ma_m)
        call h5_add(file_id, 'n', n, comment='Drive: Re[Phi_MA*exp(i*(m*theta+n*phi))]')
        call h5_add(file_id, 'phi_ma_real', real(phi_ma, dp), unit='statvolt')
        call h5_add(file_id, 'phi_ma_imag', aimag(phi_ma), unit='statvolt')
        call h5_add(file_id, 'Er', er, unit='statvolt/cm', &
                    comment='Phi_0 prime = -Er; Phi_MA = Phi_m + Er*Delta_m')
        call h5_add(file_id, 'volt_per_statvolt', volt_per_statvolt)
        call h5_add(file_id, 'ampere_per_statampere', ampere_per_statampere)
        call h5_add(file_id, 'ampere_m2_per_statampere_cm2', &
                    ampere_m2_per_statampere_cm2)
        if (jhat_valid) then
            call h5_add(file_id, 'J_hat_real', real(jhat, dp), lbound(jhat), &
                        ubound(jhat), unit='1', &
                        comment='Single species, A2=0; signed charge; physical Er; sampled phase retained')
            call h5_add(file_id, 'J_hat_imag', aimag(jhat), lbound(jhat), &
                        ubound(jhat), unit='1')
        end if
        call h5_close(file_id)
    end subroutine write_phi_ma_current
end module phi_ma_response_mod
