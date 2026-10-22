module phi1_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private

    ! Fixed namelist storage avoids allocation when Phi_1 is disabled.
    integer, parameter :: max_phi1_m = 128, unset_m = huge(0)
    integer, public :: isw_phi1 = 0, phi1_num_m = 0
    integer, public :: phi1_m(max_phi1_m) = unset_m
    real(dp), public :: phi1_qn_rtol = 1.0e-8_dp, phi1_rcond_min = 1.0e-12_dp
    public :: phi1_validate_config, phi1_basis, phi1_fieldline_weights
    public :: phi1_gram, phi1_assemble_qn, phi1_solve_qn
    public :: phi1_exb_band, phi1_density, phi1_write_output

    type, public :: phi1_result_t
        integer :: n
        real(dp) :: residual, rcond, chi, qn_flux_check, boltzmann_flux
        real(dp) :: theta_beg, phi_beg
        complex(dp), allocatable :: phi(:), gram(:,:), matrix(:,:), rhs(:)
        complex(dp), allocatable :: density_dB(:,:), response(:,:,:)
        complex(dp), allocatable :: sampled_phi(:), n1(:,:)
        real(dp), allocatable :: grid(:), weights(:), gamma_B(:), gamma_E(:)
        real(dp), allocatable :: energy_B(:), energy_E(:), torque(:)
        real(dp), allocatable :: collision_norm(:,:)
    end type phi1_result_t
    type(phi1_result_t), public, save :: phi1_result

    interface
        subroutine zgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import dp
            integer, intent(in) :: n, nrhs, lda, ldb
            integer, intent(out) :: ipiv(*), info
            complex(dp), intent(inout) :: a(lda,*), b(ldb,*)
        end subroutine zgesv

        subroutine zgecon(norm, n, a, lda, anorm, rcond, work, rwork, info)
            import dp
            character, intent(in) :: norm
            integer, intent(in) :: n, lda
            complex(dp), intent(in) :: a(lda,*)
            real(dp), intent(in) :: anorm
            real(dp), intent(out) :: rcond, rwork(*)
            complex(dp), intent(out) :: work(*)
            integer, intent(out) :: info
        end subroutine zgecon
    end interface

contains

    subroutine phi1_validate_config()
        integer :: j

        if (isw_phi1 == 0) return
        if (isw_phi1 /= 1) error stop 'isw_phi1 must be 0 or 1'
        if (phi1_num_m < 1 .or. phi1_num_m > max_phi1_m) &
            error stop 'phi1_num_m must be in 1:128'
        if (any(phi1_m(1:phi1_num_m) == unset_m)) &
            error stop 'phi1_m requires phi1_num_m entries'
        do j = 2, phi1_num_m
            if (any(phi1_m(1:j-1) == phi1_m(j))) &
                error stop 'phi1_m entries must be distinct'
        end do
        if (.not. ieee_is_finite(phi1_qn_rtol)) &
            error stop 'phi1_qn_rtol must be finite'
        if (phi1_qn_rtol <= 0.0_dp .or. phi1_qn_rtol >= 1.0_dp) &
            error stop 'phi1_qn_rtol must be in (0,1)'
        if (.not. ieee_is_finite(phi1_rcond_min)) &
            error stop 'phi1_rcond_min must be finite'
        if (phi1_rcond_min <= 0.0_dp .or. phi1_rcond_min >= 1.0_dp) &
            error stop 'phi1_rcond_min must be in (0,1)'
    end subroutine phi1_validate_config

    pure subroutine phi1_basis(phi, theta_beg, phi_beg, iota, m, n, psi)
        real(dp), intent(in) :: phi(:), theta_beg, phi_beg, iota
        integer, intent(in) :: m(:), n
        complex(dp), intent(out) :: psi(size(phi),size(m))
        real(dp) :: theta
        integer :: i, j

        do j = 1, size(m)
            do i = 1, size(phi)
                theta = theta_beg + iota*(phi(i) - phi_beg)
                psi(i,j) = exp(cmplx(0.0_dp, m(j)*theta + n*phi(i), dp))
            end do
        end do
    end subroutine phi1_basis

    ! Trapezoidal nodal weights for int(dl/Bhat) = int(dphi/|h^phi|/Bhat).
    ! The caller supplies the complete selected field-line interval, endpoints
    ! included. No velocity or streaming-branch factors belong in these weights.
    pure subroutine phi1_fieldline_weights(phi, bhat, h_phi, weights, info)
        real(dp), intent(in) :: phi(:), bhat(:), h_phi(:)
        real(dp), intent(out) :: weights(size(phi))
        integer, intent(out) :: info
        real(dp) :: step
        integer :: i

        info = -1
        weights = 0.0_dp
        if (size(phi) < 2) return
        if (size(bhat) /= size(phi) .or. size(h_phi) /= size(phi)) return
        if (.not. all(ieee_is_finite(phi))) return
        if (.not. all(ieee_is_finite(bhat))) return
        if (.not. all(ieee_is_finite(h_phi))) return
        if (any(bhat <= 0.0_dp) .or. any(abs(h_phi) <= 0.0_dp)) return
        do i = 1, size(phi)-1
            step = phi(i+1) - phi(i)
            if (step <= 0.0_dp) return
            weights(i) = weights(i) + 0.5_dp*step
            weights(i+1) = weights(i+1) + 0.5_dp*step
        end do
        weights = weights/(bhat*abs(h_phi))
        weights = weights/sum(weights)
        info = 0
    end subroutine phi1_fieldline_weights

    ! Also accepts the solver's own positive nodal surface-quadrature weights.
    pure subroutine phi1_gram(psi, weights, gram, info)
        complex(dp), intent(in) :: psi(:,:)
        real(dp), intent(in) :: weights(:)
        complex(dp), intent(out) :: gram(size(psi,2),size(psi,2))
        integer, intent(out) :: info
        real(dp) :: norm
        integer :: l, j

        info = -1
        gram = 0.0_dp
        if (size(weights) /= size(psi,1)) return
        if (.not. all(ieee_is_finite(weights))) return
        if (any(weights < 0.0_dp)) return
        norm = sum(weights)
        if (norm <= 0.0_dp) return
        do j = 1, size(psi,2)
            do l = 1, size(psi,2)
                gram(l,j) = sum(weights*conjg(psi(:,l))*psi(:,j))/norm
            end do
        end do
        info = 0
    end subroutine phi1_gram

    ! W4 density-projection hook: dB(l,a) [cm^-3], response(l,j,a)
    ! [cm^-3/statvolt], after force contraction. Charges are signed cgs charges,
    ! density is cm^-3 and temperature erg. Include adiabatic species with zero
    ! dB and response but nonzero density and temperature in chi.
    pure subroutine phi1_assemble_qn(charge, density, temperature, gram, &
                                    response, dB, matrix, rhs, chi, info)
        real(dp), intent(in) :: charge(:), density(:), temperature(:)
        complex(dp), intent(in) :: gram(:,:), response(:,:,:), dB(:,:)
        complex(dp), intent(out) :: matrix(size(gram,1),size(gram,2))
        complex(dp), intent(out) :: rhs(size(gram,1))
        real(dp), intent(out) :: chi
        integer, intent(out) :: info
        integer :: a, k, ns

        info = -1
        matrix = 0.0_dp
        rhs = 0.0_dp
        chi = 0.0_dp
        ns = size(charge)
        k = size(gram,1)
        if (ns < 1 .or. k < 1) return
        if (size(density) /= ns .or. size(temperature) /= ns) return
        if (size(gram,2) /= k) return
        if (size(response,1) /= k .or. size(response,2) /= k) return
        if (size(response,3) /= ns) return
        if (size(dB,1) /= k .or. size(dB,2) /= ns) return
        if (.not. all(ieee_is_finite(charge))) return
        if (.not. all(ieee_is_finite(density))) return
        if (.not. all(ieee_is_finite(temperature))) return
        if (any(density < 0.0_dp) .or. any(temperature <= 0.0_dp)) return
        chi = sum(charge**2*density/temperature)
        matrix = -chi*gram
        do a = 1, ns
            matrix = matrix + charge(a)*response(:,:,a)
            rhs = rhs - charge(a)*dB(:,a)
        end do
        info = 0
    end subroutine phi1_assemble_qn

    ! info: 0 success, >0 singular LU pivot, -1 invalid input,
    ! -2 ill-conditioned, -3 excessive residual. Phi is zero on rejection.
    subroutine phi1_solve_qn(matrix, rhs, rtol, rcond_min, phi, residual, &
                             rcond, info)
        complex(dp), intent(in) :: matrix(:,:), rhs(:)
        real(dp), intent(in) :: rtol, rcond_min
        complex(dp), intent(out) :: phi(size(rhs))
        real(dp), intent(out) :: residual, rcond
        integer, intent(out) :: info
        complex(dp), allocatable :: lu(:,:), solution(:,:), work(:), defect(:)
        real(dp), allocatable :: rwork(:)
        integer, allocatable :: pivots(:)
        real(dp) :: anorm, scale
        integer :: k, l, j

        info = -1
        phi = 0.0_dp
        residual = huge(1.0_dp)
        rcond = 0.0_dp
        k = size(rhs)
        if (k < 1) return
        if (size(matrix,1) /= k .or. size(matrix,2) /= k) return
        if (.not. all(ieee_is_finite(real(matrix,dp)))) return
        if (.not. all(ieee_is_finite(aimag(matrix)))) return
        if (.not. all(ieee_is_finite(real(rhs,dp)))) return
        if (.not. all(ieee_is_finite(aimag(rhs)))) return
        if (.not. ieee_is_finite(rtol)) return
        if (.not. ieee_is_finite(rcond_min)) return
        if (rtol <= 0.0_dp .or. rtol >= 1.0_dp) return
        if (rcond_min <= 0.0_dp .or. rcond_min >= 1.0_dp) return
        allocate(lu(k,k), solution(k,1), pivots(k), work(2*k), rwork(2*k), &
            defect(k))
        lu = matrix
        solution(:,1) = rhs
        ! LAPACK's complex 1-norm uses |Re(z)| + |Im(z)|.
        anorm = maxval(sum(abs(real(matrix,dp)) + abs(aimag(matrix)), dim=1))
        call zgesv(k, 1, lu, k, pivots, solution, k, info)
        if (info /= 0) return
        call zgecon('1', k, lu, k, anorm, rcond, work, rwork, info)
        if (info /= 0) return
        if (.not. ieee_is_finite(rcond) .or. rcond < rcond_min) then
            info = -2
            return
        end if
        ! Normwise backward residual, including a well-defined zero RHS.
        scale = anorm*sum(abs(solution(:,1))) + sum(abs(rhs))
        defect = -rhs
        do j = 1, k
            do l = 1, k
                defect(l) = defect(l) + matrix(l,j)*solution(j,1)
            end do
        end do
        residual = sum(abs(defect))
        if (scale > 0.0_dp) residual = residual/scale
        if (.not. ieee_is_finite(residual) .or. residual > rtol) then
            info = -3
            return
        end if
        phi = solution(:,1)
    end subroutine phi1_solve_qn

    ! D7: both streaming branches have this sign. rho is signed, in cm.
    pure subroutine phi1_exb_band(psi, m, n, btheta, bphi, iota, bref, &
                                 psi_prime, vthermal, rho, bhat, hphi, &
                                 lambda_hi, lambda_lo, us, source)
        complex(dp), intent(in) :: psi
        integer, intent(in) :: m, n
        real(dp), intent(in) :: btheta, bphi, iota, bref, psi_prime
        real(dp), intent(in) :: vthermal, rho, bhat, hphi
        real(dp), intent(in) :: lambda_hi(:), lambda_lo(:)
        complex(dp), intent(out) :: us, source(size(lambda_hi))

        us = cmplx(0.0_dp,2.9979e10_dp,dp)*(n*btheta-m*bphi)*psi &
            /(bref*psi_prime*(iota*btheta+bphi))
        source = -2.0_dp*us*(lambda_hi-lambda_lo) &
            /(vthermal*rho*bhat*hphi)
    end subroutine phi1_exb_band

    ! D8 local density or energy moment. Negative-sigma bands are reversed.
    pure function phi1_density(y, pitch_bra, speed_weight, hphi) result(value)
        complex(dp), intent(in) :: y(:,:)
        real(dp), intent(in) :: pitch_bra(:), speed_weight(:), hphi
        complex(dp) :: value
        integer :: m, nb

        nb = size(pitch_bra)
        value = 0.0_dp
        do m = 1, size(speed_weight)
            value = value + speed_weight(m)*sum(pitch_bra &
                *(y(1:nb,m)+y(2*nb:nb+1:-1,m)))
        end do
        value = 4.0_dp*hphi*value
    end function phi1_density

    subroutine phi1_write_output(file_id)
        use hdf5_tools, only: hid_t, h5_define_group, h5_close_group, h5_add
        use collisionality_mod, only: species_tag
        integer(hid_t), intent(in) :: file_id
        integer(hid_t) :: group

        if (isw_phi1 == 0) return
        if (.not. allocated(phi1_result%phi)) &
            error stop 'Phi_1 output requested without a closed solution'
        call h5_define_group(file_id, 'phi1', group)
        call h5_add(group, 'm', phi1_m(1:phi1_num_m), [1], [phi1_num_m])
        call h5_add(group, 'n', phi1_result%n)
        call h5_add(group, 'normalization_version', 1)
        call h5_add(group, 'species_tag', species_tag, &
            lbound(species_tag), ubound(species_tag))
        call h5_add(group, 'theta_beg', phi1_result%theta_beg)
        call h5_add(group, 'phi_beg', phi1_result%phi_beg)
        call vector('Phi_re', real(phi1_result%phi,dp), 'statvolt')
        call vector('Phi_im', aimag(phi1_result%phi), 'statvolt')
        call matrix('Gram_re', real(phi1_result%gram,dp))
        call matrix('Gram_im', aimag(phi1_result%gram))
        call matrix('M_re', real(phi1_result%matrix,dp))
        call matrix('M_im', aimag(phi1_result%matrix))
        call vector('rhs_re', real(phi1_result%rhs,dp), 'statC cm^-3')
        call vector('rhs_im', aimag(phi1_result%rhs), 'statC cm^-3')
        call h5_add(group, 'residual', phi1_result%residual)
        call h5_add(group, 'rcond', phi1_result%rcond)
        call h5_add(group, 'chi', phi1_result%chi)
        call vector('phi_mfl', phi1_result%grid, 'rad')
        call vector('weights', phi1_result%weights, '1')
        call vector('Phi1_re', real(phi1_result%sampled_phi,dp), 'statvolt')
        call vector('Phi1_im', aimag(phi1_result%sampled_phi), 'statvolt')
        call matrix('n1_re', real(phi1_result%n1,dp), 'cm^-3')
        call matrix('n1_im', aimag(phi1_result%n1), 'cm^-3')
        call matrix('density_dB_re', real(phi1_result%density_dB,dp), 'cm^-3')
        call matrix('density_dB_im', aimag(phi1_result%density_dB), 'cm^-3')
        call h5_add(group, 'density_response_re', real(phi1_result%response,dp), &
            lbound(phi1_result%response), ubound(phi1_result%response), &
            unit='cm^-3/statvolt')
        call h5_add(group, 'density_response_im', aimag(phi1_result%response), &
            lbound(phi1_result%response), ubound(phi1_result%response), &
            unit='cm^-3/statvolt')
        call vector('Gamma_magnetic', phi1_result%gamma_B, 'cm^-2 s^-1')
        call vector('Gamma_ExB', phi1_result%gamma_E, 'cm^-2 s^-1')
        call vector('Gamma_total', phi1_result%gamma_B+phi1_result%gamma_E, &
            'cm^-2 s^-1')
        call vector('Q_magnetic', phi1_result%energy_B, 'erg cm^-2 s^-1')
        call vector('Q_ExB', phi1_result%energy_E, 'erg cm^-2 s^-1')
        call vector('Q_total', phi1_result%energy_B+phi1_result%energy_E, &
            'erg cm^-2 s^-1')
        call vector('torque_species', phi1_result%torque, 'dyn cm^-2')
        call h5_add(group, 'torque_total', sum(phi1_result%torque), &
            comment='Torque from the requested ftilde*v_g1 product; D9 ordering applies')
        call h5_add(group, 'qn_flux_check', phi1_result%qn_flux_check)
        call h5_add(group, 'chi_Phi1_Ur', phi1_result%boltzmann_flux)
        call matrix('collision_source_norm', phi1_result%collision_norm)
        call h5_close_group(group)

    contains

        subroutine vector(name, value, unit)
            character(*), intent(in) :: name, unit
            real(dp), intent(in) :: value(:)
            call h5_add(group, name, value, [1], [size(value)], unit=unit)
        end subroutine vector

        subroutine matrix(name, value, unit)
            character(*), intent(in) :: name
            real(dp), intent(in) :: value(:,:)
            character(*), optional, intent(in) :: unit
            call h5_add(group, name, value, [1,1], shape(value), unit=unit)
        end subroutine matrix
    end subroutine phi1_write_output

end module phi1_mod
