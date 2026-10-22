module phi1_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    integer, parameter :: max_phi1_m = 128, unset_m = huge(0)
    integer, public :: isw_phi1 = 0, phi1_num_m = 0
    integer, public :: phi1_m(max_phi1_m) = unset_m
    integer, public :: phi1_file_nper = 1
    real(dp), parameter, public :: phi1_qn_rtol = 1e-8_dp, phi1_rcond_min = 1e-12_dp
    public :: phi1_assemble_qn, phi1_solve_qn, phi1_write_bc
    public :: phi1_validate_config, phi1_basis, phi1_density, phi1_write_output
    type, public :: phi1_response_t
        real(dp) :: residual, rcond, chi, qn_flux_check, boltzmann_flux
        complex(dp), allocatable :: phi(:), matrix(:,:), rhs(:)
        complex(dp), allocatable :: sampled_phi(:), n1(:,:), drift(:,:)
        real(dp), allocatable :: gamma_B(:), gamma_E(:), energy_B(:), energy_E(:)
        real(dp), allocatable :: torque(:)
        integer :: n
        real(dp) :: theta_beg, phi_beg
        complex(dp), allocatable :: basis(:,:), gram(:,:)
        complex(dp), allocatable :: density_dB(:,:), response(:,:,:)
        real(dp), allocatable :: grid(:), weights(:), collision_norm(:,:)
    end type
    type(phi1_response_t), public, save :: phi1_result
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
        if (isw_phi1 /= 1 .and. isw_phi1 /= 2) &
            error stop 'isw_phi1 must be 0 (off), 1 (response), or 2 (QN)'
        if (phi1_num_m < 1 .or. phi1_num_m > max_phi1_m) &
            error stop 'phi1_num_m must be in 1:128'
        if (any(phi1_m(1:phi1_num_m) == unset_m)) &
            error stop 'phi1_m requires phi1_num_m entries'
        do j = 2, phi1_num_m
            if (any(phi1_m(1:j-1) == phi1_m(j))) &
                error stop 'phi1_m entries must be distinct'
        end do
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

    ! D4 density-projection contract: dB(l,a) [cm^-3], response(l,j,a)
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

    subroutine phi1_write_output(file_id)
        use hdf5_tools, only: hid_t, h5_define_group, h5_close_group, h5_add
        use collisionality_mod, only: species_tag
        integer(hid_t), intent(in) :: file_id
        integer(hid_t) :: group

        if (isw_phi1 == 0) return
        call h5_define_group(file_id, 'phi1', group)
        call h5_add(group, 'm', phi1_m(1:phi1_num_m), [1], [phi1_num_m])
        call h5_add(group, 'n', phi1_result%n)
        call h5_add(group, 'normalization_version', 1)
        call h5_add(group, 'species_tag', species_tag, &
            lbound(species_tag), ubound(species_tag))
        call h5_add(group, 'theta_beg', phi1_result%theta_beg)
        call h5_add(group, 'phi_beg', phi1_result%phi_beg)
        call matrix('Gram_re', real(phi1_result%gram,dp))
        call matrix('Gram_im', aimag(phi1_result%gram))
        call vector('phi_mfl', phi1_result%grid, 'rad')
        call vector('weights', phi1_result%weights, '1')
        call matrix('density_dB_re', real(phi1_result%density_dB,dp), 'cm^-3')
        call matrix('density_dB_im', aimag(phi1_result%density_dB), 'cm^-3')
        call h5_add(group, 'density_response_re', real(phi1_result%response,dp), &
            lbound(phi1_result%response), ubound(phi1_result%response), &
            unit='cm^-3/statvolt')
        call h5_add(group, 'density_response_im', aimag(phi1_result%response), &
            lbound(phi1_result%response), ubound(phi1_result%response), &
            unit='cm^-3/statvolt')
        call matrix('collision_source_norm', phi1_result%collision_norm)
        if (isw_phi1 == 2) then
            call vector('Phi_re', real(phi1_result%phi,dp), 'statvolt')
            call vector('Phi_im', aimag(phi1_result%phi), 'statvolt')
            call matrix('M_re', real(phi1_result%matrix,dp))
            call matrix('M_im', aimag(phi1_result%matrix))
            call vector('rhs_re', real(phi1_result%rhs,dp), 'statC cm^-3')
            call vector('rhs_im', aimag(phi1_result%rhs), 'statC cm^-3')
            call h5_add(group, 'residual', phi1_result%residual)
            call h5_add(group, 'rcond', phi1_result%rcond)
            call h5_add(group, 'chi', phi1_result%chi)
            call vector('Phi1_re', real(phi1_result%sampled_phi,dp), 'statvolt')
            call vector('Phi1_im', aimag(phi1_result%sampled_phi), 'statvolt')
            call matrix('n1_re', real(phi1_result%n1,dp), 'cm^-3')
            call matrix('n1_im', aimag(phi1_result%n1), 'cm^-3')
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
        end if
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

    ! The ten-column exchange layout is shared with NEO-RT. The background
    ! supplies header geometry, and the perturbation supplies nper/angle rows.
    subroutine phi1_write_bc()
        use neo_control, only: in_file, inp_swi, lab_swi
        use neo_magfie, only: curr_pol_array, curr_tor_array, pprime_array, &
            sqrtg00_array, boozer_iota
        use mag_interface_mod, only: boozer_s
        integer :: unit, j, m0, n0, ns, nper
        real(dp) :: flux, minor, major, current_sign
        character(1024) :: line

        if (inp_swi /= 9 .or. (lab_swi /= 9 .and. lab_swi /= 10)) &
            error stop 'Phi_1 exchange requires the ten-column Boozer input'
        open(newunit=unit,file=trim(in_file),status='old')
        do j = 1, 5
            read(unit,'(a)') line
        end do
        read(unit,*) m0, n0, ns, nper, flux, minor, major
        close(unit)
        current_sign = 1.0_dp
        if (lab_swi == 10) current_sign = -1.0_dp
        nper = phi1_file_nper
        if (mod(phi1_result%n,nper) /= 0) &
            error stop 'Phi_1 toroidal mode is incompatible with perturbation nper'
        open(newunit=unit,file='phi1.bc',status='replace')
        write(unit,'(a)') 'CC quantity: Phi_1 [statvolt] in bmnc/bmns columns', &
            'CC Re[Phi_j exp(i*(m theta+n*nper phi))]; bmns=-Im(Phi_j)', &
            'CC single surface; geometry header copied/interpolated from Boozer input', &
            'CC synthetic or published provenance follows the run inputs', &
            'm0b n0b nsurf nper flux a R'
        write(unit,'(4(i8,1x),3(es25.17,1x))') &
            phi1_num_m-1, 0, 1, nper, flux, minor, major
        write(unit,'(a)') 's iota Jpol/nper Itor pprime sqrtg00', 'units as Boozer input'
        write(unit,'(6(es25.17,1x))') boozer_s, boozer_iota, &
            current_sign*curr_pol_array(1)/(2e-7_dp*nper), &
            current_sign*curr_tor_array(1)/2e-7_dp, &
            pprime_array(1), sqrtg00_array(1)
        write(unit,'(a)') 'm n rmnc rmns zmnc zmns vmnc vmns bmnc bmns'
        do j = 1, phi1_num_m
            write(unit,'(2(i8,1x),8(es25.17,1x))') phi1_m(j), &
                phi1_result%n/nper, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                0.0_dp, 0.0_dp, real(phi1_result%phi(j),dp), -aimag(phi1_result%phi(j))
        end do
        close(unit)
    end subroutine phi1_write_bc
end module phi1_mod
