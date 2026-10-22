module phi1_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private
    integer, parameter :: max_phi1_m = 128, unset_m = huge(0)
    integer, public :: isw_phi1 = 0, phi1_num_m = 0
    integer, public :: phi1_m(max_phi1_m) = unset_m
    public :: phi1_validate_config, phi1_basis, phi1_density, phi1_write_output
    type, public :: phi1_response_t
        integer :: n
        real(dp) :: theta_beg, phi_beg
        complex(dp), allocatable :: basis(:,:), gram(:,:)
        complex(dp), allocatable :: density_dB(:,:), response(:,:,:)
        real(dp), allocatable :: grid(:), weights(:), collision_norm(:,:)
    end type
    type(phi1_response_t), public, save :: phi1_result
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
