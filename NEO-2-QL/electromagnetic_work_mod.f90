! Direct magnetic work for the native physical-velocity f_11 response.
! This observable does not use a collision-momentum conservation identity.
module electromagnetic_work_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private
    public :: work_back_geometry, work_pitch_convolution
    public :: compute_electromagnetic_work
    public :: work_geometry_supported

contains

    ! Native lab 8/9 returns (s,phi,theta), with the opposite Jacobian
    ! sign to (s,theta,phi). Positive toroidal field has sqrtg < 0 here.
    pure logical function work_geometry_supported(coordinates, laboratory, &
            sqrtg, hphi, bphi_hat, iota, psi_hat, g, dl_over_b) result(supported)
        integer, intent(in) :: coordinates, laboratory
        real(dp), intent(in) :: sqrtg, hphi, bphi_hat, iota, psi_hat
        real(dp), intent(in) :: g, dl_over_b

        supported = coordinates /= 0 .and. &
            (laboratory == 8 .or. laboratory == 9) .and. &
            sqrtg < 0.0_dp .and. hphi > 0.0_dp .and. &
            bphi_hat > 0.0_dp .and. psi_hat > 0.0_dp .and. &
            abs(iota) > 0.0_dp .and. g > 0.0_dp .and. dl_over_b > 0.0_dp
    end function work_geometry_supported

    pure elemental function work_back_geometry(scale, m_phi, bn, bhat) result(back)
        real(dp), intent(in) :: scale, bhat
        integer, intent(in) :: m_phi
        complex(dp), intent(in) :: bn
        complex(dp) :: back

        back = scale*cmplx(0.0_dp, 1.0_dp, dp)*m_phi*bn*bhat
    end function work_back_geometry

    pure elemental function work_pitch_convolution(polpow0, polpow2, &
                                                   bhat, hphi) result(bra)
        real(dp), intent(in) :: polpow0, polpow2, bhat, hphi
        complex(dp) :: bra
        complex(dp) :: coefficient

        coefficient = cmplx(1.0_dp, 0.0_dp, dp)/bhat**2/hphi
        bra = (polpow0 + polpow2)*coefficient
    end function work_pitch_convolution

    ! Raw qflux already includes Re(back_bra^H SOL)/2, native energy moments,
    ! pitch reconstruction and field-line quadrature. Its source channels are
    ! (A1,A3,A2). Reconstruct f_11/F_M with source_weight, not beta_out from
    ! the transport-output normalization. CGS input gives erg/cm^3 torque.
    pure subroutine compute_electromagnetic_work(qflux, mass, charge_number, &
            density, temperature, dn_ds, dT_ds, er, epar_over_b2, bref, g, &
            dl_over_b, scale, torque, source_weight)
        real(dp), intent(in) :: qflux(1:,1:,0:,0:)
        real(dp), intent(in) :: mass(0:), charge_number(0:), density(0:)
        real(dp), intent(in) :: temperature(0:), dn_ds(0:), dT_ds(0:)
        real(dp), intent(in) :: er, epar_over_b2, bref, g, dl_over_b, scale
        real(dp), intent(out) :: torque(0:), source_weight(1:,0:)
        real(dp), parameter :: c = 2.9979e10_dp, e = 4.8032e-10_dp
        real(dp) :: vth, rho, a1, a2, a3, response
        integer :: a, b

        do b = 0, size(mass) - 1
            vth = sqrt(2.0_dp*temperature(b)/mass(b))
            rho = mass(b)*vth*c/(charge_number(b)*e*bref)
            a2 = g*dT_ds(b)/temperature(b)
            a1 = g*dn_ds(b)/density(b) - 1.5_dp*a2 &
                - charge_number(b)*e*er/temperature(b)
            a3 = charge_number(b)*e*epar_over_b2/temperature(b)
            source_weight(:,b) = [rho*a1/g, bref*a3, rho*a2/g]
        end do

        do a = 0, size(mass) - 1
            response = 0.0_dp
            do b = 0, size(mass) - 1
                response = response + sum(source_weight(:,b)*qflux(1,:,a,b))
            end do
            torque(a) = -2.0_dp*density(a)*temperature(a)*response &
                /(scale*dl_over_b)
        end do
    end subroutine compute_electromagnetic_work

end module electromagnetic_work_mod
