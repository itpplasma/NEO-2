module helical_response_mod
    use nrtype, only: dp
    use ntv_mod, only: c, e
    implicit none
    private

    character(len=*), parameter, public :: D_CURRENT_UNIT = 'cm^2 G/s'
    public :: dimensional_current_coefficients
    public :: reconstruct_complex_response

contains

    subroutine dimensional_current_coefficients(gamma, temperature, charge_number, &
            gamma_current, d_current)
        real(dp), intent(in) :: gamma(3, 3), temperature, charge_number
        real(dp), intent(out) :: gamma_current(2), d_current(2)
        real(dp) :: scale

        if (temperature <= 0.0_dp) error stop 'helical response temperature must be positive'
        if (charge_number == 0.0_dp) error stop 'helical response charge must be nonzero'

        gamma_current = gamma(3, 1:2)
        scale = 2.0_dp*temperature*c/(charge_number*e)
        d_current = gamma_current*scale
    end subroutine dimensional_current_coefficients

    pure subroutine reconstruct_complex_response(real_run, imaginary_run, baseline, response)
        real(dp), intent(in) :: real_run(:), imaginary_run(:), baseline(:)
        complex(dp), intent(out) :: response(:)

        response = cmplx(real_run - baseline, -(imaginary_run - baseline), kind=dp)
    end subroutine reconstruct_complex_response
end module helical_response_mod
