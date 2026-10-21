program test_covariant_curl_leaf
    use iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use neo_covariant_curl_generated, only: covariant_curl, &
        covariant_curl_jvp, covariant_curl_vjp
    implicit none
    real(dp) :: x(3), jet(12), direction(12), weight(3), bar(12)
    real(dp) :: actual(3), expected(3), jac, hc_plus(3), hc_minus(3)
    real(dp) :: hc_half_plus(3), hc_half_minus(3), dh(3,3), h, dx(3)
    real(dp) :: tangent(3), plus(3), minus(3), error
    integer :: sign_j, axis, k, step

    x = [0.43_dp, 0.71_dp, 1.13_dp]
    direction = [(sin(real(k,dp))/10, k=1,12)]
    weight = [0.3_dp, -0.2_dp, 0.7_dp]
    do sign_j = -1, 1, 2
        jac = sign_j*2.3_dp
        call point_jet(x, jac, jet)
        call evaluate(jet, actual)
        do axis = 1, 3
            dx = 0.0_dp
            h = 1.0e-3_dp
            dx(axis) = h
            call independent_field(x+dx, hc_plus)
            call independent_field(x-dx, hc_minus)
            call independent_field(x+dx/2, hc_half_plus)
            call independent_field(x-dx/2, hc_half_minus)
            dh(:,axis) = (4*(hc_half_plus-hc_half_minus)/h &
                - (hc_plus-hc_minus)/(2*h))/3
        end do
        expected = [dh(3,2)-dh(2,3), dh(1,3)-dh(3,1), &
            dh(2,1)-dh(1,2)]/jac
        if (any(abs(actual-expected) > &
            1.0e-8_dp*max(abs(expected),1.0e-8_dp))) &
            error stop 'independent covariant-field curl leaf'
        if (any(abs(actual) < 1.0e-8_dp)) error stop 'degenerate curl fixture'
        call covariant_curl_jvp(jet(1),jet(2),jet(3),jet(4),jet(5),jet(6), &
            jet(7),jet(8),jet(9),jet(10),jet(11),jet(12), &
            direction(1),direction(2),direction(3),direction(4), &
            direction(5),direction(6),direction(7),direction(8), &
            direction(9),direction(10),direction(11),direction(12), &
            tangent(1),tangent(2),tangent(3))
        if (.not. all(ieee_is_finite(tangent))) error stop 'nonfinite curl JVP'
        do step = 1, 2
            h = 1.0e-5_dp/2**(step-1)
            call evaluate(jet+h*direction, plus)
            call evaluate(jet-h*direction, minus)
            expected = (plus-minus)/(2*h)
            error = maxval(abs(tangent-expected)/max(abs(tangent),1.0e-8_dp))
            if (error > 1.0e-8_dp) error stop 'independent curl JVP'
        end do
        call covariant_curl_vjp(jet(1),jet(2),jet(3),jet(4),jet(5),jet(6), &
            jet(7),jet(8),jet(9),jet(10),jet(11),jet(12), &
            weight(1),weight(2),weight(3),bar(1),bar(2),bar(3),bar(4), &
            bar(5),bar(6),bar(7),bar(8),bar(9),bar(10),bar(11),bar(12))
        if (.not. all(ieee_is_finite(bar))) error stop 'nonfinite curl VJP'
        if (abs(dot_product(weight,tangent)-dot_product(bar,direction)) > &
            1.0e-10_dp*max(1.0_dp,abs(dot_product(weight,tangent)))) &
            error stop 'curl JVP/VJP pairing'
    end do
    print *, 'All covariant curl leaf tests passed!'
contains
    subroutine independent_field(position, hc)
        real(dp), intent(in) :: position(3)
        real(dp), intent(out) :: hc(3)
        real(dp) :: s, p, t, magnitude
        s=position(1); p=position(2); t=position(3)
        magnitude = 2 + 0.2_dp*s + 0.1_dp*sin(t) + 0.08_dp*cos(p)
        hc(1) = 100*(0.3_dp*s+0.2_dp*sin(t)+0.15_dp*cos(p))/magnitude
        hc(2) = 100*(3+0.1_dp*s+0.02_dp*s**2)/magnitude
        hc(3) = 100*(0.5_dp+0.07_dp*s+0.03_dp*s**2)/magnitude
    end subroutine independent_field
    subroutine point_jet(position, jacobian, v)
        real(dp), intent(in) :: position(3), jacobian
        real(dp), intent(out) :: v(12)
        real(dp) :: s, p, t
        s=position(1); p=position(2); t=position(3)
        v = [2+0.2_dp*s+0.1_dp*sin(t)+0.08_dp*cos(p), &
            0.3_dp*s+0.2_dp*sin(t)+0.15_dp*cos(p), &
            3+0.1_dp*s+0.02_dp*s**2, 0.5_dp+0.07_dp*s+0.03_dp*s**2, &
            0.2_dp, 0.1_dp*cos(t), -0.08_dp*sin(p), 0.2_dp*cos(t), &
            -0.15_dp*sin(p), 0.07_dp+0.06_dp*s, 0.1_dp+0.04_dp*s, jacobian]
    end subroutine point_jet
    subroutine evaluate(v, curl)
        real(dp), intent(in) :: v(12)
        real(dp), intent(out) :: curl(3)
        call covariant_curl(v(1),v(2),v(3),v(4),v(5),v(6),v(7),v(8), &
            v(9),v(10),v(11),v(12),curl(1),curl(2),curl(3))
        if (.not. all(ieee_is_finite(curl))) error stop 'nonfinite curl leaf'
    end subroutine evaluate
end program test_covariant_curl_leaf
