program test_spline_value_derivatives
    use iso_fortran_env, only: dp => real64
    use spline_mod, only: eva2d, eva2d_fd
    implicit none
    real(dp) :: coefficients(4,4,2,2), x, y, value, gradient(2), expected(2)
    integer :: j
    coefficients=0.0_dp
    coefficients(1,1,1,1)=3.0_dp
    coefficients(2,1,1,1)=2.0_dp
    coefficients(1,2,1,1)=-4.0_dp
    coefficients(2,2,1,1)=0.7_dp
    coefficients(4,1,1,1)=0.11_dp
    coefficients(1,4,1,1)=-0.09_dp
    coefficients(3,3,1,1)=0.13_dp
    do j=1,3
        x=0.23_dp*j; y=0.61_dp/j
        call eva2d(2,2,1,1,x,y,coefficients,value)
        call eva2d_fd(2,2,1,1,x,y,coefficients,gradient)
        if (abs(value-(3+2*x-4*y+0.7_dp*x*y+0.11_dp*x**3 &
            -0.09_dp*y**3+0.13_dp*x**2*y**2))>1.0e-12_dp) &
            error stop 'same coefficient polynomial value'
        expected=[2+0.7_dp*y+0.33_dp*x**2+0.26_dp*x*y**2, &
            -4+0.7_dp*x-0.27_dp*y**2+0.26_dp*x**2*y]
        if(any(abs(gradient-expected)>1.0e-12_dp*max(abs(expected),1.0_dp))) &
            error stop 'same polynomial theta/phi derivative ordering'
    end do
    print *, 'All coefficient-polynomial derivative tests passed!'
end program test_spline_value_derivatives
