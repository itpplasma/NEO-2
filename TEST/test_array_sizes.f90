program test_array_sizes
    use nrtype, only: I4B, DP
    implicit none
    
    interface
        subroutine splinecof3_original_dense(x, y, c1, cn, lambda1, indx, sw1, sw2, &
             a, b, c, d, m, f)
            use nrtype, only: I4B, DP
            real(DP),                   intent(inout) :: c1, cn
            real(DP),     dimension(:), intent(in)    :: x, y, lambda1
            integer(I4B), dimension(:), intent(in)    :: indx
            real(DP),     dimension(:), intent(out)   :: a, b, c, d
            integer(I4B),               intent(in)    :: sw1, sw2
            real(DP),                   intent(in)    :: m
            interface
               function f(x,m)
                 use nrtype, only : DP
                 implicit none
                 real(DP), intent(in) :: x, m
                 real(DP)             :: f
               end function f
            end interface
        end subroutine splinecof3_original_dense
    end interface
    
    integer(I4B), parameter :: n = 5
    real(DP) :: x(n), y(n)
    integer(I4B) :: indx(n)
    real(DP) :: lambda1(n)
    real(DP) :: a_test(n), b_test(n), c_test(n), d_test(n)
    real(DP) :: c1, cn, m
    integer(I4B) :: sw1, sw2, i
    
    ! Initialize test data
    x = [0.0_DP, 1.0_DP, 2.0_DP, 3.0_DP, 4.0_DP]
    y = [0.0_DP, 1.0_DP, 4.0_DP, 9.0_DP, 16.0_DP]
    indx = [(i, i=1,n)]
    lambda1 = 1.0_DP
    c1 = 0.0_DP
    cn = 0.0_DP
    sw1 = 2
    sw2 = 4
    m = 0.0_DP
    
    ! Initialize coefficient arrays with sentinel values
    a_test = -999.0_DP
    b_test = -999.0_DP
    c_test = -999.0_DP
    d_test = -999.0_DP
    
    write(*,'(A)') 'Testing array sizes for original implementation'
    write(*,'(A,I0)') 'Number of data points (n): ', n
    write(*,'(A,I0)') 'Number of intervals (n-1): ', n-1
    write(*,'(A)') ''
    
    ! Call original implementation
    call splinecof3_original_dense(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                                  a_test, b_test, c_test, d_test, m, test_function)
    
    ! Check which elements were written
    write(*,'(A)') 'Coefficient array values after spline calculation:'
    do i = 1, n
        write(*,'(A,I0,A,4F10.3)') 'i=', i, ': a,b,c,d = ', a_test(i), b_test(i), c_test(i), d_test(i)
    end do
    
    write(*,'(A)') ''
    write(*,'(A)') 'Analysis:'
    if (abs(a_test(n) + 999.0_DP) < 1.0e-10) then
        write(*,'(A)') 'Original implementation outputs n-1 coefficients (correct)'
    else
        write(*,'(A)') 'Original implementation outputs n coefficients (includes extra element)'
    end if
    
contains

    function test_function(x, m) result(f_val)
        real(DP), intent(in) :: x, m
        real(DP) :: f_val
        f_val = 1.0_DP
    end function test_function

end program test_array_sizes