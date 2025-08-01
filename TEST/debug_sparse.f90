program debug_sparse
    use nrtype, only: I4B, DP
    implicit none
    
    interface
        subroutine splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
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
        end subroutine splinecof3_a
        
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
    
    ! Simple test
    integer(I4B), parameter :: n = 3
    real(DP) :: x(n), y(n)
    integer(I4B) :: indx(2)
    real(DP) :: lambda1(2)
    real(DP) :: a_new(2), b_new(2), c_new(2), d_new(2)
    real(DP) :: a_orig(2), b_orig(2), c_orig(2), d_orig(2)
    real(DP) :: c1, cn, m
    integer(I4B) :: sw1, sw2
    
    ! Simple test case
    x = [0.0_DP, 1.0_DP, 2.0_DP]
    y = [0.0_DP, 1.0_DP, 4.0_DP]  ! x^2
    indx = [1, 3]
    lambda1 = [1.0_DP, 1.0_DP]
    c1 = 0.0_DP
    cn = 0.0_DP
    sw1 = 2  ! Natural
    sw2 = 4  ! Natural
    m = 0.0_DP
    
    write(*,'(A)') 'Testing simple case:'
    write(*,'(A,3F8.3)') 'x = ', x
    write(*,'(A,3F8.3)') 'y = ', y
    
    ! Original
    call splinecof3_original_dense(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                                  a_orig, b_orig, c_orig, d_orig, m, test_function)
    
    ! New
    c1 = 0.0_DP; cn = 0.0_DP  ! Reset
    call splinecof3_a(x, y, c1, cn, lambda1, indx, sw1, sw2, &
                      a_new, b_new, c_new, d_new, m, test_function)
    
    write(*,'(A)') ''
    write(*,'(A)') 'Original coefficients:'
    write(*,'(A,2F12.6)') 'a = ', a_orig
    write(*,'(A,2F12.6)') 'b = ', b_orig
    write(*,'(A,2F12.6)') 'c = ', c_orig
    write(*,'(A,2F12.6)') 'd = ', d_orig
    
    write(*,'(A)') ''
    write(*,'(A)') 'New coefficients:'
    write(*,'(A,2F12.6)') 'a = ', a_new
    write(*,'(A,2F12.6)') 'b = ', b_new
    write(*,'(A,2F12.6)') 'c = ', c_new
    write(*,'(A,2F12.6)') 'd = ', d_new
    
    write(*,'(A)') ''
    write(*,'(A)') 'Differences:'
    write(*,'(A,2E12.5)') 'a diff = ', abs(a_new - a_orig)
    write(*,'(A,2E12.5)') 'b diff = ', abs(b_new - b_orig)
    write(*,'(A,2E12.5)') 'c diff = ', abs(c_new - c_orig)
    write(*,'(A,2E12.5)') 'd diff = ', abs(d_new - d_orig)
    
contains

    function test_function(x, m) result(f_val)
        real(DP), intent(in) :: x, m
        real(DP) :: f_val
        f_val = 1.0_DP  ! Simple weight function
    end function test_function

end program debug_sparse