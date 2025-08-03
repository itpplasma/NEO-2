program test_nrutil
    use nrtype
    use nrutil
    implicit none
    
    integer :: test_status
    
    test_status = 0
    
    ! Test arithmetic progression functions
    call test_arth_functions(test_status)
    
    ! Test swap functions
    call test_swap_functions(test_status)
    
    ! Test assert_eq functions
    call test_assert_eq_functions(test_status)
    
    if (test_status == 0) then
        print *, "All tests passed!"
    else
        print *, "Some tests failed. Status:", test_status
        error stop
    end if
    
contains

    subroutine test_arth_functions(status)
        integer, intent(inout) :: status
        real(SP), dimension(5) :: result_r
        real(DP), dimension(5) :: result_d
        integer, dimension(5) :: result_i
        integer :: i
        
        print *, "Testing arth functions..."
        
        ! Test integer arithmetic progression
        result_i = arth(1, 2, 5)  ! Start at 1, increment by 2, 5 elements
        if (any(result_i /= [1, 3, 5, 7, 9])) then
            print *, "FAIL: arth_i test failed"
            print *, "Expected: [1, 3, 5, 7, 9]"
            print *, "Got:", result_i
            status = status + 1
        else
            print *, "PASS: arth_i test"
        end if
        
        ! Test single precision arithmetic progression
        result_r = arth(0.0_SP, 0.5_SP, 5)  ! Start at 0, increment by 0.5, 5 elements
        if (any(abs(result_r - [0.0, 0.5, 1.0, 1.5, 2.0]) > epsilon(result_r) * 10.0)) then
            print *, "FAIL: arth_r test failed"
            print *, "Expected: [0.0, 0.5, 1.0, 1.5, 2.0]"
            print *, "Got:", result_r
            status = status + 1
        else
            print *, "PASS: arth_r test"
        end if
        
        ! Test double precision arithmetic progression
        result_d = arth(1.0_DP, -0.25_DP, 5)  ! Start at 1, decrement by 0.25, 5 elements
        if (any(abs(result_d - [1.0, 0.75, 0.5, 0.25, 0.0]) > epsilon(result_d) * 10.0)) then
            print *, "FAIL: arth_d test failed"
            print *, "Expected: [1.0, 0.75, 0.5, 0.25, 0.0]"
            print *, "Got:", result_d
            status = status + 1
        else
            print *, "PASS: arth_d test"
        end if
        
    end subroutine test_arth_functions
    
    subroutine test_swap_functions(status)
        integer, intent(inout) :: status
        integer :: a_i, b_i
        real(SP) :: a_r, b_r
        real(DP) :: a_d, b_d
        
        print *, "Testing swap functions..."
        
        ! Test integer swap
        a_i = 10
        b_i = 20
        call swap(a_i, b_i)
        if (a_i /= 20 .or. b_i /= 10) then
            print *, "FAIL: swap_i test failed"
            print *, "Expected: a=20, b=10"
            print *, "Got: a=", a_i, ", b=", b_i
            status = status + 1
        else
            print *, "PASS: swap_i test"
        end if
        
        ! Test single precision swap
        a_r = 1.5_SP
        b_r = 2.5_SP
        call swap(a_r, b_r)
        if (abs(a_r - 2.5_SP) > epsilon(a_r) * 10.0 .or. abs(b_r - 1.5_SP) > epsilon(b_r) * 10.0) then
            print *, "FAIL: swap_r test failed"
            print *, "Expected: a=2.5, b=1.5"
            print *, "Got: a=", a_r, ", b=", b_r
            status = status + 1
        else
            print *, "PASS: swap_r test"
        end if
        
        ! Test double precision swap
        a_d = 3.14159_DP
        b_d = 2.71828_DP
        call swap(a_d, b_d)
        if (abs(a_d - 2.71828_DP) > epsilon(a_d) * 10.0 .or. abs(b_d - 3.14159_DP) > epsilon(b_d) * 10.0) then
            print *, "FAIL: swap_d test failed"
            print *, "Expected: a=2.71828, b=3.14159"
            print *, "Got: a=", a_d, ", b=", b_d
            status = status + 1
        else
            print *, "PASS: swap_d test"
        end if
        
    end subroutine test_swap_functions
    
    subroutine test_assert_eq_functions(status)
        integer, intent(inout) :: status
        integer :: result
        
        print *, "Testing assert_eq functions..."
        
        ! Test assert_eq2
        result = assert_eq2(5, 5, "test_assert_eq2")
        if (result /= 5) then
            print *, "FAIL: assert_eq2 test failed for equal values"
            status = status + 1
        else
            print *, "PASS: assert_eq2 test for equal values"
        end if
        
        ! Test assert_eq3
        result = assert_eq3(10, 10, 10, "test_assert_eq3")
        if (result /= 10) then
            print *, "FAIL: assert_eq3 test failed for equal values"
            status = status + 1
        else
            print *, "PASS: assert_eq3 test for equal values"
        end if
        
        ! Test assert_eq4
        result = assert_eq4(7, 7, 7, 7, "test_assert_eq4")
        if (result /= 7) then
            print *, "FAIL: assert_eq4 test failed for equal values"
            status = status + 1
        else
            print *, "PASS: assert_eq4 test for equal values"
        end if
        
    end subroutine test_assert_eq_functions

end program test_nrutil
