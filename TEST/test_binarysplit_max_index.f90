!> Calling eval_bslimit with an index beyond bsfunc_max_index must terminate
!> with the bsfunc_max_index error; CMake matches that message, so another
!> failure does not pass.
program test_binarysplit_max_index
    use nrtype, only: dp
    use binarysplit_int, only: create_bslimit, eval_bslimit
    implicit none

    logical :: ok

    call create_bslimit(max_index=10)
    ok = eval_bslimit(1.0_dp, 0, 10)
    if (.not. ok) error stop "unexpected: limit reached at maxind == max_index"
    ok = eval_bslimit(1.0_dp, 0, 11)
    print *, "FAIL: returned instead of stopping"
end program test_binarysplit_max_index
