program test_sparse_memory_config
    use, intrinsic :: iso_fortran_env, only: int64
    use sparse_mod, only: sparse_host_memory_limit_bytes, sparse_low_memory, &
                          sparse_validate_memory_options, sparse_solve
    implicit none
    character(32) :: argument
    character(128) :: record
    integer :: ios
    real(kind(1d0)) :: matrix(2, 2), rhs(2)
    namelist /settings/ sparse_host_memory_limit_bytes, sparse_low_memory
    call get_command_argument(1, argument)
    matrix = reshape([3d0, 1d0, 2d0, 4d0], [2, 2])
    rhs = [7d0, 9d0]
    select case (trim(argument))
    case ('--negative')
        sparse_host_memory_limit_bytes = -1_int64
    case ('--legacy-cap')
        sparse_host_memory_limit_bytes = 4096_int64
    case ('--legacy-low')
        sparse_low_memory = .true.
    case default
        if (sparse_host_memory_limit_bytes /= 0_int64 .or. sparse_low_memory) &
            error stop 'memory defaults changed'
        record = '&settings sparse_host_memory_limit_bytes=4294967296, '// &
                 'sparse_low_memory=.true. /'
        read (record, nml=settings, iostat=ios)
        if (ios /= 0) error stop 'memory namelist parse failed'
        if (sparse_host_memory_limit_bytes /= 4294967296_int64) &
            error stop '64-bit byte cap truncated'
        if (.not. sparse_low_memory) error stop 'LOW_MEMORY setting lost'
        sparse_host_memory_limit_bytes = 0_int64
        sparse_low_memory = .false.
    end select
    if (len_trim(argument) > 0) then
        ! Verify actual public solve entry rejects controls, rather than only
        ! testing the input-validation routine used by QL and PAR startup.
        call sparse_solve(matrix, rhs)
        error stop 'invalid configuration accepted by solve'
    end if
    call sparse_validate_memory_options()
    call sparse_solve(matrix, rhs)
    if (maxval(abs(rhs - [1d0, 2d0])) > 1d-13) error stop 'default solve changed'
    print '(a)', 'All tests passed!'
end program test_sparse_memory_config
