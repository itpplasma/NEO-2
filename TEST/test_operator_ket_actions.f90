program test_operator_ket_actions
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use sparse_mod, only: sparse_coo_action
    implicit none
    complex(dp) :: z(3), m(3), c(3), r(3), k(3), e0(3), e1(3)
    complex(dp) :: actions(3,3), gram(3,3)
    real(dp) :: entry_sum
    integer :: source
    ! Analytic exchange Laplacian: common temperature is null, each species is not.
    do source = 0, 3
        z = (0.0_dp,0.0_dp)
        if (source == 0) z = (1.0_dp,0.0_dp)
        if (source > 0) z(source) = (1.0_dp,0.0_dp)
        call sparse_coo_action([1,2,3],[1,2,3],[2.0_dp,2.0_dp,2.0_dp],z,m)
        call sparse_coo_action([1,1,2,2,3,3],[2,3,1,3,1,2], &
            [1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp],z,c)
        if (source == 0) then
            if (maxval(abs(m-c)) > 1e-14_dp) error stop 'common null action'
        else
            if (abs(sum(abs(m-c))-4.0_dp) > 1e-14_dp) error stop 'non-null species'
            actions(:,source) = m-c
        end if
    end do
    gram = matmul(transpose(conjg(actions)),actions)
    if (maxval(abs(matmul(gram,[(1.0_dp,0.0_dp),(1.0_dp,0.0_dp), &
        (1.0_dp,0.0_dp)]))) > 1e-14_dp) error stop 'full action Gram null'
    if (abs(gram(1,1)-6.0_dp) > 1e-14_dp) error stop 'full action Gram scale'
    z = [(1.0_dp,1.0_dp),(2.0_dp,-1.0_dp),(0.0_dp,0.0_dp)]
    call sparse_coo_action([1,2],[2,1],[3.0_dp,-2.0_dp],z,m)
    if (maxval(abs(m-[(6.0_dp,-3.0_dp),(-2.0_dp,-2.0_dp), &
        (0.0_dp,0.0_dp)])) > 1e-14_dp) error stop 'coordinate orientation'
    ! Analytic speed nulls: R fixes m0, but leaves m1 null in K=M+R.
    e0 = [(1.0_dp,0.0_dp),(0.0_dp,0.0_dp),(0.0_dp,0.0_dp)]
    e1 = cshift(e0,-1)
    do source = 0, 1
        z = e0
        if (source == 1) z = e1
        call sparse_coo_action([3],[3],[1.0_dp],z,m)
        call sparse_coo_action([1],[1],[1.0_dp],z,r)
        call sparse_coo_action([3,1],[3,1],[1.0_dp,1.0_dp],z,k)
        if (maxval(abs(k-m-r)) > 1e-14_dp) error stop 'K/M/R identity'
        if (source == 0 .and. abs(k(1)-1.0_dp) > 1e-14_dp) error stop 'fixed m0'
        if (source == 1 .and. maxval(abs(k)) > 1e-14_dp) error stop 'remaining m1 null'
    end do
    z = (1.0_dp,0.0_dp)
    call sparse_coo_action([1,1],[1,2],[1.0_dp,-1.0_dp],z,m,entry_sum)
    if (maxval(abs(m)) /= 0.0_dp .or. entry_sum /= 2.0_dp) &
        error stop 'cancelled action needs nonzero entry scale'
    z = (0.0_dp,0.0_dp)
    call sparse_coo_action([1,1],[1,2],[1.0_dp,-1.0_dp],z,m,entry_sum)
    if (maxval(abs(m)) /= 0.0_dp .or. entry_sum /= 0.0_dp) &
        error stop 'zero entry scale'
    print *, 'All tests passed!'
end program test_operator_ket_actions
