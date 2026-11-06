! A singular join system must report the LAPACK status, not only ierr=2.
program test_join_ripples_info
    use propagator_mod, only: prop_c
    implicit none

    integer :: ierr

    allocate(prop_c)
    allocate(prop_c%prev)

    ! One passing class, no velocity moments: the join matrix is 1 - 1*1 = 0.
    prop_c%p%npass_l = 1
    prop_c%p%nvelocity = 0
    allocate(prop_c%p%cmat(1, 1), prop_c%p%amat_p_m(1, 1), prop_c%p%source_m(1, 3))
    prop_c%p%cmat = 1.0d0
    prop_c%p%amat_p_m = 1.0d0
    prop_c%p%source_m = 0.0d0
    prop_c%fieldpropagator_tag_s = 13
    prop_c%fieldpropagator_tag_e = 13

    prop_c%prev%p%npass_r = 1
    allocate(prop_c%prev%p%cmat(1, 1), prop_c%prev%p%amat_m_p(1, 1), &
        prop_c%prev%p%source_p(1, 3))
    prop_c%prev%p%cmat = 1.0d0
    prop_c%prev%p%amat_m_p = 1.0d0
    prop_c%prev%p%source_p = 0.0d0
    prop_c%prev%fieldpropagator_tag_s = 12
    prop_c%prev%fieldpropagator_tag_e = 12

    call join_ripples(ierr)

    if (ierr /= 2) then
        print *, 'FAIL: singular join returned ierr =', ierr
    end if
end program test_join_ripples_info
