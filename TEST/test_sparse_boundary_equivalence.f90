program test_sparse_boundary_equivalence
    ! Boundary before/after safety check. Solves representative nonsymmetric
    ! real and complex sparse systems through sparse_mod's public sparse_solve
    ! and independently via a dense LAPACK reference (dgesv/zgesv on the
    ! densified matrix), then asserts the two solutions agree to < 1e-10.
    use sparse_mod, only: sparse_solve, sparse_solve_method
    implicit none

    integer, parameter :: dp = kind(1.0d0)
    integer :: status

    status = 0
    sparse_solve_method = 3

    call real_case(1, status)
    call real_case(2, status)
    call complex_case(1, status)
    call complex_case(2, status)

    if (status == 0) then
        print *, "All tests passed!"
    else
        print *, "Some tests failed. Status:", status
        error stop
    end if

contains

    ! Build a representative nonsymmetric real CSC system and its dense form.
    subroutine make_real(which, nrow, ncol, nz, irow, pcol, val, dense, b)
        integer, intent(in) :: which
        integer, intent(out) :: nrow, ncol, nz
        integer, allocatable, intent(out) :: irow(:), pcol(:)
        real(dp), allocatable, intent(out) :: val(:), b(:)
        real(dp), allocatable, intent(out) :: dense(:,:)

        integer :: i, j, k
        nrow = 5; ncol = 5
        allocate (dense(5,5), b(5))
        dense = 0.0_dp
        if (which == 1) then
            dense(1,:) = [ 4.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 2.0_dp]
            dense(2,:) = [ 0.0_dp, 3.0_dp, 1.0_dp, 0.0_dp, 0.0_dp]
            dense(3,:) = [ 1.0_dp, 0.0_dp, 5.0_dp, 2.0_dp, 0.0_dp]
            dense(4,:) = [ 0.0_dp, 0.0_dp, 0.0_dp, 6.0_dp, 1.0_dp]
            dense(5,:) = [ 3.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 7.0_dp]
        else
            dense(1,:) = [ 2.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, 0.0_dp]
            dense(2,:) = [-1.0_dp, 4.0_dp, 0.0_dp, 1.0_dp, 0.0_dp]
            dense(3,:) = [ 0.0_dp, -2.0_dp, 3.0_dp, 0.0_dp, 1.0_dp]
            dense(4,:) = [ 0.0_dp, 0.0_dp, -1.0_dp, 5.0_dp, 0.0_dp]
            dense(5,:) = [ 1.0_dp, 0.0_dp, 0.0_dp, -3.0_dp, 4.0_dp]
        end if
        do i = 1, 5
            b(i) = real(i, dp)
        end do
        call densify_real(dense, nrow, ncol, nz, irow, pcol, val)
    end subroutine make_real

    ! Convert a dense matrix to 1-based CSC (column-major).
    subroutine densify_real(dense, nrow, ncol, nz, irow, pcol, val)
        real(dp), intent(in) :: dense(:,:)
        integer, intent(in) :: nrow, ncol
        integer, intent(out) :: nz
        integer, allocatable, intent(out) :: irow(:), pcol(:)
        real(dp), allocatable, intent(out) :: val(:)
        integer :: i, j, k

        nz = count(dense /= 0.0_dp)
        allocate (irow(nz), pcol(ncol+1), val(nz))
        k = 0
        pcol(1) = 1
        do j = 1, ncol
            do i = 1, nrow
                if (dense(i,j) /= 0.0_dp) then
                    k = k + 1
                    irow(k) = i
                    val(k) = dense(i,j)
                end if
            end do
            pcol(j+1) = k + 1
        end do
    end subroutine densify_real

    subroutine real_case(which, status)
        integer, intent(in) :: which
        integer, intent(inout) :: status
        integer :: nrow, ncol, nz, info
        integer, allocatable :: irow(:), pcol(:), ipiv(:)
        real(dp), allocatable :: val(:), b(:), dense(:,:)
        real(dp), allocatable :: bsparse(:), bdense(:)
        real(dp) :: r

        call make_real(which, nrow, ncol, nz, irow, pcol, val, dense, b)
        allocate (bsparse(nrow), bdense(nrow), ipiv(nrow))
        bsparse = b
        bdense = b
        call sparse_solve(nrow, ncol, nz, irow, pcol, val, bsparse, 0)
        call dgesv(nrow, 1, dense, nrow, ipiv, bdense, nrow, info)
        r = maxval(abs(bsparse - bdense))
        if (info == 0 .and. r < 1.0e-10_dp) then
            print *, "PASS: real boundary case", which
        else
            print *, "FAIL: real boundary case", which, " diff", r, " info", info
            status = status + 1
        end if
    end subroutine real_case

    ! Build a representative nonsymmetric complex CSC system and dense form.
    subroutine make_complex(which, nrow, ncol, nz, irow, pcol, val, dense, b)
        integer, intent(in) :: which
        integer, intent(out) :: nrow, ncol, nz
        integer, allocatable, intent(out) :: irow(:), pcol(:)
        complex(dp), allocatable, intent(out) :: val(:), b(:)
        complex(dp), allocatable, intent(out) :: dense(:,:)

        integer :: i
        nrow = 4; ncol = 4
        allocate (dense(4,4), b(4))
        dense = (0.0_dp, 0.0_dp)
        if (which == 1) then
            dense(1,1) = (3.0_dp, 1.0_dp); dense(1,3) = (1.0_dp, 0.0_dp)
            dense(2,2) = (4.0_dp, -1.0_dp); dense(2,1) = (0.0_dp, 1.0_dp)
            dense(3,3) = (5.0_dp, 2.0_dp); dense(3,4) = (1.0_dp, 1.0_dp)
            dense(4,4) = (2.0_dp, -1.0_dp); dense(4,2) = (1.0_dp, 0.0_dp)
        else
            dense(1,1) = (2.0_dp, 0.0_dp); dense(1,2) = (-1.0_dp, 1.0_dp)
            dense(2,2) = (3.0_dp, 1.0_dp); dense(2,4) = (0.0_dp, 2.0_dp)
            dense(3,1) = (1.0_dp, -1.0_dp); dense(3,3) = (4.0_dp, 0.0_dp)
            dense(4,3) = (-2.0_dp, 0.0_dp); dense(4,4) = (5.0_dp, 1.0_dp)
        end if
        do i = 1, 4
            b(i) = cmplx(real(i, dp), 1.0_dp, dp)
        end do
        call densify_complex(dense, nrow, ncol, nz, irow, pcol, val)
    end subroutine make_complex

    subroutine densify_complex(dense, nrow, ncol, nz, irow, pcol, val)
        complex(dp), intent(in) :: dense(:,:)
        integer, intent(in) :: nrow, ncol
        integer, intent(out) :: nz
        integer, allocatable, intent(out) :: irow(:), pcol(:)
        complex(dp), allocatable, intent(out) :: val(:)
        integer :: i, j, k

        nz = count(dense /= (0.0_dp, 0.0_dp))
        allocate (irow(nz), pcol(ncol+1), val(nz))
        k = 0
        pcol(1) = 1
        do j = 1, ncol
            do i = 1, nrow
                if (dense(i,j) /= (0.0_dp, 0.0_dp)) then
                    k = k + 1
                    irow(k) = i
                    val(k) = dense(i,j)
                end if
            end do
            pcol(j+1) = k + 1
        end do
    end subroutine densify_complex

    subroutine complex_case(which, status)
        integer, intent(in) :: which
        integer, intent(inout) :: status
        integer :: nrow, ncol, nz, info
        integer, allocatable :: irow(:), pcol(:), ipiv(:)
        complex(dp), allocatable :: val(:), b(:), dense(:,:)
        complex(dp), allocatable :: bsparse(:), bdense(:)
        real(dp) :: r

        call make_complex(which, nrow, ncol, nz, irow, pcol, val, dense, b)
        allocate (bsparse(nrow), bdense(nrow), ipiv(nrow))
        bsparse = b
        bdense = b
        call sparse_solve(nrow, ncol, nz, irow, pcol, val, bsparse, 0)
        call zgesv(nrow, 1, dense, nrow, ipiv, bdense, nrow, info)
        r = maxval(abs(bsparse - bdense))
        if (info == 0 .and. r < 1.0e-10_dp) then
            print *, "PASS: complex boundary case", which
        else
            print *, "FAIL: complex boundary case", which, " diff", r, " info", info
            status = status + 1
        end if
    end subroutine complex_case

end program test_sparse_boundary_equivalence
