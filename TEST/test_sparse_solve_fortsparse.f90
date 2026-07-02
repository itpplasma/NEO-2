program test_sparse_solve_fortsparse
    ! Exercises every path of sparse_solve_fortsparse and the rewritten
    ! sparse_mod leaves: real and complex, iopt = 0 and the 1->2->3 reuse
    ! sequence, refine on (sparse_solve_method 2) and off (3), and a 2-D RHS
    ! through the b2_loop path. Each solution is checked against a known
    ! answer and a residual bound.
    use sparse_mod, only: sparse_solve, sparse_solve_method
    implicit none

    integer, parameter :: dp = kind(1.0d0)
    integer :: status

    status = 0

    call test_real_iopt0(3, status)
    call test_real_iopt0(2, status)
    call test_real_reuse(3, status)
    call test_real_reuse(2, status)
    call test_complex_iopt0(3, status)
    call test_complex_iopt0(2, status)
    call test_complex_reuse(3, status)
    call test_real_b2(3, status)
    call test_real_b2(2, status)
    call test_complex_b2(3, status)
    call test_poisson(3, status)
    call test_poisson(2, status)

    if (status == 0) then
        print *, "All tests passed!"
    else
        print *, "Some tests failed. Status:", status
        error stop
    end if

contains

    ! Small 3x3 nonsymmetric system with a known solution.
    subroutine small_real(nrow, ncol, nz, irow, pcol, val, b, xref)
        integer, intent(out) :: nrow, ncol, nz
        integer, allocatable, intent(out) :: irow(:), pcol(:)
        real(dp), allocatable, intent(out) :: val(:), b(:), xref(:)

        ! A = [ 2 0 1 ; 0 3 0 ; 1 0 4 ], x = [1;2;3]
        nrow = 3; ncol = 3; nz = 5
        allocate (pcol(4), irow(5), val(5), b(3), xref(3))
        pcol = [1, 3, 4, 6]
        irow = [1, 3, 2, 1, 3]
        val  = [2.0_dp, 1.0_dp, 3.0_dp, 1.0_dp, 4.0_dp]
        xref = [1.0_dp, 2.0_dp, 3.0_dp]
        b = [5.0_dp, 6.0_dp, 13.0_dp]
    end subroutine small_real

    ! Small 3x3 complex nonsymmetric system with a known solution.
    subroutine small_complex(nrow, ncol, nz, irow, pcol, val, b, xref)
        integer, intent(out) :: nrow, ncol, nz
        integer, allocatable, intent(out) :: irow(:), pcol(:)
        complex(dp), allocatable, intent(out) :: val(:), b(:), xref(:)

        complex(dp) :: a11, a31, a22, a13, a33
        nrow = 3; ncol = 3; nz = 5
        allocate (pcol(4), irow(5), val(5), b(3), xref(3))
        pcol = [1, 3, 4, 6]
        irow = [1, 3, 2, 1, 3]
        a11 = (2.0_dp, 1.0_dp); a31 = (1.0_dp, 0.0_dp); a22 = (3.0_dp, -1.0_dp)
        a13 = (0.0_dp, 1.0_dp); a33 = (4.0_dp, 2.0_dp)
        val = [a11, a31, a22, a13, a33]
        xref = [(1.0_dp, 1.0_dp), (0.0_dp, 2.0_dp), (-1.0_dp, 1.0_dp)]
        ! b = A * xref, columnwise
        b(1) = a11*xref(1) + a13*xref(3)
        b(2) = a22*xref(2)
        b(3) = a31*xref(1) + a33*xref(3)
    end subroutine small_complex

    subroutine test_real_iopt0(method, status)
        integer, intent(in) :: method
        integer, intent(inout) :: status
        integer :: nrow, ncol, nz
        integer, allocatable :: irow(:), pcol(:)
        real(dp), allocatable :: val(:), b(:), xref(:)

        sparse_solve_method = method
        call small_real(nrow, ncol, nz, irow, pcol, val, b, xref)
        call sparse_solve(nrow, ncol, nz, irow, pcol, val, b, 0)
        call check_real(b, xref, 'real iopt0', method, status)
    end subroutine test_real_iopt0

    subroutine test_real_reuse(method, status)
        integer, intent(in) :: method
        integer, intent(inout) :: status
        integer :: nrow, ncol, nz
        integer, allocatable :: irow(:), pcol(:)
        real(dp), allocatable :: val(:), b(:), xref(:), b2(:)

        sparse_solve_method = method
        call small_real(nrow, ncol, nz, irow, pcol, val, b, xref)
        allocate (b2(3))
        b2 = [7.0_dp, 9.0_dp, 14.0_dp]   ! A * [2;3;3]

        call sparse_solve(nrow, ncol, nz, irow, pcol, val, b, 1)  ! factor
        call sparse_solve(nrow, ncol, nz, irow, pcol, val, b, 2)  ! solve 1
        call check_real(b, xref, 'real reuse rhs1', method, status)
        call sparse_solve(nrow, ncol, nz, irow, pcol, val, b2, 2) ! solve 2
        call check_real(b2, [2.0_dp, 3.0_dp, 3.0_dp], 'real reuse rhs2', &
            method, status)
        call sparse_solve(nrow, ncol, nz, irow, pcol, val, b, 3)  ! free
    end subroutine test_real_reuse

    subroutine test_complex_iopt0(method, status)
        integer, intent(in) :: method
        integer, intent(inout) :: status
        integer :: nrow, ncol, nz
        integer, allocatable :: irow(:), pcol(:)
        complex(dp), allocatable :: val(:), b(:), xref(:)

        sparse_solve_method = method
        call small_complex(nrow, ncol, nz, irow, pcol, val, b, xref)
        call sparse_solve(nrow, ncol, nz, irow, pcol, val, b, 0)
        call check_complex(b, xref, 'complex iopt0', method, status)
    end subroutine test_complex_iopt0

    subroutine test_complex_reuse(method, status)
        integer, intent(in) :: method
        integer, intent(inout) :: status
        integer :: nrow, ncol, nz
        integer, allocatable :: irow(:), pcol(:)
        complex(dp), allocatable :: val(:), b(:), xref(:)

        sparse_solve_method = method
        call small_complex(nrow, ncol, nz, irow, pcol, val, b, xref)
        call sparse_solve(nrow, ncol, nz, irow, pcol, val, b, 1)
        call sparse_solve(nrow, ncol, nz, irow, pcol, val, b, 2)
        call check_complex(b, xref, 'complex reuse', method, status)
        call sparse_solve(nrow, ncol, nz, irow, pcol, val, b, 3)
    end subroutine test_complex_reuse

    subroutine test_real_b2(method, status)
        integer, intent(in) :: method
        integer, intent(inout) :: status
        integer :: nrow, ncol, nz
        integer, allocatable :: irow(:), pcol(:)
        real(dp), allocatable :: val(:), b(:), xref(:)
        real(dp), allocatable :: bb(:,:), xx(:,:)

        sparse_solve_method = method
        call small_real(nrow, ncol, nz, irow, pcol, val, b, xref)
        allocate (bb(3, 2), xx(3, 2))
        bb(:,1) = b
        bb(:,2) = [7.0_dp, 9.0_dp, 14.0_dp]
        xx(:,1) = xref
        xx(:,2) = [2.0_dp, 3.0_dp, 3.0_dp]
        call sparse_solve(nrow, ncol, nz, irow, pcol, val, bb, 0)
        call check_real(bb(:,1), xx(:,1), 'real b2 col1', method, status)
        call check_real(bb(:,2), xx(:,2), 'real b2 col2', method, status)
    end subroutine test_real_b2

    subroutine test_complex_b2(method, status)
        integer, intent(in) :: method
        integer, intent(inout) :: status
        integer :: nrow, ncol, nz
        integer, allocatable :: irow(:), pcol(:)
        complex(dp), allocatable :: val(:), b(:), xref(:)
        complex(dp), allocatable :: bb(:,:)

        sparse_solve_method = method
        call small_complex(nrow, ncol, nz, irow, pcol, val, b, xref)
        allocate (bb(3, 2))
        bb(:,1) = b
        bb(:,2) = b
        call sparse_solve(nrow, ncol, nz, irow, pcol, val, bb, 0)
        call check_complex(bb(:,1), xref, 'complex b2 col1', method, status)
        call check_complex(bb(:,2), xref, 'complex b2 col2', method, status)
    end subroutine test_complex_b2

    ! Tridiagonal Poisson(n=64): A = tridiag(-1, 2, -1), known solution x.
    subroutine test_poisson(method, status)
        integer, intent(in) :: method
        integer, intent(inout) :: status
        integer, parameter :: n = 64
        integer :: nz, j, k
        integer, allocatable :: irow(:), pcol(:)
        real(dp), allocatable :: val(:), b(:), xref(:)

        sparse_solve_method = method
        nz = 3*n - 2
        allocate (pcol(n+1), irow(nz), val(nz), b(n), xref(n))
        do j = 1, n
            xref(j) = real(j, dp)
        end do
        k = 0
        pcol(1) = 1
        do j = 1, n
            if (j > 1) then
                k = k + 1; irow(k) = j - 1; val(k) = -1.0_dp
            end if
            k = k + 1; irow(k) = j; val(k) = 2.0_dp
            if (j < n) then
                k = k + 1; irow(k) = j + 1; val(k) = -1.0_dp
            end if
            pcol(j+1) = k + 1
        end do
        ! b = A * xref
        do j = 1, n
            b(j) = 2.0_dp*xref(j)
            if (j > 1) b(j) = b(j) - xref(j-1)
            if (j < n) b(j) = b(j) - xref(j+1)
        end do
        call sparse_solve(n, n, nz, irow, pcol, val, b, 0)
        call check_real(b, xref, 'poisson', method, status)
    end subroutine test_poisson

    subroutine check_real(x, xref, name, method, status)
        real(dp), intent(in) :: x(:), xref(:)
        character(*), intent(in) :: name
        integer, intent(in) :: method
        integer, intent(inout) :: status
        real(dp) :: r

        r = maxval(abs(x - xref))
        if (r < 1.0e-9_dp) then
            print *, "PASS: ", name, " method", method
        else
            print *, "FAIL: ", name, " method", method, " residual", r
            status = status + 1
        end if
    end subroutine check_real

    subroutine check_complex(x, xref, name, method, status)
        complex(dp), intent(in) :: x(:), xref(:)
        character(*), intent(in) :: name
        integer, intent(in) :: method
        integer, intent(inout) :: status
        real(dp) :: r

        r = maxval(abs(x - xref))
        if (r < 1.0e-9_dp) then
            print *, "PASS: ", name, " method", method
        else
            print *, "FAIL: ", name, " method", method, " residual", r
            status = status + 1
        end if
    end subroutine check_complex

end program test_sparse_solve_fortsparse
