program test_native_sparse_adapter
    use sparse_mod, only: sparse_solve, sparse_solve_suitesparse, &
        sparse_solve_method, sparse_talk
    use fortsparse, only: FORTSPARSE_SINGULAR
    implicit none
    integer, parameter :: dp = kind(1d0), n = 4
    real(dp) :: a(n,n), exact(n), b(n), work(n), block(n,3), expected(n,3)
    real(dp) :: small(3), small_values(4), saved(n)
    complex(dp) :: z(n,n), zx(n), zb(n), zblock(n,3), zexpected(n,3)
    complex(dp) :: real_factor_rhs(n), real_factor_answer(n)
    complex(dp) :: zvalues(32)
    integer :: method, j, i, nz, rows(32), cols(32)
    real(dp) :: values(32)
    integer, parameter :: small_rows(4) = [1,1,2,3]
    integer, parameter :: small_csc(4) = [1,2,4,5], small_coo(4) = [1,2,2,3]
    character(20) :: argument
    call get_command_argument(1,argument)
    sparse_talk = .false.
    if (trim(argument) == '--singular') then
        print '(a,i0)', 'Expected singular code=',FORTSPARSE_SINGULAR
        a = 0d0
        b = 1d0
        call sparse_solve(a,b)
        error stop 'singular matrix unexpectedly accepted'
    end if
    a = reshape([0d0,3d0,0d0,1d0, 2d0,4d0,1d0,0d0, &
                 0d0,1d0,5d0,2d0, 1d0,0d0,2d0,6d0],[n,n])
    exact = [0.5d0,-1d0,2d0,0.25d0]
    b = matmul(a,exact)
    z = cmplx(a,0d0,dp)
    z(1,2) = z(1,2)+(0d0,0.5d0)
    z(3,1) = (0d0,-0.25d0)
    zx = cmplx(exact,exact(n:1:-1),dp)
    zb = matmul(z,zx)
    do method = 2,3
        sparse_solve_method = method
        ! Dense scalar API, explicit factor, repeated solve and automatic free.
        work = b
        call sparse_solve(a,work,1)
        call close(work,b,'factor-only overwrote RHS')
        call sparse_solve(a,work,2)
        call close(work,exact,'dense real solve')
        call close(matmul(a,work),b,'original real matrix residual')
        work = 2*b
        call sparse_solve(a,work,2)
        call close(work,2*exact,'retained real factor')
        call sparse_solve(a,work,3)
        work = b
        call sparse_solve(a,work)
        call close(work,exact,'omitted operation resets to factor/solve/free')
        do j = 1,3
            expected(:,j) = j*exact
            block(:,j) = j*b
            zexpected(:,j) = j*zx
            zblock(:,j) = j*zb
        end do
        saved = block(:,1)
        call sparse_solve(a,block,1)
        call close(block(:,1),saved,'block factor-only RHS')
        call sparse_solve(a,block,2)
        do j = 1,3
            call close(block(:,j),expected(:,j),'retained many real RHS')
        end do
        call sparse_solve(a,block,3)
        work = b
        call sparse_solve(a,work,2) ! Historical solve-only auto-factors if needed.
        call close(work,exact,'solve-only auto-factor')
        ! Changing method on a retained factor preserves with/without refinement.
        sparse_solve_method = 5-method
        work = b
        call sparse_solve(a,work,2)
        call close(work,exact,'retained method switch')
        sparse_solve_method = method
        call sparse_solve(a,work,3)
        ! Arnoldi's real-factor branch splits a complex RHS into real/imaginary
        ! solves, retaining the same real factor for both components.
        real_factor_rhs = matmul(cmplx(a,0d0,dp),zx)
        work = real(real_factor_rhs,dp)
        call sparse_solve(a,work,1)
        call sparse_solve(a,work,2)
        real_factor_answer = cmplx(work,0d0,dp)
        work = aimag(real_factor_rhs)
        call sparse_solve(a,work,2)
        real_factor_answer = real_factor_answer+cmplx(0d0,work,dp)
        call zclose(real_factor_answer,zx,'Arnoldi real factor/complex RHS split')
        call sparse_solve(a,work,3)
        ! True complex scalar and many-RHS branches.
        real_factor_answer = zb
        call sparse_solve(z,real_factor_answer,1)
        call zclose(real_factor_answer,zb,'complex factor-only RHS')
        call sparse_solve(z,real_factor_answer,2)
        call zclose(real_factor_answer,zx,'complex scalar factor')
        call zclose(matmul(z,real_factor_answer),zb,'original complex matrix residual')
        call sparse_solve(z,real_factor_answer,3)
        real_factor_answer = zb
        call sparse_solve(z,real_factor_answer)
        call zclose(real_factor_answer,zx,'complex implicit lifecycle')
        call sparse_solve(z,zblock,1)
        call sparse_solve(z,zblock,2)
        do j = 1,3
            call zclose(zblock(:,j),zexpected(:,j),'complex block reuse')
        end do
        call sparse_solve(z,zblock,3)
        ! Unsorted duplicate COO exercises the library canonicalizer.
        nz = 0
        do j = n,1,-1
            do i = n,1,-1
                if (a(i,j) == 0d0) cycle
                rows(nz+1:nz+2) = i
                cols(nz+1:nz+2) = j
                values(nz+1:nz+2) = 0.5d0*a(i,j)
                nz = nz+2
            end do
        end do
        work = b
        call sparse_solve(n,n,nz,rows(:nz),cols(:nz),values(:nz),work)
        call close(work,exact,'unsorted duplicate COO')
        do j = 1,3
            block(:,j) = j*b
        end do
        call sparse_solve(n,n,nz,rows(:nz),cols(:nz),values(:nz),block,1)
        call sparse_solve(n,n,nz,rows(:nz),cols(:nz),values(:nz),block,2)
        do j = 1,3
            call close(block(:,j),expected(:,j),'raw real block RHS')
        end do
        call sparse_solve(n,n,nz,rows(:nz),cols(:nz),values(:nz),block,3)
        nz = 0
        do j = n,1,-1
            do i = n,1,-1
                if (z(i,j) == (0d0,0d0)) cycle
                nz = nz+1
                rows(nz) = i
                cols(nz) = j
                zvalues(nz) = z(i,j)
            end do
        end do
        real_factor_answer = zb
        call sparse_solve(n,n,nz,rows(:nz),cols(:nz),zvalues(:nz),real_factor_answer)
        call zclose(real_factor_answer,zx,'raw complex scalar RHS')
        do j = 1,3
            zblock(:,j) = j*zb
        end do
        call sparse_solve(n,n,nz,rows(:nz),cols(:nz),zvalues(:nz),zblock,1)
        call sparse_solve(n,n,nz,rows(:nz),cols(:nz),zvalues(:nz),zblock,2)
        do j = 1,3
            call zclose(zblock(:,j),zexpected(:,j),'raw complex block RHS')
        end do
        call sparse_solve(n,n,nz,rows(:nz),cols(:nz),zvalues(:nz),zblock,3)
        ! CSC pointer size equals nnz; distinguish it from full column indices.
        small_values = [2d0,1d0,3d0,4d0]
        small = [1d0,1d0,1d0]
        call sparse_solve(3,3,4,small_rows,small_csc,small_values,small)
        call close(small,[1d0/3d0,1d0/3d0,0.25d0],'ambiguous-size CSC')
        small = 1d0
        call sparse_solve(3,3,4,small_rows,small_coo,small_values,small)
        call close(small,[1d0/3d0,1d0/3d0,0.25d0],'same-size COO')
        small = 1d0
        call sparse_solve_suitesparse(3,3,4,small_rows,small_csc,small_values,small,0)
        call close(small,[1d0/3d0,1d0/3d0,0.25d0],'legacy low-level name routes native')
    end do
    print *, 'All tests passed! Native sparse lifecycle, real/complex and RHS contracts.'
contains
    subroutine close(actual,expected,label)
        real(dp), intent(in) :: actual(:),expected(:)
        character(*), intent(in) :: label
        if (maxval(abs(actual-expected)) > 5d-12*max(1d0,maxval(abs(expected)))) then
            print *, 'FAIL ',label,maxval(abs(actual-expected))
            stop 1
        end if
    end subroutine
    subroutine zclose(actual,expected,label)
        complex(dp), intent(in) :: actual(:),expected(:)
        character(*), intent(in) :: label
        if (maxval(abs(actual-expected)) > 5d-12*max(1d0,maxval(abs(expected)))) then
            print *, 'FAIL ',label,maxval(abs(actual-expected))
            stop 1
        end if
    end subroutine
end program
