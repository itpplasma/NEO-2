program test_solve_report
   ! The axisymmetric system is assembled as a bulk block and a boundary block
   ! that are concatenated before the solve. A residual that walks only the
   ! bulk block measures a system that was never solved: it stays at the size
   ! of the boundary contribution however tight the solver tolerance is. These
   ! cases pin that difference down.
   use nrtype, only: dp
   use solve_report_mod, only: csr_add
   implicit none

   integer :: status

   status = 0
   call test_both_blocks_reproduce_the_solved_matrix(status)
   call test_omitted_boundary_block_floors_the_residual(status)
   call test_complex_block_matches_real_block(status)

   if (status == 0) then
      print *, "All tests passed!"
   else
      print *, "Some tests failed. Status:", status
      error stop
   end if

contains

   ! A 3x3 matrix given as two blocks in compressed-row storage. The bulk block
   ! carries rows 1 and 2, the boundary block row 3, which couples both bulk
   ! columns. f is chosen so that q = A f is exact.
   subroutine setup(f, q)
      real(dp), intent(out) :: f(3), q(3)

      f = [1.5_dp, -0.75_dp, 0.0_dp]
      q = [2.0_dp*f(1), 4.0_dp*f(2), 0.5_dp*f(1) - 0.25_dp*f(2)]
   end subroutine setup

   subroutine residual(include_boundary, f, q, value)
      logical, intent(in) :: include_boundary
      real(dp), intent(in) :: f(:), q(:)
      real(dp), intent(out) :: value
      complex(dp), allocatable :: lhs(:)
      complex(dp) :: x(size(f))
      integer :: bulk_rows(2), bulk_cols(2), bound_rows(2), bound_cols(2)
      real(dp) :: bulk_amat(2), bound_amat(2)

      bulk_rows = [1, 2]
      bulk_cols = [1, 2]
      bulk_amat = [2.0_dp, 4.0_dp]
      bound_rows = [3, 3]
      bound_cols = [1, 2]
      bound_amat = [0.5_dp, -0.25_dp]

      x = cmplx(f, 0.0_dp, dp)
      allocate(lhs(size(f)))
      lhs = 0.0d0
      call csr_add(2, bulk_rows, bulk_cols, bulk_amat, x, lhs)
      if (include_boundary) then
         call csr_add(2, bound_rows, bound_cols, bound_amat, x, lhs)
      end if
      value = sum(abs(cmplx(q, 0.0_dp, dp) - lhs)) &
         & / max(sum(abs(q)) + sum(abs(lhs)), tiny(1.0_dp))
   end subroutine residual

   subroutine test_both_blocks_reproduce_the_solved_matrix(status)
      integer, intent(inout) :: status
      real(dp) :: f(3), q(3), value

      call setup(f, q)
      call residual(.true., f, q, value)
      if (value > 1000.0_dp*tiny(1.0_dp)) then
         print *, "FAIL: complete matrix residual", value
         status = status + 1
      end if
   end subroutine test_both_blocks_reproduce_the_solved_matrix

   subroutine test_omitted_boundary_block_floors_the_residual(status)
      integer, intent(inout) :: status
      real(dp) :: f(3), q(3), value

      call setup(f, q)
      ! Dropping the boundary block leaves 0.5*1.5 - 0.25*(-0.75) = 0.9375 of
      ! a total scale of 12.9375 unaccounted, so the residual sits near 0.072
      ! and cannot fall further by tightening any tolerance.
      call residual(.false., f, q, value)
      if (value < 1.0d-2) then
         print *, "FAIL: omitted boundary block should floor the residual, got", value
         status = status + 1
      end if
   end subroutine test_omitted_boundary_block_floors_the_residual

   subroutine test_complex_block_matches_real_block(status)
      integer, intent(inout) :: status
      integer :: bulk(2)
      real(dp) :: f(3), bulk_amat(2)
      complex(dp) :: camat(2), x(3), lhs_real(3), lhs_complex(3)

      bulk = [1, 2]
      bulk_amat = [2.0_dp, 4.0_dp]
      camat = cmplx(bulk_amat, 0.0_dp, dp)
      f = [1.5_dp, -0.75_dp, 0.0_dp]
      x = cmplx(f, 0.0_dp, dp)
      lhs_real = 0.0d0
      lhs_complex = 0.0d0
      call csr_add(2, bulk, bulk, bulk_amat, x, lhs_real)
      call csr_add(2, bulk, bulk, camat, x, lhs_complex)
      if (maxval(abs(lhs_real - lhs_complex)) > 0.0_dp) then
         print *, "FAIL: real and complex blocks disagree"
         status = status + 1
      end if
      if (abs(lhs_real(1) - 3.0_dp) > 100.0_dp*tiny(1.0_dp)) then
         print *, "FAIL: bulk product entry 1", lhs_real(1)
         status = status + 1
      end if
   end subroutine test_complex_block_matches_real_block

end program test_solve_report
