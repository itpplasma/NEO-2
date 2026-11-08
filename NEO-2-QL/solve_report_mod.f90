module solve_report_mod
   !> Normwise backward error of a sparse linear system as solved.
   !>
   !> The axisymmetric and nonaxisymmetric systems are assembled from several
   !> blocks that are concatenated into the matrix handed to the solver, so a
   !> residual that walks only part of the stored blocks reports the error of
   !> a different system than the one that was solved and levels off near the
   !> size of the omitted block instead of near the solver tolerance. This
   !> module takes the blocks that form the solved matrix and no others.
   use nrtype, only: dp
   implicit none

   private
   public :: csr_add

   interface csr_add
      module procedure csr_add_real
      module procedure csr_add_complex
   end interface csr_add

contains

   !> Add A*x to lhs for one compressed-row-storage block. The axisymmetric
   !> matrix is held as real, the nonaxisymmetric one as complex; the solution
   !> and the accumulator are complex in both cases.
   subroutine csr_add_real(count, irow, icol, amat, x, lhs)
      integer, intent(in) :: count
      integer, intent(in) :: irow(:), icol(:)
      real(dp), intent(in) :: amat(:)
      complex(dp), intent(in) :: x(:)
      complex(dp), intent(inout) :: lhs(:)
      integer :: entry

      do entry = 1, count
         lhs(irow(entry)) = lhs(irow(entry)) + amat(entry)*x(icol(entry))
      end do
   end subroutine csr_add_real

   subroutine csr_add_complex(count, irow, icol, amat, x, lhs)
      integer, intent(in) :: count
      integer, intent(in) :: irow(:), icol(:)
      complex(dp), intent(in) :: amat(:), x(:)
      complex(dp), intent(inout) :: lhs(:)
      integer :: entry

      do entry = 1, count
         lhs(irow(entry)) = lhs(irow(entry)) + amat(entry)*x(icol(entry))
      end do
   end subroutine csr_add_complex

end module solve_report_mod
