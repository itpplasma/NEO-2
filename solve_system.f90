
!
! AUTHOR: Bernhard Seiwald
!
! DATE:   18.07.2001
!


! simple wrapper for solvers for real system of linear
! equations  A * X = B


MODULE solve_systems

  USE inter_precision, ONLY: I4B, DP

  IMPLICIT NONE

  PUBLIC :: solve_eqsys

! --------------------------------------------------------------------


CONTAINS

  SUBROUTINE solve_eqsys(a, b, info)

!    USE inter_interfaces, ONLY: ludcmp, lubksb 
    use mpiprovider_module
    use parallelStorage_module

    IMPLICIT NONE

    REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: a
    REAL(DP), DIMENSION(:),   INTENT(INOUT) :: b
    INTEGER(I4B),             INTENT(OUT)   :: info
    INTEGER(I4B) :: i_alloc
    INTEGER(I4B) :: n, nrhs, lda, ldb
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: ipiv
    double precision :: atime

! --------------------------------------------------------------------

    lda  = SIZE(a,1)
    n    = SIZE(a,2)
    ldb  = SIZE(b,1)
    nrhs = 1

    ALLOCATE(ipiv(n),  stat = i_alloc)
    IF(i_alloc /= 0) STOP 'solve_eqsys: Allocation for array failed!'
    atime = MPI_WTime()
    CALL dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
    parallel_storage%dgesvTime = parallel_storage%dgesvTime + MPI_WTime() - atime
    info = 0

    DEALLOCATE(ipiv,  stat = i_alloc)
    IF(i_alloc /= 0) STOP 'solve_eqsys: Deallocation for array failed!'

  END SUBROUTINE solve_eqsys

END MODULE solve_systems



