MODULE sparse_mod

  USE sparse_types_mod, ONLY: dp, long
  USE sparse_conversion_mod
  USE sparse_io_mod
  USE sparse_arithmetic_mod
  USE sparse_solvers_mod
  IMPLICIT NONE

  ! Re-export sparse_solve_method and solver constants for backward compatibility
  PUBLIC :: sparse_solve_method
  PUBLIC :: SOLVER_UMFPACK, SOLVER_BICGSTAB
  
  ! Re-export conversion routines for backward compatibility
  PUBLIC :: column_pointer2full, column_full2pointer
  PUBLIC :: sparse2full, full2sparse
  
  ! Re-export I/O routines for backward compatibility
  PUBLIC :: load_mini_example, load_compressed_example
  PUBLIC :: load_standard_example, load_octave_matrices
  PUBLIC :: find_unit
  
  ! Re-export arithmetic routines for backward compatibility
  PUBLIC :: sparse_matmul, sparse_solver_test, sparse_talk
  
  ! Re-export solver routines for backward compatibility
  PUBLIC :: sparse_solve, sparse_solve_suitesparse
  PUBLIC :: factorization_exists

  PUBLIC sparse_example

CONTAINS

  !-------------------------------------------------------------------------------
  ! this will give an example of how to use the interface
  SUBROUTINE sparse_example

    INTEGER, DIMENSION(:), ALLOCATABLE :: pcol, icol
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: val

    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: A
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: b, x
    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: bb, xx

    COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE :: z_val
    COMPLEX(kind=dp), DIMENSION(:,:), ALLOCATABLE :: z_A
    COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE :: z_b, z_x
    COMPLEX(kind=dp), DIMENSION(:,:), ALLOCATABLE :: z_bb, z_xx

    INTEGER :: nrow, ncol, nz

    PRINT *, 'starting sparse_example'

    PRINT *, 'testing a real valued 5x5 mini-matrix'
    !testing a real 5x5 matrix
    CALL load_mini_example(A)
    ALLOCATE(b(SIZE(A,1)),x(SIZE(A,2)))
    b = 1.0_dp
    CALL sparse_matmul(A,b,x)
    PRINT *, 'x=',x

    x = 0.0_dp
    CALL sparse_solve(A,b,x)
    PRINT *, 'x=',x
    b = 1.0_dp
    CALL sparse_solver_test(A,b,x)
    PRINT *, 'b_out =',b

    ALLOCATE(bb(SIZE(A,1),7),xx(SIZE(A,2),7))
    bb(:,1) = 1
    bb(:,2) = 2
    bb(:,3) = 3
    bb(:,4) = 4
    bb(:,5) = 5
    bb(:,6) = 6
    bb(:,7) = 7
    CALL sparse_matmul(A,bb,xx)
    PRINT *, 'xx(:,1)=',xx(:,1)
    PRINT *, 'xx(:,7)=',xx(:,7)

    xx = 0.0_dp
    CALL sparse_solve(A,bb,xx)
    PRINT *, 'xx(:,1)=',xx(:,1)
    PRINT *, 'xx(:,7)=',xx(:,7)
    bb(:,1) = 1
    bb(:,2) = 2
    bb(:,3) = 3
    bb(:,4) = 4
    bb(:,5) = 5
    bb(:,6) = 6
    bb(:,7) = 7
    CALL sparse_solver_test(A,bb,xx)
    PRINT *, 'bb(:,1) =',bb(:,1)
    PRINT *, 'bb(:,7) =',bb(:,7)

    ! clear previous allocation
    DEALLOCATE(A,b,x,bb,xx)

    !-----------------------------------------------------------------------------------------

    PRINT *, 'testing 10x10 real octave matrix'

    CALL load_octave_matrices('../COMMON/testcases/octave_test_sparse.real', &
         nrow,ncol,nz,icol,pcol,val, &
         b,x,bb,xx)
    CALL sparse_matmul(icol,pcol,val,b,x)
    PRINT *, 'x=',x

    x = 0.0_dp
    CALL sparse_solve(icol,pcol,val,b,x)
    PRINT *, 'x=',x
    b(:) = (/ 1,2,3,4,5,6,7,8,9,10 /)
    CALL sparse_solver_test(icol,pcol,val,b,x)
    PRINT *, 'b_out =',b

    CALL sparse_matmul(icol,pcol,val,bb,xx)
    PRINT *, 'xx(:,1)=',xx(:,1)
    PRINT *, 'xx(:,10)=',xx(:,10)

    xx = 0.0_dp
    CALL sparse_solve(icol,pcol,val,bb,xx)
    PRINT *, 'xx(:,1)=',xx(:,1)
    PRINT *, 'xx(:,10)=',xx(:,10)

    bb(:,1) = (/ 1,2,3,4,5,6,7,8,9,10 /)
    bb(:,2) = (/ 2,3,4,5,6,7,8,9,10,11 /)
    bb(:,3) = (/ 3,4,5,6,7,8,9,10,11,12 /)
    bb(:,4) = (/ 4,5,6,7,8,9,10,11,12,13 /)
    bb(:,5) = (/ 5,6,7,8,9,10,11,12,13,14 /)
    bb(:,6) = (/ 6,7,8,9,10,11,12,13,14,15 /)
    bb(:,7) = (/ 7,8,9,10,11,12,13,14,15,16 /)
    bb(:,8) = (/ 8,9,10,11,12,13,14,15,16,17 /)
    bb(:,9) = (/ 9,10,11,12,13,14,15,16,17,18 /)
    bb(:,10) = (/ 10,11,12,13,14,15,16,17,18,19 /)
    CALL sparse_solver_test(icol,pcol,val,bb,xx)
    PRINT *, 'bb(:,1) =',bb(:,1)
    PRINT *, 'bb(:,10) =',bb(:,10)

    ! clear previous allocation
    DEALLOCATE(icol,pcol,val,b,x,bb,xx)

    !-----------------------------------------------------------------------------------------

    PRINT *, 'testing 10x10 complex octave matrix'

    CALL load_octave_matrices('../COMMON/testcases/octave_test_sparse.complex', &
         nrow,ncol,nz,icol,pcol,z_val, &
         z_b,z_x,z_bb,z_xx)
    CALL sparse_matmul(icol,pcol,z_val,z_b,z_x)
    PRINT *, 'z_x=',z_x

    z_x = (0.0_dp,0.0_dp)
    CALL sparse_solve(icol,pcol,z_val,z_b,z_x)
    PRINT *, 'z_x=',z_x
    z_b(:) = (/ (1.0_dp,1.0_dp),(2.0_dp,2.0_dp),(3.0_dp,3.0_dp), &
         (4.0_dp,4.0_dp),(5.0_dp,5.0_dp),(6.0_dp,6.0_dp), &
         (7.0_dp,7.0_dp),(8.0_dp,8.0_dp),(9.0_dp,9.0_dp),(10.0_dp,10.0_dp) /)
    CALL sparse_solver_test(icol,pcol,z_val,z_b,z_x)
    PRINT *, 'z_b_out =',z_b

    CALL sparse_matmul(icol,pcol,z_val,z_bb,z_xx)
    PRINT *, 'z_xx(:,1)=',z_xx(:,1)
    PRINT *, 'z_xx(:,10)=',z_xx(:,10)

    z_xx = (0.0_dp,0.0_dp)
    CALL sparse_solve(icol,pcol,z_val,z_bb,z_xx)
    PRINT *, 'z_xx(:,1)=',z_xx(:,1)
    PRINT *, 'z_xx(:,10)=',z_xx(:,10)

    z_bb(:,1) = (/ (1.0_dp,1.0_dp),(2.0_dp,2.0_dp),(3.0_dp,3.0_dp), &
         (4.0_dp,4.0_dp),(5.0_dp,5.0_dp),(6.0_dp,6.0_dp), &
         (7.0_dp,7.0_dp),(8.0_dp,8.0_dp),(9.0_dp,9.0_dp),(10.0_dp,10.0_dp) /)
    z_bb(:,2) = (/ (2.0_dp,2.0_dp),(3.0_dp,3.0_dp),(4.0_dp,4.0_dp), &
         (5.0_dp,5.0_dp),(6.0_dp,6.0_dp),(7.0_dp,7.0_dp), &
         (8.0_dp,8.0_dp),(9.0_dp,9.0_dp),(10.0_dp,10.0_dp),(11.0_dp,11.0_dp) /)
    z_bb(:,3) = (/ (3.0_dp,3.0_dp),(4.0_dp,4.0_dp),(5.0_dp,5.0_dp), &
         (6.0_dp,6.0_dp),(7.0_dp,7.0_dp),(8.0_dp,8.0_dp), &
         (9.0_dp,9.0_dp),(10.0_dp,10.0_dp),(11.0_dp,11.0_dp),(12.0_dp,12.0_dp) /)
    z_bb(:,4) = (/ (4.0_dp,4.0_dp),(5.0_dp,5.0_dp),(6.0_dp,6.0_dp), &
         (7.0_dp,7.0_dp),(8.0_dp,8.0_dp),(9.0_dp,9.0_dp), &
         (10.0_dp,10.0_dp),(11.0_dp,11.0_dp),(12.0_dp,12.0_dp),(13.0_dp,13.0_dp) /)
    z_bb(:,5) = (/ (5.0_dp,5.0_dp),(6.0_dp,6.0_dp),(7.0_dp,7.0_dp), &
         (8.0_dp,8.0_dp),(9.0_dp,9.0_dp),(10.0_dp,10.0_dp), &
         (11.0_dp,11.0_dp),(12.0_dp,12.0_dp),(13.0_dp,13.0_dp),(14.0_dp,14.0_dp) /)
    z_bb(:,6) = (/ (6.0_dp,6.0_dp),(7.0_dp,7.0_dp),(8.0_dp,8.0_dp), &
         (9.0_dp,9.0_dp),(10.0_dp,10.0_dp),(11.0_dp,11.0_dp), &
         (12.0_dp,12.0_dp),(13.0_dp,13.0_dp),(14.0_dp,14.0_dp),(15.0_dp,15.0_dp) /)
    z_bb(:,7) = (/ (7.0_dp,7.0_dp),(8.0_dp,8.0_dp),(9.0_dp,9.0_dp), &
         (10.0_dp,10.0_dp),(11.0_dp,11.0_dp),(12.0_dp,12.0_dp), &
         (13.0_dp,13.0_dp),(14.0_dp,14.0_dp),(15.0_dp,15.0_dp),(16.0_dp,16.0_dp) /)
    z_bb(:,8) = (/ (8.0_dp,8.0_dp),(9.0_dp,9.0_dp),(10.0_dp,10.0_dp), &
         (11.0_dp,11.0_dp),(12.0_dp,12.0_dp),(13.0_dp,13.0_dp), &
         (14.0_dp,14.0_dp),(15.0_dp,15.0_dp),(16.0_dp,16.0_dp),(17.0_dp,17.0_dp) /)
    z_bb(:,9) = (/ (9.0_dp,9.0_dp),(10.0_dp,10.0_dp),(11.0_dp,11.0_dp), &
         (12.0_dp,12.0_dp),(13.0_dp,13.0_dp),(14.0_dp,14.0_dp), &
         (15.0_dp,15.0_dp),(16.0_dp,16.0_dp),(17.0_dp,17.0_dp),(18.0_dp,18.0_dp) /)
    z_bb(:,10) = (/ (10.0_dp,10.0_dp),(11.0_dp,11.0_dp),(12.0_dp,12.0_dp), &
         (13.0_dp,13.0_dp),(14.0_dp,14.0_dp),(15.0_dp,15.0_dp), &
         (16.0_dp,16.0_dp),(17.0_dp,17.0_dp),(18.0_dp,18.0_dp),(19.0_dp,19.0_dp) /)

    CALL sparse_solver_test(icol,pcol,z_val,z_bb,z_xx)
    PRINT *, 'z_bb(:,1) =',z_bb(:,1)
    PRINT *, 'z_bb(:,10) =',z_bb(:,10)

    !-----------------------------------------------------------------------------------------

    PRINT *, 'testing compressed row format example (real)'
    !testing a compressed row example

    CALL load_compressed_example('../COMMON/testcases/example.dat', &
         nrow,ncol,nz,icol,pcol,val)

    ! conversion test column_pointer2full, if you have data
    ! compressed column format pointer
    CALL column_pointer2full(pcol,A)
    PRINT *, 'test full'
    PRINT '(12f6.1)', A

    ! conversion test column_full2pointer, if you have data
    ! full column format, but want the pointer
    ALLOCATE(A(nrow,ncol))
    CALL sparse2full(icol,pcol,val,A)
    PRINT *, 'test full2pointer'
    CALL column_full2pointer(A,pcol)
    PRINT *, 'pcol  ',pcol

    !-----------------------------------------------------------------------------------------

    PRINT *, 'testing loaded standard example (real)'
    !testing to load a standard example

    CALL load_standard_example('../COMMON/testcases/bcsstk14.dat', &
         nrow,ncol,nz,icol,pcol,val)

    !-----------------------------------------------------------------------------------------

    PRINT *, 'ending sparse_example successfully'

    ! final clean-up
    IF (ALLOCATED(icol)) DEALLOCATE(icol)
    IF (ALLOCATED(pcol)) DEALLOCATE(pcol)
    IF (ALLOCATED(val)) DEALLOCATE(val)
    IF (ALLOCATED(A)) DEALLOCATE(A)
    IF (ALLOCATED(b)) DEALLOCATE(b)
    IF (ALLOCATED(x)) DEALLOCATE(x)
    IF (ALLOCATED(bb)) DEALLOCATE(bb)
    IF (ALLOCATED(xx)) DEALLOCATE(xx)
    IF (ALLOCATED(z_val)) DEALLOCATE(z_val)
    IF (ALLOCATED(z_A)) DEALLOCATE(z_A)
    IF (ALLOCATED(z_b)) DEALLOCATE(z_b)
    IF (ALLOCATED(z_x)) DEALLOCATE(z_x)
    IF (ALLOCATED(z_bb)) DEALLOCATE(z_bb)
    IF (ALLOCATED(z_xx)) DEALLOCATE(z_xx)

  END SUBROUTINE sparse_example
  !-------------------------------------------------------------------------------

END MODULE sparse_mod