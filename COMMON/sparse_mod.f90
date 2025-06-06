MODULE sparse_mod

  IMPLICIT NONE

  PUBLIC sparse_solve_method
  INTEGER :: sparse_solve_method = 3

  PUBLIC sparse_talk
  LOGICAL :: sparse_talk = .FALSE.

  PRIVATE dp
  INTEGER, PARAMETER :: dp = KIND(1.0d0)

  PRIVATE long
  INTEGER, PARAMETER :: long = 8

  PRIVATE factorization_exists
  LOGICAL :: factorization_exists = .FALSE.

  !-------------------------------------------------------------------------------
  !Initialization of the parameters of Super_LU c-Routines
  PRIVATE factors
  INTEGER(kind=long) :: factors
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  !Initialization of the SuiteSparse-Solver-Routine!
  !Solver's internal data adress pointer
  INTEGER(kind=long), PRIVATE :: symbolic, numeric
  !Solves A*x=b (e.g. sys=2 -> solves (A^T)*x=b; further options manual pg. 26)
  INTEGER(kind=long), PRIVATE :: sys=0
  !default values for control pg. 22
  REAL(kind=dp), PRIVATE :: control(20), info_suitesparse(90)
  !-------------------------------------------------------------------------------

  PUBLIC load_mini_example
  PRIVATE load_mini_ex
  INTERFACE load_mini_example
     MODULE PROCEDURE load_mini_ex
  END INTERFACE load_mini_example

  PUBLIC load_compressed_example
  PRIVATE load_compressed_ex
  INTERFACE load_compressed_example
     MODULE PROCEDURE load_compressed_ex
  END INTERFACE load_compressed_example

  PUBLIC load_standard_example
  PRIVATE load_standard_ex
  INTERFACE load_standard_example
     MODULE PROCEDURE load_standard_ex
  END INTERFACE load_standard_example

  PUBLIC load_octave_matrices
  PRIVATE load_octave_mat
  INTERFACE load_octave_matrices
     MODULE PROCEDURE load_octave_mat, load_octave_matComplex
  END INTERFACE load_octave_matrices

  PUBLIC column_pointer2full
  PRIVATE col_pointer2full
  INTERFACE column_pointer2full
     MODULE PROCEDURE col_pointer2full
  END INTERFACE column_pointer2full

  PUBLIC column_full2pointer
  PRIVATE col_full2pointer
  INTERFACE column_full2pointer
     MODULE PROCEDURE col_full2pointer
  END INTERFACE column_full2pointer

  PUBLIC sparse2full
  PRIVATE sp2full
  INTERFACE sparse2full
     MODULE PROCEDURE sp2full, sp2fullComplex
  END INTERFACE sparse2full

  PUBLIC full2sparse
  PRIVATE full2sp
  INTERFACE full2sparse
     MODULE PROCEDURE full2sp,full2spComplex
  END INTERFACE full2sparse

  PUBLIC sparse_solve
  INTERFACE sparse_solve
     MODULE PROCEDURE sparse_solveReal_b1,sparse_solveReal_b2,sparse_solveReal_A_b1,sparse_solveReal_A_b2, &
          sparse_solveComplex_b1,sparse_solveComplex_b2,sparse_solveComplex_A_b1,sparse_solveComplex_A_b2
  END INTERFACE sparse_solve

  PUBLIC sparse_solve_suitesparse
  INTERFACE sparse_solve_suitesparse
     MODULE PROCEDURE sparse_solve_suitesparse_b1, sparse_solve_suitesparse_b2_loop, &
          sparse_solve_suitesparseComplex_b1, sparse_solve_suitesparseComplex_b2_loop
  END INTERFACE sparse_solve_suitesparse

  PUBLIC sparse_matmul
  INTERFACE sparse_matmul
     MODULE PROCEDURE sp_matmul_A_b1,sp_matmul_b1,sp_matmul_A_b2,sp_matmul_b2, &
          sp_matmulComplex_A_b1, sp_matmulComplex_b1, sp_matmulComplex_A_b2, sp_matmulComplex_b2
  END INTERFACE sparse_matmul

  PUBLIC sparse_solver_test
  INTERFACE sparse_solver_test
     MODULE PROCEDURE sp_test_A_b1,sp_test_b1,sp_test_A_b2,sp_test_b2, &
          sp_testComplex_A_b1, sp_testComplex_b1, sp_testComplex_A_b2, sp_testComplex_b2
  END INTERFACE sparse_solver_test

  PUBLIC sparse_example

  PUBLIC remap_rc
  INTERFACE remap_rc
     MODULE PROCEDURE remap_rc_real, remap_rc_cmplx
  END INTERFACE remap_rc

  ! helper
  PRIVATE find_unit


CONTAINS

  !-------------------------------------------------------------------------------
  ! finds free unit
  SUBROUTINE find_unit(unit)
    INTEGER, INTENT(inout) :: unit
    LOGICAL :: opened
    DO
       INQUIRE(unit=unit,opened=opened)
       IF (.NOT. opened) EXIT
       unit = unit + 1
    END DO

  END SUBROUTINE find_unit
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  ! Examples
  SUBROUTINE sparse_example(example,subexample)
    INTEGER, INTENT(in) :: example
    INTEGER, INTENT(in), OPTIONAL :: subexample

    CHARACTER(len=100) :: name
    INTEGER :: nrow,ncol,nz, nrhs
    INTEGER, DIMENSION(:), ALLOCATABLE :: irow,pcol,icol
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: val,b,x
    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: A,bb,xx
    COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE :: z_val,z_b,z_x
    COMPLEX(kind=dp), DIMENSION(:,:), ALLOCATABLE :: z_A,z_bb,z_xx

    INTEGER :: ir,ic,icmax,subex_example6,i,unit

    nrhs = 10

    subex_example6=1
    IF(PRESENT(subexample)) subex_example6=subexample

    IF (example .EQ. 1) THEN
       ! load the test-matrix for the mini_example
       CALL load_mini_example(A)

       ! construct the rhs
       IF (ALLOCATED(b)) DEALLOCATE(b)
       ALLOCATE(b(SIZE(A,2)))
       b = 1.0_dp
       ! x is only needed because b should not be overwritten
       IF (ALLOCATED(x)) DEALLOCATE(x)
       ALLOCATE(x(SIZE(b,1)))
       x = b
       ! solve
       CALL sparse_solve(A,x)
       PRINT *,x
       ! test
       CALL sparse_solver_test(A,x,b)

    ELSEIF (example .EQ. 2) THEN
       ! load the test-matrix for the mini_example (with multiple rhs)
       CALL load_mini_example(A)
       ! convert to sparse
       CALL full2sparse(A,irow,pcol,val,nrow,ncol,nz)
       IF (sparse_talk) PRINT *, 'nrow=',nrow,' ncol=',ncol,' nz=',nz

       !construct an array of rhs
       IF (ALLOCATED(bb)) DEALLOCATE(bb)
       icmax = ncol
       ALLOCATE(bb(nrow,icmax))
       DO ic = 1, icmax
          DO ir = 1, nrow
             bb(ir,ic) = ir-1 + 10*ic
          END DO
       END DO
       IF(ALLOCATED(xx)) DEALLOCATE(xx)
       ALLOCATE(xx(SIZE(bb,1),SIZE(bb,2)))
       xx = bb
       ! solve the system for multiple rhs
       !CALL sparse_solve(nrow,ncol,nz,irow,pcol,val,xx)
       ! would also work with a full column index vector
       CALL column_pointer2full(pcol,icol)
       CALL sparse_solve(nrow,ncol,nz,irow,icol,val,xx)
       ! test
       CALL sparse_solver_test(nrow,ncol,irow,pcol,val,xx,bb)

    ELSEIF (example .EQ. 3) THEN
       ! load the test-matrix g10
       name = '/proj/plasma/Libs/SuperLU/SuperLU_3.0/DATA/g10'
       CALL load_standard_example(name,nrow,ncol,nz,irow,pcol,val)
       IF (sparse_talk) PRINT *, 'nrow=',nrow,' ncol=',ncol,' nz=',nz

       ! construct the rhs
       IF (ALLOCATED(b)) DEALLOCATE(b)
       ALLOCATE(b(nrow))
       b = 1.0_dp
       IF (ALLOCATED(x)) DEALLOCATE(x)
       ALLOCATE(x(SIZE(b,1)))
       x = b
       ! solve
       CALL sparse_solve(nrow,ncol,nz,irow,pcol,val,x)
       PRINT *,x
       ! test
       CALL sparse_solver_test(nrow,ncol,irow,pcol,val,x,b)

    ELSEIF (example .EQ. 4) THEN
       ! load the test-matrix of the compressed_example
       name = '/proj/plasma/Libs/SuperLU/SuperLU_3.0/DATA/sparse_compressed_e100_s100_D0d001.dat'
       CALL load_compressed_example(name,nrow,ncol,nz,irow,pcol,val)
       IF (sparse_talk) PRINT *, 'nrow=',nrow,' ncol=',ncol,' nz=',nz

       ! construct a rhs
       IF (ALLOCATED(b)) DEALLOCATE(b)
       ALLOCATE(b(nrow))
       b = 0.0_dp
       b(1) = 1.0_dp
       IF(ALLOCATED(x)) DEALLOCATE(x)
       ALLOCATE(x(SIZE(b,1)))
       x = b
       ! solve
       CALL sparse_solve(nrow,ncol,nz,irow,pcol,val,x)
       PRINT *,x
       ! test
       CALL sparse_solver_test(nrow,ncol,irow,pcol,val,x,b)

    ELSEIF (example .EQ. 5) THEN
       ! load the test-matrix of the compressed_example
       name = "/proj/plasma/Libs/SuperLU/SuperLU_3.0/DATA/&
            &sparse_compressed_e100_s100_D0d001.dat"
       CALL load_compressed_example(name,nrow,ncol,nz,irow,pcol,val)
       IF (sparse_talk) PRINT *, 'nrow=',nrow,' ncol=',ncol,' nz=',nz
       ! construct a rhs
       IF (ALLOCATED(bb)) DEALLOCATE(bb)
       ALLOCATE(bb(nrow,ncol))
       DO ir = 1, nrow
          bb(ir,ir) = 1.0_dp
       END DO
       IF(ALLOCATED(xx)) DEALLOCATE(xx)
       ALLOCATE(xx(SIZE(bb,1),SIZE(bb,2)))
       xx = bb
       ! solve
       CALL sparse_solve(nrow,ncol,nz,irow,pcol,val,xx)
       ! test
       CALL sparse_solver_test(nrow,ncol,irow,pcol,val,xx,bb)

    ELSEIF (example .EQ. 6) THEN
       ! load the different test-matrices generated by octave
       ! for the different test-cases
       SELECT CASE (subex_example6)
       CASE (1)
          name = '/proj/plasma/Solver_Test/TestMatrices/test_matrix1.dat'
       CASE (2)
          name = '/proj/plasma/Solver_Test/TestMatrices/test_matrix2.dat'
       CASE (3)
          name = '/proj/plasma/Solver_Test/TestMatrices/test_matrix3.dat'
       CASE (4)
          name = '/proj/plasma/Solver_Test/TestMatrices/test_matrix4.dat'
       CASE (5)
          name = '/proj/plasma/Solver_Test/TestMatrices/test_matrix5.dat'
       CASE (6)
          name = '/proj/plasma/Solver_Test/TestMatrices/test_matrix6.dat'
       CASE (7)
          name = '/proj/plasma/Solver_Test/TestMatrices/test_matrix7.dat'
       CASE (8)
          name = '/proj/plasma/Solver_Test/TestMatrices/test_matrix8.dat'
       CASE (9)
          name = '/proj/plasma/Solver_Test/TestMatrices/test_matrix9.dat'
       CASE (10)
          name = '/proj/plasma/Solver_Test/TestMatrices/test_matrix10.dat'
       CASE (11)
          name = '/proj/plasma/Solver_Test/TestMatrices/test_matrix11.dat'
       CASE (12)
          name = '/proj/plasma/Solver_Test/TestMatrices/test_matrix12.dat'
       CASE (13)
          name = '/proj/plasma/Solver_Test/TestMatrices/test_matrix13.dat'
       CASE DEFAULT
          PRINT *, 'unknown file name -> select a subexample between 1 and 13'
       END SELECT
       CALL load_octave_matrices(name,nrow,ncol,nz,irow,pcol,val)
       IF (sparse_talk) PRINT *, 'nrow=',nrow,' ncol=',ncol,' nz=',nz

       ! construct a rhs
       IF (ALLOCATED(b)) DEALLOCATE(b)
       ALLOCATE(b(nrow))
       b = 0.0_dp
       b(1) = 1.0_dp
       IF(ALLOCATED(x)) DEALLOCATE(x)
       ALLOCATE(x(SIZE(b,1)))
       x = b
       ! solve
       CALL sparse_solve(nrow,ncol,nz,irow,pcol,val,x)

       ! test
       CALL sparse_solver_test(nrow,ncol,irow,pcol,val,x,b)


    ELSEIF (example .EQ. 7) THEN

       nrow=8
       ncol=8
       nz=20

       IF (ALLOCATED(pcol)) DEALLOCATE(pcol)
       ALLOCATE(pcol(nrow+1))
       IF (ALLOCATED(irow)) DEALLOCATE(irow)
       ALLOCATE(irow(nz))
       IF (ALLOCATED(z_val)) DEALLOCATE(z_val)
       ALLOCATE(z_val(nz))

       pcol= (/1,5,8,10,12,13,16,18,21/)

       irow =(/ 1,3,6,7,2,3,5,3,8,4,7,2,3,6,8,2,7,3,7,8 /)
       z_val=(/ (7.d0, 1.d0), (1.d0,1.d0), (2.d0,1.d0), (7.d0,1.d0), (-4.d0,0.d0),&
            (8.d0,1.d0), (2.d0,1.d0),(1.d0,1.d0),(5.d0,1.d0),(7.d0,0.d0),  (9.d0,1.d0),&
            (-4d0,1.d0),(7.d0,1.d0),  (3.d0,1.d0), (8.d0,0.d0),(1.d0,1.d0),&
            (11.d0,1.d0),(-3.d0,1.d0), (2.d0,1.d0), (5.d0,0.d0)/)

       IF (ALLOCATED(z_b)) DEALLOCATE(z_b)
       ALLOCATE(z_b(nrow))
       DO ir = 1, nrow
          z_b(ir) = (1.d0,1.d0)
       END DO

       IF (ALLOCATED(z_x)) DEALLOCATE(z_x)
       ALLOCATE(z_x(nrow))
       z_x=z_b

       CALL sparse_solve(nrow,ncol,nz,irow,pcol,z_val,z_x)
       PRINT *,z_x
       CALL sparse_solver_test(nrow,ncol,irow,pcol,z_val,z_x,z_b)

    ELSEIF (example .EQ. 8) THEN
       ! load the different complex test-matrices generated by octave
       ! for the different test-cases
       SELECT CASE (subex_example6)
       CASE (1)
          name = '/proj/plasma/Solver_Test/TestMatrices/test_ComplexMatrix1.dat'
       CASE (2)
          name = '/proj/plasma/Solver_Test/TestMatrices/test_ComplexMatrix2.dat'
       CASE (3)
          name= '/proj/plasma/Solver_Test/TestMatrices/test_ComplexMatrix3.dat'
       CASE (4)
          name = '/proj/plasma/Solver_Test/TestMatrices/test_ComplexMatrix4.dat'
       CASE (5)
          name = '/proj/plasma/Solver_Test/TestMatrices/test_ComplexMatrix5.dat'
       CASE (6)
          name = '/proj/plasma/Solver_Test/TestMatrices/test_ComplexMatrix6.dat'
       CASE (7)
          name = '/proj/plasma/Solver_Test/TestMatrices/test_ComplexMatrix7.dat'
       CASE DEFAULT
          PRINT *, 'unknown file name -> select a subexample between 1 and 7'
       END SELECT
       CALL load_octave_matrices(name,nrow,ncol,nz,irow,pcol,z_val)
       IF (sparse_talk) PRINT *, 'nrow=',nrow,' ncol=',ncol,' nz=',nz


       ! construct a rhs
       IF (ALLOCATED(z_b)) DEALLOCATE(z_b)
       ALLOCATE(z_b(nrow))
       z_b = (0.0_dp,0.0_dp)
       z_b(1) = (1.0_dp,1.0_dp)
       IF(ALLOCATED(z_x)) DEALLOCATE(z_x)
       ALLOCATE(z_x(SIZE(z_b,1)))
       z_x = z_b
       ! solve
       CALL sparse_solve(nrow,ncol,nz,irow,pcol,z_val,z_x)
       ! test
       CALL sparse_solver_test(nrow,ncol,irow,pcol,z_val,z_x,z_b)

    ELSEIF (example .EQ. 9) THEN
       !load test_matrix13.dat and solve the system for different numbers
       !of right-hand sides

       name = '/proj/plasma/Solver_Test/TestMatrices/test_matrix13.dat'
       CALL load_octave_matrices(name,nrow,ncol,nz,irow,pcol,val)
       IF (sparse_talk) PRINT *, 'nrow=',nrow,' ncol=',ncol,' nz=',nz

       SELECT CASE (subex_example6)
       CASE (1)
          nrhs=10
       CASE (2)
          nrhs=25
       CASE (3)
          nrhs=50
       CASE (4)
          nrhs=75
       CASE (5)
          nrhs=100
       CASE (6)
          nrhs=250
       CASE (7)
          nrhs=500
       CASE (8)
          nrhs=750
       CASE (9)
          nrhs=1000
       CASE DEFAULT
          PRINT *, 'unknown number of rhs -> select a subexample between 1 and 9'
       END SELECT

       ! construct a rhs
       IF (ALLOCATED(bb)) DEALLOCATE(bb)
       ALLOCATE(bb(nrow,nrhs))
       DO ir = 1, nrhs
          bb(ir,ir) = 1.0_dp
       END DO
       IF(ALLOCATED(xx)) DEALLOCATE(xx)
       ALLOCATE(xx(SIZE(bb,1),SIZE(bb,2)))
       xx = bb
       ! solve
       CALL sparse_solve(nrow,ncol,nz,irow,pcol,val,xx)
       ! test
       CALL sparse_solver_test(nrow,ncol,irow,pcol,val,xx,bb)

    END IF

    IF (ALLOCATED(irow)) DEALLOCATE(irow)
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

  !-------------------------------------------------------------------------------
  ! loads a mini example
  SUBROUTINE load_mini_ex(A)
    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(out) :: A

    ALLOCATE(A(5,5))
    A(:,1) = (/1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp/)
    A(:,2) = A(:,1)*5 + 2
    A(:,3) = A(:,2)*7 + 2
    A(:,4) = A(:,3)*2 + 2
    A(:,5) = A(:,4)*9 + 2

    A(2,4) = 0.0_dp
    A(3,3) = 0.0_dp
    A(4,2) = 0.0_dp

  END SUBROUTINE load_mini_ex
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  ! loads own compressed example
  SUBROUTINE load_compressed_ex(name,nrow,ncol,nz,irow,pcol,val)
    CHARACTER(LEN=*), INTENT(in) :: name
    INTEGER, INTENT(out) :: nrow,ncol,nz
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(out) :: irow,pcol
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(out) :: val

    INTEGER :: unit,i

    unit = 10;
    CALL find_unit(unit)
    OPEN(unit=unit,file=TRIM(ADJUSTL(name)),action='read')

    READ(unit,*) nrow,ncol,nz
    ALLOCATE(irow(nz),pcol(ncol+1),val(nz))
    READ(unit,*) (irow(i), i = 1, nz)
    READ(unit,*) (pcol(i), i = 1, ncol+1)
    READ(unit,*) (val(i),  i = 1, nz)

    CLOSE(unit=unit)

  END SUBROUTINE load_compressed_ex
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  ! loads standard example from SuperLU distribution
  SUBROUTINE load_standard_ex(name,nrow,ncol,nz,irow,pcol,val)
    CHARACTER(LEN=*), INTENT(in) :: name
    INTEGER, INTENT(out) :: nrow,ncol,nz
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(out) :: irow,pcol
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(out) :: val

    INTEGER :: unit,i

    CHARACTER(len=72) :: fmt1
    CHARACTER(len=72) :: title
    CHARACTER(len=8)  :: key
    CHARACTER(len=3)  :: mxtype
    CHARACTER(len=16) :: ptrfmt,indfmt
    CHARACTER(len=20) :: valfmt,rhsfmt

    INTEGER :: totcrd,ptrcrd,indcrd,valcrd,rhscrd,neltvl

    fmt1 = '( A72, A8 / 5I14 / A3, 11X, 4I14 / 2A16, 2A20 )'

    unit = 10;
    CALL find_unit(unit)
    OPEN(unit=unit,file=TRIM(ADJUSTL(name)),action='read')

    READ (unit=unit,fmt=fmt1 ) &
         title, key, totcrd, ptrcrd, indcrd, valcrd, rhscrd, &
         mxtype, nrow, ncol, nz, neltvl, &
         ptrfmt, indfmt, valfmt, rhsfmt
    ALLOCATE(irow(nz),pcol(ncol+1),val(nz))
    READ (unit=unit,fmt=ptrfmt) ( pcol(i), i = 1, ncol+1 )
    READ (unit=unit,fmt=indfmt) ( irow(i), i = 1, nz )
    READ (unit=unit,fmt=valfmt) ( val(i),  i = 1, nz )

    CLOSE(unit=unit)

  END SUBROUTINE load_standard_ex
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE load_octave_mat(name,nrow,ncol,nz,irow,pcol,val)
    CHARACTER(LEN=*), INTENT(in) :: name
    INTEGER, INTENT(out) :: nrow,ncol,nz
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(out) :: irow,pcol
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(out) :: val

    INTEGER :: unit,i,k
    INTEGER, DIMENSION(:), ALLOCATABLE :: octave_pcol

    !open the input-file ("name")
    unit = 10;
    CALL find_unit(unit)
    OPEN(unit=unit,file=TRIM(ADJUSTL(name)),action='read')

    !read nrow, ncol, nz and allocate the arrays for
    !irow, pcol val
    READ(unit,*) nrow,ncol,nz
    ALLOCATE(irow(nz),pcol(ncol+1),octave_pcol(nz),val(nz))
    !read the sparse matrix (Octave-format)
    !storage-format for sparse matrices in ocatave
    !uses the coordinates (irow, octave_pcol) of entries (val)
    !in matrix
    DO i=1,nz
       READ(unit,*) irow(i),octave_pcol(i),val(i)
    END DO
    CLOSE(unit=unit)

    !now calculate the index of the first entry (linear index)
    !of each row (pcol)
    !first step: calculate the number of entries in each row
    pcol(1)=octave_pcol(1)
    k=1
    DO i=1,ncol
       IF (k .GT. nz) EXIT
       IF (octave_pcol(k) .EQ. i) THEN
          DO WHILE (octave_pcol(k) .EQ. i)
             pcol(i+1)=pcol(i+1)+1
             k=k+1
             IF (k .GT. nz) EXIT
          END DO
          k=k-1
       ELSE
          CYCLE
       END IF
       k=k+1
    END DO
    !second step: sum over the number of entries in each row
    !to get desired the linear index
    DO i=1,ncol
       pcol(i+1)=pcol(i)+pcol(i+1)
    END DO

  END SUBROUTINE load_octave_mat
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE load_octave_matComplex(name,nrow,ncol,nz,irow,pcol,val)
    CHARACTER(LEN=*), INTENT(in) :: name
    INTEGER, INTENT(out) :: nrow,ncol,nz
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(out) :: irow,pcol
    COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(out) :: val

    INTEGER :: unit,i,k
    INTEGER, DIMENSION(:), ALLOCATABLE :: octave_pcol

    !open the input-file ("name")
    unit = 10;
    CALL find_unit(unit)
    OPEN(unit=unit,file=TRIM(ADJUSTL(name)),action='read')

    !read nrow, ncol, nz and allocate the arrays for
    !irow, pcol val
    READ(unit,*) nrow,ncol,nz
    ALLOCATE(irow(nz),pcol(ncol+1),octave_pcol(nz),val(nz))
    !read the sparse matrix (Octave-format)
    !storage-format for sparse matrices in ocatave
    !uses the coordinates (irow, octave_pcol) of entries (val)
    !in matrix
    DO i=1,nz
       READ(unit,*) irow(i),octave_pcol(i),val(i)
    END DO
    CLOSE(unit=unit)

    !now calculate the index of the first entry (linear index)
    !of each row (pcol)
    !first step: calculate the number of entries in each row
    pcol(1)=octave_pcol(1)
    k=1
    DO i=1,ncol
       IF (k .GT. nz) EXIT
       IF (octave_pcol(k) .EQ. i) THEN
          DO WHILE (octave_pcol(k) .EQ. i)
             pcol(i+1)=pcol(i+1)+1
             k=k+1
             IF (k .GT. nz) EXIT
          END DO
          k=k-1
       ELSE
          CYCLE
       END IF
       k=k+1
    END DO
    !second step: sum over the number of entries in each row
    !to get desired the linear index
    DO i=1,ncol
       pcol(i+1)=pcol(i)+pcol(i+1)
    END DO

    RETURN

  END SUBROUTINE load_octave_matComplex
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  ! solves A*x = b for sparse A and 1-D vector b
  ! A is specified through nrow,ncol,nz,irow,pcol,val
  ! results are returned in b
  SUBROUTINE sparse_solveReal_b1(nrow,ncol,nz,irow,pcol,val,b,iopt_in)
    INTEGER, INTENT(in) :: nrow,ncol,nz
    INTEGER, DIMENSION(:), INTENT(in) :: irow,pcol
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: val
    REAL(kind=dp), DIMENSION(:), INTENT(inout) :: b
    INTEGER, INTENT(in), OPTIONAL :: iopt_in

    INTEGER :: iopt = 0
    INTEGER, DIMENSION(:), ALLOCATABLE :: pcoln
    LOGICAL :: pcol_modified = .FALSE.

    !optional input
    IF (PRESENT(iopt_in)) iopt = iopt_in

    ! make sure that pcol is a pointer, otherwise create pcoln
    IF (SIZE(pcol,1) .EQ. SIZE(irow,1)) THEN
       CALL column_full2pointer(pcol,pcoln)
       pcol_modified = .TRUE.
    END IF

    ! check about existing factorization
    IF (factorization_exists .AND. iopt .EQ. 1) THEN ! free memory first
       IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. 3) ) THEN
          IF (pcol_modified) THEN
             CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcoln,val,b,3)
          ELSE
             CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,3)
          END IF
       END IF
    END IF
    IF (.NOT. factorization_exists .AND. iopt .EQ. 2) THEN ! factorize first
       IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. 3) ) THEN
          IF (pcol_modified) THEN
             CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcoln,val,b,1)
          ELSE
             CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,1)
          END IF
       END IF
       factorization_exists = .TRUE.
    END IF
    IF (iopt .EQ. 1) factorization_exists = .TRUE.
    IF (iopt .EQ. 3) factorization_exists = .FALSE.

    IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. 3) ) THEN
       IF (pcol_modified) THEN
          CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcoln,val,b,iopt)
       ELSE
          CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,iopt)
       END IF
    ELSE
       PRINT *, 'sparse_solve_method ',sparse_solve_method,'not implemented'
       STOP
    END IF

    IF (ALLOCATED(pcoln)) DEALLOCATE(pcoln)
    RETURN
  END SUBROUTINE sparse_solveReal_b1
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  ! solves A*x = b for sparse A and 1-D vector b
  ! A is specified through nrow,ncol,nz,irow,pcol,val
  ! results are returned in b
  SUBROUTINE sparse_solveComplex_b1(nrow,ncol,nz,irow,pcol,val,b,iopt_in)
    INTEGER, INTENT(in) :: nrow,ncol,nz
    INTEGER, DIMENSION(:), INTENT(in) :: irow,pcol
    COMPLEX(kind=dp), DIMENSION(:), INTENT(in) :: val
    COMPLEX(kind=dp), DIMENSION(:), INTENT(inout) :: b
    INTEGER, INTENT(in), OPTIONAL :: iopt_in

    INTEGER :: iopt = 0
    INTEGER, DIMENSION(:), ALLOCATABLE :: pcoln
    LOGICAL :: pcol_modified = .FALSE.

    !optional input
    IF (PRESENT(iopt_in)) iopt = iopt_in

    ! make sure that pcol is a pointer, otherwise create pcoln
    IF (SIZE(pcol,1) .EQ. SIZE(irow,1)) THEN
       CALL column_full2pointer(pcol,pcoln)
       pcol_modified = .TRUE.
    END IF

    ! check about existing factorization
    IF (factorization_exists .AND. iopt .EQ. 1) THEN ! free memory first
       IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. 3) ) THEN
          IF (pcol_modified) THEN
             CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcoln,val,b,3)
          ELSE
             CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,3)
          END IF
       END IF
    END IF
    IF (.NOT. factorization_exists .AND. iopt .EQ. 2) THEN ! factorize first
       IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. 3) ) THEN
          IF (pcol_modified) THEN
             CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcoln,val,b,1)
          ELSE
             CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,1)
          END IF
       END IF
       factorization_exists = .TRUE.
    END IF
    IF (iopt .EQ. 1) factorization_exists = .TRUE.
    IF (iopt .EQ. 3) factorization_exists = .FALSE.

    IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. 3) ) THEN
       IF (pcol_modified) THEN
          CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcoln,val,b,iopt)
       ELSE
          CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,iopt)
       END IF
    ELSE
       PRINT *, 'sparse_solve_method ',sparse_solve_method,'not implemented'
       STOP
    END IF

    IF (ALLOCATED(pcoln)) DEALLOCATE(pcoln)

  END SUBROUTINE sparse_solveComplex_b1
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  ! solves A*x = b for sparse A and 1-D array b
  ! A is specified through nrow,ncol,nz,irow,pcol,val
  ! results are returned in b
  SUBROUTINE sparse_solveReal_b2(nrow,ncol,nz,irow,pcol,val,b,iopt_in)
    INTEGER, INTENT(in) :: nrow,ncol,nz
    INTEGER, DIMENSION(:), INTENT(in) :: irow,pcol
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: val
    REAL(kind=dp), DIMENSION(:,:), INTENT(inout) :: b
    INTEGER, INTENT(in), OPTIONAL :: iopt_in

    INTEGER :: iopt = 0
    INTEGER, DIMENSION(:), ALLOCATABLE :: pcoln
    LOGICAL :: pcol_modified = .FALSE.

    !optional input
    IF (PRESENT(iopt_in)) iopt = iopt_in

    ! make sure that pcol is a pointer, otherwise create pcoln
    IF (SIZE(pcol,1) .EQ. SIZE(irow,1)) THEN
       CALL column_full2pointer(pcol,pcoln)
       pcol_modified = .TRUE.
    END IF

    ! check about existing factorization
    IF (factorization_exists .AND. iopt .EQ. 1) THEN ! free memory first
       IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. 3) ) THEN
          IF (pcol_modified) THEN
             CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcoln,val,b,3)
          ELSE
             CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,3)
          END IF
       END IF
    END IF
    IF (.NOT. factorization_exists .AND. iopt .EQ. 2) THEN ! factorize first
       IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. 3) ) THEN
          IF (pcol_modified) THEN
             CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcoln,val,b,1)
          ELSE
             CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,1)
          END IF
       END IF
       factorization_exists = .TRUE.
    END IF
    IF (iopt .EQ. 1) factorization_exists = .TRUE.
    IF (iopt .EQ. 3) factorization_exists = .FALSE.

    IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. 3) ) THEN
       IF (pcol_modified) THEN
          CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcoln,val,b,iopt)
       ELSE
          CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,iopt)
       END IF
    ELSE
       PRINT *, 'sparse_solve_method ',sparse_solve_method,'not implemented'
       STOP
    END IF

    IF (ALLOCATED(pcoln)) DEALLOCATE(pcoln)

  END SUBROUTINE sparse_solveReal_b2
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  ! solves A*x = b for sparse A and 1-D array b
  ! A is specified through nrow,ncol,nz,irow,pcol,val
  ! results are returned in b
  SUBROUTINE sparse_solveComplex_b2(nrow,ncol,nz,irow,pcol,val,b,iopt_in)
    INTEGER, INTENT(in) :: nrow,ncol,nz
    INTEGER, DIMENSION(:), INTENT(in) :: irow,pcol
    COMPLEX(kind=dp), DIMENSION(:), INTENT(in) :: val
    COMPLEX(kind=dp), DIMENSION(:,:), INTENT(inout) :: b
    INTEGER, INTENT(in), OPTIONAL :: iopt_in

    INTEGER :: iopt = 0
    INTEGER, DIMENSION(:), ALLOCATABLE :: pcoln
    LOGICAL :: pcol_modified = .FALSE.

    !optional input
    IF (PRESENT(iopt_in)) iopt = iopt_in

    ! make sure that pcol is a pointer, otherwise create pcoln
    IF (SIZE(pcol,1) .EQ. SIZE(irow,1)) THEN
       CALL column_full2pointer(pcol,pcoln)
       pcol_modified = .TRUE.
    END IF

    ! check about existing factorization
    IF (factorization_exists .AND. iopt .EQ. 1) THEN ! free memory first
       IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. 3) ) THEN
          IF (pcol_modified) THEN
             CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcoln,val,b,3)
          ELSE
             CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,3)
          END IF
       END IF
    END IF
    IF (.NOT. factorization_exists .AND. iopt .EQ. 2) THEN ! factorize first
       IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. 3) ) THEN
          IF (pcol_modified) THEN
             CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcoln,val,b,1)
          ELSE
             CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,1)
          END IF
       END IF
       factorization_exists = .TRUE.
    END IF
    IF (iopt .EQ. 1) factorization_exists = .TRUE.
    IF (iopt .EQ. 3) factorization_exists = .FALSE.

    IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. 3) ) THEN
       IF (pcol_modified) THEN
          CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcoln,val,b,iopt)
       ELSE
          CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,iopt)
       END IF
    ELSE
       PRINT *, 'sparse_solve_method ',sparse_solve_method,'not implemented'
       STOP
    END IF

    IF (ALLOCATED(pcoln)) DEALLOCATE(pcoln)

  END SUBROUTINE sparse_solveComplex_b2
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  ! solves A*x = b for sparse A and 1-D vector b
  ! A is given as a full matrix
  ! results are returned in b
  SUBROUTINE sparse_solveReal_A_b1(A,b,iopt_in)
    REAL(kind=dp), DIMENSION(:,:), INTENT(in) :: A
    REAL(kind=dp), DIMENSION(:), INTENT(inout) :: b
    INTEGER, INTENT(in), OPTIONAL :: iopt_in

    INTEGER :: iopt = 0
    INTEGER :: nrow,ncol,nz
    INTEGER, DIMENSION(:), ALLOCATABLE :: irow,pcol
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: val

    !optional input
    IF (PRESENT(iopt_in)) iopt = iopt_in

    CALL full2sparse(A,irow,pcol,val,nrow,ncol,nz)

    ! check about existing factorization
    IF (factorization_exists .AND. iopt .EQ. 1) THEN ! free memory first
       IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. 3) ) THEN
          CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,3)
       END IF
    END IF
    IF (.NOT. factorization_exists .AND. iopt .EQ. 2) THEN ! factorize first
       IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. 3) ) THEN
          CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,1)
       END IF
       factorization_exists = .TRUE.
    END IF
    IF (iopt .EQ. 1) factorization_exists = .TRUE.
    IF (iopt .EQ. 3) factorization_exists = .FALSE.

    IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. 3) ) THEN
       CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,iopt)
    ELSE
       PRINT *, 'sparse_solve_method ',sparse_solve_method,'not implemented'
       STOP
    END IF

    IF (ALLOCATED(irow)) DEALLOCATE(irow)
    IF (ALLOCATED(pcol)) DEALLOCATE(pcol)
    IF (ALLOCATED(val))  DEALLOCATE(val)

  END SUBROUTINE sparse_solveReal_A_b1
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  ! solves A*x = b for sparse A and 1-D vector b
  ! A is given as a full matrix
  ! results are returned in b
  SUBROUTINE sparse_solveComplex_A_b1(A,b,iopt_in)
    COMPLEX(kind=dp), DIMENSION(:,:), INTENT(in) :: A
    COMPLEX(kind=dp), DIMENSION(:), INTENT(inout) :: b
    INTEGER, INTENT(in), OPTIONAL :: iopt_in

    INTEGER :: iopt = 0
    INTEGER :: nrow,ncol,nz
    INTEGER, DIMENSION(:), ALLOCATABLE :: irow,pcol
    COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE :: val

    !optional input
    IF (PRESENT(iopt_in)) iopt = iopt_in

    CALL full2sparse(A,irow,pcol,val,nrow,ncol,nz)

    ! check about existing factorization
    IF (factorization_exists .AND. iopt .EQ. 1) THEN ! free memory first
       IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. 3) ) THEN
          CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,3)
       END IF
    END IF
    IF (.NOT. factorization_exists .AND. iopt .EQ. 2) THEN ! factorize first
       IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. 3) ) THEN
          CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,1)
       END IF
       factorization_exists = .TRUE.
    END IF
    IF (iopt .EQ. 1) factorization_exists = .TRUE.
    IF (iopt .EQ. 3) factorization_exists = .FALSE.

    IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. 3) ) THEN
       CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,iopt)
    ELSE
       PRINT *, 'sparse_solve_method ',sparse_solve_method,'not implemented'
       STOP
    END IF

    IF (ALLOCATED(irow)) DEALLOCATE(irow)
    IF (ALLOCATED(pcol)) DEALLOCATE(pcol)
    IF (ALLOCATED(val))  DEALLOCATE(val)

  END SUBROUTINE sparse_solveComplex_A_b1
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  ! solves A*x = b for sparse A and 2-D array b
  ! A is given as a full matrix
  ! results are returned in b
  SUBROUTINE sparse_solveReal_A_b2(A,b,iopt_in)
    REAL(kind=dp), DIMENSION(:,:), INTENT(in) :: A
    REAL(kind=dp), DIMENSION(:,:), INTENT(inout) :: b
    INTEGER, INTENT(in), OPTIONAL :: iopt_in

    INTEGER :: iopt = 0
    INTEGER :: nrow,ncol,nz
    INTEGER, DIMENSION(:), ALLOCATABLE :: irow,pcol
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: val

    !optional input
    IF (PRESENT(iopt_in)) iopt = iopt_in

    CALL full2sparse(A,irow,pcol,val,nrow,ncol,nz)

    ! check about existing factorization
    IF (factorization_exists .AND. iopt .EQ. 1) THEN ! free memory first
       IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. 3) ) THEN
          CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,3)
       END IF
    END IF
    IF (.NOT. factorization_exists .AND. iopt .EQ. 2) THEN ! factorize first
       IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. 3) ) THEN
          CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,1)
       END IF
       factorization_exists = .TRUE.
    END IF
    IF (iopt .EQ. 1) factorization_exists = .TRUE.
    IF (iopt .EQ. 3) factorization_exists = .FALSE.

    IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. 3) ) THEN
       CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,iopt)
    ELSE
       PRINT *, 'sparse_solve_method ',sparse_solve_method,'not implemented'
       STOP
    END IF

    IF (ALLOCATED(irow)) DEALLOCATE(irow)
    IF (ALLOCATED(pcol)) DEALLOCATE(pcol)
    IF (ALLOCATED(val))  DEALLOCATE(val)

  END SUBROUTINE sparse_solveReal_A_b2
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  ! solves A*x = b for sparse A and 2-D array b
  ! A is given as a full matrix
  ! results are returned in b
  SUBROUTINE sparse_solveComplex_A_b2(A,b,iopt_in)
    COMPLEX(kind=dp), DIMENSION(:,:), INTENT(in) :: A
    COMPLEX(kind=dp), DIMENSION(:,:), INTENT(inout) :: b
    INTEGER, INTENT(in), OPTIONAL :: iopt_in

    INTEGER :: iopt = 0
    INTEGER :: nrow,ncol,nz
    INTEGER, DIMENSION(:), ALLOCATABLE :: irow,pcol
    COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE :: val

    !optional input
    IF (PRESENT(iopt_in)) iopt = iopt_in

    CALL full2sparse(A,irow,pcol,val,nrow,ncol,nz)

    ! check about existing factorization
    IF (factorization_exists .AND. iopt .EQ. 1) THEN ! free memory first
       IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. 3) ) THEN
          CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,3)
       END IF
    END IF
    IF (.NOT. factorization_exists .AND. iopt .EQ. 2) THEN ! factorize first
       IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. 3) ) THEN
          CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,1)
       END IF
       factorization_exists = .TRUE.
    END IF
    IF (iopt .EQ. 1) factorization_exists = .TRUE.
    IF (iopt .EQ. 3) factorization_exists = .FALSE.

    IF ( (sparse_solve_method .EQ. 2) .OR. (sparse_solve_method .EQ. 3) ) THEN
       CALL sparse_solve_suitesparse(nrow,ncol,nz,irow,pcol,val,b,iopt)
    ELSE
       PRINT *, 'sparse_solve_method ',sparse_solve_method,'not implemented'
       STOP
    END IF

    IF (ALLOCATED(irow)) DEALLOCATE(irow)
    IF (ALLOCATED(pcol)) DEALLOCATE(pcol)
    IF (ALLOCATED(val))  DEALLOCATE(val)

  END SUBROUTINE sparse_solveComplex_A_b2
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  ! solves A*x = b for sparse A and 1-D vector b
  ! A is specified through nrow,ncol,nz,irow,pcol,val
  ! results are returned in b
  ! Routines from SuperLU-Distribution
  SUBROUTINE sparse_solve_suitesparse_b1(nrow,ncol,nz,irow,pcol,val,b,iopt_in)
    INTEGER, INTENT(in) :: nrow,ncol,nz
    INTEGER, DIMENSION(:), INTENT(in) :: irow,pcol
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: val
    REAL(kind=dp), DIMENSION(:), INTENT(inout) :: b
    INTEGER, INTENT(in) :: iopt_in

    INTEGER(kind=long) :: n
    INTEGER(kind=long), ALLOCATABLE, DIMENSION(:) :: Ai, Ap  !row-index Ai, column-pointer Ap
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: x !vector to store the solution

    ALLOCATE( x(SIZE(b)) )
    ALLOCATE( Ai(SIZE(irow)) )
    ALLOCATE( Ap(SIZE(pcol)) )

    IF (SIZE(pcol,1) .NE. ncol+1) THEN
       PRINT *, 'Wrong pcol'
       STOP
    END IF

    !   set default parameters
    CALL umf4def (control)

    n = nrow !convert from 1 to 0-based indexing
    Ai=irow-1 !convert from 1 to 0-based indexing
    Ap=pcol-1 !convert from 1 to 0-based indexing

    ! First, factorize the matrix. The factors are stored in *numeric* handle.
    IF (iopt_in .EQ. 0 .OR. iopt_in .EQ. 1) THEN
       !pre-order and symbolic analysis
       CALL umf4sym (n, n, Ap, Ai, val, symbolic, control, info_suitesparse)
       IF (sparse_talk) THEN
          IF (info_suitesparse(1) .EQ. 0) THEN
          ELSE
             PRINT *, 'Error occurred in umf4sym: ', info_suitesparse (1)
          ENDIF
       ENDIF

       CALL umf4num (Ap, Ai, val, symbolic, numeric, control, info_suitesparse)

       IF (sparse_talk) THEN
          IF (info_suitesparse(1) .EQ. 0) THEN
             PRINT *, 'Factorization succeeded'

          ELSE
             PRINT *, 'INFO from factorization = ', info_suitesparse(1)
          ENDIF
       END IF
    END IF

    ! Second, solve the system using the existing factors.
    IF (iopt_in .EQ. 0 .OR. iopt_in .EQ. 2) THEN
       IF ( sparse_solve_method .EQ. 2 ) THEN ! SuiteSparse (with (=2)
          CALL umf4solr (sys, Ap, Ai, val, x, b, numeric, control, info_suitesparse) !iterative refinement
       ELSE !or without (=3)) iterative refinement
          CALL umf4sol (sys, x, b, numeric, control, info_suitesparse) !without iterative refinement
       END IF
       b=x !store solution under b

       IF (sparse_talk) THEN
          IF (info_suitesparse(1) .EQ. 0) THEN
             PRINT *, 'Solve succeeded'
          ELSE
             PRINT *, 'INFO from triangular solve = ', info_suitesparse(1)
          ENDIF
       END IF
    END IF

    ! Last, free the storage allocated inside SuiteSparse
    IF (iopt_in .EQ. 0 .OR. iopt_in .EQ. 3) THEN
       CALL umf4fnum (numeric)
       CALL umf4fsym (symbolic)
    END IF

    IF (ALLOCATED(Ai)) DEALLOCATE(Ai)
    IF (ALLOCATED(Ap)) DEALLOCATE(Ap)
    IF (ALLOCATED(x))  DEALLOCATE(x)

  END SUBROUTINE sparse_solve_suitesparse_b1
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  ! solves A*x = b for sparse A and 1-D vector b
  ! A is specified through nrow,ncol,nz,irow,pcol,val
  ! results are returned in b
  ! Routines from SuperLU-Distribution
  SUBROUTINE sparse_solve_suitesparseComplex_b1(nrow,ncol,nz,irow,pcol,val,b,iopt_in)
    INTEGER, INTENT(in) :: nrow,ncol,nz
    INTEGER, DIMENSION(:), INTENT(in) :: irow,pcol
    COMPLEX(kind=dp), DIMENSION(:), INTENT(in) :: val
    COMPLEX(kind=dp), DIMENSION(:), INTENT(inout) :: b
    INTEGER, INTENT(in) :: iopt_in

    INTEGER :: k
    INTEGER(kind=long) :: n
    INTEGER(kind=long), ALLOCATABLE, DIMENSION(:) :: Ai, Ap  !row-index Ai, column-pointer Ap
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: xx,xz !vector to store the solution (real and imag. part)
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: valx, valz !val of matrix (real and imag. part)
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: bx, bz !rhs (real and imag part)

    ALLOCATE( xx(nrow) )
    ALLOCATE( xz(nrow) )
    ALLOCATE( bx(nrow) )
    ALLOCATE( bz(nrow) )
    ALLOCATE( valx(nz) )
    ALLOCATE( valz(nz) )


    bx=DBLE(b)
    bz=AIMAG(b)

    valx=DBLE(val)
    valz=AIMAG(val)

    ALLOCATE( Ai(SIZE(irow)) )
    ALLOCATE( Ap(SIZE(pcol)) )

    IF (SIZE(pcol,1) .NE. ncol+1) THEN
       PRINT *, 'Wrong pcol'
       STOP
    END IF

    !   set default parameters
    CALL umf4zdef (control)

    n = nrow
    Ai=irow-1 !convert from 1 to 0-based indexing
    Ap=pcol-1 !convert from 1 to 0-based indexing

    ! First, factorize the matrix. The factors are stored in *numeric* handle.
    IF (iopt_in .EQ. 0 .OR. iopt_in .EQ. 1) THEN
       !pre-order and symbolic analysis
       CALL umf4zsym (n, n, Ap, Ai, valx, valz, symbolic, control, info_suitesparse)

       IF (sparse_talk) THEN
          IF (info_suitesparse(1) .EQ. 0) THEN
             WRITE(*,80)  info_suitesparse (1), info_suitesparse (16), &
                  (info_suitesparse (21) * info_suitesparse (4)) / 2**20, &
                  (info_suitesparse (22) * info_suitesparse (4)) / 2**20, &
                  info_suitesparse (23), info_suitesparse (24), &
                  info_suitesparse (25)
80           FORMAT ('symbolic analysis:',/,&
                  '   status:  ', f5.0,/, &
                  '   time:    ', e10.4, ' (sec)',/, &
                  '   estimates (upper bound) for numeric LU:',/, &
                  '   size of LU:    ', f10.2, ' (MB)',/, &
                  '   memory needed: ', f10.2, ' (MB)',/, &
                  '   flop count:    ', e10.2,/, &
                  '   nnz (L):       ', f10.0,/, &
                  '   nnz (U):       ', f10.0)
          ELSE
             PRINT *, 'Error occurred in umf4sym: ', info_suitesparse (1)
          ENDIF
       ENDIF

       CALL umf4znum (Ap, Ai, valx, valz, symbolic, numeric, control, info_suitesparse)

       IF (sparse_talk) THEN
          IF (info_suitesparse(1) .EQ. 0) THEN
             PRINT *, 'Factorization succeeded'
             WRITE(*,90) info_suitesparse (1), info_suitesparse (66),&
                  (info_suitesparse (41) * info_suitesparse (4)) / 2**20, &
                  (info_suitesparse (42) * info_suitesparse (4)) / 2**20,&
                  info_suitesparse (43), info_suitesparse (44),&
                  info_suitesparse (45)
90           FORMAT ('numeric factorization:',/, &
                  '   status:  ', f5.0, /, &
                  '   time:    ', e10.4, /, &
                  '   actual numeric LU statistics:', /, &
                  '   size of LU:    ', f10.2, ' (MB)', /, &
                  '   memory needed: ', f10.2, ' (MB)', /, &
                  '   flop count:    ', e10.2, / &
                  '   nnz (L):       ', f10.0, / &
                  '   nnz (U):       ', f10.0)
          ELSE
             PRINT *, 'INFO from factorization = ', info_suitesparse(1)
          ENDIF
       END IF
    END IF

    ! Second, solve the system using the existing factors.
    IF (iopt_in .EQ. 0 .OR. iopt_in .EQ. 2) THEN
       IF ( sparse_solve_method .EQ. 2 ) THEN ! SuiteSparse (with (=2)
          CALL umf4zsolr (sys, Ap, Ai, valx, valz, xx, xz, bx, bz, numeric, &
               control, info_suitesparse) !iterative refinement
       ELSE !or without (=3)) iterative refinement
          CALL umf4zsol (sys, xx, xz, bx, bz, numeric, control, &
               info_suitesparse) !without iterative refinement
       END IF

       b=CMPLX(xx,xz, kind=kind(0d0)) !store solution under b

       IF (sparse_talk) THEN
          IF (info_suitesparse(1) .EQ. 0) THEN
             PRINT *, 'Solve succeeded'
          ELSE
             PRINT *, 'INFO from triangular solve = ', info_suitesparse(1)
          ENDIF
       END IF
    END IF

    ! Last, free the storage allocated inside SuiteSparse
    IF (iopt_in .EQ. 0 .OR. iopt_in .EQ. 3) THEN
       CALL umf4zfnum (numeric)
       CALL umf4zfsym (symbolic)
    END IF

    IF (ALLOCATED(Ai)) DEALLOCATE(Ai)
    IF (ALLOCATED(Ap)) DEALLOCATE(Ap)
    IF (ALLOCATED(xx))  DEALLOCATE(xx)
    IF (ALLOCATED(xz))  DEALLOCATE(xz)
    IF (ALLOCATED(bx))  DEALLOCATE(bx)
    IF (ALLOCATED(bz))  DEALLOCATE(bz)
    IF (ALLOCATED(valx))  DEALLOCATE(valx)
    IF (ALLOCATED(valz))  DEALLOCATE(valz)

  END SUBROUTINE sparse_solve_suitesparseComplex_b1
  !-------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------
  ! solves A*x = b for sparse A and 2-D array b
  ! A is specified through nrow,ncol,nz,irow,pcol,val
  ! results are returned in b
  ! Routines from SuiteSparse-Distribution
  SUBROUTINE sparse_solve_suitesparse_b2_loop(nrow,ncol,nz,irow,pcol,val,b,iopt_in)
    INTEGER, INTENT(in) :: nrow,ncol,nz
    INTEGER, DIMENSION(:), INTENT(in) :: irow,pcol
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: val
    REAL(kind=dp), DIMENSION(:,:), INTENT(inout) :: b
    INTEGER, INTENT(in) :: iopt_in

    INTEGER(kind=long) :: n, i
    INTEGER(kind=long), ALLOCATABLE, DIMENSION(:) :: Ai, Ap  !row-index Ai, column-pointer Ap
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: x !vector to store the solution
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: bloc

    !**********************************************************
    ! Patch from TU Graz ITPcp Plasma - 01.09.2015
    ! Wrong allocation size of x fixed
    !**********************************************************
    ALLOCATE( x(nrow) )
    ALLOCATE( Ai(SIZE(irow)) )
    ALLOCATE( Ap(SIZE(pcol)) )
    ALLOCATE(bloc(nrow))

    IF (SIZE(pcol,1) .NE. ncol+1) THEN
       PRINT *, 'Wrong pcol'
       STOP
    END IF

    !   set default parameters
    CALL umf4def (control)

    n = nrow
    bloc = 0.0_dp
    Ai=irow-1 !convert from 1 to 0-based indexing
    Ap=pcol-1 !convert from 1 to 0-based indexing

    IF (SIZE(pcol,1) .NE. ncol+1) THEN
       PRINT *, 'Wrong pcol'
       STOP
    END IF


    ! First, factorize the matrix. The factors are stored in *numeric* handle.
    IF (iopt_in .EQ. 0 .OR. iopt_in .EQ. 1) THEN
       !pre-order and symbolic analysis
       CALL umf4sym (n, n, Ap, Ai, val, symbolic, control, info_suitesparse)
       IF (sparse_talk) THEN
          IF (info_suitesparse(1) .EQ. 0) THEN
          ELSE
             PRINT *, 'Error occurred in umf4sym: ', info_suitesparse (1)
          ENDIF
       ENDIF

       CALL umf4num (Ap, Ai, val, symbolic, numeric, control, info_suitesparse)

       IF (sparse_talk) THEN
          IF (info_suitesparse(1) .EQ. 0) THEN
             PRINT *, 'Factorization succeeded'
          ELSE
             PRINT *, 'INFO from factorization = ', info_suitesparse(1)
          ENDIF
       END IF
    END IF

    ! Second, solve the system using the existing factors.
    IF (iopt_in .EQ. 0 .OR. iopt_in .EQ. 2) THEN
       DO i = 1,SIZE(b,2)
          bloc = b(:,i)
          IF ( sparse_solve_method .EQ. 2 ) THEN ! SuiteSparse (with (=2)
             CALL umf4solr (sys, Ap, Ai, val, x, bloc, numeric, control, info_suitesparse) !iterative refinement
          ELSE !or without (=3)) iterative refinement
             CALL umf4sol (sys, x, bloc, numeric, control, info_suitesparse) !without iterative refinement
          END IF

          IF (sparse_talk) THEN
             IF (info_suitesparse(1) .EQ. 0) THEN
             ELSE
                PRINT *, 'INFO from solve = ', info_suitesparse(1)
             ENDIF
          END IF
          b(:,i) = x
       END DO
    END IF

    ! Last, free the storage allocated inside SuiteSparse
    IF (iopt_in .EQ. 0 .OR. iopt_in .EQ. 3) THEN
       CALL umf4fnum (numeric)
       CALL umf4fsym (symbolic)
    END IF

    IF (ALLOCATED(bloc)) DEALLOCATE(bloc)
    IF (ALLOCATED(Ai)) DEALLOCATE(Ai)
    IF (ALLOCATED(Ap)) DEALLOCATE(Ap)
    IF (ALLOCATED(x))  DEALLOCATE(x)

    RETURN
  END SUBROUTINE sparse_solve_suitesparse_b2_loop
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  ! solves A*x = b for sparse A and 2-D array b
  ! A is specified through nrow,ncol,nz,irow,pcol,val
  ! results are returned in b
  ! Routines from SuiteSparse-Distribution
  SUBROUTINE sparse_solve_suitesparseComplex_b2_loop(nrow,ncol,nz,irow,pcol,val,b,iopt_in)
    INTEGER, INTENT(in) :: nrow,ncol,nz
    INTEGER, DIMENSION(:), INTENT(in) :: irow,pcol
    COMPLEX(kind=dp), DIMENSION(:), INTENT(in) :: val
    COMPLEX(kind=dp), DIMENSION(:,:), INTENT(inout) :: b
    INTEGER, INTENT(in) :: iopt_in

    INTEGER(kind=long) :: n, i
    INTEGER(kind=long), ALLOCATABLE, DIMENSION(:) :: Ai, Ap  !row-index Ai, column-pointer Ap
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: xx,xz !vector to store the solution (real and imag. part)
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: valx, valz !val of matrix (real and imag. part)
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:) :: bx, bz !rhs (real and imag part)
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: blocx, blocz

    ALLOCATE( xx(nrow) )
    ALLOCATE( xz(nrow) )
    ALLOCATE( bx(nrow, SIZE(b,2)) )
    ALLOCATE( bz(nrow, SIZE(b,2)) )
    ALLOCATE( valx(nz) )
    ALLOCATE( valz(nz) )

    bx=DBLE(b)
    bz=AIMAG(b)
    valx=DBLE(val)
    valz=AIMAG(val)

    ALLOCATE( Ai(SIZE(irow)) )
    ALLOCATE( Ap(SIZE(pcol)) )
    ALLOCATE(blocx(nrow))
    ALLOCATE(blocz(nrow))

    n = nrow
    blocx = 0.0_dp
    blocz = 0.0_dp
    Ai=irow-1 !convert from 1 to 0-based indexing
    Ap=pcol-1 !convert from 1 to 0-based indexing

    IF (SIZE(pcol,1) .NE. ncol+1) THEN
       PRINT *, 'Wrong pcol'
       STOP
    END IF

    !   set default parameters
    CALL umf4zdef (control)


    ! First, factorize the matrix. The factors are stored in *numeric* handle.
    IF (iopt_in .EQ. 0 .OR. iopt_in .EQ. 1) THEN
       !pre-order and symbolic analysis
       CALL umf4zsym (n, n, Ap, Ai, valx, valz, symbolic, control, info_suitesparse)

       IF (sparse_talk) THEN
          IF (info_suitesparse(1) .EQ. 0) THEN
             WRITE(*,80)  info_suitesparse (1), info_suitesparse (16), &
                  (info_suitesparse (21) * info_suitesparse (4)) / 2**20, &
                  (info_suitesparse (22) * info_suitesparse (4)) / 2**20, &
                  info_suitesparse (23), info_suitesparse (24), &
                  info_suitesparse (25)
80           FORMAT ('symbolic analysis:',/,&
                  '   status:  ', f5.0,/, &
                  '   time:    ', e10.4, ' (sec)',/, &
                  '   estimates (upper bound) for numeric LU:',/, &
                  '   size of LU:    ', f10.2, ' (MB)',/, &
                  '   memory needed: ', f10.2, ' (MB)',/, &
                  '   flop count:    ', e10.2,/, &
                  '   nnz (L):       ', f10.0,/, &
                  '   nnz (U):       ', f10.0)

          ELSE
             PRINT *, 'Error occurred in umf4sym: ', info_suitesparse (1)
          ENDIF
       ENDIF

       CALL umf4znum (Ap, Ai, valx, valz, symbolic, numeric, control, info_suitesparse)

       IF (sparse_talk) THEN
          IF (info_suitesparse(1) .EQ. 0) THEN
             PRINT *, 'Factorization succeeded'
             WRITE(*,90) info_suitesparse (1), info_suitesparse (66),&
                  (info_suitesparse (41) * info_suitesparse (4)) / 2**20, &
                  (info_suitesparse (42) * info_suitesparse (4)) / 2**20,&
                  info_suitesparse (43), info_suitesparse (44),&
                  info_suitesparse (45)
90           FORMAT ('numeric factorization:',/, &
                  '   status:  ', f5.0, /, &
                  '   time:    ', e10.4, /, &
                  '   actual numeric LU statistics:', /, &
                  '   size of LU:    ', f10.2, ' (MB)', /, &
                  '   memory needed: ', f10.2, ' (MB)', /, &
                  '   flop count:    ', e10.2, / &
                  '   nnz (L):       ', f10.0, / &
                  '   nnz (U):       ', f10.0)
          ELSE
             PRINT *, 'INFO from factorization = ', info_suitesparse(1)
          ENDIF
       END IF
    END IF

    ! Second, solve the system using the existing factors.
    IF (iopt_in .EQ. 0 .OR. iopt_in .EQ. 2) THEN
       DO i = 1,SIZE(b,2)
          blocx = bx(:,i)
          blocz = bz(:,i)
          IF ( sparse_solve_method .EQ. 2 ) THEN ! SuiteSparse (with (=2)
             CALL umf4zsolr (sys, Ap, Ai, valx, valz, xx, xz, blocx, blocz, numeric,&
                  control, info_suitesparse) !iterative refinement
          ELSE !or without (=3)) iterative refinement
             CALL umf4zsol (sys, xx, xz, blocx, blocz, numeric,&
                  control, info_suitesparse) !without iterative refinement
          END IF

          IF (sparse_talk) THEN
             IF (info_suitesparse(1) .EQ. 0) THEN
             ELSE
                PRINT *, 'INFO from solve = ', info_suitesparse(1)
             ENDIF
          END IF
          b(:,i)=CMPLX(xx,xz, kind=kind(0d0))
       END DO
    END IF

    ! Last, free the storage allocated inside SuiteSparse
    IF (iopt_in .EQ. 0 .OR. iopt_in .EQ. 3) THEN
       CALL umf4zfnum (numeric)
       CALL umf4zfsym (symbolic)
    END IF

    IF (ALLOCATED(blocx)) DEALLOCATE(blocx)
    IF (ALLOCATED(blocz)) DEALLOCATE(blocz)
    IF (ALLOCATED(Ai)) DEALLOCATE(Ai)
    IF (ALLOCATED(Ap)) DEALLOCATE(Ap)
    IF (ALLOCATED(xx))  DEALLOCATE(xx)
    IF (ALLOCATED(xz))  DEALLOCATE(xz)
    IF (ALLOCATED(bx))  DEALLOCATE(bx)
    IF (ALLOCATED(bz))  DEALLOCATE(bz)
    IF (ALLOCATED(valx))  DEALLOCATE(valx)
    IF (ALLOCATED(valz))  DEALLOCATE(valz)

  END SUBROUTINE sparse_solve_suitesparseComplex_b2_loop
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  ! column pointer pcol to full column index icol
  SUBROUTINE col_pointer2full(pcol,icol)

    INTEGER, DIMENSION(:), INTENT(in) :: pcol
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(inout) :: icol

    INTEGER :: nz
    INTEGER :: nc_old,c,nc,ncol

    ncol = SIZE(pcol,1)-1
    nz = pcol(ncol+1) - 1
    IF (ALLOCATED(icol)) DEALLOCATE(icol)
    ALLOCATE(icol(nz))
    nc_old = 0
    DO c = 1,ncol
       nc = pcol(c+1) - pcol(c)
       icol(nc_old+1:nc_old+nc) = c;
       nc_old = nc_old + nc;
    END DO

  END SUBROUTINE col_pointer2full
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  ! full column index icol to column pointer pcol
  SUBROUTINE col_full2pointer(icol,pcol)

    INTEGER, DIMENSION(:), INTENT(in) :: icol
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(inout) :: pcol

    INTEGER :: ncol,nz
    INTEGER :: c_c,c_old,k,c,kc

    ncol = MAXVAL(icol)
    nz = SIZE(icol,1)

    IF (ALLOCATED(pcol)) DEALLOCATE(pcol)
    ALLOCATE(pcol(ncol+1))

    c_c = 1
    pcol(c_c) = 1
    c_old = 0
    DO k = 1,nz
       c = icol(k)
       IF (c .NE. c_old) THEN
          IF (c .GT. c_old + 1) THEN
             DO kc = c_old+1,c
                c_c = c_c + 1
                pcol(c_c) = k
             END DO
          ELSE
             c_c = c_c + 1
             pcol(c_c) = k+1
          END IF
          c_old = c
       ELSE
          pcol(c_c) = k+1;
       END IF
    END DO
    IF (c_c .LT. ncol+1) pcol(c_c+1:ncol+1) = pcol(c_c)

  END SUBROUTINE col_full2pointer
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  ! sparse to full conversion
  SUBROUTINE sp2full(irow,pcol,val,nrow,ncol,A)

    INTEGER, DIMENSION(:), INTENT(in) :: irow,pcol
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: val
    INTEGER, INTENT(in) :: nrow,ncol
    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: A

    INTEGER :: nz,n,ir,ic
    INTEGER, DIMENSION(:), ALLOCATABLE :: icol

    nz = SIZE(val,1)
    IF (SIZE(pcol,1) .NE. nz) THEN
       IF (SIZE(pcol,1) .EQ. ncol+1) THEN
          CALL column_pointer2full(pcol,icol)
       ELSE
          PRINT *, 'Error in sparse2full: icol is not correct'
          STOP
       END IF
    ELSE
       ALLOCATE(icol(SIZE(pcol)))
       icol = pcol
    END IF

    IF (ALLOCATED(A)) DEALLOCATE(A)
    ALLOCATE(A(nrow,ncol))
    A = 0.0_dp
    DO n = 1,nz
       ir = irow(n)
       ic = icol(n)
       A(ir,ic) = A(ir,ic) + val(n)
    END DO
    IF (ALLOCATED(icol)) DEALLOCATE(icol)

  END SUBROUTINE sp2full
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  ! sparse to full conversion for complex matrices
  SUBROUTINE sp2fullComplex(irow,pcol,val,nrow,ncol,A)

    INTEGER, DIMENSION(:), INTENT(in) :: irow,pcol
    COMPLEX(kind=dp), DIMENSION(:), INTENT(in) :: val
    INTEGER, INTENT(in) :: nrow,ncol
    COMPLEX(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: A

    INTEGER :: nz,n,ir,ic
    INTEGER, DIMENSION(:), ALLOCATABLE :: icol

    nz = SIZE(val,1)
    IF (SIZE(pcol,1) .NE. nz) THEN
       IF (SIZE(pcol,1) .EQ. ncol+1) THEN
          CALL column_pointer2full(pcol,icol)
       ELSE
          PRINT *, 'Error in sparse2full: icol is not correct'
          STOP
       END IF
    ELSE
       ALLOCATE(icol(SIZE(pcol)))
       icol = pcol
    END IF

    IF (ALLOCATED(A)) DEALLOCATE(A)
    ALLOCATE(A(nrow,ncol))
    A = 0.0_dp
    DO n = 1,nz
       ir = irow(n)
       ic = icol(n)
       A(ir,ic) = A(ir,ic) + val(n)
    END DO
    IF (ALLOCATED(icol)) DEALLOCATE(icol)

  END SUBROUTINE sp2fullComplex
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  !> \brief Full to sparse conversion for real matrix.
  !>
  !> input:
  !> A: matrix (intent in) to convert to a sparse format. Note that this
  !>   is formally the only input of this subroutine.
  !> irow: allocatable vector (intent inout), where the row indices will
  !>   be stored.
  !>   Note that the original content is not read, the intent
  !>   is because the array is deallocated if it is already allocated.
  !> pcol: allocatable vector (intent inout), stores ?
  !>   Note that the original content is not read, the intent
  !>   is because the array is deallocated if it is already allocated.
  !> values: allocatable vector (intent inout), stores the nonzero
  !>   values of the matrix.
  !>   Note that the original content is not read, the intent
  !>   is because the array is deallocated if it is already allocated.
  !> nrow, ncol: integers (intent out), giving the size of the original
  !>   matrix.
  !> nz_out: integer (intent out, optional), giving the number of
  !>   nonzero elements of matrix A, i.e. the size of the vectors irow,
  !>   pcol(?) and values on output.
  SUBROUTINE full2sp(A,irow,pcol,values,nrow,ncol,nz_out)

    REAL(kind=dp), DIMENSION(:,:), INTENT(in) :: A
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(inout) :: irow,pcol
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(inout) :: values
    INTEGER, INTENT(out) :: nrow,ncol
    INTEGER, OPTIONAL, INTENT(out) :: nz_out

    INTEGER, DIMENSION(:), ALLOCATABLE :: icol
    INTEGER :: nz,nc,nr,n

    nrow = SIZE(A,1)
    ncol = SIZE(A,2)

    nz = 0
    DO nc = 1,ncol
       DO nr = 1,nrow
          IF (A(nr,nc) .NE. 0.0_dp) nz = nz + 1
       END DO
    END DO

    IF (ALLOCATED(irow)) DEALLOCATE(irow)
    ALLOCATE(irow(nz))
    IF (ALLOCATED(icol)) DEALLOCATE(icol)
    ALLOCATE(icol(nz))
    IF (ALLOCATED(values)) DEALLOCATE(values)
    ALLOCATE(values(nz))

    n = 0
    DO nc = 1,ncol
       DO nr = 1,nrow
          IF (A(nr,nc) .NE. 0.0_dp) THEN
             n = n + 1
             irow(n) = nr
             icol(n) = nc
             values(n) = A(nr,nc)
          END IF
       END DO
    END DO

    CALL column_full2pointer(icol,pcol)

    IF (PRESENT(nz_out)) nz_out = nz
    IF (ALLOCATED(icol)) DEALLOCATE(icol)

  END SUBROUTINE full2sp
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  ! full to sparse conversion for complex matrices
  SUBROUTINE full2spComplex(A,irow,pcol,val,nrow,ncol,nz_out)

    COMPLEX(kind=dp), DIMENSION(:,:), INTENT(in) :: A
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(inout) :: irow,pcol
    COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(inout) :: val
    INTEGER, INTENT(out) :: nrow,ncol
    INTEGER, OPTIONAL, INTENT(out) :: nz_out

    INTEGER, DIMENSION(:), ALLOCATABLE :: icol
    INTEGER :: nz,nc,nr,n

    nrow = SIZE(A,1)
    ncol = SIZE(A,2)

    nz = 0
    DO nc = 1,ncol
       DO nr = 1,nrow
          IF (A(nr,nc) .NE. 0.0_dp) nz = nz + 1
       END DO
    END DO

    IF (ALLOCATED(irow)) DEALLOCATE(irow)
    ALLOCATE(irow(nz))
    IF (ALLOCATED(icol)) DEALLOCATE(icol)
    ALLOCATE(icol(nz))
    IF (ALLOCATED(val)) DEALLOCATE(val)
    ALLOCATE(val(nz))

    n = 0
    DO nc = 1,ncol
       DO nr = 1,nrow
          IF (A(nr,nc) .NE. 0.0_dp) THEN
             n = n + 1
             irow(n) = nr
             icol(n) = nc
             val(n)  = A(nr,nc)
          END IF
       END DO
    END DO

    CALL column_full2pointer(icol,pcol)

    IF (PRESENT(nz_out)) nz_out = nz
    IF (ALLOCATED(icol)) DEALLOCATE(icol)

  END SUBROUTINE full2spComplex
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  ! computes A*x for sparse A and 1-D array x
  ! A is specified through nrow,ncol,nz,irow,pcol,val
  ! results are returned in r
  SUBROUTINE sp_matmul_b1(nrow,ncol,irow,pcol,val,x,r)
    INTEGER, INTENT(in) :: nrow,ncol
    INTEGER, DIMENSION(:), INTENT(in) :: irow,pcol
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: val
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: x
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(inout) :: r

    INTEGER :: nz,n,ic,ir
    INTEGER, DIMENSION(:), ALLOCATABLE :: icol



    nz = SIZE(val,1)
    IF (SIZE(pcol,1) .NE. nz) THEN
       IF (SIZE(pcol,1) .EQ. ncol+1) THEN
          CALL column_pointer2full(pcol,icol)
       ELSE
          PRINT *, 'Error in sparse_matmul: icol is not correct'
          STOP
       END IF
    ELSE
       ALLOCATE(icol(SIZE(pcol)))
       icol = pcol
    END IF

    IF (ncol .NE. SIZE(x,1)) THEN
       PRINT *, 'Error in sparse_matmul: sizes do not fit'
       STOP
    END IF
    IF (ALLOCATED(r)) DEALLOCATE(r)

    ALLOCATE(r(nrow))
    r = 0.0_dp

    DO n = 1,nz
       ic = icol(n)
       ir = irow(n)
       r(ir) = r(ir) + val(n)*x(ic)
    END DO
    IF (ALLOCATED(icol)) DEALLOCATE(icol)

  END SUBROUTINE sp_matmul_b1
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  ! computes A*x for sparse A and 1-D array x
  ! A is specified through nrow,ncol,nz,irow,pcol,val
  ! results are returned in r
  SUBROUTINE sp_matmulComplex_b1(nrow,ncol,irow,pcol,val,x,r)
    INTEGER, INTENT(in) :: nrow,ncol
    INTEGER, DIMENSION(:), INTENT(in) :: irow,pcol
    COMPLEX(kind=dp), DIMENSION(:), INTENT(in) :: val
    COMPLEX(kind=dp), DIMENSION(:), INTENT(in) :: x
    COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(inout) :: r

    INTEGER :: nz,n,ic,ir
    INTEGER, DIMENSION(:), ALLOCATABLE :: icol



    nz = SIZE(val,1)
    IF (SIZE(pcol,1) .NE. nz) THEN
       IF (SIZE(pcol,1) .EQ. ncol+1) THEN
          CALL column_pointer2full(pcol,icol)
       ELSE
          PRINT *, 'Error in sparse_matmul: icol is not correct'
          STOP
       END IF
    ELSE
       ALLOCATE(icol(SIZE(pcol)))
       icol = pcol
    END IF

    IF (ncol .NE. SIZE(x,1)) THEN
       PRINT *, 'Error in sparse_matmul: sizes do not fit'
       STOP
    END IF
    IF (ALLOCATED(r)) DEALLOCATE(r)

    ALLOCATE(r(nrow))
    r = 0.0_dp

    DO n = 1,nz
       ic = icol(n)
       ir = irow(n)
       r(ir) = r(ir) + val(n)*x(ic)
    END DO
    IF (ALLOCATED(icol)) DEALLOCATE(icol)

  END SUBROUTINE sp_matmulComplex_b1
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  ! computes A*x for sparse A and 2-D array x
  ! A is specified through nrow,ncol,nz,irow,pcol,val
  ! results are returned in r
  SUBROUTINE sp_matmul_b2(nrow,ncol,irow,pcol,val,x,r)
    INTEGER, INTENT(in) :: nrow,ncol
    INTEGER, DIMENSION(:), INTENT(in) :: irow,pcol
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: val
    REAL(kind=dp), DIMENSION(:,:), INTENT(in) :: x
    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: r

    INTEGER :: nz,n,ic,ir
    INTEGER, DIMENSION(:), ALLOCATABLE :: icol

    nz = SIZE(val,1)
    IF (SIZE(pcol,1) .NE. nz) THEN
       IF (SIZE(pcol,1) .EQ. ncol+1) THEN
          CALL column_pointer2full(pcol,icol)
       ELSE
          PRINT *, 'Error in sparse_matmul: icol is not correct'
          STOP
       END IF
    ELSE
       ALLOCATE(icol(SIZE(pcol)))
       icol = pcol
    END IF

    IF (ncol .NE. SIZE(x,1)) THEN
       PRINT *, 'Error in sparse_matmul: sizes do not fit'
       STOP
    END IF
    IF (ALLOCATED(r)) DEALLOCATE(r)
    !ALLOCATE(r(SIZE(x,1),SIZE(x,2)))
    ALLOCATE(r(nrow,SIZE(x,2)))
    r = 0.0_dp

    DO n = 1,nz
       ic = icol(n)
       ir = irow(n)
       r(ir,:) = r(ir,:) + val(n)*x(ic,:)
    END DO
    IF (ALLOCATED(icol)) DEALLOCATE(icol)

  END SUBROUTINE sp_matmul_b2
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  ! computes A*x for sparse A and 2-D array x
  ! A is specified through nrow,ncol,nz,irow,pcol,val
  ! results are returned in r
  SUBROUTINE sp_matmulComplex_b2(nrow,ncol,irow,pcol,val,x,r)
    INTEGER, INTENT(in) :: nrow,ncol
    INTEGER, DIMENSION(:), INTENT(in) :: irow,pcol
    COMPLEX(kind=dp), DIMENSION(:), INTENT(in) :: val
    COMPLEX(kind=dp), DIMENSION(:,:), INTENT(in) :: x
    COMPLEX(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: r

    INTEGER :: nz,n,ic,ir
    INTEGER, DIMENSION(:), ALLOCATABLE :: icol

    nz = SIZE(val,1)
    IF (SIZE(pcol,1) .NE. nz) THEN
       IF (SIZE(pcol,1) .EQ. ncol+1) THEN
          CALL column_pointer2full(pcol,icol)
       ELSE
          PRINT *, 'Error in sparse_matmul: icol is not correct'
          STOP
       END IF
    ELSE
       ALLOCATE(icol(SIZE(pcol)))
       icol = pcol
    END IF

    IF (ncol .NE. SIZE(x,1)) THEN
       PRINT *, 'Error in sparse_matmul: sizes do not fit'
       STOP
    END IF
    IF (ALLOCATED(r)) DEALLOCATE(r)

    ALLOCATE(r(nrow,SIZE(x,2)))
    r = 0.0_dp

    DO n = 1,nz
       ic = icol(n)
       ir = irow(n)
       r(ir,:) = r(ir,:) + val(n)*x(ic,:)
    END DO
    IF (ALLOCATED(icol)) DEALLOCATE(icol)

  END SUBROUTINE sp_matmulComplex_b2
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  ! computes A*x for sparse A and 1-D array x
  ! results are returned in r
  SUBROUTINE sp_matmul_A_b1(A,x,r)
    REAL(kind=dp), DIMENSION(:,:), INTENT(in) :: A
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: x
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(inout) :: r

    INTEGER :: nrow,ncol
    INTEGER, DIMENSION(:), ALLOCATABLE :: irow,pcol
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: val

    CALL full2sparse(A,irow,pcol,val,nrow,ncol)
    CALL sparse_matmul(nrow,ncol,irow,pcol,val,x,r)

    IF (ALLOCATED(irow)) DEALLOCATE(irow)
    IF (ALLOCATED(pcol)) DEALLOCATE(pcol)
    IF (ALLOCATED(val))  DEALLOCATE(val)

  END SUBROUTINE sp_matmul_A_b1
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  ! computes A*x for sparse A and 1-D array x
  ! results are returned in r
  SUBROUTINE sp_matmulComplex_A_b1(A,x,r)
    COMPLEX(kind=dp), DIMENSION(:,:), INTENT(in) :: A
    COMPLEX(kind=dp), DIMENSION(:), INTENT(in) :: x
    COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(inout) :: r

    INTEGER :: nrow,ncol
    INTEGER, DIMENSION(:), ALLOCATABLE :: irow,pcol
    COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE :: val

    CALL full2sparse(A,irow,pcol,val,nrow,ncol)
    CALL sparse_matmul(nrow,ncol,irow,pcol,val,x,r)

    IF (ALLOCATED(irow)) DEALLOCATE(irow)
    IF (ALLOCATED(pcol)) DEALLOCATE(pcol)
    IF (ALLOCATED(val))  DEALLOCATE(val)

  END SUBROUTINE sp_matmulComplex_A_b1
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  ! computes A*x for sparse A and 2-D array x
  ! results are returned in r
  SUBROUTINE sp_matmul_A_b2(A,x,r)
    REAL(kind=dp), DIMENSION(:,:), INTENT(in) :: A
    REAL(kind=dp), DIMENSION(:,:), INTENT(in) :: x
    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: r

    INTEGER :: nrow,ncol
    INTEGER, DIMENSION(:), ALLOCATABLE :: irow,pcol
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: val

    CALL full2sparse(A,irow,pcol,val,nrow,ncol)
    CALL sparse_matmul(nrow,ncol,irow,pcol,val,x,r)

    IF (ALLOCATED(irow)) DEALLOCATE(irow)
    IF (ALLOCATED(pcol)) DEALLOCATE(pcol)
    IF (ALLOCATED(val))  DEALLOCATE(val)

  END SUBROUTINE sp_matmul_A_b2
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  ! computes A*x for sparse A and 2-D array x
  ! results are returned in r
  SUBROUTINE sp_matmulComplex_A_b2(A,x,r)
    COMPLEX(kind=dp), DIMENSION(:,:), INTENT(in) :: A
    COMPLEX(kind=dp), DIMENSION(:,:), INTENT(in) :: x
    COMPLEX(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: r

    INTEGER :: nrow,ncol
    INTEGER, DIMENSION(:), ALLOCATABLE :: irow,pcol
    COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE :: val

    CALL full2sparse(A,irow,pcol,val,nrow,ncol)
    CALL sparse_matmul(nrow,ncol,irow,pcol,val,x,r)

    IF (ALLOCATED(irow)) DEALLOCATE(irow)
    IF (ALLOCATED(pcol)) DEALLOCATE(pcol)
    IF (ALLOCATED(val))  DEALLOCATE(val)

  END SUBROUTINE sp_matmulComplex_A_b2
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  ! tests A*x-b and returns errors
  SUBROUTINE sp_test_b1(nrow,ncol,irow,pcol,val,x,b,max_abs_err_out,max_rel_err_out)
    INTEGER, INTENT(in) :: nrow,ncol
    INTEGER, DIMENSION(:), INTENT(in) :: irow,pcol
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: val
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: x
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: b
    REAL(kind=dp), OPTIONAL, INTENT(out) :: max_abs_err_out,max_rel_err_out

    REAL(kind=dp) :: max_abs_err,max_rel_err
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: r

    CALL sparse_matmul(nrow,ncol,irow,pcol,val,x,r)
    r = ABS(r - b)
    max_abs_err = MAXVAL(r)
    max_rel_err = max_abs_err / MAXVAL(ABS(b))
    IF (sparse_talk) THEN
       PRINT *, 'max_abs_err=',max_abs_err,' max_rel_err=',max_rel_err
    END IF

    IF (PRESENT(max_abs_err_out)) max_abs_err_out = max_abs_err
    IF (PRESENT(max_rel_err_out)) max_rel_err_out = max_rel_err
    IF (ALLOCATED(r)) DEALLOCATE(r)

  END SUBROUTINE sp_test_b1
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  ! tests A*x-b and returns errors
  SUBROUTINE sp_testComplex_b1(nrow,ncol,irow,pcol,val,x,b,max_abs_err_out,max_rel_err_out)
    INTEGER, INTENT(in) :: nrow,ncol
    INTEGER, DIMENSION(:), INTENT(in) :: irow,pcol
    COMPLEX(kind=dp), DIMENSION(:), INTENT(in) :: val
    COMPLEX(kind=dp), DIMENSION(:), INTENT(in) :: x
    COMPLEX(kind=dp), DIMENSION(:), INTENT(in) :: b
    REAL(kind=dp), OPTIONAL, INTENT(out) :: max_abs_err_out,max_rel_err_out

    REAL(kind=dp) :: max_abs_err,max_rel_err
    COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE :: r

    CALL sparse_matmul(nrow,ncol,irow,pcol,val,x,r)
    r = ABS(r - b)
    max_abs_err = MAXVAL(SQRT(REAL(r)**2+AIMAG(r)**2))
    max_rel_err = max_abs_err / MAXVAL(SQRT(REAL(b)**2+AIMAG(b)**2))
    IF (sparse_talk) THEN
       PRINT *, 'max_abs_err=',max_abs_err,' max_rel_err=',max_rel_err
    END IF

    IF (PRESENT(max_abs_err_out)) max_abs_err_out = max_abs_err
    IF (PRESENT(max_rel_err_out)) max_rel_err_out = max_rel_err
    IF (ALLOCATED(r)) DEALLOCATE(r)

  END SUBROUTINE sp_testComplex_b1
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  ! tests A*x-b and returns errors
  SUBROUTINE sp_test_b2(nrow,ncol,irow,pcol,val,x,b,max_abs_err_out,max_rel_err_out)
    INTEGER, INTENT(in) :: nrow,ncol
    INTEGER, DIMENSION(:), INTENT(in) :: irow,pcol
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: val
    REAL(kind=dp), DIMENSION(:,:), INTENT(in) :: x
    REAL(kind=dp), DIMENSION(:,:), INTENT(in) :: b
    REAL(kind=dp), OPTIONAL, INTENT(out) :: max_abs_err_out,max_rel_err_out

    REAL(kind=dp) :: max_abs_err,max_rel_err
    REAL(kind=dp) :: abs_err,rel_err
    INTEGER :: ic

    max_abs_err = 0.0_dp
    max_rel_err = 0.0_dp

    DO ic = 1,SIZE(x,2)
       CALL sparse_solver_test(nrow,ncol,irow,pcol,val,x(:,ic),b(:,ic),abs_err,rel_err)
       max_abs_err = MAX(max_abs_err,abs_err)
       max_rel_err = MAX(max_rel_err,rel_err)
    END DO
    IF (sparse_talk) THEN
       PRINT *, 'max_abs_err=',max_abs_err,' max_rel_err=',max_rel_err,' total'
    END IF

    IF (PRESENT(max_abs_err_out)) max_abs_err_out = max_abs_err
    IF (PRESENT(max_rel_err_out)) max_rel_err_out = max_rel_err

  END SUBROUTINE sp_test_b2
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  ! tests A*x-b and returns errors
  SUBROUTINE sp_testComplex_b2(nrow,ncol,irow,pcol,val,x,b,max_abs_err_out,max_rel_err_out)
    INTEGER, INTENT(in) :: nrow,ncol
    INTEGER, DIMENSION(:), INTENT(in) :: irow,pcol
    COMPLEX(kind=dp), DIMENSION(:), INTENT(in) :: val
    COMPLEX(kind=dp), DIMENSION(:,:), INTENT(in) :: x
    COMPLEX(kind=dp), DIMENSION(:,:), INTENT(in) :: b
    REAL(kind=dp), OPTIONAL, INTENT(out) :: max_abs_err_out,max_rel_err_out

    REAL(kind=dp) :: max_abs_err,max_rel_err
    REAL(kind=dp) :: abs_err,rel_err
    INTEGER :: ic

    max_abs_err = 0.0_dp
    max_rel_err = 0.0_dp

    DO ic = 1,SIZE(x,2)
       CALL sparse_solver_test(nrow,ncol,irow,pcol,val,x(:,ic),b(:,ic),abs_err,rel_err)
       max_abs_err = MAX(max_abs_err,abs_err)
       max_rel_err = MAX(max_rel_err,rel_err)
    END DO
    IF (sparse_talk) THEN
       PRINT *, 'max_abs_err=',max_abs_err,' max_rel_err=',max_rel_err,' total'
    END IF

    IF (PRESENT(max_abs_err_out)) max_abs_err_out = max_abs_err
    IF (PRESENT(max_rel_err_out)) max_rel_err_out = max_rel_err

  END SUBROUTINE sp_testComplex_b2
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  ! tests A*x-b and returns errors
  SUBROUTINE sp_test_A_b1(A,x,b,max_abs_err_out,max_rel_err_out)
    REAL(kind=dp), DIMENSION(:,:), INTENT(in) :: A
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: x
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(in) :: b
    REAL(kind=dp), OPTIONAL, INTENT(out) :: max_abs_err_out,max_rel_err_out

    REAL(kind=dp) :: max_abs_err,max_rel_err
    INTEGER :: nrow,ncol
    INTEGER, DIMENSION(:), ALLOCATABLE :: irow,pcol
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: val

    CALL full2sparse(A,irow,pcol,val,nrow,ncol)
    CALL sparse_solver_test(nrow,ncol,irow,pcol,val,x,b,max_abs_err,max_rel_err)

    IF (PRESENT(max_abs_err_out)) max_abs_err_out = max_abs_err
    IF (PRESENT(max_rel_err_out)) max_rel_err_out = max_rel_err
    IF (ALLOCATED(irow)) DEALLOCATE(irow)
    IF (ALLOCATED(pcol)) DEALLOCATE(pcol)
    IF (ALLOCATED(val))  DEALLOCATE(val)

  END SUBROUTINE sp_test_A_b1
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  ! tests A*x-b and returns errors
  SUBROUTINE sp_testComplex_A_b1(A,x,b,max_abs_err_out,max_rel_err_out)
    COMPLEX(kind=dp), DIMENSION(:,:), INTENT(in) :: A
    COMPLEX(kind=dp), DIMENSION(:), INTENT(in) :: x
    COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(in) :: b
    REAL(kind=dp), OPTIONAL, INTENT(out) :: max_abs_err_out,max_rel_err_out

    REAL(kind=dp) :: max_abs_err,max_rel_err
    INTEGER :: nrow,ncol
    INTEGER, DIMENSION(:), ALLOCATABLE :: irow,pcol
    COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE :: val

    CALL full2sparse(A,irow,pcol,val,nrow,ncol)
    CALL sparse_solver_test(nrow,ncol,irow,pcol,val,x,b,max_abs_err,max_rel_err)

    IF (PRESENT(max_abs_err_out)) max_abs_err_out = max_abs_err
    IF (PRESENT(max_rel_err_out)) max_rel_err_out = max_rel_err
    IF (ALLOCATED(irow)) DEALLOCATE(irow)
    IF (ALLOCATED(pcol)) DEALLOCATE(pcol)
    IF (ALLOCATED(val))  DEALLOCATE(val)

  END SUBROUTINE sp_testComplex_A_b1
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  ! tests A*x-b and returns errors
  SUBROUTINE sp_test_A_b2(A,x,b,max_abs_err_out,max_rel_err_out)
    REAL(kind=dp), DIMENSION(:,:), INTENT(in) :: A
    REAL(kind=dp), DIMENSION(:,:), INTENT(in) :: x
    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(in) :: b
    REAL(kind=dp), OPTIONAL, INTENT(out) :: max_abs_err_out,max_rel_err_out

    REAL(kind=dp) :: max_abs_err,max_rel_err
    INTEGER :: nrow,ncol
    INTEGER, DIMENSION(:), ALLOCATABLE :: irow,pcol
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: val

    CALL full2sparse(A,irow,pcol,val,nrow,ncol)
    CALL sparse_solver_test(nrow,ncol,irow,pcol,val,x,b,max_abs_err,max_rel_err)

    IF (PRESENT(max_abs_err_out)) max_abs_err_out = max_abs_err
    IF (PRESENT(max_rel_err_out)) max_rel_err_out = max_rel_err
    IF (ALLOCATED(irow)) DEALLOCATE(irow)
    IF (ALLOCATED(pcol)) DEALLOCATE(pcol)
    IF (ALLOCATED(val))  DEALLOCATE(val)

  END SUBROUTINE sp_test_A_b2
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  ! tests A*x-b and returns errors
  SUBROUTINE sp_testComplex_A_b2(A,x,b,max_abs_err_out,max_rel_err_out)
    COMPLEX(kind=dp), DIMENSION(:,:), INTENT(in) :: A
    COMPLEX(kind=dp), DIMENSION(:,:), INTENT(in) :: x
    COMPLEX(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(in) :: b
    REAL(kind=dp), OPTIONAL, INTENT(out) :: max_abs_err_out,max_rel_err_out

    REAL(kind=dp) :: max_abs_err,max_rel_err
    INTEGER :: nrow,ncol
    INTEGER, DIMENSION(:), ALLOCATABLE :: irow,pcol
    COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE :: val

    CALL full2sparse(A,irow,pcol,val,nrow,ncol)
    CALL sparse_solver_test(nrow,ncol,irow,pcol,val,x,b,max_abs_err,max_rel_err)

    IF (PRESENT(max_abs_err_out)) max_abs_err_out = max_abs_err
    IF (PRESENT(max_rel_err_out)) max_rel_err_out = max_rel_err
    IF (ALLOCATED(irow)) DEALLOCATE(irow)
    IF (ALLOCATED(pcol)) DEALLOCATE(pcol)
    IF (ALLOCATED(val))  DEALLOCATE(val)

  END SUBROUTINE sp_testComplex_A_b2
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE remap_rc_real(nz,nz_sqeezed,irow,icol,amat)
    ! Re-arranges matrix elements which may be unordered and may have
    ! different elements with the same row and column indices is such
    ! a way that column index, icol, forms a non-decreasing sequence
    ! and row index, irow, forms increasing sub-sequences for itervals
    ! with a fixed column index. Sums up elements of the matrix which
    ! have the same row and column indices to one element with these
    ! indices
    !
    ! Arguments:
    ! nz          - (input)  number of elements in irow,icol,amat
    ! nz_sqeezed  - (output) number of elements with different (irow(k),icol(k))
    ! irow        - (inout)  row indices
    ! icol        - (inout)  column indices
    ! amat        - (inout)  matrix values

    INTEGER, INTENT(in)                          :: nz
    INTEGER, INTENT(out)                         :: nz_sqeezed
    INTEGER, DIMENSION(nz), INTENT(inout)        :: irow,icol
    REAL(kind=dp), DIMENSION(nz), INTENT(inout)  :: amat

    INTEGER                            :: ncol,i,j,k,kbeg,kend,ips,iflag,ksq
    INTEGER, DIMENSION(:), ALLOCATABLE :: nrows,icount,ipoi
    INTEGER                            :: ksq_ne0
    INTEGER, DIMENSION(:), ALLOCATABLE :: kne0

    ncol=MAXVAL(icol)
    ALLOCATE(nrows(ncol),icount(ncol),ipoi(nz))
    nrows=0

    ! count number of rows in a given column:
    DO k=1,nz
       j=icol(k)
       nrows(j)=nrows(j)+1
    ENDDO

    ! compute starting index - 1 of rows in a general list for each column:

    icount(1)=0

    DO i=1,ncol-1
       icount(i+1)=icount(i)+nrows(i)
    ENDDO

    ! compute the pointer from the list ordered by columns to a general list

    DO k=1,nz
       j=icol(k)
       icount(j)=icount(j)+1
       ipoi(icount(j))=k
    ENDDO

    ! re-order row indices to non-decreasing sub-sequences

    DO i=1,ncol
       kend=icount(i)
       kbeg=kend-nrows(i)+1
       DO j=1,kend-kbeg
          iflag=0
          DO k=kbeg+1,kend
             IF(irow(ipoi(k)).LT.irow(ipoi(k-1))) THEN
                iflag=1
                ips=ipoi(k)
                ipoi(k)=ipoi(k-1)
                ipoi(k-1)=ips
             ENDIF
          ENDDO
          IF(iflag.EQ.0) EXIT
       ENDDO
    ENDDO

    irow=irow(ipoi)
    icol=icol(ipoi)
    amat=amat(ipoi)

    ! squeese the data - sum up matrix elements with the same indices

    ksq=1
    !
    DO k=2,nz
       IF(irow(k).EQ.irow(k-1).AND.icol(k).EQ.icol(k-1)) THEN
          amat(ksq)=amat(ksq)+amat(k)
       ELSE
          ksq=ksq+1
          irow(ksq)=irow(k)
          icol(ksq)=icol(k)
          amat(ksq)=amat(k)
       ENDIF
    ENDDO

    ! remove zeros from the sparse vector

    ALLOCATE(kne0(ksq))
    ksq_ne0=0
    DO k=1,ksq
       IF(amat(k) .NE. 0.0d0) THEN
          ksq_ne0=ksq_ne0+1
          kne0(ksq_ne0)=k
       ENDIF
    ENDDO
    IF(ksq_ne0 .EQ. 0) THEN
       PRINT *,'sparse_mod.f90/remap_rc: All entries of the sparse vector are zero!'
    ELSE
       irow(1:ksq_ne0)=irow(kne0(1:ksq_ne0))
       icol(1:ksq_ne0)=icol(kne0(1:ksq_ne0))
       amat(1:ksq_ne0)=amat(kne0(1:ksq_ne0))
    ENDIF

    nz_sqeezed=ksq_ne0
    DEALLOCATE(nrows,icount,ipoi,kne0)

  END SUBROUTINE remap_rc_real
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE remap_rc_cmplx(nz,nz_sqeezed,irow,icol,amat)
    ! Re-arranges matrix elements which may be unordered and may have
    ! different elements with the same row and column indices is such
    ! a way that column index, icol, forms a non-decreasing sequence
    ! and row index, irow, forms increasing sub-sequences for itervals
    ! with a fixed column index. Sums up elements of the matrix which
    ! have the same row and column indices to one element with these
    ! indices
    !
    ! Arguments:
    ! nz          - (input)  number of elements in irow,icol,amat
    ! nz_sqeezed  - (output) number of elements with different (irow(k),icol(k))
    ! irow        - (inout)  row indices
    ! icol        - (inout)  column indices
    ! amat        - (inout)  matrix values

    INTEGER, INTENT(in)                          :: nz
    INTEGER, INTENT(out)                         :: nz_sqeezed
    INTEGER, DIMENSION(nz), INTENT(inout)        :: irow,icol
    COMPLEX(kind=kind(1d0)), DIMENSION(nz), INTENT(inout) :: amat

    INTEGER                            :: ncol,i,j,k,kbeg,kend,ips,iflag,ksq
    INTEGER, DIMENSION(:), ALLOCATABLE :: nrows,icount,ipoi
    INTEGER                            :: ksq_ne0
    INTEGER, DIMENSION(:), ALLOCATABLE :: kne0

    ncol=MAXVAL(icol)
    ALLOCATE(nrows(ncol),icount(ncol),ipoi(nz))
    nrows=0

    ! count number of rows in a given column:

    DO k=1,nz
       j=icol(k)
       nrows(j)=nrows(j)+1
    ENDDO

    ! compute starting index - 1 of rows in a general list for each column:

    icount(1)=0
    !
    DO i=1,ncol-1
       icount(i+1)=icount(i)+nrows(i)
    ENDDO

    ! compute the pointer from the list ordered by columns to a general list

    DO k=1,nz
       j=icol(k)
       icount(j)=icount(j)+1
       ipoi(icount(j))=k
    ENDDO

    ! re-order row indices to non-decreasing sub-sequences

    DO i=1,ncol
       kend=icount(i)
       kbeg=kend-nrows(i)+1
       DO j=1,kend-kbeg
          iflag=0
          DO k=kbeg+1,kend
             IF(irow(ipoi(k)).LT.irow(ipoi(k-1))) THEN
                iflag=1
                ips=ipoi(k)
                ipoi(k)=ipoi(k-1)
                ipoi(k-1)=ips
             ENDIF
          ENDDO
          IF(iflag.EQ.0) EXIT
       ENDDO
    ENDDO

    irow=irow(ipoi)
    icol=icol(ipoi)
    amat=amat(ipoi)

    ! squeese the data - sum up matrix elements with the same indices

    ksq=1

    DO k=2,nz
       IF(irow(k).EQ.irow(k-1).AND.icol(k).EQ.icol(k-1)) THEN
          amat(ksq)=amat(ksq)+amat(k)
       ELSE
          ksq=ksq+1
          irow(ksq)=irow(k)
          icol(ksq)=icol(k)
          amat(ksq)=amat(k)
       ENDIF
    ENDDO

    ! remove zeros from the sparse vector

    ALLOCATE(kne0(ksq))
    ksq_ne0=0
    DO k=1,ksq
       IF(amat(k) .NE. (0.d0,0.d0)) THEN
          ksq_ne0=ksq_ne0+1
          kne0(ksq_ne0)=k
       ENDIF
    ENDDO
    IF(ksq_ne0 .EQ. 0) THEN
       PRINT *,'sparse_mod.f90/remap_rc: All entries of the sparse vector are zero!'
    ELSE
       irow(1:ksq_ne0)=irow(kne0(1:ksq_ne0))
       icol(1:ksq_ne0)=icol(kne0(1:ksq_ne0))
       amat(1:ksq_ne0)=amat(kne0(1:ksq_ne0))
    ENDIF

    nz_sqeezed=ksq_ne0
    DEALLOCATE(nrows,icount,ipoi,kne0)

  END SUBROUTINE remap_rc_cmplx
  !-------------------------------------------------------------------------------

END MODULE sparse_mod
