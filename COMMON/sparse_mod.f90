MODULE sparse_mod

  USE sparse_types_mod, ONLY: dp, long
  USE sparse_conversion_mod
  USE sparse_io_mod
  USE sparse_arithmetic_mod
  USE sparse_solvers_mod
  IMPLICIT NONE

  ! Re-export sparse_solve_method for backward compatibility
  PUBLIC :: sparse_solve_method
  
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

  PUBLIC remap_rc
  INTERFACE remap_rc
     MODULE PROCEDURE remap_rc_real, remap_rc_cmplx
  END INTERFACE remap_rc

  ! helper


CONTAINS


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
  ! All solver operations have been moved to sparse_solvers_mod
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