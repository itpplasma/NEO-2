MODULE sparse_arithmetic_mod
  ! Module containing sparse matrix arithmetic operations
  ! Extracted from sparse_mod.f90 for better modularity
  
  USE sparse_types_mod, ONLY: dp
  USE sparse_conversion_mod, ONLY: column_pointer2full, full2sparse
  IMPLICIT NONE
  
  PUBLIC :: sparse_matmul
  PUBLIC :: sparse_solver_test
  PUBLIC :: sparse_talk
  
  LOGICAL :: sparse_talk = .FALSE.
  
  INTERFACE sparse_matmul
    MODULE PROCEDURE sp_matmul_A_b1, sp_matmul_b1, sp_matmul_A_b2, sp_matmul_b2, &
                     sp_matmulComplex_A_b1, sp_matmulComplex_b1, sp_matmulComplex_A_b2, sp_matmulComplex_b2
  END INTERFACE sparse_matmul
  
  INTERFACE sparse_solver_test
    MODULE PROCEDURE sp_test_A_b1, sp_test_b1, sp_test_A_b2, sp_test_b2, &
                     sp_testComplex_A_b1, sp_testComplex_b1, sp_testComplex_A_b2, sp_testComplex_b2
  END INTERFACE sparse_solver_test
  
CONTAINS
  
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
    IF (SIZE(pcol,1) .EQ. ncol+1) THEN
       ! pcol is column pointers, convert to column indices
       CALL column_pointer2full(pcol,icol)
    ELSE IF (SIZE(pcol,1) .EQ. nz) THEN
       ! pcol is already column indices
       ALLOCATE(icol(SIZE(pcol)))
       icol = pcol
    ELSE
       PRINT *, 'Error in sparse_matmul: pcol size is not correct'
       STOP
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
    IF (SIZE(pcol,1) .EQ. ncol+1) THEN
       ! pcol is column pointers, convert to column indices
       CALL column_pointer2full(pcol,icol)
    ELSE IF (SIZE(pcol,1) .EQ. nz) THEN
       ! pcol is already column indices
       ALLOCATE(icol(SIZE(pcol)))
       icol = pcol
    ELSE
       PRINT *, 'Error in sparse_matmul: pcol size is not correct'
       STOP
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
    
    INTEGER :: nz,n,ic,ir,i_rhs,n_rhs
    INTEGER, DIMENSION(:), ALLOCATABLE :: icol
    
    
    
    nz = SIZE(val,1)
    IF (SIZE(pcol,1) .EQ. ncol+1) THEN
       ! pcol is column pointers, convert to column indices
       CALL column_pointer2full(pcol,icol)
    ELSE IF (SIZE(pcol,1) .EQ. nz) THEN
       ! pcol is already column indices
       ALLOCATE(icol(SIZE(pcol)))
       icol = pcol
    ELSE
       PRINT *, 'Error in sparse_matmul: pcol size is not correct'
       STOP
    END IF
    
    IF (ncol .NE. SIZE(x,1)) THEN
       PRINT *, 'Error in sparse_matmul: sizes do not fit'
       STOP
    END IF
    IF (ALLOCATED(r)) DEALLOCATE(r)
    n_rhs = SIZE(x,2)
    ALLOCATE(r(nrow,n_rhs))
    r = 0.0_dp
    
    DO i_rhs = 1,n_rhs
       DO n = 1,nz
          ic = icol(n)
          ir = irow(n)
          r(ir,i_rhs) = r(ir,i_rhs) + val(n)*x(ic,i_rhs)
       END DO
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
    
    INTEGER :: nz,n,ic,ir,i_rhs,n_rhs
    INTEGER, DIMENSION(:), ALLOCATABLE :: icol
    
    
    
    nz = SIZE(val,1)
    IF (SIZE(pcol,1) .EQ. ncol+1) THEN
       ! pcol is column pointers, convert to column indices
       CALL column_pointer2full(pcol,icol)
    ELSE IF (SIZE(pcol,1) .EQ. nz) THEN
       ! pcol is already column indices
       ALLOCATE(icol(SIZE(pcol)))
       icol = pcol
    ELSE
       PRINT *, 'Error in sparse_matmul: pcol size is not correct'
       STOP
    END IF
    
    IF (ncol .NE. SIZE(x,1)) THEN
       PRINT *, 'Error in sparse_matmul: sizes do not fit'
       STOP
    END IF
    IF (ALLOCATED(r)) DEALLOCATE(r)
    n_rhs = SIZE(x,2)
    ALLOCATE(r(nrow,n_rhs))
    r = 0.0_dp
    
    DO i_rhs = 1,n_rhs
       DO n = 1,nz
          ic = icol(n)
          ir = irow(n)
          r(ir,i_rhs) = r(ir,i_rhs) + val(n)*x(ic,i_rhs)
       END DO
    END DO
    IF (ALLOCATED(icol)) DEALLOCATE(icol)
    
  END SUBROUTINE sp_matmulComplex_b2
  !-------------------------------------------------------------------------------
  
  ! computes A*x for full A and 1-D array x
  SUBROUTINE sp_matmul_A_b1(A,x,r)
    REAL(kind=dp), DIMENSION(:,:), INTENT(in) :: A
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: x
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(inout) :: r
    
    INTEGER :: nrow,ncol,nz
    INTEGER, DIMENSION(:), ALLOCATABLE :: irow,pcol
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: val
    
    CALL full2sparse(A,irow,pcol,val,nrow,ncol,nz)
    CALL sparse_matmul(nrow,ncol,irow,pcol,val,x,r)
    IF (ALLOCATED(irow)) DEALLOCATE(irow)
    IF (ALLOCATED(pcol)) DEALLOCATE(pcol)
    IF (ALLOCATED(val)) DEALLOCATE(val)
    
  END SUBROUTINE sp_matmul_A_b1
  
  ! computes A*x for full A and 1-D array x
  SUBROUTINE sp_matmulComplex_A_b1(A,x,r)
    COMPLEX(kind=dp), DIMENSION(:,:), INTENT(in) :: A
    COMPLEX(kind=dp), DIMENSION(:), INTENT(in) :: x
    COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE, INTENT(inout) :: r
    
    INTEGER :: nrow,ncol,nz
    INTEGER, DIMENSION(:), ALLOCATABLE :: irow,pcol
    COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE :: val
    
    CALL full2sparse(A,irow,pcol,val,nrow,ncol,nz)
    CALL sparse_matmul(nrow,ncol,irow,pcol,val,x,r)
    IF (ALLOCATED(irow)) DEALLOCATE(irow)
    IF (ALLOCATED(pcol)) DEALLOCATE(pcol)
    IF (ALLOCATED(val)) DEALLOCATE(val)
    
  END SUBROUTINE sp_matmulComplex_A_b1
  
  ! computes A*x for full A and 2-D array x
  SUBROUTINE sp_matmul_A_b2(A,x,r)
    REAL(kind=dp), DIMENSION(:,:), INTENT(in) :: A
    REAL(kind=dp), DIMENSION(:,:), INTENT(in) :: x
    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: r
    
    INTEGER :: nrow,ncol,nz
    INTEGER, DIMENSION(:), ALLOCATABLE :: irow,pcol
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: val
    
    CALL full2sparse(A,irow,pcol,val,nrow,ncol,nz)
    CALL sparse_matmul(nrow,ncol,irow,pcol,val,x,r)
    IF (ALLOCATED(irow)) DEALLOCATE(irow)
    IF (ALLOCATED(pcol)) DEALLOCATE(pcol)
    IF (ALLOCATED(val)) DEALLOCATE(val)
    
  END SUBROUTINE sp_matmul_A_b2
  
  ! computes A*x for full A and 2-D array x
  SUBROUTINE sp_matmulComplex_A_b2(A,x,r)
    COMPLEX(kind=dp), DIMENSION(:,:), INTENT(in) :: A
    COMPLEX(kind=dp), DIMENSION(:,:), INTENT(in) :: x
    COMPLEX(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: r
    
    INTEGER :: nrow,ncol,nz
    INTEGER, DIMENSION(:), ALLOCATABLE :: irow,pcol
    COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE :: val
    
    CALL full2sparse(A,irow,pcol,val,nrow,ncol,nz)
    CALL sparse_matmul(nrow,ncol,irow,pcol,val,x,r)
    IF (ALLOCATED(irow)) DEALLOCATE(irow)
    IF (ALLOCATED(pcol)) DEALLOCATE(pcol)
    IF (ALLOCATED(val)) DEALLOCATE(val)
    
  END SUBROUTINE sp_matmulComplex_A_b2
  
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
    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: r
    
    INTEGER :: n_rhs,i_rhs
    
    n_rhs = SIZE(x,2)
    CALL sparse_matmul(nrow,ncol,irow,pcol,val,x,r)
    max_abs_err = 0.0_dp
    max_rel_err = 0.0_dp
    DO i_rhs = 1,n_rhs
       r(:,i_rhs) = ABS(r(:,i_rhs) - b(:,i_rhs))
       max_abs_err = MAX(max_abs_err,MAXVAL(r(:,i_rhs)))
       max_rel_err = MAX(max_rel_err,max_abs_err / MAXVAL(ABS(b(:,i_rhs))))
    END DO
    IF (sparse_talk) THEN
       PRINT *, 'max_abs_err=',max_abs_err,' max_rel_err=',max_rel_err
    END IF
    
    IF (PRESENT(max_abs_err_out)) max_abs_err_out = max_abs_err
    IF (PRESENT(max_rel_err_out)) max_rel_err_out = max_rel_err
    IF (ALLOCATED(r)) DEALLOCATE(r)
    
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
    COMPLEX(kind=dp), DIMENSION(:,:), ALLOCATABLE :: r
    
    INTEGER :: n_rhs,i_rhs
    
    n_rhs = SIZE(x,2)
    CALL sparse_matmul(nrow,ncol,irow,pcol,val,x,r)
    max_abs_err = 0.0_dp
    max_rel_err = 0.0_dp
    DO i_rhs = 1,n_rhs
       r(:,i_rhs) = ABS(r(:,i_rhs) - b(:,i_rhs))
       max_abs_err = MAX(max_abs_err,MAXVAL(SQRT(REAL(r(:,i_rhs))**2+AIMAG(r(:,i_rhs))**2)))
       max_rel_err = MAX(max_rel_err,max_abs_err / MAXVAL(SQRT(REAL(b(:,i_rhs))**2+AIMAG(b(:,i_rhs))**2)))
    END DO
    IF (sparse_talk) THEN
       PRINT *, 'max_abs_err=',max_abs_err,' max_rel_err=',max_rel_err
    END IF
    
    IF (PRESENT(max_abs_err_out)) max_abs_err_out = max_abs_err
    IF (PRESENT(max_rel_err_out)) max_rel_err_out = max_rel_err
    IF (ALLOCATED(r)) DEALLOCATE(r)
    
  END SUBROUTINE sp_testComplex_b2
  !-------------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------------
  ! tests A*x-b and returns errors
  SUBROUTINE sp_test_A_b1(A,x,b,max_abs_err_out,max_rel_err_out)
    REAL(kind=dp), DIMENSION(:,:), INTENT(in) :: A
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: x
    REAL(kind=dp), DIMENSION(:), INTENT(in) :: b
    REAL(kind=dp), OPTIONAL, INTENT(out) :: max_abs_err_out,max_rel_err_out
    
    INTEGER :: nrow,ncol,nz
    INTEGER, DIMENSION(:), ALLOCATABLE :: irow,pcol
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: val
    
    CALL full2sparse(A,irow,pcol,val,nrow,ncol,nz)
    CALL sparse_solver_test(nrow,ncol,irow,pcol,val,x,b,max_abs_err_out,max_rel_err_out)
    IF (ALLOCATED(irow)) DEALLOCATE(irow)
    IF (ALLOCATED(pcol)) DEALLOCATE(pcol)
    IF (ALLOCATED(val)) DEALLOCATE(val)
    
  END SUBROUTINE sp_test_A_b1
  !-------------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------------
  ! tests A*x-b and returns errors
  SUBROUTINE sp_testComplex_A_b1(A,x,b,max_abs_err_out,max_rel_err_out)
    COMPLEX(kind=dp), DIMENSION(:,:), INTENT(in) :: A
    COMPLEX(kind=dp), DIMENSION(:), INTENT(in) :: x
    COMPLEX(kind=dp), DIMENSION(:), INTENT(in) :: b
    REAL(kind=dp), OPTIONAL, INTENT(out) :: max_abs_err_out,max_rel_err_out
    
    INTEGER :: nrow,ncol,nz
    INTEGER, DIMENSION(:), ALLOCATABLE :: irow,pcol
    COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE :: val
    
    CALL full2sparse(A,irow,pcol,val,nrow,ncol,nz)
    CALL sparse_solver_test(nrow,ncol,irow,pcol,val,x,b,max_abs_err_out,max_rel_err_out)
    IF (ALLOCATED(irow)) DEALLOCATE(irow)
    IF (ALLOCATED(pcol)) DEALLOCATE(pcol)
    IF (ALLOCATED(val)) DEALLOCATE(val)
    
  END SUBROUTINE sp_testComplex_A_b1
  !-------------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------------
  ! tests A*x-b and returns errors
  SUBROUTINE sp_test_A_b2(A,x,b,max_abs_err_out,max_rel_err_out)
    REAL(kind=dp), DIMENSION(:,:), INTENT(in) :: A
    REAL(kind=dp), DIMENSION(:,:), INTENT(in) :: x
    REAL(kind=dp), DIMENSION(:,:), INTENT(in) :: b
    REAL(kind=dp), OPTIONAL, INTENT(out) :: max_abs_err_out,max_rel_err_out
    
    INTEGER :: nrow,ncol,nz
    INTEGER, DIMENSION(:), ALLOCATABLE :: irow,pcol
    REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: val
    
    CALL full2sparse(A,irow,pcol,val,nrow,ncol,nz)
    CALL sparse_solver_test(nrow,ncol,irow,pcol,val,x,b,max_abs_err_out,max_rel_err_out)
    IF (ALLOCATED(irow)) DEALLOCATE(irow)
    IF (ALLOCATED(pcol)) DEALLOCATE(pcol)
    IF (ALLOCATED(val)) DEALLOCATE(val)
    
  END SUBROUTINE sp_test_A_b2
  !-------------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------------
  ! tests A*x-b and returns errors
  SUBROUTINE sp_testComplex_A_b2(A,x,b,max_abs_err_out,max_rel_err_out)
    COMPLEX(kind=dp), DIMENSION(:,:), INTENT(in) :: A
    COMPLEX(kind=dp), DIMENSION(:,:), INTENT(in) :: x
    COMPLEX(kind=dp), DIMENSION(:,:), INTENT(in) :: b
    REAL(kind=dp), OPTIONAL, INTENT(out) :: max_abs_err_out,max_rel_err_out
    
    INTEGER :: nrow,ncol,nz
    INTEGER, DIMENSION(:), ALLOCATABLE :: irow,pcol
    COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE :: val
    
    CALL full2sparse(A,irow,pcol,val,nrow,ncol,nz)
    CALL sparse_solver_test(nrow,ncol,irow,pcol,val,x,b,max_abs_err_out,max_rel_err_out)
    IF (ALLOCATED(irow)) DEALLOCATE(irow)
    IF (ALLOCATED(pcol)) DEALLOCATE(pcol)
    IF (ALLOCATED(val)) DEALLOCATE(val)
    
  END SUBROUTINE sp_testComplex_A_b2
  !-------------------------------------------------------------------------------
  
END MODULE sparse_arithmetic_mod