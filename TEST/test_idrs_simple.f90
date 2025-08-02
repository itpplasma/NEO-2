PROGRAM test_idrs_simple
  !> Simple test of IDR(s) solver
  !! Tests IDR(s) on a simple matrix to verify implementation
  
  USE nrtype, ONLY: DP
  USE idrs_mod, ONLY: idrs_solve_real
  IMPLICIT NONE
  
  INTEGER, PARAMETER :: n = 10
  REAL(DP) :: A(n, n), b(n), x(n), x_exact(n)
  REAL(DP) :: residual, tolerance
  INTEGER :: i, j, stats
  LOGICAL :: test_passed
  
  WRITE(*,*) '=== IDR(s) Simple Test ==='
  
  ! Create a simple test matrix: A = I + 0.1 * J (where J is all ones)
  A = 0.0_DP
  DO i = 1, n
    A(i, i) = 1.0_DP
    DO j = 1, n
      A(i, j) = A(i, j) + 0.1_DP
    END DO
  END DO
  
  ! Create exact solution
  DO i = 1, n
    x_exact(i) = REAL(i, DP)
  END DO
  
  ! Compute right-hand side: b = A * x_exact
  DO i = 1, n
    b(i) = 0.0_DP
    DO j = 1, n
      b(i) = b(i) + A(i, j) * x_exact(j)
    END DO
  END DO
  
  ! Initial guess
  x = 0.0_DP
  
  ! Solve with IDR(s)
  CALL idrs_solve_real(A, b, x, 4, stats)
  
  ! Check residual
  residual = 0.0_DP
  DO i = 1, n
    residual = residual + (x(i) - x_exact(i))**2
  END DO
  residual = SQRT(residual)
  
  tolerance = 1.0E-10_DP
  test_passed = (residual < tolerance .AND. stats == 0)
  
  WRITE(*,'(A,ES12.5)') 'Solution error (2-norm): ', residual
  WRITE(*,'(A,I0)') 'IDR(s) status: ', stats
  WRITE(*,'(A,L1)') 'Test passed: ', test_passed
  
  IF (.NOT. test_passed) THEN
    WRITE(*,*) 'ERROR: IDR(s) simple test failed!'
    WRITE(*,*) 'Expected solution:'
    DO i = 1, n
      WRITE(*,'(I3,F12.6)') i, x_exact(i)
    END DO
    WRITE(*,*) 'IDR(s) solution:'
    DO i = 1, n
      WRITE(*,'(I3,F12.6)') i, x(i)
    END DO
    STOP 1
  END IF
  
  WRITE(*,*) 'IDR(s) simple test passed!'
  
END PROGRAM test_idrs_simple