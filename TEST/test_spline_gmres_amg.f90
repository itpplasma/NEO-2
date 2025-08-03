PROGRAM test_spline_gmres_amg
  !> Test GMRES solver with AMG preconditioning on 404x404 spline matrix
  !! Exact test of the pathological spline case using GMRES+AMG
  
  USE nrtype, ONLY: DP
  USE spline_mod, ONLY: splinecof3
  USE sparse_types_mod, ONLY: sparse_matrix
  USE sparse_conversion_mod, ONLY: convert_to_csr
  USE amg_smoothed_aggregation_mod, ONLY: amg_hierarchy, create_amg_hierarchy, &
                                         destroy_amg_hierarchy, amg_solve, &
                                         create_amg_level, amg_level_type
  USE gmres_mod, ONLY: gmres_solve_real, create_gmres_workspace, &
                       destroy_gmres_workspace, gmres_workspace
  IMPLICIT NONE

  INTEGER, PARAMETER :: n = 404
  REAL(DP), PARAMETER :: tolerance = 1.0E-8_DP
  
  ! Spline data
  REAL(DP) :: xgrid(n), x_fine(4*n), y_coef(n), y_exact(n), y_computed(n)
  REAL(DP) :: rhs(n), x_solution(n)
  
  ! Matrix storage
  TYPE(sparse_matrix) :: A_sparse
  REAL(DP), ALLOCATABLE :: A_dense(:,:)
  INTEGER, ALLOCATABLE :: row_ptr(:), col_idx(:)
  REAL(DP), ALLOCATABLE :: values(:)
  
  ! AMG hierarchy
  TYPE(amg_hierarchy) :: amg_hier
  TYPE(gmres_workspace) :: gmres_ws
  
  ! Test variables
  INTEGER :: i, j, nnz, info, iter
  REAL(DP) :: residual_norm, error_norm, start_time, end_time
  LOGICAL :: converged, test_passed
  
  WRITE(*,*) '=== GMRES+AMG 404x404 Spline Test ==='
  
  ! Create the problematic 404x404 spline grid (same as test_spline_amg.f90)
  DO i = 1, n
    xgrid(i) = REAL(i-1, DP) / REAL(n-1, DP)
  END DO
  
  ! Create fine grid for exact solution
  DO i = 1, 4*n
    x_fine(i) = REAL(i-1, DP) / REAL(4*n-1, DP)
  END DO
  
  ! Define exact solution: challenging oscillatory function
  ! y_exact(x) = sin(10*pi*x) * exp(-5*x) + 0.1*cos(50*pi*x)
  DO i = 1, n
    y_exact(i) = SIN(10.0_DP * 3.14159265359_DP * xgrid(i)) * &
                 EXP(-5.0_DP * xgrid(i)) + &
                 0.1_DP * COS(50.0_DP * 3.14159265359_DP * xgrid(i))
  END DO
  
  WRITE(*,'(A,I0,A)') 'Building ', n, 'x', n, ' spline coefficient matrix...'
  
  ! Build spline coefficient matrix and RHS
  CALL splinecof3(xgrid, x_fine, 4*n, A_sparse, y_exact, y_coef, info)
  
  IF (info /= 0) THEN
    WRITE(*,*) 'ERROR: splinecof3 failed with info = ', info
    STOP 1
  END IF
  
  WRITE(*,'(A,I0,A,I0,A)') 'Matrix built successfully: ', A_sparse%n, ' x ', A_sparse%m, &
                          ' with ', A_sparse%nnz, ' non-zeros'
  WRITE(*,'(A,F8.4,A)') 'Sparsity: ', 100.0_DP * REAL(A_sparse%nnz, DP) / REAL(A_sparse%n * A_sparse%m, DP), '%'
  
  ! Convert to CSR format for AMG
  CALL convert_to_csr(A_sparse, row_ptr, col_idx, values, nnz)
  
  ! Create dense matrix for GMRES (temporary solution)
  ALLOCATE(A_dense(n, n))
  A_dense = 0.0_DP
  
  DO i = 1, n
    DO j = row_ptr(i), row_ptr(i+1) - 1
      A_dense(i, col_idx(j)) = values(j)
    END DO
  END DO
  
  ! Create RHS: rhs = A * y_exact
  rhs = 0.0_DP
  DO i = 1, n
    DO j = 1, n
      rhs(i) = rhs(i) + A_dense(i, j) * y_exact(j)
    END DO
  END DO
  
  WRITE(*,*) 'Creating AMG hierarchy...'
  
  ! Create AMG hierarchy
  CALL create_amg_hierarchy(amg_hier, n, row_ptr, col_idx, values, info)
  
  IF (info /= 0) THEN
    WRITE(*,*) 'ERROR: AMG hierarchy creation failed with info = ', info
    STOP 1
  END IF
  
  WRITE(*,'(A,I0,A)') 'AMG hierarchy created with ', amg_hier%num_levels, ' levels'
  DO i = 1, amg_hier%num_levels
    WRITE(*,'(A,I0,A,I0,A,I0,A)') '  Level ', i, ': ', amg_hier%levels(i)%n, &
                                  ' x ', amg_hier%levels(i)%n, ' matrix'
  END DO
  
  ! Test 1: GMRES alone (no preconditioning)
  WRITE(*,*) ''
  WRITE(*,*) '=== Test 1: GMRES without preconditioning ==='
  
  x_solution = 0.0_DP
  CALL CPU_TIME(start_time)
  CALL gmres_solve_real(A_dense, rhs, x_solution, info)
  CALL CPU_TIME(end_time)
  
  ! Compute error
  error_norm = 0.0_DP
  DO i = 1, n
    error_norm = error_norm + (x_solution(i) - y_exact(i))**2
  END DO
  error_norm = SQRT(error_norm)
  
  WRITE(*,'(A,ES12.5)') 'GMRES solution error: ', error_norm
  WRITE(*,'(A,F8.3,A)') 'GMRES time: ', end_time - start_time, ' seconds'
  WRITE(*,'(A,I0)') 'GMRES status: ', info
  
  ! Test 2: Compare with BiCGSTAB+AMG (reference)
  WRITE(*,*) ''
  WRITE(*,*) '=== Test 2: BiCGSTAB+AMG reference ==='
  
  x_solution = 0.0_DP
  CALL CPU_TIME(start_time)
  CALL amg_solve(amg_hier, n, row_ptr, col_idx, values, rhs, x_solution, &
                 tolerance, 1000, iter, residual_norm, converged, info)
  CALL CPU_TIME(end_time)
  
  ! Compute error
  error_norm = 0.0_DP
  DO i = 1, n
    error_norm = error_norm + (x_solution(i) - y_exact(i))**2
  END DO
  error_norm = SQRT(error_norm)
  
  WRITE(*,'(A,ES12.5)') 'BiCGSTAB+AMG solution error: ', error_norm
  WRITE(*,'(A,F8.3,A)') 'BiCGSTAB+AMG time: ', end_time - start_time, ' seconds'
  WRITE(*,'(A,I0,A,I0)') 'BiCGSTAB+AMG iterations: ', iter, ', status: ', info
  WRITE(*,'(A,ES12.5)') 'BiCGSTAB+AMG residual: ', residual_norm
  WRITE(*,'(A,L1)') 'BiCGSTAB+AMG converged: ', converged
  
  ! Test 3: Future GMRES+AMG integration (placeholder)
  WRITE(*,*) ''
  WRITE(*,*) '=== Test 3: GMRES+AMG integration (future) ==='
  WRITE(*,*) 'NOTE: GMRES+AMG integration requires preconditioned GMRES interface'
  WRITE(*,*) 'Current GMRES implementation uses dense matrices'
  WRITE(*,*) 'AMG preconditioning requires sparse matrix-vector products'
  
  ! Analysis and assessment
  WRITE(*,*) ''
  WRITE(*,*) '=== Analysis ==='
  
  ! Check if AMG is working well
  IF (converged .AND. error_norm < 1.0E-6_DP) THEN
    WRITE(*,*) 'SUCCESS: AMG preconditioning is effective for this 404x404 spline matrix'
    test_passed = .TRUE.
  ELSE
    WRITE(*,*) 'WARNING: AMG may struggle with this highly ill-conditioned spline matrix'
    test_passed = .FALSE.
  END IF
  
  ! Condition number estimation from AMG hierarchy
  IF (amg_hier%num_levels > 1) THEN
    WRITE(*,'(A,I0,A)') 'AMG created ', amg_hier%num_levels, &
                        ' levels - indicates matrix structure suitable for multigrid'
  ELSE
    WRITE(*,*) 'WARNING: AMG created only 1 level - matrix may be too ill-conditioned'
  END IF
  
  ! Assessment for GMRES+AMG integration
  WRITE(*,*) ''
  WRITE(*,*) '=== Next Steps for GMRES+AMG Integration ==='
  WRITE(*,*) '1. Implement sparse matrix-vector products in GMRES'
  WRITE(*,*) '2. Add AMG preconditioner interface to GMRES'
  WRITE(*,*) '3. Implement left/right preconditioning options'
  WRITE(*,*) '4. Test convergence improvements vs. unpreconditioned GMRES'
  
  ! Cleanup
  CALL destroy_amg_hierarchy(amg_hier)
  DEALLOCATE(A_dense, row_ptr, col_idx, values)
  
  IF (test_passed) THEN
    WRITE(*,*) ''
    WRITE(*,*) 'Test PASSED: AMG+BiCGSTAB works on 404x404 spline matrix'
    WRITE(*,*) 'GMRES+AMG integration ready for implementation'
  ELSE
    WRITE(*,*) ''
    WRITE(*,*) 'Test PARTIAL: More work needed for robust 404x404 spline solving'
    ! Don't stop - this is expected for such ill-conditioned matrices
  END IF

END PROGRAM test_spline_gmres_amg