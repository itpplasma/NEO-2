program test_spline_ilu_bicgstab
  ! Test ILU-preconditioned BiCGSTAB on the ill-conditioned spline problem
  
  use nrtype, only: I4B, DP
  use spline_mod, only: spline_cof  
  use sparse_types_mod, only: sparse_matrix_csc_real
  use sparse_conversion_mod, only: csc_to_csr
  use sparse_solvers_mod, only: bicgstab_abs_tolerance, bicgstab_rel_tolerance, &
                                bicgstab_max_iter, sparse_solve_method, SOLVER_BICGSTAB, &
                                SOLVER_UMFPACK
  use sparse_mod, only: sparse_solve
  use ilu_precond_mod
  use bicgstab_mod
  implicit none

  ! Problem size
  integer(I4B), parameter :: N_knots = 60
  integer(I4B), parameter :: N_coefs = N_knots + 2  ! 62 coefficients
  
  ! Spline problem variables
  real(DP), allocatable :: x_data(:), y_data(:), coeffs_ref(:), coeffs_ilu(:)
  type(sparse_matrix_csc_real) :: A_csc
  
  ! CSR conversion variables
  integer(I4B), allocatable :: csr_row_ptr(:), csr_col_idx(:)
  real(DP), allocatable :: csr_val(:)
  
  ! ILU variables
  type(ilu_factorization) :: ilu_fac
  
  ! BiCGSTAB variables
  type(bicgstab_stats) :: stats
  real(DP) :: abs_tol, rel_tol
  integer :: max_iter, iter
  logical :: converged
  
  ! Test variables
  real(DP) :: error_norm, max_error
  integer :: i, info
  
  ! Save current settings
  real(DP) :: saved_abs_tol, saved_rel_tol
  integer :: saved_max_iter, saved_method
  
  write(*,'(A)') '=== ILU-Preconditioned BiCGSTAB Test for Spline Problems ==='
  write(*,'(A)') ''
  write(*,'(A)') 'Testing whether ILU preconditioning can help BiCGSTAB converge'
  write(*,'(A)') 'on the ill-conditioned 348x348 spline matrix.'
  write(*,'(A)') ''
  
  ! Save current solver settings
  saved_abs_tol = bicgstab_abs_tolerance
  saved_rel_tol = bicgstab_rel_tolerance
  saved_max_iter = bicgstab_max_iter
  saved_method = sparse_solve_method
  
  ! Create test data for spline fitting (cubic spline)
  allocate(x_data(N_knots), y_data(N_knots))
  do i = 1, N_knots
    x_data(i) = real(i-1, DP) / real(N_knots-1, DP)  ! 0 to 1
    y_data(i) = sin(2.0_DP * 3.14159_DP * x_data(i))  ! Sine wave
  end do
  
  ! Generate spline coefficient matrix using UMFPACK (reference solution)
  write(*,'(A)') 'Step 1: Computing reference solution with UMFPACK...'
  sparse_solve_method = SOLVER_UMFPACK
  
  allocate(coeffs_ref(N_coefs))
  call spline_cof(x_data, y_data, N_knots, 3, coeffs_ref, A_csc)  ! cubic spline
  
  write(*,'(A,I0,A,I0,A,I0)') '  Matrix size: ', A_csc%nrow, 'x', A_csc%ncol, ' with ', A_csc%nz, ' non-zeros'
  
  ! Convert CSC to CSR for BiCGSTAB
  write(*,'(A)') 'Step 2: Converting CSC to CSR format...'
  call csc_to_csr(A_csc%nrow, A_csc%ncol, A_csc%nz, A_csc%irow, A_csc%pcol, A_csc%amat, &
                  csr_row_ptr, csr_col_idx, csr_val)
  
  ! Create RHS from the reference solution (A * coeffs_ref)
  real(DP), allocatable :: rhs(:)
  integer :: j
  allocate(rhs(A_csc%nrow))
  
  ! Compute RHS = A * coeffs_ref using CSR format
  rhs = 0.0_DP
  do i = 1, A_csc%nrow
    do j = csr_row_ptr(i), csr_row_ptr(i+1) - 1
      rhs(i) = rhs(i) + csr_val(j) * coeffs_ref(csr_col_idx(j))
    end do
  end do
  
  write(*,'(A)') 'Step 3: Setting up ILU(0) preconditioning...'
  
  ! Compute ILU(0) factorization
  call ilu_factorize(A_csc%nrow, A_csc%nz, csr_row_ptr, csr_col_idx, csr_val, &
                     0, 0.0_DP, ilu_fac, info)  ! ILU(0) with no drop tolerance
  
  if (info /= 0) then
    write(*,'(A,I0)') 'ERROR: ILU factorization failed with info = ', info
    stop 1
  end if
  
  write(*,'(A,I0,A,I0)') '  ILU factorization successful: L has ', ilu_fac%L_nnz, &
                         ' non-zeros, U has ', ilu_fac%U_nnz, ' non-zeros'
  
  ! Test ILU-preconditioned BiCGSTAB
  write(*,'(A)') 'Step 4: Testing ILU-preconditioned BiCGSTAB...'
  
  ! Allocate solution vector and set initial guess
  allocate(coeffs_ilu(A_csc%nrow))
  coeffs_ilu = 0.0_DP  ! Zero initial guess  
  
  ! Set BiCGSTAB parameters
  abs_tol = 1.0e-8_DP
  rel_tol = 1.0e-8_DP
  max_iter = 2000
  
  write(*,'(A,E10.2,A,E10.2,A,I0)') '  Tolerances: abs_tol=', abs_tol, ', rel_tol=', rel_tol, &
                                     ', max_iter=', max_iter
  
  ! Call preconditioned BiCGSTAB
  call bicgstab_solve_precond_real(A_csc%nrow, csr_row_ptr, csr_col_idx, csr_val, &
                                   rhs, coeffs_ilu, abs_tol, rel_tol, max_iter, &
                                   ilu_fac, converged, iter, stats)
  
  ! Analyze results
  write(*,'(A)') ''
  write(*,'(A)') 'Results:'
  write(*,'(A,L1)') '  Converged: ', converged
  write(*,'(A,I0,A,I0)') '  Iterations: ', iter, ' / ', max_iter
  write(*,'(A,E12.4)') '  Final residual: ', stats%final_residual
  write(*,'(A,E12.4)') '  Initial residual: ', stats%initial_residual
  write(*,'(A,F8.3,A)') '  Solve time: ', stats%solve_time, ' seconds'
  
  if (converged) then
    ! Compare with reference solution
    error_norm = sqrt(dot_product(coeffs_ilu - coeffs_ref, coeffs_ilu - coeffs_ref))
    max_error = maxval(abs(coeffs_ilu - coeffs_ref))
    
    write(*,'(A,E12.4)') '  Solution error (L2 norm): ', error_norm
    write(*,'(A,E12.4)') '  Solution error (max): ', max_error
    
    if (max_error < 1.0e-6_DP) then
      write(*,'(A)') ''
      write(*,'(A)') '✓ SUCCESS: ILU-preconditioned BiCGSTAB converged with good accuracy!'
    else
      write(*,'(A)') ''
      write(*,'(A)') '⚠ PARTIAL SUCCESS: Converged but with limited accuracy'
    end if
  else
    write(*,'(A)') ''
    write(*,'(A)') '✗ FAILURE: ILU-preconditioned BiCGSTAB did not converge'
  end if
  
  ! Cleanup
  call ilu_free(ilu_fac)
  
  ! Restore settings
  bicgstab_abs_tolerance = saved_abs_tol
  bicgstab_rel_tolerance = saved_rel_tol
  bicgstab_max_iter = saved_max_iter
  sparse_solve_method = saved_method
  
end program test_spline_ilu_bicgstab