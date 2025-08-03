program test_simple_rap
  use sparse_types_mod
  use amg_types_mod
  use amg_smoothed_aggregation_mod
  implicit none

  ! Test RAP computation on a simple 4x4 matrix with known 2x2 result
  integer, parameter :: n_fine = 4, n_coarse = 2
  type(amg_level) :: level_fine, level_coarse
  integer :: i, j
  
  print *, "=== Testing Simple RAP Computation ==="
  
  ! Create simple fine matrix: 4x4 identity matrix
  level_fine%n = n_fine
  level_fine%nnz = n_fine
  allocate(level_fine%row_ptr(n_fine + 1))
  allocate(level_fine%col_idx(n_fine))
  allocate(level_fine%values(n_fine))
  allocate(level_fine%diagonal(n_fine))
  
  ! Identity matrix: A = diag(1,1,1,1)
  level_fine%row_ptr = [1, 2, 3, 4, 5]
  level_fine%col_idx = [1, 2, 3, 4]
  level_fine%values = [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp]
  level_fine%diagonal = [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp]
  
  ! Simple aggregation: nodes 1,2 -> aggregate 1; nodes 3,4 -> aggregate 2
  allocate(level_fine%aggregates(n_fine))
  level_fine%aggregates = [1, 1, 2, 2]
  level_fine%n_aggregates = n_coarse
  
  ! Build tentative prolongation manually for this test
  level_fine%n_fine = n_fine
  level_fine%n_coarse = n_coarse
  level_fine%nnz_p = n_fine
  allocate(level_fine%P_row_ptr(n_fine + 1))
  allocate(level_fine%P_col_idx(n_fine))
  allocate(level_fine%P_values(n_fine))
  
  ! P matrix: nodes 1,2 -> coarse node 1; nodes 3,4 -> coarse node 2
  ! With normalization: P[i,j] = 1/sqrt(2) for each aggregate
  level_fine%P_row_ptr = [1, 2, 3, 4, 5]
  level_fine%P_col_idx = [1, 1, 2, 2]
  level_fine%P_values = [1.0_dp/sqrt(2.0_dp), 1.0_dp/sqrt(2.0_dp), &
                        1.0_dp/sqrt(2.0_dp), 1.0_dp/sqrt(2.0_dp)]
  
  print *, "Fine matrix (4x4 identity):"
  do i = 1, n_fine
    print '(4F8.3)', (level_fine%values(level_fine%row_ptr(i):level_fine%row_ptr(i+1)-1))
  end do
  
  print *, "Aggregation:", level_fine%aggregates
  print *, "Prolongation values:", level_fine%P_values
  
  ! Skip strength matrix and smoothing for this simple test
  ! Just build restriction as P^T
  call build_restriction(level_fine)
  
  ! Compute RAP
  level_coarse%n = n_coarse
  call sparse_rap_product(level_fine, level_coarse)
  call extract_diagonal(level_coarse)
  
  print *, "Coarse matrix size:", level_coarse%n
  print *, "Coarse matrix nnz:", level_coarse%nnz
  print *, "Coarse diagonal:", level_coarse%diagonal
  
  ! Expected result: 2x2 matrix with diagonal entries = 1/2
  ! Since R*A*P should give: (1/sqrt(2))^2 * 1 * 2_nodes = 1/2 for each diagonal
  print *, "Expected diagonal: [0.5, 0.5]"
  
  if (all(abs(level_coarse%diagonal - 0.5_dp) < 1.0e-12_dp)) then
    print *, "PASS: RAP computation correct"
  else
    print *, "FAIL: RAP computation incorrect"
  end if
  
  ! Cleanup
  deallocate(level_fine%row_ptr, level_fine%col_idx, level_fine%values)
  deallocate(level_fine%diagonal, level_fine%aggregates)
  deallocate(level_fine%P_row_ptr, level_fine%P_col_idx, level_fine%P_values)
  if (allocated(level_fine%R_row_ptr)) deallocate(level_fine%R_row_ptr)
  if (allocated(level_fine%R_col_idx)) deallocate(level_fine%R_col_idx)
  if (allocated(level_fine%R_values)) deallocate(level_fine%R_values)
  if (allocated(level_coarse%row_ptr)) deallocate(level_coarse%row_ptr)
  if (allocated(level_coarse%col_idx)) deallocate(level_coarse%col_idx)
  if (allocated(level_coarse%values)) deallocate(level_coarse%values)
  if (allocated(level_coarse%diagonal)) deallocate(level_coarse%diagonal)

end program test_simple_rap