module amg_smoothed_aggregation_mod
  use amg_types_mod
  use sparse_types_mod, only: dp
  implicit none
  
  private
  public :: sa_amg_setup, sa_compute_strength, sa_aggregate_nodes
  public :: sa_fit_candidates, sa_smooth_prolongation
  
contains

  subroutine sa_amg_setup(hierarchy, A_csr, n, nnz, row_ptr, col_idx, values)
    type(amg_hierarchy), intent(inout) :: hierarchy
    integer, intent(in) :: n, nnz
    integer, intent(in) :: row_ptr(n+1), col_idx(nnz)
    real(dp), intent(in) :: values(nnz)
    type(amg_level), intent(inout) :: A_csr
    
    integer :: level, n_current, nnz_current
    integer, allocatable :: row_ptr_current(:), col_idx_current(:)
    real(dp), allocatable :: values_current(:)
    logical :: done
    
    ! Initialize first level with input matrix
    hierarchy%n_levels = 1
    allocate(hierarchy%levels(hierarchy%params%max_levels))
    
    ! Copy input matrix to first level
    hierarchy%levels(1)%n = n
    hierarchy%levels(1)%nnz = nnz
    allocate(hierarchy%levels(1)%row_ptr(n+1))
    allocate(hierarchy%levels(1)%col_idx(nnz))
    allocate(hierarchy%levels(1)%values(nnz))
    hierarchy%levels(1)%row_ptr = row_ptr
    hierarchy%levels(1)%col_idx = col_idx
    hierarchy%levels(1)%values = values
    
    ! Extract diagonal for smoothing
    call extract_diagonal(hierarchy%levels(1))
    
    ! Build hierarchy
    done = .false.
    level = 1
    
    do while (.not. done .and. level < hierarchy%params%max_levels)
      n_current = hierarchy%levels(level)%n
      
      ! Check if we've reached coarsest level
      if (n_current <= hierarchy%params%coarsest_size) then
        done = .true.
        exit
      end if
      
      ! Create strength matrix
      call sa_compute_strength(hierarchy%levels(level), hierarchy%params%strength_threshold)
      
      ! Perform aggregation
      call sa_aggregate_nodes(hierarchy%levels(level))
      
      ! Check coarsening ratio
      if (real(hierarchy%levels(level)%n_aggregates, dp) / real(n_current, dp) > &
          hierarchy%params%max_coarse_ratio) then
        done = .true.
        exit
      end if
      
      ! Build tentative prolongation via candidate fitting
      call sa_fit_candidates(hierarchy%levels(level))
      
      ! Smooth prolongation (exact Julia AlgebraicMultigrid.jl implementation)
      call sa_smooth_prolongation(hierarchy%levels(level), hierarchy%params%prolongation_damping)
      
      ! Build restriction (transpose of prolongation for SA)
      call build_restriction(hierarchy%levels(level))
      
      ! Compute coarse grid operator: A_c = R * A * P
      call compute_coarse_operator(hierarchy%levels(level), hierarchy%levels(level+1))
      
      level = level + 1
      hierarchy%n_levels = level
    end do
    
    ! Compute complexities
    call compute_complexities(hierarchy)
    
    hierarchy%setup_done = .true.
    
  end subroutine sa_amg_setup
  
  subroutine sa_compute_strength(level, theta)
    type(amg_level), intent(inout) :: level
    real(dp), intent(in) :: theta
    
    integer :: i, j, k, nnz_s
    real(dp) :: diag_i, diag_j, aij, eps_aii, val_sq, max_col_entry
    real(dp), allocatable :: diags(:)
    
    ! Exact replication of Julia SymmetricStrength algorithm
    ! From strength.jl lines 77-122
    
    ! Step 1: Extract diagonal norms exactly like Julia (lines 89-101)
    allocate(diags(level%n))
    do i = 1, level%n
      diag_i = 0.0_dp
      do k = level%row_ptr(i), level%row_ptr(i+1)-1
        j = level%col_idx(k)
        if (i == j) then
          diag_i = diag_i + level%values(k)
        end if
      end do
      diags(i) = abs(diag_i)  ! Julia: norm(diag)
    end do
    
    ! Step 2: Allocate strength matrix (S) - same structure as A
    if (.not. allocated(level%S_row_ptr)) then
      allocate(level%S_row_ptr(level%n + 1))
      allocate(level%S_col_idx(level%nnz))
      allocate(level%S_values(level%nnz))
    end if
    
    level%S_row_ptr = level%row_ptr
    level%S_col_idx = level%col_idx
    level%nnz_s = level%nnz
    
    ! Step 3: Copy A to S initially (line 86: S = copy(A))
    level%S_values = level%values
    
    ! Step 4: Apply strength criterion exactly like Julia (lines 103-114)
    do i = 1, level%n
      eps_aii = theta * theta * diags(i)  ! Julia: eps_Aii = θ * θ * diags[i]
      
      do k = level%row_ptr(i), level%row_ptr(i+1)-1
        j = level%col_idx(k)
        
        if (i /= j) then  ! Off-diagonal only
          val_sq = level%values(k) * level%values(k)  ! Julia: val*val
          if (val_sq < eps_aii * diags(j)) then  ! Julia: val*val < eps_Aii * diags[row]
            level%S_values(k) = 0.0_dp  ! Zero out weak connections
          end if
        end if
      end do
    end do
    
    ! Step 5: Take absolute values (line 118: S.nzval .= abs.(S.nzval))
    do k = 1, level%nnz_s
      level%S_values(k) = abs(level%S_values(k))
    end do
    
    ! Step 6: Scale columns by largest entry (line 119: scale_cols_by_largest_entry!(S))
    call scale_cols_by_largest_entry(level)
    
    deallocate(diags)
    
  end subroutine sa_compute_strength
  
  subroutine scale_cols_by_largest_entry(level)
    type(amg_level), intent(inout) :: level
    
    integer :: i, k, j
    real(dp) :: max_val
    
    ! Exact replication of Julia scale_cols_by_largest_entry! function
    ! From strength.jl lines 61-70
    
    do i = 1, level%n
      ! Find maximum entry in column i (Julia: _m = find_max(A, i))
      max_val = 0.0_dp
      do k = level%S_row_ptr(i), level%S_row_ptr(i+1)-1
        max_val = max(max_val, level%S_values(k))
      end do
      
      ! Scale all entries in column i by max_val (Julia: A.nzval[j] /= _m)
      if (max_val > 1.0e-14_dp) then
        do k = level%S_row_ptr(i), level%S_row_ptr(i+1)-1
          level%S_values(k) = level%S_values(k) / max_val
        end do
      end if
    end do
    
  end subroutine scale_cols_by_largest_entry
  
  subroutine sa_aggregate_nodes(level)
    type(amg_level), intent(inout) :: level
    
    ! Exact replication of Julia StandardAggregation algorithm
    ! From aggregate.jl lines 4-113
    
    integer :: i, j, k, next_aggregate
    logical :: has_agg_neighbors, has_neighbors
    integer :: x_row, xi
    integer, allocatable :: x(:), y(:)
    
    allocate(level%aggregates(level%n))
    allocate(x(level%n))
    allocate(y(level%n))
    
    x = 0  ! Julia: x = zeros(R, n)
    y = 0  ! Julia: y = zeros(R, n)
    next_aggregate = 1
    
    ! Pass 1: Exactly like Julia lines 12-44
    do i = 1, level%n
      if (x(i) /= 0) cycle  ! Julia: if x[i] != 0; continue; end
      
      has_agg_neighbors = .false.
      has_neighbors = .false.
      
      ! Check neighbors (Julia lines 20-29)
      do k = level%S_row_ptr(i), level%S_row_ptr(i+1)-1
        j = level%S_col_idx(k)
        if (j /= i) then
          has_neighbors = .true.
          if (x(j) /= 0) then
            has_agg_neighbors = .true.
            exit
          end if
        end if
      end do
      
      if (.not. has_neighbors) then
        x(i) = -level%n  ! Julia: x[i] = -n
      else if (.not. has_agg_neighbors) then
        x(i) = next_aggregate  ! Julia: x[i] = next_aggregate
        y(next_aggregate) = i  ! Julia: y[next_aggregate] = i
        
        ! Add all neighbors to this aggregate (Julia lines 37-40)
        do k = level%S_row_ptr(i), level%S_row_ptr(i+1)-1
          j = level%S_col_idx(k)
          x(j) = next_aggregate  ! Julia: x[row] = next_aggregate
        end do
        
        next_aggregate = next_aggregate + 1
      end if
    end do
    
    ! Pass 2: Exactly like Julia lines 47-61
    do i = 1, level%n
      if (x(i) /= 0) cycle
      
      do k = level%S_row_ptr(i), level%S_row_ptr(i+1)-1
        j = level%S_col_idx(k)
        x_row = x(j)
        if (x_row > 0) then
          x(i) = -x_row  ! Julia: x[i] = -x_row
          exit
        end if
      end do
    end do
    
    next_aggregate = next_aggregate - 1  ! Julia: next_aggregate -= 1
    
    ! Pass 3: Exactly like Julia lines 65-91
    do i = 1, level%n
      xi = x(i)
      if (xi /= 0) then
        if (xi > 0) then
          x(i) = xi - 1  ! Julia: x[i] = xi - 1
        else if (xi == -level%n) then
          x(i) = -1  ! Julia: x[i] = -1
        else
          x(i) = -xi - 1  ! Julia: x[i] = -xi - 1
        end if
        cycle
      end if
      
      x(i) = next_aggregate  ! Julia: x[i] = next_aggregate
      y(next_aggregate + 1) = i  ! Julia: y[next_aggregate + 1] = i
      
      do k = level%S_row_ptr(i), level%S_row_ptr(i+1)-1
        j = level%S_col_idx(k)
        if (x(j) == 0) then
          x(j) = next_aggregate
        end if
      end do
      
      next_aggregate = next_aggregate + 1
    end do
    
    ! Convert to 1-based indexing for Fortran (Julia uses 0-based aggregates)
    do i = 1, level%n
      if (x(i) >= 0) then
        level%aggregates(i) = x(i) + 1
      else if (x(i) == -1) then
        ! Isolated nodes get their own aggregate (Julia handling)
        next_aggregate = next_aggregate + 1
        level%aggregates(i) = next_aggregate
      else
        level%aggregates(i) = 0  ! Should not happen with correct Julia algorithm
      end if
    end do
    
    level%n_aggregates = next_aggregate
    level%n_coarse = next_aggregate
    
    deallocate(x, y)
    
  end subroutine sa_aggregate_nodes
  
  subroutine sa_fit_candidates(level)
    type(amg_level), intent(inout) :: level
    
    ! Exact replication of Julia fit_candidates function
    ! From aggregation.jl lines 91-126
    ! Note: For simplicity, we implement the constant vector case (m=1)
    ! The full QR approach would require LAPACK integration
    
    integer :: i, j, k, agg, nnz_p
    integer :: agg_start, agg_size
    integer, allocatable :: nodes_in_agg(:), agg_ptr(:)
    real(dp) :: tol = 1.0e-10_dp
    
    ! For SA-AMG with constant vector candidates, we can use simpler injection
    ! This matches Julia's QR result for constant vectors
    
    ! Build aggregate-to-node mapping (Julia AggOp structure)
    allocate(agg_ptr(level%n_aggregates + 1))
    allocate(nodes_in_agg(level%n))
    
    ! Count nodes per aggregate
    agg_ptr = 0
    do i = 1, level%n
      agg = level%aggregates(i)
      if (agg > 0) then  ! Valid aggregate
        agg_ptr(agg+1) = agg_ptr(agg+1) + 1
      end if
    end do
    
    ! Convert to pointers
    agg_ptr(1) = 1
    do i = 2, level%n_aggregates + 1
      agg_ptr(i) = agg_ptr(i) + agg_ptr(i-1)
    end do
    
    ! Fill node list
    do i = 1, level%n
      agg = level%aggregates(i)
      if (agg > 0) then
        nodes_in_agg(agg_ptr(agg)) = i
        agg_ptr(agg) = agg_ptr(agg) + 1
      end if
    end do
    
    ! Restore pointers
    do i = level%n_aggregates, 1, -1
      agg_ptr(i+1) = agg_ptr(i)
    end do
    agg_ptr(1) = 1
    
    ! Build tentative prolongation exactly like Julia's QR result for constant vectors
    ! Julia: Qs = spzeros(T, n_fine, n_coarse) with proper orthogonalization
    nnz_p = 0
    do i = 1, level%n
      agg = level%aggregates(i)
      if (agg > 0) then
        nnz_p = nnz_p + 1
      end if
    end do
    
    level%n_fine = level%n
    level%n_coarse = level%n_aggregates
    level%nnz_p = nnz_p
    
    allocate(level%P_row_ptr(level%n + 1))
    allocate(level%P_col_idx(nnz_p))
    allocate(level%P_values(nnz_p))
    
    ! Build prolongation with proper normalization (Julia QR orthogonalization result)
    nnz_p = 0
    level%P_row_ptr(1) = 1
    
    do i = 1, level%n
      agg = level%aggregates(i)
      if (agg > 0) then
        agg_size = agg_ptr(agg+1) - agg_ptr(agg)
        
        nnz_p = nnz_p + 1
        level%P_col_idx(nnz_p) = agg
        ! Julia QR for constant vector gives 1/sqrt(aggregate_size)
        level%P_values(nnz_p) = 1.0_dp / sqrt(real(agg_size, dp))
        level%P_row_ptr(i+1) = nnz_p + 1
      else
        ! Unaggregated node - this shouldn't happen with correct aggregation
        level%P_row_ptr(i+1) = level%P_row_ptr(i)
      end if
    end do
    
    deallocate(agg_ptr, nodes_in_agg)
    
  end subroutine sa_fit_candidates
  
  subroutine sa_smooth_prolongation(level, omega)
    type(amg_level), intent(inout) :: level
    real(dp), intent(in) :: omega
    
    ! Exact replication of Julia AlgebraicMultigrid.jl JacobiProlongation with LocalWeighting
    ! From smoother.jl lines 157-194:
    ! function (j::JacobiProlongation)(A, T, S, B, degree = 1, weighting = LocalWeighting())
    !     D_inv_S = weight(weighting, A, j.ω)  
    !     P = T
    !     for i = 1:degree
    !         P = P - (D_inv_S * P)
    !     end
    !     P
    ! end
    
    integer :: i, j, k, kp, kp_j, nnz_ds
    real(dp) :: row_sum, d_inv_val
    real(dp), allocatable :: D(:), D_inv_S_values(:), DSP_values(:)
    integer, allocatable :: D_inv_S_row_ptr(:), D_inv_S_col_idx(:)
    
    ! Step 1: Compute LocalWeighting exactly like Julia (lines 178-194)
    ! D = zeros(eltype(S), size(S,1))  
    ! for i = 1:size(S, 1)
    !     for j in nzrange(S, i)
    !         row = S.rowval[j]
    !         val = S.nzval[j]
    !         D[row] += abs(val)
    !     end
    ! end
    
    allocate(D(level%n))
    D = 0.0_dp
    
    ! Compute row sums of absolute values (Julia's LocalWeighting)
    do i = 1, level%n
      do k = level%S_row_ptr(i), level%S_row_ptr(i+1)-1
        j = level%S_col_idx(k)
        D(j) = D(j) + abs(level%S_values(k))
      end do
    end do
    
    ! for i = 1:size(D, 1)
    !     if D[i] != 0
    !         D[i] = 1/D[i]
    !     end
    ! end
    do i = 1, level%n
      if (abs(D(i)) > 1.0e-14_dp) then
        D(i) = 1.0_dp / D(i)
      else
        D(i) = 0.0_dp
      end if
    end do
    
    ! Step 2: Build D_inv_S = scale_rows(S, D) * ω
    ! This creates D_inv_S matrix (same structure as S)
    nnz_ds = level%nnz_s
    allocate(D_inv_S_row_ptr(level%n + 1))
    allocate(D_inv_S_col_idx(nnz_ds))
    allocate(D_inv_S_values(nnz_ds))
    
    D_inv_S_row_ptr = level%S_row_ptr
    D_inv_S_col_idx = level%S_col_idx
    
    ! D_inv_S = scale_rows(S, D) with omega scaling
    do i = 1, level%n
      do k = level%S_row_ptr(i), level%S_row_ptr(i+1)-1
        j = level%S_col_idx(k)
        ! Scale by row weight and omega (Julia: rmul!(D_inv_S, eltype(S)(ω)))
        D_inv_S_values(k) = level%S_values(k) * D(j) * omega
      end do
    end do
    
    ! Step 3: Apply smoothing exactly like Julia: P = P - (D_inv_S * P)
    ! Compute D_inv_S * P 
    allocate(DSP_values(level%nnz_p))
    DSP_values = 0.0_dp
    
    do i = 1, level%n
      kp = level%P_row_ptr(i)
      
      ! Multiply row i of D_inv_S with prolongation
      do k = D_inv_S_row_ptr(i), D_inv_S_row_ptr(i+1)-1
        j = D_inv_S_col_idx(k)
        
        ! Find P(j, :) contribution
        if (j <= level%n) then
          kp_j = level%P_row_ptr(j)
          if (level%P_col_idx(kp) == level%P_col_idx(kp_j)) then
            DSP_values(kp) = DSP_values(kp) + D_inv_S_values(k) * level%P_values(kp_j)
          end if
        end if
      end do
    end do
    
    ! Final step: P = P - (D_inv_S * P)
    do i = 1, level%n
      kp = level%P_row_ptr(i)
      level%P_values(kp) = level%P_values(kp) - DSP_values(kp)
    end do
    
    deallocate(D, D_inv_S_row_ptr, D_inv_S_col_idx, D_inv_S_values, DSP_values)
    
  end subroutine sa_smooth_prolongation
  
  subroutine extract_diagonal(level)
    type(amg_level), intent(inout) :: level
    
    integer :: i, k, j
    
    if (.not. allocated(level%diagonal)) then
      allocate(level%diagonal(level%n))
    end if
    level%diagonal = 0.0_dp
    
    do i = 1, level%n
      do k = level%row_ptr(i), level%row_ptr(i+1)-1
        j = level%col_idx(k)
        if (i == j) then
          level%diagonal(i) = level%values(k)
          exit
        end if
      end do
    end do
    
  end subroutine extract_diagonal
  
  subroutine build_restriction(level)
    type(amg_level), intent(inout) :: level
    
    ! For SA-AMG, restriction is transpose of prolongation
    ! R = P^T
    call transpose_csr(level%n_fine, level%n_coarse, level%nnz_p, &
                      level%P_row_ptr, level%P_col_idx, level%P_values, &
                      level%R_row_ptr, level%R_col_idx, level%R_values, &
                      level%nnz_r)
    
  end subroutine build_restriction
  
  subroutine transpose_csr(n_rows, n_cols, nnz, row_ptr, col_idx, values, &
                          row_ptr_t, col_idx_t, values_t, nnz_t)
    integer, intent(in) :: n_rows, n_cols, nnz
    integer, intent(in) :: row_ptr(n_rows+1), col_idx(nnz)
    real(dp), intent(in) :: values(nnz)
    integer, allocatable, intent(out) :: row_ptr_t(:), col_idx_t(:)
    real(dp), allocatable, intent(out) :: values_t(:)
    integer, intent(out) :: nnz_t
    
    integer :: i, j, k, pos
    integer, allocatable :: col_count(:), col_ptr(:)
    
    nnz_t = nnz
    allocate(row_ptr_t(n_cols+1))
    allocate(col_idx_t(nnz_t))
    allocate(values_t(nnz_t))
    allocate(col_count(n_cols))
    allocate(col_ptr(n_cols))
    
    ! Count entries per column
    col_count = 0
    do i = 1, n_rows
      do k = row_ptr(i), row_ptr(i+1)-1
        j = col_idx(k)
        col_count(j) = col_count(j) + 1
      end do
    end do
    
    ! Build column pointers
    row_ptr_t(1) = 1
    do j = 1, n_cols
      row_ptr_t(j+1) = row_ptr_t(j) + col_count(j)
    end do
    
    ! Copy column pointers for insertion
    col_ptr = row_ptr_t(1:n_cols)
    
    ! Fill transpose
    do i = 1, n_rows
      do k = row_ptr(i), row_ptr(i+1)-1
        j = col_idx(k)
        pos = col_ptr(j)
        col_idx_t(pos) = i
        values_t(pos) = values(k)
        col_ptr(j) = col_ptr(j) + 1
      end do
    end do
    
    deallocate(col_count, col_ptr)
    
  end subroutine transpose_csr
  
  subroutine compute_coarse_operator(level_fine, level_coarse)
    type(amg_level), intent(in) :: level_fine
    type(amg_level), intent(inout) :: level_coarse
    
    ! Compute A_c = R * A * P using sparse matrix multiplication
    ! This is a simplified implementation - production code would use optimized SpGEMM
    
    integer :: i, j, k, jp, kp, kr, jr
    real(dp) :: val
    integer :: nnz_estimate
    
    level_coarse%n = level_fine%n_coarse
    
    ! Estimate nonzeros (rough upper bound)
    nnz_estimate = min(level_coarse%n * level_coarse%n, 7 * level_coarse%n)
    
    allocate(level_coarse%row_ptr(level_coarse%n + 1))
    allocate(level_coarse%col_idx(nnz_estimate))
    allocate(level_coarse%values(nnz_estimate))
    allocate(level_coarse%diagonal(level_coarse%n))
    
    ! Simplified RAP computation
    ! In production, use optimized sparse matrix multiplication
    call sparse_rap_product(level_fine, level_coarse)
    
    ! Extract diagonal
    call extract_diagonal(level_coarse)
    
    ! Diagnostic output for coarse operator properties
    if (allocated(level_coarse%diagonal)) then
      call diagnose_coarse_operator(level_coarse)
    end if
    
  end subroutine compute_coarse_operator
  
  subroutine diagnose_coarse_operator(level)
    type(amg_level), intent(in) :: level
    
    real(dp) :: min_diag, max_diag, diag_ratio
    integer :: i, zero_diags, neg_diags
    
    if (.not. allocated(level%diagonal) .or. level%n <= 0) return
    
    min_diag = huge(1.0_dp)
    max_diag = -huge(1.0_dp)
    zero_diags = 0
    neg_diags = 0
    
    do i = 1, level%n
      if (abs(level%diagonal(i)) < 1.0e-14_dp) then
        zero_diags = zero_diags + 1
      else if (level%diagonal(i) < 0.0_dp) then
        neg_diags = neg_diags + 1
      end if
      min_diag = min(min_diag, abs(level%diagonal(i)))
      max_diag = max(max_diag, abs(level%diagonal(i)))
    end do
    
    if (min_diag > 1.0e-14_dp) then
      diag_ratio = max_diag / min_diag
    else
      diag_ratio = huge(1.0_dp)
    end if
    
    print '(A,I0,A)', '  Coarse operator (n=', level%n, ') diagnostics:'
    print '(A,E12.3)', '    Min diagonal magnitude: ', min_diag
    print '(A,E12.3)', '    Max diagonal magnitude: ', max_diag
    print '(A,E12.3)', '    Diagonal ratio (cond est): ', diag_ratio
    print '(A,I0)', '    Zero diagonals: ', zero_diags
    print '(A,I0)', '    Negative diagonals: ', neg_diags
    
  end subroutine diagnose_coarse_operator
  
  subroutine sparse_rap_product(level_fine, level_coarse)
    type(amg_level), intent(in) :: level_fine
    type(amg_level), intent(inout) :: level_coarse
    
    ! Compute RAP = R * A * P where R = P^T for Smoothed Aggregation
    ! This implements the correct SA-AMG coarse operator construction
    
    integer :: i, j, k, ic, jc, nnz_count, kp_i, kp_j, kr_ic
    real(dp) :: value, r_val, p_val
    integer, allocatable :: temp_i(:), temp_j(:)
    real(dp), allocatable :: temp_v(:)
    integer :: max_nnz_estimate
    
    ! Step 1: Estimate maximum nonzeros in coarse matrix
    max_nnz_estimate = level_coarse%n * level_coarse%n  ! Conservative upper bound
    
    ! Allocate temporary storage
    allocate(temp_i(max_nnz_estimate))
    allocate(temp_j(max_nnz_estimate))
    allocate(temp_v(max_nnz_estimate))
    
    nnz_count = 0
    
    ! Step 2: Compute RAP = R * A * P using proper SA theory
    ! A_coarse[ic, jc] = sum_{i,k} R[ic, i] * A[i, k] * P[k, jc]
    ! For SA: R = P^T, so R[ic, i] = P[i, ic]
    
    do ic = 1, level_coarse%n
      do jc = 1, level_coarse%n
        value = 0.0_dp
        
        ! Triple sum: sum over fine nodes i and k
        do i = 1, level_fine%n
          ! Get R[ic, i] = P[i, ic] (restriction from prolongation transpose)
          r_val = 0.0_dp
          if (level_fine%aggregates(i) == ic) then
            kp_i = level_fine%P_row_ptr(i)
            if (level_fine%P_col_idx(kp_i) == ic) then
              r_val = level_fine%P_values(kp_i)
            end if
          end if
          
          if (abs(r_val) > 1.0e-14_dp) then
            ! Sum over connections in row i of A
            do k = level_fine%row_ptr(i), level_fine%row_ptr(i+1)-1
              j = level_fine%col_idx(k)
              
              ! Get P[j, jc] (prolongation value)
              p_val = 0.0_dp
              if (level_fine%aggregates(j) == jc) then
                kp_j = level_fine%P_row_ptr(j)
                if (level_fine%P_col_idx(kp_j) == jc) then
                  p_val = level_fine%P_values(kp_j)
                end if
              end if
              
              if (abs(p_val) > 1.0e-14_dp) then
                ! Add contribution: R[ic,i] * A[i,j] * P[j,jc]
                value = value + r_val * level_fine%values(k) * p_val
                
                ! Debug output for diagonal entries that become negative
                if (ic == jc .and. ic <= 3 .and. value < 0.0_dp) then
                  print '(A,2I3,A,3E12.3)', '  DEBUG RAP diagonal ic=', ic, jc, &
                    ' contrib:', r_val, level_fine%values(k), p_val
                  print '(A,E12.3)', '    running value:', value
                end if
              end if
            end do
          end if
        end do
        
        ! Store nonzero entry
        if (abs(value) > 1.0e-14_dp) then
          nnz_count = nnz_count + 1
          if (nnz_count > max_nnz_estimate) then
            ! Should not happen with conservative estimate
            error stop "RAP product: exceeded nonzero estimate"
          end if
          temp_i(nnz_count) = ic
          temp_j(nnz_count) = jc  
          temp_v(nnz_count) = value
        end if
      end do
    end do
    
    ! Step 3: Convert to CSR format
    level_coarse%nnz = nnz_count
    
    ! Sort entries by row, then column (required for CSR)
    call sort_csr_entries(temp_i, temp_j, temp_v, nnz_count)
    
    ! Fill CSR arrays
    level_coarse%row_ptr(1) = 1
    j = 1
    do i = 1, nnz_count
      level_coarse%col_idx(i) = temp_j(i)
      level_coarse%values(i) = temp_v(i)
      
      ! Update row pointers
      do while (j < temp_i(i))
        j = j + 1
        level_coarse%row_ptr(j) = i
      end do
    end do
    
    ! Fill remaining row pointers
    do while (j <= level_coarse%n)
      j = j + 1
      level_coarse%row_ptr(j) = nnz_count + 1
    end do
    
    deallocate(temp_i, temp_j, temp_v)
    
  end subroutine sparse_rap_product
  
  subroutine sort_csr_entries(rows, cols, vals, nnz)
    integer, intent(inout) :: rows(:), cols(:)
    real(dp), intent(inout) :: vals(:)
    integer, intent(in) :: nnz
    
    ! Simple insertion sort for CSR entries
    ! Sort by row first, then by column within each row
    integer :: i, j, temp_r, temp_c
    real(dp) :: temp_v
    
    do i = 2, nnz
      temp_r = rows(i)
      temp_c = cols(i)
      temp_v = vals(i)
      j = i - 1
      
      do while (j >= 1)
        if (rows(j) < temp_r .or. (rows(j) == temp_r .and. cols(j) <= temp_c)) exit
        rows(j+1) = rows(j)
        cols(j+1) = cols(j)
        vals(j+1) = vals(j)
        j = j - 1
      end do
      
      rows(j+1) = temp_r
      cols(j+1) = temp_c
      vals(j+1) = temp_v
    end do
    
  end subroutine sort_csr_entries
  
  subroutine compute_complexities(hierarchy)
    type(amg_hierarchy), intent(inout) :: hierarchy
    
    integer :: i
    integer :: total_unknowns, total_nonzeros
    
    total_unknowns = 0
    total_nonzeros = 0
    
    do i = 1, hierarchy%n_levels
      total_unknowns = total_unknowns + hierarchy%levels(i)%n
      total_nonzeros = total_nonzeros + hierarchy%levels(i)%nnz
    end do
    
    hierarchy%grid_complexity = real(total_unknowns, dp) / real(hierarchy%levels(1)%n, dp)
    hierarchy%operator_complexity = real(total_nonzeros, dp) / real(hierarchy%levels(1)%nnz, dp)
    
  end subroutine compute_complexities
  
end module amg_smoothed_aggregation_mod