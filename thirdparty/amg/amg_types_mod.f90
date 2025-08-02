module amg_types_mod
  use sparse_types_mod, only: dp
  implicit none
  
  private
  public :: amg_hierarchy, amg_level, amg_params
  public :: AMG_CLASSICAL, AMG_SMOOTHED_AGGREGATION
  public :: AMG_VCYCLE, AMG_WCYCLE, AMG_FCYCLE
  public :: AMG_JACOBI, AMG_GAUSS_SEIDEL, AMG_SYM_GAUSS_SEIDEL
  
  ! AMG method constants
  integer, parameter :: AMG_CLASSICAL = 1
  integer, parameter :: AMG_SMOOTHED_AGGREGATION = 2
  
  ! Cycle type constants
  integer, parameter :: AMG_VCYCLE = 1
  integer, parameter :: AMG_WCYCLE = 2
  integer, parameter :: AMG_FCYCLE = 3
  
  ! Smoother constants
  integer, parameter :: AMG_JACOBI = 1
  integer, parameter :: AMG_GAUSS_SEIDEL = 2
  integer, parameter :: AMG_SYM_GAUSS_SEIDEL = 3
  
  ! AMG parameters
  type :: amg_params
    ! Method selection
    integer :: method = AMG_SMOOTHED_AGGREGATION
    integer :: cycle_type = AMG_VCYCLE
    integer :: smoother = AMG_SYM_GAUSS_SEIDEL
    
    ! Coarsening parameters
    real(dp) :: strength_threshold = 0.0_dp  ! 0.0 for SA, 0.25 for Classical
    integer :: max_levels = 20
    integer :: coarsest_size = 50
    real(dp) :: max_coarse_ratio = 0.8_dp
    
    ! Smoothing parameters
    integer :: presmoother_steps = 2
    integer :: postsmoother_steps = 2
    real(dp) :: jacobi_weight = 0.66666666666666666_dp  ! 2/3
    
    ! SA-specific parameters
    real(dp) :: prolongation_damping = 1.33333333333333333_dp  ! 4/3
    
    ! Convergence parameters
    real(dp) :: tolerance = 1.0e-10_dp
    integer :: max_iter = 100
    logical :: verbose = .false.
  end type amg_params
  
  ! Single AMG level
  type :: amg_level
    ! Sparse matrices in CSR format
    integer :: n = 0                     ! Number of rows/cols
    integer :: nnz = 0                   ! Number of nonzeros
    integer, allocatable :: row_ptr(:)   ! CSR row pointers (n+1)
    integer, allocatable :: col_idx(:)   ! CSR column indices (nnz)
    real(dp), allocatable :: values(:)   ! CSR values (nnz)
    
    ! Prolongation operator (fine-to-coarse)
    integer :: n_fine = 0
    integer :: n_coarse = 0
    integer :: nnz_p = 0
    integer, allocatable :: P_row_ptr(:)
    integer, allocatable :: P_col_idx(:)
    real(dp), allocatable :: P_values(:)
    
    ! Restriction operator (coarse-to-fine, often P^T)
    integer :: nnz_r = 0
    integer, allocatable :: R_row_ptr(:)
    integer, allocatable :: R_col_idx(:)
    real(dp), allocatable :: R_values(:)
    
    ! Strength matrix (for SA-AMG prolongation smoothing)
    integer :: nnz_s = 0
    integer, allocatable :: S_row_ptr(:)
    integer, allocatable :: S_col_idx(:)
    real(dp), allocatable :: S_values(:)
    
    ! Work arrays for smoothing
    real(dp), allocatable :: diagonal(:)  ! Diagonal of A for Jacobi/GS
    real(dp), allocatable :: res_work(:)  ! Residual workspace
    real(dp), allocatable :: sol_work(:)  ! Solution workspace
    
    ! Aggregation info (for SA-AMG)
    integer, allocatable :: aggregates(:) ! Node-to-aggregate mapping
    integer :: n_aggregates = 0
  end type amg_level
  
  ! Complete AMG hierarchy
  type :: amg_hierarchy
    type(amg_params) :: params
    type(amg_level), allocatable :: levels(:)
    integer :: n_levels = 0
    logical :: setup_done = .false.
    
    ! Statistics
    real(dp) :: grid_complexity = 0.0_dp
    real(dp) :: operator_complexity = 0.0_dp
    integer :: setup_iter = 0
    real(dp) :: setup_time = 0.0_dp
  end type amg_hierarchy
  
end module amg_types_mod