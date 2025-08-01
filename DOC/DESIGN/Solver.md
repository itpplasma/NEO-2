# NEO-2 Solver Architecture and GMRES Implementation Design

## Executive Summary

NEO-2's current solver architecture relies heavily on direct sparse solvers (UMFPACK) which creates significant memory bottlenecks, particularly for the stellarator (PAR) variant. The sparse operator inversion hits memory bounds at ~6-7 GB per MPI process, limiting scalability. This document analyzes the current solver framework and proposes GMRES integration to overcome memory limitations.

## Current Solver Architecture

### Core Components

#### 1. Sparse Matrix Infrastructure
- **Primary Module**: `COMMON/sparse_mod.f90`
  - Generic interface `sparse_solve` with multiple backend support
  - Storage: COO (assembly) → CSC (solving) format conversion
  - Current backends: SuiteSparse UMFPACK (methods 2,3)
  - Matrix size validation and runtime bounds checking

#### 2. Ripple Solver Framework
The "ripple solver" refers to the particle transport calculations in magnetic ripples:

**NEO-2-QL** (`ripple_solver_axi_test.f90`, `ripple_solver_ArnoldiOrder2_test.f90`):
- Axisymmetric tokamak geometry
- Arnoldi method for eigenvalue analysis
- Richardson iteration with eigenmode subtraction
- Smaller memory footprint due to 2D axisymmetry

**NEO-2-PAR** (`NEO-2-PAR/ripple_solver.f90`):
- Full 3D stellarator magnetic geometry
- MPI parallelization across field lines
- Direct sparse solver approach only
- Memory scaling: O(nsurf³) for surface discretization

#### 3. Linear Algebra Integration
- **Arnoldi Module** (`COMMON/arnoldi_mod.f90`): Krylov subspace eigenvalue solver (QL only)
- **LAPACK Interface** (`COMMON/lapack_band.f90`): Dense system wrapper
- **Simple Solver** (`COMMON/solve_system.f90`): Basic LAPACK integration

## Memory Bottlenecks Analysis

### 1. Sparse Matrix Assembly
```fortran
! From splinecof3_direct_sparse.f90
! Matrix size scales as O(nsurf³)
allocate(coo_row(max_nnz), coo_col(max_nnz), coo_val(max_nnz))
```
- **Issue**: Memory allocation failures in `splinecof3_direct_sparse.f90`
- **Scaling**: Cubic growth with surface resolution
- **Impact**: Limits problem size for magnetic field representation

### 2. UMFPACK Factorization
```fortran
! Current sparse_solve method = 2,3
call umf4zfac(symbolic, numeric, ...)  ! LU factorization
```
- **Memory Usage**: 5-10x matrix storage for factorization
- **PAR Bottleneck**: 6-7 GB per MPI process documented limit
- **Constraint**: Limits to ~48 threads per node

### 3. Matrix Storage Format
- **CSC Storage**: `(nnz + n + 1)` memory for values + pointers
- **Problem**: No memory reuse between different RHS vectors
- **3D Complexity**: Stellarator geometry creates much larger systems than tokamaks

## QL vs PAR Solver Differences

### Algorithmic Differences

| Aspect | NEO-2-QL | NEO-2-PAR |
|--------|----------|-----------|
| **Geometry** | Axisymmetric (2D) | Full 3D stellarator |
| **Parallelization** | OpenMP threads | MPI + OpenMP hybrid |
| **Eigenvalue Analysis** | Arnoldi + Richardson | Direct solve only |
| **Memory Scaling** | O(nsurf²) | O(nsurf³) |
| **Stability Method** | Unstable eigenmode subtraction | Full system inversion |

### Implementation Differences

**NEO-2-QL Approach**:
```fortran
! Arnoldi method for unstable eigenmodes
call arnoldi_iteration(matrix, rhs, eigenvals, eigenvecs)
! Richardson iteration with preconditioning
call richardson_preconditioned(matrix, rhs, solution, eigenvecs)
```

**NEO-2-PAR Approach**:
```fortran
! Direct sparse solver only
call sparse_solve(matrix, rhs, solution, method=3)  ! UMFPACK
```

### Computational Complexity

**QL Complexity**:
- Matrix assembly: O(nsurf² × nspecies)
- Arnoldi method: O(m × n²) for m iterations
- Richardson iteration: O(k × n²) for k iterations
- **Total**: O(nsurf² × (m + k))

**PAR Complexity**:
- Matrix assembly: O(nsurf³ × nspecies × nproc)
- UMFPACK factorization: O(n³) sparse fill-in dependent
- Multiple RHS solves: O(n²) per RHS
- **Total**: O(nsurf³ × fill-in factor)

## GMRES Implementation Design

### 1. Integration Point
**Primary Interface**: Extend `sparse_mod.f90`
```fortran
! Add new solver method
sparse_solve_method = 4  ! GMRES iterative solver

! Generic interface remains unchanged
call sparse_solve(matrix, rhs, solution, method=4)
```

### 2. GMRES Module Structure
**New Module**: `COMMON/gmres_mod.f90`
```fortran
module gmres_mod
  use iso_fortran_env
  use sparse_mod
  
  implicit none
  
  type :: gmres_solver
    integer :: restart_dim = 30     ! Krylov subspace dimension
    real(wp) :: tolerance = 1.0e-12 ! Convergence tolerance  
    integer :: max_iter = 1000      ! Maximum iterations
    logical :: use_precon = .true.  ! Preconditioning flag
  end type
  
contains
  
  subroutine gmres_solve(matrix, rhs, solution, solver_params, info)
    ! Main GMRES implementation
  end subroutine
  
  subroutine arnoldi_process(matrix, krylov_basis, hessenberg)
    ! Arnoldi iteration for orthogonal basis construction
  end subroutine
  
  subroutine apply_givens_rotations(hessenberg, rhs_vec)
    ! QR factorization via Givens rotations
  end subroutine
  
end module gmres_mod
```

### 3. Memory Comparison

**Current UMFPACK**:
- Matrix storage: `nnz × (8 + 4 + 4)` bytes (value + row + col)
- Factorization: `~5-10 × nnz × 8` bytes
- **Total**: ~50-80 GB for large PAR problems

**Proposed GMRES**:
- Matrix storage: `nnz × 16` bytes (value + indices)  
- Krylov basis: `m × n × 8` bytes (m = restart dimension)
- Hessenberg matrix: `m² × 8` bytes
- **Total**: ~5-10 GB for same problems (5-8x reduction)

### 4. Preconditioning Strategy

**Option 1: ILU(k) Preconditioning**
```fortran
! Use UMFPACK symbolic factorization for structure
call umf4zsym(symbolic_handle, matrix)
! Generate incomplete factorization  
call ilu_factorize(matrix, symbolic_handle, ilu_factors, fill_level=2)
```

**Option 2: Physics-Based Preconditioning**
```fortran
! Diagonal collision operator preconditioning
call extract_collision_diagonal(matrix, diag_preconditioner)
! Block-diagonal magnetic surface preconditioning  
call build_surface_blocks(matrix, surface_indices, block_preconditioner)
```

**Option 3: Multigrid Preconditioning**
```fortran
! Magnetic surface hierarchy for multigrid
call build_surface_hierarchy(fine_surfaces, coarse_surfaces, restriction_op)
call multigrid_vcycle(fine_matrix, coarse_matrix, rhs, correction)
```

### 5. Algorithm Implementation

**Core GMRES Algorithm**:
```fortran
subroutine gmres_solve(A, b, x, params, info)
  ! Input: matrix A, RHS b, initial guess x
  ! Output: solution x, convergence info
  
  ! Initialize
  r0 = b - A*x                    ! Initial residual
  beta = ||r0||                   ! Residual norm
  V(:,1) = r0/beta               ! First Krylov vector
  
  do restart = 1, max_restarts
    
    ! Arnoldi process  
    do j = 1, restart_dim
      w = A * V(:,j)              ! Matrix-vector product
      
      ! Gram-Schmidt orthogonalization
      do i = 1, j
        H(i,j) = dot_product(w, V(:,i))
        w = w - H(i,j) * V(:,i)
      end do
      
      H(j+1,j) = ||w||
      if (H(j+1,j) < tolerance) exit  ! Lucky breakdown
      V(:,j+1) = w / H(j+1,j)
    end do
    
    ! Solve least squares problem: min ||H*y - beta*e1||
    call qr_solve_hessenberg(H, beta, y)
    
    ! Update solution
    x = x + V(:,1:j) * y
    
    ! Check convergence  
    residual_norm = abs(H(j+1,j) * y(j))
    if (residual_norm < tolerance) exit
    
  end do
  
end subroutine
```

## Implementation Phases

### Phase 1: Basic GMRES Implementation
- [ ] Create `gmres_mod.f90` with core algorithm
- [ ] Integrate with `sparse_mod.f90` as method 4
- [ ] Basic convergence testing with simple problems
- [ ] Memory usage validation

### Phase 2: Preconditioning Integration  
- [ ] Implement ILU(k) preconditioning
- [ ] Physics-based diagonal preconditioning
- [ ] Benchmarking vs UMFPACK on medium problems

### Phase 3: PAR Integration and Optimization
- [ ] MPI-aware GMRES for distributed matrices
- [ ] Memory profiling and optimization
- [ ] Large-scale stellarator problem testing
- [ ] Performance comparison with current solver

### Phase 4: Advanced Features
- [ ] Multigrid preconditioning
- [ ] Adaptive restart strategies
- [ ] Integration with existing Arnoldi eigenvalue analysis
- [ ] Hybrid direct/iterative approach

## Benefits and Impact

### Memory Reduction
- **5-8x memory reduction** for PAR problems
- Enable larger stellarator configurations
- More MPI processes per node (better scalability)

### Algorithmic Advantages
- **Matrix-free operation**: Only need matrix-vector products
- **Restart capability**: Handle memory constraints gracefully  
- **Preconditioning flexibility**: Physics-aware acceleration
- **Convergence control**: Adaptive tolerance and restart

### Integration Benefits
- **Backward compatibility**: Existing solver methods unchanged
- **Runtime selection**: Easy switching between direct/iterative
- **Hybrid approach**: Combine with existing Arnoldi framework
- **Cross-platform**: Pure Fortran implementation

## Risk Mitigation

### Convergence Issues
- **Risk**: GMRES may not converge for ill-conditioned systems
- **Mitigation**: Robust preconditioning + fallback to UMFPACK

### Numerical Accuracy
- **Risk**: Iterative solver precision vs direct solver accuracy  
- **Mitigation**: Adaptive tolerance + residual monitoring

### Integration Complexity
- **Risk**: Disruption of existing validated solver chain
- **Mitigation**: Additive implementation (new method=4) with extensive testing

## Success Metrics

1. **Memory Usage**: <2 GB per MPI process for problems requiring 6-7 GB with UMFPACK
2. **Convergence**: <1000 iterations for 1e-12 tolerance on typical problems
3. **Accuracy**: Solutions match UMFPACK to within 1e-10 relative error
4. **Performance**: Total solve time competitive with UMFPACK for large problems
5. **Scalability**: Enable 2-3x more MPI processes per node for PAR problems

This design provides a comprehensive roadmap for implementing GMRES as a memory-efficient alternative to the current direct sparse solvers, with particular focus on addressing the memory bottlenecks in stellarator (PAR) calculations.