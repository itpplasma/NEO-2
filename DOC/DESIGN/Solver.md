# NEO-2 Solver Architecture and GMRES Implementation Design

## Executive Summary

NEO-2's current solver architecture relies heavily on direct sparse solvers (UMFPACK) which creates significant memory bottlenecks for both tokamak (QL) and stellarator (PAR) variants. The primary bottleneck is **O(lag³) scaling from collision operators** in velocity space, where `lag` controls energy/speed discretization and `leg` controls pitch angle discretization. Both variants hit memory limits when increasing velocity space resolution, limiting physics accuracy. This document analyzes the velocity space basis functions, operator structure, and proposes GMRES integration to break the cubic scaling limitation.

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

## Velocity Space Basis Functions and Memory Scaling

### Velocity Space Discretization

#### Energy/Speed Direction (`lag` parameter)
- **Physical quantity**: v/v_th (normalized velocity magnitude)  
- **Default values**: lag=10 (QL), lag=3 (PAR testing)
- **Basis function options**:
  - **Laguerre basis** (`collop_base_prj = 0`): Generalized Laguerre L^(3/2)_m(x²) [DEFAULT]
  - **B-spline basis** (`collop_base_prj = 11`): B-splines with order `collop_bspline_order` [RECOMMENDED]
  - **Polynomial basis** (`collop_base_prj = 1,2`): Standard or quadratic polynomials
  - **Cubic splines** (`collop_base_prj = 10`): Traditional cubic splines

#### Pitch Angle Direction (`leg` parameter)  
- **Physical quantity**: ξ = v_∥/v (pitch angle cosine)
- **Default values**: leg=20 (QL), leg=3 (PAR testing)
- **Basis functions**: Legendre polynomials P_l(ξ) (fixed choice for physics reasons)

### Memory Bottlenecks Analysis

#### 1. Collision Operator Assembly (Primary Bottleneck)
```fortran
! From collision_operator_mems.f90 - critical allocations:
allocate(anumm_aa(0:lag,0:lag,0:num_spec-1,0:num_spec-1))          ! O(lag² × species²)
allocate(ailmm_aa(0:lag,0:lag,0:leg,0:num_spec-1,0:num_spec-1))    ! O(lag² × leg × species²)

! Rosenbluth potential integrals
allocate(I1_mmp_s(0:lagmax, 0:lagmax, 0:legmax_local))             ! O(lag² × leg)
allocate(I2_mmp_s(0:lagmax, 0:lagmax, 0:legmax_local))             ! O(lag² × leg)
allocate(I3_mmp_s(0:lagmax, 0:lagmax, 0:legmax_local))             ! O(lag² × leg)
allocate(I4_mmp_s(0:lagmax, 0:lagmax, 0:legmax_local))             ! O(lag² × leg)
```

#### 2. Sparse Matrix Scaling
```fortran
! Main system matrix size from field line discretization:
n_2d_size = sum over field line steps: 2*(lag+1)*(npl(istep)+1)

! Matrix elements scale as: nnz ∝ lag³ × nsteps
! Documentation: "Memory scales at least with cube of lag parameter"
```

**Complete Memory Scaling**:
- **Collision operators**: O(lag² × leg × nspecies²)
- **Distribution functions**: O(lag × 4 × npart × nperiods)  
- **Sparse matrix storage**: O(lag³ × nsteps) 
- **UMFPACK factorization**: O(lag³ × nsteps × fill_factor), fill_factor ~5-10x

### 2. UMFPACK Factorization
```fortran
! Current sparse_solve method = 2,3
call umf4zfac(symbolic, numeric, ...)  ! LU factorization
```
- **Memory Usage**: 5-10x matrix storage for factorization
- **PAR Bottleneck**: 6-7 GB per MPI process documented limit
- **Constraint**: Limits to ~48 threads per node

### 3. Matrix Storage Format
- **CSC Storage**: `(nnz + n + 1)` memory for values + pointers where nnz ∝ lag³ × nsteps
- **Problem**: No memory reuse between different RHS vectors
- **Factorization dominance**: UMFPACK memory ~5-10x matrix storage

## Operator Structure and Discretization

### Differential Operators

#### Field Line Transport
```fortran
! Parallel derivative along field lines
∂f/∂s where s = arc length parameter
! Discretized via finite differences on field line grid
```
- **Magnetic drifts**: ∇B and curvature drift operators
- **Bounce/transit**: Particle classification at turning points
- **Field line following**: RK4 integration of dx/ds = B/|B|

#### Velocity Space Derivatives  
```fortran
! Collision operator derivatives in velocity space
∂/∂v (energy direction) and ∂/∂ξ (pitch angle direction)
! Discretized via basis function derivatives
```
- **Energy direction**: Laguerre/B-spline basis derivatives
- **Pitch angle**: Legendre polynomial derivatives  
- **Mixed derivatives**: Cross-terms ∂²/(∂v∂ξ)

### Integral Operators

#### Collision Integrals (Rosenbluth Potentials)
```fortran
! From collop_compute.f90 - precomputed integral matrices:
I1_mmp_s: Momentum exchange integrals    ∫∫ G(v,v') dv' terms
I2_mmp_s: Energy exchange integrals      ∫∫ H(v,v') dv' terms  
I3_mmp_s: Mixed momentum-energy terms    ∫∫ mixed dv' terms
I4_mmp_s: Angular scattering integrals   ∫∫ angular dv' terms
```
- **Rosenbluth potentials**: G(v,v') and H(v,v') computed via GSL adaptive quadrature
- **Precomputation**: Matrix elements stored and reused when possible
- **Numerical integration**: High-precision quadrature with error control

#### Field Line Integrals
```fortran
! Bounce/transit averaging integrals
∫ ds/v_∥ over particle orbits
! Numerical integration along field lines
```
- **Orbit integration**: Particle trajectory following
- **Magnetic well resolution**: Binary splitting for ripple structure
- **Boundary conditions**: Periodic/reflecting boundary treatment

## QL vs PAR Detailed Comparison

### Velocity Space Treatment (IDENTICAL)

| Component | QL Implementation | PAR Implementation |
|-----------|------------------|-------------------|
| **Energy basis** | Same options: Laguerre/B-spline/polynomial | Same options: Laguerre/B-spline/polynomial |
| **Pitch angle basis** | Legendre polynomials P_l(ξ) | Legendre polynomials P_l(ξ) |
| **Collision operators** | Same I1-I4 Rosenbluth integral assembly | Same I1-I4 Rosenbluth integral assembly |
| **Basis coefficients** | O(lag² × leg × nspecies²) | O(lag² × leg × nspecies²) |
| **Memory per collision op** | Same scaling with lag³ | Same scaling with lag³ |

### Field Line and Matrix Differences

| Aspect | NEO-2-QL | NEO-2-PAR |
|--------|----------|-----------|
| **Magnetic geometry** | Axisymmetric (2D) | Full 3D stellarator |
| **Field line complexity** | ~1-4 periods | ~100+ periods typical |
| **Parallelization** | OpenMP threads only | MPI + OpenMP hybrid |
| **Matrix assembly** | Single-process assembly | MPI-parallel with allgather |
| **Eigenvalue methods** | Arnoldi + Richardson available | Direct solve only |
| **Memory scaling** | O(lag³ × nsteps_total) | O(lag³ × nsteps_local) per process |
| **Typical nsteps** | ~1000-10000 | ~10000-100000 (distributed) |

### Algorithmic Implementation Differences

**NEO-2-QL Multi-method Approach**:
```fortran
! Option 1: Arnoldi eigenvalue analysis for stability
call arnoldi_iteration(matrix, rhs, eigenvals, eigenvecs)
call richardson_preconditioned(matrix, rhs, solution, eigenvecs)

! Option 2: Direct sparse solver  
call sparse_solve(matrix, rhs, solution, method=3)  ! UMFPACK
```

**NEO-2-PAR Direct-only Approach**:
```fortran
! MPI-parallel collision operator assembly
do a = 0, num_spec-1
   call compute_collop('a', 'b', m_spec(a), m_spec(b), anumm_aa(:,:,a,b), ...)
end do
call mpro%allgather_inplace(anumm_aa)  ! Distribute collision data

! Direct sparse solver only
call sparse_solve(matrix, rhs, solution, method=3)  ! UMFPACK
```

### Critical Memory Bottleneck Insight

**Both QL and PAR hit identical velocity space limitations**:
1. **Collision operators**: O(lag² × leg × nspecies²) scaling identical
2. **Main matrix size**: Both scale as lag³, only constants differ
3. **UMFPACK factorization**: Same ~5-10x memory penalty
4. **Practical limits**: Both reach memory walls at similar lag values (~20-30)

**Key difference**: PAR distributes field line work but **not collision operator work**
- Field line scaling: QL O(nsteps_total) vs PAR O(nsteps_local)  
- Collision scaling: **Both O(lag³) with identical computational kernels**
- Memory bottleneck: **Collision operators dominate for both variants**

### Operator Structure Comparison

**Differential Operators (Similar)**:
- Both use same finite difference schemes for field line derivatives
- Both use same basis function derivatives in velocity space
- Field line complexity differs but discretization methods identical

**Integral Operators (Identical)**:  
- Same Rosenbluth potential computation algorithms
- Same GSL numerical integration for collision integrals
- Same precomputation and storage strategies for I1-I4 matrices

**Matrix Assembly (Different parallelization)**:
- QL: Single-process assembly, full matrix storage
- PAR: MPI-parallel assembly with communication, distributed storage

### Computational Complexity

**Critical Insight: Velocity Space Parameters Dominate Memory Scaling**

Both `lag` (energy/speed) and `leg` (pitch angle) parameters create significant memory scaling:

**Complete Memory Scaling Analysis**:

**QL Complexity**:
- **Collision operator coefficients**: O(lag² × leg × nspecies²)
- **Rosenbluth integrals**: O(lag² × leg) for each of I1-I4 matrices  
- **Main sparse matrix**: O(lag³ × nsteps) where nsteps ~ field line discretization
- **Distribution functions**: O(lag × 4 × npart × nperiods)
- **UMFPACK factorization**: O(nnz^1.2-1.5 × fill_factor) where nnz ∝ lag³ × nsteps
- **Total memory**: O(lag³ × nsteps × fill_factor) + O(lag² × leg × nspecies²)

**PAR Complexity (per MPI process)**:
- **Same collision operator scaling**: O(lag² × leg × nspecies²) per process
- **Local sparse matrix**: O(lag³ × nsteps_local) where nsteps_local < nsteps_total
- **Distribution functions**: O(lag × 4 × npart_local × nperiods_local)  
- **UMFPACK factorization**: O(nnz_local^1.2-1.5 × fill_factor)
- **Total memory per process**: O(lag³ × nsteps_local × fill_factor) + O(lag² × leg × nspecies²)

**Memory Bottleneck Hierarchy**:
1. **Primary**: O(lag³) from UMFPACK factorization of collision-dominated sparse matrix
2. **Secondary**: O(lag² × leg × nspecies²) from collision operator coefficient storage
3. **Tertiary**: O(lag × nsteps) from distribution function storage
4. **Quaternary**: Field line discretization (nsteps vs nsteps_local)

**Why both QL and PAR hit similar limits**:
- Collision operator memory scaling **identical** between variants
- UMFPACK factorization penalty **identical** (5-10x matrix storage)
- Field line distribution only affects constants, not asymptotic scaling
- Practical memory limits reached at **similar lag values (~20-30)** for both

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

**Current UMFPACK (velocity space scaling)**:
- **Matrix storage**: `nnz × 16` bytes where nnz ∝ lag³ × nsteps  
- **Collision operators**: O(lag² × leg × nspecies²) coefficient matrices
- **Factorization memory**: `~5-10 × nnz × 8` bytes for LU factors
- **Total for lag=30, leg=20**: ~50-80 GB (factorization dominates)
- **Memory scaling**: O(lag³ × nsteps × fill_factor) + O(lag² × leg × nspecies²)

**Proposed GMRES**:
- **Matrix storage**: Same `nnz × 16` bytes (matrix-vector products only)
- **Collision operators**: Same O(lag² × leg × nspecies²) coefficient storage  
- **Krylov basis**: `m × n × 8` bytes where n ∝ lag × nsteps, m = restart dim
- **Hessenberg matrix**: `m² × 8` bytes (typically m = 30-50)
- **Working vectors**: O(n) temporary storage for matrix-vector products
- **Total for lag=30, leg=20**: ~5-10 GB (eliminates factorization memory)
- **Memory scaling**: O(lag × nsteps) + O(lag² × leg × nspecies²)

**Key GMRES Memory Advantage**:
- **Eliminates O(lag³) UMFPACK factorization memory** 
- **Retains O(lag² × leg) collision operator storage** (unavoidable)
- **Net scaling reduction**: From O(lag³) to O(lag²) for large lag values
- **Practical impact**: Enables lag ~50-100 instead of lag ~20-30 limit

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