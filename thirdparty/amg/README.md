# AMG Preconditioner for NEO-2

This directory contains a Fortran implementation of Algebraic Multigrid (AMG) preconditioners based on the AlgebraicMultigrid.jl package.

## Attribution

The algorithms and implementation strategies in this module are derived from:
- **AlgebraicMultigrid.jl** (https://github.com/JuliaLinearAlgebra/AlgebraicMultigrid.jl)
- Copyright (c) 2017-20: ranjanan, Julia Computing
- Licensed under MIT License (see LICENSE file)

## Overview

This implementation provides AMG preconditioners as a robust alternative to ILU for the NEO-2 sparse solver framework. AMG is particularly effective for:
- Ill-conditioned matrices with structural zeros on the diagonal
- Large sparse systems from discretized PDEs
- Problems where ILU factorization fails

## Features

### Implemented AMG Variants
1. **Classical AMG (Ruge-Stüben)**
   - Strong connection based coarsening
   - Standard interpolation
   - Suitable for M-matrices and diagonally dominant systems

2. **Smoothed Aggregation AMG**
   - Aggregation-based coarsening
   - Smoothed prolongation operators
   - Better for systems from finite element discretizations

### Multigrid Cycles
- V-cycle (default)
- W-cycle (more robust, higher cost)
- F-cycle (future implementation)

### Smoothers
- Gauss-Seidel (forward, backward, symmetric)
- Jacobi (weighted)
- Block Jacobi (future)

## Usage

The AMG preconditioner integrates seamlessly with NEO-2's solver framework:

```fortran
! In solver configuration
solver_preconditioner = PRECOND_AMG
amg_variant = AMG_CLASSICAL  ! or AMG_SMOOTHED_AGGREGATION
amg_cycle = AMG_VCYCLE       ! or AMG_WCYCLE
amg_smoother = AMG_GAUSS_SEIDEL
amg_presmoother_steps = 2
amg_postsmoother_steps = 2
amg_coarsest_size = 50
```

## Implementation Status

- [x] Basic module structure
- [ ] Classical AMG coarsening
- [ ] Interpolation operators
- [ ] Galerkin coarse grid operators
- [ ] V-cycle implementation
- [ ] Smoother implementations
- [ ] Integration with sparse_solvers_mod
- [ ] Performance optimization
- [ ] Comprehensive testing

## Technical Details

### Coarsening Strategy (Classical AMG)
- Strong connections: |a_ij| ≥ θ * max|a_ik|
- Standard Ruge-Stüben C/F splitting
- Default strength threshold θ = 0.25

### Interpolation
- Direct interpolation from C-points
- Weights based on matrix coefficients
- Preserves constant vectors (important for Poisson-like problems)

### Complexity Control
- Aggressive coarsening for memory efficiency
- Grid complexity target: < 1.5
- Operator complexity target: < 2.5

## Performance Considerations

1. **Setup Cost**: AMG has higher setup cost than ILU but often requires fewer iterations
2. **Memory Usage**: Typically 1.5-2.5x the original matrix
3. **Parallelization**: Smoother operations are naturally parallel
4. **Cache Efficiency**: Coarse grids fit in cache, improving performance

## References

1. Briggs, W. L., Henson, V. E., & McCormick, S. F. (2000). A multigrid tutorial. SIAM.
2. Stüben, K. (2001). A review of algebraic multigrid. Journal of Computational and Applied Mathematics, 128(1-2), 281-309.
3. Vaněk, P., Mandel, J., & Brezina, M. (1996). Algebraic multigrid by smoothed aggregation for second and fourth order elliptic problems. Computing, 56(3), 179-196.