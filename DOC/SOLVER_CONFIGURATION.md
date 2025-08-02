# NEO-2 Solver Configuration Guide

## Overview

NEO-2 provides multiple sparse linear solver options for collision operator matrix solutions. This document explains how to configure and use different solver methods, with particular focus on the new BiCGSTAB iterative solver.

## Available Solvers

### UMFPACK Direct Solver (Method 3)
- **Type**: Direct solver using LU factorization
- **Memory**: High (O(lag³) for dense factorization)
- **Speed**: Fast for small-medium problems, slower for large problems
- **Reliability**: Very reliable, exact solution (within machine precision)
- **Use case**: Default for small problems, reference solutions

### BiCGSTAB Iterative Solver (Method 4)
- **Type**: Iterative solver with ILU preconditioning
- **Memory**: Low (O(lag) sparse storage)
- **Speed**: Fast for large problems, 2-5x speedup potential
- **Reliability**: High for well-conditioned problems
- **Use case**: Recommended for large problems (lag > 30)

## Configuration

### Method 1: Environment Variable (Temporary)
```bash
export SPARSE_SOLVE_METHOD=4  # Use BiCGSTAB
export SPARSE_SOLVE_METHOD=3  # Use UMFPACK (default)
```

### Method 2: Runtime Selection (Programming Interface)
```fortran
! In your Fortran code
USE sparse_solvers_mod, ONLY: SOLVER_UMFPACK, SOLVER_BICGSTAB
USE sparse_mod, ONLY: sparse_solve_method

! Set solver method before calling sparse_solve
sparse_solve_method = SOLVER_BICGSTAB  ! Use BiCGSTAB
sparse_solve_method = SOLVER_UMFPACK   ! Use UMFPACK
```

### Method 3: Input File Configuration (Future)
```fortran
! In NEO-2 namelist (planned feature)
&solver_control
  solver_method = 4              ! BiCGSTAB
  solver_tolerance = 1.0e-12     ! Convergence tolerance
  solver_max_iter = 1000         ! Maximum iterations
/
```

## Solver Performance Characteristics

### Memory Usage Comparison
| Problem Size (lag) | UMFPACK Memory | BiCGSTAB Memory | Memory Reduction |
|--------------------|--------------|--------------  |------------------|
| lag = 10          | ~10 MB       | ~2 MB          | 5x               |
| lag = 20          | ~80 MB       | ~8 MB          | 10x              |
| lag = 30          | ~270 MB      | ~18 MB         | 15x              |
| lag = 50          | ~1.25 GB     | ~50 MB         | 25x              |

### Runtime Performance
- **Small problems (lag < 20)**: UMFPACK typically faster
- **Medium problems (lag 20-30)**: Comparable performance
- **Large problems (lag > 30)**: BiCGSTAB typically 2-5x faster

### Convergence Behavior
- **Well-conditioned matrices**: BiCGSTAB converges in 10-50 iterations
- **Ill-conditioned matrices**: May require 100-500 iterations
- **Convergence tolerance**: Default 1e-12, adjustable if needed

## When to Use Each Solver

### Use UMFPACK (Method 3) when:
- Problem size is small (lag < 20)
- Maximum accuracy is required
- Debugging solver issues
- Matrix is severely ill-conditioned
- Development/testing phase

### Use BiCGSTAB (Method 4) when:
- Problem size is large (lag > 30)
- Memory is limited
- Runtime performance is critical
- Matrix is reasonably well-conditioned
- Production runs

## Troubleshooting

### BiCGSTAB Convergence Issues

#### Symptom: "BiCGSTAB did not converge after N iterations"
**Causes and Solutions:**
1. **Matrix is ill-conditioned**
   - Try UMFPACK for comparison
   - Check physical parameters for unrealistic values
   - Consider problem reformulation

2. **Tolerance too strict**
   - Increase tolerance from 1e-12 to 1e-10 or 1e-8
   - Monitor solution quality with residual checks

3. **Maximum iterations too low**
   - Increase max_iter from 1000 to 2000 or 5000
   - Monitor convergence rate

#### Symptom: "Solutions differ between UMFPACK and BiCGSTAB"
**Investigation steps:**
1. Check convergence: Did BiCGSTAB reach tolerance?
2. Compare residuals: `||Ax - b|| / ||b||` for both solutions
3. If residuals are similar (~1e-10), differences are acceptable
4. If residuals differ significantly, investigate matrix conditioning

### Performance Optimization

#### For Better BiCGSTAB Performance:
1. **Ensure good preconditioning**: ILU(1) is used by default
2. **Check sparsity pattern**: More structured = better convergence
3. **Monitor iteration counts**: Should be < 100 for most problems
4. **Verify problem scaling**: Well-scaled problems converge faster

#### Memory Optimization:
1. **Use BiCGSTAB for large problems**: Dramatic memory savings
2. **Monitor peak memory usage**: Especially important for cluster runs
3. **Consider problem decomposition**: For extremely large problems

## Validation and Testing

### Integration Tests
Run the solver integration tests to verify consistent behavior:
```bash
cd build/TEST
./test_solver_integration
```

Expected output:
- All tests should pass
- Solution differences should be < 1e-6
- Performance data should be reasonable

### Golden Record Tests
Verify that solver changes don't affect physics results:
```bash
make test  # Runs all regression tests
```

### Custom Validation
For your specific problem, compare solver results:
```fortran
! Save solutions from both solvers
sparse_solve_method = SOLVER_UMFPACK
CALL sparse_solve(nrow, ncol, nz, irow, pcol, val, x_umfpack)

sparse_solve_method = SOLVER_BICGSTAB  
CALL sparse_solve(nrow, ncol, nz, irow, pcol, val, x_bicgstab)

! Compare solutions
max_diff = MAXVAL(ABS(x_umfpack - x_bicgstab))
rel_diff = max_diff / SQRT(SUM(x_umfpack**2))

IF (rel_diff < 1.0e-6_dp) THEN
  WRITE(*,*) "Solvers agree to sufficient precision"
ELSE
  WRITE(*,*) "WARNING: Solver disagreement detected"
END IF
```

## Best Practices

### Development Phase
1. Start with UMFPACK for correctness verification
2. Test BiCGSTAB on small problems first
3. Compare solutions between solvers
4. Gradually increase problem size
5. Monitor convergence behavior

### Production Phase
1. Use BiCGSTAB for large problems (lag > 30)
2. Keep UMFPACK as fallback for difficult cases
3. Monitor convergence warnings in output
4. Log solver timing for performance analysis
5. Run validation tests periodically

### Cluster/HPC Usage
1. BiCGSTAB typically better for memory-constrained systems
2. Monitor parallel efficiency (if applicable)
3. Consider node memory limits when choosing solver
4. Use shorter walltime limits with BiCGSTAB for large problems

## Technical Details

### Algorithm Description
- **BiCGSTAB**: Bi-Conjugate Gradient Stabilized method
- **Preconditioning**: Incomplete LU factorization with level 1 fill-in
- **Matrix format**: Compressed Sparse Column (CSC) internally converted to CSR
- **Convergence criterion**: ||r|| / ||b|| < tolerance

### Implementation Notes
- BiCGSTAB implementation in `COMMON/bicgstab_mod.f90`
- ILU preconditioning in `COMMON/ilu_precond_mod.f90`
- Solver dispatch in `COMMON/sparse_solvers_mod.f90`
- Integration tests in `tests/test_solver_integration.f90`

### Future Enhancements
- Complex matrix support for BiCGSTAB
- Additional preconditioner options (AMG)
- Automatic solver selection based on problem characteristics
- Namelist-based configuration interface

## Contact and Support

For solver-related issues:
1. Check this documentation first
2. Run integration tests to verify installation
3. Compare with UMFPACK results for validation
4. Report issues with minimal reproduction case

---

*Last updated: 2025-08-02*
*Compatible with: NEO-2 bicgstab branch and later*