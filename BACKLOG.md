# NEO-2 Development Backlog

## BiCGSTAB with ILU(1) Preconditioner Implementation Plan

### Overview
Replace the current UMFPACK direct solver with BiCGSTAB iterative solver using ILU(1) preconditioning to achieve:
- **2-5x memory reduction** (from O(lag³) to O(lag))
- **20-60% runtime improvement** for large problems
- **Better scalability** for high-resolution velocity space (large lag values)

### Key Decision: ILU Implementation Strategy

**Recommendation: Implement our own ILU(1) factorization**

Reasons:
1. **Full control** over memory layout and optimization
2. **Tight integration** with BiCGSTAB solver
3. **Avoid dependency issues** with SuiteSparse ILU routines
4. **Custom optimization** for NEO-2's sparse matrix structure
5. **Easier debugging** and profiling

SuiteSparse (UMFPACK) doesn't provide standalone ILU - it uses complete LU factorization. We would need to add another dependency (e.g., ILUPACK) or implement our own.

### Implementation Phases

## Phase 1: Core Infrastructure (Week 1)

### 1.1 Sparse Matrix Utilities Module
**File:** `COMMON/sparse_utils_mod.f90`
- [ ] CSR (Compressed Sparse Row) format conversion routines
- [ ] CSC ↔ CSR conversion utilities
- [ ] Matrix-vector multiplication for CSR format
- [ ] Diagonal extraction routines
- [ ] **Unit tests:** Verify conversions preserve matrix structure

### 1.2 ILU(1) Preconditioner Module
**File:** `COMMON/ilu_precond_mod.f90`
- [ ] ILU(1) factorization with level-of-fill = 1
- [ ] Forward/backward substitution solvers
- [ ] Memory-efficient storage for L and U factors
- [ ] Drop tolerance parameter support
- [ ] **Unit tests:** 
  - Small test matrices with known factorizations
  - Verify L*U ≈ A within tolerance
  - Test singular/near-singular matrix handling

### 1.3 BiCGSTAB Core Module
**File:** `COMMON/bicgstab_mod.f90`
- [ ] Basic BiCGSTAB algorithm implementation
- [ ] Convergence monitoring and criteria
- [ ] Residual norm calculation
- [ ] Iteration history tracking
- [ ] **Unit tests:**
  - Solve Ax=b for diagonal matrices
  - Solve small SPD systems
  - Verify convergence for well-conditioned problems

## Phase 2: Algorithm Implementation (Week 2)

### 2.1 Enhanced BiCGSTAB Features
- [ ] Preconditioned BiCGSTAB with ILU(1)
- [ ] Restart capability for stability
- [ ] Adaptive tolerance adjustment
- [ ] **Integration tests:**
  - Test with ILU(1) preconditioner
  - Compare convergence with/without preconditioning
  - Benchmark iteration counts

### 2.2 Robustness Enhancements
- [ ] Breakdown detection and recovery
- [ ] Stagnation detection
- [ ] Fallback to unpreconditioned iteration
- [ ] **Stress tests:**
  - Ill-conditioned matrices
  - Near-singular systems
  - Matrices with zero pivots

### 2.3 Performance Optimizations
- [ ] Cache-friendly data access patterns
- [ ] Vectorized dot products and axpy operations
- [ ] OpenMP parallelization for matrix-vector products
- [ ] **Performance tests:**
  - Profile hot spots
  - Measure FLOPS efficiency
  - Compare with BLAS implementations

## Phase 3: Integration (Week 3)

### 3.1 Sparse Module Integration
**File:** Update `COMMON/sparse_mod.f90`
- [ ] Add `sparse_solve_method = 4` for BiCGSTAB
- [ ] Implement wrapper routines for existing interface
- [ ] Support both real and complex systems
- [ ] **Integration tests:**
  - Verify same interface behavior
  - Test method switching (3 → 4)
  - Ensure backward compatibility

### 3.2 Test Suite Development
**File:** `COMMON/test_solvers_mod.f90`
- [ ] Generate test matrices:
  - Tridiagonal systems
  - Random sparse matrices
  - NEO-2 representative matrices
- [ ] Comparative solver framework:
  - UMFPACK (method 3)
  - BiCGSTAB (method 4)
  - Error norms and timing
- [ ] **Validation tests:**
  - Compare solutions to UMFPACK
  - Verify ||Ax - b|| / ||b|| < tolerance
  - Check conservation properties

## Phase 4: NEO-2 Specific Testing (Week 4)

### 4.1 Collision Operator Tests
- [ ] Extract collision operator matrices from NEO-2 runs
- [ ] Test BiCGSTAB convergence on real physics problems
- [ ] Compare memory usage with UMFPACK
- [ ] **Physics validation:**
  - Conservation of particles, momentum, energy
  - Compare transport coefficients
  - Verify numerical stability

### 4.2 Scalability Analysis
- [ ] Test with increasing lag values (10, 20, 30, 50)
- [ ] Memory profiling at each scale
- [ ] Runtime comparisons
- [ ] **Benchmark suite:**
  - QL small/medium/large cases
  - PAR test cases
  - Document speedup factors

### 4.3 Production Readiness
- [ ] Error handling and user messages
- [ ] Documentation and examples
- [ ] Configuration parameters in input files
- [ ] **Final validation:**
  - Run golden record tests
  - Compare with published results
  - Stress test on cluster

## Implementation Guidelines

### Code Standards
1. **Module structure:** One feature per module, clear interfaces
2. **Naming conventions:** Follow NEO-2 style (lowercase with underscores)
3. **Documentation:** Doxygen-style comments for all public routines
4. **Error handling:** Graceful degradation, informative messages

### Testing Strategy
1. **Unit tests first:** Test each routine in isolation
2. **Integration tests:** Test module interactions
3. **Validation tests:** Compare with known solutions
4. **Performance tests:** Profile and optimize
5. **Physics tests:** Ensure conservation laws

### Performance Targets
- **Memory:** < 20% of UMFPACK usage for large problems
- **Runtime:** 2-3x faster than UMFPACK for lag > 30
- **Convergence:** < 500 iterations for typical problems
- **Accuracy:** ||Ax - b|| / ||b|| < 1e-12

### Risk Mitigation
1. **Fallback option:** Keep UMFPACK as method 3
2. **Adaptive switching:** Auto-select solver based on problem size
3. **Extensive testing:** Each phase fully tested before proceeding
4. **Incremental integration:** Add features gradually

## Success Metrics

1. **Memory Reduction:** Achieve 5x reduction for lag=50 cases
2. **Performance:** 2-3x speedup on production problems
3. **Reliability:** Pass all golden record tests
4. **Maintainability:** Clean, documented, testable code

## Next Steps

1. Create feature branch: `feature/bicgstab-ilu-solver`
2. Set up test framework infrastructure
3. Begin Phase 1.1 implementation
4. Weekly progress reviews and adjustments

---

*Last updated: 2025-08-01*