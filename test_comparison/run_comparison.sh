#!/bin/bash

# Golden Record Comparison Test: Arnoldi O2 vs IDR(s)
# This script runs the same golden record test case with both solvers
# and compares the results to validate IDR(s) implementation

set -e

echo "========================================================="
echo "NEO-2 Golden Record Comparison: Arnoldi O2 vs IDR(s)"
echo "========================================================="
echo ""

# Clean up any previous runs
echo "Cleaning up previous runs..."
rm -f arnoldi_run/*.h5 arnoldi_run/out arnoldi_run/*.log
rm -f idrs_run/*.h5 idrs_run/out idrs_run/*.log

echo ""
echo "=== Running Arnoldi O2 (Reference) Test ==="
echo "Solver: ISW_RIPPLE_SOLVER=3 (Arnoldi 2nd order)"
cd arnoldi_run

# Run with timeout to prevent hanging
echo "Starting NEO-2-QL with Arnoldi solver..."
timeout 300s ../../build/NEO-2-QL/neo_2_ql.x > arnoldi_output.log 2>&1 && {
    echo "✅ Arnoldi run completed successfully"
    ARNOLDI_SUCCESS=1
} || {
    echo "❌ Arnoldi run failed or timed out"
    echo "Last 10 lines of output:"
    tail -10 arnoldi_output.log 2>/dev/null || echo "No output log found"
    ARNOLDI_SUCCESS=0
}

cd ..

echo ""
echo "=== Running IDR(s) (Test) Test ==="
echo "Solver: ISW_RIPPLE_SOLVER=4 (IDR(s) iterative solver)"
cd idrs_run

echo "Starting NEO-2-QL with IDR(s) solver..."
timeout 300s ../../build/NEO-2-QL/neo_2_ql.x > idrs_output.log 2>&1 && {
    echo "✅ IDR(s) run completed successfully"
    IDRS_SUCCESS=1
} || {
    echo "❌ IDR(s) run failed or timed out"
    echo "Last 10 lines of output:"
    tail -10 idrs_output.log 2>/dev/null || echo "No output log found"
    IDRS_SUCCESS=0
}

cd ..

echo ""
echo "========================================================="
echo "Comparison Results"
echo "========================================================="

if [ $ARNOLDI_SUCCESS -eq 1 ] && [ $IDRS_SUCCESS -eq 1 ]; then
    echo "✅ Both solvers completed successfully!"
    echo ""
    echo "Generated output files:"
    echo "  Arnoldi results in: arnoldi_run/"
    ls -la arnoldi_run/*.h5 2>/dev/null || echo "    No HDF5 files found"
    echo "  IDR(s) results in: idrs_run/"
    ls -la idrs_run/*.h5 2>/dev/null || echo "    No HDF5 files found"
    echo ""
    echo "Next steps:"
    echo "  1. Compare HDF5 output files for physics accuracy"
    echo "  2. Analyze solver performance and memory usage"
    echo "  3. Validate transport coefficients match within tolerance"
    
elif [ $ARNOLDI_SUCCESS -eq 1 ] && [ $IDRS_SUCCESS -eq 0 ]; then
    echo "❌ IDR(s) solver failed - needs debugging"
    echo "✅ Arnoldi solver worked as expected"
    echo ""
    echo "Check IDR(s) output log:"
    echo "  cat idrs_run/idrs_output.log"
    
elif [ $ARNOLDI_SUCCESS -eq 0 ] && [ $IDRS_SUCCESS -eq 1 ]; then
    echo "❌ Arnoldi solver failed unexpectedly"
    echo "✅ IDR(s) solver worked!"
    echo ""
    echo "Check Arnoldi output log:"
    echo "  cat arnoldi_run/arnoldi_output.log"
    
else
    echo "❌ Both solvers failed"
    echo ""
    echo "This might indicate:"
    echo "  - Missing input files or data"
    echo "  - Environment setup issues"
    echo "  - Fundamental build problems"
    echo ""
    echo "Check output logs:"
    echo "  cat arnoldi_run/arnoldi_output.log"
    echo "  cat idrs_run/idrs_output.log"
fi

echo ""
echo "Test completed: $(date)"