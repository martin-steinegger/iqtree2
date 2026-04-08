#!/bin/bash
# Run IQ-TREE2 tree inference with a 400-state (3Di+AA) custom rate matrix.
#
# Prerequisites:
#   1. Build iqtree2 with BLAS support (see below)
#   2. Prepare input files:
#      a. Numeric PHYLIP alignment (space-separated integers, gap = -1)
#         Use: python3 scripts/convert_foldmason_to_phylip.py <aa.fa> <ss.fa> <output.phy>
#      b. IQ-TREE rate matrix (full Q matrix + frequencies)
#         Use: python3 scripts/convert_cherryml_to_iqtree.py <cherryml.txt> <matrix_iq.txt>
#
# Build command (macOS with Apple Silicon):
#   cd build && cmake .. -DIQTREE_FLAGS=omp \
#     -DEIGEN3_INCLUDE_DIR=/opt/homebrew/Cellar/eigen/5.0.1/include/eigen3 \
#     -DOpenMP_CXX_FLAGS="-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include" \
#     -DOpenMP_CXX_LIB_NAMES="omp" \
#     -DOpenMP_omp_LIBRARY=/opt/homebrew/opt/libomp/lib/libomp.dylib \
#     -DCMAKE_CXX_FLAGS="-I/opt/homebrew/opt/libomp/include" \
#     -DCMAKE_EXE_LINKER_FLAGS="-L/opt/homebrew/opt/libomp/lib -lomp"
#   make -j8
#
# Build command (Linux with OpenBLAS):
#   cd build && cmake .. -DIQTREE_FLAGS=omp
#   make -j8
#   # OpenBLAS must be available (apt install libopenblas-dev, or via conda)
#
# Usage:
#   ./scripts/run_400state_tree.sh <alignment.phy> <matrix.txt> [output_prefix] [threads]

set -euo pipefail

IQTREE="${IQTREE:-build/iqtree2}"
ALIGNMENT="${1:?Usage: $0 <alignment.phy> <matrix.txt> [prefix] [threads]}"
MATRIX="${2:?Usage: $0 <alignment.phy> <matrix.txt> [prefix] [threads]}"
PREFIX="${3:-output_400state}"
THREADS="${4:-4}"

echo "=== IQ-TREE2 400-state tree inference ==="
echo "Alignment: ${ALIGNMENT}"
echo "Matrix:    ${MATRIX}"
echo "Prefix:    ${PREFIX}"
echo "Threads:   ${THREADS}"
echo ""

${IQTREE} \
    -s "${ALIGNMENT}" \
    -st MULTI \
    -m "${MATRIX}" \
    -B 1000 \
    -pre "${PREFIX}" \
    -nt "${THREADS}" \
    --redo

echo ""
echo "=== Done ==="
echo "Tree:      ${PREFIX}.treefile"
echo "Log:       ${PREFIX}.log"
echo "Report:    ${PREFIX}.iqtree"
