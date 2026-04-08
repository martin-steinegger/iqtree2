#!/usr/bin/env python3
"""
Convert a CherryML rate matrix to IQ-TREE2 custom model format.

CherryML format:
  - Tab-separated, first row is header with state names
  - Subsequent rows: state_label <tab> Q[i][0] <tab> Q[i][1] <tab> ...
  - Full rate matrix Q (400x400 for 3Di+AA, includes diagonal)

IQ-TREE format:
  - 400 rows of tab-separated rate values (full Q matrix, no labels)
  - 1 row of 400 equilibrium state frequencies (computed from Q)

Usage:
  python3 convert_cherryml_to_iqtree.py <cherryml_matrix.txt> <output_iqtree.txt>

Example:
  python3 convert_cherryml_to_iqtree.py matrix3diaa/matrix.txt test_data/matrix3diaa_iqtree.txt
"""

import sys
import re
import numpy as np


def read_cherryml(path):
    """Read a CherryML rate matrix, returning Q and state names."""
    with open(path, "rb") as f:
        raw = f.read()
    # Remove any ANSI escape sequences that may be embedded
    clean = re.sub(rb"\x1b\[[^\x40-\x7E]*[\x40-\x7E]", b"", raw).decode()
    lines = clean.strip().split("\n")

    # First line: header with state names
    state_names = lines[0].strip().split("\t")
    # Remove empty first element if present (leading tab)
    if state_names[0] == "":
        state_names = state_names[1:]

    # Data rows: label + values
    Q = []
    for line in lines[1:]:
        parts = line.strip().split("\t")
        Q.append([float(x) for x in parts[1:]])  # skip row label

    return np.array(Q), state_names


def equilibrium_frequencies(Q):
    """Compute equilibrium frequencies pi such that pi @ Q = 0."""
    eigenvalues, eigenvectors = np.linalg.eig(Q.T)
    idx = np.argmin(np.abs(eigenvalues))
    pi = np.real(eigenvectors[:, idx])
    pi = pi / pi.sum()
    return pi


def write_iqtree(path, Q, pi):
    """Write IQ-TREE custom model file: full rate matrix + frequencies."""
    n = Q.shape[0]
    with open(path, "w") as f:
        for i in range(n):
            f.write("\t".join(f"{Q[i][j]}" for j in range(n)) + "\n")
        f.write("\t".join(f"{pi[j]}" for j in range(n)) + "\n")


def main():
    if len(sys.argv) != 3:
        print("Usage: python3 convert_cherryml_to_iqtree.py <input.txt> <output.txt>")
        sys.exit(1)

    Q, state_names = read_cherryml(sys.argv[1])
    n = Q.shape[0]
    print(f"Read {n}x{n} rate matrix with {len(state_names)} states")
    print(f"States: {state_names[0]} ... {state_names[-1]}")

    pi = equilibrium_frequencies(Q)
    print(f"Equilibrium frequencies sum: {pi.sum():.10f}")
    print(f"Min freq: {pi.min():.6e}, Max freq: {pi.max():.6e}")

    # Sanity: check detailed balance (reversibility)
    db_err = max(
        abs(Q[i][j] * pi[i] - Q[j][i] * pi[j])
        for i in range(min(n, 20))
        for j in range(i + 1, min(n, 20))
    )
    print(f"Max detailed-balance error (first 20 states): {db_err:.2e}")

    write_iqtree(sys.argv[2], Q, pi)
    print(f"Written: {sys.argv[2]} ({n} rate rows + 1 frequency row)")


if __name__ == "__main__":
    main()
