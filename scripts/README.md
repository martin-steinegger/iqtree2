# 400-State (3Di+AA) IQ-TREE2 Scripts

Scripts for running IQ-TREE2 with large-alphabet (up to 400-state) substitution models, designed for joint 3Di+AA structural-sequence phylogenetics.

## Prerequisites

- Python 3 with NumPy
- IQ-TREE2 built with BLAS support (see build instructions below)
- [FoldMason](https://github.com/steineggerlab/foldmason) alignments (paired AA + 3Di FASTA files)

## Scripts

### `convert_cherryml_to_iqtree.py`

Converts a CherryML rate matrix (tab-separated, with header row and row labels) to IQ-TREE2 custom model format (full Q matrix + equilibrium frequencies).

```bash
python3 scripts/convert_cherryml_to_iqtree.py matrix3diaa/matrix.txt my_matrix_iq.txt
```

### `convert_foldmason_to_phylip.py`

Converts paired AA and 3Di FASTA alignments (FoldMason output) into a numeric PHYLIP file with space-separated integer states for `--seqtype MULTI`.

State encoding: `state = aa_idx * 20 + ss_idx` (0..399), gaps encoded as `-1`.

```bash
python3 scripts/convert_foldmason_to_phylip.py aa.fa ss.fa alignment.phy
```

### `run_400state_tree.sh`

Runs IQ-TREE2 tree inference with a custom rate matrix and ultrafast bootstrap.

```bash
./scripts/run_400state_tree.sh <alignment.phy> <matrix.txt> [output_prefix] [threads]

# Example:
./scripts/run_400state_tree.sh example3diaa/ferritin_remapped.phy test_data/ferritin_sub.txt ferritin_out 4
```

Numeric PHYLIP files use **negative integers** (typically `-1`) for gaps. Do not use `N` or the state count as a gap value -- only negative values are treated as unknown/gap by IQ-TREE2.

## Test data

- `test_data/seqs_400state.phy` -- 5 sequences, 10 sites, 400 states (quick sanity check)
- `test_data/matrix3diaa_iqtree.txt` -- full 400x400 3Di+AA rate matrix
- `test_data/ferritin_sub.txt` -- 307x307 sub-matrix for ferritin dataset
- `example3diaa/ferritin_remapped.phy` -- 53 ferritin sequences, 307 states
