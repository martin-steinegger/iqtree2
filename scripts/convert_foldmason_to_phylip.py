#!/usr/bin/env python3
"""
Convert paired AA + 3Di FASTA alignments into a PHYLIP file
with space-separated integer states for IQ-TREE2 --seqtype MULTI.

State encoding (matching matrix3diaa):
  state = aa_idx * 20 + ss_idx
  Alphabet order (both AA and 3Di use the same): A R N D C Q E G H I L K M F P S T W Y V
  Gaps (-) or unknown (X) in either channel → -1 (treated as STATE_UNKNOWN by IQ-TREE)

Usage:
  python3 convert_3diaa.py <aa.fa> <ss.fa> <output.phy>
"""

import sys

# Alphabet order matching the matrix3diaa state names (A|A, A|R, ..., V|V)
ALPHABET = list("ARNDCQEGHILKMFPSTWYV")
CHAR_TO_IDX = {c: i for i, c in enumerate(ALPHABET)}

def read_fasta(path):
    """Return ordered list of (name, sequence) from a FASTA file."""
    seqs = []
    name, seq = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    seqs.append((name, "".join(seq)))
                name = line[1:]
                seq = []
            else:
                seq.append(line)
    if name is not None:
        seqs.append((name, "".join(seq)))
    return seqs

def convert(aa_path, ss_path, out_path):
    aa_seqs = read_fasta(aa_path)
    ss_seqs = read_fasta(ss_path)

    if len(aa_seqs) != len(ss_seqs):
        raise ValueError(f"Sequence count mismatch: {len(aa_seqs)} AA vs {len(ss_seqs)} 3Di")

    # Check name correspondence
    for (aa_name, _), (ss_name, _) in zip(aa_seqs, ss_seqs):
        if aa_name != ss_name:
            raise ValueError(f"Name mismatch: '{aa_name}' vs '{ss_name}'")

    nseq = len(aa_seqs)
    nsite = len(aa_seqs[0][1])

    for name, seq in aa_seqs:
        if len(seq) != nsite:
            raise ValueError(f"Length mismatch for {name}")
    for name, seq in ss_seqs:
        if len(seq) != nsite:
            raise ValueError(f"Length mismatch for {name}")

    # Convert to integer state sequences
    int_seqs = []
    n_gaps = 0
    n_unknown_aa = 0
    n_unknown_ss = 0
    max_state = 0

    for (name, aa_seq), (_, ss_seq) in zip(aa_seqs, ss_seqs):
        states = []
        for aa_char, ss_char in zip(aa_seq, ss_seq):
            if aa_char == "-" or ss_char == "-":
                states.append(-1)
                n_gaps += 1
            elif aa_char not in CHAR_TO_IDX:
                states.append(-1)
                n_unknown_aa += 1
            elif ss_char not in CHAR_TO_IDX:
                states.append(-1)
                n_unknown_ss += 1
            else:
                s = CHAR_TO_IDX[aa_char] * 20 + CHAR_TO_IDX[ss_char]
                states.append(s)
                if s > max_state:
                    max_state = s
        int_seqs.append((name, states))

    print(f"Sequences: {nseq}, sites: {nsite}")
    print(f"Gaps/unknowns: {n_gaps} positions, unknown AA: {n_unknown_aa}, unknown 3Di: {n_unknown_ss}")
    print(f"Max state index: {max_state} (alphabet size = {max_state + 1})")
    unique_states = set(s for _, states in int_seqs for s in states if s >= 0)
    print(f"Unique states observed: {len(unique_states)}")

    # Build unique 10-char PHYLIP names.
    # For any prefix that collides, ALL sequences sharing it get a 2-digit counter
    # appended to an 8-char prefix (e.g. WP_301497 + 00..NN) to guarantee uniqueness.
    from collections import Counter
    prefix_freq = Counter(name[:10] for name, _ in int_seqs)
    prefix_counter: dict[str, int] = {}
    phylip_names = []
    used: set[str] = set()
    for name, _ in int_seqs:
        prefix = name[:10]
        if prefix_freq[prefix] == 1:
            phylip_names.append(prefix)
            used.add(prefix)
        else:
            idx = prefix_counter.get(prefix, 0)
            prefix_counter[prefix] = idx + 1
            suffix = f"{idx:02d}"
            candidate = prefix[:8] + suffix
            # If still collides (very unlikely), keep incrementing
            while candidate in used:
                idx += 1
                prefix_counter[prefix] = idx + 1
                suffix = f"{idx:02d}"
                candidate = prefix[:8] + suffix
            phylip_names.append(candidate)
            used.add(candidate)
    n_collisions = len(prefix_counter)
    if n_collisions:
        print(f"WARNING: {n_collisions} name prefix collision(s) disambiguated with 2-digit suffix")

    with open(out_path, "w") as out:
        out.write(f" {nseq} {nsite}\n")
        for (name, states), phylip_name in zip(int_seqs, phylip_names):
            padded_name = phylip_name.ljust(10)
            state_str = " ".join(str(s) for s in states)
            out.write(f"{padded_name} {state_str}\n")

    print(f"Written: {out_path}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python3 convert_3diaa.py <aa.fa> <ss.fa> <output.phy>")
        sys.exit(1)
    convert(sys.argv[1], sys.argv[2], sys.argv[3])
