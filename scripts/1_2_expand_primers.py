#!/usr/bin/env python3
import sys
from Bio.Data.IUPACData import ambiguous_dna_values
import itertools

def expand_degenerate(primer):
    bases = [ambiguous_dna_values[base] for base in primer]
    return ["".join(p) for p in itertools.product(*bases)]

input_file = sys.argv[1]
output_file = sys.argv[2]

with open(input_file) as f, open(output_file, "w") as out:
    for line in f:
        if not line.strip(): continue
        name, seq = line.strip().split()
        for i, expanded in enumerate(expand_degenerate(seq)):
            out.write(f"{name}_{i+1}\t{expanded}\n")
