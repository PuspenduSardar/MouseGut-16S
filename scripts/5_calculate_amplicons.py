#!/usr/bin/env python3

import os
import argparse
import pandas as pd
from itertools import product

# Define primer groups
FORWARD_PRIMERS = ("27F", "7F")
REVERSE_PRIMERS = ("1492R", "1510R")

def calculate_amplicons(input_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    for fname in os.listdir(input_dir):
        if not fname.endswith(".tsv"):
            continue

        in_path = os.path.join(input_dir, fname)
        out_name = os.path.splitext(fname)[0] + "_amplicons.tsv"
        out_path = os.path.join(output_dir, out_name)

        try:
            df = pd.read_csv(in_path, sep='\t', header=None)
            if df.empty:
                continue

            genome = df.iloc[0, 1] if not df.empty else "Unknown"

            forward = df[df[0].str.startswith(FORWARD_PRIMERS)]
            reverse = df[df[0].str.startswith(REVERSE_PRIMERS)]

            results = []

            for fwd, rev in product(forward.itertuples(index=False), reverse.itertuples(index=False)):
                # Raw BLAST positions
                f_raw_start, f_raw_end = int(fwd[8]), int(fwd[9])
                r_raw_start, r_raw_end = int(rev[8]), int(rev[9])

                # Normalize: smaller always = start, larger always = end
                f_start, f_end = sorted([f_raw_start, f_raw_end])
                r_start, r_end = sorted([r_raw_start, r_raw_end])

                # Amplicon boundaries
                start = min(f_start, r_start)
                end = max(f_end, r_end)
                length = end - start + 1

                orientation = "Obverse" if f_raw_start < r_raw_start else "Reverse"

                results.append([
                    genome,
                    fwd[0], rev[0],
                    f_start, f_end,
                    r_start, r_end,
                    length,
                    orientation
                ])

            if results:
                out_df = pd.DataFrame(results, columns=[
                    "Genome", "Forward_Primer", "Reverse_Primer",
                    "Forward_start", "Forward_end",
                    "Reverse_start", "Reverse_end",
                    "Amplicon_Length", "Orientation"
                ])
                out_df.to_csv(out_path, sep='\t', index=False)
        except Exception as e:
            print(f"Error processing {fname}: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_dir", required=True, help="Filtered TSV directory")
    parser.add_argument("-o", "--output_dir", required=True, help="Directory to save amplicon metadata")
    args = parser.parse_args()
    calculate_amplicons(args.input_dir, args.output_dir)
