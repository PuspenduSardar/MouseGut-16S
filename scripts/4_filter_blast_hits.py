#!/usr/bin/env python3

import os
import argparse
import pandas as pd

# Define primer groups
FORWARD_PRIMERS = ("27F")
REVERSE_PRIMERS = ("1492R", "1510R")

def filter_hits(input_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    for fname in os.listdir(input_dir):
        if not fname.endswith(".tsv"):
            continue

        in_path = os.path.join(input_dir, fname)
        out_path = os.path.join(output_dir, fname)

        try:
            df = pd.read_csv(in_path, sep='\t', header=None)
            if df.empty:
                df.to_csv(out_path, sep='\t', index=False, header=False)
                continue

            forward = df[df[0].str.startswith(FORWARD_PRIMERS)]
            reverse = df[df[0].str.startswith(REVERSE_PRIMERS)]

            forward_zero = forward[forward[4] == 0]
            reverse_zero = reverse[reverse[4] == 0]

            if not forward_zero.empty and not reverse_zero.empty:
                filtered_df = pd.concat([forward_zero, reverse_zero])
            else:
                filtered_df = df  # fallback to original if either is missing

            filtered_df.to_csv(out_path, sep='\t', index=False, header=False)
        except Exception as e:
            print(f"Error processing {fname}: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_dir", required=True, help="Input directory of raw BLAST TSVs")
    parser.add_argument("-o", "--output_dir", required=True, help="Directory to save filtered TSVs")
    args = parser.parse_args()
    filter_hits(args.input_dir, args.output_dir)

