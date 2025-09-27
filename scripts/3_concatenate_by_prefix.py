import os
import argparse
import pandas as pd
from collections import defaultdict
import re

def concatenate_by_prefix(input_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    groups = defaultdict(list)

    # Group files by prefix before the last underscore block
    for fname in os.listdir(input_dir):
        if fname.endswith('.tsv'):
            match = re.match(r'(.+?\.fasta)_.+\.tsv$', fname)
            if match:
                prefix = match.group(1)  # e.g., MGBC000001.fasta
                groups[prefix].append(os.path.join(input_dir, fname))

    for prefix, files in groups.items():
        dfs = []
        for file in files:
            try:
                df = pd.read_csv(file, sep='\t', header=None)
                dfs.append(df)
            except Exception as e:
                print(f"Skipping {file}: {e}")
        if dfs:
            result = pd.concat(dfs, ignore_index=True)
            out_file = os.path.join(output_dir, f"{prefix}.tsv")
            result.to_csv(out_file, sep='\t', index=False, header=False)
            print(f"Saved: {out_file}")
        else:
            print(f"No data found for prefix: {prefix}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Concatenate TSV files by shared filename prefix.")
    parser.add_argument("-i", "--input_dir", required=True, help="Directory containing the TSV files")
    parser.add_argument("-o", "--output_dir", required=True, help="Directory to save merged TSV files")

    args = parser.parse_args()
    concatenate_by_prefix(args.input_dir, args.output_dir)
