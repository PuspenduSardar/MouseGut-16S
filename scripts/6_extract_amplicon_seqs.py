import os
import argparse
import pandas as pd
from Bio import SeqIO

def extract_amplicons(amplicon_dir, fasta_dir, output_dir, min_len, max_len):
    os.makedirs(output_dir, exist_ok=True)

    # Build a dictionary of available FASTA files keyed by cleaned name
    fasta_map = {}
    for fasta_file in os.listdir(fasta_dir):
        if fasta_file.endswith(".fasta") or fasta_file.endswith(".fa"):
            base = fasta_file.replace("_genomic", "").replace(".fasta", "").replace(".fa", "")
            fasta_map[base] = os.path.join(fasta_dir, fasta_file)

    # Process each amplicon TSV file
    for tsv_file in os.listdir(amplicon_dir):
        if not tsv_file.endswith("_amplicons.tsv"):
            continue

        genome_id = tsv_file.replace("_amplicons.tsv", "").replace(".fasta", "")  # normalize name
        tsv_path = os.path.join(amplicon_dir, tsv_file)

        fasta_path = fasta_map.get(genome_id)
        if not fasta_path:
            print(f"FASTA file not found for {genome_id}, skipping.")
            continue

        df = pd.read_csv(tsv_path, sep='\t')

        # Load FASTA sequences
        genome_record = next(SeqIO.parse(fasta_path, "fasta"))

        # Extract and write amplicons
        for idx, row in df.iterrows():
            amp_len = row['Amplicon_Length']
            if amp_len < min_len or amp_len > max_len:
                continue

            f_start = int(row['Forward_start'])
            f_end = int(row['Forward_end'])
            r_start = int(row['Reverse_start'])
            r_end = int(row['Reverse_end'])

            orientation = row['Orientation']
            if orientation == 'Obverse':
                start = min(f_start, f_end)
                end = max(r_start, r_end)
                seq = genome_record.seq[start-1:end]  # 1-based indexing
            else:
                start = min(r_start, r_end)
                end = max(f_start, f_end)
                seq = genome_record.seq[start-1:end].reverse_complement()

            out_file = os.path.join(output_dir, f"{genome_id}_amplicon_{idx+1}.fasta")
            with open(out_file, "w") as f:
                f.write(f">{genome_id}_amplicon_{idx+1}\n{seq}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract amplicon sequences from genome FASTA files.")
    parser.add_argument("-a", "--amplicon_dir", required=True, help="Directory containing amplicon TSV files")
    parser.add_argument("-f", "--fasta_dir", required=True, help="Directory containing genome FASTA files")
    parser.add_argument("-o", "--output_dir", required=True, help="Directory to save amplicon FASTA files")
    parser.add_argument("--min_len", type=int, default=100, help="Minimum amplicon length")
    parser.add_argument("--max_len", type=int, default=2000, help="Maximum amplicon length")
    args = parser.parse_args()

    extract_amplicons(args.amplicon_dir, args.fasta_dir, args.output_dir, args.min_len, args.max_len)

