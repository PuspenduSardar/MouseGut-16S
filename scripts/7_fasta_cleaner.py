import argparse
from Bio import SeqIO
from Bio.Seq import Seq

def clean_and_filter_fasta(input_file, output_file, min_len, max_len):
    total = 0
    kept = 0
    with open(output_file, "w") as out_handle:
        for record in SeqIO.parse(input_file, "fasta"):
            total += 1
            cleaned_seq = str(record.seq).replace("N", "").replace("n", "")
            if min_len <= len(cleaned_seq) <= max_len:
                record.seq = Seq(cleaned_seq)
                SeqIO.write(record, out_handle, "fasta")
                kept += 1
    print(f"Total sequences processed: {total}")
    print(f"Sequences kept after filtering: {kept}")

def main():
    parser = argparse.ArgumentParser(description="Clean Ns and filter FASTA sequences by length.")
    parser.add_argument("-i", "--input", required=True, help="Path to input FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Path to output FASTA file")
    parser.add_argument("--min_len", type=int, default=0, help="Minimum sequence length after removing Ns (default: 0)")
    parser.add_argument("--max_len", type=int, default=1000000, help="Maximum sequence length after removing Ns (default: 1,000,000)")
    
    args = parser.parse_args()
    clean_and_filter_fasta(args.input, args.output, args.min_len, args.max_len)

if __name__ == "__main__":
    main()
