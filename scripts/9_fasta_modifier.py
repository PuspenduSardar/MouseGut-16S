#!/usr/bin/env python3
import argparse

def replace_fasta_headers(fasta_file, mapping_file, output_file):
    """
    Replace FASTA headers using a mapping file where each line contains:
    old_header|kraken:taxid|XXXX
    """
    # Read the mapping file and create a dictionary
    header_map = {}
    with open(mapping_file, 'r') as map_file:
        for line in map_file:
            parts = line.strip().split('|', 1)  # Split at first '|'
            if len(parts) > 1:
                header_map[parts[0]] = line.strip()

    # Process the FASTA file and replace headers
    replaced_count = 0
    total_headers = 0
    with open(fasta_file, 'r') as fasta, open(output_file, 'w') as output:
        for line in fasta:
            if line.startswith('>'):
                total_headers += 1
                header = line[1:].strip()  # Remove '>' and newline
                new_header = header_map.get(header, header)  # Replace if found
                if new_header != header:
                    replaced_count += 1
                output.write(f'>{new_header}\n')
            else:
                output.write(line)

    print(f"✅ Headers replaced successfully.")
    print(f"   Total headers processed: {total_headers}")
    print(f"   Headers replaced: {replaced_count}")
    print(f"   Output written to: {output_file}")

def main():
    parser = argparse.ArgumentParser(
        description="Replace FASTA headers with Kraken2-style headers using a mapping file."
    )
    parser.add_argument(
        "-f", "--fasta", required=True,
        help="Input FASTA file to modify (e.g., sequences.fasta)"
    )
    parser.add_argument(
        "-m", "--mapping", required=True,
        help="Mapping file (e.g., seqid2taxid.map)"
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Output FASTA file with modified headers (e.g., sequences_kraken.fasta)"
    )

    args = parser.parse_args()
    replace_fasta_headers(args.fasta, args.mapping, args.output)

if __name__ == "__main__":
    main()

