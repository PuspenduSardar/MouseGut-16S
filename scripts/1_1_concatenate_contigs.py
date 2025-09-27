#!/usr/bin/env python3
import os
import sys
from Bio import SeqIO

input_dir = sys.argv[1]
output_dir = sys.argv[2]
os.makedirs(output_dir, exist_ok=True)

def clean_id(name):
    return os.path.splitext(os.path.basename(name))[0]

for filename in os.listdir(input_dir):
    if filename.endswith(('.fna', '.fa', '.fasta')):
        full_path = os.path.join(input_dir, filename)
        genome_id = clean_id(filename)
        sequences = list(SeqIO.parse(full_path, "fasta"))
        concatenated_seq = "".join(str(seq.seq) for seq in sequences)
        new_record = SeqIO.SeqRecord(
            seq=concatenated_seq.upper(),
            id=genome_id,
            description=""
        )
        output_path = os.path.join(output_dir, f"{genome_id}.fasta")
        with open(output_path, "w") as out_f:
            SeqIO.write(new_record, out_f, "fasta")
