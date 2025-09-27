#!/usr/bin/env python3
import os
import sys
import subprocess
from Bio import SeqIO

genome_dir = sys.argv[1]
primer_file = sys.argv[2]
blast_dir = sys.argv[3]
threads = sys.argv[4]
identity = sys.argv[5]

os.makedirs(blast_dir, exist_ok=True)

# Read primers
primers = []
with open(primer_file) as f:
    for line in f:
        name, seq = line.strip().split()
        primers.append((name, seq))

# For each genome
for file in os.listdir(genome_dir):
    if not file.endswith(".fasta"):
        continue
    genome_path = os.path.join(genome_dir, file)
    db_name = genome_path + ".blastdb"
    if not os.path.exists(db_name + ".nin"):
        subprocess.run(["makeblastdb", "-in", genome_path, "-dbtype", "nucl"], check=True)
    for pname, pseq in primers:
        query_file = os.path.join(blast_dir, f"{pname}.fa")
        with open(query_file, "w") as qf:
            qf.write(f">{pname}\n{pseq}\n")
        out_file = os.path.join(blast_dir, f"{file}_{pname}.tsv")
        subprocess.run([
            "blastn", "-query", query_file, "-db", genome_path,
            "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
            "-num_threads", threads, "-perc_identity", identity, "-task", "blastn-short", "-evalue", "0.01",
            "-out", out_file
        ], check=True)
