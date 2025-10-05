#!/usr/bin/env python3

import pandas as pd
import hashlib
import argparse

# -------------------------------
# Command-line argument handling
# -------------------------------
parser = argparse.ArgumentParser(
    description="Generate Kraken2-compatible taxonomy files (seqid2taxid.map, names.dmp, nodes.dmp) from GTDB-style taxonomy."
)
parser.add_argument(
    "-i", "--input", required=True, help="Input taxonomy file (e.g., taxonomy.txt or silva_taxa_to_GTDB.tsv)"
)
parser.add_argument(
    "-o", "--outprefix", default="gtdb_taxonomy",
    help="Prefix for output files (default: gtdb_taxonomy). Output files will be named <prefix>_seqid2taxid.map, <prefix>_names.dmp, and <prefix>_nodes.dmp."
)
args = parser.parse_args()

input_file = args.input
output_seqid2taxid = f"{args.outprefix}_seqid2taxid.map"
output_nodes = f"{args.outprefix}_nodes.dmp"
output_names = f"{args.outprefix}_names.dmp"

# -------------------------------
# Define taxonomic structure
# -------------------------------
taxonomic_ranks = ["no rank", "domain", "phylum", "class", "order", "family", "genus", "species"]
rank_prefixes = ["root__", "d__", "p__", "c__", "o__", "f__", "g__", "s__"]

# -------------------------------
# Load taxonomy input file
# -------------------------------
taxonomy = pd.read_csv(input_file, sep="\t", header=None, names=["seq_id", "taxonomy"])

# -------------------------------
# Helper functions
# -------------------------------
def generate_taxid(taxonomy_string):
    """Generate a reproducible numeric taxid."""
    return int(hashlib.md5(taxonomy_string.encode()).hexdigest(), 16) % 1000000000

def parse_taxonomy(tax_string):
    """Parse a GTDB taxonomy string into individual ranks."""
    taxa = {prefix: "" for prefix in rank_prefixes}
    for entry in tax_string.split(";"):
        for prefix in rank_prefixes:
            if entry.startswith(prefix):
                taxa[prefix] = entry.replace(prefix, "")
    return taxa

# -------------------------------
# Initialize taxonomy containers
# -------------------------------
taxon_dict = {}   # taxid → name
parent_dict = {}  # taxid → parent_taxid
rank_dict = {}    # taxid → rank
seqid_to_taxid = []

# Add root node
taxon_dict[1] = "root"
parent_dict[1] = 0
rank_dict[1] = "no rank"

# -------------------------------
# Build taxonomy tree
# -------------------------------
for _, row in taxonomy.iterrows():
    seq_id = row["seq_id"]
    taxonomy_string = row["taxonomy"]
    taxa = parse_taxonomy(taxonomy_string)

    last_taxid = 1  # root as parent

    for i, prefix in enumerate(rank_prefixes[1:], start=1):  # skip root
        tax_name = taxa[prefix]
        if tax_name:
            taxid = generate_taxid(f"{prefix}{tax_name}")

            if taxid not in taxon_dict:
                formatted_name = f"{taxonomic_ranks[i]}__{tax_name}"
                taxon_dict[taxid] = formatted_name
                parent_dict[taxid] = last_taxid
                rank_dict[taxid] = taxonomic_ranks[i]

            last_taxid = taxid

    # Skip root-only (unclassified) sequences
    if last_taxid != 1:
        seqid_to_taxid.append([f"{seq_id}|kraken:taxid|{last_taxid}", last_taxid])

# -------------------------------
# Write seqid2taxid.map
# -------------------------------
pd.DataFrame(seqid_to_taxid).to_csv(output_seqid2taxid, sep="\t", index=False, header=False)

# -------------------------------
# Write names.dmp
# -------------------------------
with open(output_names, "w") as f:
    for taxid, name in taxon_dict.items():
        split_name = name.split("__", 1)
        formatted_name = split_name[1] if len(split_name) > 1 else name
        f.write(f"{taxid}\t|\t{formatted_name}\t|\t{name}\t|\tscientific name\t|\n")

# -------------------------------
# Write nodes.dmp
# -------------------------------
with open(output_nodes, "w") as f:
    for taxid, parent_taxid in parent_dict.items():
        rank = rank_dict[taxid]
        f.write(f"{taxid}\t|\t{parent_taxid}\t|\t{rank}\t|\n")

print(f"✅ Files generated successfully:\n- {output_seqid2taxid}\n- {output_names}\n- {output_nodes}")
