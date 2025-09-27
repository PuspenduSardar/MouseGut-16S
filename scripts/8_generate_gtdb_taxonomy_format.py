import pandas as pd
import hashlib

# Input and output filenames
input_file = "taxonomy4kraken_bracken.txt"
output_seqid2taxid = "seqid2taxid.map"
output_nodes = "nodes.dmp"
output_names = "names.dmp"

# Define taxonomic ranks and prefixes
taxonomic_ranks = ["no rank", "domain", "phylum", "class", "order", "family", "genus", "species"]
rank_prefixes = ["root__", "d__", "p__", "c__", "o__", "f__", "g__", "s__"]

# Load GTDB taxonomy file
taxonomy = pd.read_csv(input_file, sep="\t", header=None, names=["seq_id", "taxonomy"])

# Function to generate a unique taxonomic ID based on the taxonomy string
def generate_taxid(taxonomy_string):
    return int(hashlib.md5(taxonomy_string.encode()).hexdigest(), 16) % 1000000000  # Keeps it within numeric limits

# Function to parse taxonomy into structured ranks
def parse_taxonomy(tax_string):
    taxa = {prefix: "" for prefix in rank_prefixes}
    for entry in tax_string.split(";"):
        for prefix in rank_prefixes:
            if entry.startswith(prefix):
                taxa[prefix] = entry.replace(prefix, "")
    return taxa

# Dictionaries to store taxon data
taxon_dict = {}  # Maps tax_id → tax_name
parent_dict = {}  # Maps tax_id → parent_tax_id
rank_dict = {}  # Maps tax_id → rank
seqid_to_taxid = []  # Stores sequence ID mappings

# Root node (ID = 1)
taxon_dict[1] = "root"
parent_dict[1] = 0  # Root has no parent
rank_dict[1] = "no rank"

# Process each row in taxonomy file
for _, row in taxonomy.iterrows():
    seq_id = row["seq_id"]
    taxonomy_string = row["taxonomy"]
    taxa = parse_taxonomy(taxonomy_string)

    last_taxid = 1  # Start from root

    for i, prefix in enumerate(rank_prefixes[1:], start=1):  # Skip root
        tax_name = taxa[prefix]
        if tax_name:
            taxid = generate_taxid(f"{prefix}{tax_name}")
            
            if taxid not in taxon_dict:
                formatted_name = f"{taxonomic_ranks[i]}__{tax_name}"  # Add rank prefix
                taxon_dict[taxid] = formatted_name
                parent_dict[taxid] = last_taxid  # Parent = previous level
                rank_dict[taxid] = taxonomic_ranks[i]

            last_taxid = taxid  # Update last taxid

    # Store sequence ID mapping to species taxid
    seqid_to_taxid.append([seq_id + "|kraken:taxid|" + str(last_taxid), last_taxid])

# Save seqid2taxid.map in the correct format
pd.DataFrame(seqid_to_taxid).to_csv(output_seqid2taxid, sep="\t", index=False, header=False)

# Save names.dmp (correct format)
with open(output_names, "w") as f:
    for taxid, name in taxon_dict.items():
        split_name = name.split('__', 1)
        formatted_name = split_name[1] if len(split_name) > 1 else name  # Ensure no index error
        f.write(f"{taxid}\t|\t{formatted_name}\t|\t{name}\t|\tscientific name\t|\n")

# Save nodes.dmp (correct format)
with open(output_nodes, "w") as f:
    for taxid, parent_taxid in parent_dict.items():
        rank = rank_dict[taxid]
        f.write(f"{taxid}\t|\t{parent_taxid}\t|\t{rank}\t|\n")

print("Files generated: seqid2taxid.map, names.dmp, nodes.dmp")
