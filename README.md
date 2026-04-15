## Full-length 16S Extractor

### Purpose

This tool extracts full-length 16S-like regions from metagenome assembled genomes (MAGs) using degenerate primer pairs and BLAST-based mapping.

### Requirements

- Python 3 + Biopython
- BLAST+ (`makeblastdb`, `blastn`)
- Bash

### Usage

#### Step 1: Expand primers
`bash scripts/1_0_circular_16s_extraction.sh`

#### Step 2: Subset only non-zero TSV files (optional)
`bash scripts/2_copy_nonZero_tsvFiles.sh`

#### Step 3: Concatenate results from same MAG into one TSV
`python scripts/3_concatenate_by_prefix.py`

#### Step 4: Filter BLAST hits
`python scripts/4_filter_blast_hits.py`

#### Step 5: Calculate amplicon's properties (length, orientation, start-end position etc.)
`python scripts/5_calculate_amplicons.py`

#### Step 6: Extract amplicon sequences of required lengths
`python scripts/6_extract_amplicon_seqs.py`

#### Step 7: Remove any Ns from the fasta file (and filter for the lengths)
`python scripts/7_fasta_cleaner.py`



### 📜 License

MIT License
