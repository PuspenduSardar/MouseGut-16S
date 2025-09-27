#!/bin/bash
# circular_16s_extractor.sh
# Main orchestrator script for extracting full-length 16S rRNA using degenerate primers in circularized MAGs

set -e

# Default parameters
THREADS=8
MIN_LEN=700
MAX_LEN=3000
IDENTITY=70
PRIMER_FILE="primers.txt"

# Help message
usage() {
  echo "Usage: $0 --input_dir DIR --output_dir DIR [--threads N] [--min_len N] [--max_len N] [--identity N] [--primers FILE]"
  exit 1
}

# Parse arguments
while [[ "$#" -gt 0 ]]; do
  case $1 in
    --input_dir) INPUT_DIR="$2"; shift ;;
    --output_dir) OUTPUT_DIR="$2"; shift ;;
    --threads) THREADS="$2"; shift ;;
    --min_len) MIN_LEN="$2"; shift ;;
    --max_len) MAX_LEN="$2"; shift ;;
    --identity) IDENTITY="$2"; shift ;;
    --primers) PRIMER_FILE="$2"; shift ;;
    *) usage ;;
  esac
  shift
done

# Validate
if [[ -z "$INPUT_DIR" || -z "$OUTPUT_DIR" ]]; then
  usage
fi

mkdir -p "$OUTPUT_DIR"

# Step 1: Concatenate contigs
echo "🔄 Concatenating multi-contig MAGs..."
python3 1_1_concatenate_contigs.py "$INPUT_DIR" "$OUTPUT_DIR/concatenated"

# Step 2: Expand degenerate primers
echo "🧬 Expanding degenerate primers..."
python3 1_2_expand_primers.py "$PRIMER_FILE" "$OUTPUT_DIR/primers_expanded.txt"

# Step 3: BLAST primer hits
echo "🚀 Running BLAST searches..."
python3 1_3_run_blast_primers.py "$OUTPUT_DIR/concatenated" "$OUTPUT_DIR/primers_expanded.txt" "$OUTPUT_DIR/blast_hits" "$THREADS" "$IDENTITY"

echo "✅ Done! All Blast hits are stored in $OUTPUT_DIR/blast_hits"
