#!/bin.bash

input_dir="$1"
out_dir="$2"

mkdir -p "$out_dir"

find "$input_dir" -type f -name '*.tsv' ! -size 0 -exec cp {} "$out_dir" \;
