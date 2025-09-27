def replace_fasta_headers(fasta_file, mapping_file, output_file):
    # Read the mapping file and create a dictionary
    header_map = {}
    with open(mapping_file, 'r') as map_file:
        for line in map_file:
            parts = line.strip().split('|', 1)  # Split at first pipe
            if len(parts) > 1:
                header_map[parts[0]] = line.strip()
    
    # Process the FASTA file and replace headers
    with open(fasta_file, 'r') as fasta, open(output_file, 'w') as output:
        for line in fasta:
            if line.startswith('>'):
                header = line[1:].strip()  # Remove '>' and newline
                new_header = header_map.get(header, header)  # Replace if found
                output.write(f'>{new_header}\n')  # Write new header
            else:
                output.write(line)  # Write sequence lines unchanged
    
    print(f'Headers replaced successfully. Output written to {output_file}')

# Example usage
replace_fasta_headers('fasta4_krakenbracken.fasta', 'seqid2taxid_mapping.txt', 'output.fasta')
