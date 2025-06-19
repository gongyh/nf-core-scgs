#!/usr/bin/env python

import re
import argparse

def update_paths_file(fasta_path: str, paths_path: str, output_path: str) -> None:
    """
    Update contig IDs in .paths file based on corresponding IDs from .fasta file
    
    Args:
        fasta_path (str): Path to input contigs.fasta file
        paths_path (str): Path to input contigs.paths file
        output_path (str): Path for updated output file
    """
    # Step 1: Create NODE number to full ID mapping from FASTA
    node_id_map = {}
    with open(fasta_path, 'r') as f:
        for line in f:
            if line.startswith('>'):  # Process FASTA header lines
                full_id = line.strip()[1:]  # Remove '>' and newline
                components = full_id.split('_')
                node_number = components[1]  # Extract NODE number (e.g., '382')
                node_id_map[node_number] = full_id  # Map NODE number to full new ID

    # Step 2: Update paths file using the ID mapping
    with open(paths_path, 'r') as in_file, open(output_path, 'w') as out_file:
        # Regex pattern to match old-style IDs in paths file
        id_pattern = re.compile(r'NODE_(\d+)_length_\d+_cov_[\d.]+')

        for line in in_file:
            # Replacement function (removed re.Match type hint for compatibility)
            def replace_callback(match):
                matched_node = match.group(1)  # Get captured NODE number
                return node_id_map.get(matched_node, match.group(0))  # Use mapping or keep original

            # Replace all matching IDs in current line
            updated_line = id_pattern.sub(replace_callback, line)
            out_file.write(updated_line)

if __name__ == "__main__":
    # Set up command-line argument parser
    parser = argparse.ArgumentParser(description='Update contig IDs in .paths file based on .fasta file IDs')
    parser.add_argument('-f', '--fasta', required=True, help='Path to contigs.fasta file')
    parser.add_argument('-p', '--paths', required=True, help='Path to contigs.paths file')
    parser.add_argument('-o', '--output', required=True, help='Path for updated output file')
    
    args = parser.parse_args()
    update_paths_file(args.fasta, args.paths, args.output)
    print(f"Successfully updated paths file. Output saved to: {args.output}")
