#!/usr/bin/env python3
import pandas as pd
import os
import argparse
import re
from collections import defaultdict

def clean_filename(name):
    clean = re.sub(r'[^a-zA-Z0-9_\-\.]', '_', name)
    clean = clean.strip('-_')
    if not clean or not clean[0].isalpha():
        clean = "target_" + clean
    return clean

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--checkm2', required=True)
    parser.add_argument('--taxonomy', required=True)
    parser.add_argument('--bindir', required=True)
    parser.add_argument('--outdir', required=True)
    parser.add_argument('--min_comp', type=float, default=50.0)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print(f"--- 1. Reading CheckM2 report ---")
    df_checkm = pd.read_csv(args.checkm2, sep='\t', comment='#')
    name_col = 'Name' if 'Name' in df_checkm.columns else df_checkm.columns[0]
    comp_col = 'Completeness' if 'Completeness' in df_checkm.columns else df_checkm.columns[1]

    # Record quality of all bins
    bin_quality = dict(zip(df_checkm[name_col], df_checkm[comp_col]))
    low_quality_bins = df_checkm[df_checkm[comp_col] < args.min_comp][name_col].tolist()
    print(f"Found {len(low_quality_bins)} low-quality bins (Completeness < {args.min_comp}%)")

    print(f"\n--- 2. Reading Taxonomy report ---")
    df_tax = pd.read_csv(args.taxonomy, sep='\t', header=None, names=['contig', 'taxonomy'], skiprows=1)
    contig2tax = dict(zip(df_tax['contig'], df_tax['taxonomy']))

    print(f"\n--- 3. Determining dominant species for all bins ---")
    all_bin2species = {}
    for bin_name in df_checkm[name_col]:
        bin_path = os.path.join(args.bindir, f"{bin_name}.fa")
        if not os.path.exists(bin_path):
            continue

        species_counts = defaultdict(int)
        with open(bin_path, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    contig_id = line[1:].strip().split()[0]
                    tax = contig2tax.get(contig_id, "Unknown")
                    species_counts[tax] += 1

        if species_counts:
            dominant_species = max(species_counts, key=species_counts.get)
            # Fix Bug: Check if species level exists (single underscore ";s_")
            if ";s_" in dominant_species:
                all_bin2species[bin_name] = dominant_species

    print(f"\n--- 4. Smart merging of same-species fragments (including high/low quality)---")
    species2bins = defaultdict(list)
    for b, s in all_bin2species.items():
        species2bins[s].append(b)

    merged_count = 0
    for species, bins in species2bins.items():
        if len(bins) > 1:
            # Core logic: If ANY bin in this species group is low-quality, trigger merging!
            has_low_quality = any(b in low_quality_bins for b in bins)

            if has_low_quality:
                safe_species_name = clean_filename(species.split(";")[-1])
                target_file = os.path.join(args.outdir, f"{safe_species_name}.fa")

                print(f"✅ Merging: {bins} -> {target_file}")
                print(f"   (Because it contains low-quality bin(s):  {[b for b in bins if b in low_quality_bins]})")

                seen_contigs = set()
                with open(target_file, 'w') as out_f:
                    for b in bins:
                        bin_path = os.path.join(args.bindir, f"{b}.fa")
                        with open(bin_path, 'r') as in_f:
                            write_seq = False
                            for line in in_f:
                                if line.startswith('>'):
                                    contig_id = line[1:].strip().split()[0]
                                    if contig_id not in seen_contigs:
                                        seen_contigs.add(contig_id)
                                        write_seq = True
                                        out_f.write(line)
                                    else:
                                        write_seq = False
                                elif write_seq:
                                    out_f.write(line)
                merged_count += 1

    print(f"\n--- Finished ---")
    print(f"Total {merged_count} target fasta file(s) created.")

if __name__ == '__main__':
    main()
