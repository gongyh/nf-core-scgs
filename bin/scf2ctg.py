#!/usr/bin/env python

import sys
from Bio import SeqIO


def split_scaffolds_into_contigs(input_file, output_file):
    # Parse the input FASTA file
    records = list(SeqIO.parse(input_file, "fasta"))

    # Create a new file to write the contigs
    with open(output_file, "w") as output_handle:
        for record in records:
            # Split the scaffold into contigs based on 'N' gaps
            contigs = record.seq.split("N")
            # Remove any empty strings that might result from leading/trailing gaps
            contigs = [contig for contig in contigs if contig]

            if len(contigs) == 1:
                # If there is only one contig (no 'N' in the scaffold), keep the original ID
                SeqIO.write(record, output_handle, "fasta")
            else:
                for i, contig in enumerate(contigs, start=1):
                    # Create a new SeqRecord for each contig
                    new_record = SeqIO.SeqRecord(
                        seq=contig, id=f"{record.id}_contig{i}", description=record.description
                    )
                    # Write the new SeqRecord to the output file
                    SeqIO.write(new_record, output_handle, "fasta")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python scf2ctg.py input.fasta output.fasta")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    split_scaffolds_into_contigs(input_file, output_file)
