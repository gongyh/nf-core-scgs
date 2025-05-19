#!/usr/bin/env python

from Bio import SeqIO
import sys


def update_fasta_id_with_length(input_file):
    records = list(SeqIO.parse(input_file, "fasta"))
    for record in records:
        parts = record.id.split("_")
        if "length" in parts:
            actual_length = len(record)
            length_index = parts.index("length")
            parts[length_index + 1] = str(actual_length)
            new_id = "_".join(parts)
            record.id = new_id
            record.description = ""

    SeqIO.write(records, sys.stdout, "fasta")


if len(sys.argv) != 2:
    print("Usage: python fixSPAdesLen.py <input_file>")
    sys.exit(1)

input_fasta = sys.argv[1]
update_fasta_id_with_length(input_fasta)
