#!/usr/bin/env python3
import sys
import re  # Added import
from Bio import SeqIO

def find_motif_sites(fasta_path, bed_path, motif):
    """Find regex patterns in FASTA and create BED entries."""
    pattern = re.compile(motif, re.IGNORECASE)  # Regex with case insensitivity
    
    with open(bed_path, 'w') as bed_file:
        for record in SeqIO.parse(fasta_path, "fasta"):
            seq = str(record.seq)  # Remove .upper() to preserve case
            for match in pattern.finditer(seq):
                bed_entry = (record.id, match.start(), match.end())
                bed_file.write("\t".join(map(str, bed_entry)) + "\n")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python find_motif_sites.py <input.fasta> <output.bed> <motif>")
        sys.exit(1)
    fasta_path = sys.argv[1]
    bed_path = sys.argv[2]
    motif = sys.argv[3]
    find_motif_sites(fasta_path, bed_path, motif)
