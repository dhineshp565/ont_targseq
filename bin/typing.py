#!/usr/bin/env python3
"""
Sequence Typing Script

This script performs sequence typing by:
1. Cleaning FASTA headers (removing _consensus/_medaka suffixes)
2. Running Abricate to identify sequences against a custom database
3. Extracting matched sequences and creating output files
4. Combining Abricate results with sequence data for importing into LIMS

Usage: typing.py <sample_name> <input_fasta> <abricate_datadir>

Outputs:
  - {sample}_abricate.tsv: Raw Abricate results
  - {sample}_withseq.tsv: Abricate results with sequences appended
  - {sample}_filtered.fasta: FASTA containing only matched sequences
  - {sample}.tsv: TSV format of matched sequences
  - {sample}_matched_seqs.txt: List of matched sequence names
"""

import sys
import subprocess
from pathlib import Path


# Parse command-line arguments

sample = sys.argv[1]        # Sample name for output files
fasta = Path(sys.argv[2])   # Input FASTA file (consensus sequences)
datadir = sys.argv[3]       # Abricate database directory

# Define output file paths
abricate_file = Path(f"{sample}_abricate.tsv")       # Abricate raw output
filtered_fasta = Path(f"{sample}_filtered.fasta")    # FASTA with only matched sequences
withseq_file = Path(f"{sample}_withseq.tsv")         # Abricate results + sequences
cleaned_fasta = Path(f"{sample}.fasta")              # Input FASTA with cleaned headers


# Clean input FASTA headers
# Remove _consensus and _medaka suffixes from sequence headers
# This ensures consistency with Abricate database entries

with open(fasta) as infile, open(cleaned_fasta, "w") as outfile:
    for line in infile:
        if line.startswith(">"):
            # Strip pipeline-added suffixes from headers
            line = line.replace("_medaka_consensus", "")
        outfile.write(line)


# Run Abricate for sequence typing
# Search sequences against custom database with 60% identity and coverage thresholds

cmd = [
    "abricate",
    "--datadir", datadir,
    "--db", "targseq",          # Custom database name
    "-minid", "60",             # Minimum 60% sequence identity
    "-mincov", "60",            # Minimum 60% coverage
    "--quiet",                  # Suppress progress messages
    str(cleaned_fasta)
]

with open(abricate_file, "w") as outfile:
    subprocess.run(
        cmd,
        stdout=outfile,
        check=True  # Raise exception if Abricate fails
    )


# Read Abricate results

abricate_results = abricate_file.read_text().splitlines()


# Handle case with no Abricate matches
# If only header line exists in abricate output (no hits), create placeholder outputs

if len(abricate_results) < 2:

    # Write Abricate header plus placeholder row with NA values
    with open(withseq_file, "w") as f:
        f.write(f"{abricate_results[0]}\tSEQ HEADER\tSEQUENCE\n")
        f.write(
            f"{sample}.fasta\t{sample}_No_matches_found\t"
            "NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t"
            "NA\tNA\tNA\tNA\tNo_consensus\tNA\tNA\n"
        )

    # Create placeholder FASTA indicating no matches
    filtered_fasta.write_text(f">{sample}\nNo_matches_found")

    sys.exit(0)  # Exit early since no further processing needed


# Extract sequence names from Abricate results
# Field [1] in Abricate output contains the sequence ID

matched_sequences = set() # unique sequence IDs

for line in abricate_results[1:]:  # Skip header line
    fields = line.split("\t")
    matched_sequences.add(fields[1])  # Collect unique sequence IDs



# Parse FASTA and extract only matched sequences
# Only sequences that had Abricate hits are retained

sequences = []  # List of (header, sequence) tuples

current_header = None
current_sequence = []

with open(cleaned_fasta) as f:

    for line in f:
        line = line.strip()

        if not line:
            continue

        if line.startswith(">"):

            # Save previous sequence if it matched
            if current_header is not None:
                if current_header in matched_sequences:
                    sequence = "".join(current_sequence) # Concatenate sequence lines into single line
                    sequences.append(
                        (current_header, sequence) #tuple of header and sequence
                    )

            # Start new sequence
            current_header = line[1:]  # Remove '>' prefix
            current_sequence = []

        else:
            # Accumulate sequence lines
            current_sequence.append(line)


# Save the final sequence if it matched
if current_header is not None:
    if current_header in matched_sequences:
        sequences.append(
            (current_header, "".join(current_sequence))
        )

# Create filtered FASTA output
# Contains only sequences that matched the Abricate database

with open(filtered_fasta, "w") as f:

    for header, sequence in sequences:
        f.write(f">{header}\n")
        f.write(f"{sequence}\n")


# Combine Abricate results with sequence data
# Create final output file with Abricate annotations plus full sequences

# Build lookup dictionary for fast sequence retrieval
sequence_dict = {
    header: sequence
    for header, sequence in sequences
}

with open(withseq_file, "w") as output:

    # Write header: Abricate columns + SEQ HEADER + SEQUENCE
    output.write(
        f"{abricate_results[0]}\tSEQ HEADER\tSEQUENCE\n"
    )

    # Append sequence data to each Abricate result row
    for line in abricate_results[1:]:

        fields = line.split("\t")

        seq_header = fields[1]

        if seq_header in sequence_dict:

                sequence = sequence_dict[seq_header]

                # Write complete row with sequence appended
                output.write(
                    f"{line}\t{seq_header}\t{sequence}\n"
                )

