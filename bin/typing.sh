#!/usr/bin/env bash

# This script performs typing analysis using the abricate tool and processes the results.
# Usage: typing.sh <SampleName> <consensus_sequence> <abricate_database_directory>
# Arguments:
#   $1 - SampleName: The name of the sample being analyzed.
#   $2 - consensus_sequence: The path to the consensus sequence file.
#   $3 - abricate_database_directory: The directory containing the abricate database.


# Run abricate with specified database and parameters, output results to a CSV file named <SampleName>_abricate.csv
abricate --datadir $3 --db targseq -minid 60  -mincov 60 --quiet $2 1> $1_abricate.csv

# Replace '_consensus' in the output CSV file with an empty string
sed -i -e "s/_consensus//g" -e "s/_medaka//g" "$1_abricate.csv"

# Check if the output CSV file has less than 2 lines, if so, append a line indicating no consensus and none for all other fields
if [[ $(wc -l < "$1_abricate.csv") -lt 2 ]]
then
	printf "$1.fasta\t$1_No_consensus\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNo_consensus" >> $1_abricate.csv	
fi

# Extract sequence names that have abricate results (skip header line, get column 2)
tail -n +2 $1_abricate.csv | cut -f2 > $1_matched_seqs.txt

# Extract only sequences that matched in abricate results and convert to TSV
awk 'NR==FNR {matched[$1]=1; next}     # First pass: read matched sequence names and store in array
     /^>/ {                            # When we encounter a fasta header line
       if (seq && print_seq) print header "\t" seq    # Print previous sequence if it was matched
       header=substr($0,2)                            # Extract header name (remove leading >)
       gsub(/_consensus/,"",header)                   # Remove _consensus from header
       gsub(/_medaka/,"",header)                      # Remove _medaka from header
       print_seq=(header in matched)                  # Check if this header is in matched list
       seq=""                                         # Reset sequence variable
       next                                           # Move to next line
     }
     print_seq {seq=seq $0}            # If current sequence is matched, accumulate sequence lines
     END {if (seq && print_seq) print header "\t" seq}' $1_matched_seqs.txt $2 > $1.tsv

# Create a filtered FASTA file with only matched sequences
awk 'NR==FNR {matched[$1]=1; next}     # Load matched sequence names
     /^>/ {                            # When we encounter a fasta header line
       if (seq && print_seq) {print ">"header; print seq}    # Print previous sequence in FASTA format
       header=substr($0,2)              # Extract header name
       gsub(/_consensus/,"",header)     # Remove _consensus from header
       gsub(/_medaka/,"",header)        # Remove _medaka from header
       print_seq=(header in matched)    # Check if this header is in matched list
       seq=""                           # Reset sequence variable
       next
     }
     print_seq {seq=seq $0}            # Accumulate sequence lines
     END {if (seq && print_seq) {print ">"header; print seq}}' $1_matched_seqs.txt $2 > $1_filtered.fasta

# Add a header to the TSV file with sequence header and sequence columns
sed -i '1i SEQ HEADER\tSEQUENCE' $1.tsv

# Combine the abricate results CSV and the sequence TSV into a single CSV file named <SampleName>_withseq.csv
paste $1_abricate.csv $1.tsv > $1_withseq.csv