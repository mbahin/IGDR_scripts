#!/bin/bash

# Mathieu Bahin, 12/06/14

# Script to get the sequences involved in a breakpoint from a fasta file produced by chimCT
# Inputs are a fasta output file from chimCT, a breakpoint coordinates (e.g. "(15:60179225,1 / X:35772912,1 # 50,62)") and an output filename.
# The output is fasta file.

# Getting options back
while getopts "f:b:o:" OPTION
do
	case $OPTION in
    	f) fasta_file=$OPTARG;;
    	b) breakpoint=$OPTARG;;
    	o) output_file=$OPTARG;;
    esac
done

# Checking the parameters
if [[ -s "$output_file" ]]; then
	echo "Warning: There is already a non-empty file with the output filename provided."
fi

# Getting the sequence ids
feat1=$(echo $breakpoint | sed 's|^(\([^:]*\):\([^,]*\),.*|\1@-?1@\2|g')
feat2=$(echo $breakpoint | sed 's|^.*/ \([^:]*\):\([^,]*\),.*|\1@-?1@\2|g')

# Building the output file
egrep -A 1 "$feat1" $fasta_file | egrep -A 1 "$feat2" >> $output_file