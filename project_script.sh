#!/bin/bash

# Function to display script usage
display_usage() {
    echo """Usage: $0 [install | help | -f <IDs> -p <path> [-t <table>] [-m <minsize>]]
Welcome to our script. In this script, you must enter the required IDs
and the location to store the results. We will produce 5 files:
1) sequence.fasta: Contain the downloaded sequences
2) ORFs_out.txt: Contain the ORF Finder results.
3) Conserved_region.fasta: Contain the conserved nucleotides and its length.
4) Main_conserved.fasta: Contain The Longest conserved region.
5) Alignment.fasta: Contain the alignment result of the given IDs

To use the script, you must have the following packages:
- python
- Biopython
- Numpy
- pip
To install all packages, run: bash $0 install

Flags:
-h       Show help message
-f       IDs of the files, example: AB021961.1,U50395.1
-p       Path of the desired folder to store the results, example: /home/bioinformaticsnu
-t       Table used in ORF (default=0)
-m       Minimum size of ORF (default=300)"""
    exit 1
}

# Check for the correct number of arguments
if [ "$#" -lt 1 ]; then
    display_usage
fi

# Part 1: Install Dependencies
if [[ "$1" == "install" ]]; then
    sudo apt update
    sudo apt install -y python3-pip
    pip3 install numpy biopython
    sudo apt-get install -y python3-matplotlib
    echo "Congratulations! Dependencies are installed."
    exit 0
fi

# Part 2: Display Help Information
if [[ "$1" == "help" ]]; then
    display_usage
fi

while getopts ":f:p:t:m:h" flag; do
    case "${flag}" in
        f) ID=${OPTARG} ;;
        p) file_path=${OPTARG} ;;
        t) table=${OPTARG:-0} ;;  # Use default value if not provided
        m) minsize=${OPTARG:-300} ;;  # Use default value if not provided
        h) display_usage ;;
        ?) echo "Invalid option: -${OPTARG}" && display_usage ;;
    esac
done

efetch -db nucleotide -format fasta -id $ID > $file_path/sequence.fasta
cut -d ' ' -f1 $file_path/sequence.fasta > $file_path/Sequence.fasta
rm $file_path/sequence.fasta
muscle -in $file_path/Sequence.fasta -out $file_path/Alignment.fasta  \
 -clw -tree1 $file_path/Tree.phy
#! /usr/bin/python
python3 Conserved_region_linux.py $file_path
getorf -sequence $file_path/Main_conserved.fasta -outseq $file_path/ORFs_out.txt -table $table -minsize $minsize -find 3
