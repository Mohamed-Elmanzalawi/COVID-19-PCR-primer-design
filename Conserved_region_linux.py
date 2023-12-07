# Code Done by Mohamed Elsayed Elmanzalawi.
# Import necessary libraries
import numpy as np
from Bio import AlignIO, SeqIO, Phylo
import matplotlib.pyplot as plt
import textwrap
import sys

# Get the filepath from the command line arguments
filepath = sys.argv[1]

# Define file paths for input and output files
filename = f"{filepath}/Alignment.fasta"
Conserved_region_txt = f"{filepath}/Conserved_region.fasta"
Main_conserved = f"{filepath}/Main_conserved.fasta"

# Open files for writing
Main_conserved_open = open(Main_conserved, 'w')
Conserved_region_txt_open = open(Conserved_region_txt, 'w')

# Initialize lists for sequence data and conserved positions
A = []
y = []

# Read the alignment file and store sequences in 'A'
for SeqRecord in AlignIO.read(filename, 'clustal'):
    A.append(SeqRecord.seq)

# Convert 'A' to a numpy array for efficient manipulation
profile = np.array(A)

# Identify conserved regions
difference = True
for x in range(len(A[0])):
    if '-' in profile[:, x]:
        difference = True
    if len(set(profile[:, x])) == 1:
        difference = False
        y.append(x)
    if len(set(profile[:, x])) != 1:
        difference = True
    if difference or x == (len(A[0]) - 1):
        if len(y) == 1:
            # Write single conserved nucleotide
            Conserved_region_txt_open.write('>Conserved_Nucleotide(index=%s) \n%s \n\n' % (y[0], A[0][y[0]]))
            y.clear()
        if len(y) == 0:
            continue
        if len(y) != 1:
            # Write conserved region with start and end indices
            New_seq = textwrap.fill("".join(A[0][y[0]:(y[-1] + 1)]), 70)
            Conserved_region_txt_open.write('>Conserved_region(index=%s:%s)_length=%s \n%s \n\n' % (y[0], y[-1],
                                                                                                (y[-1] - y[0] + 1), New_seq))
            y.clear()

# Close the output file for conserved regions
Conserved_region_txt_open.close()

# Find the longest conserved sequence and write it to the main conserved file
max_len = 0
max_description = ""
for record in SeqIO.parse(Conserved_region_txt, "fasta"):
    if len(record.seq) > max_len:
        max_len = len(record.seq)
        Sequence = record.seq
        max_ID = record.id

# Format the sequence and write it to the main conserved file
New_conserved = textwrap.fill("".join(Sequence), 70)
Main_conserved_open.write('>' + max_ID + '\n')
Main_conserved_open.write(New_conserved)

# Read the phylogenetic tree and plot it
tree_file = f"{filepath}/Tree.phy"
tree = Phylo.read(tree_file, "newick")

# Create and save the phylogenetic tree plot
fig_2 = plt.figure(figsize=(13, 5), dpi=100)
fig_2.suptitle('The Phylogenetic Tree', fontsize=20)
axes = fig_2.add_subplot(1, 1, 1)
Phylo.draw(tree, axes=axes, branch_labels=lambda c: c.branch_length, do_show=False)
fig_2.savefig(f"{filepath}/The_Phylogenetic_Tree", figsize=(13, 5))
