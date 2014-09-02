#!/local/python/2.7/bin/python

# Mathieu Bahin, 21/05/14

# Script to check if getting the gene name from the transcript ID removing the part after the last underscrore (when present, i.e. when it's not the first transcript) is correct.
# We remove the part and then test if this gene name is present more than once.
# The input is a GFF file.

import sys,re

input = open(sys.argv[1],'r')

# To check the gene names in the GTF file
for line in input:
    tscpt_name = line.split('\t')[8].split('"')[3]
    match = re.match(r'(.*)_[0-9]*_([2-9]|[1-9][0-9]+)$',tscpt_name)
    if match:
        print 'Warning:', match.group(1)

# To check the gene names in the ENSCAFGs index
"""
for line in input:
    gene_name = line.split('\t')[1]
    match = re.match(r'(.*)_([2-9]|[1-9][0-9]+)$',gene_name)
    if match:
        print 'Warning:', match.group(1)
"""

input.close()
