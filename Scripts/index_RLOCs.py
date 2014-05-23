#!/local/python/2.7/bin/python

# Mathieu Bahin, 21/05/14

# Script to index the RLOCs into list of ENSCAFGs and gene names.
# Input is a file intersecting BROAD and ensembl annotation.
# Output is a tab-separated file with the RLOC and its ENSCAFG(s).

import sys, re

input = open(sys.argv[1],'r')
output = open(sys.argv[2],'w')

RLOCs = {}
for line in input:
    rloc = line.split('\t')[8].split('"')[1]
    enscafg = line.split('\t')[17].split('"')[1]
    if not RLOCs.has_key(rloc):
        RLOCs[rloc] = []
    if enscafg not in RLOCs[rloc]:
        RLOCs[rloc].append(enscafg)

for i in sorted(RLOCs):
    output.write(i+'\t'+'\t'.join(RLOCs[i])+'\n')

input.close()
output.close()