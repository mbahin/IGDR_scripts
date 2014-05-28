#!/local/python/2.7/bin/python

# Mathieu Bahin, 21/05/14

# Script to index the RLOCs into list of ENSCAFGs and gene names.
# Input is a file intersecting BROAD and ensembl annotation.
# Output is a tab-separated file with the RLOC and its ENSCAFG(s).

import sys, re

input = open(sys.argv[1],'r')
RLOC_output = open(sys.argv[2],'w')
ENSCAFG_output = open(sys.argv[3], 'w')

RLOCs = {}
ENSCAFGs = {}
for line in input:
    rloc = line.split('\t')[8].split('"')[1]
    name = line.split('\t')[8].split('"')[3]
    # Removing the transcript number if present
    match = re.match(r'(.*)_[0-9]*$',name)
    if match:
        name = match.group(1)
    enscafg = line.split('\t')[17].split('"')[1]
    if not RLOCs.has_key(rloc):
        RLOCs[rloc] = []
    if enscafg not in RLOCs[rloc]:
        RLOCs[rloc].append(enscafg)
    if not name.startswith('ENSCAFG') and not name.startswith('CFRNASEQ'):
        if ENSCAFGs.has_key(enscafg) and ENSCAFGs[enscafg] != name:
            sys.exit('Two names for ENSCAFG '+enscafg+'! Aborting.')
        ENSCAFGs[enscafg] = name
    else:
        if ENSCAFGs.has_key(enscafg) and ENSCAFGs[enscafg] != 'No_name':
            sys.exit('Two names for ENSCAFG '+enscafg+'! Aborting.')
        ENSCAFGs[enscafg] = 'No_name'

# Writing index files
for rloc in sorted(RLOCs):
    RLOC_output.write(rloc+'\t'+','.join(RLOCs[rloc])+'\n')

for i in sorted(ENSCAFGs):
    ENSCAFG_output.write(i+'\t'+ENSCAFGs[i]+'\n')

input.close()
RLOC_output.close()
ENSCAFG_output.close()