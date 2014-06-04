#!/local/python/2.7/bin/python

# Mathieu Bahin, 21/05/14

# Script to index the RLOCs into list of ENSCAFGs and gene names.
# Input is a file intersecting BROAD and ensembl annotation.
# Output is a tab-separated file with the RLOC and its ENSCAFG(s).

import sys, re

# Opening files
input = open(sys.argv[1],'r')
RLOC_output = open(sys.argv[2],'w')
ENSCAFG_output = open(sys.argv[3], 'w')

# Parsing the input file
RLOCs = {}
ENSCAFGs = {}
for line in input:
    # Getting the 'RLOC',the BROAD name or 'ENSCAFG' if there is no name and the ENSCAFG from Ensembl
    rloc = line.split('\t')[8].split('"')[1]
    name = line.split('\t')[8].split('"')[3]
    # Removing the transcript number if present
    match = re.match(r'(.*)_[0-9]*$',name)
    if match:
        name = match.group(1)
    biotype_BROAD = line.split('\t')[8].split('"')[5]
    biotype_Ensembl = line.split('\t')[10]
    enscafg = line.split('\t')[17].split('"')[1]
    # Filling RLOCs index
    if not RLOCs.has_key(rloc):
        RLOCs[rloc] = []
    if enscafg not in RLOCs[rloc]:
        RLOCs[rloc].append(enscafg)
    if (name.startswith('ENSCAFG')) and (name != enscafg) and (name not in RLOCs[rloc]):
        RLOCs[rloc].append(name)
    # Filling ENSCFGs index
        # Finding the cases where the BROAD and Ensmbl doesn't agree on the ENSCAFG
    #if name.startswith ('ENSCAFG') and name != enscafg:
        #print 'Oh oh... :',enscafg, name
        #exit
    if not name.startswith('ENSCAFG') and not name.startswith('CFRNASEQ') and not name.startswith('ENSCAFT'):
        if ENSCAFGs.has_key(enscafg) and ENSCAFGs[enscafg]['name'] != name:
            sys.exit('Two names for ENSCAFG '+enscafg+'! Aborting.')
        if not ENSCAFGs.has_key(enscafg):
            ENSCAFGs[enscafg] = {}
            ENSCAFGs[enscafg]['biotype'] = ''
        ENSCAFGs[enscafg]['name'] = name
    else:
        if ENSCAFGs.has_key(enscafg) and ENSCAFGs[enscafg]['name'] != 'No_name':
            sys.exit('Two names for ENSCAFG '+enscafg+'! Aborting.')
        if not ENSCAFGs.has_key(enscafg):
            ENSCAFGs[enscafg] = {}
            ENSCAFGs[enscafg]['biotype'] = ''
        ENSCAFGs[enscafg]['name'] = 'No_name'
    if ENSCAFGs[enscafg]['biotype'] != 'Ambiguous':
        if biotype_BROAD != biotype_Ensembl:
            ENSCAFGs[enscafg]['biotype'] = 'Ambiguous'
        elif ENSCAFGs[enscafg]['biotype'] == '':
            ENSCAFGs[enscafg]['biotype'] = biotype_BROAD
        elif ENSCAFGs[enscafg]['biotype'] != biotype_BROAD:
            ENSCAFGs[enscafg]['biotype'] = 'Multiple'
        

# Writing index files
for rloc in sorted(RLOCs):
    RLOC_output.write(rloc+'\t'+','.join(RLOCs[rloc])+'\n')

for enscafg in sorted(ENSCAFGs):
    ENSCAFG_output.write(enscafg+'\t'+ENSCAFGs[enscafg]['name']+'\t'+ENSCAFGs[enscafg]['biotype']+'\n')

# Closing fils
input.close()
RLOC_output.close()
ENSCAFG_output.close()