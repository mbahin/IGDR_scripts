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
    gene_name = line.split('\t')[8].split('"')[3]
    # Removing the transcript number if present
    match = re.match(r'(.*)_[0-9]*$',gene_name)
    if match:
        gene_name = match.group(1)
    enscafg = line.split('\t')[17].split('"')[1]
    if not RLOCs.has_key(rloc):
        RLOCs[rloc] = {}
        RLOCs[rloc]['enscafgs'] = []
        RLOCs[rloc]['gene_names'] = []
    if enscafg not in RLOCs[rloc]['enscafgs']:
        RLOCs[rloc]['enscafgs'].append(enscafg)
        if RLOCs[rloc]['gene_names'] == []:
            if (gene_name.startswith('ENSCAFG')) or (gene_name.startswith('CFRNASEQ')):
                RLOCs[rloc]['gene_names'].append('No_name')
            else:
                RLOCs[rloc]['gene_names'].append(gene_name)
        #elif (RLOCs[rloc]['gene_names'][0] != gene_name) and (RLOCs[rloc]['gene_names'][0] != 'No_name') and (not gene_name.startswith('CFRNASEQ')):
            #print enscafg, gene_name, RLOCs[rloc]['gene_names']
        elif gene_name not in RLOCs[rloc]['gene_names'] and not ((gene_name.startswith('ENSCAFG')) or (gene_name.startswith('CFRNASEQ'))):
            if RLOCs[rloc]['gene_names'][0] == 'No_name':
                RLOCs[rloc]['gene_names'] == [gene_name]
            else:
                RLOCs[rloc]['gene_names'].append(gene_name)

for i in sorted(RLOCs):
    output.write(i)
    #for j in range(0, len(RLOCs[i]['enscafgs'])):
    output.write('\t'+str(RLOCs[i]['enscafgs'])+'\t'+str(RLOCs[i]['gene_names'])+'\n')

input.close()
output.close()