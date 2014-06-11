#!/local/python/2.7/bin/python

# Mathieu Bahin, 21/05/14

# Script to index the RLOCs into list of ENSCAFGs and gene names.
# Input is a file intersecting BROAD and ensembl annotation.
# Output is a tab-separated file with the RLOC and its ENSCAFG(s).

import sys, re

# Opening files
intersection_file = open(sys.argv[1],'r')
gff_file = open(sys.argv[2],'r')
RLOC_output = open(sys.argv[3],'w')
ENSCAFG_output = open(sys.argv[4], 'w')

# Parsing the intersection file
RLOCs = {}
ENSCAFGs = {}
for line in intersection_file:
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
        RLOCs[rloc] = {}
        RLOCs[rloc]['enscafg'] = []
        RLOCs[rloc]['BROAD_bonus'] = ''
    if enscafg not in RLOCs[rloc]['enscafg']:
        #RLOCs[rloc].append(enscafg)
        RLOCs[rloc]['enscafg'].append(enscafg)
    if (name.startswith('ENSCAFG')) and (name != enscafg) and (name not in RLOCs[rloc]['enscafg']):
        #RLOCs[rloc].append(name)
        RLOCs[rloc]['enscafg'].append(name)
    # Filling ENSCFGs index
        # Finding the cases where the BROAD and Ensmbl doesn't agree on the ENSCAFG
    #if name.startswith ('ENSCAFG') and name != enscafg:
        #print 'Oh oh... :',enscafg, name
        #exit
    # Filling gene name information
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
    # Filling the biotype information
    if ENSCAFGs[enscafg]['biotype'] != 'Ambiguous':
        if biotype_BROAD != biotype_Ensembl:
            ENSCAFGs[enscafg]['biotype'] = 'Ambiguous'
        elif ENSCAFGs[enscafg]['biotype'] == '':
            ENSCAFGs[enscafg]['biotype'] = biotype_BROAD
        elif ENSCAFGs[enscafg]['biotype'] != biotype_BROAD:
            ENSCAFGs[enscafg]['biotype'] = 'Multiple'

# Catching more information from the GFF file (BROAD names without ENSCAFG)
###
# The intersection file only lists features for which BROAD and Ensembl annotations match at a certain percentage. The GFF file might have more information on several RLOCs (an ENSCAFG or a gene name or a human paralogous gene name). These information are added to the RLOC dictionary.
###
for line in gff_file:
    if line.split('\t')[2] == 'mRNA':
        enscafg = ''
        match = re.match(r'ID=(.*);Parent=(RLOC_.*);exons',line.split('\t')[8])
        match2 = re.match(r'(.*)_[0-9]*$',match.group(1))
        if match2:
            enscafg = match2.group(1)
        else:
            enscafg = match.group(1)
        # Checking if there is an ENSCAFG from the GFF file not found in the intersection file
        #if enscafg.startswith('ENSCAFG') and (RLOCs.has_key(match.group(2))) and (enscafg not in RLOCs[match.group(2)]['enscafg']):
            #print "Oh no !",match.group(2),enscafg,str(RLOCs[match.group(2)]['enscafg'])
        if (not RLOCs.has_key(match.group(2))) and (not match.group(1).startswith('CFRNASEQ')):
            RLOCs[match.group(2)] = {}
            if enscafg.startswith('ENSCAFG'):
                RLOCs[match.group(2)]['enscafg'] = enscafg
                RLOCs[match.group(2)]['BROAD_bonus'] = ''
            else:
                RLOCs[match.group(2)]['enscafg'] = ['No_ENSCAFG']
                RLOCs[match.group(2)]['BROAD_bonus'] = enscafg

# Writing index files
for rloc in sorted(RLOCs):
    RLOC_output.write(rloc+'\t'+','.join(RLOCs[rloc]['enscafg'])+'\t'+RLOCs[rloc]['BROAD_bonus']+'\n')

for enscafg in sorted(ENSCAFGs):
    ENSCAFG_output.write(enscafg+'\t'+ENSCAFGs[enscafg]['name']+'\t'+ENSCAFGs[enscafg]['biotype']+'\n')

# Closing fils
intersection_file.close()
gff_file.close()
RLOC_output.close()
ENSCAFG_output.close()