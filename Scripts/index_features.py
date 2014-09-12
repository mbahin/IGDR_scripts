#!/local/python/2.7/bin/python

# Mathieu Bahin, 21/05/14

# Script to index the RLOCs and ENSCAFGs, 2 index files are created.
# The RLOC index list the ENSCAFGs corresponding to each RLOC (from the current GTF file). If none is found, information can sometimes be found about human orthologous gene names. If so, the second column is 'NA' and the additive information is in the third column. Most of the times, there is the RLOC, its ENSCAFGs list and 'NA' in the third column.
# Inputs are the current GTF file, a BioMart index (for the Ensembl gene names) and the 2 output index files to produce.
# The symbolic link on the two indexes are created in the directory.

import argparse, sys, re, datetime, os

today = str('{:02d}'.format(datetime.datetime.now().month))+'-'+str('{:02d}'.format(datetime.datetime.now().day))+'-'+str(datetime.datetime.now().year)

# Getting options back
parser = argparse.ArgumentParser()
parser.add_argument('-g', dest='GTF', default='/home/genouest/umr6061/recomgen/dog/data/canFam3/annotation/MasterAnnotation/canfam3_cons_annot.gtf')
parser.add_argument('-b', dest='biomart', default='/home/genouest/umr6061/recomgen/dog/data/canFam3/annotation/Correspondence_Indexes/BioMart_paralogous.05-23-2014.txt')
parser.add_argument('-r', dest='RLOC', default='/home/genouest/umr6061/recomgen/dog/data/canFam3/annotation/Correspondence_Indexes/RLOCs_index.'+today+'.txt')
parser.add_argument('-e', dest='ENSCAFG', default='/home/genouest/umr6061/recomgen/dog/data/canFam3/annotation/Correspondence_Indexes/ENSCAFGs_index.'+today+'.txt')
options = parser.parse_args()

# Indexing the GTF file
RLOCs = {}
ENSCAFGs = {}
with open(options.GTF,'r') as GTF_file:
    for line in GTF_file:
        # Only working on the exon lines
        if line.split('\t')[2] != 'exon':
            continue

        # Getting the RLOC, ENSCAFG, gene name and biotype
        enscafg = 'NA'
        name = 'NA'
        rloc = line.split('\t')[8].split('"')[1]
        biotype = line.split('\t')[8].split('"')[9]
        
        # Getting the transcript ID (an ENSCAFG or a gene name)
        if not line.split('\t')[8].split('"')[3].startswith(('CFRNASEQ','ENSCAFT')):
            tcpt_id = re.match(r'(.*)_([2-9]|[1-9][0-9]+)$',line.split('\t')[8].split('"')[3])
            if tcpt_id:
                tcpt = tcpt_id.group(1)
            else:
                tcpt = line.split('\t')[8].split('"')[3]
        else:
            tcpt = 'NA'
        if line.split('\t')[8].split('"')[5] != 'NA':
            enscafg = line.split('\t')[8].split('"')[5].split(',')[0] # When there a several ENSCAFGs, they are identicals
            if not tcpt.startswith('ENSCAFG'):
                name = tcpt
        elif tcpt.startswith('ENSCAFG'):
            enscafg = tcpt
        else:
            name = tcpt

        # Creating the RLOC in the index if necessary
        if not RLOCs.has_key(rloc):
            RLOCs[rloc] = {}
            RLOCs[rloc]['enscafgs'] = ['NA']
            RLOCs[rloc]['orthologues'] = 'NA'
            RLOCs[rloc]['biotype'] = []

        # Creating the ENSCAFG in the index if necessary
        if (enscafg != 'NA') and (not ENSCAFGs.has_key(enscafg)):
            ENSCAFGs[enscafg] = {}
            ENSCAFGs[enscafg]['ENSEMBL_name'] = 'NA'
            ENSCAFGs[enscafg]['BROAD_name'] = 'NA'
            ENSCAFGs[enscafg]['consensus_name'] = []

        # Filling the indexes
        if enscafg != 'NA':
            if not ENSCAFGs[enscafg].has_key('RLOC'):
                ENSCAFGs[enscafg]['RLOC'] = rloc
            if enscafg not in RLOCs[rloc]['enscafgs']:
                if 'NA' in RLOCs[rloc]['enscafgs']:
                    RLOCs[rloc]['enscafgs'].remove('NA')
                RLOCs[rloc]['enscafgs'].append(enscafg)
            if name != 'NA':
                RLOCs[rloc]['orthologues'] = 'NA' # If there is an ENSCAFG and one of the transcripts has a gene name, then no name is needed in the RLOC index
                if (name != ENSCAFGs[enscafg]['BROAD_name']) and (ENSCAFGs[enscafg]['BROAD_name'] != 'NA'):
                    print 'Warning: the BROAD provided 2 names for '+enscafg+'!'
                elif ENSCAFGs[enscafg]['BROAD_name'] == 'NA':
                    ENSCAFGs[enscafg]['BROAD_name'] = name
        elif (RLOCs[rloc]['enscafgs'] == ['NA']) and (name != 'NA'):
            if (RLOCs[rloc]['orthologues'] != 'NA') and (RLOCs[rloc]['orthologues'] != name):
                print 'Warning: there are two orthologues names for '+rloc+'!'
            RLOCs[rloc]['orthologues'] = name
        if biotype not in RLOCs[rloc]['biotype']:
            RLOCs[rloc]['biotype'].append(biotype)

# Indexing the Biomart file
with open(options.biomart,'r') as biomart_file:
    biomart_file.readline()
    for line in biomart_file:
        enscafg = line.split('\t')[0]
        name = line.split('\t')[1]
        if name:
            if not ENSCAFGs.has_key(enscafg):
                print 'Warning: '+enscafg+' was not listed in the GFF file!'
            elif (ENSCAFGs[enscafg]['ENSEMBL_name'] != 'NA') and (ENSCAFGs[enscafg]['ENSEMBL_name'] != name):
                print 'Warning: Biomart provided 2 different names for '+enscafg+'!'
            else:
                ENSCAFGs[enscafg]['ENSEMBL_name'] = name

# Creating a consensus name column in the ENSCAFGs index
for enscafg in ENSCAFGs:
    if (ENSCAFGs[enscafg]['ENSEMBL_name'] == ENSCAFGs[enscafg]['BROAD_name']) or (re.match(r''+ENSCAFGs[enscafg]['ENSEMBL_name']+'_CANFA',ENSCAFGs[enscafg]['BROAD_name'])) or (ENSCAFGs[enscafg]['BROAD_name'] == 'NA'):
        # If the BROAD name is different from Ensembl name only adding '_CANFA' by the end, the consensus is Ensembl name
        ENSCAFGs[enscafg]['consensus_name'] = [ENSCAFGs[enscafg]['ENSEMBL_name']]
    elif ENSCAFGs[enscafg]['ENSEMBL_name'] == 'NA':
        ENSCAFGs[enscafg]['consensus_name'] = [ENSCAFGs[enscafg]['BROAD_name']]
    else:
        ENSCAFGs[enscafg]['consensus_name'].append(ENSCAFGs[enscafg]['ENSEMBL_name'])
        ENSCAFGs[enscafg]['consensus_name'].append(ENSCAFGs[enscafg]['BROAD_name'])

# Writing the output index files

with open(options.RLOC,'w') as RLOC_output:
    RLOC_output.write('RLOC\tENSCAFG(s)\tOrthologous_name\tBiotype(s)\n')
    for rloc in sorted(RLOCs):
        RLOC_output.write(rloc+'\t'+','.join(sorted(RLOCs[rloc]['enscafgs']))+'\t'+RLOCs[rloc]['orthologues']+'\t'+','.join(sorted(RLOCs[rloc]['biotype']))+'\n')

with open(options.ENSCAFG,'w') as ENSCAFG_output:
    ENSCAFG_output.write('ENSCAFG\tEnsembl_name\tBROAD_name\tConsensus_name(s)\tRLOC\n')
    for enscafg in sorted(ENSCAFGs):
        ENSCAFG_output.write(enscafg+'\t'+ENSCAFGs[enscafg]['ENSEMBL_name']+'\t'+ENSCAFGs[enscafg]['BROAD_name']+'\t'+'|'.join(sorted(ENSCAFGs[enscafg]['consensus_name']))+'\t'+ENSCAFGs[enscafg]['RLOC']+'\n')

# Creating the symbolic link pointing at the current indexes versions
os.remove('/home/genouest/umr6061/recomgen/dog/data/canFam3/annotation/Correspondence_Indexes/ENSCAFGs_index.txt')
os.symlink(options.ENSCAFG,'/home/genouest/umr6061/recomgen/dog/data/canFam3/annotation/Correspondence_Indexes/ENSCAFGs_index.txt')
os.remove('/home/genouest/umr6061/recomgen/dog/data/canFam3/annotation/Correspondence_Indexes/RLOCs_index.txt')
os.symlink(options.RLOC,'/home/genouest/umr6061/recomgen/dog/data/canFam3/annotation/Correspondence_Indexes/RLOCs_index.txt')