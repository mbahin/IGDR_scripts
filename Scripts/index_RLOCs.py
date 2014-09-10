#!/local/python/2.7/bin/python

# Mathieu Bahin, 21/05/14

# Script to index the RLOCs and ENSCAFGs, 2 index files are created.
# The RLOC index list the ENSCAFGs corresponding to each RLOC (from the file intersecting BROAD and Ensembl features). If none is found, information can sometimes be found in the GFF file (particularly human homologous gene names. If so, the second column is 'No_ENSCAFG' and the additive information is in the third column. Else, if no information can be found, the RLOC is not present in the index). Most of the times, there is the RLOC, its ENSCAFGs list and 'No_BROAD_bonus' in the third column (if an ENSCAFG is found in the intersection file, no information is searched in the GFF file).
# Input is a file intersecting BROAD and ensembl annotation, the GFF file corresponding and the 2 index filenames to produce.
# Output are two tab-separated files with the RLOCs and ENSCAFGs indexes.

import argparse, sys, re, datetime

today = str('{:02d}'.format(datetime.datetime.now().month))+'-'+str('{:02d}'.format(datetime.datetime.now().day))+'-'+str(datetime.datetime.now().year)

# Getting options back
parser = argparse.ArgumentParser()
parser.add_argument('-g', dest='GTF', default='/home/genouest/umr6061/recomgen/dog/data/canFam3/annotation/MasterAnnotation/BROADmRNA_lncRNA_antis.Ens75.gtfclean.09-02-2014.gtf')
parser.add_argument('-b', dest='biomart', default='/home/genouest/genouest/mbahin/Annotations/BioMart_paralogous_140523.txt')
parser.add_argument('-r', dest='RLOC', default='/home/genouest/genouest/mbahin/Annotations/RLOCs_index.'+today+'.txt')
parser.add_argument('-e', dest='ENSCAFG', default='/home/genouest/genouest/mbahin/Annotations/ENSCAFGs_index.'+today+'.txt')
options = parser.parse_args()

# Indexing the GTF file
RLOCs = {}
ENSCAFGs = {}
with open(options.GTF,'r') as GTF_file:
    for line in GTF_file:
        # Getting the RLOC, ENSCAFG, gene name and biotype
        enscafg = 'No_ENSCAFG'
        name = 'No_name'
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
            tcpt = 'No_name'
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
            RLOCs[rloc]['enscafgs'] = ['No_ENSCAFG']
            RLOCs[rloc]['orthologues'] = 'NA'
            RLOCs[rloc]['biotype'] = []

        # Creating the ENSCAFG in the index if necessary
        if enscafg and (not ENSCAFGs.has_key(enscafg)):
            ENSCAFGs[enscafg] = {}
            ENSCAFGs[enscafg]['ENSEMBL_name'] = 'No_name'
            ENSCAFGs[enscafg]['BROAD_name'] = 'No_name'

        # Filling the indexes
        if enscafg != 'No_ENSCAFG':
            if enscafg not in RLOCs[rloc]['enscafgs']:
                if 'No_ENSCAFG' in RLOCs[rloc]['enscafgs']:
                    RLOCs[rloc]['enscafgs'].remove('No_ENSCAFG')
                RLOCs[rloc]['enscafgs'].append(enscafg)
            if name != 'No_name':
                RLOCs[rloc]['orthologues'] = 'NA' # If there is an ENSCAFG and one of the transcripts has a gene name, then no name is needed in the RLOC index
                if (name != ENSCAFGs[enscafg]['BROAD_name']) and (ENSCAFGs[enscafg]['BROAD_name'] != 'No_name'):
                    print 'Warning: the BROAD provided 2 names for '+enscafg+'!'
                elif ENSCAFGs[enscafg]['BROAD_name'] == 'No_name':
                    ENSCAFGs[enscafg]['BROAD_name'] = name
        elif (RLOCs[rloc]['enscafgs'] == ['No_ENSCAFG']) and (name != 'No_name'):
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
            elif (ENSCAFGs[enscafg]['ENSEMBL_name'] != 'No_name') and (ENSCAFGs[enscafg]['ENSEMBL_name'] != name):
                print 'Warning: Biomart provided 2 different names for '+enscafg+'!'
            else:
                ENSCAFGs[enscafg]['ENSEMBL_name'] = name

# Creating a consensus name column in the ENSCAFGs index
for enscafg in ENSCAFGs:
    ENSCAFGs[enscafg]['consensus'] = []
    if (ENSCAFGs[enscafg]['ENSEMBL_name'] == ENSCAFGs[enscafg]['BROAD_name']) or (re.match(r''+ENSCAFGs[enscafg]['ENSEMBL_name']+'_CANFA',ENSCAFGs[enscafg]['BROAD_name'])) or (ENSCAFGs[enscafg]['BROAD_name'] == 'No_name'):
        # If the BROAD name is different from Ensembl name only adding '_CANFA' by the end, the consensus is Ensembl name
        ENSCAFGs[enscafg]['consensus'] = [ENSCAFGs[enscafg]['ENSEMBL_name']]
    elif ENSCAFGs[enscafg]['ENSEMBL_name'] == 'No_name':
        ENSCAFGs[enscafg]['consensus'] = [ENSCAFGs[enscafg]['BROAD_name']]
    else:
        ENSCAFGs[enscafg]['consensus'].append(ENSCAFGs[enscafg]['ENSEMBL_name'])
        ENSCAFGs[enscafg]['consensus'].append(ENSCAFGs[enscafg]['BROAD_name'])

# Writing the output index files

with open(options.RLOC,'w') as RLOC_output:
    RLOC_output.write('RLOC\tENSCAFG(s)\tOrthologous_name\tBiotype(s)\n')
    for rloc in sorted(RLOCs):
        RLOC_output.write(rloc+'\t'+','.join(sorted(RLOCs[rloc]['enscafgs']))+'\t'+RLOCs[rloc]['orthologues']+'\t'+','.join(sorted(RLOCs[rloc]['biotype']))+'\n')

with open(options.ENSCAFG,'w') as ENSCAFG_output:
    ENSCAFG_output.write('ENSCAFG\tEnsembl_name\tBROAD_name\tConsensus_name(s)\n')
    for enscafg in sorted(ENSCAFGs):
        ENSCAFG_output.write(enscafg+'\t'+ENSCAFGs[enscafg]['ENSEMBL_name']+'\t'+ENSCAFGs[enscafg]['BROAD_name']+'\t'+'|'.join(sorted(ENSCAFGs[enscafg]['consensus']))+'\n')

# Particular case: grep RLOC_00021852 $GTFv3_2
# Giving an "aberration" in the RLOC index: RLOC_00021852	ENSCAFG00000031893,ENSCAFG00000032619	SOX2