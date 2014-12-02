#!/local/python/2.7/bin/python

# Mathieu Bahin, 21/05/14

# Script to process chimCT output. The chimeras are split into 4 categories according to the number of part matching a known feature: 2 different, twice the same, only one, zero. Information are provided about the chimeras and the features matched.
# The input is the sample name and the CRAC output directory.
# Outputs are the 4 file describing the chimeras a link to the spanning reads fasta file and another link to the CRAC output directory.

import argparse, sys, shutil, re, os, glob, stat
from collections import OrderedDict

# Getting the functions from 'classics.py' and the RLOCs and ENSCAFGs indexes
sys.path.insert(1, '/home/genouest/genouest/mbahin/Scripts')
import classics, get_indexes
RLOCs = get_indexes.RLOCs
ENSCAFGs = get_indexes.ENSCAFGs
sys.path.remove('/home/genouest/genouest/mbahin/Scripts')

# Getting options back
parser = argparse.ArgumentParser()
parser.add_argument('-n', dest='sample_name', required=True)
parser.add_argument('-c', dest='crac_dir', required=True)
parser.add_argument('-s', dest='single_end', action='store_true')
options = parser.parse_args()
options.crac_dir = options.crac_dir.rstrip('/')

# Creating the directory to store result files
dir_name = options.sample_name+'_processed.dir'
classics.create_wd(dir_name)

###### Functions

def write_RLOC_info(rloc,file,part):
    #####
    # Function to write the info on a RLOC (ENSCAFGs, gene name from BioMart and the BROAD and gene biotype) into the result files.
    #####
    if RLOCs.has_key(rloc):
        if RLOCs[rloc]['enscafgs'] != ['NA']:
            file.write(','.join(RLOCs[rloc]['enscafgs']))
            names = []
            in_mutation_list = []
            in_translocation_list = []
            for enscafg in RLOCs[rloc]['enscafgs']:
                names = names + ENSCAFGs[enscafg]['consensus_names']
                if enscafg in mutations:
                    in_mutation_list.append('Yes')
                else:
                    in_mutation_list.append('No')
                if enscafg in translocations:
                    in_translocation_list.append('Yes')
                else:
                    in_translocation_list.append('No')
            biotypes = RLOCs[rloc]['biotypes']
            file.write('\t'+','.join(names)+'\t'+','.join(biotypes)+'\t'+','.join(in_mutation_list)+'\t'+','.join(in_translocation_list)+'\t')
        else:
            file.write('No_ENSCAFG\t'+RLOCs[rloc]['orthologous']+'\tNA\tNA\tNA\t')
    else:
        file.write('No_ENSCAFG\tNA\tNA\tNA\tNA\t')

def write_res(file,one_feat):
    #####
    # Function to write the information relative to a line in the chimera file from chimCT into the result files.
    #####

    if options.single_end == False:
        file.write(chim+'\t'+'\t'.join((str(orderedChim[chim]['spR']),str(orderedChim[chim]['spP']),str(orderedChim[chim]['class'])))+'\t')
    else:
        file.write(chim+'\t'+'\t'.join((str(orderedChim[chim]['spR']),'NA',str(orderedChim[chim]['class'])))+'\t')
    breakpoints = []
    for pos in orderedChim[chim]['breakpoint']:
        breakpoints.append('('+str(pos[0])+':'+str(pos[1])+','+str(pos[2])+' / '+str(pos[3])+':'+str(pos[4])+','+str(pos[5])+' # '+str(pos[6])+','+str(pos[7])+')')
    file.write(','.join(breakpoints)+'\t'+str(','.join(orderedChim[chim]['warnings']))+'\t')
    if one_feat == 'First':
        write_RLOC_info(feats[0],file,1)
        file.write('\n')
    elif one_feat == 'Second':
        write_RLOC_info(feats[1],file,2)
        file.write('\n')
    else:
        i = 1
        for rloc in feats:
            write_RLOC_info(rloc,file,i)
            i+= 1
##### Functions end

# Indexing the paralogues file (from BioMart)
with open('/home/genouest/umr6061/recomgen/dog/data/canFam3/annotation/Correspondence_Indexes/BioMart_paralogous.txt','r') as paral_file:
    paral_file.readline()
    for line in paral_file:
        enscafg = line.split('\t')[0]
        if not ENSCAFGs[enscafg].has_key('paralogous_genes'):
            ENSCAFGs[enscafg]['paralogous_genes'] = []
            ENSCAFGs[enscafg]['paralogy_type'] = []
        if (line.split('\t')[2]) and ((not ENSCAFGs[enscafg].has_key('paralogous_genes')) or (line.split('\t')[2] not in ENSCAFGs[enscafg]['paralogous_genes'])):
            paralog = line.split('\t')[2]
            paralogy_type = line.rstrip().split('\t')[3]
            ENSCAFGs[enscafg]['paralogous_genes'].append(paralog)
            ENSCAFGs[enscafg]['paralogy_type'].append(paralogy_type)

# Indexing the cancer mutation list
mutations = []
with open('/home/genouest/genouest/mbahin/Annotations/mutation_gene_list.txt','r') as mutation_file:
    for line in mutation_file:
        mutations.append(line.split('\t')[0])

# Indexing the cancer translocation list
translocations = []
with open('/home/genouest/genouest/mbahin/Annotations/translocation_gene_list.txt','r') as translocation_file:
    for line in translocation_file:
        translocations.append(line.split('\t')[0])

# Processing chimCT output
chim = {}
noFeat = []
with open(glob.glob(os.getcwd()+'/../*.chimCT.txt')[0],'r') as input:
    for line in input:
        if line.startswith('#'):
            continue
        # Catching line information
        (chimID,chr1,end1,strand1,chr2,start2,strand2) = (line.split('\t')[1:8])
        (chimClass,warnings) = (line.split('\t')[11:13])
        spR = re.match(r'Nb_spanning_reads=(.*)',line.split('\t')[21]).group(1)
        if options.single_end == False:
            spP = re.match(r'Nb_spanning_PE=(.*)',line.rstrip().split('\t')[22]).group(1)
        else:
            spP = 'NA'
        if chimClass == '4':
            continue
        match = re.match(r'FusionDistance=([^,]*)',warnings)
        if match:
            warn = match.group(1)
        else:
            warn = 'No_warning'
        # Chimeras with 2 unmatched features
        if chimID == 'N/A---N/A':
            noFeat.append((chimClass,spR,spP,warn,chr1,end1,strand1,chr2,start2,strand2))
            continue
        # Commands executed only the first time the chimera is met
        if not chim.has_key(chimID):
            chim[chimID] = {}
            chim[chimID]['spR'] = 0
            if options.single_end == False:
                chim[chimID]['spP'] = 0
            else:
                spP = 'NA'
            chim[chimID]['class'] = chimClass
            chim[chimID]['breakpoint'] = []
            chim[chimID]['warnings'] = []
        # Commands executed for each processed line
        chim[chimID]['spR'] += int(spR)
        if options.single_end == False:
            chim[chimID]['spP'] += int(spP)
        chim[chimID]['warnings'].append(warn)
        chim[chimID]['breakpoint'].append((chr1,end1,strand1,chr2,start2,strand2,spR,spP))

# Opening output files (permissi
noFeat_file = open(options.sample_name+'.noFeat.xls','w')
os.chmod(noFeat_file.name, stat.S_IRWXU | stat.S_IRGRP | stat.S_IROTH)
oneFeat_file = open(options.sample_name+'.oneFeat.xls','w')
os.chmod(oneFeat_file.name, stat.S_IRWXU | stat.S_IRGRP | stat.S_IROTH)
monoFeat_file = open(options.sample_name+'.monoFeat.xls','w')
os.chmod(monoFeat_file.name, stat.S_IRWXU | stat.S_IRGRP | stat.S_IROTH)
twoFeat_file = open(options.sample_name+'.twoFeat.xls','w')
os.chmod(twoFeat_file.name, stat.S_IRWXU | stat.S_IRGRP | stat.S_IROTH)

# Writing the file for chimera with 2 unmatched parts
noFeat_file.write('\t'.join(('nb spanning reads','nb spanning PE','class','chr1:pos1,strand1','chr2:pos2,strand2','warning'))+'\n')
noFeat.sort(key = lambda x:int(x[1
]),reverse=True)
for c in noFeat:
    noFeat_file.write(c[1]+'\t'+c[2]+'\t'+c[0]+'\t'+c[4]+':'+c[5]+','+c[6]+'\t'+c[7]+':'+c[8]+','+c[9]+'\t'+c[3]+'\n')

# Writing the files for chimeras with at least one matched part
oneFeat_file.write('\t'.join(('chimID','nb spanning reads','nb spanning PE','class','(chr1:pos1,strand1 / chr2:pos2,strand2)','warning','ENSCAFG_list','gene_name_list','biotype','in mutation list','in translocation list'))+'\n')
monoFeat_file.write('\t'.join(('chimID','nb spanning reads','nb spanning PE','class','(chr1:pos1,strand1 / chr2:pos2,strand2)','warning','ENSCAFG_list','gene_name_list','biotype','in mutation list','in translocation list'))+'\n')
twoFeat_file.write('\t'.join(('chimID','nb spanning reads','nb spanning PE','class','(chr1:pos1,strand1 / chr2:pos2,strand2)','warning','ENSCAFG1_list','gene_name1_list','biotype1','in mutation list','in translocation list','ENSCAFG2_list','gene_name2_list','biotype2','in mutation list','in translocation list','paralogy'))+'\n')
orderedChim = OrderedDict(sorted(chim.iteritems(), key = lambda x:x[1]['spR'], reverse=True))
#for i in orderedChim:
    #print i,orderedChim[i]

for chim in orderedChim:
    #RLOCs = chim.split('---')
    feats = chim.split('---')
    # Writing the file for chimera with 2 matched features
    if ('N/A' not in feats) and (feats[0] != feats[1]):
        write_res(twoFeat_file, 'Both')
        # Determining if the two features are considered paralogous
        paralogous = False
        if RLOCs.has_key(feats[0]) and RLOCs.has_key(feats[1]):
            for enscafg1 in RLOCs[feats[0]]['enscafgs']:
                for enscafg2 in RLOCs[feats[1]]['enscafgs']:
                    if ENSCAFGs.has_key(enscafg1):
                        if enscafg2 in ENSCAFGs[enscafg1]['paralogous_genes']:
                            paralogous = True
                            break
                if paralogous:
                    break
        if paralogous:
            twoFeat_file.write('Paralogous\n')
        else:
            twoFeat_file.write('No_paralogy\n')
    # Writing the file for chimera with the same feature matched twice
    elif feats[0] == feats[1]:
        write_res(monoFeat_file, 'First')
    # Writing the files for chimera with only one matched feature
    elif feats[1] == 'N/A':
        write_res(oneFeat_file, 'First')
    else:
        write_res(oneFeat_file, 'Second')

# Closing output files
noFeat_file.close()
oneFeat_file.close()
monoFeat_file.close()
twoFeat_file.close()

# Creating symbolic links to the spanning split reads fasta file and to the CRAC output directory of the analyse
os.symlink(os.path.abspath(os.getcwd()+'/../'+options.sample_name+'.fa'),options.sample_name+'.fa')
os.symlink(options.crac_dir,'link_to_BAM_files_directory.ln')