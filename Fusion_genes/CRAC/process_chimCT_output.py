#!/local/python/2.7/bin/python

# Mathieu Bahin, 21/05/14

# Script to process chimCT output. The chimeras are split into 4 categories according to the number of part matching a known feature: 2 different, twice the same, only one, zero. Information are provided about the chimeras and the features matched.
# The input is the chimera file produced by chimCT.
# Outputs are the 4 file describing the chimeras.

import sys, re
from collections import OrderedDict

###### Functions

def write_RLOC_info(rloc,file):
    #####
    # Function to write the information on a RLOC (ENSCAFGs, gene name from BioMart and the BROAD) into the result files.
    #####
    if RLOCs_index.has_key(rloc):
        file.write('\t'+','.join(RLOCs_index[rloc]))
        names = ''
        for enscafg in RLOCs_index[rloc]:
            BM_name = ENSCAFGs_BM[enscafg]['gene_name']
            TD_name = ENSCAFGs_TD[enscafg]
            if BM_name == TD_name:
                names += BM_name+','
            else:
                names += BM_name+'||'+TD_name+','
        file.write('\t'+names.rstrip(',')+'\t')
    else:
        file.write('\tNo_feat\tNA\t')

def write_res(file,one_feat):
    #####
    # Function to write the information relative to a line in the chimera file from chimCT into the result files.
    #####
    file.write(chim+'\t'+'\t'.join((str(orderedChim[chim]['spR']),str(orderedChim[chim]['spP']),str(','.join(orderedChim[chim]['warnings']))))+'\t')
    for chimClass in orderedChim[chim]['breakpoint']:
        for pos in orderedChim[chim]['breakpoint'][chimClass]:
            file.write(' ('+str(pos[0])+':'+str(pos[1])+','+str(pos[2])+' / '+str(pos[3])+':'+str(pos[4])+','+str(pos[5])+' # '+str(pos[6])+','+str(pos[7])+')')
    if one_feat == 'First':
        write_RLOC_info(RLOCs[0],file)
        file.write('\n')
    elif one_feat == 'Second':
        write_RLOC_info(RLOCs[1],file)
        file.write('\n')
    else:
        for rloc in RLOCs:
            write_RLOC_info(rloc,file)
##### Functions end

# Indexing the RLOCs index
RLOCs_file = open('/home/genouest/genouest/mbahin/Annotations/RLOCs_index.txt','r')
RLOCs_index = {}
for line in RLOCs_file:
    RLOC = line.split('\t')[0]
    RLOCs_index[RLOC] = line.rstrip().split('\t')[1].split(',')

#for i in RLOCs_index:
#    print i,RLOCs_index[i]
RLOCs_file.close()

# Indexing the ENSCAFGs index (from TD intersection between BROAD and ensembl transcript files)
ENSCAFGs_TD_file = open('/home/genouest/genouest/mbahin/Annotations/ENSCAFGs_index.txt','r')
ENSCAFGs_TD = {}
for line in ENSCAFGs_TD_file:
    enscafg = line.split('\t')[0]
    ENSCAFGs_TD[enscafg] = line.split('\t')[1].rstrip()

#for i in ENSCAFGs_TD:
#    print i,ENSCAFGs_TD[i]
ENSCAFGs_TD_file.close()

# Indexing the paralogues file (from BioMart)
paral_file = open('/home/genouest/genouest/mbahin/Fusion_genes/CRAC/paralogous_140523.txt','r')
#paral_file = open('/home/genouest/umr6061/recomgen/dog/mbahin/Fusion-genes/CRAC/chimCT/paralogous_test.txt','r')
ENSCAFGs_BM = {}
paral_file.readline()
for line in paral_file:
    enscafg,gene_name,paralog = line.split('\t')[0:3]
    paralogy_type = line.split('\t')[3].rstrip()
    if gene_name == '':
        gene_name = 'No_name'
    if not ENSCAFGs_BM.has_key(enscafg):
        ENSCAFGs_BM[enscafg] = {}
        ENSCAFGs_BM[enscafg]['gene_name'] = gene_name
        ENSCAFGs_BM[enscafg]['paralogy_type'] = []
        ENSCAFGs_BM[enscafg]['paralogous_genes'] = []
    ENSCAFGs_BM[enscafg]['paralogous_genes'].append(paralog)
    ENSCAFGs_BM[enscafg]['paralogy_type'].append(paralogy_type)

#for i in ENSCAFGs_BM:
#    print i,ENSCAFGs_BM[i]
paral_file.close()

# Processing chimCT output
input = open(sys.argv[1],'r')
chim = {}
noFeat = []
for line in input:
    if line.startswith('#'):
        continue
    # Catching line information
    (chimID,chr1,end1,strand1,chr2,start2,strand2) = (line.rstrip().split('\t')[1:8])
    (spR,spP,chimClass,warnings) = (line.rstrip().split('\t')[9:13])
    if chimClass == '4':
        continue
    match=re.match(r'FusionDistance=([^,]*)',warnings)
    if match:
        warn = match.group(1)
    else:
        warn = 'No_warning'
    # Chimeras with 2 unmatched features
    RLOCs = chimID.split('---')
    if chimID == 'N/A---N/A':
        noFeat.append((spR,spP,warn,chr1,end1,strand1,chr2,start2,strand2))
        continue
    # Commands executed only the first time the chimera is met
    if not chim.has_key(chimID):
        chim[chimID] = {}
        chim[chimID]['spR'] = 0
        chim[chimID]['spP'] = 0
        chim[chimID]['breakpoint'] = {}
        chim[chimID]['warnings'] = []
    # Commands executed for each processed line
    chim[chimID]['spR'] += int(spR)
    chim[chimID]['spP'] += int(spP)
    chim[chimID]['warnings'].append(warn)
    if not chim[chimID]['breakpoint'].has_key(chimClass):
        chim[chimID]['breakpoint'][chimClass] = []
    chim[chimID]['breakpoint'][chimClass].append((chr1,end1,strand1,chr2,start2,strand2,spR,spP))

input.close()

# Opening output files
noFeat_file = open('file.noFeat.chimera','w')
oneFeat_file = open('file.oneFeat.chimera','w')
monoFeat_file = open('file.monoFeat.chimera','w')
twoFeat_file = open('file.twoFeat.chimera','w')

# Writing the file for chimera with 2 unmatched parts
noFeat.sort(key = lambda x:int(x[0]),reverse=True)
for c in noFeat:
    noFeat_file.write('\t'.join(c[0:3])+'\t'+c[3]+':'+c[4]+','+c[5]+'\t'+c[6]+':'+c[7]+','+c[8]+'\n')

orderedChim = OrderedDict(sorted(chim.iteritems(), key = lambda x:x[1]['spR'], reverse=True))
#for i in orderedChim:
#    print i,orderedChim[i]

for chim in orderedChim:
    RLOCs = chim.split('---')
    # Writing the file for chimera with 2 matched features
    if ('N/A' not in RLOCs) and (RLOCs[0] != RLOCs[1]):
        write_res(twoFeat_file, 'Both')
        # Determining if the two features are considered paralogous
        paralogous = False
        if RLOCs_index.has_key(RLOCs[0]) and RLOCs_index.has_key(RLOCs[1]):
            for enscafg1 in RLOCs_index[RLOCs[0]]:
                for enscafg2 in RLOCs_index[RLOCs[1]]:
                    if ENSCAFGs_BM.has_key(enscafg1):
                        if enscafg2 in ENSCAFGs_BM[enscafg1]['paralogous_genes']:
                            paralogous = True
                            break
                if paralogous:
                    break
        if paralogous:
            twoFeat_file.write('Paralogous\n')
        else:
            twoFeat_file.write('No_paralogy\n')
    # Writing the file for chimera with the same feature matched twice
    elif RLOCs[0] == RLOCs[1]:
        write_res(monoFeat_file, 'First')
    # Writing the files for chimera with only one matched feature
    elif RLOCs[1] == 'N/A':
        write_res(oneFeat_file, 'First')
    else:
        write_res(oneFeat_file, 'Second')

# Closing output files
noFeat_file.close()
oneFeat_file.close()
monoFeat_file.close()
twoFeat_file.close()