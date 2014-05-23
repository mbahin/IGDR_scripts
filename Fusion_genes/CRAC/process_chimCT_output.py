#!/local/python/2.7/bin/python

# Mathieu Bahin, 21/05/14

# Script to process chimCT output

import sys, re
from collections import OrderedDict

# Indexing the RLOCs index
RLOCs_file = open ('/home/genouest/genouest/mbahin/Fusion_genes/CRAC/RLOCs_index.txt','r')
RLOCs_index = {}
for line in RLOCs_file:
    RLOC = line.split('\t')[0]
    RLOCs_index[RLOC] = line.rstrip().split('\t')[1:]

#for i in RLOCs_index:
#    print i,RLOCs_index[i]
RLOCs_file.close()

# Indexing the paralogues file
paral_file = open('/home/genouest/genouest/mbahin/Fusion_genes/CRAC/paralogous_140523.txt','r')
#paral_file = open('/home/genouest/umr6061/recomgen/dog/mbahin/Fusion-genes/CRAC/chimCT/paralogous_test.txt','r')
ENSCAFGs = {}
paral_file.readline()
for line in paral_file:
    #print line
    #print line.split('\t')
    #enscafg,gene_name,paralog,paralogy_type = line.rstrip().split('\t')[0:4]
    enscafg = line.rstrip().split('\t')[0]
    gene_name = line.rstrip().split('\t')[1]
    if gene_name == '':
        gene_name = 'No_gene_name'
    if not ENSCAFGs.has_key(enscafg):
        ENSCAFGs[enscafg] = {}
        ENSCAFGs[enscafg]['gene_name'] = gene_name
        ENSCAFGs[enscafg]['paralogy_type'] = []
        ENSCAFGs[enscafg]['paralogous_genes'] = []
    ENSCAFGs[enscafg]['paralogous_genes'].append(paralog)
    ENSCAFGs[enscafg]['paralogy_type'].append(paralogy_type)

#for i in ENSCAFGs:
#    print i,ENSCAFGs[i]
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
    # Chimeras with at least one matched feature
    	# Command executed only the first time the chimera is met
    if not chim.has_key(chimID):
        chim[chimID] = {}
        chim[chimID]['spR'] = 0
        chim[chimID]['spP'] = 0
        chim[chimID]['breakpoint'] = {}
        chim[chimID]['warnings'] = []
        chim[chimID]['paralogous'] = 'NA'
        if RLOCs[0] != 'N/A':
            if RLOCs_index.has_key(RLOCs[0]):
                chim[chimID]['ENSCAFG1'] = RLOCs_index[RLOCs[0]]
                if not chim[chimID].has_key('gene_name1'):
                    chim[chimID]['gene_name1'] = []
                for enscafg in chim[chimID]['ENSCAFG1']:
                    if ENSCAFGs.has_key(enscafg) and ENSCAFGs[enscafg]['gene_name'] != '':
                        chim[chimID]['gene_name1'].append(ENSCAFGs[enscafg]['gene_name'])
                    else:
                        chim[chimID]['gene_name1'].append('No_name')
            else:
                chim[chimID]['ENSCAFG1'] = 'No_feat'
                chim[chimID]['gene_name1'] = 'NA'
        if (RLOCs[1] != 'N/A') and (RLOCs[0] != RLOCs[1]):
            if RLOCs_index.has_key(RLOCs[1]):
                chim[chimID]['ENSCAFG2'] = RLOCs_index[RLOCs[1]]
                if not chim[chimID].has_key('gene_name2'):
                    chim[chimID]['gene_name2'] = []
                for enscafg in chim[chimID]['ENSCAFG2']:
                    if ENSCAFGs.has_key(enscafg) and ENSCAFGs[enscafg]['gene_name'] != '':
                        chim[chimID]['gene_name2'].append(ENSCAFGs[enscafg]['gene_name'])
                    else:
                        chim[chimID]['gene_name2'].append('No_name')
                    if (RLOCs[0] != 'N/A'):
                        chim[chimID]['paralogous'] = 'No_paralogy'
                        for enscafg1 in chim[chimID]['ENSCAFG1']:
                            if (ENSCAFGs.has_key(enscafg1)) and (enscafg in ENSCAFGs[enscafg1]['paralogous_genes']):
                                chim[chimID]['paralogous'] = 'Paralogous'   
            else:
                chim[chimID]['ENSCAFG2'] = 'No_feat'
                chim[chimID]['gene_name2'] = 'NA'
        # Command executed for each processed line
    chim[chimID]['spR'] += int(spR)
    chim[chimID]['spP'] += int(spP)
    chim[chimID]['warnings'].append(warn)
    if not chim[chimID]['breakpoint'].has_key(chimClass):
        chim[chimID]['breakpoint'][chimClass] = []
    chim[chimID]['breakpoint'][chimClass].append((chr1,end1,strand1,chr2,start2,strand2,spR,spP))

input.close()

# Sorting data and writing output files
noFeat_file = open('file.noFeat.chimera','w')
oneFeat_file = open('file.oneFeat.chimera','w')
monoFeat_file = open('file.monoFeat.chimera','w')
twoFeat_file = open('file.twoFeat.chimera','w')

noFeat.sort(key = lambda x:int(x[0]),reverse=True)
for c in noFeat:
    noFeat_file.write('\t'.join(c[0:3])+'\t'+c[3]+':'+c[4]+','+c[5]+'\t'+c[6]+':'+c[7]+','+c[8]+'\n')

orderedChim = OrderedDict(sorted(chim.iteritems(), key = lambda x:x[1]['spR'], reverse=True))
#for i in orderedChim:
#    print i,orderedChim[i]
for c in orderedChim:
    RLOCs = c.split('---')
    # Chimera with 2 matched features
    if (RLOCs[0] != 'N/A') and (RLOCs[1] != 'N/A') and (RLOCs[0] != RLOCs[1]):
        twoFeat_file.write(str(c)+'\t'+str(orderedChim[c]['spR'])+'\t'+str(orderedChim[c]['spP'])+'\t'+str(','.join(orderedChim[c]['warnings'])+'\t'))
        for chimClass in orderedChim[c]['breakpoint']:
            for pos in orderedChim[c]['breakpoint'][chimClass]:
                twoFeat_file.write(' ('+str(pos[0])+':'+str(pos[1])+','+str(pos[2])+' / '+str(pos[3])+':'+str(pos[4])+','+str(pos[5])+' # '+str(pos[6])+','+str(pos[7])+')')
        twoFeat_file.write('\t'+str(orderedChim[c]['ENSCAFG1'])+'\t'+str(orderedChim[c]['ENSCAFG2'])+'\t'+str(orderedChim[c]['gene_name1'])+'\t'+str(orderedChim[c]['gene_name2'])+'\t'+str(orderedChim[c]['paralogous'])+'\n')
    # Chimera with the same feature matched twice
    elif (RLOCs[0] == RLOCs[1]):
        monoFeat_file.write(str(c)+'\t'+str(orderedChim[c]['spR'])+'\t'+str(orderedChim[c]['spP'])+'\t'+str(','.join(orderedChim[c]['warnings'])+'\t'))
        for chimClass in orderedChim[c]['breakpoint']:
            for pos in orderedChim[c]['breakpoint'][chimClass]:
                monoFeat_file.write(' ('+str(pos[0])+':'+str(pos[1])+','+str(pos[2])+' / '+str(pos[3])+':'+str(pos[4])+','+str(pos[5])+' # '+str(pos[6])+','+str(pos[7])+')')
        monoFeat_file.write('\t'+str(orderedChim[c]['ENSCAFG1'])+'\t'+str(orderedChim[c]['gene_name1'])+'\n')
    # Chimera with only 1 matched feature
    else:
        oneFeat_file.write(str(c)+'\t'+str(orderedChim[c]['spR'])+'\t'+str(orderedChim[c]['spP'])+'\t'+str(','.join(orderedChim[c]['warnings'])+'\t'))
        for chimClass in orderedChim[c]['breakpoint']:
            for pos in orderedChim[c]['breakpoint'][chimClass]:
                oneFeat_file.write(' ('+str(pos[0])+':'+str(pos[1])+','+str(pos[2])+' / '+str(pos[3])+':'+str(pos[4])+','+str(pos[5])+' # '+str(pos[6])+','+str(pos[7])+')')
        if RLOCs[1] == 'N/A':
            oneFeat_file.write('\t'+str(orderedChim[c]['ENSCAFG1'])+'\t'+str(orderedChim[c]['gene_name1'])+'\n')
        else:
            oneFeat_file.write('\t'+str(orderedChim[c]['ENSCAFG2'])+'\t'+str(orderedChim[c]['gene_name2'])+'\n')
        
noFeat_file.close()
oneFeat_file.close()
monoFeat_file.close()
twoFeat_file.close()