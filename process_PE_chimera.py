#!/usr/bin/env python

# Mathieu Bahin, 06/03/14

# Script to process the chimeras detected by CRAC. There are 3 steps:
# 1) Parsing CRAC chimera output files into BED files.
# 2) Intersecting the BED files with the reference GTF (bedtools)
# 3) Parsing the bedtools output into 4 files summarizing different kind of chimera detected
# Input is a chimera output file from CRAC.
# The outputs are 4 chimera files detailing and summarizing them.

import os, re, sys
from subprocess import *
"""
NEWLINE PROBLEM
def merge_positions(pos_list):
    # Merges the overlapping positions from a list. Outputs a compact list

    # Formatting the positions list
    bedpe_file = open('file.bedpe','w')
    i = 1
    for pos in pos_list:
        details = re.match(r"\('chr([0-9X]*):([0-9]*)-([0-9]*),([+-])', 'chr([0-9X]*):([0-9]*)-([0-9]*),([+-])'\)",str(pos))
        #print details.group(0)
        #print '\t'.join((details.group(1),details.group(2),details.group(3),details.group(5),details.group(6),details.group(7),'p'+str(i),'.',details.group(4),details.group(8)))
        bedpe_file.write('Pos:\n')
        bedpe_file.write('\t'.join((details.group(1),details.group(2),details.group(3),details.group(5),details.group(6),details.group(7),'p'+str(i),'.',details.group(4),details.group(8))))#+'\n')
        i += 1
    bedpe_file.close()
    bedpe_file = open('file.bedpe','r')
    print 'BEDPE file beginning'
    for line in bedpe_file:
        print line
    print 'BEDPE file ending'
    bedpe_file.close()

    # Executing pairToPair
    p2p = open('file.pairToPair.output','w')
    for line in Popen('pairToPair -a file.bedpe -b file.bedpe', shell=True, stdout=PIPE).stdout:
        print 'L:',line
        p2p.write(line)
    p2p.close()
    
    # Processing the pairToPair results

    p2p = open('file.pairToPair.output','r')
    for line in p2p:
        print line
    p2p.close()

    return True
"""
# Defining GTF file
"""
#gtf='/omaha-beach/tderrien/DATA/canFam3/annotation/MERGED_BROAD_RENNES/merged_asm_2014-01-24/merged.withStrand.gtf'
gtf='/home/genouest/genouest/mbahin/merged.withStrand.gtf'

# Creating the output directory
if not os.path.exists('Paired-end_results'):
    os.mkdir('Paired-end_results')
    os.chdir('Paired-end_results')
else:
    sys.exit('Directory "Paired-end_chimera" already exists ! Aborting.')
"""
"""
# Opening files
chim_CRAC = open('../pe.chimera.output','r')
chim_BED = open('chimera.bed','w')

# Parsing CRAC output to get a BED file
for line in chim_CRAC:
    if not line.startswith('#'):
        pos1 = line.split(' ')[2]
        chr1 = pos1.split('|')[0]
        start1 = pos1.split('|')[1].split(',')[1]
        end1 = int(start1) + len(line.split(' ')[4]) - 1
        strand1 = pos1.split('|')[1].split(',')[0]
        if pos1.split('|')[1].split(',')[0] == '1':
            strand1 = '+'
        else:
            strand1 = '-'
        pos2 = line.split(' ')[8]
        chr2 = pos2.split('|')[0]
        start2 = pos2.split('|')[1].split(',')[1]
        end2 = int(start2) + len(line.split(' ')[10]) - 1
        strand2 = pos2.split('|')[1].split(',')[0]
        if pos2.split('|')[1].split(',')[0] == '1':
            strand2 = '+'
        else:
            strand2 = '-'
        if (pos1 != pos2):
            chim_BED.write('\t'.join((chr2,start2,str(end2),'.','.',strand2))+'\n'+'\t'.join((chr1,start1,str(end1),'.','.',strand1))+'\n')

# Closing files
chim_CRAC.close()
chim_BED.close()
"""
"""
# Intersecting BED and GTF files
intersect = open('file.intersectBed.output','w')
for line in Popen('intersectBed -a chimera.bed -b '+gtf+' -wao -s', shell=True, stdout=PIPE).stdout:
    intersect.write(line)
intersect.close()
#os.remove('chimera.bed')
"""

# Parsing the intersection results
inter = open('/home/mbahin/work/file.intersectBed.output','r')

# Parsing the intersectBed output file
# WARNING: if a first mate follows an exact same second mate just before, the script will be dephased (it is supposed to process each mate of a pair together).
index = {}
line = inter.readline().rstrip()
while True:
    # Initialiazing the mates feat matching lists
    feat = [[],[]]
    reverse = False # The 2 features order is reversed
    # Getting information from the pair's first line
    (chr1,start1,end1),strand1 = (line.split('\t')[0:3],line.split('\t')[5])
    #if line.split('\t')[-1] != '0':
        #feat[0].append(re.match(r'"(.*)";',line.split('\t')[14].split(' ')[1]).group(1))
    #else:
        #feat[0].append('No_match')
    if line.split('\t')[-1] == '0':
        feat[0].append('No_match')
    elif line.split('\t')[14].split(' ')[6] == 'gene_name':
        feat[0].append(re.match(r'"(.*)";',line.split('\t')[14].split(' ')[7]).group(1))
    else:
        #feat[0].append('No_ensembl')
        feat[0].append(re.match(r'"(.*)";',line.split('\t')[14].split(' ')[1]).group(1))
    line2 = inter.readline().rstrip()
    
    # Looping until the pair's second line (several lines to skip for a same gene)
    while (line2.split('\t')[0:3],line2.split('\t')[5]) == ([chr1,start1,end1],strand1):
        #line_feat = re.match(r'"(.*)";',line2.split('\t')[14].split(' ')[1]).group(1)
        if line2.split('\t')[14].split(' ')[6] == 'gene_name':
            line_feat = re.match(r'"(.*)";',line2.split('\t')[14].split(' ')[7]).group(1)
        else:
            #line_feat = 'No_ensembl'
            line_feat = re.match(r'"(.*)";',line2.split('\t')[14].split(' ')[1]).group(1)
        if line_feat not in feat[0]:
            feat[0].append(line_feat)
        line2 = inter.readline().rstrip()
    
    # Getting information from the pair's second line
    (chr2,start2,end2),strand2 = (line2.split('\t')[0:3],line2.split('\t')[5])
    #if line2.split('\t')[-1] != '0':
        #feat[1].append(re.match(r'"(.*)";',line2.split('\t')[14].split(' ')[1]).group(1))
    #else:
        #feat[1].append('No_match')
    if line2.split('\t')[-1] == '0':
        feat[1].append('No_match')
    elif line2.split('\t')[14].split(' ')[6] == 'gene_name':
        feat[1].append(re.match(r'"(.*)";',line2.split('\t')[14].split(' ')[7]).group(1))
    else:
        #feat[1].append('No_ensembl')
        feat[1].append(re.match(r'"(.*)";',line2.split('\t')[14].split(' ')[1]).group(1))
    
    # Looping until the next pair's first line (several lines to skip for a same gene)
    line = inter.readline().rstrip()
    if line:
        while (line.split('\t')[0] == chr2) and (line.split('\t')[1] == start2) and (line.split('\t')[2] == end2) and (line.split('\t')[5] == strand2):
            #line_feat = re.match(r'"(.*)";',line.split('\t')[14].split(' ')[1]).group(1)
            #if line_feat not in feat[1]:
                #feat[1].append(line_feat)
            if line.split('\t')[14].split(' ')[6] == 'gene_name':
                line_feat = re.match(r'"(.*)";',line.split('\t')[14].split(' ')[7]).group(1)
            else:
                #line_feat = 'No_ensembl'
                line_feat = re.match(r'"(.*)";',line.split('\t')[14].split(' ')[1]).group(1)
            line = inter.readline().rstrip()
        # Filling the index (chimID = "XLOC1%XLOC2@XLOC3%XLOC4")
        if feat[0] == ['No_match']:
            reverse = True
            chimID = '%'.join(sorted(feat[1]))+'@No_match'
        elif feat[1] == ['No_match']:
            chimID = '%'.join(sorted(feat[0]))+'@No_match'
        elif sorted(feat[0])[0] <= sorted(feat[1])[0]:
            chimID = '@'.join(('%'.join(sorted(feat[0])),'%'.join(sorted(feat[1]))))
        else:
            reverse = True
            chimID = '@'.join(('%'.join(sorted(feat[1])),'%'.join(sorted(feat[0]))))
        """
        if feat[0] == ['No_ensembl']:
            reverse = True
            chimID = '%'.join(sorted(feat[1]))+'@No_ensembl'
        elif feat[1] == ['No_ensembl']:
            chimID = '%'.join(sorted(feat[0]))+'@No_ensembl'
        elif sorted(feat[0])[0] <= sorted(feat[1])[0]:
            chimID = '@'.join(('%'.join(sorted(feat[0])),'%'.join(sorted(feat[1]))))
        else:
            reverse = True
            chimID = '@'.join(('%'.join(sorted(feat[1])),'%'.join(sorted(feat[0]))))
        #chimID = '@'.join(('%'.join(sorted(feat[0])),'%'.join(sorted(feat[1]))))
        """
        loc1 = 'chr'+chr1+':'+start1+'-'+end1+','+strand1
        loc2 = 'chr'+chr2+':'+start2+'-'+end2+','+strand2
        if index.has_key(chimID):
            index[chimID]['count'] += 1
            if not reverse:
                index[chimID]['pos'].append((loc1,loc2))
            else:
                index[chimID]['pos'].append((loc2,loc1))
        else:
            index[chimID] = {}
            index[chimID]['count'] = 1
            if not reverse:
                index[chimID]['pos'] = [(loc1,loc2)]
            else:
                index[chimID]['pos'] = [(loc2,loc1)]
        if not line:
            break
    else:
        break

#for i in index:
    #print i,index[i]

# Closing file
inter.close()
#"""
# Indexing the paralogous genes file
paralogous = open('/home/mbahin/work/paralogous.txt','r')

ENSCAFGs = {}
for line in paralogous:
    enscafg,gene_name,paralog,paralogy_type = line.rstrip().split(',')[0:4]
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
    #print i,ENSCAFGs[i]

paralogous.close()
#"""
#"""
# Writing the results
noMatch_file = open('/home/mbahin/work/Res_script/file.noMatch.chimera','w')
output = open('/home/mbahin/work/Res_script/output.chimera','w')
#oneMatch_file = open('/home/mbahin/work/Res_script/file.oneMatch.chimera','w')
#sameXLOC_file = open('/home/mbahin/work/Res_script/file.sameXLOC.chimera','w')
#twoXLOC_file = open('/home/mbahin/work/Res_script/file.twoXLOC.chimera','w')

for chimID in index:
    if (chimID == 'No_match@No_match'):
        noMatch_file.write(str(index[chimID]['count'])+'\n')
        for pos in index['No_match@No_match']['pos']:
            noMatch_file.write(str(pos)+'\n')
    else:
        if len(index[chimID]['pos']) > 1:
            merge_positions(index[chimID]['pos'])
        output.write('\t'.join((str(index[chimID]['count']),chimID.split('@')[0],chimID.split('@')[1],str(index[chimID]['pos']))))
        for feat in chimID.split('@')[0].split('%'):
            if not feat.startswith('XLOC'):
                output.write('\t'+feat+'\t'+ENSCAFGs[feat]['gene_name'])
                if len(set(ENSCAFGs[feat]['paralogy_type'])) == 1:
                    output.write('\t'+ENSCAFGs[feat]['paralogy_type'][0])
                else:
                    output.write('\t'+str(ENSCAFGs[feat]['paralogy_type']))
                output.write('\t'+str(ENSCAFGs[feat]['paralogous_genes']))
        if (chimID.split('@')[1] != 'No_match') and (chimID.split('@')[0] != chimID.split('@')[1]):
            for feat in chimID.split('@')[1].split('%'):
                if not feat.startswith('XLOC'):
                    output.write('\t'+feat+'\t'+ENSCAFGs[feat]['gene_name'])
                    if len(set(ENSCAFGs[feat]['paralogy_type'])) == 1:
                        output.write('\t'+ENSCAFGs[feat]['paralogy_type'][0])
                    else:
                        output.write('\t'+str(ENSCAFGs[feat]['paralogy_type']))
                    output.write('\t'+str(ENSCAFGs[feat]['paralogous_genes']))
        output.write('\n')

# Closing files
noMatch_file.close()
output.close()
#oneMatch_file.close()
#sameXLOC_file.close()
#twoXLOC_file.close()
#"""
"""
# Writing the results
noMatch_file = open('file.noMatch.chimera','w')
oneMatch_file = open('file.oneMatch.chimera','w')
sameXLOC_file = open('file.sameXLOC.chimera','w')
twoXLOC_file = open('file.twoXLOC.chimera','w')

for i in index:
    if (i == 'No_match@No_match'):
        #noMatch_file.write(str(index[i]['count'])+'\t'+i+'\t'+str(index[i]['pos'])+'\n')
        noMatch_file.write(str(index[i]['count'])+'\n')
        for pos in index[i]['pos']:
            noMatch_file.write(str(pos)+'\n')
    elif i.split('@')[1] == 'No_match':
        oneMatch_file.write(str(index[i]['count'])+'\t'+i+'\t'+str(index[i]['pos']))
        for xloc in i.split('@')[0].split('%'):
            oneMatch_file.write('\t'+xloc+' ->')
            #print xloc
            if XLOCs[xloc] != ['No_ensembl_annotation']:
                for enscafg in XLOCs[xloc]:
                    oneMatch_file.write(' '+enscafg)
                oneMatch_file.write('\t')
                for enscafg in XLOCs[xloc]:
                    oneMatch_file.write(' '+ENSCAFGs[enscafg]['gene_name'])
        oneMatch_file.write('\n')
    elif (i.split('@')[0] == i.split('@')[1]):
        sameXLOC_file.write(str(index[i]['count'])+'\t'+i.split('@')[0]+'\t'+str(index[i]['pos']))
        for xloc in i.split('@')[0].split('%'):
            sameXLOC_file.write('\t'+xloc+' ->')
            if XLOCs[xloc] != ['No_ensembl_annotation']: 
                for enscafg in XLOCs[xloc]:
                    sameXLOC_file.write(' '+enscafg)
                sameXLOC_file.write('\t')
                for enscafg in XLOCs[xloc]:
                    sameXLOC_file.write(' '+ENSCAFGs[enscafg]['gene_name'])
                for enscafg in XLOCs[xloc]:
                    sameXLOC_file.write('\t'+enscafg+' -> '+str(ENSCAFGs[enscafg]['paralogy_type']))
                for enscafg in XLOCs[xloc]:
                    sameXLOC_file.write('\t'+enscafg+' -> '+str(ENSCAFGs[enscafg]['paralogous_genes']))
        sameXLOC_file.write('\n')
    else:
        z = 1

noMatch_file.close()
oneMatch_file.close()
sameXLOC_file.close()
twoXLOC_file.close()
"""
"""
# Indexing the GTF file
gtf = open('/home/mbahin/work/merged.withStrand.gtf','r')

XLOCs = {}
for line in gtf:
    details = line.split('\t')[8].split(' ')
    xloc = details[1].split('"')[1]
    if not XLOCs.has_key(xloc):
        XLOCs[xloc] = ['No_ensembl_annotation']
    if details[6] == 'gene_name':
        if XLOCs[xloc] == ['No_ensembl_annotation']:
            XLOCs[xloc] = [details[7].split('"')[1]]
        elif not details[7].split('"')[1] in XLOCs[xloc]:
            XLOCs[xloc].append(details[7].split('"')[1])        

gtf.close()
"""
print 'Coucou final'
# To add:
	# - Check the two features distances ?
	# - Merging the read pairs (positions)

# Closing files
#os.remove('file.intersectBed.output')