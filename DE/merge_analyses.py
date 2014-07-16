#!/local/python/2.7/bin/python

# Mathieu Bahin, 04/04/14
# Last update: 16/07/14

import argparse

# Getting options back
parser = argparse.ArgumentParser()
parser.add_argument('-t', dest='threshold')
parser.add_argument('-d', dest='DESeq2_file')
parser.add_argument('-e', dest='edgeR_file')
options = parser.parse_args()

index = {}

# Indexing the DESeq2 results
DESeq2_file = open(options.DESeq2_file,'r')
DESeq2_file.readline()
for line in DESeq2_file:
    xloc = line.split(',')[0].split('"')[1]
    if not index.has_key(xloc):
        index[xloc] = {}
    index[xloc]['DESeq2_score'] = line.rstrip().split(',')[-1]
    if line.split(',')[2] == 'NA':
        index[xloc]['DESeq2_fc'] = 'NA'
    elif float(line.split(',')[2]) >= 0:
        index[xloc]['DESeq2_fc'] = 'up'
    else:
        index[xloc]['DESeq2_fc'] = 'down'
DESeq2_file.close()

# Indexing the edgeR results
edgeR_file = open(options.edgeR_file,'r')
edgeR_file.readline()
for line in edgeR_file:
    xloc = line.split(',')[0].split('"')[1]
    if not index.has_key(xloc):
        index[xloc] = {}
    index[xloc]['edgeR_score'] = line.rstrip().split(',')[-1]
    if float(line.split(',')[1]) >= 0:
        index[xloc]['edgeR_fc'] = 'up'
    else:
        index[xloc]['edgeR_fc'] = 'down'
edgeR_file.close()

# Producing the scores and venn files
scores_file = open('file.scores.csv','w')
venn_file = open('file.venn.csv','w')
gene_id_file = open('file.gene_id.list','w')
scores_file.write('XLOC,DESeq2,edgeR\n')
venn_file.write('XLOC,DESeq2,edgeR\n')
threshold = float(options.threshold)

DE_genes = 0
upreg = 0
downreg = 0
ambiguous = 0
cancer_known = 0
#biotypes

for xloc in index:
    scores_file.write(xloc)
    venn_file.write(xloc)
    if index[xloc].has_key('DESeq2_score') and (index[xloc]['DESeq2_score'] != 'NA') and (index[xloc]['DESeq2_fc'] != 'NA'):
        scores_file.write(','+index[xloc]['DESeq2_score'])
        venn_file.write(','+str(int(float(index[xloc]['DESeq2_score']) <= threshold)))
    else:
        scores_file.write(',NA')
        venn_file.write(',0')
    if index[xloc].has_key('edgeR_score'):
        scores_file.write(','+index[xloc]['edgeR_score']+'\n')
        venn_file.write(','+str(int(float(index[xloc]['edgeR_score']) <= threshold))+'\n')
        # Producing the DE gene list
        if (float(index[xloc]['edgeR_score']) <= threshold) and (index[xloc].has_key('DESeq2_score')) and (index[xloc]['DESeq2_score'] != 'NA') and (float(index[xloc]['DESeq2_score']) <= threshold):
            DE_genes += 1
            if (index[xloc]['DESeq2_fc'] == index[xloc]['edgeR_fc']) and (index[xloc]['DESeq2_fc'] == 'up'):
                upreg += 1
                gene_id_file.write(xloc+'\tUp-regulated\n')
            elif (index[xloc]['DESeq2_fc'] == index[xloc]['edgeR_fc']) and (index[xloc]['DESeq2_fc'] == 'down'):
                downreg += 1
                gene_id_file.write(xloc+'\tDown-regulated\n')
            else:
                ambiguous += 1
                gene_id_file.write(xloc+'\tambiguously regulated\n')
    else:
        scores_file.write(',NA\n')
        venn_file.write(',0\n')
scores_file.close()
venn_file.close()

# Printing some statistics
print 'DE gene(s) found:', DE_genes
print '\tUp-regulated gene(s):', upreg
print '\tDown-regulated gene(s):', downreg
print '\tAmbiguously regulated gene(s):', ambiguous